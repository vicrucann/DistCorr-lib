/*----------------------------------------------------------------------------

  Image Gaussian filtering functions.

  Copyright 2010-2011 rafael grompone von gioi (grompone@gmail.com)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  ----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------- Gaussian filter ------------------------------*/
/*----------------------------------------------------------------------------*/
/** @file gauss.c
    Image Gaussian filtering and subsampling.
    @author rafael grompone von gioi (grompone@gmail.com)
 */
/*----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "misc.h"
#include "image.h"
#include "ntuple.h"
#include "gauss.h"

/*----------------------------------------------------------------------------*/
/** Compute a Gaussian kernel of length 'kernel->dim',
    standard deviation 'sigma', and centered at value 'mean'.

    For example, if mean=0.5, the Gaussian will be centered
    in the middle point between values 'kernel->values[0]'
    and 'kernel->values[1]'.
 */
void gaussian_kernel(ntuple_list kernel, double sigma, double mean)
{
  double sum = 0.0;
  double val;
  unsigned int i;

  /* check parameters */
  if( kernel == NULL || kernel->values == NULL )
    error("gaussian_kernel: invalid n-tuple 'kernel'.");
  if( sigma <= 0.0 ) error("gaussian_kernel: 'sigma' must be positive.");

  /* compute Gaussian kernel */
  if( kernel->max_size < 1 ) enlarge_ntuple_list(kernel);
  kernel->size = 1;
  for(i=0;i<kernel->dim;i++)
    {
      val = ( (double) i - mean ) / sigma;
      kernel->values[i] = exp( -0.5 * val * val );
      sum += kernel->values[i];
    }

  /* normalization */
  if( sum >= 0.0 ) for(i=0;i<kernel->dim;i++) kernel->values[i] /= sum;
}

/*----------------------------------------------------------------------------*/
void gaussian_filter(image_double image, double sigma)
{
  unsigned int x,y;
  int offset,i,j,nx2,ny2,n;
  ntuple_list kernel;
  image_double tmp;
  double val,prec;

  if( sigma <= 0.0 ) error("gaussian_filter: 'sigma; must be positive.");
  if( image == NULL || image->data == NULL ||
      image->xsize <= 0 || image->ysize <= 0 )
    error("gaussian_filter: invalid image.");

  /* create temporary image */
  tmp = new_image_double(image->xsize,image->ysize);

  /* compute gaussian kernel */
  /*
     The size of the kernel is selected to guarantee that the
     the first discarted term is at least 10^prec times smaller
     than the central value. For that, h should be larger than x, with
       e^(-x^2/2sigma^2) = 1/10^prec.
     Then,
       x = sigma * sqrt( 2 * prec * ln(10) ).
   */
  prec = 3.0;
  offset = ceil( sigma * sqrt( 2.0 * prec * log(10.0) ) );
  n = 1 + 2 * offset; /* kernel size */
  kernel = new_ntuple_list(n);
  gaussian_kernel( kernel, sigma, (double) offset );

  /* auxiliary double image size variables */
  nx2 = 2*image->xsize;
  ny2 = 2*image->ysize;

  /* x axis convolution */
  for(x=0;x<image->xsize;x++)
    for(y=0;y<image->ysize;y++)
      {
        val = 0.0;
        for(i=0;i<n;i++)
          {
            j = x - offset + i;

            /* symmetry boundary condition */
            while(j<0) j+=nx2;
            while(j>=nx2) j-=nx2;
            if( j >= (int) image->xsize ) j = nx2-1-j;

            val += image->data[ j + y * image->xsize ] * kernel->values[i];
          }
        tmp->data[ x + y * tmp->xsize ] = val;
      }

  /* y axis convolution */
  for(x=0;x<image->xsize;x++)
    for(y=0;y<image->ysize;y++)
      {
        val = 0.0;
        for(i=0;i<n;i++)
          {
            j = y - offset + i;

            /* symmetry boundary condition */
            while(j<0) j+=ny2;
            while(j>=ny2) j-=ny2;
            if( j >= (int) image->ysize ) j = ny2-1-j;

            val += tmp->data[ x + j * tmp->xsize ] * kernel->values[i];
          }
        image->data[ x + y * image->xsize ] = val;
      }

  /* free memory */
  free_ntuple_list(kernel);
  free_image_double(tmp);
}

/*----------------------------------------------------------------------------*/
/** Scale the input image 'in' by a factor 'scale' by Gaussian sub-sampling.

    For example, scale=0.8 will give a result at 80% of the original size.

    The image is convolved with a Gaussian kernel
    @f[
        G(x,y) = \frac{1}{2\pi\sigma^2} e^{-\frac{x^2+y^2}{2\sigma^2}}
    @f]
    before the sub-sampling to prevent aliasing.

    The standard deviation sigma given by:
    -  sigma = sigma_scale / scale,   if scale <  1.0
    -  sigma = sigma_scale,           if scale >= 1.0

    To be able to sub-sample at non-integer steps, some interpolation
    is needed. In this implementation, the interpolation is done by
    the Gaussian kernel, so both operations (filtering and sampling)
    are done at the same time. The Gaussian kernel is computed
    centered on the coordinates of the required sample. In this way,
    when applied, it gives directly the result of convolving the image
    with the kernel and interpolated to that particular position.

    A fast algorithm is done using the separability of the Gaussian
    kernel. Applying the 2D Gaussian kernel is equivalent to applying
    first a horizontal 1D Gaussian kernel and then a vertical 1D
    Gaussian kernel (or the other way round). The reason is that
    @f[
        G(x,y) = G(x) * G(y)
    @f]
    where
    @f[
        G(x) = \frac{1}{\sqrt{2\pi}\sigma} e^{-\frac{x^2}{2\sigma^2}}.
    @f]
    The algorithm first apply a combined Gaussian kernel and sampling
    in the x axis, and then the combined Gaussian kernel and sampling
    in the y axis.
 */
image_double gaussian_sampler(image_double in, double scale, double sigma_scale)
{
  image_double aux,out;
  ntuple_list kernel;
  unsigned int N,M,h,n,x,y,i;
  int xc,yc,j,double_x_size,double_y_size;
  double sigma,xx,yy,sum,prec;

  /* check parameters */
  if( in == NULL || in->data == NULL || in->xsize == 0 || in->ysize == 0 )
    error("gaussian_sampler: invalid image.");
  if( scale <= 0.0 ) error("gaussian_sampler: 'scale' must be positive.");
  if( sigma_scale <= 0.0 )
    error("gaussian_sampler: 'sigma_scale' must be positive.");

  /* get memory for images */
  if( in->xsize * scale > (double) UINT_MAX ||
      in->ysize * scale > (double) UINT_MAX )
    error("gaussian_sampler: the output image size exceeds the handled size.");
  N = (unsigned int) floor( in->xsize * scale );
  M = (unsigned int) floor( in->ysize * scale );
  aux = new_image_double(N,in->ysize);
  out = new_image_double(N,M);

  /* sigma, kernel size and memory for the kernel */
  sigma = scale < 1.0 ? sigma_scale / scale : sigma_scale;
  /*
     The size of the kernel is selected to guarantee that the
     the first discarded term is at least 10^prec times smaller
     than the central value. For that, h should be larger than x, with
       e^(-x^2/2sigma^2) = 1/10^prec.
     Then,
       x = sigma * sqrt( 2 * prec * ln(10) ).
   */
  prec = 3.0;
  h = (unsigned int) ceil( sigma * sqrt( 2.0 * prec * log(10.0) ) );
  n = 1+2*h; /* kernel size */
  kernel = new_ntuple_list(n);

  /* auxiliary double image size variables */
  double_x_size = (int) (2 * in->xsize);
  double_y_size = (int) (2 * in->ysize);

  /* First subsampling: x axis */
  for(x=0;x<aux->xsize;x++)
    {
      /*
         x   is the coordinate in the new image.
         xx  is the corresponding x-value in the original size image.
         xc  is the integer value, the pixel coordinate of xx.
       */
      /* pixel (0,0) corresponds to coordinate (0.5,0.5)
         because the pixel goes from 0.0 to 1.0 */
      xx = ( (double) x + 0.5 ) / scale;
      xc = (int) floor( xx );
      gaussian_kernel( kernel, sigma, (double) h + xx - (double) xc - 0.5 );
      /* the kernel must be computed for each x because the fine
         offset xx-xc is different in each case */

      for(y=0;y<aux->ysize;y++)
        {
          sum = 0.0;
          for(i=0;i<kernel->dim;i++)
            {
              j = xc - h + i;

              /* symmetry boundary condition */
              while( j < 0 ) j += double_x_size;
              while( j >= double_x_size ) j -= double_x_size;
              if( j >= (int) in->xsize ) j = double_x_size-1-j;

              sum += in->data[ j + y * in->xsize ] * kernel->values[i];
            }
          aux->data[ x + y * aux->xsize ] = sum;
        }
    }

  /* Second subsampling: y axis */
  for(y=0;y<out->ysize;y++)
    {
      /*
         y   is the coordinate in the new image.
         yy  is the corresponding x-value in the original size image.
         yc  is the integer value, the pixel coordinate of xx.
       */
      /* pixel (0,0) corresponds to coordinate (0.5,0.5)
         because the pixel goes from 0.0 to 1.0 */
      yy = ( (double) y + 0.5 ) / scale;
      yc = (int) floor( yy );
      gaussian_kernel( kernel, sigma, (double) h + yy - (double) yc - 0.5 );
      /* the kernel must be computed for each y because the fine
         offset yy-yc is different in each case */

      for(x=0;x<out->xsize;x++)
        {
          sum = 0.0;
          for(i=0;i<kernel->dim;i++)
            {
              j = yc - h + i;

              /* symmetry boundary condition */
              while( j < 0 ) j += double_y_size;
              while( j >= double_y_size ) j -= double_y_size;
              if( j >= (int) in->ysize ) j = double_y_size-1-j;

              sum += aux->data[ x + j * aux->xsize ] * kernel->values[i];
            }
          out->data[ x + y * out->xsize ] = sum;
        }
    }

  /* free memory */
  free_ntuple_list(kernel);
  free_image_double(aux);

  return out;
}

/*----------------------------------------------------------------------------*/
