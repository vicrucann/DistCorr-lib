#include <stdio.h>
#include <math.h>
#include "misc.h"
#include "ntuple.h"
#include "ntuple_ll.h"

/*----------------------------------------------------------------------------*/
/** Add a 2-tuple to an 2-tuple list.
 */
static void add_2tuple( ntuple_list out, double v1, double v2 )
{
  /* check parameters */
  if( out == NULL ) error("add_2tuple: invalid n-tuple input.");
  if( out->dim != 2 ) error("add_2tuple: the n-tuple must be a 2-tuple.");

  /* if needed, alloc more tuples to 'out' */
  if( out->size == out->max_size ) enlarge_ntuple_list(out);
  if( out->values == NULL ) error("add_2tuple: invalid n-tuple input.");

  /* add new 2-tuple */
  out->values[ out->size * out->dim + 0 ] = v1;
  out->values[ out->size * out->dim + 1 ] = v2;

  /* update number of tuples counter */
  out->size++;
}

/*----------------------------------------------------------------------------*/
/** Add a 1-tuple to an 1-tuple list.
 */
static void add_1tuple( ntuple_list out, double v1)
{
  /* check parameters */
  if( out == NULL ) error("add_1tuple: invalid n-tuple input.");
  if( out->dim != 1 ) error("add_1tuple: the n-tuple must be a 1-tuple.");

  /* if needed, alloc more tuples to 'out' */
  if( out->size == out->max_size ) enlarge_ntuple_list(out);
  if( out->values == NULL ) error("add_2tuple: invalid n-tuple input.");

  /* add new 2-tuple */
  out->values[ out->size * out->dim] = v1;

  /* update number of tuples counter */
  out->size++;
}

/*----------------------------------------------------------------------------*/
/** Copy ntuple_list
 */
static void copy_ntuple(ntuple_list src, ntuple_list des)
{
  /* check parameters */
  if( src == NULL ) error("copy_ntuple: invalid n-tuple input.");
  if( des == NULL ) error("copy_ntuple: invalid n-tuple input.");
  if( src->dim != des->dim ) error("copy_ntuple: different dimension of input."); 

  int i, j, nb, dim;
  nb = src->size;
  dim = src->dim;
  
  des->size = 0;

  for(i = 0; i < nb; i++) {
    if(des->size == des->max_size) 
      enlarge_ntuple_list(des);
    for(j = 0; j < dim; j++)       
      des->values[i*dim+j] = src->values[i*dim+j];
    des->size++;
  }
}


/*--------------------------------------------------------------------------*/
ntuple_list gaussian_convol_on_curve(double unit_sigma, double Nsigma, bool resampling, bool eliminate_border, double up_factor, double down_factor, ntuple_list orig_edges)
{
  double x1, y1, x2, y2, sumx, sumy, d, norm, w;
  int j, k;
  
  double line_length, sample_step = 0.0;
  ntuple_list upsampled_edges, convolved_upsampled_edges, each_seg_length, output_pts;
  each_seg_length = new_ntuple_list(1);

  double residual_d = 0.0, lambda;
  int low_bound, high_bound;
  int size_lambda, begin_lambda;
  double small_seg_x, small_seg_y;
  double sigma;
  int flag, nb_sigma, not_border_flag;
  double real_down_factor;
  unsigned int nb_pts;

  nb_pts = orig_edges->size;

  upsampled_edges = new_ntuple_list(2);

  /*length of the whole line*/
  line_length = 0.0;

  /*the total length of the line and the lenth of each segment*/
  for(j = 0; j < (int)nb_pts-1; j++) {
    x1 = orig_edges->values[j*orig_edges->dim + 0];
    y1 = orig_edges->values[j*orig_edges->dim + 1];
    x2 = orig_edges->values[(j+1)*orig_edges->dim + 0];
    y2 = orig_edges->values[(j+1)*orig_edges->dim + 1];
    d = dist(x1, y1, x2, y2);
    
    add_1tuple(each_seg_length, d);
    line_length += d;
  }

  /*the average distance between two points for a line (original line)*/
  sample_step = line_length / (nb_pts-1);


  if (resampling) {
    /*TANG: the edge points are not regularly sampled, so first do an interpolation to get regularly sampled points*/
    if (up_factor < 1.0)
      error("up_factor must be greater than 1!\n");

    /*the line is upsampled, the average distance becomes smaller*/
    sample_step /= up_factor;

    /* Resampling */
    for(j = 0; j < (int)nb_pts-1; j++) {
      flag = 0;

      x1 = orig_edges->values[j*orig_edges->dim + 0];
      y1 = orig_edges->values[j*orig_edges->dim + 1];
      x2 = orig_edges->values[(j+1)*orig_edges->dim + 0];
      y2 = orig_edges->values[(j+1)*orig_edges->dim + 1];

      if (j == 0)
	residual_d = 0.0;

      size_lambda = (int)floor((each_seg_length->values[j]+residual_d)/sample_step);
      
      begin_lambda = 0;
      if (j != 0)
	begin_lambda = 1;

      for(k = begin_lambda; k <= size_lambda; k++) {
	lambda = (-residual_d + k*sample_step) / each_seg_length->values[j];
	small_seg_x = (1-lambda)*x1 + lambda*x2;
	small_seg_y = (1-lambda)*y1 + lambda*y2;

	add_2tuple(upsampled_edges, small_seg_x, small_seg_y);
	
	flag = 1;
      }
      
      if(flag == 1)
	residual_d = dist(x2, y2, small_seg_x, small_seg_y);
      else
	residual_d += each_seg_length->values[j];
    }

  }
  else /* no resampling, copy 'orig_edges' to 'upsampled_edges' */
    copy_ntuple(orig_edges, upsampled_edges);
  

  /* convolution */
  convolved_upsampled_edges = new_ntuple_list(2);

  for(j=0; j<(int)upsampled_edges->size; j++)
  {
    norm = 0.0; sumx = 0.0; sumy = 0.0;

    sigma = unit_sigma * sqrt(down_factor*down_factor - 1.0);
    nb_sigma = (int)floor( Nsigma * sigma / sample_step );

    low_bound = j - nb_sigma;
    if(low_bound < 0)
      low_bound = 0;
    high_bound = j + nb_sigma;
    if(high_bound > (int)upsampled_edges->size-1)
      high_bound = upsampled_edges->size-1;

    if(eliminate_border) {
      if(high_bound - low_bound == 2*nb_sigma) 
	not_border_flag = 1;
      else
	not_border_flag = 0;
    }
    else
      not_border_flag = 1;

    if(not_border_flag == 1) {
      for (k = low_bound; k <= high_bound; k++) {
	/*x1 = upsampled_edges->list[i]->values[upsampled_edges->list[i]->dim*j];
	y1 = upsampled_edges->list[i]->values[upsampled_edges->list[i]->dim*j+1];
	x2 = upsampled_edges->list[i]->values[upsampled_edges->list[i]->dim*k];
	y2 = upsampled_edges->list[i]->values[upsampled_edges->list[i]->dim*k+1];
	d = dist(x1, y1, x2, y2);*/
	x2 = upsampled_edges->values[upsampled_edges->dim*k];
	y2 = upsampled_edges->values[upsampled_edges->dim*k+1];
	d = sample_step * (k-j);
	w = exp( -d*d / 2.0 / sigma / sigma );
	sumx += w*x2;
	sumy += w*y2;
	norm += w;
      }
  
      sumx /= norm;
      sumy /= norm;

      add_2tuple(convolved_upsampled_edges, sumx, sumy);
    }
  }

  /* if resampling is enabled */
  /*
  output_pts = new_ntuple_list(2);
  if(resampling)
  {
    real_down_factor = (int)(down_factor * up_factor);
    
    nb = convolved_upsampled_edges->size / real_down_factor + 1;
    for(j = 0; j < nb; j++)
      add_2tuple(output_pts, convolved_upsampled_edges->values[j*real_down_factor*convolved_upsampled_edges->dim], convolved_upsampled_edges->values[j*real_down_factor*convolved_upsampled_edges->dim+1]);
  }
  else 
    copy_ntuple(upsampled_edges, output_pts);
  */

  output_pts = new_ntuple_list(2);
  if(resampling)
    real_down_factor = down_factor * up_factor;
  else
    real_down_factor = down_factor;
  
  /*
  nb = (int)((double)convolved_upsampled_edges->size / real_down_factor);
  
  for(j = 0; j < nb; j++)
    add_2tuple(output_pts, convolved_upsampled_edges->values[(int)(j*real_down_factor*convolved_upsampled_edges->dim)], convolved_upsampled_edges->values[(int)(j*real_down_factor*convolved_upsampled_edges->dim)+1]);
  */

  j = 0;
  while(j < (int)convolved_upsampled_edges->size) 
  { 
    add_2tuple(output_pts, convolved_upsampled_edges->values[(int)(j*convolved_upsampled_edges->dim)], convolved_upsampled_edges->values[(int)(j*convolved_upsampled_edges->dim)+1]);
    j = floor(j + real_down_factor);
  }

  free_ntuple_list(upsampled_edges);
  free_ntuple_list(convolved_upsampled_edges);
  free_ntuple_list(each_seg_length);

  return output_pts;
}


