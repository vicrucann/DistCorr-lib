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
/** @file gauss.h
    Image Gaussian filtering and subsampling.
    @author rafael grompone von gioi (grompone@gmail.com)
 */
/*----------------------------------------------------------------------------*/
#ifndef GAUSS_HEADER
#define GAUSS_HEADER

#include "image.h"
#include "ntuple.h"

void gaussian_kernel(ntuple_list kernel, double sigma, double mean);
void gaussian_filter(image_double image, double sigma);
image_double gaussian_sampler(image_double in,double scale, double sigma_scale);

#endif /* !GAUSS_HEADER */
/*----------------------------------------------------------------------------*/
