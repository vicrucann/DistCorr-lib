#ifndef GAUSSIAN_CONVOL_ON_CURVE_HEADER
#define GAUSSIAN_CONVOL_ON_CURVE_HEADER

#include "ntuple.h"

/*----------------------------------------------------------------------------*/
/** Gaussian Convolution on curve
  unit_sigma: Gaussain standard deviation for 1-sampling, typical value is 0.8
  Nsigma: the truncation of Gaussian, typical value is 3.0
  resampling: flag for resampling
  eliminate_border: flag for eliminate points on border not receving a complete Gaussian convolution
  up_factor: the upsampling factor if resampling is enabled
  down_factor: the downsampling factor (with respect to the original points)
  orig_edges: original input edge points detected on a curve (typical detected and grouped by devernay detector and LSD)
 */

ntuple_list gaussian_convol_on_curve(double unit_sigma, double Nsigma, bool resampling, bool eliminate_border, double up_factor, double down_factor, ntuple_list orig_edges);

#endif
