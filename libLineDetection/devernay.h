/*----------------------------------------------------------------------------

  Implementation of Devernay's sub-pixel edge detector.

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
/*
   Implementation of the algorithm described in
   "A Non-Maxima Suppression Method for Edge Detection with Sub-Pixel Accuracy"
   by Frederic Devernay. Rapport de recherche INRIA No.2724, November 1995.
*/
/*----------------------------------------------------------------------------*/

#ifndef DEVERNAY_HEADER
#define DEVERNAY_HEADER

#include "image.h"
#include "ntuple.h"

/*----------------------------------------------------------------------------*/
/** Devernay sub-pixel edge detector.
 */
ntuple_list devernay( image_double image, double sigma,
                      double th_low, double th_hi );

#endif /* !DEVERNAY_HEADER */
/*----------------------------------------------------------------------------*/
