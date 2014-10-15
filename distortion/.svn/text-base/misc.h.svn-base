/*----------------------------------------------------------------------------

  Miscellanea functions.

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
/*-------------------------------- Miscellanea -------------------------------*/
/*----------------------------------------------------------------------------*/
/** @file misc.h
    Miscellaneous functions.
    @author rafael grompone von gioi (grompone@gmail.com)
 */
/*----------------------------------------------------------------------------*/
#ifndef MISC_HEADER
#define MISC_HEADER

#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

/*----------------------------------------------------------------------------*/
/** Fatal error, print a formated message to standard-error output and exit.

    The syntax is exactly as the one of printf, but it adds "error: "
    before the message and then it end the program with an error number.
 */
void error(const char * msg, ...);

/*----------------------------------------------------------------------------*/
/** Doubles relative error factor
 */
#ifndef RELATIVE_ERROR_FACTOR
#define RELATIVE_ERROR_FACTOR 100.0
#endif /* !RELATIVE_ERROR_FACTOR */

/*----------------------------------------------------------------------------*/
/** Compare doubles by relative error.

    The resulting rounding error after floating point computations
    depend on the specific operations done. The same number computed by
    different algorithms could present different rounding errors. For a
    useful comparison, an estimation of the relative rounding error
    should be considered and compared to a factor times EPS. The factor
    should be related to the cumulated rounding error in the chain of
    computation. Here, as a simplification, a fixed factor is used.
 */
int double_equal(double a, double b);

/*----------------------------------------------------------------------------*/
/** Computes Euclidean distance between point (x1,y1) and point (x2,y2).
 */
double dist(double x1, double y1, double x2, double y2);

#endif /* !MISC_HEADER */
/*----------------------------------------------------------------------------*/
