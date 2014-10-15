/*----------------------------------------------------------------------------

  Image data type and basic functions.

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
/*----------------------------- Image Data Types -----------------------------*/
/*----------------------------------------------------------------------------*/
/** @file image.h
    Image data types.
    @author rafael grompone von gioi (grompone@gmail.com)
 */
/*----------------------------------------------------------------------------*/
#ifndef IMAGE_HEADER
#define IMAGE_HEADER

/*----------------------------------------------------------------------------*/
/** char image data type

    The pixel value at (x,y) is accessed by:

      image->data[ x + y * image->xsize ]

    with x and y integer.
 */
typedef struct image_char_s
{
  unsigned char * data;
  unsigned int xsize,ysize;
} * image_char;

void free_image_char(image_char i);
image_char new_image_char(unsigned int xsize, unsigned int ysize);
image_char new_image_char_ini( unsigned int xsize, unsigned int ysize,
                               unsigned char fill_value );
image_char new_image_char_copy(image_char in);

/*----------------------------------------------------------------------------*/
/** int image data type

    The pixel value at (x,y) is accessed by:

      image->data[ x + y * image->xsize ]

    with x and y integer.
 */
typedef struct image_int_s
{
  int * data;
  unsigned int xsize,ysize;
} * image_int;

void free_image_int(image_int i);
image_int new_image_int(unsigned int xsize, unsigned int ysize);
image_int new_image_int_ini( unsigned int xsize, unsigned int ysize,
                             int fill_value );
image_int new_image_int_copy(image_int in);

/*----------------------------------------------------------------------------*/
/** double image data type

    The pixel value at (x,y) is accessed by:

      image->data[ x + y * image->xsize ]

    with x and y integer.
 */
typedef struct image_double_s
{
  double * data;
  unsigned int xsize,ysize;
} * image_double;

void free_image_double(image_double i);
image_double new_image_double(unsigned int xsize, unsigned int ysize);
image_double new_image_double_ini( unsigned int xsize, unsigned int ysize,
                                   double fill_value );
image_double new_image_double_copy(image_double in);
double interpolate_image_double(image_double& in, int order, double u, double v);

bool valid_image_double(image_double& in, int x, int y);

#endif /* !IMAGE_HEADER */
/*----------------------------------------------------------------------------*/
