/*----------------------------------------------------------------------------

  PGM image file format I/O functions.

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
/*------------------------------ PGM image I/O -------------------------------*/
/*----------------------------------------------------------------------------*/
/** @file pgm_io.h
    PGM image format read and write functions.
    @author rafael grompone von gioi (grompone@gmail.com)
 */
/*----------------------------------------------------------------------------*/
#ifndef PGMIO_HEADER
#define PGMIO_HEADER

#include "image.h"

/*----------------------------------------------------------------------------*/
/** Read a PGM file into an "image_char".
    If the name is "-" the file is read from standard input.
 */
image_char read_pgm_image_char(char * name);

/*----------------------------------------------------------------------------*/
/** Read a PGM file into an "image_int".
    If the name is "-" the file is read from standard input.
 */
image_int read_pgm_image_int(char * name);

/*----------------------------------------------------------------------------*/
/** Read a PGM file into an "image_double".
    If the name is "-" the file is read from standard input.
 */
image_double read_pgm_image_double(char * name);

/*----------------------------------------------------------------------------*/
/** Write an "image_char" into a PGM file.
    If the name is "-" the file is written to standard output.
 */
void write_pgm_image_char(image_char image, char * name);

/*----------------------------------------------------------------------------*/
/** Write an "image_int" into a PGM file.
    If the name is "-" the file is written to standard output.
 */
void write_pgm_image_int(image_int image, char * name);

/*----------------------------------------------------------------------------*/
/** Write an "image_int" normalized to [0,255] into a PGM file.
    If the name is "-" the file is written to standard output.
 */
void write_pgm_image_int_normalized(image_int image, char * name);

/*----------------------------------------------------------------------------*/
/** Write an "image_double" into a PGM file.
    If the name is "-" the file is written to standard output.
 */
void write_pgm_image_double(image_double image, char * name);

/*----------------------------------------------------------------------------*/
/** Write an "image_double" normalized to [0,255] into a PGM file.
    If the name is "-" the file is written to standard output.
 */
void write_pgm_image_double_normalized(image_double image, char * name);

#endif /* !PGMIO_HEADER */
/*----------------------------------------------------------------------------*/
