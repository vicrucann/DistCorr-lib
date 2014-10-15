/*----------------------------------------------------------------------------

  n-tuple list list data type and basic functions.

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
/*------------------- 'list of list of n-tuple' data type --------------------*/
/*----------------------------------------------------------------------------*/
/** @file ntuple_ll.h
    'list of list of n-tuple' data type.
    @author rafael grompone von gioi (grompone@gmail.com)
 */
/*----------------------------------------------------------------------------*/
#ifndef NTUPLE_LL_HEADER
#define NTUPLE_LL_HEADER

#include "ntuple.h"

/*----------------------------------------------------------------------------*/
/** 'list of n-tuple' data type

    The n-tuple list number i is accessed with:

      ntll->list[i]

    The number of number of n-tuples lists in the list is:

      ntl->size

    The maximum number of n-tuples lists that can be stored in the
    list with the allocated memory at a given time is given by:

      ntl->max_size
 */
typedef struct ntuple_ll_s
{
  unsigned int size;
  unsigned int max_size;
  ntuple_list * list;
} * ntuple_ll;

/*----------------------------------------------------------------------------*/
/** Free memory used in n-tuple list list 'in'.
 */
void free_ntuple_ll(ntuple_ll in);

/*----------------------------------------------------------------------------*/
/** Create an n-tuple list list and allocate memory for 'size' n-tuple list.
    @param size number of n-tuples list to allocate.
 */
ntuple_ll new_ntuple_ll(unsigned int size);

/*----------------------------------------------------------------------------*/
/** Enlarge the allocated memory of an n-tuple list list.
 */
void enlarge_ntuple_ll(ntuple_ll ll);


/*----------------------------------------------------------------------------*/
/** Add an n-tuple list to an n-tuple list list.
 */
void add_ntuple_list(ntuple_ll ll, ntuple_list l);

#endif /* !NTUPLE_LL_HEADER */
/*----------------------------------------------------------------------------*/
