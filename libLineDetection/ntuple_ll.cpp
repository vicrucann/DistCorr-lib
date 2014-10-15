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
/** @file ntuple_ll.c
    'list of list of n-tuple' data type.
    @author rafael grompone von gioi (grompone@gmail.com)
 */
/*----------------------------------------------------------------------------*/
#include <stdlib.h>
#include "misc.h"
#include "ntuple_ll.h"

/*----------------------------------------------------------------------------*/
/** Free memory used in n-tuple list list 'in'.
 */
void free_ntuple_ll(ntuple_ll in)
{
  unsigned int i;

  /* check parameters */
  if( in == NULL || in->list == NULL )
    error("free_ntuple_ll: invalid n-tuple list list input.");

  /* free elements in the list */
  for(i=0; i<in->size; i++)
    if( in->list[i] != NULL )
      free_ntuple_list(in->list[i]);

  /* free memory */
  free( (void *) in->list );
  free( (void *) in );
}

/*----------------------------------------------------------------------------*/
/** Create an n-tuple list list and allocate memory for 'size' n-tuple list.
    @param size number of n-tuples list to allocate.
 */
ntuple_ll new_ntuple_ll(unsigned int size)
{
  ntuple_ll ll;

  /* check parameters */
  if( size == 0 ) error("new_ntuple_ll: initial size must be positive.");

  /* get memory for list structure */
  ll = (ntuple_ll) malloc( sizeof(struct ntuple_ll_s) );
  if( ll == NULL ) error("not enough memory.");

  /* initialize list */
  ll->size = 0;
  ll->max_size = size;

  /* get memory for n-tuple list list */
  ll->list = (ntuple_list *) malloc( ll->max_size * sizeof(ntuple_list) );
  if( ll->list == NULL ) error("not enough memory.");

  return ll;
}

/*----------------------------------------------------------------------------*/
/** Enlarge the allocated memory of an n-tuple list list.
 */
void enlarge_ntuple_ll(ntuple_ll ll)
{
  /* check parameters */
  if( ll == NULL || ll->list == NULL || ll->max_size == 0 )
    error("enlarge_ntuple_ll: invalid n-tuple list list.");

  /* duplicate number of tuples list size */
  ll->max_size *= 2;

  /* realloc memory */
  ll->list = (ntuple_list *) realloc( (void *) ll->list,
                                      ll->max_size * sizeof(ntuple_list) );
  if( ll->list == NULL ) error("not enough memory.");
}

/*----------------------------------------------------------------------------*/
/** Add an n-tuple list to the end of an n-tuple list list.
 */
void add_ntuple_list(ntuple_ll ll, ntuple_list l)
{
  /* check parameters */
  if( ll == NULL || ll->list == NULL )
    error("add_ntuple_list: invalid n-tuple list list input.");

  /* if needed, alloc enlarge tuples list list */
  if( ll->size == ll->max_size ) enlarge_ntuple_ll(ll);

  /* add new n-tuple list */
  ll->list[ll->size] = l;
  ll->size++;
}

/*----------------------------------------------------------------------------*/
