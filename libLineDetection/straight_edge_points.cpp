/*----------------------------------------------------------------------------

  Algorithm that produces a list of straight edges; for each one,
  a list of edge points is produced.

  Copyright 2011 rafael grompone von gioi (grompone@gmail.com),
                 Zhongwei Tang (tang@cmla.ens-cachan.fr).

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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "misc.h"
#include "image.h"
#include "ntuple.h"
#include "ntuple_ll.h"
#include "devernay.h"
#include "lsd.h"

#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif

/*----------------------------------------------------------------------------*/
struct point
{
  float x,y;
};

/*----------------------------------------------------------------------------*/
static int comp_point_x(const void * a, const void * b)
{
  if( ((struct point *)a)->x <  ((struct point *)b)->x ) return -1;
  /* if( ((struct point *)a)->x == ((struct point *)b)->x ) return  0; */
  if( ((struct point *)a)->x >  ((struct point *)b)->x ) return  1;
  return 0; /* so c compiler knows the functions does return a value */
}

/*----------------------------------------------------------------------------*/
static int comp_point_y(const void * a, const void * b)
{
  if( ((struct point *)a)->y <  ((struct point *)b)->y ) return -1;
  /* if( ((struct point *)a)->y == ((struct point *)b)->y ) return  0; */
  if( ((struct point *)a)->y >  ((struct point *)b)->y ) return  1;
  return 0; /* so c compiler knows the functions does return a value */
}

static int comp_point_y_decrease(const void *a, const void *b)
{
  if( ((struct point *)a)->y <  ((struct point *)b)->y ) return 1;
  /* if( ((struct point *)a)->y == ((struct point *)b)->y ) return 0; */
  if( ((struct point *)a)->y >  ((struct point *)b)->y ) return -1;
  return 0; /* so c compiler knows the functions does return a value */
}

/*----------------------------------------------------------------------------*/
static void sort_edge_points(ntuple_ll se, ntuple_list ls)
{
  unsigned int i,j,n;
  struct point * p;
  double theta;

  for(i=0; i<se->size; i++)
    {
      n = se->list[i]->size;

      /* get memory */
      p = (struct point *) calloc( (size_t) n, sizeof(struct point) );
      if( p == NULL ) error("not enough memory!");

      /* copy points to list */
      for(j=0;j<n;j++)
        {
          p[j].x = se->list[i]->values[ j*2 + 0 ];
          p[j].y = se->list[i]->values[ j*2 + 1 ];
        }

      /* sort */
      /* if line segment angle < PI/4, sort according to x value,
         else according to y value. */
      /*
      if( fabs(atan( ( ls->values[ i * ls->dim + 3 ]
                     - ls->values[ i * ls->dim + 1 ] ) /
                     ( ls->values[ i * ls->dim + 2 ]
                     - ls->values[ i * ls->dim + 0 ] ) )) < (M_PI / 4.0) )
        qsort( (void *)p, n, sizeof(struct point), &comp_point_x );
      else
        qsort( (void *)p, n, sizeof(struct point), &comp_point_y );
      */ 

      theta = atan( ( ls->values[ i * ls->dim + 3 ]
	      - ls->values[ i * ls->dim + 1 ] ) /
	    ( ls->values[ i * ls->dim + 2 ]
	      - ls->values[ i * ls->dim + 0 ] ) );
      if( theta < (M_PI / 4.0)  && theta > (-M_PI / 4.0))
        qsort( (void *)p, n, sizeof(struct point), &comp_point_x );
      else if (theta >= (M_PI / 4.0))
        qsort( (void *)p, n, sizeof(struct point), &comp_point_y );
      else if (theta <= -(M_PI / 4.0))
        qsort( (void *)p, n, sizeof(struct point), &comp_point_y_decrease );
      else
	error("[straight_edge_points]: theta error!\n");


      /* copy back sorted points */
      for(j=0;j<n;j++)
        {
          se->list[i]->values[ j*2 + 0 ] = p[j].x;
          se->list[i]->values[ j*2 + 1 ] = p[j].y;
        }

      /* free memory */
      free( (void *) p );
    }
}

/*----------------------------------------------------------------------------*/
static void remove_short_edges(ntuple_ll se, ntuple_list ls, double threshold)
{
  unsigned int i;
  double length;

  /* mark short edges to be removed
    and remove its list of points */
  for(i=0; i<se->size; i++)
    {
      length = dist( ls->values[i*ls->dim + 0], ls->values[i*ls->dim + 1],
                     ls->values[i*ls->dim + 2], ls->values[i*ls->dim + 3] );

      if( length < threshold )
        {
          free_ntuple_list(se->list[i]);  /* remove list of points */
          se->list[i] = NULL;             /* mark the edge to be removed */
        }
    }

  /* remove marked edges */
  for(i=0; i<se->size; i++)
    if( se->list[i] == NULL )
      {
        se->list[i] = se->list[se->size-1];
        se->size--;
        i--;
      }
}

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
ntuple_ll straight_edge_points( image_double image, double sigma,
                                double th_low, double th_hi, double min_length )
{
  ntuple_list edges;   /* edge points */
  ntuple_list ls;      /* line segments */
  image_int region; 
  ntuple_ll se;        /* straight edges */
  unsigned int i,x,y,seg;
  double xx,yy;

  /* call Devernay sub-pixel edge detector */
  edges = devernay(image,sigma,th_low,th_hi);

  /* call LSD line segment detector */
  ls = lsd_scale_region( image, 1.0, &region );

  /* create list of straight edges and list of its points */
  se = new_ntuple_ll(ls->size);
  for(i=0; i<ls->size; i++)
    add_ntuple_list(se,new_ntuple_list(2)); /* create all the n-tuple lists */ //<- memory leak

  /* assign edge points to line segments */
  for(i=0; i<edges->size; i++)
    {
      xx = edges->values[ i * edges->dim + 0 ];
      yy = edges->values[ i * edges->dim + 1 ];
      x = xx; /* interger part of the edge coordinates */
      y = yy;

      if( x<region->xsize && y<region->ysize ) /* unsigned, so x>=0 and y>=0 */
        if( ( seg = region->data[ x + y * region->xsize ] ) > 0 )
          if( (seg-1) < se->size ) /* seg corresponds to line segment seg-1 */
            add_2tuple(se->list[seg-1],xx,yy);
    }

  /* sort edge points in each line segment */
  sort_edge_points(se,ls);

  /* remove edges correspoinding to short line segments */
  remove_short_edges(se,ls,min_length);

  /* free memory */
  free_image_int(region);
  free_ntuple_list(edges);
  free_ntuple_list(ls);

  return se;
}
/*----------------------------------------------------------------------------*/
