#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

/*
*    This file is part of TEMPO2. 
* 
*    TEMPO2 is free software: you can redistribute it and/or modify 
*    it under the terms of the GNU General Public License as published by 
*    the Free Software Foundation, either version 3 of the License, or 
*    (at your option) any later version. 
*    TEMPO2 is distributed in the hope that it will be useful, 
*    but WITHOUT ANY WARRANTY; without even the implied warranty of 
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
*    GNU General Public License for more details. 
*    You should have received a copy of the GNU General Public License 
*    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
*/

/*
*    If you use TEMPO2 then please acknowledge it by citing 
*    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
*    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
*    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
*    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
*    timing model.
*/


/*-*-C-*- */

#include "dynarr.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


void DynamicArray_init(DynamicArray *a, size_t elemSize)
{
  a->data = NULL;
  a->nelem = a->nalloced = 0;
  a->elem_size = elemSize;
}

void DynamicArray_resize(DynamicArray *a, size_t nelem)
{
  size_t needed, increment;

  a->nelem = nelem;
  needed = a->nelem - a->nalloced;
  if (needed > 0)
  {
    /*  increment in minimum 20 % increments to save on reallocs */
    increment = (a->nelem < 5 ? 1 : a->nelem / 5);
    increment = (needed > increment ? needed : increment);
    a->nalloced += increment;
    if (a->data == NULL)
      {
	//	printf("MALLOCING\n");
	a->data = malloc(a->nalloced * a->elem_size);
      }
    else
      a->data = realloc(a->data, a->nalloced * a->elem_size);
    if (a->data == NULL)
    {
      fprintf(stderr, "Fatal Error: Out of memory in %s:%d!\n", 
	      __FILE__, __LINE__);
      exit(1);
    }
  }
}

void *DynamicArray_push_back(DynamicArray *a, void *elem)
{
  void *last;
  DynamicArray_resize(a, a->nelem+1);
  last = (char*)a->data+(a->nelem-1)*a->elem_size;
  memcpy(last, elem, a->elem_size);
  return last;
}

void DynamicArray_free(DynamicArray *a)
{
  free(a->data);
  a->data = NULL;
  a->nelem = a->nalloced = 0;
}
