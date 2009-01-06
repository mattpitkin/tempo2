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
