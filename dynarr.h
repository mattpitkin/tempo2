#ifndef DYNARR_H
#define DYNARR_H

#include <stdlib.h>

/* redwards dynamic array code, to help relieve C++ withdrawl symptoms */

typedef struct
{
  /* public */
  void *data;
  size_t nelem;

  /* private */
  size_t elem_size;
  size_t nalloced;
} DynamicArray;

void DynamicArray_init(DynamicArray *, size_t elemSize);
void DynamicArray_resize(DynamicArray *, size_t nelem);
void *DynamicArray_push_back(DynamicArray *, void *elem);
void DynamicArray_free(DynamicArray *);

#endif
