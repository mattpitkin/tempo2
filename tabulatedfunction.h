#ifndef TABULATEDFUNCTION_H
#define TABULATEDFUNCTION_H

#include "dynarr.h"

// redwards general ASCII table handling routines for tempo2

typedef struct
{
  double x;
  double y;
} TabulatedFunctionSample;

typedef struct
{
  char fileName[256];
  char header_line[256];
  DynamicArray samples; 
} TabulatedFunction;

extern void TabulatedFunction_load(TabulatedFunction *func, char *fileName);
extern double TabulatedFunction_getValue(TabulatedFunction *func,
					  double x);
extern double TabulatedFunction_getStartX(TabulatedFunction *func);
extern double TabulatedFunction_getEndX(TabulatedFunction *func);


#endif
