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

#include "tabulatedfunction.h"
#include "tempo2.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

void
TabulatedFunction_load(TabulatedFunction *func,
			     char *fileName)
{
  FILE *f;
  char line[1024]="", *c;
  TabulatedFunctionSample sample;
  int narg;

  DynamicArray_init(&func->samples, sizeof(TabulatedFunctionSample));

/*   sprintf(fname, "%s/clock2/%s.clk", getenv("TEMPO2"), fileName); */
  f = fopen(fileName, "r");
  if (f==NULL)
  {
    fprintf(stderr, "Fatal Error: Unable to open file %s for reading: %s\n",
	      fileName, strerror(errno));
    exit(1);
  }

  fgets(func->header_line, 1024, f); // first line is a header line
  if (func->header_line[0]!='#')
  {
    fprintf(stderr, 
	    "Error parsing file %s: first line must begin with #\n",
	    fileName);
    exit(1);
  }
  for (c=line; *c=='#'; c++)
    ;
    
  /* parse table lines */
  while (fgets(line, 1024, f)!=NULL)
  {
    if (line[0]!='#')
    {
      if (sscanf(line, "%lf %lf", &sample.x, &sample.y)==2)
	DynamicArray_push_back(&func->samples, &sample);
    }
  }

  fclose(f); 

  /* Make shortened filename and store it */
  /* find last / */
  for (c = fileName + strlen(fileName)-1; *c!='/' && c != fileName; c--)
    ;
  if (*c=='/')
    c++;
  strcpy(func->fileName, c); /* copy after / */
  //  func->fileName[strlen(func->fileName)-4] = '\0';/* wipe off .clk */
}

double
TabulatedFunction_getValue(TabulatedFunction *func,
				      double x)
{
  TabulatedFunctionSample *samp = (TabulatedFunctionSample *)func->samples.data;
  size_t isamp;

  /* check for out of bounds conditions */
  if (samp[0].x  > x)
  {
    char msg[1000],msg2[1000];
    sprintf(msg,"requested value before available data! (%s @",
	    func->fileName);
    sprintf(msg2,"%.1f)",x);
    displayMsg(1,"TAB1",msg,msg2,1);
    return samp[0].y;
  }
  if (samp[func->samples.nelem-1].x  < x)
  {
    char msg[1000],msg2[1000];
    sprintf(msg,"requested value after available data! (%s @",
	    func->fileName);
    sprintf(msg2,"%.1f)",x);
    displayMsg(1,"TAB2",msg,msg2,1);
    return samp[func->samples.nelem-1].y;
  }

  /* find first sample to fall before requested time */
#if 0
  for (isamp = 0; 
       isamp < func->samples.nelem  && samp[isamp].x <= x; 
       isamp++)
    ;  
#else  /* binary search code, aka the "The Price Is Right" method */
  int imin=0, imax=func->samples.nelem-1;
  int imid;
  do
  {
    imid = (imin+imax)/2;
    if (samp[imid].x > x)
      imax = imid;
    else
      imin = imid;
  } while (imax > imin+1);
  isamp = imax;
#endif

  /* interpolate */
  return (x - samp[isamp-1].x) / (samp[isamp].x - samp[isamp-1].x)
    * (samp[isamp].y - samp[isamp-1].y) + samp[isamp-1].y;
}

double
TabulatedFunction_getStartX(TabulatedFunction *func)
{
  return ((TabulatedFunctionSample *)func->samples.data)[0].x;
}
double
TabulatedFunction_getEndX(TabulatedFunction *func)
{
  return ((TabulatedFunctionSample *)func->samples.data)[func->samples.nelem-1].x;
}
