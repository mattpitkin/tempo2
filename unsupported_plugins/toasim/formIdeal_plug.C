//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards
//
// Plugin to form idealised site-arrival-times


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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"

using namespace std;



extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,nit,j,p;
  char fname[MAX_FILELEN];
  double globalParameter;
  FILE *fout;

  *npsr = 0;
  nit = 4;

  printf("Graphical Interface: formIdeal\n");
  printf("Author:              G. Hobbs, M. Keith\n");
  printf("Version:             1.0\n");

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[*npsr],argv[i+1]); 
	  strcpy(timFile[*npsr],argv[i+2]);
	  (*npsr)++;
	}
    }


  // Now read in all the .par and .tim files
  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

  printf("'%s' NOBS %d\n",parFile[*npsr-1],psr[0].nobs);
  for (i=0;i<nit;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);     /* Form the residuals                 */
      // Remove the residuals from the site arrival times
      for (p=0;p<*npsr;p++)
	{
	  for (j=0;j<psr[p].nobs;j++){
	    psr[p].obsn[j].sat -= (psr[p].obsn[j].residual/SECDAY);
	  }
	}
    }

  // Now ouptut new .tim files
  for (p=0;p<*npsr;p++)
    {
      sprintf(fname,"%s.sim",timFile[p]);
      writeTim(fname,&psr[p],"tempo2");
      
      // Write more information about the file
      fout = fopen(fname,"a");
      fprintf(fout,"# Filename: %s\n",fname);
      fprintf(fout,"# Original .tim file: %s\n",timFile[p]);
      fprintf(fout,"# Original .par file: %s\n",parFile[p]);
      fclose(fout);
    }

  return 0;
}

