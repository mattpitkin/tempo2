//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards
//
// Plugin to add various offsets to an idealised .tim file


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
#include "T2toolkit.h"
#include "tempo2.h"
#include "toasim.h"

using namespace std;

#define MAX_CORR 10 // MAXIMUM NUMBER OF CORRECTIONS

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,nit,j,p;
  char fname[MAX_CORR][MAX_FILELEN];
  int  ireal[MAX_CORR];
  char fname_pertf[MAX_FILELEN];
  char nname[MAX_FILELEN];
  double globalParameter;
  long double result;
  long seed = TKsetSeed();
  toasim_header_t* read_header;
  FILE* file;
  int64_t offsets[MAX_OBSN]; // Will change to doubles - should use malloc
  double offset[MAX_OBSN];
  int ncorr=0;
  char mode=0;
  
  *npsr = 0;
  nit = 1;

  printf("Graphical Interface: createRealisation\n");
  printf("Author:              G. Hobbs, M. Keith\n");
  printf("Version:             1.0\n");

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(timFile[*npsr],argv[++i]);
	  (*npsr)++;
	}
      else if (strcmp(argv[i],"-pertf")==0)
        {
	   strcpy(fname_pertf,argv[++i]);
	   if(file = fopen(fname_pertf,"r")){
		while(!feof(file)){
			fscanf(file,"%s %d\n",fname[ncorr],ireal+ncorr);
			ncorr++;
			printf("%s %d\n",fname[ncorr],ireal[ncorr]);
		}
	   }else{
		   fprintf(stderr,"Could not open pert file '%s'\n",fname_pertf);
	   }
        }
      else if (strcmp(argv[i],"-corn")==0)
        {
	  strcpy(fname[ncorr],argv[++i]);
	  ireal[ncorr]=atoi(argv[++i]);
	  ncorr++;
	}
      else if (strcmp(argv[i],"-corr")==0)
	{
	  strcpy(fname[ncorr],argv[++i]);
	  ireal[ncorr]=0;
	  ncorr++;
	}
    }
  if (ncorr==0)
    {
      printf("Must include some correction files using -corr\n");
      exit(1);
    }


  // Now read in all the .tim files
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */

  for (p=0;p<*npsr;p++)
    {
      sprintf(nname,"%s.real",timFile[p]);
      printf("Writing file: %s\n",nname);
      for (j=0;j<ncorr;j++)
	{
	  file = fopen(fname[j],"r");
	  if(file==NULL){
		  printf("ERROR: File could not be read: '%s'\n",fname[j]);
		  exit(1);
	  }
	  read_header = toasim_read_header(file);
	  if(strcmp(read_header->timfile_name,timFile[p])!=0){
		  fprintf(stderr,"\n\n*****************\nWARNING: .tim file name mismatch '%s' != '%s'\n*****************\n\n",read_header->timfile_name,timFile[p]);
	  }
	  toasim_corrections_t *read_corr= toasim_read_corrections(read_header,ireal[j],file);
	  for (i=0;i<psr[p].nobs;i++)
	    {
	      if (j==0)
		offset[i]=(double)(read_corr->offsets[i]);
	      else
		offset[i]+=(double)(read_corr->offsets[i]);
	    }
	  fclose(file);
	}
      for (i=0;i<psr[p].nobs;i++)
	psr[p].obsn[i].sat += (long double)offset[i]/SECDAY;
      writeTim(nname,&psr[p],"tempo2");
    }
  return 0;
}

char * plugVersionCheck = TEMPO2_h_VER;
