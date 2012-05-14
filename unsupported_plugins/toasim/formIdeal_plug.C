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

void copyObservation(observation *obs2,observation *obs1);

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,nit,j,p;
  char fname[MAX_FILELEN];
  double globalParameter;
  FILE *fout;
  int giveRef=0;
  long double refMJD;
  double refFreq;
  char refSite[128];

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
      else if (strcmp(argv[i],"-refmjd")==0)
	{
	  sscanf(argv[++i],"%Lf",&refMJD);
	  giveRef=1;
	}
      else if (strcmp(argv[i],"-refsite")==0)
	strcpy(refSite,argv[++i]);
      else if (strcmp(argv[i],"-reffreq")==0)
	sscanf(argv[++i],"%lf",&refFreq);
    }


  // Now read in all the .par and .tim files
  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  
  //
  // We need to put the reference MJD as the first observation
  //
  if (giveRef == 1)
    {
      // Add 1 to each observation number and then add in the reference line
      for (p=0;p<*npsr;p++)
	{
	  for (i=psr[p].nobs-1;i>=0;i--)
	    copyObservation(&(psr[p].obsn[i+1]),&(psr[p].obsn[i]));
	  psr[p].nobs++;
	  strcpy(psr[p].obsn[0].fname,"reference");
	  strcpy(psr[p].obsn[0].telID,refSite);
	  psr[p].obsn[0].sat = refMJD;
	  psr[p].obsn[0].freq = refFreq;
	  psr[p].obsn[0].toaErr = psr[p].obsn[0].origErr = 1e-4;
	  psr[p].obsn[0].nFlags = 0;
	  psr[p].obsn[0].clockCorr = 1;
	  psr[p].obsn[0].delayCorr = 1;
	  psr[p].obsn[0].deleted = 0;
	  psr[p].obsn[0].phaseOffset = 0;
	}
      writeTim("tryit.tim",&psr[0],"tempo2");
    }

  preProcess(psr,*npsr,argc,argv);

  printf("'%s' NOBS %d\n",parFile[*npsr-1],psr[0].nobs);
  for (i=0;i<nit;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,0);     /* Form the residuals                 */
      // Remove the residuals from the site arrival times
      for (p=0;p<*npsr;p++)
	{
	  for (j=0;j<psr[p].nobs;j++){
	    psr[p].obsn[j].sat -= (psr[p].obsn[j].residual/SECDAY);	    
	    //printf("Residual hello = %g (%d)\n",(double)psr[p].obsn[j].residual,j);

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

void copyObservation(observation *obs2,observation *obs1)
{
  obs2->sat = obs1->sat;
  obs2->clockCorr = obs1->clockCorr;
  obs2->delayCorr = obs1->delayCorr;
  obs2->deleted = obs1->deleted;
  obs2->freq = obs1->freq;
  obs2->toaErr = obs1->toaErr;
  obs2->origErr = obs1->origErr;
  obs2->phaseOffset = obs1->phaseOffset;
  strcpy(obs2->fname,obs1->fname);
  strcpy(obs2->telID,obs1->telID);
  obs2->nFlags = obs1->nFlags;
  for (int i=0;i<obs1->nFlags;i++)
    {
      strcpy(obs2->flagID[i],obs1->flagID[i]);
      strcpy(obs2->flagVal[i],obs1->flagVal[i]);
    }
  obs2->efac = obs1->efac;
  obs2->equad = obs1->equad;
}

char * plugVersionCheck = TEMPO2_h_VER;
