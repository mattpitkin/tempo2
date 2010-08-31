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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"

/* ******************************************** */
/* getInputs                                    */
/* Author:  G. Hobbs (02 May 2003)              */
/* Purpose: Parse the command line              */
/* Inputs:  argc, argv - command line arguments */
/* Outputs: File path for .tim and .par files   */
/*                                              */
/* Notes:                                       */
/*                                              */
/* Changes:                                     */
/* ******************************************** */

void getInputs(pulsar *psr,int argc, char *argv[],char timFile[][MAX_FILELEN],
	       char parFile[][MAX_FILELEN],int *list,int *npsr,
	       int *nGlobal,int *outRes,int *writeModel,char *outputSO,
	       int *polyco, char *polyco_args,
	       int *newpar,int *onlypre,char *dcmFile,char *covarFuncFile)
{
  int i,p;
  int gr=0;
  int parfile_num = 0;  /* Have we got a parfile? */
  int timfile_num = 0;  /* Have we got a timfile? */
  int gotTim=0;
  *list = 0;  /* Don't list parameters */
  //  *nGlobal=0; /* How many global parameters are we fitting? */

  if (argc==2) /* Just have .tim file name */
    {
      if (strcmp(argv[1],"-h")==0) /* Some help */
	{
	  printf("\n\n");
	  printf("tempo2 v1\n\n");
	  printf("examples: tempo2 mytim.tim\n");
	  printf("          tempo2 -f mypar.par mytim.tim\n");
	  printf("          tempo2 -gr plk -f mypar.par mytim.tim\n");
	  printf("\n");
	  printf("Options: \n\n");

	  printf("-epoch centre     Centres the PEPOCH in the fit\n");
	  printf("-f parFile        Selects parameter file\n");
	  printf("-dcm dcmFile      Data covariance matrix file\n");
	  printf("-gr name          Uses 'name' plugin for graphical interface\n");
	  printf("-h                This help\n");
	  printf("-list             Provides listing of clock corrections and residuals etc.\n");
	  printf("-output name      Uses 'name' plugin for output format\n");
	  printf("-pred \"args\"    Creates a predictive 2D Chebyshev polynomial.\n");
	  printf("      args = \"sitename mjd1 mjd2 freq1 freq2 ntimecoeff nfreqcoeff seg_length (s)\"\n"); 
	  printf("-polyco \"args\"  Creates a TEMPO1-style polyco file.\n");
	  printf("                   args = \"mjd1 mjd2 nspan ncoeff maxha sitename freq\"\n");
	  printf("-residuals        Outputs the residuals\n");
	  printf("-allInfo          Prints out clock, Earth orientation and similar information\n");
	  printf("-reminder         Saves the command line to T2command.input for future reference.\n");
	  printf("-norescale        Do not rescale parameter uncertainties by the sqrt(red. chisq)\n");
	  printf("\n\n");
	  printf("Available plugins\n");
	  system("ls $TEMPO2/plugins/ | grep plug | sed s/\"_\"/\" \"/ | awk '{print \"  - \" $1}' | sort | uniq");
	  printf("-----------------\n");
	  exit(1);
	}	    
      strcpy(timFile[timfile_num],argv[1]);
      strcpy(parFile[parfile_num],argv[1]);
      strcpy(parFile[parfile_num]+strlen(parFile[parfile_num])-3,"par");
      timfile_num++; parfile_num++;
      (*npsr)++;
    }
  else /* Have multiple command line arguments */
    {
      for (i=1;i<argc;i++)
	{
	  if (strcmp(argv[i],"-f")==0) /* Have .par file */
	    {
	      strcpy(parFile[parfile_num],argv[++i]);
	      if (argv[i+1][0]!='-') /* Must be tim file */
		{
		  strcpy(timFile[timfile_num],argv[++i]);
		  timfile_num++;
		  gotTim=1;
		}
	      /*	      i+=2; */
	      parfile_num++;
	    }
	  else if (strcmp(argv[i],"-fitfunc")==0)
	    strcpy(psr[0].fitFunc,argv[++i]);
	  else if (strcasecmp(argv[i],"-norescale")==0)
	    {
	      for (p=0;p<MAX_PSR;p++)
		psr[p].rescaleErrChisq=0;
	    }
	  else if (strcmp(argv[i],"-gr")==0)
	    gr=1;
	  else if (strcmp(argv[i],"-dcm")==0)
	    strcpy(dcmFile,argv[++i]);
	  else if (strcmp(argv[i],"-dcf")==0)
	    strcpy(covarFuncFile,argv[++i]);
	  else if (strcmp(argv[i],"-filter")==0)
	    {
	      strcat(psr[0].filterStr,argv[++i]);
	      strcat(psr[0].filterStr," ");
	    }
	  else if (strcmp(argv[i],"-pass")==0)
	    {
	      strcat(psr[0].passStr,argv[++i]);
	      strcat(psr[0].passStr," ");
	    }
	  else if (strcmp(argv[i],"-list")==0)
	    *list = 1;	  
	  else if (strcmp(argv[i],"-output")==0)
	      strcpy(outputSO,argv[++i]);
	  else if (strcmp(argv[i],"-residuals")==0)
	    *outRes = 1;
	  else if (strcmp(argv[i],"-del")==0)
	    strcpy(psr[0].deleteFileName,argv[++i]);	  
	  else if (strcmp(argv[i],"-model")==0)
	    *writeModel = 1;
	  else if (strcmp(argv[i],"-newpar")==0)
	    *newpar=1;
	  else if (strcmp(argv[i],"-pre")==0) /* Don't iterate to form post-fit residuals */
	    *onlypre = 1;
	  else if (strcmp(argv[i],"-pred")==0)
	    {
	      *polyco=2;
	      strcpy(polyco_args, argv[++i]);
	    }
	  else if (strcmp(argv[i],"-polyco")==0)
	    {
	      *polyco=1;
	      strcpy(polyco_args, argv[++i]);
	    }
	  else if (i==argc-1 && gotTim==0) /* Must be .tim file name */
	    {
	      strcpy(timFile[timfile_num],argv[i]);	    
	      //printf("Have %s\n",timFile[timfile_num]);
	      timfile_num++;
	    }
	}
      if (timfile_num==0 && *polyco==0 && gr==0)
	{
	  printf("ERROR [FILE1]: No .tim file given on command line\n"); 
	  exit(1);
	}      
      if (parfile_num==0 && gr==0)
	{
	  printf("ERROR [FILE2]: No .par file given on command line\n"); 
	  exit(1);
	}
      *npsr = parfile_num;
    }
}
