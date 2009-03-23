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
	       int *newpar,int *onlypre,char *dcmFile)
{
  int i;
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
	  else if (strcmp(argv[i],"-gr")==0)
	    gr=1;
	  else if (strcmp(argv[i],"-dcm")==0)
	    strcpy(dcmFile,argv[++i]);
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
	      printf("Have %s\n",timFile[timfile_num]);
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
