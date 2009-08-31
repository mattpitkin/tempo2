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

/* TEMPO2.C                                         */
/* --------                                         */
/*                                                  */
/* TEMPO2.C is a clone of TEMPO, but is written     */
/* in C/C++ and contains full descriptions of the   */
/* precision available from each routine.           */
/* Each variable is documented and can be           */
/* tabulated.                                       */
/*                                                  */
/* V1.0 G. Hobbs, R. Edwards                        */

//#include <iostream.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "tempo2.h"
#include "tempo2Util.h"
#include "tempo2pred.h"
#include "tempo2pred_int.h"
#include <dlfcn.h>

void ephemeris_routines(pulsar *psr,int npsr);
void clock_corrections(pulsar *psr,int npsr);
void extra_delays(pulsar *psr,int npsr);

int main(int argc, char *argv[])
{
  int iteration; 
  int listparms;
  int outRes=0;
  int writeModel=0;
  char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
  char outputSO[MAX_FILELEN];
  char str[MAX_FILELEN];
  longdouble coeff[MAX_COEFF]; /* For polynomial coefficients in polyco */
  char dcmFile[MAX_FILELEN]="NULL";
  int npsr;      /* The number of pulsars */
  int noWarnings=1;
  double globalParameter=0.0;
  int  displayParams,p;
  int nGlobal,i,flagPolyco=0,it,k;
  char polyco_args[128];
  int newpar=0;
  int onlypre=0;
  char tempo2MachineType[MAX_FILELEN]="";
  FILE *alias;
  char **commandLine;
  clock_t startClock,endClock;

  printf("This program comes with ABSOLUTELY NO WARRANTY.\n");
  printf("This is free software, and you are welcome to redistribute it\n");
  printf("under conditions of GPL license.\n\n");

  startClock = clock();

  pulsar *psr;

  commandLine = (char **)malloc(1000*sizeof(char *));

  for (i=0;i<1000;i++)
    commandLine[i] = (char *)malloc(sizeof(char)*1000);

  /* Parse input line for machine type */
  for (i=0;i<argc;i++)    
    {
      if (strcasecmp(argv[i],"-npsr")==0)
	sscanf(argv[i+1],"%d",&MAX_PSR);
      if (strcasecmp(argv[i],"-nobs")==0)
	sscanf(argv[i+1],"%d",&MAX_OBSN);
      else if (strcasecmp(argv[i],"-debug")==0)
	debugFlag=1;
      else if (strcasecmp(argv[i],"-veryfast")==0)
	veryFast=1;

      strcpy(commandLine[i],argv[i]);
    }
  if ((psr = (pulsar *)malloc(sizeof(pulsar)*MAX_PSR))==NULL)
    {
      printf("Not enough memory to allocate room for %d pulsars\n",MAX_PSR);
      printf("Please decrease the value of MAX_PSR_VAL in tempo2.h\n"); 
      exit(1); 
    }
  if (debugFlag==1) printf("Have allocated memory for pulsar\n");
  psr[0].jboFormat = 0;

  for (i=1;i<argc;i++)
    {
      if (strcmp(commandLine[i],"-machine")==0)
	strcpy(tempo2MachineType,commandLine[++i]);
      else if (strcasecmp(commandLine[i],"-noWarnings")==0)
	noWarnings=2;
      else if (strcasecmp(commandLine[i],"-allInfo")==0)
	noWarnings=0;
      else if (strcmp(commandLine[i],"-jbo")==0)
	psr[0].jboFormat=1;
      else if (strcmp(commandLine[i],"-test")==0) /* Use TEMPO2_TEST environment variable */
	strcpy(TEMPO2_ENVIRON,"TEMPO2_TEST");
      else
	{
	  char oldCommandLine[1000][1000];
	  int  oldArgc = argc;
	  int  oldI = i;

	  for (k=i;k<argc;k++)
	    strcpy(oldCommandLine[k-i],commandLine[k]);

	  sprintf(str,"%s/alias.dat",getenv(TEMPO2_ENVIRON));	  
	  if ((alias = fopen(str,"r")))
	    {
	      while (!feof(alias))
		{
		  char *ret;
		  fgets(str,MAX_FILELEN,alias);
		  ret = strtok(str," ");
		  if (strcmp(ret,commandLine[i])==0)
		    {
		      printf("Using alias: %s\n",ret);
		      /* Have an alias: now split up the remainder of the line */
		      while ((ret = strtok(NULL," "))!=NULL)
			{
			  strcpy(commandLine[i],ret);
			  if (commandLine[i][strlen(commandLine[i])-1]=='\n')
			    commandLine[i][strlen(commandLine[i])-1] = '\0';
			  i++;
			  argc++;
			}
		      argc--;
		      for (k=1;k<oldArgc-oldI;k++)
			strcpy(commandLine[i+k-1],oldCommandLine[k]);
		      i=oldI-1;

		    }
		}
	      fclose(alias);
	    }
	}
    }
  int ii;
  for(ii=1;ii<argc;ii++){
      if (strcasecmp(commandLine[ii],"-reminder")==0){
	  // Writing command line to log file
	  char commandfile[200] = "T2command.input";
	  FILE *fout;
	  time_t rawtime;
	  struct tm * timeinfo;
	  time ( &rawtime );
	  timeinfo = localtime ( &rawtime );
	  char timeval[200];
	  strcpy(timeval,asctime (timeinfo));
	  strcpy(&timeval[(int)strlen(timeval)-1],"");
	  fout = fopen(commandfile,"a");
	  fprintf(fout,"[%s]>> ",timeval);
	  for(i=0;i<argc;i++){
	      fprintf(fout," %s ",commandLine[i]);
	  }
	  fprintf(fout,"\n");
	  fclose(fout);
      }
  }
  /* If running from the command line ... */
  if (debugFlag==1) printf("Running initialise\n");
  initialise(psr,noWarnings); /* Initialise all */
  if (debugFlag==1) printf("Completed running initialise %d\n",psr[0].nits);
  /* Obtain login architecture */
  if (strlen(tempo2MachineType)==0)
    {
#ifdef  TEMPO2_ARCH 
      strcpy(tempo2MachineType, TEMPO2_ARCH);
#else
      if (getenv("LOGIN_ARCH")==NULL)
	{
	  printf("Unable to determine machine type: You must do one of the following:\n"
                 "Re-compile tempo2 with the standard export distrubution, or\n"
                 "Set the LOGIN_ARCH environment variable, or\n"
                 "Use -machine on the command line\n");
	  exit(1);
	}
      strcpy(tempo2MachineType, getenv("LOGIN_ARCH"));
#endif
    }

  if (sizeof(longdouble)!=16 && noWarnings<1)
    {
      printf("Warning: the size of a longdouble is only %d bytes\n",sizeof(longdouble));
      printf(" --- the size of a double is %d bytes\n",sizeof(double));
    }
  strcpy(outputSO,"");
  if (argc==1) /* No command line arguments */
    {
      printf("TEMPO2.  Usage: tempo2 XXX.tim\n");
      exit(1);
    }
  npsr = 0;   /* Initialise the number of pulsars */
  displayParams=0;
  nGlobal=0;
  /* Obtain command line arguments */
  if (debugFlag==1) printf("Running getInputs %d\n",psr[0].nits);
  getInputs(psr,argc, commandLine, timFile,parFile,&listparms,&npsr,&nGlobal,&outRes,&writeModel,
	    outputSO,&flagPolyco,polyco_args,&newpar,&onlypre,dcmFile);
  if (debugFlag==1) printf("Completed getInputs\n");

  for (i=1;i<argc;i++)
    {
      if (strcmp(commandLine[i],"-gr")==0 || strcmp(commandLine[i],"-gr2")==0) 
	/* Running from a graphical interface? */
	{      
	  char *(*entry)(int,char **,pulsar *,int *);
	  void * module;
	  if (strcmp(commandLine[i],"-gr")==0)
	    {
	      sprintf(str,"%s/plugins/%s_%s_plug.so",getenv(TEMPO2_ENVIRON),
		      commandLine[i+1],tempo2MachineType);
	    }
	  else
	    sprintf(str,"./%s_%s_plug.so",commandLine[i+1],tempo2MachineType);
	  printf("Looking for %s\n",str);
	  module = dlopen(str, RTLD_NOW); 
	  if(!module)  {
	    fprintf(stderr, "[error]: dlopen() failed while resolving symbols.\n" );
	    fprintf(stderr, "dlerror() = %s\n",dlerror());
	    return -1;
	  }
	  entry = (char*(*)(int,char **,pulsar *,int *))dlsym(module, "graphicalInterface");
	  if( entry == NULL ) {
	    dlclose(module);
	    fprintf(stderr, "[error]: dlerror() failed while  retrieving address.\n" ); 
	    fprintf(stderr, "dlerror() = %s\n",dlerror());
	    return -1;
	  }
	  entry(argc,commandLine,psr,&npsr);
	  return 0;
	}
    }
  if (debugFlag==1) printf("Reading par file\n");
  readParfile(psr,parFile,timFile,npsr); /* Read .par file to define the pulsar's initial parameters */  
  if (debugFlag==1) printf("Finished reading par file %d\n",psr[0].nits);
  if (flagPolyco==0)
    {
      if (debugFlag==1) printf("Running readTimfile\n");
      readTimfile(psr,timFile,npsr); /* Read .tim file to define the site-arrival-times */
      if (debugFlag==1) printf("Completed readTimfile %d\n",psr[0].param[param_ecc].paramSet[1]);
    }

  if (debugFlag==1) printf("Running preProcess %d\n",psr[0].nits);
  preProcess(psr,npsr,argc,commandLine);
  if (debugFlag==1) printf("Completed preProcess %d\n",psr[0].nits);
  if (flagPolyco> 0)  /* Running tempo2 in polyco mode? */
  {
    if (flagPolyco == 1)
    {
      longdouble mjd1, mjd2,maxha,freq;
      int ncoeff,nspan; 
      double seg_length; 
      char sitename[128];
      if (sscanf(polyco_args, "%Lf %Lf %d %d %Lf %s %Lf", &mjd1, &mjd2, &nspan,
		 &ncoeff, &maxha,sitename,&freq)!=7)
      {
	fprintf(stderr, "Error parsing -polyco arguments! See tempo2 -h.\n");
	printf("Have: %s\n",polyco_args);
	exit(1);
      }
      if (psr[0].tempo1 == 0)
	printf("WARNING: Should probably use -tempo1 option\n");

      if (psr[0].param[param_tzrmjd].paramSet[0]==0)
	{
	  printf("WARNING: tzrmjd not set.  Setting to %g\n",(double)psr[0].param[param_pepoch].val[0]);
	  psr[0].param[param_tzrmjd].paramSet[0]=1;
	  psr[0].param[param_tzrmjd].val[0] = psr[0].param[param_pepoch].val[0];
	}
      if (psr[0].param[param_tzrfrq].paramSet[0]==0)
	{
	  printf("WARNING: tzrfrq not set.  Setting to %g\n",(double)1400.0);
	  psr[0].param[param_tzrfrq].paramSet[0]=1;
	  psr[0].param[param_tzrfrq].val[0] = 1400.0;
	}
      if (strlen(psr[0].tzrsite)<1)
	{
	  printf("WARNING: tzrsite not set.  Setting to %s\n",sitename);
	  strcpy(psr[0].tzrsite,sitename);
	}


      polyco(psr,npsr,mjd1, mjd2,nspan,ncoeff,maxha,sitename,freq,coeff,1);
    }
    else if (flagPolyco==2)
    {
      long double seg_length;
      int ntimecoeff, nfreqcoeff;
      char sitename[64];
      long double mjd_start, mjd_end;
      long double freq_start, freq_end;

      if (sscanf(polyco_args, "%s %Lf %Lf %Lf %Lf %d %d %Lf", sitename,
		 &mjd_start, &mjd_end, &freq_start, &freq_end,
		 &ntimecoeff, &nfreqcoeff, &seg_length)!=8)
      {
	fprintf(stderr, "Error parsing -pred arguments! See tempo2 -h.\n");
	printf("Have: %s\n",polyco_args);
	exit(1);
      }
      /* Actually want seg_length in days */
      seg_length /= SECDAY;

      if (ntimecoeff%2)
      {
	ntimecoeff++;
	printf("Adjusting number of coefficients (time axis) to %d (must be even)\n", ntimecoeff);
      }
      if (nfreqcoeff%2)
      {
	nfreqcoeff++;
	printf("Adjusting number of coefficients (frequency axis) to %d (must be even)\n", nfreqcoeff);
      }
      if (psr[0].param[param_tzrmjd].paramSet[0]==0)
	{
	  printf("WARNING: tzrmjd not set.  Setting to %g\n",(double)psr[0].param[param_pepoch].val[0]);
	  psr[0].param[param_tzrmjd].paramSet[0]=1;
	  psr[0].param[param_tzrmjd].val[0] = psr[0].param[param_pepoch].val[0];
	}
      if (psr[0].param[param_tzrfrq].paramSet[0]==0)
	{
	  printf("WARNING: tzrfrq not set.  Setting to %g\n",(double)1400.0);
	  psr[0].param[param_tzrfrq].paramSet[0]=1;
	  psr[0].param[param_tzrfrq].val[0] = 1400.0;
	}
      if (strlen(psr[0].tzrsite)<1)
	{
	  printf("WARNING: tzrsite not set.  Setting to %s\n",sitename);
	  strcpy(psr[0].tzrsite,sitename);
	}

      ChebyModelSet cms;
      ChebyModelSet_Construct(&cms, psr, sitename, mjd_start, mjd_end,
			      seg_length, seg_length*0.1, 
			      freq_start, freq_end, ntimecoeff, nfreqcoeff);
      FILE *f = fopen("t2pred.dat", "w");
      if (!f)
      {
	fprintf(stderr, "Could not open t2pred.dat for writing!\n");
	exit(1);
      }
      ChebyModelSet_Write(&cms, f);
      fclose(f);
      long double rms, mav;
      ChebyModelSet_Test(&cms, psr, ntimecoeff*5*cms.nsegments, 
			 nfreqcoeff*5*cms.nsegments, &rms, &mav);
      printf("Predictive model constructed and written to t2pred.dat.\n");
      printf("RMS error = %.3Lg s MAV= %.3Lg s\n", 
	     rms/psr[0].param[param_f].val[0], mav/psr[0].param[param_f].val[0]);
      ChebyModelSet_Destroy(&cms);
    }

    return 0;
  }
  if (debugFlag==1)
    {
      printf("Number of iterations = %d\n",psr[0].nits);
      printf("Maximum number of parameters = %d\n",MAX_PARAMS);
      printf("Number of pulsars = %d\n",npsr);
    }
  for (it=0;it<psr[0].nits;it++) /* Why pulsar 0 should select the iterations? */
    {
      if (it>0) /* Copy post-fit values to pre-fit values */
	{
	  for (i=0;i<MAX_PARAMS;i++)
	    {
	      for (p=0;p<npsr;p++)
		{
		  for (k=0;k<psr[p].param[i].aSize;k++)
		    {
		      psr[p].param[i].prefit[k] = psr[p].param[i].val[k];
		      psr[p].param[i].prefitErr[k] = psr[p].param[i].err[k];
		    }
		}
	    }
	}
      for (iteration=0;iteration<2;iteration++) /* Do pre- and post- fit analysis */
	{
	  if (debugFlag==1) printf("iteration %d\n",iteration);
	  if (debugFlag==1) printf("calling formBatsAll\n");
	  printf("Formning bats\n");
	  formBatsAll(psr,npsr);                /* Form Barycentric arrival times */
	  if (debugFlag==1) printf("calling formResiduals\n");
	  printf("Formning resids\n");
	  formResiduals(psr,npsr,iteration);       /* Form residuals */
	  printf("Done resid\n");
	  if (listparms==1 && iteration==0)displayParameters(13,timFile,parFile,psr,npsr); /* List out all the parameters */  
	  if (iteration==0)          /* Only fit to pre-fit residuals */
	    {
	      if (debugFlag==1) printf("calling doFit\n");

	      if (strcmp(dcmFile,"NULL")==0)
		doFit(psr,npsr,writeModel); /* Fit to the residuals to obtain updated parameters */
	      else
		doFitDCM(psr,dcmFile,npsr,writeModel);
	      printf("Complete return\n");
	      /* doFitGlobal(psr,npsr,&globalParameter,nGlobal,writeModel);*/ /* Fit to the residuals to obtain updated parameters  */
	      if (debugFlag==1) printf("completed doFit\n");
	    }
	  if (iteration==1 || onlypre==1)
	    {
	      if (strlen(outputSO)==0)
		textOutput(psr,npsr,globalParameter,nGlobal,outRes,newpar,"new.par"); /* Output results to the screen */
	      else  /* Use a plug in for the output */
		{
		  char *(*entry)(int,char **,pulsar *,int);
		  void * module;
		  sprintf(str,"%s/plugins/%s_%s_plug.so",getenv(TEMPO2_ENVIRON),
			  outputSO,tempo2MachineType);
		  
		  /*	      sprintf(str,"%s/%s_%s_plug.so",tempo2plug,outputSO,tempo2MachineType); */
		  module = dlopen(str, RTLD_NOW); 
		  if(!module)  {
		    fprintf(stderr, "[error]: dlopen() failed while resolving symbols.\n" );
		    return -1;
		  }
		  entry = (char*(*)(int,char **,pulsar *,int))dlsym(module, "tempoOutput");
		  if( entry == NULL ) {
		    dlclose(module);
		    fprintf(stderr, "[error]: dlerror() failed while  retrieving address.\n" ); 
		    return -1;
		  }
		  entry(argc,argv,psr,npsr);
		}
	    }
	  psr[0].noWarnings=2;
	  if (onlypre==1) iteration=2;
		textOutput(psr,npsr,globalParameter,nGlobal,outRes,newpar,"new.par"); /* Output results to the screen */
	  printf("Next iteration\n");
	}
    }
  endClock = clock();
  printf("Finishing off: time taken = %.2f (s)\n",(endClock-startClock)/(float)CLOCKS_PER_SEC);
  exit(EXIT_SUCCESS);
} 

// redwards function to force linkage with library functions used by
// plugins
void
thwart_annoying_dynamic_library_stuff(int never_call_me, float or_sink)
{
  ChebyModel *cm;
  T2Predictor *t2p;
  ChebyModel_Init(cm, 0, 0);
  T2Predictor_GetPhase(t2p, 0, 0);
}
