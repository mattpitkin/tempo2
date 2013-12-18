//  Copyright (C) 2004,2006,2007,2008,2009, George Hobbs, Russell Edwards

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

/* Fake: George Hobbs (Feb 19th 2004) based on fake.f from Simon Johnston     */
/* Purpose is to produce an arrival time file that simulates a pulsar defined */
/* using a .par file                                                          */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"
#include "T2toolkit.h"

void callFit(pulsar *psr,int npsr);

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  longdouble imjd=-1.0,fmjd=-1.0,ra,grms=0.0;
  longdouble toa,ha0,almst,amjd,hnobs,mjd,trmjd,tstep=0.0;
  longdouble times[MAX_OBSN];
  longdouble out[MAX_OBSN];
  longdouble lowFreq,highFreq,alpha=-3.0,ampPL=1.0e-16;
  int setref=0;
  int nshots,npts;
  int j,count=0,i,k,ii,jj,kk;
  longdouble solsid  = 1.002737909;
  longdouble obslong = 149.0;  /* Longitude of Parkes */
  longdouble freq    = 1440.0;
  long iseed;
  int    site    = 7,p;
  int endit=0;
  long idum = 0;
  int nday = -1;
  int nit=4;
  float ngap=-1.0;
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  char str[MAX_FILELEN],str2[MAX_FILELEN];
  double hamax =-1.0;
  char random[100];
  char read_f=0;
  FILE *fout,*fin;
  char addCubic[100]="n",temp[100];
  char formstr[50]="tempo2";
  char smooth[100]="n";
  int  psrNum,giveRMS=-1,setrand=-1,setred=-1,Npsr=0;
  int timesFile=0;
  char fake_fname[100];
  char have_outfile=0;
  char outfile[100];
  char timesfname[100];
  char telID[128]="7"; // Hardcode to Parkes
  int bunching=0; // flag on whether or not observations occur in groups
  // size of gap between observing runs, and length of observing runs.
  // These defaults give 7 observations every 28 days.
  long double gapsize=21,hillsize=7,gapstartmjd;
  // Flag whether or not to ask for red noise variables.
  *npsr = 1;

  strcpy(fake_fname,"fake.rf");

  for (i=0;i<argc;i++){
    if(strcmp(argv[i],"-ndobs")==0){
      sscanf(argv[i+1],"%f",&ngap);
      printf("Have >>%f<< days between observations\n",ngap);
    }
    if(strcmp(argv[i],"-nobsd")==0){
      sscanf(argv[i+1],"%d",&nday);
      printf("Have >>%d<< observations per day\n",nday);
    }
    if (strcmp(argv[i],"-idum")==0){
	sscanf(argv[i+1],"%d",&idum);
	printf("Have idum >>%d<<\n",idum);
    }
    if(strcmp(argv[i],"-ha")==0){
      sscanf(argv[i+1],"%lf",&hamax);
      printf("Have maximum absolute HA >>%lf<<\n",hamax);
    }
    if(strcmp(argv[i],"-randha")==0){
      strcpy(&random[0],argv[i+1]);
      printf("Have random >>%s<<\n",random);
      setrand=1;
    }
    if(strcmp(argv[i],"-start")==0){
      sscanf(argv[i+1],"%Lf",&imjd);
      printf("Have initial MJD >>%lf<<\n",(double)imjd);
    }
    if(strcmp(argv[i],"-tel")==0){
      sscanf(argv[i+1],"%s",telID);
    }
    if(strcmp(argv[i],"-end")==0){
      sscanf(argv[i+1],"%Lf",&fmjd);
      printf("Have final MJD >>%lf<<\n",(double)fmjd);
    }
    if(strcmp(argv[i],"-rms")==0){
      sscanf(argv[i+1],"%Lf",&grms);
      giveRMS=1;
      printf("Have Gaussian noise rms >>%lf<<\n",(double)grms);
    }
    if(strcmp(argv[i],"-format")==0){
	sscanf(argv[i+1],"%s",&formstr);
	printf("Have output format >>%s<<\n",formstr);
    }
    if(strcmp(argv[i],"-nit")==0){
	sscanf(argv[i+1],"%d",&nit);
    }
    if (strcmp(argv[i],"-times")==0){
      timesFile=1;
      sscanf(argv[i+1],"%s",&timesfname);
      printf("Timesfile = %s\n",timesfname);
    }
    if (strcmp(argv[i],"-setref")==0)
      setref=1;

    if (strcmp(argv[i],"-o")==0){
      have_outfile=1;
      sscanf(argv[i+1],"%s",&outfile);
      printf("outfile = %s\n",outfile);
    }

    if (strcmp(argv[i],"-readtim")==0){
      read_f=1;
      printf("Read name,freq,mjd\n");
    }


    if(strcmp(argv[i],"-group")==0){
	bunching = 1;
	sscanf(argv[i+1],"%Lf",&hillsize);
	sscanf(argv[i+2],"%Lf",&gapsize);
	printf("Will simulate runs of %lg observing days long and leave a gap of %lg days between runs.\n",
	       (double)hillsize,(double)gapsize);
    }
    if(strcmp(argv[i],"-h")==0){
      printf("==========================================================================================\n");
      printf(" fake Tempo2 plugin - usage instructions.\n");
      printf(" tempo2 -gr fake -f file.par: The program will prompt you for parameters.\n");
      printf("\n Command-line arguments:\n");
      printf(" \t -f J0437-4715.par J1909-3744.par J1713+0747.par: specify a number of parfiles.\n");
      printf(" \t -ndobs xxx: specify number of days between observations.\n");
      printf(" \t -nobsd xxx: specify number of observations per day.\n");
      printf(" \t -ha xxx: specify maximal absolute Hour Angle.\n");
      printf(" \t -randha y: specify whether to use random HA coverage or not (y/n).\n");
      printf(" \t -start xxxxx: specify start MJD.\n");
      printf(" \t -end xxxxx: specify final MJD.\n");
      printf(" \t -rms xxx: specify Gaussian noise rms (in ms).\n");
      printf(" \t -times xxx: read observation times from specified file\n");
      printf(" \t           suppresses questions for red noise characteristics.\n");
      printf(" \t -format parkes : sets the tim-file format to tempo. (Default: tempo2).\n");
      printf(" \t -group 7 21 : simulate observations in groups of 7 days, with 21 days between groups.\n");
      printf(" \t               There will be nobsd observations every ndobs days during these 7 days.\n");
      printf("\n\n \t -idum xxx: specify random number seed (default = set from clock)\n");
      printf("\n The program will prompt you for the parameters not defined in the command line.\n");
      printf("\n\n Have a nice day!\n\n");
      printf("==========================================================================================\n");
      exit(0);
    }
    if(strcmp(argv[i],"-f")==0){
      Npsr=0;
      while(argv[i+Npsr+1][0]!='-'){
	strcpy(parFile[Npsr],argv[i+Npsr+1]);
	printf("Have parfile %d: >>%s<<.\n",Npsr+1,parFile[Npsr]);
	Npsr++;
	if(i+Npsr+1>=argc) break;
      }
      printf("Have %d parameter files.\n",Npsr);
    }
  }
  if (timesFile==1)
    {
      nday = 1;
      ngap = 1;
      hamax = 8;
      setrand = 1;
      imjd = 1;
      fmjd = 1;
    }
  printf("Simulate arrival times for parameter file %s\n",argv[3]);
  printf("----------------------------------------------------\n\n");
  if(ngap<0){    printf("Enter number of days between observations . "); scanf("%f",&ngap);}
  if(nday<0) {    printf("Enter number of observations/day .......... "); scanf("%d",&nday);}
  if(hamax<0 && nday!=1)  {  
      printf("Enter max absolute HA...................... "); scanf("%lf", &hamax);}
  if(setrand!=1){ printf("Random HA coverage? (y/n) ................. "); scanf("%s",&random);}
  if(imjd<0)   {  printf("Enter initial MJD ......................... "); scanf("%Lf",&imjd);}
  if(fmjd<0)   {  printf("Enter final MJD ........................... "); scanf("%Lf",&fmjd);}
  if(giveRMS<0){
    printf("Enter Gaussian noise rms  (ms/auto)........ "); scanf("%s",temp);
    giveRMS = sscanf(temp,"%Lf",&grms);
  }

  printf("GIVE RMS = %d\n",giveRMS);
  if (idum==0)
    {
      printf("Setting random number seed from the clock\n");      
      idum = TKsetSeed();
    }


  iseed = idum;
  psrNum = Npsr;
  if (timesFile==1)
    fin = fopen(timesfname,"r");

  hnobs = nday/2.0; 
  for(ii=0;ii<Npsr;ii++){
    count = 0;
    strcpy(parFile[0],parFile[ii]);
    psr[0].nJumps=0;
    psr[0].fitMode=0;
    psr[0].eclCoord=0;
    psr[0].nits=1;
    psr[0].clockFromOverride[0] = '\0';
    psr[0].nCompanion = 0;
    psr[0].bootStrap = 0;
    psr[0].units = SI_UNITS;
    psr[0].ne_sw  = NE_SW_DEFAULT; 
    psr[0].nWhite = 0;  /* No whitening by default */
    psr[0].timeEphemeris = IF99_TIMEEPH;
    psr[0].dilateFreq = 1;
    psr[0].planetShapiro = 1;
    psr[0].correctTroposphere = 1;
    psr[0].t2cMethod = T2C_IAU2000B;
    psr[0].fixedFormat=0;
    psr[0].nStorePrecision=0;
    strcpy(psr[0].deleteFileName,"NONE");
    strcpy(psr[0].tzrsite,"NULL");
    psr[0].calcShapiro=1;
    psr[0].ipm = 1;
    psr[0].swm = 0;
    psr[0].nPhaseJump=0;
    psr[0].nobs = 0;

    for (i=0;i<MAX_PARAMS;i++)
      {
	for (j=0;j<psr[0].param[i].aSize;j++)
	  {
	    psr[0].param[i].fitFlag[j] = 0;
	    psr[0].param[i].paramSet[j] = 0;
	    psr[0].param[i].err[j] = 0;
	    psr[0].param[i].val[j] = 0;
	  }
      }
    //      initialise(psr,0);              /* Initialise the structures */      
    //      printf("SET AFTER 1: %d\n",psr[0].param[param_pb].val[0]);
    readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
    //      printf("SET AFTER 2: %d %s\n",psr[0].param[param_pb].val[0],psr[0].name);
      if (psr[0].param[param_track].paramSet[0]==1)
	psr[0].param[param_track].paramSet[0] = 0;
      ra = (double)psr[0].param[param_raj].val[0]/2.0/M_PI;
      
      /* Code based on fake1.f to calculate the TOA at transit for each of these observations */
      trmjd = imjd;
      /* 47892.0 = 1990??? */
      almst = fortran_mod((trmjd-47892.0)*solsid+0.276105324+obslong/360.0,(longdouble)1.0);
      
      /* Hour angle at 00h UT */
      ha0 = almst - ra;
      /* Approximate transit time */
      if (ha0 < 0.0) ha0+=1.0;
      trmjd += 1.0-ha0;
      
      if (nday > 1)tstep = hamax/12.0/nday;  /* Was 0.4/nday */
      amjd =  trmjd;
      gapstartmjd = imjd+(hillsize)/solsid;

      do {
	  for (j=0;j<nday;j++)
	      {
		if (timesFile==1)
		  {
		    if (read_f){
		      if (fscanf(fin,"%s %Lf %Lf\n",fake_fname,&freq,&mjd)==3)
			{
			  printf("Read %g %f\n",(double)mjd, (double)freq);
			  endit=0;
			}
		      else
			endit=1;
		    } else{
		      if (fscanf(fin,"%Lf",&mjd)==1)
			{
			  printf("Read %g\n",(double)mjd);
			  endit=0;
			}
		      else
			endit=1;
		    }
		  }
		else
		  {
		    if (random[0]=='y'||random[0]=='Y')
		      amjd=trmjd + (rand()/(longdouble)RAND_MAX - 0.5)*hamax/12.0;
		    else if (nday==1)
		      amjd=trmjd;
		    else
		      amjd=trmjd + ((j+1)-hnobs)*tstep;
		    
		    mjd=amjd;
		  }
		if (endit==0)
		  {
		    if (count==0 && setref==1)
		      {
			psr[0].obsn[count].sat    = psr[0].param[param_tzrmjd].val[0];
			strcpy(psr[0].obsn[count].fname,"reference");
			psr[0].obsn[count].freq   = psr[0].param[param_tzrfrq].val[0];
			if (giveRMS!=1) grms = psr[0].param[param_tres].val[0]/1e3;
			//	    else grms=0.0;
			psr[0].obsn[count].toaErr = grms*1000.0;
			psr[0].obsn[count].origErr = grms*1000.0;
			psr[0].obsn[count].phaseOffset = 0.0;
			strcpy(psr[0].obsn[count].telID, psr[0].tzrsite);
			psr[0].obsn[count].deleted = 0;
			psr[0].obsn[count].clockCorr=1;
			psr[0].obsn[count].delayCorr=1;
			psr[0].obsn[count].efac=1;
			count++;
		      }
		    psr[0].obsn[count].sat    = mjd;
		    strcpy(psr[0].obsn[count].fname,fake_fname);
		    psr[0].obsn[count].freq   = freq;
		    if (giveRMS!=1) grms = psr[0].param[param_tres].val[0]/1e3;
		    //	    else grms=0.0;
		    psr[0].obsn[count].toaErr = grms*1000.0;
		    psr[0].obsn[count].origErr = grms*1000.0;
		    psr[0].obsn[count].phaseOffset = 0.0;
		    strcpy(psr[0].obsn[count].telID, telID);
		    psr[0].obsn[count].deleted = 0;
		    psr[0].obsn[count].clockCorr=1;
		    psr[0].obsn[count].delayCorr=1;
		    psr[0].obsn[count].efac=1;
		    count++;
		    if (count>MAX_OBSN)
		      {
			printf("Number of TOAs > MAX_OBSN.\n");
			count--;
		      }
		  }
	      }
	  if((bunching == 1) && (trmjd >= (gapstartmjd-1))){
	    trmjd += gapsize/solsid;
	    gapstartmjd += (gapsize+hillsize)/solsid;
	  }
	  else{
	    trmjd+=ngap/solsid;
	  }
      }while ((timesFile == 0 && amjd<fmjd) || (timesFile == 1 && endit==0));
      if (timesFile==1)
    fclose(fin);

  
      psr[0].nobs=count;
      
      if (have_outfile){
	      strcpy(str,outfile);
	      strcpy(timFile[0],outfile);
      }else{
	      strcpy(str,parFile[0]);
	      str[strlen(str)-4]='\0';
	      strcat(str,".simulate");
	      strcpy(timFile[0],str);
      }
      
      /* Now run the tempo2 code */
      preProcess(psr,*npsr,argc,argv);
      callFit(psr,*npsr);             /* Do all the fitting routines */

      for (j=0;j<nit;j++)
	{
	  /* Now update the site arrival times depending upon the residuals */
	  
	  for (i=0;i<psr[0].nobs;i++)  
	    {
	      psr[0].obsn[i].sat -= psr[0].obsn[i].prefitResidual/SECDAY; 
	      psr->obsn[i].nFlags = 0;
	    } 

	  writeTim(str,psr,"tempo2");
	  //	  initialise(&psr[ii],0);
	  // Reset the jumps
	  psr[0].nJumps = 0;
	  for(kk=0;kk<MAX_JUMPS;kk++){
	      psr[0].jumpVal[kk] = 0.0;
	      psr[0].jumpValErr[kk] = 0.0;
	  }
	  for(jj=0;jj<MAX_PARAMS;jj++){
	      psr[0].param[jj].nLinkTo = 0;
	      psr[0].param[jj].nLinkFrom = 0;
	  }
	  psr[0].nconstraints = 0;
	  psr[0].nobs = 0;

	  // SOMETHING TO FIX HERE!!

	  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
	  readTimfile(psr,timFile,*npsr); 
	  preProcess(psr,1,argc,argv);
	  /* Now run the tempo2 code again */
	  callFit(psr,*npsr);             /* Do all the fitting routines */
	}

      printf("Complete %d iterations\n",nit);
      for (i=0;i<psr[0].nobs;i++)  
	{
	  psr[0].obsn[i].sat -= psr[0].obsn[i].prefitResidual/SECDAY;  
	  psr->obsn[i].nFlags = 0;
	}
      
      
      for (i=0;i<psr[0].nobs;i++)
	{ 
	  times[i] = (psr[0].obsn[i].sat-psr[0].param[param_posepoch].val[0]);
	}
      npts = psr[0].nobs;
      
      /* Add Gaussian noise */
      if (giveRMS!=1) grms = psr[0].param[param_tres].val[0]/1e3;
      if (grms>0.0)
	{
	  printf("Adding Gaussian noise with rms = %f\n",(float)grms);
	  for (i=0;i<psr[0].nobs;i++)
	    psr[0].obsn[i].sat += TKgaussDev(&idum)*grms/1000.0/SECDAY;
	}
      
      printf("Output TOA file written to %s\n",str);
      writeTim(str,psr,formstr);      
  }
}

/* This function calls all of the fitting routines.             */
/* The function is looped twice, the first time for the pre-fit */
/* residuals and the second time for the post-fit residuals     */

void callFit(pulsar *psr,int npsr){
  int iteration;
  double globalParameter = 0.0;

  for (iteration=0;iteration<2;iteration++)
    {
      formBatsAll(psr,npsr);
            
      /* Form residuals */
      formResiduals(psr,npsr,0);
      
      /* Do the fitting */
      if (iteration==0) doFit(psr,npsr,0);
      else textOutput(psr,npsr,globalParameter,0,0,0,"");
    }
}

char * plugVersionCheck = TEMPO2_h_VER;
