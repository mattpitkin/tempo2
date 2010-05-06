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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "tempo2.h"
#include "tempo2Util.h"

int readValue(pulsar *psr,char *pmtr,FILE *fin,parameter *parameter,int arr);
int getValue(char *str,int v1,int v2,pulsar *psr,int l,int arr);
void removeCR(char *str);
void checkLine(pulsar *p,char *str,FILE *fin,parameter *elong,parameter *elat);
void checkAllSet(pulsar *psr,parameter elong,parameter elat,char *filename);

/* Function to set up default parameters before reading a .par file */
int setupParameterFileDefaults(pulsar *psr)
{
  strcpy(psr->clock, "TT(TAI)");
  if (getenv(TEMPO2_ENVIRON)==NULL)
    {
      printf("Environment variable >%s< not set\n",TEMPO2_ENVIRON);
      return -1;
    }
  strcpy(psr->JPL_EPHEMERIS,getenv(TEMPO2_ENVIRON));
  strcpy(psr->ephemeris,"DE405");
  strcat(psr->JPL_EPHEMERIS,"/ephemeris/DE405.1950.2050");
  return 0;
}

/* ******************************************** 
 * readSimpleParfile
 *
 * a simple read of a standard T2 parameter file 
 * mainly used for PSRCHIVE applications
 *
 * Currently not working for ecliptic coordinates
 */
int readSimpleParfile (FILE *fin, pulsar *p)
{
  char str[MAX_STRLEN];
  int nread,j;
  parameter elong,elat;

  elong.aSize = 1;
  elong.val       = (longdouble *)malloc(elong.aSize*sizeof(longdouble));
  elong.err       = (longdouble *)malloc(elong.aSize*sizeof(longdouble));
  elong.prefit    = (longdouble *)malloc(elong.aSize*sizeof(longdouble));
  elong.prefitErr = (longdouble *)malloc(elong.aSize*sizeof(longdouble));
  elong.fitFlag   = (int *)malloc(elong.aSize*sizeof(int));
  elong.paramSet  = (int *)malloc(elong.aSize*sizeof(int));
  elong.label     = (char **)malloc(elong.aSize*sizeof(char *));
  elong.shortlabel= (char **)malloc(elong.aSize*sizeof(char *));

  for (j=0;j<elong.aSize;j++)
    {
      elong.label[j] = (char *)malloc(sizeof(char)*100);
      elong.shortlabel[j] = (char *)malloc(sizeof(char)*100);
    }

  for (j=0;j<elong.aSize;j++)
    {
      elong.paramSet[j] = 0;
      elong.err[j]      = 0.0;	     
    }

  elat.aSize = 1;
  elat.val       = (longdouble *)malloc(elat.aSize*sizeof(longdouble));
  elat.err       = (longdouble *)malloc(elat.aSize*sizeof(longdouble));
  elat.prefit    = (longdouble *)malloc(elat.aSize*sizeof(longdouble));
  elat.prefitErr = (longdouble *)malloc(elat.aSize*sizeof(longdouble));
  elat.fitFlag   = (int *)malloc(elat.aSize*sizeof(int));
  elat.paramSet  = (int *)malloc(elat.aSize*sizeof(int));
  elat.label     = (char **)malloc(elat.aSize*sizeof(char *));
  elat.shortlabel= (char **)malloc(elat.aSize*sizeof(char *));

  for (j=0;j<elat.aSize;j++)
    {
      elat.label[j] = (char *)malloc(sizeof(char)*100);
      elat.shortlabel[j] = (char *)malloc(sizeof(char)*100);
    }

  for (j=0;j<elat.aSize;j++)
    {
      elat.paramSet[0] = 0;
      elat.err[j]      = 0.0;	     
    }

  while (!feof(fin))
    {
      // Read in a line from the parameter file
      nread = fscanf(fin,"%s",str);
      if (nread==1)
	checkLine(p,str,fin,&elong,&elat);
    }

  checkAllSet(p,elong,elat,"");

  return 0;
}

void checkLine(pulsar *psr,char *str,FILE *fin,parameter *elong, parameter *elat)
{
  int gval;
  if (str[0]=='#') /* Comment line */
    fgets(str,1000,fin);
  else if (strcasecmp(str,"PSR")==0 || strcasecmp(str,"PSRB")==0 || strcasecmp(str,"PSRJ")==0) /* Name of pulsar */
      fscanf(fin,"%s",psr->name);
  else if (strcasecmp(str,"CLK")==0)
    fscanf(fin,"%s",psr->clock);
  else if (strcasecmp(str,"TRES")==0)
    readValue(psr,str,fin,&(psr->param[param_tres]),0);
  else if (strcasecmp(str,"MODE")==0) /* Fitting mode */
    fscanf(fin,"%d",&(psr->fitMode));
  else if (strcasecmp(str,"NOTRACK")==0)  /* TEMPO2 uses automatic tracking */
    psr->param[param_track].paramSet[0]=0;
  else if (strcasecmp(str,"TRACK")==0)  /* TEMPO2 uses automatic tracking */
    readValue(psr,"TRACK",fin,&(psr->param[param_track]),0);
  else if (strcasecmp(str,"NO_SS_SHAPIRO")==0)
    psr->calcShapiro=-1;
  else if (strcasecmp(str,"IPM")==0)
    fscanf(fin,"%d",&psr->ipm);
  else if (strcasecmp(str,"SWM")==0)
    fscanf(fin,"%d",&psr->swm);
  else if (strcasecmp(str,"NITS")==0)
    fscanf(fin,"%d",&psr->nits);
  else if (strcasecmp(str,"TEMPO1")==0)
    {
      psr->units = TDB_UNITS;
      psr->timeEphemeris = FB90_TIMEEPH;
      psr->dilateFreq = 0;
      psr->planetShapiro=0;
      psr->t2cMethod = T2C_TEMPO;
      psr->correctTroposphere = 0;
    }
  else if (strcasecmp(str,"DILATEFREQ")==0)
    {
      char val[1000];
      fscanf(fin,"%s", val);
      psr->dilateFreq = 
	(val[0]=='1'||val[0]=='y'||val[0]=='Y');
    }
  else if (strcasecmp(str,"IBOOT")==0)
    fscanf(fin,"%d",&(psr->bootStrap));
  else if (strcasecmp(str,"PLANET_SHAPIRO")==0)
    {
      char val[100];
      fscanf(fin,"%s",val);
      psr->planetShapiro=(val[0]=='1'||val[0]=='y'||val[0]=='Y');
    }
  else if (strcasecmp(str,"CORRECT_TROPOSPHERE")==0)
    {
      char val[100];
      fscanf(fin,"%s",val);
      psr->correctTroposphere=(val[0]=='1'||val[0]=='y'||val[0]=='Y');
      printf("Setting correctTroposphere %d\n",psr->correctTroposphere);
    }
  else if (strcasecmp(str,"UNITS")==0)
    {
      char unit[1000];
      fscanf(fin,"%s", unit);
      if (strcasecmp(unit,"TDB")==0) psr->units = TDB_UNITS;
      else if (strcasecmp(unit,"SI")==0) psr->units = SI_UNITS;
    }
  else if (strcasecmp(str,"NE1AU")==0 || strcasecmp(str,"NE_SW")==0)
    fscanf(fin,"%lf",&(psr->ne_sw));
  else if (strcasecmp(str, "TIMEEPH")==0)
    {
      char unit[1000];
      fscanf(fin,"%s", unit);
      if (strcasecmp(unit,"IF99")==0) 
	psr->timeEphemeris = IF99_TIMEEPH;
      else if (strcasecmp(unit,"FB90")==0) 
	psr->timeEphemeris = FB90_TIMEEPH;
    }
  else if (strcasecmp(str, "T2CMETHOD")==0)
    {
      char unit[1000];
      fscanf(fin,"%s", unit);
      if (strcasecmp(unit,"IAU200B")==0) 
	psr->t2cMethod = T2C_IAU2000B;
      else if (strcasecmp(unit,"TEMPO")==0) 
	psr->t2cMethod = T2C_TEMPO;
    }
  else if (strcasecmp(str, "CLK_CORR_CHAIN")==0)
    {
      char rest[1024];
      fgets(rest, 1024, fin);
      defineClockCorrectionSequence(rest,psr->noWarnings);
    }
  else if (strcasecmp(str, "FJUMP")==0)
    {
      fscanf(fin,"%s", psr->fjumpID);
    }
  else if (strcasecmp(str, "JUMP")==0)
    {
      char rest[1024];
      char str1[100],str2[100],str3[100],str4[100],str5[100];
      int v5,nread;

      fgets(rest, 1024, fin);
      removeCR(rest);
      psr->nJumps++;
      strcpy(psr->jumpStr[psr->nJumps],rest);
      psr->fitJump[psr->nJumps]=1; /* Default: fit for jump */
            v5 = -1;
      nread = sscanf(psr->jumpStr[psr->nJumps],"%s %s %s %s %s",str1,str2,str3,str4,str5);
      
      if (strcasecmp(str1,"MJD")==0 || strcasecmp(str1,"FREQ")==0)
	{
	  if (nread>3)
	    {
	      sscanf(str4,"%lf",&(psr->jumpVal[psr->nJumps]));
	      if (sscanf(str5,"%d",&v5)==1)
		{
		  if (v5!=1) psr->fitJump[psr->nJumps]=0;
		}
	      else
		psr->fitJump[psr->nJumps]=0;
	    }
	}
      else if (strcasecmp(str1,"NAME")==0 || strcasecmp(str1,"TEL")==0 || str1[0]=='-')
	{
	  if (nread>2)
	    {
	      sscanf(str3,"%lf",&(psr->jumpVal[psr->nJumps]));
	      if (sscanf(str4,"%d",&v5)==1)
		{
		  if (v5!=1) psr->fitJump[psr->nJumps]=0;
		}
	      else
		psr->fitJump[psr->nJumps]=0;
	    }
	    }
      
    }
  else if (strcasecmp(str,"EPHEM")==0) 
    {
      fscanf(fin,"%s",psr->ephemeris);
      if (strcmp(psr->ephemeris, "DE414")==0)	
	sprintf(psr->JPL_EPHEMERIS,"%s/ephemeris/DE414.1960.2020",getenv("TEMPO2"));
      else if (strcmp(psr->ephemeris,"INPOP06")==0)
	sprintf(psr->JPL_EPHEMERIS,"%s/ephemeris/%s",getenv("TEMPO2"),psr->ephemeris);
      else if (strcmp(psr->ephemeris,"DE421")==0)
	sprintf(psr->JPL_EPHEMERIS,"%s/ephemeris/DE421.1950.2050",getenv("TEMPO2"));
      else
	sprintf(psr->JPL_EPHEMERIS,"%s/ephemeris/%s.1950.2050",getenv("TEMPO2"),psr->ephemeris);
    }
  else if (strcasecmp(str,"TOFFSET")==0) /* Time offset */
    {
      char str[1000];
      int k,nread;
      fgets(str,1000,fin);
      if (str[strlen(str)-1]=='\n') str[strlen(str)-1]='\0';
      nread = sscanf(str,"%lf %lf %lf %lf %s %lf",&(psr->tOffset_f1[psr->nToffset]),
		     &(psr->tOffset_f2[psr->nToffset]),
		     &(psr->tOffset_t1[psr->nToffset]),
		     &(psr->tOffset_t2[psr->nToffset]),
		     psr->tOffsetSite[psr->nToffset],
		     &(psr->tOffset[psr->nToffset]));
      for (k=0;k<(int)strlen(str);k++)
	{
	  if (str[k]=='-') /* Have flag */
	    break;
	}
      if (k==(int)strlen(str)) strcpy(psr->tOffsetFlags[psr->nToffset],"");
      else strcpy(psr->tOffsetFlags[psr->nToffset],str+k);
      psr->nToffset++;
    }
  else if (strcasecmp(str,"PMRA")==0 || strcasecmp(str,"PMLAMBDA")==0 || strcasecmp(str,"PMELONG")==0)      /* Proper motion in RA */
    {
      readValue(psr,str,fin,&(psr->param[param_pmra]),0);
      if (strcasecmp(str,"PMLAMBDA")==0 || strcasecmp(str,"PMELONG")==0)
	psr->eclCoord = 1;
    }
  else if (strcasecmp(str,"PMDEC")==0 || strcasecmp(str,"PMBETA")==0 || strcasecmp(str,"PMELAT")==0)     /* Proper motion in DECJ */
    {
      readValue(psr,str,fin,&(psr->param[param_pmdec]),0);
      if (strcasecmp(str,"PMBETA")==0 || strcasecmp(str,"PMELAT")==0)
	psr->eclCoord = 1;
    }
  else if (strcasecmp(str,"PMRV")==0)     /* Radial velocity */
    readValue(psr,str,fin,&(psr->param[param_pmrv]),0);
  else if (strcasecmp(str,"POSEPOCH")==0)  /* Position Epoch */
    readValue(psr,str,fin,&(psr->param[param_posepoch]),0);
  else if (strcasecmp(str,"WAVEEPOCH")==0)  /* Fitwaves Epoch */
    readValue(psr,str,fin,&(psr->param[param_waveepoch]),0);
  else if (strcasecmp(str,"SIFUNC")==0)  /* Set interpolation function */
    readValue(psr,str,fin,&(psr->param[param_ifunc]),0);
  else if (strcasecmp(str,"PEPOCH")==0)    /* Period Epoch */
    readValue(psr,str,fin,&(psr->param[param_pepoch]),0);
  else if (strcasecmp(str,"EPHVER")==0)    /* Ephemeris version */
    readValue(psr,str,fin,&(psr->param[param_ephver]),0);
  else if (strcasecmp(str,"DMEPOCH")==0)    /* DM Epoch */
    readValue(psr,str,fin,&(psr->param[param_dmepoch]),0);
  else if (strcasecmp(str,"RAJ")==0 || strcasecmp(str,"RA")==0)       /* Right ascension */
    readValue(psr,"RAJ",fin,&(psr->param[param_raj]),0);
  else if (strcasecmp(str,"DECJ")==0 || strcasecmp(str,"DEC")==0)      /* Declination */
    readValue(psr,"DECJ",fin,&(psr->param[param_decj]),0);
  else if (strcasecmp(str,"ELONG")==0 || strcasecmp(str,"LAMBDA")==0)
    readValue(psr,"ELONG",fin,elong,0);
  else if (strcasecmp(str,"ELAT")==0 || strcasecmp(str,"BETA")==0)
    readValue(psr,"ELAT",fin,elat,0);
  else if (strstr(str,"GLEP_")!=NULL || strstr(str,"glep_")!=NULL)
    {
      if (sscanf(str+5,"%d",&gval)==1)
	{
	  if (gval<psr->param[param_glep].aSize)
	    readValue(psr,str,fin,&(psr->param[param_glep]),gval-1);
	}
      
    }
  else if (strstr(str,"GLPH_")!=NULL || strstr(str,"glph_")!=NULL)
    {
      if (sscanf(str+5,"%d",&gval)==1)
	{
	  if (gval<psr->param[param_glph].aSize)
	    readValue(psr,str,fin,&(psr->param[param_glph]),gval-1);
	}
      
    }
  else if (strstr(str,"GLF0_")!=NULL || strstr(str,"glf0_")!=NULL)
    {
      if (sscanf(str+5,"%d",&gval)==1)
	{
	  if (gval<psr->param[param_glf0].aSize)
	    readValue(psr,str,fin,&(psr->param[param_glf0]),gval-1);
	}
      
    }
  else if (strstr(str,"GLF1_")!=NULL || strstr(str,"glf1_")!=NULL)
    {
      if (sscanf(str+5,"%d",&gval)==1)
	{
	  if (gval<psr->param[param_glf1].aSize)
	    readValue(psr,str,fin,&(psr->param[param_glf1]),gval-1);
	}
      
    }
  else if (strstr(str,"GLF2_")!=NULL || strstr(str,"glf2_")!=NULL)
    {
      if (sscanf(str+5,"%d",&gval)==1)
	{
	  if (gval<psr->param[param_glf2].aSize)
	    readValue(psr,str,fin,&(psr->param[param_glf2]),gval-1);
	}
      
    }
  else if (strstr(str,"GLF0D_")!=NULL || strstr(str,"glf0d_")!=NULL)
    {
      if (sscanf(str+6,"%d",&gval)==1)
	{
	  if (gval<psr->param[param_glf0d].aSize)
	    readValue(psr,str,fin,&(psr->param[param_glf0d]),gval-1);
	}
      
    }
  else if (strstr(str,"GLTD_")!=NULL || strstr(str,"gltd_")!=NULL)
    {
      if (sscanf(str+5,"%d",&gval)==1)
	{
	  if (gval<psr->param[param_gltd].aSize)
	    readValue(psr,str,fin,&(psr->param[param_gltd]),gval-1);
	}
      
    }
  else if (strcasecmp(str,"TZRMJD")==0)      /* TZRMJD */
    readValue(psr,str,fin,&(psr->param[param_tzrmjd]),0);
  else if (strcasecmp(str,"TZRSITE")==0)      /* TZRMJD */
    fscanf(fin,"%s",psr->tzrsite);
  else if (strcasecmp(str,"NSPAN")==0 || strcasecmp(str,"TSPAN")==0)      /* TSPAN */
    readValue(psr,str,fin,&(psr->param[param_tspan]),0);
  else if (strcasecmp(str,"TZRFRQ")==0 || strcasecmp(str,"TZRFREQ")==0)      /* TZRFRQ */
    readValue(psr,str,fin,&(psr->param[param_tzrfrq]),0);
  else if (strcasecmp(str,"FDDC")==0)/* Frequency dependent delay */
    readValue(psr,str,fin,&(psr->param[param_fddc]),0);
  else if (strcasecmp(str,"FDDI")==0) /* Frequency dependent delay */
    readValue(psr,str,fin,&(psr->param[param_fddi]),0);
  else if (strcasecmp(str,"DSHK")==0) /* Shklovskii term distance (kpc) */
    readValue(psr,str,fin,&(psr->param[param_dshk]),0);
  else if (strcasecmp(str,"START")==0)      /* START */
    readValue(psr,str,fin,&(psr->param[param_start]),0);
  else if (strcasecmp(str,"FINISH")==0)      /* START */
    readValue(psr,str,fin,&(psr->param[param_finish]),0);
  else if (strcasecmp(str,"P0")==0 || strcasecmp(str,"P")==0)      /* P0 */
    readValue(psr,str,fin,&(psr->param[param_f]),0);
  else if (strcasecmp(str,"P1")==0 || strcasecmp(str,"Pdot")==0)   /* P1 */
    readValue(psr,str,fin,&(psr->param[param_f]),1);
  else if (strcasecmp(str,"F0")==0 || strcasecmp(str,"F")==0)      /* F0 */
    readValue(psr,str,fin,&(psr->param[param_f]),0);
  else if (str[0]=='F' || str[0]=='f') /* Read higher frequency derivatives */
    {
      int fval;
      if (sscanf(str+1,"%d",&fval)==1)
	{
	  if (fval<psr->param[param_f].aSize)
	    readValue(psr,str,fin,&(psr->param[param_f]),fval);
	  else if (fval>=psr->param[param_f].aSize){
	    printf("WARNING!!! Currently only period derivatives up to order 12\n");
	    printf("WARNING!!! are available. All higher derivatives will be ignored!\n");
	  }
	}
    }
  else if (strcasecmp(str,"DM")==0) /* Dispersion measure */
    readValue(psr,str,fin,&(psr->param[param_dm]),0);
  else if ((str[0]=='D' || str[0]=='d') &&  /* Higher DM derivatives */
	   (str[1]=='M' || str[1]=='m') && isdigit(str[2]))
    {
      int dval;
      if (sscanf(str+2,"%d",&dval)==1)
	{
	  if (dval<psr->param[param_dm].aSize)
	    readValue(psr,str,fin,&(psr->param[param_dm]),dval);
	}
    }
  else if (strcasecmp(str,"PX")==0) /* Parallax */
    readValue(psr,str,fin,&(psr->param[param_px]),0);
  else if (strcasecmp(str,"D_AOP")==0) /* AOP Distance (unit: pc)*/
    readValue(psr,str,fin,&(psr->param[param_daop]),0);
  else if (strcasecmp(str,"IPERHARM")==0) /* AOP Distance */
    readValue(psr,str,fin,&(psr->param[param_iperharm]),0);
  
  /* SS Planet masses */
  else if (strncasecmp(str,"DMASSPLANET",11)==0 )
    {
      int k;
      if (sscanf(str+11,"%d",&k)==1 && k >= 1 && k <= 9)
	{
	  readValue(psr,str,fin,&(psr->param[param_dmassplanet]),k-1);
	  //	  printf("Read change of %g\n",(double)psr->param[param_dmassplanet].val[0]);
	}
    }
  
  /* ----------------- */
  /* Whitening params  */
  /* ----------------- */ 
  else if (strcasecmp(str,"WAVE_OM")==0) /* Fundamental frequency */
    readValue(psr,str,fin,&(psr->param[param_wave_om]),0);
  else if (strstr(str,"WAVE")!=NULL || strstr(str,"wave")!=NULL)
    {
      int number;
      /* Obtain parameter number */
      sscanf(str+4,"%d",&number);
      fscanf(fin,"%lf %lf",&psr->wave_sine[number-1],&psr->wave_cos[number-1]);
      if (psr->nWhite < number) psr->nWhite = number;
    }
  else if (strstr(str,"IFUNC")!=NULL || strstr(str,"ifunc")!=NULL)
    {
      int number;
      printf("Setting IFUNC\n");
      /* Obtain parameter number */
      sscanf(str+5,"%d",&number);

      fscanf(fin,"%lf %lf %lf",&psr->ifuncT[number-1],&psr->ifuncV[number-1],&psr->ifuncE[number-1]);
      if (psr->ifuncN < number) psr->ifuncN = number;
      printf("IFUNC n = %d\n",psr->ifuncN);
    }
  /* ---------------- */
  /* Phase jumps      */
  /* ---------------- */
 else if (strcasecmp(str,"PHASE")==0)
    {
      fscanf(fin,"%d %Lf",&psr->phaseJumpDir[psr->nPhaseJump],&psr->phaseJump[psr->nPhaseJump]);
      psr->phaseJumpID[psr->nPhaseJump]=-1;
      psr->nPhaseJump++;
    }  


  /* ----------------- */
  /* Binary Parameters */
  /* ----------------- */
  
  else if (strcasecmp(str,"BINARY")==0) /* Binary model */
    {
      fscanf(fin,"%s",psr->binaryModel);
      if (strcasecmp(psr->binaryModel,"BT1P")==0 ||
	  strcasecmp(psr->binaryModel,"BT2P")==0)
	{
	  displayMsg(1,"BIN3","Converting binary model to T2 model","",psr->noWarnings);
	  strcpy(psr->binaryModel,"T2");
	}
    }
  /* BT Binary parameters */
  else if (strcasecmp(str,"A1")==0)
    readValue(psr,str,fin,&(psr->param[param_a1]),0);
  else if (strcasecmp(str,"A1DOT")!=0 && (str[0]=='A' || str[0]=='a') && 
	   (str[1]=='1'))
    {
      int val;
      if (sscanf(str+3,"%d",&val)==1)
	{
	  if (val-1<psr->param[param_a1].aSize)
	    readValue(psr,str,fin,&(psr->param[param_a1]),val-1);
	}
    }       
  else if (strcasecmp(str,"EDOT")==0 || strcasecmp(str,"ECCDOT")==0)
    readValue(psr,str,fin,&(psr->param[param_edot]),0);
  else if (strcasecmp(str,"E")==0 || strcmp(str,"ECC")==0)
    readValue(psr,str,fin,&(psr->param[param_ecc]),0);
  else if (((str[0]=='E' || str[0]=='e') &&  /* Higher ECC derivatives */
	    (str[1]=='C' || str[1]=='c') &&
	    (str[2]=='C' || str[2]=='c')) || 
	   ((str[0]=='E' || str[0]=='e') && str[1]=='_'))
    {
      int eccval;
      int add=4;
      if ((str[0]=='E' || str[0]=='e') && str[1]=='_') add=2;
      if (sscanf(str+add,"%d",&eccval)==1)
	{
	  if (eccval-1<psr->param[param_ecc].aSize)
	    readValue(psr,str,fin,&(psr->param[param_ecc]),eccval-1);
	}
    }       
  else if (strcasecmp(str,"T0")==0)
    readValue(psr,str,fin,&(psr->param[param_t0]),0);
  else if ((str[0]=='T' || str[0]=='t') && 
	   (str[1]=='0'))
    {
      int val;
      if (sscanf(str+3,"%d",&val)==1)
	{
	  if (val-1<psr->param[param_t0].aSize)
	    readValue(psr,str,fin,&(psr->param[param_t0]),val-1);
	}
    }       
  else if (strcasecmp(str,"OM")==0)
    readValue(psr,str,fin,&(psr->param[param_om]),0);
  /*  else if (strcasecmp(str,"OMDOT")!=0 && (str[0]=='O' || str[0]=='o') && 
	   (str[1]=='M' || str[1]=='m'))
    {
      int val;
      if (sscanf(str+2,"%d",&val)==1)
	{
	  if (val-1<psr->param[param_om].aSize)
	    readValue(psr,str,fin,&(psr->param[param_om]),val);
	}
	}*/       
  else if (strcasecmp(str,"OMDOT")!=0 && (str[0]=='O' || str[0]=='o') && 
	   (str[1]=='M' || str[1]=='m'))
    {
      int val;
      if (sscanf(str+3,"%d",&val)==1)
	{
	  if (val-1<psr->param[param_om].aSize)
	    readValue(psr,str,fin,&(psr->param[param_om]),val-1);
	}
    }
  else if (strcasecmp(str,"PB")==0)
    {
      readValue(psr,str,fin,&(psr->param[param_pb]),0);
      if (psr->nCompanion==0) psr->nCompanion=1;
    }
  else if ((strcasecmp(str,"PBDOT")!=0) && (str[0]=='P' || str[0]=='p') &&  /* Higher Pb derivatives */
	   (str[1]=='B' || str[1]=='b'))
    {
      int pbval;
      if (sscanf(str+3,"%d",&pbval)==1)
	{
	  if (pbval-1<psr->param[param_pb].aSize)
	    {
	      readValue(psr,str,fin,&(psr->param[param_pb]),pbval-1);
	      if (pbval > psr->nCompanion) psr->nCompanion=pbval;
	    }
	}
    }       
  else if (strcasecmp(str,"GAMMA")==0)
    readValue(psr,str,fin,&(psr->param[param_gamma]),0);
  else if (strcasecmp(str,"DR")==0)
    readValue(psr,str,fin,&(psr->param[param_dr]),0);
  else if (strcasecmp(str,"DTH")==0)
    readValue(psr,str,fin,&(psr->param[param_dth]),0);
  else if (strcasecmp(str,"A0")==0)
    readValue(psr,str,fin,&(psr->param[param_a0]),0);
  else if (strcasecmp(str,"B0")==0)
    readValue(psr,str,fin,&(psr->param[param_b0]),0);
  else if (strcasecmp(str,"BP")==0)
    readValue(psr,str,fin,&(psr->param[param_bp]),0);
  else if (strcasecmp(str,"BPP")==0)
    readValue(psr,str,fin,&(psr->param[param_bpp]),0);
  else if (strcasecmp(str,"DTHETA")==0)
    readValue(psr,str,fin,&(psr->param[param_dtheta]),0);
  else if (strcasecmp(str,"PBDOT")==0)
    {
      readValue(psr,str,fin,&(psr->param[param_pbdot]),0);
      if (fabs(psr->param[param_pbdot].val[0]) > 1.0e-7) /* Check units: DO BETTER JOB */
	psr->param[param_pbdot].val[0]*=1.0e-12;
      psr->param[param_pbdot].prefit[0] = psr->param[param_pbdot].val[0];
    }
  else if (strcasecmp(str,"XPBDOT")==0)
    {
      readValue(psr,str,fin,&(psr->param[param_xpbdot]),0);
      if (fabs(psr->param[param_xpbdot].val[0]) > 1.0e-7) /* Check units: DO BETTER JOB */
	psr->param[param_xpbdot].val[0]*=1.0e-12;
      psr->param[param_xpbdot].prefit[0] = psr->param[param_xpbdot].val[0];
    }
  else if (strcasecmp(str,"OMDOT")==0)
    readValue(psr,str,fin,&(psr->param[param_omdot]),0);
  else if (strcasecmp(str,"XOMDOT")==0)
    readValue(psr,str,fin,&(psr->param[param_xomdot]),0);
  else if (strcasecmp(str,"AFAC")==0)
    readValue(psr,str,fin,&(psr->param[param_afac]),0);
  else if (strcasecmp(str,"A1DOT")==0 || strcasecmp(str,"XDOT")==0)
    {
      readValue(psr,str,fin,&(psr->param[param_a1dot]),0);
      if (psr->param[param_a1dot].val[0] > 1e-7) /* Check units: DO BETTER JOB */
	psr->param[param_a1dot].val[0]*=1.0e-12;
      psr->param[param_a1dot].prefit[0] = psr->param[param_a1dot].val[0];
    }
  else if (strcasecmp(str,"A2DOT")==0 || strcasecmp(str,"X2DOT")==0)
    {
      readValue(psr,str,fin,&(psr->param[param_a1dot]),1);
      psr->param[param_a1dot].prefit[1] = psr->param[param_a1dot].val[1];
    }
  else if (strcasecmp(str,"TASC")==0)
    readValue(psr,str,fin,&(psr->param[param_tasc]),0);
  else if (strcasecmp(str,"EPS1")==0)
    readValue(psr,str,fin,&(psr->param[param_eps1]),0);
  else if (strcasecmp(str,"EPS1DOT")==0)
    readValue(psr,str,fin,&(psr->param[param_eps1dot]),0);
  else if (strcasecmp(str,"EPS2")==0)
    readValue(psr,str,fin,&(psr->param[param_eps2]),0);
  else if (strcasecmp(str,"EPS2DOT")==0)
    readValue(psr,str,fin,&(psr->param[param_eps2dot]),0);
  else if (strcasecmp(str,"M2")==0)
    readValue(psr,str,fin,&(psr->param[param_m2]),0);
  else if (strcasecmp(str,"KOM")==0)
    readValue(psr,str,fin,&(psr->param[param_kom]),0);
  else if (strcasecmp(str,"KIN")==0)
    readValue(psr,str,fin,&(psr->param[param_kin]),0);
  else if (strcasecmp(str,"SHAPMAX")==0)
    readValue(psr,str,fin,&(psr->param[param_shapmax]),0);
  else if (strcasecmp(str,"MTOT")==0)
    readValue(psr,str,fin,&(psr->param[param_mtot]),0);
  else if (strcasecmp(str,"SINI")==0)
    readValue(psr,str,fin,&(psr->param[param_sini]),0);
  else if (strcasecmp(str,"BPJEP_1")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjep]),0);
  else if (strcasecmp(str,"BPJEP_2")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjep]),1);
  else if (strcasecmp(str,"BPJEP_3")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjep]),2);
  else if (strcasecmp(str,"BPJEP_4")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjep]),3);
  else if (strcasecmp(str,"BPJEP_5")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjep]),4);
  else if (strcasecmp(str,"BPJPH_1")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjph]),0);
  else if (strcasecmp(str,"BPJPH_2")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjph]),1);
  else if (strcasecmp(str,"BPJPH_3")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjph]),2);
  else if (strcasecmp(str,"BPJPH_4")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjph]),3);
  else if (strcasecmp(str,"BPJPH_5")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjph]),4);
  else if (strcasecmp(str,"BPJA1_1")==0)
    readValue(psr,str,fin,&(psr->param[param_bpja1]),0);
  else if (strcasecmp(str,"BPJA1_2")==0)
    readValue(psr,str,fin,&(psr->param[param_bpja1]),1);
  else if (strcasecmp(str,"BPJA1_3")==0)
    readValue(psr,str,fin,&(psr->param[param_bpja1]),2);
  else if (strcasecmp(str,"BPJA1_4")==0)
    readValue(psr,str,fin,&(psr->param[param_bpja1]),3);
  else if (strcasecmp(str,"BPJA1_5")==0)
    readValue(psr,str,fin,&(psr->param[param_bpja1]),4);
  else if (strcasecmp(str,"BPJEC_1")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjec]),0);
  else if (strcasecmp(str,"BPJEC_2")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjec]),1);
  else if (strcasecmp(str,"BPJEC_3")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjec]),2);
  else if (strcasecmp(str,"BPJEC_4")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjec]),3);
  else if (strcasecmp(str,"BPJEC_5")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjec]),4);
  else if (strcasecmp(str,"BPJOM_1")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjom]),0);
  else if (strcasecmp(str,"BPJOM_2")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjom]),1);
  else if (strcasecmp(str,"BPJOM_3")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjom]),2);
  else if (strcasecmp(str,"BPJOM_4")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjom]),3);
  else if (strcasecmp(str,"BPJOM_5")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjom]),4);
  else if (strcasecmp(str,"BPJPB_1")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjpb]),0);
  else if (strcasecmp(str,"BPJPB_2")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjpb]),1);
  else if (strcasecmp(str,"BPJPB_3")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjpb]),2);
  else if (strcasecmp(str,"BPJPB_4")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjpb]),3);
  else if (strcasecmp(str,"BPJPB_5")==0)
    readValue(psr,str,fin,&(psr->param[param_bpjpb]),4);
  else if (strcasecmp(str,"NTOA")==0)
    {
      char str[1000];
      fgets(str,1000,fin);
      //      parameter dummy[10];
      //      readValue(psr,str,fin,&dummy,0);
    }
  /* Other allowed parameters that are unused */
  else if (str[0]=='C') /* Comment line */
    fgets(str,1000,fin);
  else if (str[0]=='I') /* Comment line */
    fgets(str,1000,fin);
  else 
    {
      displayMsg(1,"MISC1","Unknown parameter in par file: ",str,psr->noWarnings);
    }
}

void checkAllSet(pulsar *psr,parameter elong,parameter elat,char *filename)
{
  /* Check if we have read a position epoch */
  if (psr->param[param_pepoch].paramSet[0] == 0)
    {
      printf("ERROR [PAR1]: Have not set a period epoch in %s\n",filename);
      exit(1);
    }
  if (psr->param[param_waveepoch].paramSet[0] == 0)
    {
      copyParam(psr->param[param_pepoch],&(psr->param[param_waveepoch]));
      strcpy(psr->param[param_waveepoch].label[0],"WAVEEPOCH (MJD)");
      strcpy(psr->param[param_waveepoch].shortlabel[0],"WAVEEPOCH"); 
   }
  if (psr->param[param_posepoch].paramSet[0] == 0)
    {
      displayMsg(1,"PAR1","Have not set a position epoch. The period epoch will be used instead.",filename,psr->noWarnings);
       copyParam(psr->param[param_pepoch],&(psr->param[param_posepoch]));
      strcpy(psr->param[param_posepoch].label[0],"POSEPOCH (MJD)");
      strcpy(psr->param[param_posepoch].shortlabel[0],"POSEPOCH");
    }
  if (psr->param[param_dmepoch].paramSet[0] == 0)
    {
      displayMsg(1,"PAR2","Have not set a DM epoch. The period epoch will be used instead.",filename,psr->noWarnings);
      copyParam(psr->param[param_pepoch],&(psr->param[param_dmepoch]));
      strcpy(psr->param[param_dmepoch].label[0],"DMEPOCH (MJD)");
      strcpy(psr->param[param_dmepoch].shortlabel[0],"DMEPOCH");
    }
  
  if (elat.paramSet[0]==1 && elong.paramSet[0]==1)
    {
      /* Just record the elat and elong into dec and ra respectively */
      /* set eclcoord flag to 1                                      */
      psr->param[param_raj].val[0]     = elong.val[0]*M_PI/180.0;
      psr->param[param_decj].val[0]    = elat.val[0]*M_PI/180.0;
      psr->param[param_raj].prefit[0]  = elong.val[0]*M_PI/180.0;
      psr->param[param_decj].prefit[0] = elat.val[0]*M_PI/180.0;
      psr->eclCoord = 1;
      
	      
      psr->param[param_raj].paramSet[0]  = 1;
      psr->param[param_decj].paramSet[0] = 1;
      if (elat.fitFlag[0]==1 || elong.fitFlag[0]==1)
	{
	  psr->param[param_raj].fitFlag[0]=1;
	  psr->param[param_decj].fitFlag[0]=1;
	}
    }
  /* correct CLK parameter if necessary */
  if (!strcmp(psr->clock, "UNCORR"))
    {
      strcpy(psr->clock, "TT(TAI)");
      strcpy(psr->clockFromOverride, "UTC");
    }
  if ((psr->clock[0]!='T' || psr->clock[1]!='T'))
    {
      char msg[1000];
      sprintf(msg,"CLK parameter '%s' is not a realization of TT!",psr->clock);
      displayMsg(1,"CLK1",msg,"",psr->noWarnings);
      
      /* try various tempo possibilities */
      char *clk;
      if (!strcmp(psr->clock, "UTC(NIST)"))
	clk = "TT(UTC(NIST))";
      else if (!strcmp(psr->clock, "UTC(BIPM)"))
	clk = "TT(TAI)";
      else if (!strcmp(psr->clock, "PTB"))
	clk = "TT(UTC(PTB))";
      else if (!strcmp(psr->clock, "AT1"))
	clk = "TT(TA(NIST))";
      else /* default to TT(TAI) */
	clk = "TT(TAI)";
      if (psr->noWarnings!=2)
	printf("-> Using %s instead.(You should set CLK to what you really mean!)\n", clk);
      strcpy(psr->clock, clk);
    }
}

/* ******************************************** */
/* readParfile                                 */
/* Author:  G. Hobbs (11 May 2003)              */
/* Purpose: Reads initial pulsar parameters     */
/*          from the .par file                  */
/* Inputs:  filename of parfile                 */
/* Outputs: Fills psr structure                 */
/*                                              */
/* Notes:                                       */
/* Changes:                                     */
/* ******************************************** */
void readParfile(pulsar *psr,char parFile[][MAX_FILELEN],char timFile[][MAX_FILELEN],int npsr)
{
  FILE *fin;
  int nread,p,j;
  char str[1000];
  parameter elong,elat;	
  int noread=0,endit;

  elong.aSize = 1;
  elong.val       = (longdouble *)malloc(elong.aSize*sizeof(longdouble));
  elong.err       = (longdouble *)malloc(elong.aSize*sizeof(longdouble));
  elong.prefit    = (longdouble *)malloc(elong.aSize*sizeof(longdouble));
  elong.prefitErr = (longdouble *)malloc(elong.aSize*sizeof(longdouble));
  elong.fitFlag   = (int *)malloc(elong.aSize*sizeof(int));
  elong.paramSet  = (int *)malloc(elong.aSize*sizeof(int));
  elong.label     = (char **)malloc(elong.aSize*sizeof(char *));
  elong.shortlabel= (char **)malloc(elong.aSize*sizeof(char *));

  for (j=0;j<elong.aSize;j++)
    {
      elong.label[j] = (char *)malloc(sizeof(char)*100);
      elong.shortlabel[j] = (char *)malloc(sizeof(char)*100);
    }
  for (j=0;j<elong.aSize;j++)
    {
      elong.paramSet[j] = 0;
      elong.err[j]      = 0.0;	     
    }

  elat.aSize = 1;
  elat.val       = (longdouble *)malloc(elat.aSize*sizeof(longdouble));
  elat.err       = (longdouble *)malloc(elat.aSize*sizeof(longdouble));
  elat.prefit    = (longdouble *)malloc(elat.aSize*sizeof(longdouble));
  elat.prefitErr = (longdouble *)malloc(elat.aSize*sizeof(longdouble));
  elat.fitFlag   = (int *)malloc(elat.aSize*sizeof(int));
  elat.paramSet  = (int *)malloc(elat.aSize*sizeof(int));
  elat.label     = (char **)malloc(elat.aSize*sizeof(char *));
  elat.shortlabel= (char **)malloc(elat.aSize*sizeof(char *));
  for (j=0;j<elat.aSize;j++)
    {
      elat.label[j] = (char *)malloc(sizeof(char)*100);
      elat.shortlabel[j] = (char *)malloc(sizeof(char)*100);
    }
  for (j=0;j<elat.aSize;j++)
    {
      elat.paramSet[0] = 0;
      elat.err[j]      = 0.0;	     
    }
  
  
  for (p=0;p<npsr;p++)
    {
      elat.paramSet[0] = 0;
      elong.paramSet[0] = 0;
      psr[p].nPhaseJump = 0;
      
      if (!(fin = fopen(parFile[p],"r")))
	{
	  /* Attempt to read fixed format header from tim file */
	  if (!(fin = fopen(timFile[p],"r")))
	    {
	      printf("ERROR [FILE3]: Unable to open parfile %s for pulsar %d\n",parFile[p],p);
	      exit(1);
	    }
	  fgets(str,1000,fin); removeCR(str);
	  if (strcasecmp(str,"HEAD")==0)
	    {
	      noread=-1;
	      printf("Reading from .tpo format file\n");
	    }
	  else
	    {
	      printf("Attempting to read parameters from fixed format header in %s\n",timFile[p]);
	      if (str[0]=='1' && (str[1]=='1' || str[1]=='0') && (str[2]=='1' || str[2]=='0'))
		{
		  if (str[1]=='1') psr[p].param[param_f].fitFlag[0] = 1;
		  if (str[2]=='1') psr[p].param[param_f].fitFlag[1] = 1;
		  if (str[3]=='1') psr[p].param[param_f].fitFlag[2] = 1;
		  if (str[4]=='1') psr[p].param[param_raj].fitFlag[0] = 1;
		  if (str[5]=='1') psr[p].param[param_decj].fitFlag[0] = 1;
		  if (str[6]=='1') psr[p].param[param_pmra].fitFlag[0] = 1;
		  if (str[7]=='1') psr[p].param[param_pmdec].fitFlag[0] = 1;
		  if (str[8]=='1') psr[p].param[param_a1].fitFlag[0] = 1;
		  if (str[9]=='1') psr[p].param[param_ecc].fitFlag[0] = 1;
		  if (str[10]=='1') psr[p].param[param_t0].fitFlag[0] = 1;
		  if (str[11]=='1') psr[p].param[param_pb].fitFlag[0] = 1;
		  if (str[12]=='1') psr[p].param[param_om].fitFlag[0] = 1;
		  if (str[13]=='1') psr[p].param[param_omdot].fitFlag[0] = 1;
		  if (str[14]=='1') psr[p].param[param_gamma].fitFlag[0] = 1;
		  if (str[15]=='1') psr[p].param[param_dm].fitFlag[0] = 1;
		  if (str[16]=='1') psr[p].param[param_px].fitFlag[0] = 1;
		  if (str[17]=='1') psr[p].param[param_pbdot].fitFlag[0] = 1;
		  /*	      if (str[18]=='1') psr[p].param[param_m1].fitFlag = 1; */
		  if (str[19]=='1') psr[p].param[param_sini].fitFlag[0] = 1;
		  if (str[20]=='1') psr[p].param[param_mtot].fitFlag[0] = 1;
		  if (str[21]=='1') psr[p].param[param_m2].fitFlag[0] = 1;
		  if (str[22]=='1') psr[p].param[param_dtheta].fitFlag[0] = 1;
		  if (str[22]=='1') psr[p].param[param_dtheta].fitFlag[0] = 1;
		  if (str[26]=='1') strcpy(psr[p].binaryModel,"BT");
		  if (str[26]=='2') strcpy(psr[p].binaryModel,"EH");
		  if (str[26]=='3') strcpy(psr[p].binaryModel,"DD");
		  if (str[26]=='4') strcpy(psr[p].binaryModel,"DDGR");
		  if (str[26]=='5') strcpy(psr[p].binaryModel,"H88");
		  if (str[26]=='6') strcpy(psr[p].binaryModel,"BT+");	      
		  if (str[26]=='7') strcpy(psr[p].binaryModel,"DDT");
		  if (str[26]=='8') strcpy(psr[p].binaryModel,"DD+");
		  if (str[26]=='9') strcpy(psr[p].binaryModel,"BT2P");
		  strcpy(psr[p].clock,"TT(TAI)");
		  strcpy(psr[p].JPL_EPHEMERIS,getenv(TEMPO2_ENVIRON));
		  strcat(psr[p].JPL_EPHEMERIS,"/ephemeris/DE200.1950.2050");
		  /* Now get parameter values */
		  fgets(str,1000,fin);	      	     
		  strcpy(psr[p].name,str); psr[p].name[12]='\0';
		  getValue(str,21,40,&psr[p],param_raj,0);
		  getValue(str,41,60,&psr[p],param_decj,0);
		  getValue(str,61,70,&psr[p],param_pmra,0);
		  getValue(str,71,80,&psr[p],param_pmdec,0);
		  
		  fgets(str,1000,fin);	      	     
		  getValue(str,2,20,&psr[p],param_f,0);
		  getValue(str,21,40,&psr[p],param_f,1);
		  getValue(str,41,60,&psr[p],param_pepoch,0);
		  getValue(str,61,70,&psr[p],param_f,2);
		  getValue(str,71,80,&psr[p],param_px,0);
		  
		  fgets(str,1000,fin);	      	     
		  getValue(str,9,20,&psr[p],param_dm,0);
		  psr[p].fixedFormat = 4;
		  
		  
		  psr[p].param[param_posepoch].val[0] = psr[p].param[param_pepoch].val[0];
		  psr[p].param[param_posepoch].prefit[0] = psr[p].param[param_pepoch].prefit[0];
		  psr[p].param[param_posepoch].paramSet[0] = 1;
		  
		  /* Now should read in binary parameters */

		  noread=1;
		}
	      else
		{
		  printf("Unable to read pulsar timing model\n");
		  exit(1);
		}
	    }
	}
    
      /* Set up some defaults             */
      /* Should use environment variables ------------- */
      if (noread!=1)
	{
	  if (psr[p].fixedFormat==0)
	    {
	      if (setupParameterFileDefaults(&psr[p]) < 0)
		exit(1);
	      endit=0;
	      while (!feof(fin) && endit==0)
		{
		  nread = fscanf(fin,"%s",str);
		  if (nread==1)
		    {
		      if (noread==-1 && strcasecmp(str,"TOAS")==0)
			endit=1;
		      checkLine(&psr[p],str,fin,&elong,&elat); 
		    }
		}
	      checkAllSet(&psr[p],elong,elat,parFile[p]);
	    }
	}
      fclose(fin);
    }
  /* Free-up the memory */
  for(j=0;j<elong.aSize;j++){
    free(elong.label[j]);
    free(elong.shortlabel[j]);
  }

  free(elong.val);  free(elong.err);  free(elong.prefit);  free(elong.prefitErr);
  free(elong.fitFlag);  free(elong.paramSet);  free(elong.label);  free(elong.shortlabel);
  
  free(elat.val);  free(elat.err);  free(elat.prefit);  free(elat.prefitErr);
  free(elat.fitFlag);  free(elat.paramSet);  free(elat.label);  free(elat.shortlabel); 
}


 
/* ******************************************** */
/* readValue                                    */
/* Author:  G. Hobbs (13 June 2003)             */
/* Purpose: Reads a parameter (value, error,    */
/*          and flag from the .par file)        */
/* Inputs:  File pointer to .par file           */
/*          The parameter to fill               */
/* Outputs: Fills the val, err and flag of      */
/*          the parameter                       */
/*                                              */
/* Notes:   Many .par file have the Fortran     */
/*          method for giving double            */
/*          exponential (e.g. 2.0D+03)          */
/*          This function changes these to e.g. */
/*          2.0e+03                             */
/*                                              */
/*          Errors on RAJ and DECJ in the .par  */
/*          are assumed to be in radians.  This */
/*          is probably not true, but TEMPO     */
/*          does not use the errors anyway!     */
/*                                              */
/* Changes:                                     */
/* ******************************************** */

int readValue(pulsar *psr,char *pmtr,FILE *fin,parameter *parameter,int arr)
{
  char str1[1000],str2[1000],str3[1000],str[1000];
  int i,nread;

  fgets(str,1000,fin);
  nread = sscanf(str,"%s %s %s",str1,str2,str3);

  if (strcasecmp(pmtr,"DSHK")==0)
    {
      if (strcasecmp(str1,"PX")==0) /* Link to parallax distance */
	{
	  parameter->linkTo[(parameter->nLinkTo)++] = param_px;
	  psr->param[param_px].linkFrom[(psr->param[param_px].nLinkFrom)++] = param_dshk;
	  //I think the next two lines are superfluous since this is now all done in preProcess.C. 
	  psr->param[param_dshk].val[0] = psr->param[param_px].val[0];
	  psr->param[param_dshk].prefit[0] = getParameterValue(psr,param_dshk,0);
	  psr->param[param_dshk].paramSet[0]=1;
	  return 0;
	}
    }
  else if (strcasecmp(pmtr,"SINI")==0)
    {
      if(strcasecmp(str1,"KIN")==0)
	{
	  parameter->linkTo[(parameter->nLinkTo)++] = param_kin;
	  psr->param[param_kin].linkFrom[(psr->param[param_kin].nLinkFrom)++] = param_sini;
	  psr->param[param_sini].paramSet[0]=1;
	  return 0;
	}
    }
  else if (strcasecmp(pmtr,"D_AOP")==0){
    if(strcasecmp(str1,"PBDOT")==0){
      parameter->linkTo[(parameter->nLinkTo)++] = param_pbdot;
      psr->param[param_pbdot].linkFrom[(psr->param[param_pbdot].nLinkFrom)++]=param_daop;

      psr->param[param_daop].paramSet[0]=1;
      return 0;
    }
  }

  /* Change e.g. D+02 to e+02 as expected in C */
  for (i=0;i<(int)strlen(str1);i++)
    {
      if (str1[i]=='D' || str1[i]=='d')
	str1[i]='e';
    }
  /* If nread = 3, then do the same for the uncertainty */
  if (nread==3)
    {
      for (i=0;i<(int)strlen(str3);i++)
	{
	  if (str3[i]=='D' || str3[i]=='d')
	    str3[i]='e';
	}
    }
  else if (nread==2)
    {
      /* If nread = 2, then the uncertainty may be the second value */
      for (i=0;i<(int)strlen(str2);i++)
	{
	  if (str2[i]=='D' || str2[i]=='d')
	    str2[i]='e';
	}
    }
  
  if (strcasecmp(pmtr,"RAJ")==0) /* For RAJ, convert from h:m:s to radians */
    {
      strcpy(psr->rajStrPre,str1);
      strcpy(psr->rajStrPost,str1);
      parameter->val[arr] = turn_deg(hms_turn(str1))*M_PI/180.0;
      if (parameter->val[arr]<0) 
	{
	  char retstr[100];
	  printf("ERROR: have negative RAJ: %.14lf\n",(double)parameter->val[arr]);
	  return -1;
	  parameter->val[arr]=2.0*M_PI+parameter->val[arr];
	  turn_hms(parameter->val[arr]/(2.0*M_PI), retstr);
	  strcpy(psr->rajStrPre,retstr);
	  strcpy(psr->rajStrPost,retstr);
	}
    }
  else if (strcasecmp(pmtr,"DECJ")==0) /* For decj convert from d:m:s to radians */
    {
      strcpy(psr->decjStrPre,str1);
      strcpy(psr->decjStrPost,str1);
      parameter->val[arr] = turn_deg(dms_turn(str1))*M_PI/180.0;
    }
  else
    parameter->val[arr] = parse_longdouble(str1);
  parameter->fitFlag[arr] = 0;
  parameter->err[arr]     = 0.0;

  if (nread==2) /* Have uncertainty as second value */
    {
      if (strcasecmp(str2,"1")==0 || strcasecmp(str2,"0")==0 || strcasecmp(str2,"2")==0) /* Have fit flag not error */
	{
	  sscanf(str2,"%d",&(parameter->fitFlag[arr]));
	  parameter->err[arr] = 0.0;
	}
      else
	parameter->prefitErr[arr] = parse_longdouble(str2);
    }
  else if (nread==3)
    {
      sscanf(str2,"%d",&(parameter->fitFlag[arr]));
      parameter->prefitErr[arr] = parse_longdouble(str3);
    }
  if (strcasecmp(pmtr,"P0")==0 || strcasecmp(pmtr,"P")==0) /* Must invert frequency */
    {
      parameter->val[arr] = 1.0/parameter->val[arr];
      parameter->err[arr] = parameter->err[arr]*parameter->val[arr]*parameter->val[arr];
    }
  else if (strcasecmp(pmtr,"P1")==0 || strcasecmp(pmtr,"Pdot")==0) /* Must invert frequency derivative */
    {
      if (psr->param[param_f].paramSet[0] == 0)
	{
	  printf("Error: Tempo2 used frequencies and frequency derivatives (instead of P0 and P1 etc.).\n");
	  printf("       You must set P0 before using P1 (or just use F0 and F1 in any order).\n");
	  return -1;
	}
      parameter->val[arr] = -1.0*psr->param[param_f].val[0]*psr->param[param_f].val[0]*parameter->val[arr]*1e-15;
      parameter->err[arr] =
	(2.0/psr->param[param_f].val[0]*parameter->val[arr]*
		       parameter->val[arr]*psr->param[param_f].err[0]) *
	(2.0/psr->param[param_f].val[0]*parameter->val[arr]*
		       parameter->val[arr]*psr->param[param_f].err[0]) +
	(psr->param[param_f].val[0]*psr->param[param_f].val[0]*parameter->err[arr]) * 
	(psr->param[param_f].val[0]*psr->param[param_f].val[0]*parameter->err[arr]);
    }
  /* Indicate that this parameter has been set */
  parameter->paramSet[arr] = 1;
  parameter->prefit[arr] = parameter->val[arr];
}



int getValue(char *str,int v1,int v2,pulsar *psr,int label,int arr)
{
  char segment[1000];
  char temp[1000],t1[1000],t2[1000],t3[1000];

  strcpy(segment,str+v1-1);
  segment[v2-v1+1]='\0';
  if (label == param_raj) /* Must convert to ':' form and to radians */
    {
      strcpy(t1,str+v1);   t1[2]='\0'; 
      strcpy(t2,str+v1+2); t2[2]='\0'; 
      strcpy(t3,str+v1+4); t3[v2-v1-4]='\0'; 
      sprintf(temp,"%s:%s:%s",t1,t2,t3);
      temp[v2-v1+1+2]='\0';
      strcpy(psr->rajStrPre,temp);
      psr->param[label].val[arr] = turn_deg(hms_turn(temp))*M_PI/180.0;	      
      psr->param[label].paramSet[arr] = 1;
      if (psr->param[label].val[arr]<0) 
	{
	  printf("ERROR: have negative RAJ: %.14lf\n",(double)psr->param[label].val[arr]);
	  return -1;
	}
      psr->param[label].prefit[arr] = psr->param[label].val[arr];   
    }
  else if (label == param_decj) /* Must convert to ':' form and to radians */
    {
      strcpy(t1,str+v1);   t1[2]='\0'; 
      strcpy(t2,str+v1+2); t2[2]='\0'; 
      strcpy(t3,str+v1+4); t3[v2-v1-4]='\0'; 
      psr->param[label].paramSet[arr] = 1;
      sprintf(temp,"%s:%s:%s",t1,t2,t3);
      temp[v2-v1+1+2]='\0';
      strcpy(psr->decjStrPre,temp);
      psr->param[label].val[arr] = turn_deg(dms_turn(temp))*M_PI/180.0;	      
      psr->param[label].prefit[arr] = psr->param[label].val[arr];   
    }
  else
    {
      if (sscanf(segment,"%Lf",&(psr->param[label].val[arr]))==1)
	{
	  psr->param[label].paramSet[arr]=1;
	  if (label==param_f && arr==0) /* Actually read in P0 */
	    psr->param[label].val[arr] = 1.0/psr->param[label].val[arr];
	  else if (label==param_f && arr==1) /* Actually read in P1 */
	    psr->param[label].val[arr] = -1.0*psr->param[param_f].val[0]*psr->param[param_f].val[0]*psr->param[label].val[arr]*1e-15;
	  else if (label==param_f && arr==2) /* Actually read in P2 */
	    psr->param[label].val[arr] = (2.0*psr->param[param_f].val[1]*psr->param[param_f].val[1]/psr->param[param_f].val[0] -
				     psr->param[label].val[arr]*psr->param[param_f].val[0]*psr->param[param_f].val[0])*1.0e-30;
	  else if (label==param_pepoch) /* Convert to MJD from JD */
	    psr->param[label].val[arr] -= 2400000.5;
	  psr->param[label].prefit[arr] = psr->param[label].val[arr];
	}
    }
}

/* Removes newline at end of string */
void removeCR(char *str)
{
  if (str[strlen(str)-1]=='\n') str[strlen(str)-1] = '\0';
}
