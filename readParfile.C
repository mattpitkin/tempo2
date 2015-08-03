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

  checkAllSet(p,elong,elat,(char *)"");

  return 0;
}

void readParfileGlobal(pulsar *psr,int npsr,char tpar[MAX_STRLEN][MAX_FILELEN],
		       char ttim[MAX_STRLEN][MAX_FILELEN])
{
  FILE *fin;
  char str[MAX_STRLEN];
  parameter elong,elat;
  int nread,p;

  for (p=0;p<npsr;p++)
    {
      if (!(fin = fopen(tpar[0],"r")))
	{
	  printf("ERROR: unable to open file >%s<\n",tpar[0]);
	  exit(1);
	}
      while (!feof(fin))
	{
	  // Read in a line from the parameter file
	  nread = fscanf(fin,"%s",str);
	  if (nread==1)
	    checkLine(psr+p,str,fin,&elong,&elat);
	}
      fclose(fin);
    }
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
  else if (strcasecmp(str,"EOP_FILE")==0)
    {
      fscanf(fin,"%s",psr->eopc04_file);
      printf("WARNING: All pulsars will use the EOPC04 file: %s\n",psr->eopc04_file);
    }
  else if (strcasecmp(str,"TRES")==0)
    readValue(psr,str,fin,&(psr->param[param_tres]),0);
  else if (strcasecmp(str,"MODE")==0 || strcasecmp(str,"WEIGHT")==0) /* Fitting mode */
    fscanf(fin,"%d",&(psr->fitMode));
  else if (strcasecmp(str,"ROBUST")==0) /* Robust Fitting mode */
    fscanf(fin,"%d",&(psr->robust));
  else if (strcasecmp(str,"NOTRACK")==0)  /* TEMPO2 uses automatic tracking */
    psr->param[param_track].paramSet[0]=0;
  else if (strcasecmp(str,"WHITE_NOISE_MODEL")==0) // Use a white noise model file
    fscanf(fin,"%s",psr->whiteNoiseModelFile);
  else if (strcasecmp(str,"TRACK")==0)  /* TEMPO2 uses automatic tracking */
    readValue(psr,(char *)"TRACK",fin,&(psr->param[param_track]),0);
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
      // printf("Setting correctTroposphere %d\n",psr->correctTroposphere);
    }
  else if (strcasecmp(str,"UNITS")==0)
    {
      char unit[1000];
      psr->setUnits=1;
      fscanf(fin,"%s", unit);
      if (strcasecmp(unit,"TDB")==0) psr->units = TDB_UNITS;
      else if (strcasecmp(unit,"TCB")==0) psr->units = SI_UNITS;
      else if (strcasecmp(unit,"SI")==0) psr->units = SI_UNITS;
    }
  else if (strcasecmp(str,"NE1AU")==0 || strcasecmp(str,"NE_SW")==0 ||
          strcasecmp(str,"SOLARN0")==0)
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
      if (strcasecmp(unit,"IAU2000B")==0) 
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
  else if (strcasecmp(str,"EPHEM_FILE")==0)
    {
      fscanf(fin,"%s",psr->ephemeris);
      strcpy(psr->JPL_EPHEMERIS,psr->ephemeris);
    }
  else if (strcasecmp(str,"EPH_FILE")==0)
    {
      char temp[1024];
      fscanf(fin,"%s",temp);
      sprintf(psr->ephemeris,"%s/ephemeris/%s",getenv("TEMPO2"),temp);
      psr->useCalceph = 1;
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
  else if (strcasecmp(str,"STEL_CLK_OFFS")==0)  /* Set clock offsets */
    readValue(psr,str,fin,&(psr->param[param_clk_offs]),0);
  else if (strcasecmp(str,"STEL_DX")==0)  /* Set interpolation function for telescope position offset*/
    {
      psr->setTelVelX=0;
      readValue(psr,str,fin,&(psr->param[param_tel_dx]),0);
    }
  else if (strcasecmp(str,"STEL_DY")==0)  /* Set interpolation function for telescope position offset*/
    {
      psr->setTelVelY=1;
      readValue(psr,str,fin,&(psr->param[param_tel_dy]),0);
    }
  else if (strcasecmp(str,"STEL_DZ")==0)  /* Set interpolation function for telescope position offset*/
    {
      psr->setTelVelZ=1;
      readValue(psr,str,fin,&(psr->param[param_tel_dz]),0);
    }
  else if (strcasecmp(str,"TEL_VX")==0)  /* Set telescope velocity in X */
      readValue(psr,str,fin,&(psr->param[param_tel_vx]),0);
  else if (strcasecmp(str,"TEL_VY")==0)  /* Set telescope velocity in Y */
    readValue(psr,str,fin,&(psr->param[param_tel_vy]),0);
  else if (strcasecmp(str,"TEL_VZ")==0)  /* Set telescope velocity in Z */
    readValue(psr,str,fin,&(psr->param[param_tel_vz]),0);
  else if (strcasecmp(str,"TEL_X0")==0)  /* Set telescope position in X */
      readValue(psr,str,fin,&(psr->param[param_tel_x0]),0);
  else if (strcasecmp(str,"TEL_Y0")==0)  /* Set telescope position in Y */
    readValue(psr,str,fin,&(psr->param[param_tel_y0]),0);
  else if (strcasecmp(str,"TEL_Z0")==0)  /* Set telescope position in Z */
    readValue(psr,str,fin,&(psr->param[param_tel_z0]),0);
  else if (strcasecmp(str,"SQIFUNC_p")==0)  /* Set quad interpolation function for plus*/
      readValue(psr,str,fin,&(psr->param[param_quad_ifunc_p]),0);
  else if (strcasecmp(str,"SQIFUNC_c")==0)  /* Set quad interpolation function for cross*/
    readValue(psr,str,fin,&(psr->param[param_quad_ifunc_c]),0);
  else if (strcasecmp(str,"PEPOCH")==0)    /* Period Epoch */
    readValue(psr,str,fin,&(psr->param[param_pepoch]),0);
  else if (strcasecmp(str,"TELEPOCH")==0) /* Epoch of telescope position */
    readValue(psr,str,fin,&(psr->param[param_telEpoch]),0);
  else if (strcasecmp(str,"EPHVER")==0)    /* Ephemeris version */
    readValue(psr,str,fin,&(psr->param[param_ephver]),0);
  else if (strcasecmp(str,"DMEPOCH")==0)    /* DM Epoch */
    readValue(psr,str,fin,&(psr->param[param_dmepoch]),0);
  else if (strcasecmp(str,"RAJ")==0 || strcasecmp(str,"RA")==0)       /* Right ascension */
    readValue(psr,(char *)"RAJ",fin,&(psr->param[param_raj]),0);
  else if (strcasecmp(str,"DECJ")==0 || strcasecmp(str,"DEC")==0)      /* Declination */
    readValue(psr,(char *)"DECJ",fin,&(psr->param[param_decj]),0);
  else if (strcasecmp(str,"ELONG")==0 || strcasecmp(str,"LAMBDA")==0)
    readValue(psr,(char *)"ELONG",fin,elong,0);
  else if (strcasecmp(str,"ELAT")==0 || strcasecmp(str,"BETA")==0)
    readValue(psr,(char *)"ELAT",fin,elat,0);
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
  else if (strncasecmp(str,"FD",2)==0 && (str[2]!=' '))
    {
      int fval;
      if (sscanf(str+2,"%d",&fval)==1)
        {
          if (fval<psr->param[param_f].aSize) 
            readValue(psr,str,fin,&(psr->param[param_fd]),fval-1);
        }
    }
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
  else if ((str[0]=='F' || str[0]=='f') && (str[1]!='B')) /* Read higher frequency derivatives */
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
  else if (strncasecmp(str,"DMX_",4)==0)
    {
      int dmxidx;
      if (sscanf(str+4,"%d",&dmxidx)==1)
        {
            //printf("got dmxidx=%d\n", dmxidx);
          dmxidx--;
          if (dmxidx<psr->param[param_dmx].aSize)
            readValue(psr,str,fin,&(psr->param[param_dmx]),dmxidx);
          if (psr->ndmx < dmxidx+1) psr->ndmx = dmxidx + 1;
        }
    }
  else if (strncasecmp(str,"DMXR1_",6)==0)
    {
      int dmxidx;
      if (sscanf(str+6,"%d",&dmxidx)==1)
        {
          dmxidx--;
          if (dmxidx<psr->param[param_dmxr1].aSize)
            readValue(psr,str,fin,&(psr->param[param_dmxr1]),dmxidx);
        }
    }
  else if (strncasecmp(str,"DMXR2_",6)==0)
    {
      int dmxidx;
      if (sscanf(str+6,"%d",&dmxidx)==1)
        {
          dmxidx--;
          if (dmxidx<psr->param[param_dmxr2].aSize)
            readValue(psr,str,fin,&(psr->param[param_dmxr2]),dmxidx);
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
  else if (strcasecmp(str,"TELX")==0)
    readValue(psr,str,fin,&(psr->param[param_telx]),0);
  else if ((str[0]=='T' || str[0]=='t') &&  /* Higher DM derivatives */
	   (str[1]=='E' || str[1]=='e') && 
	   (str[2]=='L' || str[2]=='l') &&
	   (str[3]=='X' || str[3]=='x') && isdigit(str[4]))
    {
      int dval;
      if (sscanf(str+4,"%d",&dval)==1)
	{
	  if (dval<psr->param[param_telx].aSize)
	    readValue(psr,str,fin,&(psr->param[param_telx]),dval);
	}
    }
  else if (strcasecmp(str,"TELY")==0)
    readValue(psr,str,fin,&(psr->param[param_tely]),0);
  else if ((str[0]=='T' || str[0]=='t') &&  /* Higher DM derivatives */
	   (str[1]=='E' || str[1]=='e') && 
	   (str[2]=='L' || str[2]=='l') &&
	   (str[3]=='Y' || str[3]=='y') && isdigit(str[4]))
    {
      int dval;
      if (sscanf(str+4,"%d",&dval)==1)
	{
	  if (dval<psr->param[param_tely].aSize)
	    readValue(psr,str,fin,&(psr->param[param_tely]),dval);
	}
    }
  else if (strcasecmp(str,"TELZ")==0)
    readValue(psr,str,fin,&(psr->param[param_telz]),0);
  else if ((str[0]=='T' || str[0]=='t') &&  /* Higher DM derivatives */
	   (str[1]=='E' || str[1]=='e') && 
	   (str[2]=='L' || str[2]=='l') &&
	   (str[3]=='Z' || str[3]=='z') && isdigit(str[4]))
    {
      int dval;
      if (sscanf(str+4,"%d",&dval)==1)
	{
	  if (dval<psr->param[param_telz].aSize)
	    readValue(psr,str,fin,&(psr->param[param_telz]),dval);
	  printf("Setting %d %d\n",dval,psr->param[param_telz].aSize);
	}
    }
  else if (strcasecmp(str,"PX")==0) /* Parallax */
    readValue(psr,str,fin,&(psr->param[param_px]),0);
  else if (strcasecmp(str,"DM_S1YR")==0) /* DM sinusoid */
    readValue(psr,str,fin,&(psr->param[param_dm_sin1yr]),0);
  else if (strcasecmp(str,"DM_C1YR")==0) /* DM sinusoid */
    readValue(psr,str,fin,&(psr->param[param_dm_cos1yr]),0);
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
    {
      readValue(psr,str,fin,&(psr->param[param_wave_om]),0);
      psr->waveScale = 0;
    }
  else if (strcasecmp(str,"WAVE_SCALE")==0)
    fscanf(fin,"%lf",&psr->waveScale);
  else if (strstr(str,"WAVE")!=NULL || strstr(str,"wave")!=NULL)
    {
      int number;
      /* Obtain parameter number */
      sscanf(str+4,"%d",&number);
      fscanf(fin,"%lf %lf",&psr->wave_sine[number-1],&psr->wave_cos[number-1]);
      if (psr->nWhite < number) psr->nWhite = number;
    }
   else if (strcasecmp(str,"WAVDM_OM")==0) /* Fundamental frequency */
    {
      readValue(psr,str,fin,&(psr->param[param_wave_dm]),0);
      //psr->waveScale = 0;
    }
  
  else if (strstr(str,"WAVDM")!=NULL || strstr(str,"wavedm")!=NULL)
    {
      int number;
      /* Obtain parameter number */
      sscanf(str+5,"%d",&number);
      fscanf(fin,"%lf %lf",&psr->wave_sine_dm[number-1],&psr->wave_cos_dm[number-1]);
      if (psr->nWhite_dm < number) psr->nWhite_dm = number;
    }
  


  /* ------------------- */
  /* Quad polar function */
  /* ------------------- */ 
  else if (strcasecmp(str,"QUAD_OM")==0) /* Fundamental frequency */
    readValue(psr,str,fin,&(psr->param[param_quad_om]),0);
  else if (strcasecmp(str,"QUAD_POS")==0) // Position of quadrupole
    fscanf(fin,"%lf %lf",&psr->quadRA,&psr->quadDEC);
  else if (strcasecmp(str,"QUAD_EPOCH")==0) // Epoch for quad function
    fscanf(fin,"%lf",&psr->quadEpoch);
  else if ((strstr(str,"QUAD")!=NULL || strstr(str,"quad")!=NULL)
      && (strcasecmp(str,"T2EQUAD")!=0)) // Avoid misinterpreting T2EQUAD
    {
      int number;
      /* Obtain parameter number */
      sscanf(str+4,"%d",&number);
      fscanf(fin,"%lf %lf %lf %lf",&psr->quad_aplus_r[number-1],
	     &psr->quad_aplus_i[number-1],
	     &psr->quad_across_r[number-1],
	     &psr->quad_across_i[number-1]);
      if (psr->nQuad < number) psr->nQuad = number;
    }
  /*
   * DMMODEL fitting.
   *
   */
  else if (strcasecmp(str,"DMMODEL")==0) 
    {
      readValue(psr,str,fin,&(psr->param[param_dmmodel]),0);
      psr->dmoffsCMnum=0;
      psr->dmoffsDMnum=0;
    }
  //  else if (strstr(str,"DMVAL")!=NULL || strstr(str,"dmval")!=NULL)
  else if (strcasecmp(str,"DMOFF")==0)
    {
      int nDM = psr->dmoffsDMnum;
      int nCM = psr->dmoffsCMnum;
	  double mjd,val,err;
      fscanf(fin,"%lf %lf %lf",&mjd,&val,&err);
	  psr->dmoffsDM_mjd[nDM]=mjd;
	  psr->dmoffsDM[nDM]=val;
      psr->dmoffsDM_error[nDM] = 0;
      psr->dmoffsDM_weight[nDM] = 1;
      psr->dmoffsCM_mjd[nCM]=mjd;
      psr->dmoffsCM[nCM]=0;
      psr->dmoffsCM_error[nCM] = 0;
      psr->dmoffsCM_weight[nCM] = 1;

	  (psr->dmoffsDMnum)++;
	  (psr->dmoffsCMnum)++;
	  if (psr->dmoffsDMnum > MAX_IFUNC){
		 fprintf(stderr,"ERROR: Too many DMMODEL DM values - need to increase MAX_IFUNC\n");
		 exit(1);
      } 
	  if (psr->dmoffsCMnum > MAX_IFUNC){
		 fprintf(stderr,"ERROR: Too many DMMODEL CM values - need to increase MAX_IFUNC\n");
		 exit(1);
      } 
    }
else if (strcasecmp(str,"_DM")==0)
    {
      int nDM = psr->dmoffsDMnum;
	  double mjd,val,err;
      fscanf(fin,"%lf %lf %lf",&mjd,&val,&err);
	  psr->dmoffsDM_mjd[nDM]=mjd;
	  psr->dmoffsDM[nDM]=val;
      psr->dmoffsDM_error[nDM] = 0;
      psr->dmoffsDM_weight[nDM] = 1;

	  (psr->dmoffsDMnum)++;
	  if (psr->dmoffsDMnum > MAX_IFUNC){
		 fprintf(stderr,"ERROR: Too many DMMODEL DM values - need to increase MAX_IFUNC\n");
		 exit(1);
      } 
    }
  else if (strcasecmp(str,"_CM")==0)
    {
      int nCM = psr->dmoffsCMnum;
	  double mjd,val,err;
	  err=0;
      fscanf(fin,"%lf %lf %lf",&mjd,&val,&err);
	  psr->dmoffsCM_mjd[nCM]=mjd;
	  psr->dmoffsCM[nCM]=val;
      psr->dmoffsCM_error[nCM] = err;
      psr->dmoffsCM_weight[nCM] = 1;

	  (psr->dmoffsCMnum)++;
	  if (psr->dmoffsCMnum > MAX_IFUNC){
		 fprintf(stderr,"ERROR: Too many CMMODEL CM values - need to increase MAX_IFUNC\n");
		 exit(1);
      } 
    }


  /*
   * Specify fitting constraints
   */
  else if (strcasecmp(str,"CONSTRAIN")==0){
      char cname[1024];
      fscanf(fin, "%s",cname);

      if((strcasecmp(cname,"AUTO")==0)){
		 psr->auto_constraints=1;
	  }
      /*
       * Constraints for DMMODEL.
       * The DMMODEL constraint affects 4 constraints.
       */
      if((strcasecmp(cname,"DMMODEL_X")==0)){
	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_mean;
	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_0;
      }
      if((strcasecmp(cname,"DMMODEL_DM1")==0)){
	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_dm1;
      }
      if((strcasecmp(cname,"DMMODEL_CUBIC")==0)){
	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_3;
      }
      if((strcasecmp(cname,"DMMODEL_PX")==0)){
	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_px;
      }

      if((strcasecmp(cname,"DMMODEL")==0)  || (strcasecmp(cname,"DMMODEL_OLD")==0)
	 || strcasecmp(cname,"DMMODEL_DM1")==0){
	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_mean;
	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_0;
	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_1;
	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_2;
      }
      if((strcasecmp(cname,"DMMODEL")==0)  || (strcasecmp(cname,"DMMODEL_YEAR")==0)
	 || strcasecmp(cname,"DMMODEL_DM1")==0){
	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_year_sin;
	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_year_cos;
	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_year_xsin;
	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_year_xcos;
	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_year_cos2;
	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_year_sin2;
//	      psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_px;
      }

      if(strcasecmp(cname,"IFUNC_ONLYF0F1")==0){
	      psr->constraints[psr->nconstraints++] = constraint_ifunc_0;
	      psr->constraints[psr->nconstraints++] = constraint_ifunc_1;
	      psr->constraints[psr->nconstraints++] = constraint_ifunc_2;
      }

      if(strcasecmp(cname,"IFUNC_ONLYPHI0")==0){
	      psr->constraints[psr->nconstraints++] = constraint_ifunc_0;
      }

      if(strcasecmp(cname,"IFUNC")==0){
	      psr->constraints[psr->nconstraints++] = constraint_ifunc_0;
	      psr->constraints[psr->nconstraints++] = constraint_ifunc_1;
	      psr->constraints[psr->nconstraints++] = constraint_ifunc_2;
	      psr->constraints[psr->nconstraints++] = constraint_ifunc_year_sin;
	      psr->constraints[psr->nconstraints++] = constraint_ifunc_year_cos;
	      psr->constraints[psr->nconstraints++] = constraint_ifunc_year_xsin;
	      psr->constraints[psr->nconstraints++] = constraint_ifunc_year_xcos;
	      psr->constraints[psr->nconstraints++] = constraint_ifunc_year_sin2;
	      psr->constraints[psr->nconstraints++] = constraint_ifunc_year_cos2;

      }
      if(strcasecmp(cname,"TEL_DX")==0){
	      psr->constraints[psr->nconstraints++] = constraint_tel_dx_0;
	      psr->constraints[psr->nconstraints++] = constraint_tel_dx_1;
	      psr->constraints[psr->nconstraints++] = constraint_tel_dx_2;
      }
      if(strcasecmp(cname,"TEL_DY")==0){
	      psr->constraints[psr->nconstraints++] = constraint_tel_dy_0;
	      psr->constraints[psr->nconstraints++] = constraint_tel_dy_1;
	      psr->constraints[psr->nconstraints++] = constraint_tel_dy_2;
      }
      if(strcasecmp(cname,"TEL_DZ")==0){
	      psr->constraints[psr->nconstraints++] = constraint_tel_dz_0;
	      psr->constraints[psr->nconstraints++] = constraint_tel_dz_1;
	      psr->constraints[psr->nconstraints++] = constraint_tel_dz_2;
      }
	  if(strcasecmp(cname,"QIFUNC_p")==0 || strcasecmp(cname,"QIFUNC_p_offset")==0){
		 psr->constraints[psr->nconstraints++] = constraint_quad_ifunc_p_0;
	  }
	  if(strcasecmp(cname,"QIFUNC_p")==0 || strcasecmp(cname,"QIFUNC_p_F0F1")==0){

		 psr->constraints[psr->nconstraints++] = constraint_quad_ifunc_p_1;
		 psr->constraints[psr->nconstraints++] = constraint_quad_ifunc_p_2;
	  }
	  if(strcasecmp(cname,"QIFUNC_p")==0|| strcasecmp(cname,"QIFUNC_p_YEAR")==0){

		 psr->constraints[psr->nconstraints++] = constraint_qifunc_p_year_sin;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_p_year_cos;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_p_year_xsin;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_p_year_xcos;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_p_year_sin2;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_p_year_cos2;
	  }
	  if(strcasecmp(cname,"QIFUNC_c")==0 || strcasecmp(cname,"QIFUNC_c_offset")==0){
		 psr->constraints[psr->nconstraints++] = constraint_quad_ifunc_c_0;
	  }
	  if(strcasecmp(cname,"QIFUNC_c")==0 || strcasecmp(cname,"QIFUNC_c_F0F1")==0){
		 psr->constraints[psr->nconstraints++] = constraint_quad_ifunc_c_1;
		 psr->constraints[psr->nconstraints++] = constraint_quad_ifunc_c_2;
	  }
	  if(strcasecmp(cname,"QIFUNC_c")==0|| strcasecmp(cname,"QIFUNC_c_YEAR")==0){
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_c_year_sin;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_c_year_cos;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_c_year_xsin;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_c_year_xcos;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_c_year_sin2;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_c_year_cos2;
	  }

  }
  /*
   * Single source graviational waves (GWs)
   */
  else if (strcasecmp(str,"CGW_FREQ")==0) 
    readValue(psr,str,fin,&(psr->param[param_cgw]),0);
  else if (strcasecmp(str,"GW_SINGLE")==0) 
	 readValue(psr,str,fin,&(psr->param[param_gwsingle]),0);
  else if (strcasecmp(str,"GW_POSITION")==0)
	 fscanf(fin,"%lf %lf",&psr->gwsrc_ra,&psr->gwsrc_dec);
  else if (strcasecmp(str,"GW_APLUS")==0)
	 fscanf(fin,"%lf %lf",&psr->gwsrc_aplus_r,&psr->gwsrc_aplus_i);
  else if (strcasecmp(str,"GW_ACROSS")==0)
	 fscanf(fin,"%lf %lf",&psr->gwsrc_across_r,&psr->gwsrc_across_i);
  else if (strcasecmp(str,"GW_EPOCH")==0)
	 fscanf(fin,"%lf",&psr->gwsrc_epoch);
  else if (strcasecmp(str,"GW_PSR_DIST")==0)
	 fscanf(fin,"%lf",&psr->gwsrc_psrdist);
  else if (strcasecmp(str,"CGW_H0")==0)
	 fscanf(fin,"%lf",&psr->cgw_h0);
  else if (strcasecmp(str,"CGW_COSINC")==0)
	 fscanf(fin,"%lf",&psr->cgw_cosinc);
  else if (strcasecmp(str,"CGW_ANGPOL")==0)
	 fscanf(fin,"%lf",&psr->cgw_angpol);
  else if (strcasecmp(str,"CGW_MC")==0)
    fscanf(fin,"%lf",&psr->cgw_mc);

  // Eccentric binary source graviational waves -- Vikram Ravi
  else if (strcasecmp(str,"GWECC_AMP")==0) 
    readValue(psr,str,fin,&(psr->param[param_gwecc]),0);
  else if (strcasecmp(str,"GWECC_POSITION")==0)
    fscanf(fin,"%lf %lf",&psr->gwecc_ra,&psr->gwecc_dec);
  else if (strcasecmp(str,"GWECC_E")==0)
    fscanf(fin,"%lf",&psr->gwecc_e);
  else if (strcasecmp(str,"GWECC_MASS")==0)
    fscanf(fin,"%lf %lf",&psr->gwecc_m1,&psr->gwecc_m2);
  else if (strcasecmp(str,"GWECC_INC")==0)
    fscanf(fin,"%lf",&psr->gwecc_inc);
  else if (strcasecmp(str,"GWECC_THETA_NODES")==0)
    fscanf(fin,"%lf",&psr->gwecc_theta_nodes);
  else if (strcasecmp(str,"GWECC_NODES_ORIENTATION")==0)
    fscanf(fin,"%lf",&psr->gwecc_nodes_orientation);
  else if (strcasecmp(str,"GWECC_THETA_0")==0)
    fscanf(fin,"%lf",&psr->gwecc_theta_0);
  else if (strcasecmp(str,"GWECC_ORBITAL_PERIOD")==0)
    fscanf(fin,"%lf",&psr->gwecc_orbital_period);
  else if (strcasecmp(str,"GWECC_Z")==0)
    fscanf(fin,"%lf",&psr->gwecc_redshift);
  else if (strcasecmp(str,"GWECC_DIST")==0)
    fscanf(fin,"%lf",&psr->gwecc_distance);
  else if (strcasecmp(str,"GWECC_EPOCH")==0)
    fscanf(fin,"%lf",&psr->gwecc_epoch);
  else if (strcasecmp(str,"GWECC_PSR_DIST")==0)
    fscanf(fin,"%lf",&psr->gwecc_psrdist);
  else if (strcasecmp(str,"GWECC_PSRTERM")==0)
    {
      fscanf(fin,"%d",&psr->gwecc_pulsarTermOn);
      // 0 - no, 1 - yes, 2 - only
      if (psr->param[param_gwecc].paramSet[0] == 1) 
	{
	  strcpy(psr->param[param_gwecc].label[0],"GWECC_AMP");
	  strcpy(psr->param[param_gwecc].shortlabel[0],"GWECC_AMP");
	}
    }
  //* Gravitational wave memory
  else if (strcasecmp(str,"GWM_AMP")==0)
	 readValue(psr,str,fin,&(psr->param[param_gwm_amp]),0);
  else if (strcasecmp(str,"GWM_A1")==0)
    readValue(psr,str,fin,&(psr->param[param_gwm_amp]),0);
  else if (strcasecmp(str,"GWM_A2")==0)
    readValue(psr,str,fin,&(psr->param[param_gwm_amp]),1);
  else if (strcasecmp(str,"GWM_POSITION")==0)
	 fscanf(fin,"%lf %lf",&psr->gwm_raj,&psr->gwm_decj);
  else if (strcasecmp(str,"GWM_EPOCH")==0)
	 fscanf(fin,"%lf",&psr->gwm_epoch);
  else if (strcasecmp(str,"GWM_PHI")==0)
	 fscanf(fin,"%lf",&psr->gwm_phi);
  else if (strcasecmp(str,"GWM_DPHASE")==0)
	 fscanf(fin,"%lf",&psr->gwm_dphase);
   // Ryan's gw bursts
    else if (strcasecmp(str,"GWB_AMP")==0)
	 readValue(psr,str,fin,&(psr->param[param_gwb_amp]),0);
  else if (strcasecmp(str,"GWB_A1")==0)
    readValue(psr,str,fin,&(psr->param[param_gwb_amp]),0);
  else if (strcasecmp(str,"GWB_A2")==0)
    readValue(psr,str,fin,&(psr->param[param_gwb_amp]),1);
  else if (strcasecmp(str,"GWB_POSITION")==0)
	 fscanf(fin,"%lf %lf",&psr->gwb_raj,&psr->gwb_decj);
  else if (strcasecmp(str,"GWB_EPOCH")==0)
	 fscanf(fin,"%lf",&psr->gwb_epoch);
  else if (strcasecmp(str,"GWB_WIDTH")==0)
	 fscanf(fin,"%lf",&psr->gwb_width);
  


  else if ((strstr(str,"IFUNC")!=NULL || strstr(str,"ifunc")!=NULL)
		&& strstr(str,"QIFUNC")==NULL)
  {
	 int number;
	 /* Obtain parameter number */
	 sscanf(str+5,"%d",&number);

	 fscanf(fin,"%lf %lf %lf",&psr->ifuncT[number-1],&psr->ifuncV[number-1],&psr->ifuncE[number-1]);
	 psr->ifunc_weights[number-1]=1.0;
	 if (psr->ifuncN < number) psr->ifuncN = number;
  }
  else if ((strstr(str,"TEL_CLK_OFFS")!=NULL || strstr(str,"tel_clk_offs")!=NULL) &&
		strcmp(str,"STEL_CLK_OFFS")!=0)
  {
	 int number;
	 /* Obtain parameter number */
	 sscanf(str+12,"%d",&number);
	 printf("Reading clock with %s %d\n",str,number);

	 fscanf(fin,"%lf %lf %lf",&psr->clk_offsT[number-1],&psr->clk_offsV[number-1],&psr->clk_offsE[number-1]);
	 if (psr->clkOffsN < number) psr->clkOffsN = number;
  }
  else if ((strstr(str,"TEL_DX")!=NULL || strstr(str,"tel_dx")!=NULL) && strcmp(str,"STEL_DX")!=0)
  {
	 int number;
	 /* Obtain parameter number */
	 sscanf(str+6,"%d",&number);
	 fscanf(fin,"%lf %lf %lf",&psr->telDX_t[number-1],&psr->telDX_v[number-1],&psr->telDX_e[number-1]);
	 if (psr->nTelDX < number) psr->nTelDX = number;
  }
  else if ((strstr(str,"TEL_VX")!=NULL || strstr(str,"tel_vx")!=NULL))
  {
	 int number;
	 psr->setTelVelX=1;
	 /* Obtain parameter number */
	 sscanf(str+6,"%d",&number);
	 fscanf(fin,"%lf %lf",&psr->telDX_vel[number-1],&psr->telDX_vel_e[number-1]);
  }
  else if ((strstr(str,"TEL_DY")!=NULL || strstr(str,"tel_dy")!=NULL) && strcmp(str,"STEL_DY")!=0)
  {
	 int number;
	 /* Obtain parameter number */
	 sscanf(str+6,"%d",&number);
	 if (number < 1) {
	   printf("ERROR loading TEL_DX. The value (%d) must be > 0\n",number);
	   exit(1);
	 }

	 fscanf(fin,"%lf %lf %lf",&psr->telDY_t[number-1],&psr->telDY_v[number-1],&psr->telDY_e[number-1]);
	 if (psr->nTelDY < number) psr->nTelDY = number;
  }
  else if ((strstr(str,"TEL_VY")!=NULL || strstr(str,"tel_vy")!=NULL))
  {
	 int number;
	 psr->setTelVelY=1;
	 /* Obtain parameter number */
	 sscanf(str+6,"%d",&number);
	 if (number < 1) {
	   printf("ERROR loading TEL_DY. The value (%d) must be > 0\n",number);
	   exit(1);
	 }

	 fscanf(fin,"%lf %lf",&psr->telDY_vel[number-1],&psr->telDY_vel_e[number-1]);
  }
  else if ((strstr(str,"TEL_DZ")!=NULL || strstr(str,"tel_dz")!=NULL) && strcmp(str,"STEL_DZ")!=0)
  {
	 int number;
	 /* Obtain parameter number */
	 sscanf(str+6,"%d",&number);
	 if (number < 1) {
	   printf("ERROR loading TEL_DZ. The value (%d) must be > 0\n",number);
	   exit(1);
	 }
	 fscanf(fin,"%lf %lf %lf",&psr->telDZ_t[number-1],&psr->telDZ_v[number-1],&psr->telDZ_e[number-1]);
	 if (psr->nTelDZ < number) psr->nTelDZ = number;
  }
  else if ((strstr(str,"TEL_VZ")!=NULL || strstr(str,"tel_vz")!=NULL))
  {
	 int number;
	 psr->setTelVelZ=1;
	 /* Obtain parameter number */
	 sscanf(str+6,"%d",&number);
	 fscanf(fin,"%lf %lf",&psr->telDZ_vel[number-1],&psr->telDZ_vel_e[number-1]);
  }
  else if (strstr(str,"QIFUNC_p")!=NULL || strstr(str,"qifunc_p")!=NULL)
  {
	 int number;
	 /* Obtain parameter number */
	 sscanf(str+8,"%d",&number);

	 fscanf(fin,"%lf %lf %lf",&psr->quad_ifuncT_p[number-1],&psr->quad_ifuncV_p[number-1],&psr->quad_ifuncE_p[number-1]);
	 if (psr->quad_ifuncN_p < number) psr->quad_ifuncN_p = number;
  }
  else if (strcasecmp(str,"QIFUNC_POS_p") == 0)
	 fscanf(fin,"%lf %lf",&psr->quad_ifunc_p_RA,&psr->quad_ifunc_p_DEC);
  else if (strcasecmp(str,"QIFUNC_POS_c") == 0)
	 fscanf(fin,"%lf %lf",&psr->quad_ifunc_c_RA,&psr->quad_ifunc_c_DEC);
  else if (strstr(str,"QIFUNC_c")!=NULL || strstr(str,"qifunc_c")!=NULL)
  {
	 int number;
	 /* Obtain parameter number */
	 sscanf(str+8,"%d",&number);

	 fscanf(fin,"%lf %lf %lf",&psr->quad_ifuncT_c[number-1],&psr->quad_ifuncV_c[number-1],&psr->quad_ifuncE_c[number-1]);
	 if (psr->quad_ifuncN_c < number) psr->quad_ifuncN_c = number;
  }
  // Braking index
  else if (strcasecmp(str,"BRAKE")==0)
    readValue(psr,str,fin,&(psr->param[param_brake]),0);

  /* ---------------- */
  /* Phase jumps      */
  /* ---------------- */
  else if (strcasecmp(str,"PHASE")==0)
  {
	 fscanf(fin,"%d %Lf",&psr->phaseJumpDir[psr->nPhaseJump],&psr->phaseJump[psr->nPhaseJump]);
	 psr->phaseJumpID[psr->nPhaseJump]=-1;
	 psr->nPhaseJump++;
  }  

  /* /---------\
     | T2 EFAC |
     \---------/ */
  else if( strcasecmp( str, "T2EFAC") == 0 ) // EFAC for given flag
    { 
      int nefacFlag = psr->nT2efac;
      fscanf( fin, "%s %s %lf", psr->T2efacFlagID[nefacFlag], 
              psr->T2efacFlagVal[nefacFlag], 
              &psr->T2efacVal[nefacFlag] );
      ( psr->nT2efac )++;
    }

  /* /----------\
     | T2 EQUAD |
     \----------/ */
  else if( strcasecmp( str, "T2EQUAD") == 0 ) // EQUAD for given flag
    { 
      int nequadFlag = psr->nT2equad;
      fscanf( fin, "%s %s %lf", psr->T2equadFlagID[nequadFlag], 
              psr->T2equadFlagVal[nequadFlag], 
              &psr->T2equadVal[nequadFlag] );
      ( psr->nT2equad )++;
    }


  /* /---------\
     | TN EF |
     \---------/ */
  else if( strcasecmp( str, "TNEF") == 0 ) // EFAC for given flag
    {
      int nefacFlag = psr->nTNEF;
      fscanf( fin, "%s %s %lf", psr->TNEFFlagID[nefacFlag],
              psr->TNEFFlagVal[nefacFlag],
              &psr->TNEFVal[nefacFlag] );
      ( psr->nTNEF )++;
    }
  else if( strcasecmp( str, "TNGlobalEF") == 0 ) // EFAC for given flag
    {
      fscanf( fin, "%lf", &(psr->TNGlobalEF));
    }


  /* /---------\
     | TN EQ |
     \---------/ */
  else if( strcasecmp( str, "TNEQ") == 0 ) // EFAC for given flag
    {
      int nequadFlag = psr->nTNEQ;
      fscanf( fin, "%s %s %lf", psr->TNEQFlagID[nequadFlag],
              psr->TNEQFlagVal[nequadFlag],
              &psr->TNEQVal[nequadFlag] );
      ( psr->nTNEQ )++;
    }

  else if( strcasecmp( str, "TNGlobalEQ") == 0 ) // EFAC for given flag
    {
      fscanf( fin, "%lf", &(psr->TNGlobalEQ));
    }

  else if( strcasecmp( str, "addTNGlobalEQ") == 0 ) // EFAC for given flag
    {
      fscanf( fin, "%lf", &(psr->addTNGlobalEQ));
    }


  /* /---------\
     | TN SQ |
     \---------/ */

   
  else if( strcasecmp( str, "TNSQ") == 0 ) // EFAC for given flag
    {
      int nequadFlag = psr->nTNSQ;
      fscanf( fin, "%s %s %lf", psr->TNSQFlagID[nequadFlag],
              psr->TNSQFlagVal[nequadFlag],
              &psr->TNSQVal[nequadFlag] );
      ( psr->nTNSQ )++;
    }
  
  /* /---------\
     | TN ECORR |
     \---------/ */
  else if( strcasecmp( str, "ECORR") == 0 || strcasecmp( str, "TNECORR") == 0)  
	  // ECORR for given flag
    {
      int necorrFlag = psr->nTNECORR;
      fscanf( fin, "%s %s %lf", psr->TNECORRFlagID[necorrFlag],
              psr->TNECORRFlagVal[necorrFlag],
              &psr->TNECORRVal[necorrFlag] );
      ( psr->nTNECORR )++;
    }

   /* /---------\
     | TN Noise |
     \---------/ */
     
  else if (strcasecmp(str,"TNRedAmp")==0) /* TempoNest Red noise power law amplitude */
    fscanf(fin,"%lf",&(psr->TNRedAmp));
  else if (strcasecmp(str,"TNRedGam")==0) /* TempoNest Red noise spectral index */
	fscanf(fin,"%lf",&(psr->TNRedGam));
  else if (strcasecmp(str,"TNRedC")==0) /* TempoNest Red noise spectral index */
        fscanf(fin,"%d",&(psr->TNRedC));
  else if (strcasecmp(str,"TNRedFLow")==0) /* TempoNest Red noise power law amplitude */
    fscanf(fin,"%lf",&(psr->TNRedFLow));
  else if (strcasecmp(str,"TNRedCorner")==0) /* TempoNest Red noise spectral index */
        fscanf(fin,"%lf",&(psr->TNRedCorner));
  else if(strcasecmp(str,"TNsubtractRed")==0)
	fscanf(fin,"%d",&(psr->TNsubtractRed));
  else if (strcasecmp(str,"TNDMAmp")==0) /* TempoNest Red noise power law amplitude */
	fscanf(fin,"%lf",&(psr->TNDMAmp));
  else if (strcasecmp(str,"TNDMGam")==0) /* TempoNest Red noise spectral index */
	fscanf(fin,"%lf",&(psr->TNDMGam));
  else if (strcasecmp(str,"TNDMC")==0) /* TempoNest Red noise spectral index */
        fscanf(fin,"%d",&(psr->TNDMC));
  else if(strcasecmp(str,"TNsubtractDM")==0)
        fscanf(fin,"%d",&(psr->TNsubtractDM));
  else if(strcasecmp(str,"RNAMP")==0){ /* compatibility with tempo RN notation */
	printf("\nWARNING: Using tempo RNAMP parameter: setting TNRedC to 100!\n");
	
	char tmpstr[1000];
	fscanf(fin,"%s", tmpstr);
	for (int i=0;i<(int)strlen(tmpstr);i++)
	{
	 if (tmpstr[i]=='D' || tmpstr[i]=='d')
		tmpstr[i]='e';
	}
	sscanf(tmpstr,"%lf",&(psr->TNRedAmp));
	psr->TNRedAmp = log10(2.0*M_PI*pow(3.0,0.5)/(86400.0*365.25*1e6)*psr->TNRedAmp);
	psr->TNRedC = 100; /* Since tempo doesn't have this just hard-code to 100 */
  }
  else if (strcasecmp(str,"RNIDX")==0){ /* Tempo Red noise spectral index */
	fscanf(fin,"%lf",&(psr->TNRedGam));
	psr->TNRedGam *= -1.0; /* Flip sign convention */
  }


   /* /---------\
     | TNBandDM |
     \---------/ */
     
  else if (strcasecmp(str,"TNBandDM")==0){ /* TempoNest Band DM */   
	fscanf( fin, "%lf %lf %d", 
		&psr->TNBandDMAmp,
		&psr->TNBandDMGam,
		&psr->TNBandDMC);
	}    

   /* /-------------\
      | TNBandNoise |
      \-------------/ */

  else if (strcasecmp(str,"TNBandNoise")==0){ /* TempoNest Band DM */
	int nTNBandNoiseFlag = psr->nTNBandNoise;
        fscanf( fin, "%lf %lf %lf %lf %d",
                &psr->TNBandNoiseLF[nTNBandNoiseFlag],
                &psr->TNBandNoiseHF[nTNBandNoiseFlag],
                &psr->TNBandNoiseAmp[nTNBandNoiseFlag],	
		&psr->TNBandNoiseGam[nTNBandNoiseFlag],
		&psr->TNBandNoiseC[nTNBandNoiseFlag]);
		( psr->nTNBandNoise )++;

        }


   /* /-------------\
     | TNGroupNoise |
     \------------ -/ */
     
  else if (strcasecmp(str,"TNGroupNoise")==0){ /* TempoNest Group Noise */   
		int nTNGroupNoiseFlag = psr->nTNGroupNoise;
		fscanf( fin, "%s %s %lf %lf %d", psr->TNGroupNoiseFlagID[nTNGroupNoiseFlag],
			psr->TNGroupNoiseFlagVal[nTNGroupNoiseFlag],
			&psr->TNGroupNoiseAmp[nTNGroupNoiseFlag],
			&psr->TNGroupNoiseGam[nTNGroupNoiseFlag],
			&psr->TNGroupNoiseC[nTNGroupNoiseFlag]);
		( psr->nTNGroupNoise )++;
	}    



   /* /---------\
     | TNDMEvents |
     \---------/ */
     
  else if (strcasecmp(str,"TNDMEvent")==0){ /* TempoNest DM Event Start position */   
	int nTNDMEv = psr->nDMEvents;
	fscanf( fin, "%lf %lf %lf %lf %d %d %d", 
		&psr->TNDMEvStart[nTNDMEv],
		&psr->TNDMEvLength[nTNDMEv],
		&psr->TNDMEvAmp[nTNDMEv],
		&psr->TNDMEvGam[nTNDMEv],
		&psr->TNDMEvOff[nTNDMEv],
		&psr->TNDMEvLin[nTNDMEv],
		&psr->TNDMEvQuad[nTNDMEv]);

	( psr->nDMEvents )++;
	}    




   /* /----------------\
     | TNShapeletEvents |
     \-----------------/ */

     
  else if (strcasecmp(str,"TNShapeletEvent")==0){ /* TempoNest DM Shapelet Event Start position */   
	int nTNEv = psr->nTNShapeletEvents;
	fscanf( fin, "%d %lf %lf %lf", 
		&psr->TNShapeletEvN[nTNEv],
		&psr->TNShapeletEvPos[nTNEv],
		&psr->TNShapeletEvWidth[nTNEv],
		&psr->TNShapeletEvFScale[nTNEv]);

	( psr->nTNShapeletEvents )++;
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
		displayMsg(1,(char *)"BIN3",(char *)"Converting binary model to T2 model",
			  (char *)"",psr->noWarnings);
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
  {
	 readValue(psr,str,fin,&(psr->param[param_edot]),0);
         if (fabs(psr->param[param_edot].val[0]) > 1e-7) 
	   psr->param[param_edot].val[0] *= 1.0e-12;
  }
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
  else if  (strcasecmp(str,"OM2DOT")==0)
    {
      readValue(psr,str,fin,&(psr->param[param_om2dot]),0);
      
    }
	
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
  else if  (strcasecmp(str,"ORBPX")==0)
    readValue(psr,str,fin,&(psr->param[param_orbpx]),0);
  
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
  else if (str[0]=='F' && str[1]=='B')
    {
      int fbval;
      if (strlen(str)==2)
	fbval=0;
      else
	sscanf(str+2,"%d",&fbval);
      
      readValue(psr,str,fin,&(psr->param[param_fb]),fbval);
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
      if (fabs(psr->param[param_a1dot].val[0]) > 1e-7) /* Check units: DO BETTER JOB */
	psr->param[param_a1dot].val[0]*=1.0e-12;
      psr->param[param_a1dot].prefit[0] = psr->param[param_a1dot].val[0];
    }  
  else if ( strcasecmp(str,"X2DOT")==0) 
    {
      // Ryan: set this value which is used in THE MSS model plugin
      // These shouldn't be used in the same place, so hopefully this doesn't cause anything to crash!
      readValue(psr,str,fin,&(psr->param[param_a2dot]),0);
      
      psr->param[param_a2dot].prefit[0] = psr->param[param_a2dot].val[0];
    }
  
  else if ((strcasecmp(str,"A2DOT")==0) ||  (strcasecmp(str,"X2DOT")==0))
    {
      readValue(psr,str,fin,&(psr->param[param_a1dot]),1);
      psr->param[param_a1dot].prefit[1] = psr->param[param_a1dot].val[1];
      
      
    }
  
  else if (strcasecmp(str,"TASC")==0)
    readValue(psr,str,fin,&(psr->param[param_tasc]),0);
  else if (strcasecmp(str,"EPS1")==0)
    readValue(psr,str,fin,&(psr->param[param_eps1]),0);
  else if (strcasecmp(str,"EPS1DOT")==0)
    {
      
      readValue(psr,str,fin,&(psr->param[param_eps1dot]),0);
      if (fabs(psr->param[param_eps1dot].val[0]) > 1e-7) 
	psr->param[param_eps1dot].val[0] *= 1.0e-12;
      psr->param[param_eps1dot].prefit[0] = psr->param[param_eps1dot].val[0];
    }
  else if (strcasecmp(str,"EPS2")==0)
    readValue(psr,str,fin,&(psr->param[param_eps2]),0);
  else if (strcasecmp(str,"EPS2DOT")==0)
    {
      readValue(psr,str,fin,&(psr->param[param_eps2dot]),0);
      if (fabs(psr->param[param_eps2dot].val[0]) > 1e-7) 
	psr->param[param_eps2dot].val[0] *= 1.0e-12;
      psr->param[param_eps2dot].prefit[0] = psr->param[param_eps2dot].val[0];
    }
  else if (strcasecmp(str,"M2")==0)
    readValue(psr,str,fin,&(psr->param[param_m2]),0);
  else if (strcasecmp(str,"KOM")==0)
    readValue(psr,str,fin,&(psr->param[param_kom]),0);
  else if (strcasecmp(str,"KIN")==0)
    readValue(psr,str,fin,&(psr->param[param_kin]),0);
  else if (strcasecmp(str,"SHAPMAX")==0)
    readValue(psr,str,fin,&(psr->param[param_shapmax]),0);
  else if( strcasecmp( str, "H3" ) == 0 ){
    // h3 harmonic Shapiro delay parameter for DDH model (FW10)
    readValue( psr, str, fin, &( psr->param[param_h3] ), 0 );
  }else if( strcasecmp( str, "H4" ) == 0 ){
    // h4 harmonic Shapiro delay parameter for DDH model (FW10)
    readValue( psr, str, fin, &( psr->param[param_h4] ), 0 );
  }else if( strcasecmp( str, "STIG" ) == 0 ){
    // Stigma Shapiro delay harmonic ratio for DDH model (FW10)
    readValue( psr, str, fin, &( psr->param[param_stig] ), 0 );
  }else if( strcasecmp( str, "NHARM" ) == 0 ){
    // Number of Shapiro delay harmonics to be used for DDH model (FW10)
    readValue( psr, str, fin, &( psr->param[param_nharm] ), 0 );
  }
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
      displayMsg(1,(char *)"MISC1",(char *)"Unknown parameter in par file: ",str,psr->noWarnings);
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
	  // Do we need a value for waveepoch (and does waveepoch.paramSet
	  // need to be unity), even when we're not using waves for
	  // prewhitening?
	  //
	  // No we don't.
	  if(psr->param[param_wave_om].paramSet[0] == 1){
		 copyParam(psr->param[param_pepoch],&(psr->param[param_waveepoch]));
		 strcpy(psr->param[param_waveepoch].label[0],"WAVEEPOCH (MJD)");
		 strcpy(psr->param[param_waveepoch].shortlabel[0],"WAVEEPOCH"); 
		 printf("Setting waveepoch to %g\n",(double)psr->param[param_waveepoch].val[0]);
	  }else{
		 // waveepoch isn't set, but neither is wave_om. Ergo: we're
		 // not using waves and don't need to set this epoch.
		 //
		 // Do nothing.
	  }
   }
       if (psr->param[param_waveepoch_dm].paramSet[0] == 0)
   {
	  // Do we need a value for waveepoch (and does waveepoch.paramSet
	  // need to be unity), even when we're not using waves for
	  // prewhitening?
	  //
	  // No we don't.
	  if(psr->param[param_wave_dm].paramSet[0] == 1){
		 copyParam(psr->param[param_pepoch],&(psr->param[param_waveepoch_dm]));
		 strcpy(psr->param[param_waveepoch_dm].label[0],"WAVEEPOCHDM (MJD)");
		 strcpy(psr->param[param_waveepoch_dm].shortlabel[0],"WAVEEPOCHDM"); 
		 printf("Setting waveepoch to %g\n",(double)psr->param[param_waveepoch_dm].val[0]);
	  }else{
		 // waveepoch isn't set, but neither is wave_om. Ergo: we're
		 // not using waves and don't need to set this epoch.
		 //
		 // Do nothing.
	  }
   }


   if (psr->param[param_posepoch].paramSet[0] == 0)
   {
	  displayMsg(1,(char *)"PAR1",(char *)"Have not set a position epoch. The period epoch will be used instead.",filename,psr->noWarnings);
	  copyParam(psr->param[param_pepoch],&(psr->param[param_posepoch]));
	  strcpy(psr->param[param_posepoch].label[0],"POSEPOCH (MJD)");
	  strcpy(psr->param[param_posepoch].shortlabel[0],"POSEPOCH");
   }
   if (psr->param[param_dmepoch].paramSet[0] == 0)
   {
	  displayMsg(1,(char *)"PAR2",(char *)"Have not set a DM epoch. The period epoch will be used instead.",filename,psr->noWarnings);
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
	  strcpy(psr->clock, (char *)"TT(TAI)");
	  strcpy(psr->clockFromOverride, (char *)"UTC");
   }
   if ((psr->clock[0]!='T' || psr->clock[1]!='T'))
   {
	  char msg[1000];
	  sprintf(msg,"CLK parameter '%s' is not a realization of TT!",psr->clock);
	  displayMsg(1,(char *)"CLK1",msg,(char *)"",psr->noWarnings);

	  /* try various tempo possibilities */
	  char *clk;
	  if (!strcmp(psr->clock, "UTC(NIST)"))
		 clk = (char *)"TT(UTC(NIST))";
	  else if (!strcmp(psr->clock, "UTC(BIPM)"))
		 clk = (char *)"TT(TAI)";
	  else if (!strcmp(psr->clock, "PTB"))
		 clk = (char *)"TT(UTC(PTB))";
	  else if (!strcmp(psr->clock, "AT1"))
		 clk = (char *)"TT(TA(NIST))";
	  else /* default to TT(TAI) */
		 clk = (char *)"TT(TAI)";
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
   int nread,p,j,k;
   char str[1000];
   parameter elong,elat;	
   int noread=0,endit;
   const char *CVS_verNum = "$Revision: 1.90 $";

   if (displayCVSversion == 1) CVSdisplayVersion((char *)"readParfile.C",(char *)"readParfile()",CVS_verNum);

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
	  for (k=0;k<MAX_PARAMS;k++)
	  {
		 psr[p].param[k].nLinkTo = 0;
		 psr[p].param[k].nLinkFrom = 0;
	  }
	  psr[p].nJumps=0;

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
	  free(elat.label[j]);
	  free(elat.shortlabel[j]);
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
   else if (strcasecmp(pmtr,"KIN")==0)
   {
	  if(strcasecmp(str1,"STIG")==0)
	  {
		 parameter->linkTo[(parameter->nLinkTo)++] = param_stig;
		 psr->param[param_stig].linkFrom[(psr->param[param_stig].nLinkFrom)++] = param_kin;
		 psr->param[param_kin].paramSet[0]=1;
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
   } else if (strcasecmp(pmtr,"DMMODEL")==0){
	  if(strcasecmp(str1,"DM")==0){
		 parameter->linkTo[(parameter->nLinkTo)++] = param_dm;
		 psr->param[param_dm].linkFrom[(psr->param[param_dm].nLinkFrom)++]=param_dmmodel;
		 psr->param[param_dmmodel].paramSet[0]=1;
		 if(nread == 2){
			if (strcasecmp(str2,"1")==0 || strcasecmp(str2,"0")==0 || strcasecmp(str2,"2")==0) /* Have fit flag not error */
			{
			   sscanf(str2,"%d",&(parameter->fitFlag[arr]));
			}
		 }

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
   else if (strcasecmp(pmtr,"EPHVER")==0) // Check strings for EPHVER
   {
	  if (strcasecmp(str1,"2")==0 || strcasecmp(str1,"5")==0) 
		 parameter->val[arr] = parse_longdouble(str1);
	  else if (strcasecmp(str1,"TEMPO1")==0)
		 parameter->val[arr] = 2;
	  else if (strcasecmp(str1,"TEMPO2")==0)
		 parameter->val[arr]=5;
	  else
	  {
		 printf("ERROR: Unknown EPHVER >%s<\n",str1);
		 exit(1);
	  }
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
	  parameter->err[arr] = 0.0;
	  parameter->prefitErr[arr] = parse_longdouble(str3);
   }
   if (strcasecmp(pmtr,"RAJ")==0){
	  // convert to radians
	  parameter->prefitErr[arr]/=(12.0*60.0*60.0/M_PI);
   } else if (strcasecmp(pmtr,"DECJ")==0){
	  // convert to radians
	  parameter->prefitErr[arr]/=(180*60.0*60.0/M_PI);
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
