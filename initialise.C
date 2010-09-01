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

/* Initialise all the parameters */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"

void initialise(pulsar *psr,int noWarnings)
{
  int p;

  // allocate space for covar arrays
  int fullSetup = 1;

  for (p=0;p<MAX_PSR;p++)
    //    if (psr[p].obsn == NULL)
      initialiseOne (psr+p, noWarnings, fullSetup);
}

void initialiseOne (pulsar *psr, int noWarnings, int fullSetup)
{
  int fail = 0;
  int i,j,k;
  char temp[100];

  psr->nobs = 0;
  //  psr->obsn = NULL;
  psr->covar = NULL;

  //  if (psr->obsn == NULL)
    psr->obsn = (observation *)malloc(sizeof(observation)*MAX_OBSN);

  if (psr->obsn == NULL)
    fail = 1;

  if (fullSetup && !fail) {

    psr->covar = (double **)malloc(sizeof(double *)*MAX_PARAMS);
    
    if (psr->covar != NULL) {
      for(i=0;i<MAX_PARAMS;i++) {
        psr->covar[i] = (double *)malloc(sizeof(double)*MAX_PARAMS);
        if (psr->covar[i] == NULL) {
          fail = 1;
          break;
        }
      }
    }
  }
  // Initialise the barycentre vectors:
  //
  // To the best of my knowledge, this only ever happens in
  // preProcessSimple.C, which does not get called when running Tempo2
  // through PSRchive. For reasons unknown to me, these uninitialised
  // values most commonly end up being zero (as they should be) and
  // are subsequently set later on, but once in a statistical while,
  // they will go rogue, crashing everything. For my purposes,
  // initialisation right here seems to do the trick, though a
  // substantial reanalysis of the preprocess and initialise functions
  // (and which bits belong where) is in order to make sure these
  // "glitches" don't happen with any other parameters.
  //
  // JPWV, MPIfR, 22 July 2010.
  for(int obsct = 0; obsct < MAX_OBSN; obsct++ ){
    for (int vecct = 0; vecct < 6; vecct++ ){
      psr->obsn[obsct].earthMoonBary_ssb[vecct] = 0.0;
      psr->obsn[obsct].earthMoonBary_earth[vecct] = 0.0;
      psr->obsn[obsct].observatory_earth[vecct] = 0.0;
      psr->obsn[obsct].earth_ssb[vecct] = 0.0;
    }
  }
  if (fail)
    {
      printf("Not enough memory to allocate room for %d observations\n",MAX_OBSN);
      printf("Please decrease the value of MAX_OBSN_VAL in tempo2.h or use -nobs on the command line\n"); 
      printf("You can also decrease the number of pulsars being stored in memory using -npsr\n");
      printf("Note: 1 observation requires %f kBytes, you request %f MBytes\n",(double)sizeof(observation)/1024.0,(double)sizeof(observation)/1024.0/1024.0*MAX_OBSN);
      exit(1);
    }  /* This memory gets deallocated by destroyOne */


  strcpy(psr->filterStr,"");
  strcpy(psr->passStr,"");
  strcpy(psr->fitFunc,"default");
  psr->nits=1;
  psr->clockFromOverride[0] = '\0';
  psr->nCompanion = 0;
  psr->bootStrap = 0;
  psr->units = SI_UNITS;
  psr->ne_sw  = NE_SW_DEFAULT; 
  psr->nWhite = 0;  /* No whitening by default */
  psr->ifuncN = 0;  /* No interpolation functions by default */
  psr->timeEphemeris = IF99_TIMEEPH;
  psr->dilateFreq = 1;
  psr->planetShapiro = 1;
  psr->correctTroposphere = 1;
  psr->t2cMethod = T2C_IAU2000B;
  psr->fixedFormat=0;
  psr->nStorePrecision=0;
  strcpy(psr->fjumpID,"");
  strcpy(psr->deleteFileName,"NONE");
  strcpy(psr->tzrsite,"NULL");
  psr->calcShapiro=1;
  psr->dmOffset = 0;
  psr->ipm = 1;
  psr->swm = 0;
  psr->nPhaseJump=0;
  psr->eclCoord=0;
  psr->noWarnings=noWarnings;
  psr->fitMode = 0;         /* Don't fit with errors by default (MODE 0) */
  psr->rescaleErrChisq = 1; /* Rescale parameter errors by reduced chisq  */
  strcpy(psr->name,"NOT SET");
  strcpy(psr->binaryModel,"NONE"); 
  psr->nJumps=0;
  psr->nToffset = 0;
  psr->ndmx = 0;
  psr->jboFormat=0;
  for (i=0;i<MAX_JUMPS;i++)
    {
      psr->jumpVal[i] = 0.0;
      psr->jumpValErr[i] = 0.0;
    }
  
  allocateMemory(psr,0);
  /*  psr->param[param_track].paramSet[0]=1;
  psr->param[param_track].val[0]=0.0;
  psr->param[param_track].prefit[0]=0.0;*/
  
  /* Spin-frequency parameters */
  for (j=0;j<psr->param[param_f].aSize;j++)
    {
      psr->param[param_f].val[j] = 0.0;
      sprintf(temp,"F%d (s^-%d)",j,j+1);
      strcpy(psr->param[param_f].label[j],temp);   
      sprintf(temp,"F%d",j);
      strcpy(psr->param[param_f].shortlabel[j],temp);
    }
  strcpy(psr->param[param_raj].label[0],"RAJ (rad)");
  strcpy(psr->param[param_raj].shortlabel[0],"RAJ");
  strcpy(psr->param[param_decj].label[0],"DECJ (rad)");
  strcpy(psr->param[param_decj].shortlabel[0],"DECJ");
  strcpy(psr->param[param_fddi].label[0],"FDDI"); /* Frequency dependent delay */
  strcpy(psr->param[param_fddi].shortlabel[0],"FDDI"); /* Frequency dependent delay */
  strcpy(psr->param[param_fddc].label[0],"FDDC"); /* Frequency dependent delay */
  strcpy(psr->param[param_fddc].shortlabel[0],"FDDC"); /* Frequency dependent delay */
  /* Dispersion measure and its derivative */
  for (k=0;k<psr->param[param_dm].aSize;k++)
    {
      if (k>0){
	sprintf(temp,"DM%d (cm^-3 pc yr^-%d)",k,k);
	strcpy(psr->param[param_dm].label[k],temp);
	sprintf(temp,"DM%d",k);
	strcpy(psr->param[param_dm].shortlabel[k],temp);
      }
      else
	{
	  strcpy(psr->param[param_dm].label[0],"DM (cm^-3 pc)");
	  strcpy(psr->param[param_dm].shortlabel[0],"DM");
	}
    }
  strcpy(psr->param[param_px].label[0],"PX (mas)");
  strcpy(psr->param[param_px].shortlabel[0],"PX");
  strcpy(psr->param[param_daop].label[0],"AOP dist. (kpc)");
  strcpy(psr->param[param_daop].shortlabel[0],"D_AOP");
  strcpy(psr->param[param_daop].label[0],"IPERHARM");
  strcpy(psr->param[param_daop].shortlabel[0],"IPERHARM");
  strcpy(psr->param[param_pmrv].label[0],"PMRV (mas/yr)"); strcpy(psr->param[param_pmrv].shortlabel[0],"PMRV");
  for (k=0;k<psr->param[param_dmassplanet].aSize;k++)
    {
      sprintf(psr->param[param_dmassplanet].label[k], "DMASSPLANET%d (Msun)", k+1);
      sprintf(psr->param[param_dmassplanet].shortlabel[k], "DMASSPLANET%d", k+1);
    }
  strcpy(psr->param[param_tres].label[0],"TRES");
  strcpy(psr->param[param_tres].shortlabel[0],"TRES");
  strcpy(psr->param[param_ephver].label[0],"EPHVER");
  strcpy(psr->param[param_ephver].shortlabel[0],"EPHVER");
  strcpy(psr->param[param_pmra].label[0],"PMRA (mas/yr)");
  strcpy(psr->param[param_pmra].shortlabel[0],"PMRA");
  strcpy(psr->param[param_pmdec].label[0],"PMDEC (mas/yr)");
  strcpy(psr->param[param_pmdec].shortlabel[0],"PMDEC");
  strcpy(psr->param[param_posepoch].label[0],"POSEPOCH (MJD)");
  strcpy(psr->param[param_posepoch].shortlabel[0],"POSEPOCH");
  strcpy(psr->param[param_waveepoch].label[0],"WAVEEPOCH (MJD)");
  strcpy(psr->param[param_waveepoch].shortlabel[0],"WAVEEPOCH");
  strcpy(psr->param[param_ifunc].label[0],"IFUNC");
  strcpy(psr->param[param_ifunc].shortlabel[0],"IFUNC");
  strcpy(psr->param[param_pepoch].label[0],"PEPOCH (MJD)");
  strcpy(psr->param[param_pepoch].shortlabel[0],"PEPOCH");
  strcpy(psr->param[param_dmepoch].label[0],"DMEPOCH (MJD)");
  strcpy(psr->param[param_dmepoch].shortlabel[0],"DMEPOCH");
  strcpy(psr->param[param_start].label[0],"START (MJD)");
  strcpy(psr->param[param_start].shortlabel[0],"START");
  strcpy(psr->param[param_finish].label[0],"FINISH (MJD)");
  strcpy(psr->param[param_finish].shortlabel[0],"FINISH");
  strcpy(psr->param[param_track].label[0],"TRACK (MJD)");
  strcpy(psr->param[param_track].shortlabel[0],"TRACK");
  strcpy(psr->param[param_dshk].label[0],"DSHK (kpc)");
  strcpy(psr->param[param_dshk].shortlabel[0],"DSHK");
  /* Glitch parameters */
  for (k=0;k<psr->param[param_glep].aSize;k++)
    {
      sprintf(temp,"GLEP_%d",k+1);
      strcpy(psr->param[param_glep].label[k],temp);
      strcpy(psr->param[param_glep].shortlabel[k],temp);
      sprintf(temp,"GLPH_%d",k+1);
      strcpy(psr->param[param_glph].label[k],temp);
      strcpy(psr->param[param_glph].shortlabel[k],temp);
      sprintf(temp,"GLF0_%d",k+1);
      strcpy(psr->param[param_glf0].label[k],temp);
      strcpy(psr->param[param_glf0].shortlabel[k],temp);
      sprintf(temp,"GLF1_%d",k+1);
      strcpy(psr->param[param_glf1].label[k],temp);
      strcpy(psr->param[param_glf1].shortlabel[k],temp);
      sprintf(temp,"GLF2_%d",k+1);
      strcpy(psr->param[param_glf2].label[k],temp);
      strcpy(psr->param[param_glf2].shortlabel[k],temp);
      sprintf(temp,"GLF0D_%d",k+1);
      strcpy(psr->param[param_glf0d].label[k],temp);
      strcpy(psr->param[param_glf0d].shortlabel[k],temp);
      sprintf(temp,"GLTD_%d",k+1);
      strcpy(psr->param[param_gltd].label[k],temp);
      strcpy(psr->param[param_gltd].shortlabel[k],temp);
    }
  /* Binary parameters */
  strcpy(psr->param[param_t0].label[0],"T0 (MJD)");
  strcpy(psr->param[param_t0].shortlabel[0],"T0");
  /* Dispersion measure and its derivative */
  for (k=1;k<psr->param[param_pb].aSize;k++)
    {
      sprintf(temp,"PB_%d (d)",k); strcpy(psr->param[param_pb].label[k],temp);
      sprintf(temp,"PB_%d",k); strcpy(psr->param[param_pb].shortlabel[k],temp);
	  
      sprintf(temp,"ECC_%d",k); strcpy(psr->param[param_ecc].label[k],temp);
      sprintf(temp,"ECC_%d",k); strcpy(psr->param[param_ecc].shortlabel[k],temp);

      sprintf(temp,"OM_%d (deg)",k); strcpy(psr->param[param_om].label[k],temp);
      sprintf(temp,"OM_%d",k); strcpy(psr->param[param_om].shortlabel[k],temp);

      sprintf(temp,"A1_%d (lt-s)",k); strcpy(psr->param[param_a1].label[k],temp);
      sprintf(temp,"A1_%d",k); strcpy(psr->param[param_a1].shortlabel[k],temp);

      sprintf(temp,"T0_%d (mjd)",k); strcpy(psr->param[param_t0].label[k],temp);
      sprintf(temp,"T0_%d",k); strcpy(psr->param[param_t0].shortlabel[k],temp);
    }
  strcpy(psr->param[param_pb].label[0],"PB (d)");
  strcpy(psr->param[param_pb].shortlabel[0],"PB");
  strcpy(psr->param[param_a1].label[0],"A1 (lt-s)");
  strcpy(psr->param[param_a1].shortlabel[0],"A1");
  strcpy(psr->param[param_om].label[0],"OM (deg)");
  strcpy(psr->param[param_om].shortlabel[0],"OM");
  strcpy(psr->param[param_ecc].label[0],"ECC");
  strcpy(psr->param[param_edot].shortlabel[0],"EDOT");
  strcpy(psr->param[param_edot].label[0],"EDOT");
  strcpy(psr->param[param_ecc].shortlabel[0],"ECC");
  strcpy(psr->param[param_kom].label[0],"KOM"); 
  strcpy(psr->param[param_kom].shortlabel[0],"KOM");
  strcpy(psr->param[param_kin].label[0],"KIN"); 
  strcpy(psr->param[param_kin].shortlabel[0],"KIN");
  strcpy(psr->param[param_shapmax].label[0],"SHAPMAX"); 
  strcpy(psr->param[param_shapmax].shortlabel[0],"SHAPMAX");
  strcpy(psr->param[param_m2].label[0],"M2");
  strcpy(psr->param[param_m2].shortlabel[0],"M2");
  strcpy(psr->param[param_mtot].label[0],"MTOT");
  strcpy(psr->param[param_mtot].shortlabel[0],"MTOT");
  strcpy(psr->param[param_dr].label[0],"DR"); strcpy(psr->param[param_dr].shortlabel[0],"DR");
  strcpy(psr->param[param_dth].label[0],"DTH"); strcpy(psr->param[param_dth].shortlabel[0],"DTH");
  strcpy(psr->param[param_a0].label[0],"A0"); strcpy(psr->param[param_a0].shortlabel[0],"A0");
  strcpy(psr->param[param_b0].label[0],"B0"); strcpy(psr->param[param_b0].shortlabel[0],"B0");
  strcpy(psr->param[param_bp].label[0],"BP"); strcpy(psr->param[param_bp].shortlabel[0],"BP");
  strcpy(psr->param[param_bpp].label[0],"BPP"); strcpy(psr->param[param_bpp].shortlabel[0],"BPP");
  strcpy(psr->param[param_dtheta].label[0],"DTHETA"); 
  strcpy(psr->param[param_dtheta].shortlabel[0],"DTHETA");
  strcpy(psr->param[param_sini].label[0],"SINI");
  strcpy(psr->param[param_sini].shortlabel[0],"SINI");
  strcpy(psr->param[param_gamma].label[0],"GAMMA");
  strcpy(psr->param[param_gamma].shortlabel[0],"GAMMA");
  strcpy(psr->param[param_pbdot].label[0],"PBDOT");
  strcpy(psr->param[param_pbdot].shortlabel[0],"PBDOT");
  strcpy(psr->param[param_xpbdot].label[0],"XPBDOT");
  strcpy(psr->param[param_xpbdot].shortlabel[0],"XPBDOT");
  strcpy(psr->param[param_a1dot].label[0],"XDOT");
  strcpy(psr->param[param_a1dot].shortlabel[0],"XDOT");
  strcpy(psr->param[param_xomdot].label[0],"XOMDOT");
  strcpy(psr->param[param_xomdot].shortlabel[0],"XOMDOT");
  strcpy(psr->param[param_afac].label[0],"AFAC");
  strcpy(psr->param[param_afac].shortlabel[0],"AFAC");
  strcpy(psr->param[param_omdot].label[0],"OMDOT (deg/yr)");
  strcpy(psr->param[param_omdot].shortlabel[0],"OMDOT");
  strcpy(psr->param[param_tasc].label[0],"TASC (MJD)");
  strcpy(psr->param[param_tasc].shortlabel[0],"TASC");
  strcpy(psr->param[param_eps1].label[0],"EPS1");
  strcpy(psr->param[param_eps1].shortlabel[0],"EPS1");
  strcpy(psr->param[param_eps1dot].label[0],"EPS1DOT");
  strcpy(psr->param[param_eps1dot].shortlabel[0],"EPS1DOT");
  strcpy(psr->param[param_eps2].label[0],"EPS2");
  strcpy(psr->param[param_eps2].shortlabel[0],"EPS2");
  strcpy(psr->param[param_eps2dot].label[0],"EPS2DOT");
  strcpy(psr->param[param_eps2dot].shortlabel[0],"EPS2DOT");
  strcpy(psr->param[param_tzrmjd].label[0],"TZRMJD");
  strcpy(psr->param[param_tzrmjd].shortlabel[0],"TZRMJD");
  strcpy(psr->param[param_tzrfrq].label[0],"TZRFRQ (MHz)");
  strcpy(psr->param[param_tzrfrq].shortlabel[0],"TZRFRQ");
  strcpy(psr->param[param_tspan].label[0],"TSPAN (min)"); 
  strcpy(psr->param[param_tspan].shortlabel[0],"TSPAN");

  for (k=0;k<psr->param[param_bpjep].aSize;k++)
    {
      sprintf(temp,"BPJEP_%d",k+1);
      strcpy(psr->param[param_bpjep].label[k],temp); 
      strcpy(psr->param[param_bpjep].shortlabel[k],temp);

      sprintf(temp,"BPJPH_%d",k+1);
      strcpy(psr->param[param_bpjph].label[k],temp); 
      strcpy(psr->param[param_bpjph].shortlabel[k],temp);

      sprintf(temp,"BPJA1_%d",k+1);
      strcpy(psr->param[param_bpja1].label[k],temp); 
      strcpy(psr->param[param_bpja1].shortlabel[k],temp);

      sprintf(temp,"BPJEC_%d",k+1);
      strcpy(psr->param[param_bpjec].label[k],temp); 
      strcpy(psr->param[param_bpjec].shortlabel[k],temp);

      sprintf(temp,"BPJOM_%d",k+1);
      strcpy(psr->param[param_bpjom].label[k],temp); 
      strcpy(psr->param[param_bpjom].shortlabel[k],temp);

      sprintf(temp,"BPJPB_%d",k+1);
      strcpy(psr->param[param_bpjpb].label[k],temp); 
      strcpy(psr->param[param_bpjpb].shortlabel[k],temp);


    }       

  strcpy(psr->param[param_wave_om].label[0],"WAVE_OM"); strcpy(psr->param[param_wave_om].shortlabel[0],"WAVE_OM");
  strcpy(psr->param[param_dmval].label[0],"DMVAL"); strcpy(psr->param[param_dmval].shortlabel[0],"DMVAL");

  /* Piecewise-constant DM variation (DMX) */
  for (k=0;k<psr->param[param_dmx].aSize;k++) 
    {
      sprintf(temp,"DMX_%04d (cm^-3 pc)",k+1);
      strcpy(psr->param[param_dmx].label[k],temp);
      sprintf(temp,"DMX_%04d",k+1);
      strcpy(psr->param[param_dmx].shortlabel[k],temp);

      sprintf(temp,"DMXR1_%04d (MJD)",k+1);
      strcpy(psr->param[param_dmxr1].label[k],temp);
      sprintf(temp,"DMXR1_%04d",k+1);
      strcpy(psr->param[param_dmxr1].shortlabel[k],temp);

      sprintf(temp,"DMXR2_%04d (MJD)",k+1);
      strcpy(psr->param[param_dmxr2].label[k],temp);
      sprintf(temp,"DMXR2_%04d",k+1);
      strcpy(psr->param[param_dmxr2].shortlabel[k],temp);
    }
}


void allocateMemory(pulsar *psr, int realloc)
{
  int i,j;
  for (i=0;i<MAX_PARAMS;i++)    
    {
      psr->param[i].nLinkTo   = 0;
      psr->param[i].nLinkFrom = 0;      

      if (i==param_dm)      psr->param[i].aSize = MAX_DM_DERIVATIVES;
      else if (i==param_f)  psr->param[i].aSize = MAX_FREQ_DERIVATIVES;
      else if (i==param_pb || i==param_ecc || i==param_om || i==param_t0 || i==param_a1) 
	psr->param[i].aSize = MAX_COMPANIONS;
      else if (i==param_bpjep || i==param_bpjph || i==param_bpja1 || i==param_bpjec || i==param_bpjom
	       || i==param_bpjpb)  psr->param[i].aSize = MAX_BPJ_JUMPS;
      else if (i==param_glep || i==param_glph || i==param_glf0 || i==param_glf1 || i==param_glf2 || 
	       i==param_glf0d || i==param_gltd) psr->param[i].aSize = 20;
      else if (i==param_dmassplanet)
	psr->param[i].aSize = 9;
      else if (i==param_dmx || i==param_dmxr1 || i==param_dmxr2)
        psr->param[i].aSize = MAX_DMX;
      else psr->param[i].aSize = 1;
      
      psr->param[i].val       = (longdouble *)malloc(psr->param[i].aSize*sizeof(longdouble));
      psr->param[i].err       = (longdouble *)malloc(psr->param[i].aSize*sizeof(longdouble));
      psr->param[i].prefit    = (longdouble *)malloc(psr->param[i].aSize*sizeof(longdouble));
      psr->param[i].prefitErr = (longdouble *)malloc(psr->param[i].aSize*sizeof(longdouble));
      psr->param[i].fitFlag   = (int *)malloc(psr->param[i].aSize*sizeof(int));
      psr->param[i].paramSet  = (int *)malloc(psr->param[i].aSize*sizeof(int));
      psr->param[i].label     = (char **)malloc(psr->param[i].aSize*sizeof(char *));
      psr->param[i].shortlabel= (char **)malloc(psr->param[i].aSize*sizeof(char *));

      for (j=0;j<psr->param[i].aSize;j++)
	{
	  psr->param[i].label[j] = (char *)malloc(sizeof(char)*100);
	  psr->param[i].shortlabel[j] = (char *)malloc(sizeof(char)*100);

	  psr->param[i].fitFlag[j]  = 0;
	  psr->param[i].paramSet[j] = 0;
	  psr->param[i].err[j]      = 0.0;	     
	  psr->param[i].val[j]      = 0.0;
	}
    }
}


void destroyOne (pulsar *psr)
{
  int i = 0;

  if (psr->obsn)
    free (psr->obsn);

  if (psr->covar) {
    for(i=0;i<MAX_PARAMS;i++)
      free (psr->covar[i]);
    free (psr->covar);
  }
  
  destroyMemory(psr);
}

void destroyMemory (pulsar *psr)
{
  int i,j;
  for (i=0;i<MAX_PARAMS;i++)    
    {
      for (j=0;j<psr->param[i].aSize;j++)
	{
	  free( psr->param[i].label[j] );
	  free( psr->param[i].shortlabel[j] );
	}
      
      free( psr->param[i].val );
      free( psr->param[i].err );
      free( psr->param[i].prefit );
      free( psr->param[i].prefitErr );
      free( psr->param[i].fitFlag );
      free( psr->param[i].paramSet );
      free( psr->param[i].label );
      free( psr->param[i].shortlabel );

    }
}
