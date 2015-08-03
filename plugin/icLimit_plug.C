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
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "TKspectrum.h"
#include "T2toolkit.h"
#include "GWsim.h"
#include "fftw3.h"
#include "TKfit.h"
#include <cpgplot.h>

using namespace std;


//double getStatPS(pulsar *psr,int npsr,double gwAmp,double gwAlpha,int it,char *covarFuncFile,double noise,int plot);
double getSpectra(pulsar *psr,int npsr,char *covarFuncFile,double **specX,double **specY,int *nSpec);
double getStatPS(pulsar *psr,int npsr,double gwAmp,double gwAlpha,int it,char *covarFuncFile,double noise,int plot,double *specX,double *specY,int *nSpec);
void calculateWeighting(double *avSpecY,double *specX,int nSpec,double noiseLevel,double *weighting,double gwAmp,double gwAlpha);
double calculateStatistic(double **specY,double **weighting,int *nSpec,int npsr);
long double getTspan(pulsar *psr,int npsr);
void formCholeskyMatrixPlugin(double *c,double *resx,double *resy,double *rese,int np,double **uinv);
void calculateGWCholesky(double modelAlpha,double modelFc,double fitVar,double *covFunc, double dspan);
void createGWcovarianceFunction(char *file,double gwAmp,double gwAlpha,pulsar *psr,int npsr,double *gwVar);
  
void help() /* Display help */
{
  printf("-f mypar.par mytim.tim   Parameter and arrival time file\n");
  printf("-dcf                     Name of Cholesky covariance function file\n");
  printf("-gwamp                   Amplitude of GWB\n");
  printf("-gwalpha                 Spectral exponent of GWB\n");
  printf("-ahigh                   Used in a grid search: maximum possible amplitude of GWB\n");
  printf("-alow                    Used in a grid search: minimum possible amplitude of GWB\n");
  printf("-getlimit                Do a grid search\n");
  printf("-nit                     Number of iterations\n");
  printf("-ngw                     Number of gravitational waves to simulate\n");
  printf("-fast                    Do not recalculate clock files nor ephemeris positions\n");
  printf("-noise                   Noise level\n");
  printf("-plot                    Plot the results\n");
  printf("-h                       This help file\n");
  printf("-onlyF0F1                Re-do the fit for only F0 and F1 without using the Cholesky\n");
  printf("\n\n");
  printf("Typical usage: tempo2 -gr icLimit -gwamp 5e-15 -ngw 100 -nit 100 -dcf cfunc_1 -plot -f 0437_10cm.par try.tim -noise 2e-31\n");
  exit(1);
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,p,j,k;
  double globalParameter;
  double realStatistic;
  double simStatistic;
  const char *CVS_verNum = "$Revision: 1.10 $";
  double gwAmp,gwAlpha,aHigh,aLow;
  double whiteNoise;
  int getLimit=0;
  int plot=0;
  double noiseLevel[MAX_PSR];

  // For GW simulation
  long seed = TKsetSeed();
  long double timeOffset;
  long double ra_p,dec_p;
  long double flo=0.0,fhi=0.0;
  long double kp[3];            /* Vector pointing to pulsar           */
  long double tspan;
  long double time;
  long double scale;
  long double dist[MAX_PSR];
  long double mean;
  double gwRes[MAX_OBSN];
  long double satIdeal[MAX_PSR][MAX_OBSN]; // Should use malloc

  int ngw;
  gwSrc *gw;
  // For simulations
  int it,nit;
  int ndetect;
  char covarFuncFile[128];
  int endit=0,timeThrough=0;
  int fast=0;
  double gwVar[MAX_PSR];

  char fname[128];
  FILE *fout;

  // Store spectra
  double **specX_act,**specY_act;
  double **weighting;
  float  real_fx[MAX_OBSN],real_fy[MAX_OBSN];
  double **avSpecY;
  double ***specX_sim,***specY_sim;
  int    *nSpec;
  int    onlyF0F1=0;
  double x[MAX_OBSN];

  float  **actDataX,**actDataY,**actDataE1,**actDataE2;
  float  **simDataX,**simDataY,**simDataE1,**simDataE2;

  fout = fopen("icLimit_result.dat","w");

  gwAmp   = 1e-15; // Default value
  gwAlpha = -2.0/3.0; // Default value
  nit     = 1; // Default value
  ngw = 1000;

  if (displayCVSversion == 1) CVSdisplayVersion("icLimit.C","plugin",CVS_verNum);

  *npsr = 0;  

  printf("Graphical Interface: icLimit\n");
  printf("Author:              R. Shannon, G. Hobbs, M. Keith\n");
  printf("CVS Version:         $Revision: 1.10 $\n");
  printf(" --- type 'h' for help information\n");


  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[*npsr],argv[i+1]); 
	  strcpy(timFile[*npsr],argv[i+2]);
	  (*npsr)++;
	}
      else if (strcmp(argv[i],"-dcf")==0)
	strcpy(covarFuncFile,argv[++i]); // Must be in order with fname_1, fname_2 etc.
      else if (strcmp(argv[i],"-gwamp")==0)
	sscanf(argv[++i],"%lf",&gwAmp);
      else if (strcmp(argv[i],"-gwalpha")==0)
	sscanf(argv[++i],"%lf",&gwAlpha);
      else if (strcmp(argv[i],"-ahigh")==0)
	sscanf(argv[++i],"%lf",&aHigh);
      else if (strcmp(argv[i],"-alow")==0)
	sscanf(argv[++i],"%lf",&aLow);
      else if (strcmp(argv[i],"-getlimit")==0)
	getLimit=1;
      else if (strcmp(argv[i],"-nit")==0)
	sscanf(argv[++i],"%d",&nit);
      else if (strcmp(argv[i],"-ngw")==0)
	sscanf(argv[++i],"%d",&ngw);
      else if (strcmp(argv[i],"-fast")==0)
	fast=1;
      else if (strcmp(argv[i],"-noise")==0)
	sscanf(argv[++i],"%lf",&noiseLevel[(*npsr)-1]);
      else if (strcmp(argv[i],"-plot")==0)
	plot=1;
      else if (strcmp(argv[i],"-h")==0)
	help();
      else if (strcmp(argv[i],"-onlyF0F1")==0)
	onlyF0F1=1;
    }

  // Allocate memory for spectra
  printf("Allocating memory\n");
  nSpec = (int *)malloc(sizeof(int)*(*npsr));
  specX_sim = (double ***)malloc(sizeof(double **)*(nit));
  specY_sim = (double ***)malloc(sizeof(double **)*(nit));
  specX_act = (double **)malloc(sizeof(double *)*(*npsr));
  specY_act = (double **)malloc(sizeof(double *)*(*npsr));

  actDataX = (float **)malloc(sizeof(float *)*(*npsr));
  actDataY = (float **)malloc(sizeof(float *)*(*npsr));
  actDataE1 = (float **)malloc(sizeof(float *)*(*npsr));
  actDataE2 = (float **)malloc(sizeof(float *)*(*npsr));

  simDataX = (float **)malloc(sizeof(float *)*(*npsr));
  simDataY = (float **)malloc(sizeof(float *)*(*npsr));
  simDataE1 = (float **)malloc(sizeof(float *)*(*npsr));
  simDataE2 = (float **)malloc(sizeof(float *)*(*npsr));

  weighting = (double **)malloc(sizeof(double *)*(*npsr));
  avSpecY = (double **)malloc(sizeof(double *)*(*npsr));
  for (i=0;i<nit;i++)
    {
      specX_sim[i] = (double **)malloc(sizeof(double *)*(*npsr));
      specY_sim[i] = (double **)malloc(sizeof(double *)*(*npsr));
      for (j=0;j<*npsr;j++)
	{
	  specX_sim[i][j] = (double *)malloc(sizeof(double)*MAX_OBSN);
	  specY_sim[i][j] = (double *)malloc(sizeof(double)*MAX_OBSN);
	}
    }

  for (i=0;i<*npsr;i++)
    {
      specX_act[i] = (double *)malloc(sizeof(double)*MAX_OBSN);
      specY_act[i] = (double *)malloc(sizeof(double)*MAX_OBSN);
      weighting[i] = (double *)malloc(sizeof(double)*MAX_OBSN);
      avSpecY[i] = (double *)malloc(sizeof(double)*MAX_OBSN);
      actDataX[i] = (float *)malloc(sizeof(float)*MAX_OBSN);
      actDataY[i] = (float *)malloc(sizeof(float)*MAX_OBSN);
      actDataE1[i] = (float *)malloc(sizeof(float)*MAX_OBSN);
      actDataE2[i] = (float *)malloc(sizeof(float)*MAX_OBSN);
      simDataX[i] = (float *)malloc(sizeof(float)*MAX_OBSN);
      simDataY[i] = (float *)malloc(sizeof(float)*MAX_OBSN);
      simDataE1[i] = (float *)malloc(sizeof(float)*MAX_OBSN);
      simDataE2[i] = (float *)malloc(sizeof(float)*MAX_OBSN);
    }
  printf("Complete allocating memory\n");

  scale = pow(86400.0*365.25,gwAlpha);
  gwAmp *= scale;
  aLow *= scale;
  aHigh *= scale;

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

  // Check that the data are time sorted
  for (p=0;p<*npsr;p++)
    {
      for (i=0;i<psr[p].nobs-1;i++)
	{
	  if (psr[p].obsn[i].sat > psr[p].obsn[i+1].sat)
	    {
	      printf("ERROR: Data not time sorted for psr %s\n",psr[p].name);
	      exit(1);
	    }
	}
    }

  // Should update this to take the par file that the user gives .... (don't refit with the Cholesky)

  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      if (i==0) doFitDCM(psr,"NULL",covarFuncFile,*npsr,0);   /* Do the fitting     */
      else textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
    }

  getSpectra(psr,*npsr,covarFuncFile,specX_act,specY_act,nSpec);
  //  realStatistic = getStatPS(psr,*npsr,gwAmp,gwAlpha,-1,covarFuncFile,noiseLevel,plot,specX_act,specY_act,&nSpec);

  // Form idealised site arrival times
  for (p=0;p<*npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)
	{
	  actDataX[p][i] = (float)psr[p].obsn[i].sat;
	  actDataY[p][i] = (float)psr[p].obsn[i].residual/1e-6;
	  actDataE1[p][i] = actDataY[p][i]-(float)psr[p].obsn[i].toaErr;
	  actDataE2[p][i] = actDataY[p][i]+(float)psr[p].obsn[i].toaErr;
	  satIdeal[p][i] = psr[p].obsn[i].sat-psr[p].obsn[i].residual/86400.0L;
	}
    }

  if ((gw = (gwSrc *)malloc(sizeof(gwSrc)*ngw))==NULL)
    {
      printf("Unable to allocate memory for %d GW sources\n",ngw);
      exit(1);
    }
  tspan=getTspan(psr,*npsr)*SECDAY; // Gets data span for all pulsars
  //
  // Should choose something sensible for these
  //
  flo=0.01/tspan;
  fhi = 1.0/(long double)SECDAY;
  timeOffset = psr[0].param[param_pepoch].val[0];

  do {
    ndetect = 0;
    if (getLimit==1)
      {
	//	gwAmp = (aHigh + aLow)/2.0;
	gwAmp = pow(10,(log10(aHigh) + log10(aLow))/2.0);
	printf("Searching amplitude = %g, attempt number = %d\n",gwAmp/(double)scale,timeThrough+1);
      }

    for (it=0;it < nit; it++)
      {
	GWbackground(gw,ngw,&seed,flo,fhi,gwAmp,gwAlpha,1);
	for (i=0;i<ngw;i++) setupGW(&gw[i]);

	for (p=0;p<*npsr;p++)
	  {
	    // Should do this properly - should do a random selection for each realisation
	    dist[p] = 1e19;
	    
	    setupPulsar_GWsim(psr[p].param[param_raj].val[0],psr[p].param[param_decj].val[0],kp);
	    mean = 0.0;
	    for (i=0;i<psr[p].nobs;i++)
	      {
		time = (psr[p].obsn[i].sat - timeOffset)*SECDAY;
		gwRes[i] = 0.0;
		for (j=0;j<ngw;j++)
		  gwRes[i] += calculateResidualGW(kp,&gw[j],time,dist[p]);
		mean += gwRes[i];
	      }
	    for (i=0;i<psr[p].nobs;i++)
	      {
		//		psr[p].obsn[i].toaErr = psr[p].obsn[i].origErr = 0.01;
		//		psr[p].obsn[i].efac = 1;
		//		psr[p].nT2efac = 0;
		gwRes[i] -= (mean/(long double)psr[p].nobs);
		whiteNoise = TKgaussDev(&seed)*psr[p].obsn[i].toaErr*1.0e-6;
		psr[p].obsn[i].sat = satIdeal[p][i] + (gwRes[i] + whiteNoise)/86400.0L;
	      }
	    for (i=0;i<psr[p].nobs;i++)
		x[i] = actDataX[p][i];

	    TKremovePoly_d(x,gwRes,psr[p].nobs,3);
	    gwVar[p] = TKvariance_d(gwRes,psr[p].nobs); 
	    //	    	    sprintf(fname,"psr%s.sim.tim.%d",psr[p].name,it);
	    //	    	    writeTim(fname,psr+p,"tempo2");
	    //	    exit(1);
	  }
	// Do some fitting
	createGWcovarianceFunction("simFuncFile",gwAmp,gwAlpha,psr,*npsr,gwVar);
	
	if (fast==1)  // Can speed up by not re-calculating clocks etc.
	  {
	    vectorPulsar(psr,*npsr);   
	    calculate_bclt(psr,*npsr); 
	    formBats(psr,*npsr);       		
	  }
	else
	  formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
	
	formResiduals(psr,*npsr,1);    /* Form the residuals                 */
	doFitDCM(psr,"NULL","simFuncFile",*npsr,0);   /* Do the fitting     */
      
	// Form post-fit residuals
	if (fast==1)  // Can speed up by not re-calculating clocks etc.
	  {
	    vectorPulsar(psr,*npsr);   
	    calculate_bclt(psr,*npsr); 
	    formBats(psr,*npsr);       		
	  }
	else
	  formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
	
	formResiduals(psr,*npsr,1);    /* Form the residuals                 */
	if (onlyF0F1==1)
	  {
	    // Turn off fitting for everything except F0 and F1
	    for (i=0;i<MAX_PARAMS;i++)
	      {
		for (p=0;p<*npsr;p++)
		  {
		    if (i==0) // Turn off jumps
		      {
			for (j=0;j<psr[p].nJumps;j++)
			  psr[p].fitJump[j] = 0;		    
		      }
		    for (j=0;j<psr[p].param[i].aSize;j++)
		      {
			if (psr[p].param[i].paramSet[j] == 1 &&
			    psr[p].param[i].fitFlag[j] == 1)
			  {
			    if (i != param_f)
			      psr[p].param[i].fitFlag[j] = 0;
			  }
		      }
		  }
	      }
	    // Re-do the fit without the Cholesky
	    doFit(psr,*npsr,0);     
	  }	    
	// Re-form post-fit residuals
	if (fast==1)  // Can speed up by not re-calculating clocks etc.
	  {
	    vectorPulsar(psr,*npsr);   
	    calculate_bclt(psr,*npsr); 
	    formBats(psr,*npsr);       		
	  }
	else
	  formBatsAll(psr,*npsr);       
	
	formResiduals(psr,*npsr,1);   
	



	// Must re-set the parameter values to their initial value
	for (p=0;p<*npsr;p++)
	  {
	    psr[p].nJumps=0;
	    for(j=0;j<MAX_PARAMS;j++){
	      psr[p].param[j].nLinkTo = 0;
	      psr[p].param[j].nLinkFrom = 0;
	    }
	    if (plot==1)
	      {
		for (i=0;i<psr[p].nobs;i++)
		  {
		    simDataX[p][i] = (float)psr[p].obsn[i].sat;
		    simDataY[p][i] = (float)psr[p].obsn[i].residual/1e-6;
		    simDataE1[p][i] = simDataY[p][i]-(float)psr[p].obsn[i].toaErr;
		    simDataE2[p][i] = simDataY[p][i]+(float)psr[p].obsn[i].toaErr;
		  }
	      }
	  }
	readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */

	getSpectra(psr,*npsr,"simFuncFile",specX_sim[it],specY_sim[it],nSpec);	
	for (p=0;p<*npsr;p++)
	  {
	    for (i=0;i<nSpec[p];i++)
	      {
		for (k=0;k<it+1;k++)
		  {		    
		    if (k==0)
		      avSpecY[p][i] = specY_sim[k][p][i];
		    else
		      avSpecY[p][i] += specY_sim[k][p][i];
		  }
		avSpecY[p][i] /= (double)(it+1);
	      }
	  }
	//
	// MUST Fix avSpecY for multiple pulsars
	for (p=0;p<*npsr;p++)
	  calculateWeighting(avSpecY[p],specX_act[p],nSpec[p],noiseLevel[p],weighting[p],gwAmp,gwAlpha);

	// Now calculate the statistics
	realStatistic  = calculateStatistic(specY_act,weighting,nSpec,*npsr);
	ndetect=0;
	for (i=0;i<it+1;i++)
	  {
	    simStatistic = calculateStatistic(specY_sim[i],weighting,nSpec,*npsr);
	    if (simStatistic > realStatistic) ndetect++;
	  }
	printf("%d Real statistic = %g, detected = %g percent\n",it+1,realStatistic,(double)100*ndetect/(double)(it+1));


	//	if ((it+1)%100==0) printf("Detected = %d, iteration = %d, %g percent, %g\n",ndetect,it+1,(double)100*ndetect/(double)(it+1),(double)gwAmp/(double)scale);
	if ((it+1)%10==0) printf("Detected = %d, iteration = %d, %g percent, %g\n",ndetect,it+1,(double)100*ndetect/(double)(it+1),(double)gwAmp/(double)scale);



	if (plot==1)
	  {
	    float sv[it+1];
	    float av95[nSpec[0]];
	    float mean[nSpec[0]];
	    float fy2[nSpec[0]];
	    float sig[nSpec[0]];
	    float ftx[2],fty[2];
	    double gwAmpYr;
	    char grDev[128];
	    
	    for (p=0;p<*npsr;p++)
	      {
		sprintf(grDev,"%d/xs",30+p);
		

		// Find 95% level for simulated spectra
		for (i=0;i<nSpec[p];i++)
		  {
		    for (j=0;j<it+1;j++)
		      sv[j] = specY_sim[j][p][i];
		    TKsort_f(sv,it+1);
		    av95[i] = (float)log10(sv[(int)(5.0/100.0*(it+1)+0.5)]);
		    mean[i] = (float)log10(avSpecY[p][i]);
		  }

		cpgbeg(0,grDev,1,1);
		cpgask(0);
		
		//		cpgeras();
		cpgsch(1.0);
		cpgsvp(0.1,0.95,0.75,0.90);
		cpgswin(-4,0,0,1.1);
		cpgbox("BCTSL",0,0,"BCTSN",0,0);
		cpglab("","",psr[p].name);
		gwAmpYr = gwAmp/pow(86400.0*365.25,gwAlpha);


		for (i=0;i<nSpec[p];i++)
		  {
		    fy2[i] = (float)weighting[p][i];
		    sig[i] = log10(gwAmpYr*gwAmpYr/12.0/M_PI/M_PI*pow(specX_act[p][i]*365.25,2*gwAlpha-3.0));
		    real_fx[i] = (float)log10(specX_act[p][i]);
		    real_fy[i] = (float)log10(specY_act[p][i]);
		  }
		cpgsls(1); cpgsci(7); cpgline(nSpec[p],real_fx,fy2); cpgsci(1); cpgsls(1);
		
		cpgsvp(0.1,0.95,0.4,0.75);
		cpgswin(-4,0,-35,-25);
		cpgbox("BNCTSL",0,0,"BCTSLN",0,0);
		cpglab("Frequency (d\\u-1\\d)","PSD (yr\\u-3\\d)","");
		
		
		cpgline(nSpec[p],real_fx,real_fy);
		cpgsci(4); cpgline(nSpec[p],real_fx,sig); cpgsci(1);

		cpgsci(8); cpgsls(4); cpgline(nSpec[p],real_fx,av95); cpgsls(1); cpgsci(1);
		cpgsci(8); cpgsls(1); cpgline(nSpec[p],real_fx,mean); cpgsls(1); cpgsci(1);
		// Plot noise level
		ftx[0] = -4; ftx[1] = 1;
		fty[0] = fty[1] = log10(noiseLevel[p]);
		cpgsci(2); cpgsls(4); cpgline(2,ftx,fty); cpgsls(1); cpgsci(1);
		// Plot 1/1yr and 1/6month
		fty[0] = -35; fty[1] = -25;
		ftx[0] = ftx[1] = log10(1.0/365.25);
		cpgsci(3); cpgsls(4); cpgline(2,ftx,fty); cpgsls(1); cpgsci(1);
		ftx[0] = ftx[1] = log10(2.0/365.25);
		cpgsci(3); cpgsls(4); cpgline(2,ftx,fty); cpgsls(1); cpgsci(1);

		cpgsvp(0.1,0.95,0.10,0.30);
		cpgswin(psr[p].obsn[0].sat,psr[p].obsn[psr[p].nobs-1].sat,-1,1);
		cpgbox("BNCTNS",0,0,"BCTSN",0,0);
		cpglab("Day","Residual (us)","");
		cpgerry(psr[p].nobs,actDataX[p],actDataE1[p],actDataE2[p],1);
		cpgsch(0.5); cpgpt(psr[p].nobs,actDataX[p],actDataY[p],9); cpgsch(1);
		cpgsci(7);
		cpgerry(psr[p].nobs,simDataX[p],simDataE1[p],simDataE2[p],1);
		cpgsch(0.5); cpgpt(psr[p].nobs,simDataX[p],simDataY[p],9); cpgsch(1);
		cpgsci(1);
		cpgend();		
	      }

	  }


      }
    printf("Detected = %d, %d, %g percent, %g\n",ndetect,it,(double)100*ndetect/(double)(it),(double)gwAmp/(double)scale);
    timeThrough++;
    if (getLimit==1)
      {
	if (100*ndetect/(double)it > 96)
	  aHigh = gwAmp;
	else if (100*ndetect/(double)it < 94)
	  aLow = gwAmp;
	else
	  endit=1;

	if (timeThrough==10) endit=1; // Break if not converging
      }

  } while (getLimit==1 && endit==0);
  if (getLimit==1)
    {
      if (timeThrough==10) printf("Result = NO CONVERGENCE\n");
      else printf("Result = %g\n",gwAmp/(double)scale);
      if (timeThrough==10) fprintf(fout,"Result = NO CONVERGENCE\n");
      else fprintf(fout,"Result = %g ntrial = %d\n",gwAmp/(double)scale,timeThrough);
    }
  free(gw);
  fclose(fout);

  for (i=0;i<nit;i++)
    {
      for (j=0;j<*npsr;j++)
	{
	  free(specX_sim[i][j]);
	  free(specY_sim[i][j]);
	}
      free(specX_sim[i]);
      free(specY_sim[i]);
    }
  for (i=0;i<*npsr;i++)
    {
      free(specX_act[i]);
      free(specY_act[i]);      
      free(avSpecY[i]);
      free(weighting[i]);      
      free(actDataX[i]);
      free(actDataY[i]);
      free(actDataE1[i]);
      free(actDataE2[i]);
      free(simDataX[i]);
      free(simDataY[i]);  
      free(simDataE1[i]);
      free(simDataE2[i]);
    }
  free(specX_sim);
  free(specY_sim);
  free(specX_act);
  free(specY_act);
  free(nSpec);
  free(weighting);
  free(avSpecY);
  free(actDataX);
  free(actDataY);
  free(actDataE1);
  free(actDataE2);
  free(simDataX);
  free(simDataY);
  free(simDataE1);
  free(simDataE2);

  return 0;
}

double calculateStatistic(double **specY,double **weighting,int *nSpec,int npsr)
{
  int i,p;
  double stat=0;

  for (p=0;p<npsr;p++)
    {
      for (i=0;i<nSpec[p];i++)
	{
	  stat += specY[p][i]*weighting[p][i];
	}
    }
  return stat;
}

void calculateWeighting(double *avSpecY,double *specX,int nSpec,double noiseLevel,double *weighting,double gwAmp,double gwAlpha)
{
  int i;
  double exSig;
  double gwAmpYr;

  gwAmpYr = gwAmp/pow(86400.0*365.25,gwAlpha);

  for (i=0;i<nSpec;i++)
    {
      weighting[i] = pow(avSpecY[i],2)/(pow(avSpecY[i],2)+pow(noiseLevel,2));
      exSig = gwAmpYr*gwAmpYr/12.0/M_PI/M_PI*pow(specX[i]*365.25,2*gwAlpha-3.0);
      weighting[i] *= pow(exSig,2)/(pow(exSig,2)+pow(noiseLevel,2));
    }
}

double getSpectra(pulsar *psr,int npsr,char *covarFuncFile,double **specX,double **specY,int *nSpec)
{
  double **uinv;
  double ndays;
  int i,j,k,p;
  char fname[128];
  double x[MAX_OBSN],y[MAX_OBSN],e[MAX_OBSN],escaleFactor;
  double covarFunc[MAX_OBSN];
  FILE *fin;

  //
  // Should allocate this to the correct size required
  //
  uinv = (double **)malloc(sizeof(double *)*(MAX_OBSN));
  for (i=0;i<MAX_OBSN;i++)
    uinv[i] = (double *)malloc(sizeof(double)*MAX_OBSN);

  for (p=0;p<npsr;p++)
    {
      // Obtain a power spectrum
      ndays = (int)(psr[p].obsn[psr[p].nobs-1].sat - psr[p].obsn[0].sat)+2;
      //      if (it==-1) // Using real data
	{
	  if (npsr>1)
	    sprintf(fname,"%s_%d",covarFuncFile,p+1);
	  else
	    strcpy(fname,covarFuncFile);

	  if (!(fin = fopen(fname,"r")))
	    {
	      printf("Unable to open covariance function file: %s\n",fname);
	      exit(1);
	    }
	  if (debugFlag==1) printf("ndays = %d\n",ndays);
	  fscanf(fin,"%lf",&escaleFactor);
	  for (i=0;i<ndays;i++)
	    fscanf(fin,"%lf",&covarFunc[i]);
	  fclose(fin);
	}

	for (i=0;i<psr[p].nobs;i++)
	  {
	    //
	    // Should check for deleted points ****
	    //
	    x[i] = (double)(psr[p].obsn[i].sat-psr[p].param[param_pepoch].val[0]); // In days
	    y[i] = (double)(psr[p].obsn[i].residual); // In seconds
	    e[i] = (double)(psr[p].obsn[i].toaErr*1.0e-6); // In seconds	 
	  }
	// Form the data covariance matrix
	//
	// Really only need to do this when the gwamp changes
	formCholeskyMatrixPlugin(covarFunc,x,y,e,psr[p].nobs,uinv);
	nSpec[p] = calcSpectra(uinv,x,y,psr[p].nobs,specX[p],specY[p],-1);  // 10 = number of frequency channels
    }	
  
  for (i=0;i<MAX_OBSN;i++)
    free(uinv[i]);
  free(uinv);  
}

double getStatPS(pulsar *psr,int npsr,double gwAmp,double gwAlpha,int it,char *covarFuncFile,double noise,int plot,double *specX,double *specY,int *nSpec)
{
  int i,p;
  double x[MAX_OBSN],y[MAX_OBSN],e[MAX_OBSN];
  double outY_re[MAX_OBSN],outY_im[MAX_OBSN];
  int    specType,specOut;
  FILE   *fout;
  char   fname[128];
  double statistic=0.0,s_indv,signal;
  // For Cholesky
  double **uinv;
  int ndays;
  float fx[MAX_OBSN],fy[MAX_OBSN],fy2[MAX_OBSN];
  static float real_fx[MAX_OBSN_VAL],real_fy[MAX_OBSN_VAL];
  float fnx[2],fny[2];
  static float avSpec[MAX_OBSN_VAL];
  float avSpecY[MAX_OBSN];
  float fysignal[MAX_OBSN];
  double covarFunc[MAX_OBSN];
  double escaleFactor = 1.0;
  double fc,fitVar;
  double gwAmpYr;
  FILE *fin;

  specType = 2; // Unweighted Lomb-Scargle
  specOut  = 1; // Power spectral density
  
  //
  // Should allocate this to the correct size required
  //
  uinv = (double **)malloc(sizeof(double *)*(MAX_OBSN));
  for (i=0;i<MAX_OBSN;i++)
    uinv[i] = (double *)malloc(sizeof(double)*MAX_OBSN);

  for (p=0;p<npsr;p++)
    {
      // Obtain a power spectrum
      ndays = (int)(psr[p].obsn[psr[p].nobs-1].sat - psr[p].obsn[0].sat)+2;
      //      if (it==-1) // Using real data
	{
	  if (npsr>1)
	    sprintf(fname,"%s_%d",covarFuncFile,p+1);
	  else
	    strcpy(fname,covarFuncFile);

	  //	  printf("Opening >%s<\n",fname);
	  if (!(fin = fopen(fname,"r")))
	    {
	      printf("Unable to open covariance function file: %s\n",fname);
	      exit(1);
	    }
	  if (debugFlag==1) printf("ndays = %d\n",ndays);
	  fscanf(fin,"%lf",&escaleFactor);
	  for (i=0;i<ndays;i++)
	    fscanf(fin,"%lf",&covarFunc[i]);
	  fclose(fin);
	  //	  printf("Read covariance function\n");
	}

      //      printf("WARNING: scaling all errors by: %g\n",escaleFactor);
      //      for (i=0;i<psr[p].nobs;i++)
      //	sig[i]*=escaleFactor;
	  
	  for (i=0;i<psr[p].nobs;i++)
	    {
	  //
	  // Should check for deleted points ****
	  //
	  x[i] = (double)(psr[p].obsn[i].sat-psr[p].param[param_pepoch].val[0]); // In days
	  y[i] = (double)(psr[p].obsn[i].residual); // In seconds
	  e[i] = (double)(psr[p].obsn[i].toaErr*1.0e-6); // In seconds	 
	}
      // Form the data covariance matrix
      //
      // Really only need to do this when the gwamp changes
      formCholeskyMatrixPlugin(covarFunc,x,y,e,psr[p].nobs,uinv);
      *nSpec = calcSpectra(uinv,x,y,psr[p].nobs,specX,specY,-1);  // 10 = number of frequency channels

      //      TKspectrum(x,y,e,psr[p].nobs,0,0,0,0,0,specType,1,1,specOut,specX,specY,&nSpec,0,0,outY_re,outY_im);
      //      sprintf(fname,"psr%s.spec.%d",psr[p].name,it);
      //      fout = fopen(fname,"w");
      //      for (i=0;i<nSpec;i++)
      //	fprintf(fout,"%g %g\n",specX[i],specY[i]);
      //      fclose(fout);

      s_indv = 0.0; // Individual statistic
      gwAmpYr = gwAmp/pow(86400.0*365.25,gwAlpha);
      for (i=0;i<*nSpec;i++)
	{
	  // Should use the correct units for gwAmp
	  // Should think about corner frequency or fitting
	  signal = gwAmpYr*gwAmpYr/12.0/M_PI/M_PI*pow(specX[i]*365.25,2*gwAlpha-3.0);
	  fysignal[i] = log10(signal);
	  //
	  //	  s_indv += (specY[i]*signal*signal/(signal*signal+noise*noise));
	  // Ryan's statistic
	  s_indv += (specY[i]*signal/(signal+noise));
	  if (plot==1)
	    {
	      fx[i] = log10(specX[i]);
	      fy[i] = log10(specY[i]);
	      fy2[i] = signal/(signal+noise);
	      //	      printf("y = %g\n",fy2[i]);
	    }

	  //	  printf("individual statistic%d = %d %d %g %g %g %g %g %g\n",it,p,i,signal*signal/(signal*signal+noise*noise),signal,specX[i],specY[i],s_indv,specY[i]*signal*signal/(signal*signal+noise*noise));
	}
      //      printf("statistic for pulsar %d = %g\n",p,s_indv);
      statistic += s_indv;
      if (plot==1)
	{
	  cpgeras();
	  cpgsch(1.0);
	  cpgsvp(0.1,0.95,0.8,0.95);
	  cpgswin(-4,0,0,1.1);
	  cpgbox("BCTSL",0,0,"BCTSN",0,0);
	  cpgsls(4); cpgsci(7); cpgline(*nSpec,fx,fy2); cpgsci(1); cpgsls(1);
	  
	  cpgsvp(0.1,0.95,0.1,0.8);
	  cpgswin(-4,0,-35,-27);
	  cpgbox("BNCTSL",0,0,"BCTSLN",0,0);
	  cpglab("Frequency (d\\u-1\\d)","PSD (yr\\u-3\\d)","");

	  fnx[0] = -4; fnx[1] = 1;
	  fny[0] = log10(noise); fny[1] = log10(noise);
	  cpgsls(2); cpgline(2,fnx,fny); cpgsls(1);
	  fny[0] = -35; fny[1] = -27;
	  fnx[0] = fnx[1] = log10(1.0/(365.25));
	  cpgsls(2); cpgline(2,fnx,fny); cpgsls(1);
	  fnx[0] = fnx[1] = log10(2.0/(365.25));
	  cpgsls(2); cpgline(2,fnx,fny); cpgsls(1);
	  
	  // Plot expected GWB signal based on amplitude
	  cpgsci(7); cpgline(*nSpec,fx,fysignal); cpgsci(1);	  
	}
    }
  
  
  //  printf("Total statistic = %g\n",statistic);
  
  for (i=0;i<MAX_OBSN;i++)
    free(uinv[i]);
  free(uinv);

  return statistic;
}


long double getTspan(pulsar *psr,int npsr)
{
  long double first,last;
  int i,p;
    
  
  first = psr[0].obsn[0].sat;
  last = psr[0].obsn[0].sat;

  for (p=0;p<npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)
	{
	  if (first > psr[p].obsn[i].sat)
	    first = psr[p].obsn[i].sat;
	  if (last < psr[p].obsn[i].sat)
	    last = psr[p].obsn[i].sat;
	}
    }

  return last-first;
}

void formCholeskyMatrixPlugin(double *c,double *resx,double *resy,double *rese,int np,double **uinv)
{
  double **m,**u,sum;
  double *cholp;
  int i,j,k,ix,iy;
  double t0,cint,t;
  int t1,t2;
  int debug=0;

  //  printf("Getting the covariance matrix in doFit\n");
  m = (double **)malloc(sizeof(double *)*(np+1));
  u= (double **)malloc(sizeof(double *)*(np+1));
  cholp  = (double *)malloc(sizeof(double)*(np+1));  // Was ndays

  for (i=0;i<np+1;i++)
    {
      m[i] = (double *)malloc(sizeof(double)*(np+1));
      u[i] = (double *)malloc(sizeof(double)*(np+1));
    }
  
  for (ix=0;ix<np;ix++)
    {
      for (iy=0;iy<np;iy++)
	m[ix][iy] = fabs(resx[ix]-resx[iy]);
    }
  if (debug==1)
    {
      printf("First m = \n");
      for (i=0;i<5;i++)
	{ 
	  for (j=0;j<5;j++) printf("%10g ",m[i][j]); 
	  printf("\n");
	}

    }
  // Insert the covariance which depends only on the time difference.
  // Linearly interpolate between elements on the covariance function because
  // valid covariance matrix must have decreasing off diagonal elements.
  //  printf("Inserting into the covariance matrix\n");
  for (ix=0;ix<np;ix++)
    {
      for (iy=0;iy<np;iy++)
	{
	  t0 = m[ix][iy];
	  t1 = (int)floor(t0);
	  t2 = t1+1;
	  t  = t0-t1;
	  cint = c[t1]*(1-t)+c[t2]*t; // Linear interpolation
	  m[ix][iy] = cint;
	}
    }
  //  printf("Multiplying by errors\n");
  for (ix=0;ix<np;ix++)
    m[ix][ix]+=rese[ix]*rese[ix];
  if (debug==1)
    {
      printf("m = \n\n");
      for (i=0;i<5;i++)
	{ 
	  for (j=0;j<5;j++) printf("%10g ",m[i][j]); 
	  printf("\n");
	}
    }

  // Do the Cholesky
  TKcholDecomposition(m,np,cholp);
  // Now calculate uinv
  for (i=0;i<np;i++)
    {
      m[i][i] = 1.0/cholp[i];
      uinv[i][i] = m[i][i];
      for (j=0;j<i;j++)
      	uinv[i][j] = 0.0;
      for (j=i+1;j<np;j++)
	{
	  sum=0.0;
	  for (k=i;k<j;k++) sum-=m[j][k]*m[k][i];
	  m[j][i]=sum/cholp[j];
	  uinv[i][j] = m[j][i];
	}
    } 

  if (debug==1)
    {
      printf("uinv = \n\n");
      for (i=0;i<5;i++)
	{ 
	  for (j=0;j<5;j++) printf("%10g ",uinv[i][j]); 
	  printf("\n");
	}
      printf("Completed inverting the matrix\n");
    }



  // Should free memory not required
  // (note: not freeing uinv)

  for (i=0;i<np+1;i++)
    {
      free(m[i]);
      free(u[i]);
    }
  free(m);
  free(u);
  free(cholp);
}

void calculateGWCholesky(double modelAlpha,double modelFc,double fitVar,double *covFunc, double dspan)
{
  int i,j;
  double *f; // Frequency vector
  double *p; // Model of pulsar power spectrum
  double *pe;
  double *pf; // Periodic spectrum model
  double *opf; 
  double weightVarRes,weightMeanRes;
  double weightVarHighFreqRes,weightMeanHighFreqRes;
  double mean,escale,actVar,tt,tl,bl,tl2,bl2;
  int ndays;
  int debug=0;
  double varScaleFactor = 0.6;

  ndays = (int)(dspan+1)+2; // Add two extra days for interpolation
  f  = (double *)malloc(sizeof(double)*ndays*2); // Frequency vector
  p  = (double *)malloc(sizeof(double)*ndays); // Model of pulsar power spectrum
  pf = (double *)malloc(sizeof(double)*ndays*2); // Periodic spectrum model
  opf = (double *)malloc(sizeof(double)*ndays*2); // Periodic spectrum model
  pe = (double *)malloc(sizeof(double)*ndays*2); // Periodic spectrum model

  //  printf("Number of days = %d\n",ndays);

  // Get rms of normalised high freq residuals

  // Weighted variance of residuals

  for (i=0;i<ndays;i++)
    {
      f[i] = i*1.0/(dspan)*365.25;
      // Choose a sensible model
      p[i] = 1.0/(pow((1.0+pow(f[i]/(modelFc*2),2)),(modelAlpha)/2.0));
      //p[i] = pow(f[i],-13.0/3.0);
      pf[i] = p[i];
    }
  j = ndays;
  for (i=ndays-1;i>0;i--)
    {
      f[j] = -f[i];
      // Choose a sensible model
      pf[j] = 1.0/(pow((1.0+pow(f[j]/(modelFc*2),2)),modelAlpha/2.0));
      //pf[j] = pow(f[i],-13.0/3.0);
      j++;
    }
  if (debug==1)
    {
      FILE *fout;
      fout = fopen("specModel","w");
      for (i=0;i<j;i++)
	fprintf(fout,"%d %g\n",i,pf[i]);
      fclose(fout);     
    }
  //  printf("Obtaining covariance function from analytic model\n");
  {
    fftw_complex* output;
    fftw_plan transform_plan;
    
    output = (fftw_complex*)opf;
    transform_plan = fftw_plan_dft_r2c_1d(j, pf, output, FFTW_ESTIMATE);
    fftw_execute(transform_plan);    
    fftw_destroy_plan(transform_plan);  
    for (i=0;i<=j/2;i++) 
      {
	covFunc[i] = opf[2*i];
	//	printf("covFunc: %g %g\n",opf[2*i],opf[2*i+1]);
      }
  }

  tt = covFunc[0];

  for (i=0;i<=j/2;i++)
    covFunc[i] = covFunc[i]*fitVar*varScaleFactor/tt;

    //    exit(1);
  // Now get the covariance matrix ...
  // First put the abs(time difference) in each matrix element
  //  formCholeskyMatrixPlugin(covFunc,resx,resy,rese,np,uinv);

  free(p);              
  free(f);
  free(pf);
  free(opf);
  free(pe);
}

void createGWcovarianceFunction(char *file,double gwAmp,double gwAlpha,pulsar *psr,int npsr,double *gwVar)
{
  double fitVar,fc;
  double x[MAX_OBSN],y[MAX_OBSN];
  double covarFunc[MAX_OBSN],dspan;
  int i,p;
  char fnameuse[128];
  FILE *fout;

  for (p=0;p<npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)
	{
	  //
	  // Should check for deleted points ****
	  //
	  x[i] = (double)(psr[p].obsn[i].sat-psr[p].param[param_pepoch].val[0]); // In days
	  y[i] = (double)(psr[p].obsn[i].residual); // In seconds
	}
      
      fitVar = gwVar[p]; // TKvariance_d(y,psr[p].nobs); 
      fc = 1.0/(x[psr[p].nobs-1]-x[0])*365.25;
      dspan = ceil(x[psr[p].nobs-1])-floor(x[0]);
      calculateGWCholesky((3.0-2.0*gwAlpha),fc,fitVar,covarFunc, dspan);
      
      if (npsr == 1)
	sprintf(fnameuse,"%s",file,p);
      else
	sprintf(fnameuse,"%s_%d",file,p+1);
      fout = fopen(fnameuse,"w");
      fprintf(fout,"1\n");
      for (i=0;i<=(dspan+3);i++)
	fprintf(fout,"%g\n",covarFunc[i]);
      fclose(fout);
      //    exit(1);
    }
}
char * plugVersionCheck = TEMPO2_h_VER;
