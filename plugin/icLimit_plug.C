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

using namespace std;


double getStatPS(pulsar *psr,int npsr,double gwAmp,double gwAlpha,int it,char *covarFuncFile);
long double getTspan(pulsar *psr,int npsr);
void formCholeskyMatrixPlugin(double *c,double *resx,double *resy,double *rese,int np,double **uinv);
void calculateGWCholesky(double modelAlpha,double modelFc,double modelScale,double fitVar,double **uinv,double *covFunc,
			 double *resx,double *resy,double *rese,int np);
  
void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,p,j;
  double globalParameter;
  double realStatistic;
  double simStatistic;
  const char *CVS_verNum = "$Revision$";
  double gwAmp,gwAlpha;
  double whiteNoise;

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
  long double mean,gwRes[MAX_OBSN];
  long double satIdeal[MAX_PSR][MAX_OBSN]; // Should use malloc

  int ngw;
  gwSrc *gw;
  // For simulations
  int it,nit;
  int ndetect;
  char covarFuncFile[128];

  
  char fname[128];

  gwAmp   = 1e-15; // Default value
  gwAlpha = -2.0/3.0; // Default value
  nit     = 1; // Default value

  if (displayCVSversion == 1) CVSdisplayVersion("icLimit.C","plugin",CVS_verNum);

  *npsr = 0;  

  printf("Graphical Interface: icLimit\n");
  printf("Author:              R. Shannon, G. Hobbs, M. Keith\n");
  printf("CVS Version:         $Revision$\n");
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
      else if (strcmp(argv[i],"-nit")==0)
	sscanf(argv[++i],"%d",&nit);
    }


  scale = pow(86400.0*365.25,gwAlpha);
  gwAmp *= scale;



  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      if (i==0) doFitDCM(psr,"NULL",covarFuncFile,*npsr,0);   /* Do the fitting     */
      else textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
    }

  realStatistic = getStatPS(psr,*npsr,gwAmp,gwAlpha,-1,covarFuncFile);

  // Form idealised site arrival times
  for (p=0;p<*npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)
	satIdeal[p][i] = psr[p].obsn[i].sat-psr[p].obsn[i].residual/86400.0L;
    }

  ngw = 1000;
  if ((gw = (gwSrc *)malloc(sizeof(gwSrc)*ngw))==NULL)
    {
      printf("Unable to allocate memory for %d GW sources\n",ngw);
      exit(1);
    }
  tspan=getTspan(psr,*npsr)*SECDAY; // Gets data span for all pulsars
  // Should choose something sensible for these
  flo=0.01/tspan;
  fhi = 1.0/(long double)SECDAY;
  timeOffset = psr[0].param[param_pepoch].val[0];

  ndetect = 0;
  for (it=0;it < nit; it++)
    {
      GWbackground(gw,ngw,&seed,flo,fhi,gwAmp,gwAlpha,1);
      for (i=0;i<ngw;i++)
	setupGW(&gw[i]);
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
	      whiteNoise = TKgaussDev(&seed)*psr[p].obsn[i].toaErr*1.0e-6;
	      psr[p].obsn[i].sat = satIdeal[p][i] + (gwRes[i]-mean/(long double)psr[p].nobs + whiteNoise)/86400.0L;
	    }
	  sprintf(fname,"psr%s.sim.tim.%d",psr[p].name,it);
	  writeTim(fname,psr+p,"tempo2");
	}
      // Do some fitting
      
      //
      // Must re-set the parameter values to their postfit value for the real data
      //
      for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
	{
	  formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
	  formResiduals(psr,*npsr,1);    /* Form the residuals                 */
	  if (i==0) doFit(psr,*npsr,0);   /* Do the fitting     */
	}
      // Must re-set the parameter values to their initial value
      printf("Simulation ...\n");
      simStatistic = getStatPS(psr,*npsr,gwAmp,gwAlpha,it,covarFuncFile);
      if (simStatistic > realStatistic) ndetect++;
      printf("Detected = %d, %d, %g percent\n",ndetect,it+1,(double)100*ndetect/(double)(it+1));
    }

  free(gw);

  return 0;
}

double getStatPS(pulsar *psr,int npsr,double gwAmp,double gwAlpha,int it,char *covarFuncFile)
{
  int i,p;
  double x[MAX_OBSN],y[MAX_OBSN],e[MAX_OBSN];
  double specX[MAX_OBSN],specY[MAX_OBSN],outY_re[MAX_OBSN],outY_im[MAX_OBSN];
  int    nSpec;
  int    specType,specOut;
  FILE   *fout;
  char   fname[128];
  double statistic=0.0,s_indv,signal,noise;
  // For Cholesky
  double **uinv;
  int ndays;
  double covarFunc[MAX_OBSN];
  double escaleFactor = 1.0;
  double fc,fitVar;
  double gwAmpYr;
  FILE *fin;

  specType = 2; // Unweighted Lomb-Scargle
  specOut  = 1; // Power spectral density

  uinv = (double **)malloc(sizeof(double *)*(MAX_OBSN));
  for (i=0;i<MAX_OBSN;i++)
    uinv[i] = (double *)malloc(sizeof(double)*MAX_OBSN);

  for (p=0;p<npsr;p++)
    {
      // Obtain a power spectrum
      // Assumes the data are time sorted ... should fix
      ndays = (int)(psr[p].obsn[psr[p].nobs-1].sat - psr[p].obsn[0].sat)+2;
      if (it==-1) // Using real data
	{
	  if (npsr>1)
	    sprintf(fname,"%s_%d",covarFuncFile,p+1);
	  else
	    strcpy(fname,covarFuncFile);

	  printf("Opening >%s<\n",fname);
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
	  printf("Read covariance function\n");
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
      if (it != -1) // Create a sensible covariance function for simulated data -- get sensible gwAmp
	{
	  fitVar = TKvariance_d(y,psr[p].nobs); // Check this
	  fc = 1.0/(x[psr[p].nobs-1]-x[0])*365.25; // Check this
	  //	  fc = 0.3; // CHANGE THIS
	  //	  gwAmp = 1.0e-30;
	  printf("fc = %g\n",fc);
	  calculateGWCholesky((3.0-2.0*gwAlpha)/2.0,fc,gwAmp,fitVar,uinv,covarFunc,
			      x,y,e,psr[p].nobs);
	    }
      // Form the data covariance matrix
      formCholeskyMatrixPlugin(covarFunc,x,y,e,psr[p].nobs,uinv);
      nSpec = calcSpectra(uinv,x,y,psr[p].nobs,specX,specY);
      //      TKspectrum(x,y,e,psr[p].nobs,0,0,0,0,0,specType,1,1,specOut,specX,specY,&nSpec,0,0,outY_re,outY_im);
      //      sprintf(fname,"psr%s.spec.%d",psr[p].name,it);
      //      fout = fopen(fname,"w");
      //      for (i=0;i<nSpec;i++)
      //	fprintf(fout,"%g %g\n",specX[i],specY[i]);
      //      fclose(fout);

      s_indv = 0.0; // Individual statistic
      gwAmpYr = gwAmp/pow(86400.0*365.25,gwAlpha);
      for (i=0;i<nSpec;i++)
	{
	  // Should use the correct units for gwAmp
	  // Should think about corner frequency or fitting
	  signal = gwAmpYr*gwAmpYr/12.0/M_PI/M_PI*pow(specX[i]*365.25,2*gwAlpha-3.0);
	  // Must calculate this properly
	  noise  = 5e-31;
	  s_indv += (specY[i]*signal*signal/(signal*signal+noise*noise));
	  //	  printf("individual statistic%d = %d %d %g %g %g %g %g %g\n",it,p,i,signal*signal/(signal*signal+noise*noise),signal,specX[i],specY[i],s_indv,specY[i]*signal*signal/(signal*signal+noise*noise));
	}
      printf("statistic for pulsar %d = %g\n",p,s_indv);
      statistic += s_indv;
    }
  printf("Total statistic = %g\n",statistic);

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
  int debug=1;

  printf("Getting the covariance matrix in doFit\n");
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
  printf("Inserting into the covariance matrix\n");
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
  printf("Multiplying by errors\n");
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
    }

  printf("Completed inverting the matrix\n");

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

void calculateGWCholesky(double modelAlpha,double modelFc,double modelScale,double fitVar,double **uinv,double *covFunc,
		       double *resx,double *resy,double *rese,int np)
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
  int debug=1;
  double varScaleFactor = 0.6;

  ndays = (int)(ceil(resx[np-1])-floor(resx[0])+1)+2; // Add two extra days for interpolation
  f  = (double *)malloc(sizeof(double)*ndays*2); // Frequency vector
  p  = (double *)malloc(sizeof(double)*ndays); // Model of pulsar power spectrum
  pf = (double *)malloc(sizeof(double)*ndays*2); // Periodic spectrum model
  opf = (double *)malloc(sizeof(double)*ndays*2); // Periodic spectrum model
  pe = (double *)malloc(sizeof(double)*ndays*2); // Periodic spectrum model

  printf("Number of days = %d\n",ndays);

  // Get rms of normalised high freq residuals
  /*  mean=0.0;
  for (i=0;i<np;i++)
    mean+=highFreqRes[i]/rese[i];
  mean/=(double)np;
  escale=0;
  for (i=0;i<np;i++)
    escale += pow(highFreqRes[i]/rese[i]-mean,2);
  escale/=(double)(np-1);
  escale = sqrt(escale);
  printf("Error scaling factor = %g\n",escale);
  // Scale error bars
  // NOTE: Have actually changed the rese[] array
      printf("ERROR NOT SCALING ERRORS\n");
      //printf("WARNING: scaling errors\n");
       escale = 1.0;
  
  for (i=0;i<np;i++)
    rese[i] *= escale;
  *errorScaleFactor *= escale;
  */

  // Weighted variance of residuals
  /* weightVarRes=0.0;
  weightMeanRes=0.0;
  tl=0.0; bl=0.0;
  tl2=0.0; bl2=0.0;
  for (i=0;i<np;i++)
    {
      tl+=resy[i]/(rese[i]*rese[i]);
      bl+=1.0/(rese[i]*rese[i]);
      tl2+=highFreqRes[i]/(rese[i]*rese[i]);
    }
  weightMeanRes = tl/bl;
  weightMeanHighFreqRes = tl2/bl;
  tl = bl = tl2 = 0.0;
  for (i=0;i<np;i++)
    {
      tl += pow(resy[i]-weightMeanRes,2)/(rese[i]*rese[i]);
      tl2 += pow(highFreqRes[i]-weightMeanHighFreqRes,2)/(rese[i]*rese[i]);
      bl += 1.0/(rese[i]*rese[i]);
    }
  weightVarRes = tl/bl;
  weightVarHighFreqRes = tl2/bl; */

  for (i=0;i<ndays;i++)
    {
      f[i] = i*1.0/(resx[np-1]-resx[0])*365.25;
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
  printf("Obtaining covariance function from analytic model\n");
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
	printf("covFunc: %g %g\n",opf[2*i],opf[2*i+1]);
      }
  }
  // Rescale
  printf("Rescaling %d\n",j/2);
  //  for (i=0;i<np/2;i++)
  //  fy2[i] = nmodelScale-log10(pow((1.0+pow(cholSpecX[i]*365.25/modelFc,2)),modelAlpha/2.0));
  actVar = pow(10,modelScale);
  tt = covFunc[0];
  printf("actvar = %g, weightVarRes = %g, weightVarHighFreqRes = %g, scale = %g, fitVar = %g, tt = %g\n",actVar*pow(86400.0*365.25,2),weightVarRes,weightVarHighFreqRes,weightVarRes-weightVarHighFreqRes,fitVar,tt);

  printf("WARNING: varScaleFactor = %g (used to deal with quadratic removal)\n",varScaleFactor);
  for (i=0;i<=j/2;i++)
    covFunc[i] = covFunc[i]*fitVar*varScaleFactor/tt;

  if (debug==1)
    {
      FILE *fout;
      fout = fopen("scaleCovar","w");
      for (i=0;i<=j/2;i++)
	fprintf(fout,"%d %g\n",i,covFunc[i]);
      fclose(fout);
    }
  // Now get the covariance matrix ...
  // First put the abs(time difference) in each matrix element
  //  formCholeskyMatrixPlugin(covFunc,resx,resy,rese,np,uinv);

  free(p);              
  free(f);
  free(pf);
  free(opf);
  free(pe);
}
