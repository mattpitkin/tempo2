// 1. Try and get rid of the need for a corner frequency
// 2. Look at how white the residuals actually are
// 3. It seems as if most of the specX are okay, but the average is high.  Look at this!
// 4. Try not as steep red noise

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
#include "fftw3.h"
#include "TKfit.h"
#include "T2toolkit.h"
#include <cpgplot.h>

using namespace std;

#define MAX_SPEC_BINS 1024
double OMEGA0=0; 

void calculateSpectrum(pulsar *psr,double *px,double *py_r,double *py_i,int *nSpec);
void formCholeskyMatrixPlugin(double *c,double *resx,double *resy,double *rese,int np,double **uinv);
int calcSpectra_plugin(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY_R,double *specY_I,int nfit);
void calculateCovarFunc(double modelAlpha,double modelFc,double modelScale,double *covFunc,double *resx,double *resy,double *rese,int np);

void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  double px[MAX_SPEC_BINS];
  double py_r[MAX_SPEC_BINS];
  double py_i[MAX_SPEC_BINS];
  double tspan,minx,maxx;
  int nSpec=0;
  int i;
  double globalParameter;
  const char *CVS_verNum = "$Revision$";

  if (displayCVSversion == 1) CVSdisplayVersion("grTemplate.C","plugin",CVS_verNum);

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: cholSpectra\n");
  printf("Author:              G. Hobbs\n");
  printf("CVS Version:         $Revision$\n");
  printf(" --- type 'h' for help information\n");

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[0],argv[i+1]); 
	  strcpy(timFile[0],argv[i+2]);
	}
    }

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      if (i==0) doFit(psr,*npsr,0);   /* Do the fitting     */
      else textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
    }
  printf("Producing spectra\n");
  for (i=0;i<psr[0].nobs;i++)
    {
      if (i==0)
	{
	  minx = maxx = (double)psr[0].obsn[i].sat;
	}
      else
	{
	  if (minx > (double)psr[0].obsn[i].sat) minx = (double)psr[0].obsn[i].sat;
	  if (maxx < (double)psr[0].obsn[i].sat) maxx = (double)psr[0].obsn[i].sat;
	}
    }
  tspan = maxx-minx;
  OMEGA0 = (double)(2*M_PI/tspan);
  printf("Doing the calculation\n");
  calculateSpectrum(psr,px,py_r,py_i,&nSpec);
  printf("Finished\n");
  return 0;
}

char * plugVersionCheck = TEMPO2_h_VER;

void calculateSpectrum(pulsar *psr,double *px,double *py_r,double *py_i,int *nSpec)
{
  int i;
  double pe[MAX_SPEC_BINS];
  double covarFunc[10000];
  double **uinv;
  FILE *fin;
  char fname[128];
  double escaleFactor;
  int ndays=0;
  double resx[psr->nobs],resy[psr->nobs],rese[psr->nobs];
  FILE *fout;
  long double toffset = 52601.0L;

  //  printf("Calculating the spectrum\n");
  uinv = (double **)malloc(sizeof(double *)*psr->nobs);
  for (i=0;i<psr->nobs;i++)
    uinv[i] = (double *)malloc(sizeof(double)*psr->nobs);

  sprintf(fname,"%s.res",psr->name);
  if (!(fout = fopen(fname,"w")))
    {
      printf("Unable to open output file\n");
      exit(1);
    }
  for (i=0;i<psr->nobs;i++)
    {
      if (psr->obsn[i].deleted != 0)
	{
	  printf("GET RID OF THE DELETED POINTS\n");
	  exit(1);
	}
      resx[i] = (double)(psr->obsn[i].sat-toffset);
      resy[i] = (double)(psr->obsn[i].residual);
      rese[i] = (double)(psr->obsn[i].toaErr*1.0e-6);
      fprintf(fout,"%g %g %g\n",resx[i],resy[i],rese[i]);
    }
  fclose(fout);


  //  printf("Allocated memory\n");

  // Read in the covariance function
  sprintf(fname,"covarFunc.dat_%s",psr->name);
  if (!(fin = fopen(fname,"r")))
    {
      printf("Unable to open data file >%s<\n",psr->name);
      exit(1);
    }
  fscanf(fin,"%lf",&escaleFactor);
  while (!feof(fin))
    {
      if (fscanf(fin,"%lf",&covarFunc[ndays])==1)
	ndays++;
    }
    fclose(fin);

    //  calculateCovarFunc(double modelAlpha,double modelFc,double modelScale,double *covFunc,double *resx,double *resy,double *rese,int np)


  //  calculateCovarFunc(13.0/3.0,0.09,1.61528e-24/2,covarFunc,resx,resy,rese,psr[0].nobs);
  //      calculateCovarFunc(13.0/3.0,0.09,2.8716e-26,covarFunc,resx,resy,rese,psr[0].nobs);
  //        calculateCovarFunc(13.0/3.0,0.001,8.4434e-18,covarFunc,resx,resy,rese,psr[0].nobs);
  //    calculateCovarFunc(13.0/3.0,0.001,2.8716e-26,covarFunc,resx,resy,rese,psr[0].nobs);


  //  printf("Read in the covariance function\n");

  // Form the uinv matrix
  //  printf("Forming uinv\n");
  formCholeskyMatrixPlugin(covarFunc,resx,resy,rese,psr->nobs,uinv);
  //  printf("Formed uinv\n");
  //  printf("Read the covariance function file\n");
  *nSpec = 124;

  // Must calculate uinv for the pulsar
  calcSpectra_plugin(uinv,resx,resy,psr->nobs,px,py_r,py_i,*nSpec);
  sprintf(fname,"%s.spec",psr->name);
  fout = fopen(fname,"w");
  for (i=0;i<*nSpec;i++)
    fprintf(fout,"%g %g\n",px[i],py_r[i]*py_r[i]+py_i[i]*py_i[i]);
  fclose(fout);
  // Free uinv
  for (i=0;i<psr->nobs;i++)
    free(uinv[i]);
  free(uinv);
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
  //  printf("Allocated memory\n");
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
  //  if (debug==1)
    {
      printf("np = %d, m = \n\n",np);
      for (i=np-5;i<np;i++)
	{ 
	  for (j=np-5;j<np;j++) printf("%10g ",m[i][j]);
	  printf("\n");
	}
    }

  // Do the Cholesky
  //  printf("Cholesky decomposition\n");
  TKcholDecomposition(m,np,cholp);
  //  printf("Complete cholesky decomposition\n");
  // Now calculate uinv
  for (i=0;i<np;i++)
    {
      //      printf("i = %d ... %d\n",i,np);
      m[i][i] = 1.0/cholp[i];
      //      printf("s1\n");
      uinv[i][i] = m[i][i];
      //      printf("s2\n");
      for (j=0;j<i;j++)
      	uinv[i][j] = 0.0;
      //      printf("s3\n");
      for (j=i+1;j<np;j++)
	{
	  sum=0.0;
	  for (k=i;k<j;k++) sum-=m[j][k]*m[k][i];
	  m[j][i]=sum/cholp[j];
	  uinv[i][j] = m[j][i];
	}
      //      printf("s4\n");
    } 
  //  printf("Complete cholesky\n");
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


int calcSpectra_plugin(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY_R,double *specY_I,int nfit)
{
  int i,j,k;
  //  int nfit=nres/2-1;
  int nSpec;

  if (nfit < 0)
    nfit=nres/2-1;

  double v[nfit];
  double sig[nres];
  double **newUinv;
  double **cvm;
  double chisq;
  pulsar *psr;
  int ip[nres];
  double param[nfit],error[nfit];

  cvm = (double **)alloca(sizeof(double *)*nfit);
  for (i=0;i<nfit;i++)
    cvm[i] = (double *)alloca(sizeof(double)*nfit);

  // Should fit independently to all frequencies
  for (i=0;i<nres;i++)
    {
      sig[i] = 1.0; // The errors are built into the uinv matrix
      ip[i] = 0;
    }
  for (k=0;k<nfit;k++)
    {
      //      printf("k = %d\n",k);
      //      printf("%5.2g\%\r",(double)k/(double)nfit*100.0);
      //      fflush(stdout);
      GLOBAL_OMEGA = OMEGA0*(k+1);
      //      printf("Doing leastSquares fit\n");
      //      TKleastSquares_svd_psr_dcm(resx,resy,sig,nres,param,error,3,cvm,&chisq,fitMeanSineFunc,0,psr,1.0e-40,ip,uinv);
      TKleastSquares_svd_psr_dcm(resx,resy,sig,nres,param,error,2,cvm,&chisq,fitCosSineFunc,0,psr,1.0e-40,ip,uinv);
      //      printf("Done leastSquares fit\n");
      v[k] = (resx[nres-1]-resx[0])/365.25/2.0/pow(365.25*86400.0,2); 
      specX[k] = GLOBAL_OMEGA/2.0/M_PI;
      //      specY_R[k] = sqrt(v[k])*param[1];
      //      specY_I[k] = sqrt(v[k])*param[2];
            specY_R[k] = sqrt(v[k])*param[0];
            specY_I[k] = sqrt(v[k])*param[1];
    }

  //  for (i=0;i<nfit;i++)
  //    free(cvm[i]);
  //  free(cvm);
  //  printf("Complete spectra\n");
  return nfit;
}

void calculateCovarFunc(double modelAlpha,double modelFc,double modelA,double *covFunc,double *resx,double *resy,double *rese,int np)
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

  ndays = (int)((resx[np-1])-(resx[0])+0.5); 
  f  = (double *)malloc(sizeof(double)*(ndays*2+2)); // Frequency vector
  p  = (double *)malloc(sizeof(double)*(ndays*2+2)); // Model of pulsar power spectrum
  pf = (double *)malloc(sizeof(double)*(ndays*2+2)); // Periodic spectrum model
  opf = (double *)malloc(sizeof(double)*(ndays*2+2)); // Periodic spectrum model
  pe = (double *)malloc(sizeof(double)*(ndays*2+2)); // Periodic spectrum model

  printf("Number of days = %d\n",ndays);


  //  for (i=0;i<ndays+1;i++)
  //  for (i=0;i<ndays;i++)
  for (i=0;i<2*ndays;i++)
    {
      //      f[i] = i*1.0/(resx[np-1]-resx[0])*365.25;
      f[i] = (double)(i+1)/(double)(2*ndays)*365.25;
      //      f[i] = (double)(i+1)/(double)(2*ndays)*365.25;
      
      // CHANGED TO THIS ...
      // BIG CHANGE HERE ........
      p[i] = modelA/pow(1.0+pow(f[i]/modelFc,modelAlpha/2.0),2);
	     //p[i] = 1e-14*1e-14/12.0/M_PI/M_PI*pow(f[i]*365.25,-13.0/3.0);
      pf[i] = p[i];
    }
  //  j = ndays+1;
      ndays *= 2; j = ndays;
  /*  j = ndays;
  //  for (i=ndays;i>0;i--)
  for (i=ndays-1;i>=0;i--)
    {
      // Changed from minus sign
      //f[j] = -f[i];
            f[j] = f[i];
      //      pf[j] = 1.0/(pow((1.0+pow(f[j]/(modelFc*2),2)),modelAlpha/2.0));
      pf[j] = modelA/pow(1.0+pow(f[j]/modelFc,modelAlpha/2.0),2);
      //pf[j] = 1e-14*1e-14/12.0/M_PI/M_PI*pow(f[j]*365.25,-13.0/3.0);
      j++;
    }
  //  ndays = j;
  ndays=j-1; */

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
    double tt;

    output = (fftw_complex*)opf;
    //    transform_plan = fftw_plan_dft_r2c_1d(ndays, pf, output, FFTW_ESTIMATE);
    transform_plan = fftw_plan_dft_r2c_1d(ndays, pf, output, FFTW_ESTIMATE);
    fftw_execute(transform_plan);    
    fftw_destroy_plan(transform_plan);  
    for (i=0;i<ndays;i++) 
      {
	// Xinping does not have the varScale factor
	covFunc[i] = opf[2*i]/ndays*pow(86400.0*365.25,2)*365.25*varScaleFactor;
	//		printf("covFunc: %d %g %g %d\n",i,opf[2*i],opf[2*i+1],j);
      }
    //    tt = covFunc[0];
    //    for (i=0;i<ndays;i++)
    //      covFunc[i] = 0.0;
    //    covFunc[0] = 1e-12;
  }

  if (debug==1)
    {
      FILE *fout;
      fout = fopen("scaleCovar","w");
      for (i=0;i<=j/2;i++)
	fprintf(fout,"%d %g\n",i,covFunc[i]);
      fclose(fout);
    }

  free(p);              
  free(f);
  free(pf);
  free(opf);
  free(pe);
}
