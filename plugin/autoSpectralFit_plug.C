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
#include "fftw3.h"
#include <TKspectrum.h>
#include "T2toolkit.h"
#include "TKfit.h"

using namespace std;

#define MAX_FREQ 10000


void help() /* Display help */
{
}

void get_covFunc_automatic(pulsar *psr, double expSmooth);
int obtainTimingResiduals(pulsar *psr,double *resx,double *resy,double *rese);
void cubicFit(double *resx,double *resy,double *rese,int nres,double *cubicVal,double *cubicErr);
void findSmoothCurve(double *resx,double *resy,double *rese,int nres,double *cubicVal,double *smoothModel,double expSmooth);
void getHighFreqRes(double *resy,double *smoothModel,int nres,double *highFreqRes);
void interpolate(double *resx,double *resy,double *rese,int nres,double *cubicVal,double *interpX,double *interpY,int *nInterp,int interpTime,double expSmooth);
int calculateSpectra(double *x,double *y,double *e,int n,int useErr,int preWhite,int specType,double *specX,double *specY);
void getWhiteRes(double *resx,double *resy,double *rese,int nres,double **uinv,double *cholWhiteY);
int guess_vals(double *x, double *y, int n, double *alpha,  double *fc,  int  *nfit, double wn);
int fitSpectra(double *preWhiteSpecX,double *preWhiteSpecY,int nPreWhiteSpec,double *modelAlpha,double *modelFc,int *modelNfit,double *modelScale,double *fitVar,int aval,int ipw,double ifc,double iexp,int inpt);
void calculateCholesky(double modelAlpha,double modelFc,double modelScale,double fitVar,double **uinv,double *covFunc,double *resx,double *resy,double *rese,int np,double *highFreqRes,double *errorScaleFactor);
void calculateDailyCovariance(double *x,double *y,double *e,int n,double *cv,int *in,double *zl,int usewt);
void formCholeskyMatrix_pl(double *c,double *resx,double *resy,double *rese,int np,double **uinv);
void getWhiteNoiseLevel(int n, double *y, int nlast, double *av);

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i;
  double globalParameter;
  double tsmooth;

  const char *CVS_verNum = "$Revision$";

  if (displayCVSversion == 1) CVSdisplayVersion("grTemplate.C","plugin",CVS_verNum);

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: autoSpectralFit\n");
  printf("Author:              R. Shannon, G. Hobbs\n");
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


  tsmooth = 14;
  get_covFunc_automatic(psr, tsmooth);

  return 0;
}

//This is based on Ryan Shannon's function of the same name.
void get_covFunc_automatic(pulsar *psr, double expSmooth)
{
  int nres;
  double resx[MAX_OBSN], resy[MAX_OBSN], rese[MAX_OBSN];
  char fname[MAX_FILELEN];
  char outname[MAX_FILELEN];

  // allocate the memory for the covariance fucntion
  nres = obtainTimingResiduals(psr,resx,resy,rese);

  int i,j;
  double **uinv;
  double *covFunc;
  int dspan=(int)(resx[nres-1]-resx[0]+0.5);

  

  uinv= (double **)malloc(sizeof(double *)*(nres+1));
  covFunc = (double *)malloc(sizeof(double)*(dspan + 5));
  for (i=0;i<nres+1;i++)uinv[i] = (double *)malloc(sizeof(double)*(nres+1));      


  // start with residuals, calculate Cholesky matrix uinv and covariance function
  // fit cubic to residuals (no interaction required)
  double cubicVal[4], cubicErr[4];
  cubicFit(resx,resy,rese,nres,cubicVal,cubicErr);
  

  // obtain a smooth curve that models the residuals well
  // NEED to choose smoothing time expSmooth! (will this be the same in the simulated and real data?)
  double *smoothModel;
  smoothModel = (double*) malloc(nres*sizeof(double));
  findSmoothCurve(resx,resy,rese,nres,cubicVal,smoothModel,expSmooth);


  // obtain the high frequency residuals (no interaction required)
  double *highFreqRes;
  highFreqRes= (double*) malloc(nres*sizeof(double));
  getHighFreqRes(resy,smoothModel,nres,highFreqRes);


  // calculate power spectra of high frequency residuals
  for (i=0;i<nres;i++)
    {
      for (j=0;j<nres;j++)	
	{
	  if (i==j)
	    uinv[i][j]=1.0/(rese[i]);
	  else
	    uinv[i][j]=0.0;
	}
    }
  int nHighFreqSpec;
  double highFreqSpecX[MAX_FREQ], highFreqSpecY[MAX_FREQ];
  nHighFreqSpec = calcSpectra(uinv,resx,highFreqRes,nres,highFreqSpecX,highFreqSpecY,-1);

    
  // interpolate the smooth curve
  double interpTime = 14.0;
  int nInterp = (int)((resx[nres-1]-resx[0])/interpTime+0.5);
  double *interpX, *interpY;
  interpX = (double*) malloc(nInterp*sizeof(double));
  interpY = (double*) malloc(nInterp*sizeof(double));


  interpolate(resx,resy,rese,nres,cubicVal,interpX,interpY,
	       &nInterp,interpTime,expSmooth);


   // obtain spectra with the correct level of prewhitening
   int usePreWhitening, nPreWhiteSpec;
   double *preWhiteSpecX, *preWhiteSpecY;
   preWhiteSpecX = (double*) malloc(MAX_FREQ*sizeof(double));
   preWhiteSpecY = (double*) malloc(MAX_FREQ*sizeof(double));

   usePreWhitening = 1;  //ask Bill - this seems to still give sensible results - after all, the model is iterated to find a better model, and the model only needs to be good enough to avoid leakage in the spectrum, doesn't have to be exact.
   nPreWhiteSpec = calculateSpectra(interpX,interpY,rese,nInterp,0,usePreWhitening,2,preWhiteSpecX,preWhiteSpecY);

   FILE *specfile;
   /*
   specfile=fopen("comp.dat", "w");
   for(i=0;i<nPreWhiteSpec;i++)
     {
       fprintf(specfile, "%.3le %.3le %.3le %.3le\n", preWhiteSpecX[i], preWhiteSpecY[i] ,highFreqSpecX[i],highFreqSpecY[i]);
     }
   fclose(specfile);
   */

   int nLast=20;
   if (nLast > nHighFreqSpec)
     {
       printf("Error: nLast > nHighFreqSpec because %d > %d, setting nLast = 13\n",nLast,nHighFreqSpec);
       nLast=13;
       if (nLast > nHighFreqSpec)
	 {
	   printf("Error: STILL nLast > nHighFreqSpec because %d > %d\n",nLast,nHighFreqSpec);       
	   exit(0);
	 }
     }


   // Estimate the level of white noise using highest frequency points in power spectrum. This appears to be working
   double whiteNoiseLevel;
   getWhiteNoiseLevel(nHighFreqSpec,highFreqSpecY, nLast, &whiteNoiseLevel);
   printf("In getcovFun_automatic: whiteNoiseLevel = %.3e\n", whiteNoiseLevel);



   // guess correct values for model fit
   double modelFc;
   double modelAlpha;
   int modelNfit;
   guess_vals(preWhiteSpecX, preWhiteSpecY,nPreWhiteSpec, &modelAlpha, &modelFc, &modelNfit, whiteNoiseLevel);   
   modelAlpha = -1.*modelAlpha;
   modelFc =  modelFc*365.25;  //convert modelFc to 1/years


   /*
   //For using a fixed model that resembles the GWB power spectrum:
   modelAlpha = 13.0/3.0;
   modelFc = 365.2425 * 0.7 * preWhiteSpecX[0];
   modelNfit = 5;
   */


   printf("Alpha %.3le Fc  %.3le Nfit %d\n", modelAlpha, modelFc, modelNfit);   
   

   // fit to spectra:  use aval=1 to do this automatically
   // RMS:  these don't do anything because aval=1, but set them to some arbitrary values to be pedantic
   double modelScale, fitVar;
   double ifc=1.0, iexp=1.0;
   int inpt=1, cont;
   cont = fitSpectra(preWhiteSpecX,preWhiteSpecY,nPreWhiteSpec,&modelAlpha,&modelFc,&modelNfit,&modelScale,&fitVar,1,usePreWhitening,ifc, iexp, inpt);

   //For using a fixed spectral model
   //modelScale = gwamp * gwamp / 12.0 / M_PI / M_PI / pow(modelFc,modelAlpha);

   printf("Alpha %.3lg Fc  %.3lg Nfit %d modelScale %g\n", modelAlpha, modelFc, modelNfit, modelScale);


   //PRINT OUT MODELS TO FILE!!!!!!!!!!!!!
   sprintf(fname,"specModel_%s",psr[0].name);
   specfile=fopen(fname, "w");
   for (i=0;i<nPreWhiteSpec;i++)
     {
       fprintf(specfile,"%d %g %g\n", i, preWhiteSpecX[i]*365.25, 1.0 / (pow((1.0 + pow(preWhiteSpecX[i] * 365.25 / (modelFc),2.0)), (modelAlpha) / 2.0)) );   //consider having modelScale as the numerator - NO modelScale hasn't been set properly yet!!! See below "fitSpectra" command...
     }
   fclose(specfile);

   sprintf(fname,"shortSpecModel_%s",psr[0].name);  //this file contains just the model parameters.
   specfile=fopen(fname, "w");
   fprintf(specfile,"Alpha %.3lg Fc  %.3lg Nfit %d modelScale %g\n", modelAlpha, modelFc, modelNfit, modelScale);
   fclose(specfile);

   sprintf(fname,"%s.model",psr[0].name);
   specfile = fopen(fname,"w");
   fprintf(specfile,"MODEL 1\n");
   fprintf(specfile,"ALPHA %g\n",modelAlpha);
   fprintf(specfile,"FC %g\n",modelFc);
   fprintf(specfile,"AMP %g\n",modelScale);
   fprintf(specfile,"WHITENOISE %g\n",whiteNoiseLevel);
   fclose(specfile);




   // calculate the Cholesky matrix
   double errorScaleFactor=1.0;
   calculateCholesky(modelAlpha,modelFc,modelScale,fitVar,uinv,covFunc,resx,resy,rese,nres,highFreqRes,&errorScaleFactor);


   // get the white residuals using the Cholesky matrix
   double cholWhiteY[MAX_OBSN];
   getWhiteRes(resx,resy,rese,nres,uinv,cholWhiteY);

   // get spectrum of whitened data
   int nCholWspec;
   double cholWspecX[MAX_FREQ],cholWspecY[MAX_FREQ]; 
   nCholWspec = calculateSpectra(resx,cholWhiteY,rese,nres,0,0,2,cholWspecX,cholWspecY);

   // get covariance of white residuals
   int whiteCovarNpts[dspan];
   double whiteCovar[dspan];
   double zerolagWhiteCovar;
   calculateDailyCovariance(resx,cholWhiteY,rese,nres,whiteCovar,whiteCovarNpts,&zerolagWhiteCovar,0);

   // improve the spectral estimate
   int nCholSpec;
   double cholSpecX[MAX_FREQ],cholSpecY[MAX_FREQ];
   nCholSpec = calcSpectra(uinv,resx,resy,nres,cholSpecX,cholSpecY,-1);

   double nmodelScale;
   fitSpectra(cholSpecX,cholSpecY,nCholSpec,&modelAlpha,&modelFc,&modelNfit,&nmodelScale,&fitVar,1,usePreWhitening,ifc, iexp, inpt);

   //For using a fixed spectral model that resembles the GWB:
   //modelScale = gwamp * gwamp / 12.0 / M_PI / M_PI / pow(modelFc,modelAlpha);

   // recalculate Cholesky matrix
   calculateCholesky(modelAlpha,modelFc,nmodelScale,fitVar,uinv,covFunc,resx,resy,rese,nres,highFreqRes,&errorScaleFactor);
   

  
   
   // write out covariance function for diagnostic purposes
   FILE *covarfile;
   sprintf(outname,"covarFunc.dat_%s",psr[0].name);
   covarfile=fopen(outname, "w");
   fprintf(covarfile, "%.15g\n", errorScaleFactor);
    for (i=0;i<=(dspan+3);i++)
      {
	fprintf(covarfile,"%.15g\n",covFunc[i]);
      }
   fclose(covarfile);

   // free memory (don't need to keep covFunc)
   free(smoothModel);
   free(highFreqRes);
   for (i=0;i<nres+1;i++) free(uinv[i]);
   free(uinv);
   free(covFunc);
   return;
}


// Based on Ryan Shannon's function of the same name. Fill resx with the SATs (in days), resy with post-fit residuals (s) and rese
// with TOA errors (s)
//
int obtainTimingResiduals(pulsar *psr,double *resx,double *resy,double *rese)
{
  int i;
  int nres=0;

  for (i=0;i<psr[0].nobs;i++)
    {
      if (psr[0].obsn[i].deleted==0)
	{
	  resx[nres] = (double)(psr[0].obsn[i].sat - psr[0].param[param_pepoch].val[0]);
	  // Check to make sure that the data are time sorted
	  if (nres>0 && resx[nres] < resx[nres-1])
	    {
	      printf("ERROR: Data are not time sorted\n");
	      exit(1);
	    }
	 
	  resy[nres] = (double)(psr[0].obsn[i].residual);
	  rese[nres] = (double)(psr[0].obsn[i].toaErr*1.0e-6);
	  //fprintf(stderr, "%.13le %.13le %.13le\n", resx[nres],(double) (psr[1].obsn[i].residual) , rese[nres]);
	  nres++;
	}
    }
  return nres;
}


// Based on Ryan Shannon's function
// NOTE: Bill does an unweighted fit
//
void cubicFit(double *resx,double *resy,double *rese,int nres,double *cubicVal,double *cubicErr)
{
  double **cvm;
  int i;
  int nfit=4; // a+b*x+c*x^2 + d*x^3 
  double chisq;
  int useWeight=0;

  // Allocate memory for fit covariance matrix
  cvm = (double **)malloc(sizeof(double *)*nfit);
  for (i=0;i<nfit;i++) cvm[i] = (double *)malloc(sizeof(double)*nfit);

  TKleastSquares_svd(resx,resy,rese,nres,cubicVal,cubicErr,nfit,cvm,&chisq,TKfitPoly,useWeight);

 
  

  // Free memory allocation
  for (i=0;i<nfit;i++) free(cvm[i]);
  free(cvm);
  
}

//Based on Ryan Shannon's function
void getHighFreqRes(double *resy,double *smoothModel,int nres,double *highFreqRes)
{
  int i;
  for (i=0;i<nres;i++)
    highFreqRes[i] = resy[i] - smoothModel[i];
}

// Based on Ryan Shannon's function
// Exponential smoothing and interpolation
void interpolate(double *resx,double *resy,double *rese,int nres,double *cubicVal,double *interpX,double *interpY,int *nInterp,int interpTime,double expSmooth)
{
  int i,j;
  double tl,bl,dt;
  double yval;
  int writeFiles=0;

  *nInterp = (int)((resx[nres-1]-resx[0])/interpTime+0.5); // Number of interpolated points
  for (i=0;i<*nInterp;i++)
    {
      interpX[i] = resx[0]+i*interpTime;
      tl = 0.0;
      bl = 0.0;
      // Interpolate difference between points and cubic
      for (j=0;j<nres;j++)
	{
	  dt = interpX[i] - resx[j];
	  yval = resy[j]-(cubicVal[0] + cubicVal[1]*resx[j] + 
		  cubicVal[2]*pow(resx[j],2) + cubicVal[3]*pow(resx[j],3));


	  tl += yval/rese[j]/rese[j]*exp(-fabs(dt/expSmooth));
	  bl += 1.0/rese[j]/rese[j]*exp(-fabs(dt/expSmooth));
	}	    
      interpY[i] = tl/bl;      
      // Add on cubic
      interpY[i] += (cubicVal[0] + cubicVal[1]*interpX[i] + 
		     cubicVal[2]*pow(interpX[i],2) + cubicVal[3]*pow(interpX[i],3));
    }
  if(writeFiles)
  {
    FILE *fout;
    fout = fopen("interp.dat","w");
    for (i=0;i<*nInterp;i++)
      {
	fprintf(fout,"%g %g\n",interpX[i],interpY[i]);
      }
    fclose(fout);
    fout = fopen("res.dat","w");
    for (i=0;i<nres;i++)
      fprintf(fout,"%g %g\n",resx[i],resy[i]);
    fclose(fout);
  }
}

int calculateSpectra(double *x,double *y,double *e,int n,int useErr,int preWhite,int specType,double *specX,double *specY)
{
  int nSpec;
  int i;
  int passSpecType;
  double *outY_im, *outY_re;
  outY_im = (double*) malloc(MAX_FREQ*sizeof(double));
  outY_re = (double*) malloc(MAX_FREQ*sizeof(double));

  if (specType==1) // Weighted Lomb Scargle
    passSpecType=4;
  else if (specType==2)
    passSpecType=2;

  //  TKspectrum(x,y,e,n,0,0,0,0,preWhite,passSpecType,1,1,1,specX,specY,&nSpec,0,0,outY_re,outY_im,0);
  TKspectrum(x,y,e,n,0,0,0,0,preWhite,passSpecType,1,1,1,specX,specY,&nSpec,0,0,outY_re,outY_im);
  free(outY_im);
  free(outY_re);
  return nSpec;
}

//Based on Ryan Shannon's function
void getWhiteNoiseLevel(int n, double *y, int nlast, double *av)
{
  int i;
  
  *av =0;
  
  for(i=1;i<=nlast;i++)
    {
      *av += y[n-i];
    }

  *av /=(double) nlast;
}

//Based on Ryan Shannon's function
int guess_vals(double *x, double *y, int n, double *alpha,  double *fc,  int  *nfit, double wn)
{
  double *lx;
  double *ly;
  
  //find how many points in the spectrum are above the white noise level.
  int i=0;

  //IMPORTANTLY this loop will stop AS SOON AS the spectrum dips below the white noise level, meaning we don't include any random noisy points above the white noise level at high frequencies. Rather, we only include the low frequency power excess.
  while(y[i] > wn && i < n)
    {
      i++;
    }

  //and then take the first point of the white noise. UNLESS every point in the spectrum is above the white noise level.
  if (i < n) i++;

  //if there are less than 2 points above the white noise level
  if(i <2)
    {
      fprintf(stderr, "Error, data isn't **red**[?? Ryan had white written here...] enough to calculate covariance function this way. Therefore am auto-fitting to the first 10 points of the spectrum.\n");
      //exit(0);
      i = 10;
    }
  printf("There are %d points above the white noise level in guess_vals\n",i);
  *nfit = i; //DY added the "-1" to avoid problems with nan's where the number of fitted parameters exceeds the length of the string.

  // for all the relevant values get the logs
  lx = (double*) malloc((*nfit)*sizeof(double));
  ly = (double*) malloc((*nfit)*sizeof(double));

  for(i=0;i<*nfit;i++)
    {
      lx[i] = log(x[i]);
      ly[i] = log(y[i]);
    }

  // linear regression;
  double alx, aly;
  double vlx,vly,clxy;
  double tx, ty;

  alx=0;
  aly=0;
  vlx=0;
  vly=0;
  clxy=0;
  for(i=0;i<*nfit;i++)
    {
      alx += lx[i]/((double) *nfit);
      aly += ly[i]/((double) *nfit);
    }
  
  for(i=0;i<*nfit;i++)
    {
      tx = lx[i]-alx;
      ty = ly[i]-aly;
      vlx += tx*tx/((double) *nfit);
      vly += ty*ty/((double) *nfit);
      clxy += tx*ty/((double) *nfit);
    }
  
  // according to wikipedia if y = a+b*x
  // then b = cov(x,y)/var(x)
  // alpha = b is the slope of the line
  

  *alpha = clxy/vlx;
  printf("HERE alpha = %g, clxy = %g, vlx = %g, but it seems that the Cholesky plugin CAN cope with a positive slope (i.e. a negative alpha - recall that alpha goes in the bottom line of the spec Model in the Cholesky plugin).\n",*alpha,clxy,vlx);

  // set fc to the lowest frequency bin
  *fc = x[0];
  
  return 1;
}

// From Ryan Shannon, probably from George Hobbs before that
int fitSpectra(double *preWhiteSpecX,double *preWhiteSpecY,int nPreWhiteSpec,double *modelAlpha,double *modelFc,int *modelNfit,double *modelScale,double *fitVar,int aval,int ipw,double ifc,double iexp,int inpt)
{
  static int time=1;
  double v1,v2,m;
  double df;
  int i;

  if (time==2 && aval==0)
    {
      int redo;
      if (inpt == -1)
	{
	  printf("Redo fit (1 = yes, 0 = no)\n");
	  scanf("%d",&redo);
	  if (redo==0)
	    return 0;
	}
      else
	return 0;
    }
  time=2;
  if (aval==0)
    {
      if (ifc == -1) {printf("Enter corner freq (yr-1) "); scanf("%lf",modelFc);}
      else *modelFc = ifc;
      if (iexp == 0) {printf("Enter power law exponential (should be positive) "); scanf("%lf",modelAlpha);}
      else *modelAlpha = iexp;	
      if (inpt == -1) {printf("Enter nfit "); scanf("%d",modelNfit);}
      else *modelNfit = inpt;
    }
  // Do the fit
  // This fit is useful for fitting spectra where the error is proportional to the mean
  // The error on each point is simply taken as the model value squared and so we solve for
  // chisq = sum (P_d(f) - aM(f))^2/M^2(f)
  // which simplifies to a simple formula
  v1 = 0.0;
  for (i=0;i<*modelNfit;i++)
    {
      m = 1.0/pow((1.0+pow(preWhiteSpecX[i]*365.25/(*modelFc),2)),(*modelAlpha)/2.0);
      //m = pow(preWhiteSpecX[i]*365.25,(-1.* *modelAlpha));
      v1 += preWhiteSpecY[i]/m;
    }
  //  *modelScale = log10(v1/(double)(*modelNfit));
  *modelScale = (v1/(double)(*modelNfit));
  printf("Model scale = %g\n",*modelScale);

  // Get area under the spectra
  *fitVar=0.0;
  df = preWhiteSpecX[0]*365.25;
  // Changed based on Bill's advice
  for (i=0;i<nPreWhiteSpec;i++)
    {
      //  for (i=0;i<*modelNfit;i++)
      //      *fitVar+=preWhiteSpecY[i]*df;
      m = 1.0/pow((1.0+pow(preWhiteSpecX[i]*365.25/(*modelFc),2)),(*modelAlpha)/2.0);
      *fitVar+=(*modelScale*m*df);
    }
  (*fitVar)*=pow(86400.0*365.25,2);
  return 1;
}


// From Ryan Shannon. Prob based on George Hobbs's code.
void calculateCholesky(double modelAlpha,double modelFc,double modelScale,double fitVar,double **uinv,double *covFunc,double *resx,double *resy,double *rese,int np,double *highFreqRes,double *errorScaleFactor)
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
  mean=0.0;
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

  // Weighted variance of residuals
  weightVarRes=0.0;
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
  weightVarHighFreqRes = tl2/bl;

  for (i=0;i<ndays;i++)
    {
      // Update based on Bill's advice
      //      f[i] = i*1.0/(resx[np-1]-resx[0])*365.25;
      f[i] = i*0.5/(resx[np-1]-resx[0])*365.25;
      
      p[i] = 1.0/(pow((1.0+pow(f[i]/(modelFc*2),2)),(modelAlpha)/2.0));
      pf[i] = p[i];
    }
  j = ndays;
  for (i=ndays-1;i>0;i--)
    {
      f[j] = -f[i];
      pf[j] = 1.0/(pow((1.0+pow(f[j]/(modelFc*2),2)),modelAlpha/2.0));
      j++;
    }
  if (debug==1)
    {
      FILE *fout;
      fout = fopen("specModel","w");
      for (i=0;i<j;i++)
	fprintf(fout,"%d %g %g\n",i,f[i],pf[i]);
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
	//printf("covFunc: %g %g\n",opf[2*i],opf[2*i+1]);
      }
  }
  // Rescale
  printf("Rescaling %d\n",j/2);
  //  for (i=0;i<np/2;i++)
  //  fy2[i] = nmodelScale-log10(pow((1.0+pow(cholSpecX[i]*365.25/modelFc,2)),modelAlpha/2.0));
  actVar = pow(10,modelScale);
  printf("modelScale = %g, actVar = %g\n",modelScale,actVar);
  tt = covFunc[0];
  printf("actVar*(numsec in year)^2 = %g, weightVarRes = %g, weightVarHighFreqRes = %g, scale = %g, fitVar = %g, tt = %g\n",actVar*pow(86400.0*365.25,2),weightVarRes,weightVarHighFreqRes,weightVarRes-weightVarHighFreqRes,fitVar,tt);

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
  formCholeskyMatrix_pl(covFunc,resx,resy,rese,np,uinv);

  free(p);              
  free(f);
  free(pf);
  free(opf);
  free(pe);
}


void getWhiteRes(double *resx,double *resy,double *rese,int nres,double **uinv,double *cholWhiteY)
{
  int i,j;
  double sum;
  printf("Getting white residuals\n");
  for (i=0;i<nres;i++)
    {
      sum=0.0;
      for (j=0;j<nres;j++)
	{
	  sum+=uinv[j][i]*resy[j];
	  //	  if (i==0) printf("uinv = %g\n",uinv[j+1][i+1]);
	}
      cholWhiteY[i]=sum;
    }
  //  exit(1);
}

void calculateDailyCovariance(double *x,double *y,double *e,int n,double *cv,int *in,double *zl,int usewt)
{
  int i,j;
  int nd = (int)(x[n-1]-x[0]+0.5);
  int nc=0;
  int nzl=0;
  int npts[nd];
  double dt;
  double wt[nd];
  double nzerolag=0.0;

  for (i=0;i<nd;i++)
    {
      cv[i] = 0.0;
      in[i] = 0;
      wt[i] = 0.0;
    }
  *zl = 0.0;

 // Bin in 1 day intervals
  for (i=0;i<n;i++)
    {
      for (j=i;j<n;j++)
	{
	  dt = fabs(x[i]-x[j]);
	  if (dt == 0) // Absolute zero time delay
	    {
	      (*zl)+=(y[i]*y[j]);
	      nzl++;
	    }
	  else 
	    {
	      if (usewt==0)
		cv[(int)(dt+0.5)]+=(y[i]*y[j]);
	      else
		cv[(int)(dt+0.5)]+=(y[i]*y[j])/(e[i]*e[j]);
	      in[(int)(dt+0.5)]++;
	      wt[(int)(dt+0.5)]+=1.0/(e[i]*e[j]);
	    }
	}
    }
  printf("zl = %g\n",*zl);
  (*zl)/=(double)nzl;
  for (i=0;i<nd;i++)
    {
      if (in[i] > 0)
	{
	  if (usewt==0)
	    cv[i]/=(double)in[i];
	  else
	    cv[i]/=(double)wt[i];
	}
    }

}


void findSmoothCurve(double *resx,double *resy,double *rese,int nres,double *cubicVal,double *smoothModel,double expSmooth)
{
  int i,j;
  double tl,bl;
  double diffy[nres];
  double dt;

  // Obtain difference between residuals and the cubic fit
  for (i=0;i<nres;i++)
    {
      diffy[i] = resy[i] - (cubicVal[0] + cubicVal[1]*resx[i] + 
			    cubicVal[2]*pow(resx[i],2) + cubicVal[3]*pow(resx[i],3));
    }
  // Smooth to get model of the timing residuals
  for (i=0;i<nres;i++)
    {
      tl = 0.0;
      bl = 0.0;
      for (j=0;j<nres;j++)
	{
	  dt = resx[i] - resx[j];
	  tl += diffy[j]*exp(-fabs(dt/expSmooth))/rese[j]/rese[j];
	  bl += 1.0*exp(-fabs(dt/expSmooth))/rese[j]/rese[j];
	  //	  tl += diffy[j]*exp(-fabs(dt/expSmooth));
	  //	  bl += 1.0*exp(-fabs(dt/expSmooth));
	}
      smoothModel[i] = tl/bl + (cubicVal[0] + cubicVal[1]*resx[i] + 
				cubicVal[2]*pow(resx[i],2) + 
				cubicVal[3]*pow(resx[i],3));
    }
}


char * plugVersionCheck = TEMPO2_h_VER;

void formCholeskyMatrix_pl(double *c,double *resx,double *resy,double *rese,int np,double **uinv)
{
  double **m,**u,sum;
  double *cholp;
  int i,j,k,ix,iy;
  double t0,cint,t;
  int t1,t2;
  int debug;

  printf("Getting the covariance matrix\n");
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
  // Insert the covariance which depends only on the time difference. Linearly interpolate between elements on the covariance function because valid covariance matrix must have decreasing off diagonal elements.
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
