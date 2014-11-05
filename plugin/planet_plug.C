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

/* Template for a tempo2 plugin */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include <cpgplot.h>
#include "T2toolkit.h"
#include "TKspectrum.h"
#include "TKfit.h"
#include "choleskyRoutines.h"
#include "cholesky.h"

using namespace std;


char pgdevice[80];
double G_OMEGA;
void plot6(double *cholSpecX,double *cholSpecY,int nCholSpec,double *cholWspecX,
	   double *cholWspecY,int nCholWspec,double *highFreqSpecX,
	   double *highFreqSpecY,int nHighFreqSpec,int makeps);
void doPlugin(pulsar *psr,double idt,int ipw,double ifc,double iexp,int inpt,int makeps,double amp,char* dcf_file, int *npsr,  char *argv[], int argc, char parFile[][MAX_FILELEN], char timFile[][MAX_FILELEN], int nit );
int obtainTimingResiduals(pulsar *psr,double *resx,double *resy,double *rese,int *ip);
void fitSineFunc(double x,double *v,int nfit,pulsar *psr,int ival);
void plot1(double *resx,double *resy,double *rese,int nres,double *cubicVal,double *smoothModel,double *highFreqRes,double *hfNormCovar,int *hfNormCovarNpts,double hfZerolagNormCovar);
void removeMean(double *resx,double *resy,int n);
void fileOutput3(char *fname,double *x,double *y,double *z,int n);
void fileOutput2(char *fname,double *x,double *y,int n);
void findSmoothCurve(double *resx,double *resy,double *rese,
		     int nres,double *cubicVal,double *smoothModel,double expSmooth);
void getHighFreqRes(double *resy,double *smoothModel,int nres,double *highFreqRes);
void getHighFreqCovar(double *resx,double *rese,double *highFreqRes,int nres,double *hfNormCovar,int *hfNormCovarNpts,double *hfZerolagNormCovar);
void calculateDailyCovariance(double *x,double *y,double *e,int n,double *cv,int *in,double *zl,int usew);
int calculateSpectra(double *x,double *y,double *e,int n,int useErr,int preWhite,
		     int specType,double *specX,double *specY);
void plot2(double *origSpecX,double *origSpecY,int nOrigSpec,double *smoothSpecX0,
	   double *smoothSpecY0,int nSmoothSpec0,double *smoothSpecX1,
	   double *smoothSpecY1,int nSmoothSpec1,double *smoothSpecX2,
	   double *smoothSpecY2,int nSmoothSpec2,double *highFreqSpecX,
	   double *highFreqSpecY,int nHighFreqSpec,int makeps);
void plot3(double *preWhiteSpecX,double *preWhiteSpecY,int nPreWhiteSpec,
	   int usePreWhitening,double *highFreqSpecX,double *highFreqSpecY,
	   int nHighFreqSpec,double modelAlpha,double modelFc,int modelNfit,double modelScale,int closeit,float *minx,float *maxx, float wn);
void plot3a(double *resx,double *resy,int nres,double *rawCovar,int *rawCovarNpts,
	    double zerolagRawCovar,double *ampFit,double *chisqFit,int nGridFit,
	    double bestAmp,double bestLag,double bestChisq,int makeps);

void plot4(double *resx,double *resy,double *rese,int nres,double *cholWhiteY,double *whiteCovar,int *whiteCovarNpts,double zerolagWhiteCovar);
void plot5(double *preWhiteSpecX,double *preWhiteSpecY,int nPreWhiteSpec,
	   int usePreWhitening,double *highFreqSpecX,double *highFreqSpecY,int nHighFreqSpec,
	   double modelAlpha,double modelFc,int modelNfit,double modelScale,
	   double nmodelScale,double *cholSpecX,double *cholSpecY,int nCholSpec,
	   double *cholWspecX,double *cholWspecY,int nCholWspec,int makeps, double wn, double pb);
void outputMatrix(double **uinv,int nres);
void fitExponential(double *resx,int nres,double *rawCovar,int *rawCovarNpts,double *ampFit,double *chisqFit,double *bestAmp,double *bestLag,double *bestChisq,int *nGridFit);
void calculateCholeskyCovarFunc(double bestAmp,double bestLag,int nGridFit,double **uinv,double *resx,
				double *resy,double *rese,int nres,double *covarFunc);
void outputCovarianceFunction(double *covFunc,int n,double errorScaleFactor,pulsar *psr);





double modelfcn(double freq, double nmodelScale, double modelFc, double modelAlpha, double wn);


int T2fitSpectraRMS(double *preWhiteSpecX,double *preWhiteSpecY,int nPreWhiteSpec,double *modelAlpha,double *modelFc,int *modelNfit,double *modelScale,double *fitVar,int aval,int ipw,double ifc,double iexp,int inpt,double amp, double *wn)
{
  static int time=1;
  double v1,v2,m;
  double df;
  
  int i;
  printf("choleskyRoutines: fitSpectra RMS\n");
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
      
      printf("Enter white noise level "); scanf("%lf",wn);
 
       *wn = pow(10., *wn);
  }

 
  // Do the fit
  // This fit is useful for fitting spectra where the error is proportional to the mean
  // The error on each point is simply taken as the model value squared and so we solve for
  // chisq = sum (P_d(f) - aM(f))^2/M^2(f)
  // which simplifies to a simple formula
  printf("Got here with amp = %g\n",amp);
  if (amp == -1)
    {
      v1 = 0.0;
      for (i=0;i<*modelNfit;i++)
	{
	  //	  m = 1.0/pow((1.0+pow(preWhiteSpecX[i]*365.25/(*modelFc),2)),(*modelAlpha)/2.0);
//	  m = 1.0/pow(1.0+pow(preWhiteSpecX[i]*365.25/(*modelFc),*modelAlpha/2.0),2); // MJK OLD CODE

	  m=1.0/pow(1.0+pow(fabs(preWhiteSpecX[i]*365.25)/(*modelFc),2),*modelAlpha/2.0);

	  v1 += preWhiteSpecY[i]/m;
	  printf("Here with %g %g %g %d\n",v1,preWhiteSpecY[i],m,*modelNfit);

	}
      //  *modelScale = log10(v1/(double)(*modelNfit));
      *modelScale = (v1/(double)(*modelNfit));
    }
  else
    *modelScale = amp;
  printf("Model scale = %g %g %g %d\n",*modelScale,v1,m,*modelNfit);

  // Get area under the spectra
  *fitVar=0.0;
  df = preWhiteSpecX[0]*365.25;
    for (i=0;i<*modelNfit;i++)
      //for (i=0;i<nPreWhiteSpec;i++)
    {
      *fitVar+=preWhiteSpecY[i]*df;
      // CHANGED HERE ...
      //      m = 1.0/pow((1.0+pow(preWhiteSpecX[i]*365.25/(*modelFc),2)),(*modelAlpha)/2.0);
      //*fitVar+=(*modelScale*m*df);
    }
  (*fitVar)*=pow(86400.0*365.25,2);
  return 1;
}


void help() /* Display help */
{
  /* This function should contain usage information about the plugin which should (in general) be accessed */
  /* by the user pressing 'h'                                                                              */
}

char skipstep2=0; // test to skip step 2 added by MJK 2011-07.
bool writeFiles=true;
int skipprocess=0;

/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  char dcf_file[MAX_FILELEN];
  char covarFuncFile[MAX_FILELEN];
  int i;
  double globalParameter;
  double idt=0;
  double ifc=-1;
  double iexp=0;
  double amp=-1;
  int inpt=-1;
  int ipw=-1;
  int makeps=0;
  int nit=20;

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: planet plugin\n");
  printf("Author:              Ryan Shannon\n");
  printf("Version:             1.1\n");
  printf("The techniques used here based on spectralModel by  Coles and Hobbs\n");
  printf(" --- type 'h' for help information\n");

  dcf_file[0]='\0';
  strcpy(covarFuncFile,"NULL");

  /* Obtain the .par and the .tim file from the command line */
  if (argc==4) /* Only provided .tim name */
    {
      strcpy(timFile[0],argv[3]);
      strcpy(parFile[0],argv[3]);
      parFile[0][strlen(parFile[0])-3] = '\0';
      strcat(parFile[0],"par");
    }

  strcpy(pgdevice,"97/xs");
  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[0],argv[i+1]); 
	  strcpy(timFile[0],argv[i+2]);
	}
      else if (strcmp(argv[i],"-t")==0)
	sscanf(argv[++i],"%lf",&idt);
      else if (strcmp(argv[i],"-g")==0)
	strcpy(pgdevice,argv[++i]);
      else if (strcmp(argv[i],"-dcf")==0){
	strcpy(dcf_file,argv[++i]);
	strcpy(covarFuncFile,dcf_file);
	  }
      else if (strcmp(argv[i],"-fc")==0)
	sscanf(argv[++i],"%lf",&ifc);
      else if (strcmp(argv[i],"-exp")==0)
	sscanf(argv[++i],"%lf",&iexp);
      else if (strcmp(argv[i],"-pw")==0)
	sscanf(argv[++i],"%d",&ipw);
      else if (strcmp(argv[i],"-setamp")==0)
	sscanf(argv[++i],"%lf",&amp);
      else if (strcmp(argv[i],"-nfit")==0)
	sscanf(argv[++i],"%d",&inpt);
      else if (strcmp(argv[i],"-skipstep2")==0)
	skipstep2=1;
      else if (strcmp(argv[i],"-skipprocess")==0)
	skipprocess=1;
      else if (strcmp(argv[i],"-makeps")==0)
	makeps=1;
      else if (strcmp(argv[i],"-nofiles")==0)
	writeFiles=false;
      else if (strcmp(argv[i], "-nit")==0)
	{
	  sscanf(argv[++i],"%d",&nit);
	}

    }

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      if (i==0) {
		 doFitAll(psr,*npsr,covarFuncFile);
	  }
      else textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
    }

  //
  printf("Plugin to obtain a spectral model of pulsar timing residuals\n");
  doPlugin(psr,idt,ipw,ifc,iexp,inpt,makeps,amp,dcf_file, npsr, argv, argc, parFile, timFile, nit);

  return 0; 
} 

void doPlugin(pulsar *psr,double idt,int ipw,double ifc,double iexp,int inpt,int makeps,double amp,char* dcf_file, int *npsr, char *argv[], int argc, char parFile[][MAX_FILELEN], char timFile[][MAX_FILELEN], int nit)
{
  int i,j;

  double resx[MAX_OBSN],resy[MAX_OBSN],rese[MAX_OBSN]; // Timing residuals
  int ip[MAX_OBSN];
  int    nres;                                         // Number of timing residuals
  double cubicVal[4],cubicErr[4];                      // Cubic fit to timing residuals
  double smoothModel[MAX_OBSN];                        // Smooth curve fitting the residuals
  double highFreqRes[MAX_OBSN];                        // High frequency residuals
  double expSmooth;                                    // Smoothing timescale
  double hfNormCovar[MAX_OBSN];                        // High frequency normalised covariance
  double hfZerolagNormCovar;                           // Zerolag covariance of HF residuals
  int hfNormCovarNpts[MAX_OBSN];                       // Number of covariances measured
  //
  double interpX[MAX_OBSN],interpY[MAX_OBSN];          // Interpolation of smooth curve
  int    nInterp;
  double origSpecX[MAX_OBSN],origSpecY[MAX_OBSN];      // Spectrum of raw data
  int nOrigSpec;
  double smoothSpecX0[MAX_OBSN],smoothSpecY0[MAX_OBSN];// 0-whitening of smooth curve
  int nSmoothSpec0;
  double smoothSpecX1[MAX_OBSN],smoothSpecY1[MAX_OBSN];// 1-whitening of smooth curve
  int nSmoothSpec1;
  double smoothSpecX2[MAX_OBSN],smoothSpecY2[MAX_OBSN];// 2-whitening of smooth curve
  int nSmoothSpec2;
  double highFreqSpecX[MAX_OBSN],highFreqSpecY[MAX_OBSN];// Spectrum of HF residuals
  int nHighFreqSpec;
  int usePreWhitening=0;                               // Prewhitening to use
  double preWhiteSpecX[MAX_OBSN],preWhiteSpecY[MAX_OBSN]; // Spectrum with final prewhitening
  double *covFunc;                                     // Covariance function
  int nPreWhiteSpec;
  double modelAlpha,modelFc,modelScale,nmodelScale;    // Model parameters
  double wn;
  int modelNfit=-1;
  double **uinv;                                       // Whitening matrix
  double **uinvI;
  double cholWhiteY[MAX_OBSN];                         // Cholesky white residuals
  double cholWspecX[MAX_OBSN],cholWspecY[MAX_OBSN];    // Cholesky whitened spectrum
  int    nCholWspec;
  double fitVar;

  int interpTime=14;                                   // Grid for interpolated values (days)
  int cont=1;
  double nSmooth;
  double zerolagRawCovar;
  int rawCovarNpts[MAX_OBSN];
  double rawCovar[MAX_OBSN];
  double zerolagWhiteCovar;
  int whiteCovarNpts[MAX_OBSN];
  double whiteCovar[MAX_OBSN];
  double cholSpecX[MAX_OBSN],cholSpecY[MAX_OBSN],cholSpecE[MAX_OBSN];
  double cholSpecXp[MAX_OBSN],cholSpecYp[MAX_OBSN],cholSpecEp[MAX_OBSN];
  int nCholSpec;
  double ampFit[MAX_OBSN],chisqFit[MAX_OBSN],bestAmp,bestLag,bestChisq;
  int nGridFit;
  double errorScaleFactor = 1;
  int tempTime=1;
  char dummy[100];

  

  verbose_calc_spectra=true;

  // Set up some defaults
  if (idt!=0) expSmooth = idt; 
  else expSmooth = 20;
  // Step 1a: obtain time sorted post-fit residuals
  nres = obtainTimingResiduals(psr,resx,resy,rese,ip);
  // Step 1b: remove mean from residuals (x and y)
  removeMean(resx,resy,nres);
  if(writeFiles)fileOutput3("tresiduals.dat",resx,resy,rese,nres);

	 // Step 1
	 do {
	   
	   if (tempTime==1)
		 {
	   uinv=malloc_uinv(nres);
	   covFunc = (double *)malloc(sizeof(double)*((int)(resx[nres-1]-resx[0])+5));
	   tempTime=2;
		 }
	   // Step 1c: fit a cubic to the timing residuals
	   T2cubicFit(resx,resy,rese,nres,cubicVal,cubicErr);
	   // Step 1d: obtain a smooth curve that models the residuals well
	   T2findSmoothCurve(resx,resy,rese,nres,cubicVal,smoothModel,expSmooth);
	   if(writeFiles)fileOutput3("smoothCurve.dat",resx,smoothModel,rese,nres);
	   // Step 1e: obtain high-freq. residuals
	   T2getHighFreqRes(resy,smoothModel,nres,highFreqRes);
	   if(writeFiles)fileOutput3("highFreqRes.dat",resx,highFreqRes,rese,nres);
	   // Step 1f: obtain covariance of high-freq. residuals
	   getHighFreqCovar(resx,rese,highFreqRes,nres,hfNormCovar,hfNormCovarNpts,&hfZerolagNormCovar);
	   // Step 1g: plot the results
	   if (skipprocess==0)
		 {
	   plot1(resx,resy,rese,nres,cubicVal,smoothModel,highFreqRes,hfNormCovar,hfNormCovarNpts,hfZerolagNormCovar);
		 }
	   if (idt!=0) cont=0;
	   else {
		 printf("smooth = %g days. Press 0 to continue, or another number to update the smooth\n",expSmooth);
		 scanf("%lf",&nSmooth);
		 if (nSmooth==0) cont=0;
		 else expSmooth = nSmooth;
	   }
	 } while (cont==1);

  if(dcf_file[0]=='\0'){

	 if(!skipstep2){
	 // Step 2:

	 // Put errors into uinv matrix
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

	 // Step 2a: Obtain spectra of original residuals without any prewhitening
	 //  nOrigSpec = calculateSpectra(resx,resy,rese,nres,1,0,1,origSpecX,origSpecY);
	 if (skipprocess==0)
	   {
		 nOrigSpec = calcSpectra(uinv,resx,resy,nres,origSpecX,origSpecY,-1);
		 if(writeFiles)fileOutput2("origSpectra.dat",origSpecX,origSpecY,nOrigSpec);

		 // Step 2b: interpolate the smooth curve
		 T2interpolate(resx,resy,rese,nres,cubicVal,interpX,interpY,
			 &nInterp,interpTime,expSmooth);
	   }
	 // Step 2c: Obtain spectra of smooth interpolated model without any prewhitening
	 // Errors are ignored
	 //  nSmoothSpec0 = calculateSpectra(interpX,interpY,rese,nInterp,0,0,2,
	 //  				  smoothSpecX0,smoothSpecY0);
/*	 uinvI = malloc_uinv(nInterp);
	 for (i=0;i<nInterp;i++)
	   {
		 for (j=0;j<nInterp;j++)
	   {
		 if (i==j) uinvI[i][j]=1.0;
		 else	    uinvI[i][j]=0.0;
	   }
	   }*/
	 if (skipprocess==0)
	   {
		 printf("Calculating spectra without prewhitening\n");
		 nSmoothSpec0 = T2calculateSpectra(interpX,interpY,rese,nInterp,0,0,2,
						 smoothSpecX0,smoothSpecY0);
		 if(writeFiles)fileOutput2("zeroprewhite.dat",smoothSpecX0,smoothSpecY0,nSmoothSpec0);

		 //  nSmoothSpec0 = calcSpectra(uinvI,interpX,interpY,nInterp,smoothSpecX0,smoothSpecY0,-1);
		 printf("Done calculating spectra\n");
	   }
		 if(writeFiles)fileOutput2("zeroprewhite.dat",smoothSpecX0,smoothSpecY0,nSmoothSpec0);
	 // TESTING
	 //  nSmoothSpec0 = calcSpectra(uinv,resx,resy,rese,nres,
	 //			     smoothSpecX0,smoothSpecY0);

	 // Step 2d: Obtain spectra of smooth interpolated model with 1st order prewhitening
		 if (skipprocess==0)
	   {
		 nSmoothSpec1 = T2calculateSpectra(interpX,interpY,rese,nInterp,0,1,2,
						 smoothSpecX1,smoothSpecY1);
		 if(writeFiles)fileOutput2("oneprewhite.dat",smoothSpecX1,smoothSpecY1,nSmoothSpec1);
		 // Step 2e: Obtain spectra of smooth interpolated model with 2nd order prewhitening
		 nSmoothSpec2 = T2calculateSpectra(interpX,interpY,rese,nInterp,0,2,2,
						 smoothSpecX2,smoothSpecY2);
		 if(writeFiles)fileOutput2("twoprewhite.dat",smoothSpecX2,smoothSpecY2,nSmoothSpec2);
		 // Step 2f: Obtain spectra of high frequency residuals
	   }
	 if(writeFiles)
	 {
	   long seed= -123;
	   FILE *fout = fopen("highfreqresiduals.dat","w");
	   for (i=0;i<nres;i++)
		 {
	   //highFreqRes[i] = TKgaussDev(&seed)*rese[i];
	   fprintf(fout,"%g %g %g\n",resx[i],highFreqRes[i],rese[i]);
		 }
	   fclose(fout);
	 }
	 //  nHighFreqSpec = calculateSpectra(resx,highFreqRes,rese,nres,1,0,1,highFreqSpecX,
	 //				   highFreqSpecY);
	 if (skipprocess==0)
	   {
		 nHighFreqSpec = calcSpectra(uinv,resx,highFreqRes,nres,highFreqSpecX,highFreqSpecY,-1);
		 if(writeFiles)fileOutput2("highfreqspec.dat",highFreqSpecX,highFreqSpecY,nHighFreqSpec);
		 
		 // Step 2g: make the plot
		 plot2(origSpecX,origSpecY,nOrigSpec,smoothSpecX0,smoothSpecY0,nSmoothSpec0,
		   smoothSpecX1,smoothSpecY1,nSmoothSpec1,smoothSpecX1,smoothSpecY2,nSmoothSpec2,
		   highFreqSpecX,highFreqSpecY,nHighFreqSpec,makeps);
		 
	   }
	 }
	 // Step 3: select a prewhitening
	 if (ipw==-1){
	   printf("Select prewhitening required (0,1,2) (type -1 to obtain covariance function from the data) ");
	   scanf("%d",&usePreWhitening);
	 }
	 else
	   usePreWhitening=ipw;

	 if (usePreWhitening!=-1) // If we're calculating the covariance function from a spectrum
	   {
		 float mx,my;
		 // Step 3a: obtain spectra with this prewhitening
		 if (skipprocess==0)
	   {
		 nPreWhiteSpec = T2calculateSpectra(interpX,interpY,rese,nInterp,0,usePreWhitening,2,
						  preWhiteSpecX,preWhiteSpecY);
	   }
		 // Step 3b: plot spectra
		 cont=1;
		 //      if (amp==-1)
		 {
	   do {
		 if (skipprocess==0)
		   {
			 plot3(preWhiteSpecX,preWhiteSpecY,nPreWhiteSpec,usePreWhitening,
			       highFreqSpecX,highFreqSpecY,nHighFreqSpec,modelAlpha,modelFc,modelNfit,modelScale,1,&mx,&my, wn);
		   }
		 // Step 3c: Fit to spectra
		 cont = T2fitSpectraRMS(preWhiteSpecX,preWhiteSpecY,nPreWhiteSpec,&modelAlpha,&modelFc,&modelNfit,&modelScale,&fitVar,0,ipw,ifc, iexp, inpt,amp, &wn);
		 printf("modelScale = %g\n",modelScale);
	   } while (cont==1);
		 }
	   //      else
	   //	{
	   //	  modelScale = amp;
	   //	}

		 // Step 4a: calculate the Cholesky whitening matrix (uinv)
		 T2calculateCholesky(modelAlpha,modelFc,modelScale,fitVar,uinv,covFunc,resx,resy,rese,nres,highFreqRes,&errorScaleFactor,0,0,0);
		 for (i=0;i<100;i++)
	   printf("cov: %lg %g %g %g %d %g %g\n",covFunc[i],resx[i],resy[i],rese[i],nres,highFreqRes[i],errorScaleFactor);
	   }
	 else // Calculate covariance function from the data
	   {
		 double tt;
		 // Step 3a: calculate covariance function of the raw data
		 printf("Step 1\n");
		 calculateDailyCovariance(resx,resy,rese,nres,rawCovar,rawCovarNpts,&zerolagRawCovar,1);
		 //      getHighFreqCovar(resx,rese,resy,nres,rawCovar,rawCovarNpts,&zerolagRawCovar);
		 // Step 3b: fit for an exponential function
		 printf("Step 2 %g\n",zerolagRawCovar);
		 fitExponential(resx,nres,rawCovar,rawCovarNpts,ampFit,chisqFit,&bestAmp,&bestLag,&bestChisq,&nGridFit);
		 printf("Do plot\n");
		 do {
	   plot3a(resx,resy,nres,rawCovar,rawCovarNpts,zerolagRawCovar,ampFit,chisqFit,nGridFit,bestAmp,bestLag,bestChisq,makeps);
	   printf("Chosen lag = %g (press '-1' to continue or type in a new lag) ",bestLag); scanf("%lf",&tt);
	   if (tt!=-1) {
		 bestLag=tt;
		 // Find closest
		 bestAmp = ampFit[(int)(tt+0.5)];
		 bestChisq = chisqFit[(int)(tt+0.5)];
	   }
		 } while (tt!=-1);
		 calculateCholeskyCovarFunc(bestAmp,bestLag,nGridFit,uinv,resx,resy,rese,nres,covFunc);
		 //      exit(1);
	   }

   } else {
	// Alternaitve to first steps - use a DCF to get the spectrum.
	
		 float mx,my;
	  uinv=malloc_uinv(psr->nobs);
	  getCholeskyMatrix(uinv,covarFuncFile,psr,resx,resy,rese,psr->nobs,0,ip);
		 nOrigSpec = calcSpectra(uinv,resx,resy,nres,origSpecX,origSpecY,-1);
		 usePreWhitening=4;

		 nHighFreqSpec = calcSpectra(uinv,resx,highFreqRes,nres,highFreqSpecX,highFreqSpecY,-1);
		 if(writeFiles)fileOutput2("highfreqspec.dat",highFreqSpecX,highFreqSpecY,nHighFreqSpec);

	   memcpy(preWhiteSpecX,origSpecX,nOrigSpec*sizeof(double));
	   memcpy(preWhiteSpecY,origSpecY,nOrigSpec*sizeof(double));
	   nPreWhiteSpec=nOrigSpec;
	   do {
			  plot3(preWhiteSpecX,preWhiteSpecY,nPreWhiteSpec,usePreWhitening,
				highFreqSpecX,highFreqSpecY,nHighFreqSpec,modelAlpha,modelFc,modelNfit,modelScale,1,&mx,&my, wn);
		 // Step 3c: Fit to spectra
			  cont = T2fitSpectraRMS(origSpecX,origSpecY,nOrigSpec,&modelAlpha,&modelFc,&modelNfit,&modelScale,&fitVar,0,ipw,ifc, iexp, inpt,amp, &wn);
		 printf("modelScale = %g\n",modelScale);
	   } while (cont==1);

		 // Step 4a: calculate the Cholesky whitening matrix (uinv)
	   T2calculateCholesky(modelAlpha,modelFc,modelScale,fitVar,uinv,covFunc,resx,resy,rese,nres,highFreqRes,&errorScaleFactor,0,0,0);

   }

 
  //  fileOutput2("cholWhiteSpec.dat",cholWspecX,cholWspecY,nCholWspec);
  // Step 4e: get covariance of white residuals
  T2calculateDailyCovariance(resx,cholWhiteY,rese,nres,whiteCovar,whiteCovarNpts,&zerolagWhiteCovar,0);
  // Step 4f: plot
  plot4(resx,resy,rese,nres,cholWhiteY,whiteCovar,whiteCovarNpts,zerolagWhiteCovar);
  if (inpt==-1) {printf("Continue (press '1') "); scanf("%s",dummy);}

  // Step 5: now improve the spectral estimate
  // Step 5a: Obtain a spectrum with the whitening routine
  printf("Forming spectrum of whitened data\n");

  /*  {
    FILE *fin;
    fin = fopen("/u/hob044/uinv1.dat","r");
    for (i=0;i<nres;i++)
      {
	for (j=0;j<nres;j++)
	  fscanf(fin,"%lf",&uinv[i][j]);
      }
    
    fclose(fin);
    
    }*/

  nCholSpec = calcSpectraErr(uinv,resx,resy,nres,cholSpecX,cholSpecY,cholSpecE,-1);
  fileOutput3("cholSpectra.dat",cholSpecX,cholSpecY,cholSpecE,nCholSpec);

  // Step 5b: refit the model
  if (usePreWhitening!=-1) // If we're calculating the covariance function from a spectrum
    {
      FILE *fout;
      char fname[128];

      T2fitSpectraRMS(cholSpecX,cholSpecY,nCholSpec,&modelAlpha,&modelFc,&modelNfit,&nmodelScale,&fitVar,1,ipw,ifc, iexp, inpt,amp, &wn);
      sprintf(fname,"%s.model",psr[0].name);
      fout = fopen(fname,"w");
      fprintf(fout,"MODEL 1\n");
      fprintf(fout,"ALPHA %g\n",modelAlpha);
      fprintf(fout,"FC %g\n",modelFc);
      fprintf(fout,"AMP %g\n",nmodelScale);
      fclose(fout);
      
      /*plot5(preWhiteSpecX,preWhiteSpecY,nPreWhiteSpec,usePreWhitening,highFreqSpecX,highFreqSpecY,nHighFreqSpec,
	modelAlpha,modelFc,modelNfit,modelScale,nmodelScale,cholSpecX,cholSpecY,nCholSpec,cholWspecX,cholWspecY,nCholWspec,makeps,wn,0);
      */

      // Step 5c: recalculate the Cholesky matrix
      T2calculateCholesky(modelAlpha,modelFc,nmodelScale,fitVar,uinv,covFunc,resx,resy,rese,nres,highFreqRes,&errorScaleFactor,0,0,0);
    }
  else
  // Step 6: output the covariance function
  outputCovarianceFunction(covFunc,(int)(resx[nres-1]-resx[0])+2,errorScaleFactor,psr);
  //    outputMatrix(uinv,nres);
  
  
  T2getWhiteRes(resx,resy,rese,nres,uinv,cholWhiteY);
  if(writeFiles)fileOutput3("cholWhiteRes.dat",resx,cholWhiteY,rese,nres);


  // put in negative x // this should simulate the absence of a planet?

  long double x;
 

  int it;
 
  long int idum;
  idum=-10;

  double mass,Gn,Pb;
  Gn = 6.67e-8;

  Pb = psr[0].param[param_pb].val[0]*24*3600;

  x = 1e-4;
  int ispec;
  
  ispec=10;

 

  float x_old, mass_old, fdet_old;
  float x_int, mass_int, fdet_int;
 

  FILE *detfile;
  detfile=fopen("detect.dat", "w");
  fclose(detfile);
  detfile=fopen("det_spec.dat", "w");
  fclose(detfile);
  
  for(ispec=1;ispec<nCholSpec;ispec++)
    {
      
      

      float fdet;
      fdet=0.;

      x=1e-5;

      Pb =  1./origSpecX[ispec]*24*3600;
      //fprintf(stderr, "%.3le %.3le\n",origSpecX[ispec], Pb);
      //exit(0);
      int firstit=0;

      x=1e-5;
       mass = (double) pow(4*M_PI*M_PI/Gn/Pb/Pb, 0.33)*x*3e10*pow(2e33,0.67)/6e27;

      do
	{

	  x_old=x;
	  mass_old=mass;
	  fdet_old =fdet;


	  if (firstit!=0)
	    {

	      if (fdet < 0.1)
		{
		  
		  x *=2;
		}
	      else if (fdet < 0.5)
		{
		  x *=1.5;
		}
	      else if (fdet < 0.94)
		{
		  x *= 1.1;
		}
	    }
	  firstit=1;


	  mass = (double) pow(4*M_PI*M_PI/Gn/Pb/Pb, 0.33)*x*3e10*pow(2e33,0.67)/6e27;

	 
	 
	  
 
	  fdet=0.;
    
	  //  psr[0].param[param_om].val[0] = 0;
	  //  psr[0].param[param_ecc].val[0] =0;
	  //psr[0].param[param_t0].val[0]=50000;
	  //psr[0].param[param_pb].val[0] = 500;
	  
      
	  for(it=0;it<nit;it++)
	    {
	      
	  
	  
	      int jj,kk;
	      psr[0].nconstraints =0;
	      psr[0].nJumps = 0;
	      for(kk=0;kk<MAX_JUMPS;kk++)
		{
		  psr[0].jumpVal[kk] = 0.0;
		  psr[0].jumpValErr[kk] = 0.0;
		}
	      for(jj=0;jj<MAX_PARAMS;jj++){
		psr[0].param[jj].nLinkTo = 0;
		psr[0].param[jj].nLinkFrom = 0;
	      }
	    
	  
 

	
	      

	      readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
	     


	      readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
	      preProcess(psr,*npsr,argc,argv);
	      psr[0].param[param_pb].val[0] = (long double) Pb/24./3600.;
	      psr[0].param[param_a1].val[0] = x;
	      psr[0].param[param_t0].val[0] = 50000 + Pb/24./3600.*TKranDev(&idum);
	      psr[0].param[param_ecc].val[0] = 0.8*TKranDev(&idum);
	      psr[0].param[param_om].val[0]= 2*M_PI*TKranDev(&idum);
	      
	      
	      for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
		{
		  formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
		  formResiduals(psr,*npsr,1);    /* Form the residuals                 */
		  if (i==0) {
		    doFitAll(psr,*npsr,"J0401-7608.model");
		  }
		  //else textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
		}
	      
	      
	      nres = obtainTimingResiduals(psr,resx,resy,rese,ip);
	      
	      
      

	      
	      fprintf(stderr, "MASS: %.3le\n", mass);
	      
	      //if (inpt==-1) {printf("Continue (press '1') "); scanf("%s",dummy);}
	      
	      nCholSpec = calcSpectraErr(uinv,resx,resy,nres,cholSpecXp,cholSpecYp,cholSpecEp,-1);
	      
	  
	  
	      
	      T2calculateDailyCovariance(resx,cholWhiteY,rese,nres,whiteCovar,whiteCovarNpts,&zerolagWhiteCovar,0);
	      
	      double exsig;
	      exsig=  modelfcn(cholSpecX[ispec], nmodelScale, modelFc, modelAlpha, wn);  


	      
	      if(cholSpecYp[ispec] > cholSpecY[ispec])
		{
		  fdet += 1./((float) nit);
		}
	      

	      T2getWhiteRes(resx,resy,rese,nres,uinv,cholWhiteY);
	      
	  
	      plot4(resx,resy,rese,nres,cholWhiteY,whiteCovar,whiteCovarNpts,zerolagWhiteCovar);
	      
	      // Step 4c: get white residuals using the Cholesky matrix
	      
	      
	      nCholWspec = T2calculateSpectra(resx,cholWhiteY,rese,nres,0,0,2,
					      cholWspecX,cholWspecY);
	  
	      
	      
	      
	      //if (inpt==-1) {printf("Continue (press '1') "); scanf("%s",dummy);}
	      
	      plot5(cholSpecXp,cholSpecYp,nCholSpec,usePreWhitening,highFreqSpecX,highFreqSpecY,nHighFreqSpec,
		modelAlpha,modelFc,modelNfit,modelScale,nmodelScale,cholSpecX,cholSpecY,nCholSpec,cholWspecX,cholWspecY,nCholWspec,makeps,wn,Pb);
	      
	    }

	  detfile = fopen("det_spec.dat", "a");
      
	  fprintf(detfile,"%.3Le %.3le %.3le %.3e\n",  x , mass, cholSpecX[ispec], fdet);
	  fclose(detfile);


	  


	  
	  


	}while(fdet < 0.9 );
	  
      detfile = fopen("detect.dat", "a");
      // should probably interpoltae to get mass and x
     

      float fac;
      fac=(fdet-0.95)/(fdet-fdet_old);

      x_int = (1-fac)*x + fac*x_old;
      mass_int = (1-fac)*mass +fac*mass_old;



      fprintf(detfile,"%.3e %.3e %.3e %.3e\n",  cholSpecX[ispec],x_int,mass_int, x_old);
      fclose(detfile);
     
    }

 
  



  // Deallocate memory 
  free_uinv(uinv);
  //  free_uinv(uinvI);
  free(covFunc);
}


double modelfcn(double freq, double nmodelScale, double modelFc, double modelAlpha, double wn)
{
  double val;

  val = nmodelScale*(1.0/pow((1.0+pow(freq*365.25/modelFc,2)),modelAlpha/2.0)) + wn;

  return val;
}

void outputCovarianceFunction(double *covFunc,int n,double errorScaleFactor,pulsar* psr)
{
  FILE *fout;
  char fname[100];
  int i;

  fprintf(stderr, "here\n");
  exit(0);



  sprintf(fname,"covarFunc.dat_%s",psr[0].name);

  if (!(fout = fopen(fname,"w")))
    {
      printf("Unable to open output file: %s\n",fname);
      exit(1);
    }
  fprintf(fout,"%.15g\n",errorScaleFactor);
  for (i=0;i<n;i++)
    {
      fprintf(fout,"%.15g\n",covFunc[i]);
    }
  fclose(fout);
}

void removeMean(double *resx,double *resy,int n)
{
  double meanx,meany;
  double maxx,minx;
  int i;
  maxx = TKfindMax_d(resx,n);
  minx = TKfindMin_d(resx,n);
  meany = TKmean_d(resy,n);
  for (i=0;i<n;i++)
    {
      resx[i] -= (maxx+minx)*0.5;
      resy[i] -= meany;
    }
}

void fileOutput3(char *fname,double *x,double *y,double *z,int n)
{
  FILE *fout;
  int i;
  if (!(fout = fopen(fname,"w")))
    {
      printf("Unable to open file >%s< for output\n",fname);
      exit(1);
    }
  for (i=0;i<n;i++)
    fprintf(fout,"%g %g %g\n",x[i],y[i],z[i]);
  fclose(fout);
}

void fileOutput2(char *fname,double *x,double *y,int n)
{
  FILE *fout;
  int i;
  if (!(fout = fopen(fname,"w")))
    {
      printf("Unable to open file >%s< for output\n",fname);
      exit(1);
    }
  for (i=0;i<n;i++)
    fprintf(fout,"%g %g\n",x[i],y[i]);
  fclose(fout);
}

void calculateCholeskyCovarFunc(double bestAmp,double bestLag,int nGridFit,double **uinv,double *resx,
				double *resy,double *rese,int nres,double *covarFunc)
{
  int ndays = (int)(ceil(resx[nres-1])-floor(resx[0])+1)+2; // Add two extra days for interpolation  
  int i,j;
  double tt;
  int debug=1;

  for (i=0;i<ndays;i++)
    covarFunc[i] = bestAmp*exp(-(double)i/bestLag);
  tt = covarFunc[0];
  //  for (i=0;i<ndays;i++)
  //    covarFunc[i] = covarFunc[i]*pow(86400.0*365.25,2)/tt*1.0e-12;
  if (debug==1)
    {
      FILE *fout;
      fout = fopen("covarfunc","w");
      for (i=0;i<ndays;i++)
	fprintf(fout,"%g %g\n",(double)(i+1),covarFunc[i]);
      fclose(fout);
    }
  double** m=malloc_uinv(nres);
  cholesky_covarFunc2matrix(m,covarFunc,ndays,resx,resy,rese,nres,0);
  for(i=0;i<nres;i++){
	 m[i][i]+=rese[i]*rese[i];
  }

  cholesky_formUinv(uinv,m,nres);
  free_uinv(m);

}

// Do a grid search for minimum in chisq for fitting an exponential function
void fitExponential(double *resx,int nres,
		    double *rawCovar,int *rawCovarNpts,double *ampFit,double *chisqFit,double *bestAmp,
		    double *bestLag,double *bestChisq,int *nGridFit)
{
  int i,j;
  int nc=0;
  int dspan=(int)(resx[nres-1]-resx[0]+0.5);
  double lag;
  double tl,bl;
  double x[dspan],y[dspan],e[dspan];
  int foundMin=0;

  for (i=1;i<dspan;i++)
    {
      if (rawCovarNpts[i] > 0)
	{
	  x[nc]=i; // Fitting linearly with x
	  y[nc]=rawCovar[i];
	  e[nc]=1;
	  nc++;
	}
    }
  
  // Do the grid search
  for (i=0;i<nc;i++)
    {
      lag = (double)x[i];
      tl = 0;
      bl = 0;
      for (j=0;j<nc;j++)
	{
	  tl += y[j]*exp(-x[j]/lag)/e[j];
	  bl += exp(-2*x[j]/lag)/e[j];
	}
      ampFit[i] = tl/bl;
      // Calculate chisq
      chisqFit[i] = 0.0;
      for (j=0;j<nc;j++)
	chisqFit[i]+=pow((y[j]-ampFit[i]*exp(-x[j]/lag))/e[j],2);
      //      printf("ampFit = %g %g %g\n",ampFit[i],chisqFit[i],lag);
      //
      // Find first minimum
      //
      //            printf("Have %g %g\n",lag,chisqFit[i]);
      if (i==0)
	{
	  *bestLag=lag;
	  *bestAmp=ampFit[i];
	  *bestChisq = chisqFit[i];
	}
      else if (chisqFit[i] < *bestChisq && foundMin==0)
	{
	  *bestChisq = chisqFit[i];
	  *bestLag=lag;
	  *bestAmp=ampFit[i];	  
	}
      else if (chisqFit[i] > *bestChisq && foundMin==0)
	foundMin=1;
    }
  *nGridFit = nc;
  //  printf("Returning %g %g %g %g\n",*bestLag,*bestAmp,*bestChisq,(double)(*nGridFit));
}

void plot3a(double *resx,double *resy,int nres,double *rawCovar,int *rawCovarNpts,double zerolagRawCovar,
	    double *ampFit,double *fitChisq,int nGridFit,double bestAmp,double bestLag,double bestChisq,
	    int makeps)
{
  int dspan=(int)(resx[nres-1]-resx[0]+0.5);
  float fx[dspan],fy[dspan];
  float fx2[10],fy2[10];
  float fx3[nGridFit],fy3[nGridFit],fy4[nGridFit];
  float fy5[dspan];
  float fx6[dspan],fy6[dspan];
  float minx,maxx,miny,maxy;
  int i,j,np,ncovar;
  int nc=0,nc2=0;

  for (i=1;i<dspan;i++)
    {
      if (rawCovarNpts[i] > 0)
	{
	  fx[nc]=log10(i);
	  fy[nc]=rawCovar[i];
	  fy5[nc] = bestAmp*exp(-pow(10,fx[nc])/bestLag);
	  nc++;
	}
    }
  nc2=0;

  for (i=1;i<dspan;i++)
    {
      if (rawCovarNpts[i] > 0)
	{
	  fx3[nc2] = log10(i);
	  fy3[nc2] = ampFit[nc2];
	  fy4[nc2] = fitChisq[nc2];      
	  nc2++;
	}
    }
  minx = TKfindMin_f(fx,nc);
  maxx = TKfindMax_f(fx,nc);


  // WARNING REMOVE
  //  zerolagRawCovar = 8.5e-11;

  miny = -2*fabs(zerolagRawCovar);//TKfindMin_f(fy,nc);
  maxy = 2*fabs(zerolagRawCovar);//TKfindMax_f(fy,nc);
  printf("Here with %g %g\n",miny,maxy);
  if (makeps==1)
    {
      int addi;
      cpgbeg(0,"spectralPlot3a.ps/ps",1,1);
      cpgsch(2);
      cpgpap(0,0.400);
      cpgsch(2);
      cpgsfs(2);
      cpgslw(2);
      //      cpgsch(1.4);  cpgsfs(2);  cpgslw(2);
      printf("MIN/MAX = %g %g %g %g\n",minx,maxx,miny,maxy);
      cpgenv(minx,maxx,miny/1.3,maxy/1.3,0,10);
      cpglab("Lag (d)","Covariance","");
      cpgtext(0.1,-1e-10,"(b)");
      //      cpgpt(nc,fx,fy,1);

      ncovar=0;
      addi=1;
      for (i=1;i<dspan;i+=addi)
	{
	  np = 0;
	  fx6[ncovar] = 0.0;
	  fy6[ncovar] = 0.0;
	  printf("Have %d %d\n",i,addi);
	  for (j=i;j<i+addi;j++)
	    {
	      if (rawCovarNpts[j] > 0)
		{
		  fx6[ncovar] += j*rawCovarNpts[j];
		  fy6[ncovar] += (rawCovar[j]*rawCovarNpts[j]);
		  np+=rawCovarNpts[j];
		}
	    }
	  if (np>0)
	    {
	      fx6[ncovar]=log10(fx6[ncovar]/(double)np);
	      fy6[ncovar]/=(double)np;
	      ncovar++;
	    }
	  if (i==10 && addi==1) {i=0; addi=10;}
	  if (i==100 && addi==10) {i=0; addi=100;}
	  if (i==1000 && addi==100) {i=0; addi=1000;}
	}
      cpgsci(1); cpgpt(ncovar,fx6,fy6,5);   cpgsci(1);

            cpgsci(1); cpgsls(1); cpgline(nc,fx,fy5); cpgsci(1); cpgsls(1);
      fx2[0] = minx;
      fx2[1] = maxx;
      fy2[0] = 0;
      fy2[1] = 0;
      cpgsls(1); cpgline(2,fx2,fy2);
      fx2[0] = minx;
      fy2[0] = zerolagRawCovar;
       cpgsch(3); cpgsci(1); cpgpt(1,fx2,fy2,23); cpgsch(1.4); cpgsci(1);            
      cpgend();
    }

  cpgbeg(0,pgdevice,1,1);
  cpgsch(1.4);
  cpgsfs(2);
  cpgslw(2);
  //  cpgenv(minx,maxx,miny,maxy,0,10);
  cpgsvp(0.1,0.95,0.6,0.95);
  cpgswin(minx,maxx,miny,maxy);
  cpgsch(1);
  cpgbox("BCTLS",0,0,"BCNTS",0,0);
  cpgsci(14);  cpgbox("G",0,0,"G",0,0);  cpgsci(1);
  cpglab("","Covariance","");
  cpgsch(1.4);
  cpgpt(nc,fx,fy,20);
  cpgsci(3); cpgsls(3); cpgline(nc,fx,fy5); cpgsci(1); cpgsls(1);
  fx2[0] = minx;
  fy2[0] = zerolagRawCovar;
  cpgsch(2); cpgsci(2); cpgpt(1,fx2,fy2,15); cpgsch(1.4); cpgsci(1);

  // Plot chisq
  cpgsvp(0.1,0.95,0.35,0.6);
  miny = TKfindMin_f(fy4,nc2);
  maxy = TKfindMax_f(fy4,nc2);
  cpgswin(minx,maxx,miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny));
  cpgsch(1);
  cpgbox("BCTLS",0,0,"BCNTS",0,0);
  cpgsci(14);  cpgbox("G",0,0,"G",0,0);  cpgsci(1);
  cpglab("","Chisq","");
  cpgsch(1.4);
  cpgsci(2); cpgline(nc2,fx3,fy4); cpgsci(1);
  fx2[0] = log10(bestLag); fy2[0] = bestChisq;
  printf("Plotting at %g %g %g\n",bestLag,bestChisq,bestAmp);
  cpgsch(2); cpgsci(2); cpgpt(1,fx2,fy2,15); cpgsch(1.4); cpgsci(1);

  cpgsvp(0.1,0.95,0.1,0.35);
  miny = TKfindMin_f(fy3,nc2);
  maxy = TKfindMax_f(fy3,nc2);
  cpgswin(minx,maxx,miny,maxy);
  cpgsch(1);
  cpgbox("BCTNLS",0,0,"BCNTS",0,0);
  cpgsci(14);  cpgbox("G",0,0,"G",0,0);  cpgsci(1);
  cpglab("Lag (d)","A","");
  cpgsch(1.4);
  cpgsci(2); cpgline(nc2,fx3,fy3); cpgsci(1);
  fx2[0] = log10(bestLag); fy2[0] = bestAmp;
  cpgsch(2); cpgsci(2); cpgpt(1,fx2,fy2,15); cpgsch(1.4); cpgsci(1);

  cpgend();
}

void outputMatrix(double **uinv,int nres)
{
  int i,j;
  FILE *fout;

  fout = fopen("/DATA/BRAHE_1/hob044/idcm.dat","w");
  for (i=0;i<nres;i++)
    {
      for (j=0;j<nres;j++)
	fprintf(fout,"%.15g ",uinv[i][j]);
    }
  fclose(fout);


}

void plot6(double *cholSpecX,double *cholSpecY,int nCholSpec,double *cholWspecX,
	   double *cholWspecY,int nCholWspec,double *highFreqSpecX,
	   double *highFreqSpecY,int nHighFreqSpec,int makeps)
{
  int i,j;
  float fx[nCholSpec],fy[nCholSpec];
  float fx2[nCholSpec],fy2[nCholSpec];
  float fx3[nHighFreqSpec],fy3[nHighFreqSpec];
  float minx,maxx,miny,maxy;
  FILE *fout;

  fout = fopen("cholWhiteSpec.dat","w");
  printf("In plot 6\n");
  for (i=0;i<nCholSpec;i++)
    {
      fx[i] = log10(cholSpecX[i]*365.25);
      fy[i] = log10(cholSpecY[i]);
    }
  for (i=0;i<nCholWspec;i++)
    {
      fx2[i] = log10(cholWspecX[i]*365.25);
      fy2[i] = log10(cholWspecY[i]*pow(365.25*86400,2));
      fprintf(fout,"%g %g\n",cholWspecX[i],pow(10,fy2[i]));
    }
  fclose(fout);
  for (i=0;i<nHighFreqSpec;i++)
    {
      fx3[i] = log10(highFreqSpecX[i]*365.25);
      fy3[i] = log10(highFreqSpecY[i]);
    }
  minx = TKfindMin_f(fx,nCholSpec);
  maxx = TKfindMax_f(fx,nCholSpec);
  miny = TKfindMin_f(fy,nCholSpec);
  maxy = TKfindMax_f(fy,nCholSpec);

  cpgbeg(0,pgdevice,1,1);
  cpgsvp(0.1,0.95,0.45,0.95);
  cpgswin(minx-0.1*(maxx-minx),maxx+0.1*(maxx-minx),miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny));
  cpgbox("BCLTS",0,0,"BCNTS",0,0);
  cpglab("","Power Spectral Density (yr\\u3\\d)","");
  cpgsci(7); cpgline(nCholSpec,fx,fy); cpgpt(nCholSpec,fx,fy,20); cpgsci(1);
  cpgsci(1); cpgline(nHighFreqSpec,fx3,fy3); cpgpt(nHighFreqSpec,fx3,fy3,20); cpgsci(1);
  for (i=(int)(miny-0.1*(maxy-miny))-1;i<(int)maxy+1;i++)
    {
      fx[0] = minx-0.1*(maxx-minx); fx[1] = maxx+0.1*(maxx-minx);
      for (j=0;j<10;j++)
	{
	  fy[0] = fy[1] = log10(pow(10,i)*(j+1));
	  cpgsls(4); cpgsci(14); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
	}
    }
  for (i=(int)(minx-0.1*(maxx-minx))-1;i<(int)maxx+1;i++)
    {
      fy[0] = miny-0.1*(maxy-miny); fy[1] = maxy+0.1*(maxy-miny);
      for (j=0;j<10;j++)
	{
	  fx[0] = fx[1] = log10(pow(10,i)*(j+1));
	  cpgsls(4); cpgsci(14); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
	}
    }



  minx = TKfindMin_f(fx2,nCholWspec);
  maxx = TKfindMax_f(fx2,nCholWspec);
  miny = TKfindMin_f(fy2,nCholWspec);
  maxy = TKfindMax_f(fy2,nCholWspec);
  cpgsvp(0.1,0.95,0.10,0.45);
  cpgswin(minx-0.1*(maxx-minx),maxx+0.1*(maxx-minx),miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny));
  cpgbox("BCNLTS",0,0,"BCNTS",0,0);
  cpglab("Frequency (yr\\u-1\\d)","Power Spectral Density (yr)","");
  cpgsci(2); cpgline(nCholWspec,fx2,fy2); cpgpt(nCholWspec,fx2,fy2,20); cpgsci(1);
  for (i=(int)(miny-0.1*(maxy-miny))-1;i<(int)maxy+1;i++)
    {
      fx[0] = minx-0.1*(maxx-minx); fx[1] = maxx+0.1*(maxx-minx);
      for (j=0;j<10;j++)
	{
	  fy[0] = fy[1] = log10(pow(10,i)*(j+1));
	  cpgsls(4); cpgsci(14); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
	}
    }
  for (i=(int)(minx-0.1*(maxx-minx))-1;i<(int)maxx+1;i++)
    {
      fy[0] = miny-0.1*(maxy-miny); fy[1] = maxy+0.1*(maxy-miny);
      for (j=0;j<10;j++)
	{
	  fx[0] = fx[1] = log10(pow(10,i)*(j+1));
	  cpgsls(4); cpgsci(14); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
	}
    }
  fx[0] = minx-0.1*(maxx-minx); fx[1] = maxx+0.1*(maxx-minx);
  fy[0] = fy[1] = log10(1.0/pow(10,fx2[nCholWspec-1]));
  cpgsci(3); cpgline(2,fx,fy);
  fy[0] = fy[1] = log10(3.0/pow(10,fx2[nCholWspec-1]));
  cpgsci(3); cpgline(2,fx,fy);
  fy[0] = fy[1] = log10(0.05/pow(10,fx2[nCholWspec-1]));
  cpgsci(3); cpgline(2,fx,fy);
  cpgsci(1);
  cpgend();

  if (makeps==1)
    {
      cpgbeg(1,"spectralPlot6a.ps/ps",1,1);
      cpgsch(2.0);
      cpgpap(0,0.400);
      cpgsch(2); cpgsfs(2); cpgslw(2);
      // Plot the white data spectrum
      cpgenv(minx-0.1*(maxx-minx),maxx+0.1*(maxx-minx),miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny),0,10);
      cpglab("Frequency (yr\\u-1\\d)","log\\d10\\u[Spectral density (yr)]","");
      cpgtext(-1.1,-3,"(c)");
      cpgsci(7); cpgline(nCholWspec,fx2,fy2); 
      cpgpt(nCholWspec,fx2,fy2,20); cpgsci(1);
      
      fx[0] = minx-0.1*(maxx-minx); fx[1] = maxx+0.1*(maxx-minx);
      fy[0] = fy[1] = log10(1.0/pow(10,fx2[nCholWspec-1]));
      cpgsls(4); cpgline(2,fx,fy);
      fy[0] = fy[1] = log10(3.0/pow(10,fx2[nCholWspec-1]));
      cpgsls(4); cpgline(2,fx,fy);
      fy[0] = fy[1] = log10(0.05/pow(10,fx2[nCholWspec-1]));
      cpgsls(4); cpgline(2,fx,fy);
      cpgsls(1);
      fx[0] = fx[1] = 0;
      fy[0] = miny-0.1*(maxy-miny); fy[1] = maxy+0.1*(maxy-miny);
      cpgsls(4); cpgline(2,fx,fy); cpgsls(1);
      cpgend();
      
    }

}

void plot5(double *preWhiteSpecX,double *preWhiteSpecY,int nPreWhiteSpec,
	   int usePreWhitening,double *highFreqSpecX,double *highFreqSpecY,int nHighFreqSpec,
	   double modelAlpha,double modelFc,int modelNfit,double modelScale,
	   double nmodelScale,double *cholSpecX,double *cholSpecY,int nCholSpec,
	   double *cholWspecX,double *cholWspecY,int nCholWspec,int makeps, double wn, double pb)
{
  int i,j;
  float fx[MAX_OBSN],fy[MAX_OBSN];
  float fx1[MAX_OBSN],fy1[MAX_OBSN];
  float fx2[nCholSpec],fy2[nCholSpec];
  float fx4[nCholSpec], fy4[nCholSpec];
  float pbx[2], pby[2];
  float fx3[nCholWspec],fy3[nCholWspec];
  float minx,maxx,miny,maxy;
  FILE *fout1,*fout2,*fout;

  // First do exactly the same as plot3
  cpgbeg(0,pgdevice,1,1);
  cpgask(0);
  cpgsch(1);
  cpgsfs(2);
  cpgslw(2);

  plot3(preWhiteSpecX,preWhiteSpecY,nPreWhiteSpec,usePreWhitening,highFreqSpecX,highFreqSpecY,
	nHighFreqSpec, modelAlpha, modelFc, modelNfit, modelScale,0,&minx,&maxx, wn);
  // Overlay spectrum obtained using Cholesky
  
  // Overlay new model
  for (i=0;i<nCholSpec;i++)
    {
      fx1[i] = log10(cholSpecX[i]*365.25);
      fy1[i] = log10(cholSpecY[i]);
    }
  cpgsci(7); cpgline(nCholSpec,fx1,fy1); 
  cpgpt(nCholSpec,fx1,fy1,20); cpgsci(1);

  // Now overplot the new model
    {
      for (i=0;i<nCholSpec;i++)
	{
	  fx2[i] = log10(cholSpecX[i]*365.25);
	  //	  fy2[i] = nmodelScale-log10(pow((1.0+pow(cholSpecX[i]*365.25/modelFc,2)),modelAlpha/2.0));
	  fy2[i] = log10(nmodelScale*(1.0/pow((1.0+pow(cholSpecX[i]*365.25/modelFc,2)),modelAlpha/2.0)));
	  //	  fy2[i] = log10(pow(10,nmodelScale)/365.25)-log10(pow((1.0+pow(cholSpecX[i]*365.25/modelFc,2)),modelAlpha/2.0));
	  fx4[i] = fx2[i];
	  fy4[i] = log10( sqrt((float) nCholSpec)*modelfcn(cholSpecX[i],  nmodelScale, modelFc,  modelAlpha, wn));


			 //fy4[i] = log10(log((double) nCholSpec)*(pow(10., fy2[i]) + wn));
	  //fprintf(stderr, "%.3e %.3e %d\n", fx4[i], fy4[i], nCholSpec);

	}
      cpgsls(2); cpgsci(2); cpgline(nCholSpec,fx4,fy4); cpgsci(1); cpgsls(1);
      pbx[0] = log10(365.25*24*3600/pb);
      pbx[1] = log10(365.25*24*3600/pb);
      pby[0] = -30;
      pby[1] = -10;
      
      cpgline(2, pbx, pby); 
    }
    printf("New model spectrum = %g\n",nmodelScale);
    //
    fout = fopen("cholWhiteSpec.dat","w");
  for (i=0;i<nCholWspec;i++)
    {
      fx3[i] = log10(cholWspecX[i]*365.25);
      fy3[i] = log10(cholWspecY[i]*pow(365.25*86400,2));
      fprintf(fout,"%g %g\n",cholWspecX[i],pow(10,fy3[i]));
    }
  fclose(fout);
  //  minx = TKfindMin_f(fx3,nCholWspec);
  //  maxx = TKfindMax_f(fx3,nCholWspec);
  miny = TKfindMin_f(fy3,nCholWspec);
  maxy = TKfindMax_f(fy3,nCholWspec);

  cpgsvp(0.1,0.95,0.15,0.45);
  cpgswin(minx-0.1*(maxx-minx),maxx+0.1*(maxx-minx),miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny));
  cpgbox("BCNLTS",0,0,"BCNTS",0,0);

  cpglab("Frequency (yr\\u-1\\d)","Power Spectral Density (yr)","");
  for (i=(int)(miny-0.1*(maxy-miny))-1;i<(int)maxy+1;i++)
    {
      fx[0] = minx-0.1*(maxx-minx); fx[1] = maxx+0.1*(maxx-minx);
      for (j=0;j<10;j++)
	{
	  fy[0] = fy[1] = log10(pow(10,i)*(j+1));
	  cpgsls(4); cpgsci(14); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
	}
    }
  for (i=(int)(minx-0.1*(maxx-minx))-1;i<(int)maxx+1;i++)
    {
      fy[0] = miny-0.1*(maxy-miny); fy[1] = maxy+0.1*(maxy-miny);
      for (j=0;j<10;j++)
	{
	  fx[0] = fx[1] = log10(pow(10,i)*(j+1));
	  cpgsls(4); cpgsci(14); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
	}
    }

  cpgsci(7); cpgline(nCholWspec,fx3,fy3); 
  cpgpt(nCholWspec,fx3,fy3,20); cpgsci(1);

  fx[0] = minx; fx[1] = maxx;
  fy[0] = fy[1] = log10(1.0/pow(10,fx3[nCholWspec-1]));
  cpgsci(3); cpgline(2,fx,fy);
  fy[0] = fy[1] = log10(3.0/pow(10,fx3[nCholWspec-1]));
  cpgsci(3); cpgline(2,fx,fy);
  fy[0] = fy[1] = log10(0.05/pow(10,fx3[nCholWspec-1]));
  cpgsci(3); cpgline(2,fx,fy);
  cpgsci(1);
  cpgend();

  if (makeps==1) // Make plot for paper
    {
      float nx1[MAX_OBSN],ny1[MAX_OBSN];
      float nx2[MAX_OBSN],ny2[MAX_OBSN];

      fout1 = fopen("model1.dat","w");
      fout2 = fopen("model2.dat","w");
      cpgbeg(1,"spectralPlot5a.ps/ps",1,1);
      cpgsch(2.0);
      cpgpap(0,0.400);
      cpgsch(2); cpgsfs(2); cpgslw(2);
      // Plot the white data spectrum
      cpgenv(-1.2,1,miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny),0,10);
      //      cpgenv(minx-0.1*(maxx-minx),maxx+0.1*(maxx-minx),miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny),0,10);
      cpglab("Frequency (yr\\u-1\\d)","log\\d10\\u[Power spectral density (yr)]","");
      cpgtext(-1.1,-0.2,"(c)");
      cpgsci(7); cpgline(nCholWspec,fx3,fy3); 
      cpgpt(nCholWspec,fx3,fy3,20); cpgsci(1);
      
      fx[0] = minx-0.1*(maxx-minx); fx[1] = maxx+0.1*(maxx-minx);
      fy[0] = fy[1] = log10(1.0/pow(10,fx3[nCholWspec-1]));
      cpgsls(4); cpgline(2,fx,fy);
      fy[0] = fy[1] = log10(3.0/pow(10,fx3[nCholWspec-1]));
      cpgsls(4); cpgline(2,fx,fy);
      fy[0] = fy[1] = log10(0.05/pow(10,fx3[nCholWspec-1]));
      cpgsls(4); cpgline(2,fx,fy);
      cpgsls(1);
      fx[0] = fx[1] = 0;
      fy[0] = miny-0.1*(maxy-miny); fy[1] = maxy+0.1*(maxy-miny);
      cpgsls(4); cpgline(2,fx,fy); cpgsls(1);
      cpgend();
	
      // Plot the high freq res. and Cholesky spectrum
      cpgbeg(1,"spectralPlot5b.ps/ps",1,1);
      cpgsch(2.0);
      cpgpap(0,0.400);
      cpgsch(2); cpgsfs(2); cpgslw(2);
      miny = TKfindMin_f(fy2,nCholSpec);
      maxy = TKfindMax_f(fy2,nCholSpec);
      for (i=0;i<nHighFreqSpec;i++)
	{
	  nx1[i] = log10(highFreqSpecX[i]*365.25);
	  ny1[i] = log10(highFreqSpecY[i]);
	  if (ny1[i] > maxy) maxy = ny1[i];
	  if (ny1[i] < miny) miny = ny1[i];
	}
      // Also have fx1, fy1;
      //      cpgenv(minx,maxx,miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny),0,10);
      cpgenv(-1.2,1,miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny),0,10);
      cpglab("Frequency (yr\\u-1\\d)","log\\d10\\u[Power spectral density (yr\\u3\\d)]","");
      cpgtext(-1.1,-19,"(d)");
      cpgsls(4); cpgline(nHighFreqSpec,nx1,ny1); // cpgpt(nHighFreqSpec,nx1,ny1,16);
      cpgsls(1); cpgline(nCholSpec,fx1,fy1); cpgsls(1);
      cpgline(nCholSpec,fx2,fy2);
      // old model
      for (i=0;i<nCholSpec;i++)
	{
	  fx2[i] = log10(cholSpecX[i]*365.25);
	  //	  fy2[i] = nmodelScale-log10(pow((1.0+pow(cholSpecX[i]*365.25/modelFc,2)),modelAlpha/2.0));

	  fy2[i] = log10(nmodelScale*(1.0/pow((1.0+pow(cholSpecX[i]*365.25/modelFc,2)),modelAlpha/2.0)));
	  fprintf(fout2,"%g %g\n",cholSpecX[i],pow(10,fy2[i]));

	  fy2[i] = log10(modelScale*(1.0/pow((1.0+pow(cholSpecX[i]*365.25/modelFc,2)),modelAlpha/2.0)));
	  fprintf(fout1,"%g %g\n",cholSpecX[i],pow(10,fy2[i]));
	}
      cpgsls(2); cpgsci(2); cpgline(nCholSpec,fx2,fy2); cpgsci(1); cpgsls(1);
      fx[0] = fx[1] = 0;
      fy[0] = miny-0.1*(maxy-miny); fy[1] = maxy+0.1*(maxy-miny);
      cpgsls(4); cpgline(2,fx,fy); cpgsls(1);

      cpgend();
      fclose(fout1); fclose(fout2);

      //      plot3(preWhiteSpecX,preWhiteSpecY,nPreWhiteSpec,usePreWhitening,highFreqSpecX,highFreqSpecY,
      //	    nHighFreqSpec, modelAlpha, modelFc, modelNfit, modelScale,0);
  // Overlay spectrum obtained using Cholesky
  
  // Overlay new model
  //  for (i=0;i<nCholSpec;i++)
  //    {
  //      fx[i] = log10(cholSpecX[i]*365.25);
  //      fy[i] = log10(cholSpecY[i]);
  //    }


      
    }
  
}


void plot4(double *resx,double *resy,double *rese,int nres,double *cholWhiteY,
	   double *whiteCovar,int *whiteCovarNpts,double zerolagWhiteCovar)
{
  int dspan=(int)(resx[nres-1]-resx[0]+0.5);
  float fx1[nres],fy1[nres],fy1e1[nres],fy1e2[nres];
  float fx2[nres],fy2[nres];
  float fx3[dspan],fy3[dspan];
  float minx,maxx,miny,maxy;
  int nc=0;
  int i;

  for (i=0;i<nres;i++)
    {
      fx1[i] = fx2[i] = (float)resx[i];
      fy1[i] = (float)resy[i];
      fy2[i] = (float)cholWhiteY[i];
      fy1e1[i] = fy1[i] - (float)rese[i];
      fy1e2[i] = fy1[i] + (float)rese[i];
    }
  for (i=1;i<dspan;i++)
    {
      if (whiteCovarNpts[i] > 0)
	{
	  fx3[nc]=log10(i);
	  fy3[nc]=whiteCovar[i];
	  nc++;
	}
    }

  cpgbeg(0,"98/xs",1,1);
  cpgsch(1.4);
  cpgsfs(2);
  cpgslw(2);
  // Plot the original timing residuals
  minx = TKfindMin_f(fx1,nres);
  maxx = TKfindMax_f(fx1,nres);
  miny = TKfindMin_f(fy1,nres);
  maxy = TKfindMax_f(fy1,nres);
  
  cpgsvp(0.1,0.95,0.71,0.95);
  cpgswin(minx-0.1*(maxx-minx),maxx+0.1*(maxx-minx),miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny));
  cpgsch(1);
  cpgbox("BCTS",0,0,"BCNTS",0,0);
  cpgsci(14);  cpgbox("G",0,0,"G",0,0);  cpgsci(1);
  cpglab("","Residual (s)","");
  cpgsch(1.4);
  cpgpt(nres,fx1,fy1,20);
  cpgerry(nres,fx1,fy1e1,fy1e2,1);

  // Plot the Cholesky whitened timing residuals
  minx = TKfindMin_f(fx2,nres);
  maxx = TKfindMax_f(fx2,nres);
  miny = TKfindMin_f(fy2,nres);
  maxy = TKfindMax_f(fy2,nres);
  cpgsvp(0.1,0.95,0.47,0.71);
  cpgswin(minx-0.1*(maxx-minx),maxx+0.1*(maxx-minx),miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny));
  cpgsch(1);
  cpgbox("BCNTS",0,0,"BCNTS",0,0);
  cpgsci(14);  cpgbox("G",0,0,"G",0,0);  cpgsci(1);
  cpglab("Day","White Res. (norm.)","");
  cpgsch(1.4);
  cpgpt(nres,fx2,fy2,20);

  minx = TKfindMin_f(fx3,nc);
  maxx = TKfindMax_f(fx3,nc);
  miny = -2*fabs(zerolagWhiteCovar); //TKfindMin_f(fy3,nc);
  maxy = 2*fabs(zerolagWhiteCovar); //TKfindMax_f(fy3,nc);
  
  cpgsvp(0.1,0.95,0.10,0.35);
  cpgswin(minx,maxx,miny,maxy);
  cpgsch(1);
  cpgbox("BCNLTS",0,0,"BCNTS",0,0);
  cpgsci(14);  cpgbox("G",0,0,"G",0,0);  cpgsci(1);
  cpglab("Lag","Cov.","");
  cpgpt(nc,fx3,fy3,20);
  fx3[0] = minx;
  fy3[0] = zerolagWhiteCovar;
  cpgsch(2);
  cpgsci(2); cpgpt(1,fx3,fy3,18); cpgsci(1); cpgsch(1.4);
  //


  cpgend();
}







void plot3(double *preWhiteSpecX,double *preWhiteSpecY,int nPreWhiteSpec,
	   int usePreWhitening,double *highFreqSpecX,double *highFreqSpecY,
	   int nHighFreqSpec,double modelAlpha,double modelFc,int modelNfit,double modelScale,
	   int closeit,float *minx,float *maxx, float wn)
{
  int i,j;
  float fx1[nPreWhiteSpec],fy1[nPreWhiteSpec];
  float fx2[nHighFreqSpec],fy2[nHighFreqSpec];
  float fx3[nPreWhiteSpec],fy3[nPreWhiteSpec];
   float fx4[nPreWhiteSpec],fy4[nPreWhiteSpec];
  float fx[2],fy[2];
  float miny,maxy;

  *minx=*maxx=(float)log10(preWhiteSpecX[0]*365.25);
  miny=maxy=(float)log10(preWhiteSpecY[0]);

  for (i=0;i<nPreWhiteSpec;i++)
    {
      fx1[i] = log10(preWhiteSpecX[i]*365.25);
      fy1[i] = log10(preWhiteSpecY[i]);
      if (*minx > fx1[i]) *minx = fx1[i];
      if (*maxx < fx1[i]) *maxx = fx1[i];
      if (miny > fy1[i]) miny = fy1[i];
      if (maxy < fy1[i]) maxy = fy1[i];      
    }

  for (i=0;i<nHighFreqSpec;i++)
    {
      fx2[i] = log10(highFreqSpecX[i]*365.25);
      fy2[i] = log10(highFreqSpecY[i]);
      if (*minx > fx2[i]) *minx = fx2[i];
      if (*maxx < fx2[i]) *maxx = fx2[i];
      if (miny > fy2[i]) miny = fy2[i];
      if (maxy < fy2[i]) maxy = fy2[i];      
    }
  if (closeit==1)
    {
      cpgbeg(0,pgdevice,1,1);
      cpgask(0);
      cpgsch(1.4);
      cpgsfs(2);
      cpgslw(2);
      
      cpgenv(*minx-0.1*((*maxx)-(*minx)),(*maxx)+0.1*((*maxx)-(*minx)),
	     miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny),0,10);
      cpglab("Frequency (yr\\u-1\\d)","Power Spectral Density (yr\\u3\\d)","");

    }
  else
    {
      cpgsvp(0.1,0.95,0.45,0.95);
      cpgswin((*minx)-0.1*((*maxx)-(*minx)),(*maxx)+0.1*((*maxx)-(*minx)),miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny));
      cpgbox("BCLTS",0,0,"BCNTS",0,0);
      cpglab("","Power Spectral Density (yr\\u3\\d)","");

    }
  for (i=(int)(miny-0.1*(maxy-miny))-1;i<(int)maxy+1;i++)
    {
      fx[0] = (*minx)-0.1*((*maxx)-(*minx)); fx[1] = (*maxx)+0.1*((*maxx)-(*minx));
      for (j=0;j<10;j++)
	{
	  fy[0] = fy[1] = log10(pow(10,i)*(j+1));
	  cpgsls(4); cpgsci(14); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
	}
    }
  for (i=(int)((*minx)-0.1*((*maxx)-(*minx)))-1;i<(int)(*maxx)+1;i++)
    {
      fy[0] = miny-0.1*(maxy-miny); fy[1] = maxy+0.1*(maxy-miny);
      for (j=0;j<10;j++)
	{
	  fx[0] = fx[1] = log10(pow(10,i)*(j+1));
	  cpgsls(4); cpgsci(14); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
	}
    }
  cpgsci(14);  cpgbox("",0,0,"G",0,0);  cpgsci(1);

  cpgsci(2+usePreWhitening);
  cpgline(nPreWhiteSpec,fx1,fy1);
  cpgpt(nPreWhiteSpec,fx1,fy1,20);

  cpgsci(5);
  cpgline(nHighFreqSpec,fx2,fy2);
  cpgpt(nHighFreqSpec,fx2,fy2,20);
  cpgsci(1);

  // Draw model if required
  if (modelNfit!=-1)
    {
      float fx3[nPreWhiteSpec],fy3[nPreWhiteSpec];
      for (i=0;i<nPreWhiteSpec;i++)
	{
	  fx3[i] = log10(preWhiteSpecX[i]*365.25);
	  //	  fy3[i] = modelScale-log10(pow((1.0+pow(preWhiteSpecX[i]*365.25/modelFc,2)),modelAlpha/2.0));
	  fy3[i] = log10(modelScale*(1.0/pow((1.0+pow(preWhiteSpecX[i]*365.25/modelFc,2)),modelAlpha/2.0)));
	  //	  fy3[i] = log10(pow(10,modelScale)/365.25)-log10(pow((1.0+pow(preWhiteSpecX[i]*365.25/modelFc,2)),modelAlpha/2.0));
	  
	  fx4[i] = fx3[i];
	  fy4[i] = log10(powf(10., fy3[i])+ wn);
	  
	}
      cpgsls(3); cpgsci(1); cpgline(nPreWhiteSpec,fx4,fy4);  cpgsci(1); cpgsls(1);
    }
  if (closeit==1) cpgend();
}


void plot2(double *origSpecX,double *origSpecY,int nOrigSpec,double *smoothSpecX0,
	   double *smoothSpecY0,int nSmoothSpec0,double *smoothSpecX1,
	   double *smoothSpecY1,int nSmoothSpec1,double *smoothSpecX2,
	   double *smoothSpecY2,int nSmoothSpec2,double *highFreqSpecX,
	   double *highFreqSpecY,int nHighFreqSpec,int makeps)
{
  float fx[2],fy[2];
  float fx1[nOrigSpec],fy1[nOrigSpec];
  float fx2[nSmoothSpec0],fy2[nSmoothSpec0];
  float fx3[nSmoothSpec1],fy3[nSmoothSpec1];
  float fx4[nSmoothSpec2],fy4[nSmoothSpec2];
  float fx5[nHighFreqSpec],fy5[nHighFreqSpec];
  float minx,maxx,miny,maxy;
  int i;

  for (i=0;i<nOrigSpec;i++)
    {
      fx1[i] = log10(origSpecX[i]*365.25);
      fy1[i] = log10(origSpecY[i]);
      if (i==0)
	{
	  minx = maxx = fx1[i];
	  miny = maxy = fy1[i];
	}
      else
	{
	  if (minx > fx1[i]) minx = fx1[i];
	  if (maxx < fx1[i]) maxx = fx1[i];
	  if (miny > fy1[i]) miny = fy1[i];
	  if (maxy < fy1[i]) maxy = fy1[i];
	}

    }
  // Smooth spec
  for (i=0;i<nSmoothSpec0;i++)
    {
      fx2[i] = log10(smoothSpecX0[i]*365.25);
      fy2[i] = log10(smoothSpecY0[i]);
      if (minx > fx2[i]) minx = fx2[i];
      if (maxx < fx2[i]) maxx = fx2[i];
      if (miny > fy2[i]) miny = fy2[i];
      if (maxy < fy2[i]) maxy = fy2[i];

    }

  // Smooth spec prewhite 1
  for (i=0;i<nSmoothSpec1;i++)
    {
      fx3[i] = log10(smoothSpecX1[i]*365.25);
      fy3[i] = log10(smoothSpecY1[i]);
      if (minx > fx3[i]) minx = fx3[i];
      if (maxx < fx3[i]) maxx = fx3[i];
      if (miny > fy3[i]) miny = fy3[i];
      if (maxy < fy3[i]) maxy = fy3[i];
    }
  // Smooth spec prewhite 2
  for (i=0;i<nSmoothSpec2;i++)
    {
      fx4[i] = log10(smoothSpecX2[i]*365.25);
      fy4[i] = log10(smoothSpecY2[i]);
      if (minx > fx4[i]) minx = fx4[i];
      if (maxx < fx4[i]) maxx = fx4[i];
      if (miny > fy4[i]) miny = fy4[i];
      if (maxy < fy4[i]) maxy = fy4[i];
    }
  // High frequency residuals
  for (i=0;i<nHighFreqSpec;i++)
    {
      fx5[i] = log10(highFreqSpecX[i]*365.25);
      fy5[i] = log10(highFreqSpecY[i]);
      if (minx > fx5[i]) minx = fx5[i];
      if (maxx < fx5[i]) maxx = fx5[i];
      if (miny > fy5[i]) miny = fy5[i];
      if (maxy < fy5[i]) maxy = fy5[i];
    }

  if (makeps==1) // For paper
    {
      float sx[MAX_OBSN],sy[MAX_OBSN];
      cpgbeg(0,"spectralPlot2.ps/ps",1,1);
      cpgsch(2.0);
      cpgpap(0,0.400);
      cpgsch(2.0); cpgsfs(2); cpgslw(2);
      //      cpgenv(minx,1,miny+1,maxy,0,10);
      cpgenv(-1.2,1,miny+1,maxy,0,10);
      cpglab("Frequency (yr\\u-1\\d)","log\\d10\\u[Power Spectral Density (yr\\u3\\d)]","");
      cpgtext(-1.1,-19,"(b)");
      //      cpgline(nOrigSpec,fx1,fy1);
      //      cpgsls(2);  cpgline(nHighFreqSpec,fx2,fy2);  cpgpt(nHighFreqSpec,fx2,fy2,20);  cpgsls(1);
      //      cpgsls(2);  cpgline(nHighFreqSpec,fx2,fy2);  cpgsls(1);
      cpgsls(1);  cpgline(nHighFreqSpec,fx3,fy3);    cpgsls(1);
      //      cpgsci(7);  cpgline(nHighFreqSpec,fx5,fy5); cpgpt(nHighFreqSpec,fx5,fy5,16);cpgsci(1);
      cpgsls(4);  cpgline(nHighFreqSpec,fx5,fy5); cpgsci(1);
      sx[0] = sx[1] = 0;
      sy[0] = miny; sy[1] = maxy;
      cpgsls(4); cpgline(2,sx,sy); cpgsls(1);
      // Overlay model
      for (i=0;i<nHighFreqSpec;i++)
	{
	  sx[i] = fx2[i];
	  	  sy[i] = log10(3.27306e-18*(1.0/pow((1.0+pow(smoothSpecX1[i]*365.25/0.2,2)),5.0/2.0)));
		  // sy[i] = log10(7.88532e-23*(1.0/pow((1.0+pow(smoothSpecX1[i]*365.25/0.04,2)),5.0/2.0)));
	}
      cpgslw(3); cpgline(nHighFreqSpec,sx,sy); cpgslw(2);
      cpgend();
    }


  cpgbeg(0,pgdevice,1,1);
  cpgask(0);
  cpgsch(1.4);
  cpgsfs(2);
  cpgslw(2);
  cpgenv(minx,maxx,miny,maxy,0,10);
  cpgsci(14);  cpgbox("G",0,0,"G",0,0);  cpgsci(1);
  cpglab("Frequency (yr\\u-1\\d)","Power Spectral Density (yr\\u3\\d)","");
  cpgline(nOrigSpec,fx1,fy1);
  cpgpt(nOrigSpec,fx1,fy1,20);
  cpgsci(2);  cpgline(nSmoothSpec0,fx2,fy2);  cpgpt(nSmoothSpec0,fx2,fy2,20);  cpgsci(1);
  cpgsci(3);  cpgline(nSmoothSpec1,fx3,fy3);  cpgpt(nSmoothSpec1,fx3,fy3,20);  cpgsci(1);
  cpgsci(4);  cpgline(nSmoothSpec2,fx4,fy4);  cpgpt(nSmoothSpec2,fx4,fy4,20);  cpgsci(1);
  cpgsci(5);  cpgline(nHighFreqSpec,fx5,fy5); cpgpt(nHighFreqSpec,fx5,fy5,20);  cpgsci(1);

  // Overplot slopes of -2,-4,-6
  fx[0] = minx;
  fy[0] = maxy;
  fx[1] = maxx;
  fy[1] = maxy + (maxx-minx)*(-2);
  cpgsls(4); cpgline(2,fx,fy); cpgsls(1);
  fy[1] = maxy + (maxx-minx)*(-4);
  cpgsci(7); cpgsls(4); cpgline(2,fx,fy); cpgsls(1); cpgsci(1);
  fy[1] = maxy + (maxx-minx)*(-6);
  cpgsci(8); cpgsls(4); cpgline(2,fx,fy); cpgsls(1); cpgsci(1);


}

int calculateSpectra(double *x,double *y,double *e,int n,int useErr,int preWhite,
		     int specType,double *specX,double *specY)
{
  int nSpec;
  int i;
  int passSpecType;
  double outY_im[n],outY_re[n];

  if (specType==1) // Weighted Lomb Scargle
    passSpecType=4;
  else if (specType==2)
    passSpecType=2;

  TKspectrum(x,y,e,n,0,0,0,0,preWhite,passSpecType,1,1,1,specX,specY,&nSpec,0,0,outY_re,outY_im);
  
  return nSpec;
}
void getHighFreqCovar(double *resx,double *rese,double *highFreqRes,int nres,double *hfNormCovar,int *hfNormCovarNpts,double *hfZerolagNormCovar)
{
  int i,j;
  double yv[nres];
  double mean;

  // Should remove normal mean or weighted mean??
  mean = TKmean_d(highFreqRes,nres);
  for (i=0;i<nres;i++)
    yv[i] = (highFreqRes[i]-mean)/rese[i];

  calculateDailyCovariance(resx,yv,rese,nres,hfNormCovar,hfNormCovarNpts,hfZerolagNormCovar,0);

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

void getHighFreqRes(double *resy,double *smoothModel,int nres,double *highFreqRes)
{
  int i;
  for (i=0;i<nres;i++)
    highFreqRes[i] = resy[i] - smoothModel[i];
}

// Note: do a weighted exponential smoothing.  Bill originally used an unweighted smoother
void findSmoothCurve(double *resx,double *resy,double *rese,
		     int nres,double *cubicVal,double *smoothModel,double expSmooth)
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

void plot1(double *resx,double *resy,double *rese,int nres,double *cubicVal,double *smoothModel,double *highFreqRes,double *hfNormCovar,int *hfNormCovarNpts,double hfZerolagNormCovar)
{
  int dspan=(int)(resx[nres-1]-resx[0]+0.5);
  float fx[nres],fy[nres],fe1[nres],fe2[nres];
  float fx2[dspan],fy2[dspan];
  float fy3[nres],fy4[nres];
  float fx5[dspan],fy5[dspan];
  float minx,maxx,miny,maxy;
  int i,j,ncovar,np;

  cpgbeg(0,pgdevice,1,1); cpgsch(1.4); cpgsfs(2); cpgslw(2);
  // Plot timing residuals
  for (i=0;i<nres;i++)
    {
      fx[i] = (float)resx[i];
      fy[i] = (float)resy[i];
      fe1[i] = (float)(resy[i]-rese[i]);
      fe2[i] = (float)(resy[i]+rese[i]);
    }
  minx = TKfindMin_f(fx,nres);  maxx = TKfindMax_f(fx,nres);
  miny = TKfindMin_f(fy,nres);  maxy = TKfindMax_f(fy,nres);
  cpgsvp(0.1,0.95,0.71,0.95);
  cpgswin(minx-0.1*(maxx-minx),maxx+0.1*(maxx-minx),miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny));
  cpgsch(1);
  cpgbox("BCTS",0,0,"BCNTS",0,0); cpgsci(14);  cpgbox("G",0,0,"G",0,0);  cpgsci(1);
  cpglab("","Residual (s)","");
  cpgsch(1.4);
  cpgpt(nres,fx,fy,1);
  cpgerry(nres,fx,fe1,fe2,1);
  
  // Overplot the cubic fit
  for (i=0;i<dspan;i++)
    {
      fx2[i] = resx[0]+i;
      fy2[i] = cubicVal[0] + cubicVal[1]*fx2[i] + cubicVal[2]*pow(fx2[i],2) + 
	cubicVal[3]*pow(fx2[i],3);
    }
  cpgsci(2); cpgline(dspan,fx2,fy2); cpgsci(1);

  // Overplot the smoothed model of the residuals
  for (i=0;i<nres;i++) fy3[i] = (float)smoothModel[i];
  cpgsci(3); cpgsls(2); cpgline(nres,fx,fy3); cpgsls(1); cpgsci(1);

  // Plot the high frequency residuals
  for (i=0;i<nres;i++) fy4[i] = (float)highFreqRes[i];
  miny = TKfindMin_f(fy4,nres);
  maxy = TKfindMax_f(fy4,nres);

  cpgsvp(0.1,0.95,0.47,0.71);
  cpgswin(minx-0.1*(maxx-minx),maxx+0.1*(maxx-minx),miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny));
  cpgsch(1);
  cpgbox("BNCTS",0,0,"BCNTS",0,0);
  cpgsci(14);  cpgbox("G",0,0,"G",0,0);  cpgsci(1);
  cpglab("Day","HF Residual (s)","");
  cpgpt(nres,fx,fy4,20);
  
  // Plot the normalised covariances
  ncovar=0;
  for (i=1;i<dspan;i++)
    {
      if (hfNormCovarNpts[i] > 0)
	{
	  fx5[ncovar] = log10(i);
	  fy5[ncovar] = hfNormCovar[i];
	  ncovar++;
	}
    }
  minx = TKfindMin_f(fx5,ncovar);
  maxx = TKfindMax_f(fx5,ncovar);
  miny = -2.0*fabs(hfZerolagNormCovar); //TKfindMin_f(fy5,ncovar);
  maxy = 2.0*fabs(hfZerolagNormCovar); //TKfindMax_f(fy5,ncovar);
  cpgsvp(0.1,0.95,0.10,0.35);
  cpgswin(minx,maxx,miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny));
  
  cpglab("Lag (d)","Norm. Covariance","");
  //    cpgenv(minx,maxx,miny,maxy,0,10);
  cpgsci(14);  cpgbox("G",0,0,"G",0,0);  cpgsci(1);
  cpgbox("BNCTLS",0,0,"BCNTS",0,0);
  cpgpt(ncovar,fx5,fy5,20);
  fx5[0] = 0;fy5[0] = (float)hfZerolagNormCovar;
  cpgsch(2); cpgsci(2); cpgpt(1,fx5,fy5,18); cpgsci(1); cpgsch(1.4);

  // Now overplot the normalised covariance binned in 10 day intervals
  ncovar=0;
    for (i=1;i<dspan;i+=10)
    {
      np = 0;
      fx5[ncovar] = 0.0;
      fy5[ncovar] = 0.0;

      for (j=i;j<i+10;j++)
	{
	  if (hfNormCovarNpts[j] > 0)
	    {
	      fx5[ncovar] += j*hfNormCovarNpts[j];
	      fy5[ncovar] += (hfNormCovar[j]*hfNormCovarNpts[j]);
	      np+=hfNormCovarNpts[j];
	    }
	}
      if (np>0)
	{
	  fx5[ncovar]=log10(fx5[ncovar]/(double)np);
	  fy5[ncovar]/=(double)np;
	  ncovar++;
	}
    }
  cpgsci(3); cpgpt(ncovar,fx5,fy5,7);   cpgsci(1);

  // Now overplot the normalised covariance binned in 100 day intervals
  ncovar=0;
    for (i=1;i<dspan;i+=100)
    {
      np = 0;
      fx5[ncovar] = 0.0;
      fy5[ncovar] = 0.0;

      for (j=i;j<i+100;j++)
	{
	  if (hfNormCovarNpts[j] > 0)
	    {
	      fx5[ncovar] += j*hfNormCovarNpts[j];
	      fy5[ncovar] += (hfNormCovar[j]*hfNormCovarNpts[j]);
	      np+=hfNormCovarNpts[j];
	    }
	}
      if (np>0)
	{
	  fx5[ncovar]=log10(fx5[ncovar]/(double)np);
	  fy5[ncovar]/=(double)np;
	  ncovar++;
	}
    }
  cpgsci(4); cpgpt(ncovar,fx5,fy5,16);   cpgsci(1);
  
  cpgend();
}



// Fill resx with the SATs (in days), resy with post-fit residuals (s) and rese
// with TOA errors (s)
//
int obtainTimingResiduals(pulsar *psr,double *resx,double *resy,double *rese,int* ip)
{
  int i;
  int nres=0;
  sortToAs(psr);

  for (i=0;i<psr[0].nobs;i++)
    {
      if (psr[0].obsn[i].deleted==0)
	{
	  resx[nres] = (double)(psr[0].obsn[i].sat - psr[0].param[param_pepoch].val[0]);
	  // Check to make sure that the data are time sorted
	  if (nres>0 && resx[nres] < resx[nres-1])
	    {
		   logmsg("i=%d x[i]=%lg, x[i-1]=%lg",nres,resx[nres],resx[nres-1]);
	      logerr("ERROR: Data are not time sorted");
	      exit(1);
	    }
	  resy[nres] = (double)(psr[0].obsn[i].residual);
	  rese[nres] = (double)(psr[0].obsn[i].toaErr*1.0e-6);
	  ip[nres]=nres;
	  nres++;
	}
    }

  
  

  return nres;
}







char * plugVersionCheck = TEMPO2_h_VER;
