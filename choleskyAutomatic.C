/* Routines based around the TEMPO2 use of the Cholesky algorithm    */
/* These have been put together by G. Hobbs, W. Coles and R. Shannon */

#include "choleskyRoutines.h"
#include "tempo2.h"
#include "math.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "TKfit.h"
#include "TKspectrum.h"
#include "T2toolkit.h"
#include "cholesky.h"

// Automatic determination of the covariance function
// R. Shannon

void T2get_covFunc_automatic(pulsar *psr, double expSmooth, char *outname, double *fc_w,double *fc_r, double *modelAlpha_out, double *modelVal, double *whiteNoiseLevel, int realflag, int dcmflag)
{

  int nres;


  double *resx, *resy, *rese;
  double *specX, *specY, *specY_err;

  resx = (double*) malloc(MAX_OBSN* sizeof(float));
  resy = (double*) malloc(MAX_OBSN* sizeof(float));
  rese = (double*) malloc(MAX_OBSN* sizeof(float));
  specX = (double*) malloc(MAX_OBSN* sizeof(float));
  specY = (double*) malloc(MAX_OBSN* sizeof(float));
  specY_err = (double*) malloc(MAX_OBSN* sizeof(float));
  
  //double resx[MAX_OBSN], resy[MAX_OBSN], rese[MAX_OBSN];
  //double specX[MAX_OBSN], specY[MAX_OBSN], specY_err[MAX_OBSN];
  double modelAmp;

  // allocate the memory for the covariance fucntino


  nres = T2obtainTimingResiduals(psr,resx,resy,rese);
    
  int i,j;
  double **uinv;
  double *covFunc;
  int dspan=ceil(resx[nres-1]-resx[0]);
  
  uinv=malloc_uinv(nres);
  covFunc = (double *)malloc(sizeof(double)*((int)(resx[nres-1]-resx[0])+5));
  
  int nHighFreqSpec;
  double *highFreqSpecX, *highFreqSpecY;
  highFreqSpecX = (double*) malloc(MAX_OBSN*sizeof(double));
  highFreqSpecY = (double*) malloc(MAX_OBSN*sizeof(double));
  
  
  
  // start with residuals, calculate Cholesky matrix uinv and covariance function
 
  // fit cubic to residuals (no interaction required)
  
  double cubicVal[4], cubicErr[4];

  fprintf(stderr, "starting cubic fit\n");
  T2cubicFit(resx,resy,rese,nres,cubicVal,cubicErr);
  
  // obitain a smooth curve that models the residuals well
  // NEED to choose smoothing time expSmooth! (will this be the same in the simulated and real data?)
 
  double *smoothModel;
  smoothModel = (double*) malloc(nres*sizeof(double));


  double *highFreqRes;
  highFreqRes= (double*) malloc(nres*sizeof(double));
  //  double interpTime = 28;
  double interpTime = 28;
  int nInterp;

  double *interpX, *interpY;
  //double interpX[MAX_OBSN], interpY[MAX_OBSN];
  interpX = (double*) malloc(MAX_OBSN*sizeof(double));
  interpY = (double*) malloc(MAX_OBSN*sizeof(double));


  // obtain spectra with the correct level of prewhitening
   int usePreWhitening, nPreWhiteSpec;
   
   double* preWhiteSpecX, *preWhiteSpecY;
   double dprewhiteSpecX, dcholSpecX;
   preWhiteSpecX = (double*) malloc(MAX_OBSN*sizeof(double));
   preWhiteSpecY = (double*) malloc(MAX_OBSN*sizeof(double));
   //double preWhiteSpecX[MAX_OBSN], preWhiteSpecY[MAX_OBSN];

   int nSpec;

   FILE *specfile;
   int nLast;// =15;
   //double whiteNoiseLevel;

      double modelFc;
   double modelAlpha;
   int modelNfit;
   double fc_white;
   
   double  modelScale,fitVar;

   int ismooth, nsmooth;
   nsmooth=2;


   int prewhiteFlag;

   FILE *covfile;
   char *covname;
   covname =(char*)calloc(80, sizeof(char));
 
  // start loop to determine optimal level of smoothing and white noise level

   expSmooth=40;
   //exit(0);

   if (realflag > 0)
     {
       sprintf(covname,"%d.cov",realflag);
       //sprintf(covname,"%s.cov", psr->name);
       //fprintf(stderr, "%s", covname);
       covfile=fopen(covname, "r");
       fscanf(covfile, "%lf %d %lf %lf %d", &expSmooth, &usePreWhitening, &modelFc, &modelAlpha, &modelNfit);
       //fprintf(stderr, "%f %d %g %f %d\n", expSmooth, usePreWhitening, modelFc, modelAlpha, modelNfit);
       fclose(covfile);
       //exit(0);
     }
   else
     {
       expSmooth =40;
       usePreWhitening=0;
     }
   free(covname);

   expSmooth =40;
   for (ismooth=0;ismooth<nsmooth;ismooth++)
     {


       fprintf(stderr, "smoothing curve\n");
       T2findSmoothCurve(resx,resy,rese,nres,cubicVal,smoothModel,expSmooth);

       
       // obtain the high frequency residuals (no interaction required)
       
       fprintf(stderr, "high frequency residuals\n");
       T2getHighFreqRes(resy,smoothModel,nres,highFreqRes);

       // cput errors into univ matrix
       
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
  
 
       //double highFreqSpecX[MAX_OBSN], highFreqSpecY[MAX_OBSN];
       fprintf(stderr, "calculating spectra\n");
       nHighFreqSpec = calcSpectra(uinv,resx,highFreqRes,nres,highFreqSpecX,highFreqSpecY,-1);
       
       
       // interpolate the smooth curve
       
       fprintf(stderr, "interpolating\n");
       T2interpolate(resx,resy,rese,nres,cubicVal,interpX,interpY,
		   &nInterp,interpTime,expSmooth);
       

       
       fprintf(stderr, "calculating prewhitened spectra\n");
       usePreWhitening=0;
       prewhiteFlag=0;


     


    
       
       while(prewhiteFlag !=1 )
	 {
	   
	   
	   nPreWhiteSpec = T2calculateSpectra(interpX,interpY,rese,nInterp,0,usePreWhitening,2,
				       preWhiteSpecX,preWhiteSpecY);

	   printf("Here nPreWhiteSpec = %d\n",nPreWhiteSpec);
       
	   dprewhiteSpecX = preWhiteSpecX[1]-preWhiteSpecX[0];

	   // Estimate the level of white noise using highest frequency points in power spectrum
	   // This appears to be working
	   if(ismooth ==0)
	     {
	       nLast = nHighFreqSpec/2;
	     }
	   else
	     {
	         nLast = nHighFreqSpec/2;
		 //nLast = nHighFreqSpec*highFreqSpecX[0]/fc_white/2;///highFreqSpecX
	     }

	   fprintf(stderr, "getting white noise level\n");
	   //getWhiteNoiseLevel(nPreWhiteSpec,preWhiteSpecY, nLast, whiteNoiseLevel);
	   T2getWhiteNoiseLevel(nHighFreqSpec, highFreqSpecY, nLast, whiteNoiseLevel);

	   specfile=fopen("comp.dat", "w");
	   for(i=0;i<nPreWhiteSpec;i++)
	     {
	       fprintf(specfile, "%.3le %.3le %.3le %.3le\n",preWhiteSpecX[i], preWhiteSpecY[i] ,highFreqSpecX[i],highFreqSpecY[i]);
	     }
	   fclose(specfile);
	   


	   // get the correct spectra:  need to figure out way to automatically fit for the spectrum (e.g., get reasonable values for modelAlpha, modelFc, modelNfit;
	   
	   // guess correct values for model fit
	   //NEED to estimate high frequency noise level:
	   
	   

	   if (realflag ==0)
	     {
	       T2guess_vals(preWhiteSpecX, preWhiteSpecY,nPreWhiteSpec, &modelAlpha, &modelAmp, &modelFc, &modelNfit, *whiteNoiseLevel, &fc_white,  usePreWhitening);
	     }
	   
	   if(modelAlpha  < -0.1)
	     {
	       //modelAlpha *=1.2;  
	     }

	   //modelAlpha = -5;
	   //modelNfit =5;

	   //exit(0);

	   fprintf(stderr, "Using %d order prewhitening alpha = %.2f fc_white = %.3e white noise level %.3e\n", usePreWhitening, modelAlpha, fc_white, *whiteNoiseLevel);
	   //exit(0);

	   //exit(0);
	   
	   // use some logic to determine what order prewhitening is necessary

	   if( (usePreWhitening ==0) && (modelAlpha > -1))
	     {
	       fprintf(stderr, "here \n");
	       
	       prewhiteFlag=1;
	     }
	   else if( (usePreWhitening ==1) && (modelAlpha > -3.))
	     {
	       fprintf(stderr, "made it here\n");
	       //exit(0);
	       prewhiteFlag=1;
	     }
	   else if ((usePreWhitening ==2))
	     {
	       prewhiteFlag=1;
	       fprintf(stderr, "made it here\n");
	       //exit(0);
	     }
	   usePreWhitening++;
	   
	 }
   
       

       expSmooth=0.2/fc_white;
      EXPSMOOTH = expSmooth;
     
fprintf(stderr,"fc_white %.3e expsmooth %.3e\n", fc_white, expSmooth);
       //exit(0);
       
       usePreWhitening--;
     
     }
  
   UPW = usePreWhitening;

   free(interpX);
   free(interpY);

   
   free(highFreqSpecX);
   free(highFreqSpecY);

   modelAlpha = -1.*modelAlpha;
   modelFc =  modelFc*365.25;

  
   *fc_w = fc_white;

   // RMS:  these don't do anything because acor=1, but set them to some arbitrary values to be pedantic
   double ifc=1, iexp=1;
   int inpt=1;
   int cont;

   
   cont = T2fitSpectra(preWhiteSpecX,preWhiteSpecY,nPreWhiteSpec,&modelAlpha,&modelFc,&modelNfit,&modelScale,&fitVar,1,usePreWhitening,modelFc, modelAlpha, modelNfit, -1,0,NULL,-1);

   NFIT = modelNfit;

   fprintf(stderr, "npreWhiteSpec %d\n", nPreWhiteSpec);
   fprintf(stderr, "modelScale %.3e modelAlpha  %.3e modelFc  %.3e fitvar %.3e nfit %d\n", modelScale, modelAlpha, modelFc, fitVar, modelNfit);
   

   free(preWhiteSpecX);
   free(preWhiteSpecY);

  
   

   double errorScaleFactor=1;

   fprintf(stderr, "calculating Cholesky\n");
   T2calculateCholesky(modelAlpha,modelFc,modelScale,fitVar,uinv,covFunc,resx,resy,rese,nres,highFreqRes,&errorScaleFactor, 0,0,0);

   
   
   

   // get the white residuals using the Cholesky matrix

   double *cholWhiteY;
   cholWhiteY = (double*) malloc(MAX_OBSN*sizeof(double)); 
     //double cholWhiteY[MAX_OBSN];
   

   fprintf(stderr, "get white residuals 2\n");
   T2getWhiteRes(resx,resy,rese,nres,uinv,cholWhiteY);

   // get spectrum of whitened data
   int nCholWspec;
   double *cholWspecX, *cholWspecY;
   //  double cholWspecX[MAX_OBSN],cholWspecY[MAX_OBSN]; 
   
   cholWspecX = (double*)malloc(MAX_OBSN*sizeof(double));
   cholWspecY = (double*)malloc(MAX_OBSN*sizeof(double));


   fprintf(stderr, "calc spectral 2\n");
   nCholWspec = T2calculateSpectra(resx,cholWhiteY,rese,nres,0,0,2,
				  cholWspecX,cholWspecY);

   printf("Have nCholWspec = %d\n",nCholWspec);
   
   specfile=fopen("compchol.dat", "w");
   for(i=0;i<nCholWspec;i++)
     {
       fprintf(specfile, "%.3le %.3le \n",cholWspecX[i], cholWspecY[i]);
     }
   fclose(specfile);
  

   

   free(cholWspecX);
   free(cholWspecY);

   // get covariance of white residuals

   int *whiteCovarNpts;//whiteCovarNpts[MAX_OBSN];
   double *whiteCovar;//whiteCovar[MAX_OBSN];
   double zerolagWhiteCovar;//zerolagWhiteCovar;
   
   whiteCovarNpts = (int*) malloc(MAX_OBSN*sizeof(int));
   whiteCovar = (double*) malloc(MAX_OBSN*sizeof(double));




   fprintf(stderr, "recal covariance function\n");
   T2calculateDailyCovariance(resx,cholWhiteY,rese,nres,whiteCovar,whiteCovarNpts,&zerolagWhiteCovar,0);
 
   free(cholWhiteY);
   free(whiteCovarNpts);
   free(whiteCovar);
   

   fprintf(stderr, "finished recalculating  covariance\n");

   // improve the spectral estimate
   int nCholSpec;
  double *cholSpecX, *cholSpecY;
  cholSpecX =  (double*) malloc(MAX_OBSN*sizeof(double));
  cholSpecY = (double*) malloc(MAX_OBSN*sizeof(double));
   
   //double cholSpecX[MAX_OBSN],cholSpecY[MAX_OBSN];

   fprintf(stderr, "calc spectra 3\n");

   FCFINAL=modelFc;
   FCALPHA=modelAlpha;
   WNLEVEL=*whiteNoiseLevel;

   nCholSpec = calcSpectra(uinv,resx,resy,nres,cholSpecX,cholSpecY,-1);


   
   specfile=fopen("cholSpec.dat", "w");
   for(i=0;i<nCholSpec;i++)
     {
       fprintf(specfile, "%.3le %.3le \n",cholSpecX[i], cholSpecY[i]);
     }
   fclose(specfile);
   
   double nmodelScale;


   fprintf(stderr, "fit spectra \n");
   //guess_vals(cholSpecX, cholSpecY,nCholSpec, &modelAlpha, &modelAmp,  &modelFc, &modelNfit, *whiteNoiseLevel, &fc_white);
   //guess_vals(preWhiteSpecX, preWhiteSpecY,nPreWhiteSpec, &modelAlpha, &modelAmp, &modelFc, &modelNfit, *whiteNoiseLevel, &fc_white);
   
   nLast  = nCholSpec/3;
   // guess_vals(cholSpecX, cholSpecY,nCholSpec, &modelAlpha, &modelAmp,  &modelFc, &modelNfit, *whiteNoiseLevel, &fc_white);
   
   
   T2getWhiteNoiseLevel(nCholSpec, cholSpecY, nLast, whiteNoiseLevel);
   
   //guess_vals(cholSpecX, cholSpecY,nCholSpec, &modelAlpha, &modelAmp,  &modelFc, &modelNfit, *whiteNoiseLevel, &fc_white);

   // need to get new modelNFit
   // not the same for nonCholesky spectra
   int modelNfit_chol;
   
   dcholSpecX =  cholSpecX[1]- cholSpecX[0];

   modelNfit_chol = modelNfit*dprewhiteSpecX/dcholSpecX;

   if (modelAmp != 1000)
     {
      
       T2fitSpectra(cholSpecX,cholSpecY,nCholSpec,&modelAlpha, &modelFc,&modelNfit_chol,&nmodelScale,&fitVar,1,usePreWhitening,modelFc, modelAlpha, modelNfit, -1,0,NULL,-1);
     }
   else
     {
       //fprintf(stderr, "made it here\n");
       // exit(0);
       T2fitSpectra(cholSpecX,cholSpecY,nCholSpec,&modelAlpha, &modelFc,&modelNfit_chol,&nmodelScale,&fitVar,1,usePreWhitening,modelFc, modelAlpha, modelNfit, -1,0,NULL,-1);
     }

   T2fitSpectra(cholSpecX,cholSpecY,nCholSpec,&modelAlpha, &modelFc,&modelNfit_chol,&nmodelScale,&fitVar,1,usePreWhitening,modelFc, modelAlpha, modelNfit, -1,0,NULL,-1);

   free(cholSpecX);
   free(cholSpecY);

   //guess_vals(cholSpecX, cholSpecY,nCholSpec, &modelAlpha, &modelFc, &modelNfit, whiteNoiseLevel, &fc_white);
   *modelVal = nmodelScale; 
   *modelAlpha_out = modelAlpha;

   // recalculate Cholesky matrix

   //aaaa

    fprintf(stderr, "re calc cholesky \n");
    // does this need to be one beforehand?
    errorScaleFactor=1;

    T2calculateCholesky(modelAlpha,modelFc,nmodelScale,fitVar,uinv,covFunc,resx,resy,rese,nres,highFreqRes,&errorScaleFactor,0,0,0);
   nSpec =  calcSpectra(uinv,resx,resy, nres,specX,specY, -1);
   //  nSpec = calcSpectra_with_err(uinv,resx,resy, nres,specX,specY,specY_err,-1);  
   
   
   
   
 
   // nLast = 10;
   //fprintf(stderr, "get whitenoise 3 \n");
   nLast = nSpec/3;
   T2getWhiteNoiseLevel(nSpec, specY, nLast, whiteNoiseLevel);
   //fprintf(stderr, "guess vals 2\n");
   T2guess_vals(specX, specY,nSpec, &modelAlpha, &modelAmp, &modelFc, &modelNfit, *whiteNoiseLevel, &fc_white, usePreWhitening);
   modelAlpha = -1.*modelAlpha;
   modelFc =  modelFc*365.25;
   *fc_r = modelFc;

   //fitSpectra(cholSpecX,cholSpecY,nCholSpec,&modelAlpha,&modelFc,&modelNfit,&nmodelScale,&fitVar,1,usePreWhitening,ifc, iexp, inpt);

   T2fitSpectra(specX,specY,nSpec,&modelAlpha,&modelFc,&modelNfit,&nmodelScale,&fitVar,1,usePreWhitening,ifc, iexp, inpt, -1,0,NULL,-1);
   *modelVal = modelAmp; 
   *modelAlpha_out = modelAlpha;


   //fprintf(stderr, "%.3e %.3e %.3e\n", modelAmp, whiteNoiseLevel, *modelVal);
   //exit(0);
   fprintf(stderr, "%.3e %.3e\n", modelAlpha, modelFc);
   
  
   //fprintf(stderr, "%.3e\n", errorScaleFactor);
   //exit(0);
   
   // write out covariance function for diagnostic purposes
   //exit(0);


   FILE *covarfile;
   
   covarfile=fopen(outname, "w");

  
   fprintf(covarfile, "%.15g\n", errorScaleFactor);
   //   for (i=0;i<=(dspan+3);i++)
   for (i=0;i<(dspan+1);i++)
     {
       fprintf(covarfile,"%.15g\n",covFunc[i]);
     }
     
  
   
   fclose(covarfile);

   //fprintf(stderr, "%.10le\n", modelAlpha);
   
   //exit(0);

   // free memory (don't need to keep covFunc)

   free(smoothModel);
   free(highFreqRes);

    free_uinv(uinv);
    
    free(covFunc);

    

  free(resx);
  free(resy);
  free(rese);
  free(specX);
  free(specY);
  free(specY_err);


 return;
}
