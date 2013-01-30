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

   
   cont = T2fitSpectra(preWhiteSpecX,preWhiteSpecY,nPreWhiteSpec,&modelAlpha,&modelFc,&modelNfit,&modelScale,&fitVar,1,usePreWhitening,modelFc, modelAlpha, modelNfit, -1);

   NFIT = modelNfit;

   fprintf(stderr, "npreWhiteSpec %d\n", nPreWhiteSpec);
   fprintf(stderr, "modelScale %.3e modelAlpha  %.3e modelFc  %.3e fitvar %.3e nfit %d\n", modelScale, modelAlpha, modelFc, fitVar, modelNfit);
   

   free(preWhiteSpecX);
   free(preWhiteSpecY);

  
   

   double errorScaleFactor=1;

   fprintf(stderr, "calculating Cholesky\n");
   T2calculateCholesky(modelAlpha,modelFc,modelScale,fitVar,uinv,covFunc,resx,resy,rese,nres,highFreqRes,&errorScaleFactor, 0);

   
   
   

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
      
       T2fitSpectra(cholSpecX,cholSpecY,nCholSpec,&modelAlpha, &modelFc,&modelNfit_chol,&nmodelScale,&fitVar,1,usePreWhitening,modelFc, modelAlpha, modelNfit,-1);
     }
   else
     {
       //fprintf(stderr, "made it here\n");
       // exit(0);
       T2fitSpectra(cholSpecX,cholSpecY,nCholSpec,&modelAlpha, &modelFc,&modelNfit_chol,&nmodelScale,&fitVar,1,usePreWhitening,modelFc, modelAlpha, modelNfit,-1);
     }

    T2fitSpectra(cholSpecX,cholSpecY,nCholSpec,&modelAlpha, &modelFc,&modelNfit_chol,&nmodelScale,&fitVar,1,usePreWhitening,modelFc, modelAlpha, modelNfit,-1);

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

    T2calculateCholesky(modelAlpha,modelFc,nmodelScale,fitVar,uinv,covFunc,resx,resy,rese,nres,highFreqRes,&errorScaleFactor,0);
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

   T2fitSpectra(specX,specY,nSpec,&modelAlpha,&modelFc,&modelNfit,&nmodelScale,&fitVar,1,usePreWhitening,ifc, iexp, inpt,-1);
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

// Used in the automatic determination of the covariance function
// R. Shannon
int T2guess_vals(double *x, double *y, int n, double *alpha, double *amp,  double *fc,  int  *nfit, double wn, double *fc_white, int prewhite)
{
 
  double *lx;
  double *ly;
  
  int i;
  int j;
  

  double dgw;
  dgw = x[1]-x[0];


  //exit(0);

  i=1;


 

  while((y[i] > wn) /* &&  (1./365.25-x[i]  > dgw )*/)
    {
      i++;
    }

  if((prewhite !=0) /*&& (1./365.25-x[i+1]  > dgw )*/)
    {
      i++;
    }
  //fprintf(stderr, "made it here %d\n", i);
  
  
  //i=i+4;
  //i = 5;
  
  
 

  j=i;

  if(i < 2)
    {
      //fprintf(stderr, "Error, data isn't white enough to calculate covariance function this way");
      //exit(0);
      *alpha = -0.01 ;
      *amp = 0;

      j=0;
      *fc = x[0];
      *nfit = n/2;

  
      *fc_white = x[0];
      //return 1;
    }
  else
    {
      
      
      *nfit=i;
      //*nfit =i-1;
      
      // for all the relevant values get the logs
      
 

      lx = (double*) malloc((*nfit)*sizeof(double));
      ly = (double*) malloc((*nfit)*sizeof(double));
      
      lx[0] = log10(x[0]);
      ly[0] = log10(y[0]);

      for(i=1;i<*nfit;i++)
	{
	  lx[i] = log10(x[i]);
	  ly[i] = log10(y[i]);//-log10(wn);
	  fprintf(stderr, "%.3e %.3e\n", lx[i], ly[i]);
	}
      
      

      //exit(0);
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

  //alpha =b is the slope of the line
  
    
      *alpha = clxy/vlx;
      

      *amp =  aly -(*alpha)*alx;
      
      // need to convert amplitude to physical units
     
      

      *amp = powf(10., *amp)*powf(x[0], *alpha);
      
      
      
      if(!(*alpha < 0))
	{
	  *alpha = -0.01; 
	  *amp =0.;
	  *fc_white = x[0];
	}
      if (*amp <0)
	{
	  *alpha = -0.01;
	  *amp =0.;
	  *fc_white = x[0];
	}
      else
	{
	  *fc_white = x[j];
	}
	  
	  
      *fc = x[0];

  
      

      

     
      //*fc_white = x[j];
  

      free(lx);
      free(ly);
    

      *fc =  x[0];
      /*
      if (*amp >0)
	{
	  // *fc_white=x[j];
	  *fc_white = pow(wn/(*amp), -1./(*alpha));
	}
      else
	{

	  *fc_white= x[0];
	  
	}
      */ 
     
}
  //fprintf(stderr, "%.3e %.3le %.3le %.3le %.3le\n", wn, *amp, wn/(*amp), *alpha, 1/(*fc_white));  
  //exit(0);
    
		      

  return 1;
}

// Used in the automatic determination of the covariance function
// R. Shannon
void T2getWhiteNoiseLevel(int n, double *y, int nlast, double *av)
{
  int i;
  
  *av =0;
  
  for(i=1;i<=nlast;i++)
    {
     
      *av += y[n-i];
    }

  *av /=(double) nlast;
}

void T2getWhiteRes(double *resx,double *resy,double *rese,int nres,double **uinv,double *cholWhiteY)
{
  int i,j;
  double sum;
  printf("choleskyRoutines: Getting white residuals\n");
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

// NOTE: Bill does an unweighted fit
//
void T2cubicFit(double *resx,double *resy,double *rese,int nres,double *cubicVal,double *cubicErr)
{
  double **cvm;
  int i;
  int nfit=4; // a+b*x+c*x^2 + d*x^3 
  double chisq;
  int useWeight=0;
  printf("choleskyRoutines: cubic\n");
  // Allocate memory for fit covariance matrix
  cvm = (double **)malloc(sizeof(double *)*nfit);
  for (i=0;i<nfit;i++) cvm[i] = (double *)malloc(sizeof(double)*nfit);

  TKleastSquares_svd(resx,resy,rese,nres,cubicVal,cubicErr,nfit,cvm,&chisq,TKfitPoly,useWeight);

  // Provide information about the fit
  printf("---------------------------------------------------------\n");
  printf("Cubic fit to timing residuals:\n");
  printf("y = %.5g + %.5g * x + %.5g * x^2 + %.5g * x^3\n",cubicVal[0],cubicVal[1],
	 cubicVal[2],cubicVal[3]);
  printf("chisq of fit = %g, reduced chisq = %g\n",chisq,chisq/(nres-4.0));
  printf("\n");
  printf("This curve is plotted as a red line through the timing residuals\n");
  printf("---------------------------------------------------------\n");
  

  // Free memory allocation
  for (i=0;i<nfit;i++) free(cvm[i]);
  free(cvm);
  
}

// Note: do a weighted exponential smoothing.  Bill originally used an unweighted smoother
void T2findSmoothCurve(double *resx,double *resy,double *rese,
		     int nres,double *cubicVal,double *smoothModel,double expSmooth)
{
  int i,j;
  double tl,bl;
  double diffy[nres];
  double dt;
  printf("choleskyRoutines: findSmoothCurve\n");
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

// Exponential smoothing and interpolation
void T2interpolate(double *resx,double *resy,double *rese,
		 int nres,double *cubicVal,double *interpX,
		 double *interpY,int *nInterp,int interpTime,double expSmooth)
{
  int i,j;
  double tl,bl,dt;
  double yval;
  printf("choleskyRoutines: interpolate\n");  
  *nInterp = ceil((resx[nres-1]-resx[0])/interpTime); // Number of interpolated points
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
  /*  if(writeFiles)
  {
    FILE *fout;
    fout = fopen("interp.dat","w");
    for (i=0;i<*nInterp;i++)
      fprintf(fout,"%g %g\n",interpX[i],interpY[i]);
    fclose(fout);
    }*/
}

void T2getHighFreqRes(double *resy,double *smoothModel,int nres,double *highFreqRes)
{
  int i;
  for (i=0;i<nres;i++)
    highFreqRes[i] = resy[i] - smoothModel[i];
}

int T2calculateSpectra(double *x,double *y,double *e,int n,int useErr,int preWhite,
		     int specType,double *specX,double *specY)
{
  int nSpec;
  int i;
  int passSpecType;
  double outY_im[n],outY_re[n];
  printf("choleskyRoutines: calculateSpectra\n");

  if (specType==1) // Weighted Lomb Scargle
    passSpecType=4;
  else if (specType==2)
    passSpecType=2;

  TKspectrum(x,y,e,n,0,0,0,0,preWhite,passSpecType,1,1,1,specX,specY,&nSpec,0,0,outY_re,outY_im);
  
  return nSpec;
}

int T2fitSpectra(double *preWhiteSpecX,double *preWhiteSpecY,int nPreWhiteSpec,double *modelAlpha,double *modelFc,int *modelNfit,double *modelScale,double *fitVar,int aval,int ipw,double ifc,double iexp,int inpt,double amp)
{
  static int time=1;
  double v1,v2,m;
  double df;
  int i;
  printf("choleskyRoutines: fitSpectra\n");
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



/*
 * New algorithm for computing the 
 *
 *
 */
int T2calculateCovarFunc(double modelAlpha,double modelFc,double modelA,double *covFunc,double *resx,double *resy,double *rese,int np){
   int ndays;
   int npts,i;
   double *p_r,*p_i;
   double freq;
   double P;
   double varScaleFactor=0.6;

   ndays=ceil((resx[np-1])-(resx[0])+1e-10);
   npts=128;
   while(npts<(ndays+1)*2)npts*=2;

   p_r=(double*)malloc(sizeof(double)*npts);
   p_i=(double*)malloc(sizeof(double)*npts);

   logmsg("Generating powerlaw model covariance function. nday=%d npts=%d",ndays,npts);


   for (i=1;i<npts;i++){
	  p_r[i]=0;
	  p_i[i]=0;
   }

   double delta=1.0/365.25;
   double N=(double)npts;

   p_r[0]=modelA/pow(1.0+pow(fabs(0)/modelFc,2),modelAlpha/2.0);
   for (i=1;i<=npts/2;i++){
	  freq=double(i)/(N*delta);
	  P=modelA/pow(1.0+pow(fabs(freq)/modelFc,2),modelAlpha/2.0);
	  p_r[i]=P;
	  p_r[npts-i]=P;
	  p_i[i]=0;
	  p_i[npts-i]=0;
   }

   TK_fft(1,npts,p_r,p_i);

   for (i=0; i <= ndays; i++){
	covFunc[i]=p_r[i]*pow(86400.0*365.25,2)*365.25*varScaleFactor;
   }

   free(p_r);
   free(p_i);
   return ndays;

}
   /*
	* I am pretty sure this routine is broken - i.e. computes the covariance function incorrectly
	* because it uses the 1-sided PSD rather than a symetric function.
	*
	* I have also re-written it to use TKspectrum rather than fftw so that it can be moved to be
	* part of libtempo2.
	*
	* M.Keith 2012. Code checked by W.Coles.
	*
	*
int T2calculateCovarFunc(double modelAlpha,double modelFc,double modelA,double *covFunc,double *resx,double *resy,double *rese,int np,double fitVar)
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

  printf("choleskyRoutines: calcCovarFunc\n");
  ndays = ceil((resx[np-1])-(resx[0])); 
  f  = (double *)malloc(sizeof(double)*(ndays*2+2)); // Frequency vector
  p  = (double *)malloc(sizeof(double)*(ndays*2+2)); // Model of pulsar power spectrum
  pf = (double *)malloc(sizeof(double)*(ndays*2+2)); // Periodic spectrum model
  opf = (double *)malloc(sizeof(double)*(ndays*2+2)); // Periodic spectrum model
  pe = (double *)malloc(sizeof(double)*(ndays*2+2)); // Periodic spectrum model

  printf("Number of days = %d\n",ndays);


  for (i=0;i<=2*ndays;i++)
    {
      f[i] = (double)(i)/(double)(2*ndays)*365.25;
            pf[i] = modelA/pow(1.0+pow(f[i]/modelFc,modelAlpha/2.0),2);
    }

  ndays *= 2;  
  
  printf("Obtaining covariance function from analytic model\n");
  {
    fftw_complex* output;
    fftw_plan transform_plan;
    double tt;

    output = (fftw_complex*)opf;
    transform_plan = fftw_plan_dft_r2c_1d(ndays, pf, output, FFTW_ESTIMATE);
    //   transform_plan = fftw_plan_dft_r2c_1d(j, pf, output, FFTW_ESTIMATE);
    fftw_execute(transform_plan);    
    fftw_destroy_plan(transform_plan);  
    for (i=0;i<=ceil((resx[np-1])-(resx[0]));i++) 
      {
	// Xinping does not have the varScale factor
	covFunc[i] = opf[2*i]/ndays*pow(86400.0*365.25,2)*365.25*varScaleFactor;
	printf("AAAA %lg %lg\n",opf[2*i],opf[2*i+1]);
      }

  printf("WARNING: varScaleFactor = %g (used to deal with quadratic removal)\n",varScaleFactor);
  }

  free(p);              
  free(f);
  free(pf);
  free(opf);
  free(pe);

  return ceil((resx[np-1])-(resx[0])); 
}*/

void T2calculateCholesky(double modelAlpha,double modelFc,double modelA,
double fitVar,double **uinv,double *covFunc,double *resx,double *resy,
double *rese,int np,double *highFreqRes,double *errorScaleFactor, 
int dcmflag)
{
  int i,j,nc;
  printf("choleskyRoutines: calculateCholesky\n");
  nc = T2calculateCovarFunc(modelAlpha,modelFc,modelA,covFunc,resx,resy,rese,np);
  T2formCholeskyMatrix_pl(covFunc,nc,resx,resy,rese,np,uinv);
  printf("Complete calculateCholesky\n");
}

void T2calculateDailyCovariance(double *x,double *y,double *e,int n,double *cv,int *in,double *zl,int usewt)
{
  int i,j;
  int nd = (int)(x[n-1]-x[0]+0.5);
  int nc=0;
  int nzl=0;
  int npts[nd];
  double dt;
  double wt[nd];
  double nzerolag=0.0;
  printf("choleskyRoutines: calculateDailyCovar\n");
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

void T2formCholeskyMatrix_pl(double *c,int cSize,double *resx,double *resy,double *rese,int np,double **uinv)
{
  double **m,**u,sum;
  double *cholp;
  int i,j,k,ix,iy;
  double t0,cint,t;
  int t1,t2;
  int debug=1;

  printf("choleskyRoutines: formCholeskyMatrix_pl\n");
  printf("Getting the covariance matrix, np = %d\n",np);
  m = (double **)malloc(sizeof(double *)*(np));
  u= (double **)malloc(sizeof(double *)*(np));
  cholp  = (double *)malloc(sizeof(double)*(np));  // Was ndays

  for (i=0;i<np;i++)
    {
      m[i] = (double *)malloc(sizeof(double)*(np));
      u[i] = (double *)malloc(sizeof(double)*(np));
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
	  if (t1 > cSize || t2 > cSize)
	    {
	      printf("ERROR: Something wrong in the Cholesky - requesting covariance function value greater than array size\n");
	      exit(1);
	    }
	  cint = c[t1]*(1-t)+c[t2]*t; // Linear interpolation
	  m[ix][iy] = cint;
	  //	  printf("searching for t1 = %d, t2 = %d, cint = %g\n",t1,t2,cint);
	}
    }
  printf("Multiplying by errors\n");
  for (ix=0;ix<np;ix++)
    m[ix][ix]+=(rese[ix]*rese[ix]);

  if (debug==1)
    {
      printf("m = \n\n");
      //      for (i=125;i<131;i++)
      //	{ 
      //	  for (j=125;j<131;j++) printf("%10g ",m[i][j]); 
      //	  printf("\n");
      //	}
    }
  printf("Doing the Cholesky Decomposition\n");
  // Do the Cholesky
  TKcholDecomposition(m,np,cholp);
  printf("Done the Cholesky decomposition\n");
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

  for (i=0;i<np;i++)
    {
      free(m[i]);
      free(u[i]);
    }
  free(m);
  free(u);
  free(cholp);
}

// Fill resx with the SATs (in days), resy with post-fit residuals (s) and rese
// with TOA errors (s)
//
int T2obtainTimingResiduals(pulsar *psr,double *resx,double *resy,double *rese)
{
  int i;
  int nres=0;
  printf("choleskyRoutines: obtainTimingResiduaks\n");
  for (i=0;i<psr[0].nobs;i++)
    {
      if (psr[0].obsn[i].deleted==0)
	{
	  resx[nres] = (double)(psr[0].obsn[i].sat - psr[0].param[param_pepoch].val[0]);
	  // Check to make sure that the data are time sorted
	  /*if (nres>0 && resx[nres] < resx[nres-1])
	    {
	      printf("ERROR: Data are not time sorted\n");
	      exit(1);
	    }
	  */	 

	  resy[nres] = (double)(psr[0].obsn[i].residual);
	  rese[nres] = (double)(psr[0].obsn[i].toaErr*1.0e-6);
	  //fprintf(stderr, "%.13le %.13le %.13le\n", resx[nres],(double) (psr[1].obsn[i].residual) , rese[nres]);
	  nres++;
	}
    }

  
  

  return nres;
}

void T2writeCovarFuncModel(double alpha,double fc,double val,double white,char *fname)
{
  FILE *fout;
  
  fout = fopen(fname,"w");
  fprintf(fout,"MODEL 1\n");
  fprintf(fout,"ALPHA %g\n",alpha);
  fprintf(fout,"FC %g\n",fc);
  fprintf(fout,"AMP %g\n",val);
  fprintf(fout,"WHITENOISE %g\n",white);
  fclose(fout);
}
