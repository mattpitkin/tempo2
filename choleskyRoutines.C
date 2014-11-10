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
  logmsg("Getting white residuals\n");
  logmsg("Updated version after transposing uinv\n");
  for (i=0;i<nres;i++)
    {
      sum=0.0;
      for (j=0;j<nres;j++)
	{
	  //	  sum+=uinv[j][i]*resy[j];
	  sum+=uinv[i][j]*resy[j];
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

int T2fitSpectra(double *preWhiteSpecX,double *preWhiteSpecY,int nPreWhiteSpec,double *modelAlpha,double *modelFc,int *modelNfit,double *modelScale,double *fitVar,int aval,int ipw,double ifc,double iexp,int inpt,double amp,int useBeta,double *betaVal,double cutoff)
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
      else if (ifc == -2) { // Automatically determine the corner frequency
	*modelFc = preWhiteSpecX[0]*365.25;
	printf("Automatically determining corner frequency to be %g\n",*modelFc);
      }
      else *modelFc = ifc;
      if (useBeta==1)
	{
	  printf("Enter beta "); scanf("%lf",betaVal);
	}

      if (iexp == 0) {printf("Enter power law exponential (should be positive) "); scanf("%lf",modelAlpha);}
      else if (iexp!=-2) *modelAlpha = iexp;	

      if (inpt == -1) {printf("Enter nfit "); scanf("%d",modelNfit);}
      else if (inpt == -2) // Automatically determine this
	{
	  int i;
	  *modelNfit = -1;

	  for (i=0;i<nPreWhiteSpec;i++)
	    {
	      //	      printf("testing %g %g\n",preWhiteSpecY[i],mean+rms);
	      if (preWhiteSpecY[i] < cutoff)
		{
		  //		  printf("STOPPING %d\n",i);
		  *modelNfit = i+1;
		  break;
		}
	    }
	  if (*modelNfit == -1)
	    {
	      printf("Cannot fit Nfit\n");
	      exit(1);
	    }
	  printf("Automatically found nfit = %d\n",*modelNfit);
	}
      else *modelNfit = inpt;

      if (iexp == -2)
	{
	  int i;
	  int n = *modelNfit;	 
	  double fx[n],fy[n],vals[2];

	  if (n==1) // Be careful
	    {
	      double ratio = preWhiteSpecY[0]/cutoff;
	      if (ratio < 2) // Probably all white
		{
		  *modelAlpha = 0;
		  *modelNfit = 6; // This should be set from the position of the filter and not hardcoded
		}
	      else // Probably the first point is really red noise and so alpha can be calculated. Probably should use the model
		*modelAlpha = (log10(preWhiteSpecY[0])-log10(preWhiteSpecY[1]))/(log10(preWhiteSpecX[0])-log10(preWhiteSpecX[1]));
	      printf("ratio = %g\n",ratio);
	    }
	  else
	    {
	      for (i=0;i<n;i++)
		{
		  fx[i] = log10(preWhiteSpecX[i]);
		  fy[i] = log10(preWhiteSpecY[i]);
		}
	      // Fit a straight line to these points
	      TKfindPoly_d(fx,fy,n,2,vals);
	      printf("Fit result = %g %g\n",vals[0],vals[1]);
	      *modelAlpha = -vals[1];
	    }
	  printf("Automatically found alpha = %g\n",*modelAlpha);
	}


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

	  if (useBeta==0)
	    m=1.0/pow(1.0+pow(fabs(preWhiteSpecX[i]*365.25)/(*modelFc),2),*modelAlpha/2.0);
	  else
	    m=pow(preWhiteSpecX[i]*365.25/(*modelFc),*betaVal)/pow(1.0+pow(fabs(preWhiteSpecX[i]*365.25)/(*modelFc),2),*modelAlpha/2.0);
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
int T2calculateCovarFunc(double modelAlpha,double modelFc,double modelA,int useBeta,double betaVal,double *covFunc,double *resx,double *resy,double *rese,int np){
   int ndays;
   int npts,i;
   double *p_r,*p_i;
   double freq;
   double P;
   double varScaleFactor=0.6;
   varScaleFactor=1.0; // test mjk 2013

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
   int noTF=0;
   FILE* tf = fopen("testf","w");
   if (!tf) {printf("Warning: unable to open output file: testf\n"); noTF=1;}

   if (useBeta==0)
     p_r[0]=modelA/pow(1.0+pow(fabs(0)/modelFc,2),modelAlpha/2.0);
   else
     p_r[0]=modelA*pow(fabs(0)/modelFc,betaVal)/pow(1.0+pow(fabs(0)/modelFc,2),modelAlpha/2.0);

   if (noTF==0) fprintf(tf,"%g %g\n",0,p_r[0]);
   for (i=1;i<=npts/2;i++){
	  freq=double(i)/(N*delta);
	  if (useBeta==0)
	    P=modelA/pow(1.0+pow(fabs(freq)/modelFc,2),modelAlpha/2.0);
	  else
	    P=modelA*pow(fabs(freq)/modelFc,betaVal)/pow(1.0+pow(fabs(freq)/modelFc,2),modelAlpha/2.0);
	  if (noTF==0) fprintf(tf,"%g %g\n",freq,P);
	  p_r[i]=P;
	  p_r[npts-i]=P;
	  p_i[i]=0;
	  p_i[npts-i]=0;
   }
   if (noTF == 0) fclose(tf);

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
double fitVar,double **uinv,double *covarFunc,double *resx,double *resy,
double *rese,int np,double *highFreqRes,double *errorScaleFactor, 
			 int dcmflag,int useBeta,double betaVal)
{
   int i,j,ndays;
   double** m=malloc_uinv(np);
   
   printf("choleskyRoutines: calculateCholesky\n");
   ndays=ceil((resx[np-1])-(resx[0])+1e-10);
   int ndays_out = T2calculateCovarFunc(modelAlpha,modelFc,modelA,useBeta,betaVal,covarFunc,resx,resy,rese,np);

   if(ndays!=ndays_out){
	  logerr("Ndays in != Ndays out!");
   }
   cholesky_covarFunc2matrix(m,covarFunc,ndays,resx,resy,rese,np,0);
   for(i=0;i<np;i++){
	  m[i][i]+=rese[i]*rese[i];
   }

   cholesky_formUinv(uinv,m,np);

   free_uinv(m);
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




void T2cholDecomposition(double **a, int n, double *p)
{
   int i,j,k;
   long double sum;
   // float sum;
   /*  for (i=0;i<n;i++)
	   {
	   for (j=0;j<n;j++)
	   printf("i=%d, j=%d, a=%g\n",i,j,a[i][j]);
	   }*/
   for (i=0;i<n;i++)
   {
	  for (j=i;j<n;j++)
	  {
		 for (sum=a[i][j],k=i-1;k>=0;k--) sum-=a[i][k]*a[j][k]; 
		 if (i==j)
		 {
			//	      printf("Currently have %d %d %Lg\n",i,j,sum);
			if (sum <= 0.0)
			{
			   printf("Here with %d %d %g %Lg\n",i,j,a[i][j],sum);
			   for (sum=a[i][j],k=i-1;k>=0;k--)
			   {
				  sum-=a[i][k]*a[j][k];
				  printf("Failed: %d %d %d %g %g %Lg\n",i,j,k,a[i][k],a[j][k],sum);
			   }
			   printf("Failed - the matrix is not positive definite\n");
			   exit(1);
			}
			p[i] = sqrt(sum);
			//	      printf("Currently have %d %d %Lg %g\n",i,j,sum,p[i]);
		 }
		 else
		 {
			a[j][i] = (double)(sum/p[i]);
			/*	      if (j==120)
					  printf("j=120, setting %g %Lg %g %d\n",a[j][i],sum,p[i],i);
					  if (j==130)
					  printf("j=130, setting %g %Lg %g %d\n",a[j][i],sum,p[i],i);*/
		 }
	  }
   }
}
