//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russel Edwards

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

/* Routines useful for spectral analysis 
 * Based on work carried out by G. Hobbs
 * and W. Coles
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tempo2.h"
#include "T2toolkit.h"
#include "TKfit.h"
#include "TKspectrum.h"

double GLOBAL_OMEGA = 0;
/* Global variables used in the fortran code for sigmaz*/
int npt,nusewt,nxunits,ntunits,nformat,nwriteres,nbintype;
int npt1last,npt2last,ncubic,ncubics,ntau,linfile,indx[90000],ndim;
double data[90000],utjd[90000],taumin,sigmai[90000],permax,root2;
double utjd1,utjd2,tmin,tmax,xmin,xmax,utjdlast,tausec,taumax,tauday;
double prtl[5],utmean,secyear,taulog,addvar,tauyear,tauensure,tdiffmin;
double utfirst,utlast;


double globalOmega;
bool verbose_calc_spectra=false;

/* ************************************************************** *
 * General routines:
 *
 * TKspectrum() - spectral analysis of a data set
 * x = array of arrival times
 * y = array of residuals 
 * n = number of observations
 * averageTime > 0 => average to sampling of given number of days
 *             = 0 => no averaging
 * smoothWidth = smoothing width (days; odd, integer) 
 * smoothType   = 0 => No smoothing
 *              = 1 => Hann smoothing
 * preWhite = 0 => no pre-whitening
 *          = 1 => first differencing
 *          = 2 => second differencing
 * specType = 1 => DFT
 *          = 2 => Lomb-Scargle periodogram
 *          = 3 => FFT 
 *          = 4 => Weighted L-S periodogram
 * ofac = oversampling factor
 * hifac = the highest frequency to which we want to calculate, expressed as a multiple of the Nyquist frequency
 * output   = 0 => no text output
 *          = 1 => print sorted input data set, then stop
 *          = 2 => print averaged data points, then stop
 * 
 * Output:
 *
 * outX = frequency channels
 * outY = power 
 * nout = number of frequency channels
 *
 * NOTE: All arrays start from zero
 * ************************************************************** */

double TKspectrum(double *x,double *y,double *e,int n,int averageTime,int smoothWidth,int smoothType,
		int fitSpline,int preWhite,int specType,double ofac,double hifac,int specOut,double *outX,
		double *outY,int *nout,int calcWhite,int output, double *outY_re, double *outY_im)
{
  double avPtsX[n],avPtsY[n];
  double interpX[n],interpY[n];
  double smoothX[n],smoothY[n];  
  double specX[(int)ceil(1 + ofac * hifac * MAX_OBSN / 2.0)],specY[(int)ceil(1 + ofac * hifac * MAX_OBSN / 2.0)];
  int splineDaily=0;
  double mean;
  int    i,j;
  int    nAv=0;
  int    nInterp=0;
  int    nSmooth=0;
  int    nSpec=0;

  // Step 1: Sort the data in increasing time order
  //NB THIS SORTING CAN CHANGE THE ORDER OF SATS AND RESIDUALS!!!! Causes many strange problems
  TKsortit(x,y,n);
  
  if (output==1)
    {
      for (i=0;i<n;i++)
	printf("sortdata: %d %g %g\n",i,x[i],y[i]);
      exit(1);
    }
  // Step 2: Average the points within a given time of each other
  // CURRENTLY NOT USING ERROR BARS
  if (averageTime==0) // No averaging
    {
      for (i=0;i<n;i++)
	{
	  avPtsX[i] = x[i];
	  avPtsY[i] = y[i];
	}
      nAv = n;
    }
  else
    TKaveragePts(x,y,n,averageTime,avPtsX,avPtsY,&nAv);
  //  for(i=0;i<n;i++){
      //printf("Joris: %.10lg %.10lg %.10lg %.10lg.\n",x[i],avPtsX[i],y[i],avPtsY[i]);
  //}
  //  printf("Joris2: %d %d.\n",n,nAv);
  if (output==2)
    {
      for (i=0;i<nAv;i++)
	printf("averagePoints: %d %g %g\n",i,avPtsX[i],avPtsY[i]);
      exit(1);
    }

  // Step 3: Remove mean from the timing residuals
  mean=TKmean_d(avPtsY,nAv);
  for (i=0;i<nAv;i++)
    avPtsY[i]-=mean;
  if (output==3)
    {
      for (i=0;i<nAv;i++)
	printf("meanRemoved: %d %g %g\n",i,avPtsX[i],avPtsY[i]);
      exit(1);
    }

  // Step 4: Fit a spline to the data and obtain equally spaced sampling
  if (fitSpline==1)
    {
      double yd[MAX_OBSN][4],h;
      
      TKcmonot(nAv,avPtsX,avPtsY,yd);
      if (splineDaily==1)
	{
	  nInterp=0;
	  do{
	    interpX[nInterp] = avPtsX[0]+nInterp;
	    nInterp++;
	  }while (interpX[nInterp-1]<avPtsX[nAv-1]);
	  nInterp--;
	}
      else
	{
	  //	nInterp = nAv;
	  nInterp = 4096;
	  for (i=0;i<nInterp;i++)
	    interpX[i] = avPtsX[0]+i*(avPtsX[nAv-1]-avPtsX[0])/(double)nInterp;
	  //	  interpX[i] = avPtsX[i];
	  
	}
      TKspline_interpolate(nAv,avPtsX,avPtsY,yd,interpX,interpY,nInterp);
      if (output==4)
	{
	  for (i=0;i<nInterp;i++)
	    printf("spline: %d %g %g\n",i,interpX[i],interpY[i]);
	  exit(1);
	}
    }
  else
    {
      for (i=0;i<nAv;i++)
	{
	  interpX[i] = avPtsX[i];
	  interpY[i] = avPtsY[i];
	}
      nInterp = nAv;
    }
  // Step 5: Smooth using Hann filter
  if (smoothType==0) // No smoothing
    {
      for (i=0;i<nInterp;i++)
	{
	  smoothX[i] = interpX[i];
	  smoothY[i] = interpY[i];
	}
      nSmooth = nInterp;
    }
  else if (smoothType==1) // Hann smoother
    {
      int smoothWidth2;
      // Convert to days
      //      smoothWidth2 = (int)(smoothWidth*nInterp/(interpX[nInterp-1]-interpX[0]));
      smoothWidth2=smoothWidth;
      if (smoothWidth < 1) smoothWidth = 1;
      printf("SmoothWidth = %d\n",smoothWidth2);
      TKhann(interpX,interpY,nInterp,smoothX,smoothY,&nSmooth,smoothWidth2);
    }
  else
    {
      printf("SmoothType = %d unknown\n",smoothType);
      exit(1);
    }
  if (output==5)
    {
      for (i=0;i<nSmooth;i++)
	printf("smooth: %d %g %g\n",i,smoothX[i],smoothY[i]);
      exit(1);
    }

  // Do we want to calculate the white noise component?
  if (calcWhite==1)
    {
      double reInterpX[MAX_OBSN],reInterpY[MAX_OBSN];
      double yd[MAX_OBSN][4],h;
      int nReInterp=n,c=0;
      double xv=0.0,xvsq=0.0;
      
      TKcmonot(nSmooth,smoothX,smoothY,yd);
      for (i=0;i<n;i++)
	reInterpX[i] = x[i];
      TKspline_interpolate(nSmooth,smoothX,smoothY,yd,reInterpX,reInterpY,nReInterp);
      for (i=0;i<nReInterp;i++)
	{
	  //	  printf("reinterp: %d %g %g\n",i,reInterpX[i],reInterpY[i]);
	  y[i]-=reInterpY[i];
	  //	  printf("white: %d %g %g\n",i,x[i],y[i]);
	  if (x[i]>x[0]+smoothWidth/2.0 && x[i]<x[n-1]-smoothWidth/2.0)
	    {
	      xv+=y[i];
	      xvsq+=(y[i]*y[i]);
	      c++;
	    }
	}
      printf("rms = %g\n",sqrt(xvsq/(double)c-pow(xv/(double)c,2)));
      return sqrt(xvsq/(double)c-pow(xv/(double)c,2));
    }



  // Pre-whiten if requested
  if (preWhite==1)
    {
      TKfirstDifference(smoothX,smoothY,nSmooth);
      smoothY[0] = 0.0;      
      //      for (i=0;i<nSmooth;i++)
      //	printf("firstdiff %g %g\n",smoothX[i],smoothY[i]);
    }
  else if (preWhite==2)
    {
      TKfirstDifference(smoothX,smoothY,nSmooth);
      TKfirstDifference(smoothX,smoothY,nSmooth);
      smoothY[0] = 0.0;      
      smoothY[1] = 0.0;
    }
  // Step 6: Calculate the spectrum
  {
    int jmax;
    double prob,var,tspan;
    tspan = TKrange_d(smoothX,nSmooth);
    
    if (specType==1) // DFT
      {
	TK_dft(smoothX,smoothY,nSmooth,specX,specY,&nSpec,outY_re, outY_im);
	// The following code to form the PSD in units of yr^3 and 
	// frequency in 1/days (required for the whitening/postdarkening)
	// has been checked by G. Hobbs (30/10/08).
	//
	// Do not modify without leaving a comment here.
	//
	for (i=0;i<nSpec;i++)
	  {
	    if (specOut==1)
	      specY[i] = (specY[i]/pow(365.25*86400.0,2))*2*(tspan/365.25)/(double)nSmooth/(double)nSmooth;
	    else if (specOut==2) // Amplitude
	      specY[i]=sqrt(specY[i])/nSpec;
	    else if (specOut==3) // Power
	      specY[i]=specY[i]/nSpec/nSpec;


	  }
	// End of checked section
      }
    else if (specType==2) // Lomb-Scargle
      {
 	// The following code to form the PSD in units of yr^3 and 
	// frequency in 1/days (required for the whitening/postdarkening)
	// has been checked by G. Hobbs (30/10/08).
	//
	// Do not modify without leaving a comment here.
	//
	// 05 Jan 09: G Hobbs: changed to use TKlomb_d instead of numerical recipes
	// 05 Jan 09: G Hobbs: added PSD, amp and pow outputs (originally just PSD)
	// 06 May 09: D Yardley: added "*ofac/hifac" in PSD calculation. Now should be consistent measured power regardless of oversampling.
	// 06 May 09: D Yardley: added "*ofac*hifac" in power calculation. Now should be consistent measured power regardless of oversampling.
	// 11 Aug 09: D Yardley: added specOut = 4 option for normalized power (important when adding power spectra together)

	TKlomb_d(smoothX,smoothY,nSmooth,ofac,hifac,specX,specY,&nSpec,&var);

	for (i=0;i<nSpec;i++)
	  {
	    if (specOut==1) // PSD
	      {
		specY[i]*=(2.0*var*nSpec*ofac/hifac);
		specY[i] = (specY[i]/pow(365.25*86400.0,2))*2*(tspan/365.25)/(double)nSmooth/(double)nSmooth;
	      }
	    else if (specOut==2) // Amplitude
	      specY[i]=sqrt(specY[i]*2.0*var/nSpec);
	    else if (specOut==3) // Power
	      specY[i]=specY[i]*2.0*var/nSpec*ofac*hifac;
	    else if (specOut==4) // Power Normalised by variance.
	      specY[i]=specY[i]*2.0/nSpec*ofac*hifac;
	  }
	// End of checked section
      }
    else if (specType==3)  //fast fourier transform
      {
	double imag[nSmooth];

 	// The following code to form the PSD in units of yr^3 and 
	// frequency in 1/days (required for the whitening/postdarkening)
	// has been checked by G. Hobbs (30/10/08).
	// Note that a 2**N number of points is required for an FFT 
	// (this is not checked)
	//
	// Do not modify without leaving a comment here.
	// G. Hobbs 05 Jan 2009: Use TK_FFT instead of TK_four1

	for (i=0;i<nSmooth;i++)imag[i] = 0.0;
	TK_fft(1,nSmooth,smoothY,imag);
	nSpec = nSmooth/2;
	for (i=1;i<nSpec;i++)
	  {
	    specX[i-1] = i/(tspan);
	    specY[i-1] = pow(smoothY[i],2)+pow(imag[i],2);
	    if (specOut==1)
	      specY[i-1] = (specY[i-1]/pow(365.25*86400.0,2))*2*(tspan/365.25);///(double)nSmooth/(double)nSmooth;
	    else if (specOut==2) // Amplitude
	      specY[i-1]=sqrt(specY[i-1]*4.0);
	    else if (specOut==3) // Power
	      specY[i-1]=specY[i-1]*4.0;
	  }
	nSpec--;
	// End of checked version
      }
    else if (specType==4) // Weighted Lomb-Scargle
      {
	TK_weightLS(smoothX,smoothY,e,nSmooth,specX,specY,&nSpec,outY_re,outY_im);
	var=1;
	for (i=0;i<nSpec;i++)
	  {
	    if (specOut==1) // PSD
	      {
		//		specY[i]*=(2.0*var*nSpec);
		specY[i] = (specY[i]/pow(365.25*86400.0,2))*(tspan/365.25)/2.0; ///(double)nSpec;///(double)nSmooth/(double)nSmooth;
	      }
	    else if (specOut==2) // Amplitude
	      specY[i]=sqrt(specY[i]);
	    else if (specOut==3) // Power
	      specY[i]=specY[i]; //*2.0*var/nSpec;
	  }
      }
    else if (specType==5) // Fit sinusoids
      {
	TK_fitSinusoids(smoothX,smoothY,e,nSmooth,specX,specY,&nSpec);
	var=1;
	for (i=0;i<nSpec;i++)
	  {
	    if (specOut==1) // PSD
	      {
		//		specY[i]*=(2.0*var*nSpec);
		specY[i] = (specY[i]/pow(365.25*86400.0,2))*(tspan/365.25)/2.0; ///(double)nSpec;///(double)nSmooth/(double)nSmooth;
	      }
	    else if (specOut==2) // Amplitude
	      specY[i]=sqrt(specY[i]);
	    else if (specOut==3) // Power
	      specY[i]=specY[i]; //*2.0*var/nSpec;
	  }
      }
    // Post-darken if requested
    if (preWhite==1) 
      {
	for (i=0;i<nSpec;i++)
	  {
	    if (specX[i]!=0)
	      specY[i]/=(4.0*pow(sin(M_PI*specX[i]*(smoothX[1]-smoothX[0])),2));
	  }
      }
    else if (preWhite==2)
      {
	for (i=0;i<nSpec;i++)
	  {
	    if (specX[i]!=0.0)
	      specY[i]/=pow(4.0*pow(sin(M_PI*specX[i]*(smoothX[1]-smoothX[0])),2),2);
	  }
      }
      
    if (output==6)
      {
	for (i=0;i<nSpec;i++)
	  printf("spectrum: %d %g %g\n",i,specX[i],specY[i]);
	exit(1);
      }
  }

  // Return the spectrum
  for (i=0;i<nSpec;i++)
    {
      outX[i] = specX[i];
      outY[i] = specY[i];
    }
  *nout = nSpec-1;  

  return 0.0;
}

/* PUTBACK: removed TK_fitSine */


void TKfirstDifference(double *x,double *y,int n)
{
  int i;
  float sy[n];
  for (i=0;i<n;i++)
    sy[i] = y[i];

  for (i=0;i<n;i++)
    {
      if (i!=0)
	y[i] = sy[i]-sy[i-1];
      else
	y[i] = 0.0;
    }
}

/* TK_fitSinusoids attempts to fit a bunch of harmonically related sinusoids to the data. */

void TK_fitSinusoids(double *x,double *y,double *sig,int n,double *outX,double *outY,int *outN)
{
  double omega0,chisq;
  int i,j;
  double **cvm,p[2],e[2];
  double recC[MAX_OBSN],recS[MAX_OBSN];
  double outE[MAX_OBSN];
  double pred;

  cvm = (double **)malloc(sizeof(double *)*2);
  for (i=0;i<2;i++) cvm[i] = (double *)malloc(sizeof(double)*2);

  omega0 = 2.0*M_PI/(TKrange_d(x,n)*(double)n/(double)(n-1));
  *outN = n/2;

  for (j=0;j<*outN;j++)
    {
      GLOBAL_OMEGA = omega0*(j+1);
      TKleastSquares_svd(x,y,sig,n,p,e,2,cvm,&chisq,sineFunc,1);
      outX[j] = GLOBAL_OMEGA/2.0/M_PI;
      outY[j] = p[0]*p[0]+p[1]*p[1];
      outE[j] = sqrt(4*p[0]*p[0]*e[0]*e[0]+4*p[1]*p[1]*e[1]*e[1]);
      recC[j] = p[0];
      recS[j] = p[1];
    }
  /*  for (i=0;i<n;i++)
    {
      printf("data: %g %g %g\n",x[i],y[i],sig[i]);
      } */

  /*  for (i=0;i<1000;i++)
    {
      pred=0.0;    
      for (j=0;j<*outN;j++)
	{
	  GLOBAL_OMEGA = omega0*(j+1);
	  pred += (recC[j]*cos(GLOBAL_OMEGA*(x[n-1]-x[0])/1000.0*i+x[0])+recS[j]*sin(GLOBAL_OMEGA*(x[n-1]-x[0])/1000.0*i+x[0]));
	}
      printf("series: %g %g\n",(x[n-1]-x[0])/1000.0*i+x[0],pred);
      } */
}

void sineFunc(double x,double *v,int ma)
{
  v[0] = sin(GLOBAL_OMEGA*x);
  v[1] = cos(GLOBAL_OMEGA*x);
}

//GEORGE'S ALGORITHM!!!! I think mine is better.     Calculates a weighted Lomb-Scargle periodogram
//void TK_weightLS(double *x,double *y,double *sig,int n,double *outX,double *outY,int *outN)
void TK_weightLS(double *x,double *y,double *sig,int n,double *outX,double *outY,int *outN, double *outY_re, double *outY_im)
{
  long double s1,s2,s3,s4,s5;
  //double recA[MAX_OBSN],recB[MAX_OBSN];
  double pred;
  double omega=0.0;
  double si,ci;
  double a,b;
  double omega0;
  double var;
  int i,j;
  
  //omega0 = 2.0*M_PI/TKrange_d(x,n);
  omega0 = 2.0L*M_PI/(TKrange_d(x,n) * (double)n / (double)(n-1));
  *outN = n/2 - 1;

  
  //USE THIS FOR AN UNWEIGHTED LEAST SQUARES FIT
  //printf("IGNORING ERROR BARS --> unweighted LSq Fit in TK_weightLS, all input errors are now equal to 1.0 (even in your own code - this is a side effect!!!\n");
  //for (i=0;i<n;i++)
  //  sig[i]=1.0;
  

  for (j=0;j<*outN;j++)
    {
      omega = omega0*(j+1);  //the DC term will be zero always since we explicitly removed a mean from the input data set in Step 3 above.
      s1 = s2 = s3 = s4 = s5 = 0.0;
      for (i=0;i<n;i++)
	{
	  si = sin(omega*x[i]);
	  ci = cos(omega*x[i]);
	  
	  s1 += (long double)(y[i]*si/sig[i]/sig[i]);
	  s2 += (long double)(si*si/sig[i]/sig[i]);
	  s3 += (long double)(si*ci/sig[i]/sig[i]);
	  s4 += (long double)(y[i]*ci/sig[i]/sig[i]);
	  s5 += (long double)ci*ci/sig[i]/sig[i];    //NB!!!!!!!!!!!!!!!! THIS IS NOT TYPECAST PROPERLY!!!! - DY
	}
      b = (double)((s4-s1/s2)/(s5-s3/s2));  //the real Fourier component (since it has y*cos) George's solution
      a = (double)((s1-b*s3)/s2);           //the imaginary Fourier component (since it has y*sin) George
      //recA[j] = a;
      outY_im[j] = a;
      //recB[j] = b;
      outY_re[j] = b;
      outX[j] = omega/2.0/M_PI;
      outY[j] = a*a+b*b;
    }

    //for (i=0;i<1000;i++)
    //{
    //  pred=0.0;
    //  for (j=0;j<*outN;j++)
    //{
    //  GLOBAL_OMEGA = omega0*(j+1);
    //  pred += (recB[j]*cos(GLOBAL_OMEGA*(x[n-1]-x[0])/1000.0*i+x[0])+recA[j]*sin(GLOBAL_OMEGA*(x[n-1]-x[0])/1000.0*i+x[0]));
    //}
    //printf("series2: %g %g\n",(x[n-1]-x[0])/1000.0*i+x[0],pred);
    //}
}



/* Calculates the discrete Fourier transform of a data-set */
/* Assumes that the x-values are in units of days and the  */
/* y-values are in seconds                                 */
/* This code was checked by G. Hobbs - 30/10/08            */
/* DO NOT MODIFY WITHOUT LEAVING COMMENTS HERE             */
/* 13th November: DY added capability for retrieval of imaginary part into main code */
void TK_dft(double *x,double *y,int n,double *outX,double *outY,int *outN, double *outY_re, double *outY_im)
{
  complexVal spec[n];
  double xo;
  int nfreq;
  int i,j;
  double f0;
  double range,t,f;

  range = TKrange_d(x,n)*(double)n/(double)(n-1);
  nfreq = n/2;

  for (i=0;i<n;i++)
    outX[i]=1.0/range*(i+1);

  for (i=0;i<n;i++)
    {
      spec[i].real = 0.0;
      spec[i].imag = 0.0;
      for (j=0;j<n;j++)
	{
	  t = x[j]*86400.0; // In seconds
	  f = outX[i]/86400.0; // In Hz
	  spec[i].real+=y[j]*cos(2*M_PI*f*t);
	  spec[i].imag+=y[j]*sin(2*M_PI*f*t);
	}
    }
  // Form spectrum
  for (i=0;i<n;i++)
    {
      outY_re[i] = spec[i].real;
      outY_im[i] = spec[i].imag;
      outY[i] = (pow(spec[i].real,2)+pow(spec[i].imag,2));
    }
  *outN = n/2;
}

void TKaveragePts(double *x,double *y,int n,int width,double *meanX,double *meanY,int *nMean)
{
  double max,min;
  int i;
  meanX[*nMean] = 0.0;
  meanY[*nMean] = 0.0;

  max = x[0],min=x[0];
  for (i=0;i<n;i++)
    {
      if (max < x[i]) max=x[i];
      if (min > x[i]) min=x[i];
      meanY[*nMean]+=y[i];
      meanX[*nMean]+=x[i];
    }
  meanY[*nMean]/=(double)n;
  meanX[*nMean]/=(double)n;
  if (max-min < (double)width)
    {
      (*nMean)++;
      return;
    }  
  else
    {
      // Find largest gap in the data-points
      //      printf("Finding largest gap\n");
      double gap=-1;
      int   igap=-1;

      for (i=1;i<n;i++)
	{
	  if (gap < fabs(x[i]-x[i-1]))
	    {
	      gap = fabs(x[i]-x[i-1]);
	      igap=i;
	    }
	}
      if (igap==-1)
	{
	  printf("Error when averaging points x[0] = %g, x[1] = %g\n",x[0],x[1]);
	  exit(1);
	}
      //      printf("LEFT SIDE %d %d %g [%g %g]\n",n,igap,gap,x[igap],x[igap-1]);
      // Consider left hand side
      TKaveragePts(x,y,igap,width,meanX,meanY,nMean);
      // Consider right hand side
      //      printf("RIGHT SIDE %d %d %g [%g %g]\n",n,igap,gap,x[igap],x[igap-1]);
      TKaveragePts(x+igap,y+igap,(n-igap),width,meanX,meanY,nMean);            
    }
}

/* This sorting routine should be replaced by something better in toolkit.h  */
void TKsortit(double *x,double *y,int n)
{
  int sort=0,i;
  double t1,t2;
  do {
    sort=0;
    for (i=0;i<n-1;i++)
      {
	if (x[i]>x[i+1])
	  {
	    sort=1;
	    t1=x[i];
	    t2=y[i];
	    x[i]=x[i+1];
	    y[i]=y[i+1];
	    x[i+1] = t1;
	    y[i+1] = t2;
	    break;
	  }
      }
  } while (sort==1);
}


/* ************************************************************** *
 * Smoothing routines:                                            *
 *                                                                * 
 * These routines input an irregularly sampled time series        * 
 * and output a smoothed version of the time series               *
 *                                                                * 
 * boxcar                                                         * 
 * hann - also called raised cosine smoothing                     *
 *                                                                * 
 * ************************************************************** */


/* Boxcar smoothing 
 * The width of the box-car must be an odd integer
 */
void TKboxcar(double *x,double *y,int n,double *ox,double *oy,int *on,int width)
{
  int i,j;
  double sy;
  int nStore;
  *on = 0;
  if (width/2.0 == (int)width/2)
    printf("WARNING in TKboxcar - should use odd number for the width\n");

  for (i=(width-1)/2;i<n-(width-1)/2;i++)
    {
      sy = 0.0;
      nStore=0;
      for (j=i-(width-1)/2;j<i+(width-1)/2;j++)
	{
	  sy += y[j];
	  nStore++;
	}
      ox[*on] = x[i];
      oy[*on] = sy/(double)nStore;
      (*on)++;
    }
}

/* 
 * Hann (or raised cosine) smoother
 */
void TKhann(double *x,double *y,int n,double *ox,double *oy,int *on,int width)
{
  int i,j;
  double sy;
  int nStore;
  *on = 0;
  for (i=(width-1)/2;i<n-(width-1)/2;i++)
    {
      sy = 0.0;
      nStore=0;
      for (j=i-(width-1)/2;j<i+(width-1)/2;j++)
	{
	  sy += 0.5*(1.0-cos(2*M_PI*(j-(i-(width-1)/2))/(double)width))*y[j];
	  nStore++;
	}
      ox[*on] = x[i];
      // Normalise so that the area under the Hann smoother = 1
      oy[*on] = sy/((double)width/2.0);
      (*on)++;
    }
}

/* ************************************************************** *
 * Spline routines                                                *
 * ************************************************************** */
#define ABS(x) ((x) < 0 ? -(x) : (x))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))


void TKcmonot (int n, double x[], double y[], double yd[][4])
/*****************************************************************************
compute cubic interpolation coefficients via the Fritsch-Carlson method,
which preserves monotonicity
******************************************************************************
Input:
n		number of samples
x  		array[n] of monotonically increasing or decreasing abscissae
y		array[n] of ordinates

Output:
yd		array[n][4] of cubic interpolation coefficients (see notes)
******************************************************************************
Notes:
The computed cubic spline coefficients are as follows:
yd[i][0] = y(x[i])    (the value of y at x = x[i])
yd[i][1] = y'(x[i])   (the 1st derivative of y at x = x[i])
yd[i][2] = y''(x[i])  (the 2nd derivative of y at x = x[i])
yd[i][3] = y'''(x[i]) (the 3rd derivative of y at x = x[i])

To evaluate y(x) for x between x[i] and x[i+1] and h = x-x[i],
use the computed coefficients as follows:
y(x) = yd[i][0]+h*(yd[i][1]+h*(yd[i][2]/2.0+h*yd[i][3]/6.0))

The Fritsch-Carlson method yields continuous 1st derivatives, but 2nd
and 3rd derivatives are discontinuous.  The method will yield a
monotonic interpolant for monotonic data.  1st derivatives are set
to zero wherever first divided differences change sign.

For more information, see Fritsch, F. N., and Carlson, R. E., 1980,
Monotone piecewise cubic interpolation:  SIAM J. Numer. Anal., v. 17,
n. 2, p. 238-246.

Also, see the book by Kahaner, D., Moler, C., and Nash, S., 1989,
Numerical Methods and Software, Prentice Hall.  This function was
derived from SUBROUTINE PCHEZ contained on the diskette that comes
with the book.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 09/30/89
Modified:  Dave Hale, Colorado School of Mines, 02/28/91
	changed to work for n=1.
Modified:  Dave Hale, Colorado School of Mines, 08/04/91
	fixed bug in computation of left end derivative
*****************************************************************************/
{
  int i;
  double h1,h2,del1,del2,dmin,dmax,hsum,hsum3,w1,w2,drat1,drat2,divdf3;

  /* copy ordinates into output array */
  for (i=0; i<n; i++)
    yd[i][0] = y[i];
  
  /* if n=1, then use constant interpolation */
  if (n==1) {
    yd[0][1] = 0.0;
    yd[0][2] = 0.0;
    yd[0][3] = 0.0;
    return;
    
    /* else, if n=2, then use linear interpolation */
  } else if (n==2) {
    yd[0][1] = yd[1][1] = (y[1]-y[0])/(x[1]-x[0]);
    yd[0][2] = yd[1][2] = 0.0;
    yd[0][3] = yd[1][3] = 0.0;
    return;
  }
  
  /* set left end derivative via shape-preserving 3-point formula */
  h1 = x[1]-x[0];
  h2 = x[2]-x[1];
  hsum = h1+h2;
  del1 = (y[1]-y[0])/h1;
  del2 = (y[2]-y[1])/h2;
  w1 = (h1+hsum)/hsum;
  w2 = -h1/hsum;
  yd[0][1] = w1*del1+w2*del2;
  if (yd[0][1]*del1<=0.0)
    yd[0][1] = 0.0;
  else if (del1*del2<0.0) {
    dmax = 3.0*del1;
    if (ABS(yd[0][1])>ABS(dmax)) yd[0][1] = dmax;
  }

  /* loop over interior points */
  for (i=1; i<n-1; i++) {
    
    /* compute intervals and slopes */
    h1 = x[i]-x[i-1];
    h2 = x[i+1]-x[i];
    hsum = h1+h2;
    del1 = (y[i]-y[i-1])/h1;
    del2 = (y[i+1]-y[i])/h2;
    
    /* if not strictly monotonic, zero derivative */
    if (del1*del2<=0.0) {
      yd[i][1] = 0.0;
      
      /*
       * else, if strictly monotonic, use Butland's formula:
       *      3*(h1+h2)*del1*del2
       * -------------------------------
       * ((2*h1+h2)*del1+(h1+2*h2)*del2)
       * computed as follows to avoid roundoff error
       */
    } else {
      hsum3 = hsum+hsum+hsum;
      w1 = (hsum+h1)/hsum3;
      w2 = (hsum+h2)/hsum3;
      dmin = MIN(ABS(del1),ABS(del2));
      dmax = MAX(ABS(del1),ABS(del2));
      drat1 = del1/dmax;
      drat2 = del2/dmax;

      yd[i][1] = dmin/(w1*drat1+w2*drat2);
    }
  }
  
  /* set right end derivative via shape-preserving 3-point formula */
  w1 = -h2/hsum;
  w2 = (h2+hsum)/hsum;
  yd[n-1][1] = w1*del1+w2*del2;
  if (yd[n-1][1]*del2<=0.0)
    yd[n-1][1] = 0.0;
  else if (del1*del2<0.0) {
    dmax = 3.0*del2;
    if (ABS(yd[n-1][1])>ABS(dmax)) yd[n-1][1] = dmax;
  }
  
  /* compute 2nd and 3rd derivatives of cubic polynomials */
  for (i=0; i<n-1; i++) {
    h2 = x[i+1]-x[i];
    del2 = (y[i+1]-y[i])/h2;
    divdf3 = yd[i][1]+yd[i+1][1]-2.0*del2;
    yd[i][2] = 2.0*(del2-yd[i][1]-divdf3)/h2;
    yd[i][3] = (divdf3/h2)*(6.0/h2);
  }
  yd[n-1][2] = yd[n-2][2]+(x[n-1]-x[n-2])*yd[n-2][3];
  yd[n-1][3] = yd[n-2][3];
}

// Interpolate the spline fit on to a given set of x-values
void TKspline_interpolate(int n,double *x,double *y,double yd[][4],double *interpX,
			  double *interpY,int nInterp)
{
  double h;
  int jpos;
  int i,j;

  //To evaluate y(x) for x between x[i] and x[i+1] and h = x-x[i],
  //use the computed coefficients as follows:
  //y(x) = yd[i][0]+h*(yd[i][1]+h*(yd[i][2]/2.0+h*yd[i][3]/6.0))
  
  // Assume the data are sorted so x[i] < x[i+1]
  for (i=0;i<nInterp;i++)
    {
      jpos=-1;
      for (j=0;j<n-1;j++)
	{
	  if (interpX[i]>=x[j] && interpX[i]<x[j+1])
	    {
	      jpos=j;
	      break;
	    }
	}
      if (jpos!=-1)
	{
	  h = interpX[i]-x[jpos];
	  interpY[i] = yd[jpos][0]+h*(yd[jpos][1]+h*(yd[jpos][2]/2.0+h*yd[jpos][3]/6.0));
	}
      else
	interpY[i] = 0.0;
    }
}

void TKlomb_d(double *x,double *y,int n,double ofac,double hifac,double *ox,double *oy,int *outN,double *var)
{ 
  double store;
  double freq0;         /* Initial frequency */
  int i,j;
  double arg;
  double sum1,sum2,sum3,sum4;
  double tau;
  double omegaTau,omega;
  double aveX;
  double ybar;
  double alpha[n],delta,beta[n]; /* See equation 5.5.7 */
  double ctheta,stheta,theta,omega0;
  double cosTau,sinTau,st,ct;
  double cosTheta[n],sinTheta[n];
  double scale;
  double y2;

  ybar  = TKmean_d(y,n);
  *var  = TKvariance_d(y,n);

  // Use to get the L-S periodogram to give the same results as the DFT
  scale = (double)n/(double)(n-1);
  //  scale = 1;
  freq0 = 1.0/(TKrange_d(x,n)*ofac*scale);
  aveX  = (TKfindMin_d(x,n)+TKfindMax_d(x,n))/2.0;
  *outN = (int)(0.5*ofac*hifac*n);

  /* Find initial value of cos and sin theta for the trigonometric recurrence */
  for (i=0;i<n;i++)
    {
      omega0      = 2.0*M_PI*freq0;
      theta       = omega0*(x[i]-aveX);
      alpha[i]    = -2.0*sin(theta/2.0)*sin(theta/2.0);
      sinTheta[i] = beta[i] = sin(theta);
      cosTheta[i] = cos(theta);
    }

  for (j=0;j<*outN;j++)
    {
      ox[j] = freq0*(j+1);
      sum1 = 0.0; sum2 = 0.0;
      for (i=0;i<n;i++)
	{
	  st = sinTheta[i];
	  ct = cosTheta[i];
	  sum1  += st*ct;
	  sum2  += (ct-st)*(ct+st);
	}

      omegaTau = 0.5*atan2(2.0*sum1,sum2);
      cosTau = cos(omegaTau);
      sinTau = sin(omegaTau);
      sum1 = sum2 = sum3 = sum4 = 0.0;
      for (i=0;i<n;i++)
	{
	  y2 = y[i]-ybar;
	  ct = cosTheta[i];
	  st = sinTheta[i];

	  sum1 += y2*(ct*cosTau+st*sinTau);
	  sum2 += y2*(st*cosTau-ct*sinTau); 
	  sum3 += pow((ct*cosTau)+st*sinTau,2); 
	  sum4 += pow((ct*sinTau)-st*cosTau,2); 

	  /* Now update the recurrence relations */
	  cosTheta[i]+=(alpha[i]*ct-beta[i]*sinTheta[i]);
	  sinTheta[i]+=(alpha[i]*sinTheta[i]+beta[i]*ct);
	}
      oy[j] = 0.5/(*var)*(sum1*sum1/sum3 + sum2*sum2/sum4);
    }
}

// Fourier transform
// Requires npts = 2**M value - does not check for this
// Routine modified from http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/dft
/*
   This computes an in-place complex-to-complex FFT 
   x and y are the real and imaginary arrays of n = 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform 
*/
int TK_fft(short int dir,long n,double *x,double *y)
{
   long i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;
   int m;
   
   m = (int)(log(n)/log(2)+0.1);
   printf("m = %d, n = %d %g %g %g %d\n",m,n,log(n),log(2),log(n)/log(2),(int)(log(n)/log(2)+0.1));
   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   //   printf("Here\n");
   for (i=0;i<n-1;i++) {
     //     printf("I = %d, n = %d\n",i,n);
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
	//	printf("k = %d %d\n",k,j);
         j -= k;
         k >>= 1;
      }
      j += k;
   }
   //   printf("Here 2\n");
   /* Compute the FFT */
   c1 = -1.0; 
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0; 
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1; 
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1) 
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }
   printf("Here 3\n");
   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
   }
   
   return 0;
}


// Routines for sigmaz calculation
/* Direct copy of the original Fortran */
/* using identical fitting routines etc. */
/* Note that all arrays start from 1 */
//void calcSigmaz(pulsar *psr,int npsr,int time,int psrfile,double mintau,
//		int nSel,int *selPsr,int weights,int style,int average,
//		float *avX,float *avY,int *nav)
void TKcalcSigmaz(pulsar psr,int weights,double *ret_tau,double *ret_szbias,double *ret_e1,
		double *ret_e2,int *ret_nval,double mintau)
{
  int n,nbad,nfit,nread,idiag,niter,ntable,ios,nsigs,nuncor,ntot;
  int ndostats,nptot,nuncorbins,ntries;
  double x,szlog,sz,cov4,rootsz,p4,sumwt,sumwtuncor,error,scale,bias;
  double utjd2nuncor,wt,dplus,dminus,errlog,chidf,avewt,scalelog;
  double taulog1,addvsq,szbias,dp,dm,tabledata[101];
  int endit;
  int   j;

  *ret_nval=0;
  nxunits = 0;
  nusewt = weights; 
  tmin =-1e50;
  tmax = 1.0e50;
  xmin =-1.0e50;
  xmax = 1.0e50;
  root2 = sqrt(2.0);
  secyear = 86400.0*365.25;
  utfirst = 1.0e50;
  utlast = -1.0e50;
  npt1last = 0;
  npt2last = 0;
  utjd1 = 0.0;
  ndim = 0;
  if (nxunits >= 0) ndim=4;
  if (nxunits <= -1 && nxunits >= -3) ndim=3;
  if (nxunits <= -4 && nxunits >= -6) ndim=1;

  int nobs = psr.nobs;
  double t2times[MAX_OBSN],t2resid[MAX_OBSN],t2error[MAX_OBSN],t0;
  int i;
  scale = 2.0*sqrt(5.0);
  scalelog = log10(scale);
  for (n=1;n<=100;n++)
    tabledata[n] = 0.0;
  
  readin(psr);
  tausec = taumax;
  tauday = tausec/86400.0;
  ntau = 0;
  ncubics=0;
  if (nusewt >= 0)
    {
      tauyear = tauday/365.25;
      utjd1 = utjd[indx[1]]-1.0e-2;
      utjd2 = utjd[indx[npt]]+1.0e-2;
      utmean = (utjd1+utjd2)/2.0;
      for (niter=1;niter<=nusewt+1;niter++)
	{
	  fit4(&nfit,&p4,&cov4,ndostats,&chidf,&avewt);
	  addvsq=addvar*addvar+(chidf-1.0)/avewt;
	  if (addvsq > 0.0)
	    addvar=sqrt(addvsq);
	  else
	    {
	      if (nusewt != 0 && nusewt != -3) addvar = 0.0;
	    }
	}
    }
  // Should check nbintype
  ndostats=0;

  // 110 continue statement
  do {
    ntau++;
    ncubic = 0;
    tauday = tausec/86400.0;
    tauyear = tauday/365.25;
    taulog = log10(tausec);
    sumwt = 0.0;
    sz=0.0;
    utjd2nuncor = 0.0;
    sumwtuncor = 0.0;
    nuncorbins = 0;
    ntot = 0;
    ntries = 0;
    nptot = 0;
    npt1last = 0;
    npt2last = 0;
    
    if (nbintype==0) utjd1=utjd[indx[1]]-tausec;
    if (nbintype==1) utjd1=utjd[indx[1]]-taumin;
    if (nbintype==2) utjd1=utjd[indx[1]]-tausec;
    endit=0;

    do {
      //      printf("utjd1 = %g\n",utjd1);
      if (nbintype==0) utjd1+=tausec;
      if (nbintype==1) utjd1+=taumin;
      if (nbintype==2) utjd1+=tausec;
      if (utjd1 < utjd[indx[npt]])
	{
	  utjd2=utjd1+tausec;
	  utmean=(utjd1+utjd2)/2.0;
	  if (utjd1 >= utjd2nuncor)
	    {
	      nuncor=1;
	      utjd2nuncor=utjd2;
	    }
	  else
	    {
	      // Should check overlap
	      nuncor=0;
	    }
	  ncubics++;		  
	  fit4(&nfit,&p4,&cov4,ndostats,&chidf,&avewt);
	  ntries++;
	  if (nfit==1)
	    {
	      wt=1.0/cov4;
	      sz=sz+wt*p4*p4;
	      sumwt=sumwt+wt;
	      ntot++;
	      nptot=nptot+npt2last-npt1last+1;
	      if (nuncor==1)
		{
		  sumwtuncor=sumwtuncor+wt;
		  nuncorbins=nuncorbins+1;
		}	     
	      //	      printf("Here %g %g %d %d\n",(double)sumwtuncor,(double)wt,(int)nuncorbins,nuncor);
	    }
	}
      else
	endit=1;
    } while (endit==0); // 130 loop
    //    printf("Add or not (Exit): %d %g %g\n",ntot,sz,sumwt);
    if (ntot==0 || sz <=0.0 || sumwt<=0.0)
      {
      }
    else
      {
	rootsz = tausec*tausec*sqrt(sz/sumwt)/scale;
	szlog = log10(rootsz);
	errlog=0.0;
	bias=0.0;
	dp=0.0;
	dm=0.0;
	if (sumwtuncor > 1.0e-50)
	  {
	    error=nuncorbins/sumwtuncor;
	    if (error > 0.0) errlog=2.0*taulog+0.5*log10(error)-scalelog;
	    x=nuncorbins;
	    bias=0.17/x;
	    dm = 0.31/sqrt(x);
	    if (nuncorbins > 1)
	      dp=0.31/sqrt(x-1.0);
	    else
	      dp=0.52;
	    
	  }
	szbias = szlog-bias;
	dplus = szbias + dp;
	dminus = szbias - dm;

	ret_szbias[*ret_nval] = szbias;
	ret_tau[*ret_nval]    = log10(tauyear);
	ret_e1[*ret_nval]     = dp;
	ret_e2[*ret_nval]     = -dm;
	//	printf("Have %g %g %g %g\n",szbias,log10(tauyear),dp,-dm);

	//	exit(1);
	(*ret_nval)++;
      }
    if (nbintype==0) tausec=tausec/2.0;
    if (nbintype==1) tausec=tausec-taumin;
    if (nbintype==2) tausec=tausec+taumin;
    //    exit(1);
    //    printf("Exited (testing) with %d %g %g %g %g %d\n",ntot,tausec,tauensure,tausec/86400.0/365.25,mintau,*ret_nval);
    //	fflush(stdout);
  } while (ntot > 0 || (tausec >= tauensure && tauensure >= 0.0) || (tausec/86400.0/365.25 > mintau));
  //  printf("Exited with %d %g %g %g %g %d\n",ntot,tausec,tauensure,tausec/86400.0/365.25,mintau,*ret_nval);
  //	fflush(stdout);
}



/* WHAT AM I RETURNING IN AND OUT OF THIS FUNCTION? */
void fit4(int *nfit,double *p4,double *cov4,int ndostats,double *chidf,double *avewt)
{
  int n,i,j,nbad,nok,nf,npt1,npt2,m,k;
  static int npass=0;
  static int ntaulast=-1;
  double am[21][21],aminv[21][21],fvec[21],determ,y,ylast,x,ave,rms,alv,yf;
  double res,sumwt,par[5],wt,chisqu,degfree,utave,err4;
  static double scalefact=1.0;
  static double scalecon=1.0;

  npass++;
  ncubic++;
  ncubics++;
  if (ndim > 1) scalefact = pow(86400.0,ndim-1);
  alv = 0.0;
  for (i=1;i<=ndim;i++)
    {
      fvec[i] = 0.0;
      par[i] = 0.0;
      for (j=1;j<=ndim;j++)
	am[i][j]=0.0;
    }
  sumwt=0.0;
  npt1=0;
  npt2=0;
  *nfit=0;
  if (ndim != 4)
    {
      if (ntunits >= -2)
	scalecon=3.0;
      else
	{
	  if (ndim==3) scalecon = 3.0*tdiffmin;
	  if (ndim==2) scalecon = 6.0*tdiffmin*tdiffmin;
	  if (ndim==1) scalecon = 6.0*tdiffmin*tdiffmin*tdiffmin;
	}
    }
  for (n=1;n<=npt;n++)
    {
      if (utjd[indx[n]] >= utjd1-1.0e-5)
	{
	  if (utjd[indx[n]]+(4-ndim)*tdiffmin > utjd2+1.0e-5) goto pos90;
	  if (npt2 == 0) npt1=n;
	  npt2=n;
	}
    }
 pos90:
  if (npt1 == 0 || npt2 == 0 || npt2-npt1+1 < ndim || npt2 == 0 || (npt1 == npt1last && npt2 == npt2last) 
      || (utjd[indx[npt2]]+(4.0-ndim)*tdiffmin-utjd[indx[npt1]] < tausec/root2))
    {
      goto pos1000;
    }
  npt1last = npt1;
  npt2last = npt2;
  for (n=npt1;n<=npt2;n++)
    {
      getweights(n,&wt);
      sumwt=sumwt+wt;
      getprtj(n);
      y = data[indx[n]];
      for (i=1;i<=ndim;i++)
	{
	  fvec[i] = fvec[i] + y*prtl[i]*wt;
	  for (j=1;j<=ndim;j++)
	    {
	      am[i][j]=am[i][j]+prtl[i]*prtl[j]*wt;
	    }
	}
    }
  mat20(am,aminv,ndim,&determ,&nbad);
  if (nbad != 0.0)
    {
      goto pos1000;     
    }
  for (i=1;i<=ndim;i++)
    {
      for (j=1;j<=ndim;j++)
	par[i]=par[i]+aminv[i][j]*fvec[j];
    }
  *nfit=1;
  *p4=par[ndim]/scalefact;
  *cov4=aminv[ndim][ndim]/scalefact/scalefact;
  if (ndim!=4)
    {
      *p4=(*p4)/scalecon;
      *cov4=(*cov4)/scalecon/scalecon;
    }
  utave=(utjd1/86400.0+utjd2/86400.0)/2.0;
  err4=sqrt(*cov4);
  if (ntau!=ntaulast) ntaulast=ntau;
  ave=0.0;
  nok=0;
  sumwt=0.0;
  chisqu=0.0;
  for (n=npt1;n<=npt2;n++)
    {
      yf = 0.0;
      getweights(n,&wt);
      getprtj(n);
      for (i=1;i<=ndim;i++)
	yf=yf+prtl[i]*par[i];
      res=data[indx[n]]-yf;
      ave=ave+wt*res;
      chisqu=chisqu+wt*res*res;
      sumwt=sumwt+wt;
      nok++;
      if (nok>1) alv=alv+(res-ylast)*(res-ylast);
      ylast=res;
    }
  // 401 Continue
  if (nok > 1 && sumwt > 1.0e-50)
    {
      ave = ave/sumwt;
      rms = chisqu/sumwt-ave*ave;
      if (rms > 0.0) rms=sqrt(rms);
      if (alv > 0.0) alv=sqrt(alv/2.0/(nok-1.0));
      *avewt=sumwt/nok;
      degfree=nok-ndim;
      if (nok > 4)
	*chidf = chisqu/degfree;
      else
	*chidf=1.0;
    }
 pos1000:
  return;
}

void mat20(double sam[21][21],double a[21][21],int n,double *determ,int *nbad)
{
  int index[21][3],ipivot[21];
  double pivot[21],amax,swap,x,t;
  int i,j,k,l,l1,irow,icolum,jrow,jcolum;

  *nbad=0;
  for (i=1;i<=n;i++)
    {
      for (j=1;j<=n;j++)
	a[i][j]=sam[i][j];
    }
  *determ=0.0;
  for (j=1;j<=n;j++)
    ipivot[j]=0;
  for (i=1;i<=n;i++) // End 550
    {
      amax=0.0;
      for (j=1;j<=n;j++) // End 105
	{
	  if (ipivot[j]-1 != 0) // CHECK THIS <<<
	    {
	      for (k=1;k<=n;k++)
		{
		  if (ipivot[k]-1 < 0) goto pos_mat80;
		  if (ipivot[k]-1 == 0) goto pos_mat100;
		  if (ipivot[k]-1 > 0) goto pos_mat750;
		pos_mat80:
		  if (fabs(amax)-fabs(a[j][k]) < 0)
		    {
		      irow=j;
		      icolum=k;
		      amax=a[j][k];
		    }
		pos_mat100:
		  {
		  }
		} // 100
	    }
	} // 105
      ipivot[icolum]=ipivot[icolum]+1;
      if (irow-icolum != 0)
	{
	  for (l=1;l<=n;l++)
	    {
	      swap = a[irow][l];
	      a[irow][l]=a[icolum][l];
	      a[icolum][l]=swap;
	    }
	} // 260
      index[i][1]=irow;
      index[i][2]=icolum;
      pivot[i]=a[icolum][icolum];
      x=fabs(pivot[i]);
      if (x<=1.0e-70 || x>=1e70)
	{
	  *nbad=1;
	  return;
	}
      *determ=*determ+log10(fabs(pivot[i]));
      a[icolum][icolum]=1.0;
      for (l=1;l<=n;l++)
	a[icolum][l]=a[icolum][l]/pivot[i];
      for (l1=1;l1<=n;l1++)
	{
	  if (l1-icolum!=0)
	    {
	      t=a[l1][icolum];
	      a[l1][icolum]=0.0;
	      for (l=1;l<=n;l++)
		a[l1][l]=a[l1][l]-a[icolum][l]*t;
	    }
	}
    }// End 550

  for (i=1;i<=n;i++)
    {
      l=n+1-i;
      if ((index[l][1]-index[l][2])!=0.0)
	{
	  jrow=index[l][1];
	  jcolum=index[l][2];
	  for (k=1;k<=n;k++)
	    {
	      swap=a[k][jrow];
	      a[k][jrow]=a[k][jcolum];
	      a[k][jcolum]=swap;
	    }
	}
    }
 pos_mat750:
  return;
}

void getprtj(int n)
{
  double x;
  x = (utjd[indx[n]]-utmean)/86400.0;
  prtl[1] = 1.0;
  prtl[2] = x;
  prtl[3] = x*x;
  prtl[4] = x*x*x;
}

void getweights(int n, double *wt)
{
  if (nusewt == 0) *wt = 1.0/addvar/addvar;
  if (nusewt == -1) *wt = 1.0/pow(sigmai[indx[n]],2);
  if (nusewt == -2 || nusewt > 0) *wt = 1.0/(pow(sigmai[indx[n]],2)+pow(addvar,2));
  if (nusewt == -3) *wt = 1.0/addvar/addvar;
}

void readin(pulsar psr)
{
  int nalv,nread,n;
  static int idiag = 0;
  double tlast,ave,rms,var,sumwt,alv,t,x,sig,wt,t1,xr,xi,tautest;
  double xold,rootvar,permin;
  static double unitfact = 1.0;
  double mean=0.0;
  int i;
  int loop=0;

  for (i=0;i<psr.nobs;i++)
    mean+=psr.obsn[i].residual;
  mean/=(double)psr.nobs;

  npt = 0;
  utjdlast=-1.0; 
  ave=0.0;
  rms=0.0;
  var=0.0;
  sumwt=0.0;
  alv=0.0;
  nalv=0;
  tlast=0.0;
  sig=1.0;
  wt=1.0;

  // Should check ntunits and nxunits - assume both are zero
  ntunits=0;

  tauensure = -1.0;
  tauensure = tauensure*86400.0;
  unitfact = 1.0e-6;

  // 10 continue statement
  int ic=0;
  int breakit=0;

  do {
    do {
      do {
	// Assuming that nformat = 0
	
	if (nusewt == 0 || nusewt == -3)
	  {
	    t = (double)psr.obsn[ic].sat;
	    x = (double)(psr.obsn[ic].residual-mean)/1e-6;
	    if (ic>0) {
	      if (psr.obsn[ic-1].sat > psr.obsn[ic].sat)
		{
		  printf("1: ERROR: The TOAs need to be time sorted, sorry! - sats %g and %g\n",(double)psr.obsn[ic-1].sat,(double)psr.obsn[ic].sat);
		  //		  exit(1);
		}
	    }
	    ic++;
	    if (ic == psr.nobs) breakit=1;
	  }
	else
	  {
	    t = (double)psr.obsn[ic].sat;
	    x = (double)(psr.obsn[ic].residual-mean)/1e-6;
	    sig = (double)(psr.obsn[ic].toaErr);
	    if (ic>0) {
	      if (psr.obsn[ic-1].sat > psr.obsn[ic].sat)
		{
		  printf("2: ERROR: The TOAs need to be time sorted, sorry! sats: %.15g and %.15g\n",(double)psr.obsn[ic-1].sat,(double)psr.obsn[ic].sat);
		  //		  exit(1);
		}
	    }
	    ic++;
	    if (ic == psr.nobs) breakit=1;
	  }
	nread++;
      } while ((t < tmin || t > tmax) && breakit==0);
      if (npt == 0) t1=t;
      if (t < utfirst) utfirst=t;
      if (t > utlast) utlast=t;
      t=t-t1; // to avoid losing precision
      if (ntunits==0) t=t*86400.0;
    } while ((x < xmin || x > xmax) && breakit==0);
    x = unitfact*x;
    if (nusewt != 0 && nusewt != -3) sig = unitfact*sig;
    npt++;
    indx[npt] = npt;
    utjd[indx[npt]]=t;
    data[indx[npt]]=x;
    sigmai[indx[npt]]=sig;
    utjdlast = t;
    wt = 1.0/sig/sig;
    ave=ave+wt*x;
    var=var+wt*x*x;
    sumwt=sumwt+wt;
    if (npt > 1)
      {
	nalv++;
	alv = alv+pow(x-xold,2);
      }
  } while (breakit==0);
  // Ave different because I have a better mean removal or something similar
  if (npt <= 0) exit(1);
  indexx8(npt,utjd,indx);
  tdiffmin = utjd[indx[2]]-utjd[indx[1]];
  for (n=1;n<=npt-1;n++)
    {
      if (tdiffmin > utjd[indx[n+1]]-utjd[indx[n]])
	tdiffmin = utjd[indx[n+1]]-utjd[indx[n]];
    }
  ave=ave/sumwt;
  var=var/sumwt;
  rootvar=sqrt(var);
  rms=var-ave*ave;
  rms = sqrt(rms);
  if (nalv > 0 && alv > 0.0) alv=sqrt(alv/2.0/nalv);
  permax = utjd[indx[npt]]-utjd[indx[1]];
  permin = permax/(npt-1.0);
  taumin = permin;
  taumax = permax;
  if (tauensure > 0.0) {
    taumax = tauensure;
    //105 CONTINUE
    loop=0;
    do
      {
	if (taumax < permax)
	  {
	    loop=1;
	    taumax=2.0*taumax;
	  }
      }while (loop==1);
  }
  // Should check nbintype
  tautest = taumax;
  do {
    tautest=tautest/2.0;
  } while (tautest > taumin);
  taumin=tautest;

  if (nusewt >= 0.0) addvar = rms;
}


// Check conversion from Fortran carefully
void indexx8(int n,double *arrin,int *indx)
{
  //     from numerical recipes
  double q;
  int j,l,ir,i=0,indxt;

  for (j=1;j<=n;j++)
    indx[j] = j;
  l = n/2+1;
  ir = n;

 pos10:
  if (l > 1){
    l--;
    indxt = indx[l];
    q=arrin[indxt];
  } else {
    indxt=indx[ir];
    q = arrin[indxt];
    indx[ir]=indx[1];
    ir--;
    if (ir == 1){
      indx[1] = indxt;
      return;
    }
  }
  i=l;
  j=l+l;
 pos20:
 if (j <= ir) {
    if (j < ir) {
      if (arrin[indx[j]] < arrin[indx[j+1]])j++;
    }
    if (q < arrin[indx[j]]){
      indx[i] = indx[j];
      i=j;
      j=j+j;
    } else {
      j=ir+1;
    }
    goto pos20; // Sorry about this - copying the fortran!
 }
 indx[i]=indxt;
 goto pos10;
}
//Adapted from Stefan / George's plugin
//interpolation (spline): this function interpolates a data set using constrained spline onto an input set of interpX and nInterp values
void TKinterpolateSplineSmoothFixedXPts(double *inX, double *inY, int inN, double *interpX, double *interpY, int nInterp)
{
    //array needed by TKcmonot
    double yd[MAX_OBSN][4];

    //auxilary 'i'
    int i;

    double tempX[MAX_OBSN];
    int nTemp = nInterp;

    for (i=0;i<nTemp;i++)
      {
	tempX[i] = interpX[i];
      }

    TKcmonot(inN, inX, inY, yd);

    //only need to determine what interpY is
    TKspline_interpolate(inN, inX, inY, yd, tempX, interpY, nTemp);
} //interpolateSplineSmoothFixedXPts

// Spectral analysis using covariance matrix
// note: uinv array must start from 0, not 1
// NEW FEATURE:
// set nfit < 0 to automatically set it to nres/2-1
int calcSpectraErr(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY,double* specE,int nfit)
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

  psr = (pulsar *)malloc(sizeof(pulsar));
  initialiseOne(psr,0,1);

  cvm = (double **)malloc(sizeof(double *)*nfit);
  for (i=0;i<nfit;i++)
    cvm[i] = (double *)malloc(sizeof(double)*nfit);

  printf("Entering calcSpectra\n");

  // Should fit independently to all frequencies
  for (i=0;i<nres;i++)
    {
      sig[i] = 1.0; // The errors are built into the uinv matrix
      ip[i] = 0;
    }
  for (k=0;k<nfit;k++)
    {
	    if(verbose_calc_spectra){
		    printf("%5.1f%%   \r",(double)k/(double)nfit*100.0);
		    fflush(stdout);
	    }
      GLOBAL_OMEGA = 2.0*M_PI/((resx[nres-1]-resx[0])*(double)nres/(double)(nres-1))*(k+1);
      //      TKleastSquares_svd_psr_dcm(resx,resy,sig,nres,param,error,3,cvm,&chisq,fitMeanSineFunc,0,psr,1.0e-40,ip,uinv);
      //      v[k] = (resx[nres-1]-resx[0])/365.25/2.0*(pow(param[1],2)+pow(param[2],2))/pow(365.25*86400.0,2); 
      //      specX[k] = GLOBAL_OMEGA/2.0/M_PI;
      //      specY[k] = v[k];
      //      if(specE!=NULL){
      //	      specE[k] = (resx[nres-1]-resx[0])/365.25/2.0*(pow(error[1],2)+pow(error[2],2))/pow(365.25*86400.0,2);
      //      }
      TKleastSquares_svd_psr_dcm(resx,resy,sig,nres,param,error,2,cvm,&chisq,fitCosSineFunc,0,psr,1.0e-40,ip,uinv);



      v[k] = (resx[nres-1]-resx[0])/365.25/2.0*(pow(param[0],2)+pow(param[1],2))/pow(365.25*86400.0,2); 
      specX[k] = GLOBAL_OMEGA/2.0/M_PI;
      specY[k] = v[k];
      if(specE!=NULL){
	      specE[k] = (resx[nres-1]-resx[0])/365.25/2.0*(pow(error[0],2)+pow(error[1],2))/pow(365.25*86400.0,2);
      }

    }

  for (i=0;i<nfit;i++)
	  free(cvm[i]);
  free(cvm);
  if(verbose_calc_spectra){
	  printf("\b\b\b\b\b\b\b\b100.0%%\n");
  }
  free(psr);
  printf("Leaving calcSpectra\n");
  return nfit;
}


int calcSpectra(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY,int nfit) {
	return calcSpectraErr(uinv,resx,resy,nres,specX,specY,NULL,nfit);
}

// Spectral analysis using covariance matrix
// note: uinv array must start from 0, not 1
// NEW FEATURE:
// set nfit < 0 to automatically set it to nres/2-1
int calcSpectra_ri(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY_R,double *specY_I,int nfit,pulsar* psr) {

  //  pulsar *psr=NULL;
  return calcSpectra_ri_T(uinv,resx,resy,nres,specX,specY_R,specY_I,nfit,(resx[nres-1]-resx[0]),'N',psr);
}

int calcSpectra_ri_T(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY_R,double *specY_I,int nfit,double T,char fitfuncMode, pulsar* psr) {
  int i,j,k;
  int nSpec;

  if (nfit < 0)
    nfit=nres/2-1;

  double v[nfit];
  double chisq;
  int ip[nres];
  double param[3];

  // to allow for computing the spectrum of IFUNC/CM etc we need different fit functions.
  void (*FIT_FUNC)(double, double [], int,pulsar *,int,int); // the fit function we will use
  FIT_FUNC=fitMeanSineFunc; // this is the standard periodogram/
  if(fitfuncMode=='I'){
	 if (psr==NULL){
		logerr("ERROR: cannot use IFUNC spectral estimation without a real pulsar\nCheck call to calcSpectra_ri_T()");
		return -1;
	 }
	 FIT_FUNC=fitMeanSineFunc_IFUNC; // this accounts for smoothing of the CM/IFUNC
  }
  double binfactor = (double)nres/(double)(nres-1);
  if(fitfuncMode=='T'){
	 binfactor=1.0;
  }

  // Should fit independently to all frequencies
  for (i=0;i<nres;i++)
  {
	 ip[i] = i;
  }
  logmsg("Computing %d spectral channels",nfit);
  for (k=0;k<nfit;k++)
  {
	 GLOBAL_OMEGA = 2.0*M_PI/(T*binfactor)*(k+1);
	 printf("In here k = %d\n",k);
	 TKleastSquares_single_pulsar(resx,resy,nres,param,NULL,3,NULL,&chisq,FIT_FUNC,psr,1.0e-40,ip,1,uinv);
	 printf("Out here k = %d\n",k);
	 v[k] = (resx[nres-1]-resx[0])/365.25/2.0/pow(365.25*86400.0,2); 
	 specX[k] = GLOBAL_OMEGA/2.0/M_PI;
	 specY_R[k] = sqrt(v[k])*param[1];
	 specY_I[k] = sqrt(v[k])*param[2];
  }
  return nfit;
}

// Fit for mean and sine and cosine terms at a specified frequency (G_OMEGA)
// The psr and ival parameters are ignored
void fitMeanSineFunc(double x,double *v,int nfit,pulsar *psr,int ival,int ipsr)
{
   int i;
   v[0] = 1; // Fit for mean
   v[1] = cos(GLOBAL_OMEGA*x);
   v[2] = sin(GLOBAL_OMEGA*x);

}

/*
 * This only works for residuals faked from an IFUNC.
 * Note that ival must be in order.
 */
void fitMeanSineFunc_IFUNC(double x,double *v,int nfit,pulsar *psr,int ival,int ipsr)
{
   double m,c; // for the straight line.
   double x0,x1;
   double i0,i1;
   double w = GLOBAL_OMEGA;
   double A=0;

   //   logmsg("%lf %lf %d %d",x,w,ival,psr->nobs);
   v[0] = 0; // Fit for mean
   v[1]=0; // will be cos
   v[2]=0; // will be sin

   x0 = x;
   if (ival< psr->nobs-1){
	  x1 = x0 + (psr->obsn[ival+1].sat-psr->obsn[ival].sat); // this is the next x-val, under the assumptions.
	  m=-1/pow(x0-x1,2);
	  c=x1/pow(x0-x1,2);

	  //logmsg("R) x0=%lf x1=%lf m=%lg c=%lg",x0,x1,m,c);
	  i0=m*x0*x0/2.0 + c*x0;
	  i1=m*x1*x1/2.0 + c*x1;
	  v[0]+=i1-i0;

	  // integrate the right side of the triangle.
	  i0=(m*cos(w*x0) - w*(c+m*x0)*sin(w*x0))/w/w; // integral of (m*x+c).cos(w*x)
	  i1=(m*cos(w*x1) - w*(c+m*x1)*sin(w*x1))/w/w; // integral of (m*x+c).sin(w*x)
	  v[1]+=i1-i0;

	  i0=(m*sin(w*x0) - w*(c+m*x0)*cos(w*x0))/w/w;
	  i1=(m*sin(w*x1) - w*(c+m*x1)*cos(w*x1))/w/w;
	  v[2]+=i1-i0;
	  A+=0.5; // area
   }
   if (ival > 0){
	  x1 = x0 + (psr->obsn[ival-1].sat-psr->obsn[ival].sat); // this is the next x-val, under the assumptions.
	  m=1/pow(x0-x1,2);
	  c=-x1/pow(x0-x1,2);

	  //	  logmsg("L) x0=%lf x1=%lf m=%lg c=%lg",x0,x1,m,c);
	  i0=m*x0*x0/2.0 + c*x0;
	  i1=m*x1*x1/2.0 + c*x1;
	  v[0]+=i0-i1;

	  // integrate the right side of the triangle.
	  i0=(m*cos(w*x0) - w*(c+m*x0)*sin(w*x0))/w/w;
	  i1=(m*cos(w*x1) - w*(c+m*x1)*sin(w*x1))/w/w;
	  v[1]+=i0-i1;

	  i0=(m*sin(w*x0) - w*(c+m*x0)*cos(w*x0))/w/w;
	  i1=(m*sin(w*x1) - w*(c+m*x1)*cos(w*x1))/w/w;
	  v[2]+=i0-i1;
	  A+=0.5; // area
   }
   v[0]/=A; // normalise to 1 in case we only had half the triangle.
   v[1]/=A; // normalise to 1 in case we only had half the triangle.
   v[2]/=A;

   //  double old[3];
   //  fitMeanSineFunc(x, old,nfit,psr,ival);
   //   logmsg("%lg %lg --  %lg %lg -- %lg %lg\n",v[0],old[0],v[1],old[1],v[2],old[2]);
}


// Fit for mean and sine and cosine terms at a specified frequency (G_OMEGA)
// The psr and ival parameters are ignored
void fitCosSineFunc(double x,double *v,int nfit,pulsar *psr,int ival,int ipsr)
{
   int i;
   v[0] = cos(GLOBAL_OMEGA*x);
   v[1] = sin(GLOBAL_OMEGA*x);

}
