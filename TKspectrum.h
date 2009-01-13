/* Routines useful for spectral analysis 
 * Based on work carried out by G. Hobbs
 * and W. Coles
 */

#include <stdio.h>
#include <math.h>
#include "T2toolkit.h"
#include "TKfit.h"

double GLOBAL_OMEGA = 0;

void sineFunc(double x,double *v,int ma);
void TKsortit(double *x,double *y,int n);
void TKaveragePts(double *x,double *y,int n,int width,double *meanX,double *meanY,int *nMean);
void TKcmonot (int n, double x[], double y[], double yd[][4]);
void TKspline_interpolate(int n,double *x,double *y,double yd[][4],double *interpX,
			  double *interpY,int nInterp);
void TKhann(double *x,double *y,int n,double *ox,double *oy,int *on,int width);
void TKfirstDifference(double *x,double *y,int n);
void TK_fitSine(double *x,double *y,double *e,int n,int wErr,double *outX,double *outY,int *outN);

void TKlomb_d(double *x,double *y,int n,double ofac,double hifac,double *ox,double *oy,int *outN,double *var);
int TK_fft(short int dir,long n,double *x,double *y);
void TK_dft(double *x,double *y,int n,double *outX,double *outY,int *outN);
void TK_weightLS(double *x,double *y,double *sig,int n,double *outX,double *outY,int *outN);
void TK_fitSinusoids(double *x,double *y,double *sig,int n,double *outX,double *outY,int *outN);

double globalOmega;

typedef struct complexVal {
  double real;
  double imag;
}complexVal;

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
		int fitSpline,int preWhite,int specType,int ofac,int specOut,double *outX,
		double *outY,int *nout,int calcWhite,int output)
{
  double avPtsX[n],avPtsY[n];
  double interpX[MAX_OBSN],interpY[MAX_OBSN];
  double smoothX[MAX_OBSN],smoothY[MAX_OBSN];  
  double specX[MAX_OBSN],specY[MAX_OBSN];
  int splineDaily=0;
  double mean;
  int    i,j;
  int    nAv=0;
  int    nInterp=0;
  int    nSmooth=0;
  int    nSpec=0;

  // Step 1: Sort the data in increasing time order
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
	TK_dft(smoothX,smoothY,nSmooth,specX,specY,&nSpec);
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
	//
	TKlomb_d(smoothX,smoothY,nSmooth,ofac,1,specX,specY,&nSpec,&var);
	for (i=0;i<nSpec;i++)
	  {
	    if (specOut==1) // PSD
	      {
		specY[i]*=(2.0*var*nSpec);
		specY[i] = (specY[i]/pow(365.25*86400.0,2))*2*(tspan/365.25)/(double)nSmooth/(double)nSmooth;
	      }
	    else if (specOut==2) // Amplitude
	      specY[i]=sqrt(specY[i]*2.0*var/nSpec);
	    else if (specOut==3) // Power
	      specY[i]=specY[i]*2.0*var/nSpec;
	  }
	// End of checked section
      }
    else if (specType==3)
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
	TK_weightLS(smoothX,smoothY,e,nSmooth,specX,specY,&nSpec);
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
  *nout = nSpec;  

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

  omega0 = 2.0*M_PI/TKrange_d(x,n);
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

/* Calculates a weighted Lomb-Scargle periodogram */
void TK_weightLS(double *x,double *y,double *sig,int n,double *outX,double *outY,int *outN)
{
  long double s1,s2,s3,s4,s5;
  double recA[MAX_OBSN],recB[MAX_OBSN];
  double pred;
  double omega=0.0;
  double si,ci;
  double a,b;
  double omega0;
  double var;
  int i,j;
  
  omega0 = 2.0*M_PI/TKrange_d(x,n);
  *outN = n/2;


  //    for (i=0;i<n;i++)
  //      sig[i]=1.0;

  for (j=0;j<*outN;j++)
    {
      omega = omega0*(j+1);
      s1 = s2 = s3 = s4 = s5 = 0.0;
      for (i=0;i<n;i++)
	{
	  si = sin(omega*x[i]);
	  ci = cos(omega*x[i]);
	  
	  s1 += (long double)(y[i]*si/sig[i]/sig[i]);
	  s2 += (long double)(si*si/sig[i]/sig[i]);
	  s3 += (long double)(si*ci/sig[i]/sig[i]);
	  s4 += (long double)(y[i]*ci/sig[i]/sig[i]);
	  s5 += (long double)ci*ci/sig[i]/sig[i];
	}
      b = (double)((s4-s1/s2)/(s5-s3/s2));
      a = (double)((s1-b*s3)/s2);
      recA[j] = a;
      recB[j] = b;
      outX[j] = omega/2.0/M_PI;
      outY[j] = a*a+b*b;
    }
  /*  for (i=0;i<1000;i++)
    {
      pred=0.0;
      for (j=0;j<*outN;j++)
	{
	  GLOBAL_OMEGA = omega0*(j+1);
	  pred += (recB[j]*cos(GLOBAL_OMEGA*(x[n-1]-x[0])/1000.0*i+x[0])+recA[j]*sin(GLOBAL_OMEGA*(x[n-1]-x[0])/1000.0*i+x[0]));
	}
      printf("series2: %g %g\n",(x[n-1]-x[0])/1000.0*i+x[0],pred);
      } */

}

/* Calculates the discrete Fourier transform of a data-set */
/* Assumes that the x-values are in units of days and the  */
/* y-values are in seconds                                 */
/* This code was checked by G. Hobbs - 30/10/08            */
/* DO NOT MODIFY WITHOUT LEAVING COMMENTS HERE             */
void TK_dft(double *x,double *y,int n,double *outX,double *outY,int *outN)
{
  complexVal spec[n];
  double xo;
  int nfreq;
  int i,j;
  double f0;
  double range,t,f;

  range = TKrange_d(x,n);
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
    outY[i] = (pow(spec[i].real,2)+pow(spec[i].imag,2));

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
  //  scale = (double)n/(double)(n-1);
  scale = 1;
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

   m = (int)(log(n)/log(2));
   printf("m = %d\n",m);
   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
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
         j -= k;
         k >>= 1;
      }
      j += k;
   }

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

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
   }
   
   return 0;
}
