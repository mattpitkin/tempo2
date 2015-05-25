//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russel Edwards

#ifndef __TKspectrum_h
#define __TKspectrum_h
#include "tempo2.h"
#include "TKfit.h"
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

extern double GLOBAL_OMEGA;
extern bool verbose_calc_spectra;
/* Global variables used in the fortran code for sigmaz*/

void readin(pulsar psr);
void getprtj(int n);
void indexx8(int n,double *arrin,int *indx);
void getweights(int n, double *wt);
void fit4(int *nfit,double *p4,double *cov4,int ndostats,double *chidf,double *avewt);
void mat20(double sam[21][21],double a[21][21],int n,double *determ,int *nbad);

void sineFunc(double x,double *v,int ma);
void TKsortit(double *x,double *y,int n);
void TKaveragePts(double *x,double *y,int n,int width,double *meanX,double *meanY,int *nMean);
void TKcmonot (int n, double x[], double y[], double yd[][4]);
void TKspline_interpolate(int n,double *x,double *y,double yd[][4],double *interpX,
			  double *interpY,int nInterp);
void TKinterpolateSplineSmoothFixedXPts(double *inX, double *inY, int inN, double *interpX, double *interpY, int nInterp);
void TKhann(double *x,double *y,int n,double *ox,double *oy,int *on,int width);
void TKfirstDifference(double *x,double *y,int n);
void TK_fitSine(double *x,double *y,double *e,int n,int wErr,double *outX,double *outY,int *outN);

void TKlomb_d(double *x,double *y,int n,double ofac,double hifac,double *ox,double *oy,int *outN,double *var);
int TK_fft(short int dir,long n,double *x,double *y);
void TK_dft(double *x,double *y,int n,double *outX,double *outY,int *outN, double *outY_re, double *outY_im);
void TK_weightLS(double *x,double *y,double *sig,int n,double *outX,double *outY,int *outN, double *outY_re, double *outY_im);
void TK_fitSinusoids(double *x,double *y,double *sig,int n,double *outX,double *outY,int *outN);
void fitMeanSineFunc(double x,double *v,int nfit,pulsar *psr,int ival,int ipsr);
void fitCosSineFunc(double x,double *v,int nfit,pulsar *psr,int ival,int ipsr);

int calcSpectraErr(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY,double* specE,int nfit);


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
		double *outY,int *nout,int calcWhite,int output, double *outY_re, double *outY_im);

void TKfirstDifference(double *x,double *y,int n);


/* TK_fitSinusoids attempts to fit a bunch of harmonically related sinusoids to the data. */

void TK_fitSinusoids(double *x,double *y,double *sig,int n,double *outX,double *outY,int *outN);
void sineFunc(double x,double *v,int ma);

//GEORGE'S ALGORITHM!!!! I think mine is better.     Calculates a weighted Lomb-Scargle periodogram
//void TK_weightLS(double *x,double *y,double *sig,int n,double *outX,double *outY,int *outN)
void TK_weightLS(double *x,double *y,double *sig,int n,double *outX,double *outY,int *outN, double *outY_re, double *outY_im);


/* Calculates the discrete Fourier transform of a data-set */
/* Assumes that the x-values are in units of days and the  */
/* y-values are in seconds                                 */
/* This code was checked by G. Hobbs - 30/10/08            */
/* DO NOT MODIFY WITHOUT LEAVING COMMENTS HERE             */
/* 13th November: DY added capability for retrieval of imaginary part into main code */
void TK_dft(double *x,double *y,int n,double *outX,double *outY,int *outN, double *outY_re, double *outY_im);


void TKaveragePts(double *x,double *y,int n,int width,double *meanX,double *meanY,int *nMean);
/* This sorting routine should be replaced by something better in toolkit.h  */
void TKsortit(double *x,double *y,int n);

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
void TKboxcar(double *x,double *y,int n,double *ox,double *oy,int *on,int width);

/* 
 * Hann (or raised cosine) smoother
 */
void TKhann(double *x,double *y,int n,double *ox,double *oy,int *on,int width);

/* ************************************************************** *
 * Spline routines                                                *
 * ************************************************************** */
#define ABS(x) ((x) < 0 ? -(x) : (x))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))


void TKcmonot (int n, double x[], double y[], double yd[][4]);
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
// Interpolate the spline fit on to a given set of x-values
void TKspline_interpolate(int n,double *x,double *y,double yd[][4],double *interpX,
			  double *interpY,int nInterp);

void TKlomb_d(double *x,double *y,int n,double ofac,double hifac,double *ox,double *oy,int *outN,double *var);
// Fourier transform
// Requires npts = 2**M value - does not check for this
// Routine modified from http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/dft
/*
   This computes an in-place complex-to-complex FFT 
   x and y are the real and imaginary arrays of n = 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform 
*/
int TK_fft(short int dir,long n,double *x,double *y);


// Routines for sigmaz calculation
/* Direct copy of the original Fortran */
/* using identical fitting routines etc. */
/* Note that all arrays start from 1 */
//void calcSigmaz(pulsar *psr,int npsr,int time,int psrfile,double mintau,
//		int nSel,int *selPsr,int weights,int style,int average,
//		float *avX,float *avY,int *nav)
void TKcalcSigmaz(pulsar psr,int weights,double *ret_tau,double *ret_szbias,double *ret_e1,
		double *ret_e2,int *ret_nval,double mintau);



/* WHAT AM I RETURNING IN AND OUT OF THIS FUNCTION? */
void fit4(int *nfit,double *p4,double *cov4,int ndostats,double *chidf,double *avewt);
void readin(pulsar psr);

// Check conversion from Fortran carefully
void indexx8(int n,double *arrin,int *indx);
//Adapted from Stefan / George's plugin
//interpolation (spline): this function interpolates a data set using constrained spline onto an input set of interpX and nInterp values
void TKinterpolateSplineSmoothFixedXPts(double *inX, double *inY, int inN, double *interpX, double *interpY, int nInterp);


int calcSpectra(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY,int nfit);

// Spectral analysis using covariance matrix
// note: uinv array must start from 0, not 1
// NEW FEATURE:
// set nfit < 0 to automatically set it to nres/2-1
int calcSpectra_ri(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY_R,double *specY_I,int nfit,pulsar *psr);
int calcSpectra_ri_T(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY_R,double *specY_I,int nfit,double T,char useCM,pulsar* psr);


// Fit for mean and sine and cosine terms at a specified frequency (G_OMEGA)
// The psr and ival parameters are ignored
void fitMeanSineFunc(double x,double *v,int nfit,pulsar *psr,int ival,int ipsr);
void fitMeanSineFunc_IFUNC(double x,double *v,int nfit,pulsar *psr,int ival,int ipsr);

// Fit for mean and sine and cosine terms at a specified frequency (G_OMEGA)
// The psr and ival parameters are ignored
void fitCosSineFunc(double x,double *v,int nfit,pulsar *psr,int ival);

#endif
