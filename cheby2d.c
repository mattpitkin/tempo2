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

// redwards code for finding and evaluating 2D Chebyshev polynomials.
// 
// Note, 2D (Cartesian) Chebyshev polynomials actually follow quite trivially
// from the 1D case.
//
// Defining 
//
//   Tnm(x) = Tn(x)Tm(y)  , where Tn(x) is the Chebyshev polnomial of deg n,
//
//  the discrete orthogonality condition,
//
//    sum_a=1^n sum_b=1^m Tij(xa,yb)Tkl(xa,yb) = 0
//
//  is satisfied by the standard 1D choice of xa, yb, i.e. the zeros of
//  the 1D basis function, i.e. x = cos(pi[k-1/2]/n) for Tn . Then,
//  rewriting we have
//
//  sum_a=1^n Tij(xa,yb) sum_b=1^m Tkl(xa,yb) = 0 
//
//  which implies  ! (j==l && i==k)
//
//  i.e. the bases are orthonormal
//
// this means we can define
//
// c_ij = sum_a=1^n sum_b=1^m f(xa,yb)Tij(xa,yb)
//
// and then have f(x,y) ~= sum_i=1^n sum_j=1^m  c_ij Tij(x,y)
//
//  where the relationship is exact for x=xa , y=ya and pretty good elsewhere

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef sun
#include <sunmath.h>
#endif

#include "tempo2pred.h"
#include "tempo2pred_int.h"

#ifndef M_PIl
#define M_PIl 3.14159265358979323846264338327950288L
#endif



void Cheby2D_Init(Cheby2D *cheby, int nx, int ny)
{
  cheby->nx = nx;
  cheby->ny = ny;
  if (nx && ny)
    cheby->coeff = (long double*)malloc(nx*ny*sizeof(long double));
  else
    cheby->coeff = 0;
}

void Cheby2D_Destroy(Cheby2D *cheby)
{
  free(cheby->coeff);
}

void Cheby2D_Copy(Cheby2D *cheby, const Cheby2D* from)
{
  Cheby2D_Destroy (cheby);
  Cheby2D_Init (cheby, from->nx, from->ny);
  memcpy (cheby->coeff, from->coeff, from->nx*from->ny*sizeof(long double));
}

void
Cheby2D_Construct(Cheby2D *cheby,
		  void
		  (*func)(long double *x, long double *y, 
			  int nx, int ny, long double *z, void *info),
		  void *info)
{
  int i, j, a, b;
  int nx = cheby->nx;
  int ny = cheby->ny;

  // first of all pre-compute all the x's and y's
  long double *x = (long double *)malloc(nx*sizeof(long double));
  long double *y = (long double *)malloc(ny*sizeof(long double));

  long double *f = (long double*)malloc(nx*ny*sizeof(long double));
  long double *fcurr = f;

  long double *Tx = (long double*)malloc(nx*nx*sizeof(long double));
  long double *Txcurr = Tx;

  long double *Ty = (long double*)malloc(ny*ny*sizeof(long double));
  long double *Tycurr = Ty;
  for (i=0; i < nx; i++)
    x[i] = cosl(M_PIl * (i+0.5L) / nx);
  for (i=0; i < ny; i++)
    y[i] = cosl(M_PIl * (i+0.5L) / ny);
  // and pre-compute the function on the resultant grid
  func(x, y, nx, ny, f, info);
//   for (b=0; b < ny; b++)
//     for (a=0; a < nx; a++)
//       *(fcurr++) = func(x[a], y[b], info);

  // and, precompute the x and y basis functions
  for (i=0 ; i < nx; i++)
    for (a=0; a < nx; a++)
      *(Txcurr++) = cosl(M_PIl * (i+0) * (a+0.5L) / nx); //=Ti(xa)

  for (j=0 ; j < ny; j++)
    for (b=0; b < ny; b++)
      *(Tycurr++) = cosl(M_PIl * (j+0) * (b+0.5L) / ny);//=Tj(yb)

  // now compute the inner product of the function and the basis, for
  // each basis.
  // for efficiency we can write 
  // c_ij = sum_b=1^m Tj(yb) sum_a=1^n Ti(xa) f(xa,yb)

  for (i=0 ; i < nx; i++)
  {
    Tycurr = Ty;
    for (j=0; j < ny; j++)
    {
      long double sum = 0.0L;
      fcurr = f;
      for (b=0; b < ny; b++) // for each row
      {
	long double sum2 = 0.0L;
	Txcurr = Tx+i*nx;
	for (a=0; a < nx; a++) // get inner prod of f with x-basis
	  sum2 += *(fcurr++) * *(Txcurr++);  //=f(xa,xb)*Ti(xa)
	sum += *(Tycurr++) * sum2; // = Tj(yb)*sum_a[f(xa,xb)*Ti(xa)]
      }
      cheby->coeff[j*nx+i] = 4.0*sum/(nx*ny); // store result, c_ij = sum_a=1^n sum_b=1^m f(xa,yb)Tij(xa,yb)
    }
  }
  free(Ty);
  free(Tx);
  free(f);
  free(y);
  free(x);
}


// Here is a function to construct a function that's the partial derivative
// wrt x of the a second chebyshev polynomial
void
Cheby2D_Construct_x_Derivative(Cheby2D *dcheby, const Cheby2D *cheby)
{
  int ix, iy, i=0, nx=cheby->nx;

  if (dcheby->nx != cheby->nx || dcheby->ny != cheby->ny)
  {
    // if this happens, the derivative cheby hasn't been initialized properly
    fprintf(stderr, "Programming error, %s:%d\n", __FILE__, __LINE__);
    exit(1);
  }

  for (iy=0; iy < cheby->ny; iy++)
  {
    dcheby->coeff[i+nx-1] = 0.0;
    dcheby->coeff[i+nx-2] = 2*(nx-1)*cheby->coeff[i+nx-1];    
    for (ix=nx-3; ix >= 0; ix--)
      dcheby->coeff[i+ix]= dcheby->coeff[i+ix+2]+2*(ix+1)*cheby->coeff[i+ix+1];
    i+=nx;
  }
}

// here is the function to evaulate a 2D Chebyshev polynomial,
// f(x,y) ~= sum_i=1^n sum_j=1^m  c_ij Tij(x,y)
//         = sum_j=1^m Tj(yb) sum_i=1^n c_ij Ti(xa)
//  This can be tackled in much the same way as the 1D case:
//  the inner loop is a simple 1D Chebyshev with c_j as the coefficients,
// while the outer loop is a 1D Cheby with the result of the inner as the
// coefficients.
// These are evaluated using Clenshaw's recurrence, as in Numerical Recipes,
//   f(x) = xd_1 - d_2 + c_0/2
//   d_j = 2xd_(j+1) - d_(j+2) + c_j
//   d_(m+1) = d_m = 0

long double 
Cheby2D_Evaluate(const Cheby2D *cheby, long double x, long double y)
{
  long double di, di1, di2; // inner loop values d_i, d_(i+1), d_(i+2)
  long double dj, dj1, dj2; // as above, for outer loop
  int i, j;
  long double *coeffcurr = cheby->coeff;

  dj = dj1 = 0.0L;
  for (j=cheby->ny-1; j >= 0; j--)
  {
    di = di1 = 0.0L;
    coeffcurr = cheby->coeff+(j+1)*cheby->nx-1;
    for (i=cheby->nx-1; i > 0; i--)
    {
      di2 = di1;
      di1 = di;
      di = 2.0L*x*di1 - di2 + *(coeffcurr--); // i.e. + c_ij
    }
    // compute di0, special case
    di = x*di - di1 + 0.5L*(*coeffcurr--);
    dj2 = dj1;
    dj1 = dj;
    if (j > 0)
      dj = 2.0L*y*dj1 - dj2 + di;
    else // do special case inside loop since we still need to know di0 (di)
      dj = y*dj1 - dj2 + 0.5L*di;
  }
  return dj;
}


void
Cheby2D_Test(Cheby2D *cheby, int nx_test, int ny_test,
	     void
	     (*func)(long double *x, long double *y, 
		     int nx, int ny, long double *z, void *info),
	     void *info,
	     long double *residualRMS, long double *residualMAV)
{
  int ix_test, iy_test;
  long double  fprime;
  long double sdiff, ssdiff, maxdiff, diff;
  long double *f = (long double*)malloc(nx_test*ny_test*sizeof(long double));
  long double *x = (long double *)malloc(nx_test*sizeof(long double));
  long double *y = (long double *)malloc(ny_test*sizeof(long double));

  sdiff=ssdiff=maxdiff=0.0L;

  // call func to get actual values on grid
  for (iy_test=0; iy_test<ny_test; iy_test++)
    y[iy_test] = -1.0L + 2.0L*(iy_test+0.5L)/ny_test; 
  for (ix_test=0; ix_test<nx_test; ix_test++)
    x[ix_test] = -1.0L + 2.0L*(ix_test+0.5L)/nx_test; 
  func(x, y, nx_test, ny_test, f, info);

  // compare
  for (iy_test=0; iy_test<ny_test; iy_test++)
  {
    for (ix_test=0; ix_test<nx_test; ix_test++)
    {
      fprime = Cheby2D_Evaluate(cheby, x[ix_test], y[iy_test]);
      diff = f[iy_test*nx_test+ix_test]-fprime;

#if 0
      {
      // redwards hack to simulate 1-D polyco.. first get prediction from
      // centre freq
      long double phase1d = Cheby2D_Evaluate(cheby, x[ix_test], 0.0);
      // second, get instantaneous freq at band center
      long double dx = 1.0e-5L;
      long double freq1 = 
	(Cheby2D_Evaluate(cheby, x[ix_test]+dx*0.5L, y[iy_test])-
	 Cheby2D_Evaluate(cheby, x[ix_test]-dx*0.5L, y[iy_test])) / dx;
      long double freq=
	(Cheby2D_Evaluate(cheby, x[ix_test]+dx*0.5L, 0.0L)-
      	 Cheby2D_Evaluate(cheby, x[ix_test]-dx*0.5L, 0.0L)) / dx;
      long double skyfreq = 653.0L+y[iy_test]*32.0L; // MHzx
      long double dmdelay = 1.0L/2.41e-4L *(1.0L/ (skyfreq*skyfreq) -1.0L/(653.0L*653.0L))
	*48.901787L; // s
      long double dphase = -freq*dmdelay;
      long double diff1d = f[iy_test*nx_test+ix_test]-(phase1d+dphase);

      freq *= 2.0L/0.03L / 86400.0L; // normalized->day^-1->Hz
      freq1 *= 2.0L/0.03L / 86400.0L;

//             printf("YYY %Lf %Lf %Lf\n", freq, dmdelay, dphase);
//       printf("%Lf %Lf %Lg %Lg %Lg %Lg %Lg XXX\n", x[ix_test], y[iy_test], 
//   	     f[iy_test*nx_test+ix_test], fprime, diff, diff1d, 
// 	     (freq1-freq)/freq);
      }
#endif
      if (fabs(diff) > maxdiff)
	maxdiff = fabs(diff);
      sdiff += diff;
      ssdiff += diff*diff;
    }
    //      printf("\n");
  }
  sdiff /=nx_test*ny_test;
  *residualRMS = sqrtl(ssdiff/(nx_test*ny_test));
  *residualMAV = maxdiff;

  free(f);
  free(x);
  free(y);
}

void testFunc(long double *x, long double *y, 
	      int nx, int ny, long double *z, void *info)
{
  int ix, iy;
  for (iy=0; iy < ny; iy++)
    for (ix=0; ix < nx; ix++)
      *(z++) =  0.5*x[ix] + 0.1*sinl(x[ix]) + 0.7*y[iy]+0.1*y[iy]*y[iy];
}

void testCheby2D()
{
  int nx=30, ny=30;
  long double rms, mav;
  Cheby2D cheby;
  Cheby2D_Init(&cheby, nx, ny);
  printf("Constructing..."); fflush(stdout);
  Cheby2D_Construct(&cheby, testFunc, NULL);
  printf("\nTesting..."); fflush(stdout);
  Cheby2D_Test(&cheby, nx*3, ny*3,testFunc, NULL, &rms, &mav);
  printf("\nRMS= %Lg MAV= %Lg\n", rms, mav);
  Cheby2D_Destroy(&cheby);
}


void
ChebyModel_Init(ChebyModel *cm, int nmjdcoeff, int nfreqcoeff)
{
  Cheby2D_Init(&cm->cheby, nmjdcoeff, nfreqcoeff);
  Cheby2D_Init(&cm->frequency_cheby, nmjdcoeff, nfreqcoeff);
}

void
ChebyModel_Copy(ChebyModel *cm, ChebyModel *from)
{
  strcpy (cm->psrname, from->psrname);
  strcpy (cm->sitename, from->sitename);
  cm->mjd_start = from->mjd_start;
  cm->mjd_end = from->mjd_end;
  cm->freq_start = from->freq_start;
  cm->freq_end = from->freq_end;
  cm->dispersion_constant = from->dispersion_constant;

  Cheby2D_Copy(&cm->cheby, &from->cheby);
  Cheby2D_Copy(&cm->frequency_cheby, &from->frequency_cheby);
}

void
ChebyModel_Destroy(ChebyModel *cm)
{
  Cheby2D_Destroy(&cm->cheby);
  Cheby2D_Destroy(&cm->frequency_cheby);
}

long double
ChebyModel_GetPhase(const ChebyModel *cm, long double mjd, long double freq)
{
  if (!cm)
    return -1;

  return Cheby2D_Evaluate
    (&cm->cheby, 
     -1.0L+2.0L*(mjd-cm->mjd_start)/(cm->mjd_end-cm->mjd_start),
     -1.0L+2.0L*(freq-cm->freq_start)/(cm->freq_end-cm->freq_start))
    +      cm->dispersion_constant / (freq*freq);
}

long double
ChebyModel_GetFrequency(const ChebyModel *cm, long double mjd, long double freq)
{
  if (!cm)
    return -1;

  return Cheby2D_Evaluate
    (&cm->frequency_cheby, 
     -1.0L+2.0L*(mjd-cm->mjd_start)/(cm->mjd_end-cm->mjd_start),
     -1.0L+2.0L*(freq-cm->freq_start)/(cm->freq_end-cm->freq_start))
    // this gives cycles per half of the whole MJD interval.. scale to per day
    * 2 / (cm->mjd_end-cm->mjd_start)
    // then scale to per second
    / 86400.0;
}



void ChebyModel_Write(const ChebyModel *cm, FILE *f)
{
  int ix, iy;
  fprintf(f, "ChebyModel BEGIN\n");
  fprintf(f, "PSRNAME %s\n", cm->psrname);
  fprintf(f, "SITENAME %s\n", cm->sitename);
  fprintf(f, "TIME_RANGE %.34Lg %.34Lg\n", cm->mjd_start, cm->mjd_end);
  fprintf(f, "FREQ_RANGE %.34Lg %.34Lg\n", cm->freq_start, cm->freq_end);
  fprintf(f, "DISPERSION_CONSTANT %.34Lg\n", cm->dispersion_constant);
  fprintf(f, "NCOEFF_TIME %d\n", cm->cheby.nx);
  fprintf(f, "NCOEFF_FREQ %d\n", cm->cheby.ny);

  for (ix=0; ix < cm->cheby.nx; ix++)
  {
    fprintf(f, "COEFFS");
    for (iy=0; iy < cm->cheby.ny; iy++)
      {
	//	fprintf(f, " %.34Lg", cm->cheby.coeff[iy*cm->cheby.nx+ix]);
	fprintf(f, " %.25Lg", cm->cheby.coeff[iy*cm->cheby.nx+ix]);
	if ((iy+1)%3==0) fprintf(f, "\n");  // Every 3 coefficients put a new line
      }
    fprintf(f, "\n");
  }
  fprintf(f, "ChebyModel END\n");
}

int ChebyModel_Read(ChebyModel *cm, FILE *f)
{
  int done = 0;
  int first = 1;
  char line[1024], keyword[64], arg[64], junk[1024];
  int nx=-1, ny=-1, ix=0, iy;
  int ichar, nread;

  cm->cheby.coeff=NULL;

  do
  {
    if (fgets(line, 1024, f)!=line)
      return -1;
    if (sscanf(line, "%s", keyword)!=1)
      continue; // skip blank lines
    if (sscanf(line, "%s %s", keyword, arg)!=2)
      return -2;
    if (line[0]=='#')
      continue; // skip comment lines
    // check first line
    if (first && (strcasecmp(keyword, "ChebyModel")||strcasecmp(arg, "BEGIN")))
      return -3;
    // parse based on keyword
    if (!strcasecmp(keyword, "PSRNAME"))
      strcpy(cm->psrname, arg);
    else if (!strcasecmp(keyword, "SITENAME"))
      strcpy(cm->sitename, arg);
    else if (!strcasecmp(keyword, "TIME_RANGE"))
    {
      if (sscanf(line, "%*s %Lf %Lf", &cm->mjd_start, &cm->mjd_end)!=2)
	return -4;
    }
    else if (!strcasecmp(keyword, "FREQ_RANGE"))
    {
      if (sscanf(line, "%*s %Lf %Lf", &cm->freq_start, &cm->freq_end)!=2)
	return -5;
    }
    else if (!strcasecmp(keyword, "DISPERSION_CONSTANT"))
    {
      if (sscanf(arg, "%Lf", &cm->dispersion_constant)!=1)
	return -6;
    }
    else if (!strcasecmp(keyword, "NCOEFF_TIME"))
    {
      if (sscanf(arg, "%d", &nx)!=1) 
	return -7;
    }
    else if (!strcasecmp(keyword, "NCOEFF_FREQ"))
    {
      if (sscanf(arg, "%d", &ny)!=1)
	return -8;
    }
    else if (!strcasecmp(keyword, "COEFFS"))
    {
      if (cm->cheby.coeff==NULL) // first instance of COEFF keyword
      {
	if (nx < 0 && ny < 0) // oops, these should come first!
	  return -8;
	ChebyModel_Init(cm, nx, ny);
      }
      if (ix >= nx)
	return -9; // too many coefficient lines!!

      sscanf(line, "%*s %n", &ichar);
      if (ny<4) // All on one line
	{
	  for (iy=0; iy < cm->cheby.ny; iy++)
	    {
	      if (sscanf(line+ichar, "%Lf %n", 
			 &cm->cheby.coeff[iy*cm->cheby.nx+ix], &nread)!=1)
		return -10;
	      ichar += nread;
	    }
	}
      else   // Code added by G. Hobbs for multiple lines in the predictor file
	{
	  for (iy=0; iy < cm->cheby.ny; iy++)
	    {
	      if (sscanf(line+ichar, "%Lf %n", 
			 &cm->cheby.coeff[iy*cm->cheby.nx+ix], &nread)!=1)
		return -10;
	      ichar += nread;
	      if ((iy+1)%3==0)
		{
		  if (sscanf(line+ichar, "%s", junk)==1)
		    return -11; // excess stuff at end of line		  
		  ichar = 0;
		  if (fgets(line, 1024, f)!=line)
		    return -1;
		}
	    }
	}
      if (sscanf(line+ichar, "%s", junk)==1)
	return -11; // excess stuff at end of line
      ix++;
    }
    else if (!strcasecmp(keyword, "ChebyModel"))
    {
      if ((!first) && !strcasecmp(arg, "BEGIN"))
	return -12;
      else if (!strcasecmp(arg, "END"))
      {
	if (cm->cheby.coeff==NULL || ix!=nx)
	  return -13; // haven't read enough coefficients yet!!
	else
	{
	  Cheby2D_Construct_x_Derivative(&cm->frequency_cheby, &cm->cheby);

	  done = 1;
	}
      }
	       
    }
    else
      return -14; // unrecognized keyword!! 
    first = 0;    
  } while (!done);

  return 0;
}


int ChebyModelSet_GetNearestIndex(const ChebyModelSet *cms, long double mjd)
{
  int inearest=-1;

  long double best_offset=1e6, offset;
  int iseg;

  ChebyModelSet_OutOfRange = 0;

  for (iseg=0; iseg < cms->nsegments ; iseg++)
  {
    if (mjd < cms->segments[iseg].mjd_start)
      continue;
    if (mjd > cms->segments[iseg].mjd_end)
      continue;

    offset = fabs((cms->segments[iseg].mjd_start+cms->segments[iseg].mjd_end)/2
		  - mjd);
    if (offset < best_offset)
    {
      inearest = iseg;
      best_offset = offset;
    }
  }

  return inearest;
}


int ChebyModelSet_OutOfRange = 0;


ChebyModel *ChebyModelSet_GetNearest(const ChebyModelSet *cms, long double mjd)
{
  int inearest = ChebyModelSet_GetNearestIndex(cms,mjd);

  if (inearest < 0) {
    ChebyModelSet_OutOfRange = 1;
    return 0;
  }

  return &cms->segments[inearest];
}


long double ChebyModelSet_GetPhase(const ChebyModelSet *cms, long double mjd, long double freq)
{
  return ChebyModel_GetPhase(ChebyModelSet_GetNearest(cms, mjd), mjd, freq);
}
long double ChebyModelSet_GetFrequency(const ChebyModelSet *cms, long double mjd, long double freq)
{
  return ChebyModel_GetFrequency(ChebyModelSet_GetNearest(cms, mjd), mjd, freq);
}

 

void ChebyModelSet_Write(const ChebyModelSet *cms, FILE *f)
{
  int iseg;

  fprintf(f, "ChebyModelSet %d segments\n", cms->nsegments);
  for (iseg=0; iseg < cms->nsegments ; iseg++)
    ChebyModel_Write(&cms->segments[iseg], f);
}

int ChebyModelSet_Read(ChebyModelSet *cms, FILE *f)
{
  char line[1024], keyword[64];
  int iseg;
  int ret;
  if (fgets(line, 1024, f)!=line)
    return -1;
  if (sscanf(line, "%s %d", keyword, &cms->nsegments)!=2)
    return -1;
  if (strcasecmp(keyword, "ChebyModelSet"))
    return -1;
  cms->segments = (ChebyModel *)malloc(cms->nsegments*sizeof(ChebyModel));

  for (iseg=0; iseg < cms->nsegments ; iseg++)
    if ((ret=ChebyModel_Read(&cms->segments[iseg], f)) != 0)
      return ret;
  return 0;
}

void ChebyModelSet_Init(ChebyModelSet *cms)
{
  cms->segments = 0;
  cms->nsegments = 0;
}

int ChebyModelSet_Insert(ChebyModelSet *cms, const ChebyModelSet *from)
{
  int old_nseg = cms->nsegments;
  int iseg;

  cms->nsegments += from->nsegments;
  cms->segments = (ChebyModel *) realloc (cms->segments, 
					  cms->nsegments*sizeof(ChebyModel));

  for (iseg=old_nseg; iseg < cms->nsegments ; iseg++)
  {
    ChebyModel_Init(&cms->segments[iseg], 0, 0);
    ChebyModel_Copy(&cms->segments[iseg], &from->segments[iseg-old_nseg]);
  }
}

/*
  This method destroys all ChebyModel elements that are no longer required
  The MJD array defines what is required; only the nearest ChebyModel
  to each MJD is kept.
*/ 
void
ChebyModelSet_Keep(ChebyModelSet *cms, unsigned nmjd, const long double* mjd)
{
  unsigned nseg = cms->nsegments;  // number of segments in input
  unsigned new_nseg = nseg;        // number of segments kept

  unsigned iseg = 0;               // current segment index
  unsigned rem_iseg = 0;           // remaining segment index
  
  unsigned i = 0;                  // counter
  char* keep = malloc (nseg);      // array of "ChebyModel to be kept" flags

  memset (keep, 0, nseg);          // set all keep flags to false

  for (i=0; i<nmjd; i++)
  {
    int index = ChebyModelSet_GetNearestIndex (cms, mjd[i]);
    if (index >= 0)
      keep[index] = 1;         // set flag to true
  }

  for (i=0; i<nseg; i++)
    if ( keep[i] )
      iseg ++;
    else
      {
	ChebyModel_Destroy(&cms->segments[iseg]);
	// shift the remaining ChebyModel elements "to the left"
	for (rem_iseg = iseg; rem_iseg < new_nseg-1; rem_iseg++)
	  cms->segments[rem_iseg] = cms->segments[rem_iseg+1];

	new_nseg --;
      }

  cms->nsegments = new_nseg;

  free (keep);
}


void ChebyModelSet_Destroy(ChebyModelSet *cms)
{
  int iseg;
  for (iseg=0; iseg < cms->nsegments ; iseg++)
    ChebyModel_Destroy(&cms->segments[iseg]);

  free(cms->segments);
}
 

#if 0

int
main(int argc, char *argv[])
{
  testCheby2D();

  return 0;
}
#endif
