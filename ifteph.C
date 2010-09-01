#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
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
#include <errno.h>
#include <stdlib.h>
#include <math.h>
#include "ifteph.h"


struct IFTE_interpolation_info
{
   double pc[18],vc[18], twot;
   int np, nv;
};

static void IFTEinterp( struct IFTE_interpolation_info *iinfo,
                 const double coef[], const double t[2], const int ncf,
                 const int ncm, const int na, const int ifl, double posvel[]);

typedef struct
{
  char title[256];
  double startJD, endJD, stepJD;
  int ephver;
  double L_C;
  int swap_endian;
  int reclen;
  int irec;
  double buf[322];
  FILE *f;
  struct IFTE_interpolation_info iinfo;
  int ipt[2][3];
} IFTEphemeris;



static IFTEphemeris ifte;

/* helper functions for endianness */
void
IFTswap4(char *word)
{
  char tmp;
  tmp = word[0]; word[0] = word[3]; word[3] = tmp;
  tmp = word[1]; word[1] = word[2]; word[2] = tmp;
} 

void
IFTswapInt(int *word)
{
  IFTswap4((char *)word);
}

void
IFTswapInts(int *word, int n)
{
  int i;
  for (i=0; i < n; i++)
    IFTswap4((char *)(word+i));
}

void
IFTswap8(char *dword)
{
  char tmp;
  tmp = dword[0]; dword[0] = dword[7]; dword[7] = tmp;
  tmp = dword[1]; dword[1] = dword[6]; dword[6] = tmp;
  tmp = dword[2]; dword[2] = dword[5]; dword[5] = tmp;
  tmp = dword[3]; dword[3] = dword[4]; dword[4] = tmp;
}

void
IFTswapDouble(double *dbl)
{
  IFTswap8((char *)dbl);
}

void
IFTswap8N(char *dwords, int n) /* a bit of code duplication for speed */
{
  char tmp;
  int i;
  for (i=0; i < n; i++)
  {
    tmp = dwords[0]; dwords[0] = dwords[7]; dwords[7] = tmp;
    tmp = dwords[1]; dwords[1] = dwords[6]; dwords[6] = tmp;
    tmp = dwords[2]; dwords[2] = dwords[5]; dwords[5] = tmp;
    tmp = dwords[3]; dwords[3] = dwords[4]; dwords[4] = tmp;
    dwords += 8;
  }
}

void
IFTswapDoubles(double *dbl, int N)
{
  IFTswap8N((char *)dbl, N);
}

void
IFTE_init(const char *fname)
{
  FILE *f;
  static int time=0;
  char buf[1024];
  int ncon;
  double double_in;
  /* open file */
  //  if (time==1 && ifte.f != NULL)
  //    fclose(ifte.f);  

  time=1;
  f = fopen(fname, "r");
  if (!f)
  {
    fprintf(stderr, "Error opening time ephemeris file '%s' : %s\n",
	    fname, strerror(errno));
    exit(1);
  }

  /* read in header info */
  fread(buf, 1, 252, f); /* read CHARACTER*6 TTL(14,3) */
  fread(buf, 1, 12, f); /* read CHARACTER*6 CNAM(2) */
  fread(&ifte.startJD, 1, 8, f);
  fread(&ifte.endJD, 1, 8, f);
  fread(&ifte.stepJD, 1, 8, f);
  fread(&ncon, 1, 4, f);
  
  ifte.swap_endian = (ncon!=2); /* check for endianness */

  if (ifte.swap_endian)  IFTswapInt(&ncon);
  if (ncon!=2) /* check that we can decode the file */
  {
    fprintf(stderr, "Cannot understan format of time ephemeris file '%s'!\n",
	    fname);
    exit(1);
  }
  if (ifte.swap_endian) 
  {
    IFTswapDouble(&ifte.startJD); 
    IFTswapDouble(&ifte.endJD); 
    IFTswapDouble(&ifte.stepJD); 
  }

  fread(ifte.ipt, 8, 3, f);
  if (ifte.swap_endian) IFTswapInts(&ifte.ipt[0][0], 6); 

  /* figure out the record length */
  ifte.reclen = 4 * 2*(ifte.ipt[1][0]-1 + 3*ifte.ipt[1][1]*ifte.ipt[1][2]);
  
  /* get the constants from record "2" */
  fseek(f, ifte.reclen, SEEK_SET);
  fread(&double_in, 8, 1, f);
  if (ifte.swap_endian) IFTswapDouble(&double_in); 
  ifte.ephver = (int)floor(double_in);
  fread(&ifte.L_C, 8, 1, f);
  if (ifte.swap_endian) IFTswapDouble(&ifte.L_C); 

  ifte.f = f;
  ifte.iinfo.np = 2;
  ifte.iinfo.nv = 3;
  ifte.iinfo.pc[0] = 1.0;
  ifte.iinfo.pc[1] = 0.0;
  ifte.iinfo.vc[1] = 1.0;
  ifte.irec = -1;
  /* Note: file is not closed as it is used by other routines */
}
 
void IFTE_close_file()
{
  if (ifte.f != NULL)
    fclose(ifte.f);
}

/* general purpose value-getter */
void
IFTE_get_Vals(double JDeph0, double JDeph1, int kind,
	      double *res)
{
  /* Get dates into one part that is exactly an integer + 0.5, plus
   * a fractional part. Here is the fortran code:
      S=ET2(1)-.5D0
      CALL SPLIT(S,PJD(1))
      CALL SPLIT(ET2(2),PJD(3))
      PJD(1)=PJD(1)+PJD(3)+.5D0
      PJD(2)=PJD(2)+PJD(4)
      CALL SPLIT(PJD(2),PJD(3))
      PJD(1)=PJD(1)+PJD(3)
  */
  double whole0, whole1, frac0, frac1;
  int irec;
  double t[2];
  int ncoeff = ifte.reclen/8;
  size_t nread;

  //  printf("Starting ITFE_get_vals with JDeph0 = %g JDeph1 = %g kind = %d\n",JDeph0,JDeph1,kind);
  
  whole0 = floor(JDeph0-0.5);
  frac0 = JDeph0-0.5-whole0;
  whole1 = floor(JDeph1);
  frac1 = JDeph1-whole1;
  whole0 += whole1 + 0.5;
  frac0 += frac1;
  whole1 = floor(frac0);
  frac1 = frac0-whole1;
  whole0 += whole1;

  JDeph0 = whole0;
  JDeph1 = frac1;

  if (JDeph0 < ifte.startJD) {
    fprintf (stderr, "Error: Requested JD=%lf is less than start JD=%lf\n",
	     JDeph0, ifte.startJD);
    exit(1);
  }

  /* CALCULATE RECORD # AND RELATIVE TIME IN INTERVAL */
  irec = (int)floor((JDeph0-ifte.startJD)/ifte.stepJD)+2;
  //  printf("Calc irec = %d %g %g %g\n",irec,JDeph0,ifte.startJD,ifte.stepJD);
  if (JDeph0 == ifte.endJD)
    irec--;
  //  printf("Calc irec2 = %d %g %g %g\n",irec,JDeph0,ifte.startJD,ifte.stepJD);
  t[0] = (JDeph0-(ifte.startJD+ifte.stepJD*(irec-2))+JDeph1)/ifte.stepJD;
  t[1] = ifte.stepJD;
  //  printf("In ifteph: tempo2 -- with ncoeff = %d\n",ncoeff);
  //  printf("JDeph0 = %g, JDeph1 = %g, ifte.startJD = %g ifte.endJD = %g\n",JDeph0,JDeph1,
  //	 ifte.startJD,ifte.endJD);
  //  printf("irec = %d. ifte.irec = %d, ifte.reclen = %d\n",irec,ifte.irec,ifte.reclen);
  //  printf("ifte.stepJD = %g\n",ifte.stepJD);
  /* READ CORRECT RECORD IF NOT IN CORE */
  if (irec != ifte.irec)
  {
    if (fseek(ifte.f, ifte.reclen*irec, SEEK_SET) < 0)  
    {
      perror("Error reading time ephemeris");
      exit(1);
    }
    nread = fread(ifte.buf, 8, ncoeff, ifte.f);
    if ((int)nread < ncoeff)
    {
      fprintf(stderr, "Error reading time ephemeris: Only read %d coefficients, wanted %d!\n", 
	      nread, ncoeff);
      exit(1);
    }
    if (ifte.swap_endian) IFTswapDoubles(ifte.buf, ncoeff); 
  }
  
  /*  INTERPOLATE time ephemeris */
  if (kind==1)
    IFTEinterp(&ifte.iinfo, ifte.buf+ifte.ipt[0][0]-1, t, 
	       ifte.ipt[0][1], 1, ifte.ipt[0][2], 2, res);
  else
    IFTEinterp(&ifte.iinfo, ifte.buf+ifte.ipt[1][0]-1, t, 
	       ifte.ipt[1][1], 3, ifte.ipt[1][2], 2, res);

}

/* convenience interfaces */
void IFTE_get_DeltaT_DeltaTDot(double Teph0, double Teph1,
				   double *DeltaT, double *DeltaTDot)
{
  double res[2];
  IFTE_get_Vals(Teph0, Teph1, 1, res);
  *DeltaT = res[0];
  *DeltaTDot = res[1];
}

double IFTE_DeltaT(double Teph0, double Teph1)
{
  double DeltaT, DeltaTDot;
  IFTE_get_DeltaT_DeltaTDot(Teph0, Teph1, &DeltaT, &DeltaTDot);
  return DeltaT;
}

double IFTE_DeltaTDot(double Teph0, double Teph1)
{
  double DeltaT, DeltaTDot;
  IFTE_get_DeltaT_DeltaTDot(Teph0, Teph1, &DeltaT, &DeltaTDot);
  return DeltaTDot;
}

void IFTE_get_vE_vEDot(double Teph0, double Teph1,
				   double *vE, double *vEDot)
{
  double res[6];
  int i;

  IFTE_get_Vals(Teph0, Teph1, 2, res);

  for (i=0; i < 3; i++)
  {
    vE[i] = res[i];
    vEDot[i] = res[i+3];
  }
}

void IFTE_get_vE(double Teph0, double Teph1, double *vE) 
{
  double vEDot[3];
  IFTE_get_vE_vEDot(Teph0, Teph1, vE, vEDot);
}

void IFTE_get_vEDot(double Teph0, double Teph1, double *vEDot)
{
  double vE[3];
  IFTE_get_vE_vEDot(Teph0, Teph1, vE, vEDot);
}

/*  The following routine is borrowed from the JPL ephemeris C code */
/*****************************************************************************
*        *****    jpl planetary and lunar ephemerides    *****     C ver.1.2 *
******************************************************************************
*                                                                            *
*  This program was written in standard fortran-77 and it was manually       *
*  translated to C language by Piotr A. Dybczynski (dybol@phys.amu.edu.pl),  *
*  subsequently revised heavily by Bill J Gray (pluto@gwi.net).              *
*                                                                            *
******************************************************************************/

/*****************************************************************************
**                     interp(buf,t,ncf,ncm,na,ifl,pv)                      **
******************************************************************************
**                                                                          **
**    this subroutine differentiates and interpolates a                     **
**    set of chebyshev coefficients to give position and velocity           **
**                                                                          **
**    calling sequence parameters:                                          **
**                                                                          **
**      input:                                                              **
**                                                                          **
**      iinfo   stores certain chunks of interpolation info,  in hopes      **
**              that if you call again with similar parameters,  the        **
**              function won't have to re-compute all coefficients/data.    **
**                                                                          **
**       coef   1st location of array of d.p. chebyshev coefficients        **
**              of position                                                 **
**                                                                          **
**          t   t[0] is double fractional time in interval covered by       **
**              coefficients at which interpolation is wanted               **
**              (0 <= t[0] <= 1).  t[1] is dp length of whole               **
**              interval in input time units.                               **
**                                                                          **
**        ncf   # of coefficients per component                             **
**                                                                          **
**        ncm   # of components per set of coefficients                     **
**                                                                          **
**         na   # of sets of coefficients in full array                     **
**              (i.e., # of sub-intervals in full interval)                 **
**                                                                          **
**         ifl  integer flag: =1 for positions only                         **
**                            =2 for pos and vel                            **
**                                                                          **
**                                                                          **
**      output:                                                             **
**                                                                          **
**    posvel   interpolated quantities requested.  dimension                **
**              expected is posvel[ncm*ifl], double precision.              **
**                                                                          **
*****************************************************************************/
static void IFTEinterp( struct IFTE_interpolation_info *iinfo,
                 const double coef[], const double t[2], const int ncf,
                 const int ncm, const int na, const int ifl, double posvel[])
{
  double dna, dt1, temp, tc, vfac, temp1;
  double *pc_ptr;
  int l, i, j; 

/*  entry point. get correct sub-interval number for this set
    of coefficients and then get normalized chebyshev time
    within that subinterval.                                             */

  dna = (double)na;
  modf( t[0], &dt1);
  temp = dna * t[0];
  l = (int)(temp - dt1);

/*  tc is the normalized chebyshev time (-1 <= tc <= 1)    */

  tc = 2.0 * (modf( temp, &temp1) + dt1) - 1.0;

/*  check to see whether chebyshev time has changed,
    and compute new polynomial values if it has.
    (the element iinfo->pc[1] is the value of t1[tc] and hence
    contains the value of tc on the previous call.)     */

  if(tc != iinfo->pc[1])
    {
      iinfo->np = 2;
      iinfo->nv = 3;
      iinfo->pc[1] = tc;
      iinfo->twot = tc+tc;
    }

/*  be sure that at least 'ncf' polynomials have been evaluated
    and are stored in the array 'iinfo->pc'.    */

  if(iinfo->np < ncf)
    {
    pc_ptr = iinfo->pc + iinfo->np;

    for(i=ncf - iinfo->np; i; i--, pc_ptr++)
       *pc_ptr = iinfo->twot * pc_ptr[-1] - pc_ptr[-2];
    iinfo->np=ncf;
    }

/*  interpolate to get position for each component  */

  for( i = 0; i < ncm; ++i)        /* ncm is a number of coordinates */
     {
     const double *coeff_ptr = coef + ncf * (i + l * ncm + 1);
     const double *pc_ptr = iinfo->pc + ncf;

     posvel[i]=0.0;
     for( j = ncf; j; j--)
        posvel[i] += (*--pc_ptr) * (*--coeff_ptr);
     }

   if(ifl <= 1) return;

/*  if velocity interpolation is wanted, be sure enough
    derivative polynomials have been generated and stored.    */

  vfac=(dna+dna)/t[1];
  iinfo->vc[2] = iinfo->twot + iinfo->twot;
  if( iinfo->nv < ncf)
     {
     double *vc_ptr = iinfo->vc + iinfo->nv;
     const double *pc_ptr = iinfo->pc + iinfo->nv - 1;

     for( i = ncf - iinfo->nv; i; i--, vc_ptr++, pc_ptr++)
        *vc_ptr = iinfo->twot * vc_ptr[-1] + *pc_ptr + *pc_ptr - vc_ptr[-2];
     iinfo->nv = ncf;
     }

/*  interpolate to get velocity for each component    */

   for( i = 0; i < ncm; ++i)
      {
      double tval = 0.;
      const double *coeff_ptr = coef + ncf * (i + l * ncm + 1);
      const double *vc_ptr = iinfo->vc + ncf;

      for( j = ncf; j; j--)
         tval += (*--vc_ptr) * (*--coeff_ptr);
      posvel[ i + ncm] = tval * vfac;
      }
   return;
}
