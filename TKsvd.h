//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

#ifndef __TKsvd_h
#define __TKsvd_h

#include "TKmatrix.h"
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

#ifdef __cplusplus
namespace TK{
// class definition
template<typename D>
class SVD {
    public:
        SVD(const TK::matrix<D> &in);
        SVD(TK::matrix<D> &in);

    private:
       TK::matrix<D> in;
    public:
       TK::matrix<D> U;
       TK::matrix<D> V;
       TK::diagonal<D> S;
};

}




void TKsingularValueDecomposition_lsq(longdouble **designMatrix,int n,int nf,longdouble **v,longdouble *w,longdouble **u);
void TKbacksubstitution_svd(longdouble **V, longdouble *w,longdouble **U,longdouble *b,longdouble *x,int n,int nf);
longdouble TKpythag(longdouble a,longdouble b);
void TKbidiagonal(longdouble **a,longdouble *anorm,int ndata,int nfit,longdouble **v,longdouble *w,longdouble **u,longdouble *rv1);

void TKsingularValueDecomposition_lsq(double **designMatrix,int n,int nf,double **v,double *w,double **u);
void TKbacksubstitution_svd(double **V, double *w,double **U,double *b,double *x,int n,int nf);
double TKpythag(double a,double b);
void TKbidiagonal(double **a,double *anorm,int ndata,int nfit,double **v,double *w,double **u,double *rv1);

void TKsingularValueDecomposition_lsq(float **designMatrix,int n,int nf,float **v,float *w,float **u);
void TKbacksubstitution_svd(float **V, float *w,float **U,float *b,float *x,int n,int nf);
float TKpythag(float a,float b);
void TKbidiagonal(float **a,float *anorm,int ndata,int nfit,float **v,float *w,float **u,float *rv1);

#endif
#endif
