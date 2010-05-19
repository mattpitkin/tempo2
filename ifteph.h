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

#ifndef IFTEPH_H
#define IFTEPH_H

#include "tempo2.h"

/* C Code for reading Irwin & Fukushima's time ephemeris

   Russell Edwards March 2005 */
#ifdef __cplusplus
extern "C" {
#endif


void IFTE_init(const char *fname);

void IFTE_get_DeltaT_DeltaTDot(double Teph0, double Teph1,
				   double *DeltaT, double *DeltaTDot);

double IFTE_DeltaT(double Teph0, double Teph1);
double IFTE_DeltaTDot(double Teph0, double Teph1);
  void IFTE_close_file();
void IFTE_get_vE_vEDot(double Teph0, double Teph1,
				   double *ve, double *vEDot);
void IFTE_get_vE(double Teph0, double Teph1, double *vE); 
void IFTE_get_vEDot(double Teph0, double Teph1, double *vEDot); 

/* constants determined by Irwin & Fukushima */

#define IFTE_JD0  2443144.5003725    /* Epoch of TCB, TCG and TT */
#define IFTE_MJD0 43144.0003725
#define IFTE_TEPH0 -65.564518e-6
#define IFTE_LC   1.48082686742e-8

//This is the value used by if99 : #define IFTE_KM1  1.55051979154e-8 */
// However we should use the IAU value of L_B that follows from
// their definition of L_G: L_B = 1.55051976772e-8, K=1/(1-L_B)
#define IFTE_KM1 1.55051979176e-8
#define IFTE_K    (((longdouble)1.0) + ((longdouble)IFTE_KM1))  /* needs quad precision */

#ifdef __cplusplus
}
#endif

#endif
