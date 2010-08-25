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
#include "tempo2.h"
//#include "tempo2Util.h"
//#include "tempo2pred.h"
//#include "tempo2pred_int.h"

char TEMPO2_ENVIRON[MAX_STRLEN]="TEMPO2";
char TEMPO2_ERROR[MAX_STRLEN]="";

int MAX_PSR   = MAX_PSR_VAL;    /* Maximum number of pulsars to fit simultaneously  */
int MAX_OBSN  = MAX_OBSN_VAL;
double ECLIPTIC_OBLIQUITY = ECLIPTIC_OBLIQUITY_VAL;
int debugFlag = 0;
int veryFast = 0;


void extra_delays(pulsar *psr,int npsr)
{  
  calculate_bclt(psr,npsr);/* 3. Calculate bclt  */
  /*  shapiro_delay(psr,npsr); */ /* 1. Calculate the Shapiro delay */
  /* dm_delays(psr,npsr); */    /* 2. Extra dispersion measure delays */  
}

void clock_corrections(pulsar *psr,int npsr)
{  
  if (debugFlag==1) printf("Calling toa2utc\n");
  toa2utc(psr,npsr);        /* 1. UTC(Observatory) -> UTC(NIST) */
  if (debugFlag==1) printf("Calling tai2ut1\n");
//   utc2tai(psr,npsr);     /* 2. UTC(NIST) -> TAI              */
  tai2ut1(psr,npsr);        /* 3. TAI -> UT1                    */
//   tai2tt(psr,npsr);      /* 4. TAI -> TT                     */
  if (debugFlag==1) printf("Calling tt2tb\n");
  tt2tb(psr,npsr);          /* 5. Rough estimate of TT-TB (+-2.2 microsec) */
  if (debugFlag==1) printf("Done clock corrections\n");
}

void ephemeris_routines(pulsar *psr,int npsr)
{ 
  vectorPulsar(psr,npsr);   /* 1. Form a vector pointing at the pulsar */
  readEphemeris(psr,npsr,0);/* 2. Read the ephemeris */
  get_obsCoord(psr,npsr);   /* 3. Get Coordinate of observatory relative to Earth's centre */
  tt2tb(psr,npsr);          /* Observatory/time-dependent part of TT-TB */
  readEphemeris(psr,npsr,0);  /* Re-evaluate ephemeris with correct TB */ 
}

void formBatsAll(pulsar *psr,int npsr)
{
  if (debugFlag==1) printf("Calling clock corrections\n");
  clock_corrections(psr,npsr);          /* Clock corrections  ... */  
  if (debugFlag==1) printf("Reading ephemeris routines\n");
  ephemeris_routines(psr,npsr);         /* Ephemeris routines ... */
  if (debugFlag==1) printf("Reading extra delays\n");
  extra_delays(psr,npsr);               /* Other time delays  ... */
  formBats(psr,npsr);                   /* Form Barycentric arrival times */
  secularMotion(psr,npsr); 
}

// Only recalculate that which is likely
// to change if psr position has been altered.
void updateBatsAll(pulsar *psr, int npsr)
{
  vectorPulsar(psr, npsr);
  calculate_bclt(psr, npsr);
  formBats(psr, npsr);
  secularMotion(psr,npsr); 
}
