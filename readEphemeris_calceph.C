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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tempo2.h"
#include "ifteph.h"

void convertUnits(double *val,int units);

#ifdef HAVE_CALCEPH
#include <calceph.h>

/*
 * readEphemeris_calceph
 *
 * ephemeris reading routine making use of the calceph library
 *
 * This routine sets:
 *  psr[p].obsn[i].sun_ssb[0->5]
 *  psr[p].obsn[i].sun_earth[0->5]
 *  psr[p].obsn[i].planet_ssb[0->9][0->5]
 *  psr[p].obsn[i].jupiter_earth[0->5]
 *  psr[p].obsn[i].saturn_earth[0->5]
 *  psr[p].obsn[i].venus_earth[0->5]
 *  psr[p].obsn[i].uranus_earth[0->5]
 *  psr[p].obsn[i].neptune_earth[0->5]
 *  psr[p].obsn[i].earth_ssb[0->5]
 * 
 * Also
 *  update SSB position for reflex motion due to error in planetary masses
 *  convert to ecliptic coordinates if requested
 *  set earth_ssb with TEL_DX, TEL_DY, TEL_DZ parameters are required
 */

void readEphemeris_calceph(pulsar *psr,int npsr)
{
  t_calcephbin *eph;
  int i,p;
  long double jd;
  double jd0,jd1;
 

  for (p=0;p<npsr;p++)
    {
      eph = calceph_open(psr[p].ephemeris);
      if (eph) {
	printf("Successfully opened ephemeris >%s<\n",psr[p].ephemeris);
      } else {
	printf("Error: unable to open ephemeris >%s< for pulsar >%s<\n",psr[p].ephemeris,psr[p].name);
	exit(1);
      }
      // Now read the ephemeris for each observation
      for (i=0;i<psr[p].nobs;i++)
	{
	  jd = psr[p].obsn[i].sat + getCorrectionTT(psr[p].obsn+i)/SECDAY + 
	    psr[p].obsn[i].correctionTT_Teph/SECDAY+2400000.5; 
	  jd0 = (double)((int)jd);
	  jd1 = (double)(jd-(int)jd);

	  // Calculate the Earth to SSB vector
	  calceph_compute_unit(eph,jd0,jd1,3,12,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,psr[p].obsn[i].earth_ssb);
	  convertUnits(psr[p].obsn[i].earth_ssb,psr[p].units);

	  // Calculate the Sun to SSB vector
	  calceph_compute_unit(eph,jd0,jd1,11,12,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,psr[p].obsn[i].sun_ssb);
	  convertUnits(psr[p].obsn[i].sun_ssb,psr[p].units);
	}	
      calceph_close(eph);
  }
}
#else
void readEphemeris_calceph(pulsar *psr,int npsr)
{
  printf("ERROR: unable to use calceph library routines as library not installed\n");
  exit(1);
}
#endif

void convertUnits(double *val,int units)
{
  double scale=1;
  if (units == SI_UNITS)
    scale = IFTE_K;
  val[0] = val[0]*1000.0*scale/SPEED_LIGHT;
  val[1] = val[1]*1000.0*scale/SPEED_LIGHT;
  val[2] = val[2]*1000.0*scale/SPEED_LIGHT;
  val[3] = val[3]*1000.0*scale/SPEED_LIGHT;
  val[4] = val[4]*1000.0*scale/SPEED_LIGHT;
  val[5] = val[5]*1000.0*scale/SPEED_LIGHT;
}
