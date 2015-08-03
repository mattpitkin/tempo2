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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"


// NOTE: DOES NOT MAKE USE OF DT_SSB
void shapiro_delay(pulsar *psr,int npsr,int p,int i,double delt,double dt_SSB)
{
  double delay,ctheta,r,rsa[3],pospos;
  int j,k;
  const char *CVS_verNum = "$Revision: 1.9 $";

  if (displayCVSversion == 1) CVSdisplayVersion("shapiro_delay.C","shapiro_delay()",CVS_verNum);
  
  if ((strcmp(psr[p].obsn[i].telID,"STL_BAT")==0) || psr[p].obsn[i].delayCorr==0) /* No correction */
  //if ( psr[p].obsn[i].delayCorr==0) /* No correction */
    {
      psr[p].obsn[i].shapiroDelaySun = 0.0;
      psr[p].obsn[i].shapiroDelayJupiter = 0.0;
    }
  else
    {
      for (k=0;k<3;k++)
	psr[p].obsn[i].psrPos[k] = psr[p].posPulsar[k]+delt*psr[p].velPulsar[k];      
      
      pospos = sqrt(psr[p].obsn[i].psrPos[0]*psr[p].obsn[i].psrPos[0]+
		    psr[p].obsn[i].psrPos[1]*psr[p].obsn[i].psrPos[1]+
		    psr[p].obsn[i].psrPos[2]*psr[p].obsn[i].psrPos[2]);
      for (k=0;k<3;k++)
	psr[p].obsn[i].psrPos[k] /= pospos;

      /* Note rsa = vector from observatory to centre of Sun */
      for (j=0;j<3;j++)
	{
	  /* Require vector from solar system body (e.g. Sun) to observatory at the time of the 
	   * closest approach of the photon to that body */
	  if (strcmp(psr[p].obsn[i].telID,"STL_FBAT")==0)
	    rsa[j] = -psr[p].obsn[i].sun_ssb[j] + psr[p].obsn[i].observatory_earth[j]; 
	  else
	    rsa[j] = -psr[p].obsn[i].sun_ssb[j] + psr[p].obsn[i].earth_ssb[j] + psr[p].obsn[i].observatory_earth[j]; 

	  /*	  rsa[j] = -psr[p].obsn[i].sun_ssb[j] + psr[p].obsn[i].earthMoonBary_ssb[j] -
		  psr[p].obsn[i].earthMoonBary_earth[j] + psr[p].obsn[i].observatory_earth[j]; */
	}
      r = sqrt(dotproduct(rsa,rsa));
      ctheta = dotproduct(psr[p].obsn[i].psrPos,rsa)/r;  
      delay = -2.0*GM_C3 * log(r/AULTSC*(1.0+ctheta));/* Note log is `ln' not `log10' */
      psr[p].obsn[i].shapiroDelaySun = delay;
      //      printf("Have %g [%g %g %g; %g %g %g]\n",(double)delay,rsa[0],rsa[1],rsa[2],
      //	     psr[p].obsn[i].observatory_earth[0],psr[p].obsn[i].observatory_earth[1],
      //	     psr[p].obsn[i].observatory_earth[2]);

      if (strcmp(psr[p].obsn[i].telID,"STL_FBAT")!=0)
	{
	  /* Now calculate the Shapiro delay due to Jupiter */
	  for (j=0;j<3;j++)
	    rsa[j] = psr[p].obsn[i].jupiter_earth[j] + psr[p].obsn[i].observatory_earth[j];
	  r = sqrt(dotproduct(rsa,rsa));
	  ctheta = dotproduct(psr[p].obsn[i].psrPos,rsa)/r;  
	  delay = -2.0*GMJ_C3 * log(r/AULTSC*(1.0+ctheta));  
	  psr[p].obsn[i].shapiroDelayJupiter = psr[p].planetShapiro*delay;
	  //	  printf("ctheta = %g\n",ctheta);
	  /* Now calculate the Shapiro delay due to Saturn */
	  for (j=0;j<3;j++)
	    rsa[j] = psr[p].obsn[i].saturn_earth[j] + psr[p].obsn[i].observatory_earth[j];
	  r = sqrt(dotproduct(rsa,rsa));
	  ctheta = dotproduct(psr[p].obsn[i].psrPos,rsa)/r;  
	  delay = -2.0*GMS_C3 * log(r/AULTSC*(1.0+ctheta));  
	  psr[p].obsn[i].shapiroDelaySaturn = psr[p].planetShapiro*delay;
	  
	  /* Now calculate the Shapiro delay due to Venus */
	  for (j=0;j<3;j++)
	    rsa[j] = psr[p].obsn[i].venus_earth[j] + psr[p].obsn[i].observatory_earth[j];
	  r = sqrt(dotproduct(rsa,rsa));
	  ctheta = dotproduct(psr[p].obsn[i].psrPos,rsa)/r;  
	  delay = -2.0*GMV_C3 * log(r/AULTSC*(1.0+ctheta));  
	  psr[p].obsn[i].shapiroDelayVenus = psr[p].planetShapiro*delay;
	  
	  /* Now calculate the Shapiro delay due to Uranus */
	  for (j=0;j<3;j++)
	    rsa[j] = psr[p].obsn[i].uranus_earth[j] + psr[p].obsn[i].observatory_earth[j];
	  r = sqrt(dotproduct(rsa,rsa));
	  ctheta = dotproduct(psr[p].obsn[i].psrPos,rsa)/r;  
	  delay = -2.0*GMU_C3 * log(r/AULTSC*(1.0+ctheta));  
	  psr[p].obsn[i].shapiroDelayUranus = psr[p].planetShapiro*delay;
	  
	  /* Now calculate the Shapiro delay due to Neptune */
	  for (j=0;j<3;j++)
	    rsa[j] = psr[p].obsn[i].neptune_earth[j] + psr[p].obsn[i].observatory_earth[j];
	  r = sqrt(dotproduct(rsa,rsa));
	  ctheta = dotproduct(psr[p].obsn[i].psrPos,rsa)/r;  
	  delay = -2.0*GMN_C3 * log(r/AULTSC*(1.0+ctheta));  
	  psr[p].obsn[i].shapiroDelayNeptune = psr[p].planetShapiro*delay;
	}
      else
	{
	  psr[p].obsn[i].shapiroDelayNeptune = 0.0;
	  psr[p].obsn[i].shapiroDelayJupiter = 0.0;
	  psr[p].obsn[i].shapiroDelaySaturn = 0.0;
	  psr[p].obsn[i].shapiroDelayUranus = 0.0;
	  psr[p].obsn[i].shapiroDelayVenus = 0.0;
	}
    }
  recordPrecision(&psr[p],1.0e-10,"SS. Shapiro delay","Total guess!");
}


