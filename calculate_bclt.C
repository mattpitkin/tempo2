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

void calculate_bclt(pulsar *psr,int npsr)
{
  int i,j,p;
  double dt_px,dt_pm,dt_pmtt,dt_pmtr,pmtrans_rcos2,dt_SSB,dt_SSB_old,rcos1,pmtrans,delt;
  double pmrvrad;
  double rr,pxConv;
  double rca[3];
  int loop;
  const char *CVS_verNum = "$Revision: 1.14 $";

  if (displayCVSversion == 1) CVSdisplayVersion("calculate_bclt.C","calculate_bclt()",CVS_verNum);

  logdbg("In calculate_bclt with number of psr = %d, nobs (psr[0]) = %d",npsr,psr[0].nobs);

  /* Conversion of mas to radians for parallax */
  pxConv = 1.74532925199432958E-2/3600.0e3;

  for (p=0;p<npsr;p++)
    {
      if (psr[p].correctTroposphere==1) {
	//	printf("Correcting troposphere\n");
	logdbg("Correcting troposphere");
	compute_tropospheric_delays(&psr[p], 1);
	logdbg("Complete correcting troposphere");
      }

      for (i=0;i<psr[p].nobs;i++)	
	{
	  if (psr[p].correctTroposphere==0) psr[0].obsn[i].troposphericDelay=0;
	  logdbg("In tdis2 calculate_bclt with observation %d %d",i,psr[p].obsn[i].delayCorr);
	  if (psr[p].obsn[i].delayCorr==0) /* No correction */
	    {
	      psr[p].obsn[i].freqSSB = psr[p].obsn[i].freq;
	      psr[p].obsn[i].roemer = 0.0;
	    }
	  else
	    {
	      /* Calculate vector from SSB to observatory */
	      /*	      for (j=0;j<3;j++)
		rca[j] = psr[p].obsn[i].earthMoonBary_ssb[j] - psr[p].obsn[i].earthMoonBary_earth[j] + 
		psr[p].obsn[i].observatory_earth[j]; */
	      //	      printf("roemer: %g %g %g %g %g %g\n",psr[p].obsn[i].earth_ssb[0],psr[p].obsn[i].earth_ssb[1],psr[p].obsn[i].earth_ssb[2],psr[p].obsn[i].observatory_earth[0],psr[p].obsn[i].observatory_earth[1],psr[p].obsn[i].observatory_earth[2]);
	      
	      if (strcmp(psr[p].obsn[i].telID,"STL_FBAT")==0)
		{ // Observatory_earth is actually observatory to ssb
		  for (j=0;j<3;j++)
		    rca[j] = psr[p].obsn[i].observatory_earth[j]; 
		}
	      else
		{
		  for (j=0;j<3;j++)
		    rca[j] = psr[p].obsn[i].earth_ssb[j] + psr[p].obsn[i].observatory_earth[j]; 
		}
	      /* Radial velocity */
	      if (psr[p].param[param_pmrv].paramSet[0]==1)
		pmrvrad = psr[p].param[param_pmrv].val[0]*(2*M_PI/360.0)/36000.0; /* In radians/century - CHECK 36000*/
	      else
		pmrvrad = 0.0;
	      
	      dt_SSB = 0.0;
	      loop=0;
	      do {
		dt_SSB_old = dt_SSB;
		
		rcos1 = dotproduct(psr[p].posPulsar,rca);
		rr    = dotproduct(rca,rca);
		delt = (psr[p].obsn[i].sat-psr[p].param[param_posepoch].val[0] + 
			(getCorrectionTT(psr[p].obsn+i)+psr[p].obsn[i].correctionTT_TB+dt_SSB)/SECDAY)/36525.0;
		
		pmtrans_rcos2 = dotproduct(psr[p].velPulsar,rca);      
		pmtrans       = sqrt(dotproduct(psr[p].velPulsar,psr[p].velPulsar));
		dt_pm         = (delt)*pmtrans_rcos2;
		dt_pmtt       = -0.5*pmtrans*pmtrans*(delt)*(delt)*rcos1;
			       		
		/* Parallax */
		dt_px = -0.5*psr[p].param[param_px].val[0]*pxConv*(rr-rcos1*rcos1)/AULTSC;
		/* dt_pmtr = transverse velocity */
		dt_pmtr = -pow((delt),2)*pmrvrad*pmtrans_rcos2;

		logdbg("Calculating roemer using %g %g %g %g %g (delt = %g; dt_SSB = %g)",(double)rcos1,(double)dt_pm, (double)dt_pmtt,(double)dt_px,(double)dt_pmtr,(double)delt,(double)dt_SSB);		
		psr[p].obsn[i].roemer = rcos1 + dt_pm + dt_pmtt + dt_px + dt_pmtr;		
		shapiro_delay(psr,npsr,p,i,delt,dt_SSB); /* Now calculate the Shapiro delay */
		logdbg("In tdis2 calculate_bclt with observation %d %d calling dmdelays",i,psr[p].obsn[i].delayCorr);
		dm_delays(psr,npsr,p,i,delt,dt_SSB);     /* Now calculate the dispersion measure delays */

		logdbg("Complete dm_delays");
		dt_SSB = psr[p].obsn[i].roemer-(psr[p].obsn[i].tdis1+psr[p].obsn[i].tdis2)-
		  (psr[p].obsn[i].shapiroDelaySun+psr[p].planetShapiro*psr[p].obsn[i].shapiroDelayJupiter); 
				logdbg("dt_SSB: %g %g %g %g %g",(double)psr[p].obsn[i].roemer,(double)psr[p].obsn[i].tdis1,(double)psr[p].obsn[i].tdis2,(double)psr[p].obsn[i].shapiroDelaySun,(double)psr[p].planetShapiro*psr[p].obsn[i].shapiroDelayJupiter);
		if (veryFast==1) break;
		loop++;
	      } while (fabs(dt_SSB-dt_SSB_old)>1.0e-10 && psr[p].obsn[i].deleted!=1 && loop < 100);
	      if (loop==100)
		{
		  printf("Warning: the loop to obtain dt_SSB reached 100 iterations\n");
		}
	      //	      printf("roemer %g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %d\n",(double)psr[p].obsn[i].sat,(double)psr[p].obsn[i].roemer,(double)rcos1,(double)dt_pm,(double)dt_pmtt,(double)dt_px,(double)dt_pmtr,(double)rr,loop);
	      //	      printf("roemer %g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",(double)psr[p].obsn[i].sat,psr[p].posPulsar[0],psr[p].posPulsar[1],psr[p].posPulsar[2],rca[0],rca[1],rca[2],rr,(double)psr[p].obsn[i].roemer);
	    }
	}
    }
  logdbg("Leaving calculate_bclt");
}
