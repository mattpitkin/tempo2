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

  if (debugFlag==1) printf("In calculate_bclt with number of psr = %d, nobs (psr[0]) = %d\n",npsr,psr[0].nobs);

  /* Conversion of mas to radians for parallax */
  pxConv = 1.74532925199432958E-2/3600.0e3;

  for (p=0;p<npsr;p++)
    {
      if (psr[p].correctTroposphere==1) compute_tropospheric_delays(&psr[p], 1);


      for (i=0;i<psr[p].nobs;i++)	
	{
	  if (debugFlag==1)printf("In tdis2 calculate_bclt with observation %d %d\n",i,psr[p].obsn[i].delayCorr);
	  if (psr[p].obsn[i].delayCorr==0) /* No correction */
	    {
	      psr[p].obsn[i].freqSSB = 0.0; //psr[p].obsn[i].freq;
	      psr[p].obsn[i].roemer = 0.0;
	    }
	  else
	    {
	      /* Calculate vector from SSB to observatory */
	      /*	      for (j=0;j<3;j++)
		rca[j] = psr[p].obsn[i].earthMoonBary_ssb[j] - psr[p].obsn[i].earthMoonBary_earth[j] + 
		psr[p].obsn[i].observatory_earth[j]; */
	      for (j=0;j<3;j++)
		rca[j] = psr[p].obsn[i].earth_ssb[j] + psr[p].obsn[i].observatory_earth[j]; 
	      /* Radial velocity */
	      if (psr[p].param[param_pmrv].paramSet[0]==1)
		pmrvrad = psr[p].param[param_pmrv].val[0]*(2*M_PI/360.0)/36000.0; /* In radians/century - CHECK 36000*/
	      else
		pmrvrad = 0.0;
	      
	      dt_SSB = 0.0;
	      
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

		//		printf("Calculating roemer using %g %g %g %g %g (delt = %g; dt_SSB = %g)\n",(double)rcos1,(double)dt_pm, (double)dt_pmtt,(double)dt_px,(double)dt_pmtr,(double)delt,(double)dt_SSB);
		psr[p].obsn[i].roemer = rcos1 + dt_pm + dt_pmtt + dt_px + dt_pmtr;		
		shapiro_delay(psr,npsr,p,i,delt,dt_SSB); /* Now calculate the Shapiro delay */
		if (debugFlag==1)printf("In tdis2 calculate_bclt with observation %d %d calling dmdelays\n",i,psr[p].obsn[i].delayCorr);
		dm_delays(psr,npsr,p,i,delt,dt_SSB);     /* Now calculate the dispersion measure delays */
		dt_SSB = psr[p].obsn[i].roemer-(psr[p].obsn[i].tdis1+psr[p].obsn[i].tdis2)-
		  (psr[p].obsn[i].shapiroDelaySun+psr[p].planetShapiro*psr[p].obsn[i].shapiroDelayJupiter); 
		//		printf("dt_SSB: %g %g %g %g %g \n",(double)psr[p].obsn[i].roemer,(double)psr[p].obsn[i].tdis1,(double)psr[p].obsn[i].tdis2,(double)psr[p].obsn[i].shapiroDelaySun,(double)psr[p].planetShapiro*psr[p].obsn[i].shapiroDelayJupiter);
		if (veryFast==1) break;
	      } while (fabs(dt_SSB-dt_SSB_old)>1.0e-10 && psr[p].obsn[i].deleted!=1);
	    }
	}
    }
}
