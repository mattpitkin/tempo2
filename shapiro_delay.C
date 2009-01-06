#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"

void shapiro_delay(pulsar *psr,int npsr,int p,int i,double delt,double dt_SSB)
{
  double delay,ctheta,r,rsa[3],pospos;
  int j,k;
  
  if (psr[p].obsn[i].delayCorr==0) /* No correction */
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
	  rsa[j] = -psr[p].obsn[i].sun_ssb[j] + psr[p].obsn[i].earth_ssb[j] + psr[p].obsn[i].observatory_earth[j]; 

	  /*	  rsa[j] = -psr[p].obsn[i].sun_ssb[j] + psr[p].obsn[i].earthMoonBary_ssb[j] -
		  psr[p].obsn[i].earthMoonBary_earth[j] + psr[p].obsn[i].observatory_earth[j]; */
	}
      r = sqrt(dotproduct(rsa,rsa));
      ctheta = dotproduct(psr[p].obsn[i].psrPos,rsa)/r;  
      delay = -2.0*GM_C3 * log(r/AULTSC*(1.0+ctheta));/* Note log is `ln' not `log10' */
      psr[p].obsn[i].shapiroDelaySun = delay;

      /* Now calculate the Shapiro delay due to Jupiter */
      for (j=0;j<3;j++)
	rsa[j] = psr[p].obsn[i].jupiter_earth[j] + psr[p].obsn[i].observatory_earth[j];
      r = sqrt(dotproduct(rsa,rsa));
      ctheta = dotproduct(psr[p].obsn[i].psrPos,rsa)/r;  
      delay = -2.0*GMJ_C3 * log(r/AULTSC*(1.0+ctheta));  
      psr[p].obsn[i].shapiroDelayJupiter = psr[p].planetShapiro*delay;

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
  recordPrecision(&psr[p],1.0e-10,"SS. Shapiro delay","Total guess!");
}


