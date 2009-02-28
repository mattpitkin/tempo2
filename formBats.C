#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"

void formBats(pulsar *psr,int npsr)
{
  int i,p;
  //  double etatdm;
  double shapiroDelay;

  if (debugFlag==1)printf("In formBats\n");

  for (p=0;p<npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)
	{
	  //	  if (psr[p].obsn[i].clockCorr==0 || psr[p].obsn[i].delayCorr==0)
	  if (psr[p].obsn[i].delayCorr==0)
	    psr[p].obsn[i].bat = psr[p].obsn[i].sat;
	  else
	    {
	      /* ORIGINAL TEMPO calculation */
	      if (psr[p].calcShapiro==-1)
		shapiroDelay = 0.0;
	      else
		shapiroDelay = psr[p].obsn[i].shapiroDelaySun + psr[p].planetShapiro*
		  (psr[p].obsn[i].shapiroDelayVenus+
		   psr[p].obsn[i].shapiroDelayJupiter+psr[p].obsn[i].shapiroDelaySaturn
		   +psr[p].obsn[i].shapiroDelayUranus + psr[p].obsn[i].shapiroDelayNeptune);
	      //	      	      printf("Forming bat with sat = %g TT = %g TT_TB = %g trop = %g roemer = %g shap = %g tdis1 = %g tdis2 = %g\n",(double)psr[p].obsn[i].sat,(double)getCorrectionTT(psr[p].obsn+i),(double)psr[p].obsn[i].correctionTT_TB,(double)psr[p].obsn[i].troposphericDelay,(double)psr[p].obsn[i].roemer,(double)shapiroDelay,(double)psr[p].obsn[i].tdis1,(double)psr[p].obsn[i].tdis2);
	      psr[p].obsn[i].bat = psr[p].obsn[i].sat + 
		getCorrectionTT(psr[p].obsn+i)/SECDAY 
		+ (psr[p].obsn[i].correctionTT_TB
		   -psr[p].obsn[i].troposphericDelay
		   +psr[p].obsn[i].roemer - 
		   shapiroDelay - psr[p].obsn[i].tdis1 - psr[p].obsn[i].tdis2)/SECDAY;

	      /* Now calculate binary barycentric arrival time (bbat) */
	      psr[p].obsn[i].bbat = psr[p].obsn[i].bat - (psr[p].obsn[i].shklovskii)/SECDAY;
	    }
	}
    }
  if (debugFlag==1)printf("Leaving formBats\n");
}

