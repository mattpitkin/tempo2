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
#include <string.h>
#include "tempo2.h"

void formBats(pulsar *psr,int npsr)
{
  int i,p;
  //  double etatdm;
  double shapiroDelay;
  const char *CVS_verNum = "$Revision: 1.8 $";

  if (displayCVSversion == 1) CVSdisplayVersion("formBats.C","formBats()",CVS_verNum);

  logdbg("In formBats");

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
	     	   //if(i==0){printf("Forming bat with sat = %.15g TT = %.15g TT_TB = %.15g trop = %.15g roemer = %.15g shap = %.15g tdis1 = %.15g tdis2 = %.15g\n",(double)psr[p].obsn[i].sat,(double)getCorrectionTT(psr[p].obsn+i),(double)psr[p].obsn[i].correctionTT_TB,(double)psr[p].obsn[i].troposphericDelay,(double)psr[p].obsn[i].roemer,(double)shapiroDelay,(double)psr[p].obsn[i].tdis1,(double)psr[p].obsn[i].tdis2);}
		     

		long double batcorr=getCorrectionTT(psr[p].obsn+i)/SECDAY
                + (psr[p].obsn[i].correctionTT_TB
                   -psr[p].obsn[i].troposphericDelay
                   +psr[p].obsn[i].roemer -
                   shapiroDelay - psr[p].obsn[i].tdis1 - psr[p].obsn[i].tdis2)/SECDAY;

                psr[p].obsn[i].batCorr = batcorr;

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
  logdbg("Leaving formBats");
}

