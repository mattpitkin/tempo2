#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"

void secularMotion(pulsar *psr,int npsr)
{
  int i,j,p;
  longdouble t0;
  longdouble kpc2m = 3.08568025e19L;           /* 1 kpc in m        */
  longdouble mas_yr2rad_s = 1.536281850e-16L;  /* 1 mas/yr in rad/s */
  for (p=0;p<npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)	
	{
	  if (psr[0].param[param_dshk].paramSet[0]==1 && 
	      (psr[0].param[param_pmra].paramSet[0]==1 || psr[0].param[param_pmdec].paramSet[0]==1))
	    {
	      t0 = (psr[p].obsn[i].bat - psr[p].param[param_posepoch].val[0])*86400.0L;	      
	      psr[p].obsn[i].shklovskii = t0*t0/2.0L/299792458.0L*
		(getParameterValue(&psr[p],param_dshk,0)*kpc2m)*
		(psr[p].param[param_pmra].val[0]*psr[p].param[param_pmra].val[0]*
		 mas_yr2rad_s*mas_yr2rad_s+
		 psr[p].param[param_pmdec].val[0]*psr[p].param[param_pmdec].val[0]*
		 mas_yr2rad_s*mas_yr2rad_s);
	    }
	  else
	    psr[p].obsn[i].shklovskii = 0.0;
	  
	  psr[p].obsn[i].bbat = psr[p].obsn[i].bat - (psr[p].obsn[i].shklovskii)/86400.0L;
	}
    }
}
