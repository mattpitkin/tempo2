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

void secularMotion(pulsar *psr,int npsr)
{
  int i,j,p;
  longdouble t0;
  longdouble kpc2m = 3.08568025e19L;           /* 1 kpc in m        */
  longdouble mas_yr2rad_s = 1.536281850e-16L;  /* 1 mas/yr in rad/s */
  const char *CVS_verNum = "$Revision: 1.5 $";

  if (displayCVSversion == 1) CVSdisplayVersion("secularMotion.C","secularMotion()",CVS_verNum);

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
