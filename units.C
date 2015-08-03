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


#include "tempo2.h"
#include "ifteph.h"

#include <math.h>

/* redwards function to transform units of paramaters */

/* helpers */
void scale_param(struct parameter *p, int arr,longdouble f)
{
  p->val[arr] *= f;
  p->err[arr] *= f;
  p->prefit[arr] *= f;
  p->prefitErr[arr] *=f;
}

void xform_mjd(struct parameter *p,int arr, longdouble f)
{
  p->val[arr] = (p->val[arr] - IFTE_MJD0)*f + IFTE_MJD0;
  p->err[arr] *= f;
  p->prefit[arr] = (p->prefit[0] - IFTE_MJD0)*f + IFTE_MJD0;
  p->prefitErr[arr] *= f;
}
 
void
transform_units(struct pulsar *psr, int from, int to)
{
  const char *CVS_verNum = "$Revision: 1.6 $";

  if (displayCVSversion == 1) CVSdisplayVersion("units.C","transform_units()",CVS_verNum);

  longdouble f = 
    (to == TDB_UNITS ? ((longdouble)1.0) : ((longdouble)IFTE_K))
     / (from == TDB_UNITS ? ((longdouble)1.0) : ((longdouble)IFTE_K)) ;
  longdouble one = LONGDOUBLE_ONE;
  longdouble val;
  int k;

  // all length and time dimensions get scaled by f ....
  // these are done in same order as initialise.C
  scale_param(&psr->param[param_f],0, one/f);
  scale_param(&psr->param[param_f],1, one/(f*f));
  scale_param(&psr->param[param_f],2, one/(f*f*f));
  scale_param(&psr->param[param_f],3, one/(f*f*f*f));
  scale_param(&psr->param[param_f],4, one/(f*f*f*f*f));
  scale_param(&psr->param[param_f],5, one/(f*f*f*f*f*f));
  scale_param(&psr->param[param_f],6, one/(f*f*f*f*f*f*f));
  scale_param(&psr->param[param_f],7, one/(f*f*f*f*f*f*f*f));
  scale_param(&psr->param[param_f],8, one/(f*f*f*f*f*f*f*f*f));
  scale_param(&psr->param[param_f],9, one/(f*f*f*f*f*f*f*f*f*f));
  scale_param(&psr->param[param_fddc],0, pow(f, one+(longdouble)psr->param[param_fddi].val[0]));

  // DM is physically Length^-2 . However, it is actually defined in terms of
  // a published constant, not physical constants. If this constant is left 
  // the same then really DM absorbs the dimensions of the constant, and becomes
  // time delay per square frequency, i.e. Time^3
  // BUT! These are barycentric frequencies! Tempo neglected the "Einstein rate"
  // correction to the frequency; this is not surprising as it amounts to
  // only 4x10^-10 in rate, or about 10^-9 in DM delay ...
  // OTOH using TCB there is a mean difference of ~1.5x10^-8 in rate,
  // or 3x10^-8 in DM delay, which could become significant
  // XXX Check this! My new code gives ~1.5x10^-8 for TDB itself!

  val=f*f*f;
  for (k=0;k<psr->param[param_dm].aSize;k++)
    {
      scale_param(&psr->param[param_dm],k,  val);
      val/=f;  /* Is this correct: note dm1 and dm2 in the old method */
      /*      scale_param(&psr->param[param_dm],k, f*f);
	      scale_param(&psr->param[param_dm],k, f);
	      scale_param(&psr->param[param_dm],k, one/f);
	      scale_param(&psr->param[param_dm],k, one/(f*f));
	      scale_param(&psr->param[param_dm],k, one/(f*f*f));
	      scale_param(&psr->param[param_dm],k, one/(f*f*f*f));
	      scale_param(&psr->param[param_dm],k, one/(f*f*f*f*f));
	      scale_param(&psr->param[param_dm],k, one/(f*f*f*f*f*f));
	      scale_param(&psr->param[param_dm],k, one/(f*f*f*f*f*f*f)); */
    }
  // Parallax is quoted as mas but actually it enters the timing model
  // with units seconds of delay per square distance / time/distance , 
  // ie length^-1 , which makes sense as it is 1/distance to source
  //scale_param(&psr->param[param_px], one/f);
  // XXXXXXXXX NOTE!!! 1 AU will have a different definition in the two
  // frames!! Check that this is the case and the old value / tempo value
  // is actually correct for Ephemeris units
  // XXXX Is that right??? The definition of a "parsec" should be the 
  // same in the two systems so there should be no transformation (!?)
  // commented out!

  // Radial velocity in mas/yr (???) 
  scale_param(&psr->param[param_pmrv],0, one/f);
  scale_param(&psr->param[param_pmra],0, one/f);
  scale_param(&psr->param[param_pmdec],0, one/f);
  
  // Epoch transformations ... these have to be translated to the common
  // origin before scaling
  xform_mjd(&psr->param[param_posepoch],0, f);
  xform_mjd(&psr->param[param_pepoch],0, f);
  xform_mjd(&psr->param[param_start],0, f);
  xform_mjd(&psr->param[param_finish],0, f);

  // Glitch parameters
  // XXX Need to know wtf these glitch parameters are before xforming them

  // Binary parameters
  xform_mjd(&psr->param[param_t0],0, f);
  scale_param(&psr->param[param_pb],0, f);
  scale_param(&psr->param[param_fb],0, one/f);
  scale_param(&psr->param[param_a1],0, f);
  scale_param(&psr->param[param_edot],0, one/f);
  // XXX need to investigate DR, DTHETA, GAMMA
  scale_param(&psr->param[param_omdot],0, one/f);
  xform_mjd(&psr->param[param_tasc],0, f);

  xform_mjd(&psr->param[param_tzrmjd],0, f); // XXX barycentric (TB) MJD?
  scale_param(&psr->param[param_tzrfrq],0, one/f); // XXX Barycentric freq??

  for (k=0;k<psr->param[param_bpjep].aSize;k++)
    {
      xform_mjd(&psr->param[param_bpjep],k, f);
      scale_param(&psr->param[param_bpja1],k, f); 
      scale_param(&psr->param[param_bpjpb],k, f); 
    }
}
