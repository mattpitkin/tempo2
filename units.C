
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
