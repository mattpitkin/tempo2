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
#include "tempo2.h"

/* This binary model is based on the DD model, but also includes
 * the effects of orbital parallax (Kopeikin, ApJ 439, L5 1995)
 * and secular variations of x, w and i due to proper motion of the binary
 * (Kopeikin, ApJ 467, L94, 1996
 *
 * Partly based on the work by Willem van Straten in bnryddk.f
 */

double DDKmodel(pulsar *psr,int p,int ipos,int param)
{
  double an;
  double pb,k;
  double rad2deg = 180.0/M_PI;
  double SUNMASS = 4.925490947e-6;
  double m2,tt0,t0,ecc,er,xdot,edot,dr,dth,eth,am2,ct;
  double pbdot,xpbdot,phase,u,du,gamma;
  double orbits;
  int norbits;
  double  cu,onemecu,cae,sae,ae,omega,omz,sw,cw,alpha,beta,bg,dre,drep,drepp,anhat,su;
  double sqr1me2,cume,brace,si,dlogbr,ds,da,a0,b0,d2bar,torb,ci,ti;
  double csigma,ce,cx,comega,cgamma,cdth,cm2,csi;
  double sin_omega,cos_omega,ki;
  double pmra,pmdec;
  double ki_dot,sini,cosi,tani,delta_i0,delta_j0,asi;
  double cos_alpha,sin_alpha,cos_delta,sin_delta,xpr,ypr,dpara;
  double pxConv = 1.74532925199432958E-2/3600.0e3;
  const char *CVS_verNum = "$Revision: 1.5 $";

  if (displayCVSversion == 1) CVSdisplayVersion("DDKmodel.C","DDKmodel()",CVS_verNum);

  dr  = 0.0; /* WHAT SHOULD THESE BE SET TO? */
  dth = 0.0; 

  /* Obtain KOM: the longitude of the ascending node */
  sin_omega = sin(psr[p].param[param_kom].val[0]*M_PI/180.0);
  cos_omega = cos(psr[p].param[param_kom].val[0]*M_PI/180.0);

  /* ... and KIN, the inclination angle */
  ki = psr[p].param[param_kin].val[0]*M_PI/180.0;       
  si = sin(ki);
  ci = cos(ki);
  ti = si/ci;

  /* Obtain proper motion in radians/sec */
  pmra  = psr[p].param[param_pmra].val[0]*M_PI/(180.0*3600.0e3)/(365.25*86400.0);
  pmdec = psr[p].param[param_pmdec].val[0]*M_PI/(180.0*3600.0e3)/(365.25*86400.0);

  /* Obtain parallax */
  if (psr[p].param[param_px].paramSet[0]==1)
    dpara = psr[p].param[param_px].val[0]*pxConv;
  else dpara = 0.0;

  /* Obtain vector pointing at the pulsar */
  sin_delta = psr[p].obsn[ipos].psrPos[2];
  cos_delta = cos(asin(sin_delta));
  sin_alpha = psr[p].obsn[ipos].psrPos[1]/cos_delta;
  cos_alpha = psr[p].obsn[ipos].psrPos[0]/cos_delta;

  if (psr[p].param[param_m2].paramSet[0]==1) am2 = psr[p].param[param_m2].val[0];
  else am2 = 0.0;

  pb = psr[p].param[param_pb].val[0]*SECDAY;
  an = 2.0*M_PI/pb;
  k = psr[p].param[param_omdot].val[0]/(rad2deg*365.25*86400.0*an);

  m2 = am2*SUNMASS;
  t0 = psr[p].param[param_t0].val[0];
  ct = psr[p].obsn[ipos].bbat;    
  
  tt0 = (ct-t0)*SECDAY;

  if (psr[p].param[param_gamma].paramSet[0]==1)
    gamma = psr[p].param[param_gamma].val[0];
  else
    gamma = 0.0;
  a0    = 0.0; /* WHAT SHOULD THIS BE SET TO? */
  b0    = 0.0; /* WHAT SHOULD THIS BE SET TO? */

  if (psr[p].param[param_om].paramSet[0]==1) omz = psr[p].param[param_om].val[0]*M_PI/180.0;
  else omz = 0.0;

  if (psr[p].param[param_a1dot].paramSet[0]==1) xdot  = psr[p].param[param_a1dot].val[0];
  else xdot  = 0.0;

  if (psr[p].param[param_pbdot].paramSet[0] == 1) pbdot = psr[p].param[param_pbdot].val[0];
  else pbdot = 0.0;

  if (psr[p].param[param_edot].paramSet[0] == 1) edot = psr[p].param[param_edot].val[0];
  else edot = 0.0;

  if (psr[p].param[param_xpbdot].paramSet[0] == 1) xpbdot = psr[p].param[param_xpbdot].val[0];
  else xpbdot = 0.0;


  asi = psr[p].param[param_a1].val[0]+xdot*tt0;
  ecc = psr[p].param[param_ecc].val[0]+edot*tt0;
  er = ecc*(1.0+dr);
  eth = ecc*(1.0+dth);

  /* Equation 10 in Kopeikin 1996 */
  ki_dot = -pmra * sin_omega + pmdec*cos_omega;
  ki    += ki_dot*tt0;
  sini = sin(ki);
  cosi = cos(ki);
  tani = sini/cosi;

  /* Equation 8 in Kopeikin 1996 */
  asi += (asi*ki_dot/tani)*tt0;
  /* Equation 9 in Kopeikin 1996 */
  omz += (pmra*cos_omega+pmdec*sin_omega)/sini*tt0;

  
  /* Now modify x and omega due to the annual-orbital parallax term 
   * as described in Kopeikin 1995 
   *
   * Require knowledge of the barycentric earth position vector - earth_ssb
   */
  
  /* Equation 15 in Kopeikin 1995 */
  delta_i0 = -psr[p].obsn[ipos].earth_ssb[0]/AULTSC*sin_alpha+
    psr[p].obsn[ipos].earth_ssb[1]/AULTSC*cos_alpha;
  /* Equation 16 in Kopeikin 1995 */
  delta_j0 = -psr[p].obsn[ipos].earth_ssb[0]/AULTSC*sin_delta*cos_alpha-
    psr[p].obsn[ipos].earth_ssb[1]/AULTSC*sin_delta*sin_alpha+
    psr[p].obsn[ipos].earth_ssb[2]/AULTSC*cos_delta;
  
  xpr = delta_i0*sin_omega - delta_j0*cos_omega;
  ypr = delta_i0*cos_omega + delta_j0*sin_omega;

  /* Equations 18 and 19 in Kopeikin 1995 */
  asi += asi/tani * dpara * xpr;
  si  += ci*dpara*xpr;
  omz -= 1.0/si*dpara*ypr;

  /* Continue as for the DD model now that asi and omz have been updated */
  orbits  = tt0/pb - 0.5*(pbdot+xpbdot)*pow(tt0/pb,2);
  norbits = (int)orbits;
  if (orbits<0.0) norbits--;
  phase=2.0*M_PI*(orbits-norbits);
  /*  Compute eccentric anomaly u by iterating Kepler's equation. */
  u=phase+ecc*sin(phase)*(1.0+ecc*cos(phase));
  do {
    du=(phase-(u-ecc*sin(u)))/(1.0-ecc*cos(u));
    u=u+du;
  } while (fabs(du)>1.0e-12);
  
  /*  DD equations 17b, 17c, 29, and 46 through 52 */
  su=sin(u);
  cu=cos(u);
  onemecu=1.0-ecc*cu;
  cae=(cu-ecc)/onemecu;
  sae=sqrt(1.0-pow(ecc,2))*su/onemecu;
  ae=atan2(sae,cae);
  if(ae<0.0) ae=ae+2.0*M_PI;
  ae=2.0*M_PI*orbits + ae - phase;
  omega=omz + k*ae;
  sw=sin(omega);
  cw=cos(omega);
  alpha=asi*sw;
  beta=asi*sqrt(1-pow(eth,2))*cw;
  bg=beta+gamma;
  dre=alpha*(cu-er) + bg*su;
  drep=-alpha*su + bg*cu;
  drepp=-alpha*cu - bg*su;
  anhat=an/onemecu;

  /* DD equations 26, 27, 57: */
  sqr1me2=sqrt(1-pow(ecc,2));
  cume=cu-ecc;
  brace=onemecu-si*(sw*cume+sqr1me2*cw*su);
  dlogbr=log(brace);
  ds=-2*m2*dlogbr;
  da=a0*(sin(omega+ae) + ecc*sw) + b0*(cos(omega+ae) + ecc*cw);

  /*  Now compute d2bar, the orbital time correction in DD equation 42. */
  d2bar=dre*(1-anhat*drep+(pow(anhat,2))*(pow(drep,2) + 0.5*dre*drepp -
					  0.5*ecc*su*dre*drep/onemecu)) + ds + da;
  torb=-d2bar;

  if (param==-1) return torb;
  
  /*  Now we need the partial derivatives. Use DD equations 62a - 62k. */
  csigma=asi*(-sw*su+sqr1me2*cw*cu)/onemecu;
  ce=su*csigma-asi*sw-ecc*asi*cw*su/sqr1me2;
  cx=sw*cume+sqr1me2*cw*su;
  comega=asi*(cw*cume-sqr1me2*sw*su);
  cgamma=su;
  cdth=-ecc*ecc*asi*cw*su/sqr1me2;
  cm2=-2*dlogbr;
  csi=2*m2*(sw*cume+sqr1me2*cw*su)/brace; 

  if (param==param_pb)
    return -csigma*an*SECDAY*tt0/(pb*SECDAY); 
  else if (param==param_a1)
    return cx;
  else if (param==param_ecc)
    return ce;
  else if (param==param_om)
    return comega;
  else if (param==param_omdot)
    return ae*comega/(an*360.0/(2.0*M_PI)*365.25*SECDAY);
  else if (param==param_t0)
    return -csigma*an*SECDAY;
  else if (param==param_pbdot)
    return 0.5*tt0*(-csigma*an*SECDAY*tt0/(pb*SECDAY));
  else if (param==param_sini)
    return csi;
  else if (param==param_gamma)
    return cgamma;
  else if (param==param_m2)
    return cm2*SUNMASS;
  else if (param==param_a1dot) /* Also known as xdot */
    return cx*tt0;

  /* Should calculate and return Kopeikin parameters */
  return 0;
}


void updateDDK(pulsar *psr,double val,double err,int pos)
{
  if (pos==param_pb)
    {
      psr->param[param_pb].val[0] += val/SECDAY;
      psr->param[param_pb].err[0]  = err/SECDAY;
    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_sini || pos==param_m2
	   || pos == param_gamma)
    {
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
    }
  else if (pos==param_om)
    {
      psr->param[pos].val[0] += val*180.0/M_PI;
      psr->param[pos].err[0]  = err*180.0/M_PI;
    }
  else if (pos==param_pbdot)
    {
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
    }
  else if (pos==param_a1dot)
    {
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
    }
  else if (pos==param_omdot)
    {
      psr->param[pos].val[0] += val; /* *(SECDAY*365.25)*180.0/M_PI; */
      psr->param[pos].err[0]  = err; /* *(SECDAY*365.25)*180.0/M_PI; */
    }
}
