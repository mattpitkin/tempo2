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
#include <math.h>
#include "tempo2.h"

/* Based on bnrybtx.f */

double BTXmodel(pulsar *psr,int p,int ipos,int param,int k)
{
  double torb;
  double tt0;
  double orbits;
  int i,j;
  //  double pb;     /* Orbital period (sec) */
  //  double pbdot;
  //  double xpbdot;
  double ecc;    /* Orbital eccentricity */
  double edot;
  double asini;
  double xdot;
  double omdot;
  double omega;
  double gamma;
  int    norbits;
  double phase;
  double ep,dep,bige,tt,som,com;
  double alpha,beta,sbe,cbe,q,r,s,fac;
  const char *CVS_verNum = "$Revision: 1.2 $";

  if (displayCVSversion == 1) CVSdisplayVersion("BTXmodel.C","BTXmodel()",CVS_verNum);

  tt0 = (psr[p].obsn[ipos].bbat - psr[p].param[param_t0].val[0])*SECDAY;

  //  pb     = psr[p].param[param_pb].val[0] * SECDAY;
  edot   = 0.0;
  ecc    = psr[p].param[param_ecc].val[0] + edot*tt0;

  if (ecc < 0.0 || ecc > 1.0)
    {
      printf("BTXmodel: problem with eccentricity = %Lg\n",psr[p].param[param_ecc].val[0]);
      exit(1);
    }

  //  if (psr[p].param[param_pbdot].paramSet[0] == 1) pbdot  = psr[p].param[param_pbdot].val[0];
  //  else pbdot=0.0;
  //  xpbdot = 0.0;
  if (psr[p].param[param_a1dot].paramSet[0] == 1) xdot = psr[p].param[param_a1dot].val[0];
  else xdot = 0.0;
  asini  = psr[p].param[param_a1].val[0] + xdot*tt0;
  if (psr[p].param[param_omdot].paramSet[0] == 1) omdot = psr[p].param[param_omdot].val[0];
  else omdot  = 0.0;
  omega  = (psr[p].param[param_om].val[0] + omdot*tt0/(SECDAY * 365.25))/(180.0/M_PI);
  if (psr[p].param[param_gamma].paramSet[0]==1) gamma = psr[p].param[param_gamma].val[0];
  else gamma  = 0.0;

  torb = 0.0;

  // Note for the BTX model we use the FB parameters
  orbits = psr[p].param[param_fb].val[0]*tt0;
  // DO CHECK FOR NFBJ
  //  printf("Not checking for NFBJ %g\n",orbits);
  fac = 1.0;
  for (j=1;j<psr[p].param[param_fb].aSize;j++)
    {
      fac = fac/(double)(j+1);
      if (psr[p].param[param_fb].paramSet[j]==1)
	orbits += fac*psr[p].param[param_fb].val[j]*powl(tt0,j+1);
    }
  // UPDATE FOR EDOT, XDOT etc.
  //  printf("Not updating EDOT, XDOT, OMEGA etc. %g\n",orbits);

  //  orbits = tt0/pb - 0.5*(pbdot+xpbdot)*pow(tt0/pb,2); 
  norbits = (int)orbits;
  if (orbits < 0.0) norbits--;
  
  phase = 2.0*M_PI * (orbits-norbits);
 
  /* Using Pat Wallace's method of solving Kepler's equation -- code based on bnrybt.f */
  ep = phase + ecc*sin(phase)*(1.0+ecc*cos(phase));

  /* This line is wrong in the original tempo: should be inside the do loop */
  /*  denom = 1.0 - ecc*cos(ep);*/
  
  do {
    dep = (phase - (ep-ecc*sin(ep)))/(1.0 - ecc*cos(ep));
    ep += dep;
  } while (fabs(dep) > 1.0e-12);
  bige = ep;

  tt = 1.0-ecc*ecc;
  som = sin(omega);
  com = cos(omega);

  alpha = asini*som;
  beta = asini*com*sqrt(tt);
  sbe = sin(bige);
  cbe = cos(bige);
  q = alpha * (cbe-ecc) + (beta+gamma)*sbe;
  r = -alpha*sbe + beta*cbe;
  s = 1.0/(1.0-ecc*cbe);
  //  printf("Have torb = %g %g %g %g\n",torb,q,r,s);
  torb = -q+(2*M_PI*(double)psr[p].param[param_fb].val[0])*q*r*s + torb;
  //  printf("Have torb2 = %g %g %g %g %g\n",torb,q,r,s,(double)psr[p].param[param_fb].val[0]);
  if (param==-1) return torb;
  //  if (param==param_pb)
  //    return -2.0*M_PI*r*s/pb*SECDAY*tt0/(SECDAY*pb) * SECDAY;  /* fctn(12+j) */
  if (param==param_fb)
    {
      double diff;
      tt0/=1.0e7;
      diff = 2.0*M_PI*r*s*tt0;
      for (j=1;j<=k;j++)
	{
	  //	  tt0/=1.0e7;
	  diff = (1.0/(double)(j+1))*tt0*diff;      
	}
      //      printf("Returning: param_fb %d %g %g %g %g\n",k,diff,(double)r,(double)s,(double)tt0);
      return diff;
    }
  else if (param==param_a1)
    {
      //      printf("Returning A1: %g\n",(som*(cbe-ecc) + com*sbe*sqrt(tt)));
      return (som*(cbe-ecc) + com*sbe*sqrt(tt));                /* fctn(9+j) */
    }
  else if (param==param_ecc)
    return -(alpha*(1.0+sbe*sbe-ecc*cbe)*tt - beta*(cbe-ecc)*sbe)*s/tt; /* fctn(10+j) */
  else if (param==param_om)
   return asini*(com*(cbe-ecc) - som*sqrt(tt)*sbe);          /* fctn(13+j) */
  else if (param==param_t0)
    {
      //      printf("Returning T0 %g\n",-2.0*M_PI*psr[p].param[param_fb].val[0]*r*s*SECDAY);
      return -2.0*M_PI*psr[p].param[param_fb].val[0]*r*s*SECDAY;                           /* fctn(11+j) */
    }
  //  else if (param==param_pbdot)
  //    return 0.5*(-2.0*M_PI*r*s/pb*SECDAY*tt0/(SECDAY*pb))*tt0; /* fctn(18+j) */
  else if (param==param_a1dot)
    return (som*(cbe-ecc) + com*sbe*sqrt(tt))*tt0;            /* fctn(24+j) */
  else if (param==param_omdot)
    return asini*(com*(cbe-ecc) - som*sqrt(tt)*sbe)*tt0;      /* fctn(14+j) */
  else if (param==param_edot)                            
    return (-(alpha*(1.0+sbe*sbe-ecc*cbe)*tt - beta*(cbe-ecc)*sbe)*s/tt)*tt0; /* fctn(25+j) */
  else if (param==param_gamma) 
    return sbe;                                               /* fctn(15+j) */
  return 0.0;
}

void updateBTX(pulsar *psr,double val,double err,int pos,int k)
{
  if (pos==param_fb)
    {
      psr->param[param_fb].val[k] += (val/powl(1.0e7,k+1));
      psr->param[param_fb].err[k]  = err/powl(1.0e7,k+1);
    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_gamma || pos==param_edot)
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
  else if (pos==param_omdot)
    {
      psr->param[pos].val[0] += val*(SECDAY*365.25)*180.0/M_PI;
      psr->param[pos].err[0]  = err*(SECDAY*365.25)*180.0/M_PI;
    }
  else if (pos==param_a1dot)
    {
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
    }
}
