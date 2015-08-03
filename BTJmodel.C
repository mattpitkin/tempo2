#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russel Edwards

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

/* 
 * This is the same as the BT model, but modified to include step changes in orbital 
 * parameters
 *
 */

double BTJmodel(pulsar *psr,int p,int ipos,int param,int arr)
{
  double torb;
  double tt0;
  double orbits;
  double pb;     /* Orbital period (sec) */
  double pbdot;
  double xpbdot;
  double ecc;    /* Orbital eccentricity */
  double edot;
  double asini;
  double xdot;
  double omdot;
  double omega;
  double gamma;
  int    norbits,i;
  double phase;
  double ep,dep,bige,tt,som,com;
  double alpha,beta,sbe,cbe,q,r,s,W;
  const char *CVS_verNum = "$Revision: 1.4 $";

  if (displayCVSversion == 1) CVSdisplayVersion("BTJmodel.C","BTJmodel()",CVS_verNum);

  tt0 = (psr[p].obsn[ipos].bbat - psr[p].param[param_t0].val[0])*SECDAY;

  pb     = psr[p].param[param_pb].val[0] * SECDAY;
  edot   = 0.0;
  ecc    = psr[p].param[param_ecc].val[0] + edot*tt0;
  if (psr[p].param[param_pbdot].paramSet[0] == 1) pbdot  = psr[p].param[param_pbdot].val[0];
  else pbdot=0.0;
  xpbdot = 0.0;
  if (psr[p].param[param_a1dot].paramSet[0] == 1) xdot = psr[p].param[param_a1dot].val[0];
  else xdot = 0.0;
  asini  = psr[p].param[param_a1].val[0] + xdot*tt0;
  if (psr[p].param[param_omdot].paramSet[0] == 1) omdot = psr[p].param[param_omdot].val[0];
  else omdot  = 0.0;
  omega  = (psr[p].param[param_om].val[0] + omdot*tt0/(SECDAY * 365.25))/(180.0/M_PI);
  if (psr[p].param[param_gamma].paramSet[0]==1) gamma = psr[p].param[param_gamma].val[0];
  else gamma  = 0.0;

  torb = 0.0;

  /* Now add in the jumps */
  for (i=0;i<psr[p].param[param_bpjep].aSize;i++)
    {
      if (psr[p].param[param_bpjep].paramSet[i]==1 && 
	  psr[p].obsn[ipos].bbat > psr[p].param[param_bpjep].val[i])
	{
	  torb -= (double)(psr[p].param[param_bpjph].val[i]/psr[p].param[param_f].val[0]);  
	  asini += (double)psr[p].param[param_bpja1].val[i];
	  ecc += (double)psr[p].param[param_bpjec].val[i];
	  omega += (double)psr[p].param[param_bpjom].val[i]; /* Check units */
	  pb += (double)psr[p].param[param_bpjpb].val[i];    /* Check units */ 
	}
    }
  /* Should ct be the barycentric arrival time? -- YES */
  orbits = tt0/pb - 0.5*(pbdot+xpbdot)*pow(tt0/pb,2); 
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



  /* SOME CODE HERE IN BNRYBT about nbin -- what is this ?? */
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
  torb = -q+(2*M_PI/pb)*q*r*s + torb;

  if (param==-1) return torb;

  W = asini*(sin(omega)*sbe - sqrt(1.0-ecc*ecc)*cos(omega)*cbe)/(1.0-ecc*cbe);

  if (param==param_pb)
    return -2.0*M_PI*r*s/pb*SECDAY*tt0/(SECDAY*pb) * SECDAY;  /* fctn(12+j) */
  else if (param==param_a1)
    return (som*(cbe-ecc) + com*sbe*sqrt(tt));                /* fctn(9+j) */
  else if (param==param_ecc)
    return -(alpha*(1.0+sbe*sbe-ecc*cbe)*tt - beta*(cbe-ecc)*sbe)*s/tt; /* fctn(10+j) */
  else if (param==param_om)
    return asini*(com*(cbe-ecc) - som*sqrt(tt)*sbe);          /* fctn(13+j) */
  else if (param==param_t0)
    return -2.0*M_PI/pb*r*s*SECDAY;                           /* fctn(11+j) */
  else if (param==param_pbdot)
    return 0.5*(-2.0*M_PI*r*s/pb*SECDAY*tt0/(SECDAY*pb))*tt0; /* fctn(18+j) */
  else if (param==param_a1dot)
    return (som*(cbe-ecc) + com*sbe*sqrt(tt))*tt0;            /* fctn(24+j) */
  else if (param==param_omdot)
    return asini*(com*(cbe-ecc) - som*sqrt(tt)*sbe)*tt0;      /* fctn(14+j) */
  else if (param==param_edot)                            
    return (-(alpha*(1.0+sbe*sbe-ecc*cbe)*tt - beta*(cbe-ecc)*sbe)*s/tt)*tt0; /* fctn(25+j) */
  else if (param==param_gamma) 
    return sbe;                                               /* fctn(15+j) */

  if (psr[p].param[param_bpjep].paramSet[arr]==1 && 
      psr[p].obsn[ipos].bbat > psr[p].param[param_bpjep].val[arr])
    {
      if (param==param_bpjph)      return 1.0/psr[p].param[param_f].val[0];
      else if (param==param_bpja1) return (som*(cbe-ecc) + com*sqrt(tt)*sbe);
      else if (param==param_bpjec) return -(alpha*(1.0+sbe*sbe-ecc*cbe)*tt - beta*(cbe-ecc)*sbe)*s/tt;
      else if (param==param_bpjom) return asini*(com*(cbe-ecc) - som*sqrt(tt)*sbe); 
      else if (param==param_bpjpb) return -2.0*M_PI*r*s/pb*SECDAY*tt0/(SECDAY*pb) * SECDAY;
    }
  return 0.0;
}

void updateBTJ(pulsar *psr,double val,double err,int pos,int arr)
{
 if (pos==param_pb || pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_gamma || pos==param_edot
     || pos==param_bpjph
     || pos==param_bpja1
     || pos==param_bpjec
     || pos==param_bpjom
     || pos==param_bpjpb
     )
    {
      psr->param[pos].val[arr] += val;
      psr->param[pos].err[arr]  = err;
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
