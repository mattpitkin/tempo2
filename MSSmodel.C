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

/* Based on bnrymss.f */
/* Model for main-sequence star binary pulsars (Wex 1998, astro-ph/9706086).
  -> changes in x proportonal to Ae(u)
  -> use of second time derivative in omega and x: om2dot, x2dot

  Computes pulsar orbit time, torb, at time of observation t=ct(n)-pepoch.
  Pulsar proper time is then TP=T+TORB.
  Units are such that c=g=1. Thus masses have units of seconds, with
  one solar mass = 4.925490947 usec.

  Also computes the binary orbit-related values of fctn: partial
  derivatives of each arrival time residual with respect to the model
  parameters. */

double MSSmodel(pulsar *psr,int p,int obs,int param)
{
  double SUNMASS=4.925490947e-6;
  double RAD = 180.0/M_PI;
  double an,ecc0,x0,omega0,k,xi,m2,orbits,phase,ecc,er,eth,du,om2dot,x2dot,si,b0,a0;
  double u,su,cu,onemecu,cae,sae,ae,x,sw,cw,alpha,beta,bg,dre,drep,drepp,anhat;
  double sqr1me2,cume,brace,dlogbr,ds,da,d2bar,torb,edot,dr,dth,gamma;
  double csigma,ce,cx,comega,cgamma,cdth,cm2,csi;
  int norbits;  
  double pb,eccentricity,a1,omega,omdot,xdot,am2,pbdot,tt0;
  double shapmax;
  double e2dot;
  double orbpx;
  const char *CVS_verNum = "$Revision: 1.7 $";

  if (displayCVSversion == 1) CVSdisplayVersion("MSSmodel.C","MSSmodel()",CVS_verNum);

  tt0 = (psr[p].obsn[obs].bbat - psr[p].param[param_t0].val[0])*SECDAY;

  if (psr[p].param[param_pbdot].paramSet[0] == 1) pbdot  = psr[p].param[param_pbdot].val[0];
  else pbdot = 0.0;

  if (psr[p].param[param_m2].paramSet[0] == 1) am2  = psr[p].param[param_m2].val[0];
  else am2 = 0.0;

  if (psr[p].param[param_a1dot].paramSet[0] == 1) xdot = psr[p].param[param_a1dot].val[0];
  else xdot = 0.0;
  if (psr[p].param[param_a2dot].paramSet[0] == 1) x2dot = psr[p].param[param_a2dot].val[0];
  else x2dot = 0.0;

  if (psr[p].param[param_e2dot].paramSet[0] == 1) e2dot = 1e-20*psr[p].param[param_e2dot].val[0];
  else e2dot = 0.0;
    
  if(psr[p].param[param_orbpx].paramSet[0] ==1) orbpx = 1./3.086e21*psr[p].param[param_orbpx].val[0];
  else orbpx = 0.0;


  //fprintf(stderr, "AAAAAA %.3e\n", x2dot);
  //exit(0);

  if (psr[p].param[param_dtheta].paramSet[0] == 1) dth = psr[p].param[param_dtheta].val[0];
  else dth=0.0;

  if (psr[p].param[param_sini].paramSet[0] == 1) si = psr[p].param[param_sini].val[0];
  else si=0.0;

  if (psr[p].param[param_a0].paramSet[0] == 1) a0 = psr[p].param[param_a0].val[0];
  else a0=0.0;

  if (psr[p].param[param_b0].paramSet[0] == 1) b0 = psr[p].param[param_b0].val[0];
  else b0=0.0;

  if (psr[p].param[param_dr].paramSet[0] == 1) dr = psr[p].param[param_dr].val[0];
  else dr=0.0;

  if (psr[p].param[param_omdot].paramSet[0]==1) omdot = psr[p].param[param_omdot].val[0];
  else omdot  = 0.0;

  if (psr[p].param[param_gamma].paramSet[0]==1) gamma = psr[p].param[param_gamma].val[0];
  else gamma  = 0.0;

  if (psr[p].param[param_om2dot].paramSet[0]==1) om2dot = psr[p].param[param_om2dot].val[0];
  else om2dot  = 0.0;

  if (psr[p].param[param_edot].paramSet[0]==1) edot = psr[p].param[param_edot].val[0];
  else edot  = 0.0;

  if (psr[p].param[param_shapmax].paramSet[0]==1) shapmax = psr[p].param[param_shapmax].val[0];
   else shapmax  = 0.0;


  pb  = psr[p].param[param_pb].val[0] * SECDAY;
  eccentricity    = psr[p].param[param_ecc].val[0]; 
  a1  = psr[p].param[param_a1].val[0];
  omega  = (psr[p].param[param_om].val[0]);


  an=2.0*M_PI/pb;
  ecc0=eccentricity;
  x0=a1;
  omega0=omega/RAD;
  k=omdot/an/(RAD*365.25*SECDAY);
  xi=xdot/an;
  m2=am2*SUNMASS;

  tt0 = ((double)psr[p].obsn[obs].bbat - (double)psr[p].param[param_t0].val[0])*SECDAY;
  orbits=tt0/pb - 0.5*pbdot*pow(tt0/pb,2);
  norbits=(int)orbits;
  if(orbits<0.0) norbits=norbits-1;
  phase=2.0*M_PI*(orbits-norbits);
  ecc=ecc0 + edot*tt0  +0.5*e2dot*tt0*tt0;
  er =ecc*(1.0+dr);
  eth=ecc*(1.0+dth);
  /*  Compute eccentric anomaly u by iterating Kepler's equation.*/
  u=phase+ecc*sin(phase)*(1+ecc*cos(phase));
  //  printf("params: %g %g %g %d\n",(double)tt0,(double)si,(double)am2,psr[p].param[param_om].paramSet[2]);

  do
    {
      du=(phase-(u-ecc*sin(u)))/(1.0-ecc*cos(u));
      u=u+du;
    } while (fabs(du)>1.0e-14);

  /* DD equations 17b, 17c, 29, and 46 through 52 */
  su=sin(u);
  cu=cos(u);
  onemecu=1.0-ecc*cu;
  cae=(cu-ecc)/onemecu;
  sae=sqrt(1.0-pow(ecc,2))*su/onemecu;
  ae=atan2(sae,cae);
  if(ae < 0.0) ae=ae+2.0*M_PI;
  ae=2.0*M_PI*orbits + ae - phase;

  //  omega = omega0 +  k*ae + 0.5*om2dot*pow(tt0,2); /* Wex 1998 */
  //x = x0+xi*ae+0.5*x2dot*pow(tt0,2);

  double xii =  1e-20*x2dot/an/an;
  double oii =  1e-20*om2dot/an/an; ///RAD/pow(365.25*SECDAY,2.);

  x = x0+xi*ae+0.5*xii*pow(ae,2); 
  omega = omega0+k*ae +0.5*oii*pow(ae,2.);

  //
  //omega = omega0 +  k*ae + 0.5*om2dot*tt0*tt0; /* Wex 1998 */
  //x = x0+xi*ae+0.5*x2dot*tt0*tt0;
  



  sw=sin(omega);
  cw=cos(omega);
  alpha=x*sw;
  beta=x*sqrt(1-pow(eth,2))*cw;
  bg=beta+gamma;
  dre=alpha*(cu-er) + bg*su;
  drep=-alpha*su + bg*cu;
  drepp=-alpha*cu - bg*su;
  anhat=an/onemecu;

  /*  DD equations 26, 27, 57 */
  sqr1me2=sqrt(1-pow(ecc,2));
  cume=cu-ecc;
  brace=onemecu-si*(sw*cume+sqr1me2*cw*su);
  dlogbr=log(brace);
  ds=-2.0*m2*dlogbr;
  da=a0*(sin(omega+ae) + ecc*sw) + b0*(cos(omega+ae) + ecc*cw);
  
  /* Now compute d2bar, the orbital time correction in DD equation 42 */
  d2bar=dre*(1-anhat*drep+(pow(anhat,2))*(pow(drep,2) + 0.5*dre*drepp -
					  +    0.5*ecc*su*dre*drep/onemecu)) + ds + da;
  torb=-d2bar;

  double shapparam;
  shapparam = -log(1-ecc*cu-(sin(omega)*(cu-ecc)+ sqrt(1-ecc*ecc)*cos(omega)*su)*si);

  torb -= shapmax*shapparam;


 

  // add in orbital parallax

 
  
  double so = sin(omega);
  double co  = cos(omega);
  double cpx; 

  

  cpx= 1e2*SPEED_LIGHT*x*x/2.*(  1./si/si-0.5+0.5*ecc*ecc*(1 + so*so -3/si/si) 
			       -2*ecc*(1./si/si- so*so)*(cu  -ecc)
			       + sqrt(1-ecc*ecc)*sin(2*omega)*(ecc*su-0.5*sin(2*u))
			       +0.5*( cos(2*omega) +  ecc*ecc*(1./si/si+ co*co))*cos(2*u));

  //csi = 1e2*SPEED_LIGHT*x*x/2.*orbpx*( -2*ecc*(-2./si/si/si)*cu + 0.5*ecc*ecc*(-2./si/si/si)*cos(2*u));

  //cpx = 1e2*SPEED_LIGHT*x0*x0/2.*(0.5*(cos(2*omega) + ecc*ecc*(1./si/si+ co*co))*cos(2*u));
 

  if( orbpx != 0)
    {
      
      torb -= cpx*orbpx;
    }      
  
 
  
  
  /*  printf("MSS here: %.20g %.20g %.20g %.20g\n",dlogbr,ds,da,torb);*/

  if (param==-1) return torb;

  /*  Partial derivatives, DD equations 62a - 62k*/
  csigma=x*(-sw*su+sqr1me2*cw*cu)/onemecu;
  ce=su*csigma-x*sw-ecc*x*cw*su/sqr1me2;
  cx=sw*cume+sqr1me2*cw*su;
  comega=x*(cw*cume-sqr1me2*sw*su);
  cgamma=su;
  cdth=-ecc*ecc*x*cw*su/sqr1me2;
  cm2=-2*dlogbr;
  //csi=2*m2*(sw*cume+sqr1me2*cw*su)/brace;


  /* Otherwise here for fitting */
  if (param==param_a1) return cx;
  else if (param==param_ecc) return ce;
  else if (param==param_om) return comega;
  else if (param==param_omdot) return comega*ae/an;
  else if (param==param_pb) return -csigma*an*SECDAY*tt0/(pb*SECDAY);
  else if (param==param_t0) return -csigma*an*SECDAY;  
  else if (param==param_a1dot) return cx*ae/an;
  else if (param==param_pbdot) return  0.5*tt0*(-csigma*an*SECDAY*tt0/(pb*SECDAY));
  else if (param==param_a2dot) return 1e-20*0.5*cx*ae/an*ae/an;   
  else if (param==param_om2dot) return  1e-20*0.5*comega*ae/an*ae/an;    
  else if (param==param_shapmax) return shapparam;
  else if (param==param_m2) return cm2*SUNMASS;

  else if (param==param_sini) return csi;
  else if (param==param_edot) return ce*tt0;
  else if (param==param_e2dot) return 1e-20*0.5*ce*tt0*tt0;
  else if (param==param_orbpx) return 1./3.086e21*cpx;
  return 0.0;
}

void updateMSS(pulsar *psr,double val,double err,int pos)
{
  
 
  if (pos==param_pb)
    {
      psr->param[param_pb].val[0] += val/SECDAY;
      psr->param[param_pb].err[0]  = err/SECDAY;
    }
    else if (pos==param_a1 || pos==param_ecc || pos==param_t0)
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
      psr->param[pos].val[0] += val*180/M_PI*365.25*SECDAY;
      psr->param[pos].err[0]  = err*180/M_PI*365.25*SECDAY;
    }
  else if (pos==param_a1dot)
    {
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
    }
   else if (pos==param_a2dot)
    {
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
    }
   else if (pos==param_om2dot)
    {
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
    }
   else if(pos==param_shapmax)
     {
       psr->param[pos].val[0] += val;
       psr->param[pos].err[0]  = err;
     }
   else if (pos==param_m2)
     {
       psr->param[pos].val[0] += val;
       psr->param[pos].err[0] = err;
     }

  else if (pos==param_sini)
    {
     psr->param[pos].val[0] += val;
    psr->param[pos].err[0]  = err;
    }
    else if (pos==param_edot)
    {
     psr->param[pos].val[0] += val;
    psr->param[pos].err[0]  = err;
    }
   else if (pos==param_e2dot)
    {
     psr->param[pos].val[0] += val;
    psr->param[pos].err[0]  = err;
    }
else if (pos==param_orbpx)
    {
     psr->param[pos].val[0] += val;
    psr->param[pos].err[0]  = err;
    }
  


}
