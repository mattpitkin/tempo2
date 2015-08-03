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


/* Timing model      */
/* Based on bnryddgr.f */

/* Damour & Deruelle model assuming the general theory of relativity */
/* Computes the pulsar orbit time, torb, at the time of observation, t =ct(n)-pepoch */
/* Pulsar proper time is TP = T + TORB */
/* Units are c=g=1. */

void mass2dd(double am,double am2,double x,double ecc,double an,double *arr,double *ar,
	     double *xk,double *si,double *gamma,double *pbdot);


double DDGRmodel(pulsar *psr,int p,int ipos,int param)
{
  double an,afac;
  double pb,k;
  double rad2deg = 180.0/M_PI;
  double SUNMASS = 4.925490947e-6;
  double m2,tt0,t0,x,ecc,er,xdot,edot,dr,dth,eth,am2,ct,f0,cm;
  double pbdot,xpbdot,phase,u,du,gamma,m,m1,arr,ar,xk,fac,xomdot;
  double dtdm2,dgmdm2,dsidm2,dkdm2,dthdm2,dpbdm2,darrdm2,dnum,denom;
  double orbits,a0aligned;
  int norbits;
  double  cu,onemecu,ae,omega,omz,sw,cw,alpha,beta,bg,dre,drep,drepp,anhat,su;
  double sqr1me2,cume,brace,si,dlogbr,ds,da,a0,b0,d2bar,torb,am;
  double csigma,ce,cx,comega,cgamma,cdth,cm2,csi;
  double fact1,fact2,fact3,fact4,fact5,fact6,fact7,fact8,denumm,denomm,darrdm,ck,dkdm,cdr;
  double ddrdm,cpbdot,dpbdm,csini,dsidm,an0,dgamdm,dthdm,ddrdm2;
  const char *CVS_verNum = "$Revision: 1.5 $";

  if (displayCVSversion == 1) CVSdisplayVersion("DDGRmodel.C","DDGRmodel()",CVS_verNum);

  t0 = psr[p].param[param_t0].val[0];
  ct = psr[p].obsn[ipos].bbat;    
  
  tt0 = (ct-t0)*SECDAY;

  f0 = psr[p].param[param_f].val[0];

  xomdot = 0.0;  /* WHAT SHOULD THIS BE??? */
  afac = 0.0;    /* WHAT SHOULD THIS BE??? */

  if (psr[p].param[param_sini].paramSet[0]==1) si = getParameterValue(&psr[p],param_sini,0);
  else si = 0.0;

  if (psr[p].param[param_m2].paramSet[0]==1) am2 = psr[p].param[param_m2].val[0];
  else am2 = 0.0;

  pb = psr[p].param[param_pb].val[0]*SECDAY;
  an = 2.0*M_PI/pb;

  if (psr[p].param[param_mtot].paramSet[0]==1) am = psr[p].param[param_mtot].val[0];
  else am = 0.0;

  m  = am*SUNMASS;
  m2 = am2*SUNMASS;
  m1 = m-m2;

  if (psr[p].param[param_om].paramSet[0]==1) omz = psr[p].param[param_om].val[0];
  else omz = 0.0;

  if (psr[p].param[param_a1dot].paramSet[0]==1) xdot  = psr[p].param[param_a1dot].val[0];
  else xdot  = 0.0;

  if (psr[p].param[param_pbdot].paramSet[0] == 1) pbdot = psr[p].param[param_pbdot].val[0];
  else pbdot = 0.0;

  if (psr[p].param[param_edot].paramSet[0] == 1) edot = psr[p].param[param_edot].val[0];
  else edot = 0.0;

  if (psr[p].param[param_xpbdot].paramSet[0] == 1) xpbdot = psr[p].param[param_xpbdot].val[0];
  else xpbdot = 0.0;
  
  x = psr[p].param[param_a1].val[0]+xdot*tt0;
  ecc = psr[p].param[param_ecc].val[0]+edot*tt0;

  /* Given system masses m,m2 and keplerian parameters x,ecc,an, calculate the values
   * of arr,ar,si,gamma,pbdot under general relativity */

  mass2dd(am,am2,x,ecc,an,&arr,&ar,&xk,&si,&gamma,&pbdot);

  k=xk;

  dr = (3.0*pow(m1,2) + 6.0*m1*m2 + 2.0*pow(m2,2))/(arr*m);
  er = ecc*(1.0+dr);
  dth = (3.5*m1*m1 + 6*m1*m2 + 2*m2*m2)/(arr*m);
  eth = ecc*(1.0+dth);

  orbits = tt0/pb - 0.5*(pbdot+xpbdot)*pow(tt0/pb,2);
  //  printf("xpbdot = %.14g %.14g\n",(double)xpbdot/1e-12,(double)orbits);
  norbits = (int)orbits;
  if (orbits<0.0) norbits--;
  phase=2.0*M_PI*(orbits-norbits);
  /*  Compute eccentric anomaly u by iterating Kepler's equation. */
  u=phase+ecc*sin(phase)*(1.0+ecc*cos(phase));
  fac = 1.0/(1.0-ecc*cos(u));  /* NOTE COULD BE WRONG IN DDmodel - SEE USE OF FAC !!!! */
  do {
    du=(phase-(u-ecc*sin(u)))*fac; 
    u=u+du;
  } while (fabs(du)>1.0e-14);  /* 1e-12 in DDmodel */
  
  /*  DD equations 17a, 29 */
  ae = 2.0*atan(sqrt((1+ecc)/(1-ecc))*tan(0.5*u));
  if(ae<0.0) ae=ae+2.0*M_PI;
  ae = 2.0*M_PI*orbits + ae-phase;
  omega=omz/rad2deg + (k+xomdot/(an*rad2deg*365.25*86400.0))*ae;
  /* DD equations 46 through 52 */
  su=sin(u);
  cu=cos(u);
  sw=sin(omega);
  cw=cos(omega);
  alpha=x*sw;
  beta=x*sqrt(1-pow(eth,2))*cw;
  bg=beta+gamma;
  dre=alpha*(cu-er) + bg*su;
  drep=-alpha*su + bg*cu;
  drepp=-alpha*cu - bg*su;
  onemecu=1.0-ecc*cu;
  anhat=an/onemecu;

  /* DD equations 26,27,57 */

  cume=cu-ecc;
  sqr1me2=sqrt(1-pow(ecc,2));
  brace=onemecu-si*(sw*cume+sqr1me2*cw*su);
  if (brace<=0)
    {
      printf("ERROR: In DDGR model, brace < 0\n");
      exit(1);
    }
  dlogbr=log(brace);
  ds=-2*m2*dlogbr;

  /* These will be different if spin axis not aligned -- IS THIS AN ASSUMPTION OF THE MODEL? */
  a0aligned = an*ar/(2.0*M_PI*f0*si*sqr1me2);
  a0 = afac*a0aligned;
  b0 = 0.0;
  da = a0*(sin(omega+ae)+ecc*sw) + b0*(cos(omega+ae) + ecc*cw);


  /*  Now compute d2bar, the orbital time correction in DD equation 42. */
  d2bar=dre*(1-anhat*drep+(pow(anhat,2))*(pow(drep,2) + 0.5*dre*drepp -
					  +    0.5*ecc*su*dre*drep/onemecu)) + ds + da;
  torb=-d2bar;

  if (param==-1)  return torb;
  
  /* Now get partial derivatives */
  an0 = sqrt(m/pow(arr,3));
  
  csigma=x*(-sw*su+sqr1me2*cw*cu)/onemecu;
  ce=su*csigma-x*sw-ecc*x*cw*su/sqr1me2;
  cx=sw*cume+sqr1me2*cw*su;
  comega=x*(cw*cume-sqr1me2*sw*su);
  cgamma=su;
  cm2=-2*dlogbr;
  csi=2*m2*(sw*cume+sqr1me2*cw*su)/brace; 

  fact1=(m/(2*arr)) * ((m-m2)*m2/pow(m,2) - 9);
  fact2=(3*m/(2*pow(arr,4))) * (1.0 + fact1);
  fact3=(m/2*arr) * (m2/pow(m,2)-2*(m-m2)*m2/pow(m,3));
  fact4=(1+fact1)*3*m/(2*pow(arr,4)*an0);
  fact5=an0*fact1/arr;

  denumm=(1+fact1)/(2*pow(arr,3)*an0) + an0*(fact1/m+fact3);
  denomm=fact4+fact5;
  darrdm=denumm/denomm;
  
  dnum = an0*(m-2*m2)/(2*arr*m);
  denom = an0*fact1/arr + fact2/an0;

  darrdm2 = dnum/denom;

  dgmdm2 = ((m+2*m2)/arr - (m2*(m+m2)*darrdm2/pow(arr,2)))*ecc/(an*m);
  cdth=-ecc*ecc*x*cw*su/sqr1me2;
  dthdm2 = -dth*darrdm2/arr - (m+m2)/(arr*m);

  dkdm = k/m - k*darrdm/arr;
  dsidm2 = -(m*x/(arr*m2))*(1.0/m2+darrdm2/arr);
  ck = ae*comega;
  dkdm2 = -k*darrdm2/arr;
  cdr = -ecc*x*sw;
  ddrdm2 = -dr*darrdm2/arr - 2*m2/(arr*m);
  dtdm2 = -2*dlogbr;

  csini = 2*m2*(sw*cume+sqr1me2*cw*su)/brace;

  dsidm=-(m*x/(arr*m2))*(-1.0/m+darrdm/arr);
  dpbdm = pbdot/(m-m2) - pbdot/(3*m);
  cpbdot = -csigma*an*pow(tt0,2)/(2*pb);
  ddrdm = -dr/m - dr*darrdm/arr + 6/arr;

  dpbdm2 = pbdot/m2 - pbdot/(m-m2);

  cm2 = dtdm2+cgamma*dgmdm2+csini*dsidm2+ck*dkdm2+cdr*ddrdm2+cdth*dthdm2+cpbdot*dpbdm2;



  fact6=1.0/(arr*m);
  fact7=-(m+m2)/(arr*pow(m,2));
  fact8=-(m+m2)*darrdm/(pow(arr,2)*m);
  dgamdm = (ecc*m2/an)*(fact6+fact7+fact8);

  dthdm=-dth/m - dth*darrdm/arr + (7*m-m2)/(arr*m);

  cm = ck*dkdm+cgamma*dgamdm+cdr*ddrdm+cdth*dthdm+cpbdot*dpbdm+csini*dsidm;

  if (param==-2)  /* Set derived parameters */
    {
      /* calculated values (assuming GR) */
      psr[p].param[param_sini].paramSet[0]=1;
      psr[p].param[param_sini].val[0]=si;
      
      /* Should be xomdot??? */
      psr[p].param[param_omdot].paramSet[0]=1;
      psr[p].param[param_omdot].val[0]=360.0*365.25*xk/(pb/SECDAY);
      
      psr[p].param[param_gamma].paramSet[0]=1;
      psr[p].param[param_gamma].val[0]=gamma;

      psr[p].param[param_pbdot].paramSet[0]=1;
      psr[p].param[param_pbdot].val[0]=pbdot;
  
      psr[p].param[param_dtheta].paramSet[0]=1;
      psr[p].param[param_dtheta].val[0]=dth;

      psr[p].param[param_dr].paramSet[0]=1;
      psr[p].param[param_dr].val[0]=dr;

      return 0;

    }

  if (param==param_pb)
    return -csigma*an*SECDAY*tt0/(pb*SECDAY); 
  else if (param==param_a1)
    return cx;
  else if (param==param_ecc)
    return ce;
  else if (param==param_om)
    return comega;
  else if (param==param_t0)
    return -csigma*an*SECDAY;
  else if (param==param_pbdot)
    return 0.5*tt0*(-csigma*an*SECDAY*tt0/(pb*SECDAY));
  else if (param==param_xpbdot)
    return 0.5*tt0*(-csigma*an*SECDAY*tt0/(pb*SECDAY));
  else if (param==param_sini)
    return csi;
  else if (param==param_m2)
    return cm2*SUNMASS;
  else if (param==param_mtot)
      return cm*SUNMASS;
  else if (param==param_a1dot) /* Also known as xdot */
    return cx*tt0;

  return 0.0;
}

void updateDDGR(pulsar *psr,double val,double err,int pos)
{
  if (pos==param_pb)
    {
      psr->param[param_pb].val[0] += val/SECDAY;
      psr->param[param_pb].err[0]  = err/SECDAY;
    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_sini || pos==param_m2
	   || pos==param_mtot)
    {
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
    }
  else if (pos==param_om)
    {
      psr->param[pos].val[0] += val*180.0/M_PI;
      psr->param[pos].err[0]  = err*180.0/M_PI;
    }
  else if (pos==param_pbdot || pos==param_xpbdot)
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
      psr->param[pos].val[0] += val*(SECDAY*365.25)*180.0/M_PI;
      psr->param[pos].err[0]  = err*(SECDAY*365.25)*180.0/M_PI;
    }
}

/* Given system masses of m,m2 and Keplerian parameters x,ecc and an, this 
 * routine calculates values of arr, ar, si, gamma and pbdot under GR */

void mass2dd(double am,double am2,double x,double ecc,double an,double *arr,double *ar,
	     double *xk,double *si,double *gamma,double *pbdot)
{
  double SUNMASS = 4.925490947e-6;
  double ARRTOL = 1.0e-10;

  double m,m2,m1,arr0,arrold;

  m=am*SUNMASS;
  m2 = am2*SUNMASS;
  m1 = m-m2;
  if (m<0)
    {
      printf("ERROR: problem in DDGR model (mtot < 0)\n");
      exit(1);
    }
  arr0 = pow(m/(an*an),1.0/3.0);
  *arr = arr0;
  do {
    arrold = *arr;
    *arr = arr0*pow(1.0+(m1*m2/pow(m,2) - 9.0)*0.5*m/(*arr),2.0/3.0);
  } while (fabs(((*arr)-arrold)/(*arr)) > ARRTOL);


  *arr = arr0*pow(1.0+(m1*m2/pow(m,2) - 9.0)*0.5*m/(*arr),2.0/3.0);
  *ar = (*arr)*m2/m;


  *si=x/(*ar);
  *xk=3.0*m/((*arr)*(1.0-ecc*ecc));
  *gamma = ecc*m2*(m1+2*m2)/(an*(*arr)*m);
  *pbdot = -(96.0*2.0*M_PI/5.0)*pow(an,5.0/3.0)*pow(1.0-pow(ecc,2),-3.5)
    * (1+(73.0/24)*pow(ecc,2) + (37.0/96)*pow(ecc,4)) * m1*m2*pow(m,-1.0/3.0);
}
