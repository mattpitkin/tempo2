//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards
/* 
   This file implements the harmonic Shapiro delay fitting as
   described by Freire & Wex (2010). It is a modification of the
   Tempo2 file DDmodel.C, which in turn was based on bnrydd.f, which
   was part of the tempo pulsar timing package.

   DDHmodel.C was created by J.P.W. Verbiest at Parkes observatory in May 2011.
 */
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
#include <stdlib.h>

/* Timing model      */
/* Based on bnrydd.f */

double DDHmodel(pulsar *psr,int p,int ipos,int param){
  long double an; // angular velocity
  long double pb,k; // orbital period and unitless orbital precession
  // conversion factors
  long double rad2deg = 180.0/M_PI; 
  long double SUNMASS = 4.925490947e-6;
  // binary parameters
  long double m2,tt0,t0,x,ecc,er,xdot,edot,dr,dth,eth,am2,ct;
  // more obscure orbital parameters
  long double pbdot,xpbdot,phase,u,du,gamma;
  long double orbits;
  int norbits;
  // DDH specific orbital parameters:
  long double h3, stig;// 3rd, 4th harmonics and harmonic ratio
  // needed for intermediate steps in calculations
  long double  cu,onemecu,cae,sae,ae,omega,omz,sw,cw,alpha,beta,bg,dre,
    drep,drepp,anhat,su;
  long double sqr1me2,cume,brace,si,dlogbr,ds,da,a0,b0,d2bar,torb;
  long double csigma,ce,cx,comega,cgamma,cdth,cm2,csi;
  // Aberration parameters (JPWV thinks so, at least).
  dr = 0.0; /* WHAT SHOULD THESE BE SET TO? */
  dth = 0.0; 

  // Sanity check
  if( psr[p].param[param_h3].paramSet[0] * psr[p].param[param_stig].paramSet[0] != 1 ){
    printf("ERROR! For the DDH model, both h3 and stig must be set.\n Exiting.\n");
    exit( 1 );
  }
  h3 = getParameterValue( &psr[p], param_h3, 0 );
  stig = getParameterValue( &psr[p], param_stig, 0 );

  // Define sin(i) and m2 for calculation of the orbital phases etc.
  // FW10, Eq. 22:
  si = 2.0 * stig / ( 1.0 + pow( stig, 2.0 ) );
  
  // FW10, Eq. 20:
  m2 = h3 / pow( stig, 3.0 ); // Shapiro r, not just M2.

  if (si > 1.0){
    displayMsg(1,"BIN1",
               "SIN I > 1.0, setting to 1: should probably use DDS model",
               "",psr[p].noWarnings);
    si = 1.0;
    psr[p].param[param_sini].val[0] = 1.0L;
  }

  pb = psr[p].param[param_pb].val[0]*SECDAY;
  an = 2.0*M_PI/pb;
  // k is the orbital precession expressed in orbits per orbit
  k = psr[p].param[param_omdot].val[0]/(rad2deg*365.25*86400.0*an);
  // this equation is effectively the same as (omdot/360) * (pb/365.25)

  t0 = psr[p].param[param_t0].val[0];
  ct = psr[p].obsn[ipos].bbat;    
  
  tt0 = (ct-t0)*SECDAY;

  if (psr[p].param[param_gamma].paramSet[0]==1)
    gamma = psr[p].param[param_gamma].val[0];
  else
    gamma = 0.0;

  // JPWV thinks a0 and b0 are the aberration delays and are at
  // present not implemented. No way to be sure though, short of
  // asking George. Or, more likely, to spend more time to try and get
  // to the bottom of this. I cannot be bothered either way.
  a0    = 0.0;
  b0    = 0.0;

  if (psr[p].param[param_om].paramSet[0]==1) 
    omz = psr[p].param[param_om].val[0];
  else omz = 0.0;

  if (psr[p].param[param_a1dot].paramSet[0]==1) 
    xdot  = psr[p].param[param_a1dot].val[0];
  else xdot  = 0.0;

  if (psr[p].param[param_pbdot].paramSet[0] == 1) 
    pbdot = psr[p].param[param_pbdot].val[0];
  else pbdot = 0.0;

  if (psr[p].param[param_edot].paramSet[0] == 1) 
    edot = psr[p].param[param_edot].val[0];
  else edot = 0.0;

  if (psr[p].param[param_xpbdot].paramSet[0] == 1) 
    xpbdot = psr[p].param[param_xpbdot].val[0];
  else xpbdot = 0.0;


  x = psr[p].param[param_a1].val[0]+xdot*tt0;
  ecc = psr[p].param[param_ecc].val[0]+edot*tt0;
  er = ecc*(1.0+dr);
  eth = ecc*(1.0+dth);
  
  if (ecc < 0.0 || ecc > 1.0){
    printf("DDHmodel: problem with eccentricity = %Lg [%s]\n",
           psr[p].param[param_ecc].val[0],psr[p].name);
    exit(1);
  }
  
  orbits = tt0/pb - 0.5L*(pbdot+xpbdot)*(tt0/pb)*(tt0/pb);
  norbits = (int)orbits;
  if (orbits<0.0) norbits--;
  phase=2.0*M_PI*(orbits-norbits);

  //  Compute eccentric anomaly u by iterating Kepler's equation.
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
  omega=omz/rad2deg + k*ae;
  sw=sin(omega);
  cw=cos(omega);
  alpha=x*sw;
  beta=x*sqrt(1-pow(eth,2))*cw;
  bg=beta+gamma;
  dre=alpha*(cu-er) + bg*su;
  drep=-alpha*su + bg*cu;
  drepp=-alpha*cu - bg*su;
  anhat=an/onemecu;

  // DD equations 26, 27, 57:
  sqr1me2=sqrt(1-pow(ecc,2));
  cume=cu-ecc;
  brace=onemecu-si*(sw*cume+sqr1me2*cw*su);
  dlogbr=log(brace);

  // Shapiro delay, ehm06, Eq. 79, identical to DD86, Eq. 26.
  ds=-2*m2*dlogbr;

  // Aberration delay, ehm06, Eq. 76, identical to DD86, Eq. 27.
  da=a0*(sin(omega+ae) + ecc*sw) + b0*(cos(omega+ae) + ecc*cw); // Zero by design

  // Now compute d2bar, the orbital time correction in DD equation 42.
  d2bar=dre*(1-anhat*drep+(pow(anhat,2))*(pow(drep,2) + 0.5*dre*drepp -
					  0.5*ecc*su*dre*drep/onemecu)) + ds + da;
  torb=-d2bar;

  if (param==-1) return torb;
  
  // Now we need the partial derivatives. Use DD equations 62a - 62k.
  csigma=x*(-sw*su+sqr1me2*cw*cu)/onemecu;
  ce=su*csigma-x*sw-ecc*x*cw*su/sqr1me2;
  cx=sw*cume+sqr1me2*cw*su;
  comega=x*(cw*cume-sqr1me2*sw*su);
  cgamma=su;
  cdth=-ecc*ecc*x*cw*su/sqr1me2;
  cm2=-2*dlogbr;
  csi=2*m2*(sw*cume+sqr1me2*cw*su)/brace; 

  long double lgf = log( ( 1.0 + ecc * cos( ae ) ) / 
                         (1.0 + pow( stig, 2.0 ) - 2.0 * stig * sin( ae + omega ) ) );

  if( param == param_pb )
    return -csigma * an * SECDAY * tt0 / ( pb * SECDAY ); 
  
  else if( param == param_a1 )
    return cx;

  else if( param == param_ecc )
    return ce;

  else if( param == param_om )
    return comega;

  else if( param == param_omdot )
    return ae*comega/(an*360.0/(2.0*M_PI)*365.25*SECDAY);

  else if( param == param_t0 )
    return -csigma*an*SECDAY;

  else if( param == param_pbdot )
    return 0.5*tt0*(-csigma*an*SECDAY*tt0/(pb*SECDAY));

  else if( param == param_sini )
    return csi;

  else if( param == param_gamma )
    return cgamma;

  else if( param == param_m2 )
    return cm2*SUNMASS;

  else if( param == param_a1dot ) // Also known as xdot
    return cx*tt0;

  else if( param == param_h3 ) 
    return( 2.0 / pow( stig, 3.0 ) * lgf );

  else if( param == param_stig )
    return( -2.0 * m2 / stig * 
            ( 1.0 + 3.0 * lgf - ( 1.0 - stig * stig ) / 
              ( 1.0 + stig * stig - 2.0 * stig * sin( ae + omega ) ) ) );

  return 0.0;
}


void updateDDH(pulsar *psr,double val,double err,int pos){
  if (pos==param_pb){
    psr->param[param_pb].val[0] += val/SECDAY;
    psr->param[param_pb].err[0]  = err/SECDAY;
  }else if( pos == param_a1 || pos == param_ecc || pos == param_t0 || 
            pos == param_sini || pos == param_m2 || pos == param_gamma || 
            pos == param_h3 || pos == param_stig ){
    psr->param[pos].val[0] += val;
    psr->param[pos].err[0]  = err;
  }else if (pos==param_om){
    psr->param[pos].val[0] += val*180.0/M_PI;
    psr->param[pos].err[0]  = err*180.0/M_PI;
  }else if (pos==param_pbdot){
    psr->param[pos].val[0] += val;
    psr->param[pos].err[0]  = err;
  }else if (pos==param_a1dot){
    psr->param[pos].val[0] += val;
    psr->param[pos].err[0]  = err;
  }else if (pos==param_omdot){
    psr->param[pos].val[0] += val; /* *(SECDAY*365.25)*180.0/M_PI; */
    psr->param[pos].err[0]  = err; /* *(SECDAY*365.25)*180.0/M_PI; */
  }
}
