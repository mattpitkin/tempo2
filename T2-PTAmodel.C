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

/* Timing model for tempo2 for PTA-pulsars                */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"

long double computeU(long double phase,long double e);

double T2_PTAmodel(pulsar *psr,int p,int ipos,int param,int arr)
{
  const char *CVS_verNum = "$Revision: 1.1 $";
  long double SUNMASS = 4.925490947e-6;
  long double T0,tbbat,dt;
  long double e,e0,edot;
  long double pb,nval,orbits;
  long double u;
  int norbits;
  long double phase;
  long double Ae,kval;
  long double omega;
  long double Sval,basicRoemer,x;
  long double pbdot;
  long double deltaB=0.0L,dDB;
  long double deltaSB=0.0L,model,r,s,deltaSR=0.0L,deltaAOP=0.0;
  long double kom,kin,pmra,pmdec,dt0,Cval;
	

  // MUST CHECK IMPLEMENTATION OF OMDOT -- DIFFERENT FROM DD MODEL

  if (displayCVSversion == 1) 
    CVSdisplayVersion((char *)"T2_PTAmodel.C",(char *)"T2_PTAmodel()",CVS_verNum);

  if (psr[p].param[param_t0].paramSet[0] == 0) {
    printf("Must use T0 in the binary model\n"); exit(1);
  }
  T0 = psr[p].param[param_t0].val[0];
  do {
    tbbat = psr[p].obsn[ipos].bbat - deltaB/86400.0L;
    
    dt = (tbbat - T0)*86400.0L;
    
    // Eccentricity
    if (psr[p].param[param_ecc].paramSet[0] == 0) {
      printf("Must use ECC in the binary model\n"); exit(1);
    }
    e0 = psr[p].param[param_ecc].val[0];
    if (psr[p].param[param_edot].paramSet[0] == 1) 
      edot = psr[p].param[param_edot].val[0];
    else
      edot = 0.0L;
    
    // Equation 64 in Edwards et al. 2006
    e = e0 + edot*dt;
    
    // Equation 61 in Edwards et al. 2006
    if (psr[p].param[param_pb].paramSet[0] == 0) {
      printf("Must use PB in the binary model\n"); exit(1);
    }
    pb = psr[p].param[param_pb].val[0]*86400.0L;
    
    nval = 2.0L*M_PI/pb; 
    pbdot = 0.0L;
    if (psr[p].param[param_pbdot].paramSet[0]==1)
      {
	pbdot = psr[p].param[param_pbdot].val[0];
	// This is a negative sign in the DD model
	// THIS NEEDS CAREFUL CHECKING
	nval += M_PI*pbdot*dt/pb/pb;
      }

    // HOW DOES THIS RELATE TO THE EDWARDS PAPER?
    orbits = nval/2.0L/M_PI*dt;
    norbits = (int)(orbits);
    if (orbits < 0.0L) norbits --;  // SHOULD CHECK THIS
    phase = 2.0L*M_PI*(orbits-norbits);
    
    //    phase = nval*dt;
    //    printf("Phase = %g %g %g %g\n",(double)phase,(double)pb,(double)dt,(double)norbits);
    // Equation 56: must solve for u
    u = computeU(phase,e);
    //    printf("U = %g\n",(double)u);

    // Equation 58 -- should check how to calculate this atan - USE LONG DOUBLE MATHS?
    //Ae = 2.0L*atanl(powl((1.0L+e)/(1.0L-e),2)*tanl(u/2.0L));
    Ae = 2.0L*atan2((1.0L+e)*(1.0L+e)*tan(u/2.0),(1.0L-e)*(1.0L-e));
    
    //    Ae = atan2(sqrt(1.0-e*e)*sin(u)/(1.0-e*cos(u)),(cos(u)-e)/(1.0-e*cos(u)));
        if(Ae<0.0) Ae=Ae+2.0*M_PI;
    // THIS IS DIFFERENT TO THE EDWARD'S ET AL PAPER -- MUST CHECK
    //Ae = Ae + (int)phase;
    //    Ae += 2*M_PI*orbits - phase;
    
    // Equation 63
    if (psr[p].param[param_omdot].paramSet[0]==1){
      kval = psr[p].param[param_omdot].val[0]*M_PI/180.0/(365.25*86400.0L)/nval;
    }
    else
      kval = 0.0L;
    
    // Equation 62
    if (psr[p].param[param_om].paramSet[0] == 0) {
      printf("Must use OM in the binary model\n"); exit(1);
    }
    omega = (psr[p].param[param_om].val[0]*M_PI/180.0L+kval*Ae);
    //    printf("omega %g %g %g %g %g\n",(double)psr[0].obsn[ipos].bbat,(double)omega,(double)Ae,(double)kval,(double)(kval*Ae));
    // Equation 65
    Sval = sinl(omega)*(cosl(u)-e)+cosl(omega)*sqrt(1.0L-e*e)*sinl(u);
    
    
    // Equation 71
    if (psr[p].param[param_a1].paramSet[0] == 0) {
      printf("Must use A1 in the binary model\n"); exit(1);
    }
    x = psr[p].param[param_a1].val[0];
    if (psr[p].param[param_a1dot].paramSet[0] == 1) 
      x += (psr[p].param[param_a1dot].val[0]*dt);
    
    // Equation 70
    basicRoemer = x*Sval;

    // Calculate Shapiro delay
    // Note that we use KIN instead of sin(i)
    if (psr[p].param[param_kin].paramSet[0]==1 &&
	psr[p].param[param_m2].paramSet[0]==1)
      {
	kin = psr[p].param[param_kin].val[0]*M_PI/180.0L;
	s = sin(kin);
	r = psr[p].param[param_m2].val[0]*SUNMASS;
	deltaSB = -2*r*log(1.0L-e*cosl(u) - s*(sinl(omega)*(cosl(u)-e) + sqrt(1.0-e*e)*cos(omega)*sin(u)));
      }

    // Calculate Kopeikin terms
    if (psr[p].param[param_pmra].paramSet[0]==1 &&
	psr[p].param[param_pmdec].paramSet[0]==1 &&
	psr[p].param[param_kin].paramSet[0]==1 &&
	psr[p].param[param_kom].paramSet[0]==1)
      {
	dt0 = (psr[p].obsn[ipos].bbat - psr[p].param[param_posepoch].val[0])*86400.0L;

	kin = psr[p].param[param_kin].val[0]*M_PI/180.0;
	kom = psr[p].param[param_kom].val[0]*M_PI/180.0;
	pmra = psr[p].param[param_pmra].val[0]*M_PI/(180.0*3600.0e3)/(365.25*86400.0L);
	pmra = psr[p].param[param_pmdec].val[0]*M_PI/(180.0*3600.0e3)/(365.25*86400.0L);

	Cval = cosl(omega)*(cosl(u)-e)-sinl(omega)*sqrtl(1.0L-e*e)*sinl(u);
	
	// Equation 73
	deltaSR = x*dt0*((pmra*sin(kom)+pmdec*cos(kom))*Cval/sin(kin) + (pmra*cos(kom)-pmdec*sin(kom))*Sval/tan(kin));

	// Equation 74
	{
	  long double pi_daop,rdote1,rdote2,r1,r2,r3,sin_delta,cos_delta,sin_alpha,cos_alpha;
	  long double px;
	  
	  r1 = psr->obsn[ipos].earth_ssb[0]/AULTSC;
	  r2 = psr->obsn[ipos].earth_ssb[1]/AULTSC;
	  r3 = psr->obsn[ipos].earth_ssb[2]/AULTSC;

	  /* Obtain vector pointing at the pulsar */
	  sin_delta = psr->obsn[ipos].psrPos[2];
	  cos_delta = cos(asin(sin_delta));
	  sin_alpha = psr->obsn[ipos].psrPos[1]/cos_delta;
	  cos_alpha = psr->obsn[ipos].psrPos[0]/cos_delta;

	  // These have been copied from T2model.C -- MUST CHECK
	  rdote1 = -psr->obsn[ipos].earth_ssb[0]/AULTSC*sin_alpha+
	    psr->obsn[ipos].earth_ssb[1]/AULTSC*cos_alpha;
	  rdote2 = -psr->obsn[ipos].earth_ssb[0]/AULTSC*sin_delta*cos_alpha-
	    psr->obsn[ipos].earth_ssb[1]/AULTSC*sin_delta*sin_alpha+
	    psr->obsn[ipos].earth_ssb[2]/AULTSC*cos_delta;


	  px = psr[p].param[param_px].val[0];
	  pi_daop = px*M_PI/180.0/3600.0*1e-3;
	  deltaAOP = -x*pi_daop*((rdote1*sin(kom)+rdote2*cos(kom))*Cval/sin(kin) + (rdote1*cos(kom)-rdote2*sin(kom))*Sval/tan(kin));
	}
      }

    model = basicRoemer + deltaSB + deltaSR + deltaAOP;
    dDB = model - deltaB;
    deltaB = model;
    //    printf("Here with %g %g\n",(double)deltaB,(double)dDB);
  } while (fabs(dDB) > 1e-14);
  printf("basicRoemer = %g\n",(double)basicRoemer);
  // Why negative sign?
  if (param==-1) return -model;

  // Calculate differentials
  long double du_dpb,dD_dpb,dD_da1,dD_domega,du_dT0,dD_dT0;
  long double du_de,dD_de,du_dpdot,dD_dpdot,dD_domegaDot;
  long double domega_domegaDot,dD_dkom,dD_di,dD_dm2;

  // Not including time derivatives of anything yet
  //  du_dpb = (-2.0L*M_PI*dt/pb/pb - 2.0*M_PI*dt*dt*pbdot/powl(pb,3))/(1.0L-e*cosl(u));
  du_dpb = (-2.0L*M_PI*dt/pb/pb)/(1.0L-e*cosl(u));
  du_dT0 = -2.0L*M_PI/pb/(1.0L-e*cosl(u));
  du_de  = sinl(u)/(1.0L-e*cosl(u));
  // -ve sign in DDmodel.C
  du_dpdot = M_PI*dt*dt/pb/pb/(1.0L-e*cosl(u));

  // This seems odd
  domega_domegaDot = Ae/(nval*180.0/M_PI*365.25*SECDAY);

  dD_da1 = sinl(omega)*(cosl(u)-e)+cosl(omega)*sqrt(1.0L-e*e)*sinl(u);
  dD_de = -x*sinl(omega)*sinl(u)*du_de - x*sinl(omega) 
    - x*cosl(omega)/sqrt(1.0L-e*e)*e*sinl(u) + x*cosl(omega)*sqrt(1-e*e)*cosl(u)*du_de;
  dD_dT0 = (-x*sinl(omega)*sinl(u) + x*cosl(omega)*sqrt(1.0L-e*e)*cosl(u))*du_dT0;
  dD_dpb = (-x*sinl(omega)*sinl(u) + x*cosl(omega)*sqrt(1.0L-e*e)*cosl(u))*du_dpb;
  dD_dpdot = (-x*sinl(omega)*sinl(u) + x*cosl(omega)*sqrt(1.0L-e*e)*cosl(u))*du_dpdot;
  dD_domega = x*(cosl(u)-e)*cosl(omega) - x*sqrt(1.0L-e*e)*sinl(u)*sinl(omega);
  dD_domegaDot = dD_domega*domega_domegaDot;


  // Currently only including Shapiro term
  dD_dm2 = deltaSB/r;
  dD_di = -2*r/(1.0L-e*cosl(u) - s*(sinl(omega)*(cosl(u)-e) + sqrt(1.0-e*e)*cos(omega)*sin(u)))*(-cos(kin)*(sinl(omega)*(cosl(u)-e) + sqrt(1.0-e*e)*cos(omega)*sin(u)));

  // 
  if (param==param_kom)
    dD_dkom = x*dt0*(Cval/sin(kin)*(pmra*cos(kom)-pmdec*sin(kom)) +
		     (Sval/tan(kin)*(-pmra*sin(kom) - pmdec*cos(kom))));

  printf("output: %g %g %g %g %g %g\n",(double)psr[p].obsn[ipos].bbat,(double)u,(double)Ae,(double)nval,(double)deltaSR,(double)deltaAOP);
  if (param==param_pb)
    return dD_dpb;
  else if (param==param_pbdot)
    return dD_dpdot;
  else if (param==param_a1)
    return dD_da1;
  else if (param==param_a1dot)
    return dD_da1*dt;
  else if (param==param_om)
    return dD_domega;
  else if (param==param_omdot)
    return dD_domegaDot;
  else if (param==param_kom)
    return dD_dkom;
  else if (param==param_kin)
    return dD_di;
  else if (param==param_m2)
    return dD_dm2;
  else if (param==param_t0)
    {
      printf("dT0 = %g\n",(double)dD_dT0);
      return dD_dT0*86400.0L;
    }
  else if (param==param_ecc)
    return dD_de;
}

long double computeU(long double phase,long double ecc)
{
  long double du;
  long double u;

  // DD uses this. T2 uses the one below.
  u=phase+ecc*sin(phase)*(1.0+ecc*cos(phase));
  //  u = phase+ecc*sin(phase)/sqrt(1.0-2*ecc*cos(phase)+ecc*ecc);
  do {
    du=(phase-(u-ecc*sin(u)))/(1.0-ecc*cos(u));
    (u)+=du;
  } while (fabs(du)>1.0e-14);
  return u;
}

void updateT2_PTA(pulsar *psr,double val,double err,int pos,int arr){
  long double SUNMASS = 4.925490947e-6;
  if (pos==param_pb) // || pos==param_t0)
    {
      printf("In here with param_pb: %g %g %g\n",(double)psr->param[pos].val[0],(double)val,(double)err);
      psr->param[pos].val[0] += val/SECDAY;
      psr->param[pos].err[0]  = err/SECDAY;
    }
  else if (pos==param_om) // || pos==param_t0)
    {
      psr->param[pos].val[0] += val*180.0L/M_PI;
      psr->param[pos].err[0]  = err*180.0L/M_PI;
    }
  else if (pos==param_kom) // || pos==param_t0)
    {
      psr->param[pos].val[0] += val*180.0L/M_PI;
      psr->param[pos].err[0]  = err*180.0L/M_PI;
    }
  else if (pos==param_kin) 
    {
      psr->param[pos].val[0] += val*180.0L/M_PI;
      psr->param[pos].err[0]  = err*180.0L/M_PI;
    }
  else if (pos==param_m2) 
    {
      psr->param[pos].val[0] += val/SUNMASS;
      psr->param[pos].err[0]  = err/SUNMASS;
    }
  else if (pos==param_omdot) 
    {
      psr->param[pos].val[0] += val;//*180.0L/M_PI*86400.0*365.25;
      psr->param[pos].err[0]  = err;//*180.0L/M_PI*86400.0*365.25;
    }
  else 
    {
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
    }
}
