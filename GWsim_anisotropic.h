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


/*
 * Set of routines necessary for simulating gravitational wave sources on pulsar
 * timing data
 *
 * Most of these routines are described in Hobbs, Jenet, Verbiest, Yardley, Manchester,
 * Lommen, Coles, Edwards & Shettigara (2009) MNRAS: "TEMPO2, a new pulsar timing 
 * package. III: Gravitational wave simulation
 *
 * These routines are used in the following plugins:
 *
 * GWbkgrd: simulation of a GW background
 *
 * This code has been developed by G. Hobbs and D. Yardley with help from F. Jenet and K. J. Lee.
 */

#include "T2toolkit.h"
#include "GWsim.h"

/* Routine to generate a normalised associated Legendre polynomial suitable for construction of spherical harmonics on the sky. */
double sphharm(int l, int m, double x) {
	//printf("In sphharm l = %i, m = %i, x = %lf\n",l,m,x);
	double fact,oldfact,pll,pmm,pmmp1,pmmp2,omx2;
	int i,ll;
	if (m < 0 || m > l || abs(x) > 1.) {
		printf("Bad arguments in routine sphharm!\n");
		exit(0);
	}
	pmm=1.;
	if (m > 0) {
		omx2=(1.-x)*(1.+x);
		fact=1.;
		for (i=1;i<=m;i++) {
			pmm*=omx2*fact/(fact+1.);
			fact+=2.;
		}
	}
	pmm=sqrt(((double)(2*m+1))*pmm/(4.*M_PI));
	if (m & 1)
		pmm=-pmm;
	if (l == m)
		return(pmm);
	else {
		pmmp1=x*sqrt((double)(2*m+3))*pmm;
		if (l == (m+1))
			return(pmmp1);
		else {
			oldfact=sqrt((double)(2*m+3));
			for (ll=m+2;ll<=l;ll++) {
				fact=sqrt(((double)(4*ll*ll-1))/((double)(ll*ll-m*m)));
				pll=(x*pmmp1-pmm/oldfact)*fact;
				oldfact=fact;
				pmm=pmmp1;
				pmmp1=pll;
			}
			return(pll);
		}	
	}
}

int Ngrid=1000;

/* Produce an anisotropic background of gravitational waves */
void GWanisotropicbackground(gwSrc *gw,int numberGW,long *idum,long double flo,long double fhi,
		  double gwAmp,double alpha,int loglin, double **harmlist, int nharms)
{
  int k;
  while (Ngrid*Ngrid < numberGW)
	Ngrid*=2;
  double *probs,min,norm,dcth,dphi,cth,phi,isocomp,*plm,sum;
  probs=(double *)malloc(2.*Ngrid*Ngrid*sizeof(double));
  plm=(double *)malloc(nharms*sizeof(double));
  norm=0.;
  dcth=2./((double)Ngrid);
  dphi=2.*M_PI/((double)Ngrid);
  cth=-1.+0.5*dcth;
  isocomp=0.;
  for (int hi=0;hi<nharms;hi++) {
	if (!(harmlist[hi][0]) && !(harmlist[hi][1]))
		isocomp=harmlist[hi][2];
  }
  k=0;
  for (int i=0;i<Ngrid;i++) {
	phi=0.5*dphi;
	for (int hi=0;hi<nharms;hi++) {
		plm[hi]=sphharm((int)(harmlist[hi][0]),abs((int)(harmlist[hi][1])),cth);
		//printf("%lf\n",plm[hi]);
	}
	for (int j=0;j<Ngrid;j++) {
		probs[k]=0.;
		for (int hi=0;hi<nharms;hi++) {
			if ((int)(harmlist[hi][1]) >= 0) 
				probs[k]+=harmlist[hi][2]*plm[hi]*cos(harmlist[hi][1]*phi);
			else
				probs[k]+=harmlist[hi][2]*plm[hi]*sin(-harmlist[hi][1]*phi);
		}
		if (i+j) {
			if (probs[k]<min)
				min=probs[k];
		} else
			min=probs[k];
		norm+=probs[k];
		k++;
		phi+=dphi;
	}
	cth+=dcth;
  }
  free(plm);
  if (min < 0.) {
	printf("Warning: distribution computed from harmonic file has negative probabilities at cetain sky positions. Increasing isotropic component amplitude to %6.4lf to renormalise this to a non-negative pdf.\n",isocomp-min);
	norm+=Ngrid*Ngrid*(-min);
  } else
	min=0.;
  sum=0.;
  for (int i=0;i<Ngrid*Ngrid;i++) {
	sum+=(probs[i]-min);
	probs[i]=sum/norm;
	//printf("%lf\n",probs[i]);
  }
  double thisprob;
  int highindx,lowindx,newindx;

  for (k=0;k<numberGW;k++)
    {
	thisprob=TKranDev(idum);
	lowindx=-1;
	highindx=Ngrid*Ngrid-1;
	while (highindx-lowindx > 1) {
		newindx=(int)((highindx+lowindx)/2);
		if (probs[newindx]>thisprob)
			highindx=newindx;
		else
			lowindx=newindx;
	}
      gw[k].theta_g=acos(-1.+((double)((int)(highindx/Ngrid)))*2./((double)(Ngrid))+dcth*TKranDev(idum));
      gw[k].phi_g=(((double)(highindx%Ngrid))+TKranDev(idum))*dphi;
      //printf("OUTPUT %Lf %Lf %i %i\n",cos(gw[k].theta_g),gw[k].phi_g,(int)(highindx/Ngrid),(highindx%Ngrid));
      gw[k].phi_polar_g = 0.0;  
      gw[k].phase_g     = TKranDev(idum)*2*M_PI;   
      if (loglin==1)  /* Use equal sampling in log */
	{
	  gw[k].omega_g  = 2*M_PI*exp(log(flo)+TKranDev(idum)*(log(fhi/flo))); 
	  gw[k].aplus_g  = gwAmp*pow(gw[k].omega_g/(2.0*M_PI), alpha)/sqrt((double)(numberGW)/log(fhi/flo))*TKgaussDev(idum);
	  gw[k].across_g = gwAmp*pow(gw[k].omega_g/(2.0*M_PI), alpha)/sqrt((double)(numberGW)/log(fhi/flo))*TKgaussDev(idum);
	  gw[k].aplus_im_g = 0.0;
	  gw[k].across_im_g = 0.0;
	}      
      else
	{
	  gw[k].omega_g = 2*M_PI*(TKranDev(idum)*(fhi-flo)+flo); 
	  gw[k].aplus_g = gwAmp*pow(gw[k].omega_g/(2.0*M_PI), alpha)*sqrt(1.0/(double)numberGW)*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
	  gw[k].across_g = gwAmp*pow(gw[k].omega_g/(2.0*M_PI), alpha)*sqrt(1.0/(double)numberGW)*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum); 
	  gw[k].aplus_im_g = 0.0;
	  gw[k].across_im_g = 0.0;
	}
      //      printf("omega = %g\n",(double)(1.0/(gw[k].omega_g/2.0/M_PI)/86400.0));

    }
    free(probs);
}

