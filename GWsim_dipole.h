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

/* Find numerical root for phi with tolerance TOL to the equation 

prob = phi/(2 pi) + amp * [sin(phi + phase) - sin(phase)]

with phi in the range [0, 2 pi]. */
double Findphi(double prob, double amp, double phase) {
	//printf("%lf %lf %lf\n",prob,amp,phase);
	double TOL=1.e-10;
	double low,high,lowval,highval,phi,phival,twopi;
	low=0.;
	twopi=2.*M_PI;
	high=twopi;
	lowval=0.;
	highval=1.;
	while (high-low > TOL) {
		phi=0.5*(low+high);
		phival=phi/twopi+amp*sin(phi+phase)-amp*sin(phase)-prob;
		if ((phival*highval) > 0.)
			high=phi;
		else
			low=phi;
		//printf("%6.4e %6.4e %6.4e\n",lowval,highval,phival);
	}
	return(phi);
}

/* Produce a dipole background of gravitational waves */
void GWdipolebackground(gwSrc *gw,int numberGW,long *idum,long double flo,long double fhi,
		  double gwAmp,double alpha,int loglin, double *dipoleamps)
{
  int k;
  printf("DIPOLE AMPS: %lf %lf %lf\n",dipoleamps[0],dipoleamps[1],dipoleamps[2]);

  double cth,thisprob,dipamp,dipphs,twopi;
  double sign=1.;
  if (dipoleamps[0]<0.)
	sign=-1.;
  twopi=2.*M_PI;
  dipamp=sqrt(dipoleamps[1]*dipoleamps[1]+dipoleamps[2]*dipoleamps[2]);
  if (dipoleamps[1])
  	dipphs=atan(-dipoleamps[2]/dipoleamps[1]);
  else
	dipphs=M_PI/2.;
  if (dipoleamps[1] < 0.)
	dipphs+=M_PI;
  for (k=0;k<numberGW;k++)
    {
      //gw[k].theta_g     = acos((TKranDev(idum)-0.5)*2);
      thisprob=TKranDev(idum);
      cth=(-0.5+sqrt(0.25-2.*twopi*dipoleamps[0]*(0.5-thisprob-M_PI*dipoleamps[0])))/(twopi*dipoleamps[0]);
      
      gw[k].theta_g     = acos(cth);
      gw[k].phi_g       = Findphi(TKranDev(idum),dipamp*sqrt(1.-cth*cth)/(0.5+twopi*cth*dipoleamps[0]),dipphs);   
      //printf("OUTPUT %6.4e %Lf\n",cth,gw[k].phi_g);

      gw[k].phi_polar_g = 0.0; 
      gw[k].phase_g     = TKranDev(idum)*twopi;
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

}

