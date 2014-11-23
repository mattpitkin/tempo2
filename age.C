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
*    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issuse 2, 
*    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
*    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
*    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
*    timing model.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include "age.h"
//#include "/u/sha355/include/nr_c.h"




void rk4rms(double y[], double  dydx[], int n, double x, double h, double yout[],
	void (*derivs)(double, double [], double []))
{
	int i;
	double xh,hh,h6;

        double dym[n];
	double  dyt[n];
	double yt[n];
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
	(*derivs)(xh,yt,dyt);
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
	(*derivs)(xh,yt,dym);
	for (i=0;i<n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt);
	for (i=0;i<n;i++)
	  yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	//free_vector(yt,1,n);
	//free_vector(dyt,1,n);
	//free_vector(dym,1,n);
	return;
}


void derivs(double  x,double y[], double dydx[])
{

  //double lambda_tot=3;

  //y[0]= frequency
  // first deriv of freq is freq derive

 

  dydx[0]=y[1];
  
  //y[1]= frequency derivative
  dydx[1]= 3*y[1]*y[1]/y[0]-y[1]*LAMBDA_EXCESS;
  
  //fprintf(stderr, "%.3e %.3e %.3e %.3e\n", dydx[1], dydx[2], y[1], y[2]);

  
  return;
}


void calc_age(double f0, double f1,double f2, double *age)
{
  
  double freq[2];
  double fderiv[2];
  double freqout[2];

  //double f0;
  //double f1;
  //double f2;
  double tyr=1000;
  double tscale=86400.*365.25*tyr;
 
  
  //Vela
  //f0=11.2;
  //f1=-1.5e-11;
  //f2=7.4e-22;
  //f2=3*f1*f1/f0;
  
  // J0857-4424
  //f0=3.06;
  //f1=-2e-13;
  //f2=0.;
  //f2=3*f1*f1/f0;
  //f2=6e-23;
  
  fprintf(stderr, "%.3e\n", f2);

  
  //f2=0;
  

  freq[0]=f0;
  freq[1]=f1*tscale;
  fderiv[0]=f1*tscale;
  fderiv[1]=f2*tscale*tscale;
  freqout[0]=0;
  freqout[1]=0;

  LAMBDA_EXCESS = (double) (f2/f1-3.*f1/f0)*tscale;
  //LAMBDA_EXCESS=0; 

  //fprintf(stderr, "excess: %.3le %.5e\n", LAMBDA_EXCESS, (f2/f1 -3.*f1/f0)/(3*f1/f0  ));


  double h;
  // h is in units of tscale
  h=0.05;
  int i;
  i=0;
  
    while(freqout[0] < 1e2)
      {
	//	fprintf(stderr, "%.3e\n", freqout[0]);
      
      

      
      rk4rms(freq,fderiv,2,0.,-h*i,freqout,derivs);
      
      //update derivatives
      //derivs(0,freqout,

      //freq[0]=freqout[0];
      //freq[1]=freqout[1];
      

      //fderiv[0]=freqout[1];

      //fprintf(stdout, "%.3e %.3e %.3e %.3e %.3e\n", i*h*tyr, freqout[0], freqout[1],fderiv[0], fderiv [1]);
      i++;
    }
    
    //fprintf(stderr, "value here %.3e\n", h*i*tyr);
    //exit(0);
    
    *age= (double) h*i*tyr;
    


  return;
}

