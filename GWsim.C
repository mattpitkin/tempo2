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

#include "tempo2.h"
#include "GWsim.h"
#include "T2toolkit.h"
#include <math.h>
#include <stdlib.h>


/* Set up GW:                                                              
 * Sets up the vector pointing at the GW source                            
 * Sets up the strain matrix
 */


int gwsim_Ngrid=1000;


#ifdef __cplusplus
extern "C"
#endif
void setupGW(gwSrc *gw)
{
  long double deg2rad = M_PI/180.0;
  long double convert[3][3],trans[3][3];
  long double out[3][3];
  long double out_im[3][3];
  int i,j;
  long double st,ct,cp,sp,cpp,spp;

  st = sin(gw->theta_g);
  ct = cos(gw->theta_g);
  sp = sin(gw->phi_g);
  cp = cos(gw->phi_g);
  cpp = cos(gw->phi_polar_g);
  spp = sin(gw->phi_polar_g);

  gw->kg[0] = st*cp;
  gw->kg[1] = st*sp;
  gw->kg[2] = ct;

  /* Set up the GW strain: assume general relativity */
  gw->h[0][0] = gw->aplus_g;
  gw->h[0][1] = gw->across_g;
  gw->h[1][0] = gw->across_g;
  gw->h[1][1] = -gw->aplus_g;
  gw->h[2][0] = 0.0;
  gw->h[2][1] = 0.0;
  gw->h[2][2] = 0.0;
  gw->h[0][2] = 0.0;
  gw->h[1][2] = 0.0;
  
  gw->h_im[0][0] = gw->aplus_im_g; // ALL SHOULD HAVE sqrt(-1) *  out the front
  gw->h_im[0][1] = gw->across_im_g; // only used for calcuateResidualGW and setupGW
  gw->h_im[1][0] = gw->across_im_g;
  gw->h_im[1][1] = -gw->aplus_im_g;
  gw->h_im[2][0] = 0.0;
  gw->h_im[2][1] = 0.0;
  gw->h_im[2][2] = 0.0;
  gw->h_im[0][2] = 0.0;
  gw->h_im[1][2] = 0.0;
  
  /* Now carry out a coordinate conversion (to convert to same coordinates 
     as the pulsar and GW source coordinates?)                                */

  /* These are from Euler's angles with pitch-roll-yaw convention */
  /* Also different definition of theta in MathWorld (mw) => sin(theta_mw) = -cos(theta_g),
     cos(theta_mw) = sin(theta_g) */
  /* Should check as there seems to be a minus sign difference in the top line
     with http://mathworld.wolfram.com/EulerAngles.html. Also the top line
     of this matrix is the bottom line of the wolfram matrix.  Probably doesn't
     matter, but should check */

  convert[0][0] = cp*ct*cpp-sp*spp; 
  convert[0][1] = sp*ct*cpp+cp*spp;
  convert[0][2] = -st*cpp;
       
  convert[1][0]= -cp*ct*spp-sp*cpp;
  convert[1][1]= -sp*ct*spp+cp*cpp;
  convert[1][2]= st*spp;
  
  convert[2][0]= cp*st;
  convert[2][1]= sp*st;
  convert[2][2]= ct;

  /* Calculate h.convert */
  matrixMult(gw->h,convert,out);
  matrixMult(gw->h_im,convert,out_im);

  /* Calculate transpose of 'convert' */
  for (i=0;i<3;i++)
    {
      for (j=0;j<3;j++)
	{
	  //	  printf("%g (%g, %g) ",(double)convert[j][i],(double)gw->h[j][i],(double)out[j][i]);
	  trans[i][j] = convert[j][i];
	}
      //      printf("\n");
    }
  /* Pre-multiply "out" by transpose */
  /* why multiply here -- we have already rotated it???? */
  matrixMult(trans,out,gw->h);
  matrixMult(trans,out_im,gw->h_im);
  //  matrixMult(trans,convert,out);
  /*  for (i=0;i<3;i++)
    {
      for (j=0;j<3;j++)
	{
	  printf("%g ",(double)gw->h[j][i]);
	}
      printf("\n");
      }*/

}

/* Multiply 3x3 matrices and put output in "out" */
#ifdef __cplusplus
extern "C"
#endif
long double matrixMult(long double m1[3][3],long double m2[3][3],long double out[3][3])
{
  int i,j,k;
  for (i=0;i<3;i++)
    {
      for (j=0;j<3;j++)
	{
	  out[i][j] = 0.0;
	  for (k=0;k<3;k++)
	    out[i][j]+=m1[i][k]*m2[k][j];
	}
    }
}
/* Calculates the dot product between two vectors */
#ifdef __cplusplus
extern "C"
#endif
long double dotProduct(long double *m1,long double *m2)
{
  long double dprod;

  dprod = m1[0]*m2[0] + m1[1]*m2[1] + m1[2]*m2[2];
  return dprod;
}

/* Produce a background of gravitational waves */
#ifdef __cplusplus
extern "C"
#endif
void GWbackground(gwSrc *gw,int numberGW,long *idum,long double flo,long double fhi,
		  double gwAmp,double alpha,int loglin)
{
  int k;
  for (k=0;k<numberGW;k++)
    {
      gw[k].theta_g     = acos((TKranDev(idum)-0.5)*2);  
      gw[k].phi_g       = TKranDev(idum)*2*M_PI;   
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

}


/* Calculate the pulsar timing residual induced by the gravitational wave */
#ifdef __cplusplus
extern "C"
#endif
long double calculateResidualGW(long double *kp,gwSrc *gw,long double time,long double dist)
{
  long double cosMu;
  long double psrVal1,psrVal2,psrVal1_im;
  long double earthVal1;
  long double tempVal[3];
  long double tempVal_im[3];
  long double geo;
  long double residual;
  long double deg2rad = M_PI/180.0L;
  //  long double VC = 299792456.200L;
  long double VC = 299792458.0L; /* Speed of light (m/s)                       */
  int i,k;

  for (i=0;i<3;i++)
    {
      tempVal[i] = 0.0;
      tempVal_im[i] = 0.0;
      for (k=0;k<3;k++)
	{
	  tempVal[i]    += gw->h[i][k]*kp[k];
	  tempVal_im[i] += gw->h_im[i][k]*kp[k];
	}
      //      printf("tempVal[%d] = %g ",i,(double)tempVal[i]);
    }
  //  printf("\n");
  //  printf("kp = %g %g %g\n",(double)kp[0],(double)kp[1],(double)kp[2]);

		
  cosMu = dotProduct(kp,gw->kg);   /* Angle between pulsar and GW source */
  //  printf(" cosmu = %g\n",(double)cosMu);
  if ((1+cosMu)!=0) 
  {
    psrVal1    = 0.5L/(1.0L+cosMu)*dotProduct(kp,tempVal); 
    psrVal1_im = 0.5L/(1.0L+cosMu)*dotProduct(kp,tempVal_im);
    //    printf("psrval1 = %g %g %g\n",(double)psrVal1,(double)(1+cosMu),(double)(1-cosMu));
  }
  else psrVal1 = 0.0L;
  
  geo      = psrVal1;
  //  printf("psrVal = %g %g\n",(double)psrVal1,(double)psrVal1_im);
  earthVal1 = psrVal1*sinl(gw->phase_g+time*gw->omega_g)+ 
    psrVal1_im*cosl(gw->phase_g+time*gw->omega_g); // sin -> cos

  //    earthVal1 = 0;
  //  printf("dist = %g\n",(double)dist);
  if (dist==0) /* No pulsar term */
    psrVal2=0.0L;
  else
    psrVal2 = psrVal1*sinl(gw->phase_g-(1+cosMu)*dist/VC*gw->omega_g+time*gw->omega_g) + psrVal1_im*cosl(gw->phase_g-(1+cosMu)*dist/VC*gw->omega_g+time*gw->omega_g); 

    residual = (earthVal1-psrVal2)/gw->omega_g; 
  

  //    printf("GWsim %g %g %g %g %g %g %g %g %g %g %g %g %g\n",(double)residual,(double)earthVal1,(double)psrVal2,(double)psrVal1,(double)cosMu,(double)psrVal1_im,(double)kp[0],(double)kp[1],(double)kp[2],(double)cosMu,(double)gw->h_im[0][0],(double)gw->h_im[0][1],(double)gw->h_im[1][0]);
  return residual;
}

/* Set up pulsar: note: in KJ this is based in Gwave.cpp: Load_Pulsar_Data */
/* Sets up the vector pointing at the pulsar                               */
#ifdef __cplusplus
extern "C"
#endif
void setupPulsar_GWsim(long double ra_p,long double dec_p,long double *kp)
{
  kp[0] = cos(dec_p)*cos(ra_p);
  kp[1] = cos(dec_p)*sin(ra_p);
  kp[2] = sin(dec_p);
}



#ifdef __cplusplus
extern "C"
#endif
int GWbackground_read(gwSrc *gw, FILE *file, int ireal){
	char key[13];
	int nreal;
	int ngw,id,igw,i;
	const unsigned int gwsize = 9*sizeof(long double);

	for (i=0;i<=ireal;i++){
		fread(&id,sizeof(int),1,file);
		printf("%d\n",id);
		if(id==ireal)break;
		fread(&ngw,sizeof(int),1,file);
		fseek(file,ngw*gwsize,SEEK_CUR);
		if(feof(file)){
			fprintf(stderr,"Could not read enough file to find realiation required\n");
			return -1;
		}
	}

	fread(&ngw,sizeof(int),1,file);

	printf("Reading %d GW sources from real %d (%d)\n",ngw,id,ireal);
	
	for (igw=0;igw<ngw;igw++){
		fread(&(gw[igw]),gwsize,1,file);
	}

	return ngw;
}



#ifdef __cplusplus
extern "C"
#endif
void GWbackground_write(gwSrc *gw, FILE *file,int ngw, int ireal){
	int igw;
	const unsigned int gwsize = 9*sizeof(long double);
	fwrite(&ireal,sizeof(int),1,file);
	fwrite(&ngw,sizeof(int),1,file);
	for (igw=0;igw<ngw;igw++){
		fwrite(&(gw[igw]),gwsize,1,file);
	}
	fflush(file);
}


// Vikram Ravi's code to include an eccentric, binary black hole system into the 
// tempo2 timing model
#ifdef __cplusplus
extern "C"
#endif
long double eccRes(pulsar *psr, int i, int *coalesceFlag, double *prev_p, double *prev_e, double *prev_a, double *prev_epoch, double *prev_theta) 
{


  long double pert;

  if (*coalesceFlag==1) {
    pert=0.;
    return pert;
  }
  else 
    {

      double ra_p, dec_p, ra_g, dec_g;
      double m1, m2, e_0, inc, theta_n, theta_0, phi, z, D, t_0, p_0;
      int ih;
      
      

      const double G = 6.67e-11;     // in m^3/kg/s^2
      const double c = 2.998e8;      // in m/s
      const double msun = 1.9891e30; // in kg
      const double Mpc_to_m = 3.08568e22; // in m
      const double yr_to_sec = 365.25*86400.; // in s
      
      ra_p = (double)psr->param[param_raj].val[0];
      dec_p   = (double)psr->param[param_decj].val[0];
      ra_g  = psr->gwecc_ra;
      dec_g     = psr->gwecc_dec;
      m1 = psr->gwecc_m1*msun; // READ IN Msun UNITS!
      m2 = psr->gwecc_m2*msun; // READ IN Msun UNITS!
      e_0 = psr->gwecc_e;
      inc = psr->gwecc_inc;
      theta_n = psr->gwecc_theta_nodes;
      theta_0 = psr->gwecc_theta_0;
      phi = psr->gwecc_nodes_orientation;
      z = psr->gwecc_redshift;
      D = psr->gwecc_distance*Mpc_to_m; // READ IN Mpc UNITS!
      t_0 = psr->gwecc_epoch; // READ IN MJD!
      p_0 = (psr->gwecc_orbital_period)*yr_to_sec; // READ IN YEARS AT SOURCE!   

      // schwarzschild radius
      double rs_m1 = Rs(m1);
      
      // current time
      // this is assuming gwecc_epoch defined in earth frame
      double time    = ((double)(psr->obsn[i].sat-psr->obsn[i].residual/86400.L) - *prev_epoch)/(1.+z);
      
      // find integration interval (roughly less than 0.01 days)
      double tint = 86400.*(time)/(fabs(100.*floor(time))); // this is in days here
      int nt = (int)fabs((100.*floor(time)));

      if (nt>1000000) {
	tint=86400.*(time)/1000000.;
	nt=1000000;
	printf("running %d ints with tint %g\n",nt,tint/86400.);
      }
      

      // find theta and ecc using 4th order runge-kutta
      double theta = *prev_theta;
      double ec = *prev_e;
      double ax = *prev_a;
      double M_c = (m1+m2)*pow(m1*m2/pow(m1+m2,2.),3./5.);
      double k1, k2, k3, k4;
      double l1, l2, l3, l4;
      double o1, o2, o3, o4;
		 		 
      ih=0;
      while (ih<nt && *coalesceFlag==0)
	{
	  
	  k1=tint*dadt(ec,ax,m1,m2);
	  l1=tint*dedt(ec,ax,m1,m2);
	  o1=tint*dtdt(ec,theta,*prev_p);
	  
	  if (ax+k1/2.<1.5*rs_m1) {
	    *coalesceFlag=1;
	    printf("***********SOURCE HAS COALESCED***********\n");
	  }

	  k2=tint*dadt(ec+l1/2.,ax+k1/2.,m1,m2);
	  l2=tint*dedt(ec+l1/2.,ax+k1/2.,m1,m2);
	  o2=tint*dtdt(ec+l1/2.,theta+o1/2.,*prev_p);
	  
	  if (ax+k2/2.<1.5*rs_m1) {
	    *coalesceFlag=1;
	    printf("***********SOURCE HAS COALESCED***********\n");
	  }

	  k3=tint*dadt(ec+l2/2.,ax+k2/2.,m1,m2);
	  l3=tint*dedt(ec+l2/2.,ax+k2/2.,m1,m2);
	  o3=tint*dtdt(ec+l2/2.,theta+o2/2.,*prev_p);
	  
	  if (ax+k3<1.5*rs_m1) {
	    *coalesceFlag=1;
	    printf("***********SOURCE HAS COALESCED***********\n");
	  }

	  k4=tint*dadt(ec+l3,ax+k3,m1,m2);
	  l4=tint*dedt(ec+l3,ax+k3,m1,m2);
	  o4=tint*dtdt(ec+l3,theta+o3,*prev_p);
	  
	  
	  ax+=(k1+2.*k2+2.*k3+k4)/6.;
	  ec+=(l1+2.*l2+2.*l3+l4)/6.;
	  theta+=(o1+2.*o2+2.*o3+o4)/6.;
	  

	  if (ax<1.5*rs_m1) {
	    *coalesceFlag=1;
	    printf("***********SOURCE HAS COALESCED***********\n");
	  }
		     
	  *prev_p=(1./(sqrt(G*(m1+m2)*pow(ax,-3.))/(2.*M_PI)));
	  
		     
	  ih++;
	  
	}

      printf("Have time %g theta %g\n",((double)(psr->obsn[i].sat-psr->obsn[i].residual/86400.L-t_0)/365.25)/(*prev_p/(86400.*365.25)),theta/(2.*M_PI));

      if (*coalesceFlag==1) {
	pert=0.;
      }
      else {

	*prev_e = ec;
	*prev_epoch = (double)(psr->obsn[i].sat-psr->obsn[i].residual/86400.L);
	*prev_theta = theta;
	*prev_a = ax;
	
	
	// Now put in the residuals
	

	double Ag = 4*pow(G,2.)*m1*m2/(pow(c,4.)*D*(1-ec*ec)*ax);
	double frontTerm = -Ag*pow(1-ec*ec,3./2.)*(*prev_p)*sin(theta)/(32.*M_PI*(1+ec*cos(theta)));
	double rplus = frontTerm*(ec*cos(2.*(inc-phi))+ec*cos(2.*(inc+phi))-16.*cos(inc)*sin(2.*phi)*sin(theta-2.*theta_n)+8.*ec*cos(inc)*sin(2.*phi)*sin(2.*theta_n)+2.*cos(2.*phi)*(-ec+(3.+cos(2.*inc))*(ec+2.*cos(theta))*cos(2.*theta_n)+2.*(3.+cos(2.*inc))*sin(theta)*sin(2.*theta_n)));
	double rcross = frontTerm*(ec*sin(2.*(inc-phi))-ec*sin(2.*(inc+phi))+8.*cos(inc)*cos(2.*phi)*(-2.*sin(theta-2.*theta_n)+ec*sin(2.*theta_n))-2.*sin(2.*phi)*(-ec+(3.+cos(2.*inc))*(ec+2.*cos(theta))*cos(2.*theta_n)+2.*(3.+cos(2.*inc))*sin(theta)*sin(2.*theta_n)));

	double dalpha = ra_g-ra_p;
	double gplus = (-(1.+3.*cos(2.*dec_p))*pow(cos(dec_g),2.)+cos(2.*dalpha)*(-3.+cos(2.*dec_g))*pow(sin(dec_p),2.)+2.*cos(dalpha)*sin(2.*dec_g)*sin(2.*dec_p))/(8.*(-1.+cos(dec_p)*sin(dec_g)+cos(dalpha)*cos(dec_g)*sin(dec_p)));
	double gcross = ((-cos(dec_g)*cos(dec_p)+cos(dalpha)*sin(dec_g)*sin(dec_p))*sin(dalpha)*sin(dec_p))/(-1.+cos(dec_p)*sin(dec_g)+cos(dalpha)*cos(dec_g)*sin(dec_p));

	pert = 1e9*(rplus*gplus+rcross*gcross);
	
      }

      printf("Have pert %Lg\n",pert);
      return pert;
      

    }

}

#ifdef __cplusplus
extern "C"
#endif
long double eccResWithEnergy(pulsar *psr, int i, int *coalesceFlag, double *prev_p, double *prev_e, double *prev_a, double *prev_epoch, double *prev_theta, float *eOut) 
{


  long double pert;

  if (*coalesceFlag==1) {
    pert=0.;
    *eOut=0.;
    return pert;
  }
  else 
    {

      double ra_p, dec_p, ra_g, dec_g;
      double m1, m2, e_0, inc, theta_n, theta_0, phi, z, D, t_0, p_0;
      int ih;
      
      

      const double G = 6.67e-11;     // in m^3/kg/s^2
      const double c = 2.998e8;      // in m/s
      const double msun = 1.9891e30; // in kg
      const double Mpc_to_m = 3.08568e22; // in m
      const double yr_to_sec = 365.25*86400.; // in s
      
      ra_p = (double)psr->param[param_raj].val[0];
      dec_p   = (double)psr->param[param_decj].val[0];
      ra_g  = psr->gwecc_ra;
      dec_g     = psr->gwecc_dec;
      m1 = psr->gwecc_m1*msun; // READ IN Msun UNITS!
      m2 = psr->gwecc_m2*msun; // READ IN Msun UNITS!
      e_0 = psr->gwecc_e;
      inc = psr->gwecc_inc;
      theta_n = psr->gwecc_theta_nodes;
      theta_0 = psr->gwecc_theta_0;
      phi = psr->gwecc_nodes_orientation;
      z = psr->gwecc_redshift;
      D = psr->gwecc_distance*Mpc_to_m; // READ IN Mpc UNITS!
      t_0 = psr->gwecc_epoch; // READ IN MJD!
      p_0 = (psr->gwecc_orbital_period)*yr_to_sec; // READ IN YEARS AT SOURCE!   

      // schwarzschild radius
      double rs_m1 = Rs(m1);
      
      // current time
      // this is assuming gwecc_epoch defined in earth frame
      double time    = ((double)(psr->obsn[i].sat-psr->obsn[i].residual/86400.L) - *prev_epoch)/(1.+z);
      
      // find integration interval (roughly less than 0.01 days)
      double tint = 86400.*(time)/(fabs(100.*floor(time))); // this is in days here
      int nt = (int)fabs((100.*floor(time)));

      if (nt>1000000) {
	tint=86400.*(time)/1000000.;
	nt=1000000;
	printf("running %d ints with tint %g\n",nt,tint/86400.);
      }
      

      // find theta and ecc using 4th order runge-kutta
      double theta = *prev_theta;
      double ec = *prev_e;
      double ax = *prev_a;
      double M_c = (m1+m2)*pow(m1*m2/pow(m1+m2,2.),3./5.);
      double k1, k2, k3, k4;
      double l1, l2, l3, l4;
      double o1, o2, o3, o4;
		 		 
      ih=0;
      while (ih<nt && *coalesceFlag==0)
	{
	  
	  k1=tint*dadt(ec,ax,m1,m2);
	  l1=tint*dedt(ec,ax,m1,m2);
	  o1=tint*dtdt(ec,theta,*prev_p);
	  
	  if (ax+k1/2.<1.5*rs_m1) {
	    *coalesceFlag=1;
	    printf("***********SOURCE HAS COALESCED***********\n");
	  }

	  k2=tint*dadt(ec+l1/2.,ax+k1/2.,m1,m2);
	  l2=tint*dedt(ec+l1/2.,ax+k1/2.,m1,m2);
	  o2=tint*dtdt(ec+l1/2.,theta+o1/2.,*prev_p);
	  
	  if (ax+k2/2.<1.5*rs_m1) {
	    *coalesceFlag=1;
	    printf("***********SOURCE HAS COALESCED***********\n");
	  }

	  k3=tint*dadt(ec+l2/2.,ax+k2/2.,m1,m2);
	  l3=tint*dedt(ec+l2/2.,ax+k2/2.,m1,m2);
	  o3=tint*dtdt(ec+l2/2.,theta+o2/2.,*prev_p);
	  
	  if (ax+k3<1.5*rs_m1) {
	    *coalesceFlag=1;
	    printf("***********SOURCE HAS COALESCED***********\n");
	  }

	  k4=tint*dadt(ec+l3,ax+k3,m1,m2);
	  l4=tint*dedt(ec+l3,ax+k3,m1,m2);
	  o4=tint*dtdt(ec+l3,theta+o3,*prev_p);
	  
	  
	  ax+=(k1+2.*k2+2.*k3+k4)/6.;
	  ec+=(l1+2.*l2+2.*l3+l4)/6.;
	  theta+=(o1+2.*o2+2.*o3+o4)/6.;
	  

	  if (ax<1.5*rs_m1) {
	    *coalesceFlag=1;
	    printf("***********SOURCE HAS COALESCED***********\n");
	  }
		     
	  *prev_p=(1./(sqrt(G*(m1+m2)*pow(ax,-3.))/(2.*M_PI)));
	  
		     
	  ih++;
	  
	}

      printf("Have time %g theta %g\n",((double)(psr->obsn[i].sat-psr->obsn[i].residual/86400.L-t_0)/365.25)/(*prev_p/(86400.*365.25)),theta/(2.*M_PI));

      if (*coalesceFlag==1) {
	pert=0.;
	*eOut=0.;
      }
      else {

	*prev_e = ec;
	*prev_epoch = (double)(psr->obsn[i].sat-psr->obsn[i].residual/86400.L);
	*prev_theta = theta;
	*prev_a = ax;
	
	
	// Now put in the residuals
	

	double Ag = 4*pow(G,2.)*m1*m2/(pow(c,4.)*D*(1-ec*ec)*ax);
	double frontTerm = -Ag*pow(1-ec*ec,3./2.)*(*prev_p)*sin(theta)/(32.*M_PI*(1+ec*cos(theta)));
	double rplus = frontTerm*(ec*cos(2.*(inc-phi))+ec*cos(2.*(inc+phi))-16.*cos(inc)*sin(2.*phi)*sin(theta-2.*theta_n)+8.*ec*cos(inc)*sin(2.*phi)*sin(2.*theta_n)+2.*cos(2.*phi)*(-ec+(3.+cos(2.*inc))*(ec+2.*cos(theta))*cos(2.*theta_n)+2.*(3.+cos(2.*inc))*sin(theta)*sin(2.*theta_n)));
	double rcross = frontTerm*(ec*sin(2.*(inc-phi))-ec*sin(2.*(inc+phi))+8.*cos(inc)*cos(2.*phi)*(-2.*sin(theta-2.*theta_n)+ec*sin(2.*theta_n))-2.*sin(2.*phi)*(-ec+(3.+cos(2.*inc))*(ec+2.*cos(theta))*cos(2.*theta_n)+2.*(3.+cos(2.*inc))*sin(theta)*sin(2.*theta_n)));

	double dalpha = ra_g-ra_p;
	double gplus = (-(1.+3.*cos(2.*dec_p))*pow(cos(dec_g),2.)+cos(2.*dalpha)*(-3.+cos(2.*dec_g))*pow(sin(dec_p),2.)+2.*cos(dalpha)*sin(2.*dec_g)*sin(2.*dec_p))/(8.*(-1.+cos(dec_p)*sin(dec_g)+cos(dalpha)*cos(dec_g)*sin(dec_p)));
	double gcross = ((-cos(dec_g)*cos(dec_p)+cos(dalpha)*sin(dec_g)*sin(dec_p))*sin(dalpha)*sin(dec_p))/(-1.+cos(dec_p)*sin(dec_g)+cos(dalpha)*cos(dec_g)*sin(dec_p));

	pert = 1e9*(rplus*gplus+rcross*gcross);
	
	
	// do energy estimate

	(*eOut) = (float)(pow(c,3.)/(16.*M_PI*G))*(((Ag*pow(M_PI,3.)*pow(1.+ec*cos(theta),2.))/(pow(*prev_p,2.)*pow(ec*ec-1.,3.)))*(-1536.-6184.*ec*ec-1847.*pow(ec,4.)-96.*pow(ec,6.)+256.*cos(4.*theta)+ec*(-4.*(1488.+1689.*ec*ec+116.*pow(ec,4.))*cos(theta)+2.*ec*(-1270.-243.*ec*ec+48.*pow(ec,4.))*cos(2.*theta)+4.*(32.+267.*ec*ec+76.*pow(ec,4.))*cos(3.*theta)+2.*ec*(200.+479.*ec*ec)*cos(4.*theta)+4.*(176.+87.*ec*ec+40.*pow(ec,4.))*cos(5.*theta)+14.*ec*(46.+5.*ec*ec)*cos(6.*theta)+220.*ec*ec*cos(7.*theta)+25.*pow(ec,3.)*cos(8.*theta))));


      }

      printf("Have pert %Lg energy %f\n",pert,*eOut);
      return pert;
      

    }

}


#ifdef __cplusplus
extern "C"
#endif
double Rs(double m1) {

  double G = 6.67e-11;
  double c = 2.998e8;
  double retval = 2*G*m1/(c*c);

  return retval;

}

#ifdef __cplusplus
extern "C"
#endif
double dadt(double ec, double a, double m1, double m2) {

  double G = 6.67e-11;
  double c = 2.998e8;
  double retval = (-64./5.)*pow(G/a,3.)*m1*m2*(m1+m2)*pow(c,-5.)*Fe(ec);
  return retval;

}

#ifdef __cplusplus
extern "C"
#endif
double dedt(double ec, double a, double m1, double m2) {

  double G = 6.67e-11;
  double c = 2.998e8;
  double retval = (-304./15.)*pow(G,3.)*m1*m2*(m1+m2)*pow(c,-5.)*pow(a,-4.)*pow(1.-ec*ec,-5./2.)*ec*(1.+121.*ec*ec/304.);
  return retval;

}

#ifdef __cplusplus
extern "C"
#endif
double dtdt(double ec, double t, double p) {

  double G = 6.67e-11;
  //double retval = (2.*M_PI/p)*pow(1.+ec*cos(t),2.)/(pow(1.-pow(ec,2.),3./2.));
  double retval = (2.*M_PI/p);
  return retval;

}


#ifdef __cplusplus
extern "C"
#endif
double Fe(double ec) {

  double retval = pow(1.-ec*ec,-7./2.)*(1.+73.*pow(ec,2.)/24.+37.*pow(ec,4.)/96.);
  return retval;

}


// returns radians
#ifdef __cplusplus
extern "C"
#endif
double psrangle(double centre_long,double centre_lat,double psr_long,double psr_lat)
{
  double dlon,dlat,a,c;
  double deg2rad = M_PI/180.0;
  
  /* Apply the Haversine formula */
  dlon = (psr_long - centre_long);
  dlat = (psr_lat  - centre_lat);
  a = pow(sin(dlat/2.0),2) + cos(centre_lat) * 
    cos(psr_lat)*pow(sin(dlon/2.0),2);
  if (a==1)
    c = M_PI;
  else
    c = 2.0 * atan2(sqrt(a),sqrt(1.0-a));  
  return c;
}




/* Routine to generate a normalised associated Legendre polynomial suitable for construction of spherical harmonics on the sky. */
#ifdef __cplusplus
extern "C"
#endif
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


/* Find numerical root for phi with tolerance TOL to the equation 

prob = phi/(2 pi) + amp * [sin(phi + phase) - sin(phase)]

with phi in the range [0, 2 pi]. */
#ifdef __cplusplus
extern "C"
#endif
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





#ifdef __cplusplus
extern "C"
#endif
void setupgeneralGW(gwgeneralSrc *gw)
{
  long double deg2rad = M_PI/180.0;
  long double convert[3][3],trans[3][3];
  long double out[3][3];
  long double out_im[3][3];
  int i,j;
  long double st,ct,cp,sp,cpp,spp;

  st = sin(gw->theta_g);
  ct = cos(gw->theta_g);
  sp = sin(gw->phi_g);
  cp = cos(gw->phi_g);
  cpp = cos(gw->phi_polar_g);
  spp = sin(gw->phi_polar_g);

  gw->kg[0] = st*cp;
  gw->kg[1] = st*sp;
  gw->kg[2] = ct;

  /* Set up the GW strain: do not assume general relativity */
  gw->h[0][0] = gw->aplus_g+gw->ast_g;
  gw->h[0][1] = gw->across_g;
  gw->h[1][0] = gw->across_g;
  gw->h[1][1] = -gw->aplus_g+gw->ast_g+gw->asl_g;
  gw->h[2][0] = gw->avx_g;
  gw->h[2][1] = gw->avy_g;
  gw->h[2][2] = gw->asl_g;
  gw->h[0][2] = gw->avx_g;
  gw->h[1][2] = gw->avy_g;
  
  gw->h_im[0][0] = gw->aplus_im_g+gw->ast_im_g; // ALL SHOULD HAVE sqrt(-1) *  out the front
  gw->h_im[0][1] = gw->across_im_g; // only used for calcuateResidualGW and setupGW
  gw->h_im[1][0] = gw->across_im_g;
  gw->h_im[1][1] = -gw->aplus_im_g+gw->ast_im_g+gw->asl_im_g;
  gw->h_im[2][0] = gw->avx_im_g;
  gw->h_im[2][1] = gw->avy_im_g;
  gw->h_im[2][2] = gw->asl_im_g;
  gw->h_im[0][2] = gw->avx_im_g;
  gw->h_im[1][2] = gw->avy_im_g;
 
  /*printf("%10.8Le %10.8Le %10.8Le\n",gw->aplus_g,gw->ast_g,gw->asl_g);
  for (int jj=0;jj<3;jj++) {
                for (int kk=0;kk<3;kk++) {
                        printf("%10.8Le ",gw->h[jj][kk]);
                }
        }
        printf("\n");*/

  /* Now carry out a coordinate conversion (to convert to same coordinates 
     as the pulsar and GW source coordinates?)                                */

  /* These are from Euler's angles with pitch-roll-yaw convention */
  /* Also different definition of theta in MathWorld (mw) => sin(theta_mw) = -cos(theta_g),
     cos(theta_mw) = sin(theta_g) */
  /* Should check as there seems to be a minus sign difference in the top line
     with http://mathworld.wolfram.com/EulerAngles.html. Also the top line
     of this matrix is the bottom line of the wolfram matrix.  Probably doesn't
     matter, but should check */

  convert[0][0] = cp*ct*cpp-sp*spp; 
  convert[0][1] = sp*ct*cpp+cp*spp;
  convert[0][2] = -st*cpp;
       
  convert[1][0]= -cp*ct*spp-sp*cpp;
  convert[1][1]= -sp*ct*spp+cp*cpp;
  convert[1][2]= st*spp;
  
  convert[2][0]= cp*st;
  convert[2][1]= sp*st;
  convert[2][2]= ct;

  /* Calculate h.convert */
  matrixMult(gw->h,convert,out);
  matrixMult(gw->h_im,convert,out_im);

  /* Calculate transpose of 'convert' */
  for (i=0;i<3;i++)
    {
      for (j=0;j<3;j++)
	{
	  //	  printf("%g (%g, %g) ",(double)convert[j][i],(double)gw->h[j][i],(double)out[j][i]);
	  trans[i][j] = convert[j][i];
	}
      //      printf("\n");
    }
  /* Pre-multiply "out" by transpose */
  /* why multiply here -- we have already rotated it???? */
  matrixMult(trans,out,gw->h);
  matrixMult(trans,out_im,gw->h_im);
  //  matrixMult(trans,convert,out);
  /*  for (i=0;i<3;i++)
    {
      for (j=0;j<3;j++)
	{
	  printf("%g ",(double)gw->h[j][i]);
	}
      printf("\n");
      }*/

}






/* Produce a background of gravitational waves */
#ifdef __cplusplus
extern "C"
#endif
void GWgeneralbackground(gwgeneralSrc *gw,int *numberGW,long *idum,long double flo,long double fhi,
		  gwgenSpec gwAmps, int loglin)
{
  int j,k,l;
  k=0;
  for (j=0;j<4;j++) {
  	for (l=0;l<numberGW[j];l++)
    	{
      	gw[k].theta_g     = acos((TKranDev(idum)-0.5)*2);  
      	gw[k].phi_g       = TKranDev(idum)*2*M_PI;   
      	gw[k].phi_polar_g = 0.0; 
      	gw[k].phase_g     = TKranDev(idum)*2*M_PI;   
      	if (loglin==1)  /* Use equal sampling in log */
		{
	  	gw[k].omega_g  = 2*M_PI*exp(log(flo)+TKranDev(idum)*(log(fhi/flo))); 
		if (j == 0) {
	  		gw[k].aplus_g  = gwAmps.tensor_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.tensor_alpha)/sqrt((double)(numberGW[j])/log(fhi/flo))*TKgaussDev(idum);
	  		gw[k].across_g = gwAmps.tensor_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.tensor_alpha)/sqrt((double)(numberGW[j])/log(fhi/flo))*TKgaussDev(idum);
	  		gw[k].ast_g=gw[k].asl_g=gw[k].avx_g=gw[k].avy_g=0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		} else if (j == 1) {
			gw[k].ast_g=gwAmps.st_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.st_alpha)/sqrt((double)(numberGW[j])/log(fhi/flo))*TKgaussDev(idum);
	  		gw[k].aplus_g=gw[k].across_g=gw[k].asl_g=gw[k].avx_g=gw[k].avy_g=0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		} else if (j == 2) {
			gw[k].asl_g=gwAmps.sl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.sl_alpha)/sqrt((double)(numberGW[j])/log(fhi/flo))*TKgaussDev(idum);
	  		gw[k].aplus_g=gw[k].across_g=gw[k].ast_g=gw[k].avx_g=gw[k].avy_g=0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		} else if (j == 3) {
			gw[k].avx_g=gwAmps.vl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.vl_alpha)/sqrt((double)(numberGW[j])/log(fhi/flo))*TKgaussDev(idum);
			gw[k].avy_g=gwAmps.vl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.vl_alpha)/sqrt((double)(numberGW[j])/log(fhi/flo))*TKgaussDev(idum);
	  		gw[k].ast_g=gw[k].asl_g=gw[k].aplus_g=gw[k].across_g=0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		}
	} 
      else
	{
	  gw[k].omega_g = 2*M_PI*(TKranDev(idum)*(fhi-flo)+flo); 
	  	if (j == 0) {
	  		gw[k].aplus_g  = gwAmps.tensor_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.tensor_alpha)*sqrt(1.0/(double)numberGW[j])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
	  		gw[k].across_g  = gwAmps.tensor_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.tensor_alpha)*sqrt(1.0/(double)numberGW[j])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
	  		gw[k].ast_g=gw[k].asl_g=gw[k].avx_g=gw[k].avy_g=0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		} else if (j == 1) {
			gw[k].ast_g=gwAmps.st_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.st_alpha)*sqrt(1.0/(double)numberGW[j])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
	  		gw[k].aplus_g=gw[k].across_g=gw[k].asl_g=gw[k].avx_g=gw[k].avy_g=0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		} else if (j == 2) {
			gw[k].asl_g=gwAmps.sl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.sl_alpha)*sqrt(1.0/(double)numberGW[j])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
	  		gw[k].aplus_g=gw[k].across_g=gw[k].ast_g=gw[k].avx_g=gw[k].avy_g=0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		} else if (j == 3) {
			gw[k].avx_g=gwAmps.vl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.vl_alpha)*sqrt(1.0/(double)numberGW[j])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
			gw[k].avy_g=gwAmps.vl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.vl_alpha)*sqrt(1.0/(double)numberGW[j])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
	  		gw[k].ast_g=gw[k].asl_g=gw[k].aplus_g=gw[k].across_g==0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		}
	}
      //      printf("omega = %g\n",(double)(1.0/(gw[k].omega_g/2.0/M_PI)/86400.0));
	k++;
    }
  }
}

/* Calculate the pulsar timing residual induced by the gravitational wave */
#ifdef __cplusplus
extern "C"
#endif
long double calculateResidualgeneralGW(long double *kp,gwgeneralSrc *gw,long double time,long double dist)
{
  long double cosMu;
  long double psrVal1,psrVal2,psrVal1_im;
  long double earthVal1;
  long double tempVal[3];
  long double tempVal_im[3];
  long double geo;
  long double residual;
  long double deg2rad = M_PI/180.0L;
  long double VC = 299792456.200L;
  int i,k;

  for (i=0;i<3;i++)
    {
      tempVal[i] = 0.0;
      tempVal_im[i] = 0.0;
      for (k=0;k<3;k++)
	{
	  tempVal[i]    += gw->h[i][k]*kp[k];
	  tempVal_im[i] += gw->h_im[i][k]*kp[k];
	  //printf("%10.8Le ",gw->h[i][k]);
	}
            //printf("tempVal[%d] = %g ",i,(double)tempVal[i]);
    }
   //printf("\n");
   // printf("kp = %g %g %g\n",(double)kp[0],(double)kp[1],(double)kp[2]);

		
  cosMu = dotProduct(kp,gw->kg);   /* Angle between pulsar and GW source */
  //    printf(" cosmu = %g\n",(double)cosMu);
  if ((1+cosMu)!=0) 
  {
    psrVal1    = 0.5L/(1.0L+cosMu)*dotProduct(kp,tempVal); 
    psrVal1_im = 0.5L/(1.0L+cosMu)*dotProduct(kp,tempVal_im);
    //    printf("psrval1 = %g %g %g\n",(double)psrVal1,(double)(1+cosMu),(double)(1-cosMu));
  }
  else psrVal1 = 0.0L;
  
  geo      = psrVal1;
  //  printf("psrVal = %g %g\n",(double)psrVal1,(double)psrVal1_im);
  earthVal1 = psrVal1*sinl(gw->phase_g+time*gw->omega_g)+ 
    psrVal1_im*cosl(gw->phase_g+time*gw->omega_g); // sin -> cos
  if (dist==0) /* No pulsar term */
    psrVal2=0.0L;
  else
    psrVal2 = psrVal1*sinl(gw->phase_g-(1+cosMu)*dist/VC*gw->omega_g+time*gw->omega_g) + 
      psrVal1_im*cosl(gw->phase_g-(1+cosMu)*dist/VC*gw->omega_g+time*gw->omega_g); 

  //printf("%10.8Le %10.8Le %10.8Le\n",tempVal[0],tempVal[1],tempVal[2]);
  //printf("%10.8Le %10.8Le %10.8Le %10.8Le %10.8Le\n",psrVal1, psrVal1_im,earthVal1,psrVal2,gw->omega_g);
  residual = (earthVal1-psrVal2)/gw->omega_g; 
  //    printf("GWsim %g %g %g %g %g %g %g %g %g %g %g %g %g\n",(double)residual,(double)earthVal1,(double)psrVal2,(double)psrVal1,(double)cosMu,(double)psrVal1_im,(double)kp[0],(double)kp[1],(double)kp[2],(double)cosMu,(double)gw->h_im[0][0],(double)gw->h_im[0][1],(double)gw->h_im[1][0]);
  return residual;
}


/* Produce an anisotropic background of gravitational waves with arbitrary polarisations */
#ifdef __cplusplus
extern "C"
#endif
void GWgeneralanisotropicbackground(gwgeneralSrc *gw,int *numberGW,long *idum,long double flo,long double fhi,
		  gwgenSpec gwAmps,int loglin, double ***harmlist, int *nharms)
{
  int k,polj,l;
  double *probs,min,norm,dcth,dphi,cth,phi,isocomp,*plm,sum;
  double thisprob;
  int highindx,lowindx,newindx,Ngridref;
  int Ngrid;
  Ngridref=gwsim_Ngrid;
  k=0;
  for (polj=0;polj<4;polj++) {
	if (numberGW[polj]) {
		Ngrid=Ngridref;
  		while (Ngrid*Ngrid < numberGW[polj])
			Ngrid*=2;
  		probs=(double *)malloc(2.*Ngrid*Ngrid*sizeof(double));
  		plm=(double *)malloc(nharms[polj]*sizeof(double));
  		norm=0.;
  		dcth=2./((double)Ngrid);
  		dphi=2.*M_PI/((double)Ngrid);
  		cth=-1.+0.5*dcth;
  		isocomp=0.;
  		for (int hi=0;hi<nharms[polj];hi++) {
			if (!(harmlist[polj][hi][0]) && !(harmlist[polj][hi][1]))
				isocomp=harmlist[polj][hi][2];
  		}
  		l=0;
  		for (int i=0;i<Ngrid;i++) {
			phi=0.5*dphi;
			for (int hi=0;hi<nharms[polj];hi++) {
				plm[hi]=sphharm((int)(harmlist[polj][hi][0]),abs((int)(harmlist[polj][hi][1])),cth);
				//printf("%lf\n",plm[hi]);
			}
			for (int j=0;j<Ngrid;j++) {
				probs[l]=0.;
				for (int hi=0;hi<nharms[polj];hi++) {
					if ((int)(harmlist[polj][hi][1]) >= 0) 
						probs[l]+=harmlist[polj][hi][2]*plm[hi]*cos(harmlist[polj][hi][1]*phi);
					else
						probs[l]+=harmlist[polj][hi][2]*plm[hi]*sin(-harmlist[polj][hi][1]*phi);
				}
				if (i+j) {
					if (probs[l]<min)
						min=probs[l];
				} else
					min=probs[l];
				norm+=probs[l];
				l++;
				phi+=dphi;
			}
			cth+=dcth;
  		}
  		free(plm);
  		if (min < 0.) {
			printf("Warning: distribution computed from harmonic file has negative probabilities at certain sky positions. Increasing isotropic component amplitude to %6.4lf to renormalise this to a non-negative pdf.\n",isocomp-min);
			norm+=Ngrid*Ngrid*(-min);
		  } else
			min=0.;
		  sum=0.;
		  for (int i=0;i<Ngrid*Ngrid;i++) {
			sum+=(probs[i]-min);
			probs[i]=sum/norm;
			//printf("%lf\n",probs[i]);
		  }
	}
	for (l=0;l<numberGW[polj];l++) {
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
      		//printf("OUTPUT %Lf %Lf %i %i %lf\n",cosl(gw[k].theta_g),gw[k].phi_g,(int)(highindx/Ngrid),(highindx%Ngrid),dphi);
      		gw[k].phi_polar_g = 0.0;  
      		gw[k].phase_g     = TKranDev(idum)*2*M_PI;   
		if (loglin==1)  /* Use equal sampling in log */
                {
                	gw[k].omega_g  = 2*M_PI*exp(log(flo)+TKranDev(idum)*(log(fhi/flo)));
                	if (polj == 0) {
                        	gw[k].aplus_g  = gwAmps.tensor_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.tensor_alpha)/sqrt((double)(numberGW[polj])/log(fhi/flo))*TKgaussDev(idum);
                        	gw[k].across_g = gwAmps.tensor_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.tensor_alpha)/sqrt((double)(numberGW[polj])/log(fhi/flo))*TKgaussDev(idum);
                        	gw[k].ast_g=gw[k].asl_g=gw[k].avx_g=gw[k].avy_g=(long double)0.0;
                        	gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=(long double)0.0;
                	} else if (polj == 1) {
                        	gw[k].ast_g=gwAmps.st_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.st_alpha)/sqrt((double)(numberGW[polj])/log(fhi/flo))*TKgaussDev(idum);
                        	gw[k].aplus_g=gw[k].across_g=gw[k].asl_g=gw[k].avx_g=gw[k].avy_g=(long double)0.0;
                        	gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=(long double)0.0;
                	} else if (polj == 2) {
                        	gw[k].asl_g=gwAmps.sl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.sl_alpha)/sqrt((double)(numberGW[polj])/log(fhi/flo))*TKgaussDev(idum);
                        	gw[k].aplus_g=gw[k].across_g=gw[k].ast_g=gw[k].avx_g=gw[k].avy_g=(long double)0.0;
                        	gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=(long double)0.0;
			//printf("%i %i %6.4Le %6.4Le\n",k,polj,gw[k].asl_g,gw[k].avx_g);
                	} else if (polj == 3) {
                        	gw[k].avx_g=gwAmps.vl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.vl_alpha)/sqrt((double)(numberGW[polj])/log(fhi/flo))*TKgaussDev(idum);
                        	gw[k].avy_g=gwAmps.vl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.vl_alpha)/sqrt((double)(numberGW[polj])/log(fhi/flo))*TKgaussDev(idum);
                        	gw[k].ast_g=gw[k].asl_g=gw[k].aplus_g=gw[k].across_g=(long double)0.0;
                        	gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=(long double)0.0;
			//printf("%i %i %6.4Le %6.4Le\n",k,polj,gw[k].avx_g,gw[k].avy_g);
                	}
			//printf("%i %i %6.4Le %6.4Le\n",k,polj,gw[k].aplus_g,gw[k].ast_g);
        	} else {
          		gw[k].omega_g = 2*M_PI*(TKranDev(idum)*(fhi-flo)+flo);
			if (polj == 0) {
				gw[k].aplus_g  = gwAmps.tensor_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.tensor_alpha)*sqrt(1.0/(double)numberGW[polj])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
	                        gw[k].across_g  = gwAmps.tensor_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.tensor_alpha)*sqrt(1.0/(double)numberGW[polj])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
                        	gw[k].ast_g=gw[k].asl_g=gw[k].avx_g=gw[k].avy_g=0.;
                        	gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
                	} else if (polj == 1) {
                        	gw[k].ast_g=gwAmps.st_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.st_alpha)*sqrt(1.0/(double)numberGW[polj])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
                        	gw[k].aplus_g=gw[k].across_g=gw[k].asl_g=gw[k].avx_g=gw[k].avy_g=0.;
                        	gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
                	} else if (polj == 2) {
                        	gw[k].asl_g=gwAmps.sl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.sl_alpha)*sqrt(1.0/(double)numberGW[polj])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
                        	gw[k].aplus_g=gw[k].across_g=gw[k].ast_g=gw[k].avx_g=gw[k].avy_g=0.;
                        	gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
                	} else if (polj == 3) {
                        	gw[k].avx_g=gwAmps.vl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.vl_alpha)*sqrt(1.0/(double)numberGW[polj])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
                        	gw[k].avy_g=gwAmps.vl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.vl_alpha)*sqrt(1.0/(double)numberGW[polj])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
                        	gw[k].ast_g=gw[k].asl_g=gw[k].aplus_g=gw[k].across_g==0.;
                        	gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
                	}
        	}
		k++;
            //printf("omega = %g\n",(double)(1.0/(gw[k].omega_g/2.0/M_PI)/86400.0));
    	}
	if (numberGW[polj])
    		free(probs);
    }
}



/* Produce an anisotropic background of gravitational waves */
#ifdef __cplusplus
extern "C"
#endif
void GWanisotropicbackground(gwSrc *gw,int numberGW,long *idum,long double flo,long double fhi,
		  double gwAmp,double alpha,int loglin, double **harmlist, int nharms)
{
  int k, Ngrid;
  Ngrid = gwsim_Ngrid;
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


/* Produce a dipole background of gravitational waves */
#ifdef __cplusplus
extern "C"
#endif
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










#ifdef __cplusplus
extern "C"
#endif
int GWgeneralbackground_read(gwgeneralSrc *gw, FILE *file, int ireal){
	char key[13];
	int nreal;
	int ngw,id,igw,i;
	const unsigned int gwsize = 17*sizeof(long double);

	for (i=0;i<=ireal;i++){
		fread(&id,sizeof(int),1,file);
		printf("%d\n",id);
		if(id==ireal)break;
		fread(&ngw,sizeof(int),1,file);
		fseek(file,ngw*gwsize,SEEK_CUR);
		if(feof(file)){
			fprintf(stderr,"Could not read enough file to find realiation required\n");
			return -1;
		}
	}

	fread(&ngw,sizeof(int),1,file);

	printf("Reading %d GW sources from real %d (%d)\n",ngw,id,ireal);
	
	for (igw=0;igw<ngw;igw++){
		fread(&(gw[igw]),gwsize,1,file);
	}

	return ngw;
}



void GWgeneralbackground_write(gwgeneralSrc *gw, FILE *file,int ngw, int ireal){
	int igw;
	const unsigned int gwsize = 17*sizeof(long double);
	fwrite(&ireal,sizeof(int),1,file);
	fwrite(&ngw,sizeof(int),1,file);
	for (igw=0;igw<ngw;igw++){
		fwrite(&(gw[igw]),gwsize,1,file);
	}
	fflush(file);
}
