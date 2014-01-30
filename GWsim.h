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
 *
 *
 *
 *
 * Routines for generating anisotropic GW backgrounds and general GW backgrounds
 * is written by Jonathan Gair and Steve Taylor. Integrated in tempo2 by Rutger
 * van Haasteren. These routines are used in the plugins:
 *
 * GWanisobkgrd, GWbkgrdfromfile, GWdipolebkgrd, GWdipolebkgrd,
 * GWgeneralanisobkgrd
 *
 */

#ifndef __Tempo2_GWsim_h

#define __Tempo2_GWsim_h
#ifndef __Tempo2_h
#include "tempo2.h"
#endif /* __Tempo2_h */


/* Have a structure to define a gravitational wave source */
typedef struct gwSrc
{
  long double theta_g;     /* Angle of source from "z"-direction (e.g. 90-declination)             */
  long double phi_g;       /* Azimuthal angle of source (e.g. right ascension)                     */
  long double omega_g;     /* Frequency of gravitational wave source (Hz) in observer's rest frame */   
  long double phi_polar_g; /* Polarization angle of the gravitational wave source w.r.t. the unit vector along te direction of increasing elevation*/
  long double phase_g;
  long double aplus_g;
  long double aplus_im_g; /* imagainary part of plus polarization */
  long double across_g;
  long double across_im_g; /*imag part of cross polariation */
  
  /*Important if the source is a supermassive blackhole binary. These are parameters of the binary. */
  
  long double phi_bin;    /* orientation of line of nodes on sky, recall line of nodes is intersection of plane of sky with plane of binary orbit. phi is
							 measured relative to the unit vector in the direction of increasing(?) right ascension */
  long double theta_bin;  /* orbital phase, assume time stationary */
  //long double chirp_mass; /* chirp mass of binary = (m_1+m_2)*[(m_1*m_2)/(m_1+m_2)^2]^(3/5) */
  long double inc_bin;    /* orbital inclination angle w.r.t. plane of the sky*/
  long double dist_bin;   /* proper distance to the binary system */

  /* Derived parameters                                                             */
  long double h[3][3];     /* The gravitational wave strain                         */
  long double h_im[3][3];  /* same as h but with im amplitudes. Ideally would have h = h + i*h_im */
  long double kg[3];       /* The gravitational wave unit vector                    */
}gwSrc;



#ifdef __cplusplus
extern "C" {
#endif



// Functions relating to Vikram Ravi's eccentric binary black holes
double Fe(double ec);
double dadt(double ec, double a, double m1, double m2);
double dedt(double ec, double a, double m1, double m2);
double dtdt(double ec, double t, double p);
double Rs(double m1); // in m

long double eccRes(pulsar *psr, int i, int *coalesceFlag, double *prev_p, double *prev_e, double *prev_a, double *prev_epoch, double *prev_theta);

long double eccResWithEnergy(pulsar *psr, int i, int *coalesceFlag, double *prev_p, double *prev_e, double *prev_a, double *prev_epoch, double *prev_theta, float *eOut);


/* Some useful functions */
long double matrixMult(long double m1[3][3],long double m2[3][3],long double out[3][3]);
long double dotProduct(long double *m1,long double *m2);




void setupGW(gwSrc *gw);

long double matrixMult(long double m1[3][3],long double m2[3][3],long double out[3][3]);

long double dotProduct(long double *m1,long double *m2);

void GWbackground(gwSrc *gw,int numberGW,long *idum,long double flo,long double fhi,
		  double gwAmp,double alpha,int loglin);

long double calculateResidualGW(long double *kp,gwSrc *gw,long double time,long double dist);

void setupPulsar_GWsim(long double ra_p,long double dec_p,long double *kp);

int GWbackground_read(gwSrc *gw, FILE *file, int ireal);

void GWbackground_write(gwSrc *gw, FILE *file,int ngw, int ireal);

long double eccRes(pulsar *psr, int i, int *coalesceFlag, double *prev_p, double *prev_e, double *prev_a, double *prev_epoch, double *prev_theta);

long double eccResWithEnergy(pulsar *psr, int i, int *coalesceFlag, double *prev_p, double *prev_e, double *prev_a, double *prev_epoch, double *prev_theta, float *eOut) ;

double Rs(double m1);

double dadt(double ec, double a, double m1, double m2);
double dedt(double ec, double a, double m1, double m2);
double dtdt(double ec, double t, double p);
double Fe(double ec);
double psrangle(double centre_long,double centre_lat,double psr_long,double psr_lat);

#ifdef __cplusplus
}
#endif



/* Now the anisotropy content */


/* Have a structure to define a gravitational wave source */
typedef struct gwgeneralSrc
{
  long double theta_g;     /* Angle of source from "z"-direction (e.g. 90-declination)             */
  long double phi_g;       /* Azimuthal angle of source (e.g. right ascension)                     */
  long double omega_g;     /* Frequency of gravitational wave source (Hz) in observer's rest frame */   
  long double phi_polar_g; /* Polarization angle of the gravitational wave source w.r.t. the unit vector along te direction of increasing elevation*/
  long double phase_g;
  long double aplus_g;
  long double aplus_im_g; /* imagainary part of plus polarization */
  long double across_g;
  long double across_im_g; /*imag part of cross polarization */
  long double ast_g; /* real part of scalar transverse polarization */
  long double ast_im_g; /*imag part of scalar transverse polarization */
  long double asl_g; /* real part of scalar longitudinal polarization */
  long double asl_im_g; /*imag part of scalar longitudinal polarization */
  long double avx_g; /* real part of X vector longitudinal polarization */
  long double avx_im_g; /*imag part of X vector longitudinal polarization */
  long double avy_g; /* real part of Y vector longitudinal polarization */
  long double avy_im_g; /*imag part of Y vector longitudinal polarization */
  
  /*Important if the source is a supermassive blackhole binary. These are parameters of the binary. */
  
  long double phi_bin;    /* orientation of line of nodes on sky, recall line of nodes is intersection of plane of sky with plane of binary orbit. phi is
							 measured relative to the unit vector in the direction of increasing(?) right ascension */
  long double theta_bin;  /* orbital phase, assume time stationary */
  //long double chirp_mass; /* chirp mass of binary = (m_1+m_2)*[(m_1*m_2)/(m_1+m_2)^2]^(3/5) */
  long double inc_bin;    /* orbital inclination angle w.r.t. plane of the sky*/
  long double dist_bin;   /* proper distance to the binary system */

  /* Derived parameters                                                             */
  long double h[3][3];     /* The gravitational wave strain                         */
  long double h_im[3][3];  /* same as h but with im amplitudes. Ideally would have h = h + i*h_im */
  long double kg[3];       /* The gravitational wave unit vector                    */
}gwgeneralSrc;

typedef struct gwgenSpec {
	double tensor_amp;
	double st_amp;
	double sl_amp;
	double vl_amp;
	double tensor_alpha;
	double st_alpha;
	double sl_alpha;
	double vl_alpha;
}gwgenSpec;




#ifdef __cplusplus
extern "C" {
#endif

double sphharm(int l, int m, double x);
double Findphi(double prob, double amp, double phase);


void setupgeneralGW(gwgeneralSrc *gw);

/* Produce a background of gravitational waves */
void GWgeneralbackground(gwgeneralSrc *gw,int *numberGW,long *idum,long double flo,long double fhi,
		  gwgenSpec gwAmps, int loglin);
void GWgeneralanisotropicbackground(gwgeneralSrc *gw,int *numberGW,long *idum,long double flo,long double fhi,
		  gwgenSpec gwAmps,int loglin, double ***harmlist, int *nharms);
void GWanisotropicbackground(gwSrc *gw,int numberGW,long *idum,long double flo,long double fhi,
		  double gwAmp,double alpha,int loglin, double **harmlist, int nharms);
void GWdipolebackground(gwSrc *gw,int numberGW,long *idum,long double flo,long double fhi,
		  double gwAmp,double alpha,int loglin, double *dipoleamps);

long double calculateResidualgeneralGW(long double *kp,gwgeneralSrc *gw,long double time,long double dist);

int GWgeneralbackground_read(gwgeneralSrc *gw, FILE *file, int ireal);
void GWgeneralbackground_write(gwgeneralSrc *gw, FILE *file,int ngw, int ireal);

#ifdef __cplusplus
}
#endif






#endif /* Defined __Tempo2_GWsim_h */

