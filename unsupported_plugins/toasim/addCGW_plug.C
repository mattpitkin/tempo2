//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards
//
// Plugin to read a set of arrival time files and produce a list of Gaussian random numbers based on the TOA uncertainties 


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
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "T2toolkit.h"
#include "tempo2.h"
#include "toasim.h"

double calculateOffset(pulsar *psr,int p,int obs,double cgw_freq,double cgw_ra,double cgw_dec,double cgw_h0,double cgw_epoch,double cgw_cosinc,double cgw_angpol,double cgw_mc);


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
	char parFile[MAX_PSR][MAX_FILELEN];
	char timFile[MAX_PSR][MAX_FILELEN];
	int i,nit,j,p;
	char fname[MAX_FILELEN];
	double globalParameter;
	long double result;
	long seed = TKsetSeed();

	//
	// For the output file
	//
	toasim_header_t* header;
	toasim_header_t* read_header;
	FILE* file;
	double offsets[MAX_OBSN]; // Will change to doubles - should use malloc
	// Create a set of corrections.
	toasim_corrections_t* corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));

	// User defined parameters
	double cgw_freq;  int set_freq=0;
	double cgw_ra;    int set_ra=0;
	double cgw_dec;   int set_dec=0;
	double cgw_h0;    int set_h0=0;
	double cgw_epoch; int set_epoch=0;
	double cgw_cosinc;int set_cosinc=0;
	double cgw_angpol;int set_angpol=0;
	double cgw_mc;    int set_mc=0;

	corr->offsets=offsets;
	corr->params=""; // Normally leave as NULL. Can store this along with each realisation. 
	// Same length string in every iteration - defined in r_param_length see below
	corr->a0=0; // constant
	corr->a1=0; // a1*x
	corr->a2=0; // a2*x*X

	*npsr = 0;
	nit = 1;

	printf("Graphical Interface: addCGW\n");
	printf("Author:              G. Hobbs\n");
	printf("Version:             1.0\n");

	/* Obtain all parameters from the command line */
	for (i=2;i<argc;i++)
	{

	  if (strcmp(argv[i],"-nreal")==0){ // Note deterministic signal so this just returns the same signal multiple times
	    nit=atoi(argv[++i]);
	  }
	  if (strcmp(argv[i],"-f")==0)
	    {
	      strcpy(parFile[*npsr],argv[++i]); 
	      strcpy(timFile[*npsr],argv[++i]);
	      (*npsr)++;
	    }
	  else if (strcmp(argv[i],"-freq")==0)
	    {
	      sscanf(argv[++i],"%lf",&cgw_freq); set_freq=1;
	    }
	  else if (strcmp(argv[i],"-ra")==0)
	    { sscanf(argv[++i],"%lf",&cgw_ra); set_ra=1;}
	  else if (strcmp(argv[i],"-dec")==0)
	    { sscanf(argv[++i],"%lf",&cgw_dec); set_dec=1;}
	  else if (strcmp(argv[i],"-h0")==0)
	    { sscanf(argv[++i],"%lf",&cgw_h0); set_h0=1;}
	  else if (strcmp(argv[i],"-ep")==0)
	    { sscanf(argv[++i],"%lf",&cgw_epoch); set_epoch=1;}
	  else if (strcmp(argv[i],"-cosinc")==0)
	    { sscanf(argv[++i],"%lf",&cgw_cosinc); set_cosinc=1;}
	  else if (strcmp(argv[i],"-angpol")==0)
	    { sscanf(argv[++i],"%lf",&cgw_angpol); set_angpol=1;}
	  else if (strcmp(argv[i],"-mc")==0) // Solar masses
	    { sscanf(argv[++i],"%lf",&cgw_mc); set_mc=1;}
	}
	if (set_freq==0 || set_ra==0 || set_dec==0 || set_h0==0 || set_epoch==0 || set_cosinc==0 || set_angpol==0){
	  printf("Have not set all required parameters:\n");
	  printf("-freq, -ra, -dec, -h0, -ep, -cosinc, -angpol\n");
	  exit(1);
	}

	readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
	// Now read in all the .tim files
	readTimfile(psr,timFile,*npsr); /* Load the arrival times    */

	preProcess(psr,*npsr,argc,argv);

	for (p=0;p<*npsr;p++)
	{
		printf("NTOA = %d\n",psr[p].nobs);
		header = toasim_init_header();
		strcpy(header->short_desc,"addCGW");
		strcpy(header->invocation,argv[0]);
		strcpy(header->timfile_name,timFile[p]);
		strcpy(header->parfile_name,"Unknown");
		header->idealised_toas="NotSet"; // What should this be
		header->orig_parfile="NA";
		header->gparam_desc=""; // Global parameters
		header->gparam_vals="";
		header->rparam_desc=""; // Desciprtion of the parameters
		header->rparam_len=0; // Size of the string
		header->seed = seed;

		header->ntoa = psr[p].nobs;
		header->nrealisations = nit;

		// First we write the header...
		sprintf(fname,"%s.addCGW",timFile[p]);
		file = toasim_write_header(header,fname);


		for (i=0;i<nit;i++)
		{
			if(i%10 == 0){
				printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
				printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
				printf("Iteration %d/%d",i+1,nit);
				fflush(stdout);
			}
			for (j=0;j<psr[p].nobs;j++){
			  offsets[j] = calculateOffset(psr,p,j,cgw_freq,cgw_ra,cgw_dec,cgw_h0,cgw_epoch,cgw_cosinc,cgw_angpol,cgw_mc);
			}
			toasim_write_corrections(corr,header,file);
		}
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("Iteration %d/%d\n",i,nit);

		printf("Close file\n");
		fclose(file);
	}
	return 0;
}

double calculateOffset(pulsar *psr,int p,int obs,double cgw_freq,double cgw_ra,double cgw_dec,double cgw_h0,double cgw_epoch,double cgw_cosinc,double cgw_angpol,double cgw_mc)
{
  double kp_theta,kp_phi,kp_kg,p_plus,p_cross,gamma,omega_g;
  //	       double res_e,res_i;
  long double resp,resc,res_r,res_i;
  double theta_p,theta_g,phi_p,phi_g;
  double lambda_p,beta_p,lambda,beta;
  long double time;	      
  double n1,n2,n3;
  double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
  double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;
  double cosTheta;

  time    = (psr[p].obsn[obs].sat - cgw_epoch)*86400.0L;

  lambda_p = (double)psr[p].param[param_raj].val[0];
  beta_p   = (double)psr[p].param[param_decj].val[0];
  lambda   = cgw_ra;
  beta     = cgw_dec;
  // Pulsar vector
  n1 = cosl(lambda_p)*cosl(beta_p);
  n2 = sinl(lambda_p)*cosl(beta_p);
  n3 = sinl(beta_p);
  cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
    sinl(beta)*sinl(beta_p);

  // From KJ's paper
  // Gravitational wave matrix
  e11p = pow(sinl(lambda),2)-pow(cosl(lambda),2)*pow(sinl(beta),2);
  e21p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
  e31p = cosl(lambda)*sinl(beta)*cosl(beta);
  //	       printf("ex1p = %g %g %g (%g %g)\n",e11p,e21p,e31p,lambda,beta);
  
  e12p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
  e22p = pow(cosl(lambda),2)-pow(sinl(lambda),2)*pow(sinl(beta),2);
  e32p = sinl(lambda)*sinl(beta)*cosl(beta);
  //	       printf("ex2p = %g %g %g\n",e12p,e22p,e32p);
  
  e13p = cosl(lambda)*sinl(beta)*cosl(beta);
  e23p = sinl(lambda)*sinl(beta)*cosl(beta);
  e33p = -powl(cosl(beta),2);
  //	       printf("ex3p = %g %g %g\n",e13p,e23p,e33p);
  //	       exit(1);
  
  
  omega_g = cgw_freq;
  
  resp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
	  n2*(n1*e21p+n2*e22p+n3*e23p)+
	  n3*(n1*e31p+n2*e32p+n3*e33p));
  //	       printf("resp = %Lg\n",resp);
  
  // Determine cross term
  e11c = sin(2*lambda)*sin(beta);
  e21c = -cos(2*lambda)*sin(beta);
  e31c = -sin(lambda)*cos(beta);
  
  e12c = -cos(2*lambda)*sin(beta);
  e22c = -sin(2*lambda)*sin(beta);
  e32c = cos(lambda)*cos(beta);
  
  e13c = -sin(lambda)*cos(beta);
  e23c = cos(lambda)*cos(beta);
  e33c  = 0;
  
  
  resc = (n1*(n1*e11c+n2*e12c+n3*e13c)+
	  n2*(n1*e21c+n2*e22c+n3*e23c)+
	  n3*(n1*e31c+n2*e32c+n3*e33c));
  

  res_r = (cgw_h0/omega_g*((1+pow(cgw_cosinc,2))*cos(2*cgw_angpol)*sin(omega_g*time)+2*cgw_cosinc*sin(2*cgw_angpol)*cos(omega_g*time)))*resp 
    + (cgw_h0/omega_g*((1+pow(cgw_cosinc,2))*sin(2*cgw_angpol)*sin(omega_g*time)-2*cgw_cosinc*cos(2*cgw_angpol)*cos(omega_g*time)))*resc; 
  //  printf("Have res_r = %g %g %g %g\n",(double)res_r,sin(omega_g*time),cos(omega_g*time),(double)time);
  res_i = 0.0;

  if (psr[p].gwsrc_psrdist>0) // Add in the pulsar term  (NOTE: using subtraction here)
    {
      double omega_prime_g;
      double h0_prime;
      
      if (cgw_mc == 0) {omega_prime_g = omega_g; h0_prime = cgw_h0;}
      else {
	//			 omega_prime_g = omega_g - 2*M_PI*2.77e-8*pow(psr[p].cgw_mc/1e8,5.0/3.0)*pow(omega_g/2.0/M_PI/1e-7,11.0/3.0)*(psr[p].gwsrc_psrdist/PCM/1000.0)*(1-cosTheta);
	
	omega_prime_g = 2*M_PI*pow((1-cosTheta)*psr[p].gwsrc_psrdist/SPEED_LIGHT*256.0/5.0/pow(SPEED_LIGHT,5)*pow(M_PI,8.0/3.0)*pow(GM*cgw_mc,5.0/3.0)+pow(omega_g/2.0/M_PI,-8.0/3.0),-3.0/8.0);

	h0_prime = cgw_h0*pow(omega_prime_g/omega_g,2.0/3.0);
      }
      res_r -= ((h0_prime/omega_prime_g*((1+pow(cgw_cosinc,2))*cos(2*cgw_angpol)*sin(omega_prime_g*time-(1-cosTheta)*psr[p].gwsrc_psrdist/SPEED_LIGHT*omega_prime_g)+2*cgw_cosinc*sin(2*cgw_angpol)*cos(omega_prime_g*time-(1-cosTheta)*psr[p].gwsrc_psrdist/SPEED_LIGHT*omega_prime_g)))*resp 
		+ (h0_prime/omega_prime_g*((1+pow(cgw_cosinc,2))*sin(2*cgw_angpol)*sin(omega_prime_g*time-(1-cosTheta)*psr[p].gwsrc_psrdist/SPEED_LIGHT*omega_prime_g)-2*cgw_cosinc*cos(2*cgw_angpol)*cos(omega_prime_g*time-(1-cosTheta)*psr[p].gwsrc_psrdist/SPEED_LIGHT*omega_prime_g)))*resc); 
    }
  
  if ((1-cosTheta)==0.0)
    {
      res_r = 0.0;
      res_i = 0.0;
    }
  else
    {
      res_r = -1.0L/(2.0L*(1.0L-cosTheta))*(res_r); 
      res_i = -1.0L/(2.0L*(1.0L-cosTheta))*(res_i); 
    }

  return -(res_r+res_i); // Xingjiang Zhu argued on 31st March 2014 that this should be negative
}


char * plugVersionCheck = TEMPO2_h_VER;
