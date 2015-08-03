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
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "T2toolkit.h"
#include "../unsupported_plugins/toasim/toasim.h"

void calculateAngularFactors(pulsar *psr);
double calculateD(pulsar *psr);

using namespace std;

#define MAX_CORR 3


void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,j,k,p;
  double globalParameter;
  const char *CVS_verNum = "$Revision: 1.1 $";
  long seed = TKsetSeed();
  double a1,a1e,a2,a2e,epoch,raj,decj;
  long double ep,ep0,ep1,eps;
  double ra,dec,d;
  toasim_header_t* read_header;
  char fname[MAX_CORR][MAX_FILELEN];
  int  ireal[MAX_CORR];
  int ncorr=0;
  FILE* file;
  FILE *fin;
  FILE *fout;
  FILE *output_fp;
  char outputName[128];
  char useName[128];
  char outputFile[128] = "result.dat";
  double cvma1a1;
  double cvma2a2;
  double cvma1a2;
  double c1;
  double c2;
  int    npts;
  int    nfit;
  int    doSimulation=1;
  char covarFuncFile[MAX_FILELEN];

  // Set up epoch range
   ep0 = 54000.0;
    ep1 = 55500.0;
    eps = 100.0;

    //ep0 = 54356.0; 
    //  ep1 = 54357.0;
    //  eps = 1;

  if (displayCVSversion == 1) CVSdisplayVersion((char *)"gwmStats_plug.C",(char *)"plugin",CVS_verNum);

  *npsr = 0;
  strcpy(covarFuncFile,"NULL");

  printf("Graphical Interface: gwmStats\n");
  printf("Author:              G. Hobbs\n");
  printf("CVS Version:         $Revision: 1.1 $\n");
  printf(" --- type 'h' for help information\n");

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[*npsr],argv[++i]); 
	  strcpy(timFile[*npsr],argv[++i]);
	  (*npsr)++;
	}
      else if (strcmp(argv[i],"-useReal")==0)
	doSimulation=0;
      else if (strcmp(argv[i],"-dcf")==0)
	strcpy(covarFuncFile,argv[++i]);
      else if (strcmp(argv[i],"-o")==0)
	strcpy(outputFile,argv[++i]);
      else if (strcmp(argv[i],"-corn")==0)
        {
	  strcpy(fname[ncorr],argv[++i]);
	  ireal[ncorr]=atoi(argv[++i]);
	  ncorr++;
	}
    }
  output_fp = fopen(outputFile,"w");
  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

  // Start with postfit residuals
  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      if (i==0) doFitAll(psr,*npsr,covarFuncFile);   /* Do the fitting     */
      else textOutput(psr,*npsr,globalParameter,0,0,0,(char *)"");  /* Display the output */
    }

  // Run multiple iterations
  printf("simulating residuals\n");
  if (!(fin = fopen("34grid.dat","r")))
    {
      printf("Unable to open grid file 34grid.dat\n");
      exit(1);
    }

  // Create a data set
  if (doSimulation==1)
    {
      for (p=0;p<*npsr;p++)
	{
	  for (j=0;j<ncorr;j++)
	    {
	      strcpy(useName,timFile[p]);
	      strcat(useName,".");
	      strcat(useName,fname[j]);
	      printf("Opening file >%s<\n",useName);
	      if (!(file = fopen(useName,"r"))){
		printf("ERROR: File could not be read: '%s'\n",useName);
	      } else {
		read_header = toasim_read_header(file);
		printf("timfile name = >%s<\n",read_header->timfile_name);
		if(strcmp(read_header->timfile_name,timFile[p])!=0){
		  fprintf(stderr,"\n\n*****************\nWARNING: .tim file name mismatch '%s' != '%s'\n*****************\n\n",read_header->timfile_name,timFile[p]);
		}
		toasim_corrections_t *read_corr= toasim_read_corrections(read_header,ireal[j],file);
		for (i=0;i<psr[p].nobs;i++)
		  {
		    if (j==0)
		      psr[p].obsn[i].residual = (double)(read_corr->offsets[i]);
		    else
		      psr[p].obsn[i].residual += (double)(read_corr->offsets[i]);
		  }
		fclose(file);
	      }
	  //	  exit(1);
	    }
	  
	}
    }  
  for (p=0;p<*npsr;p++)
    {
      // ADD IN STUFF FOR CREATEREALISATION to read correction files
      //      for (i=0;i<psr[p].nobs;i++)
      //	psr[p].obsn[i].residual = TKgaussDev(&seed)*psr[p].obsn[i].toaErr*1.0e-6;
      // Set up the fitting
      psr[p].param[param_gwm_amp].val[0] = 0;
      psr[p].param[param_gwm_amp].paramSet[0] = 1;
      psr[p].param[param_gwm_amp].fitFlag[0] = 2;
      
      psr[p].param[param_gwm_amp].val[1] = 0;
      psr[p].param[param_gwm_amp].paramSet[1] = 1;
      psr[p].param[param_gwm_amp].fitFlag[1] = 2;
      
    }
  
  // Plot the residuals
  for (p=0;p<*npsr;p++)
    {
      sprintf(outputName,"%s.res",psr[p].name);
      fout = fopen(outputName,"w");
      for (i=0;i<psr[p].nobs;i++)
	fprintf(fout,"%.15Lf %.5g %.5g\n",psr[p].obsn[i].sat,(double)psr[p].obsn[i].residual,(double)psr[p].obsn[i].toaErr*1e-6);
      fclose(fout);
    }

  printf(" ... complete\n");
  
  while (!feof(fin))
    {
      if (fscanf(fin,"%lf %lf",&ra,&dec)==2)
	{
	  for (ep=ep0;ep<=ep1;ep+=eps)
	    {
	      for (p=0;p<*npsr;p++)
		{
		  psr[p].gwm_epoch = ep;
		  psr[p].gwm_raj = ra;
		  psr[p].gwm_decj = dec;
		  
		  calculateAngularFactors(&psr[p]);
		}	      
	      printf("Fitting\n");
	      doFitAll(psr,*npsr,covarFuncFile);   
	      printf(" ... complete\n");
	      // Extract results
	      a1 = (double)psr[0].param[param_gwm_amp].val[0];
	      a1e = (double)psr[0].param[param_gwm_amp].err[0];
	      a2 = (double)psr[0].param[param_gwm_amp].val[1];
	      a2e = (double)psr[0].param[param_gwm_amp].err[1];
	      epoch = (double)psr[0].gwm_epoch;
	      raj = (double)psr[0].gwm_raj;
	      decj = (double)psr[0].gwm_decj;
	      d = calculateD(&psr[0]);
	      
	      printf("Results: %g %g %g %g %g %g %g %g\n",epoch,raj,decj,a1,a1e,a2,a2e,d);
 	      cvma1a1=psr[0].covar[0][0];
	      cvma2a2=psr[0].covar[1][1];
	      cvma1a2=psr[0].covar[0][1];
	      npts = psr[0].globalNoConstrain;
	      nfit = psr[0].globalNfit;

	      c1 = sqrt(psr[0].fitChisq/(npts-nfit));
	      c2 = psr[0].fitChisq;
	      fprintf(output_fp,"%g %g %g %g %g %g %g %g %g %g %g %g %d %d %g\n",raj,decj,epoch,a1,a1e,a2,a2e,cvma1a1,cvma2a2,cvma1a2,c1,c2,npts,nfit,d);
	      fflush(output_fp);
	      // Should restore all the parameters
	      for (p=0;p<*npsr;p++)
		{
		  for (i=0;i<MAX_PARAMS;i++)
		    {
		      for (k=0;k<psr[p].param[i].aSize;k++)
			{
			  if (psr[p].param[i].fitFlag[k] > 0)
			    psr[p].param[i].val[k] = psr[p].param[i].prefit[k];
			}
		    }
		}
	    }
	}
    }
  fclose(fin);

  fclose(output_fp);
  return 0;
}

void calculateAngularFactors(pulsar *psr)
{
  long double wi,t1,t2;
  long double dt,speriod,tt;
  int k;
  double m,c,ival;
  double kp_theta,kp_phi,kp_kg,p_plus,p_cross,gamma,omega_g;
  //	       double res_e,res_i;
  long double resp,resc,res_r,res_i;
  double theta_p,theta_g,phi_p,phi_g;
  double lambda_p,beta_p,lambda,beta;
  double n1,n2,n3;
  double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
  double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;
  double cosTheta;
  
  lambda_p = (double)psr->param[param_raj].val[0];
  beta_p   = (double)psr->param[param_decj].val[0];
  lambda   = psr->gwm_raj;
  beta     = psr->gwm_decj;

		   // Pulsar vector
  n1 = cosl(lambda_p)*cosl(beta_p);
  n2 = sinl(lambda_p)*cosl(beta_p);
  n3 = sinl(beta_p);
  
  cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
    sinl(beta)*sinl(beta_p);
  
  // From KJ's paper
  // Gravitational wave matrix
  
  // NOTE: This is for the plus terms.  For cross should use different terms
  e11p = pow(sinl(lambda),2)-pow(cosl(lambda),2)*pow(sinl(beta),2);
  e21p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
  e31p = cosl(lambda)*sinl(beta)*cosl(beta);
  
  e12p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
  e22p = pow(cosl(lambda),2)-pow(sinl(lambda),2)*pow(sinl(beta),2);
  e32p = sinl(lambda)*sinl(beta)*cosl(beta);
  
  e13p = cosl(lambda)*sinl(beta)*cosl(beta);
  e23p = sinl(lambda)*sinl(beta)*cosl(beta);
  e33p = -powl(cosl(beta),2);
  
  resp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
	  n2*(n1*e21p+n2*e22p+n3*e23p)+
	  n3*(n1*e31p+n2*e32p+n3*e33p));
  
  //		   printf("Resp = %s %g %g\n",psr[p].name,(double)resp,(double)cosTheta);
  if ((1-cosTheta)==0.0)
    resp = 0.0;  // Check if this is sensible
  else
    resp = 1.0L/(2.0L*(1.0L-cosTheta))*(resp); 
  
  psr->quad_ifunc_geom_p = resp;
  
  lambda   = psr->gwm_raj;
  beta     = psr->gwm_decj;
  
  
  e11c = sinl(2*lambda)*sinl(beta);
  e21c = -cosl(2*lambda)*sinl(beta);
  e31c = -sinl(lambda)*cosl(beta);
  
  e12c = -cosl(2*lambda)*sinl(beta);
  e22c = -sinl(2*lambda)*sinl(beta);
  e32c = cosl(lambda)*cosl(beta);
  
  e13c = -sinl(lambda)*cosl(beta);
  e23c = cosl(lambda)*cosl(beta);
  e33c  = 0;
  
  resc = (n1*(n1*e11c+n2*e12c+n3*e13c)+
	  n2*(n1*e21c+n2*e22c+n3*e23c)+
	  n3*(n1*e31c+n2*e32c+n3*e33c));
  //		   printf("Resc = %s %g %g\n",psr[p].name,(double)resc,(double)cosTheta);
  if ((1-cosTheta)==0.0)
    resc = 0.0;  // Check if this is sensible
  else
    resc = 1.0L/(2.0L*(1.0L-cosTheta))*(resc); 
  psr->quad_ifunc_geom_c = resc;
  

  
}

double calculateD(pulsar *psr)
{
  double c11,c22,c12,a1,a2,det,s11,s22,s12,d;
  double chisq;
  int npts,nfit;
  int i0,i1;

  // These should be checked if more than 1 global parameter is being fitted for
  i0 = 0;
  i1 = 1;

  c11 = psr->covar[i0][i0];
  c22 = psr->covar[i1][i1];
  c12 = psr->covar[i0][i1];

  a1 = (double)psr->param[param_gwm_amp].val[0];
  a2 = (double)psr->param[param_gwm_amp].val[1];
  npts = psr->globalNoConstrain;
  nfit = psr->globalNfit;

  printf("npts,nfit = %d %d\n",npts,nfit);
  chisq = psr->fitChisq; // CHECK THIS --- JINGBO HAS CHISQR

  c11 = c11*chisq/(double)(npts-nfit);
  c22 = c22*chisq/(double)(npts-nfit);
  c12 = c12*chisq/(double)(npts-nfit);

  det = c11*c22-c12*c12;
  s11 = c22/det;
  s22 = c11/det;
  s12 = -c12/det;
  d=a1*a1*s11+a2*a2*s22+2*a1*a2*s12;
  return d;

}

char * plugVersionCheck = (char *)TEMPO2_h_VER;
