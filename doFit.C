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
#include <string.h>
#include <dlfcn.h>
#include "tempo2.h"
#include "TKfit.h"
#include "constraints.h"



int getNparams(pulsar psr);
void formCholeskyMatrix(double *c,double *resx,double *resy,double *rese,int np, int nc,double **uinv);
double getConstraintDeriv(pulsar *psr,int ipos,int i,int k);

/* Main routines for fitting in TEMPO2               */
void doFit(pulsar *psr,int npsr,int writeModel) 
{
  int p,npol,i,j,okay,ip[MAX_OBSN],k;
  double *x,*y,*sig,*val,chisq;
  double *error;
  double tol = 1.0e-27;  /* Tolerence for singular value decomposition routine */
  double newStart=-1.0,newFinish=-1.0;
  const char *CVS_verNum = "$Revision$";

  if (displayCVSversion == 1) CVSdisplayVersion("doFit.C","doFit()",CVS_verNum);
  if (debugFlag==1) printf("Entering doFit\n");
  if (debugFlag==1) printf("Fitting with function: %s\n",psr[0].fitFunc);

  if (strcmp(psr[0].fitFunc,"default")!=0)
    {
      char *(*entry)(pulsar *,int,int);
      void * module;
      char str[100];
      char tempo2MachineType[MAX_FILELEN]="";

      printf("Calling fitting plugin: %s\n",psr[0].fitFunc);
      strcpy(tempo2MachineType, getenv("LOGIN_ARCH"));

      for (int iplug=0; iplug < tempo2_plug_path_len; iplug++){
	      sprintf(str,"%s/%s_fitFunc_%s_plug.t2",tempo2_plug_path[iplug],
			      psr[0].fitFunc,tempo2MachineType);
	      printf("Looking for %s\n",str);
	      module = dlopen(str, RTLD_NOW); 
	      if(module==NULL){	  
		      printf("dlerror() = %s\n",dlerror());
	      } else break;
      }

      module = dlopen(str, RTLD_NOW); 
      if(!module)  {
	fprintf(stderr, "[error]: dlopen() unable to open plugin: %s.\n",str);
	exit(1);
      }
      /*
       * Check that the plugin is compiled against the same version of tempo2.h
       */
      char ** pv  = (char**)dlsym(module, "plugVersionCheck");
      if(pv!=NULL){
	      // there is a version check for this plugin
	      if(strcmp(TEMPO2_h_VER,*pv)){
		      fprintf(stderr, "[error]: Plugin version mismatch\n");
		      fprintf(stderr, " '%s' != '%s'\n",TEMPO2_h_VER,*pv);
		      fprintf(stderr, " Please recompile plugin against same tempo2 version!\n");
		      dlclose(module);
		      exit(1);
	      }
      }


      entry = (char*(*)(pulsar *,int,int))dlsym(module, "pluginFitFunc");
      if( entry == NULL ) {
	dlclose(module);
	fprintf(stderr, "[error]: dlerror() failed while retrieving address.\n" ); 
	exit(1);
      }
      entry(psr,npsr,writeModel);
      printf("Returning\n");
      
      return;
    }

  for (p=0;p<npsr;p++)  /* Loop over all the pulsars */
    {

	/*
	 * Check for "broken" combinations of fit parameters
	 */
	if (psr->param[param_dmmodel].paramSet[0] &&
			psr->param[param_dm].fitFlag[0]!= 0 &&
			psr->param[param_dmmodel].linkTo[0] != param_dm){
		// if we are using DMMODEL, and we want to fit for DM then
		// DMMODEL must be linked to DM... otherwise badness will occur!
		printf("WARNING: DM cannot be fit with DMMODEL\n         unless you set 'DMMODEL DM 1' in par file.\n");
		psr->param[param_dm].fitFlag[0]=0;
	}


      strcpy(psr[p].rajStrPre,psr[p].rajStrPost);
      strcpy(psr[p].decjStrPre,psr[p].decjStrPost);
      /* How many parameters are we fitting for */
      npol = getNparams(psr[p]);
      x     = (double *)malloc((psr[p].nobs+psr[p].nconstraints)*sizeof(double)); // max fit data size is nobs+nconstraints
      y     = (double *)malloc((psr[p].nobs+psr[p].nconstraints)*sizeof(double));
      sig   = (double *)malloc((psr[p].nobs+psr[p].nconstraints)*sizeof(double));
      val   = (double *)malloc(npol*sizeof(double));
      error = (double *)malloc(npol*sizeof(double));
      int count;
      count=0;
      for (i=0;i<psr[p].nobs;i++)
	{	  
	  psr[p].obsn[i].prefitResidual = psr[p].obsn[i].residual;
	  if (psr[p].obsn[i].deleted==0)
	    {
	      okay=1;
	      
	      /* Check for START and FINISH flags */
	      if (psr[p].param[param_start].paramSet[0]==1 && psr[p].param[param_start].fitFlag[0]==1 &&
		  (psr[p].param[param_start].val[0] > psr[p].obsn[i].sat))
		okay=0;
	      if (psr[p].param[param_finish].paramSet[0]==1 && psr[p].param[param_finish].fitFlag[0]==1 &&
		  psr[p].param[param_finish].val[0] < psr[p].obsn[i].sat)
		okay=0;
	      if (okay==1)
		{
		  x[count]   = (double)(psr[p].obsn[i].bbat-psr[p].param[param_pepoch].val[0]);
		  if (count==0)
		    {
		      newStart = (double)psr[p].obsn[i].sat;
		      newFinish = (double)psr[p].obsn[i].sat;
		    }
		  else
		    {
		      if (newStart > psr[p].obsn[i].sat)  newStart = (double)psr[p].obsn[i].sat;
		      if (newFinish < psr[p].obsn[i].sat) newFinish = (double)psr[p].obsn[i].sat;
		    } 
		  y[count]   = (double)psr[p].obsn[i].prefitResidual;
		  ip[count]  = i;
		  if (psr[p].fitMode==0) sig[count] = 1.0; 
		  else sig[count] = psr[p].obsn[i].toaErr*1e-6; /* Error in seconds */
		  count++;
		}
	    }
	}
     // add constraints as extra pseudo observations
     // These point to non-existant observations after the nobs array
     // These are later caught by getParamDeriv.
     for (i=0; i < psr[p].nconstraints; i++){
	ip[count] = psr->nobs+i;
	x[count]=0;
	y[count]=0;
	sig[count]=1e-12;
	count++;
      }

      psr[p].nFit = count;
      psr[p].param[param_start].val[0] = newStart-0.001; 
      psr[p].param[param_finish].val[0] = newFinish+0.001;
      psr[p].param[param_start].paramSet[0] = 1;
      psr[p].param[param_finish].paramSet[0] = 1; 
      /* Do the fit */
      if (npol!=0) /* Are we actually  doing any fitting? */ 
	{ 

	  if (debugFlag==1) printf("Doing the fit\n");
	  TKleastSquares_svd_psr(x,y,sig,psr[p].nFit,val,error,npol,psr[p].covar,&chisq,
				 FITfuncs,psr[p].fitMode,&psr[p],tol,ip);
	  //	  svdfit(x,y,sig,psr[p].nFit,val,npol,u,v,w,&chisq,FITfuncs,&psr[p],tol,ip);
	  if (debugFlag==1) printf("Complete fit: chisq = %f\n",(double)chisq);
	  //
	  // Display the covariance matrix 
	  //	  for (i=0;i<npol;i++)
	  //	    {
	  //	      for (j=0;j<npol;j++)
	  //		printf("%+.4f ",psr[0].covar[i][j]/sqrt((psr[0].covar[j][j]*psr[0].covar[i][i])));
	  //	      printf("\n");
	  //	    }


	  psr[p].fitChisq = chisq; 
	  psr[p].fitNfree = psr[p].nFit-npol;
	  //	  printf("Chisq = %g, reduced chisq = %g\n",(double)psr[p].fitChisq,(double)(psr[p].fitChisq/psr[p].fitNfree));
	  /* Now update the parameters */
	  if (debugFlag==1) printf("Updating the parameters\n");
	  //	  for (i=0;i<npol;i++)
	  //	    printf("Fit values and errors are %g %g\n",val[i],error[i]);
	  updateParameters(psr,p,val,error);
	  if (debugFlag==1) printf("Completed updating the parameters\n");
	}
      /* Free the vectors and matrices */
      free(error);
      free(val);
      free(sig);
      free(y);      
      free(x);
      
      if (psr[p].param[param_track].paramSet[0]==1 && psr[p].param[param_track].val[0]==20)
	psr[p].param[param_track].val[0]=40;

      /* Update errors if using a bootstrap error analysis.
       * See bootmc.f in tempo1.  Also description in Numerical Recipes in C for
       * description of bootstrap techniques
       */

      if (psr[p].bootStrap > 0)
	{
	  printf("Calculating uncertainties on fitted parameters using a Monte-Carlo bootstrap method (%d)\n",psr[p].bootStrap);
	  bootstrap(psr,p,npsr);
	}
    }
  if (debugFlag==1) printf("Leaving doFit\n");
}

/* Fitting routine with input data covariance matrix */
void doFitDCM(pulsar *psr,char *dcmFile,char *covarFuncFile,int npsr,int writeModel) 
{
  int i,j,k;
  double **uinv;
  double whiteres[MAX_OBSN],sum;
  FILE *fin,*fout;
  char fname[100],temp[100];
  int p,npol,okay,ip[MAX_OBSN];
  double *x,*y,*sig,*val,chisq;
  double *error;
  double tol = 1.0e-40;  /* Tolerence for singular value decomposition routine */
  double newStart=-1.0,newFinish=-1.0;
  long double meanRes=0.0;

  clock_t clk;

  //  printf("WARNING: Switching weighting off for the fit\n");
  //  printf("WARNING: THE .TIM FILE MUST BE SORTED - not checked for\n");
  clk=clock();  
  printf("Tcheck: Starting Cholesky fit (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
  if (strcmp(psr[0].fitFunc,"default")!=0)
    {
      char *(*entry)(pulsar *,int,int);
      void * module;
      char str[100];
      char tempo2MachineType[MAX_FILELEN]="";

      printf("Calling fitting plugin: %s\n",psr[0].fitFunc);
      strcpy(tempo2MachineType, getenv("LOGIN_ARCH"));

      for (int iplug=0; iplug < tempo2_plug_path_len; iplug++){
	      sprintf(str,"%s/%s_fitFunc_%s_plug.t2",tempo2_plug_path[iplug],
			      psr[0].fitFunc,tempo2MachineType);
	      printf("Looking for %s\n",str);
	      module = dlopen(str, RTLD_NOW); 
	      if(module==NULL){	  
		      printf("dlerror() = %s\n",dlerror());
	      } else break;
      }

      module = dlopen(str, RTLD_NOW); 
      if(!module)  {
	fprintf(stderr, "[error]: dlopen() unable to open plugin: %s.\n",str);
	exit(1);
      }
      /*
       * Check that the plugin is compiled against the same version of tempo2.h
       */
      char ** pv  = (char**)dlsym(module, "plugVersionCheck");
      if(pv!=NULL){
	      // there is a version check for this plugin
	      if(strcmp(TEMPO2_h_VER,*pv)){
		      fprintf(stderr, "[error]: Plugin version mismatch\n");
		      fprintf(stderr, " '%s' != '%s'\n",TEMPO2_h_VER,*pv);
		      fprintf(stderr, " Please recompile plugin against same tempo2 version!\n");
		      dlclose(module);
		      exit(1);
	      }
      }


      entry = (char*(*)(pulsar *,int,int))dlsym(module, "pluginFitFunc");
      if( entry == NULL ) {
	dlclose(module);
	fprintf(stderr, "[error]: dlerror() failed while retrieving address.\n" ); 
	exit(1);
      }
      entry(psr,npsr,writeModel);
      printf("Returning\n");
      
      return;
    }


  for (p=0;p<npsr;p++)  /* Loop over all the pulsars */
    {
      int nobs_and_constraints = psr[p].nobs + psr[p].nconstraints;
      printf("Tcheck: Processing pulsar %d\n",p);
      psr[p].fitMode = 0;

      /*
       * Check for "broken" combinations of fit parameters
       */
      if (psr->param[param_dmmodel].paramSet[0] &&
		      psr->param[param_dm].fitFlag[0]!= 0 &&
		      psr->param[param_dmmodel].linkTo[0] != param_dm){
	      // if we are using DMMODEL, and we want to fit for DM then
	      // DMMODEL must be linked to DM... otherwise badness will occur!
	      printf("WARNING: DM cannot be fit with DMMODEL\n         unless you set 'DMMODEL DM 1' in par file.\n");
	      psr->param[param_dm].fitFlag[0]=0;
      }

      strcpy(psr[p].rajStrPre,psr[p].rajStrPost);
      strcpy(psr[p].decjStrPre,psr[p].decjStrPost);
      /* How many parameters are we fitting for */
      printf("Tcheck: Determining which parameters we are fitting for  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
      npol = getNparams(psr[p]);
      x     = (double *)malloc(nobs_and_constraints*sizeof(double));
      y     = (double *)malloc(nobs_and_constraints*sizeof(double));
      sig   = (double *)malloc(nobs_and_constraints*sizeof(double));
      val   = (double *)malloc(npol*sizeof(double));
      error = (double *)malloc(npol*sizeof(double));
      int count;
      count=0;
      for (i=0;i<psr[p].nobs;i++)
	{	  
	  psr[p].obsn[i].prefitResidual = psr[p].obsn[i].residual;
	  if (psr[p].obsn[i].deleted==0)
	    {
	      okay=1;
	      
	      /* Check for START and FINISH flags */
	      if (psr[p].param[param_start].paramSet[0]==1 && psr[p].param[param_start].fitFlag[0]==1 &&
		  (psr[p].param[param_start].val[0] > psr[p].obsn[i].sat))
		okay=0;
	      if (psr[p].param[param_finish].paramSet[0]==1 && psr[p].param[param_finish].fitFlag[0]==1 &&
		  psr[p].param[param_finish].val[0] < psr[p].obsn[i].sat)
		okay=0;
	      if (okay==1)
		{
		  x[count]   = (double)(psr[p].obsn[i].bbat-psr[p].param[param_pepoch].val[0]);
		  if (count==0)
		    {
		      newStart = (double)psr[p].obsn[i].sat;
		      newFinish = (double)psr[p].obsn[i].sat;
		    }
		  else
		    {
		      if (newStart > psr[p].obsn[i].sat)  newStart = (double)psr[p].obsn[i].sat;
		      if (newFinish < psr[p].obsn[i].sat) newFinish = (double)psr[p].obsn[i].sat;
		    } 
		  y[count]   = (double)psr[p].obsn[i].prefitResidual;
		  ip[count]  = i;
		  sig[count] = psr[p].obsn[i].toaErr*1e-6; /* Error in seconds */
		  count++;
		}
	    }
	}
      // add constraints as extra pseudo observations
     // These point to non-existant observations after the nobs array
     // These are later caught by getParamDeriv.
     for (i=0; i < psr[p].nconstraints; i++){
	ip[count] = psr->nobs+i;
	x[count]=0;
	y[count]=0;
	sig[count]=1e-12;
	count++;
      }

     printf("Tcheck: Complete determining which parameters we are fitting for  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
      //      printf("Count = %d\n",count);
      psr[p].nFit = count;
      psr[p].param[param_start].val[0] = newStart-0.001; 
      psr[p].param[param_finish].val[0] = newFinish+0.001;
      psr[p].param[param_start].paramSet[0] = 1;
      psr[p].param[param_finish].paramSet[0] = 1; 
      
      printf("Tcheck: removing mean from the residuals??  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
      for (i=0;i<psr[p].nobs;i++)
	meanRes+=(long double)psr[p].obsn[i].residual;
      meanRes/=(long double)psr[p].nobs;
      for (i=0;i<psr[p].nobs;i++)
	psr[p].obsn[i].residual-=meanRes;
      printf("Tcheck: complete removing mean from the residuals??  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
      printf("Tcheck: allocating memory for uinv  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
      uinv = (double **)malloc(sizeof(double *)*nobs_and_constraints);

      for (i=0;i<nobs_and_constraints;i++)
	uinv[i] = (double *)malloc(sizeof(double)*nobs_and_constraints);
      printf("Tcheck: complete allocating memory for uinv  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
      // If we have the data covariance matrix on disk
      if (strcmp(dcmFile,"NULL")!=0)
	{
	  printf("Outputing dcm file\n");
	  sprintf(fname,"dcm_original_%d.dat",p+1);
	  fout = fopen(fname,"w");
	  for (i=0;i<psr[p].nobs;i++)
	    fprintf(fout,"%g %g %g\n",(double)(psr[p].obsn[i].sat-psr[p].param[param_pepoch].val[0]),(double)psr[p].obsn[i].residual,(double)psr[p].obsn[i].toaErr);
	  fclose(fout);
	  printf("Done outputing dcm file\n");
	  strcpy(fname,dcmFile);
	  if (npsr>1)
	    {
	      sprintf(temp,"%s_%d",fname,p+1);
	      strcpy(fname,temp);
	    }
	  //	  printf("Opening >%s<\n",fname);
	  if (!(fin = fopen(fname,"r")))
	    {
	      printf("Unable to open inverse cholesky matrix: %s\n",fname);
	      exit(1);
	    }
	  i=0;
	  j=0;
	  while (!feof(fin))
	    {
	      if (fscanf(fin,"%lf",&uinv[j][i])==1)
		{
		  i++;
		  if (i==psr[p].nobs)
		    {
		      i=0;
		      j++;
		      if (j==psr[p].nobs+1)
			{
			  printf("The matrix file is the wrong size - too large\n");
			  printf("N_obs = %d\n",psr[p].nobs);
			  exit(1);
			}
		    }
		}      
	    }
	  if (j!=psr[p].nobs && i!=0)
	    {
	      printf("The matrix file is the wrong size - too short\n");
	      printf("j = %d, i = %d, nobs = %d\n",j,i,psr[p].nobs);
	      exit(1);
	    }
	}
      else // Use data covariance function and calculate the covariance matrix
	{
	  int ndays = (int)(x[count-1-psr[p].nconstraints]-x[0])+2;
	  double covarFunc[ndays];
	  double escaleFactor = 1.0;
	  
	  printf("Tcheck: reading covariance function from disk  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
	  strcpy(fname,covarFuncFile);
	  if (npsr>1)
	    {
	      sprintf(temp,"%s_%d",fname,p+1);
	      strcpy(fname,temp);
	    }
	  //	  printf("Opening >%s<\n",fname);
	  if (!(fin = fopen(fname,"r")))
	    {
	      printf("Unable to open covariance function file: %s\n",fname);
	      exit(1);
	    }
	  if (debugFlag==1) printf("ndays = %d\n",ndays);
	  fscanf(fin,"%lf",&escaleFactor);
	  for (i=0;i<ndays;i++)
	    {
	      if (fscanf(fin,"%lf",&covarFunc[i])!=1)
		{
		  printf("ERROR: Incorrect number of days in the Cholesky matrix\n");
		  exit(1);
		}
	    }
	  fclose(fin);
	  //	  printf("Read covariance function\n");
	  //	  printf("WARNING: scaling all errors by: %g\n",escaleFactor);
	  for (i=0;i<count;i++)
	    sig[i]*=escaleFactor;
	  printf("Tcheck: complete reading covariance function from disk  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
	  // Form the data covariance matrix
	  printf("Tcheck: forming Cholesky matrix  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
	  formCholeskyMatrix(covarFunc,x,y,sig,count,psr[p].nconstraints,uinv);
	  printf("Tcheck: complete forming Cholesky matrix  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
	}
      printf("Tcheck: writing debug info to disk\n");
      sprintf(fname,"whitedata_%d.dat",p+1);
      fout = fopen(fname,"w");
      for (i=0;i<psr[p].nobs;i++)
	{
	  sum=0.0;
	  for (j=0;j<psr[p].nobs;j++)
	    sum+=uinv[j][i]*psr[p].obsn[j].residual;
	  whiteres[i] = sum;
	  //            fprintf(fout,"%g %g %g\n",(double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]),(double)psr[0].obsn[i].residual,(double)psr[0].obsn[i].toaErr);
	  fprintf(fout,"%g %g %g %g\n",(double)(psr[p].obsn[i].sat-psr[p].param[param_pepoch].val[0]),
		  whiteres[i],(double)psr[p].obsn[i].residual,(double)psr[p].obsn[i].toaErr);
	}
      fclose(fout);
      printf("Tcheck: complete writing debug info to disk (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);

      /* Do the fit */
      if (npol!=0) /* Are we actually  doing any fitting? */ 
	{ 
	  if (debugFlag==1) printf("Doing the fit\n");
	  printf("Tcheck: doing the fit  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
	  TKleastSquares_svd_psr_dcm(x,y,sig,psr[p].nFit,val,error,npol,psr[p].covar,&chisq,
				 FITfuncs,psr[p].fitMode,&psr[p],tol,ip,uinv);
	  printf("Tcheck: complete doing the fit  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
	  //	  svdfit(x,y,sig,psr[p].nFit,val,npol,u,v,w,&chisq,FITfuncs,&psr[p],tol,ip);
	  if (debugFlag==1) printf("Complete fit: chisq = %f\n",(double)chisq);
	  //	  printf("chisq = %g\n",chisq);
	  psr[p].fitChisq = chisq; 
	  psr[p].fitNfree = psr[p].nFit-npol;
	  
	  /* Now update the parameters */
	  if (debugFlag==1) printf("Updating the parameters\n");
	  printf("Tcheck: updating the parameter values  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
	  updateParameters(psr,p,val,error);
	  printf("Tcheck: complete updating the parameter values  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
	  if (debugFlag==1) printf("Completed updating the parameters\n");
	}    
      /* Free the vectors and matrices */
      free(error);
      free(val);
      free(sig);
      free(y);      
      free(x);
      
      if (psr[p].param[param_track].paramSet[0]==1 && psr[p].param[param_track].val[0]==20)
	psr[p].param[param_track].val[0]=40;

      /* Update errors if using a bootstrap error analysis.
       * See bootmc.f in tempo1.  Also description in Numerical Recipes in C for
       * description of bootstrap techniques
       */

      if (psr[p].bootStrap > 0)
	{
	  printf("Calculating uncertainties on fitted parameters using a Monte-Carlo bootstrap method (%d)\n",psr[p].bootStrap);
	  bootstrap(psr,p,npsr);
	}
	  printf("Tcheck: freeing memory  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
  for (i=0;i<psr[p].nobs;i++)
    free(uinv[i]);
  free(uinv);
	  printf("Tcheck: complete freeing memory  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
    }
}


int getNparams(pulsar psr)
{
  int npol;
  int i,k;

  npol = 1;
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr.param[i].aSize;k++)
	{
	  if (psr.param[i].fitFlag[k]==1) {
	    if (i!=param_start && i!=param_finish && i!=param_dmmodel && i!=param_gwsingle)
	      npol++;
	  }
	}
    }
  /* Add extra parameters for jumps */
  for (i=1;i<=psr.nJumps;i++)
    {
      if (psr.fitJump[i]==1)
	npol++;
    }
  /* Add extra parameters for sinusoidal whitening */
  if (psr.param[param_wave_om].fitFlag[0]==1)
    {
      printf("waveScale at this point = %d\n",psr.waveScale);
      if (psr.waveScale==1)
	npol+=psr.nWhite*2-1;
      else if (psr.waveScale==2)
	npol+=psr.nWhite*4-1;
      else
	npol+=psr.nWhite*2-1;      
    }
  if (psr.param[param_quad_om].fitFlag[0]==1)
    npol+=(psr.nQuad*4)-1;
  if (psr.param[param_ifunc].fitFlag[0]==1)
      npol+=(psr.ifuncN-1);
  if (psr.param[param_tel_dx].fitFlag[0]==1)
      npol+=(psr.nTelDX-1);
  if (psr.param[param_tel_dy].fitFlag[0]==1)
      npol+=(psr.nTelDY-1);
  if (psr.param[param_tel_dz].fitFlag[0]==1)
      npol+=(psr.nTelDZ-1);
  if (psr.param[param_quad_ifunc_p].fitFlag[0]==1)
      npol+=(psr.quad_ifuncN_p-1);
  if (psr.param[param_quad_ifunc_c].fitFlag[0]==1)
      npol+=(psr.quad_ifuncN_c-1);

  /* Add extra parameters for DMMODEL fitting */
  if (psr.param[param_dmmodel].fitFlag[0]==1)
    npol+=(int)((psr.dmoffsNum)*2); // *2 because we fit for the DM and constant offset
  /* Add extra parameters for GW single source fitting */
  if (psr.param[param_gwsingle].fitFlag[0]==1)
    npol+=4; 
  return npol;
}

void FITfuncs(double x,double afunc[],int ma,pulsar *psr,int ipos)
{
  int i,n=0,k,j,l,found;
  /*
   * IMPORTANT NOTICE: To allow for constraints, ipos may be larger than psr->nobs!
   * 
   * Therefore, if you modify this function, make sure to check (ipos < psr->nobs) before
   * using ipos, and call getParamDeriv to get the derivative since this function checks
   * for constraints and calls the appropriate constraint function/
   * 
   * M. Keith August 2011
   *
   */
  
  if(ipos < psr->nobs)afunc[n++] = 1;  /* Always fit for an arbitrary offset (unless this obs is a constraint!)*/
  else afunc[n++] = 0;
  /* See what we are fitting for */
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr->param[i].aSize;k++)
	{
	  if (psr->param[i].fitFlag[k]==1) /* If we are fitting for this parameter */
	    {
	      if (i!=param_start && i!=param_finish)
		{

		  if (debugFlag==1) printf("Fitting for %d (%s)\n",i,psr->param[i].label[k]);
		  if (i==param_wave_om)
		    {
		      if (psr->waveScale==2)
			{
			  // 			  for (j=0;j<psr->nWhite*2;j++) 			  
			  for (j=0;j<psr->nWhite*4;j++)
			    afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
			}
		      else
			{
 			  for (j=0;j<psr->nWhite*2;j++)
			    afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
			}
		    }
		  else if (i==param_quad_om)
		    {
		      for (j=0;j<psr->nQuad*4;j++)
			afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
		    }
		  else if (i==param_ifunc)
		    {
		      for (j=0;j<psr->ifuncN;j++)
			afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
		    }
		  else if (i==param_tel_dx)
		    {
		      for (j=0;j<psr->nTelDX;j++)
			afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
		    }
		  else if (i==param_tel_dy)
		    {
		      for (j=0;j<psr->nTelDY;j++)
			afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
		    }
		  else if (i==param_tel_dz)
		    {
		      for (j=0;j<psr->nTelDZ;j++)
			afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
		    }
		  else if (i==param_quad_ifunc_p)
		    {
		      for (j=0;j<psr->quad_ifuncN_p;j++)
			afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
		    }
		  else if (i==param_quad_ifunc_c)
		    {
		      for (j=0;j<psr->quad_ifuncN_c;j++)
			afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
		    }
		  else if (i==param_gwsingle)
		    {
		      afunc[n++] = getParamDeriv(psr,ipos,x,i,0);
		      afunc[n++] = getParamDeriv(psr,ipos,x,i,1);
		      afunc[n++] = getParamDeriv(psr,ipos,x,i,2);
		      afunc[n++] = getParamDeriv(psr,ipos,x,i,3);
		    }
		  else if (i==param_dmmodel)
		    {		      

		    double dmf = 1;
		    if (ipos < psr->nobs)dmf = 1.0/(DM_CONST*powl(psr->obsn[ipos].freqSSB/1.0e6,2));

		      for (j=0;j<(int)psr->dmoffsNum;j++)
			{
			  // This is for the actual frequency-dep fit
			  afunc[n++] = dmf*getParamDeriv(psr,ipos,x,i,j);
			  afunc[n++] =     getParamDeriv(psr,ipos,x,i,j+psr->dmoffsNum);
			}
		    }
		  else
		    afunc[n++] = getParamDeriv(psr,ipos,x,i,k);
		}
	    }
	}
    }
  /* JUMPS */
  
  for (i=1;i<=psr->nJumps;i++)
    {
      if (psr->fitJump[i]==1)
	{
	  found = 0;
	  if(ipos < psr->nobs){
	  for (l=0;l<psr->obsn[ipos].obsNjump;l++)
	    {
	      if (psr->obsn[ipos].jump[l]==i)
		{
		  found = 1;
		  break;
		}
	      //	      else
	      //		found = 0.0;
	    }
	  }
	  afunc[n++] = found;

	}
    } 
  if (n!=ma) { 
    printf("Problem in fitting routine n = %d, ma = %d\n",n,ma);
  }
}


double getParamDeriv(pulsar *psr,int ipos,double x,int i,int k)
{

  double afunc=0.0;
  /*
   * If we have ipos > nobs, then this is not a real observation but a 
   * constraint "pseudo observation". Therefore we call the constraint
   * function.
   *
   * Otherwise, we check what parameter it is and return the appropriate derivative.
   */
  //  printf("In here with %s %d %d %d\n",psr->name,psr->nobs,ipos,i);
  if(ipos >= psr->nobs)
    {
      //      printf("In here with ipos = %d, psr->nobs=%d, i=%d, str = %s, k= %d\n",ipos,psr->nobs,i,psr->param[i].shortlabel[0],k);
      afunc= getConstraintDeriv(psr,ipos-psr->nobs,i,k); // this is a constraint pseudo observation.
      //      printf("In here afunc = %g\n",afunc);
    }
//	  printf("MJK -- get deriv for constraint %d i=%d x=%lf k=%d af=%f\n",ipos-psr->nobs,i,x,k,afunc);
  else if (i==param_f)         /* Rotational frequency */
    {
      if (k==0)
	afunc = x*24.0L*3600.0L/psr->param[param_f].val[0];
      else if (k==1)    /* Rotational frequency derivative */
	afunc = 0.5L*x*x;
      else if (k==2)    /* Rotational frequency second derivative */
	afunc = 1.0L/6.0L*x*x*x/1.0e9L;
      else if (k==3)
	afunc = 1.0L/24.0L*x/1.0e18L*x*x*x;
      else if (k==4)
	afunc = 1.0L/120.0L*x*x*x*x*x/1.0e18L;
      else if (k==5)
	afunc = 1.0L/720.0L*powl(x,6.0L)/1.0e18L;
      else if (k==6)
	afunc = 1.0L/5040.0L*powl(x,7.0L)/1.0e18L;
      else if (k==7)
	afunc = 1.0L/40320.0L*powl(x,8.0L)/1.0e18L;
      else if (k==8)
	afunc = 1.0L/362880.0L*powl(x,9.0L)/1.0e18L;
      else if (k==9)
	afunc = 1.0L/3628800.0L*powl(x,10.0L)/1.0e18L;
      else if (k==10)
	afunc = 1.0L/3628800.0L/11.0L*powl(x,11.0L)/1.0e23L;
      else if (k==11)
	afunc = 1.0L/3628800.0L/11.0L/12.0L*powl(x,12.0L)/1.0e23L;
      else if (k==12)
	afunc = 1.0L/3628800.0L/11.0L/12.0L/13.0L*powl(x,13.0L)/1.0e23L;
    }
  else if (i==param_dshk)
    {
      longdouble kpc2m = 3.08568025e19L;           /* 1 kpc in m        */
      longdouble mas_yr2rad_s = 1.536281850e-16L;  /* 1 mas/yr in rad/s */
      longdouble t0;
      t0 = ((x + psr->param[param_pepoch].val[0])
	    - psr->param[param_posepoch].val[0])*SECDAY; 
      afunc = t0*t0/2.0L/SPEED_LIGHT*(psr->param[param_pmra].val[0]*psr->param[param_pmra].val[0]*
				     mas_yr2rad_s*mas_yr2rad_s+
				     psr->param[param_pmdec].val[0]*psr->param[param_pmdec].val[0]*
				     mas_yr2rad_s*mas_yr2rad_s)*kpc2m;
    }
  else if (i==param_glph)
    {	      
      if (psr->obsn[ipos].bbat >= psr->param[param_glep].val[k])
	afunc = 1.0/psr->param[param_f].val[0];
      else
	afunc = 0.0;
    }
  else if (i==param_glf0d)
    {
      longdouble dt1,expf,tp,tgl;
      
      tp = (psr->obsn[ipos].bbat-psr->param[param_pepoch].val[0])*86400.0;
      tgl = (psr->param[param_glep].val[k] - psr->param[param_pepoch].val[0])*86400.0;
      
      dt1 = tp-tgl;
      
      if (psr->param[param_gltd].val[k]!=0.0)
	expf = exp(-dt1/86400.0/psr->param[param_gltd].val[k]);
      else
	expf = 1.0;
      
      if (psr->obsn[ipos].bbat >= psr->param[param_glep].val[k])
	{
	  afunc = psr->param[param_gltd].val[k]*SECDAY*(1.0-expf)/psr->param[param_f].val[0]; ///psr->param[param_f].val[0];
	  //	  printf("Glitch diff = %d %.10f %.10f %.10f\n",k+1,afunc,(double)tp,(double)tgl,(double)psr->param[param_gltd].val[0]);

	}
      else
	afunc = 0.0;
    }
  else if (i==param_gltd)
    {
      longdouble dt1,expf,tp,tgl;
      
      tp = (psr->obsn[ipos].bbat-psr->param[param_pepoch].val[0])*86400.0L;
      tgl = (psr->param[param_glep].val[k] - psr->param[param_pepoch].val[0])*86400.0L;
      
      dt1 = tp-tgl;
      
      if (psr->param[param_gltd].val[k]!=0.0)
	expf = exp(-dt1/86400.0L/psr->param[param_gltd].val[k]);
      else
	expf = 1.0;
      
      if (psr->obsn[ipos].bbat >= psr->param[param_glep].val[k])
	afunc = psr->param[param_glf0d].val[k]*
	  (1.0-(1.0+dt1/SECDAY/(psr->param[param_gltd].val[k]))*expf)/psr->param[param_f].val[0]*SECDAY;
      else
	afunc = 0.0;
    }
  else if (i==param_glf0)
    {
      if (psr->obsn[ipos].bbat >= psr->param[param_glep].val[k])	
	  afunc = (psr->obsn[ipos].bbat-psr->param[param_glep].val[k])*86400.0/psr->param[param_f].val[0];	
      else
	afunc = 0.0;	      
    }
  else if (i==param_glf1)
    {
      if (psr->obsn[ipos].bbat >= psr->param[param_glep].val[k])
	afunc = 0.5*pow((psr->obsn[ipos].bbat-psr->param[param_glep].val[k])*86400.0,2)/psr->param[param_f].val[0];
      else
	afunc = 0.0;	      
    }
  else if (i==param_glf2)
    {
      if (psr->obsn[ipos].bbat >= psr->param[param_glep].val[k])
	{
	  afunc = (double)(1.0L/6.0L*powl((psr->obsn[ipos].bbat-psr->param[param_glep].val[k])*86400.0L/1.0e9,3)/psr->param[param_f].val[0]);
	  //	  printf("Fit = %.15g\n",afunc);
	}
      else
	afunc = 0.0;	      
    }
  else if (i==param_dm)    /* Dispersion measure */
    {
      double yrs;
      /* What about Doppler effect for the frequency -- change to barycentre?? */
      /* look at Blanford(?) paper */
      /* Should have a check to see if only one frequency exists in data
	 in which case fitting for DM does not make sense */
      if (k==0) 
	afunc = 1.0/(DM_CONST*powl(psr->obsn[ipos].freqSSB/1.0e6,2));
      else
	{
	  yrs = (psr->obsn[ipos].sat - psr->param[param_pepoch].val[0])/365.25;
	  afunc = 1.0/(DM_CONST*pow(psr->obsn[ipos].freqSSB/1.0e6,2))*pow(yrs,k);
	}
    }
  else if (i==param_dmx)
    {
        if ((psr->obsn[ipos].sat > psr->param[param_dmxr1].val[k])
                && (psr->obsn[ipos].sat < psr->param[param_dmxr2].val[k]))
	  afunc = 1.0/(DM_CONST*powl(psr->obsn[ipos].freqSSB/1.0e6,2));
        else 
          afunc = 0.0;
    }
  else if (i==param_fddc)    /* Frequency dependent term */
    afunc = 1.0/(pow(psr->obsn[ipos].freqSSB/1.0e6,psr->param[param_fddi].val[0]));
  else if (i==param_fddi)    /* Frequency dependent index */
    {
      printf("ERROR: Cannot currently fit for FDDI with a linear least-squares fit\n");
      exit(1);
    }
  /*	  else if (i==param_jump)
	  {
	  if (ipos>5)
	  {
	  printf("HERE with %.14Lf\n",psr->obsn[ipos].sat);
	  afunc = 1.0;
	  }
	  else
	  afunc = 0.0;
	  } */
  else if (i==param_telx)
    {
      long double dt,arg;
      if (psr->param[param_telEpoch].paramSet[0]==1)
	dt = (x + psr->param[param_pepoch].val[0]) - psr->param[param_telEpoch].val[0];
      else
	dt = 0.0;
      dt *= 86400.0;

      printf("k = %d\n",k);
      if (strcmp(psr->obsn[ipos].telID,"STL_FBAT")==0)
	{
	  if (k==0) afunc = psr->posPulsar[0];
	  arg = dt; if (k==1) afunc = psr->posPulsar[0]*arg;
	  arg *= dt; if (k==2) afunc = 0.5*psr->posPulsar[0]*arg;
	  arg *= dt; if (k==3) afunc = 1.0/6.0*psr->posPulsar[0]*arg;
	}
      else
	afunc = 0;
    }
  else if (i==param_tely)
    {
      long double dt,arg;
      if (psr->param[param_telEpoch].paramSet[0]==1)
	dt = (x + psr->param[param_pepoch].val[0]) - psr->param[param_telEpoch].val[0];
      else
	dt = 0.0;
      dt*=86400.0L;

      printf("k = %d\n",k);
      if (strcmp(psr->obsn[ipos].telID,"STL_FBAT")==0)
	{
	  if (k==0) afunc = psr->posPulsar[1];
	  arg = dt; if (k==1) afunc = psr->posPulsar[1]*arg;
	  arg *= dt; if (k==2) afunc = 0.5*psr->posPulsar[1]*arg;
	  arg *= dt; if (k==3) afunc = 1.0/6.0*psr->posPulsar[1]*arg;
	}
      else
	afunc = 0;
    }
  else if (i==param_telz)
    {
      long double dt,arg;
      if (psr->param[param_telEpoch].paramSet[0]==1)
	dt = (x + psr->param[param_pepoch].val[0]) - psr->param[param_telEpoch].val[0];
      else
	dt = 0.0;
      dt *= 86400.0L;

      printf("k = %d\n",k);
      if (strcmp(psr->obsn[ipos].telID,"STL_FBAT")==0)
	{
	  if (k==0) afunc = psr->posPulsar[2];
	  arg = dt; if (k==1) afunc = psr->posPulsar[2]*arg;
	  arg *= dt; if (k==2) afunc = 0.5*psr->posPulsar[2]*arg;
	  arg *= dt; if (k==3) afunc = 1.0/6.0*psr->posPulsar[2]*arg;
	}
      else
	afunc = 0;
    }
  else if (i==param_raj || i==param_decj || i==param_pmra || i==param_pmdec ||
	   i==param_px || i==param_pmrv)
    {
      double rce[3],re,deltae,alphae,psrra,psrdec,axy,s;
      int l;
      
      /* What about observatory to Earth ???? */
      /* Calculated centre of Earth from SSB */
      for (l=0;l<3;l++) {
	rce[l] = psr->obsn[ipos].earth_ssb[l];
	/*		rce[k] = -psr->obsn[ipos].earthMoonBary_earth[k] + 
			psr->obsn[ipos].earthMoonBary_ssb[k]; */
      }
      re = sqrt(dotproduct(rce,rce));
      
      /* Calculate position of Earth w.r.t SSB */
      axy = rce[2]/AULTSC;
      s = axy/(re/AULTSC);  /* Why is this AULT and not AULTSC in TEMPO ??? */
      deltae = atan2(s,sqrt(1.0-s*s));
      alphae = atan2(rce[1],rce[0]);
      
      /* Calculate position of pulsar w.r.t SSB */
      /* IS THIS JUST THE RA AND DEC OF THE PULSAR? */
      psrra  = psr->param[param_raj].val[0];
      psrdec = psr->param[param_decj].val[0];
      if (i==param_raj)
	{
	  afunc = re*cos(deltae)*cos(psrdec)*sin(psrra - alphae);
	}
      else if (i==param_decj)
	{
	afunc = re*(cos(deltae)*sin(psrdec)*cos(psrra - alphae) - sin(deltae)*cos(psrdec));
	//	  printf("dofit: %g %.10g 0.0\n",(double)(psr->obsn[ipos].sat-psr->param[param_pepoch].val[0]),(double)afunc);

	}
      else if (i==param_pmra)
	afunc = re*cos(deltae)*cos(psrdec)*sin(psrra - alphae) * x;
      else if (i==param_pmdec) /* pmdec */
	afunc = re*(cos(deltae)*sin(psrdec)*cos(psrra - alphae) - sin(deltae)*cos(psrdec))*x;
      else if (i==param_px)
	{
	  int l;
	  double pxConv = 1.74532925199432958E-2/3600.0e3,rca[3],rr,rcos1;
	  for (l=0;l<3;l++)
	    rca[l] = psr->obsn[ipos].earth_ssb[l] + psr->obsn[ipos].observatory_earth[l];
	  /*		    rca[k] = psr->obsn[ipos].earthMoonBary_ssb[k] -
			    psr->obsn[ipos].earthMoonBary_earth[k] + 
			    psr->obsn[ipos].observatory_earth[k]; */
	  rr    = dotproduct(rca,rca);
	  rcos1 = dotproduct(psr->posPulsar,rca);
	  afunc = 0.5*pxConv*(rr-rcos1*rcos1)/AULTSC;
	  if (debugFlag==1) printf("output fitting %g %g\n",(double)psr->obsn[ipos].bat,(double)afunc);
	  /* Now consider adding the effects of other distance determinations */
	  if (psr->param[i].nLinkFrom == 1 &&
	      psr->param[i].linkFrom[0] == param_dshk)
	    {
	      longdouble kpc2m = 3.08568025e19L;           /* 1 kpc in m        */
	      longdouble mas_yr2rad_s = 1.536281850e-16L;  /* 1 mas/yr in rad/s */
	      longdouble t0,afunc2;

	      t0 = ((x + psr->param[param_pepoch].val[0])
		    - psr->param[param_posepoch].val[0])*SECDAY; 
	      afunc2 = (t0*t0/2.0L/SPEED_LIGHT*
		(psr->param[param_pmra].val[0]*psr->param[param_pmra].val[0]*
		 mas_yr2rad_s*mas_yr2rad_s+
		 psr->param[param_pmdec].val[0]*psr->param[param_pmdec].val[0]*
		 mas_yr2rad_s*mas_yr2rad_s)*kpc2m);	     
	      afunc += (-afunc2/psr->param[param_px].val[0]/psr->param[param_px].val[0]);
	    }
	}
      else if (i==param_pmrv)
	{
	  int j;
	  double delt,etat,dt_SSB,pmtrans_rcos2,rca[4];
	  
	  for (j=0;j<3;j++)
	    rca[j] = psr->obsn[ipos].earth_ssb[j] + psr->obsn[ipos].observatory_earth[j];
	  /*		    
			   rca[j] = psr->obsn[ipos].earthMoonBary_ssb[j] -
			   psr->obsn[ipos].earthMoonBary_earth[j] + psr->obsn[ipos].observatory_earth[j]; */
	  
	  pmtrans_rcos2 = dotproduct(psr->velPulsar,rca);      
	  dt_SSB = 0.0; /* What should this be? */
	  etat = 32.149618300000;/* NOT TRUE WHAT SHOULD THIS BE? (see dm_delays.c) */
	  delt = (psr->obsn[ipos].sat-psr->param[param_posepoch].val[0] + (etat+dt_SSB)/SECDAY)/36525.0; 
	  afunc = pow(delt,2)*pmtrans_rcos2;
	}
    }
  else if (i==param_start)
    {
    }
  else if (i==param_finish)
    {
    }
  else if (i==param_wave_om) /* Whitening procedure using sinusoids */
    {
      double      om    = psr->param[param_wave_om].val[0];
      if (psr->waveScale==0)
	{
	  if (k%2==0) afunc = cos(om*(floor(k/2.0)+1)*x); 
	  else        afunc = sin(om*(floor(k/2.0)+1)*x); 
	  //	  printf("Value = %d %f %f %f %g\n",k,floor(k/2.0)+1,x,om,afunc);

	}
      else if (psr->waveScale==1)
	{
	  double freq = psr->obsn[ipos].freqSSB/1.0e6;
	  if (k%2==0) afunc = cos(om*(floor(k/2.0)+1)*x)/(DM_CONST*freq*freq); 
	  else afunc = sin(om*(floor(k/2.0)+1)*x)/(DM_CONST*freq*freq); 
	}
      else if (psr->waveScale==2)
	{
	  double freq = psr->obsn[ipos].freqSSB/1.0e6;
	  //	  printf("Have %d %d %d %d %g\n",k,psr->nWhite*2,k%2,(k-(psr->nWhite*2))%2,floor((k-(psr->nWhite*2))/2.0)+1);
	  if (k < psr->nWhite*2)
	    {
	      if (k%2==0) afunc = cos(om*(floor(k/2.0)+1)*x)/(DM_CONST*freq*freq); 
	      else afunc = sin(om*(floor(k/2.0)+1)*x)/(DM_CONST*freq*freq); 
	    }
	  else
	    {
	      if ((k-(psr->nWhite*2))%2==0) afunc = cos(om*(floor((k-(psr->nWhite*2))/2.0)+1)*x);
	      else afunc = sin(om*(floor((k-(psr->nWhite*2))/2.0)+1)*x);
	    }
	}
    }
  else if (i==param_tel_dx)
    {
      double yoffs[MAX_TEL_DX];
      double sat = (double)psr->obsn[ipos].sat;
      for (int ioff =0;ioff<psr->nTelDX;ioff++){
	yoffs[ioff]=0;
      }
      yoffs[k] = 1;
      
      if (sat < (double)psr->telDX_t[0]){
	// we are before the first jump
	// so our gradient is just the zeroth offset.
	afunc = yoffs[0]*psr->posPulsar[0];
;
      } else if(sat > (double)psr->telDX_t[(int)psr->nTelDX-1]){
	afunc = yoffs[(int)psr->nTelDX-1]*psr->posPulsar[0];
      } else{
	// find the pair we are between...
	for (int ioff =0;ioff<psr->nTelDX;ioff++){
	  if(sat >= psr->telDX_t[ioff] && sat < psr->telDX_t[ioff+1]){
	    double x1 = psr->telDX_t[ioff];
	    double x2 = psr->telDX_t[ioff+1];
	    double x = (sat-x1)/(x2-x1);
	    double y1=yoffs[ioff];
	    double y2=yoffs[ioff+1];
	    afunc = ((y2-y1)*x + y1)*psr->posPulsar[0];
	    break;
	  }
	}
      }
      //      printf("afunc = %g\n",afunc);
    }
  else if (i==param_tel_dy)
    {
      double yoffs[MAX_TEL_DY];
      double sat = (double)psr->obsn[ipos].sat;
      for (int ioff =0;ioff<psr->nTelDY;ioff++){
	yoffs[ioff]=0;
      }
      yoffs[k] = 1;
      
      if (sat < (double)psr->telDY_t[0]){
	// we are before the first jump
	// so our gradient is just the zeroth offset.
	afunc = yoffs[0]*psr->posPulsar[1];
;
      } else if(sat > (double)psr->telDY_t[(int)psr->nTelDY-1]){
	afunc = yoffs[(int)psr->nTelDY-1]*psr->posPulsar[1];
      } else{
	// find the pair we are between...
	for (int ioff =0;ioff<psr->nTelDY;ioff++){
	  if(sat >= psr->telDY_t[ioff] && sat < psr->telDY_t[ioff+1]){
	    double x1 = psr->telDY_t[ioff];
	    double x2 = psr->telDY_t[ioff+1];
	    double x = (sat-x1)/(x2-x1);
	    double y1=yoffs[ioff];
	    double y2=yoffs[ioff+1];
	    afunc = ((y2-y1)*x + y1)*psr->posPulsar[1];
	    break;
	  }
	}
      }
      //      printf("afunc = %g\n",afunc);
    }
  else if (i==param_tel_dz)
    {
      double yoffs[MAX_TEL_DZ];
      double sat = (double)psr->obsn[ipos].sat;
      for (int ioff =0;ioff<psr->nTelDZ;ioff++){
	yoffs[ioff]=0;
      }
      yoffs[k] = 1;
      
      if (sat < (double)psr->telDZ_t[0]){
	// we are before the first jump
	// so our gradient is just the zeroth offset.
	afunc = yoffs[0]*psr->posPulsar[2];
;
      } else if(sat > (double)psr->telDZ_t[(int)psr->nTelDZ-1]){
	afunc = yoffs[(int)psr->nTelDZ-1]*psr->posPulsar[2];
      } else{
	// find the pair we are between...
	for (int ioff =0;ioff<psr->nTelDZ;ioff++){
	  if(sat >= psr->telDZ_t[ioff] && sat < psr->telDZ_t[ioff+1]){
	    double x1 = psr->telDZ_t[ioff];
	    double x2 = psr->telDZ_t[ioff+1];
	    double x = (sat-x1)/(x2-x1);
	    double y1=yoffs[ioff];
	    double y2=yoffs[ioff+1];
	    afunc = ((y2-y1)*x + y1)*psr->posPulsar[2];
	    break;
	  }
	}
      }
      //      printf("afunc = %g\n",afunc);
    }
  else if (i==param_ifunc) /* Whitening procedure using interpolated function */
    {
      if (psr->param[param_ifunc].val[0]==1)
	{
	  double dt = (x + psr->param[param_pepoch].val[0]) - psr->ifuncT[k];
	  double tt = M_PI/(psr->ifuncT[1] - psr->ifuncT[0])*dt;
	  double t1;
	  double t2=0.0;
	  int l;
	  
	  t1 = sin(tt)/tt;
	  afunc = t1;
	}
      else if (psr->param[param_ifunc].val[0]==2) // Linear interpolation
	{
	  double yoffs[MAX_IFUNC];
	  double sat = (double)psr->obsn[ipos].sat;
	  for (int ioff =0;ioff<psr->ifuncN;ioff++){
	      yoffs[ioff]=0;
	  }
	  yoffs[k] = 1;

	  if (sat < (double)psr->ifuncT[0]){
	    // we are before the first jump
	    // so our gradient is just the zeroth offset.
	    afunc = yoffs[0];
	  } else if(sat > (double)psr->ifuncT[(int)psr->ifuncN-1]){
	    afunc = yoffs[(int)psr->ifuncN-1];
	  } else{
	    // find the pair we are between...
	    for (int ioff =0;ioff<psr->ifuncN;ioff++){
	      if(sat >= psr->ifuncT[ioff] && sat < psr->ifuncT[ioff+1]){
		double x1 = psr->ifuncT[ioff];
		double x2 = psr->ifuncT[ioff+1];
		double x = (sat-x1)/(x2-x1);
		double y1=yoffs[ioff];
		double y2=yoffs[ioff+1];
		afunc = (y2-y1)*x + y1;
		break;
	      }
	    }
	  }
	}
    }
    else if (i==param_quad_ifunc_p) /* Whitening procedure using interpolated function */
    {
      double yoffs[100];
      double sat = (double)psr->obsn[ipos].sat;
      for (int ioff =0;ioff<psr->quad_ifuncN_p;ioff++){
	yoffs[ioff]=0;
      }
      yoffs[k] = 1;
      
      if (sat < (double)psr->quad_ifuncT_p[0]){
	// we are before the first jump
	// so our gradient is just the zeroth offset.
	afunc = yoffs[0];
      } else if(sat > (double)psr->quad_ifuncT_p[(int)psr->quad_ifuncN_p-1]){
	afunc = yoffs[(int)psr->quad_ifuncN_p-1];
      } else{
	// find the pair we are between...
	for (int ioff =0;ioff<psr->quad_ifuncN_p;ioff++){
	  if(sat >= psr->quad_ifuncT_p[ioff] && sat < psr->quad_ifuncT_p[ioff+1]){
	    double x1 = psr->quad_ifuncT_p[ioff];
	    double x2 = psr->quad_ifuncT_p[ioff+1];
	    double x = (sat-x1)/(x2-x1);
	    double y1=yoffs[ioff];
	    double y2=yoffs[ioff+1];
	    afunc = (y2-y1)*x + y1;
	    break;
	  }
	}
      }
      afunc *= psr->quad_ifunc_geom_p;
    }
    else if (i==param_quad_ifunc_c) /* Whitening procedure using interpolated function */
    {
      double yoffs[100];
      double sat = (double)psr->obsn[ipos].sat;
      for (int ioff =0;ioff<psr->quad_ifuncN_c;ioff++){
	yoffs[ioff]=0;
      }
      yoffs[k] = 1;
      
      if (sat < (double)psr->quad_ifuncT_c[0]){
	// we are before the first jump
	// so our gradient is just the zeroth offset.
	afunc = yoffs[0];
      } else if(sat > (double)psr->quad_ifuncT_c[(int)psr->quad_ifuncN_c-1]){
	afunc = yoffs[(int)psr->quad_ifuncN_c-1];
      } else{
	// find the pair we are between...
	for (int ioff =0;ioff<psr->quad_ifuncN_c;ioff++){
	  if(sat >= psr->quad_ifuncT_c[ioff] && sat < psr->quad_ifuncT_c[ioff+1]){
	    double x1 = psr->quad_ifuncT_c[ioff];
	    double x2 = psr->quad_ifuncT_c[ioff+1];
	    double x = (sat-x1)/(x2-x1);
	    double y1=yoffs[ioff];
	    double y2=yoffs[ioff+1];
	    afunc = (y2-y1)*x + y1;
	    break;
	  }
	}
      }
      afunc *= psr->quad_ifunc_geom_c;
    }
  else if (i==param_quad_om)
    {
      double n1,n2,n3;
      double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
      double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;
      double cosTheta,omega_g;
      long double resp,resc,res_r,res_i;
      double theta_p,theta_g,phi_p,phi_g;
      double lambda_p,beta_p,lambda,beta;
      long double time;
      time    = (psr->obsn[ipos].bbat - psr->quadEpoch)*86400.0L;
      lambda_p = (double)psr->param[param_raj].val[0];
      beta_p   = (double)psr->param[param_decj].val[0];
      lambda   = psr->quadRA;
      beta     = psr->quadDEC;
      // Pulsar vector
      n1 = cosl(lambda_p)*cosl(beta_p);
      n2 = sinl(lambda_p)*cosl(beta_p);
      n3 = sinl(beta_p);
      cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
	sinl(beta)*sinl(beta_p);

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

      omega_g = (double)psr->param[param_quad_om].val[0]*((int)(k/4.0)+1);

      //            printf("In fitting with k = %d %d\n",k,k%4);
      if ((1.0L-cosTheta)==0)
	afunc=0;
      else
	{
	  if (k%4==0)      afunc = resp*sin(omega_g*time)/(2.0L*omega_g*(1.0L-cosTheta)); // aplus_re
	  else if (k%4==1) afunc = resp*cos(omega_g*time)/(2.0L*omega_g*(1.0L-cosTheta)); // aplus_im
	  else if (k%4==2) afunc = resc*sin(omega_g*time)/(2.0L*omega_g*(1.0L-cosTheta)); // across_re
	  else if (k%4==3) afunc = resc*cos(omega_g*time)/(2.0L*omega_g*(1.0L-cosTheta)); // across_im
	  //	  if (k%4==0)      afunc = resp*sin(omega_g*time)/(2.0L*omega_g*(1.0L-cosTheta)); // aplus_re
	  // else if (k%4==1) afunc = resp*cos(omega_g*time)/(2.0L*omega_g*(1.0L-cosTheta)); // aplus_im
	  //else if (k%4==2) afunc = resc*sin(omega_g*time)/(2.0L*omega_g*(1.0L-cosTheta)); // across_re
	  //else if (k%4==3) afunc = resc*cos(omega_g*time)/(2.0L*omega_g*(1.0L-cosTheta)); // across_im


	  //	  printf("Here with %d %d %g %g %g %g %g %g %g\n",k,k%4,afunc,sin(omega_g*time),cos(omega_g*time),(double)(omega_g*time),(double)resp,(double)resc,(double)cosTheta);
	  //	  else if (k%4==2) afunc = resp*(cos(omega_g*time)-1)/(2.0L*omega_g*(1.0L-cosTheta)); // aplus_im
	  // else if (k%4==3) afunc = resc*(cos(omega_g*time)-1)/(2.0L*omega_g*(1.0L-cosTheta)); // across_im
	}
      //      printf("Returning: %g %g %g\n",(double)afunc,(double)omega_g,(double)time);
    }
  else if (i==param_gwsingle)
    {
      double n1,n2,n3;
      double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
      double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;
      double cosTheta,omega_g;
      long double resp,resc,res_r,res_i;
      double theta_p,theta_g,phi_p,phi_g;
      double lambda_p,beta_p,lambda,beta;
      long double time;
      
      time    = (psr->obsn[ipos].bbat - psr->gwsrc_epoch)*86400.0L;
      lambda_p = (double)psr->param[param_raj].val[0];
      beta_p   = (double)psr->param[param_decj].val[0];
      lambda   = psr->gwsrc_ra;
      beta     = psr->gwsrc_dec;
      // Pulsar vector
      n1 = cosl(lambda_p)*cosl(beta_p);
      n2 = sinl(lambda_p)*cosl(beta_p);
      n3 = sinl(beta_p);
      cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
	sinl(beta)*sinl(beta_p);

      e11p = pow(sinl(lambda),2)-pow(cosl(lambda),2)*pow(sinl(beta),2);
      e21p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
      e31p = cosl(lambda)*sinl(beta)*cosl(beta);
      
      e12p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
      e22p = pow(cosl(lambda),2)-pow(sinl(lambda),2)*pow(sinl(beta),2);
      e32p = sinl(lambda)*sinl(beta)*cosl(beta);
      
      e13p = cosl(lambda)*sinl(beta)*cosl(beta);
      e23p = sinl(lambda)*sinl(beta)*cosl(beta);
      e33p = -powl(cosl(beta),2);

      omega_g = (double)psr->param[param_gwsingle].val[0];

      resp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
	      n2*(n1*e21p+n2*e22p+n3*e23p)+
	      n3*(n1*e31p+n2*e32p+n3*e33p));

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

      if ((1.0L-cosTheta)==0)
	afunc=0;
      else
	{
	  if (k==0)      afunc = resp*sin(omega_g*time)/(2.0L*omega_g*(1.0L-cosTheta)); // aplus_re
	  else if (k==1) afunc = resc*sin(omega_g*time)/(2.0L*omega_g*(1.0L-cosTheta)); // across_re
	  else if (k==2) afunc = resp*(cos(omega_g*time))/(2.0L*omega_g*(1.0L-cosTheta)); // aplus_im
	  else if (k==3) afunc = resc*(cos(omega_g*time))/(2.0L*omega_g*(1.0L-cosTheta)); // across_im

	  //	  else if (k==2) afunc = resp*(cos(omega_g*time)-1)/(2.0L*omega_g*(1.0L-cosTheta)); // aplus_im
	  //	  else if (k==3) afunc = resc*(cos(omega_g*time)-1)/(2.0L*omega_g*(1.0L-cosTheta)); // across_im
	}
    }
  else if (i==param_dmmodel)
    {
      double sat = (double)psr->obsn[ipos].sat;
      double yoffs[100];
      for (int ioff =0;ioff<psr->dmoffsNum;ioff++){
	      if (ioff==(k%psr->dmoffsNum)){ // we use modulo to get both the DM and CM signal.
		      yoffs[ioff]=1;
	      } else {
		      yoffs[ioff]=0;//-1.0/(double)(psr->dmoffsNum-1);
	      }
      }

      if (sat < (double)psr->dmoffsMJD[0]){
	      // we are before the first jump
	      // so our gradient is just the zeroth offset.
	      afunc = yoffs[0];
      } else if(sat > (double)psr->dmoffsMJD[(int)psr->dmoffsNum-1]){
	      afunc = yoffs[(int)psr->dmoffsNum-1];
      } else{
	      // find the pair we are between...
	      for (int ioff =0;ioff<psr->dmoffsNum;ioff++){
		      if(sat >= psr->dmoffsMJD[ioff] && sat < psr->dmoffsMJD[ioff+1]){
			      double x1 = psr->dmoffsMJD[ioff];
			      double x2 = psr->dmoffsMJD[ioff+1];
			      double x = (sat-x1)/(x2-x1);
			      double y1=yoffs[ioff];
			      double y2=yoffs[ioff+1];
			      afunc = (y2-y1)*x + y1;
			      break;
		      }
	      }
      }
      
/*
      double ti = (double)psr->obsn[ipos].sat;
      double tm1 = (double)psr->dmoffsMJD[k-1];
      double t0 = (double)psr->dmoffsMJD[k];
      double t1 = (double)psr->dmoffsMJD[k+1];
      if (ti >= t0 && ti < t1)
	{
	  afunc=1.0-(1.0/(t1-t0))*(ti-t0);
	}
      else if (ti >= tm1 && ti < t0)
	{
	  afunc=(1.0/(t0-tm1))*(ti-tm1);
	}
      else
	afunc = 0;
	*/
//      if(fabs(af2-afunc) > 1e-2)printf("XX %lf\t%lf\t%lf\n",sat,af2,afunc);
    }
  else if (i==param_dmassplanet)
    {
      afunc = dotproduct(psr->posPulsar,psr->obsn[ipos].planet_ssb[k]);
      //      printf("afunc = %.15g %g\n",(double)psr[0].obsn[ipos].sat,(double)afunc);
    }
  else if (strcmp(psr->binaryModel,"BT")==0) /* Must be below other parameters */
    afunc = BTmodel(psr,0,ipos,i);
  else if (strcmp(psr->binaryModel,"BTJ")==0) 
    afunc = BTJmodel(psr,0,ipos,i,k);
  else if (strcmp(psr->binaryModel,"ELL1")==0) 
    afunc = ELL1model(psr,0,ipos,i);	  
  else if (strcmp(psr->binaryModel,"DD")==0) 
    afunc = DDmodel(psr,0,ipos,i);	  
  else if (strcmp(psr->binaryModel,"DDK")==0) 
    afunc = DDKmodel(psr,0,ipos,i);	  
  else if (strcmp(psr->binaryModel,"DDS")==0) 
    afunc = DDSmodel(psr,0,ipos,i);	  
  else if (strcmp(psr->binaryModel,"DDGR")==0) 
    afunc = DDGRmodel(psr,0,ipos,i);	  
  else if (strcmp(psr->binaryModel,"MSS")==0) 
    afunc = MSSmodel(psr,0,ipos,i);
  else if (strcmp(psr->binaryModel,"T2")==0) 
    afunc = T2model(psr,0,ipos,i,k);	  
  

  return afunc;
}


void updateParameters(pulsar *psr,int p,double *val,double *error)
{
  int i,j,k;

  if (debugFlag==1) printf("Updating parameters\n");
  psr[p].offset = val[0];
  psr[p].offset_e = error[0];
  j=1;
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[p].param[i].aSize;k++)
	{
	  if (psr[p].param[i].fitFlag[k]==1 && (i!=param_start && i!=param_finish))
	    {
	      if (i==param_f) 
		{
		  if (k==0)
		    {
		      psr[p].param[param_f].val[k] *= (1.0-val[j]/psr[p].param[param_f].val[0]);	
		      psr[p].param[param_f].err[k]  = error[j];
		    }
		  else 
		    {
		      longdouble scale;
		      scale=1.0L;
		      if (k==2)      scale=1.0e9L;
		      else if (k>2 && k<10)  scale=1.0e18L;
		      else if (k>9) scale=1.0e23L;
		      
		      psr[p].param[param_f].val[k] = psr[p].param[param_f].val[k] - 
			(psr[p].param[param_f].val[0]*(val[j]/pow(24.0*3600.0,k+1))/scale);
		      psr[p].param[param_f].err[k] = error[j]/(pow(24.0*3600.0,k+1))/scale*
			psr[p].param[param_f].val[0];
		    }
		}
	      else if (i==param_dm || i==param_px || i==param_fddc || i==param_fddi || i==param_dmassplanet || i==param_dmx)
		{
		  psr[p].param[i].val[k] += val[j];
		  psr[p].param[i].err[k]  = error[j];
		  // This is slow - should be a better approach
		  if (i==param_dm){
		    psr[p].dmOffset+=val[j];
		  }
		}
	      else if (i==param_dshk)
		{
		  psr[p].param[i].val[k] += val[j];
		  psr[p].param[i].err[k]  = error[j];
		}
	      else if (i==param_pmrv)
		{
		  psr[p].param[i].val[k] += 10.0*val[j]*360.0*60.0*60.0/(2.0*M_PI);
		  psr[p].param[i].err[k]  = 10.0*error[j]*360.0*60.0*60.0/(2.0*M_PI);
		}
	      else if (i==param_glph) /* Glitch phase */
		{
		  psr[p].param[i].val[k] -= val[j];     
		  psr[p].param[i].err[k]  = error[j];   
		}
	      else if (i==param_glf0d) /* Glitch */
		{
		  psr[p].param[i].val[k] -= val[j]; 
		  psr[p].param[i].err[k]  = error[j]; 
		}
	      else if (i==param_gltd) /* Glitch time delay */
		{
		  psr[p].param[i].val[k] -= val[j]; 
		  psr[p].param[i].err[k]  = error[j]; 
		}
	      else if (i==param_glf0) /* Glitch permanent pulse frequency increment */
		{
		  psr[p].param[i].val[k] -= val[j]; 
		  psr[p].param[i].err[k]  = error[j];                        
		}
	      else if (i==param_glf1) /* Glitch permanent pulse frequency deriv. increment */
		{
		  psr[p].param[i].val[k] -= val[j]; //*psr[p].param[param_f].val[0]; 
		  psr[p].param[i].err[k]  = error[j];                        
		}
	      else if (i==param_glf2) /* Glitch permanent pulse frequency second deriv. increment */
		{
		  psr[p].param[i].val[k] -= val[j]*1.0e-27; //*psr[p].param[param_f].val[0]; 
		  psr[p].param[i].err[k]  = error[j]*1.0e-27;                        
		}
 	      else if (i==param_telx || i==param_tely || i==param_telz)
		{
		  psr[p].param[i].val[k] -= val[j];
		  psr[p].param[i].err[k] = error[j];
		}	      
	      else if (i==param_raj)
		{
		  char retstr[100];
		  psr[p].param[param_raj].val[k] += val[j];
		  psr[p].param[param_raj].err[k] = error[j];
		  
		  /* Must obtain this in hms form */
		  turn_hms(psr[p].param[param_raj].val[k]/(2.0*M_PI), retstr);
		  strcpy(psr[p].rajStrPost,retstr);
		}
	      else if (i==param_decj)
		{
		  char retstr[100];
		  psr[p].param[param_decj].val[k] += val[j];
		  psr[p].param[param_decj].err[k] = error[j];
		  /* Must obtain this in dms form */
		  turn_dms(psr[p].param[param_decj].val[k]/(2.0*M_PI), retstr);
		  strcpy(psr[p].decjStrPost,retstr);
		}
	      else if (i==param_pmra) /* Return in radian/sec */
		{
		  psr[p].param[param_pmra].val[k] += val[j]*180.0/M_PI*60.0*60.0*
		    1000.0*SECDAY*365.25/24.0/3600.0*cos(psr[p].param[param_decj].val[0]);
		  psr[p].param[param_pmra].err[k] = error[j]*180.0/M_PI*60.0*60.0*
		    1000.0*SECDAY*365.25/24.0/3600.0*cos(psr[p].param[param_decj].val[0]);
		}
	      else if (i==param_pmdec) /* Return in radian/sec */
		{
		  psr[p].param[param_pmdec].val[k] += val[j]*180.0/M_PI*60.0*60.0*1000.0*
		    SECDAY*365.25/24.0/3600.0;
		  psr[p].param[param_pmdec].err[k] = error[j]*180.0/M_PI*60.0*60.0*1000.0*
		    SECDAY*365.25/24.0/3600.0;
		}	       
	      else if (i==param_wave_om) /* Whitening procedure using sinusoids */
		{
		  int k;
		  for (k=0;k<psr->nWhite;k++)
		    {
		      psr[p].wave_cos[k]  -= val[j]; 
		      psr[p].wave_cos_err[k] = error[j]; j++;
		      psr[p].wave_sine[k] -= val[j]; 
		      psr[p].wave_sine_err[k] = error[j]; j++;	      
		    }
		  if (psr->waveScale==2) // Ignore the non-frequency derivative terms
		    {
		      for (k=0;k<psr->nWhite;k++)
			{
			  printf("Ignoring cos %g %g\n",val[j],error[j]); j++;
			  printf("Ignoring sin %g %g\n",val[j],error[j]); j++;
			}
		      //		      j+=psr->nWhite*2; 
		    }
		  j--;
		}		  
	      else if (i==param_ifunc) 
		{
		  int k;
		  for (k=0;k<psr->ifuncN;k++)
		    {
		      printf("val ret = %g\n",val[j]);
		      psr[p].ifuncV[k] -= val[j];
		      psr[p].ifuncE[k] = error[j];
		      j++;
		    }
		}		  
	      else if (i==param_tel_dx) 
		{
		  int k;
		  for (k=0;k<psr->nTelDX;k++)
		    {
		      printf("val ret = %g\n",val[j]);
		      psr[p].telDX_v[k] -= val[j];
		      psr[p].telDX_e[k] = error[j];
		      j++;
		    }
		}		  
	      else if (i==param_tel_dy) 
		{
		  int k;
		  for (k=0;k<psr->nTelDY;k++)
		    {
		      printf("val ret = %g\n",val[j]);
		      psr[p].telDY_v[k] -= val[j];
		      psr[p].telDY_e[k] = error[j];
		      j++;
		    }
		}		  
	      else if (i==param_tel_dz) 
		{
		  int k;
		  for (k=0;k<psr->nTelDZ;k++)
		    {
		      printf("val ret = %g\n",val[j]);
		      psr[p].telDZ_v[k] -= val[j];
		      psr[p].telDZ_e[k] = error[j];
		      j++;
		    }
		}		  
	      else if (i==param_quad_ifunc_p) 
		{
		  int k;
		  for (k=0;k<psr->quad_ifuncN_p;k++)
		    {
		      printf("val ret = %g\n",val[j]);
		      psr[p].quad_ifuncV_p[k] -= val[j];
		      psr[p].quad_ifuncE_p[k] = error[j];
		      j++;
		    }
		}		  
	      else if (i==param_quad_ifunc_c) 
		{
		  int k;
		  for (k=0;k<psr->quad_ifuncN_c;k++)
		    {
		      printf("val ret = %g\n",val[j]);
		      psr[p].quad_ifuncV_c[k] -= val[j];
		      psr[p].quad_ifuncE_c[k] = error[j];
		      j++;
		    }
		}		  
	      else if (i==param_gwsingle)
		{
		  printf("%d GW SINGLE -- in here %g %g %g %g\n",j,val[j],error[j],val[j+1],error[j+1]);
		  psr[p].gwsrc_aplus_r -= val[j];
		  j++;
		  psr[p].gwsrc_across_r -= val[j];
		  j++;
		  psr[p].gwsrc_aplus_i -= val[j];
		  j++;
		  psr[p].gwsrc_across_i -= val[j];
		  printf("Now have: %g %g\n",psr[p].gwsrc_aplus_r,psr[p].gwsrc_across_r);
		}
	      else if (i==param_dmmodel)
		{
		  for (k=0;k<psr[p].dmoffsNum;k++)
		    {
		      //		      printf("Updating %g by %g\n",psr[p].dmvalsDM[k],val[j]);
		      //
		      psr[p].dmoffsDM[k] += val[j];
		      psr[p].dmoffsDMe[k] = error[j];
		      j++;

		      psr[p].dmoffsOffset[k] = val[j];
		      psr[p].dmoffsError[k] =  error[j];
		      j++;
		    }
		  j--;
		}
	      else if (i==param_start)
		{
		}
	      else if (i==param_finish)
		{
		}
	      else if (strcmp(psr[p].binaryModel,"BT")==0)
		updateBT(&psr[p],val[j],error[j],i);
	      else if (strcmp(psr[p].binaryModel,"BTJ")==0)
		updateBTJ(&psr[p],val[j],error[j],i,k);
	      else if (strcmp(psr[p].binaryModel,"ELL1")==0)
		updateELL1(&psr[p],val[j],error[j],i);
	      else if (strcmp(psr[p].binaryModel,"DD")==0)
		updateDD(&psr[p],val[j],error[j],i);
	      else if (strcmp(psr[p].binaryModel,"DDK")==0)
		updateDDK(&psr[p],val[j],error[j],i);
	      else if (strcmp(psr[p].binaryModel,"DDS")==0)
		updateDDS(&psr[p],val[j],error[j],i);
	      else if (strcmp(psr[p].binaryModel,"DDGR")==0)
		updateDDGR(&psr[p],val[j],error[j],i);
	      else if (strcmp(psr[p].binaryModel,"MSS")==0)
		updateMSS(&psr[p],val[j],error[j],i);
	      else if (strcmp(psr[p].binaryModel,"T2")==0)
		updateT2(&psr[p],val[j],error[j],i,k);
	      j++; /* Increment position in fit list */
	    }
	}
    }
  if (strcmp(psr[p].binaryModel,"DDGR")==0) 
    DDGRmodel(psr,0,0,-2);  /* Update GR parameters */	  
  
  if (debugFlag==1) printf("Updating jumps; nJumps = %d\n",psr[p].nJumps);
  /* Now check jumps */
  for (i=1;i<=psr[p].nJumps;i++)
    {
      if (debugFlag==1) printf("%d fitJump = %d\n",i,psr[p].fitJump[i]);
      if (psr[p].fitJump[i]==1)
	{
	  if (debugFlag==1) printf("%d Jump changes\n",i);
	  if (debugFlag==1) printf("value = %g\n",(double)val[j]);
	  if (debugFlag==1) printf("error = %g\n",(double)error[j]);
	  psr[p].jumpVal[i] += -val[j];
	  psr[p].jumpValErr[i] = error[j];
	  j++;
	}
      /*	      printf("Have jumps %g %g\n",(double)val[j],error[j][j]); */
    }
  if (debugFlag==1) printf("Complete updating parameters\n");
}


void formCholeskyMatrix(double *c,double *resx,double *resy,double *rese,int np,int nc,double **uinv)
{
  double **m,**u,sum;
  double *cholp;
  int i,j,k,ix,iy;
  double t0,cint,t;
  int t1,t2;
  int debug=0;
  clock_t clk;

  clk = clock();

  printf("Tcheck: forming Cholesky matrix ... allocating memory  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
//  printf("Getting the covariance matrix in doFit\n");
  m = (double **)malloc(sizeof(double *)*(np+1));
  u= (double **)malloc(sizeof(double *)*(np+1));
  cholp  = (double *)malloc(sizeof(double)*(np+1));  // Was ndays

  
  for (i=0;i<np+1;i++)
    {
      m[i] = (double *)malloc(sizeof(double)*(np+1));
      u[i] = (double *)malloc(sizeof(double)*(np+1));
    }
  printf("Tcheck: forming Cholesky matrix ... complete allocating memory  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
  printf("Tcheck: forming Cholesky matrix ... determing m[ix][iy] = fabs(resx[ix]-resx[iy])  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
  for (ix=0;ix<np;ix++)
    {
      for (iy=0;iy<np;iy++)
	m[ix][iy] = fabs(resx[ix]-resx[iy]);
    }
  printf("Tcheck: forming Cholesky matrix ... complete determing m[ix][iy] = fabs(resx[ix]-resx[iy])  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
  if (debug==1)
    {
      printf("First m = \n");
      for (i=0;i<5;i++)
	{ 
	  for (j=0;j<5;j++) printf("%10g ",m[i][j]); 
	  printf("\n");
	}
    }
  
  // Insert the covariance which depends only on the time difference.
  // Linearly interpolate between elements on the covariance function because
  // valid covariance matrix must have decreasing off diagonal elements.
  //  printf("Inserting into the covariance matrix\n");
  printf("Tcheck: forming Cholesky matrix ... determing covariance based on time difference  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
  for (ix=0;ix<np;ix++)
    {
      for (iy=0;iy<np;iy++)
	{
	  t0 = m[ix][iy];
	  t1 = (int)floor(t0);
	  t2 = t1+1;
	  t  = t0-t1;
	  cint = c[t1]*(1-t)+c[t2]*t; // Linear interpolation
	  m[ix][iy] = cint;
	}
    }
  printf("Tcheck: forming Cholesky matrix ... complete determing covariance based on time difference  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);

  // add the values for the constraints
  // Constraints are not covariant with anything so it's all zero!
  for (i=np-nc; i < np; i++){
	  for (j=0; j < np; j++){
		  m[i][j]=0;
		  m[j][i]=0;
	  }
  }


  //  printf("Multiplying by errors\n");
  printf("Tcheck: forming Cholesky matrix ... multiplying by errors  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
  for (ix=0;ix<np;ix++)
    m[ix][ix]+=rese[ix]*rese[ix];
  printf("Tcheck: forming Cholesky matrix ... complete multiplying by errors  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
  if (debug==1)
    {
      printf("m = \n\n");
      for (i=np-5;i<np;i++)
	{ 
	  for (j=np-5;j<np;j++) printf("%10g ",m[i][j]); 
	  printf("\n");
	}
    }

  
  // Do the Cholesky
  printf("Tcheck: forming Cholesky matrix ... do Cholesky decomposition  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
  TKcholDecomposition(m,np,cholp);
  printf("Tcheck: forming Cholesky matrix ... complete do Cholesky decomposition  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
  // Now calculate uinv
  printf("Tcheck: forming Cholesky matrix ... calculate uinv  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
    for (i=0;i<np;i++)
    {
      m[i][i] = 1.0/cholp[i];
      uinv[i][i] = m[i][i];
      for (j=0;j<i;j++)
      	uinv[i][j] = 0.0;
      for (j=i+1;j<np;j++)
	{
	  sum=0.0;
	  for (k=i;k<j;k++) sum-=m[j][k]*m[k][i];
	  m[j][i]=sum/cholp[j];
	  uinv[i][j] = m[j][i];
	}
	} 
    // Note we can use something like the following - see p98 in nrec, but it doesn't speed things up much 
    /*
      for (i=0;i<np;i++)
      {
      uinv[i][i] = 1.0/cholp[i];
      for (j=i+1;j<np;j++)
      {
      sum=0.0;
      for (k=i;k<j;k++) sum-=m[j][k]*m[k][i];
      uinv[i][j] = sum/cholp[j];
      }
      }*/
 
  printf("Tcheck: forming Cholesky matrix ... complete calculate uinv  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
  if (debug==1)
    {
      printf("uinv = \n\n");
      for (i=0;i<5;i++)
	{ 
	  for (j=0;j<5;j++) printf("%10g ",uinv[i][j]); 
	  printf("\n");
	}
    }

  //  printf("Completed inverting the matrix\n");

  // Should free memory not required
  // (note: not freeing uinv)
  printf("Tcheck: forming Cholesky matrix ... free memory  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
  for (i=0;i<np+1;i++)
    {
      free(m[i]);
      free(u[i]);
    }
  free(m);
  free(u);
  free(cholp);
  printf("Tcheck: forming Cholesky matrix ... complete free memory  (%.2f)\n",(clock()-clk)/(float)CLOCKS_PER_SEC);
}


/*
 *
 * Constraint functions are in constraints.C
 *
 */
double getConstraintDeriv(pulsar *psr,int iconstraint,int i,int k){
	int order=0;
		switch(psr->constraints[iconstraint]){
		case constraint_dmmodel_mean:
			return consFunc_dmmodel_mean(psr,i,k);

		// Notice that these case statements fall through, so order is defined
		// based on which constraint name is used.
		case constraint_dmmodel_cw_2:
			order++;
		case constraint_dmmodel_cw_1:
			order++;
		case constraint_dmmodel_cw_0:
			return consFunc_dmmodel_cw(psr,i,k,order);
		case constraint_ifunc_2:
			order++;
		case constraint_ifunc_1:
			order++;
		case constraint_ifunc_0:
			return consFunc_ifunc(psr,i,k,order);
		case constraint_tel_dx_2:
		  order++;
		case constraint_tel_dx_1:
		  order++;
		case constraint_tel_dx_0:
		  return consFunc_tel_dx(psr,i,k,order);
		case constraint_tel_dy_2:
		  order++;
		case constraint_tel_dy_1:
		  order++;
		case constraint_tel_dy_0:
		  return consFunc_tel_dy(psr,i,k,order);
		case constraint_tel_dz_2:
		  order++;
		case constraint_tel_dz_1:
		  order++;
		case constraint_tel_dz_0:
		  return consFunc_tel_dz(psr,i,k,order);
		case constraint_quad_ifunc_p_2:
			order++;
		case constraint_quad_ifunc_p_1:
			order++;
		case constraint_quad_ifunc_p_0:
			return consFunc_quad_ifunc_p(psr,i,k,order);
		case constraint_quad_ifunc_c_2:
			order++;
		case constraint_quad_ifunc_c_1:
			order++;
		case constraint_quad_ifunc_c_0:
			return consFunc_quad_ifunc_c(psr,i,k,order);

		default:
			return 0;
			}
}


