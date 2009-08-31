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

int getNparams(pulsar psr);

/* Main routines for fitting in TEMPO2               */
void doFit(pulsar *psr,int npsr,int writeModel) 
{
  int p,npol,i,j,okay,ip[MAX_OBSN],k;
  double *x,*y,*sig,*val,chisq;
  double *error;
  double tol = 1.0e-40;  /* Tolerence for singular value decomposition routine */
  double newStart=-1.0,newFinish=-1.0;

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
      sprintf(str,"%s/plugins/%s_fitFunc_%s_plug.so",getenv(TEMPO2_ENVIRON),
	      psr[0].fitFunc,tempo2MachineType);
      module = dlopen(str, RTLD_NOW); 
      if(!module)  {
	fprintf(stderr, "[error]: dlopen() unable to open plugin: %s.\n",str);
	exit(1);
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
      //      strcpy(psr[p].rajStrPost,psr[p].rajStrPre);
      //      strcpy(psr[p].decjStrPost,psr[p].decjStrPre);

      strcpy(psr[p].rajStrPre,psr[p].rajStrPost);
      strcpy(psr[p].decjStrPre,psr[p].decjStrPost);
      /* How many parameters are we fitting for */
      npol = getNparams(psr[p]);
      x     = (double *)malloc(psr[p].nobs*sizeof(double));
      y     = (double *)malloc(psr[p].nobs*sizeof(double));
      sig   = (double *)malloc(psr[p].nobs*sizeof(double));
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
	  psr[p].fitChisq = chisq; 
	  psr[p].fitNfree = psr[p].nFit-npol;
	  //	  printf("Chisq = %g, reduced chisq = %g\n",(double)psr[p].fitChisq,(double)(psr[p].fitChisq/psr[p].fitNfree));
	  /* Now update the parameters */
	  if (debugFlag==1) printf("Updating the parameters\n");
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
void doFitDCM(pulsar *psr,char *dcmFile,int npsr,int writeModel) 
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

  printf("WARNING: Switching weighting off for the fit\n");
  printf("WARNING: THE .TIM FILE MUST BE SORTED - not checked for\n");
  for (p=0;p<npsr;p++)  /* Loop over all the pulsars */
    {
      psr[p].fitMode = 0;
      for (i=0;i<psr[p].nobs;i++)
	meanRes+=(long double)psr[p].obsn[i].residual;
      meanRes/=(long double)psr[p].nobs;
      for (i=0;i<psr[p].nobs;i++)
	psr[p].obsn[i].residual-=meanRes;
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
      printf("Opening >%s<\n",fname);
      if (!(fin = fopen(fname,"r")))
	{
	  printf("Unable to open inverse cholesky matrix: %s\n",fname);
	  exit(1);
	}
      
      uinv = (double **)malloc(sizeof(double *)*psr[p].nobs);
      for (i=0;i<psr[p].nobs;i++)
	uinv[i] = (double *)malloc(sizeof(double)*psr[p].nobs);
      
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
      sprintf(fname,"whitedata_%d.dat",p+1);
      fout = fopen(fname,"w");
      for (i=0;i<psr[p].nobs;i++)
	{
	  sum=0.0;
	  for (j=0;j<psr[p].nobs;j++)
	    sum+=uinv[j][i]*psr[p].obsn[j].residual;
	  printf("Have: %d %g %g %g\n",i,sum,uinv[0][i],(double)psr[p].obsn[i].residual);
	  whiteres[i] = sum;
	  //            fprintf(fout,"%g %g %g\n",(double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]),(double)psr[0].obsn[i].residual,(double)psr[0].obsn[i].toaErr);
	  fprintf(fout,"%g %g %g %g\n",(double)(psr[p].obsn[i].sat-psr[p].param[param_pepoch].val[0]),
		  (double)psr[p].obsn[i].residual,whiteres[i],(double)psr[p].obsn[i].toaErr);
	}
      fclose(fout);
      
      //  for (p=0;p<npsr;p++)  /* Loop over all the pulsars */
      //    {
      //      strcpy(psr[p].rajStrPost,psr[p].rajStrPre);
      //      strcpy(psr[p].decjStrPost,psr[p].decjStrPre);
      
      strcpy(psr[p].rajStrPre,psr[p].rajStrPost);
      strcpy(psr[p].decjStrPre,psr[p].decjStrPost);
      /* How many parameters are we fitting for */
      npol = getNparams(psr[p]);
      x     = (double *)malloc(psr[p].nobs*sizeof(double));
      y     = (double *)malloc(psr[p].nobs*sizeof(double));
      sig   = (double *)malloc(psr[p].nobs*sizeof(double));
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
      

      psr[p].nFit = count;
      psr[p].param[param_start].val[0] = newStart-0.001; 
      psr[p].param[param_finish].val[0] = newFinish+0.001;
      psr[p].param[param_start].paramSet[0] = 1;
      psr[p].param[param_finish].paramSet[0] = 1; 
      /* Do the fit */
      if (npol!=0) /* Are we actually  doing any fitting? */ 
	{ 
	  if (debugFlag==1) printf("Doing the fit\n");
	  TKleastSquares_svd_psr_dcm(x,y,sig,psr[p].nFit,val,error,npol,psr[p].covar,&chisq,
				 FITfuncs,psr[p].fitMode,&psr[p],tol,ip,uinv);
	  //	  svdfit(x,y,sig,psr[p].nFit,val,npol,u,v,w,&chisq,FITfuncs,&psr[p],tol,ip);
	  if (debugFlag==1) printf("Complete fit: chisq = %f\n",(double)chisq);
	  printf("chisq = %g\n",chisq);
	  psr[p].fitChisq = chisq; 
	  psr[p].fitNfree = psr[p].nFit-npol;
	  
	  /* Now update the parameters */
	  if (debugFlag==1) printf("Updating the parameters\n");
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
    



  for (i=0;i<psr[p].nobs;i++)
    free(uinv[i]);
  free(uinv);
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
	    if (i!=param_start && i!=param_finish)
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
    npol+=psr.nWhite*2-1;
  
  return npol;
}

void FITfuncs(double x,double afunc[],int ma,pulsar *psr,int ipos)
{
  int i,n=0,k,j;
  
  afunc[n++] = 1;  /* Always fit for an arbitrary offset */
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
		      for (j=0;j<psr->nWhite*2;j++)
			afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
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
	  if (psr->obsn[ipos].jump==i)
	    afunc[n++] = 1.0;
	  else
	    afunc[n++] = 0.0;
	}
    } 
  if (n!=ma) { 
    printf("Problem in fitting routine n = %d, ma = %d\n",n,ma);
  }
}


double getParamDeriv(pulsar *psr,int ipos,double x,int i,int k)
{
  double afunc=0.0;
  if (i==param_f)         /* Rotational frequency */
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
	afunc = psr->param[param_gltd].val[0]*SECDAY*1.0e-9*(1.0-expf)*1.0e9/psr->param[param_f].val[0];
      else
	afunc = 0.0;
    }
  else if (i==param_gltd)
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
	afunc = psr->param[param_glf0d].val[k]*
	  (1.0-(1.0+dt1/1.0e9/(psr->param[param_gltd].val[0]*SECDAY*1.0e-9))*expf)/psr->param[param_f].val[0]*SECDAY;
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
	afunc = 1.0/6.0*pow((psr->obsn[ipos].bbat-psr->param[param_glep].val[k])*86400.0,3)/psr->param[param_f].val[0];
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
	afunc = 1.0/(DM_CONST*pow(psr->obsn[ipos].freqSSB/1.0e6,2));
      else
	{
	  yrs = (psr->obsn[ipos].sat - psr->param[param_pepoch].val[0])/365.25;
	  afunc = 1.0/(DM_CONST*pow(psr->obsn[ipos].freqSSB/1.0e6,2))*pow(yrs,k);
	}
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
      if (k%2==0) afunc = cos(om*(floor(k/2.0)+1)*x); 
      else        afunc = sin(om*(floor(k/2.0)+1)*x); 
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
	      else if (i==param_dm || i==param_px || i==param_fddc || i==param_fddi || i==param_dmassplanet)
		{
		  psr[p].param[i].val[k] += val[j];
		  psr[p].param[i].err[k]  = error[j];
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
		  psr[p].param[i].val[k] -= val[j];  /* Why - sign ??? */
		  psr[p].param[i].err[k]  = error[j];                        /* ERROR IS WRONG */
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
		  psr[p].param[i].val[k] -= val[j]; //*psr[p].param[param_f].val[0];
		  psr[p].param[i].err[k]  = error[j];                        
		}
	      else if (i==param_glf1) /* Glitch permanent pulse frequency deriv. increment */
		{
		  psr[p].param[i].val[k] -= val[j]; //*psr[p].param[param_f].val[0]; 
		  psr[p].param[i].err[k]  = error[j];                        
		}
	      else if (i==param_glf2) /* Glitch permanent pulse frequency second deriv. increment */
		{
		  psr[p].param[i].val[k] -= val[j]; //*psr[p].param[param_f].val[0]; 
		  psr[p].param[i].err[k]  = error[j];                        
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
