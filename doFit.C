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
#include <algorithm>
#include <string.h>
#include <dlfcn.h>
#include "tempo2.h"
#include "TKfit.h"
#include "constraints.h"
#include "TKsvd.h"
#include "T2accel.h"
void globalFITfuncs(double x,double afunc[],int ma,pulsar *psr,int ipos,int ipsr);

void updateGlobalParameters(pulsar* psr,int npsr, double* val,double* error);

int getNparams(pulsar *psr,int offset);
int getNglobal(pulsar *psr,int npsr);
double getConstraintDeriv(pulsar *psr,int ipos,int i,int k);


///////////////////Functions for TempoNest maximum likelihood Fitting////////////////

#ifdef HAVE_LAPACK
#ifdef HAVE_BLAS
void getTempoNestMaxLike(pulsar *pulse, int npsr);

void dgesvd_ftoc(double *in, double **out, int rows, int cols);
double* dgesvd_ctof(double **in, int rows, int cols);
void dgesvd(double **A, int m, int n, double *S, double **U, double **VT);
void dgemv(double **A, double *vecin,double *vecout,int rowa, int cola, char AT);
double *dgemv_ctof(double **in, int rows, int cols);
void dgemv_ftoc(double *in, double **out, int rows, int cols);
void dgemm(double **A, double **B,double **C,int rowa, int cola, int rowb, int colb, char AT, char BT);
double *dgemm_ctof(double **in, int rows, int cols);
void dgemm_ftoc(double *in, double **out, int rows, int cols);
void dpotri(double **A, int msize);
double *dpotri_ctof(double **in, int rows, int cols);
void dpotri_ftoc(double *in, double **out, int rows, int cols);
void dpotrf(double **A, int msize, double &det);
double *dpotrf_ctof(double **in, int rows, int cols);
void dpotrf_ftoc(double *in, double **out, int rows, int cols);

extern "C" void dpotrf_(char *UPLO, int *msize, double *a, int *lda, int *info);
extern "C" void dpotri_(char *UPLO, int *msize, double *a, int *lda, int *info);
extern "C" void dgesvd_(char *jobu, char *jobvt, int *m, int *n,
		double *a, int *lda, double *s, double *u,
		int *ldu, double *vt, int *ldvt, double *work,
		int *lwork, int *info);
extern "C" void dgemv_(char *jobu, int *m, int *n,
		double *alpha, double *a, int *lda,
		double *x, int *incx, double *beta, double *y, int *incy);
extern "C" void dgemm_(char *jobu, char *jobvt, int *m, int *n,
		int *k, double *alpha, double *a, int *lda,
		double *b, int *ldb, double *beta, double *c,
		int *ldc);
#endif
#endif
/////////////////End Functions for TempoNest maximum likelihood Fitting/////////////////////

// Backwards compatibility
void doFit(pulsar *psr,int npsr,int writeModel) 
{
	logmsg("Deprecated call to doFit() -> better to use doFitAll()");
	doFitAll(psr,npsr,NULL);
}


// Backwards compatibility
void doFitDCM(pulsar *psr,char *dcmFile,char *covarFuncFile,int npsr,int writeModel) {
	logmsg("Deprecated call to doFitDCM() -> better to use doFitAll()");
	if (strcmp(dcmFile,"NULL")){
		doFitAll(psr,npsr,dcmFile);
	} else {
		doFitAll(psr,npsr,covarFuncFile);
	}
}

/**
 * Master fitting routine with or without cholesky, global or not.
 */
void doFitAll(pulsar *psr,int npsr, char *covarFuncFile) {
	int i,j,k;
	double whiteres[MAX_OBSN],sum;
	FILE *fin,*fout;
	char fname[100],temp[100];
	int p,okay;
	int nobs_noconstrain=0;
	double *x,*y,*sig,*val,chisq;
	double *error;
	double tol = 1.0e-27;  /* Tolerence for singular value decomposition routine */
	double newStart=-1.0,newFinish=-1.0;
	long double meanRes=0.0;
	int count;
	int nobs_and_constraints;
	int offsetNp = 0;
	bool DO_GLOBAL_FIT=false;

	double **xx = (double**)malloc(sizeof(double*)*npsr);
	double **yy = (double**)malloc(sizeof(double*)*npsr);
	double ***uinvs = (double***)malloc(sizeof(double**)*npsr);
	int npol;
	int nglobal;
	int **ip=(int**)malloc(sizeof(int*)*npsr);


	//  printf("WARNING: Switching weighting off for the fit\n");
	//  printf("WARNING: THE .TIM FILE MUST BE SORTED - not checked for\n");


	//If TNRed or TNDM or ECORR parameters defined, call TempoNest maxlike function and return
	if((psr[0].TNRedAmp != 0 && psr[0].TNRedGam != 0) || (psr[0].TNDMAmp != 0 && psr[0].TNDMGam != 0 || (psr[0].nTNECORR != 0))){


#ifdef HAVE_LAPACK
#ifdef HAVE_BLAS
		printf("\n\nIncluding red noise and/or ECORR parameters, calling new dofit\n\n");
		getTempoNestMaxLike(psr, npsr);
		return;
#else
		printf("LAPACK required for use of temponest noise parameters\n");
#endif
#endif
	}




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
		entry(psr,npsr,0);
		printf("Returning\n");

		return;
	}

	nglobal=getNglobal(psr,npsr);
	offsetNp += nglobal;
	if (offsetNp > MAX_FIT)
	{
		printf("ERROR1: number of fitted parameters > MAX_FIT (%d)\n",MAX_FIT);
		exit(1);
	}
	if(nglobal > 0 || forceGlobalFit){
		DO_GLOBAL_FIT=true;
		logmsg("GLOBAL fit enabled. Number of global fit parameters = %d",nglobal);
	}

	for (p=0;p<npsr;p++)  /* Loop over all the pulsars */
	{
		// if we are doing a Cholesky fit, we need to sort the data or it will fail!
		if (covarFuncFile!=NULL && strcmp(covarFuncFile,"NULL"))sortToAs(psr+p);
		if (psr[p].auto_constraints)autoConstraints(psr,p,npsr);
		nobs_and_constraints = psr[p].nobs + psr[p].nconstraints;
		nobs_noconstrain += psr[p].nobs;
		ip[p]=(int*)malloc(sizeof(int)*nobs_and_constraints);
		if (ip[p]==NULL) {printf("Unable to allocate memory in doFit.C (ip[p])\n"); exit(1);}
		logtchk("Processing pulsar %d",p);

		/*
		 * Check for "broken" combinations of fit parameters
		 */
		if (psr[p].param[param_dmmodel].paramSet[0] &&
				psr[p].param[param_dm].fitFlag[0]!= 0 &&
				psr[p].param[param_dmmodel].linkTo[0] != param_dm){
			// if we are using DMMODEL, and we want to fit for DM then
			// DMMODEL must be linked to DM... otherwise badness will occur!
			printf("WARNING: DM cannot be fit with DMMODEL\n         unless you set 'DMMODEL DM 1' in par file.\n");
			psr[p].param[param_dm].fitFlag[0]=0;
		}
		strcpy(psr[p].rajStrPre,psr[p].rajStrPost);
		strcpy(psr[p].decjStrPre,psr[p].decjStrPost);
		/* How many parameters are we fitting for */
		logtchk("Determining which parameters we are fitting for");
		npol = getNparams(psr+p,offsetNp);
		offsetNp += npol;
		if (offsetNp > MAX_FIT)
		{
			printf("ERROR2: number of fitted parameters > MAX_FIT (%d)\n",MAX_FIT);
			exit(1);
		}

		logtchk("Complete determining which parameters we are fitting for");
		// Note that x gets free'd as xx
		x     = (double *)malloc(nobs_and_constraints*sizeof(double));
		if (x == NULL) {printf("Unable to allocate memory in doFit.C (x)\n"); exit(1);}
		y     = (double *)malloc(nobs_and_constraints*sizeof(double));
		if (y == NULL) {printf("Unable to allocate memory in doFit.C (y)\n"); exit(1);}
		xx[p]=x;
		yy[p]=y;
		sig   = (double *)malloc(nobs_and_constraints*sizeof(double));
		if (sig == NULL) {printf("Unable to allocate memory in doFit.C (sig)\n"); exit(1);}
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
					ip[p][count]  = i;
					// note that for an unweighted fit we later set this to 1.
					sig[count] = psr[p].obsn[i].toaErr*1e-6; /* Error in seconds */

					count++;
				}
			}
		}
		if (count == 0){
			printf("-----------------------------------------------\n");
			printf("ERROR: no observations to process for pulsar %s\n",psr[p].name);
			printf("Perhaps you have filtered out all the observations?\n");
			exit(1);
		}
		int last=count-1;
		// add constraints as extra pseudo observations
		// These point to non-existant observations after the nobs array
		// These are later caught by getParamDeriv.
		double constraint_sigma=1e-12;
		if (psr[p].fitMode == 0) constraint_sigma=1e-6;
		for (i=0; i < psr[p].nconstraints; i++){
			ip[p][count] = psr[p].nobs+i;
			x[count]=x[last];
			y[count]=0;
			sig[count]=constraint_sigma;
			count++;
		}

		

		if (covarFuncFile!=NULL && strcmp(covarFuncFile,"NULL")){
			// fit with a covariance function.
			logtchk("allocating memory for uinv ");
			uinvs[p]=malloc_uinv(count);
			logtchk("complete allocating memory for uinv");
			psr[p].fitMode = 1; // Note: forcing this to 1 as the Cholesky fit is a weighted fit
			logmsg("Doing a FULL COVARIANCE MATRIX fit");
		} else {
			// fit without covariance function.
			logtchk("allocating memory for uinv");
			uinvs[p]=malloc_blas(1,count); // store diagonal matrix as a 1xN
			logtchk("complete allocating memory for uinv");
			if(psr[p].fitMode == 0){
				logdbg("Doing an UNWEIGHTED fit");
				// unweighted fit - set sigma to 1 (except for constraints)
				for (i=0; i <= last; i++){
					sig[i]=1.0;
				}
			} else {
				logdbg("Doing a WEIGHTED fit");
			}
		}
		logtchk("Compute uinv");
		bool uselongdouble=false;
		// note that this works even for a non-cholesky fit.
		if(!uselongdouble)
		  {
		    getCholeskyMatrix(uinvs[p],covarFuncFile,psr+p,x,y,sig,count,psr[p].nconstraints,ip[p]);
		    logtchk("Completed computing uinv usign doubles");
		  }
		else
		  {
		    getCholeskyMatrix(uinvs[p],covarFuncFile,psr+p,x,y,sig,count,psr[p].nconstraints,ip[p]);
		    logtchk("Completed computing uinv usign long doubles");
		  }
		psr[p].nFit = count;
		psr[p].param[param_start].val[0] = newStart-0.001; 
		psr[p].param[param_finish].val[0] = newFinish+0.001;
		psr[p].param[param_start].paramSet[0] = 1;
		psr[p].param[param_finish].paramSet[0] = 1; 
		/*
		   logtchk("removing mean from the residuals??  (%.2f)",(clock()-clk)/(float)CLOCKS_PER_SEC);
		   for (i=0;i<psr[p].nobs;i++)
		   meanRes+=(long double)psr[p].obsn[i].residual;
		   meanRes/=(long double)psr[p].nobs;
		   for (i=0;i<psr[p].nobs;i++)
		   psr[p].obsn[i].residual-=meanRes;
		   logtchk("complete removing mean from the residuals??  (%.2f)",(clock()-clk)/(float)CLOCKS_PER_SEC);
		   */
		logdbg("Get constraint weghts");
		logtchk("Get Constraint weights");
		//		 printf("COMPUTING constraint weights for PSR %d %s\n",p,psr[p].name);
		computeConstraintWeights(psr+p);
		psr[p].nParam = npol;
		psr[p].nGlobal = nglobal;
		free(sig);
	}

	if (DO_GLOBAL_FIT){
		int n[npsr],nf[npsr];
		int ntot=0;
		int nobs=0;


		for (p=0;p<npsr;p++) {
			n[p]=psr[p].nFit;
			nobs+=n[p];
			nf[p]=psr[p].nParam;
			ntot+=nf[p];
		}
		ntot+=nglobal;
		logmsg("GLOBAL FIT MODE %d",ntot);
		double** cvm=malloc_uinv(ntot);
		double *val=(double*)malloc(sizeof(double)*ntot);
		double *error=(double*)malloc(sizeof(double)*ntot);

		TKleastSquares_global_pulsar(xx,yy,n,val,error,nf,nglobal,cvm,&chisq,globalFITfuncs,psr,tol,ip,1,uinvs,npsr);

		// update parameters

		if(nglobal > 0){
			logmsg("Update global parameters");
			//		 for (int i=0;i<nglobal;i++)
			//		   printf("PARAM: %g %g\n",val[i],error[i]);
			updateGlobalParameters(psr,npsr,val,error);
		}

		int offset=nglobal;
		for (p=0;p<npsr;p++) {
		    psr[p].robust=psr[0].robust; // robust fitting is controled by psr 0.
			psr[p].fitChisq = chisq; 
			// G. Hobbs, update so that the global fitting 
			// gives the same reduced chisq for each pulsar
			psr[p].fitNfree = nobs-ntot;//n[p]-nf[p]-nglobal;
			psr[p].globalNfit = ntot;
			psr[p].globalNoConstrain = nobs_noconstrain;
			logmsg("Update normal parameters for %s",psr->name);
			updateParameters(psr,p,val+offset,error+offset);
			offset+=nf[p];
		}

		// Ryan:  trying to output covariance matrix etc

		FILE *paramout;
		paramout = fopen("cvm.param", "w");

		int ii,jj;

		// global variables are always assigned first

		int iglobal, kglobal;
		int pglobal, qglobal;
		// print out covariance matrix of qifunc?

		//fprintf(stderr, "%d %d\n",nglobal, psr[0].nGlobal);
		for (ii=0;ii<nglobal;ii++)
		{

			iglobal = psr[0].fitParamGlobalI[ii];
			kglobal = psr[0].fitParamGlobalK[ii];
			//fprintf(stderr, "%d %d %d\n", ii,param_quad_ifunc_p, param_quad_ifunc_c);

			if ((iglobal == param_quad_ifunc_p) || (iglobal == param_quad_ifunc_c))
			{

				for(jj=0;jj<nglobal;jj++)
				{
					pglobal =  psr[0].fitParamGlobalI[jj];
					qglobal = psr[0].fitParamGlobalK[jj];

					if ((pglobal == param_quad_ifunc_p) || (pglobal  == param_quad_ifunc_c))
					{
						// let's not normalize it for now
						//fprintf(paramout, "%.3e ", cvm[ii][jj]/sqrt(cvm[ii][ii]*cvm[jj][jj]));
						if (paramout) fprintf(paramout, "%.3e \n", cvm[ii][jj]);

						//fprintf(stderr, "%d %d %d %s %.3e %.3e\n",ii,iglobal, kglobal,psr[0].param[iglobal].label[0], (float) psr[0].quad_ifuncV_p[kglobal], val[ii]);
					}
				}
				//		  fprintf(paramout, "\n");
			}

		}



		// all of the non-global parameters are named


		if (paramout) fclose(paramout);
		//exit(0);




		// Record the covariance matrix
		for (i=0;i<nglobal;i++)
		{
			for (j=0;j<nglobal;j++)
			{
				psr[0].covar[i][j] = cvm[i][j];
				//		  printf("COVAR STUFF %g\n",psr[0].covar[i][j]);
			}
		}

		/* Free the arrays created inside this section */
		free_uinv(cvm);
		free(val);
		free(error);
	} else {
		for (p=0;p<npsr;p++) { /* Loop over all the pulsars */
			/* Do the fit */
			npol=psr[p].nParam;

			val   = (double *)malloc(npol*sizeof(double));
			error = (double *)malloc(npol*sizeof(double));
			if (npol!=0) /* Are we actually  doing any fitting? */ 
			{ 

				logdbg("Doing the fit, npol = %d, nfit = %d",npol,psr[p].nFit);
				logtchk("doing the fit");
				TKleastSquares_single_pulsar(xx[p],yy[p],psr[p].nFit,val,error,npol,psr[p].covar,&chisq,
						FITfuncs,psr+p,tol,ip[p],1,uinvs[p]);
				logtchk("complete doing the fit");
				//	  svdfit(x,y,sig,psr[p].nFit,val,npol,u,v,w,&chisq,FITfuncs,&psr[p],tol,ip);
				logdbg("Complete fit: chisq = %f",(double)chisq);
				//	  printf("chisq = %g\n",chisq);
				psr[p].fitChisq = chisq; 
				psr[p].fitNfree = (psr[p].nFit-psr[p].nconstraints)-npol;

				/* Now update the parameters */
				logdbg("Updating the parameters");
				logtchk("updating the parameter values");
				updateParameters(psr,p,val,error);
				logtchk("complete updating the parameter values");
				logdbg("Completed updating the parameters");
			}
			/* Free the arrays created inside this section */
			free(error);
			free(val);

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
	}
	logtchk("freeing memory");
	for (p=0;p<npsr;p++) {
		free(yy[p]);      
		free(xx[p]);
		free(ip[p]);
		free_uinv(uinvs[p]);
	}
	free(ip);
	free(xx);
	free(yy);
	free(uinvs);

	logtchk("complete freeing memory");
}


int getNglobal(pulsar *psr,int npsr){
	int nGlobal=0;
	int i,k,j;
	// Add global parameters
	for (i=0;i<MAX_PARAMS;i++)
	{
		for (k=0;k<psr[0].param[i].aSize;k++)
		{
			if (psr[0].param[i].fitFlag[k]==2) {
				{
					if (i!=param_wave_om && i!= param_ifunc && i!=param_quad_om &&
							i!=param_tel_dx && i!= param_tel_dy && i!=param_tel_dz &&
					    i!=param_quad_ifunc_p && i!=param_quad_ifunc_c && i!=param_gwsingle && i!=param_wave_dm)
					{
						psr->fitParamI[nGlobal]  = i;
						psr->fitParamK[nGlobal]  = k;

						nGlobal++;
					}
				}
			}
		}
	}
	/* Add extra parameters for sinusoidal whitening */
	if (psr[0].param[param_wave_om].fitFlag[0]==2)
	{
		for (i=0;i<psr->nWhite*2;i++)
		{psr->fitParamI[nGlobal+i]  = param_wave_om; psr->fitParamK[nGlobal+i]  = i;}
		nGlobal+=psr[0].nWhite*2;
	}
	if (psr[0].param[param_wave_dm].fitFlag[0]==2)
	  {
	    for (i=0;i<psr->nWhite_dm*2;i++)
	      {psr->fitParamI[nGlobal+i]  = param_wave_dm; psr->fitParamK[nGlobal+i]  = i;}
	    nGlobal+=psr[0].nWhite_dm*2;
	  }

	



	if (psr[0].param[param_ifunc].fitFlag[0]==2){
	  if (psr[0].param[param_ifunc].val[0] == 0)
	    {
	      for (i=0;i<psr->ifuncN-1;i++)
		{psr->fitParamI[nGlobal+i]  = param_ifunc; psr->fitParamK[nGlobal+i]  = i;}
	      nGlobal+=psr[0].ifuncN-1;
	    }
	  else
	    {
	      for (i=0;i<psr->ifuncN;i++)
		{psr->fitParamI[nGlobal+i]  = param_ifunc; psr->fitParamK[nGlobal+i]  = i;}
	      nGlobal+=psr[0].ifuncN;
	    }
	}
	if (psr[0].param[param_quad_om].fitFlag[0]==2){
		for (i=0;i<psr->nQuad*4;i++)
		{psr->fitParamI[nGlobal+i]  = param_quad_om; psr->fitParamK[nGlobal+i]  = i;}
		nGlobal+=psr[0].nQuad*4;
	}

	if (psr[0].param[param_tel_dx].fitFlag[0]==2)
	{
		if (psr->param[param_tel_dx].val[0] < 2){
			for (i=0;i<psr->nTelDX;i++)
			{psr->fitParamI[nGlobal+i]  = param_tel_dx; psr->fitParamK[nGlobal+i]  = i;}
			nGlobal+=(psr[0].nTelDX);
		} else if (psr->param[param_tel_dx].val[0] == 2){
			for (i=0;i<psr->nTelDX-1;i++)
			{psr->fitParamI[nGlobal+i]  = param_tel_dx; psr->fitParamK[nGlobal+i]  = i;}
			nGlobal+=(psr[0].nTelDX-1);
		}
	}
	if (psr[0].param[param_tel_dy].fitFlag[0]==2)
	{     
		if (psr->param[param_tel_dy].val[0] < 2){
			for (i=0;i<psr->nTelDY;i++)
			{psr->fitParamI[nGlobal+i]  = param_tel_dy; psr->fitParamK[nGlobal+i]  = i;}
			nGlobal+=(psr[0].nTelDY);
		} else if (psr->param[param_tel_dy].val[0] < 2) {
			for (i=0;i<psr->nTelDY-1;i++)
			{psr->fitParamI[nGlobal+i]  = param_tel_dy; psr->fitParamK[nGlobal+i]  = i;}
			nGlobal+=(psr[0].nTelDY-1);
		}
	}
	if (psr[0].param[param_tel_dz].fitFlag[0]==2)
	{
		if (psr->param[param_tel_dz].val[0] < 2){
			for (i=0;i<psr->nTelDZ;i++)
			{psr->fitParamI[nGlobal+i]  = param_tel_dz; psr->fitParamK[nGlobal+i]  = i;}

			nGlobal+=(psr[0].nTelDZ);
		} else if (psr->param[param_tel_dz].val[0] == 2){
			for (i=0;i<psr->nTelDZ-1;i++)
			{psr->fitParamI[nGlobal+i]  = param_tel_dz; psr->fitParamK[nGlobal+i]  = i;}

			nGlobal+=(psr[0].nTelDZ-1);
		}
	}
	if (psr->param[param_quad_ifunc_p].fitFlag[0]==2)
	{
		for (i=0;i<psr->quad_ifuncN_p;i++)
		{psr->fitParamI[nGlobal+i]  = param_quad_ifunc_p; psr->fitParamK[nGlobal+i]  = i;}

		nGlobal+=(psr->quad_ifuncN_p);
	}
	if (psr->param[param_quad_ifunc_c].fitFlag[0]==2)
	{
		for (i=0;i<psr->quad_ifuncN_c;i++)
		{psr->fitParamI[nGlobal+i]  = param_quad_ifunc_c; psr->fitParamK[nGlobal+i]  = i;}
		nGlobal+=(psr->quad_ifuncN_c);
	}
	if (psr[0].param[param_gwsingle].fitFlag[0]==2)
	{
		for (i=0;i<4;i++)
		{psr->fitParamI[nGlobal+i]  = param_gwsingle; psr->fitParamK[nGlobal+i]  = i;}

		nGlobal+=(4);
	}
	return nGlobal;

}

int getNparams(pulsar *psr,int offset)
{
	int npol;
	int i,k;

	npol = 1;
	for (i=0;i<MAX_PARAMS;i++)
	{
		for (k=0;k<psr->param[i].aSize;k++)
		{
			if (psr->param[i].paramSet[k]==1 && psr->param[i].fitFlag[k]==1) {
				if (i!=param_start && i!=param_finish && i!=param_dmmodel && i!=param_gwsingle)
				{
					psr->fitParamI[npol+offset]  = i;
					psr->fitParamK[npol+offset]  = k;
					npol++;
				}
			}
		}
	}
	/* Add extra parameters for jumps */
	for (i=1;i<=psr->nJumps;i++)
	{
		if (psr->fitJump[i]==1)
		{
			psr->fitParamI[npol+offset]  = -1;
			psr->fitParamK[npol+offset]  = 0;
			npol++;
		}
	}
	/* Add extra parameters for sinusoidal whitening */
	if (psr->param[param_wave_om].fitFlag[0]==1)
	{
		printf("waveScale at this point = %d\n",psr->waveScale);
		if (psr->waveScale==1)
		{
			for (i=0;i<psr->nWhite*2-1;i++)
			{psr->fitParamI[npol+i+offset]  = param_wave_om; psr->fitParamK[npol+i+offset]  = i;}
			npol+=psr->nWhite*2-1;

		}
		else if (psr->waveScale==2)
		{
			for (i=0;i<psr->nWhite*4-1;i++)
			{psr->fitParamI[npol+i+offset]  = param_wave_om; psr->fitParamK[npol+i+offset]  = i;}

			npol+=psr->nWhite*4-1;
		}
		else
		{
			for (i=0;i<psr->nWhite*2-1;i++)
			{psr->fitParamI[npol+i+offset]  = param_wave_om; psr->fitParamK[npol+i+offset]  = i;}
			npol+=psr->nWhite*2-1;      
		}
	}
	 if (psr->param[param_wave_dm].fitFlag[0]==1)
	   {
	
	for (i=0;i<psr->nWhite_dm*2-1;i++)
	  {psr->fitParamI[npol+i+offset]  = param_wave_dm; psr->fitParamK[npol+i+offset]  = i;}
	npol+=psr->nWhite_dm*2-1;
	
      }


	if (psr->param[param_quad_om].fitFlag[0]==1)
	{
		for (i=0;i<psr->nQuad*4-1;i++)
		{psr->fitParamI[npol+i+offset]  = param_quad_om; psr->fitParamK[npol+i+offset]  = i;}

		npol+=(psr->nQuad*4)-1;
	}
	if (psr->param[param_ifunc].fitFlag[0]==1)
	{
	  if (psr->param[param_ifunc].val[0] == 0)
	    {
	      for (i=0;i<psr->ifuncN-1;i++)
		{psr->fitParamI[npol+i+offset]  = param_ifunc; psr->fitParamK[npol+i+offset]  = i;}
	      
	      npol+=(psr->ifuncN-1);
	    }
	  else
	    {
	      for (i=0;i<psr->ifuncN-1;i++)
		{psr->fitParamI[npol+i+offset]  = param_ifunc; psr->fitParamK[npol+i+offset]  = i;}
	      
	      npol+=(psr->ifuncN-1);
	    }
	}
	if (psr->param[param_clk_offs].fitFlag[0]==1)
	{
		for (i=0;i<psr->clkOffsN-1;i++)
		{psr->fitParamI[npol+i+offset]  = param_clk_offs; psr->fitParamK[npol+i+offset]  = i;}

		npol+=(psr->clkOffsN-1);
	}
	if (psr->param[param_tel_dx].fitFlag[0]==1 && psr->param[param_tel_dx].val[0] < 2)
	{
		for (i=0;i<psr->nTelDX-1;i++)
		{psr->fitParamI[npol+i+offset]  = param_tel_dx; psr->fitParamK[npol+i+offset]  = i;}

		npol+=(psr->nTelDX-1);
	}
	else if (psr->param[param_tel_dx].fitFlag[0]==1 && psr->param[param_tel_dx].val[0] == 2)
	{
		for (i=0;i<psr->nTelDX-2;i++)
		{psr->fitParamI[npol+i+offset]  = param_tel_dx; psr->fitParamK[npol+i+offset]  = i;}
		npol+=(psr->nTelDX-2);
	}
	if (psr->param[param_tel_dy].fitFlag[0]==1 && psr->param[param_tel_dy].val[0] < 2)
	{
		for (i=0;i<psr->nTelDY-1;i++)
		{psr->fitParamI[npol+i+offset]  = param_tel_dy; psr->fitParamK[npol+i+offset]  = i;}

		npol+=(psr->nTelDY-1);
	}
	else if (psr->param[param_tel_dy].fitFlag[0]==1 && psr->param[param_tel_dy].val[0] == 2)
	{
		for (i=0;i<psr->nTelDY-2;i++)
		{psr->fitParamI[npol+i+offset]  = param_tel_dy; psr->fitParamK[npol+i+offset]  = i;}


		npol+=(psr->nTelDY-2);
	}
	if (psr->param[param_tel_dz].fitFlag[0]==1 && psr->param[param_tel_dz].val[0] < 2)
	{
		for (i=0;i<psr->nTelDZ-1;i++)
		{psr->fitParamI[npol+i+offset]  = param_tel_dz; psr->fitParamK[npol+i+offset]  = i;}

		npol+=(psr->nTelDZ-1);
	}
	else if (psr->param[param_tel_dz].fitFlag[0]==1 && psr->param[param_tel_dz].val[0] == 2)
	{
		for (i=0;i<psr->nTelDZ-2;i++)
		{psr->fitParamI[npol+i+offset]  = param_tel_dz; psr->fitParamK[npol+i+offset]  = i;}

		npol+=(psr->nTelDZ-2);
	}
	if (psr->param[param_quad_ifunc_p].fitFlag[0]==1)
	{
		for (i=0;i<psr->quad_ifuncN_p-1;i++)
		{psr->fitParamI[npol+i+offset]  = param_quad_ifunc_p; psr->fitParamK[npol+i+offset]  = i;}

		npol+=(psr->quad_ifuncN_p-1);
	}
	if (psr->param[param_quad_ifunc_c].fitFlag[0]==1)
	{
		for (i=0;i<psr->quad_ifuncN_c-1;i++)
		{psr->fitParamI[npol+i+offset]  = param_quad_ifunc_c; psr->fitParamK[npol+i+offset]  = i;}

		npol+=(psr->quad_ifuncN_c-1);
	}
	/* Add extra parameters for DMMODEL fitting */
	if (psr->param[param_dmmodel].fitFlag[0]==1){
		for (i=0;i<psr->dmoffsDMnum;i++)
		{psr->fitParamI[npol+i+offset]  = -2; psr->fitParamK[npol+i+offset]  = i;}

		npol+=psr->dmoffsDMnum;
		for (i=0;i<psr->dmoffsCMnum;i++)
		{psr->fitParamI[npol+i+offset]  = -3; psr->fitParamK[npol+i+offset]  = i;}

		npol+=psr->dmoffsCMnum;
	}
	/* Add extra parameters for GW single source fitting */
	if (psr->param[param_gwsingle].fitFlag[0]==1)
	{
		for (i=0;i<4;i++)
		{psr->fitParamI[npol+i+offset]  = param_gwsingle; psr->fitParamK[npol+i+offset]  = i;}
		npol+=4; 
	}
	return npol;
}


void globalFITfuncs(double x,double afunc[],int ma,pulsar *psr,int ipos,int p){
	int i,j,k;
	int nglobal = psr[p].nGlobal;
	for (i=0;i<ma;i++) afunc[i]=0.0;

	// add global parameters.
	// Global fit
	int c=0;
	int kk;
	for (i=0;i<MAX_PARAMS;i++)
	{
		for (k=0;k<psr[p].param[i].aSize;k++)
		{
			if (psr[p].param[i].fitFlag[k] == 2)
			{
				//          afunc[c] = dotproduct(psr[p].posPulsar,psr[p].obsn[ipos].planet_ssb[4]);
				if (i==param_wave_om)
				{
					for (kk=0;kk<2*psr[0].nWhite;kk++)
					{
						afunc[c]= getParamDeriv(&psr[p],ipos,x+(double)psr[p].param[param_pepoch].val[0] - (double)psr[0].param[param_waveepoch].val[0],i,kk);
						c++;
					}
				}

				else if (i==param_wave_dm)
				  {
				    for (kk=0;kk<2*psr[0].nWhite_dm;kk++)
				      {
					afunc[c]= getParamDeriv(&psr[p],ipos,x+(double)psr[p].param[param_pepoch].val[0] - (double)psr[0].param[param_waveepoch_dm].val[0],i,kk);
					c++;
				      }
				  }

				else if(i==param_ifunc){
				  if (psr[p].param[param_ifunc].val[0] == 0)
				    {
				      for (j=0;j<psr[p].ifuncN-1;j++)
					{
					  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
					  //                printf("ifc=%d %d %g\n",counter,c,afunc[c]);
					  c++;
					}
				    }
				  else
				    {
				      for (j=0;j<psr[p].ifuncN;j++)
					{
					  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
					  //                printf("ifc=%d %d %g\n",counter,c,afunc[c]);
					  c++;
					}
				    }

				}
				else if(i==param_quad_ifunc_p){
					for (j=0;j<psr[p].quad_ifuncN_p;j++)
					{
						afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
						c++;
					}

				}
				else if(i==param_quad_ifunc_c){
					for (j=0;j<psr[p].quad_ifuncN_c;j++)
					{
						afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
						c++;
					}

				}
				else if(i==param_quad_om){
					for (j=0;j<psr[p].nQuad*4;j++)
					{
						afunc[c] = getParamDeriv(&psr[p],ipos,x+(double)psr[p].param[param_pepoch].val[0],i,j);
						//                printf("ifc=%d %d %g\n",counter,c,afunc[c]);
						c++;
					}
				}
				else if(i==param_tel_dx){
					if (psr[p].param[i].val[0]<2){
						for (j=0;j<psr[p].nTelDX;j++)
						{
							afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
							c++;
						}
					} else {
						for (j=0;j<psr[p].nTelDX-1;j++)
						{
							afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
							c++;
						}
					}
				}
				else if(i==param_tel_dy){
					if (psr[p].param[i].val[0]<2){
						for (j=0;j<psr[p].nTelDY;j++)
						{
							afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
							c++;
						}
					} else {
						for (j=0;j<psr[p].nTelDY-1;j++)
						{
							afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
							c++;
						}
					}
				}
				else if(i==param_tel_dz){
					if (psr[p].param[i].val[0]<2){
						for (j=0;j<psr[p].nTelDZ;j++)
						{
							afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
							c++;
						}
					}else{
						for (j=0;j<psr[p].nTelDZ-1;j++)
						{
							afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
							c++;
						}
					}
				}
				else if (i==param_gwsingle)
				{
					afunc[c] = getParamDeriv(&psr[p],ipos,x,i,0);
					c++;
					afunc[c] = getParamDeriv(&psr[p],ipos,x,i,1);
					c++;
					afunc[c] = getParamDeriv(&psr[p],ipos,x,i,2);
					c++;
					afunc[c] = getParamDeriv(&psr[p],ipos,x,i,3);
					c++;
				}
				else
				{
					afunc[c] = getParamDeriv(&psr[p],ipos,x,i,k);
					c++;
				}
			}
		}

	}

	FITfuncs(x,afunc+nglobal,ma-nglobal,psr,ipos,p); // the non-global parameters.
}


void FITfuncs(double x,double afunc[],int ma,pulsar *psr,int ipos,int ipsr)
{

	int i,n=0,k,j,l,found;
	psr+=ipsr; // this avoids having to change the code now we have unified with global fitting.

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
			if (psr->param[i].paramSet[k]==1 && psr->param[i].fitFlag[k]==1) /* If we are fitting for this parameter */
			{
			  if (i!=param_start && i!=param_finish)
			    {
				  
					logdbg("Fitting for %d (%s)",i,psr->param[i].label[k]);
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
					else if (i==param_wave_dm)
					  {
					    

					    fprintf(stderr, "here\n");
					      exit(0);
					    for (j=0;j<psr->nWhite_dm*2;j++)
					      afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
					    
					  }
					else if (i==param_quad_om)
					{
						for (j=0;j<psr->nQuad*4;j++)
							afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
					}
					else if (i==param_ifunc)
					{
					  if (psr->param[param_ifunc].val[0] == 0)
					    {
					      for (j=0;j<psr->ifuncN-1;j++)
						afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
					    }
					  else
					    {
					      for (j=0;j<psr->ifuncN;j++)
						afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
					    }
					}
					else if (i==param_ifunc)
					{
						for (j=0;j<psr->clkOffsN;j++)
							afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
					}
					else if (i==param_tel_dx)
					{
						if (psr->param[param_tel_dx].val[0]<2)
						{
							for (j=0;j<psr->nTelDX;j++)
								afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
						}
						else if (psr->param[param_tel_dx].val[0]==2)
						{
							for (j=0;j<psr->nTelDX-1;j++)
								afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
						}
					}
					else if (i==param_tel_dy)
					{
						if (psr->param[param_tel_dy].val[0]<2)
						{
							for (j=0;j<psr->nTelDY;j++)
								afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
						}
						else if (psr->param[param_tel_dy].val[0] == 2)
						{
							for (j=0;j<psr->nTelDY-1;j++)
								afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
						}
					}
					else if (i==param_tel_dz)
					{
						if (psr->param[param_tel_dz].val[0]<2)
						{
							for (j=0;j<psr->nTelDZ;j++)
								afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
						}
						else if (psr->param[param_tel_dz].val[0] == 2)
						{
							for (j=0;j<psr->nTelDZ-1;j++)
								afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
						}
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

						for (j=0;j<(int)psr->dmoffsDMnum;j++)
							afunc[n++] = dmf*getParamDeriv(psr,ipos,x,i,j);
						for (j=0;j<(int)psr->dmoffsCMnum;j++)
							afunc[n++] =     getParamDeriv(psr,ipos,x,i,j+psr->dmoffsDMnum);
					}
					else
					{
						afunc[n++] = getParamDeriv(psr,ipos,x,i,k);
						//		      printf("getParamDeriv: n = %d, result = %g, %d %d %d %g %s\n",n-1,afunc[n-1],ipos,i,k,x,psr->param[i].shortlabel[k]);
					}
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

	long double arg3,arg4;
	
	long double f1,f0,bindex;
	long double t2;
	
	t2 = x*86400.L;
	    
	arg3=t2*t2*t2;
	arg4=t2*t2*t2*t2;
	
	f0 = psr->param[param_f].val[0];
	f1 = psr->param[param_f].val[1];
	bindex = 0;//psr->param[param_brake].val[0];
	
	//fprintf(stderr, "%.3e %.3Le %.3Le\n", x, f0,f1);

	

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
		  {
		    afunc = x*24.0L*3600.0L/psr->param[param_f].val[0];
		    if (psr->param[param_brake].paramSet[0] ==1)
		      {
			afunc +=   (-bindex*f1*f1/f0/f0*arg3/6.L -2*bindex*(2*bindex-1)*f1*f1*f1/f0/f0/f0*arg4/24.)
			  *84000.L/f0;

			//printf("Braking part %Lg %Lg %Lg \n", bindex, arg3, arg4);
		      }
		    
		    
		  }
		else if (k==1)    /* Rotational frequency derivative */
		  {
		    afunc = 0.5L*x*x;
		    if (psr->param[param_brake].paramSet[0] ==1)
		      {
			afunc += (2*bindex*f1/f0*arg3/6.L + 3*bindex*(2*bindex-1)*f1*f1/f0/f0*arg4/24.L)
			  *86400.L*86400.L/f0;
		      }

		    
		  }
		else if (k==2)    /* Rotational frequency second derivative */
		  {
		    afunc = 1.0L/6.0L*x*x*x/1.0e9L;
		 
		  }

		else if (k==3)
		  {
		    afunc = 1.0L/24.0L*x/1.0e18L*x*x*x;
		  }
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
		else if (k==13)
			afunc = 1.0L/3628800.0L/11.0L/12.0L/13.0L/14.0L*powl(x,14.0L)/1.0e23L;
	}
	else if (i==param_brake)
	  {
	    
	    

	    
	    afunc = f1*f1/f0*arg3/6.0L + (4*bindex-1)*f1*f1*f1/f0/f0*arg4/24.L;



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
		if(psr->obsn[ipos].freq==0){
			// no DM fit for infinite frequency TOAS
			afunc=0;
		} else {
			double yrs;
			/* What about Doppler effect for the frequency -- change to barycentre?? */
			/* look at Blanford(?) paper */
			/* Should have a check to see if only one frequency exists in data
			   in which case fitting for DM does not make sense */
			if (k==0) 
				afunc = 1.0/(DM_CONST*powl(psr->obsn[ipos].freqSSB/1.0e6,2));
			else
			{
				yrs = (psr->obsn[ipos].sat - psr->param[param_dmepoch].val[0])/365.25;
				afunc = 1.0/(DM_CONST*pow(psr->obsn[ipos].freqSSB/1.0e6,2))*pow(yrs,k);
			}
		}
	}
	else if (i==param_dm_sin1yr)
		afunc = 1.0/(DM_CONST*powl(psr->obsn[ipos].freqSSB/1.0e6,2))*sin(2*M_PI/(365.25)*(psr->obsn[ipos].sat - psr->param[param_dmepoch].val[0]));
	else if (i==param_dm_cos1yr)
		afunc = 1.0/(DM_CONST*powl(psr->obsn[ipos].freqSSB/1.0e6,2))*cos(2*M_PI/(365.25)*(psr->obsn[ipos].sat - psr->param[param_dmepoch].val[0]));
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
	else if (i==param_fd)
	{
		afunc = pow(log(psr->obsn[ipos].freqSSB/1e9),k+1);
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
		else if (i==param_pmra) {
			longdouble t0 = ((x + psr->param[param_pepoch].val[0])
					- psr->param[param_posepoch].val[0]);
			afunc = re*cos(deltae)*cos(psrdec)*sin(psrra - alphae) * t0;
		}
		else if (i==param_pmdec) { /* pmdec */
			longdouble t0 = ((x + psr->param[param_pepoch].val[0])
					- psr->param[param_posepoch].val[0]);
			afunc = re*(cos(deltae)*sin(psrdec)*cos(psrra - alphae) - sin(deltae)*cos(psrdec))*t0;
		}
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
	 
	  double  Xoff = 0; // this Xoff is because x is referenced to pepoch not waveepoch!
	  
	  if (psr->param[param_waveepoch].paramSet[0]){
	    Xoff = psr->param[param_waveepoch].val[0] - psr->param[param_pepoch].val[0];
	  }
	  
	  //fprintf(stderr, "Xoff in paramderiv %.3e\n", Xoff);
	  
	  
	  double      om    = psr->param[param_wave_om].val[0];
	  if (psr->waveScale==0)
	    {
	      if (k%2==0) afunc = cos(om*(floor(k/2.0)+1)*(x-Xoff)); 
	      else        afunc = sin(om*(floor(k/2.0)+1)*(x-Xoff)); 
	      	  // printf("Value = %d %f %f %f %g\n",k,floor(k/2.0)+1,x,om,afunc);
	      
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

	   else if (i==param_wave_dm) /* Whitening procedure using sinusoids */
	     {
	  double  Xoff = 0; // this Xoff is because x is referenced to pepoch not waveepoch!
	  if (psr->param[param_waveepoch_dm].paramSet[0]){
	    Xoff = psr->param[param_waveepoch_dm].val[0] - psr->param[param_pepoch].val[0];
	  }
	  
	  double      om    = psr->param[param_wave_dm].val[0];
	  double freq = psr->obsn[ipos].freqSSB/1.0e6;
	  if (k%2==0) afunc = cos(om*(floor(k/2.0)+1)*(x-Xoff))/(DM_CONST*freq*freq); 
	  else        afunc = sin(om*(floor(k/2.0)+1)*(x-Xoff))/(DM_CONST*freq*freq); 
	
	  
	  
	  //  printf("Value = %d %f %f %f %g\n",k,floor(k/2.0)+1,x,om,afunc);

	  
	 
	     }


	else if (i==param_tel_dx)
	{
		double yoffs[MAX_TEL_DX];
		double sat = (double)psr->obsn[ipos].sat;
		if (psr->param[param_tel_dx].val[0] == -1)
		{
			printf("Here with %g %g\n",psr->posPulsar[0],psr->telDX_v[0]);
			if (strcmp(psr->obsn[ipos].telID,"STL_FBAT")==0)
				return psr->posPulsar[0]/SPEED_LIGHT*1000.0;
			else
				return 0;
		}
		if (psr->param[param_tel_dx].val[0] == 2)
		{
			if (psr->telDX_t[k] <=  psr->obsn[ipos].sat &&
					psr->telDX_t[k+1] > psr->obsn[ipos].sat)
				return psr->posPulsar[0];
			else 
				return 0;
		}
		else
		{

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
		}
		//      printf("afunc = %g\n",afunc);
	}
	else if (i==param_tel_dy)
	{
		double yoffs[MAX_TEL_DY];
		double sat = (double)psr->obsn[ipos].sat;


		if (psr->param[param_tel_dy].val[0] == -1)
		{
			if (strcmp(psr->obsn[ipos].telID,"STL_FBAT")==0)
				return psr->posPulsar[1]/SPEED_LIGHT*1000.0;
			else 
				return 0;
		}
		if (psr->param[param_tel_dy].val[0] == 2)
		{
			if (psr->telDY_t[k] <=  psr->obsn[ipos].sat &&
					psr->telDY_t[k+1] > psr->obsn[ipos].sat)
				return psr->posPulsar[1];
			else 
				return 0;
		}
		else
		{
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
		}
		//      printf("afunc = %g\n",afunc);
	}
	else if (i==param_tel_dz)
	{
		double yoffs[MAX_TEL_DZ];
		double sat = (double)psr->obsn[ipos].sat;
		if (psr->param[param_tel_dz].val[0] == -1)
		{
			if (strcmp(psr->obsn[ipos].telID,"STL_FBAT")==0)
				return psr->posPulsar[2]/SPEED_LIGHT*1000.0;
			else
				return 0;
		}

		if (psr->param[param_tel_dz].val[0] == 2)
		{
			// MUST SET SOME OF THESE TO ZERO!!
			if (psr->telDZ_t[k] <=  psr->obsn[ipos].sat &&
					psr->telDZ_t[k+1] > psr->obsn[ipos].sat)
				return psr->posPulsar[2];
			else 
				return 0;
		}
		else
		{
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
		}
		//      printf("afunc = %g\n",afunc);
	}
	else if (i==param_tel_x0) // satellite position
	{
		if (strcmp(psr->obsn[ipos].telID,"STL_FBAT")==0)
			return psr->posPulsar[0]/SPEED_LIGHT*1000.0;
		else
			return 0;
	}
	else if (i==param_tel_y0) // satellite position
	{
		if (strcmp(psr->obsn[ipos].telID,"STL_FBAT")==0)
			return psr->posPulsar[1]/SPEED_LIGHT*1000.0;
		else
			return 0;
	}
	else if (i==param_tel_z0) // satellite position
	{
		if (strcmp(psr->obsn[ipos].telID,"STL_FBAT")==0)
			return psr->posPulsar[2]/SPEED_LIGHT*1000.0;
		else
			return 0;
	}
	else if (i==param_tel_vx) // Velocity of satellite
	{
		long double dt = (psr->obsn[ipos].sat - psr->param[param_telEpoch].val[0])*SECDAY;
		if (strcmp(psr->obsn[ipos].telID,"STL_FBAT")==0)
			return psr->posPulsar[0]/SPEED_LIGHT*1000.0*dt;
		else
			return 0;
	}
	else if (i==param_tel_vy) // Velocity of satellite
	{
		long double dt = (psr->obsn[ipos].sat - psr->param[param_telEpoch].val[0])*SECDAY;
		if (strcmp(psr->obsn[ipos].telID,"STL_FBAT")==0)
			return psr->posPulsar[1]/SPEED_LIGHT*1000.0*dt;
		else
			return 0;
	}
	else if (i==param_tel_vz) // Velocity of satellite
	{
		long double dt = (psr->obsn[ipos].sat - psr->param[param_telEpoch].val[0])*SECDAY;
		if (strcmp(psr->obsn[ipos].telID,"STL_FBAT")==0)
			return psr->posPulsar[2]/SPEED_LIGHT*1000.0*dt;
		else
			return 0;
	}
	else if (i==param_clk_offs) /* Whitening procedure using interpolated function */
	{
		if (psr->param[param_clk_offs].val[0]==2)
		{
			int j;
			for (j=0;j<psr->clkOffsN-1;j++)
			{
				if (psr->obsn[ipos].sat >= psr->clk_offsT[j] &&
						psr->obsn[ipos].sat < psr->clk_offsT[j+1])
					return 1;
			}
			return 0;
		}
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
				int ioff;
				// find the pair we are between...
				for (ioff =0;ioff<psr->ifuncN;ioff++){
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
				if (ioff == psr->ifuncN)
				{
					printf("We should never reach this bit of code in doFit.C -- something wrong\n");
					exit(1);
				}
			}
		}
		else if (psr->param[param_ifunc].val[0]==0) // No interpolation
		{
		  if (psr->ifuncT[k] <= psr->obsn[ipos].sat &&
		      psr->ifuncT[k+1] > psr->obsn[ipos].sat)
		    return 1;
		  else
		    return 0;
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
		if (psr->param[param_raj].paramSet[1] == 1)
			lambda_p = (double)psr->param[param_raj].val[1];
		else
			lambda_p = (double)psr->param[param_raj].val[0];

		if (psr->param[param_raj].paramSet[1] == 1)
			beta_p   = (double)psr->param[param_decj].val[1];
		else
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

			if (psr->gwsrc_psrdist > 0) // Add in the pulsar term
			{
				if (k%4==0) afunc      -=  resp*sinl(omega_g*time-(1-cosTheta)*psr->gwsrc_psrdist/SPEED_LIGHT*omega_g)/(2.0L*omega_g*(1.0L-cosTheta)); // aplus_im
				else if (k%4==1) afunc -=  resp*cosl(omega_g*time-(1-cosTheta)*psr->gwsrc_psrdist/SPEED_LIGHT*omega_g)/(2.0L*omega_g*(1.0L-cosTheta)); // aplus_im
				else if (k%4==2) afunc -=  resc*sinl(omega_g*time-(1-cosTheta)*psr->gwsrc_psrdist/SPEED_LIGHT*omega_g)/(2.0L*omega_g*(1.0L-cosTheta)); // aplus_im
				else if (k%4==3) afunc -=  resc*cosl(omega_g*time-(1-cosTheta)*psr->gwsrc_psrdist/SPEED_LIGHT*omega_g)/(2.0L*omega_g*(1.0L-cosTheta)); // aplus_im
			}


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

	 else if (i== param_gwb_amp)
     {
       long double dt;
       long double prefac;
       if (psr->param[param_gwb_amp].paramSet[1]==1)
	 {
	   dt = (psr->obsn[ipos].bbat - psr->gwb_epoch)/psr->gwb_width;
	   prefac = dt*exp( (double) -dt*dt/2.);
	   
	   if (k==0)
	     {
		 afunc = psr->gwb_geom_p*prefac;
	       }
	     else if (k==1)
	       {
		 afunc = psr->gwb_geom_c*prefac;
	       }
	     else
	       { 
		 afunc = 0;
	       }
	 }
     }
	

	

	else if (i==param_gwm_amp)
	{
		long double dt;

		if (psr->param[param_gwm_amp].paramSet[1]==1){
			dt = (psr->obsn[ipos].bbat - psr->gwm_epoch)*86400.0L;
			if (dt > 0)
			{
				if (k==0)
					afunc = psr->quad_ifunc_geom_p*dt;
				else if (k==1)
					afunc = psr->quad_ifunc_geom_c*dt;
			}
			else afunc = 0;
		}
		else
		{
			double n1,n2,n3;
			double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
			double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;
			double cosTheta,omega_g;
			long double resp,resc,res_r,res_i;
			double theta_p,theta_g,phi_p,phi_g;
			double lambda_p,beta_p,lambda,beta;
			long double time;
			double g1,g2,g3;

			time    = (psr->obsn[ipos].bbat - psr->gwm_epoch)*86400.0L;

			if (psr->param[param_raj].paramSet[1] == 1)
				lambda_p = (double)psr->param[param_raj].val[1];
			else
				lambda_p = (double)psr->param[param_raj].val[0];

			if (psr->param[param_decj].paramSet[1] == 1)
				beta_p   = (double)psr->param[param_decj].val[1];
			else
				beta_p   = (double)psr->param[param_decj].val[0];

			lambda   = psr->gwm_raj;
			beta     = psr->gwm_decj;

			// GW vector
			g1 = -cosl(lambda)*cosl(beta);
			g2 = -sinl(lambda)*cosl(beta);
			g3 = -sinl(beta);

			// Pulsar vector
			n1 = cosl(lambda_p)*cosl(beta_p);
			n2 = sinl(lambda_p)*cosl(beta_p);
			n3 = sinl(beta_p);
			cosTheta = -(cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
					sinl(beta)*sinl(beta_p));

			/* Only has effect after the glitch epoch */
			if (psr->obsn[ipos].sat >= psr->gwm_epoch)
			{
				long double dt,scale;
				double cos2Phi;
				double cosPhi;
				double l1,l2,l3,n5,m1,m2,m3;
				double beta_m;
				double d1,d2,d3,md;
				double a1,a2,a3,ma;

				//   if  (g3 != 0) 
				//	   {beta_m = atan2(-cos(beta)*cos(lambda-psr->gwm_phi),sin(beta));}
				//  else  
				//      {beta_m = atan2(sinl(psr->gwm_phi),cosl(psr->gwm_phi));
				//       psr->gwm_phi = lambda + 1.5708;}
				//  m1 = cosl(psr->gwm_phi)*cosl(beta_m);
				//  m2 = sinl(psr->gwm_phi)*cosl(beta_m);
				//  m3 = sinl(beta_m);

				if (beta == 0.0 )
				{
					d1 = 0.0;
					d2 = 0.0;
					d3 = 1.0;
				}

				if ( beta > 0)
				{
					d1 = g1*cosl(0.5*M_PI - beta);
					d2 = g2*cosl(0.5*M_PI - beta);
					d3 = 1.0 + g3*cos(0.5*M_PI - beta);
					md = sqrt(d1*d1 + d2*d2 + d3*d3);
					d1 = d1/md;
					d2 = d2/md;
					d3 = d3/md;
					/*covert d to unit vector */
				} 
				else if (beta < 0) 
				{
					d1 = g1*cosl(-0.5*M_PI - beta);
					d2 = g2*cosl(-0.5*M_PI - beta);
					d3 = -1.0 + g3*cos(-0.5*M_PI - beta);
					md = sqrt(d1*d1 + d2*d2 + d3*d3);
					d1 = d1/md;
					d2 = d2/md;
					d3 = d3/md;
				} 

				//if (g2*d3-d2*g3 != 0)
				// {
				//  a1 = 1.0; 
				//  a2 = (d1*g3-g1*d3)/(g2*d3-d2*g3);
				//  a3 = (g2*d1-g1*d2)/(g3*d2-g2*d3); 
				// }
				//else if (g1*d3-d1*g3 != 0)
				// {
				//  a1 = (g3*d2-d3*g2)/(g1*d3-g3*d1); 
				//  a2 = 1.0;
				//  a3 = (g1*d2-d1*g2)/(g3*d1-d1*d3);
				// }
				//else if (d2*g1-g2*d1 != 0)			
				// {
				//  a1 = (g2*d3-d2*g3)/(d2*g1-g2*d1); 
				//  a2 = (g1*d3-d1*g3)/(d1*g2-g1*d2);
				//  a3 =1.0; 
				// }
				a1 =  (d2*g3-d3*g2);
				a2 =  (d3*g1-d1*g3);
				a3 =  (d1*g2-d2*g1);

				/* conver it to unit vector */
				ma = sqrt(a1*a1 +a2*a2 + a3*a3);
				a1 = a1/ma;
				a2 = a2/ma;
				a3 = a3/ma;

				/* polarisation vector of GW source */
				m1 = d1*cosl(psr->gwm_phi)	+ a1*sinl(psr->gwm_phi);   
				m2 = d2*cosl(psr->gwm_phi)	+ a2*sinl(psr->gwm_phi);
				m3 = d3*cosl(psr->gwm_phi)	+ a3*sinl(psr->gwm_phi);

				if  (cosTheta != 1.0 && cosTheta != -1.0)
				{g1 = g1*cosTheta; 
					g2 = g2*cosTheta;
					g3 = g3*cosTheta;
					l1 = n1 - g1;
					l2 = n2 - g2;
					l3 = n3 - g3;
					cosPhi = (l1*m1 + l2*m2 + l3*m3)/sqrt(l1*l1 + l2*l2 + l3*l3);
					//		 if  (cosPhi >= 1.0/sqrt(2.0))
					cos2Phi = 2*cosPhi*cosPhi - 1.0;
					//		 else
					//		     cos2Phi = 2*sqrt(1.0 - cosPhi*cosPhi)*sqrt(1.0 - cosPhi*cosPhi) - 1.0;
				}
				else 
				{cos2Phi = 0;}

				dt = (psr->obsn[ipos].sat - psr->gwm_epoch)*86400.0;
				scale = -0.5*cos2Phi*(1-cosTheta);
				//	    scale=1;
				afunc = scale*dt;
			}
			else
				afunc = 0;
		}
	}
	else if (i==param_dmmodel)
	{
		int N=psr->dmoffsDMnum;
		double sat = (double)psr->obsn[ipos].sat;
		double yoffs[MAX_IFUNC];
		double* mjd=psr->dmoffsDM_mjd;
		if (k >= N){
			// if we are at k >= N, then we are in the CM, not the DM
			// so swap the variables over.
			mjd=psr->dmoffsCM_mjd;
			k-=N;
			N=psr->dmoffsCMnum;
		}


		for (int ioff =0;ioff<N;ioff++){
			if (ioff==k){
				yoffs[ioff]=1;
			} else {
				yoffs[ioff]=0;
			}
		}

		if (sat < mjd[0]){
			// we are before the first jump
			// so our gradient is just the zeroth offset.
			afunc = yoffs[0];
		} else if(sat > mjd[N-1]){
			afunc = yoffs[N-1];
		} else{
			// find the pair we are between...
			for (int ioff =0;ioff<N;ioff++){
				if(sat >= mjd[ioff] && sat < mjd[ioff+1]){
					double x1 = mjd[ioff];
					double x2 = mjd[ioff+1];
					double x = (sat-x1)/(x2-x1);
					double y1=yoffs[ioff];
					double y2=yoffs[ioff+1];
					afunc = (y2-y1)*x + y1;
					break;
				}
			}
		}

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
	else if (strcmp(psr->binaryModel,"BTX")==0) 
		afunc = BTXmodel(psr,0,ipos,i,k);
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
	else if (strcmp(psr->binaryModel,"T2-PTA")==0) 
		afunc = T2_PTAmodel(psr,0,ipos,i,k);	
	else if (strcmp(psr->binaryModel,"DDH")==0)
		afunc = DDHmodel(psr,0,ipos,i);
	else if (strcmp(psr->binaryModel,"ELL1H")==0)
		afunc = ELL1Hmodel(psr,0,ipos,i);


	return afunc;
}


void updateGlobalParameters(pulsar* psr,int npsr, double* val,double* error){
	int offset=0;
	int i,j,k,p;
	// update global parameters
	for (i=0;i<MAX_PARAMS;i++)
	{
		for (k=0;k<psr[0].param[i].aSize;k++)
		{
			if (psr[0].param[i].fitFlag[k] == 2)
			{
				//			printf("Have global parameter %d %d %d %g %g\n",i,param_wave_om,param_ifunc,val[offset],error[offset]);
				if (i==param_wave_om)
				{
					int kk;
					for (kk=0;kk<psr[0].nWhite;kk++)
					{
						for (p=0;p<npsr;p++)
						{
							psr[p].wave_cos[kk]  -= val[offset];
							psr[p].wave_cos_err[kk] = error[offset];
							if (p==0) printf("Have wave %d %d %g %g\n",offset,kk,val[offset],error[offset]);
						}
						offset++;
						for (p=0;p<npsr;p++)
						{
							psr[p].wave_sine[kk] -= val[offset];
							psr[p].wave_sine_err[kk] = error[offset];
						}
						offset++;
					}
					offset--;
				}
				else if (i==param_wave_dm) /* Whitening procedure using sinusoids */
				  {
				    int k;
				    //exit(0);
				    
				    
				    for (k=0;k<psr[p].nWhite_dm;k++)
				      {
					//fprintf(stderr, "%.3e\n", val[j]);
					psr[p].wave_cos_dm[k]  -= val[j]; 
					psr[p].wave_cos_dm_err[k] = error[j]; j++;
					psr[p].wave_sine_dm[k] -= val[j]; 
					psr[p].wave_sine_dm_err[k] = error[j]; j++;	      
				      }
				    j--;  
				    
				    //exit(0);
				  }


			

				else if (i==param_gwm_amp)
				{
					//			  printf("In here with offset = %d\n",offset);
					for (p=0;p<npsr;p++)
					{			      
						if (psr[p].param[param_gwm_amp].paramSet[1]==1)
						{
							psr[p].param[i].val[k] -= val[offset];
							psr[p].param[i].err[k] = error[offset];
						}
						else
						{
							psr[p].param[i].val[0] -= val[offset];
							psr[p].param[i].err[0] = error[offset];
						}
					}
					//			  printf("Setting %d %g %g\n",offset,val[offset],error[offset]);
					//			  offset++;

				}
				
				else if (i==param_gwb_amp)
					  {
					    //			  printf("In here with offset = %d\n",offset);
					    for (p=0;p<npsr;p++)
					      {			      
						if (psr[p].param[param_gwb_amp].paramSet[1]==1)
						  {
						    psr[p].param[i].val[k] -= val[offset];
						    psr[p].param[i].err[k] = error[offset];
						  }
						else
						  {
						    psr[p].param[i].val[0] -= val[offset];
						    psr[p].param[i].err[0] = error[offset];
						  }
					      }
			  //			  printf("Setting %d %g %g\n",offset,val[offset],error[offset]);
			  //			  offset++;

			}


				else if(i==param_ifunc) {
					printf("Updating %d point\n",psr[0].ifuncN);
					if (psr[0].param[param_ifunc].val[0]==0)
					  {
					    for (j=0;j<psr[0].ifuncN-1;j++)
					      {
						printf("Updating %d %g\n",offset,val[offset]);
						for (p=0;p<npsr;p++)
						  {
						    psr[p].ifuncV[j]-=val[offset];
						    psr[p].ifuncE[j]=error[offset];
						  }
						offset++;
					      }
					  }
					else
					  {
					    for (j=0;j<psr[0].ifuncN;j++)
					      {
						printf("Updating %d %g\n",offset,val[offset]);
						for (p=0;p<npsr;p++)
						  {
						    psr[p].ifuncV[j]-=val[offset];
						    psr[p].ifuncE[j]=error[offset];
						  }
						offset++;
					      }
					  }
					offset--;
				}
				else if (i==param_quad_ifunc_p)
				{
					for (j=0;j<psr[0].quad_ifuncN_p;j++)
					{
						printf("Updating %g\n",val[offset]);
						for (p=0;p<npsr;p++)
						{
							psr[p].quad_ifuncV_p[j]-=val[offset];
							psr[p].quad_ifuncE_p[j]=error[offset];
						}
						offset++;
					}
					offset--;
				}
				else if (i==param_quad_ifunc_c)
				{
					for (j=0;j<psr[0].quad_ifuncN_c;j++)
					{
						printf("Updating %g\n",val[offset]);
						for (p=0;p<npsr;p++)
						{
							psr[p].quad_ifuncV_c[j]-=val[offset];
							psr[p].quad_ifuncE_c[j]=error[offset];
						}
						offset++;
					}
					offset--;
				}
				else if(i==param_quad_om) {
					printf("Updating %d point\n",psr[0].ifuncN);
					for (j=0;j<psr[0].nQuad;j++)
					{
						printf("Updating %d %g\n",offset,val[offset]);
						for (p=0;p<npsr;p++)
						{
							psr[p].quad_aplus_r[j]    -= val[offset];
							psr[p].quad_aplus_i[j]    -= val[offset+1];
							psr[p].quad_across_r[j]   -= val[offset+2];
							psr[p].quad_across_i[j]   -= val[offset+3];
							psr[p].quad_aplus_r_e[j]   = error[offset];
							psr[p].quad_aplus_i_e[j]   = error[offset+1];
							psr[p].quad_across_r_e[j]  = error[offset+2];
							psr[p].quad_across_i_e[j]  = error[offset+3];
						}
						offset+=4;
					}
					offset--;
				}
				else if (i==param_tel_dx)
				{
					if (psr[0].param[i].val[0]<2){
						for (j=0;j<psr[0].nTelDX;j++)
						{
							printf("Setting a: %d %g\n",j,val[offset]);
							for (p=0;p<npsr;p++)
							{
								psr[p].telDX_v[j]-=val[offset];
								psr[p].telDX_e[j]=error[offset];
							}
							offset++;
						}
						offset--;
					}  else {
						for (j=0;j<psr[0].nTelDX-1;j++)
						{
							printf("Setting b: %d %g\n",j,val[offset]);
							for (p=0;p<npsr;p++)
							{
								psr[p].telDX_v[j]-=val[offset];
								psr[p].telDX_e[j]=error[offset];
							}
							offset++;
						}
						offset--;

					}
				}
				else if (i==param_tel_dy)
				{
					if (psr[0].param[i].val[0]<2){
						for (j=0;j<psr[0].nTelDY;j++)
						{
							printf("Setting c: %d %g\n",j,val[offset]);
							for (p=0;p<npsr;p++)
							{
								psr[p].telDY_v[j]-=val[offset];
								psr[p].telDY_e[j]=error[offset];
							}
							offset++;
						}
						offset--;
					} else {
						for (j=0;j<psr[0].nTelDY-1;j++)
						{
							printf("Setting d: %d %g\n",j,val[offset]);
							for (p=0;p<npsr;p++)
							{
								psr[p].telDY_v[j]-=val[offset];
								psr[p].telDY_e[j]=error[offset];
							}
							offset++;
						}
						offset--;
					}
				}
				else if (i==param_tel_dz)
				{
					if (psr[0].param[i].val[0]<2){
						for (j=0;j<psr[0].nTelDZ;j++)
						{
							printf("Setting e: %d %g\n",j,val[offset]);
							for (p=0;p<npsr;p++)
							{
								psr[p].telDZ_v[j]-=val[offset];
								psr[p].telDZ_e[j]=error[offset];
							}
							offset++;
						}
						offset--;
					} else {
						for (j=0;j<psr[0].nTelDZ-1;j++)
						{
							for (p=0;p<npsr;p++)
							{
								psr[p].telDZ_v[j]-=val[offset];
								psr[p].telDZ_e[j]=error[offset];
							}
							offset++;
						}
						offset--;

					}
				}
				else if (i==param_gwsingle)
				{
					for (p=0;p<npsr;p++)
					{
						psr[p].gwsrc_aplus_r -= val[offset];
						psr[p].gwsrc_across_r -= val[offset+1];
						psr[p].gwsrc_aplus_r_e = error[offset];
						psr[p].gwsrc_across_r_e = error[offset+1];
						psr[p].gwsrc_aplus_i -= val[offset+2];
						psr[p].gwsrc_across_i -= val[offset+3];
						psr[p].gwsrc_aplus_i_e = error[offset+2];
						psr[p].gwsrc_across_i_e = error[offset+3];
					}
					offset+=3;
				}
				else
				{
					for (p=0;p<npsr;p++)
					{
						if (i==param_telx || i==param_tely || i==param_telz)
							psr[p].param[i].val[k] -= val[offset];
						else
							psr[p].param[i].val[k] += val[offset];
						psr[p].param[i].err[k] = error[offset];
					}
				}
				offset++;
			}
		}
	}


}
void updateParameters(pulsar *psr,int p,double *val,double *error)
{
	int i,j,k;
	logdbg("Updating parameters");
	psr[p].offset = val[0];
	psr[p].offset_e = error[0];
	j=1;
	for (i=0;i<MAX_PARAMS;i++)
	{
		for (k=0;k<psr[p].param[i].aSize;k++)
		{
			if (psr[p].param[i].paramSet[k]==1 && psr[p].param[i].fitFlag[k]==1 && (i!=param_start && i!=param_finish))
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
				else if (i==param_dm || i==param_px || i==param_fddc || i==param_fddi || i==param_dmassplanet || i==param_dmx || i==param_fd || i==param_dm_sin1yr || i==param_dm_cos1yr)
				{
					psr[p].param[i].val[k] += val[j];
					psr[p].param[i].err[k]  = error[j];
					// The following lines break the -dmo tim-file option and have therefore been disabled. 
					// As far as I know, they don't really have an effect anyway.
					//                                               JPWV, 08.05.2014
					// This is slow - should be a better approach
					// if (i==param_dm){
					//   psr[p].dmOffset+=val[j];
					// }
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
				
				  /*
				  int k;
					for (k=0;k<psr[p].nWhite;k++)
					{
						psr[p].wave_cos[k]  -= val[j]; 
						psr[p].wave_cos_err[k] = error[j]; j++;
						psr[p].wave_sine[k] -= val[j]; 
						psr[p].wave_sine_err[k] = error[j]; j++;	      
					}
					if (psr->waveScale==2) // Ignore the non-frequency derivative terms
					{
						for (k=0;k<psr[p].nWhite;k++)
						{
							printf("Ignoring cos %g %g\n",val[j],error[j]); j++;
							printf("Ignoring sin %g %g\n",val[j],error[j]); j++;
						}
						//		      j+=psr->nWhite*2; 
					}
					j--;

				  */


				  
				  if (psr->waveScale != 2)
				    {
				      
				      for (k=0;k<psr[p].nWhite;k++)
					{
					  fprintf(stderr, "%d %.3e \n", j, val[j]);
					  psr[p].wave_cos[k]  -= val[j]; 
					  psr[p].wave_cos_err[k] = error[j]; j++;
					  psr[p].wave_sine[k] -= val[j]; 
					  psr[p].wave_sine_err[k] = error[j]; j++;	      
					}
				    }
				  else if(psr->waveScale==2) // Ignore the non-frequency derivative terms
				    {
				      
				      
				      for (k=0;k<psr[p].nWhite;k++)
					{
					  //printf("Ignoring cos %g %g\n",val[j],error[j]); j++;
					  //printf("Ignoring sin %g %g\n",val[j],error[j]); j++;
					  psr[p].wave_cos_dm[k]-= val[j];
					  psr[p].wave_cos_dm_err[k]= error[j]; j++;
					  psr[p].wave_sine_dm[k] -= val[j];
					  psr[p].wave_sine_dm_err[k]= error[j]; j++;
					}
				      // now do timing noise terms
				      for (k=0;k<psr[p].nWhite;k++)
					{
					  psr[p].wave_cos[k]  -= val[j]; 
					  psr[p].wave_cos_err[k] = error[j]; j++;
					  psr[p].wave_sine[k] -= val[j]; 
					  psr[p].wave_sine_err[k] = error[j]; j++;	      
					}
				      
				      //		      j+=psr->nWhite*2; 
				    }
				  j--;
				  
				}	


				else if (i==param_wave_dm) /* Whitening procedure using sinusoids */
				  {
				    int k;
				    //exit(0);
				    
				    
				    for (k=0;k<psr[p].nWhite_dm;k++)
				      {
					//fprintf(stderr, "%.3e\n", val[j]);
					psr[p].wave_cos_dm[k]  -= val[j]; 
					psr[p].wave_cos_dm_err[k] = error[j]; j++;
					psr[p].wave_sine_dm[k] -= val[j]; 
					psr[p].wave_sine_dm_err[k] = error[j]; j++;	      
				      }
				    j--;  
				    
				    //exit(0);
				  }

	  
				else if (i==param_ifunc) 
				{
					int k;
					for (k=0;k<psr[p].ifuncN;k++)
					{
						psr[p].ifuncV[k] -= val[j];
						psr[p].ifuncE[k] = error[j];
						j++;
					}
				}		  
				else if (i==param_clk_offs) 
				{
					int k;
					for (k=0;k<psr[p].clkOffsN;k++)
					{
						psr[p].clk_offsV[k] -= val[j];
						psr[p].clk_offsE[k] = error[j];
						j++;
					}
				}		  
				else if (i==param_tel_dx) 
				{
					int k;
					if (psr[p].param[param_tel_dx].val[0] == -1)
					{
						printf("Updating with %g %g\n",val[j],error[j]);
						psr[p].telDX_v[0] -= val[j]; //*SPEED_LIGHT/1000.0;
						psr[p].telDX_e[0] = error[j]; //*SPEED_LIGHT/1000.0;
					}
					else if (psr[p].param[param_tel_dx].val[0] < 2)
					{
						for (k=0;k<psr[p].nTelDX;k++)
						{
							psr[p].telDX_v[k] -= val[j];
							psr[p].telDX_e[k] = error[j];
							j++;
						}
					}
					else
					{
						for (k=0;k<psr[p].nTelDX-1;k++)
						{
							psr[p].telDX_v[k] -= val[j];
							psr[p].telDX_e[k] = error[j];
							j++;
						}
					}
				}
				else if (i==param_tel_dy) 
				{
					int k;
					if (psr[p].param[param_tel_dy].val[0] == -1)
					{
						psr[p].telDY_v[0] -= val[j]; //*SPEED_LIGHT/1000.0;
						psr[p].telDY_e[0] = error[j]; //*SPEED_LIGHT/1000.0;
					}
					else if (psr[p].param[param_tel_dy].val[0] < 2)
					{
						for (k=0;k<psr[p].nTelDY;k++)
						{
							psr[p].telDY_v[k] -= val[j];
							psr[p].telDY_e[k] = error[j];
							j++;
						}
					}
					else
					{
						for (k=0;k<psr[p].nTelDY-1;k++)
						{
							psr[p].telDY_v[k] -= val[j];
							psr[p].telDY_e[k] = error[j];
							j++;
						}
					}


				}		  
				else if (i==param_tel_dz) 
				{
					int k;
					if (psr[p].param[param_tel_dz].val[0] == -1)
					{
						psr[p].telDZ_v[0] -= val[j]; //*SPEED_LIGHT/1000.0;
						psr[p].telDZ_e[0] = error[j]; //*SPEED_LIGHT/1000.0;
					}
					else if (psr[p].param[param_tel_dz].val[0] < 2)
					{
						for (k=0;k<psr[p].nTelDZ;k++)
						{
							psr[p].telDZ_v[k] -= val[j];
							psr[p].telDZ_e[k] = error[j];
							j++;
						}
					}
					else
					{
						for (k=0;k<psr[p].nTelDZ-1;k++)
						{
							psr[p].telDZ_v[k] -= val[j];
							psr[p].telDZ_e[k] = error[j];
							j++;
						}
					}
				}		  
				else if (i==param_quad_ifunc_p) 
				{
					int k;
					for (k=0;k<psr[p].quad_ifuncN_p;k++)
					{
						psr[p].quad_ifuncV_p[k] -= val[j];
						psr[p].quad_ifuncE_p[k] = error[j];
						j++;
					}
				}		  
				else if (i==param_quad_ifunc_c) 
				{
					int k;
					for (k=0;k<psr[p].quad_ifuncN_c;k++)
					{
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
				else if (i==param_gwm_amp)
				{
					printf("Should I be in here???\n");
					if (psr[p].param[param_gwm_amp].paramSet[1]==1)
					{
						psr[p].param[i].val[k] -= val[j];
						psr[p].param[i].err[k] = error[j];
						j++;
					}
					else
					{
						psr[p].param[i].val[0] -= val[j];
						psr[p].param[i].err[0] = error[j];
						j++;
					}
				}

				else if (i==param_gwb_amp)
				  {
				    printf("Should I be in here???\n");
				    if (psr[p].param[param_gwb_amp].paramSet[1]==1)
				      {
					psr[p].param[i].val[k] -= val[j];
					psr[p].param[i].err[k] = error[j];
					j++;
				      }
				    else
				      {
					psr[p].param[i].val[0] -= val[j];
					psr[p].param[i].err[0] = error[j];
					j++;
				      }
				  }
				else if (i==param_brake)
				  {
				    psr[p].param[i].val[0] -= val[j];
				    psr[p].param[i].err[0] = error[j];
				    
				  }


				else if (i==param_dmmodel)
				{
					for (int k=0;k<psr[p].dmoffsDMnum;k++)
					{
						psr[p].dmoffsDM[k] += val[j];
						psr[p].dmoffsDM_error[k] = error[j];
						j++;
					}
					for (int k=0;k<psr[p].dmoffsCMnum;k++){

						psr[p].dmoffsCM[k] = val[j];
						psr[p].dmoffsCM_error[k] =  error[j];
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
				else if (strcmp(psr[p].binaryModel,"BTX")==0)
					updateBTX(&psr[p],val[j],error[j],i,k);
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
				else if (strcmp(psr[p].binaryModel,"T2-PTA")==0)
					updateT2_PTA(&psr[p],val[j],error[j],i,k);
				else if (strcmp(psr[p].binaryModel,"DDH")==0)
					updateDDH(&psr[p],val[j],error[j],i);
				else if (strcmp(psr[p].binaryModel,"ELL1H")==0)
				  updateELL1H(&psr[p],val[j],error[j],i);
			
				j++; /* Increment position in fit list */
	
		}
		}
	}
	if (strcmp(psr[p].binaryModel,"DDGR")==0) 
		DDGRmodel(psr,0,0,-2);  /* Update GR parameters */	  

	logdbg("Updating jumps; nJumps = %d",psr[p].nJumps);
	/* Now check jumps */
	for (i=1;i<=psr[p].nJumps;i++)
	{
		logdbg("%d fitJump = %d",i,psr[p].fitJump[i]);
		if (psr[p].fitJump[i]==1)
		{
			logdbg("%d Jump changes",i);
			logdbg("value = %g",(double)val[j]);
			logdbg("error = %g",(double)error[j]);
			psr[p].jumpVal[i] += -val[j];
			psr[p].jumpValErr[i] = error[j];
			j++;
		}
		/*	      printf("Have jumps %g %g\n",(double)val[j],error[j][j]); */
	}
	logdbg("Complete updating parameters");
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
                case constraint_dmmodel_dm1:
                       return consFunc_dmmodel_dm1(psr,i,k);
		case constraint_dmmodel_cw_3:
			order++;
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

			// DMMODEL annual terms
		case constraint_dmmodel_cw_px:
			order++;
		case constraint_dmmodel_cw_year_cos2:
			order++;
		case constraint_dmmodel_cw_year_sin2:
			order++;
		case constraint_dmmodel_cw_year_xcos:
			order++;
		case constraint_dmmodel_cw_year_xsin:
			order++;
		case constraint_dmmodel_cw_year_cos:
			order++;
		case constraint_dmmodel_cw_year_sin:
			return consFunc_dmmodel_cw_year(psr,i,k,order);

			// IFUNC annual terms
		case constraint_ifunc_year_cos2:
			order++;
		case constraint_ifunc_year_sin2:
			order++;
		case constraint_ifunc_year_xcos:
			order++;
		case constraint_ifunc_year_xsin:
			order++;
		case constraint_ifunc_year_cos:
			order++;
		case constraint_ifunc_year_sin:
			return consFunc_ifunc_year(psr,i,k,order);

			// QIFUNC_p annual terms
		case constraint_qifunc_p_year_cos2:
			order++;
		case constraint_qifunc_p_year_sin2:
			order++;
		case constraint_qifunc_p_year_xcos:
			order++;
		case constraint_qifunc_p_year_xsin:
			order++;
		case constraint_qifunc_p_year_cos:
			order++;
		case constraint_qifunc_p_year_sin:
			return consFunc_qifunc_p_year(psr,i,k,order);

			// QIFUNC_c annual terms
		case constraint_qifunc_c_year_cos2:
			order++;
		case constraint_qifunc_c_year_sin2:
			order++;
		case constraint_qifunc_c_year_xcos:
			order++;
		case constraint_qifunc_c_year_xsin:
			order++;
		case constraint_qifunc_c_year_cos:
			order++;
		case constraint_qifunc_c_year_sin:
			return consFunc_qifunc_c_year(psr,i,k,order);

		default:
			return 0;
	}
}

#ifdef HAVE_LAPACK
#ifdef HAVE_BLAS



void othpl(int n,double x,double *pl){


	double a=2.0;
	double b=0.0;
	double c=1.0;
	double y0=1.0;
	double y1=2.0*x;
	pl[0]=1.0;
	pl[1]=2.0*x;



	for(int k=2;k<n;k++){

		double c=2.0*(k-1.0);
		double yn=(a*x+b)*y1-c*y0;
		pl[k]=yn;
		y0=y1;
		y1=yn;

	}



}




void getTempoNestMaxLike(pulsar *pulse, int npsr){

	int subDM=pulse->TNsubtractDM;
	int subRed=pulse->TNsubtractRed;

	pulse->TNsubtractDM=0;
	pulse->TNsubtractRed=0;


	formBatsAll(pulse,npsr);       /* Form Barycentric arrival times */
	formResiduals(pulse,npsr,1);       /* Form residuals */


	pulse->TNsubtractDM=subDM;
	pulse->TNsubtractRed=subRed;

	/////////////////////////////////////////////////////////////////////////////////////////////  
	/////////////////////////Form the Design Matrix////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////  

	int numtofit = 1;
	int numTime=0;
	int numJumps=0;
	for(int i=0; i<MAX_PARAMS; i++) {
		for(int k=0; k<pulse->param[i].aSize; k++) {
			if(pulse->param[i].fitFlag[k]==1) {
				if(i!=param_start && i!=param_finish) {
					numtofit++;
					numTime++;
				}
			} // if fitFlag
		} // for k
	} // for i

	/* Add extra parameters for jumps */
	for(int i=1; i<=pulse->nJumps; i++) {
		if(pulse->fitJump[i]==1){
			numtofit++;
			numJumps++;
		}
	}


	double **TNDM=new double*[pulse->nobs];
	for(int i=0;i<pulse->nobs;i++){
		TNDM[i]=new double[numtofit];
	}


	double pdParamDeriv[MAX_PARAMS];
	double *TNDMScale=new double[numtofit];
	for(int j=0; j<numtofit; j++) {
		TNDMScale[j]=0;
	}

	for(int i=0; i < pulse->nobs; i++) {
		FITfuncs(pulse[0].obsn[i].bat - pulse[0].param[param_pepoch].val[0], pdParamDeriv, numtofit, pulse, i,0);
		for(int j=0; j<numtofit; j++) {
			TNDM[i][j]=pdParamDeriv[j];
//			printf("PD: %i %i %g \n", i, j, pdParamDeriv[j]);
			TNDMScale[j]+=TNDM[i][j]*TNDM[i][j];
		} 
	} 


	for(int j=0; j<numtofit; j++) {
		TNDMScale[j]=sqrt(TNDMScale[j]);
	}


	for(int i=0; i < pulse->nobs; i++) {
		for(int j=0; j<numtofit; j++) {
			TNDM[i][j] = TNDM[i][j]/TNDMScale[j];
		}
	}


	int useOrthogonal=0;	

	/*


	   double* S = new double[numtofit];
	   double** U = new double*[pulse->nobs];
	   for(int k=0; k < pulse->nobs; k++){
	   U[k] = new double[pulse->nobs];
	   }
	   double** VT = new double*[numtofit]; 
	   for (int k=0; k<numtofit; k++) VT[k] = new double[numtofit];

	   dgesvd(TNDM,pulse->nobs, numtofit, S, U, VT);


	   double **V=new double*[numtofit];


	   for(int i=0;i<numtofit;i++){
	   V[i]=new double[numtofit];
	//	printf("DVD %i %g \n", i, S[i]);
	if(S[i] < 5*pow(10.0,-8)){
	useOrthogonal=1;
	printf("SVD Element %i below 5E-8: %g \n", i, S[i]);
	printf("Design matrix numerically unstable, using orthogonal representation\n");
	}
	}


	for(int j=0;j < numtofit;j++){
	for(int k=0;k < numtofit;k++){
	V[j][k]=VT[k][j];
	}
	}


	if(useOrthogonal==1){
	for(int j=0;j<pulse->nobs;j++){
	for(int k=0;k < numtofit;k++){
	TNDM[j][k]=U[j][k];
	}
	}
	}
	*/
	//////////////////////////////////////////////////////////////////////////////////////////  
	///////////////////////Noise Hyperparameters//////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////


	double start,end;
	int go=0;
	for (int i=0;i<pulse->nobs;i++)
	{
		if (pulse->obsn[i].deleted==0)
		{
			if (go==0)
			{
				go = 1;
				start = (double)pulse->obsn[i].bat;
				end  = start;
			}
			else
			{
				if (start > (double)pulse->obsn[i].bat)
					start = (double)pulse->obsn[i].bat;
				if (end < (double)pulse->obsn[i].bat)
					end = (double)pulse->obsn[i].bat;
			}
		}
	}

	double maxtspan=1*(end-start);
	double averageTSamp=2*maxtspan/pulse->nobs;

	double **DMEventInfo;


	double RedAmp=0;
	double RedIndex=0;
        double RedFLow=0;
        double RedCorner=0;
	double DMAmp=0;
	double DMIndex=0;
	int FitRedCoeff=0;
	int FitDMCoeff=0;
	int DMEventCoeff=0;
	int totCoeff=0;
	int DMEventQuadTerms=0;

	if(pulse->TNRedAmp != 0 && pulse->TNRedGam != 0){
		RedAmp=pulse->TNRedAmp;
		RedIndex=pulse->TNRedGam;
		FitRedCoeff=2*pulse->TNRedC;
                RedFLow=pow(10.0, pulse->TNRedFLow);
                RedCorner=pulse->TNRedCorner/maxtspan;

		for(int i = 0; i < 200; i++){
			pulse->TNRedCoeffs[i] = 0;
		}
	

		printf("\nIncluding Red noise with %i Frequencies, %g Log_10 Amplitude, %g Spectral Index\n", FitRedCoeff/2,RedAmp,RedIndex);
	}
	if(pulse->TNDMAmp != 0 && pulse->TNDMGam != 0){
		DMAmp=pulse->TNDMAmp;
		DMIndex=pulse->TNDMGam;
		FitDMCoeff=2*pulse->TNDMC;


                for(int i = 0; i < 200; i++){
                        pulse->TNDMCoeffs[i] = 0;
                }



		printf("\nIncluding DM Variations with %i Frequencies, %g Log_10 Amplitude, %g Spectral Index\n", FitDMCoeff/2,DMAmp,DMIndex);
	}

	if(pulse->nDMEvents > 0){
		DMEventInfo=new double*[pulse->nDMEvents];
		for(int i=0; i < pulse->nDMEvents; i++){


			printf("\nIncluding DM Event %i : %g Start, %g  Length, %g Log_10 Amp, %g Spectral Index\n", pulse->TNDMEvStart[i], pulse->TNDMEvLength[i], pulse->TNDMEvAmp[i], pulse->TNDMEvGam[i]);

			DMEventInfo[i]=new double[4];
			DMEventInfo[i][0]=pulse->TNDMEvStart[i]; //Start time
			DMEventInfo[i][1]=pulse->TNDMEvLength[i]; //Stop Time
			DMEventInfo[i][2]=pow(10.0, pulse->TNDMEvAmp[i]); //Amplitude
			DMEventInfo[i][3]=pulse->TNDMEvGam[i]; //SpectralIndex
			DMEventCoeff+=2*int(DMEventInfo[i][1]/averageTSamp);

			if(pulse->TNDMEvOff[i]==1){printf(" Including DMEvent offset \n"); DMEventQuadTerms++;}
			if(pulse->TNDMEvLin[i]==1){printf(" Including DMEvent Linear Term \n");DMEventQuadTerms++;}
			if(pulse->TNDMEvQuad[i]==1){printf(" Including DMEvent Quadratic Term \n");DMEventQuadTerms++;}

		}
	}


	// count ECORR values
	if(pulse->nTNECORR > 0){
		for(int i=0; i<pulse->nTNECORR; i++){
			printf("\nIncluding ECORR value for backend %s: %g mus", \
					pulse->TNECORRFlagVal[i], pulse->TNECORRVal[i]);

		}
	}


	// find number of epochs (default dt= 10 s)
	int *Processed = new int[pulse->nobs];

	// initialize processed flags
	for (int i=0;i<pulse->nobs;i++){
		Processed[i] = 1;
	}

	// make sure we only process the epochs with the chosen flags
	for (int i=0;i<pulse->nobs;i++){
		for (int j=0;j<pulse->obsn[i].nFlags;j++){
			for (int k=0;k<pulse->nTNECORR;k++){
				if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k])==0){
					if (strcmp(pulse->obsn[i].flagVal[j],pulse->TNECORRFlagVal[k])==0){
						Processed[i] = 0;
					}
				}
			}
		}
	}

	double dt = 10.0 / SECDAY;
	double satmin;
	double satmax;
	int nepoch = 0;
	int in = 0;
	int allProcessed = 0;
	while (!allProcessed){
		for (int i=0;i<pulse->nobs;i++){
			if (Processed[i]==0){
				satmin = (double)pulse->obsn[i].bat - dt;
				satmax = (double)pulse->obsn[i].bat + dt;
				break;
			}
		}
		for (int i=0;i<pulse->nobs;i++){
			for (int j=0;j<pulse->obsn[i].nFlags;j++){
				for (int k=0;k<pulse->nTNECORR;k++){
					if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k])==0){
						if (strcmp(pulse->obsn[i].flagVal[j],pulse->TNECORRFlagVal[k])==0){
							if ((double)pulse->obsn[i].bat > satmin && \
									(double)pulse->obsn[i].bat < satmax){
								Processed[i] = 1;
								in++;
							}
						}
					}
				}
			}
		}
		if (in != 0){
			nepoch++;
			in = 0;
		}
		allProcessed = 1;
		for (int i=0;i<pulse->nobs;i++){
			if (Processed[i]==0){
				allProcessed = 0;
				break;
			}
		}
	}


	if (nepoch > 0){
		printf("\n\nUsing %d epochs for PSR %s\n\n", nepoch, pulse->name);
	}
	// Total coefficients in F-matrix (include Jitter matrix here if present)
	totCoeff += FitRedCoeff;
	totCoeff += FitDMCoeff;
	totCoeff += DMEventCoeff;
	totCoeff += nepoch;

	if(pulse->TNBandDMAmp != 0 && pulse->TNBandDMGam != 0){
		printf("Including Band DM Noise: Amp %g   Index %g \n", pulse->TNBandDMAmp, pulse->TNBandDMGam);
		totCoeff += 6*pulse->TNBandDMC;
	}

        for(int i =0; i < pulse->nTNBandNoise; i++){
                printf("Including Band Noise between %g and %g MHz: Amp %g   Index %g \n", pulse->TNBandNoiseLF[i], pulse->TNBandNoiseHF[i], pulse->TNBandNoiseAmp[i], pulse->TNBandNoiseGam[i]);
                totCoeff += 2*pulse->TNBandNoiseC[i];
        }

	for(int i =0; i < pulse->nTNGroupNoise; i++){
		printf("Including Group Noise : Amp %g   Index %g \n", pulse->TNGroupNoiseAmp[i], pulse->TNGroupNoiseGam[i]);
		totCoeff += 2*pulse->TNGroupNoiseC[i];
	}

	printf("\nTotal number of coefficients: %d\n",totCoeff);

	double **FMatrix=new double*[pulse->nobs];
	for(int i=0;i<pulse->nobs;i++){
		FMatrix[i]=new double[totCoeff];
		for(int j=0;j<totCoeff;j++){
			FMatrix[i][j] = 0;
		}
	}


	double *freqs = new double[totCoeff];

	double *DMVec=new double[pulse->nobs];
	double DMKappa = 2.410*pow(10.0,-16);
	int startpos=0;
	double freqdet=0;

	double *powercoeff=new double[totCoeff];
	for(int o=0;o<totCoeff; o++){
		powercoeff[o]=0;
	}

	double Tspan = maxtspan;
	double f1yr = 1.0/3.16e7;


	//////////////////////////////////////////////////////////////////////////////////////////  
	///////////////////////Red Noise///////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	RedAmp=pow(10.0, RedAmp);

	for (int i=0; i<FitRedCoeff/2; i++){

		freqs[startpos+i]=RedFLow*((double)(i+1.0))/(maxtspan);
		freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

		double rho = pow((1+(pow((1.0/365.25)/RedCorner,RedIndex/2))),2)*(RedAmp*RedAmp/12.0/(M_PI*M_PI))/pow((1+(pow(freqs[i]/RedCorner,RedIndex/2))),2)/(maxtspan*24*60*60)*pow(f1yr,-3.0);
		powercoeff[i]+= rho;
		powercoeff[i+FitRedCoeff/2]+= rho;
		//printf("T2 RedC: %i %g %g %g %g %g \n", i, rho, RedAmp, RedIndex, freqs[i], maxtspan);
	}


	for(int i=0;i<FitRedCoeff/2;i++){
		for(int k=0;k<pulse->nobs;k++){
			double time=(double)pulse->obsn[k].bat;
			FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);

		}
	}

	for(int i=0;i<FitRedCoeff/2;i++){
		for(int k=0;k<pulse->nobs;k++){
			double time=(double)pulse->obsn[k].bat;
			FMatrix[k][i+FitRedCoeff/2]=sin(2*M_PI*freqs[i]*time);
		}
	}


	//////////////////////////////////////////////////////////////////////////////////////////  
	///////////////////////DM Variations//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////


	startpos=FitRedCoeff;
	DMAmp=pow(10.0, DMAmp);

	for (int i=0; i<FitDMCoeff/2; i++){

		freqs[startpos+i]=((double)(i+1.0))/maxtspan;
		freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];

		double rho = (DMAmp*DMAmp)*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-DMIndex))/(maxtspan*24*60*60);	
		powercoeff[startpos+i]+=rho;
		powercoeff[startpos+i+FitDMCoeff/2]+=rho;
	}



	for(int o=0;o<pulse->nobs; o++){
		DMVec[o]=1.0/(DMKappa*pow((double)pulse->obsn[o].freqSSB,2));
	}

	for(int i=0;i<FitDMCoeff/2;i++){
		for(int k=0;k<pulse->nobs;k++){
			double time=(double)pulse->obsn[k].bat;
			FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
		}
	}

	for(int i=0;i<FitDMCoeff/2;i++){
		for(int k=0;k<pulse->nobs;k++){
			double time=(double)pulse->obsn[k].bat;
			FMatrix[k][startpos+i+FitDMCoeff/2]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
		}
	}



	//////////////////////////////////////////////////////////////////////////////////////////  
	///////////////////////DM Events//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////	

	startpos+=FitDMCoeff;

	if(pulse->nDMEvents > 0){

		for(int i =0; i < pulse->nDMEvents; i++){

			double DMamp=DMEventInfo[i][2];
			double DMindex=DMEventInfo[i][3];

			double Tspan = DMEventInfo[i][1];
			double f1yr = 1.0/3.16e7;
			int DMEventnumCoeff=int(Tspan/averageTSamp);

			for (int c=0; c<DMEventnumCoeff; c++){
				freqs[startpos+c]=((double)(c+1))/Tspan;
				freqs[startpos+c+DMEventnumCoeff]=freqs[startpos+c];

				double rho = (DMamp*DMamp)*pow(f1yr,(-3)) * pow(freqs[startpos+c]*365.25,(-DMindex))/(maxtspan*24*60*60);
				powercoeff[startpos+c]+=rho;
				powercoeff[startpos+c+DMEventnumCoeff]+=rho;
				freqdet=freqdet+2*log(powercoeff[startpos+c]);
			}


			for(int c=0;c<DMEventnumCoeff;c++){
				for(int k=0;k<pulse->nobs;k++){
					double time=(double)pulse->obsn[k].bat;
					if(time < DMEventInfo[i][0]+Tspan && time > DMEventInfo[i][0]){
						FMatrix[k][startpos+c]=cos(2*M_PI*freqs[startpos+c]*time)*DMVec[k];
						FMatrix[k][startpos+c+DMEventnumCoeff]=sin(2*M_PI*freqs[startpos+c]*time)*DMVec[k];
					}
					else{
						FMatrix[k][startpos+c]=0;
						FMatrix[k][startpos+c+DMEventnumCoeff]=0;
					}
				}
			}

			startpos+=2*DMEventnumCoeff;

		}
	}

	// shapelet events
	int ShapeEventTerms=0;
	if(pulse->nTNShapeletEvents > 0){
		for(int i =0; i < pulse->nTNShapeletEvents; i++){
			ShapeEventTerms += pulse->TNShapeletEvN[i];
		}	

	}
	//printf("Total event coeffs: %i \n", ShapeEventTerms);


	startpos += DMEventCoeff;

	//////////////////////////////////////////////////////////////////////////////////////////  
	///////////////////////DM Band Noise//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	if(pulse->TNBandDMAmp != 0 && pulse->TNBandDMGam != 0){

		double BandDMAmp=pow(10.0, pulse->TNBandDMAmp);



		double startfreq=0;
		double stopfreq=0;

		/////////////////////////////////50CM/////////////////////////////////////////////////////////////////

		startfreq = 0;
		stopfreq=1000;

		for (int i=0; i<pulse->TNBandDMC; i++){

			freqs[startpos+i]=((double)(i+1.0))/maxtspan;
			freqs[startpos+i+pulse->TNBandDMC]=freqs[startpos+i];
			
			double rho = (BandDMAmp*BandDMAmp)*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-pulse->TNBandDMGam))/(maxtspan*24*60*60);	
			powercoeff[startpos+i]+=rho;
			powercoeff[startpos+i+pulse->TNBandDMC]+=rho;
		}
		
		


		for(int i=0;i<pulse->TNBandDMC;i++){
		        for(int k=0;k<pulse->nobs;k++){
				if(pulse->obsn[k].freq > startfreq && pulse->obsn[k].freq < stopfreq){
		               		double time=(double)pulse->obsn[k].bat;
		                	FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
				}
				else{	
					FMatrix[k][startpos+i]=0;
				}
		        }
		}

		for(int i=0;i<pulse->TNBandDMC;i++){
		        for(int k=0;k<pulse->nobs;k++){
				if(pulse->obsn[k].freq > startfreq && pulse->obsn[k].freq < stopfreq){		
				        double time=(double)pulse->obsn[k].bat;
				        FMatrix[k][startpos+i+pulse->TNBandDMC]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
				}
				else{
					FMatrix[k][startpos+i+pulse->TNBandDMC]=0;
				}
		        }
		}


		startpos=startpos+2*pulse->TNBandDMC;

		/////////////////////////////////20CM/////////////////////////////////////////////////////////////////

		startfreq = 1000;
		stopfreq=1800;

		for (int i=0; i<pulse->TNBandDMC; i++){

			freqs[startpos+i]=((double)(i+1.0))/maxtspan;
			freqs[startpos+i+pulse->TNBandDMC]=freqs[startpos+i];
			
			double rho = (BandDMAmp*BandDMAmp)*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-pulse->TNBandDMGam))/(maxtspan*24*60*60);	
			powercoeff[startpos+i]+=rho;
			powercoeff[startpos+i+pulse->TNBandDMC]+=rho;
		}
		
		


		for(int i=0;i<pulse->TNBandDMC;i++){
		        for(int k=0;k<pulse->nobs;k++){
				if(pulse->obsn[k].freq > startfreq && pulse->obsn[k].freq < stopfreq){
		               		double time=(double)pulse->obsn[k].bat;
		                	FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
				}
				else{	
					FMatrix[k][startpos+i]=0;
				}
		        }
		}

		for(int i=0;i<pulse->TNBandDMC;i++){
		        for(int k=0;k<pulse->nobs;k++){
				if(pulse->obsn[k].freq > startfreq && pulse->obsn[k].freq < stopfreq){		
				        double time=(double)pulse->obsn[k].bat;
				        FMatrix[k][startpos+i+pulse->TNBandDMC]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
				}
				else{
					FMatrix[k][startpos+i+pulse->TNBandDMC]=0;
				}
		        }
		}


		startpos=startpos+2*pulse->TNBandDMC;

		/////////////////////////////////10CM/////////////////////////////////////////////////////////////////


		startfreq = 1800;
		stopfreq=10000;

		for (int i=0; i<pulse->TNBandDMC; i++){

			freqs[startpos+i]=((double)(i+1.0))/maxtspan;
			freqs[startpos+i+pulse->TNBandDMC]=freqs[startpos+i];
			
			double rho = (BandDMAmp*BandDMAmp)*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-pulse->TNBandDMGam))/(maxtspan*24*60*60);	
			powercoeff[startpos+i]+=rho;
			powercoeff[startpos+i+pulse->TNBandDMC]+=rho;
		}
		
		


		for(int i=0;i<pulse->TNBandDMC;i++){
		        for(int k=0;k<pulse->nobs;k++){
				if(pulse->obsn[k].freq > startfreq && pulse->obsn[k].freq < stopfreq){
		               		double time=(double)pulse->obsn[k].bat;
		                	FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
				}
				else{	
					FMatrix[k][startpos+i]=0;
				}
		        }
		}

		for(int i=0;i<pulse->TNBandDMC;i++){
		        for(int k=0;k<pulse->nobs;k++){
				if(pulse->obsn[k].freq > startfreq && pulse->obsn[k].freq < stopfreq){		
				        double time=(double)pulse->obsn[k].bat;
				        FMatrix[k][startpos+i+pulse->TNBandDMC]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
				}
				else{
					FMatrix[k][startpos+i+pulse->TNBandDMC]=0;
				}
		        }
		}


		startpos=startpos+2*pulse->TNBandDMC;
	}


	//////////////////////////////////////////////////////////////////////////////////////////  
	///////////////////////Band Noise//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	int totalBandNoiseCoeff=0;
	for(int g =0; g < pulse->nTNBandNoise; g++){

		double BandLF = pulse->TNBandNoiseLF[g];
		double BandHF = pulse->TNBandNoiseHF[g];
		double BandAmp=pow(10.0, pulse->TNBandNoiseAmp[g]);
		double BandSpec=pulse->TNBandNoiseGam[g];
		int BandC=pulse->TNBandNoiseC[g];

		totalBandNoiseCoeff+=2*BandC;



		for (int i=0; i<BandC; i++){

			freqs[startpos+i]=((double)(i+1.0))/maxtspan;
			freqs[startpos+i+BandC]=freqs[startpos+i];
			
			double rho = (BandAmp*BandAmp/12.0/(M_PI*M_PI))*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-BandSpec))/(maxtspan*24*60*60);	
			powercoeff[startpos+i]+=rho;
			powercoeff[startpos+i+BandC]+=rho;
		}
		
		


		for(int i=0;i<BandC;i++){
		        for(int k=0; k<pulse->nobs; k++){
				if(pulse->obsn[k].freq > BandLF && pulse->obsn[k].freq < BandHF){
					double time=(double)pulse->obsn[k].bat;
					FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time);
					FMatrix[k][startpos+i+BandC]=sin(2*M_PI*freqs[startpos+i]*time);

				}
				else{
						FMatrix[k][startpos+i]=0;
						FMatrix[k][startpos+i+BandC]=0;
				}
					
			}
		}




		startpos=startpos+2*BandC;
	}


	//////////////////////////////////////////////////////////////////////////////////////////  
	///////////////////////Group Noise//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	int totalGroupCoeff=0;
	for(int g =0; g < pulse->nTNGroupNoise; g++){

		double GroupAmp=pow(10.0, pulse->TNGroupNoiseAmp[g]);
		double GroupSpec=pulse->TNGroupNoiseGam[g];
		int GroupC=pulse->TNGroupNoiseC[g];

		totalGroupCoeff+=2*GroupC;



		for (int i=0; i<GroupC; i++){

			freqs[startpos+i]=((double)(i+1.0))/maxtspan;
			freqs[startpos+i+GroupC]=freqs[startpos+i];
			
			double rho = (GroupAmp*GroupAmp/12.0/(M_PI*M_PI))*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-GroupSpec))/(maxtspan*24*60*60);	
			powercoeff[startpos+i]+=rho;
			powercoeff[startpos+i+GroupC]+=rho;
		}
		
		


		for(int i=0;i<GroupC;i++){
		        for(int k=0; k<pulse->nobs; k++){
				 int set=0;
				 for (int j=0; j < pulse->obsn[k].nFlags; j++){
                
					   //Check Group Noise Flag
	
			         	 if (strcmp(pulse->obsn[k].flagID[j],pulse->TNGroupNoiseFlagID[g])==0){
			    	          	if (strcmp(pulse->obsn[k].flagVal[j],pulse->TNGroupNoiseFlagVal[g])==0){
					       		double time=(double)pulse->obsn[k].bat;
							FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time);
							FMatrix[k][startpos+i+GroupC]=sin(2*M_PI*freqs[startpos+i]*time);
							set = 1;

						}
						else{
							if(set == 0){	
								FMatrix[k][startpos+i]=0;
								FMatrix[k][startpos+i+GroupC]=0;
							}
						}
					}
					else{
						if(set == 0){
							FMatrix[k][startpos+i]=0;
							FMatrix[k][startpos+i+GroupC]=0;
						}
					}	
					
				}
			}
		}




		startpos=startpos+2*GroupC;
	}



	//////////////////////////////////////////////////////////////////////////////////////////  
	///////////////////////U Jitter Matrix////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////


	//// construct jitter matrix //

	
	//printf("\nStarting Coefficients at %d for Jitter:\n",startpos);

	// initialize processed flags
	for (int i=0;i<pulse->nobs;i++){
		Processed[i] = 1;
	}

	// make sure we only process the epochs with the chosen flags
	for (int i=0;i<pulse->nobs;i++){
		for (int j=0;j<pulse->obsn[i].nFlags;j++){
			for (int k=0;k<pulse->nTNECORR;k++){
				if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k])==0){
					if (strcmp(pulse->obsn[i].flagVal[j],pulse->TNECORRFlagVal[k])==0){
						Processed[i] = 0;
					}
				}
			}
		}
	}


	nepoch = 0;
	in = 0;
	allProcessed = 0;
	while (!allProcessed){
		for (int i=0;i<pulse->nobs;i++){
			if (Processed[i]==0){
				satmin = (double)pulse->obsn[i].bat - dt;
				satmax = (double)pulse->obsn[i].bat + dt;
				break;
			}
		}
		for (int i=0;i<pulse->nobs;i++){
			for (int j=0;j<pulse->obsn[i].nFlags;j++){
				for (int k=0;k<pulse->nTNECORR;k++){
					if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k])==0){
						if (strcmp(pulse->obsn[i].flagVal[j],pulse->TNECORRFlagVal[k])==0){
							if ((double)pulse->obsn[i].bat > satmin && \
									(double)pulse->obsn[i].bat < satmax){
								Processed[i] = 1;
								FMatrix[i][startpos+nepoch] = 1.0;
								powercoeff[startpos+nepoch] = pow(pulse->TNECORRVal[k], \
										2.0)*1e-12;
								in++;
							}
							else{
								FMatrix[i][startpos+nepoch] = 0.0;
							}
						}
					}
				}
			}
		}
		if (in != 0){
			nepoch++;
			in = 0;
		}
		allProcessed = 1;
		for (int i=0;i<pulse->nobs;i++){
			if (Processed[i]==0){
				allProcessed = 0;
				break;
			}
		}
	}


	/*for (int i;i<totCoeff;i++){
	  printf("powercoeff %d: %g\n",i,powercoeff[i]);
	  }*/

	// get Phi matrix
	double **PPFM=new double*[totCoeff];
	for(int i=0;i<totCoeff;i++){
		PPFM[i]=new double[totCoeff];
		for(int j=0;j<totCoeff;j++){
			PPFM[i][j]=0;
		}
	}


	for(int c1=0; c1<totCoeff; c1++){
//		printf("PPFM %i %g \n", c1, powercoeff[c1]);
		PPFM[c1][c1]=1.0/powercoeff[c1];
	}

	//////////////////////////////////////////////////////////////////////////////////////////  
	//////////////////////Form Total Matrix///////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	int totalsize=numtofit+totCoeff+DMEventQuadTerms+ShapeEventTerms;
	int totalnoisesize=numtofit+FitRedCoeff+FitDMCoeff+totalBandNoiseCoeff+6*pulse->TNBandDMC+totalGroupCoeff;
	int totalDMsize = DMEventQuadTerms+ShapeEventTerms+FitDMCoeff;
	
	double **TotalMatrix=new double*[pulse->nobs];
	for(int i =0;i<pulse->nobs;i++){
		TotalMatrix[i]=new double[totalsize];
		for(int j =0;j<totalsize; j++){
			TotalMatrix[i][j]=0;
		}
	}


	for(int i =0;i<pulse->nobs;i++){
		int startpoint=0;
		double time=((double)pulse->obsn[i].bat);

		for(int j =0;j<numtofit; j++){
			TotalMatrix[i][j]=TNDM[i][j];
		}

		startpoint+=numtofit;		
		for(int j =0;j<totCoeff; j++){
			TotalMatrix[i][j+startpoint]=FMatrix[i][j];
		}

		startpoint += totCoeff;
		int DMEvterms=0;
		for(int e =0; e < pulse->nDMEvents; e++){
			if(time < DMEventInfo[e][0]+DMEventInfo[e][1] && time > DMEventInfo[e][0]){
				if(pulse->TNDMEvOff[e]==1){	
					TotalMatrix[i][startpoint+DMEvterms]=DMVec[i];
					DMEvterms++;
				}
				if(pulse->TNDMEvLin[e]==1){	
					TotalMatrix[i][startpoint+DMEvterms]=(time-DMEventInfo[e][0])*DMVec[i];
					DMEvterms++;
				}
				if(pulse->TNDMEvQuad[e]==1){	
					TotalMatrix[i][startpoint+DMEvterms]=pow((time-DMEventInfo[e][0]),2)*DMVec[i];
					DMEvterms++;
				}
			}
		}

		startpoint += DMEventQuadTerms;
		int ShapeEvterms=0;
		for(int e =0; e < pulse->nTNShapeletEvents; e++){

			double *shapeVec =  new double[pulse->TNShapeletEvN[e]];
			double HVal=(time-pulse->TNShapeletEvPos[e])/(sqrt(2.0)*pulse->TNShapeletEvWidth[e]);
			othpl(pulse->TNShapeletEvN[e],HVal,shapeVec);

			for(int s=0; s < pulse->TNShapeletEvN[e]; s++){

				double NormTerm=1.0/sqrt(sqrt(2.0*M_PI)*pow(2.0,s));
				TotalMatrix[i][startpoint+ShapeEvterms] = NormTerm*shapeVec[s]*exp(-0.5*pow((time-pulse->TNShapeletEvPos[e])/pulse->TNShapeletEvWidth[e], 2))*pow(DMVec[i], pulse->TNShapeletEvFScale[e]/2.0);

				ShapeEvterms++;

			}
			delete[] shapeVec;
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////  
	//////////////////////Get Residuals Vector////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	double *Resvec=new double[pulse->nobs];
	for(int o=0;o<pulse->nobs; o++){
		Resvec[o]=(double)pulse->obsn[o].residual;
		pulse->obsn[o].prefitResidual = pulse->obsn[o].residual;
	}




	//////////////////////////////////////////////////////////////////////////////////////////  
	//////////////////////Get White Noise Vector//////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	double *Noise=new double[pulse->nobs];
	for(int o=0;o<pulse->nobs; o++){
		Noise[o]=pow(((pulse->obsn[o].toaErr)*pow(10.0,-6)),2);
	}


	//////////////////////////////////////////////////////////////////////////////////////////  
	//////////////////////Do Algebra//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	double **NG = new double*[pulse->nobs]; for (int k=0; k<pulse->nobs; k++) NG[k] = new double[totalsize];
	double **GNG = new double*[totalsize]; for (int k=0; k<totalsize; k++) GNG[k] = new double[totalsize];



	for(int i=0;i<pulse->nobs;i++){
		for(int j=0;j<totalsize; j++){

			NG[i][j]=TotalMatrix[i][j]/Noise[i];

		}
	}


	dgemm(TotalMatrix, NG,GNG,pulse->nobs, totalsize,pulse->nobs, totalsize, 'T','N');


	for(int j =0;j<totCoeff; j++){
		GNG[numtofit+j][numtofit+j]+=PPFM[j][j];
	}


	double tdet=0;
	dpotrf(GNG, totalsize, tdet);
	dpotri(GNG,totalsize);
//	printf("logdet GNG = %g\n",tdet);

	double *dG=new double[totalsize];
	dgemv(NG,Resvec,dG,pulse->nobs,totalsize,'T');

	double *maxcoeff=new double[totalsize];
	dgemv(GNG,dG,maxcoeff,totalsize,totalsize,'N');

	long double *Errorvec=new long double[totalsize];

	for(int i =0; i < totalsize; i++){
//		printf("GNG %i %g \n", i, GNG[i][i]);
		Errorvec[i]=pow(GNG[i][i], 0.5);
	}

	double *Scoeff=new double[numtofit];
	long double *Serr=new long double[numtofit];

	/*	for(int i =0; i < numtofit; i++){
		if(S[i] >= 0){
		Scoeff[i]=maxcoeff[i]/S[i];
		Serr[i]=Errorvec[i]/S[i];
		}
		else{
		Scoeff[i]=0;
		Serr[i]=0;
		}
		}
		*/
	double *TempoCoeff = new double[numtofit];
	double *TempoErr =  new double[numtofit];
	//	dgemv(V,Scoeff,TempoCoeff,numtofit,numtofit, 'N');

	for(int i=0;i<numtofit; i++){

		//		long double errsum=0;
		//          	for(int j=0;j<numtofit; j++){
		//			errsum += pow(V[i][j]*Serr[j],2);
		//               }
		//          	TempoErr[i]=pow(errsum,0.5);///TNDMScale[i];
		//		TempoCoeff[i]=TempoCoeff[i];///TNDMScale[i];

		if(useOrthogonal==0){
			TempoCoeff[i]=maxcoeff[i]/TNDMScale[i];
			TempoErr[i]=((double)Errorvec[i])/TNDMScale[i];
		}
	}
	updateParameters(pulse,0,TempoCoeff,TempoErr);


	//////////////////Get Red noise and DM Coeffs and errors/////////////////////
	
	int redcounter = 0;
	int dmcounter = 0;
	for(int j=0;j<totalsize; j++){
		if(j>=numtofit && j < numtofit+FitRedCoeff){
			//pulse->TNRedCoeffs[redcounter] = maxcoeff[j];
			//pulse->TNRedCoeffs[redcounter+100] = Errorvec[j];
			//printf("Red Coeff: %i %g +/- %g \n", redcounter, pulse->TNRedCoeffs[redcounter], pulse->TNRedCoeffs[redcounter+100]);
			redcounter++;
		}

		if(j>=FitRedCoeff+numtofit && j < numtofit+FitRedCoeff+FitDMCoeff){
		        //pulse->TNDMCoeffs[dmcounter] = maxcoeff[j];
                        //pulse->TNDMCoeffs[dmcounter+100] = Errorvec[j];
                        //printf("DM Coeff: %i %g +/- %g \n", dmcounter, pulse->TNDMCoeffs[dmcounter], pulse->TNDMCoeffs[dmcounter+100]);
			dmcounter++;
                }
	}

	double chisq=0;
	printf("subtractRed is %i \n", pulse->TNsubtractRed);
	printf("subtractDM is %i \n", pulse->TNsubtractDM);

	for(int i=0;i<pulse->nobs;i++){
		double dsum=0;
		double redsum=0;
		double rederr=0;
		double dmsum=0;
		double dmerr=0;
		double shapesum=0;
		double shapeerr=0;
		for(int j=0;j<totalsize; j++){
			//			if(i==0)printf("Max coeff: %i %g \n", j,maxcoeff[j]);		
			dsum=dsum+TotalMatrix[i][j]*maxcoeff[j];

			if(j>=numtofit && j < numtofit+FitRedCoeff){
				//if(i==20){
				//printf("TM: %i %g %g \n", j-numtofit, TotalMatrix[i][j], maxcoeff[j]);
				//}
				redsum+=TotalMatrix[i][j]*maxcoeff[j];
				rederr+=pow(TotalMatrix[i][j]*Errorvec[j],2);
			}
			if(j>=FitRedCoeff+numtofit && j < numtofit+FitRedCoeff+FitDMCoeff){
				//if(i==20){
                                //printf("TM: %i %g %g \n", j-numtofit, TotalMatrix[i][j], maxcoeff[j]);
                                //}

				dmsum+=TotalMatrix[i][j]*maxcoeff[j];
				dmerr+=pow(TotalMatrix[i][j]*Errorvec[j],2);
			}

			if(j>=numtofit+FitRedCoeff+FitDMCoeff && j < totalnoisesize){
					redsum+=TotalMatrix[i][j]*maxcoeff[j];
					rederr+=pow(TotalMatrix[i][j]*Errorvec[j],2);
			}		
                        if(j>=numtofit+totCoeff && j < numtofit+totCoeff+ShapeEventTerms){
//					printf("Shapeevent terms %i %i %g %g \n", i, j, TotalMatrix[i][j]*maxcoeff[j], pow(TotalMatrix[i][j]*Errorvec[j],2));
                                        dmsum+=TotalMatrix[i][j]*maxcoeff[j];
                                        dmerr+=pow(TotalMatrix[i][j]*Errorvec[j],2);
					shapesum+=TotalMatrix[i][j]*maxcoeff[j];
					shapeerr+=pow(TotalMatrix[i][j]*Errorvec[j],2);
                        }		
			

		}

//		if(fabs(shapesum) > pow(10.0, -10)){printf("Shapeevent terms %i %.10g %g %g \n", i, (double)pulse->obsn[i].bat, shapesum/DMVec[i], sqrt(shapeerr)/DMVec[i]);}


		double freq=(double)pulse->obsn[i].freqSSB;
		long double yrs = (pulse->obsn[i].bat - pulse->param[param_dmepoch].val[0])/365.25;
		long double arg = 1.0;
		double dmDot=0;
		double dmDotErr=0;
		for (int d=0;d<9;d++){
			if(d>0){
				arg *= yrs;
			}
			if (pulse->param[param_dm].paramSet[d]==1){
				if(d>0){
					dmDot+=(double)(pulse->param[param_dm].val[d]*arg);
				}
				dmDotErr+=pow((double)(pulse->param[param_dm].err[d]*arg),2);
			}
		}

		double pDotErr=0;
		pDotErr+=pow(TempoErr[0],2);
		if(pulse->param[param_f].paramSet[0]==1){
			arg=((pulse->obsn[i].bat - pulse->param[param_pepoch].val[0])/ \
					pulse[0].param[param_f].val[0])*86400.0;
		        pDotErr+=pow((double)(pulse->param[param_f].err[0]*arg),2);

		}

		if(pulse->param[param_f].paramSet[1]==1){
			arg=0.5*pow((double)(pulse->obsn[i].bat - pulse->param[param_pepoch].val[0]), 2);
			long double argerr = (pulse->param[param_f].err[1]/ \
					pulse[0].param[param_f].val[0])*86400.0*86400;
			pDotErr+=pow((double)(argerr*arg),2);
		}


		chisq+=(Resvec[i]-dsum)*(Resvec[i]-dsum)/(Noise[i]);
		pulse->obsn[i].TNRedSignal=redsum;
		pulse->obsn[i].TNRedErr=pow(rederr+pDotErr,0.5);

		pulse->obsn[i].TNDMSignal=dmsum;
		pulse->obsn[i].TNDMErr=pow(dmerr/pow(DMVec[i],2) + dmDotErr,0.5);
	}

	//////////////////////////////////////////////////////////////////////////////////////////  
	//////////////////////Calculate GLS Chi-squared///////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	for (int i=0; i<pulse->nobs;i++){
		for (int j=0; j<numtofit;j++){
			Resvec[i] -= TotalMatrix[i][j]*maxcoeff[j];
		}
	}

	int nother = totalsize-numtofit;
	double **TMat = new double*[pulse->nobs];
	for (int k=0;k<pulse->nobs;k++){
		TMat[k] = new double[nother];
	}

	for(int i=0;i<pulse->nobs;i++){
		for(int j=0;j<nother; j++){
			TMat[i][j] = TotalMatrix[i][j+numtofit];
		}
	}

	double **NT = new double*[pulse->nobs];
	for (int k=0;k<pulse->nobs;k++){
		NT[k] = new double[nother];
	}

	double **TNT = new double*[nother]; 
	for (int k=0; k<nother; k++){
		TNT[k] = new double[nother];
	}

	for(int i=0;i<pulse->nobs;i++){
		for(int j=0;j<nother; j++){
			NT[i][j] = TMat[i][j]/Noise[i];
		}
	}

	dgemm(TMat, NT, TNT, pulse->nobs, nother, pulse->nobs, nother, 'T', 'N');

	for(int j =0;j<totCoeff; j++){
		TNT[j][j] += PPFM[j][j];
	}

	tdet = 0.0;
	dpotrf(TNT, nother, tdet);
	dpotri(TNT, nother);

	double *dT = new double[nother];
	dgemv(NT, Resvec, dT, pulse->nobs, nother, 'T');

	double *Sigmad = new double[nother];
	dgemv(TNT, dT, Sigmad, nother, nother, 'N');

	chisq = 0.0;
	double timesq=0;
	for (int i=0; i<pulse->nobs; i++){
		timesq+=Resvec[i]*Resvec[i]/Noise[i];
		chisq += Resvec[i]*Resvec[i]/Noise[i];
	}
		

	double freqsq=0;
	for (int k=0; k<nother; k++){
		freqsq+=dT[k]*Sigmad[k];
		chisq -= dT[k]*Sigmad[k];
	}

	pulse->fitChisq = chisq; 
	pulse->fitNfree = pulse->nobs-numtofit;
	pulse->nFit=pulse->nobs;
	pulse->nParam = numtofit;


	/*	for (int j = 0; j < pulse->nobs; j++){
		delete[]U[j];
		}
		delete[]U;

		delete[]S;
		delete[]Scoeff;	
		delete[]Serr;	

		for (int j = 0; j < numtofit; j++){
		delete[]VT[j];
		delete[]V[j];
		}

		delete[]V; 
		delete[]VT;*/
	delete[] TempoErr;
	delete[] TempoCoeff;
	delete[] DMVec;
	delete[] dT;
	delete[] Sigmad;
	for(int i=0;i<pulse->nobs;i++){
		delete[] TNDM[i];
	}
	delete[] TNDM;
	for(int i=0; i<pulse->nobs; i++){
		delete[] TMat[i];
		delete[] NT[i];
	}
	delete[] TMat;
	delete[] NT;
	for(int i=0; i<nother; i++){
		delete[] TNT[i];
	}
	delete[] TNT;

	delete[] dG;
	delete[] maxcoeff;
	delete[] Errorvec;
	delete[] Resvec;
	delete[] Noise;

	for (int k=0; k<pulse->nobs; k++){
		delete[] NG[k];
	}
	delete[] NG;
	for (int k=0; k<totalsize; k++){
		delete[] GNG[k];
	}
	delete[] GNG;

	for(int i =0;i<pulse->nobs;i++){
		delete[] TotalMatrix[i];
	}
	delete[] TotalMatrix;

	for(int i=0;i<totCoeff;i++){
		delete[] PPFM[i];
	}
	delete[] PPFM;

	delete[] powercoeff;
	delete[] freqs;
	for(int i=0;i<pulse->nobs;i++){
		delete[] FMatrix[i];
	}
	delete[] FMatrix;
	delete[]Processed;

}


void dgesvd(double **A, int m, int n, double *S, double **U, double **VT)
{
	char jobu, jobvt;
	int lda, ldu, ldvt, lwork, info;
	double *a, *u, *vt, *work;

	int minmn, maxmn;

	jobu = 'A'; /* Specifies options for computing U.
A: all M columns of U are returned in array U;
S: the first min(m,n) columns of U (the left
singular vectors) are returned in the array U;
O: the first min(m,n) columns of U (the left
singular vectors) are overwritten on the array A;
N: no columns of U (no left singular vectors) are
computed. */

	jobvt = 'A'; /* Specifies options for computing VT.
A: all N rows of V**T are returned in the array
VT;
S: the first min(m,n) rows of V**T (the right
singular vectors) are returned in the array VT;
O: the first min(m,n) rows of V**T (the right
singular vectors) are overwritten on the array A;
N: no rows of V**T (no right singular vectors) are
computed. */

	lda = m; // The leading dimension of the matrix a.
	a = dgesvd_ctof(A, lda, n); /* Convert the matrix A from double pointer
								   C form to single pointer Fortran form. */

	ldu = m;


	maxmn = m;
	minmn = n;

	ldu = m; // Left singular vector matrix
	u = new double[ldu*ldu];

	ldvt = n; // Right singular vector matrix
	vt = new double[ldvt*n];

	int LMAX=100000;

	work = new double[LMAX];
	lwork = -1; // Set up the work array, larger than needed.

	// 	printf("parm 11 %i %i\n",ldu,ldvt);
	dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, S, u,&ldu, vt, &ldvt, work, &lwork, &info);

	lwork = std::min(LMAX,int(work[0]));

	dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, S, u,&ldu, vt, &ldvt, work, &lwork, &info);
	// 	printf("parm 11 out %i %i\n",ldu,ldvt);
	dgesvd_ftoc(u, U, ldu, ldu);
	dgesvd_ftoc(vt, VT, ldvt, n);

	delete a;
	delete u;
	delete vt;
	delete work;
}


double* dgesvd_ctof(double **in, int rows, int cols)
{
	double *out;
	int i, j;

	out = new double[rows*cols];
	for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
	return(out);
}


void dgesvd_ftoc(double *in, double **out, int rows, int cols)
{
	int i, j;

	for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}


void dgemv(double **A, double *vecin,double *vecout,int rowa, int cola, char AT)
{

	int M,N,K;
	double *a;

	double alpha=1;
	double beta=0;
	int incX=1;
	int incY=1;

	a = dgemv_ctof(A, rowa, cola); 

	dgemv_(&AT, &rowa, &cola, &alpha, a, &rowa, vecin, &incX, &beta, vecout, &incY);

	delete a;
}


double* dgemv_ctof(double **in, int rows, int cols)
{
	double *out;
	int i, j;

	out = new double[rows*cols];
	for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
	return(out);
}


void dgemv_ftoc(double *in, double **out, int rows, int cols)
{
	int i, j;

	for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}

void dgemm(double **A, double **B,double **C,int rowa, int cola, int rowb, int colb, char AT, char BT)
{

	int M,N,K;
	double *a, *b, *c;

	double alpha=1;
	double beta=0;


	if(AT == 'N'){
		M=rowa;
		K=cola;
	}
	else if(AT == 'T'){
		M=cola;
		K=rowa;
	}

	if(BT == 'N'){
		N=colb;
	}
	else if(BT == 'T'){
		N=rowa;
	}	
	/*
	   (TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)*/


	a = dgemm_ctof(A, rowa, cola); 
	b = dgemm_ctof(B, rowb, colb);
	c = new double[M*N];

	dgemm_(&AT, &BT, &M, &N, &K, &alpha, a, &rowa,b, &rowb, &beta, c, &M);
	dgemm_ftoc(c, C, M, N);


	delete a;
	delete b;
	delete c;

}


double* dgemm_ctof(double **in, int rows, int cols)
{
	double *out;
	int i, j;

	out = new double[rows*cols];
	for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
	return(out);
}


void dgemm_ftoc(double *in, double **out, int rows, int cols)
{
	int i, j;

	for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}

void dpotri(double **A, int msize)
{

	int info;
	double *a;
	char UPLO='L';

	a = dpotri_ctof(A, msize, msize); 

	dpotri_(&UPLO, &msize, a, &msize, &info);
	dpotri_ftoc(a, A, msize, msize);

	for(int i=0;i<msize;i++){
		for(int j=0;j<i;j++){
			A[j][i]=A[i][j];
		}
	}


	delete a;

}


double* dpotri_ctof(double **in, int rows, int cols)
{
	double *out;
	int i, j;

	out = new double[rows*cols];
	for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
	return(out);
}


void dpotri_ftoc(double *in, double **out, int rows, int cols)
{
	int i, j;

	for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}

void dpotrf(double **A, int msize, double &det)
{

	int info;
	double *a;
	char UPLO='L';

	a = dpotrf_ctof(A, msize, msize); 

	dpotrf_(&UPLO, &msize, a, &msize, &info);
	dpotrf_ftoc(a, A, msize, msize);

	det=0;
	for(int i=0;i<msize;i++){
		det+=log(A[i][i]);
	}

	det=det*2;

	//printf("info: %i \n", info);
	delete a;

}


double* dpotrf_ctof(double **in, int rows, int cols)
{
	double *out;
	int i, j;

	out = new double[rows*cols];
	for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
	return(out);
}


void dpotrf_ftoc(double *in, double **out, int rows, int cols)
{
	int i, j;

	for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}

#endif
#endif
