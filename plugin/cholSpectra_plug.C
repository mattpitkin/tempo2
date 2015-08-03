
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

/**
 * 2012-11-01 (M.Keith) Updated to use standard libtempo routines
 *                      rather than custom 'plugin' versions.
 *                      No longer requires fftw or pgplot.
 */

/* Not sure who's comments these are...*/
// 1. Try and get rid of the need for a corner frequency
// 2. Look at how white the residuals actually are
// 3. It seems as if most of the specX are okay, but the average is high.  Look at this!
// 4. Try not as steep red noise


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "TKspectrum.h"
#include "TKfit.h"
#include "T2toolkit.h"
#include "constraints.h"
#include "choleskyRoutines.C"

using namespace std;

double OMEGA0=0; 

long double toffset = 52601.0L;
void calculateSpectrum(pulsar *psr, double T, int nSpec, double *px, double *py_r, double *py_i,int outWhite,int outUinv);

void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
   char parFile[MAX_PSR][MAX_FILELEN];
   char timFile[MAX_PSR][MAX_FILELEN];
   double tspan = -1;
   double minx,maxx;
   int nSpec=0;
   int newpar=0;
   char newparname[MAX_FILELEN];
   char xspec=0;
   char autoNspec=1;
   char fname[MAX_FILELEN];
   FILE *fout;
   double xfactor=1.0;
   char use_ccache=0;
   int cache_hits=0;
   char roundem = 0;
   int outWhite=0;
   int outUinv=0;

   int i,p;
   double globalParameter;
   const char *CVS_verNum = "$Revision: 1.18 $";

   if (displayCVSversion == 1) CVSdisplayVersion("cholSpecra_plug.C","plugin",CVS_verNum);

   *npsr = 0;  /* Now we do more than one pulsar */

   printf("Graphical Interface: cholSpectra\n");
   printf("Author:              M. Keith, G. Hobbs\n");
   printf("CVS Version:         $Revision: 1.18 $\n");
   printf(" --- type 'h' for help information\n");

   /* Obtain all parameters from the command line */
   for (i=2;i<argc;i++)
   {
	  if (strcmp(argv[i],"-f")==0) {
		 strcpy(parFile[*npsr],argv[++i]); 
		 strcpy(timFile[*npsr],argv[++i]);
		 (*npsr)++;
	  }
	  if ((strcmp(argv[i],"-dcf")==0) || (strcmp(argv[i],"-chol")==0)){
		 strcpy(covarFuncFile,argv[++i]);
	  }
	  if (strcmp(argv[i],"-nspec")==0) {
		 nSpec=atoi(argv[++i]);
		 autoNspec=0;
	  }
	  if (strcmp(argv[i],"-xspec")==0) {
		 xspec=1;
	  }

	  if (strcmp(argv[i],"-cache")==0) {
		 use_ccache=1;
	  }
	  if (strcmp(argv[i],"-round")==0) {
		 roundem=1;
	  }
	  if (strcmp(argv[i],"-outWhite")==0) {
		 outWhite=1;
	  }
	  if (strcmp(argv[i],"-outUinv")==0) {
		 outUinv=1;
	  }

	  if (strcmp(argv[i],"-yr")==0) {
		 xfactor=365.25;
	  }

	  if (strcmp(argv[i],"-tspan")==0) {
		 tspan=atof(argv[++i]);
	  }

	  else if (strcmp(argv[i],"-outpar")==0){
		 newpar=1;
		 strcpy(newparname,argv[++i]);
	  }

   }
   if (*npsr < 1){
	  logerr("Need at least 1 pulsar!");
	  exit(1);
   }

   readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
   readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
   preProcess(psr,*npsr,argc,argv);

   for (p=0; p < *npsr ; p++){

	  sortToAs(psr+p); // sort the ToAs
   }
   logmsg("Do the fitting");


   for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
   {
	  formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
	  formResiduals(psr,*npsr,1);    /* Form the residuals                 */
	  if (i==0) doFitAll(psr,*npsr,covarFuncFile);   /* Do the fitting     */
	  else textOutput(psr,*npsr,globalParameter,0,0,newpar,newparname);  /* Display the output */
   }

   for (int p = 0; p < *npsr; p++){
	  for (int param=0; param < MAX_PARAMS; param++){
		 psr[p].param[param].fitFlag[0]=0;
	  }
   }

   if(xfactor!=1.0){
	  logmsg("Multiplying X axis by %lf",xfactor);
   }

   if(xspec){
	  if(autoNspec)nSpec=16;
	  int p1,p2;
	  double crossR,crossI,P1,P2;
	  double *p1_x,*p2_x,*p1_yr,*p2_yr,*p1_yi,*p2_yi;
	  p1_x = (double*)malloc(nSpec*sizeof(double));
	  p2_x = (double*)malloc(nSpec*sizeof(double));
	  p1_yr = (double*)malloc(nSpec*sizeof(double));
	  p1_yi = (double*)malloc(nSpec*sizeof(double));
	  p2_yr = (double*)malloc(nSpec*sizeof(double));
	  p2_yi = (double*)malloc(nSpec*sizeof(double));

	  double **cc_x, **cc_yr, **cc_yi;
	  if (use_ccache){
		 cc_x = (double**)malloc(sizeof(double*)*(*npsr));
		 cc_yr= (double**)malloc(sizeof(double*)*(*npsr));
		 cc_yi= (double**)malloc(sizeof(double*)*(*npsr));
		 for (p1 = 0; p1 < *npsr; p1++){
			cc_x[p1]=NULL;
		 }
	  }
	  if(roundem){
		 for (p1 = 0; p1 < *npsr; p1++){
			psr[p1].param[param_start].val[0] = double(int(psr[p1].param[param_start].val[0]/50.0)*50.0);
			psr[p1].param[param_finish].val[0] = double(int(psr[p1].param[param_finish].val[0]/50.0+1)*50.0);

		 }
	  }

	  logmsg("Producing Cross Spectra");
	  for (p1 = 0; p1 < *npsr; p1++){
		 for (p2 = p1+1; p2 < *npsr; p2++){
			double p1_start_o = psr[p1].param[param_start].val[0];
			double p2_start_o = psr[p2].param[param_start].val[0];
			double p1_finish_o = psr[p1].param[param_finish].val[0];
			double p2_finish_o = psr[p2].param[param_finish].val[0];

			double maxx = min(p1_finish_o,p2_finish_o);
			double minx = max(p1_start_o,p2_start_o);
			logmsg("%s X %s -- Overlap: %lf -> %lf",psr[p1].name,psr[p2].name,minx,maxx);
			if ((maxx - minx < 1)){
			   logmsg("Less than 1 day overlap, cannot continue");
			   continue;
			}
			// Here we set the start/finish so that it only computes spectra of overlap data.
			psr[p1].param[param_start].paramSet[0]=1;
			psr[p1].param[param_start].fitFlag[0]=1;
			psr[p1].param[param_finish].paramSet[0]=1;
			psr[p1].param[param_finish].fitFlag[0]=1;

			psr[p1].param[param_f].fitFlag[0]=1;
			psr[p1].param[param_f].fitFlag[1]=1;

			psr[p2].param[param_start].paramSet[0]=1;
			psr[p2].param[param_start].fitFlag[0]=1;
			psr[p2].param[param_finish].paramSet[0]=1;
			psr[p2].param[param_finish].fitFlag[0]=1;

			psr[p2].param[param_f].fitFlag[0]=1;
			psr[p2].param[param_f].fitFlag[1]=1;

			psr[p1].param[param_start].val[0]=minx;
			psr[p2].param[param_start].val[0]=minx;
			psr[p1].param[param_finish].val[0]=maxx;
			psr[p2].param[param_finish].val[0]=maxx;

			double realTspan=0;
			if (tspan < 0)
			   realTspan = maxx-minx;
			else{
			   if (1.2*(maxx-minx) < realTspan){
				  logerr("Refusing to make tspan larger than 1.2 times actual data span");
				  realTspan=maxx-minx;
			   }
			   printf("NOTE: using tspan=%lf, default would have been %lf\n",tspan,maxx-minx);
			   realTspan=tspan;
			}
			realTspan=-realTspan; // this tells it to use this to set omega, rather than scaling by nsamples/(nsamples-1)
			logmsg("Tspan=%lf",realTspan);

			verbose_calc_spectra=true;
			char dop1=1;
			char dop2=1;
			if (use_ccache){
			   if (cc_x[p1] != NULL && p1_start_o == minx && p1_finish_o == maxx){
				  // can use cached p1
				  memcpy(p1_x,cc_x[p1],sizeof(double)*nSpec);
				  memcpy(p1_yr,cc_yr[p1],sizeof(double)*nSpec);
				  memcpy(p1_yi,cc_yi[p1],sizeof(double)*nSpec);
				  logmsg("CACHED %s",psr[p1].name);
				  cache_hits++;
				  dop1=0;
			   }
			   if (cc_x[p2] != NULL && p2_start_o == minx && p2_finish_o == maxx){
				  // can use cached p2
				  memcpy(p2_x,cc_x[p2],sizeof(double)*nSpec);
				  memcpy(p2_yr,cc_yr[p2],sizeof(double)*nSpec);
				  memcpy(p2_yi,cc_yi[p2],sizeof(double)*nSpec);
				  logmsg("CACHED %s",psr[p2].name);
				  cache_hits++;
				  dop2=0;
			   }
			}
			if (dop1){
			   formBatsAll(psr+p1,1);        
			   formResiduals(psr+p1,1,1);   
			   doFitAll(psr+p1,1,covarFuncFile);   
			   formBatsAll(psr+p1,1);        
			   formResiduals(psr+p1,1,1);   
			   textOutput(psr+p1,1,globalParameter,0,0,newpar,newparname);
			   calculateSpectrum(psr+p1,realTspan,nSpec,p1_x,p1_yr,p1_yi,outWhite,outUinv);
			}
			if (dop2){
			   doFitAll(psr+p2,1,covarFuncFile);   
			   formBatsAll(psr+p2,1);        
			   formResiduals(psr+p2,1,1);   
			   textOutput(psr+p2,1,globalParameter,0,0,newpar,newparname);
			   calculateSpectrum(psr+p2,realTspan,nSpec,p2_x,p2_yr,p2_yi,outWhite,outUinv);
			}
			if (use_ccache){
			   if (cc_x[p1] == NULL && p1_start_o == minx && p1_finish_o == maxx){
				  // can cache p1
				  cc_x[p1] = (double*)malloc(sizeof(double)*nSpec);
				  cc_yr[p1] = (double*)malloc(sizeof(double)*nSpec);
				  cc_yi[p1] = (double*)malloc(sizeof(double)*nSpec);
				  memcpy(cc_x[p1],p1_x,sizeof(double)*nSpec);
				  memcpy(cc_yr[p1],p1_yr,sizeof(double)*nSpec);
				  memcpy(cc_yi[p1],p1_yi,sizeof(double)*nSpec);
			   }
			   if (cc_x[p2] == NULL && p2_start_o == minx && p2_finish_o == maxx){
				  // can cache p2
				  cc_x[p2] = (double*)malloc(sizeof(double)*nSpec);
				  cc_yr[p2] = (double*)malloc(sizeof(double)*nSpec);
				  cc_yi[p2] = (double*)malloc(sizeof(double)*nSpec);
				  memcpy(cc_x[p2],p2_x,sizeof(double)*nSpec);
				  memcpy(cc_yr[p2],p2_yr,sizeof(double)*nSpec);
				  memcpy(cc_yi[p2],p2_yi,sizeof(double)*nSpec);
			   }
			}
			sprintf(fname,"%s.%s.xspec",psr[p1].name,psr[p2].name);
			fout=fopen(fname,"w");
			for (i=0; i < nSpec; i++){
			   //365.25*crossX[i],crossY_r[i],crossY_i[i],py_r1[i]*py_r1[i]+py_i1[i]*py_i1[i],py_r2[i]*py_r2[i]+py_i2[i]*py_i2[i])
			   // X = 1 x 2*
			   if (p1_x[i]!=p2_x[i]){
				  logerr("X scale for %s not same as for %s!!! %lf %lf",psr[p1].name,psr[p2].name,p1_x[i],p2_x[i]);
				  exit(2);
			   }
			   crossR = p1_yr[i]*p2_yr[i] + p1_yi[i]*p2_yi[i];
			   crossI = p1_yi[i]*p2_yr[i] - p1_yr[i]*p2_yi[i];
			   P1     = p1_yr[i]*p1_yr[i] + p1_yi[i]*p1_yi[i];
			   P2     = p2_yr[i]*p2_yr[i] + p2_yi[i]*p2_yi[i];
			   fprintf(fout,"%.8g\t%.4g\t%.4g\t%.4g\t%.4g\n",xfactor*p1_x[i],crossR,crossI,P1,P2);
			}
			fclose(fout);

			// Set the start/finish back to orig values for the next iteration.
			psr[p1].param[param_start].val[0]=p1_start_o;
			psr[p2].param[param_start].val[0]=p2_start_o;
			psr[p1].param[param_finish].val[0]=p1_finish_o;
			psr[p2].param[param_finish].val[0]=p2_finish_o;



		 }
	  }

	  free(p1_x);
	  free(p2_x);
	  free(p1_yr);
	  free(p1_yi);
	  free(p2_yr);
	  free(p2_yi);
	  if(use_ccache){
		 logmsg("CACHE HITS = %d",cache_hits);
		 for (p1 = 0; p1 < *npsr; p1++){
			if(cc_x[p1]!=NULL){
			   logmsg("Free cc: %s",psr[p1].name);
			   free(cc_x[p1]);
			   free(cc_yi[p1]);
			   free(cc_yr[p1]);
			}
		 }
		 free(cc_x);
		 free(cc_yr);
		 free(cc_yi);


	  }
   }else{
	  logmsg("Producing spectra");
	  double* px;
	  double* py_r;
	  double* py_i;
	  double origTspan=tspan;
	  for (p=0; p < *npsr ; p++){
		 pulsar* thepulsar = psr+p;
		 tspan=origTspan;

		 maxx=thepulsar->param[param_finish].val[0];
		 minx=thepulsar->param[param_start].val[0];
		 if (tspan < 0)
			tspan = maxx-minx;
		 else{
			if (1.2*(maxx-minx) < tspan){
			   printf("tspan=%lg, maxx-minx=%lg",maxx-minx,tspan);
			   logerr("Refusing to make tspan larger than 1.2 times actual data span");
			   exit(1);
			}
			printf("NOTE: using tspan=%lf, default would have been %lf\n",tspan,maxx-minx);
			tspan=-tspan;
		 }

		 if (autoNspec){
			nSpec = 128;
			while (nSpec > psr->nFit && nSpec > 8){
			   nSpec /=2;
			}
		 }
		 px = (double*)malloc(sizeof(double)*nSpec);
		 py_r = (double*)malloc(sizeof(double)*nSpec);
		 py_i = (double*)malloc(sizeof(double)*nSpec);

		 logmsg("Doing the calculation %s Nspec=%d\n",thepulsar->name,nSpec);
		 verbose_calc_spectra=true;
		 calculateSpectrum(thepulsar,tspan,nSpec,px,py_r,py_i,outWhite,outUinv);
		 logmsg("Write files");
		 sprintf(fname,"%s.spec",thepulsar->name);
		 fout = fopen(fname,"w");
		 for (i=0;i<nSpec;i++)
			fprintf(fout,"%g %g %g %g\n",xfactor*px[i],py_r[i]*py_r[i]+py_i[i]*py_i[i],py_r[i],py_i[i]);
		 fclose(fout);
		 free(px);
		 free(py_r);
		 free(py_i);

		 logmsg("Finished\n");
	  }
   }
   return 0;
}

char * plugVersionCheck = TEMPO2_h_VER;

void calculateSpectrum(pulsar *psr, double T, int nSpec, double *px, double *py_r, double *py_i,int outWhite,int outUinv) {
   int i;
   double **uinv;
   FILE *fin;
   int ndays=0;
   double resx[psr->nobs],resy[psr->nobs],rese[psr->nobs];
   int ip[psr->nobs];


   logmsg("calculateSpectrum '%s' %f %d",psr->name,T,nSpec);
// only do between start and finish

   int nobs=0;
   for (i=0;i<psr->nobs;i++){
	  if(psr->obsn[i].deleted !=0)continue;
	  if (psr->param[param_start].paramSet[0]==1 && psr->param[param_start].fitFlag[0]==1 &&
			(psr->param[param_start].val[0] > psr->obsn[i].sat))
		 continue;
	  if (psr->param[param_finish].paramSet[0]==1 && psr->param[param_finish].fitFlag[0]==1 &&
			psr->param[param_finish].val[0] < psr->obsn[i].sat)
		 continue;
	  nobs++;
   }

   logmsg("Total obs = %d, Good obs = %d",psr->nobs,nobs);
   //  printf("Calculating the spectrum\n");
   uinv = malloc_uinv(nobs);

   int ir=0;
   for (i=0;i<psr->nobs;i++)
   {
	  if(psr->obsn[i].deleted !=0)continue;
	  if (psr->param[param_start].paramSet[0]==1 && psr->param[param_start].fitFlag[0]==1 &&
			(psr->param[param_start].val[0] > psr->obsn[i].sat))
		 continue;
	  if (psr->param[param_finish].paramSet[0]==1 && psr->param[param_finish].fitFlag[0]==1 &&
			psr->param[param_finish].val[0] < psr->obsn[i].sat)
		 continue;

	  resx[ir] = (double)(psr->obsn[i].sat-toffset);
	  resy[ir] = (double)(psr->obsn[i].residual);
	  rese[ir] = (double)(psr->obsn[i].toaErr*1.0e-6);
	  ip[ir]=i;
	  ir++;
   }


   //  printf("Allocated memory\n");

   // Read in the covariance function
   if (strcmp(covarFuncFile,"NULL")==0){
	  sprintf(covarFuncFile,"covarFunc.dat_%s",psr->name);
	  logmsg("Warning: Assuming -dcf %s",covarFuncFile);
   }

   logmsg("Get Cholesky 'uinv' matrix from '%s'",covarFuncFile);
   getCholeskyMatrix(uinv,covarFuncFile,psr,resx,resy,rese,nobs,0,ip);
   logdbg("Got uinv, now compute spectrum.");

   // Must calculate uinv for the pulsar
   //
   //   int calcSpectra_ri_T(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY_R,double *specY_I,int nfit,double T,char fitfuncMode, pulsar* psr)
   char mde = 'N';
   if (T < 0){
	  mde='T';
	  T=-T;
   }
   calcSpectra_ri_T(uinv,resx,resy,nobs,px,py_r,py_i,nSpec,T,mde,psr);

   if (outUinv==1)
     {
       FILE *fout;
       int i,j;

       fout = fopen("cholSpectra.uinv","w");
       for (i=0;i<nobs;i++){
	 for (j=0;j<nobs;j++)
	   {
	     //	     printf("%d %d %g\n",i,j,uinv[i][j]);
	     fprintf(fout,"%d %d %g\n",i,j,uinv[i][j]);
	   }
       }
       fclose(fout);
     }
   if (outWhite==1)
     {
       double outRes[nobs];
       FILE *fout;

       printf("Outputting whitened residuals\n");
       T2getWhiteRes(resx,resy,rese,nobs,uinv,outRes);
       if (!(fout = fopen("whiteRes.dat","w"))){
	 printf("Unable to open file whiteRes.dat\n");
	 exit(1);
       }
       for (i=0;i<nobs;i++)
	 {
	   printf("%g %g %g %g %.5Lf\n",resx[i],resy[i],rese[i],outRes[i],psr->obsn[i].sat);
	   fprintf(fout,"%g %g %g %g %.5Lf\n",resx[i],resy[i],rese[i],outRes[i],psr->obsn[i].sat);
	 }
       fclose(fout);
     }

   // Free uinv
   free_uinv(uinv);

}


