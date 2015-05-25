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

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "tempo2.h"
#include "TKsvd.h"
#include "TKmatrix.h"
#include "T2toolkit.h"
#include "TKfit.h"

void TKremovePoly_f(float *px,float *py,int n,int m)
{
   int i,j;
   double x[n],y[n];
   double p[m];
   double chisq;
   double v[m];

   for (i=0;i<n;i++)
   {
	  x[i] = (float)px[i];
	  y[i] = (float)py[i];
   }
   TKleastSquares_svd_noErr(x,y,n, p, m, TKfitPoly);  

   for (i=0;i<n;i++)
   {
	  TKfitPoly(x[i],v,m);
	  for (j=0;j<m;j++)
		 py[i] -= v[j]*p[j];
   }
}

void TKremovePoly_d(double *x,double *y,int n,int m)
{
   int i,j;
   double p[m];
   double chisq;
   double v[m];

   logdbg("Remove polynomial n=%d m=%d",n,m);
   TKleastSquares_svd_noErr(x,y,n, p, m, TKfitPoly);  

   for (i=0;i<n;i++)
   {
	  TKfitPoly(x[i],v,m);
	  for (j=0;j<m;j++)
		 y[i] -= v[j]*p[j];
   }
}

void TKfindPoly_d(double *x,double *y,int n,int m,double* p){
   int i,j;
   double chisq;
   double v[m];

   TKleastSquares_svd_noErr(x,y,n, p, m, TKfitPoly);  

}

void TKfitPoly(double x,double *v,int m)
{
   int i;
   double t=1;
   for (i=0;i<m;i++)
   {
	  v[i] = t;
	  t*=x;
   }
}


/* Least squares fitting routines */


/**
 * TKleastSquares performs a least squares fit.
 *
 *	double* b: Array of Y values.
 *	double* white_b: Array of whitened Y values. (Uinv.Y)
 *	double** designMatrix: Fit matrix
 *	double** white_designMatrix: Whitened fit matrix
 *	int n: size of "b"
 *	int nf: number of fit parameters (i.e. columns of designMatrix)
 *	double tol:  filter to remove small values of the SVD
 *	char rescale_errors: boolean to say if resultant errors should be scaled by chisq
 *	double* outP: output fit parameters
 *	double* e: output error in fit parameters
 *	double **cvm: nf*nf output covariance matrix for fit parameters.
 *
 */
double TKleastSquares(double* b, double* white_b,
	  double** designMatrix, double** white_designMatrix,
	  int n,int nf, double tol, char rescale_errors,
	  double* outP, double* e, double** cvm){
   return TKrobustLeastSquares(b,white_b,
		 designMatrix,white_designMatrix,
		 n,nf,tol,rescale_errors,outP, e,cvm,0);
}

double TKrobustLeastSquares(double* b, double* white_b,
	  double** designMatrix, double** white_designMatrix,
	  int n,int nf, double tol, char rescale_errors,
	  double* outP, double* e, double** cvm, int robust){

   double chisq = 0;
   int i,j,k;
   if (nf > MAX_PARAMS)
   {
	  printf("Number of fitted parameters, %d, is greater than MAX_PARAMS. Please update MAX_PARAMS and reinstall\n",nf);
	  exit(1);
   }

   logdbg("TKleastSquares n=%d nf=%d",n,nf);
   if(nf > n){
	  logerr("Number of fit parameters exceeds number of data points\nFit will crash");
	  for (k=0;k<nf;k++)outP[k]=0;
	  for (k=0;k<nf;k++)e[k]=0;
	  for (i=0;i<nf;i++){
		 for (k=0;k<nf;k++)cvm[i][k]=0;
	  }
	  return 0;
   }
   // quad precision arrays for fitting.
   longdouble svd_V[n];
   longdouble **svd_M=malloc_2dLL(n,nf);
   longdouble **v=malloc_2dLL(nf,nf);
   longdouble **u=malloc_2dLL(n,nf);
   longdouble w[nf],wt[nf];
   longdouble p[nf];

   logdbg("TKleastSquares()");
   // other variables
   longdouble sum,wmax,sum_w;

   bool computeErrors = (e!=NULL);
   bool computeCVM = (cvm!=NULL);
   bool computeParam = (outP!=NULL && b!=NULL);
   bool needToFreeCVM=false;

   if(computeErrors && ! computeCVM){
	  // we can't easily compute the errors without the CVM matrix
	  // so we have to create one.
	  cvm=malloc_uinv(nf);
	  computeCVM=true;
	  needToFreeCVM=true;
   }

   if(writeResiduals==1 && white_b!=NULL && b!=NULL){
	  logdbg("Writing out whitened residuals");
	  FILE* wFile=fopen("prefit.res","w");
	  if (!wFile){
		 printf("Unable to write out whitened residuals: cannot open file prefit.res\n");
	  }
	  else
	  {
		 for (i=0;i<n;i++)
			fprintf(wFile,"%d %lg %lg\n",i,b[i],white_b[i]);
		 fclose(wFile);
	  }

   }
   if(writeResiduals==1){
	  logdbg("Writing out design matrix");
	  FILE * wFile=fopen("design.matrix","w");
	  if (!wFile){
		 printf("Unable to write out design matrix: cannot open file design.matrix\n");
	  }
	  else
	  {
		 for (i=0;i<n;i++) {
			for (j=0;j<nf;j++){
			   fprintf(wFile,"%d %d %lg %lg\n",i,j,designMatrix[i][j],white_designMatrix[i][j]);
			}
			fprintf(wFile,"\n");
		 }
		 fclose(wFile);
	  }
   }

   // Now go to longdouble precision!
   for (i=0;i<n;i++){
	  for (j=0;j<nf;j++) svd_M[i][j] = white_designMatrix[i][j];
   }

   /* Now carry out the singular value decomposition */
   // note that this modifies svd_M
   logdbg("Do SVD");
   TKsingularValueDecomposition_lsq(svd_M,n,nf,v,w,u);

   wmax = TKfindMax_Ld(w,nf);
   longdouble sensible_wmax=pow(2,sizeof(longdouble)*8-17);
   if (wmax > sensible_wmax){
	  logerr("Warning: wmax very large. Precision issues likely to break fit\nwmax=%Lf\ngood=%Lf",wmax,sensible_wmax);
   }

   for (i=0;i<nf;i++)
   {
	  if (w[i] < tol*wmax) w[i]=0.0;
   }
   /* Back substitution */


   /* Now form the covariance matrix */
   if(computeCVM){
	  logdbg("Compute CVM");

	  for (i=0;i<nf;i++)
	  {
		 if (w[i]!=0) wt[i] = 1.0/w[i]/w[i];
		 else wt[i] = 0.0;     
	  }

	  for (i=0;i<nf;i++)
	  {
		 for (j=0;j<=i;j++)
		 {
			sum=0.0;
			for (k=0;k<nf;k++)
			   sum+=v[i][k]*v[j][k]*wt[k];
			cvm[i][j] = cvm[j][i] = (double)sum;
		 }
	  } 
	  if(debugFlag==1) {
		 FILE *fout;
		 fout = fopen("cvm.matrix","w");
		 if (!fout){
			printf("Unable to open file cvm.matrix for writing\n");
		 }
		 else{
			for (i=0;i<nf;i++)
			{
			   for (j=0;j<=i;j++)
			   {
				  fprintf(fout,"%+.8f ",cvm[i][j]/sqrt(cvm[i][i]*cvm[j][j]));
			   }
			   fprintf(fout,"\n");
			}
			fclose(fout);
		 }
	  }
	  if(computeErrors){
		 logdbg("Compute Errors");
		 for (i=0;i<nf;i++){e[i]=sqrt(cvm[i][i]);}
	  }
   }


   if (computeParam){

	  logdbg("Compute Params");
	  for (i=0;i<n;i++){
		 svd_V[i] = white_b[i];
	  }

	  logdbg("Do backsubstitution");
	  TKbacksubstitution_svd(v, w, svd_M, svd_V, p, n, nf);

	  for (k=0;k<nf;k++)outP[k]=(double)(p[k]);
	  // compute chisq
	  chisq = 0.0;
	  for (j=0;j<n;j++)
	  {
		 sum = 0.0;
		 for (k=0;k<nf;k++)
		 {
			sum+=p[k]*white_designMatrix[j][k];
		 }
		 chisq += pow((white_b[j]-sum),2);
	  }

	  if(computeErrors && rescale_errors){
		 //	    printf("Error scaling = %g [chisq = %g] [n = %d] [nf = %d]\n",sqrt(chisq/(n-nf)),(double)chisq,n,nf);
		 // This is not the place for this message: this is the only thing one
		 // sees when a fit is done. Perhaps move to the results
		 // summary?    -- Rutger van Haasteren & Michele Vallisneri
		 // printf("Error scaling = %g\n",sqrt(chisq/(n-nf)));
		 for (j=0;j<nf;j++)
			e[j] *= sqrt(chisq/(n-nf));
	  }

	  if (writeResiduals){
		 FILE* wFile=fopen("postfit.res","w");
		 if (!wFile){
			printf("Unable to open file postfit.res for writing\n");
		 }
		 else
		 {
			for (i=0;i<n;i++)
			{
			   sum=0;
			   sum_w=0;
			   for (j=0;j<nf;j++){
				  sum += designMatrix[i][j]*p[j];
				  sum_w += white_designMatrix[i][j]*p[j];
			   }
			   fprintf(wFile,"%d %lg %lg\n",i,(double)(b[i]-sum),(double)(white_b[i]-sum_w));
			}
			fclose(wFile);
		 }
	  }
   } // computeParam

   // this funny method of freeing is because of the BLAS style matricies. M.Keith 2012
   free_2dLL(v);     // free-TKleastSquares_svd_psr_dcm-v**
   free_2dLL(u);     // free-TKleastSquares_svd_psr_dcm-u**
   free_2dLL(svd_M);

   if(needToFreeCVM){
	  // we created CVM, so free it
	  free_uinv(cvm);
   }

   /** Robust Estimator code by Wang YiDi, Univ. Manchester 2015 **/
   if(robust > 0){
	  longdouble resid[n];//residual after the calculation of leastsquare
	  longdouble Weight[n];//reweights

	  longdouble sigma; // Median of the residuals
	  double c0 = 1.345; // robust tunning parameter
	  //robust estimation process-----wyd    
	  for (j=0;j<n;j++)
	  {
		 sum = 0.0;
		 for (k=0;k<nf;k++)
		 {
			sum += p[k]*white_designMatrix[j][k];
		 }

		 resid[j]=white_b[j]-sum;
	  }

	  //Finding the median of residuals
	  longdouble median = 0.0;
	  longdouble newval[n];
	  longdouble resid_new[n];
	  longdouble store;
	  int changed;
	  int count = n;

	  for (i=0;i<count;i++)
	  {
		 newval[i]=fabs(resid[i]);
	  }

	  do
	  {
		 changed=0;
		 for (j=0;j<n-1;j++)
		 {
			if (fabs(resid[i]) > fabs(resid[i+1]))
			{
			   store=newval[i+1];
			   newval[i+1]=newval[i];
			   newval[i]=store;
			   changed=1;
			}
		 }
	  }while (changed==1);


	  if (count%2==0)
	  {
		 median = (newval[count/2-1]+newval[count/2])/2.0;
	  }
	  else
	  {
		 median = newval[(count-1)/2];
	  }

	  if(writeResiduals){
		 FILE *fpp = fopen("Wmedian.txt", "w");
		 fprintf(fpp, "%lg\n", (double)median);
		 fclose(fpp);
	  }

	  //Reweight 
	  c0 = 1.345;

	  for (j=0;j<n;j++)
	  {
		 resid_new[j] = resid[j]/median;

		 if (fabs(resid_new[j])<c0)
		 {
			Weight[j] = 1.0;
		 }
		 else
		 {
			Weight[j] = c0/fabs(resid_new[j]);
		 }
	  }

	  for (i=0; i<n; i++)
	  {
		 white_b[i] *= Weight[i];

		 for (j=0; j<nf; j++)  white_designMatrix[i][j] *= Weight[i];
	  }

	  if(writeResiduals){
		 FILE *pp = fopen("Weight.txt", "w");
		 for (i=0; i<n; i++)
		 {
			fprintf(pp, "%lg\n", (double)Weight[i]);
		 }
		 fclose(pp);
	  }

	  //Calcluation via robust estimation algorithm
	  double newchisq = TKrobustLeastSquares(b, white_b, designMatrix, white_designMatrix, n, nf, tol, rescale_errors, outP, e, cvm, robust-1);

	  if(writeResiduals){
		 FILE *pp = fopen("Parameter_RLS.txt", "w");
		 for (i=0; i<nf; i++)
		 {
			fprintf(pp, "%lg\n", outP[i]);
		 }
		 fclose(pp);
		 FILE *fp_result = fopen("Wchisq.txt", "w");
		 fprintf(fp_result, "WLS=%lg   RLS=%lg\n", chisq, newchisq);
		 fclose(fp_result);
	  }
	  chisq=newchisq;

   }
   // End of Robust estimator code

   return chisq;
}





// routine for pulsar fitting.
void TKleastSquares_single_pulsar(double *x,double *y,int n,double *outP,double *e,int nf,double **cvm, double *chisq, void (*fitFuncs)(double, double [], int,pulsar *,int, int),pulsar *psr,double tol, int *ip,char rescale_errors, double **uinv) {

   double **designMatrix, **white_designMatrix;
   double *b,*white_b;

   TKfit_getPulsarDesignMatrix(x,y,n,nf,fitFuncs,psr,ip,uinv,0,&designMatrix,&white_designMatrix,&b,&white_b);
   *chisq = TKrobustLeastSquares(b,white_b,designMatrix,white_designMatrix,
		 n,nf,tol,rescale_errors,
		 outP,e,cvm,psr->robust);
   free_blas(designMatrix); // free-TKleastSquares_svd_psr_dcm-designMatrix**
   free_blas(white_designMatrix);  // free-TKleastSquares_svd_psr_dcm-white_designMatrix**
   free(b);
   free(white_b);

}

void TKleastSquares_global_pulsar(double **x,double **y,int *n,
	  double *outP,double *e,int* nf, int nglobal,double **cvm, double *chisq, void (*fitFuncs)(double, double [], int,pulsar *,int,int),pulsar *psr,double tol, int **ip,char rescale_errors, double ***uinv, int npsr) {

   double **designMatrix, **white_designMatrix;
   double **psr_DM, **psr_wDM;
   double *b,*white_b, *psr_b,*psr_wb;
   int ipsr;
   int totalFit=0;
   int totalObs=0;
   int i,j;
   int off_r=0;
   int off_f=0;

   for (ipsr=0; ipsr < npsr; ipsr++){
	  totalFit+=nf[ipsr];
	  totalObs+=n[ipsr];
   }
   totalFit+=nglobal;

   white_designMatrix=malloc_blas(totalObs,totalFit);
   designMatrix=malloc_blas(totalObs,totalFit);
   b=(double*)calloc(totalObs,sizeof(double));
   white_b=(double*)calloc(totalObs,sizeof(double));

   for (ipsr=0; ipsr < npsr; ipsr++){
	  logdbg("Getting design matrix / whitened residuals for psr %d    off_r=%d off_f=%d nglobal=%d",ipsr,off_r,off_f,nglobal);
	  TKfit_getPulsarDesignMatrix(x[ipsr],y[ipsr],n[ipsr],nf[ipsr]+nglobal,fitFuncs,psr,ip[ipsr],uinv[ipsr],ipsr,&psr_DM,&psr_wDM,&psr_b,&psr_wb);

	  // the global fit parameters
	  for(i=0; i < n[ipsr]; i++){
		 for(j=0; j < nglobal; j++){
			designMatrix[i+off_r][j] = psr_DM[i][j];
			white_designMatrix[i+off_r][j] = psr_wDM[i][j];
		 }
	  }
	  // the regular fit parameters
	  for(i=0; i < n[ipsr]; i++){
		 for(j=0; j < nf[ipsr]; j++){
			designMatrix[i+off_r][j+off_f+nglobal] = psr_DM[i][j+nglobal];
			white_designMatrix[i+off_r][j+off_f+nglobal] = psr_wDM[i][j+nglobal];
		 }
	  }
	  // the residuals
	  for(i=0; i < n[ipsr]; i++){
		 b[i+off_r] = psr_b[i];
		 white_b[i+off_r] = psr_wb[i];
	  }
	  // increment the offset.
	  off_r += n[ipsr];
	  off_f += nf[ipsr];

	  // free temp matricies.
	  free_blas(psr_DM);
	  free_blas(psr_wDM);
	  free(psr_b);
	  free(psr_wb);
   }


   // go ahead and do the fit!

   *chisq = TKrobustLeastSquares(b,white_b,designMatrix,white_designMatrix,
		 totalObs,totalFit,tol,rescale_errors,
		 outP,e,cvm,psr[0].robust);

   free_blas(designMatrix); // free-TKleastSquares_svd_psr_dcm-designMatrix**
   free_blas(white_designMatrix);  // free-TKleastSquares_svd_psr_dcm-white_designMatrix**
   free(b);
   free(white_b);

}



void TKfit_getPulsarDesignMatrix(double *x,double *y,int n,int nf,void (*fitFuncs)(double, double [], int,pulsar *,int,int), pulsar *psr, int* ip, double **uinv,int ipsr,double ***OUT_designMatrix,double ***OUT_white_designMatrix,double** OUT_b, double** OUT_wb){

   //double precision arrays for matrix algebra.
   double **designMatrix, **white_designMatrix;
   double basisFunc[nf];
   double *b,*white_b;
   int    i,j,k;
   int nrows=get_blas_rows(uinv);
   int ncols=get_blas_cols(uinv);
   if (ncols!=n){
	  logmsg("n=%d ncols=%d",n,ncols);
	  logerr("uinv error. Either you did not use malloc_uinv() to create uinv or np!=ncols");
	  exit(1);
   }

   if (nrows!=n && nrows != 1){
	  logmsg("n=%d nrows=%d",n,nrows);
	  logerr("uinv error. Either you did not use malloc_uinv() to create uinv or np!=nrows");
	  exit(1);
   }


   // double arrays
   white_designMatrix=malloc_blas(n,nf);
   designMatrix=malloc_blas(n,nf);
   b=(double*)malloc(sizeof(double)*n);
   white_b=(double*)malloc(sizeof(double)*n);

   /* This routine has been developed from Section 15 in Numerical Recipes */

   /* Determine the design matrix - eq 15.4.4 
	* and the vector 'b' - eq 15.4.5 
	*/
   for (i=0;i<n;i++)
   {
	  // fitFuncs is not threadsafe!
	  fitFuncs(x[i],basisFunc,nf,psr,ip[i],ipsr);
	  for (j=0;j<nf;j++) designMatrix[i][j] = basisFunc[j];
	  b[i] = y[i];
   }
   // Take into account the data covariance matrix

   if(nrows==1){
	  // we have only diagonal elements
	  for (i=0;i<n;i++){
		 white_b[i]=b[i]*uinv[0][i];
		 for (j=0;j<nf;j++){
			white_designMatrix[i][j] = designMatrix[i][j]*uinv[0][i];
		 }
	  }
   } else {
	  TKmultMatrix_sq(uinv,designMatrix,n,nf,white_designMatrix);  
	  TKmultMatrixVec_sq(uinv,b,n,white_b);
   }

   *OUT_designMatrix=designMatrix;
   *OUT_white_designMatrix=white_designMatrix;
   *OUT_b=b;
   *OUT_wb=white_b;
}


// legacy method.
void TKleastSquares_svd_psr_dcm(double *x,double *y,double *sig,int n,double *outP,double *e,int nf,double **cvm, double *chisq, void (*fitFuncs)(double, double [], int,pulsar * , int,int),int weight,pulsar *psr,double tol, int *ip,double **uinv) {
   logmsg("Warning: Deprecated method TKleastSquares_svd_psr_dcm() -> TKleastSquares_single_pulsar()");
   TKleastSquares_single_pulsar(x,y,n,outP,e,nf,cvm,chisq,fitFuncs,psr,tol,ip,(weight==0 || (weight==1 && psr->rescaleErrChisq==1)),uinv);
}

// same as above but without a uinv matrix.
void TKleastSquares_svd_psr(double *x,double *y,double *sig,int n,double *p,double *e,int nf,double **cvm, double *chisq, void (*fitFuncs)(double, double [], int,pulsar *,int,int),int weight,pulsar *psr,double tol, int *ip)
{
   logmsg("Warning: Deprecated method TKleastSquares_svd_psr() -> TKleastSquares_single_pulsar()");
   int i;
   double ** uinv=malloc_blas(1,n);
   if (weight==1){
	  for (i=0; i<n;i++){
		 uinv[0][i]=1.0/sig[i];
	  }
   } else{
	  for (i=0; i<n;i++){
		 uinv[0][i]=1.0;
	  }
   }
   TKleastSquares_single_pulsar(x,y,n,p,e,nf,cvm,chisq,fitFuncs,psr,tol,ip,(weight==0 || (weight==1 && psr->rescaleErrChisq==1)),uinv);
   free_blas(uinv);
}



void TKleastSquares_svd_noErr(double *x,double *y,int n,double *p,int nf, void (*fitFuncs)(double, double [], int))
{
   double chisq=0;
   TKleastSquares_svd(x,y,NULL,n,p,NULL,nf,NULL,&chisq,fitFuncs,0);
}

// Non-pulsar fit. No cholesky yet though...
void TKleastSquares_svd(double *x,double *y,double *sig,int n,double *p,double *e,int nf,double **cvm, double *chisq, void (*fitFuncs)(double, double [], int),int weight)
{
   double **designMatrix, **white_designMatrix;
   double basisFunc[nf];
   double *b,*white_b;
   int    i,j,k;

   // double arrays
   white_designMatrix=malloc_blas(n,nf);
   designMatrix=malloc_blas(n,nf);
   b=(double*)malloc(sizeof(double)*n);
   white_b=(double*)malloc(sizeof(double)*n);

   logdbg("Non pulsar least-squares fit. n=%d nf=%d",n,nf);
   /* Determine the design matrix - eq 15.4.4 
	* and the vector 'b' - eq 15.4.5 
	*/
   for (i=0;i<n;i++)
   {
	  // fitFuncs is not threadsafe!
	  fitFuncs(x[i],basisFunc,nf);
	  for (j=0;j<nf;j++) designMatrix[i][j] = basisFunc[j];
	  b[i] = y[i];
   }

   // deal with the weights if we are doing a weighted fit.
   if(weight==1 && sig!=NULL){
	  logdbg("Divide by errors");
	  for (i=0;i<n;i++){
		 white_b[i]=b[i]/sig[i];
		 for (j=0;j<nf;j++) white_designMatrix[i][j] = designMatrix[i][j]/sig[i];
	  }
   } else {
	  for (i=0;i<n;i++){
		 white_b[i]=b[i];
		 for (j=0;j<nf;j++) white_designMatrix[i][j] = designMatrix[i][j];
	  }

   }

   // go ahead and do the fit!

   *chisq = TKleastSquares(b,white_b,designMatrix,white_designMatrix,
		 n,nf,1e-10,1,
		 p,e,cvm);

   free_blas(designMatrix); // free-TKleastSquares_svd_psr_dcm-designMatrix**
   free_blas(white_designMatrix);  // free-TKleastSquares_svd_psr_dcm-white_designMatrix**
   free(b);
   free(white_b);


}

void TKleastSquares_svd_passN(double *x,double *y,double *sig2,int n,double *p,double *e,int nf,double **cvm, double *chisq, void (*fitFuncs)(double, double [], int,int),int weight)
{
   logerr("This method no longer implemented.");
   exit(-1);
}

longdouble TKfindMax_Ld(longdouble *x,int n)
{
   longdouble ret;
   int i;

   ret = x[0];
   for (i=0;i<n;i++)
   {
	  if (x[i] > ret) ret = x[i];
   }
   return ret;
}

