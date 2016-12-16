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
#include "TKlongdouble.h"
#include "TKlog.h"
#include "TKsvd.h"
#include "TKmatrix.h"
#include "T2toolkit.h"
#include "T2accel.h"
#include "TKfit.h"
#include "TKrobust.h"

void TKremovePoly_f(float *px,float *py,int n,int m)
{
    int i,j;
    double x[n],y[n];
    double p[m];
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
    return TKrobustConstrainedLeastSquares(b,white_b,
            designMatrix,white_designMatrix,NULL,
            n,nf,0,tol,rescale_errors,outP, e,cvm,0);

}
double TKconstrainedLeastSquares(double* b, double* white_b,
        double** designMatrix, double** white_designMatrix,
        double** constraintsMatrix,
        int n,int nf, int nconstraints, double tol, char rescale_errors,
        double* outP, double* e, double** cvm){
    return TKrobustConstrainedLeastSquares(b,white_b,
            designMatrix,white_designMatrix,constraintsMatrix,
            n,nf,nconstraints,tol,rescale_errors,outP, e,cvm,0);

}



/*********************************************
 * Robust fitting. Options are
 * 0       Disable robust
 * H       Huber
 * B       Bisquare
 * R       Hampel
 * W       Welsch
 **********************************************/
double TKrobustLeastSquares(double* b, double* white_b,
        double** designMatrix, double** white_designMatrix,
        int n,int nf, double tol, char rescale_errors,
        double* outP, double* e, double** cvm, char robust){
    return TKrobustConstrainedLeastSquares(b,white_b,
            designMatrix,white_designMatrix,NULL,
            n,nf,0,tol,rescale_errors,outP, e,cvm,0);


}


double TKrobustConstrainedLeastSquares(double* data, double* white_data,
        double** designMatrix, double** white_designMatrix,
        double** constraintsMatrix, int ndata,int nparams, int nconstraints, double tol, char rescale_errors,
        double* outP, double* e, double** Ocvm, char robust){

    if (robust > 48 ){
       return TKrobust(data,white_data,
               designMatrix, white_designMatrix, constraintsMatrix, ndata, nparams, nconstraints, tol,
               rescale_errors,outP,e,Ocvm, robust);
    }//end of if robust



    double chisq = 0;
    int i,j,k;

    logdbg("TKleastSquares ndata=%d nparams=%d nconstraints=%d",ndata,nparams,nconstraints);
    if(nparams > ndata + nconstraints){
        logerr("Number of fit parameters exceeds number of data points\nFit will crash");
        return 0;
    }
    bool computeErrors = (e!=NULL);
    bool computeCVM = (Ocvm!=NULL);
    bool computeParam = (outP!=NULL && data!=NULL);
    bool needToFreeCVM=false;  

    if(computeErrors && ! computeCVM){
        // we can't easily compute the errors without the CVM matrix
        // so we have to create one.
        computeCVM=true;
    }
    double** cvm=NULL;
    if (computeCVM){
        cvm=malloc_uinv(nparams);
        needToFreeCVM=true;
    }

    if((writeResiduals&1) && white_data!=NULL && data != NULL){
        logdbg("Writing out whitened residuals");
        FILE* wFile=fopen("prefit.res","w");
        if (!wFile){
            printf("Unable to write out whitened residuals: cannot open file prefit.res\n");
        }
        else
        {
            for (i=0;i<ndata;i++){
                fprintf(wFile,"%d %lg %lg\n",i,data[i],white_data[i]);
            }
            fclose(wFile);
        }

    }
    if(writeResiduals&2){
        logdbg("Writing out design matrix");
        FILE * wFile=fopen("design.matrix","w");
        if (!wFile){
            printf("Unable to write out design matrix: cannot open file design.matrix\n");
        }
        else
        {
            for (i=0;i<ndata;i++) {
                for (j=0;j<nparams;j++){
                    fprintf(wFile,"%d %d %lg %lg\n",i,j,designMatrix[i][j],white_designMatrix[i][j]);
                }
                fprintf(wFile,"\n");
            }
            fclose(wFile);
        }
        wFile=fopen("constraints.matrix","w");
        if (!wFile){
            printf("Unable to write out constraints matrix: cannot open file constraints.matrix\n");
        }
        else
        {
            for (i=0;i<nconstraints;i++) {
                for (j=0;j<nparams;j++){
                    fprintf(wFile,"%d %d %lg\n",i,j,constraintsMatrix[i][j]);
                }
                fprintf(wFile,"\n");
            }
            fclose(wFile);
        }

    }


#ifdef ACCEL_LSQ
    if(useT2accel==2){

        double augmented_white_data[ndata+nconstraints];
        double **augmented_DM=NULL;
        if (computeParam){
            augmented_DM  = malloc_blas(ndata+nconstraints,nparams);
            for (i=0;i<ndata;i++){
                augmented_white_data[i] = white_data[i];
            }
        }
        for (i=0;i<ndata;i++){
            for (j=0;j<nparams;j++) augmented_DM[i][j] = white_designMatrix[i][j];
        }

        // add the extra equations for constraints to the end of the least-squares problem.
        logmsg("QR nparams=%d nconst=%d",nparams,nconstraints);
        for (i=0;i<nconstraints;i++){
            augmented_white_data[i+ndata] = 0;
            for (j=0;j<nparams;j++) {
                augmented_DM[i+ndata][j] = constraintsMatrix[i][j];
                //if(i==j)logmsg("Cmatrix ic=%d ip=%d %lg",i,j,constraintsMatrix[i][j]);
            }
        }
        if(writeResiduals&2){
            logdbg("Writing out augmented design matrix");
            FILE * wFile=fopen("adesign.matrix","w");
            if (!wFile){
                printf("Unable to write out augmented design matrix: cannot open file adesign.matrix\n");
            }
            else
            {
                for (i=0;i<ndata+nconstraints;i++) {
                    for (j=0;j<nparams;j++){
                        fprintf(wFile,"%d %d %lg\n",i,j,augmented_DM[i][j]);
                    }
                    fprintf(wFile,"\n");
                }
                fclose(wFile);
            }
        }


        chisq = accel_lsq_qr(augmented_DM,augmented_white_data,outP,ndata+nconstraints,nparams,cvm);
        rescale_errors=false;
        free_blas(augmented_DM);

        if (computeErrors){
            for (i=0;i<nparams;i++){e[i]=sqrt(cvm[i][i]);}
        }

    } else {
#endif
        // quad precision arrays for fitting if using SVD
        // the augmented data matrix
        longdouble augmented_white_data[ndata+nconstraints];

        // the augmented design matrix
        longdouble **augmented_DM = malloc_2dLL(ndata+nconstraints,nparams);
        longdouble **v=malloc_2dLL(nparams,nparams);
        longdouble **u=malloc_2dLL(ndata+nconstraints,nparams);
        longdouble w[nparams];
        longdouble wt[nparams];
        longdouble p[nparams];
        // other variables
        longdouble sum,wmax;

        logdbg("TKleastSquares()");

        // Now go to longdouble precision and augment the DM and data vector
        for (i=0;i<ndata;i++){
            if(computeParam) augmented_white_data[i] = white_data[i];
            for (j=0;j<nparams;j++) augmented_DM[i][j] = white_designMatrix[i][j];
        }

        // add the extra equations for constraints to the end of the least-squares problem.
        logmsg("SVD nparams=%d nconst=%d",nparams,nconstraints);
        for (i=0;i<nconstraints;i++){
            augmented_white_data[i+ndata] = 0;
            for (j=0;j<nparams;j++) {
                augmented_DM[i+ndata][j] = constraintsMatrix[i][j];
                //if(i==j)logmsg("Cmatrix ic=%d ip=%d %lg",i,j,constraintsMatrix[i][j]);
            }
        }

        /* Now carry out the singular value decomposition */
        // note that this modifies svd_M
        logdbg("Do SVD");
        TKsingularValueDecomposition_lsq(augmented_DM,ndata+nconstraints,nparams,v,w,u);

        wmax = TKfindMax_Ld(w,nparams);
        longdouble sensible_wmax=pow(2,sizeof(longdouble)*8-17);
        if (wmax > sensible_wmax){
            logerr("Warning: wmax very large. Precision issues likely to break fit\nwmax=%lf\ngood=%lf",(double)wmax,(double)sensible_wmax);
        }

        for (i=0;i<nparams;i++)
        {
            if (w[i] < tol*wmax) w[i]=0.0;
        }
        /* Back substitution */


        /* Now form the covariance matrix */

        if(computeCVM){
            logdbg("Compute CVM");
            for (i=0;i<nparams;i++)
            {
                if (w[i]!=0) wt[i] = 1.0/w[i]/w[i];
                else wt[i] = 0.0;
            }
            for (i=0;i<nparams;i++)
            {
                for (j=0;j<=i;j++)
                {
                    sum=0.0;
                    for (k=0;k<nparams;k++)
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
                    for (i=0;i<nparams;i++)
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
                for (i=0;i<nparams;i++){e[i]=sqrt(cvm[i][i]);}
            }

        }

        if (computeParam) {

            logdbg("Compute Params");

            logdbg("Do backsubstitution");
            TKbacksubstitution_svd(v, w, augmented_DM, augmented_white_data, p, ndata+nconstraints, nparams);

            for (k=0;k<nparams;k++)outP[k]=(double)(p[k]);

            // compute chisq
            chisq = 0.0;
            for (j=0;j<ndata;j++)
            {
                sum = 0.0;
                for (k=0;k<nparams;k++)
                {
                    sum+=p[k]*white_designMatrix[j][k];
                }
                chisq += pow((white_data[j]-sum),2);
            }
        } // computeParam
        // this funny method of freeing is because of the BLAS style matricies. M.Keith 2012
        free_2dLL(v);     // free-TKleastSquares_svd_psr_dcm-v**
        free_2dLL(u);     // free-TKleastSquares_svd_psr_dcm-u**
        free_2dLL(augmented_DM);


#ifdef ACCEL_LSQ
    } // accel
#endif



    if(computeErrors && rescale_errors){
        //	    printf("Error scaling = %g [chisq = %g] [n = %d] [nf = %d]\n",sqrt(chisq/(n-nf)),(double)chisq,n,nf);
        // This is not the place for this message: this is the only thing one
        // sees when a fit is done. Perhaps move to the results
        // summary?    -- Rutger van Haasteren & Michele Vallisneri
        // printf("Error scaling = %g\n",sqrt(chisq/(n-nf)));
        for (j=0;j<nparams;j++)
            e[j] *= sqrt(chisq/(ndata+nconstraints-nparams));

    }

    if (writeResiduals&4){
        double sum,sum_w;
        FILE* wFile=fopen("postfit.res","w");
        if (!wFile){
            printf("Unable to open file postfit.res for writing\n");
        }
        else
        {
            if(outP!=NULL){
                for (i=0;i<ndata;i++)
                {
                    sum=0;
                    sum_w=0;
                    for (j=0;j<nparams;j++){
                        sum += designMatrix[i][j]*outP[j];
                        sum_w += white_designMatrix[i][j]*outP[j];
                    }
                    fprintf(wFile,"%d %lg %lg\n",i,(double)(data[i]-sum),(double)(white_data[i]-sum_w));
                }
                fclose(wFile);
            }
        }
    }
    if(needToFreeCVM){
        if (Ocvm != NULL){
            for (i=0; i < nparams; i++){
                for(j=0; j < nparams; j++){
                    Ocvm[i][j] = cvm[i][j]; // deal with the fact that the cvm matrix may not be allocated properly
                }
            }
        }
        // we created CVM, so free it
        free_uinv(cvm);
    }

    /** Robust Estimator code by Wang YiDi, Univ. Manchester 2015 **/

   
    return chisq;
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
    int    i,j;

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


