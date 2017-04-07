#include "t2fit.h"
#include "t2fit_stdFitFuncs.h"
#include "t2fit_nestlike.h"
#include "constraints.h"
#include "constraints_nestlike.h"
#include <TKfit.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <enum_str.h>



/**********
 * t2Fit replaces the old "doFit" routines with a slightly more refined code.
 *
 * t2Fit does the linear fitting part of tempo2.
 *
 *
 * To add new fit parameters:
 * 1) Create a function that computes the gradients, with the following signature
 *
 *   double myFitFunc(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
 *
 * Where, psr is the pulsar array, ipsr is the current pulsar, x is bat-pepoch, ipos is the
 * obsn index, label is the fit parameter enum and k is a sub-parameter counter. This function
 * should return the derivative of the data with respect to the fit parameter.
 *
 * 2) Create a function that updates the pulsar parameters with the post-fit values,
 * with the following signature:
 *
 *   void t2UpdateFunc_planet(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);
 *
 * where psr, ipsr, label and k are the same as above, and val and err are the post-fit value
 * and uncertanty
 *
 * 3) Add your fit parameter to t2fit_fillOneParameterFitInfo().
 *
 * Best bet is to see how the other parameters are implemented and copy that!
 *
 * Have fun!
 */


void callFitFuncPlugin(pulsar *psr,int npsr, const char *covarFuncFile); // currently in doFit.C

void t2fit_postfit(pulsar* psr, int npsr);
void t2fit_prefit(pulsar* psr, int npsr);

// Remove elements from SVD sigma matrix below this value.
#define T2_SVD_TOL 1e-27

void t2Fit(pulsar *psr,unsigned int npsr, const char *covarFuncFile){

    if (strcmp(psr[0].fitFunc,"default")!=0){
        callFitFuncPlugin(psr,npsr,covarFuncFile);
        return;
    }

    t2fit_prefit(psr,npsr);

    // if we have a model for the data covariance function, then use it.
    // Otherwise we we will just whiten using the error bars.
    bool haveCovar = (covarFuncFile!=NULL && strcmp(covarFuncFile,"NULL"));

    /**
     * Find out if there are any global parameters and what they are...
     */
    FitInfo global_fitinfo;
    t2Fit_fillGlobalFitInfo(psr,npsr,global_fitinfo);
    logdbg("Nglobal parameters = %d",global_fitinfo.nParams);

    // If we had any global parameters (or constraints) then we need to do a global fit
    // otherwise we can do a fit for each pulsar individually, which is quicker
    // and saves memory.
    bool doGlobalFit = (global_fitinfo.nParams > 0) || (global_fitinfo.nConstraints > 0);

    unsigned long long totalGlobalData=0; // the number of data points across all pulsars
    unsigned int gParams=global_fitinfo.nParams; // the number of global fit parameters
    unsigned int gConstraints=global_fitinfo.nConstraints; // the number of global constraints

    unsigned long long totalGlobalParams=gParams;
    unsigned long long totalGlobalConstraints=gConstraints;

    //    double* gX[MAX_PSR]; // "x" values for each pulsar
    double* gY[MAX_PSR]; // "y" values for each pulsar
    double* gW[MAX_PSR]; // whitened "y" values for each pulsar
    double** gDM[MAX_PSR]; // design matrix for each pulsar
    double** gWDM[MAX_PSR]; // whitened design matrix for each pulsar
    double** gCM[MAX_PSR]; // constraints matrix for each pulsar
    unsigned int gNdata[MAX_PSR]; // number of data points for each pulsar (size of x and y)

    logmsg("NEW fit routine. GlobalFit=%s",doGlobalFit ? "true" : "false");

    /**
     * However we are going to do the fit, we want to loop over all the pulsars
     * to get the input data and design matricies etc.
     */
    for (size_t ipsr=0; ipsr < npsr; ipsr++) {
        double *psr_x   = (double*)malloc(sizeof(double)*psr[ipsr].nobs);
        double *psr_y   = (double*)malloc(sizeof(double)*psr[ipsr].nobs);
        double *psr_white_y   = (double*)malloc(sizeof(double)*psr[ipsr].nobs);
        double *psr_e   = (double*)malloc(sizeof(double)*psr[ipsr].nobs);
        int *psr_toaidx = (int*)malloc(sizeof(int)*psr[ipsr].nobs); // mapping from fit data to observation number
        double** uinv; // the whitening matrix.

        /**
         * Working out which data contributes to the fit is done in this routine.
         * Basically gets values for all observations within START and FINISH which are
         * not deleted.
         *
         * returns the number of data points.
         */
        const unsigned int psr_ndata = t2Fit_getFitData(psr+ipsr,psr_x,psr_y,psr_e,psr_toaidx);
        assert(psr_ndata > 0u);
        psr[ipsr].nFit = psr_ndata; // pulsar.nFit is the number of data points used in the fit.


        /**
         * Now we work out which parameters are being fit for, how many parameters,
         * and determine the gradient functions for the design matrix and the update functions
         * which update the pulsar struct.
         */
        t2Fit_fillFitInfo(psr+ipsr,psr[ipsr].fitinfo,global_fitinfo,
                psr_x,psr_toaidx, psr_ndata);

        if(writeResiduals&2){
            FILE *fp = fopen("param.labels","w");
            ///TODO
            for(int ip=0; ip < psr[ipsr].fitinfo.nParams; ++ip){
                fprintf(fp,"%d %s %d\n",ip,label_str[psr[ipsr].fitinfo.paramIndex[ip]],psr[ipsr].fitinfo.paramCounters[ip]);
            }
            fclose(fp);
        }


        /**
         * The whitening matrix behaves diferently if we have a covariance matrix.
         * If we have a covariance matrix, uinv is an ndata x ndata triangular matrix.
         * Otherwise, it only has diagonal elements, so we efficiently store it as 
         * a 1-d ndata array.
         */
        if (haveCovar) {
            // ToAs must be sorted for covariance function code
            sortToAs(psr+ipsr);

            // malloc_uinv does a blas-compatible allocation of a 2-d array.
            uinv = malloc_uinv(psr_ndata);
            psr[ipsr].fitMode=1; // Note: forcing this to 1 as the Cholesky fit is a weighted fit
            logmsg("Doing a FULL COVARIANCE MATRIX fit");
        } else {
            // Here the whitening matrix is just a diagonal
            // weighting matrix. Store diagonal matrix as 1xN
            // so that types match later.
            uinv=malloc_blas(1,psr_ndata); 
            if(psr[ipsr].fitMode == 0){
                // if we are doing an unweighted fit then we should set the errors to 1.0
                // to give uniform weighting.
                logdbg("Doing an UNWEIGHTED fit");
                for (unsigned int i=0; i < psr_ndata; i++){
                    psr_e[i]=1.0;
                }
            } else {
                logdbg("Doing a WEIGHTED fit");
            }
        }
        assert(uinv!=NULL);

        /**
         * Now we form the whitening matrix, uinv.
         * Note that getCholeskyMatrix() is clever enough to see that we 
         * have created a 1 x ndata matrix if we have only diagonal elements.
         */
        getCholeskyMatrix(uinv,covarFuncFile,psr+ipsr,
                psr_x,psr_y,psr_e,
                psr_ndata,0,psr_toaidx);

        logtchk("got Uinv");

        // define some convinience variables
        const unsigned nParams=psr[ipsr].fitinfo.nParams;
        const unsigned nConstraints=psr[ipsr].fitinfo.nConstraints;


        /**
         * The design matrix is the matrix of gradients for the least-squares.
         * If the design matrix is M, parameters p, and data d, we are solving
         * M.p = d
         * It is ndata x nparams in size. We also allocate the whitened DM here.
         */
        double** designMatrix = malloc_blas(psr_ndata,nParams);
        double** white_designMatrix = malloc_blas(psr_ndata,nParams);
        for (unsigned int idata =0; idata < psr_ndata; ++idata){
            // t2Fit_buildDesignMatrix is a replacement for the old FITfuncs routine.
            // it fills one row of the design matrix.
            t2Fit_buildDesignMatrix(psr,ipsr,psr_x[idata], psr_toaidx[idata], designMatrix[idata]);
        }


        logtchk("made design matrix");

        /**
         * The constraints matrix is similar to the design matrix, but here we are solving:
         * B.p = 0
         * Where B is the constraints matrix and p is the parameters. we solve both this
         * and the DM equation set simultaniously. TKleastSquares will do this for us.
         *
         * If there are no constraints we leave it as NULL, which is detected in TKfit as
         * no constraints anyway.
         */
        double** constraintsMatrix =NULL;
        if(psr[ipsr].fitinfo.nConstraints > 0){

            computeConstraintWeights(psr+ipsr);
            constraintsMatrix = malloc_blas(nConstraints,nParams);
            for (unsigned int iconstraint =0; iconstraint < nConstraints; ++iconstraint){
                // similar to t2Fit_buildDesignMatrix, t2Fit_buildConstraintsMatrix
                // creates one row of the constraints matrix.
                t2Fit_buildConstraintsMatrix(psr, ipsr, iconstraint, constraintsMatrix[iconstraint]);
            }
        }

        logtchk("made constraints matrix");

        /**
         * Now we multiply the design matrix and the data vector by the whitening matrix.
         * If we just have variances (uinv is diagonal) then we do it traditionally, otherwise
         * we use TKmultMatrix as this is usually backed by LAPACK and so is fast :)
         */
        if(haveCovar){
            TKmultMatrixVec(uinv,psr_y,psr_ndata,psr_ndata,psr_white_y);
            TKmultMatrix_sq(uinv,designMatrix,psr_ndata,nParams,white_designMatrix);
        } else {
            for(unsigned i=0;i<psr_ndata;++i){
                psr_white_y[i]=psr_y[i]*uinv[0][i];
                for(unsigned j=0;j<nParams;++j){
                    white_designMatrix[i][j] = designMatrix[i][j]*uinv[0][i];
                }
            }
        }

        free_blas(uinv);
        free(psr_e);
        free(psr_toaidx);

        logtchk("done whitening");
        /*
         * Now - if we are going to do a global fit, we store all the above for later
         *       otherwise
         */
        if (doGlobalFit){
            // we are going to do a global fit, so need to store the values for later
            //gX[ipsr] = psr_x;
            gY[ipsr] = psr_y;
            gW[ipsr] = psr_white_y;
            gDM[ipsr] = designMatrix;
            gWDM[ipsr] = white_designMatrix;
            gCM[ipsr] = constraintsMatrix;
            gNdata[ipsr] = psr_ndata;
            totalGlobalData += psr_ndata;
            totalGlobalParams += nParams - gParams;
            totalGlobalConstraints += nConstraints - gConstraints;
        } else {
            // NOT GLOBAL
            // so do one fit at a time...

            double chisq; // the post-fit chi-squared

            // allocate memory for the output of TKleastSquares
            double* parameterEstimates = (double*)malloc(sizeof(double)*nParams);
            double* errorEstimates = (double*)malloc(sizeof(double)*nParams);

            /*
             * Call TKleastSquares, or in fact, TKrobustConstrainedLeastSquares,
             * since we might want robust fitting and/or constraints/
             *
             * The arguments here are explained in TKfit.C
             *
             */
            chisq = TKrobustConstrainedLeastSquares(psr_y,psr_white_y,
                    designMatrix,white_designMatrix,constraintsMatrix,
                    psr_ndata,nParams,nConstraints,
                    T2_SVD_TOL,1,parameterEstimates,errorEstimates,psr[ipsr].covar,
                    psr[ipsr].robust);

            // update the pulsar struct as appropriate
            psr[ipsr].fitChisq = chisq;
            psr[ipsr].fitNfree = psr_ndata + nConstraints - nParams;


            psr[ipsr].nParam = psr[ipsr].fitinfo.nParams;
            for (int iparam = 0; iparam <  psr[ipsr].nParam; ++iparam){
                psr[ipsr].fitParamI[iparam] = psr[ipsr].fitinfo.paramIndex[iparam];
                psr[ipsr].fitParamK[iparam] = psr[ipsr].fitinfo.paramCounters[iparam];
            }

            logdbg("Updating the parameters");
            logtchk("updating the parameter values");
            /*
             * This routine calls the appropriate update functions to apply the result of the fit
             * to the origianal (non-linearised) pulsar parameters.
             */
            t2Fit_updateParameters(psr,ipsr,parameterEstimates,errorEstimates);
            logtchk("complete updating the parameter values");
            logdbg("Completed updating the parameters");

            /*
             * If we are not doing a global fit, we can clean up the memory for this pulsar.
             * Might make a difference for very large datasets.
             */
            logdbg("Free fit memory");
            free(parameterEstimates);
            free(errorEstimates);
            free_blas(designMatrix);
            free_blas(white_designMatrix);
            if (constraintsMatrix) free_blas(constraintsMatrix);
            free(psr_x);
            free(psr_y);
            free(psr_white_y);
        }
    }
    if (doGlobalFit){

        const unsigned int nobs = totalGlobalData;
        double** designMatrix = malloc_blas(nobs,totalGlobalParams);
        double** white_designMatrix = malloc_blas(nobs,totalGlobalParams);

        double** constraintsMatrix = malloc_blas(totalGlobalConstraints,totalGlobalParams);

        double *y   = (double*)malloc(sizeof(double)*nobs);
        double *white_y   = (double*)malloc(sizeof(double)*nobs);

        unsigned int off_f = gParams; // leave space for globals
        unsigned int off_r = 0;
        unsigned int off_c = gConstraints;

        logdbg("Building matricies for global fit... npsr=%u",npsr);
        logdbg("nobs=%u, totalGlobalParams=%u, totalGlobalConstraints=%u",nobs,totalGlobalParams,totalGlobalConstraints);
        logwarn("This mode is not supported yet!!!");


        for (unsigned int ipsr = 0; ipsr < npsr ; ++ipsr){
            unsigned int nLocal = psr[ipsr].fitinfo.nParams-gParams;
            logdbg("ipsr=%u, off_r = %u, off_c=%u, off_f=%u, nlocal=%u",
                    ipsr,off_r,off_c,off_f,nLocal);

            // the fit parameters
            for(unsigned int i=0; i < gNdata[ipsr]; i++){

                // the global params (they go first)
                for(unsigned int g= 0; g < gParams; g++){
                    unsigned int j = g+nLocal;
                    if(ipsr==0 && i==0 && writeResiduals&0x08){
                        logmsg("Row %d = %s %s(%d)",g,"global",label_str[global_fitinfo.paramIndex[g]],global_fitinfo.paramCounters[g]);
                    }
                    designMatrix[i+off_r][g] = gDM[ipsr][i][j];
                    white_designMatrix[i+off_r][g] = gWDM[ipsr][i][j];
                }

                for(unsigned int j=0; j < nLocal; j++){
                    if(i==0 && writeResiduals&0x08){
                        logmsg("Row %d = %s %s(%d)",j+off_f,psr[ipsr].name,label_str[psr[ipsr].fitinfo.paramIndex[j]],psr[ipsr].fitinfo.paramCounters[j]);
                    }
                    designMatrix[i+off_r][j+off_f] = gDM[ipsr][i][j];
                    white_designMatrix[i+off_r][j+off_f] = gWDM[ipsr][i][j];
                }
            }
            // the data
            for(unsigned int i=0; i < gNdata[ipsr]; ++i){
                y[i+off_r] = gY[ipsr][i];
                white_y[i+off_r] = gW[ipsr][i];
            }

            for(unsigned int i=0; i < psr[ipsr].fitinfo.nConstraints; i++){
                for(unsigned int j=0; j < nLocal; j++){
                    constraintsMatrix[i+off_c][j+off_f] = gCM[ipsr][i][j];
                }

                // the global params (they go first)
                for(unsigned int g= 0; g < gParams; g++){
                    unsigned int j = g+nLocal;
                    constraintsMatrix[i+off_c][g] = gCM[ipsr][i][j];
                }

            }

            off_r += gNdata[ipsr];
            off_f += nLocal;
            off_c += psr[ipsr].fitinfo.nConstraints;

            free(gY[ipsr]);
            free(gW[ipsr]);
            free_blas(gDM[ipsr]);
            if (gCM[ipsr]) free_blas(gCM[ipsr]);
            free_blas(gWDM[ipsr]);
        }

        double chisq; // the post-fit chi-squared

        double* parameterEstimates = (double*)malloc(sizeof(double)*totalGlobalParams);
        double* errorEstimates = (double*)malloc(sizeof(double)*totalGlobalParams);
        chisq = TKrobustConstrainedLeastSquares(y,white_y,
                designMatrix,white_designMatrix,constraintsMatrix,
                nobs,totalGlobalParams,totalGlobalConstraints,
                T2_SVD_TOL,1,parameterEstimates,errorEstimates,psr[0].covar,
                psr[0].robust);
        // for now the CVM ends up in psr[0].covar.

        int off_p = gParams;
        for (unsigned int ipsr = 0; ipsr < npsr ; ++ipsr){
            // update the pulsar struct as appropriate
            psr[ipsr].fitChisq = chisq;
            psr[ipsr].fitNfree = nobs + totalGlobalConstraints - totalGlobalParams;



            double* psr_parameterEstimates = (double*)malloc(sizeof(double)*psr[ipsr].fitinfo.nParams);
            double* psr_errorEstimates = (double*)malloc(sizeof(double)*psr[ipsr].fitinfo.nParams);
            const unsigned np = psr[ipsr].fitinfo.nParams-gParams;

            /* extract the fit output for the individual pulsars. 
             * I.e. detangle the global fit
             * Notice: Globals go at the end of the individual pulsar arrays.
             */
            for (unsigned i = 0; i < np; ++i){
                psr_parameterEstimates[i] = parameterEstimates[off_p];
                psr_errorEstimates[i] = errorEstimates[off_p];
                ++off_p;
            }

            for (unsigned i = 0; i < gParams; ++i){
                psr_parameterEstimates[i+np] = parameterEstimates[i];
                psr_errorEstimates[i+np] = errorEstimates[i];
            }

            logdbg("Updating the parameters");
            logtchk("updating the parameter values");

            /*
             * This routine calls the appropriate update functions to apply the result of the fit
             * to the origianal (non-linearised) pulsar parameters.
             */
            t2Fit_updateParameters(psr,ipsr,psr_parameterEstimates,psr_errorEstimates);
            logtchk("complete updating the parameter values");
            logdbg("Completed updating the parameters");
            free(psr_parameterEstimates);
            free(psr_errorEstimates);
        }
        free(parameterEstimates);
        free(errorEstimates);
        free(white_y);
        free(y);
        free_blas(designMatrix);
        free_blas(white_designMatrix);
        free_blas(constraintsMatrix);


    }

    t2fit_postfit(psr,npsr);
}

/**
 * Some routines need to do stuff before the fit. They go in here.
 */
void t2fit_prefit(pulsar* psr, int npsr){
    for (int ipsr=0; ipsr < npsr; ++ipsr){
        /**
         * Temponest parameters write random crap to the observation struct.
         * So we should clear them first or it will just do wierd stuff.
         */
        if (psr->TNRedAmp && psr->TNRedGam) {
            for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
                if(psr[0].TNsubtractRed==0)psr[ipsr].obsn[iobs].TNRedSignal = 0;
                psr[ipsr].obsn[iobs].TNRedErr = 0;
            }
        }
        if (psr->TNDMAmp && psr->TNDMGam) {
            for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
                psr[ipsr].obsn[iobs].TNDMErr = 0;
            }
        }
    }


}

/**
 * Some fit stuff might need to do some things after the fit. They go in here.
 */
void t2fit_postfit(pulsar* psr, int npsr){

    /**
     * Need to squareroot the TN error bars.
     */
    for (int ipsr = 0; ipsr < npsr; ++ipsr){
        if (psr->TNRedAmp && psr->TNRedGam) {
            for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
                psr[ipsr].obsn[iobs].TNRedErr = sqrt(psr[ipsr].obsn[iobs].TNRedErr);
            }
        }
        if (psr->TNDMAmp && psr->TNDMGam) {
            for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
                psr[ipsr].obsn[iobs].TNDMErr = sqrt(psr[ipsr].obsn[iobs].TNDMErr);
            }
        }
    }

}



unsigned int t2Fit_getFitData(pulsar *psr, double* x, double* y,
        double* e, int* ip){

    int iobs;
    unsigned int ndata=0;

    /*
     * in tempo2, if START and/or FINISH are marked as to fit
     * i.e. fitFlag=1, then it means that we want to only fit
     * over the time between START and FINISH (or just one bound
     * if only one is set).
     *
     * The bools below deal with that.
     */
    bool startSet = psr->param[param_start].paramSet[0]==1 
        && psr->param[param_start].fitFlag[0]==1;
    bool finishSet = psr->param[param_finish].paramSet[0]==1 
        && psr->param[param_finish].fitFlag[0]==1;

    /*
     * Variables for the start/finish. Initialise to some crazy values.
     * 1e10 in MJD is about 27 million years in the future, so not likey to
     * happen (call it the Y27000000 bug).
     */
    longdouble start = 1e10;
    longdouble finish = 0;

    // if we are fixing start/finish then use the specified values.
    if (startSet) start = psr->param[param_start].val[0];
    if (finishSet) finish = psr->param[param_finish].val[0];

    for (iobs=0; iobs < psr->nobs; ++iobs){
        // a convinience pointer for the current observation.
        observation *o = psr->obsn+iobs;
        // copy the current residual to the prefit.
        o->prefitResidual = o->residual;

        // skip deleted points
        if (o->deleted) continue;

        // if start/finish is set, skip points outside of the range
        if (startSet && o->sat < start) continue;
        if (finishSet && o->sat > finish) continue;

        // update start/finish if it isn't set.
        if (!startSet && o->sat < start) start=o->sat;
        if (!finishSet && o->sat > finish) finish=o->sat;

        x[ndata] = (double)(o->bbat - psr->param[param_pepoch].val[0]);
        y[ndata] = o->residual;
        ip[ndata] = iobs;         // index
        e[ndata] = o->toaErr*1e-6; // convert error to seconds.
        ++ndata;
    }

    // save the start/finish values.
    psr->param[param_start].val[0] = start-0.001;
    psr->param[param_finish].val[0] = finish+0.001;
    psr->param[param_start].paramSet[0] = 1;
    psr->param[param_finish].paramSet[0] = 1;

    return ndata;
}

void t2Fit_buildDesignMatrix(pulsar* psr,int ipsr,double x, int ipos, double* afunc){
    // a convinience pointer to the fitinfo struct.
    FitInfo* fitinfo = &(psr[ipsr].fitinfo);
    for (unsigned int ifit = 0; ifit < fitinfo->nParams; ifit++){
        param_label param = fitinfo->paramIndex[ifit];
        // get the fit function
        paramDerivFunc func = fitinfo->paramDerivs[ifit];
        // call the function allocated to this fit parameter
        afunc[ifit] = func(psr,ipsr,x,ipos,param,fitinfo->paramCounters[ifit]);
    }
}

void t2Fit_buildConstraintsMatrix(pulsar* psr,int ipsr, int iconstraint, double* afunc){
    // a convinience pointer to the fitinfo struct.
    FitInfo* fitinfo = &(psr[ipsr].fitinfo);
    constraint_label c_label = fitinfo->constraintIndex[iconstraint];
    constraintDerivFunc func = fitinfo->constraintDerivs[iconstraint];
    for (unsigned int ifit = 0; ifit < fitinfo->nParams; ifit++){
        param_label p_label= fitinfo->paramIndex[ifit];
        // call the function allocated to this constraint
        afunc[ifit] = func(psr,ipsr,c_label,p_label,fitinfo->constraintCounters[iconstraint],fitinfo->paramCounters[ifit]);
    }
}


void t2Fit_updateParameters(pulsar *psr,int ipsr,double *val,double *error){
    // a convinience pointer to the fitinfo struct.
    FitInfo* fitinfo = &(psr[ipsr].fitinfo);
    for (unsigned int ifit = 0; ifit < fitinfo->nParams; ifit++){
        param_label param = fitinfo->paramIndex[ifit];
        // call the function allocated to this fit parameter
        // void updateParameterFunction(pulsar* psr, int ipsr, param_label param, int subparamid, double param, double err);
        paramUpdateFunc func = fitinfo->updateFunctions[ifit];
        func(psr,ipsr,param,fitinfo->paramCounters[ifit],val[ifit],error[ifit]);
    }
}

void t2Fit_fillFitInfo_INNER(pulsar* psr, FitInfo &OUT, const int globalflag);

void t2Fit_fillGlobalFitInfo(pulsar* psr, unsigned int npsr,FitInfo &OUT){
    OUT.nParams=0;
    OUT.nConstraints=0;
    for (int i=1;i<=psr->nJumps;i++) {
        if (psr->fitJump[i]==2)
        {
            logdbg("Global Jump");
            OUT.paramIndex[OUT.nParams]=param_JUMP;
            OUT.paramCounters[OUT.nParams]=i;
            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_jump;
            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_jump;
            ++OUT.nParams;
        }
    }

    // Use pulsar zero to define the global parameters.
    t2Fit_fillFitInfo_INNER(psr,OUT,2);
}

void t2Fit_fillFitInfo(pulsar* psr, FitInfo &OUT, const FitInfo &globals, const double* psr_x, const int* psr_toaidx, const int psr_ndata){

    OUT.nParams=1;
    OUT.nConstraints=0;
    // Add the zero offset
    OUT.paramCounters[0]=0;
    OUT.paramDerivs[0]=t2FitFunc_zero;
    OUT.updateFunctions[0]=t2UpdateFunc_zero;
    OUT.paramIndex[0]=param_ZERO;

    /*
     * these are the "classical" tempo2 constraints.
     */
    for (int i=0;i<psr->nconstraints;++i) {
        OUT.constraintIndex[OUT.nConstraints]=psr->constraints[i];
        OUT.constraintCounters[OUT.nConstraints]=0;
        switch(psr->constraints[i]){
            default:
                // this is a quick fix to avoid re-writing code.
                OUT.constraintDerivs[OUT.nConstraints] = standardConstraintFunctions;
                ++OUT.nConstraints;
                break;
        }
    }


    /*
     * These are the temponest derived constraints
     * @TODO: Probably want to extract this as a method
     */
    {
        if (psr->TNRedAmp && psr->TNRedGam) {
            for (int i=0;i<psr->TNRedC;++i) {
                /**
                 * Temporary fix here!
                 * Enable the TN parameters.
                 */
                psr->param[param_red_sin].fitFlag[0]=1;
                psr->param[param_red_cos].fitFlag[0]=1;

                psr->param[param_red_sin].paramSet[0]=1;
                psr->param[param_red_cos].paramSet[0]=1;

                // End of temporary fix.

                OUT.constraintIndex[OUT.nConstraints]=constraint_red_sin;
                OUT.constraintCounters[OUT.nConstraints]=i;
                OUT.constraintDerivs[OUT.nConstraints] = constraints_nestlike_red;
                ++OUT.nConstraints;
                OUT.constraintIndex[OUT.nConstraints]=constraint_red_cos;
                OUT.constraintCounters[OUT.nConstraints]=i;
                OUT.constraintDerivs[OUT.nConstraints] = constraints_nestlike_red;
                ++OUT.nConstraints;
            }
        }
        if (psr->nTNBandNoise > 0){
            int counter=0;
            for (int iband ; iband < psr->nTNBandNoise; ++iband){
                for (int i=0;i<psr->TNBandNoiseC[iband];++i) {
                    /**
                     * Temporary fix here!
                     * Enable the TN parameters.
                     */
                    psr->param[param_band_red_sin].fitFlag[0]=1;
                    psr->param[param_band_red_cos].fitFlag[0]=1;

                    psr->param[param_band_red_sin].paramSet[0]=1;
                    psr->param[param_band_red_cos].paramSet[0]=1;

                    // End of temporary fix.

                    OUT.constraintIndex[OUT.nConstraints]=constraint_band_red_sin;
                    OUT.constraintCounters[OUT.nConstraints]=counter;
                    OUT.constraintDerivs[OUT.nConstraints] = constraints_nestlike_band;
                    ++OUT.nConstraints;
                    OUT.constraintIndex[OUT.nConstraints]=constraint_band_red_cos;
                    OUT.constraintCounters[OUT.nConstraints]=counter;
                    OUT.constraintDerivs[OUT.nConstraints] = constraints_nestlike_band;
                    ++OUT.nConstraints;
                    ++counter;
                }
            }
        }

        if (psr->TNDMAmp && psr->TNDMGam) {
            for (int i=0;i<psr->TNDMC;++i) {
                /**
                 * Temporary fix here!
                 * Enable the TN parameters.
                 */
                psr->param[param_red_dm_sin].fitFlag[0]=1;
                psr->param[param_red_dm_cos].fitFlag[0]=1;

                psr->param[param_red_dm_sin].paramSet[0]=1;
                psr->param[param_red_dm_cos].paramSet[0]=1;

                // End of temporary fix.

                OUT.constraintIndex[OUT.nConstraints]=constraint_red_dm_sin;
                OUT.constraintCounters[OUT.nConstraints]=i;
                OUT.constraintDerivs[OUT.nConstraints] = constraints_nestlike_red_dm;
                ++OUT.nConstraints;
                OUT.constraintIndex[OUT.nConstraints]=constraint_red_dm_cos;
                OUT.constraintCounters[OUT.nConstraints]=i;
                OUT.constraintDerivs[OUT.nConstraints] = constraints_nestlike_red_dm;
                ++OUT.nConstraints;
            }
        }

        // ECORR parameters
        // Fits for a constrained jump for each "epoch" (i.e. observation)
        int nepochs=0;
        if (psr->nTNECORR > 0) {

            const double dt = 10.0/SECDAY;
            double xmin=0;
            double xmax=0;
            // flag observations that have been matched to an epoch already.
            bool flag[psr_ndata];
            for (int idata = 0; idata < psr_ndata; ++idata){
                flag[idata]= false;
            }

            for (int idata = 0; idata < psr_ndata; ++idata){
                int iobs = psr_toaidx[idata];
                // if this observation was in a previous epoch, then skip it.
                if(flag[idata])continue;
                // Here we check over all flags and ecorr params to see if it matches.
                for (int iecorr=0; iecorr < psr->nTNECORR; iecorr++){
                    for (int iflag=0;iflag < psr->obsn[iobs].nFlags; iflag++){
                        if (
                                (strcmp(psr->obsn[iobs].flagID[iflag],
                                        psr->TNECORRFlagID[iecorr])==0)
                                && (strcmp(psr->obsn[iobs].flagVal[iflag],
                                        psr->TNECORRFlagVal[iecorr])==0)

                           ) {
                            xmin=psr_x[idata] - dt;
                            xmax=psr_x[idata] + dt;
                            // don't need to search any more on this epoch.
                            flag[idata]=true;

                            OUT.constraintIndex[OUT.nConstraints]=constraint_jitter;
                            // we use the constraint counter to specify specific toa that represents this epoch
                            OUT.constraintCounters[OUT.nConstraints]=iobs;
                            OUT.constraintDerivs[OUT.nConstraints] = constraints_nestlike_jitter;
                            ++OUT.nConstraints;

                            OUT.paramIndex[OUT.nParams]=param_jitter;
                            // we use the constraint counter to specify specific toa that represents this epoch
                            OUT.paramCounters[OUT.nParams]=iobs;
                            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_nestlike_jitter;
                            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_nestlike_jitter;
                            ++OUT.nParams;

                            ++nepochs;
                            break;
                        }
                        if(flag[idata])break;
                    }
                    if(flag[idata])break;
                }
                if (flag[idata]){
                    for (int idata2 = idata; idata2 < psr_ndata; ++idata2){
                        if(flag[idata2])continue;
                        if(psr_x[idata2] > xmin && psr_x[idata2] < xmax){
                            flag[idata2]=true;
                        }
                    }
                }
            }
            logmsg("Generated %d TNepochs from %d TNECORRs",nepochs,psr->nTNECORR);
        }
    }


    for (int i=1;i<=psr->nJumps;i++) {
        if (psr->fitJump[i]==1)
        {
            OUT.paramIndex[OUT.nParams]=param_JUMP;
            OUT.paramCounters[OUT.nParams]=i;
            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_jump;
            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_jump;
            ++OUT.nParams;
        }
    }
    t2Fit_fillFitInfo_INNER(psr,OUT,1);

    // copy over global parameters.
    for (size_t i=0;i<globals.nParams;i++) {
        OUT.paramIndex[OUT.nParams]=globals.paramIndex[i];
        OUT.paramCounters[OUT.nParams]=globals.paramCounters[i];
        OUT.paramDerivs[OUT.nParams]     =globals.paramDerivs[i];
        OUT.updateFunctions[OUT.nParams] =globals.updateFunctions[i];
        ++OUT.nParams;
    }

    if (psr->TNRedAmp && psr->TNRedGam) {
        /**
         * Temporary fix here!
         * Disable the TN parameters.
         */
        psr->param[param_red_sin].fitFlag[0]=0;
        psr->param[param_red_cos].fitFlag[0]=0;
        psr->param[param_red_sin].paramSet[0]=0;
        psr->param[param_red_cos].paramSet[0]=0;
    }
    if (psr->nTNBandNoise > 0 ) {
        /**
         * Temporary fix here!
         * Disable the TN parameters.
         */
        psr->param[param_band_red_sin].fitFlag[0]=0;
        psr->param[param_band_red_cos].fitFlag[0]=0;
        psr->param[param_band_red_sin].paramSet[0]=0;
        psr->param[param_band_red_cos].paramSet[0]=0;
    }

    if (psr->TNDMAmp && psr->TNDMGam) {
        /**
         * Temporary fix here!
         * Disable the TN parameters.
         */
        psr->param[param_red_dm_sin].fitFlag[0]=0;
        psr->param[param_red_dm_cos].fitFlag[0]=0;
        psr->param[param_red_dm_sin].paramSet[0]=0;
        psr->param[param_red_dm_cos].paramSet[0]=0;
    }


}


void t2fit_fillOneParameterFitInfo(pulsar* psr, param_label fit_param, const int k, FitInfo& OUT){

    unsigned N;
    bool sifunc;
    OUT.paramIndex[OUT.nParams]=fit_param;
    OUT.paramCounters[OUT.nParams]=k;

    switch(fit_param){
        case param_raj:
        case param_decj:
        case param_pmra:
        case param_pmdec:
        case param_px:
        case param_pmrv:
        case param_dshk:
            // positional parameters and parallax
            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_stdPosition;
            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_stdPosition;
            ++OUT.nParams;
            break;
        case param_start:
        case param_finish:
            // these parameters are not actually fitted for...
            break;
        case param_f:
            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_stdFreq;
            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_stdFreq;
            ++OUT.nParams;
            break;
        case param_brake:
            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_stdFreq;
            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_simpleMinus;
            ++OUT.nParams;
            break;
        case param_sini:
        case param_pb:
        case param_fb:
        case param_t0:
        case param_a1:
        case param_om:
        case param_ecc:
        case param_edot:
        case param_e2dot:
        case param_xpbdot:
        case param_pbdot:
        case param_a1dot:
        case param_a2dot:
        case param_omdot:
        case param_orbpx:
        case param_tasc:
        case param_eps1:
        case param_eps2:
        case param_eps1dot:
        case param_eps2dot:
        case param_m2:
        case param_gamma:
        case param_mtot:
        case param_bp:
        case param_bpp:
        case param_dr:
        case param_dtheta:
        case param_bpjep:
        case param_bpjph:
        case param_bpja1:
        case param_bpjec:
        case param_bpjom:
        case param_bpjpb:
        case param_h3:
        case param_h4:
        case param_stig:
        case param_kin:
        case param_kom:
            // all binary models are handled by a routine that identifies the correct binary model.
            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_binaryModels;
            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_binaryModels;
            ++OUT.nParams;
            break;
        case param_dm:
            // Dispersion measure and derivatives
            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_stdDm;
            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_simpleAdd;
            ++OUT.nParams;
            break;
        case param_fddc:
            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_fddc;
            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_simpleAdd;
            ++OUT.nParams;
            break;
        case param_fd:
            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_fd;
            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_simpleAdd;
            ++OUT.nParams;
            break;
        case param_dm_sin1yr:
        case param_dm_cos1yr:
            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_dmsinusoids;
            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_simpleAdd;
            ++OUT.nParams;
            break;
        case param_fddi:
            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_notImplemented;
            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_notImplemented;
            ++OUT.nParams;
            break;
        case param_dmx:
            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_dmx;
            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_simpleAdd;
            ++OUT.nParams;
            break;
        case param_dmassplanet:
            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_planet;
            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_planet;
            ++OUT.nParams;
            break;
        case param_wave_om:
            // fitwaves has many parameters to fit.
            N=2;
            if (psr->waveScale == 2)N=4;
            for (unsigned i = 0; i < static_cast<unsigned>(psr->nWhite*N); ++i){
                OUT.paramDerivs[OUT.nParams]     =t2FitFunc_fitwaves;
                OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_fitwaves;
                OUT.paramCounters[OUT.nParams]=i;
                OUT.paramIndex[OUT.nParams]=fit_param;
                ++OUT.nParams;
            }
            break;
        case param_wave_dm:
            // fitwaves has many parameters to fit.
            for (unsigned i = 0; i < static_cast<unsigned>(psr->nWhite_dm*2); ++i){
                OUT.paramDerivs[OUT.nParams]     =t2FitFunc_fitwaves;
                OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_fitwaves;
                OUT.paramCounters[OUT.nParams]=i;
                OUT.paramIndex[OUT.nParams]=fit_param;
                ++OUT.nParams;
            }
            break;

        case param_glep:
        case param_glph:
        case param_glf0:
        case param_glf1:
        case param_glf0d:
        case param_gltd:
        case param_glf2:
            // glitches
            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_stdGlitch;
            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_simpleMinus;
            if(fit_param==param_glf2)OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_stdGlitch;
            ++OUT.nParams;
            break;

        case param_telx:
        case param_tely:
        case param_telz:
            OUT.paramDerivs[OUT.nParams]     =t2FitFunc_telPos;
            OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_telPos;
            ++OUT.nParams;
            break;
        case param_tel_dx:
        case param_tel_dy:
        case param_tel_dz:
            // complicated to work out how many of these parameters there are
            // (blame George??)
            if(fit_param==param_tel_dx)N=psr->nTelDX;
            if(fit_param==param_tel_dy)N=psr->nTelDY;
            if(fit_param==param_tel_dz)N=psr->nTelDZ;
            if(psr->param[fit_param].val[0] == -1)N=0;
            else if (psr->param[fit_param].val[0] < 2 )N=N;
            else N=N-1;
            for (unsigned i = 0; i < N; ++i){
                OUT.paramDerivs[OUT.nParams]     =t2FitFunc_telPos;
                OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_telPos;
                OUT.paramIndex[OUT.nParams]=fit_param;
                ++OUT.nParams;
            }
            break;

        case param_ifunc:
        case param_clk_offs:
        case param_quad_ifunc_p:
        case param_quad_ifunc_c:
            // ifunc-alikes
            N=0;
            sifunc=psr->param[fit_param].val[0]==2; // use sinusoids?
            if(fit_param==param_ifunc){
                N=psr->ifuncN;
                sifunc=!sifunc;
            }
            if(fit_param==param_clk_offs)N=psr->clkOffsN;
            if(fit_param==param_quad_ifunc_p)N=psr->quad_ifuncN_p;
            if(fit_param==param_quad_ifunc_c)N=psr->quad_ifuncN_c;
            if(sifunc){
                logwarn("Sinusoidal ifuncs");
            }
            for (unsigned i = 0; i < N; ++i){
                if(sifunc) OUT.paramDerivs[OUT.nParams] = t2FitFunc_sifunc;
                else OUT.paramDerivs[OUT.nParams] = t2FitFunc_ifunc;
                OUT.updateFunctions[OUT.nParams] = t2UpdateFunc_ifunc;
                OUT.paramCounters[OUT.nParams]=i;
                OUT.paramIndex[OUT.nParams]=fit_param;
                ++OUT.nParams;
            }
            break;

        case param_gwsingle:
            for (unsigned i = 0; i < 4; ++i){
                OUT.paramDerivs[OUT.nParams]     =t2FitFunc_stdGravWav;
                OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_stdGravWav;
                OUT.paramCounters[OUT.nParams]=i;
                OUT.paramIndex[OUT.nParams]=fit_param;
                ++OUT.nParams;
            }
            break;
        case param_dmmodel:
            for (int i = 0; i < psr->dmoffsDMnum; ++i){
                OUT.paramDerivs[OUT.nParams]     =t2FitFunc_dmmodelDM;
                OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_dmmodelDM;
                OUT.paramCounters[OUT.nParams]=i;
                OUT.paramIndex[OUT.nParams]=fit_param;
                ++OUT.nParams;
            }
            for (int i = 0; i < psr->dmoffsCMnum; ++i){
                OUT.paramDerivs[OUT.nParams]     =t2FitFunc_dmmodelCM;
                OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_dmmodelCM;
                OUT.paramCounters[OUT.nParams]=i+psr->dmoffsDMnum;
                OUT.paramIndex[OUT.nParams]=fit_param;
                ++OUT.nParams;
            }
            break;

        case param_red_sin:
        case param_red_cos:
            for (int i=0; i < psr->TNRedC ; ++i){
                OUT.paramDerivs[OUT.nParams]     =t2FitFunc_nestlike_red;
                OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_nestlike_red;
                OUT.paramCounters[OUT.nParams]=i;
                OUT.paramIndex[OUT.nParams]=fit_param;
                ++OUT.nParams;
            }
            break;

        case param_band_red_sin:
        case param_band_red_cos:
            {
                int counter=0;
                for (int iband=0; iband < psr->nTNBandNoise ; ++iband){
                    for (int i=0; i < psr->TNBandNoiseC[iband] ; ++i){
                        OUT.paramDerivs[OUT.nParams]     =t2FitFunc_nestlike_band;
                        OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_nestlike_band;
                        OUT.paramCounters[OUT.nParams]=counter;
                        OUT.paramIndex[OUT.nParams]=fit_param;
                        ++OUT.nParams;
                        ++counter;
                    }
                }
            }
            break;


        case param_red_dm_sin:
        case param_red_dm_cos:
            {
                for (int i=0; i < psr->TNDMC ; ++i){
                    OUT.paramDerivs[OUT.nParams]     =t2FitFunc_nestlike_red_dm;
                    OUT.updateFunctions[OUT.nParams] =t2UpdateFunc_nestlike_red_dm;
                    OUT.paramCounters[OUT.nParams]=i;
                    OUT.paramIndex[OUT.nParams]=fit_param;
                    ++OUT.nParams;
                }

            }
            break;
        default:
            logerr("ERROR: No methods for fitting parameter %s (%d)",label_str[fit_param],fit_param);
            break;
    }
}

paramDerivFunc getDerivFunction(pulsar* psr, param_label fit_param, const int k){
    FitInfo nfo;
    nfo.nParams=0;
    nfo.nConstraints=0;

    t2fit_fillOneParameterFitInfo(psr,fit_param,k,nfo);

    unsigned ip;
    for(ip = 0; ip < nfo.nParams; ++ip){
        if (nfo.paramCounters[ip]==k) break;
    }
    return nfo.paramDerivs[ip];
}


void t2Fit_fillFitInfo_INNER(pulsar* psr, FitInfo &OUT, const int globalflag){
    for (param_label fit_param=0; fit_param < param_LAST; ++fit_param){
        if (fit_param == param_ifunc){
            logmsg("if: %d %d %d",psr->param[fit_param].paramSet[0],psr->param[fit_param].fitFlag[0],psr->param[fit_param].aSize);
        }
        for(int k=0; k < psr->param[fit_param].aSize;k++){
            if (psr->param[fit_param].paramSet[k]>0 && psr->param[fit_param].fitFlag[k]==globalflag) {
                t2fit_fillOneParameterFitInfo(psr,fit_param,k,OUT);
            }
        }
    }
}


double t2Fit_getParamDeriv(pulsar* psr, const int i, const double x, const param_label fit_param, const int k){
    paramDerivFunc f = getDerivFunction(psr,fit_param,k);
    return f(psr,0,x,i,fit_param,k);
}












//###### LEGACY ROUTINES

// routine for pulsar fitting.
void TKleastSquares_single_pulsar(double *x,double *y,int n,double *outP,double *e,int nf,double **cvm, double *chisq, void (*fitFuncs)(double, double [], int,pulsar *,int, int),pulsar *psr,double tol, int *ip,char rescale_errors, double **uinv) {

    double **designMatrix, **white_designMatrix, **constraintsMatrix;
    double *b,*white_b;
    constraintsMatrix=NULL;

    TKfit_getPulsarDesignMatrix(x,y,n,nf,fitFuncs,psr,ip,uinv,0,&designMatrix,&white_designMatrix,&b,&white_b);
    if(psr->nconstraints > 0){
        logmsg("Get constraint weights");
        computeConstraintWeights(psr);
        logmsg("fill constraints matrix");
        constraintsMatrix = malloc_blas(psr->nconstraints,nf);
        for (int ic=0; ic < psr->nconstraints; ic++){
            CONSTRAINTfuncs(psr,0,nf,psr->constraints[ic],constraintsMatrix[ic]);
        }
    }

    *chisq = TKrobustConstrainedLeastSquares(b,white_b,designMatrix,white_designMatrix,constraintsMatrix,
            n,nf,psr->nconstraints,tol,rescale_errors,
            outP,e,cvm,psr->robust);
    free_blas(designMatrix); // free-TKleastSquares_svd_psr_dcm-designMatrix**
    free_blas(white_designMatrix);  // free-TKleastSquares_svd_psr_dcm-white_designMatrix**
    if(psr->nconstraints > 0) free_blas(constraintsMatrix);
    free(b);
    free(white_b);

}

void TKleastSquares_global_pulsar(double **x,double **y,int *n,
        double *outP,double *e,int* nf, int nglobal,double **cvm, double *chisq, void (*fitFuncs)(double, double [], int,pulsar *,int,int),pulsar *psr,double tol, int **ip,char rescale_errors, double ***uinv, int npsr) {

    double **designMatrix, **white_designMatrix;
    double **constraintsMatrix;
    double **psr_DM, **psr_wDM;
    double *b,*white_b, *psr_b,*psr_wb;
    int ipsr;
    int totalFit=0;
    int totalObs=0;
    int totalConstraints=0;
    int i,j;
    int off_r=0;
    int off_f=0;
    int off_c=0;

    for (ipsr=0; ipsr < npsr; ipsr++){
        totalFit+=nf[ipsr];
        totalObs+=n[ipsr];
        totalConstraints+=psr[ipsr].nconstraints;
    }
    totalFit+=nglobal;

    white_designMatrix=malloc_blas(totalObs,totalFit);
    designMatrix=malloc_blas(totalObs,totalFit);
    constraintsMatrix=malloc_blas(totalConstraints,totalFit);
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

        if(psr[ipsr].nconstraints > 0){
            logmsg("Get constraint weights");
            computeConstraintWeights(psr+ipsr);
            logmsg("fill constraints matrix");
            for (int ic=0; ic < psr[ipsr].nconstraints; ic++){
                CONSTRAINTfuncs(psr,ipsr,nf[ipsr],psr->constraints[ic],constraintsMatrix[ic+off_c]+off_f);
            }
        }



        // increment the offset.
        off_r += n[ipsr];
        off_f += nf[ipsr];
        off_c += psr[ipsr].nconstraints;

        // free temp matricies.
        free_blas(psr_DM);
        free_blas(psr_wDM);
        free(psr_b);
        free(psr_wb);
    }


    // go ahead and do the fit!

    *chisq = TKrobustConstrainedLeastSquares(b,white_b,designMatrix,white_designMatrix,
            constraintsMatrix,
            totalObs,totalFit,totalConstraints,tol,rescale_errors,
            outP,e,cvm,psr[0].robust);

    free_blas(designMatrix); // free-TKleastSquares_svd_psr_dcm-designMatrix**
    free_blas(white_designMatrix);  // free-TKleastSquares_svd_psr_dcm-white_designMatrix**
    free_blas(constraintsMatrix);  // free-TKleastSquares_svd_psr_dcm-white_designMatrix**
    free(b);
    free(white_b);

}



void TKfit_getPulsarDesignMatrix(double *x,double *y,int n,int nf,void (*fitFuncs)(double, double [], int,pulsar *,int,int), pulsar *psr, int* ip, double **uinv,int ipsr,double ***OUT_designMatrix,double ***OUT_white_designMatrix,double** OUT_b, double** OUT_wb){

    //double precision arrays for matrix algebra.
    double **designMatrix, **white_designMatrix;
    double basisFunc[nf];
    double *b,*white_b;
    int    i,j;
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


void TKleastSquares_svd_passN(double *x,double *y,double *sig2,int n,double *p,double *e,int nf,double **cvm, double *chisq, void (*fitFuncs)(double, double [], int,int),int weight)
{
    logerr("This method no longer implemented.");
    exit(-1);
}

// END OF LEGACY ROUTINES


