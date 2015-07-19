#include "t2fit.h"
#include <stdlib.h>
#include <string.h>


void t2Fit(pulsar *psr,unsigned int npsr, char *covarFuncFile){
    unsigned int ipsr;
    const double tol=1e-27;
    bool haveCovar = (covarFuncFile!=NULL && strcmp(covarFuncFile,"NULL"));
    // "uinv" is the whitening matrix

    // get global fit params/constraints.
    FitInfo global_fitinfo;
    t2Fit_fillGlobalFitInfo(psr,npsr,global_fitinfo);
    bool doGlobalFit = global_fitinfo.nParams | global_fitinfo.nConstraints;

    unsigned long long totalGlobalData=0;
    unsigned int gParams=global_fitinfo.nParams;
    unsigned int gConstraints=global_fitinfo.nConstraints;
    double** gUinvs[MAX_PSR];
    double* gX[MAX_PSR];
    double* gY[MAX_PSR];

    for (ipsr=0; ipsr < npsr; ipsr++) {
        double *psr_x   = (double*)malloc(sizeof(double)*psr[ipsr].nobs);
        double *psr_y   = (double*)malloc(sizeof(double)*psr[ipsr].nobs);
        double *psr_e   = (double*)malloc(sizeof(double)*psr[ipsr].nobs);
        int *psr_toaidx = (int*)malloc(sizeof(int)*psr[ipsr].nobs);
        double** uinv;

        // get the fit data for this pulsar
        int psr_ndata = t2Fit_getFitData(psr+ipsr,psr_x,psr_y,psr_e,psr_toaidx);

        // get the params/constraints etc for this fit
        t2Fit_fillFitInfo(psr+ipsr,psr[ipsr].fitinfo);

        // fill the whitening matrix
        if (haveCovar) {
            // ToAs must be sorted for covariance function code
            sortToAs(psr+ipsr);

            uinv = malloc_uinv(psr_ndata);
            psr[ipsr].fitMode=1; // Note: forcing this to 1 as the Cholesky fit is a weighted fit
            logmsg("Doing a FULL COVARIANCE MATRIX fit");
        } else {
            // Here the whitening matrix is just a diagonal
            // weighting matrix. Store diagonal matrix as 1xN
            // so that types match later.
            uinv=malloc_blas(1,psr_ndata); 
            if(psr[ipsr].fitMode == 0){
                logdbg("Doing an UNWEIGHTED fit");
                // unweighted fit - set sigma to 1
                for (int i=0; i < psr_ndata; i++){
                    psr_e[i]=1.0;
                }
            } else {
                logdbg("Doing a WEIGHTED fit");
            }
        }

        getCholeskyMatrix(uinv,covarFuncFile,psr+ipsr,
                psr_x,psr_y,psr_e,
                psr_ndata,0,psr_toaidx);
        psr[ipsr].nFit = psr_ndata;

        if (doGlobalFit){
            // we are going to do a global fit, so need to store the values for later
            gX[ipsr] = psr_x;
            gY[ipsr] = psr_y;
            totalGlobalData += psr_ndata;
            // MORE HERE.
        } else {
            // we do one fit at a time...
            double chisq;
            double* parameterEstimates = (double*)malloc(sizeof(double)*psr[ipsr].fitinfo.nParams);
            double* errorEstimates = (double*)malloc(sizeof(double)*psr[ipsr].fitinfo.nParams);
            // run the fit.
            TKleastSquares_single_pulsar(psr_x,psr_y,psr_ndata,
                    parameterEstimates,errorEstimates,
                    psr[ipsr].fitinfo.nParams,
                    psr[ipsr].covar,&chisq,
                    t2Fit_buildDesignMatrix,psr+ipsr,tol,psr_toaidx,
                    1,uinv);

            psr[ipsr].fitChisq = chisq;
                psr[ipsr].fitNfree = psr_ndata
                + psr[ipsr].fitinfo.nConstraints
                - psr[ipsr].fitinfo.nParams;
            logdbg("Updating the parameters");
            logtchk("updating the parameter values");
            t2Fit_updateParameters(psr,ipsr,parameterEstimates,errorEstimates);
            logtchk("complete updating the parameter values");
            logdbg("Completed updating the parameters");
            logdbg("Free fit memory");
            free(parameterEstimates);
            free(errorEstimates);
            free_blas(uinv);
            free(psr_x);
            free(psr_y);
            free(psr_e);
            free(psr_toaidx);
        }
    }
}

unsigned int t2Fit_getFitData(pulsar *psr, double* x, double* y,
        double* e, int* ip){

    unsigned int iobs;
    unsigned int ndata=0;

    bool startSet = psr->param[param_start].paramSet[0]==1 
        && psr->param[param_start].fitFlag[0]==1;
    bool finishSet = psr->param[param_finish].paramSet[0]==1 
        && psr->param[param_finish].fitFlag[0]==1;
    longdouble start = 1e10;
    longdouble finish = 0;
    if (startSet) start = psr->param[param_start].val[0];
    if (finishSet) finish = psr->param[param_finish].val[0];

    for (iobs=0; iobs < psr->nobs; ++iobs){
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

    // save the new start/finish values.
    psr->param[param_start].val[0] = start;
    psr->param[param_finish].val[0] = finish;
    psr->param[param_start].paramSet[0] = 1;
    psr->param[param_finish].paramSet[0] = 1;

    return ndata;
}

void t2Fit_buildDesignMatrix(double x,double afunc[],int ma,pulsar *psr,int ipos,int ipsr){
    FitInfo* fitinfo = &(psr[ipsr].fitinfo);
    for (int ifit = 0; ifit < fitinfo->nParams; ifit++){
        param_label param = fitinfo->paramIndex[ifit];
        // call the function allocated to this fit parameter
        // double paramDerivFunc(pulsar* psr, int ipsr, double x, int obsnid, param_label label, int subparamid)
        paramDerivFunc func = fitinfo->paramDerivs[ifit];
        afunc[ifit] = func(psr,ipsr,x,ipos,param,fitinfo->paramCounters[ifit]);
    }
}

void t2Fit_buildConstraintsMatrix(pulsar* psr,int ipsr, int iconstraint, double* afunc){
    FitInfo* fitinfo = &(psr[ipsr].fitinfo);
    constraint_label c_label = fitinfo->constraintIndex[iconstraint];
    constraintDerivFunc func = fitinfo->constraintDerivs[iconstraint];
    for (int ifit = 0; ifit < fitinfo->nParams; ifit++){
        param_label p_label= fitinfo->paramIndex[ifit];
        // call the function allocated to this constraint
        afunc[ifit] = func(psr,ipsr,c_label,p_label,fitinfo->constraintCounters[iconstraint],fitinfo->paramCounters[ifit]);
    }
}


void t2Fit_updateParameters(pulsar *psr,int ipsr,double *val,double *error){
    FitInfo* fitinfo = &(psr[ipsr].fitinfo);
    for (int ifit = 0; ifit < fitinfo->nParams; ifit++){
        param_label param = fitinfo->paramIndex[ifit];
        // call the function allocated to this fit parameter
        // void updateParameterFunction(pulsar* psr, int ipsr, param_label param, int subparamid, double param, double err);
        paramUpdateFunc func = fitinfo->updateFunctions[ifit];
        func(psr,ipsr,param,fitinfo->paramCounters[ifit],val[ifit],error[ifit]);
    }
}

void t2Fit_fillGlobalFitInfo(pulsar* psr, unsigned int npsr,FitInfo &OUT){
    OUT.nParams=0;
    OUT.nConstraints=0;
}

void t2Fit_fillFitInfo(pulsar* psr, FitInfo &OUT){
    OUT.nParams=0;
    OUT.nConstraints=0;
}

