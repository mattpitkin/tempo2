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
    double* gW[MAX_PSR];
    double** gDM[MAX_PSR];
    double** gWDM[MAX_PSR];
    unsigned int* gNdata[MAX_PSR];

    for (ipsr=0; ipsr < npsr; ipsr++) {
        double *psr_x   = (double*)malloc(sizeof(double)*psr[ipsr].nobs);
        double *psr_y   = (double*)malloc(sizeof(double)*psr[ipsr].nobs);
        double *psr_white_y   = (double*)malloc(sizeof(double)*psr[ipsr].nobs);
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

        unsigned nParams=psr[ipsr].fitinfo.nParams;
        double** constraintsMatrix =NULL;
        if(psr[ipsr].fitinfo.nConstraints > 0){
            constraintsMatrix = malloc_blas(ndata,psr[ipsr].fitinfo.nConstraints);
            for (unsigned int iconstraint =0; idata < psr_ndata; ++idata){
                t2Fit_buildConstraintsMatrix(psr, ipsr, iconstraint, constraintsMatrix[iconstraint]);
            }
        }

        double** designMatrix = malloc_blas(ndata,nParams);
        double** white_designMatrix = malloc_blas(ndata,nParams);
        for (unsigned int idata =0; idata < psr_ndata; ++idata){
            t2Fit_buildDesignMatrix(psr,ipsr, psr_toaidx[idata],psr_x[idata], designMatrix[idata]);
        }

        // do the "pre whitening"
        if(haveCovar){
            TKmultMatrixVec(uinv,psr_y,psr_ndata,psr_ndata,psr_white_y);
            TKmultMatrix_sq(uinv,designMatrix,psr_ndata,nParams,white_designMatrix);
        } else {
            for(unsigned i=0;i<psr_ndata;++i){
                psr_white_y[i]=psr_y[i]*uinv[0][i];
                for(unsigned j=0;j<nParam;++j){
                    white_designMatrix[i] = designMatrix[i][j]*uinv[0][i];
                }
            }
        }

        if (doGlobalFit){
            // we are going to do a global fit, so need to store the values for later
            gX[ipsr] = psr_x;
            gY[ipsr] = psr_y;
            gW[ipsr] = psr_white_y;
            gDM[ipsr] = designMatrix;
            gWDM[ipsr] = white_designMatrix;
            gNdata[ipsr] = psr_ndata;
            totalGlobalData += psr_ndata;
            // MORE HERE.
        } else {
            // NOT GLOBAL
            // so do one fit at a time...
            double chisq;
            double* parameterEstimates = (double*)malloc(sizeof(double)*nParams);
            double* errorEstimates = (double*)malloc(sizeof(double)*nParams);

            chisq = TKrobustConstrainedLeastSquares(psr_y,psr_white_y,
                    designMatrix,white_designMatrix,constraintsMatrix,
                    psr_ndata,nParams,psr[ipsr].fitinfo.nConstraints,
                    tol,1,parameterEstimates,errorEstimates,psr[ipsr].covar,
                    psr[ipsr].robust);

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
            free_blas(designMatrix);
            free_blas(white_designMatrix);
            free_blas(uinv);
            free(psr_x);
            free(psr_y);
            free(psr_white_y);
            free(psr_e);
            free(psr_toaidx);
        }
    }
    if (doGlobalFit){
        // add the global fit parameters
        // form the final design matrix
        // do the global fit
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

void t2Fit_buildDesignMatrix(pulsar* psr,int ipsr,double x, int ipos, double* afunc){
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
    unsigned N;
    for (unsigned iparam=0; iparam < MAX_PARAMS; iparam++){
        for(unsigned k=0; k < psr->param[i].aSize;k++){
            if (psr->param[iparam].paramSet[k]==1 && psr->param[iparam].fitFlag[k]==1) {

                OUT.paramCounters[OUT.nparams]=k;
                switch(iparam){
                    case param_raj:
                    case param_decj:
                    case param_pmra:
                    case param_pmdec:
                    case param_px:
                    case param_pmrv:
                    case param_dshk:
                        // positional parameters and parallax
                        OUT.paramDerivs[OUT.nparams]     =t2FitFunc_stdPosition;
                        OUT.updateFunctions[OUT.nparams] =t2UpdateFunc_stdPosition;
                        ++OUT.nparams;
                        break
                    case param_start:
                    case param_finish:
                            // these parameters are not actually fitted for...
                            break;
                    case param_f:
                            OUT.paramDerivs[OUT.nparams]     =t2FitFunc_stdFreq;
                            OUT.updateFunctions[OUT.nparams] =t2UpdateFunc_stdFreq;
                            ++OUT.nparams;
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
                            // all binary models are handled by a routine that identifies the correct binary model.
                            OUT.paramDerivs[OUT.nparams]     =t2FitFunc_binaryModels;
                            OUT.updateFunctions[OUT.nparams] =t2UpdateFunc_binaryModels;
                            ++OUT.nparams;
                            break;
                    case param_dm:
                            // Dispersion measure and derivatives
                            OUT.paramDerivs[OUT.nparams]     =t2FitFunc_stdDm;
                            OUT.updateFunctions[OUT.nparams] =t2UpdateFunc_stdDm;
                            ++OUT.nparams;
                            break;
                    case param_fddc:
                    case param_fddi:
                    case param_fd:
                    case param_dm_sin1yr:
                    case param_dm_cos1yr:
                    case param_dmx:
                            OUT.paramDerivs[OUT.nparams]     =t2FitFunc_miscDm;
                            OUT.updateFunctions[OUT.nparams] =t2UpdateFunc_miscDm;
                            ++OUT.nparams;
                            break;
                    case param_dmassplanet:
                            OUT.paramDerivs[OUT.nparams]     =t2FitFunc_planet;
                            OUT.updateFunctions[OUT.nparams] =t2UpdateFunc_planet;
                            ++OUT.nparams;
                            break;
                    case param_wave_om:
                            // fitwaves has many parameters to fit.
                            for (unsigned i = 0; i < psr->nWhite; ++i){
                                OUT.paramDerivs[OUT.nparams]     =t2FitFunc_fitwaves;
                                OUT.updateFunctions[OUT.nparams] =t2UpdateFunc_fitwaves;
                                OUT.paramCounters[OUT.nparams]=i;
                                ++OUT.nparams;
                            }
                            break;
                    case param_wave_dm:
                            // fitwaves has many parameters to fit.
                            for (unsigned i = 0; i < psr->nWhite_dm; ++i){
                                OUT.paramDerivs[OUT.nparams]     =t2FitFunc_fitwaves;
                                OUT.updateFunctions[OUT.nparams] =t2UpdateFunc_fitwaves;
                                OUT.paramCounters[OUT.nparams]=i;
                                ++OUT.nparams;
                            }
                            break;

                    case param_glep:
                    case param_glph:
                    case param_glf0:
                    case param_glf1:
                    case param_glf2:
                    case param_glf0d:
                    case param_gltd:
                            // glitches
                            OUT.paramDerivs[OUT.nparams]     =t2FitFunc_stdGlitch;
                            OUT.updateFunctions[OUT.nparams] =t2UpdateFunc_stdGlitch;
                            ++OUT.nparams;
                            break;

                    case param_telx:
                    case param_tely:
                    case param_telz:
                            OUT.paramDerivs[OUT.nparams]     =t2FitFunc_telPos;
                            OUT.updateFunctions[OUT.nparams] =t2UpdateFunc_telPos;
                            ++OUT.nparams;
                            break;
                    case param_tel_dx:
                    case param_tel_dy:
                    case param_tel_dz:
                            // complicated to work out how many of these parameters there are
                            // (blame George??)
                            if(iparam==param_tel_dx)N=psr[p].nTelDX;
                            if(iparam==param_tel_dy)N=psr[p].nTelDY;
                            if(iparam==param_tel_dz)N=psr[p].nTelDZ;
                            if(psr->param[param].val[0] == -1)N=0;
                            else if (psr->param[param].val[0] < 2 )N=N;
                            else N=N-1;
                            for (unsigned i = 0; i < N; ++i){
                                OUT.paramDerivs[OUT.nparams]     =t2FitFunc_telPos;
                                OUT.updateFunctions[OUT.nparams] =t2UpdateFunc_telPos;
                                ++OUT.nparams;
                            }
                            break;

                    case param_ifunc:
                    case param_clk_offs:
                    case param_quad_ifunc_p:
                    case param_quad_ifunc_c:
                            // ifunc-alikes
                            if(iparam==param_ifunc)N=psr->ifuncN;
                            if(iparam==param_clock)N=psr->clkOffsN;
                            if(iparam==param_quad_ifunc_p)N=psr->quad_ifuncN_p;
                            if(iparam==param_quad_ifunc_c)N=psr->quad_ifuncN_c;
                            for (unsigned i = 0; i < psr->N; ++i){
                                OUT.paramDerivs[OUT.nparams]     =t2FitFunc_ifunc;
                                OUT.updateFunctions[OUT.nparams] =t2UpdateFunc_ifunc;
                                OUT.paramCounters[OUT.nparams]=i;
                                ++OUT.nparams;
                            }
                            break;

                    case param_gwsingle:
                            for (unsigned i = 0; i < 4; ++i){
                                OUT.paramDerivs[OUT.nparams]     =t2FitFunc_stdGravWav;
                                OUT.updateFunctions[OUT.nparams] =t2UpdateFunc_stdGravWav;
                                OUT.paramCounters[OUT.nparams]=i;
                                ++OUT.nparams;
                            }
                            break;
                    case param_gwb_amp:
                    case param_gwm_amp:
                            break;

                    default:
                            logerr("ERROR: No methods for fitting parameter %d",i);
                            break;
                }
            }
        }
    }
}

