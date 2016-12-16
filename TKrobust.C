#include <TKfit.h>
#include <TKrobust.h>
#include <TKlog.h>
#include <T2toolkit.h>
#include <cstring>

void TKrobust_reweight_huber(double* resid, double* data, double** M, double scale, int ndata, int nparams,
        double &sum_phi_sq, double &sum_dphi);

void TKrobust_reweight_bisquare(double* resid, double* data, double** M, double scale, int ndata, int nparams,
        double &sum_phi_sq, double &sum_dphi);

void TKrobust_reweight_welsch(double* resid, double* data, double** M, double scale, int ndata, int nparams,
        double &sum_phi_sq, double &sum_dphi);

/**
 * Perform an iterative robust least-squares by calling TKleastSquares
 *
 * Based on Wang, Keith, Stappers and Zheng (2016).
 */
double TKrobust(double* data, double* white_data,
        double** designMatrix, double** white_designMatrix,
        double** constraintsMatrix, int ndata,int nparams, int nconstraints, double tol, char rescale_errors,
        double* outP, double* e, double** Ocvm, char robust) {
    logmsg("Doing robust fit!");

    double **cvm_ls = malloc_blas(nparams,nparams);//store the covariance matrix of LS
    double **cvm = malloc_blas(nparams,nparams);//store the covariance matrix of LS

    int save_writeResiduals = writeResiduals; // we will disable parts of this later to prevent writing over the prefit residuals


    // First we need to solve for the scale factor s. This is not needed for least-squares since it is basically
    // scaling by the chi-squared. For robust fitting we cannot just scale in this way.
    //
    // We must iteratively fit and estimate s from equation 17 in the Wang paper.

    double** Mw = malloc_blas(ndata,nparams);
    double   modified_data[ndata];

    bool first =  true;

    double sum_phi_sq;
    double sum_dphi;
    double sigma=1;
    double sigma_old=0;

    double resid[ndata];
    double abs_resid[ndata];


    while(fabs(sigma_old-sigma) > 1e-5){

        memcpy(modified_data,white_data,sizeof(double)*ndata);
        memcpy(*Mw,*white_designMatrix,sizeof(double)*ndata*nparams);
        if (!first) {
            // first time around we just use normal least-squares
            TKrobust_reweight_huber(resid, modified_data, Mw, sigma, ndata, nparams, sum_phi_sq, sum_dphi);
        }

        TKrobustConstrainedLeastSquares(data, modified_data,
                designMatrix, Mw,
                constraintsMatrix,
                ndata, nparams, nconstraints,
                tol, rescale_errors, outP, e, cvm,0);
        if (first) {
            memcpy(*cvm_ls, *cvm, sizeof(double)*nparams*nparams); // save this for later.
            writeResiduals &= 0x4; // disable writing the prefit residuals for the rest of this routine.
            first=false;
        }


        for (int j=0; j<ndata; j++){
            double sum = 0.0;
            for (int k=0; k<nparams; k++){
                sum += outP[k]*white_designMatrix[j][k];
            }
            resid[j] = white_data[j]-sum;
            abs_resid[j] = fabs(resid[j]);
        }

        double median_abs_deviation = TKfindMedian_d(abs_resid, ndata);

        sigma_old = sigma;
        sigma = median_abs_deviation/0.6745;
        logdbg("Sigma = %lg",sigma);
    }


    // Once the scale has converged we can do the final fit and,
    // importantly, work out the parameter covariance matrix.

    memcpy(modified_data,white_data,sizeof(double)*ndata);
    memcpy(*Mw,*white_designMatrix,sizeof(double)*ndata*nparams);
    switch (robust) {
        case 'H':
            TKrobust_reweight_huber(resid, modified_data, Mw, sigma, ndata, nparams, sum_phi_sq, sum_dphi);
            break;
        case 'B':
            TKrobust_reweight_bisquare(resid, modified_data, Mw, sigma, ndata, nparams, sum_phi_sq, sum_dphi);
            break;
        case 'W':
            TKrobust_reweight_welsch(resid, modified_data, Mw, sigma, ndata, nparams, sum_phi_sq, sum_dphi);
            break;
        default:
            logerr("Robust type '%c' not supported. (fall back to huber!)",robust);
            TKrobust_reweight_huber(resid, modified_data, Mw, sigma, ndata, nparams, sum_phi_sq, sum_dphi);
            break;

    }
    TKrobustConstrainedLeastSquares(data, modified_data,
            designMatrix, Mw,
            constraintsMatrix,
            ndata, nparams, nconstraints,
            tol, rescale_errors, outP, e, Ocvm,0);

    double robust_chisq = sigma*sigma * ndata*ndata * sum_phi_sq / (sum_dphi*sum_dphi);

    double nfree = ndata + nconstraints - nparams;

    // Some strangeness here. I think that tempo2 will automatically assume you will correct Ocvm by sqrt(chisq/nfree)
    // so we don't correct it here.
    double ecorrect=1.0;
    if (rescale_errors)ecorrect = robust_chisq / nfree;
    for (int i=0; i< nparams; ++i){
        for (int j=0; j< nparams; ++j){
            Ocvm[i][j] = cvm_ls[i][j];
        }
        e[i] = sqrt(cvm[i][i]* ecorrect);
    }

    // Free memory
    free_blas(cvm_ls);
    free_blas(cvm);

    writeResiduals = save_writeResiduals;
    return robust_chisq;
}



void TKrobust_reweight_huber(double* resid, double* data, double** M, double scale, int ndata, int nparams,
        double &sum_phi_sq, double &sum_dphi){

    const double c0 = 1.345;

    double weight;
    sum_phi_sq = 0;
    sum_dphi   = 0;

    for (int idata=0; idata < ndata; ++idata){
        resid[idata] /= scale;
        if(fabs(resid[idata]) < c0){
            // phi = r
            weight = 1.0;
            sum_phi_sq += resid[idata]*resid[idata];
            sum_dphi   += 1.0;
        } else {
            // phi = c0
            weight = c0/fabs(resid[idata]);
            sum_phi_sq += c0*c0;
            sum_dphi   += 0.0;
        }
        data[idata] *= weight;
        for(int iparam=0; iparam < nparams; ++iparam){
            M[idata][iparam] *=weight;
        }


    }
}


void TKrobust_reweight_welsch(double* resid, double* data, double** M, double scale, int ndata, int nparams,
        double &sum_phi_sq, double &sum_dphi){

    const double c0 = 1.345;
    const double aWelsch = 2.11;

    double weight;
    sum_phi_sq = 0;
    sum_dphi   = 0;

    for (int idata=0; idata < ndata; ++idata){
        resid[idata] /= scale;
        // phi = a^2[1-exp(-(eta/a)^2)]
        weight = 2*aWelsch*exp(-pow(resid[idata],2));

        sum_phi_sq += pow(2*aWelsch*resid[idata]*exp(-pow(resid[idata]/aWelsch, 2)), 2);
        sum_dphi   += 2*aWelsch*(1-2*pow(resid[idata]/aWelsch, 2))*exp(-pow(resid[idata]/aWelsch, 2));

        data[idata] *= weight;
        for(int iparam=0; iparam < nparams; ++iparam){
            M[idata][iparam] *=weight;
        }

    }

}



void TKrobust_reweight_bisquare(double* resid, double* data, double** M, double scale, int ndata, int nparams,
        double &sum_phi_sq, double &sum_dphi){

    const double aBisquare = 4.685;

    double weight;
    sum_phi_sq = 0;
    sum_dphi   = 0;

    for (int idata=0; idata < ndata; ++idata){
        resid[idata] /= scale;
        if(fabs(resid[idata]) < aBisquare){
            // phi = bisquare
            weight = pow(1-pow(resid[idata]/aBisquare, 2),2);
            sum_phi_sq += pow(resid[idata]*pow(1-pow(resid[idata]/aBisquare, 2),2), 2);
            sum_dphi   += (1-5*pow(resid[idata]/aBisquare, 2))*pow(1-pow(resid[idata]/aBisquare, 2),2);
        } else {
            // phi = 0
            weight = 0.0;
            sum_phi_sq += 0.0;
            sum_dphi   += 0.0;
        }
        data[idata] *= weight;
        for(int iparam=0; iparam < nparams; ++iparam){
            M[idata][iparam] *=weight;
        }


    }
}
