#include <cmath>
#include <assert.h>
#include "t2fit_nestlike.h"
#include "enum_str.h"


double t2FitFunc_nestlike_red(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
    double RedFLow = pow(10., psr[ipsr].TNRedFLow);
    double freq = RedFLow*((double)(k+1.0))/(maxtspan);
    double ret=0;
    switch (label){
        case param_red_cos:
            ret = cos(2.0*M_PI*freq*x);
            break;
        case param_red_sin:
            ret = sin(2.0*M_PI*freq*x);
            break;
        default:
            assert(0);
            break;
    }
    return ret;
}

double t2FitFunc_nestlike_red_dm(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
    double freq = ((double)(k+1.0))/(maxtspan);
    double kappa = DM_CONST*1e-12;
    double ret=0;
    switch (label){
        case param_red_dm_cos:
            ret = cos(2.0*M_PI*freq*x)/(kappa*(pow((double)psr[ipsr].obsn[ipos].freqSSB,2)));
            break;
        case param_red_dm_sin:
            ret = sin(2.0*M_PI*freq*x)/(kappa*(pow((double)psr[ipsr].obsn[ipos].freqSSB,2)));
            break;
        default:
            assert(0);
            break;
    }
    return ret;
}



void t2UpdateFunc_nestlike_red(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    logmsg("%d %s %d %lg %lg",ipsr,label_str[label],k,val,err);
    for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
        double x = (double)(psr[ipsr].obsn[iobs].bbat - psr[ipsr].param[param_pepoch].val[0]);
            double y = t2FitFunc_nestlike_red(psr,ipsr,x,iobs,label,k);
        psr[ipsr].obsn[iobs].TNRedSignal  += y *val;
        psr[ipsr].obsn[iobs].TNRedErr     += pow(y*err,2);

    }
}


void t2UpdateFunc_nestlike_red_dm(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    logmsg("%d %s %d %lg %lg",ipsr,label_str[label],k,val,err);
    for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
        double x = (double)(psr[ipsr].obsn[iobs].bbat - psr[ipsr].param[param_pepoch].val[0]);
            double y = t2FitFunc_nestlike_red_dm(psr,ipsr,x,iobs,label,k);
        psr[ipsr].obsn[iobs].TNDMSignal  += y *val;
        psr[ipsr].obsn[iobs].TNDMErr     += pow(y*err,2);

    }
}



double t2FitFunc_nestlike_jitter(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    const double dt=10.0/SECDAY;
    const double epoch = psr[ipsr].obsn[k].bat;
    const double batmin = epoch-dt;
    const double batmax = epoch+dt;
    const double bat = psr[ipsr].obsn[ipos].bat;
    if (bat > batmin && bat < batmax) return 1;
    else return 0;
    

}

void t2UpdateFunc_nestlike_jitter(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    logmsg("%d %s %.2f %lg %lg",ipsr,label_str[label],(double)psr[ipsr].obsn[k].bat,val,err);
}





double t2FitFunc_nestlike_band(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
	int totalBandNoiseCoeff=0;

    int iband = 0;
    int ichan = k;
    while (ichan >= psr[ipsr].TNBandNoiseC[iband]){
        ichan -= psr[ipsr].TNBandNoiseC[iband];
        ++iband;
    }

    // we are in channel ichan, band iband.


    double ret = 0;

    double BandLF = psr[ipsr].TNBandNoiseLF[iband];
    double BandHF = psr[ipsr].TNBandNoiseHF[iband];
    double BandAmp=pow(10.0, psr[ipsr].TNBandNoiseAmp[iband]);
    double BandSpec=psr[ipsr].TNBandNoiseGam[iband];


    double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
    double freq = ((double)(ichan+1.0))/(maxtspan);

    // check that we are indeed in a band
    if (psr[ipsr].obsn[ipos].freq > BandLF && psr[ipsr].obsn[ipos].freq < BandHF) {
        switch (label){
            case param_band_red_cos:
                ret = cos(2.0*M_PI*freq*x);
                break;
            case param_band_red_sin:
                ret = sin(2.0*M_PI*freq*x);
                break;
            default:
                assert(0);
                break;
        }
    }
    return ret;
}

void t2UpdateFunc_nestlike_band(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    int iband = 0;
    int ichan = k;
    while (ichan >= psr[ipsr].TNBandNoiseC[iband]){
        ichan -= psr[ipsr].TNBandNoiseC[iband];
        ++iband;
    }

    // we are in channel ichan, band iband.
    logmsg("%d %s %d %d %d %lg %lg",ipsr,label_str[label],k,iband,ichan,val,err);
    for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
        double x = (double)(psr[ipsr].obsn[iobs].bbat - psr[ipsr].param[param_pepoch].val[0]);
        double y = t2FitFunc_nestlike_band(psr,ipsr,x,iobs,label,k);
        psr[ipsr].obsn[iobs].TNRedSignal  += y *val;
        psr[ipsr].obsn[iobs].TNRedErr     += pow(y*err,2);

    }
}
