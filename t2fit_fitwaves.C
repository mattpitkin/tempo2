#include <tempo2.h>
#include <math.h>
#include <assert.h>

double t2FitFunc_fitwaves(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int kk){

    assert(label==param_wave_dm || label==param_wave_om);
    double Xoff = 0; // this Xoff is because x is referenced to pepoch not waveepoch!
    double om   = psr[ipsr].param[label].val[0];
    double afactor=1.0;
    param_label w_epoch_par = param_waveepoch;
    int waveScale=psr[ipsr].waveScale;

    if (label==param_wave_dm){
        w_epoch_par = param_waveepoch_dm;
    }
    if (label==param_wave_dm || waveScale==1 || (waveScale==2 && kk < psr[ipsr].nWhite*2)){
        double freq = psr[ipsr].obsn[ipos].freqSSB/1.0e6;
        afactor=1.0/(DM_CONST*freq*freq);
    }
    int k = kk % (psr[ipsr].nWhite*2);

    if (psr[ipsr].param[w_epoch_par].paramSet[0]){
        Xoff = psr[ipsr].param[w_epoch_par].val[0] - psr[ipsr].param[param_pepoch].val[0];
    }

    if (k%2==0) return afactor * cos(om*(floor(k/2.0)+1)*(x-Xoff));
    else        return afactor * sin(om*(floor(k/2.0)+1)*(x-Xoff));

}

void t2UpdateFunc_fitwaves(pulsar *psr, int ipsr ,param_label label,int k, double val, double error){
    int waveScale=psr[ipsr].waveScale;
    if (label==param_wave_dm || waveScale==1 || (waveScale==2 && k < psr[ipsr].nWhite*2)){
        // it's a DM wave
        if(k%2==0){
            logmsg("dm.COS(%d) = %.3e", k/2, val);
            psr[ipsr].wave_cos_dm[k/2]  -= val;
            psr[ipsr].wave_cos_dm_err[k/2] = error;
        } else {
            logmsg("dm.SIN(%d) = %.3e", k/2, val);
            psr[ipsr].wave_sine_dm[k/2] -= val;
            psr[ipsr].wave_sine_dm_err[k/2] = error;
        }
    } else {
        k=k%(psr[ipsr].nWhite*2);
        // white noise
        if(k%2==0){
            logmsg("COS(%d) = %.3e", k/2, val);
            psr[ipsr].wave_cos[k/2]  -= val;
            psr[ipsr].wave_cos_err[k/2] = error;
        } else {
            logmsg("SIN(%d) = %.3e", k/2, val);
            psr[ipsr].wave_sine[k/2] -= val;
            psr[ipsr].wave_sine_err[k/2] = error;
        }
    }
}

