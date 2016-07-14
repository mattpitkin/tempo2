#include <t2fit_dmother.h>
#include <assert.h>
/*
   case param_fddc:
   case param_fddi:
   case param_fd:
   case param_dm_sin1yr:
   case param_dm_cos1yr:
   case param_dmx:
   */

double t2FitFunc_dmx(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    assert(label==param_dmx);

    if ((psr[ipsr].obsn[ipos].sat > psr[ipsr].param[param_dmxr1].val[k])
            && (psr[ipsr].obsn[ipos].sat < psr[ipsr].param[param_dmxr2].val[k])) {

        return 1.0/(DM_CONST*powl(psr[ipsr].obsn[ipos].freqSSB/1.0e6,2));
    } else {
        return 0.0;
    }
}

double t2FitFunc_fddc(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    assert(label==param_fddc);
    return 1.0/(pow(psr[ipsr].obsn[ipos].freqSSB/1.0e6,psr[ipsr].param[param_fddi].val[0]));
}

double t2FitFunc_fd(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    assert(label==param_fd);
    return pow(log(psr[ipsr].obsn[ipos].freqSSB/1e9),k+1);
}

double t2FitFunc_dmsinusoids(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    if (label==param_dm_sin1yr){
        return 1.0/(DM_CONST*powl(psr[ipsr].obsn[ipos].freqSSB/1.0e6,2))*sin(2*M_PI/(365.25)*(psr[ipsr].obsn[ipos].sat - psr[ipsr].param[param_dmepoch].val[0]));
    } else if (label==param_dm_cos1yr) {
        return 1.0/(DM_CONST*powl(psr[ipsr].obsn[ipos].freqSSB/1.0e6,2))*cos(2*M_PI/(365.25)*(psr[ipsr].obsn[ipos].sat - psr[ipsr].param[param_dmepoch].val[0]));
    } else {
        logerr("Called dmsinusoids without dmsin1yr or dmcos1yr");
        return 0;
    }
}
