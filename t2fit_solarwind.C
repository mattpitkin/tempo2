#include <t2fit_solarwind.h>
#include <assert.h>

#include "ifunc.h"

double t2FitFunc_ne_sw_sin(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    assert(label==param_ne_sw_sin);
    return -cos(2*M_PI*(psr[ipsr].obsn[ipos].sat - SOLAR_CYCLE_EPOCH)/SOLAR_CYCLE_PERIOD) * psr[ipsr].obsn[ipos].spherical_solar_wind;
}

double t2FitFunc_ne_sw(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    assert(label==param_ne_sw);
    return psr[ipsr].obsn[ipos].spherical_solar_wind;
}

double t2FitFunc_ne_sw_ifunc(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    assert(k < psr[ipsr].ne_sw_ifuncN);
    assert(label==param_ne_sw_ifunc);
    return psr[ipsr].obsn[ipos].spherical_solar_wind * ifunc(psr[ipsr].ne_sw_ifuncT,static_cast<double>(psr[ipsr].obsn[ipos].sat),psr[ipsr].ne_sw_ifuncN,k);
}


void t2UpdateFunc_ne_sw_ifunc(pulsar *psr, int ipsr ,param_label label,int k, double val, double error) {
    assert(k < psr[ipsr].ne_sw_ifuncN);
    psr[ipsr].ne_sw_ifuncV[k] += val;
    psr[ipsr].ne_sw_ifuncE[k] = error;
}

void t2UpdateFunc_ne_sw(pulsar *psr, int ipsr ,param_label label,int k, double val, double error) {
    assert(label==param_ne_sw);
    psr[ipsr].param[label].val[k] += val;
    psr[ipsr].param[label].err[k] = error;
    psr[ipsr].ne_sw = psr[ipsr].param[label].val[k];
}


double constraint_ne_sw_ifunc_function(pulsar *psr,int ipsr, int iconstraint,int iparam,int constraintk,int k,void* extra) {
    if(iparam == param_ne_sw_ifunc && k==constraintk){
        double err = *((double*)extra);
        return 1.0/err;
    } else {
        return 0;
    }
}

