#include <t2fit_solarwind.h>
#include <assert.h>


double t2FitFunc_ne_sw_sin(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    assert(label==param_ne_sw_sin);
    return -cos(2*M_PI*(psr[ipsr].obsn[ipos].sat - SOLAR_CYCLE_EPOCH)/SOLAR_CYCLE_PERIOD) * psr[ipsr].obsn[ipos].spherical_solar_wind;
}

double t2FitFunc_ne_sw(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    assert(label==param_ne_sw);
    return psr[ipsr].obsn[ipos].spherical_solar_wind;
}


void t2UpdateFunc_ne_sw(pulsar *psr, int ipsr ,param_label label,int k, double val, double error) {
    assert(label==param_ne_sw);
    psr[ipsr].param[label].val[k] += val;
    psr[ipsr].param[label].err[k] = error;
    psr[ipsr].ne_sw = psr[ipsr].param[label].val[k];
}
