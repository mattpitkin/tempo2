#include "t2fit_dmmodel.h"
#include "t2fit_ifunc.h"
#include "TKlog.h"
#include <cmath>
#include <assert.h>



double t2FitFunc_dmmodelDM(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    const double dmf = 1.0/(DM_CONST*pow(psr->obsn[ipos].freqSSB/1.0e6,2));
    return dmf * ifunc(psr[ipsr].dmoffsDM_mjd,static_cast<double>(psr[ipsr].obsn[ipos].sat),psr[ipsr].dmoffsDMnum,k);
}
void t2UpdateFunc_dmmodelDM(pulsar *psr, int ipsr ,param_label label,int k, double val, double err){
    psr[ipsr].dmoffsDM[k] += val;
    psr[ipsr].dmoffsDM_error[k] = err;
}
double t2FitFunc_dmmodelCM(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    int idx=k-psr[ipsr].dmoffsDMnum;
    return ifunc(psr[ipsr].dmoffsCM_mjd,static_cast<double>(psr[ipsr].obsn[ipos].sat),psr[ipsr].dmoffsCMnum,idx);
}
void t2UpdateFunc_dmmodelCM(pulsar *psr, int ipsr ,param_label label,int k, double val, double err){
    int idx=k-psr[ipsr].dmoffsDMnum;
    psr[ipsr].dmoffsCM[idx] += val;
    psr[ipsr].dmoffsCM_error[idx] = err;
}
