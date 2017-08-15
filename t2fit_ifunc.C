#include <cmath>
#include <assert.h>
#include "ifunc.h"
#include "t2fit_ifunc.h"



double t2FitFunc_sifunc(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    double ret=0;
    switch(label){
        case param_ifunc:
            ret=sinfunc(psr[ipsr].ifuncT,static_cast<double>(psr[ipsr].obsn[ipos].bat),k);
            break;
        case param_quad_ifunc_p:
            ret = psr[ipsr].quad_ifunc_geom_p *
                sinfunc(psr[ipsr].quad_ifuncT_p,static_cast<double>(psr[ipsr].obsn[ipos].bat),k);
            break;
        case param_quad_ifunc_c:
            ret = psr[ipsr].quad_ifunc_geom_c*
                sinfunc(psr[ipsr].quad_ifuncT_c,static_cast<double>(psr[ipsr].obsn[ipos].bat),k);
            break;
        default:
            assert(0);
            return 0;
    }
    return ret;
}

double t2FitFunc_ifunc(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    double ret=0;
    switch(label){
        case param_ifunc:
            ret= ifunc(psr[ipsr].ifuncT,static_cast<double>(psr[ipsr].obsn[ipos].sat),psr[ipsr].ifuncN,k);
            break;
        case param_quad_ifunc_p:
            ret = psr[ipsr].quad_ifunc_geom_p *
                ifunc(psr[ipsr].quad_ifuncT_p,static_cast<double>(psr[ipsr].obsn[ipos].sat),psr[ipsr].quad_ifuncN_p,k);
            break;
        case param_quad_ifunc_c:
            ret = psr[ipsr].quad_ifunc_geom_c *
                ifunc(psr[ipsr].quad_ifuncT_c,static_cast<double>(psr[ipsr].obsn[ipos].sat),psr[ipsr].quad_ifuncN_c,k);
            break;
        default:
            assert(0);
            break;
    }
    return ret;
}
void t2UpdateFunc_ifunc(pulsar *psr, int ipsr ,param_label label,int k, double val, double err){

    switch(label){
        case param_ifunc:
            psr[ipsr].ifuncV[k] -= val;
            psr[ipsr].ifuncE[k]  = err;
            break;
        case param_quad_ifunc_p:
            psr[ipsr].quad_ifuncV_p[k] -= val;
            psr[ipsr].quad_ifuncE_p[k]  = err;
            break;
        case param_quad_ifunc_c:
            psr[ipsr].quad_ifuncV_c[k] -= val;
            psr[ipsr].quad_ifuncE_c[k]  = err;
            break;
        default:
            assert(0);
    }

}
