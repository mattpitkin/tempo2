#include <cmath>
#include <assert.h>
#include "t2fit_ifunc.h"


double ifunc(const double* mjd, const double t,const int N, const int k){
    double yoffs[MAX_IFUNC];
    logdbg("ifunc %lg => %lg : %lg %d %d",mjd[0],mjd[N-1],t,N,k);
    assert(k<N);

    for (int ioff =0;ioff<N;ioff++){
        if (ioff==k){
            yoffs[ioff]=1;
        } else {
            yoffs[ioff]=0;
        }
    }

    if (t < mjd[0]){
        // we are before the first jump
        // so our gradient is just the zeroth offset.
        return yoffs[0];
    } else if(t > mjd[N-1]){
        return yoffs[N-1];
    } else{
        // find the pair we are between...
        for (int ioff =0;ioff<N;ioff++){
            if(t >= mjd[ioff] && t < mjd[ioff+1]){
                double x1 = mjd[ioff];
                double x2 = mjd[ioff+1];
                double x = (t-x1)/(x2-x1);
                double y1=yoffs[ioff];
                double y2=yoffs[ioff+1];
                return (y2-y1)*x + y1;
            }
        }
    }
    logerr("Shouldn't get here");
    assert(0);
    return 0;
}

double sinfunc(const double *T, const double bat, const int k){
    double dt = bat - T[k];
    double tt = M_PI/(T[1] - T[0])*dt;
    return sin(tt)/tt;
}

double t2FitFunc_sifunc(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    switch(label){
        case param_ifunc:
            return sinfunc(psr[ipsr].ifuncT,static_cast<double>(psr[ipsr].obsn[ipos].bat),k);
        case param_clk_offs:
        case param_quad_ifunc_p:
        case param_quad_ifunc_c:
        default:
            assert(0);
            return 0;
    }
}

double t2FitFunc_ifunc(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    double ret=0;
    switch(label){
        case param_ifunc:
            ret = ifunc(psr[ipsr].ifuncT,static_cast<double>(psr[ipsr].obsn[ipos].sat),psr[ipsr].ifuncN,k);
            break;
        case param_clk_offs:
        case param_quad_ifunc_p:
        case param_quad_ifunc_c:
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
        case param_clk_offs:
        case param_quad_ifunc_p:
        case param_quad_ifunc_c:
        default:
            assert(0);
    }

}
