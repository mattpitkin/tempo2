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



void t2UpdateFunc_nestlike_red(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    logmsg("%d %s %d %lg %lg",ipsr,label_str[label],k,val,err);
}
