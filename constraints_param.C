
#include "tempo2.h"
#include "enum_str.h"
#include "constraints_param.h"
/*
 * struct constraint_param_info {
    label param;
    int param_k;
    double val;
    double err;
} */

double constraint_param_function(pulsar *psr,int ipsr, int iconstraint,int iparam,int constraintk,int k,void* extra) {
    struct constraint_param_info *nfo = static_cast<struct constraint_param_info*>(extra);
    logdbg("__: %s %s %s %d %d",psr[ipsr].name,label_str[nfo->param],label_str[iparam],k,nfo->param_k);

    if(iparam == nfo->param && k==nfo->param_k){
        logdbg("!!: %s %s %d",psr[ipsr].name,label_str[iparam],k);
        return 1.0/nfo->err;
    } else {
        return 0;
    }
}

