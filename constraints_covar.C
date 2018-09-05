#include <cmath>
#include "constraints_covar.h"
#include <assert.h>
#include <string.h>


double constraints_covar_ifunc(pulsar *psr,int ipsr, int iconstraint,int iparam,int constraintk,int k,void* constraintSpecial){

    if((label)iparam == param_ifunc){
        int idx=psr[ipsr].ifuncN;
        if (idx >= 0){
            double* Binv_row = (double*)constraintSpecial;
//            logmsg("ipsr=%d iconstraint=%d iparam=%d constraintk=%d k=%d elem=%lg",ipsr,iconstraint,iparam,constraintk,idx, Binv_row[idx]);
            return Binv_row[idx];
        }
    }
    return 0;
}



double constraints_covar_dmmodel_cm(pulsar *psr,int ipsr, int iconstraint,int iparam,int constraintk,int k,void* constraintSpecial){

    if((label)iparam == param_dmmodel){
        int idx=k-psr[ipsr].dmoffsDMnum;
        if (idx >= 0){
            double* Binv_row = (double*)constraintSpecial;
//            logmsg("ipsr=%d iconstraint=%d iparam=%d constraintk=%d k=%d elem=%lg",ipsr,iconstraint,iparam,constraintk,idx, Binv_row[idx]);
            return Binv_row[idx];
        }
    }
    return 0;

}


double constraints_covar_dmmodel_dm(pulsar *psr,int ipsr, int iconstraint,int iparam,int constraintk,int k,void* constraintSpecial){

    if((label)iparam == param_dmmodel){
        int idx=k;
        if (idx < psr[ipsr].dmoffsDMnum){
            double* Binv_row = (double*)constraintSpecial;
//            logmsg("ipsr=%d iconstraint=%d iparam=%d constraintk=%d k=%d elem=%lg",ipsr,iconstraint,iparam,constraintk,idx, Binv_row[idx]);
            return Binv_row[idx];
        }
    }
    return 0;

}
