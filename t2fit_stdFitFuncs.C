#include "t2fit_stdFitFuncs.h"
#include <enum_str.h>
#include <cmath>
#include <cstdlib>
#include <assert.h>
#include <cstring>
#include "ifunc.h"

/**
 *
 * The zero offset.
 *
 */
double t2FitFunc_zero(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    return 1;
}
void t2UpdateFunc_zero(pulsar *psr, int ipsr ,param_label label,int k, double val, double err){
    psr[ipsr].offset = val;
    psr[ipsr].offset_e = err;
}

/**
 * The pulse frequency, and derivatives
 */
double t2FitFunc_stdFreq(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    assert(label==param_f || label==param_brake);

    // this seems to be something to do with the breaking index fitting??
    long double arg3,arg4;

    long double f1,f0,bindex;
    long double t2;

    t2 = x*86400.L;

    arg3=t2*t2*t2;
    arg4=t2*t2*t2*t2;

    f0 = psr[ipsr].param[param_f].val[0];
    f1 = psr[ipsr].param[param_f].val[1];

    bindex = 0;//psr[ipsr].param[param_brake].val[0];

    if (label==param_brake) return f1*f1/f0*arg3/6.0L + (4*bindex-1)*f1*f1*f1/f0/f0*arg4/24.L;


    double v=0;

    switch(k){
        case 0:
            v = x*24.0L*3600.0L/f0;
            if (psr[ipsr].param[param_brake].paramSet[0] ==1)
            {
                v += (-bindex*f1*f1/f0/f0*arg3/6.L -2*bindex*(2*bindex-1)*f1*f1*f1/f0/f0/f0*arg4/24.) *84000.L/f0;
            }
            return v;
        case 1:
            v = 0.5L*x*x;
            if (psr[ipsr].param[param_brake].paramSet[0] ==1)
            {
                v += (2*bindex*f1/f0*arg3/6.L + 3*bindex*(2*bindex-1)*f1*f1/f0/f0*arg4/24.L)
                    *86400.L*86400.L/f0;
            }
            return v;
        case 2:
            return 1.0L/6.0L*x*x*x/1.0e9L;
        case 3:
            return 1.0L/24.0L*x/1.0e18L*x*x*x;
        case 4:
            return 1.0L/120.0L*x*x*x*x*x/1.0e18L;
        case 5:
            return 1.0L/720.0L*powl(x,6.0L)/1.0e18L;
        case 6:
            return 1.0L/5040.0L*powl(x,7.0L)/1.0e18L;
        case 7:
            return 1.0L/40320.0L*powl(x,8.0L)/1.0e18L;
        case 8:
            return 1.0L/362880.0L*powl(x,9.0L)/1.0e18L;
        case 9:
            return 1.0L/3628800.0L*powl(x,10.0L)/1.0e18L;
        case 10:
            return 1.0L/3628800.0L/11.0L*powl(x,11.0L)/1.0e23L;
        case 11:
            return 1.0L/3628800.0L/11.0L/12.0L*powl(x,12.0L)/1.0e23L;
        case 12:
            return 1.0L/3628800.0L/11.0L/12.0L/13.0L*powl(x,13.0L)/1.0e23L;
        case 13:
            return 1.0L/3628800.0L/11.0L/12.0L/13.0L/14.0L*powl(x,14.0L)/1.0e23L;
        default:
            // make an attempt at bigger powers... not sure it will work though.
            return 1.0L/3628800.0L/11.0L/12.0L/13.0L/14.0L*powl(x/(longdouble)k,(k+1))/1.0e23L;
    }

}
void t2UpdateFunc_stdFreq(pulsar *psr, int ipsr ,param_label label,int k, double val, double error){
    if (k==0)
    {
        psr[ipsr].param[param_f].val[k] *= (1.0-val/psr[ipsr].param[param_f].val[0]);
        psr[ipsr].param[param_f].err[k]  = error;
    }
    else
    {
        longdouble scale;
        scale=1.0L;
        if (k==2)      scale=1.0e9L;
        else if (k>2 && k<10)  scale=1.0e18L;
        else if (k>9) scale=1.0e23L;

        //logmsg("%lg",(double)(1.0/(pow(24.0*3600.0,k+1))/scale* psr[ipsr].param[param_f].val[0]));

        psr[ipsr].param[param_f].val[k] = psr[ipsr].param[param_f].val[k] -
            (psr[ipsr].param[param_f].val[0]*(val/pow(24.0*3600.0,k+1))/scale);
        psr[ipsr].param[param_f].err[k] = error/(pow(24.0*3600.0,k+1))/scale*
            psr[ipsr].param[param_f].val[0];
    }
}





/**
 * Binary models - need to select on psr[ipsr].binaryModel
 */
double t2FitFunc_binaryModels(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    double afunc;
    psr+=ipsr; // cheap way to get pointer to current pulsar.
    if (strcmp(psr->binaryModel,"BT")==0) /* Must be below other parameters */
        afunc = BTmodel(psr,0,ipos,label);
    else if (strcmp(psr->binaryModel,"BTJ")==0)
        afunc = BTJmodel(psr,0,ipos,label,k);
    else if (strcmp(psr->binaryModel,"BTX")==0)
        afunc = BTXmodel(psr,0,ipos,label,k);
    else if (strcmp(psr->binaryModel,"ELL1")==0)
        afunc = ELL1model(psr,0,ipos,label,k);
    else if (strcmp(psr->binaryModel,"DD")==0)
        afunc = DDmodel(psr,0,ipos,label);
    else if (strcmp(psr->binaryModel,"DDK")==0)
        afunc = DDKmodel(psr,0,ipos,label);
    else if (strcmp(psr->binaryModel,"DDS")==0)
        afunc = DDSmodel(psr,0,ipos,label);
    else if (strcmp(psr->binaryModel,"DDGR")==0)
        afunc = DDGRmodel(psr,0,ipos,label);
    else if (strcmp(psr->binaryModel,"MSS")==0)
        afunc = MSSmodel(psr,0,ipos,label);
    else if (strcmp(psr->binaryModel,"T2")==0)
        afunc = T2model(psr,0,ipos,label,k);
    else if (strcmp(psr->binaryModel,"T2-PTA")==0)
        afunc = T2_PTAmodel(psr,0,ipos,label,k);
    else if (strcmp(psr->binaryModel,"DDH")==0)
        afunc = DDHmodel(psr,0,ipos,label);
    else if (strcmp(psr->binaryModel,"ELL1H")==0)
        afunc = ELL1Hmodel(psr,0,ipos,label);
    else if (strcmp(psr->binaryModel,"ELL1k")==0)
        afunc = ELL1kmodel(psr,0,ipos,label);

    logdbg("%s(%d) %s %d %s(%d) %g",psr->name,ipsr,psr->binaryModel,ipos,label_str[label],k,afunc);
    return afunc;


}
void t2UpdateFunc_binaryModels(pulsar *psr, int ipsr ,param_label label,int k, double val, double error){
    if (strcmp(psr[ipsr].binaryModel,"BT")==0)
        updateBT(&psr[ipsr],val,error,label);
    else if (strcmp(psr[ipsr].binaryModel,"BTJ")==0)
        updateBTJ(&psr[ipsr],val,error,label,k);
    else if (strcmp(psr[ipsr].binaryModel,"BTX")==0)
        updateBTX(&psr[ipsr],val,error,label,k);
    else if (strcmp(psr[ipsr].binaryModel,"ELL1")==0)
        updateELL1(&psr[ipsr],val,error,label,k);
    else if (strcmp(psr[ipsr].binaryModel,"DD")==0)
        updateDD(&psr[ipsr],val,error,label);
    else if (strcmp(psr[ipsr].binaryModel,"DDK")==0)
        updateDDK(&psr[ipsr],val,error,label);
    else if (strcmp(psr[ipsr].binaryModel,"DDS")==0)
        updateDDS(&psr[ipsr],val,error,label);
    else if (strcmp(psr[ipsr].binaryModel,"DDGR")==0)
        updateDDGR(&psr[ipsr],val,error,label);
    else if (strcmp(psr[ipsr].binaryModel,"MSS")==0)
        updateMSS(&psr[ipsr],val,error,label);
    else if (strcmp(psr[ipsr].binaryModel,"T2")==0)
        updateT2(&psr[ipsr],val,error,label,k);
    else if (strcmp(psr[ipsr].binaryModel,"T2-PTA")==0)
        updateT2_PTA(&psr[ipsr],val,error,label,k);
    else if (strcmp(psr[ipsr].binaryModel,"DDH")==0)
        updateDDH(&psr[ipsr],val,error,label);
    else if (strcmp(psr[ipsr].binaryModel,"ELL1H")==0)
        updateELL1H(&psr[ipsr],val,error,label);
    else if (strcmp(psr[ipsr].binaryModel,"ELL1k")==0)
        updateELL1k(&psr[ipsr],val,error,label);
}


double t2FitFunc_planet(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    assert(label==param_dmassplanet);
    return dotproduct(psr->posPulsar,psr->obsn[ipos].planet_ssb[k]);
}

double t2FitFunc_stdDm(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    assert(label==param_dm);
    // freq=0 is infinite frequency, so no effect.
    if(psr[ipsr].obsn[ipos].freq==0) return 0;
    else if (k==0)
        return 1.0/(DM_CONST*powl(psr[ipsr].obsn[ipos].freqSSB/1.0e6,2));
    else
    {
        double yrs = (psr[ipsr].obsn[ipos].sat - psr[ipsr].param[param_dmepoch].val[0])/365.25;
        double series_factor=1.0;
        if (psr[ipsr].dm_series_type == series_taylor_pn) {
            for (int i =1; i <= k; ++i) series_factor *= i;
        }
        return 1.0/(DM_CONST*pow(psr[ipsr].obsn[ipos].freqSSB/1.0e6,2))*pow(yrs,k)/series_factor;
    }

}

double t2FitFunc_stdCm(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    assert(label==param_cm);
    double gam=psr[ipsr].TNChromIdx;
    //fprintf(stderr, "gam=%.3e\n",gam);
    // freq=0 is infinite frequency, so no effect.
    if(psr[ipsr].obsn[ipos].freq==0) return 0;
    else if (k==0)
    {
        // shouldn't fit for CM only its derivatives
        return 0;
        //return 1.0/(DM_CONST*powl(psr[ipsr].obsn[ipos].freqSSB/1.0e6,gam));
    }
    else
    {
        double yrs = (psr[ipsr].obsn[ipos].sat - psr[ipsr].param[param_dmepoch].val[0])/365.25;
        return 1.0/(DM_CONST*powl(psr[ipsr].obsn[ipos].freqSSB/1.0e6,gam))*pow(yrs,k);
    }

}


void t2UpdateFunc_simpleAdd(pulsar *psr, int ipsr ,param_label label,int k, double val, double error){
    psr[ipsr].param[label].val[k] += val;
    psr[ipsr].param[label].err[k]  = error;
}

void t2UpdateFunc_simpleMinus(pulsar *psr, int ipsr ,param_label label,int k, double val, double error){
    psr[ipsr].param[label].val[k] -= val;
    psr[ipsr].param[label].err[k]  = error;
}





double t2FitFunc_telPos_delta(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    double ret=0;

    int xyz;
    double* T;
    int N;
    switch(label) {
        case param_tel_dx:
            xyz=0;
            N=psr[ipsr].nTelDX;
            T=psr[ipsr].telDX_t;
            break;
        case param_tel_dy:
            xyz=1;
            N=psr[ipsr].nTelDY;
            T=psr[ipsr].telDY_t;
            break;
            break;
        case param_tel_dz:
            xyz=2;
            N=psr[ipsr].nTelDZ;
            T=psr[ipsr].telDZ_t;
            break;
            break;
        default:
            assert(0);
    }

    switch((int)(psr[ipsr].param[label].val[0])){
        case -1:
            if (strcmp(psr[ipsr].obsn[ipos].telID,"STL_FBAT")==0)
                ret= psr[ipsr].posPulsar[xyz]/SPEED_LIGHT*1000.0;
            break;
        case 2:
            if (T[k] <=  psr[ipsr].obsn[ipos].sat &&
                    T[k+1] > psr[ipsr].obsn[ipos].sat){
                ret= psr[ipsr].posPulsar[xyz];
            }
            break;
        default:
            ret = psr[ipsr].posPulsar[xyz] *
                ifunc(T,static_cast<double>(psr[ipsr].obsn[ipos].sat),N,k);
            break;
    }
    return ret;
}



void t2UpdateFunc_telPos_delta(pulsar *psr, int ipsr ,param_label label,int k, double val, double err){
    switch(label) {
        case param_tel_dx:
            psr[ipsr].telDX_v[k] -= val;
            psr[ipsr].telDX_e[k]  = err;
            break;
        case param_tel_dy:
            psr[ipsr].telDY_v[k] -= val;
            psr[ipsr].telDY_e[k]  = err;
            break;
        case param_tel_dz:
            psr[ipsr].telDZ_v[k] -= val;
            psr[ipsr].telDZ_e[k]  = err;
            break;
        default:
            assert(0);
    }
}


double t2FitFunc_telPos(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){

    double maclaurin_scale_factor=1.0;
    if(k==2) maclaurin_scale_factor=0.5;
    if(k==3) maclaurin_scale_factor=1.0/6.0;
    assert(k<=3);

    longdouble dt = (x + psr[ipsr].param[param_pepoch].val[0]) - psr[ipsr].param[param_telEpoch].val[0];
    longdouble arg=dt;
    for (int i=0; i< k;++i){
        arg*=dt;
    }

    arg*=maclaurin_scale_factor;

    if (strcmp(psr[ipsr].obsn[ipos].telID,"STL_FBAT")==0) {
        switch(label){

            case param_telx:
                return  psr[ipsr].posPulsar[0]*arg;
            case param_tely:
                return  psr[ipsr].posPulsar[1]*arg;
            case param_telz:
                return  psr[ipsr].posPulsar[2]*arg;
            default:
                assert(0);
        }

    } else {
        return 0;
    }
}

double t2FitFunc_jump(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    for (int l=0;l<psr[ipsr].obsn[ipos].obsNjump;l++){
        if (psr[ipsr].obsn[ipos].jump[l]==k) {
            if (psr[ipsr].jumpSAT[l]==0) {
                return -1.0;
            } else {
                //logmsg("SATJUMP  %s %d",label_str[label],k);
                return -1.0 * ((psr[ipsr].obsn[ipos].freq * 1e6) / psr[ipsr].obsn[ipos].freqSSB);
            }
        }
    }
    return 0;
}
void t2UpdateFunc_jump(pulsar *psr, int ipsr ,param_label label,int k, double val, double err){
    psr[ipsr].jumpVal[k] += val;
    psr[ipsr].jumpValErr[k] = err;
}



double t2FitFunc_fdjump(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    for (int l=0;l<psr[ipsr].obsn[ipos].obsNfdjump;l++){
        //fprintf(stderr, "%d %d\n", psr[ipsr].obsn[ipos].fdjump[l],k);
        if (psr[ipsr].obsn[ipos].fdjump[l]==k) {
            int idx=psr[ipsr].fdjumpIdx[k];    
            
            return pow(psr[ipsr].obsn[ipos].freqSSB/1e9,idx);
            //return pow(log(psr[ipsr].obsn[ipos].freqSSB/1e9),idx);
             }
    }
    return 0;
}


void t2UpdateFunc_fdjump(pulsar *psr, int ipsr ,param_label label,int k, double val, double err){
    psr[ipsr].fdjumpVal[k] += val;
    psr[ipsr].fdjumpValErr[k] = err;
}


double t2FitFunc_notImplemented(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    logerr("Parameter not implemented");
    exit(1);
}
void t2UpdateFunc_notImplemented(pulsar *psr, int ipsr ,param_label label,int k, double val, double err){
    logerr("Parameter not implemented");
    exit(1);
}
