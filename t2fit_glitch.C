#include <tempo2.h>
#include <math.h>
#include <assert.h>

double t2FitFunc_stdGlitch(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){

    if (label==param_glph) {
        if (psr[ipsr].obsn[ipos].bbat >= psr[ipsr].param[param_glep].val[k])
            return 1.0/psr[ipsr].param[param_f].val[0];
        else
            return  0.0;
    }
    else if (label==param_glf0d)
    {
        longdouble dt1,expf,tp,tgl;

        tp = (psr[ipsr].obsn[ipos].bbat-psr[ipsr].param[param_pepoch].val[0])*86400.0;
        tgl = (psr[ipsr].param[param_glep].val[k] - psr[ipsr].param[param_pepoch].val[0])*86400.0;

        dt1 = tp-tgl;

        if (psr[ipsr].param[param_gltd].val[k]!=0.0)
            expf = exp(-dt1/86400.0/psr[ipsr].param[param_gltd].val[k]);
        else
            expf = 1.0;

        if (psr[ipsr].obsn[ipos].bbat >= psr[ipsr].param[param_glep].val[k])
        {
            return  psr[ipsr].param[param_gltd].val[k]*SECDAY*(1.0-expf)/psr[ipsr].param[param_f].val[0]; ///psr[ipsr].param[param_f].val[0];
            //	  printf("Glitch diff = %d %.10f %.10f %.10f\n",k+1,afunc,(double)tp,(double)tgl,(double)psr[ipsr].param[param_gltd].val[0]);

        }
        else
            return  0.0;
    }
    else if (label==param_gltd)
    {
        longdouble dt1,expf,tp,tgl;

        tp = (psr[ipsr].obsn[ipos].bbat-psr[ipsr].param[param_pepoch].val[0])*86400.0L;
        tgl = (psr[ipsr].param[param_glep].val[k] - psr[ipsr].param[param_pepoch].val[0])*86400.0L;

        dt1 = tp-tgl;

        if (psr[ipsr].param[param_gltd].val[k]!=0.0)
            expf = exp(-dt1/86400.0L/psr[ipsr].param[param_gltd].val[k]);
        else
            expf = 1.0;

        if (psr[ipsr].obsn[ipos].bbat >= psr[ipsr].param[param_glep].val[k])
            return  psr[ipsr].param[param_glf0d].val[k]*
                (1.0-(1.0+dt1/SECDAY/(psr[ipsr].param[param_gltd].val[k]))*expf)/psr[ipsr].param[param_f].val[0]*SECDAY;
        else
            return  0.0;
    }
    else if (label==param_glf0)
    {
        if (psr[ipsr].obsn[ipos].bbat >= psr[ipsr].param[param_glep].val[k]){
            return  (psr[ipsr].obsn[ipos].bbat-psr[ipsr].param[param_glep].val[k])*86400.0/psr[ipsr].param[param_f].val[0];
        }
        else
            return  0.0;
    }
    else if (label==param_glf1)
    {
        if (psr[ipsr].obsn[ipos].bbat >= psr[ipsr].param[param_glep].val[k])
            return  0.5*pow((psr[ipsr].obsn[ipos].bbat-psr[ipsr].param[param_glep].val[k])*86400.0,2)/psr[ipsr].param[param_f].val[0];
        else
            return  0.0;
    }
    else if (label==param_glf2)
    {
        if (psr[ipsr].obsn[ipos].bbat >= psr[ipsr].param[param_glep].val[k])
        {
            return  (double)(1.0L/6.0L*powl((psr[ipsr].obsn[ipos].bbat-psr[ipsr].param[param_glep].val[k])*86400.0L/1.0e9,3)/psr[ipsr].param[param_f].val[0]);
        }
        else
            return  0.0;
    }

    assert(false);
}

void t2UpdateFunc_stdGlitch(pulsar *psr, int ipsr ,param_label label,int k, double val, double error){
    double F=1;
    if(label==param_glf2)F=1e-27;
    psr[ipsr].param[label].val[k] -= val*F;
    psr[ipsr].param[label].err[k]  = error*F;
}
