#include <cmath>
#include <assert.h>
#include <string.h>
#include <cstdio>
#include "t2fit_nestlike.h"
#include "shapelet.h"
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

double t2FitFunc_nestlike_red_dm(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
    double freq = ((double)(k+1.0))/(maxtspan);
    double kappa = DM_CONST*1e-12;
    double ret=0;
    switch (label){
        case param_red_dm_cos:
	    ret = cos(2.0*M_PI*freq*x)/(kappa*(pow((double)psr[ipsr].obsn[ipos].freqSSB,2)));
            break;
        case param_red_dm_sin:
            ret = sin(2.0*M_PI*freq*x)/(kappa*(pow((double)psr[ipsr].obsn[ipos].freqSSB,2)));
            break;
        default:
            assert(0);
            break;
    }
    return ret;
}

double t2FitFunc_nestlike_red_chrom(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
    double freq = ((double)(k+1.0))/(maxtspan);

    double ret=0;
   
    double index=psr[ipsr].TNChromIdx;
    
    //double kappa = DM_CONST*powf(1e-6,index);

    //fprintf(stderr, "%.8le\n", psr[ipsr].obsn[ipos].freqSSB);

    double prefac=1.;

    switch (label){
        case param_red_chrom_cos:
	    ret = prefac*cos(2.0*M_PI*freq*x)/((pow((double)psr[ipsr].obsn[ipos].freqSSB/1.4e9,index)));
            break;
        case param_red_chrom_sin:
            ret = prefac*sin(2.0*M_PI*freq*x)/((pow((double)psr[ipsr].obsn[ipos].freqSSB/1.4e9,index)));
            break;
        default:
            assert(0);
            break;
    }
    return ret;
}





void t2UpdateFunc_nestlike_red(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    if (k==0 && label==param_red_sin){
        double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
        double RedFLow = pow(10., psr[ipsr].TNRedFLow);
        double freq = RedFLow*((double)(k+1.0))/(maxtspan);
        //printf("WAVE_FREQ %.18lg\n",freq*365.25);
        //printf("Temponest equivilent fitwaves parameters\n");
        //printf("========================================\n");
        //printf("\n");
        //printf("WAVE_OM %.18lg\n",freq*2.0*M_PI);
        //printf("WAVEEPOCH %.18lg\n",(double)psr[ipsr].param[param_pepoch].val[0]);
        if (writeResiduals&4){
            FILE *fout;
            fout = fopen("tnred.meta","w");
            if (!fout){
                printf("Unable to open file tnred.meta for writing\n");
            }
            fprintf(fout,"RED_OMEGA %lg\n",freq*2.0*M_PI);
            fprintf(fout,"RED_EPOCH %lg\n",(double)psr[ipsr].param[param_pepoch].val[0]);
            fclose(fout);
        }
    }
    /*
    if (label==param_red_sin){
        printf("WAVE_SIN%d %.18lg\n",k+1,-val);
    }
    if (label==param_red_cos){
        printf("WAVE_COS%d %.18lg\n",k+1,-val);
    }
    printf("\n");
    */

    //    logmsg("%d %s %d %lg %lg",ipsr,label_str[label],k,val,err);
    for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
        double x = (double)(psr[ipsr].obsn[iobs].bbat - psr[ipsr].param[param_pepoch].val[0]);
        double y = t2FitFunc_nestlike_red(psr,ipsr,x,iobs,label,k);
        psr[ipsr].obsn[iobs].TNRedSignal  += y *val;
        psr[ipsr].obsn[iobs].TNRedErr     += pow(y*err,2);
        //fprintf(stderr, "are we here????\n");

    }
}


void t2UpdateFunc_nestlike_red_dm(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {

    if (k==0 && label==param_red_dm_sin){
        if (writeResiduals&4){
            double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
            double freq = ((double)(k+1.0))/(maxtspan);
            FILE *fout;
            fout = fopen("tnreddm.meta","w");
            if (!fout){
                printf("Unable to open file tnred.meta for writing\n");
            }
            fprintf(fout,"REDDM_OMEGA %lg\n",freq*2.0*M_PI);
            fprintf(fout,"REDDM_EPOCH %lg\n",(double)psr[ipsr].param[param_pepoch].val[0]);
            fprintf(fout,"REDDM_DMEPOCH %lg\n",(double)psr[ipsr].param[param_dmepoch].val[0]);
            fprintf(fout,"REDDM_DM %lg\n",(double)psr[ipsr].param[param_dm].prefit[0]);
            fprintf(fout,"REDDM_DM1 %lg\n",(double)psr[ipsr].param[param_dm].prefit[1]);
            fprintf(fout,"REDDM_DM2 %lg\n",(double)psr[ipsr].param[param_dm].prefit[2]);
            fclose(fout);
        }
    }
    logdbg("%d %s %d %lg %lg",ipsr,label_str[label],k,val,err);
    for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
        double x = (double)(psr[ipsr].obsn[iobs].bbat - psr[ipsr].param[param_pepoch].val[0]);
        double y = t2FitFunc_nestlike_red_dm(psr,ipsr,x,iobs,label,k);
        psr[ipsr].obsn[iobs].TNDMSignal  += y *val;
        psr[ipsr].obsn[iobs].TNDMErr     += pow(y*err,2);

    }
}



void t2UpdateFunc_nestlike_red_chrom(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {

    if (k==0 && label==param_red_chrom_sin){
        if (writeResiduals&4){
            double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
            double freq = ((double)(k+1.0))/(maxtspan);
            FILE *fout; 
            fout = fopen("tnreddm.meta","w");
            if (!fout){
                printf("Unable to open file tnred.meta for writing\n");
            }
            fprintf(fout,"REDDM_OMEGA %lg\n",freq*2.0*M_PI);
            fprintf(fout,"REDDM_EPOCH %lg\n",(double)psr[ipsr].param[param_pepoch].val[0]);
            fprintf(fout,"REDDM_DMEPOCH %lg\n",(double)psr[ipsr].param[param_dmepoch].val[0]);
            fprintf(fout,"REDDM_DM %lg\n",(double)psr[ipsr].param[param_dm].prefit[0]);
            fprintf(fout,"REDDM_DM1 %lg\n",(double)psr[ipsr].param[param_dm].prefit[1]);
            fprintf(fout,"REDDM_DM2 %lg\n",(double)psr[ipsr].param[param_dm].prefit[2]);
            fclose(fout);
        }
    }
    logmsg("%d %s %d %lg %lg",ipsr,label_str[label],k,val,err);
    for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
        double x = (double)(psr[ipsr].obsn[iobs].bbat - psr[ipsr].param[param_pepoch].val[0]);
        double y = t2FitFunc_nestlike_red_chrom(psr,ipsr,x,iobs,label,k);
        psr[ipsr].obsn[iobs].TNChromSignal  += y *val;
        psr[ipsr].obsn[iobs].TNChromErr     += pow(y*err,2);

    }
}




double t2FitFunc_nestlike_jitter(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    const double dt=1./SECDAY;
    const double epoch = psr[ipsr].obsn[k].bat;
    const double batmin = epoch-dt;
    const double batmax = epoch+dt;
    const double bat = psr[ipsr].obsn[ipos].bat;
    if (bat > batmin && bat < batmax) return 1;
    else return 0;


}

void t2UpdateFunc_nestlike_jitter(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    logmsg("%d %s %.2f %lg %lg",ipsr,label_str[label],(double)psr[ipsr].obsn[k].bat,val,err);
}





double t2FitFunc_nestlike_band(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    //int totalBandNoiseCoeff=0; // unused

    int iband = 0;
    int ichan = k;
    while (ichan >= psr[ipsr].TNBandNoiseC[iband]){
        ichan -= psr[ipsr].TNBandNoiseC[iband];
        ++iband;
    }

    // we are in channel ichan, band iband.


    double ret = 0;

    double BandLF = psr[ipsr].TNBandNoiseLF[iband];
    double BandHF = psr[ipsr].TNBandNoiseHF[iband];
    //double BandAmp=pow(10.0, psr[ipsr].TNBandNoiseAmp[iband]);
    //double BandSpec=psr[ipsr].TNBandNoiseGam[iband];


    double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
    double freq = ((double)(ichan+1.0))/(maxtspan);

    // check that we are indeed in a band
    if (psr[ipsr].obsn[ipos].freq > BandLF && psr[ipsr].obsn[ipos].freq < BandHF) {
        switch (label){
            case param_band_red_cos:
                ret = cos(2.0*M_PI*freq*x);
                break;
            case param_band_red_sin:
                ret = sin(2.0*M_PI*freq*x);
                break;
            default:
                assert(0);
                break;
        }
    }
    return ret;
}


double t2FitFunc_nestlike_group(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){

    //int totalGroupCoeff=0;

    int igroup = 0;
    int ichan = k;
    while (ichan >= psr[ipsr].TNGroupNoiseC[igroup]){
        ichan -= psr[ipsr].TNGroupNoiseC[igroup];
        ++igroup;
    }

    // we are in channel ichan, group igroup.


    double ret=0.0;

    //double GroupAmp=pow(10.0, psr[ipsr].TNGroupNoiseAmp[igroup]);
    //double GroupSpec=psr[ipsr].TNGroupNoiseGam[igroup];
    //int GroupC=psr[ipsr].TNGroupNoiseC[igroup];

    double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
    double freq = ((double)(ichan+1.0))/(maxtspan);

    bool ingrp=false;

    for (int j=0; j < psr[ipsr].obsn[ipos].nFlags; j++) {
        //Check Group Noise Flag
        if (strcmp(
                    psr[ipsr].obsn[ipos].flagID[j],
                    psr[ipsr].TNGroupNoiseFlagID[igroup]
                  )==0 ) {
            if (strcmp(
                        psr[ipsr].obsn[ipos].flagVal[j],
                        psr[ipsr].TNGroupNoiseFlagVal[igroup]
                      )==0 ) {
                ingrp=true;
                break;
            }
        }
    }

    if(ingrp){ // we are in this group
        switch (label){
            case param_group_red_cos:
                ret = cos(2.0*M_PI*freq*x);
                break;
            case param_group_red_sin:
                ret = sin(2.0*M_PI*freq*x);
                break;
            default:
                assert(0);
                break;
        }

    }
    return ret;

}



void t2UpdateFunc_nestlike_band(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    int iband = 0;
    int ichan = k;
    while (ichan >= psr[ipsr].TNBandNoiseC[iband]){
        ichan -= psr[ipsr].TNBandNoiseC[iband];
        ++iband;
    }
    
    // we are in channel ichan, band iband.
    logmsg("%d %s %d %d %d %lg %lg",ipsr,label_str[label],k,iband,ichan,val,err);
    for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
        double x = (double)(psr[ipsr].obsn[iobs].bbat - psr[ipsr].param[param_pepoch].val[0]);
        double y = t2FitFunc_nestlike_band(psr,ipsr,x,iobs,label,k);
        psr[ipsr].obsn[iobs].TNRedSignal  += y *val;
        psr[ipsr].obsn[iobs].TNRedErr     += pow(y*err,2);

    }
}

void t2UpdateFunc_nestlike_group(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    int igroup = 0;
    int ichan = k;
    while (ichan >= psr[ipsr].TNGroupNoiseC[igroup]){
        ichan -= psr[ipsr].TNGroupNoiseC[igroup];
        ++igroup;
    }

    // we are in channel ichan, group igroup.
    logmsg("%d %s %d %d %d %lg %lg",ipsr,label_str[label],k,igroup,ichan,val,err);
    for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
        double x = (double)(psr[ipsr].obsn[iobs].bbat - psr[ipsr].param[param_pepoch].val[0]);
        double y = t2FitFunc_nestlike_group(psr,ipsr,x,iobs,label,k);
        psr[ipsr].obsn[iobs].TNRedSignal  += y *val;
        psr[ipsr].obsn[iobs].TNRedErr     += pow(y*err,2);

    }
}


double t2FitFunc_nestlike_shape_red(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    int ishape = k/MAX_TNShapeCoef;
    int icoef = k%MAX_TNShapeCoef;
    double pos = psr[ipsr].TNShapeletEvPos[ishape];
    double width = psr[ipsr].TNShapeletEvWidth[ishape];
    double t = psr[ipsr].obsn[ipos].bat;
    int N=icoef+1;
    double shapecoef[N];
    for (int ic=0; ic < N ;++ic)shapecoef[ic]=0;
    shapecoef[icoef]=1.0;
    double shape = evaluateShapelet(N,pos,width,shapecoef,t);
    return shape;
}

double t2FitFunc_nestlike_shape_dm(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    int ishape = k/MAX_TNShapeCoef;
    int icoef = k%MAX_TNShapeCoef;
    double pos = psr[ipsr].TNShapeletEvPos[ishape];
    double width = psr[ipsr].TNShapeletEvWidth[ishape];
    double t = psr[ipsr].obsn[ipos].bat;
    int N=icoef+1;
    double shapecoef[N];
    for (int ic=0; ic < N ;++ic)shapecoef[ic]=0;
    shapecoef[icoef]=1.0;
    double shape = evaluateShapelet(N,pos,width,shapecoef,t);
    return shape *pow((double)psr[ipsr].obsn[ipos].freqSSB,-2)*1e12/DM_CONST;
}

void t2UpdateFunc_nestlike_shape(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    int ishape = k/MAX_TNShapeCoef;
    int icoef = k%MAX_TNShapeCoef;
    psr[ipsr].TNShapeletEvCoef[ishape][icoef] += val;
    psr[ipsr].TNShapeletEvCoefErr[ishape][icoef] = err;
}
