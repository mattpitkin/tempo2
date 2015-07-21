#include <tempo2.h>
#include <math.h>
#include <string.h>
#include <assert.h>

double t2FitFunc_stdPosition(pulsar *psr, int ipsr ,double x ,int ipos ,param_label i,int k){

    double rce[3],re,deltae,alphae,psrra,psrdec,axy,s;
    int l;
    double conv_epoch = psr[ipsr].param[param_pepoch].val[0] - psr[ipsr].param[param_posepoch].val[0];

    /* What about observatory to Earth ???? */
    /* Calculated centre of Earth from SSB */
    for (l=0;l<3;l++) {
        rce[l] = psr[ipsr].obsn[ipos].earth_ssb[l];
    }
    re = sqrt(dotproduct(rce,rce));

    /* Calculate position of Earth w.r.t SSB */
    axy = rce[2]/AULTSC;
    s = axy/(re/AULTSC);  /* Why is this AULT and not AULTSC in TEMPO ??? */
    deltae = atan2(s,sqrt(1.0-s*s));
    alphae = atan2(rce[1],rce[0]);

    /* Calculate position of pulsar w.r.t SSB */
    /* IS THIS JUST THE RA AND DEC OF THE PULSAR? */
    psrra  = psr[ipsr].param[param_raj].val[0];
    psrdec = psr[ipsr].param[param_decj].val[0];
    if (i==param_raj)
    {
        return re*cos(deltae)*cos(psrdec)*sin(psrra - alphae);
    }
    else if (i==param_decj)
    {
        return re*(cos(deltae)*sin(psrdec)*cos(psrra - alphae) - sin(deltae)*cos(psrdec));
    }
    else if (i==param_pmra) {
        longdouble t0 = x + conv_epoch;
        return re*cos(deltae)*cos(psrdec)*sin(psrra - alphae) * t0;
    }
    else if (i==param_pmdec) { /* pmdec */
        longdouble t0 = x + conv_epoch;
        return re*(cos(deltae)*sin(psrdec)*cos(psrra - alphae) - sin(deltae)*cos(psrdec))*t0;
    }
    else if (i==param_px)
    {
        int l;
        double pxConv = 1.74532925199432958E-2/3600.0e3,rca[3],rr,rcos1;
        for (l=0;l<3;l++)
            rca[l] = psr[ipsr].obsn[ipos].earth_ssb[l] + psr[ipsr].obsn[ipos].observatory_earth[l];

        rr    = dotproduct(rca,rca);
        rcos1 = dotproduct(psr[ipsr].posPulsar,rca);
        longdouble afunc = 0.5*pxConv*(rr-rcos1*rcos1)/AULTSC;
        /* Now consider adding the effects of other distance determinations */
        if (psr[ipsr].param[i].nLinkFrom == 1 &&
                psr[ipsr].param[i].linkFrom[0] == param_dshk)
        {
            longdouble kpc2m = 3.08568025e19L;           /* 1 kpc in m        */
            longdouble mas_yr2rad_s = 1.536281850e-16L;  /* 1 mas/yr in rad/s */
            longdouble t0,afunc2;

            t0=x+conv_epoch;
            afunc2 = (t0*t0/2.0L/SPEED_LIGHT*
                    (psr[ipsr].param[param_pmra].val[0]*psr[ipsr].param[param_pmra].val[0]*
                     mas_yr2rad_s*mas_yr2rad_s+
                     psr[ipsr].param[param_pmdec].val[0]*psr[ipsr].param[param_pmdec].val[0]*
                     mas_yr2rad_s*mas_yr2rad_s)*kpc2m);	     
            afunc += (-afunc2/psr[ipsr].param[param_px].val[0]/psr[ipsr].param[param_px].val[0]);
        }
        return (double)afunc;
    }
    else if (i==param_pmrv)
    {
        int j;
        double delt,etat,dt_SSB,pmtrans_rcos2,rca[4];

        for (j=0;j<3;j++)
            rca[j] = psr[ipsr].obsn[ipos].earth_ssb[j] + psr[ipsr].obsn[ipos].observatory_earth[j];

        pmtrans_rcos2 = dotproduct(psr[ipsr].velPulsar,rca);      
        dt_SSB = 0.0; /* What should this be? */
        etat = 32.149618300000;/* NOT TRUE WHAT SHOULD THIS BE? (see dm_delays.c) */
        delt = (psr[ipsr].obsn[ipos].sat-psr[ipsr].param[param_posepoch].val[0] + (etat+dt_SSB)/SECDAY)/36525.0; 
        return pow(delt,2)*pmtrans_rcos2;
    }else if (i==param_dshk) {
        longdouble kpc2m = 3.08568025e19L;           /* 1 kpc in m        */
        longdouble mas_yr2rad_s = 1.536281850e-16L;  /* 1 mas/yr in rad/s */
        longdouble t0;
        t0 = ((x + psr[ipsr].param[param_pepoch].val[0])
                - psr[ipsr].param[param_posepoch].val[0])*SECDAY;
        return t0*t0/2.0L/SPEED_LIGHT*(psr[ipsr].param[param_pmra].val[0]*psr[ipsr].param[param_pmra].val[0]*
                mas_yr2rad_s*mas_yr2rad_s+
                psr[ipsr].param[param_pmdec].val[0]*psr[ipsr].param[param_pmdec].val[0]*
                mas_yr2rad_s*mas_yr2rad_s)*kpc2m;
    }

    assert(false);

}
void t2UpdateFunc_stdPosition(pulsar *psr, int ipsr ,param_label label,int k, double val, double err){
    char retstr[100];
    switch(label){
        case param_raj:
            {
                psr[ipsr].param[param_raj].val[k] += val;
                psr[ipsr].param[param_raj].err[k] = err;

                /* Must obtain this in hms form */
                turn_hms(psr[ipsr].param[param_raj].val[k]/(2.0*M_PI), retstr);
                strcpy(psr[ipsr].rajStrPost,retstr);
            }
            break;
        case param_decj:
            psr[ipsr].param[param_decj].val[k] += val;
            psr[ipsr].param[param_decj].err[k] = err;
            /* Must obtain this in dms form */
            turn_dms(psr[ipsr].param[param_decj].val[k]/(2.0*M_PI), retstr);
            strcpy(psr[ipsr].decjStrPost,retstr);
            break;
        case param_pmra:
            psr[ipsr].param[param_pmra].val[k] += val*180.0/M_PI*60.0*60.0*
                1000.0*SECDAY*365.25/24.0/3600.0*cos(psr[ipsr].param[param_decj].val[0]);
            psr[ipsr].param[param_pmra].err[k] = err*180.0/M_PI*60.0*60.0*
                1000.0*SECDAY*365.25/24.0/3600.0*cos(psr[ipsr].param[param_decj].val[0]);
            break;
        case param_pmdec:
            psr[ipsr].param[param_pmdec].val[k] += val*180.0/M_PI*60.0*60.0*1000.0*
                SECDAY*365.25/24.0/3600.0;
            psr[ipsr].param[param_pmdec].err[k] = err*180.0/M_PI*60.0*60.0*1000.0*
                SECDAY*365.25/24.0/3600.0;
            break;
        case param_pmrv:
            psr[ipsr].param[label].val[k] += 10.0*val*360.0*60.0*60.0/(2.0*M_PI);
            psr[ipsr].param[label].err[k]  = 10.0*err*360.0*60.0*60.0/(2.0*M_PI);
            break;
        default:
            // others are "simple" addition
            psr[ipsr].param[label].val[k] += val;
            psr[ipsr].param[label].err[k] = err;
            break;
    }
}
