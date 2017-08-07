
#include <tempo2.h>
#include <string.h>
void updateEpoch_str(pulsar* psr, int p, const char* newEpoch)
{
    int i;
    longdouble nMJD,dt;
    longdouble earliest=-1;
    longdouble latest=-1;
    int okay=1;

    //=psr[p].obsn[0].sat,latest=psr[p].obsn[0].sat;
    for (i=0;i<psr[p].nobs;i++)
    {
        okay=1;
        if (psr[p].obsn[i].deleted==1) okay=0;
        if (psr[p].param[param_start].paramSet[0]==1 && psr[p].param[param_start].fitFlag[0]==1 &&
                (psr[p].param[param_start].val[0] > psr[p].obsn[i].sat))
            okay=0;
        if (psr[p].param[param_finish].paramSet[0]==1 && psr[p].param[param_finish].fitFlag[0]==1 &&
                psr[p].param[param_finish].val[0] < psr[p].obsn[i].sat)
            okay=0;

        if (okay==1)
        {
            if (earliest==-1)
            {
                earliest = psr[p].obsn[i].sat;
                latest   = psr[p].obsn[i].sat;
            }
            if (earliest > psr[p].obsn[i].sat) earliest = psr[p].obsn[i].sat;
            if (latest < psr[p].obsn[i].sat) latest = psr[p].obsn[i].sat;
        }
    }


    if (strcasecmp(newEpoch,"CENTRE")==0 ||
            strcasecmp(newEpoch,"CENTER")==0) /* Find centre of data */
        nMJD = (int)((earliest+latest)/2.0);
    else if (strcasecmp(newEpoch,"LEFT")==0)
        nMJD = (int)(earliest);
    else if (strcasecmp(newEpoch,"RIGHT")==0)
        nMJD = (int)(latest);
    else
        nMJD = parse_longdouble(newEpoch);

    updateEpoch(psr,p,nMJD);
}

void updateEpoch(pulsar* psr, int p, longdouble nMJD) {
    double dt = (nMJD - psr[p].param[param_pepoch].val[0])*86400.0;
    printf("dt = %g\n",(double)dt);
    psr[p].param[param_f].val[0] = psr[p].param[param_f].val[0]+
        psr[p].param[param_f].val[1]*dt + 0.5*psr[p].param[param_f].val[2]*dt*dt;
    psr[p].param[param_f].val[1] += psr[p].param[param_f].val[2]*dt;
    if (psr[p].param[param_f].paramSet[3]==1)
    {
        psr[p].param[param_f].val[0] += 1.0/6.0*psr[p].param[param_f].val[3]*dt*dt*dt;
        psr[p].param[param_f].val[1] += 0.5*psr[p].param[param_f].val[3]*dt*dt;
        psr[p].param[param_f].val[2] += psr[p].param[param_f].val[3]*dt;
    }

    psr[p].param[param_f].prefit[0] = psr[p].param[param_f].val[0];
    psr[p].param[param_f].prefit[1] = psr[p].param[param_f].val[1];
    psr[p].param[param_pepoch].val[0] = nMJD;
    psr[p].param[param_pepoch].prefit[0] = nMJD;

    /* Update position epoch */
    dt = (nMJD - psr[p].param[param_posepoch].val[0])/365.25;
    printf("pos dt = %g\n",(double)dt);
    if (psr[p].param[param_pmra].paramSet[0]==1)
    {
        char retstr[1000];
        printf("Updating RAJ\n");
        psr[p].param[param_raj].val[0] = psr[p].param[param_raj].val[0]+psr[p].param[param_pmra].val[0]
            /cos(psr[p].param[param_decj].val[0])/1000.0*(M_PI/180.0)/60.0/60.0*dt; 
        psr[p].param[param_raj].prefit[0] = psr[p].param[param_raj].val[0];
        /* Must obtain this in hms form */
        turn_hms(psr[p].param[param_raj].val[0]/(2.0*M_PI), retstr);
        strcpy(psr[p].rajStrPost,retstr);
        strcpy(psr[p].rajStrPre,retstr);
    }

    if (psr[p].param[param_pmdec].paramSet[0]==1)
    {
        char retstr[1000];
        psr[p].param[param_decj].val[0] = psr[p].param[param_decj].val[0]+psr[p].param[param_pmdec].val[0]/1000.0*(M_PI/180.0)/60.0/60.0*dt; 
        psr[p].param[param_decj].prefit[0] = psr[p].param[param_decj].val[0];
        turn_dms(psr[p].param[param_decj].val[0]/(2.0*M_PI), retstr);
        strcpy(psr[p].decjStrPost,retstr);
        strcpy(psr[p].decjStrPre,retstr);
    }

    psr[p].param[param_posepoch].val[0] = nMJD;
    psr[p].param[param_posepoch].prefit[0] = nMJD;

    /* Update dmepoch */
    dt = (nMJD - psr[p].param[param_dmepoch].val[0])/365.25;
    psr[p].param[param_dm].val[0] = psr[p].param[param_dm].val[0]+
        psr[p].param[param_dm].val[1]*dt + 0.5*psr[p].param[param_dm].val[2]*dt*dt;
    psr[p].param[param_dmepoch].val[0] = nMJD;
    psr[p].param[param_dmepoch].prefit[0] = nMJD;

    /* Should update wave epoch */
    if (psr[p].param[param_waveepoch].paramSet[0]==1)
    {
        printf("ERROR: Not updating the FITWAVES to new EPOCH\n");
    }
    if (psr[p].param[param_waveepoch_dm].paramSet[0]==1)
    {
        printf("ERROR: Not updating the FITWAVES(DM) to new EPOCH\n");
    }



    /* Update binary parameters if necessary */
    if (psr[p].param[param_pb].paramSet[0]==1 && 0 == 1)  /* Binary pulsar */
    {
        longdouble orbits,pb,tt0,pbdot,xpbdot,t0p,t0m=0.0;
        int        norbits;

        if (psr[p].param[param_pbdot].paramSet[0]==1) pbdot = psr[p].param[param_pbdot].val[0];
        else pbdot = 0.0;

        if (psr[p].param[param_xpbdot].paramSet[0] == 1) xpbdot = psr[p].param[param_xpbdot].val[0];
        else xpbdot = 0.0;

        pb = psr[p].param[param_pb].val[0]*SECDAY;

        if (psr[p].param[param_tasc].paramSet[0]==1)
            tt0 = (nMJD-psr[p].param[param_tasc].val[0])*SECDAY;
        else
            tt0 = (nMJD-psr[p].param[param_t0].val[0])*SECDAY;
        ld_printf("tt0 = %.14Lf %.14Lf %.14Lf %.14Lf\n",tt0,nMJD,psr[p].param[param_t0].val[0],
                nMJD-psr[p].param[param_t0].val[0]);

        orbits = tt0/pb - 0.5*(pbdot+xpbdot)*(tt0/pb)*(tt0/pb);
        norbits = (int)(orbits+0.5);
        ld_printf("Orbits = %.5Lf  (%d)\n",orbits,norbits);

        if (xpbdot > 0 || pbdot > 0)
        {
            /*		  t0p = (1.0/pb + sqrtl(1.0/(pb*pb)-2.0*(pbdot+xpbdot)*norbits/pb/pb))
                      /((pbdot+xpbdot)/pb/pb);
                      t0m = (1.0/pb - sqrtl(1.0/(pb*pb)-2.0*(pbdot+xpbdot)*norbits/pb/pb))/((pbdot+xpbdot)/pb/pb); */
            t0p = pb/(pbdot+xpbdot)*(1.0+sqrtl(1.0-2.0*(pbdot+xpbdot)*(longdouble)norbits));
            t0m = pb/(pbdot+xpbdot)*(1.0-sqrtl(1.0-2.0*(pbdot+xpbdot)*(longdouble)norbits));

            if (psr[p].param[param_tasc].paramSet[0]==1)
            {
                t0p = psr[p].param[param_tasc].val[0]+t0p/SECDAY;
                t0m = psr[p].param[param_tasc].val[0]+t0m/SECDAY;     

                if (fabs(t0p-psr[p].param[param_tasc].val[0]) > fabs(t0m-psr[p].param[param_tasc].val[0]))
                    t0p = t0m;
            }
            else
            {
                t0p = psr[p].param[param_t0].val[0]+t0p/SECDAY;
                t0m = psr[p].param[param_t0].val[0]+t0m/SECDAY;     

                if (fabs(t0p-psr[p].param[param_t0].val[0]) > fabs(t0m-psr[p].param[param_t0].val[0]))
                    t0p = t0m;
            }
        }
        else
        {
            t0p = norbits*pb;
            if (psr[p].param[param_tasc].paramSet[0]==1)		  
                t0p = psr[p].param[param_tasc].val[0]+t0p/SECDAY;
            else
                t0p = psr[p].param[param_t0].val[0]+t0p/SECDAY;
        }

        if (psr[p].param[param_tasc].paramSet[0]==1)		  
        {
            psr[p].param[param_tasc].val[0] = t0p;
            psr[p].param[param_tasc].prefit[0] = psr[p].param[param_tasc].val[0];
        }
        else
        {
            psr[p].param[param_t0].val[0] = t0p;
            psr[p].param[param_t0].prefit[0] = psr[p].param[param_t0].val[0];
        }
        ld_printf("Result = %.14Lf %.14Lf\n",t0p,t0m);

        psr[p].param[param_pb].val[0] += psr[p].param[param_pbdot].val[0]*tt0/SECDAY;
        psr[p].param[param_pb].prefit[0] = psr[p].param[param_pb].val[0];

        psr[p].param[param_om].val[0] += psr[p].param[param_omdot].val[0]*tt0/SECDAY/365.25;
        psr[p].param[param_om].prefit[0] = psr[p].param[param_om].val[0];

        if (psr[p].param[param_a1dot].paramSet[0]==1)
        {
            psr[p].param[param_a1].val[0] += psr[p].param[param_a1dot].val[0]*1.0e-12*tt0;
            psr[p].param[param_a1].prefit[0] = psr[p].param[param_a1].val[0];
        }



    }
}
