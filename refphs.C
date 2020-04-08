#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "tempo2.h"

void refphs_init(pulsar* psr, int npsr){
    for (int p =0; p < npsr; ++p){
        if(psr[p].refphs==REFPHS_TZR){
            logdbg("Creating extra observation to store TZR parameters for reference phase");
            // make a fake observation for TZR parameters and add it to the end of the array
            if (psr[p].nobs==MAX_OBSN){
                logerr("Error, need at least 1 spare observation to reference to REFPHS TZR");
                exit(1);
            }
            psr[p].nobs++;
            observation* obs = &(psr[p].obsn[psr[p].nobs-1]); // convinience pointer
            memset(obs,0,sizeof(observation));
            obs->sat = psr[p].param[param_tzrmjd].val[0];
            obs->freq = psr[p].param[param_tzrfrq].val[0];
            strcpy(obs->telID,psr[p].tzrsite);
            if (obs->telID[0]=='@'            /* Have barycentric arrival time */
                    || strcasecmp(obs->telID,"bat")==0) 
            {
                obs->clockCorr=0;  /* therefore don't do clock corrections */
                obs->delayCorr=0;
            }
            else if (strcmp(obs->telID,"STL")==0)
            {
                obs->clockCorr=0;  /* don't do clock corrections */
                obs->delayCorr=1;
            }
            else if (strcmp(obs->telID,"STL_FBAT")==0)
            {
                obs->clockCorr=0;  /* don't do clock corrections */
                obs->delayCorr=1;
            }
            else
            {
                obs->clockCorr=1;
                obs->delayCorr=1;
            }
        }
    }

}

void refphs_clean(pulsar* psr, int npsr){
    for (int p =0; p < npsr; ++p){
        logdbg("Deleting extra observation to store TZR parameters for reference phase");
        if(psr->refphs==REFPHS_TZR){
            // remove our TZR observation and store it in tzrobs for later use.
            memcpy(&(psr[p].tzrobs),&(psr[p].obsn[psr[p].nobs-1]),sizeof(observation));
            psr[p].nobs--;
        }
    }

}
