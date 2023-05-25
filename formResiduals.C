#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

/*
 *    This file is part of TEMPO2. 
 * 
 *    TEMPO2 is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    TEMPO2 is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
 */

/*
 *    If you use TEMPO2 then please acknowledge it by citing 
 *    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
 *    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
 *    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
 *    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
 *    timing model.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"
#include "GWsim.h"
#include "ifunc.h"
#include "shapelet.h"
#include <vector>
#include <algorithm>
/* Form the timing residuals from the timing model and the barycentric arrival times */
void residualTracking(pulsar *psr);




void averageDMResiduals(pulsar *psr, int npsr){

    int systemcount = 0;
    std::vector<std::string>systemnames;

    double mintime=100000.0;
    double maxtime=0.0;
    float timestep = psr[0].AverageEpochWidth;
    char* flagName = psr[0].AverageFlag;
    fprintf(stderr, "width %.3e flag %s\n", timestep, flagName);

    for(int o=0;o<psr[0].nobs;o++){
        int found=0;
        if((double)psr[0].obsn[o].bat > maxtime){maxtime=(double)psr[0].obsn[o].bat;}
        if((double)psr[0].obsn[o].bat < mintime){mintime=(double)psr[0].obsn[o].bat;}
        for (int f=0;f<psr[0].obsn[o].nFlags;f++){
            if(strcasecmp(psr[0].obsn[o].flagID[f],flagName)==0){

	      //fprintf(stderr, "found flag\n");                
                if(std::find(systemnames.begin(), systemnames.end(), psr[0].obsn[o].flagVal[f]) != systemnames.end()) {
                } else {

                    systemnames.push_back(psr[0].obsn[o].flagVal[f]);
                    systemcount++;
                }
                found=1;
            }
        }

        if(found==0){
            printf("Observation %i is missing the -f flag, please check before continuing\n",o);return;
        }
    }

    mintime=floor(mintime);
    maxtime=ceil(maxtime);


    int numbins=ceil((maxtime-mintime)/timestep);
    double **AverageFreq = new double*[systemcount];
    double **AverageRes = new double*[systemcount];
    double **AverageWeight = new double*[systemcount];
    double **AverageBat = new double*[systemcount];
    for(int i = 0; i < systemcount; i++){
        AverageRes[i] = new double[numbins];
        AverageWeight[i] = new double[numbins];
        AverageBat[i] = new double[numbins];
	AverageFreq[i] = new double[numbins];
        for(int j = 0; j < numbins; j++){
            AverageRes[i][j] = 0;
            AverageWeight[i][j] = 0;
            AverageBat[i][j] = 0;
	    AverageFreq[i][j]= 0;
        }
    }




    for(int o=0;o<psr[0].nobs;o++){

        for (int f=0;f<psr[0].obsn[o].nFlags;f++){
            if(strcasecmp(psr[0].obsn[o].flagID[f],flagName)==0){

                int flagindex=-1;
                for(int s = 0; s < static_cast<int>(systemnames.size()); s++){
                    if(strcasecmp(psr[0].obsn[o].flagVal[f],systemnames[s].c_str())==0){
                        flagindex = s;
                    }
                }

                double freq=psr[0].obsn[o].freq;

                double resid=(double)psr[0].obsn[o].residual;
                //double residDM=0;
                //double residTN=0;

                if (psr[0].TNsubtractRed ==1)
                {
                    resid  -= psr[0].obsn[o].TNRedSignal;

                }

                if (psr[0].TNsubtractDM ==1)
                {
                    resid -= psr[0].obsn[o].TNDMSignal;
                }
                if (psr[0].TNsubtractChrom ==1)
                {
                    resid -= psr[0].obsn[o].TNChromSignal;
                }



                int bin = floor(((double)psr[0].obsn[o].bat-mintime)/timestep);
                double adjustedErr = pow(psr[0].obsn[o].toaErr*pow(10.0, -6), 2);
                AverageWeight[flagindex][bin] += 1.0/adjustedErr;
                AverageRes[flagindex][bin] += (double)resid*powf(freq/1400.,2.)/adjustedErr;
                AverageBat[flagindex][bin] += (double)psr[0].obsn[o].bat/adjustedErr;
                AverageFreq[flagindex][bin] += freq/adjustedErr;
            }
        }

    }

    for(int o=0;o<psr[0].nobs;o++){


        //Is there an ECORR??
        double EcorrVal=0;
        for (int j=0;j<psr->obsn[o].nFlags;j++){
            for (int k=0;k<psr->nTNECORR;k++){
                if (strcmp(psr->obsn[o].flagID[j], psr->TNECORRFlagID[k])==0){
                    if (strcmp(psr->obsn[o].flagVal[j],psr->TNECORRFlagVal[k])==0){
                        EcorrVal=psr->TNECORRVal[k]*pow(10.0, -6);
                        //printf("Ecorr: %i %g \n", o, EcorrVal);
                    }
                }
            }
        }

        // SECORR
        double SEcorrVal=0;
        for (int j=0;j<psr->obsn[o].nFlags;j++){
            for (int k=0;k<psr->nTNSECORR;k++){
                if (strcmp(psr->obsn[o].flagID[j], psr->TNSECORRFlagID[k])==0){
                    if (strcmp(psr->obsn[o].flagVal[j],psr->TNSECORRFlagVal[k])==0){
		      SEcorrVal=psr->TNSECORRVal[k]*pow(10.0, -6)/sqrt(psr->obsn[o].tobs/3600.);
		      //printf("SEcorr: %i %g \n", o, SEcorrVal);
		      //	exit(0);
                    }
                }
            }
        }



        for (int f=0;f<psr[0].obsn[o].nFlags;f++){
            if(strcasecmp(psr[0].obsn[o].flagID[f],flagName)==0){
                int flagindex=-1;
                for(int s = 0; s < static_cast<int>(systemnames.size()); s++){
                    if(strcasecmp(psr[0].obsn[o].flagVal[f],systemnames[s].c_str())==0){
                        flagindex = s;
                    }
                }
                int bin = floor(((double)psr[0].obsn[o].bat-mintime)/timestep);
		double freq = AverageFreq[flagindex][bin]/AverageWeight[flagindex][bin];
		
		fprintf(stderr, "freq: %.3lf\n", freq);
                psr[0].obsn[o].averagedmbat = AverageBat[flagindex][bin]/AverageWeight[flagindex][bin];; 
                psr[0].obsn[o].averagedmres  = AverageRes[flagindex][bin]/AverageWeight[flagindex][bin];
                psr[0].obsn[o].averagedmerr = powf(freq/1400,2.)*sqrt(1.0/AverageWeight[flagindex][bin]  + pow(EcorrVal, 2) + pow(SEcorrVal,2));
            }
        }
    }

    for(int o=0;o<psr[0].nobs;o++){
        //		printf("Average: %g %g %g \n", psr[0].obsn[o].averagebat, psr[0].obsn[o].averageres, psr[0].obsn[o].averageerr);
    }



    for(int i = 0; i < systemcount; i++){
        delete[] AverageRes[i];
        delete[] AverageWeight[i];
        delete[] AverageBat[i];
    }
    delete[] AverageRes;
    delete[] AverageWeight;
    delete[] AverageBat;


}


void averageResiduals(pulsar *psr, int npsr){

    int systemcount = 0;
    std::vector<std::string>systemnames;

    double mintime=100000.0;
    double maxtime=0.0;
    float timestep = psr[0].AverageEpochWidth;
    char* flagName = psr[0].AverageFlag;
    fprintf(stderr, "width %.3e flag %s\n", timestep, flagName);

    for(int o=0;o<psr[0].nobs;o++){
        int found=0;
        if((double)psr[0].obsn[o].bat > maxtime){maxtime=(double)psr[0].obsn[o].bat;}
        if((double)psr[0].obsn[o].bat < mintime){mintime=(double)psr[0].obsn[o].bat;}
        for (int f=0;f<psr[0].obsn[o].nFlags;f++){
            if(strcasecmp(psr[0].obsn[o].flagID[f],flagName)==0){

	      //fprintf(stderr, "found flag\n");                
                if(std::find(systemnames.begin(), systemnames.end(), psr[0].obsn[o].flagVal[f]) != systemnames.end()) {
                } else {

                    systemnames.push_back(psr[0].obsn[o].flagVal[f]);
                    systemcount++;
                }
                found=1;
            }
        }

        if(found==0){
            printf("Observation %i is missing the -f flag, please check before continuing\n",o);
        }
    }

    mintime=floor(mintime);
    maxtime=ceil(maxtime);


    int numbins=ceil((maxtime-mintime)/timestep);
    double **AverageRes = new double*[systemcount];
    double **AverageWeight = new double*[systemcount];
    double **AverageBat = new double*[systemcount];
    for(int i = 0; i < systemcount; i++){
        AverageRes[i] = new double[numbins];
        AverageWeight[i] = new double[numbins];
        AverageBat[i] = new double[numbins];
        for(int j = 0; j < numbins; j++){
            AverageRes[i][j] = 0;
            AverageWeight[i][j] = 0;
            AverageBat[i][j] = 0;
        }
    }




    for(int o=0;o<psr[0].nobs;o++){

        for (int f=0;f<psr[0].obsn[o].nFlags;f++){
            if(strcasecmp(psr[0].obsn[o].flagID[f],flagName)==0){

                int flagindex=-1;
                for(int s = 0; s < static_cast<int>(systemnames.size()); s++){
                    if(strcasecmp(psr[0].obsn[o].flagVal[f],systemnames[s].c_str())==0){
                        flagindex = s;
                    }
                }

                int bin = floor(((double)psr[0].obsn[o].bat-mintime)/timestep);
                double adjustedErr = pow(psr[0].obsn[o].toaErr*pow(10.0, -6), 2);
                AverageWeight[flagindex][bin] += 1.0/adjustedErr;
                
		double resid=psr[0].obsn[o].residual;
		//
		if (psr[0].TNsubtractRed == 1)
		  {
		    resid -= (double) psr[0].obsn[o].TNRedSignal;
		  }

		if (psr[0].TNsubtractDM ==1)
		  {
		    resid -= (double) psr[0].obsn[o].TNDMSignal;
		  }
		if (psr[0].TNsubtractChrom ==1)
		  {
		    resid -= (double) psr[0].obsn[o].TNChromSignal;
		  }

		AverageRes[flagindex][bin] += resid/adjustedErr;
                AverageBat[flagindex][bin] += (double)psr[0].obsn[o].bat/adjustedErr;
		

	    }
        }
    }

    for(int o=0;o<psr[0].nobs;o++){


        //Is there an ECORR??
        double EcorrVal=0;
        for (int j=0;j<psr->obsn[o].nFlags;j++){
            for (int k=0;k<psr->nTNECORR;k++){
                if (strcmp(psr->obsn[o].flagID[j], psr->TNECORRFlagID[k])==0){
                    if (strcmp(psr->obsn[o].flagVal[j],psr->TNECORRFlagVal[k])==0){
                        EcorrVal=psr->TNECORRVal[k]*pow(10.0, -6);
                        //printf("Ecorr: %i %g \n", o, EcorrVal);
                    }
                }
            }
        }
    
        // SECORR
        

        double SEcorrVal=0;

        for (int j=0;j<psr->obsn[o].nFlags;j++){
            for (int k=0;k<psr->nTNSECORR;k++){
                if (strcmp(psr->obsn[o].flagID[j], psr->TNSECORRFlagID[k])==0){
                    if (strcmp(psr->obsn[o].flagVal[j],psr->TNSECORRFlagVal[k])==0){
		      SEcorrVal=psr->TNSECORRVal[k]*pow(10.0, -6)/sqrt(psr->obsn[o].tobs/3600.);
		      //printf("SEcorr: %i %g \n", o, SEcorrVal);
                    }
                }
            }
        }

        for (int f=0;f<psr[0].obsn[o].nFlags;f++){
            if(strcasecmp(psr[0].obsn[o].flagID[f],flagName)==0){
                int flagindex=-1;
                for(int s = 0; s < static_cast<int>(systemnames.size()); s++){
                    if(strcasecmp(psr[0].obsn[o].flagVal[f],systemnames[s].c_str())==0){
                        flagindex = s;
                    }
                }
                int bin = floor(((double)psr[0].obsn[o].bat-mintime)/timestep);


                psr[0].obsn[o].averagebat = AverageBat[flagindex][bin]/AverageWeight[flagindex][bin];; 
                psr[0].obsn[o].averageres  = AverageRes[flagindex][bin]/AverageWeight[flagindex][bin];
                psr[0].obsn[o].averageerr = sqrt(1.0/AverageWeight[flagindex][bin]  + pow(EcorrVal, 2)+pow(SEcorrVal,2));
            }
        }
    }

    FILE *avgfile;
    avgfile=fopen("avg.dat", "w");
    for(int o=0;o<psr[0].nobs;o++){
        fprintf(avgfile,  "%.9lf %.3e %3e \n", psr[0].obsn[o].averagebat, psr[0].obsn[o].averageres, psr[0].obsn[o].averageerr);
    }
    fclose(avgfile);


    for(int i = 0; i < systemcount; i++){
        delete[] AverageRes[i];
        delete[] AverageWeight[i];
        delete[] AverageBat[i];
    }
    delete[] AverageRes;
    delete[] AverageWeight;
    delete[] AverageBat;


}


void formResiduals(pulsar *psr,int npsr,int removeMean)
{

    longdouble residual;  /* Residual in phase */
    longdouble nphase,phase2,phase3,phase4,phase2state,lastResidual=0,priorResidual=0;
    longdouble *phase5 = new longdouble[MAX_OBSN];
    longdouble lastBat=0.0,priorBat=0.0;
    longdouble phaseJ,phaseW,phaseShape;
    longdouble ftpd,fct,ff0,phaseint;
    longdouble torb,deltaT,dt00=0.0,phas1=0.0;
    longdouble mean,tnmean, ct00=0.0;
    int dtm1s=0;
    int nmean,ntnmean;
    int ntpd,nf0;
    int i,p,k,l;
    int time=0;
    int ntrk=0;
    int gotit=0;
    long long pn0=-1;
    long long pnAdd=-1;

    // for gwecc Earth
    double prev_p;
    double prev_e;
    double prev_epoch;
    double prev_theta;
    double prev_a;
    int coalesceFlag = 0;

    // for gwecc pulsar
    double prev_p_p;
    double prev_e_p;
    double prev_epoch_p;
    double prev_theta_p;
    double prev_a_p;
    int coalesceFlag_p = 0;

    const char *CVS_verNum = "$Id$";



    if (displayCVSversion == 1) CVSdisplayVersion("formResiduals.C","formResiduals()",CVS_verNum);
    logtchk("Enter formresiduals()");

    for (p=0;p<npsr;p++)
    {
        mean = longdouble(0.0);
        nmean = 0;

        tnmean=longdouble(0.0);
        ntnmean=0;

        if(psr[p].refphs==REFPHS_TZR){
            // reinstate the extra TZR observation so we can compute the reference phase
            if (psr[p].nobs==MAX_OBSN){
                logerr("Error, need at least 1 spare observation to reference to REFPHS TZR");
                exit(1);
            }
            psr[p].nobs++;
            memcpy(&(psr[p].obsn[psr[p].nobs-1]),&(psr[p].tzrobs),sizeof(observation));
        }

        for (i=0;i<psr[p].nobs;i++)
        {

            torb = 0.0; 

            /* Binary parameters */
            if (psr[p].param[param_pb].paramSet[0]==1 || psr[p].param[param_fb].paramSet[0]==1) 
            {
                if (strcmp(psr[p].binaryModel,"BT")==0)         torb = BTmodel(psr,p,i,-1);
                else if (strcmp(psr[p].binaryModel,"BTJ")==0)   torb = BTJmodel(psr,p,i,-1,0);
                else if (strcmp(psr[p].binaryModel,"BTX")==0)   torb = BTXmodel(psr,p,i,-1,0);
                else if (strcmp(psr[p].binaryModel,"ELL1")==0)  torb = ELL1model(psr,p,i,-1,0);
                else if (strcmp(psr[p].binaryModel,"DD")==0)    torb = DDmodel(psr,p,i,-1);
                else if (strcmp(psr[p].binaryModel,"DDK")==0)   torb = DDKmodel(psr,p,i,-1);
                else if (strcmp(psr[p].binaryModel,"DDS")==0)   torb = DDSmodel(psr,p,i,-1);
                else if (strcmp(psr[p].binaryModel,"MSS")==0)   torb = MSSmodel(psr,p,i,-1);
                else if (strcmp(psr[p].binaryModel,"DDGR")==0)  torb = DDGRmodel(psr,p,i,-1);
                else if (strcmp(psr[p].binaryModel,"T2")==0)    torb = T2model(psr,p,i,-1,0);
                else if (strcmp(psr[p].binaryModel,"T2-PTA")==0) torb = T2_PTAmodel(psr,p,i,-1,0);
                else if( strcmp( psr[p].binaryModel, "DDH" ) == 0) torb = DDHmodel( psr, p, i, -1 );
                else if( strcmp( psr[p].binaryModel, "ELL1H" ) == 0) torb = ELL1Hmodel( psr, p, i, -1 );
		else if (strcmp( psr[p].binaryModel, "ELL1k" ) ==0)  torb = ELL1kmodel( psr, p, i, -1);
                else {printf("Warning: Unknown binary model '%s'\n",psr[p].binaryModel); exit(1);}
            }
            psr[p].obsn[i].torb = torb; // save for plotting etc
            deltaT = (psr[p].obsn[i].bbat-psr[p].param[param_pepoch].val[0])*86400.0+torb;
            fct = psr[p].obsn[i].bbat-(int)psr[p].obsn[i].bbat - (psr[p].param[param_pepoch].val[0] - 
                    (int)psr[p].param[param_pepoch].val[0]);	   
            ftpd = fct+torb/86400.0;	   
            psr[p].obsn[i].pet = psr[p].obsn[i].bbat - torb/SECDAY;

            /* What is the extra -ntpd*(int)(psr.f0) term on the end??? */
            /* The extra (not included integral parts (int)f0.(int)bat) are not necessary! */
            /* ntpd = (int)psr[p].obsn[i].bat-(int)psr[p].param[param_pepoch].val[0]; */
            /*	   phase2 = deltaT*(psr[p].param[param_f0].val)-(ntpd*(int)(psr[p].param[param_f0].val)*86400.0); */

            nf0  = (int)psr[p].param[param_f].val[0];
            ntpd = ((int)psr[p].obsn[i].bbat-(int)psr[p].param[param_pepoch].val[0]);
            ff0  = psr[p].param[param_f].val[0]-(int)psr[p].param[param_f].val[0];
            phase2 = (nf0*ftpd+ntpd*ff0+ftpd*ff0); 
            phase2 *= 86400.0;
            /* redwards, all these calls to pow are slow & imprecise. changed */
            longdouble arg = deltaT*deltaT;	   
            phase3 = 0.5*psr[p].param[param_f].val[1]*arg; 
            arg *= deltaT; if (psr[p].param[param_f].paramSet[2]==1) phase3 += (psr[p].param[param_f].val[2]/longdouble(6.0))*arg;
            arg *= deltaT; if (psr[p].param[param_f].paramSet[3]==1) phase3 += (psr[p].param[param_f].val[3]/longdouble(24.0))*arg;
            arg *= deltaT; if (psr[p].param[param_f].paramSet[4]==1) phase3 += (psr[p].param[param_f].val[4]/longdouble(120.0))*arg;
            arg *= deltaT; if (psr[p].param[param_f].paramSet[5]==1) phase3 += (psr[p].param[param_f].val[5]/longdouble(720.0))*arg;
            arg *= deltaT; if (psr[p].param[param_f].paramSet[6]==1) phase3 += (psr[p].param[param_f].val[6]/longdouble(5040.0))*arg;
            arg *= deltaT; if (psr[p].param[param_f].paramSet[7]==1) phase3 += (psr[p].param[param_f].val[7]/longdouble(40320.0))*arg;
            arg *= deltaT; if (psr[p].param[param_f].paramSet[8]==1) phase3 += (psr[p].param[param_f].val[8]/longdouble(362880.0))*arg;
            arg *= deltaT; if (psr[p].param[param_f].paramSet[9]==1) phase3 += (psr[p].param[param_f].val[9]/longdouble(3628800.0))*arg; 
            arg *= deltaT; if (psr[p].param[param_f].paramSet[10]==1) phase3 += (psr[p].param[param_f].val[10]/longdouble(39916800.0))*arg; 
            arg *= deltaT; if (psr[p].param[param_f].paramSet[11]==1) phase3 += (psr[p].param[param_f].val[11]/longdouble(479001600.0))*arg; 
            arg *= deltaT; if (psr[p].param[param_f].paramSet[12]==1) phase3 += (psr[p].param[param_f].val[12]/longdouble(6227020800.0))*arg; 

            // Must check for two-state changes
            phase2state = 0.0;
            for (k=0;k<psr[p].param[param_stateSwitchT].aSize;k++)
            {
                if (psr[p].param[param_stateSwitchT].paramSet[k]==1)
                {
                    longdouble tp,dt1,tgl;
                    int signV=1;
                    tp = (ntpd+ftpd)*86400.0;	
                    tgl = (psr[p].param[param_stateSwitchT].val[k] - psr[p].param[param_pepoch].val[0])*86400.0; 
                    if (tp >= tgl) 
                    {
                        dt1 = tp-tgl;
                        if (k%2 == 0)
                            signV=1;
                        else
                            signV=-1;
                        phase2state+=(0.5*(psr[p].param[param_df1].val[0]*signV)*dt1*dt1); 
                        //		       printf("Using sign: %d %d %g\n",k,signV,(double)phase2state);

                    }		   
                }
            }

            // add in term from pulsar braking


            if ((psr[p].param[param_brake].paramSet[0] == 1) && (psr[p].param[param_f].paramSet[1] ==1)) 
            {
                longdouble F2brake;
                longdouble F3brake;
                longdouble F1,F0;
                longdouble bindex;
                longdouble arg3;
                longdouble arg4;
                arg3= deltaT*deltaT*deltaT;
                arg4 = arg3*deltaT;

                F1 = psr[p].param[param_f].val[1];
                F0 =  psr[p].param[param_f].val[0];
                bindex= psr[p].param[param_brake].val[0];
                // f1squared= psr[p].param[param_f].val[1];

                F2brake= bindex*F1*F1/F0;
                F3brake= bindex*(2*bindex-1)*F1*F1*F1/F0/F0;

                //arg = deltaT*deltaT*deltaT;

                phase3 += F2brake*arg3/6.L + F3brake*arg4/24.L;

            }



            /* Must check glitch parameters */
            phase4 = 0.0; 
            for (k=0;k<psr[p].param[param_glep].aSize;k++)
            {
                if (psr[p].param[param_glep].paramSet[k]==1)
                {
                    longdouble tp,dt1,expf,tgl;
                    tp = (ntpd+ftpd)*86400.0;	
                    tgl = (psr[p].param[param_glep].val[k] - psr[p].param[param_pepoch].val[0])*86400.0; // - psr[p].param[param_pepoch].val[0])*86400.0;
                    /* Only has effect after the glitch epoch */
                    if (tp >= tgl) 
                    {
                        double expf2,expf3;
                        dt1 = tp-tgl;
                        if (psr[p].param[param_gltd].val[k]!=0.0)
                            expf = exp(-dt1/86400.0/psr[p].param[param_gltd].val[k]);
                        else
                            expf = 1.0;

                        if (psr[p].param[param_gltd2].val[k]!=0.0)
                            expf2 = exp(-dt1/86400.0/psr[p].param[param_gltd2].val[k]);
                        else
                            expf2 = 1.0;

                        if (psr[p].param[param_gltd3].val[k]!=0.0)
                            expf3 = exp(-dt1/86400.0/psr[p].param[param_gltd3].val[k]);
                        else
                            expf3 = 1.0;



                        // What happens if GLF2 (or GLF1) is not set?
                        phase4+=psr[p].param[param_glph].val[k]+
                            psr[p].param[param_glf0].val[k]*dt1 + 
                            0.5*psr[p].param[param_glf1].val[k]*dt1*dt1 +
                            1.0/6.0*psr[p].param[param_glf2].val[k]*dt1*dt1*dt1 +
                            psr[p].param[param_glf0d].val[k]*  // first exponential
                            psr[p].param[param_gltd].val[k]*86400.0*(1.0-expf) + 
                            psr[p].param[param_glf0d2].val[k]*  // second exponential
                            psr[p].param[param_gltd2].val[k]*86400.0*(1.0-expf2)+
                            psr[p].param[param_glf0d3].val[k]* // third exponential
                            psr[p].param[param_gltd3].val[k]*86400.0*(1.0-expf3);
                        //		       printf("Glitch phase = %10f\n",(double)(psr[p].param[param_glph].val[k]+
                        //			 psr[p].param[param_glf0].val[k]*dt1 + 
                        //			 0.5*psr[p].param[param_glf1].val[k]*dt1*dt1 +
                        //			 1.0/6.0*psr[p].param[param_glf2].val[k]*dt1*dt1*dt1 +
                        //			 psr[p].param[param_glf0d].val[k]*
                        //			      psr[p].param[param_gltd].val[k]*86400.0*(1.0-expf)));
                    }
                }
            }
	    

    for (k=0;k<psr[p].param[param_expep].aSize;k++)
            {
                if (psr[p].param[param_expep].paramSet[k]==1)
		  {
		    long double freq, dt, dm, tau, gamma;
		    if (psr[p].param[param_expindex].paramSet[k] ==1)
		      {
			gamma=psr[p].param[param_expindex].val[k];
		      }
		    else
		      {
			gamma=-2;
		      }
		    freq=psr[p].obsn[i].freqSSB/1.4e9;
		    dt=(psr[p].obsn[i].bbat - psr[p].param[param_expep].val[k]);
		    dm=psr[p].param[param_expph].val[k];
		    tau=psr[p].param[param_exptau].val[k];
		    
		    //fprintf(stderr, "DT %.3Le TAU %.3Le\n", dt, tau);

		    if (dt  >0)
		      {
			phase4 +=  dm*powl(freq, gamma)*exp(-dt/tau)*psr[p].param[param_f].val[0];   
			//fprintf(stderr, "%.3Le\n", phase4);
		      }
		  
		      
		  }
	    }
    
    for(k=0;k<psr[p].param[param_gausep].aSize;k++)
    {
        if(psr[p].param[param_gausep].paramSet[k] ==1)
            {
                long double freq, dt, amp, sig, gamma,val;            
                // reference to 1.4 GHz  to agree with Enterprise   
                freq= psr[p].obsn[i].freqSSB/1.4e9; 
  
                dt=(psr[p].obsn[i].bbat - psr[p].param[param_gausep].val[k]);
                amp=psr[p].param[param_gausamp].val[k];
                sig=psr[p].param[param_gaussig].val[k];


                // if index is not set assume the Gaussian is achromatic
                 if (psr[p].param[param_gausindex].paramSet[k] ==1)
                    {
                        gamma=psr[p].param[param_gausindex].val[k];
                    }
                else{
                    gamma=0;
    
                    }
  
                    // model  is
                    // amp*powl(freq, gamma)*exp(-dt*dt/2./sig/sig);

                val= amp*powl(freq, gamma)*exp(-dt*dt/2./sig/sig);     
                phase4 += val*psr[p].param[param_f].val[0];            
            }

    }           



	    


            if(psr[p].param[param_glep].aSize>0){
                //printf("Glitch %i %g %Lg %Lg\n", i, (double) psr[p].obsn[i].bat, phase4/psr[p].param[param_f].val[0], (1.0/psr[p].param[param_f].val[0]) );
            } 
            /* Add in extra phase due to jumps */
            phaseJ = 0.0;
            for (k=1;k<=psr[p].nJumps;k++)	    
            {
                for (l=0;l<psr[p].obsn[i].obsNjump;l++)
                {		
                    if (psr[p].obsn[i].jump[l]==k && psr[p].jumpSAT[l]==0)
                        phaseJ+=psr[p].jumpVal[k]*psr[p].param[param_f].val[0];
                }
            }
    
            // add in the newfdjumps	    
            for (k=1;k<=psr[p].nfdJumps;k++)	    
            {
                for (l=0;l<psr[p].obsn[i].obsNfdjump;l++)
                {		
                    if (psr[p].obsn[i].fdjump[l]==k)
                    {
                        int idx;
                        idx=psr[p].fdjumpIdx[k];
                        //fprintf(stderr, "%d %d %.5le  %.5le  %.5le\n",idx,k, psr[p].fdjumpVal[k],pow(psr[p].obsn[i].freqSSB/1e9,idx),psr[p].param[param_f].val[0]);
                        //phaseJ-=psr[p].fdjumpVal[k]*pow(log( psr[p].obsn[i].freqSSB/1e9),idx)*psr[p].param[param_f].val[0];
                   
                        phaseJ-=psr[p].fdjumpVal[k]*pow( psr[p].obsn[i].freqSSB/1e9,idx)*psr[p].param[param_f].val[0];
                    }
                }
            }

            /* Red Shapelet Events (M. Keith 2023) */
            phaseShape=0;
            for (int iTNShape=0; iTNShape < psr[p].nTNShapeletEvents; ++iTNShape) {
                // We only want to do ones with a spectral index equal to zero here
                if (psr[p].TNShapeletEvFScale[iTNShape] == 0.0) {
                    // I think that the red shape is a time term...
                    phaseShape -= evaluateShapelet(psr->TNShapeletEvN[iTNShape],
                            psr->TNShapeletEvPos[iTNShape],
                            psr->TNShapeletEvWidth[iTNShape],
                            psr->TNShapeletEvCoef[iTNShape],
                            (double)psr[p].obsn[i].bat)*psr[p].param[param_f].val[0];
                }
            }





            /* Add in extra phase due to whitening procedures */
            phaseW = 0.0;
            if (psr[p].param[param_wave_om].paramSet[0] == 1)
            {
                double om;  /* Fundamental frequency */
                double dt;  /* Change in time from pepoch */
                double om_eff;

                dt = psr[p].obsn[i].bbat - psr[p].param[param_waveepoch].val[0] ;
                om = psr[p].param[param_wave_om].val[0];
                for (k=0;k<psr[p].nWhite;k++)
                {
                    if (psr[p].waveScale==0)
                    {	
                        if( k==0)
                        {
                            om_eff = om*(k+1.0);
                        }
                        else
                        {
                            om_eff = om*(k+1.0);
                        }

                        phaseW += psr[p].wave_sine[k]*sin(om_eff*dt)*psr[p].param[param_f].val[0] +
                            psr[p].wave_cos[k]*cos(om_eff*dt)*psr[p].param[param_f].val[0];
                    }
                    else if(psr[p].waveScale ==1)
                    {
                        double freq= psr[p].obsn[i].freq;
                        phaseW += (psr[p].wave_sine[k]*sin(om*(k+1)*dt)*psr[p].param[param_f].val[0] +
                                psr[p].wave_cos[k]*cos(om*(k+1)*dt)*psr[p].param[param_f].val[0])/(DM_CONST*freq*freq);
                    }
                    else if (psr[p].waveScale ==2)
                    {
                        double freq= psr[p].obsn[i].freq;
                        phaseW+= psr[p].wave_sine[k]*sin(om*(k+1)*dt)*psr[p].param[param_f].val[0] +
                            psr[p].wave_cos[k]*cos(om*(k+1)*dt)*psr[p].param[param_f].val[0];
                        phaseW += (psr[p].wave_sine_dm[k]*sin(om*(k+1)*dt)*psr[p].param[param_f].val[0] +
                                psr[p].wave_cos_dm[k]*cos(om*(k+1)*dt)*psr[p].param[param_f].val[0])/(DM_CONST*freq*freq); 

                    }


                }
            }

            if (psr[p].param[param_wave_dm].paramSet[0] == 1)
            {
                double om;  /* Fundamental frequency */
                double dt;  /* Change in time from pepoch */

                dt = psr[p].obsn[i].bbat - psr[p].param[param_waveepoch_dm].val[0] ;
                om = psr[p].param[param_wave_dm].val[0];

                double freq= psr[p].obsn[i].freq;
                for (k=0;k<psr[p].nWhite_dm;k++)
                {



                    phaseW += (psr[p].wave_sine_dm[k]*sin(om*(k+1)*dt)*psr[p].param[param_f].val[0] +
                            psr[p].wave_cos_dm[k]*cos(om*(k+1)*dt)*psr[p].param[param_f].val[0])/(DM_CONST*freq*freq); 
                }



            }



            // Add in extra phase due to quadrupolar signal 
            if (psr[p].param[param_quad_om].paramSet[0] == 1)
            {
                double omega_g;
                //	       double res_e,res_i;
                longdouble resp,resc,res_r,res_i;
// UNUSED VARIABLE //                 double phi_g;
                double lambda_p,beta_p,lambda,beta;
// UNUSED VARIABLE //                 longdouble time;
                double n1,n2,n3;
                double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
                double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;
                double cosTheta;

// UNUSED VARIABLE //                 double om;  /* Fundamental frequency */
                double dt;  /* Change in time from pepoch */

                dt = (psr[p].obsn[i].bbat - psr[p].quadEpoch)*86400.0;

                if (psr[p].param[param_raj].paramSet[1] == 1)
                    lambda_p = (double)psr[p].param[param_raj].val[1];
                else
                    lambda_p = (double)psr[p].param[param_raj].val[0];

                if (psr[p].param[param_decj].paramSet[1] == 1)
                    beta_p   = (double)psr[p].param[param_decj].val[1];
                else
                    beta_p   = (double)psr[p].param[param_decj].val[0];
                lambda   = psr[p].quadRA;
                beta     = psr[p].quadDEC;

                // Pulsar vector
                n1 = cosl(lambda_p)*cosl(beta_p);
                n2 = sinl(lambda_p)*cosl(beta_p);
                n3 = sinl(beta_p);
                cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
                    sinl(beta)*sinl(beta_p);

                // From KJ's paper
                // Gravitational wave matrix
                e11p = pow(sinl(lambda),2)-pow(cosl(lambda),2)*pow(sinl(beta),2);
                e21p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
                e31p = cosl(lambda)*sinl(beta)*cosl(beta);

                e12p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
                e22p = pow(cosl(lambda),2)-pow(sinl(lambda),2)*pow(sinl(beta),2);
                e32p = sinl(lambda)*sinl(beta)*cosl(beta);

                e13p = cosl(lambda)*sinl(beta)*cosl(beta);
                e23p = sinl(lambda)*sinl(beta)*cosl(beta);
                e33p = -powl(cosl(beta),2);

                resp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
                        n2*(n1*e21p+n2*e22p+n3*e23p)+
                        n3*(n1*e31p+n2*e32p+n3*e33p));
                //	       ld_printf("resp = %Lg\n",resp);

                // Determine cross term
                e11c = sin(2*lambda)*sin(beta);
                e21c = -cos(2*lambda)*sin(beta);
                e31c = -sin(lambda)*cos(beta);

                e12c = -cos(2*lambda)*sin(beta);
                e22c = -sin(2*lambda)*sin(beta);
                e32c = cos(lambda)*cos(beta);

                e13c = -sin(lambda)*cos(beta);
                e23c = cos(lambda)*cos(beta);
                e33c  = 0;

                resc = (n1*(n1*e11c+n2*e12c+n3*e13c)+
                        n2*(n1*e21c+n2*e22c+n3*e23c)+
                        n3*(n1*e31c+n2*e32c+n3*e33c));

                for (k=0;k<psr[p].nQuad;k++)
                {
                    omega_g = (double)psr[p].param[param_quad_om].val[0]*(k+1);		   
                    res_r = (psr[p].quad_aplus_r[k]*resp+psr[p].quad_across_r[k]*resc)*sinl(omega_g*dt);
                    //		   res_i = (psr[p].quad_aplus_i[k]*resp+psr[p].quad_across_i[k]*resc)*(cos(omega_g*dt)-1);
                    res_i = (psr[p].quad_aplus_i[k]*resp+psr[p].quad_across_i[k]*resc)*(cosl(omega_g*dt));


                    if (psr[p].gwsrc_psrdist>0) // Add in the pulsar term
                    {
                        printf("Pulsar distance = %g\n",(double)psr[p].gwsrc_psrdist);
                        res_r -= (psr[p].quad_aplus_r[k]*resp+psr[p].quad_across_r[k]*resc)*sinl(omega_g*dt-(1-cosTheta)*psr[p].gwsrc_psrdist/SPEED_LIGHT*omega_g); 
                        //
                        res_i -= (psr[p].quad_aplus_i[k]*resp+psr[p].quad_across_i[k]*resc)*(cosl(omega_g*dt-(1-cosTheta)*psr[p].gwsrc_psrdist/SPEED_LIGHT*omega_g));
                    }
                    //		   printf("cosTheta = %g\n",(double)cosTheta);
                    if ((1-cosTheta)==0.0)
                    {
                        res_r = 0.0;
                        res_i = 0.0;
                    }
                    else
                    {
                        res_r = longdouble(1.0)/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta))*(res_r); 
                        res_i = longdouble(1.0)/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta))*(res_i); 
                        //		       res_r = longdouble(1.0)/(omega_g)*(res_r);
                        //res_i = longdouble(1.0)/(omega_g)*(res_i);
                    }
                    phaseW += (res_r+res_i)*psr[p].param[param_f].val[0];

                }
            } 


            // Add in eccentric binary gravitational wave signal (Vikram Ravi)
            if (psr[p].param[param_gwecc].paramSet[0]==1)
            {

                if (i==0) {

                    double ra_p = (double)psr[p].param[param_raj].val[0];
                    double dec_p   = -(double)psr[p].param[param_decj].val[0]+M_PI/2.;
                    double ra_g  = psr[p].gwecc_ra;
                    double dec_g     = psr[p].gwecc_dec;

                    printf("Have ra_p %g ra_g %g dec_p %g dec_g %g\n",ra_p,ra_g,dec_p,dec_g);

                    // cos(separation) between source and pulsar vectors
                    double cosMu = cosl(dec_g)*cosl(dec_p)*cosl(ra_g-ra_p)+sinl(dec_g)*sinl(dec_p);

                    prev_p=(psr[p].gwecc_orbital_period)*365.25*86400.; // READ IN YEARS!
                    prev_e=psr[p].gwecc_e;
                    prev_epoch = psr[p].gwecc_epoch; // READ IN MJD!
                    prev_theta = psr[p].gwecc_theta_0;
                    prev_a = pow(6.67e-11*(psr[p].gwecc_m1+psr[p].gwecc_m2)*1.9891e30*pow(2.*M_PI/prev_p,-2.),1./3.);

                    prev_p_p=(psr[p].gwecc_orbital_period)*365.25*86400.; // READ IN YEARS!
                    prev_e_p=psr[p].gwecc_e;
                    prev_epoch_p = psr[p].gwecc_epoch+(1./86400.)*(psr[p].gwecc_psrdist*3.08568e19/2.998e8)*(1.-cosMu); // READ IN MJD! // READ IN KPC
                    prev_theta_p = psr[p].gwecc_theta_0;
                    prev_a_p = pow(6.67e-11*(psr[p].gwecc_m1+psr[p].gwecc_m2)*1.9891e30*pow(2.*M_PI/prev_p,-2.),1./3.);

                }


                if (psr[p].gwecc_pulsarTermOn==0) {

                    printf("ONLY EARTH\n");

                    phaseW+=psr[p].param[param_gwecc].val[0]*(eccRes(&psr[p],i,&coalesceFlag,&prev_p,&prev_e,&prev_a,&prev_epoch,&prev_theta))*psr[p].param[param_f].val[0];
                }

                if (psr[p].gwecc_pulsarTermOn==1) {

                    printf("PULSAR + EARTH\n");

                    phaseW+=psr[p].param[param_gwecc].val[0]*(eccRes(&psr[p],i,&coalesceFlag,&prev_p,&prev_e,&prev_a,&prev_epoch,&prev_theta)-eccRes(&psr[p],i,&coalesceFlag_p,&prev_p_p,&prev_e_p,&prev_a_p,&prev_epoch_p,&prev_theta_p))*psr[p].param[param_f].val[0];

                }

                if (psr[p].gwecc_pulsarTermOn==2) {

                    printf("ONLY PULSAR\n");

                    phaseW+=psr[p].param[param_gwecc].val[0]*(-eccRes(&psr[p],i,&coalesceFlag_p,&prev_p_p,&prev_e_p,&prev_a_p,&prev_epoch_p,&prev_theta_p))*psr[p].param[param_f].val[0];
                }

            }



            /* Add in extra phase due to gravitational wave signal */
            if (psr[p].param[param_gwsingle].paramSet[0]==1 || psr[p].param[param_cgw].paramSet[0]==1)
            {
                double omega_g;
                //	       double res_e,res_i;
                longdouble resp,resc,res_r,res_i;
// UNUSED VARIABLE //                 double phi_g;
                double lambda_p,beta_p,lambda,beta;
                longdouble time;	      
                double n1,n2,n3;
                double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
                double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;
                double cosTheta;

                time    = (psr[p].obsn[i].bbat - psr[p].gwsrc_epoch)*longdouble(86400.0);

                lambda_p = (double)psr[p].param[param_raj].val[0];
                beta_p   = (double)psr[p].param[param_decj].val[0];
                lambda   = psr[p].gwsrc_ra;
                //	       beta     = M_PI/2.0-psr[p].gwsrc_dec;
                beta     = psr[p].gwsrc_dec;
                //			       phi_g   = psr[p].gwsrc_ra;

                // Pulsar vector
                n1 = cosl(lambda_p)*cosl(beta_p);
                n2 = sinl(lambda_p)*cosl(beta_p);
                n3 = sinl(beta_p);
                //	       printf("n = %g %g %g\n",n1,n2,n3);
                cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
                    sinl(beta)*sinl(beta_p);
                //	       printf("cosTheta = %g\n",cosTheta);

                // From KJ's paper
                // Gravitational wave matrix
                e11p = pow(sinl(lambda),2)-pow(cosl(lambda),2)*pow(sinl(beta),2);
                e21p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
                e31p = cosl(lambda)*sinl(beta)*cosl(beta);
                //	       printf("ex1p = %g %g %g (%g %g)\n",e11p,e21p,e31p,lambda,beta);

                e12p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
                e22p = pow(cosl(lambda),2)-pow(sinl(lambda),2)*pow(sinl(beta),2);
                e32p = sinl(lambda)*sinl(beta)*cosl(beta);
                //	       printf("ex2p = %g %g %g\n",e12p,e22p,e32p);

                e13p = cosl(lambda)*sinl(beta)*cosl(beta);
                e23p = sinl(lambda)*sinl(beta)*cosl(beta);
                e33p = -powl(cosl(beta),2);
                //	       printf("ex3p = %g %g %g\n",e13p,e23p,e33p);
                //	       exit(1);


                if (psr[p].param[param_gwsingle].paramSet[0]==1)
                    omega_g = (double)psr[p].param[param_gwsingle].val[0];
                else if (psr[p].param[param_cgw].paramSet[0]==1)
                    omega_g = (double)psr[p].param[param_cgw].val[0];

                resp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
                        n2*(n1*e21p+n2*e22p+n3*e23p)+
                        n3*(n1*e31p+n2*e32p+n3*e33p));
                //	       ld_printf("resp = %Lg\n",resp);

                // Determine cross term
                e11c = sin(2*lambda)*sin(beta);
                e21c = -cos(2*lambda)*sin(beta);
                e31c = -sin(lambda)*cos(beta);

                e12c = -cos(2*lambda)*sin(beta);
                e22c = -sin(2*lambda)*sin(beta);
                e32c = cos(lambda)*cos(beta);

                e13c = -sin(lambda)*cos(beta);
                e23c = cos(lambda)*cos(beta);
                e33c  = 0;


                resc = (n1*(n1*e11c+n2*e12c+n3*e13c)+
                        n2*(n1*e21c+n2*e22c+n3*e23c)+
                        n3*(n1*e31c+n2*e32c+n3*e33c));

                if (psr[p].param[param_gwsingle].paramSet[0]==1)
                {
                    res_r = (psr[p].gwsrc_aplus_r*resp+psr[p].gwsrc_across_r*resc)*sinl(omega_g*time);
                    res_i = (psr[p].gwsrc_aplus_i*resp+psr[p].gwsrc_across_i*resc)*(cosl(omega_g*time));
                    //	       res_i = (psr[p].gwsrc_aplus_i*resp+psr[p].gwsrc_across_i*resc)*(cos(omega_g*time)-1);
                    if (psr[p].gwsrc_psrdist>0) // Add in the pulsar term
                    {
                        res_r += (psr[p].gwsrc_aplus_r*resp+psr[p].gwsrc_across_r*resc)*sinl(omega_g*time-(1-cosTheta)*psr[p].gwsrc_psrdist/SPEED_LIGHT*omega_g);
                        //		   res_i += (psr[p].gwsrc_aplus_i*resp+psr[p].gwsrc_across_i*resc)*(cos(omega_g*time-(1-cosTheta)*psr[p].gwsrc_psrdist/SPEED_LIGHT*omega_g)-1);
                        res_i += (psr[p].gwsrc_aplus_i*resp+psr[p].gwsrc_across_i*resc)*(cosl(omega_g*time-(1-cosTheta)*psr[p].gwsrc_psrdist/SPEED_LIGHT*omega_g));
                    }
                    if ((1-cosTheta)==0.0)
                    {
                        res_r = 0.0;
                        res_i = 0.0;
                    }
                    else
                    {
                        res_r = longdouble(1.0)/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta))*(res_r); 
                        res_i = longdouble(1.0)/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta))*(res_i); 
                    }
                }
                else if (psr[p].param[param_cgw].paramSet[0]==1)
                {
                    res_r = (psr[p].cgw_h0/omega_g*((1+pow(psr[p].cgw_cosinc,2))*cos(2*psr[p].cgw_angpol)*sin(omega_g*time)+2*psr[p].cgw_cosinc*sin(2*psr[p].cgw_angpol)*cos(omega_g*time)))*resp 
                        + (psr[p].cgw_h0/omega_g*((1+pow(psr[p].cgw_cosinc,2))*sin(2*psr[p].cgw_angpol)*sin(omega_g*time)-2*psr[p].cgw_cosinc*cos(2*psr[p].cgw_angpol)*cos(omega_g*time)))*resc; 
                    //(psr[p].gwsrc_aplus_r*resp+psr[p].gwsrc_across_r*resc)*sinl(omega_g*time);
                    res_i = 0.0;
                    printf("res_r = %g %g %g %g omega = %g %g %g %g time = %g\n",(double)res_r,(double)res_i,(double)resp,(double)resc,(double)omega_g,
                            (double)psr[p].cgw_h0,(double)psr[p].cgw_cosinc,(double)psr[p].cgw_angpol,(double)time);
                    if (psr[p].gwsrc_psrdist>0) // Add in the pulsar term  (NOTE: using subtraction here)
                    {
                        double omega_prime_g;
                        double h0_prime;

                        if (psr[p].cgw_mc == 0) {omega_prime_g = omega_g; h0_prime = psr[p].cgw_h0;}
                        else {
                            //			 omega_prime_g = omega_g - 2*M_PI*2.77e-8*pow(psr[p].cgw_mc/1e8,5.0/3.0)*pow(omega_g/2.0/M_PI/1e-7,11.0/3.0)*(psr[p].gwsrc_psrdist/PCM/1000.0)*(1-cosTheta);

                            omega_prime_g = 2*M_PI*pow((1-cosTheta)*psr[p].gwsrc_psrdist/SPEED_LIGHT*256.0/5.0/pow(SPEED_LIGHT,5)*pow(M_PI,8.0/3.0)*pow(GM*psr[p].cgw_mc,5.0/3.0)+pow(omega_g/2.0/M_PI,-8.0/3.0),-3.0/8.0);

                            h0_prime = psr[p].cgw_h0*pow(omega_prime_g/omega_g,2.0/3.0);
                            printf("Using: omega_prime_g = %g, omega_g = %g, h0_prime = %g, h0 = %g\n",omega_prime_g,omega_g,h0_prime,psr[p].cgw_h0);
                        }
                        res_r -= ((h0_prime/omega_prime_g*((1+pow(psr[p].cgw_cosinc,2))*cos(2*psr[p].cgw_angpol)*sin(omega_prime_g*time-(1-cosTheta)*psr[p].gwsrc_psrdist/SPEED_LIGHT*omega_prime_g)+2*psr[p].cgw_cosinc*sin(2*psr[p].cgw_angpol)*cos(omega_prime_g*time-(1-cosTheta)*psr[p].gwsrc_psrdist/SPEED_LIGHT*omega_prime_g)))*resp 
                                + (h0_prime/omega_prime_g*((1+pow(psr[p].cgw_cosinc,2))*sin(2*psr[p].cgw_angpol)*sin(omega_prime_g*time-(1-cosTheta)*psr[p].gwsrc_psrdist/SPEED_LIGHT*omega_prime_g)-2*psr[p].cgw_cosinc*cos(2*psr[p].cgw_angpol)*cos(omega_prime_g*time-(1-cosTheta)*psr[p].gwsrc_psrdist/SPEED_LIGHT*omega_prime_g)))*resc); 
                    }

                    if ((1-cosTheta)==0.0)
                    {
                        res_r = 0.0;
                        res_i = 0.0;
                    }
                    else
                    {
                        res_r = -longdouble(1.0)/(longdouble(2.0)*(longdouble(1.0)-cosTheta))*(res_r); 
                        res_i = -longdouble(1.0)/(longdouble(2.0)*(longdouble(1.0)-cosTheta))*(res_i); 
                    }
                }

                phaseW += (res_r+res_i)*psr[p].param[param_f].val[0];
                //	       printf("Res = %g\n",(double)res);
            }

            // Jingbo Wang's method for GWM modelling
            // Uses an amplitude and phase
            // See below for method that uses two amplitudes and no phase
            //
            if (psr[p].param[param_gwm_amp].paramSet[0]==1 &&
                    psr[p].param[param_gwm_amp].paramSet[1]==0)
            {
// UNUSED VARIABLE //                 double omega_g;
                //	       double res_e,res_i;
// UNUSED VARIABLE //                 longdouble res_i;
// UNUSED VARIABLE //                 double phi_g;
                double lambda_p,beta_p,lambda,beta;
// UNUSED VARIABLE //                 longdouble time;	      
                double n1,n2,n3;
// UNUSED VARIABLE //                 double e33p;
// UNUSED VARIABLE //                 double e33c;
                double cosTheta;
                double g1,g2,g3;
                time    = (psr[p].obsn[i].bbat - psr[p].gwm_epoch)*longdouble(86400.0);

                if (psr[p].param[param_raj].paramSet[1] == 1)
                    lambda_p = (double)psr[p].param[param_raj].val[1];
                else
                    lambda_p = (double)psr[p].param[param_raj].val[0];

                if (psr[p].param[param_decj].paramSet[1] == 1)
                    beta_p   = (double)psr[p].param[param_decj].val[1];
                else
                    beta_p   = (double)psr[p].param[param_decj].val[0];
                lambda   = psr[p].gwm_raj;
                //	       beta     = M_PI/2.0-psr[p].gwsrc_dec;
                beta     = psr[p].gwm_decj;
                //			       phi_g   = psr[p].gwsrc_ra;

                // GW vector
                g1 = -cosl(lambda)*cosl(beta);
                g2 = -sinl(lambda)*cosl(beta);
                g3 = -sinl(beta);


                // Pulsar vector
                n1 = cosl(lambda_p)*cosl(beta_p);
                n2 = sinl(lambda_p)*cosl(beta_p);
                n3 = sinl(beta_p);
                //	       printf("n = %g %g %g\n",n1,n2,n3);
                cosTheta = -(cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p) + sinl(beta)*sinl(beta_p));

                //	       printf("cosTheta = %g\n",cosTheta);

                /* Only has effect after the glitch epoch */
                if (psr[p].obsn[i].sat >= psr[p].gwm_epoch)
                {
                    longdouble dt,scale;
                    double cos2Phi;
                    double cosPhi;
                    double l1,l2,l3,m1,m2,m3;
                    //double beta_m;
                    double d1,d2,d3,md;
                    double a1,a2,a3,ma;

                    /* define the GW coordinate system see Hobbs,G. (2009)*/
                    /* the d vector point to north pole or south pole*/ 

                    if (beta == 0.0 )
                    {
                        d1 = 0.0;
                        d2 = 0.0;
                        d3 = 1.0;
                    }

                    if ( beta > 0)
                    {
                        d1 = g1*cosl(0.5*M_PI - beta);
                        d2 = g2*cosl(0.5*M_PI - beta);
                        d3 = 1.0 + g3*cos(0.5*M_PI - beta);
                        md = sqrt(d1*d1 + d2*d2 + d3*d3);
                        d1 = d1/md;
                        d2 = d2/md;
                        d3 = d3/md;
                        /*covert d to unit vector */
                    } 
                    else if (beta < 0)
                    {
                        d1 = g1*cosl(-0.5*M_PI - beta);
                        d2 = g2*cosl(-0.5*M_PI - beta);
                        d3 = -1.0 + g3*cos(-0.5*M_PI - beta);
                        md = sqrt(d1*d1 + d2*d2 + d3*d3);
                        d1 = d1/md;
                        d2 = d2/md;
                        d3 = d3/md;
                    } 

                    //if (g2*d3-d2*g3 != 0)
                    // {
                    //  a1 = 1.0; 
                    //  a2 = (d1*g3-g1*d3)/(g2*d3-d2*g3);
                    //  a3 = (g2*d1-g1*d2)/(g3*d2-g2*d3); 
                    // }
                    //else if (g1*d3-d1*g3 != 0)
                    // {
                    //  a1 = (g3*d2-d3*g2)/(g1*d3-g3*d1); 
                    //  a2 = 1.0;
                    //  a3 = (g1*d2-d1*g2)/(g3*d1-d1*d3);
                    // }
                    //else if (d2*g1-g2*d1 != 0)			
                    // {
                    //  a1 = (g2*d3-d2*g3)/(d2*g1-g2*d1); 
                    //  a2 = (g1*d3-d1*g3)/(d1*g2-g1*d2);
                    //  a3 =1.0; 
                    // }
                    a1 =  (d2*g3-d3*g2);
                    a2 =  (d3*g1-d1*g3);
                    a3 =  (d1*g2-d2*g1);
                    /* conver it to unit vector */
                    ma = sqrt(a1*a1 +a2*a2 + a3*a3);
                    a1 = a1/ma;
                    a2 = a2/ma;
                    a3 = a3/ma;

                    /* polarisation vector of GW source */
                    m1 = d1*cosl(psr[p].gwm_phi)	+ a1*sinl(psr[p].gwm_phi);   
                    m2 = d2*cosl(psr[p].gwm_phi)	+ a2*sinl(psr[p].gwm_phi);
                    m3 = d3*cosl(psr[p].gwm_phi)	+ a3*sinl(psr[p].gwm_phi);

                    //////		   if  (g3 != 0) 
                    //		     {beta_m = atan2(-cos(beta)*cos(lambda-psr[p].gwm_phi),sin(beta));}
                    //		   else  
                    //	     {
                    //	       beta_m = atan2(sinl(psr[p].gwm_phi),cosl(psr[p].gwm_phi));
                    //	       psr[p].gwm_phi = lambda + 1.5708;
                    //     }	
                    //   m1 = cosl(psr[p].gwm_phi)*cosl(beta_m);
                    //   m2 = sinl(psr[p].gwm_phi)*cosl(beta_m);
                    //   m3 = sinl(beta_m);
                    if  (cosTheta != 1.0 && cosTheta != -1.0)
                    {g1 = g1*cosTheta; 
                        g2 = g2*cosTheta;
                        g3 = g3*cosTheta;

                        /*l is  the projection of pulsar vector on the plane which pependicular to the GW source direction */
                        l1 = n1 - g1;
                        l2 = n2 - g2;
                        l3 = n3 - g3;
                        cosPhi = (l1*m1 + l2*m2 + l3*m3)/sqrt(l1*l1 + l2*l2 + l3*l3);
                        //		        if  (cosPhi >= 1.0/sqrt(2.0))
                        cos2Phi = 2*cosPhi*cosPhi - 1.0;
                        //		       else
                        //		      	   cos2Phi = 2*sqrt(1.0 - cosPhi*cosPhi)*sqrt(1.0 - cosPhi*cosPhi) - 1.0;
                    }
                    else 
                    {cos2Phi = 0;}

                    dt = (psr[p].obsn[i].sat - psr[p].gwm_epoch)*86400.0;
                    scale = -0.5*cos2Phi*(1-cosTheta);
                    //		   scale=1.0;
                    printf("%s scale = %g %g %g\n",psr[p].name,(double)scale,(double)cos2Phi,(double)cosTheta);
                    // MUST CHECK CAREFULLY HOW TO USE GWM_DPHASE
                    phaseW += scale*(psr[p].param[param_gwm_amp].val[0]*psr[p].param[param_f].val[0])*dt  + psr[p].gwm_dphase;

                }
                //	       printf("Res = %g\n",(double)res);
            }


            // GWM model that includes two amplitudes
            if (psr[p].param[param_gwm_amp].paramSet[0]==1 &&
                    psr[p].param[param_gwm_amp].paramSet[1]==1)
            {
                if (psr[p].obsn[i].bbat >= psr[p].gwm_epoch)
                {
// UNUSED VARIABLE //                     longdouble t2;
                    longdouble dt;
// UNUSED VARIABLE //                     int k;
// UNUSED VARIABLE //                     double ival;
// UNUSED VARIABLE //                     double omega_g;
                    //	       double res_e,res_i;
                    longdouble resp,resc;
// UNUSED VARIABLE //                     double phi_g;
                    double lambda_p,beta_p,lambda,beta;
                    double n1,n2,n3;
                    double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
                    double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;
                    double cosTheta;

                    lambda_p = (double)psr[p].param[param_raj].val[0];
                    beta_p   = (double)psr[p].param[param_decj].val[0];
                    lambda   = psr[p].gwm_raj;
                    beta     = psr[p].gwm_decj;

                    // Pulsar vector
                    n1 = cosl(lambda_p)*cosl(beta_p);
                    n2 = sinl(lambda_p)*cosl(beta_p);
                    n3 = sinl(beta_p);

                    cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
                        sinl(beta)*sinl(beta_p);

                    // From KJ's paper
                    // Gravitational wave matrix

                    // NOTE: This is for the plus terms.  For cross should use different terms
                    e11p = pow(sinl(lambda),2)-pow(cosl(lambda),2)*pow(sinl(beta),2);
                    e21p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
                    e31p = cosl(lambda)*sinl(beta)*cosl(beta);

                    e12p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
                    e22p = pow(cosl(lambda),2)-pow(sinl(lambda),2)*pow(sinl(beta),2);
                    e32p = sinl(lambda)*sinl(beta)*cosl(beta);

                    e13p = cosl(lambda)*sinl(beta)*cosl(beta);
                    e23p = sinl(lambda)*sinl(beta)*cosl(beta);
                    e33p = -powl(cosl(beta),2);

                    resp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
                            n2*(n1*e21p+n2*e22p+n3*e23p)+
                            n3*(n1*e31p+n2*e32p+n3*e33p));

                    //		   printf("Resp = %s %g %g\n",psr[p].name,(double)resp,(double)cosTheta);
                    if ((1-cosTheta)==0.0)
                        resp = 0.0;  // Check if this is sensible
                    else
                        resp = longdouble(1.0)/(longdouble(2.0)*(longdouble(1.0)-cosTheta))*(resp); 

                    psr[p].quad_ifunc_geom_p = resp;
                    logdbg("Resp2 = %s %g %g",psr[p].name,(double)resp,(double)psr[p].quad_ifunc_geom_p);
                    // NOTE: These are for the cross terms. 
                    lambda   = psr[p].gwm_raj;
                    beta     = psr[p].gwm_decj;


                    e11c = sinl(2*lambda)*sinl(beta);
                    e21c = -cosl(2*lambda)*sinl(beta);
                    e31c = -sinl(lambda)*cosl(beta);

                    e12c = -cosl(2*lambda)*sinl(beta);
                    e22c = -sinl(2*lambda)*sinl(beta);
                    e32c = cosl(lambda)*cosl(beta);

                    e13c = -sinl(lambda)*cosl(beta);
                    e23c = cosl(lambda)*cosl(beta);
                    e33c  = 0;

                    resc = (n1*(n1*e11c+n2*e12c+n3*e13c)+
                            n2*(n1*e21c+n2*e22c+n3*e23c)+
                            n3*(n1*e31c+n2*e32c+n3*e33c));
                    //		   printf("Resc = %s %g %g\n",psr[p].name,(double)resc,(double)cosTheta);
                    if ((1-cosTheta)==0.0)
                        resc = 0.0;  // Check if this is sensible
                    else
                        resc = longdouble(1.0)/(longdouble(2.0)*(longdouble(1.0)-cosTheta))*(resc); 
                    psr[p].quad_ifunc_geom_c = resc;
                    // printf("Resc2 = %s %g %g %g %g\n",psr[p].name,(double)resc,(double)psr[p].quad_ifunc_geom_c,(double)cosTheta,(double)(longdouble(1.0)/(longdouble(2.0)*(longdouble(1.0)-cosTheta))));
                    dt = (psr[p].obsn[i].bbat - psr[p].gwm_epoch)*longdouble(86400.0);
                    //		   scale = -0.5*cos2Phi*(1-cosTheta);
                    phaseW += (psr[p].param[param_f].val[0]*dt*psr[p].param[param_gwm_amp].val[0]*psr[p].quad_ifunc_geom_p); 				
                    phaseW += (psr[p].param[param_f].val[0]*dt*psr[p].param[param_gwm_amp].val[1]*psr[p].quad_ifunc_geom_c); 				
                }
            }	 


	    // Gravitational wave signal caused by a single cosmic string burst
	    // Added in by G. Hobbs (20th Aug) based on equations from Naoyuki
            //
            if (psr[p].param[param_gwcs_amp].paramSet[0]==1)
            {
                double lambda_p,beta_p,lambda,beta;
		long double dt;
                double n1,n2,n3;
                double cosTheta;
                //double g1,g2,g3;
		double width,width_day;
		double extra1,extra2;

		longdouble resp,resc;
		double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
		double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;


		lambda_p = (double)psr[p].param[param_raj].val[0];
		beta_p   = (double)psr[p].param[param_decj].val[0];
		lambda   = psr[p].gwcs_raj;
		beta     = psr[p].gwcs_decj;
		
                    // Pulsar vector
		n1 = cosl(lambda_p)*cosl(beta_p);
		n2 = sinl(lambda_p)*cosl(beta_p);
		n3 = sinl(beta_p);
		
		cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
		  sinl(beta)*sinl(beta_p);
		
		// From KJ's paper
		// Gravitational wave matrix
		
		// NOTE: This is for the plus terms.  For cross should use different terms
		e11p = pow(sinl(lambda),2)-pow(cosl(lambda),2)*pow(sinl(beta),2);
		e21p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
		e31p = cosl(lambda)*sinl(beta)*cosl(beta);
		
		e12p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
		e22p = pow(cosl(lambda),2)-pow(sinl(lambda),2)*pow(sinl(beta),2);
		e32p = sinl(lambda)*sinl(beta)*cosl(beta);
		
		e13p = cosl(lambda)*sinl(beta)*cosl(beta);
		e23p = sinl(lambda)*sinl(beta)*cosl(beta);
		e33p = -powl(cosl(beta),2);

		resp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
			n2*(n1*e21p+n2*e22p+n3*e23p)+
			n3*(n1*e31p+n2*e32p+n3*e33p));
		
		//		   printf("Resp = %s %g %g\n",psr[p].name,(double)resp,(double)cosTheta);
		if ((1-cosTheta)==0.0)
		  resp = 0.0;  // Check if this is sensible
		else
		  resp = longdouble(1.0)/(longdouble(2.0)*(longdouble(1.0)-cosTheta))*(resp); 
		psr[p].gwcs_geom_p = resp;
		logdbg("Resp2 = %s %g %g",psr[p].name,(double)resp,(double)psr[p].quad_ifunc_geom_p);
		// NOTE: These are for the cross terms. 
		lambda   = psr[p].gwcs_raj;
		beta     = psr[p].gwcs_decj;
		

		e11c = sinl(2*lambda)*sinl(beta);
		e21c = -cosl(2*lambda)*sinl(beta);
		e31c = -sinl(lambda)*cosl(beta);
		
		e12c = -cosl(2*lambda)*sinl(beta);
		e22c = -sinl(2*lambda)*sinl(beta);
		e32c = cosl(lambda)*cosl(beta);
		
		e13c = -sinl(lambda)*cosl(beta);
		e23c = cosl(lambda)*cosl(beta);
		e33c  = 0;
		
		resc = (n1*(n1*e11c+n2*e12c+n3*e13c)+
			n2*(n1*e21c+n2*e22c+n3*e23c)+
			n3*(n1*e31c+n2*e32c+n3*e33c));
		//		   printf("Resc = %s %g %g\n",psr[p].name,(double)resc,(double)cosTheta);
		if ((1-cosTheta)==0.0)
		  resc = 0.0;  // Check if this is sensible
		else
		  resc = longdouble(1.0)/(longdouble(2.0)*(longdouble(1.0)-cosTheta))*(resc); 
		psr[p].gwcs_geom_c = resc;

		    		

		dt = (psr[p].obsn[i].sat - psr[p].gwcs_epoch)*86400.0;
		width = psr[p].gwcs_width*86400.0;
		width_day = psr[p].gwcs_width;
		
                /* Only has effect after the event epoch */
		if (psr[p].obsn[i].sat < psr[p].gwcs_epoch-width_day/2.0)
		  {
		    extra1=0;
		    extra2=0;
		  }
		else if (psr[p].obsn[i].sat <= psr[p].gwcs_epoch)
		  {		    
		    extra1 = (psr[p].param[param_gwcs_amp].val[0]*
			     (3.0/4.0*(pow(0.5*width,4.0/3.0)-pow(fabs(dt),4.0/3.0))-
			      pow(0.5*width,1.0/3.0)*(dt+0.5*width)));
		    extra2 = (psr[p].param[param_gwcs_amp].val[1]*
			     (3.0/4.0*(pow(0.5*width,4.0/3.0)-pow(fabs(dt),4.0/3.0))-
			      pow(0.5*width,1.0/3.0)*(dt+0.5*width)));
		  }
		else if (psr[p].obsn[i].sat <= psr[p].gwcs_epoch+width_day/2.0)
		  {
		    extra1 = (psr[p].param[param_gwcs_amp].val[0]*
			     (3.0/4.0*(pow(0.5*width,4.0/3.0)+pow(fabs(dt),4.0/3.0))-
			      pow(0.5*width,1.0/3.0)*(dt+0.5*width)));
		    extra2 = (psr[p].param[param_gwcs_amp].val[1]*
			     (3.0/4.0*(pow(0.5*width,4.0/3.0)+pow(fabs(dt),4.0/3.0))-
			      pow(0.5*width,1.0/3.0)*(dt+0.5*width)));
		  }
		else
		  {
		    extra1 = -0.25*(pow(0.5,1.0/3.0)*psr[p].param[param_gwcs_amp].val[0]*pow(width,4.0/3.0));
		    extra2 = -0.25*(pow(0.5,1.0/3.0)*psr[p].param[param_gwcs_amp].val[1]*pow(width,4.0/3.0));		    
		  }
		// Noting here that Ax = 0 --- this should become a phase term??
		phaseW += (psr[p].param[param_f].val[0]*(extra1*psr[p].gwcs_geom_p+extra2*psr[p].gwcs_geom_c));
            }


            /* Add in extra phase due to clock offset */
            if (psr[p].param[param_clk_offs].paramSet[0] == 1)
            {
                int j;
                printf("IN HERE SETTING THE CLOCKS %d\n",psr[p].clkOffsN);
                for (j=0;j<psr[p].clkOffsN-1;j++)
                {
                    if (psr[p].obsn[i].bbat >= psr[p].clk_offsT[j] &&
                            psr[p].obsn[i].bbat < psr[p].clk_offsT[j+1])
                    {
                        phaseW += (psr[p].clk_offsV[j]*psr[p].param[param_f].val[0]);
                        break;
                    }
                }
            }

            /* Add in extra phase due to interpolation */
            if (psr[p].param[param_ifunc].paramSet[0] == 1)
            {
                longdouble wi,t1;
                longdouble dt,speriod,tt;

                if (psr[p].param[param_ifunc].val[0] == 1) // Sinc interpolation
                {
                    t1=longdouble(0.0);
                    speriod = (longdouble)(psr[p].ifuncT[1]-psr[p].ifuncT[0]); 
                    //	       printf("ifuncN = %d\n",psr[p].ifuncN);
                    for (k=0;k<psr[p].ifuncN;k++)
                        //	       for (k=3;k<4;k++)
                    {
                        //		   printf("Have %g %g\n",psr[p].ifuncT[k],psr[p].ifuncV[k]);
                        dt = psr[p].obsn[i].bbat - (longdouble)psr[p].ifuncT[k];
                        wi=1;
                        if (dt==0)
                        {
                            t1 += wi*(longdouble)psr[p].ifuncV[k];
                        }
                        else
                        {
                            tt = M_PI/speriod*(dt);
                            t1 += wi*(longdouble)psr[p].ifuncV[k]*sinl(tt)/(tt);
                            //		       t2 += wi*sinl(tt)/(tt);
                        }
                    }
                    phaseW += (psr[p].param[param_f].val[0]*t1);
                }
                else if (psr[p].param[param_ifunc].val[0] == 2) // Linear interpolation
                {
                    double ival = ifunc(psr[p].ifuncT,psr[p].ifuncV,(double)psr[p].obsn[i].sat,psr[p].ifuncN);
                    phaseW += (psr[p].param[param_f].val[0]*ival);
                }
                else if (psr[p].param[param_ifunc].val[0] == 0) // No interpolation
                {
                    int k;
                    double ival;
                    for (k=0;k<psr[p].ifuncN-1;k++)
                    {
                        if ((double)psr[p].obsn[i].sat >= psr[p].ifuncT[k])
                        {
                            ival = psr[p].ifuncV[k];
                            break;
                        }
                    }
                    phaseW += (psr[p].param[param_f].val[0]*ival); 				
                }
            }

            /* plus term for Quad ifuncs*/
            if (psr[p].param[param_quad_ifunc_p].paramSet[0] == 1)
            {
                longdouble resp;
                double lambda_p,beta_p,lambda,beta;
                double n1,n2,n3;
                double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
                double cosTheta;

                if ((psr[p].quad_ifunc_geom_p == 0) && (psr[p].param[param_quad_ifunc_p].val[0] == 0 || psr[p].param[param_quad_ifunc_p].val[0] == 1 || psr[p].param[param_quad_ifunc_p].val[0] == 3)) 
                {
                    if (psr[p].param[param_quad_ifunc_p].val[0] == 3)
                        psr[p].quad_ifunc_geom_p = 1;
                    else
                    {
                        lambda_p = (double)psr[p].param[param_raj].val[0];
                        beta_p   = (double)psr[p].param[param_decj].val[0];
                        lambda   = psr[p].quad_ifunc_p_RA;
                        beta     = psr[p].quad_ifunc_p_DEC;

                        // Pulsar vector
                        n1 = cosl(lambda_p)*cosl(beta_p);
                        n2 = sinl(lambda_p)*cosl(beta_p);
                        n3 = sinl(beta_p);

                        cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
                            sinl(beta)*sinl(beta_p);

                        // From KJ's paper
                        // Gravitational wave matrix

                        // NOTE: This is for the plus terms.  For cross should use different terms
                        e11p = pow(sinl(lambda),2)-pow(cosl(lambda),2)*pow(sinl(beta),2);
                        e21p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
                        e31p = cosl(lambda)*sinl(beta)*cosl(beta);

                        e12p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
                        e22p = pow(cosl(lambda),2)-pow(sinl(lambda),2)*pow(sinl(beta),2);
                        e32p = sinl(lambda)*sinl(beta)*cosl(beta);

                        e13p = cosl(lambda)*sinl(beta)*cosl(beta);
                        e23p = sinl(lambda)*sinl(beta)*cosl(beta);
                        e33p = -powl(cosl(beta),2);

                        resp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
                                n2*(n1*e21p+n2*e22p+n3*e23p)+
                                n3*(n1*e31p+n2*e32p+n3*e33p));

                        if ((1-cosTheta)==0.0)
                            resp = 0.0;  // Check if this is sensible
                        else
                            resp = longdouble(1.0)/(longdouble(2.0)*(longdouble(1.0)-cosTheta))*(resp); 

                        psr[p].quad_ifunc_geom_p = resp;
                    }
                }

                double ival = ifunc(psr[p].quad_ifuncT_p,psr[p].quad_ifuncV_p,(double)psr[p].obsn[i].sat,psr[p].quad_ifuncN_p);
                phaseW += (psr[p].param[param_f].val[0]*ival*psr[p].quad_ifunc_geom_p);
            }

            // Cross term
            if (psr[p].param[param_quad_ifunc_c].paramSet[0] == 1)
            {
                longdouble resc;
                double lambda_p,beta_p,lambda,beta;
                double n1,n2,n3;
                double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;
                double cosTheta;

                if ((psr[p].quad_ifunc_geom_c == 0) && (psr[p].param[param_quad_ifunc_c].val[0] == 0 || psr[p].param[param_quad_ifunc_c].val[0] == 1 || psr[p].param[param_quad_ifunc_c].val[0] == 3)) 
                {
                    if (psr[p].param[param_quad_ifunc_c].val[0] == 3)
                        psr[p].quad_ifunc_geom_c = 1;
                    else
                    {
                        lambda_p = (double)psr[p].param[param_raj].val[0];
                        beta_p   = (double)psr[p].param[param_decj].val[0];
                        lambda   = psr[p].quad_ifunc_c_RA;
                        beta     = psr[p].quad_ifunc_c_DEC;

                        // Pulsar vector
                        n1 = cosl(lambda_p)*cosl(beta_p);
                        n2 = sinl(lambda_p)*cosl(beta_p);
                        n3 = sinl(beta_p);

                        cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
                            sinl(beta)*sinl(beta_p);

                        // From KJ's paper
                        // Gravitational wave matrix

                        // NOTE: These are for the cross terms. 
                        e11c = sin(2*lambda)*sin(beta);
                        e21c = -cos(2*lambda)*sin(beta);
                        e31c = -sin(lambda)*cos(beta);

                        e12c = -cos(2*lambda)*sin(beta);
                        e22c = -sin(2*lambda)*sin(beta);
                        e32c = cos(lambda)*cos(beta);

                        e13c = -sin(lambda)*cos(beta);
                        e23c = cos(lambda)*cos(beta);
                        e33c  = 0;

                        resc = (n1*(n1*e11c+n2*e12c+n3*e13c)+
                                n2*(n1*e21c+n2*e22c+n3*e23c)+
                                n3*(n1*e31c+n2*e32c+n3*e33c));

                        if ((1-cosTheta)==0.0)
                            resc = 0.0;  // Check if this is sensible
                        else
                            resc = longdouble(1.0)/(longdouble(2.0)*(longdouble(1.0)-cosTheta))*(resc); 
                        psr[p].quad_ifunc_geom_c = resc;
                    }
                }

                double ival = ifunc(psr[p].quad_ifuncT_c,psr[p].quad_ifuncV_c,(double)psr[p].obsn[i].sat,psr[p].quad_ifuncN_c);
                phaseW += (psr[p].param[param_f].val[0]*ival*psr[p].quad_ifunc_geom_c); 				
            }


            // Ryan's Geometrical fitting function for plus

            if ((psr[p].param[param_quad_ifunc_p].paramSet[0] ==1) &&  (psr[p].param[param_quad_ifunc_p].val[0] == 2))
            {


                double lp, bp;
                double lg, bg;


                if (psr[p].quad_ifunc_geom_p == 0)
                {



                    if(psr[p].simflag==1)
                    {
                        lp = psr[p].rasim;
                        bp = psr[p].decsim;

                    }
                    else
                    {
                        lp = (double)psr[p].param[param_raj].val[0];
                        bp   = (double)psr[p].param[param_decj].val[0];  
                    }

                    lg= psr[p].quad_ifunc_p_RA;
                    bg = psr[p].quad_ifunc_p_DEC; 

                    double px, py,pz,gx,gy,gz;

                    px = cos(bp)*cos(lp);
                    py = cos(bp)*sin(lp);
                    pz = sin(bp);

                    gx = cos(bg)*cos(lg);
                    gy = cos(bg)*sin(lg);
                    gz = sin(bg);






                    double ctheta;

                    ctheta = px*gx + py*gy + pz*gz;

                    // phi is the angle between the direction of principle polarization
                    // the projection of  the pulsar onto the plane of polarization



                    double rho, rhox, rhoy, rhoz;

                    rhox = px - ctheta*gx;
                    rhoy = py - ctheta*gy;
                    rhoz = pz - ctheta*gz;

                    double pol,polx, poly, polz;
                    // polarization vector is projection of vector pointing to NCP on plane


                    // plus polarization in \hat theta direction? 
                    // this ought to point to SCP

                    polx = sin(bg)*cos(lg);
                    poly = sin(bg)*sin(lg);
                    polz = -cos(bg);





                    rho = sqrt(rhox*rhox + rhoy*rhoy + rhoz*rhoz);
                    pol = sqrt(polx*polx + poly*poly+ polz*polz);




                    double phi, cphi,c2phi;
                    phi = acos((rhox*polx+ rhoy*poly+rhoz*polz)/(rho*pol));

                    cphi = (rhox*polx+rhoy*poly+rhoz*polz)/(rho*pol);
                    c2phi = 2*cphi*cphi-1;


                    fprintf(stderr, "%.3e\n", phi);


                    longdouble resp;

                    resp =0.5*(1-ctheta)*c2phi;


                    psr[p].quad_ifunc_geom_p = resp;
                }



                double ival = ifunc(psr[p].quad_ifuncT_p,psr[p].quad_ifuncV_p,(double)psr[p].obsn[i].sat,psr[p].quad_ifuncN_p);
                phaseW += (psr[p].param[param_f].val[0]*ival*psr[p].quad_ifunc_geom_p); 	


            }


            // Ryan's geometrical fitting function for cross

            if ((psr[p].param[param_quad_ifunc_c].paramSet[0] ==1) &&  (psr[p].param[param_quad_ifunc_c].val[0] == 2))
            {
                double lp, bp;
                double lg, bg;


                if (psr[p].quad_ifunc_geom_c == 0)
                {



                    if(psr[p].simflag==1)
                    {
                        lp = psr[p].rasim;
                        bp = psr[p].decsim;

                    }
                    else
                    {
                        lp = (double)psr[p].param[param_raj].val[0];
                        bp   = (double)psr[p].param[param_decj].val[0];  
                    }

                    lg= psr[p].quad_ifunc_c_RA;
                    bg = psr[p].quad_ifunc_c_DEC; 

                    double px, py,pz,gx,gy,gz;

                    px = cos(bp)*cos(lp);
                    py = cos(bp)*sin(lp);
                    pz = sin(bp);

                    gx = cos(bg)*cos(lg);
                    gy = cos(bg)*sin(lg);
                    gz = sin(bg);



                    double ctheta;

                    ctheta = px*gx + py*gy + pz*gz;

                    // phi is the angle between the direction of principle polarization
                    // direction of principle polarizatino of cross???
                    // the projection of  the pulsar onto the plane of polarization


                    double rho, rhox, rhoy, rhoz;

                    rhox = px - ctheta*gx;
                    rhoy = py - ctheta*gy;
                    rhoz = pz - ctheta*gz;


                    double pol, polx, poly, polz;
                    double thetax, thetay, thetaz;
                    double phix, phiy, phiz;

                    thetax = sin(bg)*cos(lg);
                    thetay = sin(bg)*sin(lg);
                    thetaz = -cos(bg);

                    phix = -sin(lg);
                    phiy = cos(lg);
                    phiz =0;

                    // polarization 45 degrees from theta towards phi?

                    polx = thetax+phix;
                    poly = thetay+phiy;
                    polz = thetaz+phiz;



                    pol = sqrt(polx*polx + poly*poly+ polz*polz);


                    rho = sqrt(rhox*rhox + rhoy*rhoy + rhoz*rhoz);

                    //fprintf(stderr, "%.3e %.3e\n", (polx*thetax+poly*thetay+polz*thetaz)/pol, pol);


                    //exit(0);



                    double  cphi, c2phi;
                    cphi = (rhox*polx + rhoy*poly + rhoz*polz)/rho/pol;
                    c2phi = 2*cphi*cphi-1;

                    longdouble resc;

                    resc =0.5*(1-ctheta)*c2phi;


                    psr[p].quad_ifunc_geom_c = resc;
                }



                double ival = ifunc(psr[p].quad_ifuncT_c,psr[p].quad_ifuncV_c,(double)psr[p].obsn[i].sat,psr[p].quad_ifuncN_c);
                phaseW += (psr[p].param[param_f].val[0]*ival*psr[p].quad_ifunc_geom_c); 	


            }




            // Ryan's burst fitting 
            if (psr[p].param[param_gwb_amp].paramSet[0]==1 &&
                    psr[p].param[param_gwb_amp].paramSet[1]==1)
            {

                //longdouble wi,t1,t2;
                longdouble dt, prefac;
// UNUSED VARIABLE //                 int k;
// UNUSED VARIABLE //                 double ival;
// UNUSED VARIABLE //                 double omega_g;
                //	       double res_e,res_i;
                longdouble resp,resc;
// UNUSED VARIABLE //                 double phi_g;
                double lambda_p,beta_p,lambda,beta;
                double n1,n2,n3;
                double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
                double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;
                double cosTheta;

                lambda_p = (double)psr[p].param[param_raj].val[0];
                beta_p   = (double)psr[p].param[param_decj].val[0];
                lambda   = psr[p].gwb_raj;
                beta     = psr[p].gwb_decj;

                // Pulsar vector
                n1 = cosl(lambda_p)*cosl(beta_p);
                n2 = sinl(lambda_p)*cosl(beta_p);
                n3 = sinl(beta_p);

                cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
                    sinl(beta)*sinl(beta_p);

                // From KJ's paper
                // Gravitational wave matrix

                // NOTE: This is for the plus terms.  For cross should use different terms
                e11p = pow(sinl(lambda),2)-pow(cosl(lambda),2)*pow(sinl(beta),2);
                e21p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
                e31p = cosl(lambda)*sinl(beta)*cosl(beta);

                e12p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
                e22p = pow(cosl(lambda),2)-pow(sinl(lambda),2)*pow(sinl(beta),2);
                e32p = sinl(lambda)*sinl(beta)*cosl(beta);

                e13p = cosl(lambda)*sinl(beta)*cosl(beta);
                e23p = sinl(lambda)*sinl(beta)*cosl(beta);
                e33p = -powl(cosl(beta),2);

                resp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
                        n2*(n1*e21p+n2*e22p+n3*e23p)+
                        n3*(n1*e31p+n2*e32p+n3*e33p));

                //		   printf("Resp = %s %g %g\n",psr[p].name,(double)resp,(double)cosTheta);
                if ((1-cosTheta)==0.0)
                    resp = 0.0;  // Check if this is sensible
                else
                    resp = longdouble(1.0)/(longdouble(2.0)*(longdouble(1.0)-cosTheta))*(resp); 

                psr[p].gwb_geom_p = resp;
                //		   printf("Resp2 = %s %g %g\n",psr[p].name,(double)resp,(double)psr[p].quad_ifunc_geom_p);
                // NOTE: These are for the cross terms. 
                lambda   = psr[p].gwb_raj;
                beta     = psr[p].gwb_decj;


                e11c = sinl(2*lambda)*sinl(beta);
                e21c = -cosl(2*lambda)*sinl(beta);
                e31c = -sinl(lambda)*cosl(beta);

                e12c = -cosl(2*lambda)*sinl(beta);
                e22c = -sinl(2*lambda)*sinl(beta);
                e32c = cosl(lambda)*cosl(beta);

                e13c = -sinl(lambda)*cosl(beta);
                e23c = cosl(lambda)*cosl(beta);
                e33c  = 0;

                resc = (n1*(n1*e11c+n2*e12c+n3*e13c)+
                        n2*(n1*e21c+n2*e22c+n3*e23c)+
                        n3*(n1*e31c+n2*e32c+n3*e33c));
                //		   printf("Resc = %s %g %g\n",psr[p].name,(double)resc,(double)cosTheta);
                if ((1-cosTheta)==0.0)
                    resc = 0.0;  // Check if this is sensible
                else
                    resc = longdouble(1.0)/(longdouble(2.0)*(longdouble(1.0)-cosTheta))*(resc); 
                psr[p].gwb_geom_c = resc;
                //printf("Resc2 = %s %g %g %g %g\n",psr[p].name,(double)resc,(double)psr[p].quad_ifunc_geom_c,(double)cosTheta,(double)(longdouble(1.0)/(longdouble(2.0)*(longdouble(1.0)-cosTheta))));

                // exp(-(t-T0)**2/2./w**2)



                dt = (psr[p].obsn[i].bbat - psr[p].gwb_epoch)/psr[p].gwb_width;

                prefac = dt*exp( (double) -dt*dt/2.);


                fprintf(stderr, "%.3e %.3e %.3e\n",(double) psr[p].gwb_width, (double) psr[p].gwb_epoch, (double) prefac);
                //exit(0);


                //		   scale = -0.5*cos2Phi*(1-cosTheta);
                phaseW += psr[p].param[param_f].val[0]*prefac*psr[p].param[param_gwb_amp].val[0]*psr[p].gwb_geom_p; 


                phaseW +=  psr[p].param[param_f].val[0]*prefac*psr[p].param[param_gwb_amp].val[1]*psr[p].gwb_geom_c; 				

            }	


            // END of Ryan's code




            phase5[i] = phase2+phase3+phase4+phaseJ+phaseW+phase2state + phaseShape;
            //	   printf("Point 1: %.5f %.5f %.5f %.5f %.5f %.5f\n",(double)phase5[i],(double)phase2,(double)phase3,(double)phase4,(double)phaseJ,(double)phaseW);
            if (psr[p].obsn[i].nFlags>0) /* Look for extra factor to add to residuals */
            {
                int k;
                longdouble extra;

                for (k=0;k<psr[p].obsn[i].nFlags;k++)
                {
                    if (strcmp(psr[p].obsn[i].flagID[k],"-radd")==0) // Add in extra time
                    {
                        extra = parse_longdouble(psr[p].obsn[i].flagVal[k]);
                        /* psr[p].obsn[i].residual+=extra; */
                        phase5[i]+=(extra*psr[p].param[param_f].val[0]);
                    }
                    if (strcmp(psr[p].obsn[i].flagID[k],"-padd")==0) // Add in extra phase
                    {
                        extra = parse_longdouble(psr[p].obsn[i].flagVal[k]);
                        /* psr[p].obsn[i].residual+=extra; */
                        phase5[i]+=(extra);
                    }
                }
            }

            /* Set the first residual to be zero */
            /* in this case fortran_mod) returns [phase5 - (int)(phase5)] 
             * ie. returns the fractional part of phase5
             */
            if (psr[p].param[param_iperharm].paramSet[0]==1)
            {  /* This code has been added for observations taken at the wrong rotational period */
                /* M. Keith Nov 2018 - I have changed this code a bit so that pulse numbers can 
                 * still be written when using IPERHARM. The old code threw away the integar
                 * part of the phase here. Now I keep it and just fix the fractional part to
                 * lie in the correct range. To be honest, it behaves strangely with PN still
                 * since the PN refers to a real pulse, not a iperharm pulse  */
                double phaseTMP = fortran_mod(phase5[i],1.0);
                if (phaseTMP >= 1.0/2.0/psr[p].param[param_iperharm].val[0])
                {
                    while (phaseTMP >= 1.0/2.0/psr[p].param[param_iperharm].val[0]){
                        phase5[i] -= 1.0/psr[p].param[param_iperharm].val[0];
                        phaseTMP -= 1.0/psr[p].param[param_iperharm].val[0];
                    }
                }
                else if (phaseTMP < -1.0/2.0/psr[p].param[param_iperharm].val[0])
                {
                    while (phaseTMP < -1.0/2.0/psr[p].param[param_iperharm].val[0]){
                        phase5[i] += 1.0/psr[p].param[param_iperharm].val[0];
                        phaseTMP += 1.0/psr[p].param[param_iperharm].val[0];

                    }
                }
                //	       if (i==0) phas1  = phase5[i];
                phase5[i] = phase5[i]; // -phas1; 	       
            }
            else
            {
                // This is where the first residual is set to equal zero
                //	       if (i==0) phas1  = fortran_mod((phase5[i]),longdouble(1.0)); 
                phase5[i] = phase5[i]; //-phas1; 
            }

        }
        // Now calculate phas1 to set the first residual equal to zero
        // 
        bool startSet = psr[p].param[param_start].paramSet[0]==1
            && psr[p].param[param_start].fitFlag[0]==1;
        bool finishSet = psr[p].param[param_finish].paramSet[0]==1
            && psr[p].param[param_finish].fitFlag[0]==1;

        bool bat_startSet = psr[p].param[param_start].paramSet[0]==1
            && psr[p].param[param_start].fitFlag[0]==2;
        bool bat_finishSet = psr[p].param[param_finish].paramSet[0]==1
            && psr[p].param[param_finish].fitFlag[0]==2;

        longdouble start = 1e10;
        longdouble finish = 0;

        // if we are fixing start/finish then use the specified values.
        if (startSet||bat_startSet) start = psr->param[param_start].val[0];
        if (finishSet||bat_finishSet) finish = psr->param[param_finish].val[0];


        for (i=0;i<psr[p].nobs;i++)
        {
            /* MJK 2021 - update to use same logic for start/finish as t2Fit */
            observation *o = psr[0].obsn+i;
            // skip deleted points
            if (o->deleted) continue;

            // if start/finish is set, skip points outside of the range
            if (startSet && o->sat < (start-START_FINISH_DELTA)) continue;
            if (finishSet && o->sat > (finish+START_FINISH_DELTA)) continue;

            if (bat_startSet && o->bat < (start-START_FINISH_DELTA)) continue;
            if (bat_finishSet && o->bat > (finish+START_FINISH_DELTA)) continue;

            phas1 = fortran_mod((phase5[i]),longdouble(1.0));
            break;

        }
        for (i=0;i<psr[p].nobs;i++)
        {
            /* phase5 - nphase is also the fractional part of phase5
             * -- this is different to using fortran_mod() as above as the nearest integer
             *    to a negative or a positive number is calculated differently. 
             */
            phase5[i] -= phas1;
            nphase = (longdouble)fortran_nlong(phase5[i]);
            psr[p].obsn[i].nphase = nphase;

            residual = phase5[i] - nphase;
            /* residual = residual in phase */
            if (psr[p].obsn[i].deleted!=1)
                gotit = 1;
            else
                gotit = 0;

            //	   printf("%s deleted = %d\n",psr[p].obsn[i].fname,psr[p].obsn[i].deleted);
            if (psr[p].param[param_track].paramSet[0]==1 && psr[p].param[param_track].val[0] > 0 && time==1)
            {
                residual+=ntrk;
                // Note that this requires that the points be in time order
                if (gotit==1)
                {
                    if (psr[p].obsn[i].bbat - ct00 < longdouble(0.0) && fabs(psr[p].obsn[i].bbat-ct00) > 1 && i > 0)
                    {
                        printf("ERROR: Points must be in time order for tracking to work\n");
                        printf("Pulsar = %s\n",psr[p].name);
                        printf("Observation %d (%s) has BBAT = %.5g\n",i,psr[p].obsn[i].fname,(double)psr[p].obsn[i].bbat);
                        printf("Observation %d (%s) has BBAT = %.5g\n",i-1,psr[p].obsn[i-1].fname,(double)psr[p].obsn[i-1].bbat);
                        ld_printf("Difference = %Lg %Lg %Lg\n",psr[p].obsn[i].bbat, ct00,psr[p].obsn[i].bbat - ct00);
                        exit(1);
                    }	       
                    if (psr[p].obsn[i].bbat - ct00 < fabs(psr[p].param[param_track].val[0]))
                    {
                        if (fabs(residual+1.0-dt00) < fabs(residual-dt00))
                        {
                            residual+=1.0;
                            ntrk+=1;
                            printf("ADDING PHASE %s\n",psr[p].obsn[i].fname);
                        } 
                        else if (fabs(residual-1.0-dt00) < fabs(residual-dt00))
                        {
                            residual-=1.0;
                            ntrk-=1;
                            printf("SUBTRACT PHASE %s\n",psr[p].obsn[i].fname);
                        } 
                    }
                }
            }

            if ((double)psr[p].param[param_track].val[0] == -1) // Do extra tracking
            {
                if (dtm1s==0) printf("Attempting tracking via the gradient method\n");

                if (dtm1s>1 && fabs(psr[p].obsn[i].bat-lastBat) > 1 
                        && fabs(lastBat-priorBat) > 1) // Have 3 points each separated by more than 1 day
                {
                    double m1,m2,m3,m4;
                    m1 = (double)(lastResidual-priorResidual)/(lastBat-priorBat);
                    m2 = (double)(residual-lastResidual)/(psr[p].obsn[i].bat-lastBat);
                    m3 = (double)(residual-1-lastResidual)/(psr[p].obsn[i].bat-lastBat);
                    m4 = (double)(residual+1-lastResidual)/(psr[p].obsn[i].bat-lastBat);
                    //		   printf("pos %d %g %g %g %g (%g) (%g) (%g)\n",i,m1,m2,m3,m4,fabs(m2-m1),fabs(m3-m1),fabs(m4-m1));

                    if (fabs(m1-m2) < fabs(m1-m3) && fabs(m1-m2) < fabs(m1-m4))
                    {
                    }
                    else if (fabs(m1-m3) < fabs(m1-m2) && fabs(m1-m3) < fabs(m1-m4))
                    {
                        //		       printf("Updating neg\n");
                        residual-=1.0;
                        ntrk-=1;
                    }
                    else if (fabs(m1-m4) < fabs(m1-m2) && fabs(m1-m4) < fabs(m1-m3))
                    {
                        //		       printf("Updating pos\n");
                        residual+=1.0;
                        ntrk+=1;
                    } 
                }
            }
            if ((double)psr[p].param[param_track].val[0] == -2) // Track on pulse number
            {
                long long pnNew;
                long long pnAct;
                long long addPhase;

                nf0  = (int)psr[p].param[param_f].val[0];
                //ntpd = ((int)psr[p].obsn[i].bbat-(int)psr[p].param[param_pepoch].val[0]);
                // M. Keith 2013. Changed to be referenced to bat0 rather than pepoch
                ntpd = ((int)psr[p].obsn[i].bbat-(int)psr[p].obsn[0].bbat);
                phaseint = nf0*ntpd*86400.0;
                pnNew = (long long)(phaseint + fortran_nlong(phase5[i]));
                logdbg("Have %g %lld %lld %lld",(double)psr[p].obsn[i].sat,pnNew-pn0,pn0,pnNew);
                if (pn0 == -1)
                {
                    pn0 = pnNew;
                    pnNew = 0;
                }
                else
                    pnNew -= pn0;
                // Compare with flag
                for (int kk=0;kk<psr[p].obsn[i].nFlags;kk++)
                {
                    if (strcmp(psr[p].obsn[i].flagID[kk],"-pnadd")==0){
                        sscanf(psr[p].obsn[i].flagVal[kk],"%lld",&pnAct);
                        pnAdd+=pnAct;
                    }
                }
                for (int kk=0;kk<psr[p].obsn[i].nFlags;kk++)
                {
                    if (strcmp(psr[p].obsn[i].flagID[kk],"-pn")==0)
                    {
                        sscanf(psr[p].obsn[i].flagVal[kk],"%lld",&pnAct);
                        pnAct+=pnAdd;
                        addPhase = pnNew-pnAct;
                        residual += addPhase;
                        ntrk = addPhase;
                        logdbg("*** Adding phase %lld ***",addPhase);
                    }
                }
            }
            if (gotit==1)
            {
                dtm1s ++;
                dt00  = residual;

                ct00 = psr[p].obsn[i].bbat;
            }
            if (gotit==1 && time==0)
                time=1;
            psr[p].obsn[i].residual = residual/psr[p].param[param_f].val[0];
	    // setting it to fit residual tn to some default value 
	    psr[p].obsn[i].residualtn = residual/psr[p].param[param_f].val[0];
	    
            priorResidual=lastResidual;
            priorBat = lastBat;
            lastResidual=residual;
            lastBat = psr[p].obsn[i].bat;

            /* Check for phase offsets */
            for (k=0;k<psr[p].nPhaseJump;k++)
            {
                //	       if (psr[p].obsn[i].sat > psr[p].phaseJump[k])
                //	       printf("Comparing with %g\n",(double)psr[p].obsn[psr[p].phaseJumpID[k]].sat);
                if (psr[p].obsn[i].sat > psr[p].obsn[psr[p].phaseJumpID[k]].sat)
                    psr[p].obsn[i].residual += (double)psr[p].phaseJumpDir[k]/psr[p].param[param_f].val[0];
            }
            nf0  = (int)psr[p].param[param_f].val[0];
            ntpd = ((int)psr[p].obsn[i].bbat-(int)psr[p].param[param_pepoch].val[0]);

            phaseint = nf0*ntpd*86400.0;
            /*	   psr[p].obsn[i].phase = nint(phase5)+phaseint; */
            /* dt in resid.f */
            /*	   psr[p].obsn[i].phase = phaseint+nphase+phase5-nphase+dphase+ddnprd; */	   
            psr[p].obsn[i].phase = phaseint+phase5[i];  
            psr[p].obsn[i].pulseN = (long long)(phaseint + fortran_nlong(phase5[i])) - (long long)ntrk;

            // correct pulse numbering for phase wraps
            for (k=0;k<psr[p].nPhaseJump;k++)
                if (psr[p].obsn[i].sat > psr[p].obsn[psr[p].phaseJumpID[k]].sat)
                    psr[p].obsn[i].pulseN -= psr[p].phaseJumpDir[k];

        

            if (psr[p].obsn[i].deleted!=1)
            {
                mean+=psr[p].obsn[i].residual;
                nmean++;

            }
        }



        if (removeMean==1)
        {
            mean/=(longdouble)nmean;
            //tnmean/=(longdouble) ntnmean;
            for (i=0;i<psr[p].nobs;i++)
            {
                psr[p].obsn[i].residual-=mean;
                //psr[p].obsn[i].residualtn-=tnmean;
            }
            // psr[p].obsn[i].residual-=0;
        }



        double newmean = 0;
        for (i=0;i<psr[p].nobs;i++){
            newmean+=((double)psr[p].obsn[i].residual-psr[p].obsn[i].TNRedSignal);
        }
        newmean = newmean/psr[p].nobs;

        // deal with the REFPHS option before we look at red noise etc.
        if(psr[p].refphs==REFPHS_TZR){
            // remove our TZR observation and store it in tzrobs for later use.
            memcpy(&(psr[p].tzrobs),&(psr[p].obsn[psr[p].nobs-1]),sizeof(observation));
            psr[p].nobs--;

            logmsg("SHIFT RESIDUALS BY %g (%g %g)\n",(double)psr[p].tzrobs.residual,(double)psr[p].tzrobs.bbat,(double)psr[p].tzrobs.sat);
            for (i=0;i<psr[p].nobs;i++){
                psr[p].obsn[i].residual -= psr[p].tzrobs.residual;
            }
        }



        if ((psr[p].TNsubtractDM ==1) || ( psr[p].TNsubtractRed ==1) || (psr[p].TNsubtractChrom ==1))
        {
            for(i=0;i<=psr[p].nobs; i++)
            {
                psr[p].obsn[i].residualtn = psr[p].obsn[i].residual;

                if (psr[p].TNsubtractRed ==1)
                {
                    psr[p].obsn[i].residualtn-=psr[p].obsn[i].TNRedSignal;
                    // Also subtract Common Mode signal if present
                    if (psr[p].dmoffsCMnum > 0) {
                        psr[p].obsn[i].residualtn -= ifunc(psr[p].dmoffsCM_mjd,psr[p].dmoffsCM,(double)psr[p].obsn[i].sat,psr[p].dmoffsCMnum);
                    }
                }

                if (psr[p].TNsubtractDM ==1)
                {
                    psr[p].obsn[i].residualtn-=psr[p].obsn[i].TNDMSignal;
                }


                if (psr[p].TNsubtractChrom ==1)
                {
                    psr[p].obsn[i].residualtn-=psr[p].obsn[i].TNChromSignal;
                }


                if (psr[p].obsn[i].deleted!=1)
                {   
                    tnmean += psr[p].obsn[i].residualtn;
                    ntnmean++;
                }
            }

            if (removeMean ==1)
            {
                tnmean /= (long double) ntnmean;
                for(i=0;i<=psr[p].nobs;i++)
                {
                    psr[p].obsn[i].residualtn -= tnmean;
                }        

            }


	  }
	  
        if(psr[p].AverageResiduals == 1){
          
	  
	  averageResiduals(psr, 1);
	    
	}
	
        if(psr[p].AverageDMResiduals == 1){
            averageDMResiduals(psr,1);

        }
        
    }

    delete[] phase5;
    logtchk("Leave formresiduals()");
}

