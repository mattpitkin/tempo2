//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards
//


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
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "T2toolkit.h"
#include "tempo2.h"
#include "toasim.h"
#include "makeRedNoise.h"

using namespace std;

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
    char parFile[MAX_PSR][MAX_FILELEN];
    char timFile[MAX_PSR][MAX_FILELEN];
    int i,nit,j,p;
    char fname[MAX_FILELEN];
    double globalParameter;
    long double result;
    long seed = TKsetSeed();

    double secperyear=365*86400.0;
    // my parameters
    int npts=1024;
    float TNChromAmp = 0;
    float TNChromIdx = 4.0;
    float TNChromGam = 3.0;
    float ref_freq=1400; // MHz
    char writeTextFiles=0;
    double lastMJD=1e99;

    int force_nreal=0;

    //
    // For the output file
    //
    toasim_header_t* header;
    toasim_header_t* read_header;
    FILE* file;
    double offsets[MAX_OBSN]; // Will change to doubles - should use malloc
    double chrom[MAX_OBSN]; // Will change to doubles - should use malloc
    // Create a set of corrections.
    toasim_corrections_t* corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));

    corr->offsets=offsets;
    corr->params=NULL; // Normally leave as NULL. Can store this along with each realisation. 
    // Same length string in every iteration - defined in r_param_length see below
    corr->a0=0; // constant
    corr->a1=0; // a1*x
    corr->a2=0; // a2*x*X

    *npsr = 0;
    nit = 1;

    printf("Graphical Interface: addDmVar\n");
    printf("Author:              M. Keith\n");
    printf("Version:             1.0\n");

    /* Obtain all parameters from the command line */
    for (i=2;i<argc;i++)
    {

        if (strcmp(argv[i],"-nreal")==0){
            nit=atoi(argv[++i]);
        }
        if (strcmp(argv[i],"-f")==0)
        {
            strcpy(parFile[*npsr],argv[++i]); 
            strcpy(timFile[*npsr],argv[++i]);
            (*npsr)++;
        }
        if (strcmp(argv[i],"-npts")==0){
            npts=atoi(argv[++i]);
        }

        if (strcmp(argv[i],"-reffreq")==0){
            ref_freq=atof(argv[++i]);
        }

        if(strcmp(argv[i],"-lastMJD")==0){
            lastMJD=atof(argv[++i]);
        }
        if (strcmp(argv[i],"-debug")==0){
            writeTextFiles=1;
        }
        if (strcmp(argv[i],"-gam")==0){
            TNChromGam=atof(argv[++i]);
            if(TNChromGam < 0){
                logmsg("Warning: gamma should normally be positive!!");
            }
        }
        if (strcmp(argv[i],"-amp")==0){
            TNChromAmp=atof(argv[++i]);
        }
        if (strcmp(argv[i],"-idx")==0){
            TNChromIdx=atof(argv[++i]);
        }

        else if (strcmp(argv[i],"-seed")==0){
            sscanf(argv[++i],"%ld",&seed);
            if (seed > 0)seed=-seed;
        }
        if (strcmp(argv[i],"-forceperiodic")==0){
            force_nreal=1;
        }


    }


    if (TNChromAmp == 0 ){
        printf("Must set -amp to set TNChromAmp\n");
        exit(1);
    }

    readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
    // Now read in all the .tim files
    readTimfile(psr,timFile,*npsr); /* Load the arrival times    */

    preProcess(psr,*npsr,argc,argv);
    formBatsAll(psr,*npsr);

    for (p=0;p<*npsr;p++)
    {
        printf("NTOA = %d\n",psr[p].nobs);
        header = toasim_init_header();
        strcpy(header->short_desc,"addChromVar");
        strcpy(header->invocation,argv[0]);
        strcpy(header->timfile_name,timFile[p]);
        strcpy(header->parfile_name,parFile[p]);
        header->seed = seed;

        header->ntoa = psr[p].nobs;
        header->nrealisations = nit;

        // First we write the header...
        sprintf(fname,"%s.addChromVar",timFile[p]);
        file = toasim_write_header(header,fname);

        double mjd_start=1000000.0;
        double mjd_end=-10000000.0;
        for (j=0;j<psr[p].nobs;j++){
            // find the start and end times
            if(psr[p].obsn[j].bat < mjd_start)mjd_start=(double)psr[p].obsn[j].bat;
            if(psr[p].obsn[j].bat > mjd_end)mjd_end=(double)psr[p].obsn[j].bat;
        }


        printf("start    = %f (mjd)\n",mjd_start);
        printf("end      = %f (mjd)\n",mjd_end  );
        printf("npts     = %d (days)\n",npts     );
        printf("ChromGam = %f\n",TNChromGam   );
        printf("ChromAmp = %f\n",TNChromAmp    );
        printf("ChromIDX = %f\n",TNChromIdx    );

        printf("ref_freq = %f (MHz)\n",ref_freq    );

        printf("seed     = %ld\n",seed);


        double pism = pow(10,2*TNChromAmp)/12.0/M_PI/M_PI;
        pism *= secperyear*secperyear;

        printf("pism(1yr)  = %g (yr^3) @reffreq\n",pism);

        printf("\n");
        printf("Generating red noise...\n");

        rednoisemodel_t* model = setupRedNoiseModel(mjd_start,mjd_end,npts,nit,pism,-TNChromGam);
        if (force_nreal){
            printf("WARNING: FORCING PERIODIC\n");
            model->nreal = nit;
        }


        populateRedNoiseModel(model,seed);

        if (writeTextFiles){
            FILE *log_spec = fopen("chromvar.spec","w");
            float* pwr_spec=getPowerSpectrum(model);
            float ps_fres=1.0/((model->end-model->start)/365.25); //yr^-1
            for (j=0;j<(model->npt/2+1);j++){
                fprintf(log_spec,"%10.10g %10.10g\n",ps_fres*j,pwr_spec[j]);
            }

            fclose(log_spec);
        }

        int itjmp=nit/50;
        if (itjmp<1)itjmp=1;
        int dots=0;
        printf("v");
        for (i=0;i<nit/itjmp;i++){
            printf("_");
        }
        printf("v\n");
        printf("[");
        fflush(stdout);
        for (i=0;i<nit;i++)
        {
            if (i%itjmp==0){
                int v = i/itjmp;
                v-=dots;
                while (v > 0){
                    printf(".");
                    fflush(stdout);
                    v--;
                    dots++;
                }
            }

            for (j=0;j<psr[p].nobs;j++){
                double t = (double)(psr[p].obsn[j].bat);
                if(t > lastMJD)t=lastMJD;
                chrom[j]=getRedNoiseValue(model,t,i);
            }
            FILE *log_ts;
            if (writeTextFiles)
                log_ts = fopen("chromvar.ts","w");
            double sum=0;
            for (j=0;j<psr[p].nobs;j++){
                sum+=chrom[j];
            }
            sum/=psr[p].nobs;
            int mm=-1;
            for (j=0;j<psr[p].nobs;j++){
                chrom[j]-=sum;
                double ofreq=psr[p].obsn[j].freqSSB/1e6; // convert to MHz
                offsets[j] = (double)(chrom[j]*pow(ref_freq/ofreq,TNChromIdx)); // offset is scaled by chromidx
                if (writeTextFiles)
                    fprintf(log_ts,"%lg %lg %lg %lg\n",(double)psr[p].obsn[j].bat,chrom[j],offsets[j],(double)ofreq);
            }
            toasim_write_corrections(corr,header,file);
            if (writeTextFiles)
                fclose(log_ts);
        }
        int v = i/itjmp;
        v-=dots;
        while (v > 0){
            printf(".");
            v--;
            dots++;
        }

        printf("]\n");

        printf("Close file\n");
        fclose(file);
    }
    return 0;
}

const char * plugVersionCheck = TEMPO2_h_VER;
