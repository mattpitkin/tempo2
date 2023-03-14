//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards
//
// Plugin to read a set of arrival time files and produce a list of Gaussian random numbers based on the TOA uncertainties 


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
    double alpha= -8.0/3.0;
    int npts=1024;
    float D_d=0; // us
    float TNDMAmp=0;
    float ref_freq=1400; // MHz
    float d=1000;
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
    double dms[MAX_OBSN]; // Will change to doubles - should use malloc
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

        if (strcmp(argv[i],"-D")==0){
           D_d=atof(argv[++i]);
        }
        if (strcmp(argv[i],"-t")==0){
            d=atof(argv[++i]);
        }
        if (strcmp(argv[i],"-nu")==0){
            ref_freq=atof(argv[++i]);
        }

        if(strcmp(argv[i],"-lastMJD")==0){
            lastMJD=atof(argv[++i]);
        }
        if (strcmp(argv[i],"-debug")==0){
            writeTextFiles=1;
        }
        if (strcmp(argv[i],"-a")==0){
            alpha=atof(argv[++i]);
            if(alpha > 0){
                logmsg("Warning: alpha should normally be negative!!");
            }

        }
        if (strcmp(argv[i],"-TNDMAmp")==0){
            TNDMAmp=atof(argv[++i]);
        }
        else if (strcmp(argv[i],"-seed")==0){
            sscanf(argv[++i],"%ld",&seed);
            if (seed > 0)seed=-seed;
        }
        if (strcmp(argv[i],"-forceperiodic")==0){
            force_nreal=1;
        }


    }


    if (TNDMAmp == 0 && D_d == 0){
        printf("Must set either -TNDMAmp or -D option\n");
        printf("Try -D 1 for a reasonable DM variation!\n");
        exit(1);
    }

    if (D_d != 0 && alpha != -8.0/3.0) {
        printf("WARNING: -D option only works for -8/3 kolmogorov turbulance");
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
        strcpy(header->short_desc,"addDmVar");
        strcpy(header->invocation,argv[0]);
        strcpy(header->timfile_name,timFile[p]);
        strcpy(header->parfile_name,"Unknown");
        header->seed = seed;

        header->ntoa = psr[p].nobs;
        header->nrealisations = nit;

        // First we write the header...
        sprintf(fname,"%s.addDmVar",timFile[p]);
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
        if (TNDMAmp == 0) {
            printf("D_d(%f)  = %f (us^2)\n",d,D_d    );
        } else {
            printf("alpha  = %f\n",alpha    );
            printf("TNDMAmp  = %f\n",TNDMAmp    );
        }
//        printf("DM_CONST  = %f\n",DM_CONST    );


        printf("ref_freq = %f (MHz)\n",ref_freq    );

        printf("seed     = %ld\n",seed);


        d*=86400.0;  // convert days to seconds.

        double pism=0;
        if (TNDMAmp == 0) {
            D_d *=1e-12; // convert us to seconds
            pism = 0.0112 * D_d * pow(d,(-5.0/3.0)) * pow(secperyear,-1.0/3.0);
        } else {
            pism = pow(10.0,(2.0*TNDMAmp)) / (pow(DM_CONST,2.0) * pow(ref_freq,4.0));
        }
        printf("pism(1yr)  = %g (yr^3) \n",pism );

        double yr2dm = secperyear * DM_CONST*pow(ref_freq,2.0);
        printf("yr2dm   = %f \n",yr2dm    );

        pism *= pow(yr2dm,2); // convert yr^3 to yr.cm^-3.pc

        printf("pism(1yr)  = %g ((cm^-3pc)^2 yr) \n",pism );

        printf("\n");
        printf("Generating red noise...\n");

        rednoisemodel_t* model = setupRedNoiseModel(mjd_start,mjd_end,npts,nit,pism,alpha);
        if (force_nreal){
            printf("WARNING: FORCING PERIODIC\n");
            model->nreal = nit;
        }


        populateRedNoiseModel(model,seed);

        if (writeTextFiles){
            FILE *log_spec = fopen("dmvar.spec","w");
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
                dms[j]=getRedNoiseValue(model,t,i);
            }
            FILE *log_ts;
            if (writeTextFiles)
                log_ts = fopen("dmvar.ts","w");
            double sum=0;
            for (j=0;j<psr[p].nobs;j++){
                sum+=dms[j];
            }
            sum/=psr[p].nobs;
            int mm=-1;
            for (j=0;j<psr[p].nobs;j++){
                dms[j]-=sum;
                double ofreq=psr[p].obsn[j].freqSSB;
                offsets[j] = (double)(dms[j]/DM_CONST/ofreq/ofreq)*1e12;
                if (writeTextFiles)
                    fprintf(log_ts,"%lg %lg %lg %lg\n",(double)psr[p].obsn[j].bat,dms[j],offsets[j],(double)ofreq);
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
