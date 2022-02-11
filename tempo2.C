#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2006,2007,2008,2009,2010,2011 George Hobbs, Russell Edwards

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
#include <time.h>
#include "tempo2.h"
#include "tempo2Util.h"
#include "tempo2pred.h"
#include "tempo2pred_int.h"
#include "T2accel.h"
#include "t2fit.h"
#include <dlfcn.h>

// #include "T2toolkit.h"

void ephemeris_routines(pulsar *psr,int npsr);
void clock_corrections(pulsar *psr,int npsr);
void extra_delays(pulsar *psr,int npsr);

int main(int argc, char *argv[])
{
    int iteration; 
    int listparms;
    int outRes=0;
    int writeModel=0;
    int writeTimFile=0;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
    char outputSO[100];
    char str[MAX_STRLEN];
    char newparname[MAX_FILELEN];
    longdouble coeff[MAX_COEFF]; /* For polynomial coefficients in polyco */
    int npsr;      /* The number of pulsars */
    int noWarnings=1;
    double globalParameter=0.0;
    int  p;
    int nGlobal,i,flagPolyco=0,it,k;
    char polyco_args[128];
    char polyco_file[128]; /* buffer for optional polyco filenames */
    int newpar=0;
    int onlypre=0;
    FILE *alias;
    char **commandLine;
    clock_t startClock,endClock;
    const char *CVS_verNum = "$Id$";

    polyco_file[0] = '\0';

    timer_clk=clock();  
    if (argc > 1 && ((strcmp(argv[1],"-version") == 0 ) || (strcmp(argv[1],"-v") == 0))){
        printf("%s\n",VERSION);
        exit(1);
    }
    printf("This program comes with ABSOLUTELY NO WARRANTY.\n");
    printf("This is free software, and you are welcome to redistribute it\n");
    printf("under conditions of GPL license.\n\n");

    startClock = clock();

    pulsar *psr;

    commandLine = (char **)malloc(1000*sizeof(char *));

    for (i=0;i<1000;i++)
        commandLine[i] = (char *)malloc(sizeof(char)*1000);

    /* Parse input line for machine type */
    for (i=0;i<argc;i++)    
    {
        if (strcasecmp(argv[i],"-npsr")==0)
            sscanf(argv[i+1],"%d",&MAX_PSR);
        else if (strcasecmp(argv[i],"-displayVersion")==0)
            displayCVSversion = 1;
        else if (strcasecmp(argv[i],"-nobs")==0)
            sscanf(argv[i+1],"%d",&MAX_OBSN);
        else if (strcasecmp(argv[i],"-debug")==0){
            debugFlag=1;
            tcheck=1;
            writeResiduals=0xff;
        }
        else if (strcasecmp(argv[i],"-tcheck")==0)
            tcheck=1;
        else if (strcasecmp(argv[i],"-noaccel")==0)
            useT2accel=0;
        else if (strcasecmp(argv[i],"-qrfit")==0)
            useT2accel=2;
        else if (strcasecmp(argv[i],"-svdfit")==0){
            if(useT2accel) useT2accel=1;
        }
        else if (strcasecmp(argv[i],"-oldfit")==0)
            NEWFIT=0;
        else if (strcasecmp(argv[i],"-newfit")==0)
            NEWFIT=1;
        else if (strcasecmp(argv[i],"-writeres")==0)
            writeResiduals=0xff;
        else if (strcasecmp(argv[i],"-writetim")==0)
            writeTimFile=1;
        else if (strcasecmp(argv[i],"-veryfast")==0)
            veryFast=1;
        else if (strcasecmp(argv[i],"-allglitch")==0)
            forceAlwaysFitForGlitches=1;
        strcpy(commandLine[i],argv[i]);
    }
    if (displayCVSversion == 1) CVSdisplayVersion("tempo2.C","main()",CVS_verNum);

    if ((psr = (pulsar *)malloc(sizeof(pulsar)*MAX_PSR))==NULL)
    {
        printf("Not enough memory to allocate room for %d pulsars\n",MAX_PSR);
        printf("Please decrease the value of MAX_PSR_VAL in tempo2.h\n"); 
        exit(1); 
    }
    logdbg("Have allocated memory for pulsar");
    bool jboFormat = false;

    for (i=1;i<argc;i++)
    {
        if (strcmp(commandLine[i],"-machine")==0)
            strncpy(tempo2MachineType,commandLine[++i],99);
        else if (strcasecmp(commandLine[i],"-noWarnings")==0)
            noWarnings=2;
        else if (strcasecmp(commandLine[i],"-allInfo")==0)
            noWarnings=0;
        else if (strcmp(commandLine[i],"-jbo")==0){
            logmsg("Reading JBO format!!");
            jboFormat=true;
        }
        else if (strcmp(commandLine[i],"-test")==0) /* Use TEMPO2_TEST environment variable */
            strcpy(TEMPO2_ENVIRON,"TEMPO2_TEST");
        else
        {
            char oldCommandLine[1000][1000];
            int  oldArgc = argc;
            int  oldI = i;

            for (k=i;k<argc;k++)
                strcpy(oldCommandLine[k-i],commandLine[k]);

            snprintf(str,MAX_FILELEN,"%s/alias.dat",getenv(TEMPO2_ENVIRON));	  
            if ((alias = fopen(str,"r")))
            {
                while (!feof(alias))
                {
                    char *ret;
                    fgets(str,MAX_FILELEN,alias);
                    ret = strtok(str," ");
                    if (strcmp(ret,commandLine[i])==0)
                    {
                        printf("Using alias: %s\n",ret);
                        /* Have an alias: now split up the remainder of the line */
                        while ((ret = strtok(NULL," "))!=NULL)
                        {
                            strcpy(commandLine[i],ret);
                            if (commandLine[i][strlen(commandLine[i])-1]=='\n')
                                commandLine[i][strlen(commandLine[i])-1] = '\0';
                            i++;
                            argc++;
                        }
                        argc--;
                        for (k=1;k<oldArgc-oldI;k++)
                            strcpy(commandLine[i+k-1],oldCommandLine[k]);
                        i=oldI-1;

                    }
                }
                fclose(alias);
            }
        }
    }
    int ii;
    for(ii=1;ii<argc;ii++){
        if (strcasecmp(commandLine[ii],"-reminder")==0){
            // Writing command line to log file
            char commandfile[200] = "T2command.input";
            FILE *fout;
            time_t rawtime;
            struct tm * timeinfo;
            time ( &rawtime );
            timeinfo = localtime ( &rawtime );
            char timeval[200];
            strcpy(timeval,asctime (timeinfo));
            strcpy(&timeval[(int)strlen(timeval)-1],"");
            fout = fopen(commandfile,"a");
            fprintf(fout,"[%s]>> ",timeval);
            for(i=0;i<argc;i++){
                fprintf(fout," %s ",commandLine[i]);
            }
            fprintf(fout,"\n");
            fclose(fout);
        }
    }

    if (getenv(TEMPO2_ENVIRON)==NULL)
    {
        printf("Environment variable >%s< not set\n",TEMPO2_ENVIRON);
        exit(1);
    }

    /* get path to look for plugins */
    setPlugPath();


    /* If running from the command line ... */
    logdbg("Running initialise");
    initialise(psr,noWarnings); /* Initialise all */
    if(jboFormat) psr[0].jboFormat=1;
    logdbg("Completed running initialise %d",psr[0].nits);
    /* Obtain login architecture */
    if (strlen(tempo2MachineType)==0)
    {
#ifdef  TEMPO2_ARCH 
        strncpy(tempo2MachineType, TEMPO2_ARCH,99);
#else
        if (getenv("LOGIN_ARCH")==NULL)
        {
            printf("Unable to determine machine type: You must do one of the following:\n"
                    "Re-compile tempo2 with the standard export distrubution, or\n"
                    "Set the LOGIN_ARCH environment variable, or\n"
                    "Use -machine on the command line\n");
            exit(1);
        }
        strncpy(tempo2MachineType, getenv("LOGIN_ARCH"),99);
#endif
    }

    if (noWarnings<1)
    {
#ifdef LONGDOUBLE_IS_IEEE754
        printf("Warning: longdouble is an IEEE754 80-bit float.\n");
        printf(" --- the size of a double is %u bytes\n",(unsigned)sizeof(double));
#endif
    }
    strcpy(outputSO,"");
    if (argc==1) /* No command line arguments */
    {
        printf("%s\n",PACKAGE_STRING);
        printf("  Usage:     %s -f XXX.par XXX.tim\n",argv[0]);
        printf("Plugin search paths:\n");
        for (i=0; i < tempo2_plug_path_len; i++){
            printf(" -- %s/*.t2\n",tempo2_plug_path[i]);
        }
        printf("\n");
#ifdef HAVE_LAPACK
        printf("* Using LAPACK acceleration for Cholesky decomposition\n");
#else

        printf("* Using slow linear algebra code (compile with LAPACK for speed improvements)\n");
#endif
#ifdef HAVE_BLAS
        printf("* Using BLAS acceleration for matrix mulitplication\n");
#else
        printf("* Using slow matrix code (compile with BLAS for speed improvements)\n");
#endif
        printf("\nFor more help, use %s -h\n",argv[0]);
        exit(1);
    }
    npsr = 0;   /* Initialise the number of pulsars */
    nGlobal=0;



    // set the extra clock path if
    if (getenv("TEMPO2_CLOCK_DIR")!=NULL){
        strncpy(tempo2_clock_path,getenv("TEMPO2_CLOCK_DIR"), MAX_STRLEN);
    }
    /* Obtain command line arguments */
    logdbg("Running getInputs %d",psr[0].nits);
    getInputs(psr,argc, commandLine, timFile,parFile,&listparms,&npsr,&nGlobal,&outRes,&writeModel,
            outputSO,&flagPolyco,polyco_args,polyco_file,&newpar,&onlypre,dcmFile,covarFuncFile,newparname);
    logdbg("Completed getInputs");

    for (i=1;i<argc;i++)
    {
        if (strcmp(commandLine[i],"-gr")==0 || strcmp(commandLine[i],"-gr2")==0) 
            /* Running from a graphical interface? */
        {      
            char *(*entry)(int,char **,pulsar *,int *);
            void * module;

            if (strcmp(commandLine[i],"-gr2")==0){
                snprintf(str,MAX_FILELEN,"./%s_%s_plug.t2",commandLine[i+1],tempo2MachineType);
                printf("Looking for %s\n",str);
                module = dlopen(str, RTLD_NOW|RTLD_GLOBAL);
            } else{
                for (int iplug=0; iplug < tempo2_plug_path_len; iplug++){
                    snprintf(str,MAX_STRLEN,"%s/%s_%s_plug.t2",tempo2_plug_path[iplug],
                            commandLine[i+1],tempo2MachineType);
                    printf("Looking for %s\n",str);
                    module = dlopen(str, RTLD_NOW|RTLD_GLOBAL); 
                    if(module==NULL){	  
                        printf("dlerror() = %s\n",dlerror());
                    } else break;
                }
            }
            if(!module)  {
                fprintf(stderr, "[error]: dlopen() failed while resolving symbols.\n" );
                //	    fprintf(stderr, "dlerror() = %s\n",dlerror());
                return -1;
            }

            /*
             * Check that the plugin is compiled against the same version of tempo2.h
             */
            char ** pv  = (char**)dlsym(module, "plugVersionCheck");
            if(pv!=NULL){
                // there is a version check for this plugin
                if(strcmp(TEMPO2_h_VER,*pv)){
                    fprintf(stderr, "[error]: Plugin version mismatch\n");
                    fprintf(stderr, " '%s' != '%s'\n",TEMPO2_h_VER,*pv);
                    fprintf(stderr, " Please recompile plugin against same tempo2 version!\n");
                    dlclose(module);
                    return -1;
                }
            }

            entry = (char*(*)(int,char **,pulsar *,int *))dlsym(module, "graphicalInterface");
            if( entry == NULL ) {
                dlclose(module);
                fprintf(stderr, "[error]: dlerror() failed while  retrieving address.\n" ); 
                fprintf(stderr, "dlerror() = %s\n",dlerror());
                return -1;
            }
            logdbg("--ENTER GRAPHICAL PLUGIN--");
            logtchk("Start graphical plugin");
            entry(argc,commandLine,psr,&npsr);
            logtchk("End graphical plugin");
            return 0;
        }
    }
    logdbg("Reading par file");
    readParfile(psr,parFile,timFile,npsr); /* Read .par file to define the pulsar's initial parameters */  
    logdbg("Finished reading par file %d",psr[0].nits);
    if (flagPolyco==0)
    {
        logdbg("Running readTimfile");
        readTimfile(psr,timFile,npsr); /* Read .tim file to define the site-arrival-times */
        logdbg("Completed readTimfile %d",psr[0].param[param_ecc].paramSet[1]);
    }

    logdbg("Running preProcess %d",psr[0].nits);
    preProcess(psr,npsr,argc,commandLine);
    logdbg("Completed preProcess %d",psr[0].nits);
    if(writeTimFile) {
        logmsg("write out.tim");
        writeTim("out.tim",psr,"tempo2");
    }


    if (flagPolyco> 0)  /* Running tempo2 in polyco mode? */
    {
        if (flagPolyco == 1)
        {
            longdouble mjd1, mjd2,maxha,freq;
            int ncoeff,nspan; 
            char sitename[128];
            char str1[128];
            char str2[128];
            char str3[128];
            char str4[128];
            if (sscanf(polyco_args, "%s %s %d %d %s %s %s", str1, str2, &nspan,
                        &ncoeff, str3,sitename,str4)!=7)
            {
                fprintf(stderr, "Error parsing -polyco arguments! See tempo2 -h.\n");
                printf("Have: %s\n",polyco_args);
                exit(1);
            }
            mjd1  = parse_longdouble(str1);
            mjd2  = parse_longdouble(str2);
            maxha = parse_longdouble(str3);
            freq  = parse_longdouble(str4);
            if (psr[0].tempo1 == 0)
                printf("WARNING: Should probably use -tempo1 option\n");

            if (psr[0].param[param_tzrmjd].paramSet[0]==0)
            {
                printf("WARNING: tzrmjd not set.  Setting to %g\n",(double)psr[0].param[param_pepoch].val[0]);
                psr[0].param[param_tzrmjd].paramSet[0]=1;
                psr[0].param[param_tzrmjd].val[0] = psr[0].param[param_pepoch].val[0];
            }
            if (psr[0].param[param_tzrfrq].paramSet[0]==0)
            {
                printf("WARNING: tzrfrq not set.  Setting to %g\n",(double)1400.0);
                psr[0].param[param_tzrfrq].paramSet[0]=1;
                psr[0].param[param_tzrfrq].val[0] = 1400.0;
            }
            if (strcmp(psr[0].tzrsite,"NULL")==0)
            {
                printf("WARNING: tzrsite not set.  Setting to %s\n",sitename);
                strcpy(psr[0].tzrsite,sitename);
            }


            polyco(psr,npsr,mjd1, mjd2,nspan,ncoeff,maxha,sitename,freq,coeff,1,polyco_file);
        }
        else if (flagPolyco==2)
        {
            longdouble seg_length;
            int ntimecoeff, nfreqcoeff;
            char sitename[64];
            char str1[128];
            char str2[128];
            char str3[128];
            char str4[128];
            char str5[128];
            longdouble mjd_start, mjd_end;
            longdouble freq_start, freq_end;
            printf("Calculating predictor >%s<\n",psr[0].tzrsite);
            if (sscanf(polyco_args, "%s %s %s %s %s %d %d %s", sitename,
                        str1, str2, str3, str4,
                        &ntimecoeff, &nfreqcoeff, str5)!=8)
            {
                fprintf(stderr, "Error parsing -pred arguments! See tempo2 -h.\n");
                printf("Have: %s\n",polyco_args);
                exit(1);
            }
            mjd_start  = parse_longdouble(str1);
            mjd_end    = parse_longdouble(str2);
            freq_start = parse_longdouble(str3);
            freq_end   = parse_longdouble(str4);
            seg_length = parse_longdouble(str5);
            /* Actually want seg_length in days */
            seg_length /= SECDAYl;

            if (ntimecoeff%2)
            {
                ntimecoeff++;
                printf("Adjusting number of coefficients (time axis) to %d (must be even)\n", ntimecoeff);
            }
            if (nfreqcoeff%2)
            {
                nfreqcoeff++;
                printf("Adjusting number of coefficients (frequency axis) to %d (must be even)\n", nfreqcoeff);
            }
            if (psr[0].param[param_tzrmjd].paramSet[0]==0)
            {
                printf("WARNING: tzrmjd not set.  Setting to %g\n",(double)psr[0].param[param_pepoch].val[0]);
                psr[0].param[param_tzrmjd].paramSet[0]=1;
                psr[0].param[param_tzrmjd].val[0] = psr[0].param[param_pepoch].val[0];
            }
            if (psr[0].param[param_tzrfrq].paramSet[0]==0)
            {
                printf("WARNING: tzrfrq not set.  Setting to %g\n",(double)1400.0);
                psr[0].param[param_tzrfrq].paramSet[0]=1;
                psr[0].param[param_tzrfrq].val[0] = 1400.0;
            }
            if (strcmp(psr[0].tzrsite,"NULL")==0)
            {
                printf("WARNING: tzrsite not set.  Setting to %s\n",sitename);
                strcpy(psr[0].tzrsite,sitename);
            }

            ChebyModelSet cms;
            ChebyModelSet_Construct(&cms, psr, sitename, mjd_start, mjd_end,
                    seg_length, seg_length*0.1, 
                    freq_start, freq_end, ntimecoeff, nfreqcoeff);
            FILE *f = fopen("t2pred.dat", "w");
            if (!f)
            {
                fprintf(stderr, "Could not open t2pred.dat for writing!\n");
                exit(1);
            }
            ChebyModelSet_Write(&cms, f);
            fclose(f);
            long double rms, mav; // long double only for predictors
            ChebyModelSet_Test(&cms, psr, ntimecoeff*5*cms.nsegments, 
                    nfreqcoeff*5*cms.nsegments, &rms, &mav);
            printf("Predictive model constructed and written to t2pred.dat.\n");
            printf("RMS error = %.3Lg s MAV= %.3Lg s\n", 
                    rms/static_cast<long double>(psr[0].param[param_f].val[0]), static_cast<long double>(mav/psr[0].param[param_f].val[0]));
            ChebyModelSet_Destroy(&cms);
        }

        logtchk("Return from main()");
        return 0;
    }
    if (debugFlag==1)
    {
        logdbg("Number of iterations = %d",psr[0].nits);
        logdbg("Maximum number of parameters = %d",MAX_PARAMS);
        logdbg("Number of pulsars = %d",npsr);
    }
    for (it=0;it<psr[0].nits;it++) /* Why pulsar 0 should select the iterations? */
    {
        if (it>0) /* Copy post-fit values to pre-fit values */
        {
            for (i=0;i<MAX_PARAMS;i++)
            {
                for (p=0;p<npsr;p++)
                {
                    for (k=0;k<psr[p].param[i].aSize;k++)
                    {
                        psr[p].param[i].prefit[k] = psr[p].param[i].val[k];
                        psr[p].param[i].prefitErr[k] = psr[p].param[i].err[k];
                    }
                }
            }
        }
        //      long seed = TKsetSeed();
        for (iteration=0;iteration<2;iteration++) /* Do pre- and post- fit analysis */
        {
            logdbg("iteration %d",iteration);
            logdbg("calling formBatsAll");
            logtchk("call formBatsAll()");
            //	  printf("Calling formBats\n");
            formBatsAll(psr,npsr);                /* Form Barycentric arrival times */
            logdbg("calling formResiduals");
            logtchk("call formResiduals()");
            formResiduals(psr,npsr,1);       /* Form residuals */
            // 
            //	  printf("WARNING: SIMULATING GAUSSIAN NOISE\n");
            //	  {
            //	    FILE *fin;
            //	    fin = fopen("tempRes.dat","r");
            //	    for (i=0;i<psr[0].nobs;i++)
            //	      {
            //		fscanf(fin,"%Lf %Lf",&psr[0].obsn[i].sat,&psr[0].obsn[i].residual);
            //		printf("Res0 = %g\n",(double)psr[0].obsn[0].residual);
            //		psr[0].obsn[i].bat = psr[0].obsn[i].bbat = psr[0].obsn[i].sat;
            //		psr[0].obsn[i].bbat+=psr[0].param[param_pepoch].val[0];
            //		psr[0].obsn[i].toaErr = 1.0e6;
            //		      psr[0].obsn[i].sat  = i;
            //	  	      psr[0].obsn[i].bat  = i; // psr[0].obsn[0].sat+i*14;
            //	  	      psr[0].obsn[i].bbat = psr[0].obsn[i].bat;
            //	  	      psr[0].obsn[i].residual = TKgaussDev(&seed);
            //		      fprintf(fout,"%.15f %.15f\n",(double)psr[0].obsn[i].sat,(double)psr[0].obsn[i].residual);
            //	      }
            //	    fclose(fin);
            //	  }
            if (listparms==1 && iteration==0)displayParameters(13,timFile,parFile,psr,npsr); /* List out all the parameters */  
            if (iteration==0)          /* Only fit to pre-fit residuals */
            {
                logdbg("calling doFit");

                logtchk("calling t2Fit");
                t2Fit(psr,npsr,covarFuncFile);
                logmsg("Complete fit");
                /* doFitGlobal(psr,npsr,&globalParameter,nGlobal,writeModel);*/ /* Fit to the residuals to obtain updated parameters  */
                logdbg("completed doFit");
            }
            if (iteration==1 || onlypre==1)
            {
                if (strlen(outputSO)==0)
                    textOutput(psr,npsr,globalParameter,nGlobal,outRes,newpar,newparname); /* Output results to the screen */
                else  /* Use a plug in for the output */
                {
                    char *(*entry)(int,char **,pulsar *,int);
                    void * module;
                    for (int iplug=0; iplug < tempo2_plug_path_len; iplug++){
                        snprintf(str,MAX_STRLEN,"%s/%s_%s_plug.t2",tempo2_plug_path[iplug],
                                outputSO,tempo2MachineType);
                        printf("Looking for %s\n",str);
                        module = dlopen(str, RTLD_NOW|RTLD_GLOBAL); 
                        if(module==NULL){	  
                            printf("dlerror() = %s\n",dlerror());
                        } else break;
                    }
                    if(!module)  {
                        fprintf(stderr, "[error]: dlopen() failed while resolving symbols.\n" );
                        return -1;
                    }
                    /*
                     * Check that the plugin is compiled against the same version of tempo2.h
                     */
                    char ** pv  = (char**)dlsym(module, "plugVersionCheck");
                    if(pv!=NULL){
                        // there is a version check for this plugin
                        if(strcmp(TEMPO2_h_VER,*pv)){
                            fprintf(stderr, "[error]: Plugin version mismatch\n");
                            fprintf(stderr, " '%s' != '%s'\n",TEMPO2_h_VER,*pv);
                            fprintf(stderr, " Please recompile plugin against same tempo2 version!\n");
                            dlclose(module);
                            return -1;
                        }
                    }


                    entry = (char*(*)(int,char **,pulsar *,int))dlsym(module, "tempoOutput");
                    if( entry == NULL ) {
                        dlclose(module);
                        fprintf(stderr, "[error]: dlerror() failed while  retrieving address.\n" ); 
                        return -1;
                    }
                    entry(argc,argv,psr,npsr);
                }
            }
            psr[0].noWarnings=2;
            if (onlypre==1) iteration=2;
            /*	  textOutput(psr,npsr,globalParameter,nGlobal,outRes,newpar,"new.par");*/ /* Output results to the screen */
            /*	  printf("Next iteration\n");*/
        }
    }
    endClock = clock();
    printf("Finishing off: time taken = %.2f (s)\n",(endClock-startClock)/(float)CLOCKS_PER_SEC);
    logtchk("Exit from main()");
    exit(EXIT_SUCCESS);
} 






// redwards function to force linkage with library functions used by
// plugins
    void
thwart_annoying_dynamic_library_stuff(int never_call_me, float or_sink)
{
    ChebyModel *cm=0;
    T2Predictor *t2p=0;
    ChebyModel_Init(cm, 0, 0);
    T2Predictor_GetPhase(t2p, 0, 0);
}
