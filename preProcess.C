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
#include "tempo2.h"
#include <math.h>
#include <string.h>
#include <dlfcn.h>

void readWhiteNoiseModelFile(pulsar *psr,int p);

void preProcess(pulsar *psr,int npsr,int argc,char **argv)
{
    int p,i,k,fitN=0,setN=0,j;
    char fitStr[10][100];
    char setStr[10][100];
    float dmvals[10000];
    float startdmmjd = 0;
    int ndm;
    longdouble setVal[10];
    FILE *fdmin;
    char newEpoch[100]="NONE";
    char selectFname[1000]="";
    char globalFname[1000]="";
    char selectPlugName[1000]="";
    char line[MAX_STRLEN];
    char hashcheck;
    char name[100];
    char dmfile[100]="";
    int setName=0;
    double last=-1;
    int tempo1=0;
    int nojump=0;
    int nofit=0;
    int modify=0;
    //trim data to match dm correction, but don't correct for dm
    int trimonly = 0;
    char modifyFname[100];
    double simulate=0;
    const char *CVS_verNum = "$Id$";

    if (displayCVSversion == 1) CVSdisplayVersion("preProcess.C","preProcess()",CVS_verNum);

    // logmsg("PreProcess");
    logdbg("In preProcess");

    //MAX_PSR   = MAX_PSR_VAL;    /* Maximum number of pulsars to fit simultaneously  */
    //MAX_OBSN  = MAX_OBSN_VAL;
    ECLIPTIC_OBLIQUITY = ECLIPTIC_OBLIQUITY_VAL;
    //  debugFlag = 0;

    for (i=0;i<argc;i++)
    {
        if (strcmp(argv[i],"-epoch")==0)
            strcpy(newEpoch,argv[++i]);
        else if (strcmp(argv[i],"-last")==0)
            sscanf(argv[++i],"%lf",&last);
        else if (strcmp(argv[i],"-setdm")==0) 
            sscanf(argv[++i],"%s",dmfile); // Should deal with multiple pulsars
        else if (strcmp(argv[i],"-trimonly")==0)
            trimonly = 1;
        else if (strcmp(argv[i],"-tempo1")==0)
            tempo1=1;
        else if (strcmp(argv[i],"-nojump")==0)
            nojump=1;
        else if (strcmp(argv[i],"-select")==0)
            sscanf(argv[++i],"%s",selectFname);
        else if (strcmp(argv[i],"-splug")==0){
            sscanf(argv[++i],"%s",selectPlugName);
            logdbg("Splug = %s\n",selectPlugName);
        }
        else if (strcmp(argv[i],"-global")==0){
            forceGlobalFit=1;
            sscanf(argv[++i],"%s",globalFname);
        }
        else if (strcmp(argv[i],"-modify")==0)
        {
            modify=1;
            sscanf(argv[++i],"%s",modifyFname);
        }
        else if (strcmp(argv[i],"-name")==0)
        {
            setName=1;
            strcpy(name,argv[i+1]);
        }
        else if (strcmp(argv[i],"-fit")==0)
            strcpy(fitStr[fitN++],argv[i+1]);
        else if (strcmp(argv[i],"-set")==0)
        {
            strcpy(setStr[setN],argv[i+1]);
            //sscanf(argv[i+2],"%Lf",&setVal[setN]);
            setVal[setN] = parse_longdouble(argv[i+2]);
            setN++;
        }
        else if (strcmp(argv[i],"-simulate")==0)
        {
            sscanf(argv[i+1],"%lf",&simulate);
        }
        else if (strcmp(argv[i],"-nofit")==0)
            nofit = 1;
        else if (strcmp(argv[i],"-clock")==0)
        {
            for (p=0;p<npsr;p++)
                strcpy(psr[p].clock,argv[i+1]);
        }
    }
    logdbg("Parsed command line");
    for (p=0;p<npsr;p++)
    {
        for (i=0;i<MAX_PARAMS;i++){
            if(psr[p].param[i].nLinkTo>0){
                psr[p].param[i].val[0] = getParameterValue(&psr[p],i,0);
                psr[p].param[i].prefit[0] = getParameterValue(&psr[p],i,0);
            }
        }
        if (setName==1)
            strcpy(psr[p].name,name);
        if (nojump==1)
            psr[p].nJumps=0;
        // Check using white noise model file
        if (strcmp(psr[p].whiteNoiseModelFile,"NULL")!=0)
        {
            readWhiteNoiseModelFile(psr,p);
        }

        // check if dm_series mode is undefined
        if (psr[p].dm_series_type == series_undefined ){
            psr[p].dm_series_type = series_taylor_pn;
            for (int k=2; k < psr[p].param[param_dm].aSize ; ++k){
                if (psr[p].param[param_dm].paramSet[k]) {
                    logwarn("PSR %s uses DM2+ but does not define DM_SERIES. Assume Taylor. This has behaviour has changed since June 2020!\nSee https://bitbucket.org/psrsoft/tempo2/issues/27/tempo2-dm-polynomial-is-not-a-taylor\n", psr[p].name);
                    break;
                }
            }
        }

        if (nofit==1)
        {
            for (i=0;i<MAX_PARAMS;i++)
            {
                if (i!=param_start && i!=param_finish){
                    for (k=0;k<psr[p].param[i].aSize;k++)
                        psr[p].param[i].fitFlag[k] = 0;
                }
            }
            // Turn off fitting for jumps
            for (i=0;i<=psr[p].nJumps;i++)
                psr[p].fitJump[i]=0;
        }
        /* Select command line fitting */
        for (i=0;i<fitN;i++)
        {
            for (j=0;j<MAX_PARAMS;j++)
            {
                for (k=0;k<psr[p].param[j].aSize;k++)
                {
                    if (strcasecmp(fitStr[i],psr[p].param[j].shortlabel[k])==0)
                    {
                        if (psr[p].param[j].paramSet[k]!=1)
                        {
                            psr[p].param[j].paramSet[k]=1;
                            psr[p].param[j].val[k]=0.0;
                            psr[p].param[j].prefit[k]=0.0;		      
                        }
                        psr[p].param[j].fitFlag[k]=1;
                    }
                }
            }
        }
        /* Simulate global parameter */
        if (simulate!=0.0)
        {
            for (i=0;i<psr[p].nobs;i++)
            {
                psr[p].obsn[i].sat += simulate/SECDAY*sin(0.003*psr[p].obsn[i].sat);
            }
        }

        /* Set jump values if already set */
        /* MOVED FOLLOWING INTO READPARFILE.C */
        /*      for (k=1;k<=psr[p].nJumps;k++)
                {	
                v5 = -1;
                nread = sscanf(psr[p].jumpStr[k],"%s %s %s %s %s",str1,str2,str3,str4,str5);

                if (strcasecmp(str1,"MJD")==0 || strcasecmp(str1,"FREQ")==0)
                {
                if (nread>3)
                {
                sscanf(str4,"%lf",&(psr[p].jumpVal[k]));
                if (sscanf(str5,"%d",&v5)==1)
                {
                if (v5!=1) psr[p].fitJump[k]=0;
                }
                else
                psr[p].fitJump[k]=0;
                }
                }
                else if (strcasecmp(str1,"NAME")==0 || strcasecmp(str1,"TEL")==0 || str1[0]=='-')
                {
                if (nread>2)
                {
                sscanf(str3,"%lf",&(psr[p].jumpVal[k]));
                if (sscanf(str4,"%d",&v5)==1)
                {
                if (v5!=1) psr[p].fitJump[k]=0;
                }
                else
                psr[p].fitJump[k]=0;
                }
                }
                } */
        /* Select command line parameter setting */
        for (i=0;i<setN;i++)
        {
            if (strcasecmp(setStr[i],"NITS")==0)
                psr[p].nits = (int)setVal[i];

            for (j=0;j<MAX_PARAMS;j++)
            {
                for (k=0;k<psr[p].param[j].aSize;k++)
                {
                    if (strcasecmp(setStr[i],psr[p].param[j].shortlabel[k])==0)
                    {
                        psr[p].param[j].val[k]=setVal[i];
                        psr[p].param[j].prefit[k]=setVal[i];
                        psr[p].param[j].paramSet[k]=1;
                    }
                }
            }
        }
        preProcessSimple1 (psr + p, tempo1, last);
        /* Update period epoch if necessary */
        if (strcmp(newEpoch,"NONE")!=0) {
            updateEpoch_str(psr,p,newEpoch);
        }

        if (psr[p].param[param_pepoch].paramSet[0]==1 && psr[p].param[param_pepoch].fitFlag[0]==1)
        {printf("Warning: Cannot fit for pepoch\n"); psr[p].param[param_pepoch].fitFlag[0]=0;}
        if (psr[p].param[param_posepoch].paramSet[0]==1 && psr[p].param[param_posepoch].fitFlag[0]==1)
        {printf("Warning: Cannot fit for posepoch\n"); psr[p].param[param_posepoch].fitFlag[0]=0;}
        if (psr[p].param[param_dmepoch].paramSet[0]==1 && psr[p].param[param_dmepoch].fitFlag[0]==1)
        {printf("Warning: Cannot fit for dmepoch\n"); psr[p].param[param_dmepoch].fitFlag[0]=0;}
        if (psr[p].param[param_track].paramSet[0]==1 && psr[p].param[param_track].fitFlag[0]==1)
        {printf("Warning: Cannot fit for track\n"); psr[p].param[param_track].fitFlag[0]=0;}

        if (strlen(dmfile)>0)
        {
            float tt;
            fdmin = fopen(dmfile,"r");
            ndm=0;
            while (!feof(fdmin))
            {
                //check wheter a line is a comment (= starts with a hash)
                fgets(line,MAX_STRLEN,fdmin);
                sscanf(line,"%c",&hashcheck);
                if (hashcheck == '#') {
                    //do nothing, perhaps give a debug message
                    logdbg("preProces():skipping line in dmfile\n");
                }
                else if (sscanf(line,"%f %f",&tt,&dmvals[ndm])==2)
                {
                    if (ndm==0)
                        startdmmjd = tt;
                    ndm++;
                }
            }
            fclose(fdmin);
        }
        //Check efacs and equads
        if (psr[p].nT2efac > 0 || psr[p].nT2equad > 0 || psr[p].T2globalEfac!=1.0)
        {
            double err;
            printf("Updating TOA errors using T2EFAC, T2EQUAD and T2GLOBALEFAC\n");
            for (i=0;i<psr[p].nobs;i++)
            {
                err = psr[p].obsn[i].toaErr;
                for (j=0;j<psr[p].obsn[i].nFlags;j++)
                {
                    //	   Check equad
                    for (k=0;k<psr[p].nT2equad;k++)
                    {
                        if (strcmp(psr[p].obsn[i].flagID[j],psr[p].T2equadFlagID[k])==0)
                        {
                            if (strcmp(psr[p].obsn[i].flagVal[j],psr[p].T2equadFlagVal[k])==0)
                                err = (sqrt(pow(err,2)+pow(psr[p].T2equadVal[k],2)));
                        }
                    }
                    //	   Check efac
                    for (k=0;k<psr[p].nT2efac;k++)
                    {
                        if (strcmp(psr[p].obsn[i].flagID[j],psr[p].T2efacFlagID[k])==0)
                        {
                            if (strcmp(psr[p].obsn[i].flagVal[j],psr[p].T2efacFlagVal[k])==0)
                                err *= psr[p].T2efacVal[k];
                        }
                    }
                }
                err *= psr[p].T2globalEfac;
                psr[p].obsn[i].toaErr = err;
            }
        }

		

        // Check TNEF and TNEQ
        if (psr[p].nTNEF > 0 || psr[p].nTNEQ > 0 || psr[p].nTNSQ > 0  || psr[p].TNGlobalEF > 0 || psr[p].TNGlobalEQ != 0)
        {
            double err;
            printf("Updating TOA errors using TN parameters.\n");
            for (i=0;i<psr[p].nobs;i++)
            {
                err = psr[p].obsn[i].toaErr;

		if(psr[p].TNGlobalEF > 0){
                           err *= psr[p].TNGlobalEF;
                }

		if(psr[p].TNGlobalEQ != 0){
			double TNEquad = pow(10.0,psr[p].TNGlobalEQ+6)*pow(10.0,psr[p].TNGlobalEQ+6);
                        err = sqrt(err*err + TNEquad);
                }

                for (j=0;j<psr[p].obsn[i].nFlags;j++)
                {
                    //Check efac
		

                    for (k=0;k<psr[p].nTNEF;k++)
                    {
                        if (strcmp(psr[p].obsn[i].flagID[j],psr[p].TNEFFlagID[k])==0)
                        {
                            if (strcmp(psr[p].obsn[i].flagVal[j],psr[p].TNEFFlagVal[k])==0)
                                err *= psr[p].TNEFVal[k];
                        }
                    }

                    //Check equad
                    for (k=0;k<psr[p].nTNEQ;k++)
                    {
                        if (strcmp(psr[p].obsn[i].flagID[j],psr[p].TNEQFlagID[k])==0)
                        {
                            if (strcmp(psr[p].obsn[i].flagVal[j],psr[p].TNEQFlagVal[k])==0){
                                double TNEquad = pow(10.0,psr[p].TNEQVal[k]+6)*pow(10.0,psr[p].TNEQVal[k]+6);
                                err = sqrt(err*err + TNEquad);
                            }
                        }
                    }
                    //Check Squad
                    for (k=0;k<psr[p].nTNSQ;k++)
                    {
                        if (strcmp(psr[p].obsn[i].flagID[j],psr[p].TNSQFlagID[k])==0)
                        {
                            if (strcmp(psr[p].obsn[i].flagVal[j],psr[p].TNSQFlagVal[k])==0){

                                double tobsval=0;
                                for (int tf=0;tf<psr[p].obsn[i].nFlags;tf++){
                                    if(strcasecmp(psr[p].obsn[i].flagID[tf],"-tobs")==0){
                                        if(strcasecmp(psr[p].obsn[i].flagVal[tf],"UNKNOWN")==0){
                                            tobsval=1;
                                        }
                                        else{
                                            double tobs=atof(psr[p].obsn[i].flagVal[tf]);
                                            tobsval=tobs;
                                        }
                                    }		
                                }
                                double TNSquad = pow(10.0,psr[p].TNSQVal[k]+6)*pow(10.0,psr[p].TNSQVal[k]+6)/tobsval;
                                err = sqrt(err*err + TNSquad);
                            }
                        }
                    }



                }
                psr[p].obsn[i].toaErr = err;
            }
        }

	
	// populate TOBS if using TNSECORR


	if (psr[p].nTNSECORR > 0)
	  {
	    for(i=0;i<psr[p].nobs;i++)
	      {
		double tobsval;
		
		for (int tf=0;tf<psr[p].obsn[i].nFlags;tf++){
		  if(strcasecmp(psr[p].obsn[i].flagID[tf],"-tobs")==0){
		    if(strcasecmp(psr[p].obsn[i].flagVal[tf],"UNKNOWN")==0){
		      tobsval=1;
		    }
		    else{
		      double tobs=atof(psr[p].obsn[i].flagVal[tf]);
		      tobsval=tobs;
		    }
		  }	
		}
		psr[p].obsn[i].tobs = tobsval;
	      }
	  }
	

        // Modify TOA flags if required
        if (modify==1)
        {
            FILE *fin;
            double mjd1,mjd2;
            char flag1[MAX_STRLEN],flag2[MAX_STRLEN],flag3[MAX_STRLEN];
            if (!(fin = fopen(modifyFname,"r")))
            {
                printf("Unable to open >%s< to modify the flags\n",modifyFname);
                exit(1);
            }
            while (!feof(fin))
            {
                if (fscanf(fin,"%s %s %lf %lf %s",flag1,flag2,&mjd1,&mjd2,flag3)==5)
                {
                    for (i=0;i<psr[p].nobs;i++)
                    {
                        for (j=0;j<psr[p].obsn[i].nFlags;j++)
                        {
                            if (strcmp(psr[p].obsn[i].flagID[j],flag1)==0 &&
                                    strcmp(psr[p].obsn[i].flagVal[j],flag2)==0 &&
                                    (double)psr[p].obsn[i].sat > mjd1 &&
                                    (double)psr[p].obsn[i].sat < mjd2)
                                strcpy(psr[p].obsn[i].flagVal[j],flag3);
                        }
                    }
                }
            }
            fclose(fin);
        }
        preProcessSimple2 (psr + p, startdmmjd, ndm, dmvals, trimonly);
        // Check for select file
        if (strlen(selectPlugName) > 0)
        {
            char *(*entry)(int,char **,pulsar *,int *);
            void * module;
            logdbg("Dealing with the select plugin"); 


            char *str = new char[2600];
            for (int iplug=0; iplug < tempo2_plug_path_len; iplug++){
                sprintf(str,"%s/%s_%s_splug.t2",tempo2_plug_path[iplug],
                        selectPlugName,tempo2MachineType);
                logmsg("Looking for %s",str);
                module = dlopen(str, RTLD_NOW|RTLD_GLOBAL); 
                if(module==NULL){	  
                    printf("dlerror() = %s\n",dlerror());
                } else break;
            }
            delete[] str;

            //	  strcpy(tempo2MachineType, getenv("LOGIN_ARCH"));
            //	  sprintf(str,"%s/plugins/%s_%s_splugt2",getenv(TEMPO2_ENVIRON),
            //		  selectPlugName,tempo2MachineType);
            //	  printf("Looking for %s\n",str);
            //	  module = dlopen(str, RTLD_NOW|RTLD_GLOBAL); 
            if(!module)  {
                fprintf(stderr, "[error]: dlopen() failed while resolving symbols.\n" );
                fprintf(stderr, "dlerror() = %s\n",dlerror());
                exit(1);
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
                    exit(1);
                }
            }


            entry = (char*(*)(int,char **,pulsar *,int *))dlsym(module, "selectInterface");
            if( entry == NULL ) {
                dlclose(module);
                fprintf(stderr, "[error]: dlerror() failed while  retrieving address.\n" ); 
                fprintf(stderr, "dlerror() = %s\n",dlerror());
                exit(1);
            }
            entry(argc,argv,psr,&npsr);
        }
        else
        {
            if (strlen(selectFname) > 0)
            {
                logdbg("Using select file");
                useSelectFile(selectFname,psr,npsr);
                logdbg("Complete using select file");
            }
        }
        preProcessSimple3 (psr + p);
    }

    // Now check for global parameters
    if (strlen(globalFname) > 0)
    {
        char tpar[MAX_STRLEN][MAX_FILELEN];
        char ttim[MAX_STRLEN][MAX_FILELEN];

        sprintf(tpar[0],"%s",globalFname);
        printf("Setting global parameters\n");
        readParfileGlobal(psr,npsr,tpar,ttim);
        printf("Complete setting global parameters\n");
    }
} 

void useSelectFile(char *fname,pulsar *psr,int npsr)
{
    int i,j,k,p;
    char line[1000];
    char first[1000];
    char second[1000];
    double v1,v2;
    int nread;
    FILE *fin;

    if (!(fin = fopen(fname,"r")))
    {
        printf("Unable to open select file: >%s<\n",fname);
        exit(1);
    }
    while (!feof(fin))
    {
        fgets(line,1000,fin);
        nread = sscanf(line,"%s %s",first,second);
        if (nread==2)
        {
            if (strcmp(first,"PASS")==0 && strcmp(second,"MJD")==0)
            {
                sscanf(line,"%s %s %lf %lf",first,second,&v1,&v2);
                for (p=0;p<npsr;p++)
                {
                    for (i=0;i<psr[p].nobs;i++)
                    {
                        if (psr[p].obsn[i].deleted==0 &&
                                (
                                 (double)psr[p].obsn[i].sat < v1 ||
                                 (double)psr[p].obsn[i].sat > v2 )
                           )
                            psr[p].obsn[i].deleted=1;
                    }
                }
            }
            else if (strcmp(first,"REJECT")==0 && strcmp(second,"MJD")==0)
            {
                sscanf(line,"%s %s %lf %lf",first,second,&v1,&v2);
                for (p=0;p<npsr;p++)
                {
                    for (i=0;i<psr[p].nobs;i++)
                    {
                        if ((double)psr[p].obsn[i].sat >= v1 &&
                                (double)psr[p].obsn[i].sat <= v2 &&
                                psr[p].obsn[i].deleted==0)
                            psr[p].obsn[i].deleted=1;			
                    }
                }
            }
            else if (strcmp(first,"PASS")==0 && strcmp(second,"TOAERR")==0)
            {
                sscanf(line,"%s %s %lf %lf",first,second,&v1,&v2);
                for (p=0;p<npsr;p++)
                {
                    for (i=0;i<psr[p].nobs;i++)
                    {
                        if ( ( (double)psr[p].obsn[i].toaErr < v1 ||
                                    (double)psr[p].obsn[i].toaErr > v2 ) &&
                                psr[p].obsn[i].deleted==0)
                            psr[p].obsn[i].deleted=1;			
                    }
                }
            }
            else if (strcmp(first,"REJECT")==0 && strcmp(second,"FREQ")==0)
            {
                sscanf(line,"%s %s %lf %lf",first,second,&v1,&v2);
                for (p=0;p<npsr;p++)
                {
                    for (i=0;i<psr[p].nobs;i++)
                    {
                        if ((double)psr[p].obsn[i].freq >= v1 &&
                                (double)psr[p].obsn[i].freq <= v2 &&
                                psr[p].obsn[i].deleted==0)
                            psr[p].obsn[i].deleted=1;			
                    }
                }
            }
            else if (strcmp(first,"PASS")==0 && strcmp(second,"FREQ")==0)
            {
                sscanf(line,"%s %s %lf %lf",first,second,&v1,&v2);
                for (p=0;p<npsr;p++)
                {
                    for (i=0;i<psr[p].nobs;i++)
                    {
                        if ( ( (double)psr[p].obsn[i].freq < v1 ||
                                    (double)psr[p].obsn[i].freq > v2 ) &&
                                psr[p].obsn[i].deleted==0)
                            psr[p].obsn[i].deleted=1;			
                    }
                }
            }
            else if (strcmp(first,"PASS")==0 && second[0]=='-')
            {
                char pos[1000];
                char *res;
                int nf=0;
                int found;
                char fi[20][100];
                strcpy(pos,line);
                strcpy(line,strtok(pos," "));
                strcpy(second,strtok(NULL," "));
                do
                {
                    res = strtok(NULL," ");
                    if (res!=NULL)
                    {
                        strcpy(fi[nf++],res);
                        if (fi[nf-1][strlen(fi[nf-1])-1] == '\n')
                            fi[nf-1][strlen(fi[nf-1])-1] = '\0';
                        if (nf == 20) 
                        {
                            printf("ERROR: too many flags in the select file\n");
                            exit(1);
                        }
                    }
                } while (res != NULL);
                //	      sscanf(line,"%s %s",first,second);
                for (p=0;p<npsr;p++)
                {
                    for (i=0;i<psr[p].nobs;i++)
                    {
                        if (psr[p].obsn[i].deleted==0)
                        {
                            found=0;
                            for (j=0;j<psr[p].obsn[i].nFlags;j++)
                            {
                                if (strcmp(psr[p].obsn[i].flagID[j],second)==0)
                                {
                                    for (k=0;k<nf;k++)
                                    {
                                        if (strcmp(psr[p].obsn[i].flagVal[j],fi[k])==0)
                                        {found=1; j=psr[p].obsn[i].nFlags; break;}
                                    }
                                }			       
                            }
                            if (found==0)
                                psr[p].obsn[i].deleted=1;			
                        }
                    }
                }
            }
            else if (strcmp(first,"REJECT")==0 && strcmp(second,"TOAERR")==0)
            {
                sscanf(line,"%s %s %lf %lf",first,second,&v1,&v2);
                for (p=0;p<npsr;p++)
                {
                    for (i=0;i<psr[p].nobs;i++)
                    {
                        if ((double)psr[p].obsn[i].toaErr >= v1 &&
                                (double)psr[p].obsn[i].toaErr <= v2 &&
                                psr[p].obsn[i].deleted==0)
                            psr[p].obsn[i].deleted=1;			
                    }
                }
            }
            else if (strcmp(first,"PROCESS")==0 && strcmp(second,"SIMUL")==0)
                processSimultaneous(line,psr,npsr);
            else if (strcmp(first,"PROCESS")==0 && second[0]=='-')
                processFlag(line,psr,npsr);
            else if (strcasecmp(first,"LOGIC")==0 && second[0]=='-')
                logicFlag(line,psr,npsr);
        }      
    }
    fclose(fin);
}

// Function to process specified flags
void logicFlag(char *line,pulsar *psr,int npsr)
{
    char t1[100],flagID[100],s1[100],s2[100],s3[100];
    int p,i,k,found;
    int fVal=0;

    sscanf(line,"%s %s %s %s %s",t1,flagID,s1,s2,s3);
    for (p=0;p<npsr;p++)
    {
        for (i=0;i<psr[p].nobs;i++)
        {
            found=0;
            for (k=0;k<psr[p].obsn[i].nFlags;k++)
            {
                if (strcmp(psr[p].obsn[i].flagID[k],flagID)==0 &&
                        psr[p].obsn[i].deleted==0)
                {found=1; fVal=k; break;}
            }
            if (found==0) // Haven't found an observation
            {
                if (strcasecmp(s1,"EXIST")==0)
                {
                    if (strcasecmp(s2,"PASS")==0)
                        psr[p].obsn[i].deleted=1;
                }	      
            }
            else if (found==1) // Got a point to process
            {
                double val;
                if (strcasecmp(s1,"EXIST")==0)
                {
                    if (strcasecmp(s2,"REJECT")==0)
                        psr[p].obsn[i].deleted=1;
                }	      
                else if (sscanf(psr[p].obsn[i].flagVal[fVal],"%lf",&val)==1)
                {
                    if (strcasecmp(s3,"REJECT")==0)
                    {
                        if (strcmp(s1,">")==0)
                        {
                            double val2;
                            sscanf(s2,"%lf",&val2);
                            if (val > val2)
                                psr[p].obsn[i].deleted=1;
                        }
                        else if (strcmp(s1,"<")==0)
                        {
                            double val2;
                            sscanf(s2,"%lf",&val2);
                            if (val < val2)
                                psr[p].obsn[i].deleted=1;
                        }
                        else if (strcmp(s1,"=")==0)
                        {
                            double val2;
                            sscanf(s2,"%lf",&val2);
                            if (val == val2)
                                psr[p].obsn[i].deleted=1;
                        }
                    }
                }
            }
        }
    }
}

// Function to process specified flags
void processFlag(char *line,pulsar *psr,int npsr)
{
    char t1[100],flagID[100],flagV[100],s1[100];
    double v1,v2;
    int p,i,k,found;

    sscanf(line,"%s %s %s %s %lf %lf",t1,flagID,flagV,s1,&v1,&v2);
    for (p=0;p<npsr;p++)
    {
        for (i=0;i<psr[p].nobs;i++)
        {
            found=0;
            for (k=0;k<psr[p].obsn[i].nFlags;k++)
            {
                if (strcmp(psr[p].obsn[i].flagID[k],flagID)==0 &&
                        strcmp(psr[p].obsn[i].flagVal[k],flagV)==0 &&
                        psr[p].obsn[i].deleted==0)
                {found=1; break;}
            }
            if (found==1) // Got a point to process
            {
                if (strcasecmp(s1,"mjdpass")==0 && ((double)psr[p].obsn[i].sat < v1 ||
                            (double)psr[p].obsn[i].sat > v2))
                    psr[p].obsn[i].deleted=1;
                if (strcasecmp(s1,"reject")==0)
                    psr[p].obsn[i].deleted=1;
                if (strcasecmp(s1,"mjdreject")==0 && ((double)psr[p].obsn[i].sat >= v1 &&
                            (double)psr[p].obsn[i].sat <= v2))
                    psr[p].obsn[i].deleted=1;
                if (strcasecmp(s1,"errpass")==0 && ((double)psr[p].obsn[i].toaErr < v1 ||
                            (double)psr[p].obsn[i].toaErr > v2))
                    psr[p].obsn[i].deleted=1;
                if (strcasecmp(s1,"errreject")==0 && ((double)psr[p].obsn[i].toaErr >= v1 &&
                            (double)psr[p].obsn[i].toaErr <= v2))
                    psr[p].obsn[i].deleted=1;
                if (strcasecmp(s1,"freqpass")==0 && ((double)psr[p].obsn[i].freq < v1 ||
                            (double)psr[p].obsn[i].freq > v2))
                    psr[p].obsn[i].deleted=1;
                if (strcasecmp(s1,"freqreject")==0 && ((double)psr[p].obsn[i].freq >= v1 &&
                            (double)psr[p].obsn[i].freq <= v2))
                    psr[p].obsn[i].deleted=1;
            }
        }
    }

}

//
// Function to process simultaneous points
//
void processSimultaneous(char *line,pulsar *psr, int npsr)
{
    int i,j,k,p,nread;
    int found1,found2,simul;
    char s1[100],s2[100],s3[100],s4[100],t1[100],t2[100];
    double tdiff;
    nread = sscanf(line,"%s %s %s %s %s %s %lf",t1,t2,s1,s2,s3,s4,&tdiff);
    if (nread!=6 && nread!=7)
    {
        printf("Problem in processing line: %s\n",line);
        return;
    }
    if (nread==6)
        tdiff=60;

    if (strcmp(t1,"PROCESS")==0 && strcmp(t2,"SIMUL")==0 && strcasecmp(s4,"only")==0) // Only identify simultaneous points
    {
        for (p=0;p<npsr;p++)
        {
            for (i=0;i<psr[p].nobs;i++)
            {
                simul=0;
                for (j=0;j<psr[p].nobs;j++)
                {
                    if (i!=j && (fabs(psr[p].obsn[i].sat - psr[p].obsn[j].sat)*SECDAY < tdiff))
                    {
                        found1=0;
                        found2=0;
                        for (k=0;k<psr[p].obsn[i].nFlags;k++)
                        {
                            if (strcmp(psr[p].obsn[i].flagID[k],s1)==0 &&
                                    strcmp(psr[p].obsn[i].flagVal[k],s2)==0 &&
                                    psr[p].obsn[i].deleted==0)
                            {found1=1; break;}
                        }
                        for (k=0;k<psr[p].obsn[j].nFlags;k++)
                        {
                            if (strcmp(psr[p].obsn[j].flagID[k],s1)==0 &&
                                    strcmp(psr[p].obsn[j].flagVal[k],s3)==0 &&
                                    psr[p].obsn[j].deleted==0)
                            {found2=1; break;}
                        }
                        if (found1==1 && found2==1)
                        {
                            simul=1;
                            break;
                        }

                        for (k=0;k<psr[p].obsn[i].nFlags;k++)
                        {
                            if (strcmp(psr[p].obsn[i].flagID[k],s1)==0 &&
                                    strcmp(psr[p].obsn[i].flagVal[k],s3)==0 &&
                                    psr[p].obsn[i].deleted==0)
                            {found1=1; break;}
                        }
                        for (k=0;k<psr[p].obsn[j].nFlags;k++)
                        {
                            if (strcmp(psr[p].obsn[j].flagID[k],s1)==0 &&
                                    strcmp(psr[p].obsn[j].flagVal[k],s2)==0 &&
                                    psr[p].obsn[j].deleted==0)
                            {found2=1; break;}
                        }
                        if (found1==1 && found2==1)
                        {
                            simul=1;
                            break;
                        }

                    }
                }
                if (simul==0)
                    psr[p].obsn[i].deleted=1;
            }
        }
    }
    else
    {
        for (p=0;p<npsr;p++)
        {
            for (i=0;i<psr[p].nobs;i++)
            {
                for (j=0;j<psr[p].nobs;j++)
                {
                    if (i!=j && (fabs(psr[p].obsn[i].sat - psr[p].obsn[j].sat)*SECDAY < tdiff))
                    {		    
                        found1=0;
                        found2=0;
                        for (k=0;k<psr[p].obsn[i].nFlags;k++)
                        {
                            if (strcmp(psr[p].obsn[i].flagID[k],s1)==0 &&
                                    strcmp(psr[p].obsn[i].flagVal[k],s2)==0 &&
                                    psr[p].obsn[i].deleted==0)
                            {found1=1; break;}
                        }
                        for (k=0;k<psr[p].obsn[j].nFlags;k++)
                        {
                            if (strcmp(psr[p].obsn[j].flagID[k],s1)==0 &&
                                    strcmp(psr[p].obsn[j].flagVal[k],s3)==0 &&
                                    psr[p].obsn[j].deleted==0)
                            {found2=1; break;}
                        }
                        if (found1==1 && found2==1)
                        {
                            if (strcasecmp(s4,"toaerr")==0)
                            {
                                if (psr[p].obsn[i].toaErr < psr[p].obsn[j].toaErr)
                                    psr[p].obsn[j].deleted=1;
                                else
                                    psr[p].obsn[i].deleted=1;
                            }
                            else if (strcasecmp(s4,"1")==0)
                                psr[p].obsn[j].deleted=1;
                            else if (strcasecmp(s4,"2")==0)
                                psr[p].obsn[i].deleted=1;
                        }
                    }
                }
            }
        }
    }
}

void readWhiteNoiseModelFile(pulsar *psr,int p)
{
    int i,j;
    FILE *fin;
    char str[128];
    char f1_id[128];
    char f1_val[128];
    char f2_id[128];
    double e0;
    double tobs;

    if (!(fin = fopen(psr[p].whiteNoiseModelFile,"r")))
    {
        printf("ERROR: Unable to open file white noise model file >%s<\n",psr[p].whiteNoiseModelFile);
        exit(1);
    }
    printf("Reading white noise model from >%s<\n",psr[p].whiteNoiseModelFile);
    while (!feof(fin))
    {
        if (fscanf(fin,"%s",str)==1)
        {
            if (strcmp(str,"JITTER")==0)
            {
                double t0,w0,wn;
                int process;
                fscanf(fin,"%s",f1_id);
                fscanf(fin,"%s",f1_val);
                fscanf(fin,"%s",f2_id);
                fscanf(fin,"%lf",&w0);
                fscanf(fin,"%lf",&t0);
                for (i=0;i<psr[p].nobs;i++)
                {
                    // Check if we should process
                    process = 0;
                    for (j=0;j<psr[p].obsn[i].nFlags;j++)
                    {
                        if (strcmp(psr[p].obsn[i].flagID[j],f1_id)==0 &&
                                strcmp(psr[p].obsn[i].flagVal[j],f1_val)==0)
                        {
                            process=1; 
                            break;
                        }
                    }
                    if (process==1)
                    {
                        // Calculate TOBS
                        tobs = -1;
                        for (j=0;j<psr[p].obsn[i].nFlags;j++)
                        {
                            if (strcmp(psr[p].obsn[i].flagID[j],f2_id)==0)
                                sscanf(psr[p].obsn[i].flagVal[j],"%lf",&tobs); 
                        }
                        if (tobs > -1)
                        {
                            e0 = psr[p].obsn[i].toaErr*1e-6;
                            wn = w0*sqrt(t0/tobs);
                            psr[p].obsn[i].toaErr = sqrt(pow(e0,2)+pow(wn,2))/1e-6; 
                        }
                    }
                }
            }
            else if (strcmp(str,"START_TIME_JITTER")==0)
            {
                double w0;
                int process;
                fscanf(fin,"%s",f1_id);
                fscanf(fin,"%s",f1_val);
                fscanf(fin,"%lf",&w0);
                for (i=0;i<psr[p].nobs;i++)
                {
                    // Check if we should process
                    process = 0;
                    for (j=0;j<psr[p].obsn[i].nFlags;j++)
                    {
                        if (strcmp(psr[p].obsn[i].flagID[j],f1_id)==0 &&
                                strcmp(psr[p].obsn[i].flagVal[j],f1_val)==0)
                        {
                            process=1; 
                            break;
                        }
                    }
                    if (process==1)
                    {
                        e0 = psr[p].obsn[i].toaErr*1e-6;
                        psr[p].obsn[i].toaErr = sqrt(pow(e0,2)+pow(w0,2))/1e-6; 
                    }
                }
            }
            else if (strcmp(str,"TOAERR_ZERO")==0)
            {
                for (i=0;i<psr[p].nobs;i++)
                    psr[p].obsn[i].toaErr = 0.0;
            }
        }
    }
    fclose(fin);

}
