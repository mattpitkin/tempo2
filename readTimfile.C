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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "tempo2.h"

/* ******************************************** */
/* read_timfile                                 */
/* Author:  G. Hobbs (02 May 2003)              */
/* Purpose: Reads a .tim file and fills         */
/*          observation structure               */
/* Inputs:  timFile - filename, obsn - structure*/
/*          of observations, nObs number of     */
/*          lines in .tim file                  */
/* Outputs: Fills obsn and nObs                 */
/*                                              */
/* Notes:   Default is to use the new .tim      */
/*          format file, but can read the old   */
/*          Parkes and Jodrell tempo format.    */ 
/*                                              */
/* Changes:                                     */
/* ******************************************** */
void readTim(char *timname,pulsar *psr,int *jumpVal, int *fdjumpVal);
void removeCR2(char *str);

void readTimfile(pulsar *psr,char timFile[][MAX_FILELEN],int npsr)
{
    int p,i;
    int jumpVal=0;
    int fdjumpVal=0;
    FILE *fin;
    const char *CVS_verNum = "$Id$";

    if (displayCVSversion == 1) CVSdisplayVersion("readTimfile.C","readTimfile()",CVS_verNum);

    logdbg("In reading tim file");

    for (p=0;p<npsr;p++)
    {

        // Next five lines have been moved to initialise.C
        //psr[p].nobs=0;
        //psr[p].dmOffset=0;
        //psr[p].nT2efac = 0;
        //psr[p].nT2equad = 0;
        //psr[p].T2globalEfac = 1;
        /*      jumpVal=psr[p].nJumps;
                for (i=0;i<jumpVal;i++)
                psr[p].jumpVal[i]=0.0; */
        jumpVal=0;
        fdjumpVal=0;

        if (psr[0].jboFormat!=0){
            logmsg("Reading JBO format!?");
            readJBO_bat(timFile[p],&psr[p],p);
        }
        else
            readTim(timFile[p],&psr[p],&jumpVal,&fdjumpVal);
        logdbg("Checking for deleted points >%s<",psr[0].deleteFileName);

        /* Check for deleted points in separate file */
        if (strcmp(psr[0].deleteFileName,"NONE")!=0)
        {
            logdbg("In checking for deleted points");
            if (!(fin = fopen(psr[0].deleteFileName,"r")))
                printf("Warning: unable to open %s\n",psr[0].deleteFileName);
            else
            {
                longdouble input;
                while (!feof(fin))
                {
                    char inputstr[1024];
                    if (fscanf(fin,"%1023s",inputstr)==1)
                    {
                        input = parse_longdouble(inputstr);
                        for (i=0;i<psr[p].nobs;i++)
                        {
                            if (fabs(psr[p].obsn[i].sat-input)<0.001)
                                psr[p].obsn[i].deleted=1;
                        }
                    }
                }
            }
        }
    }
    logdbg("Leaving readTimfile");  
}

void readTim(char *timname,pulsar *psr,int *jumpVal, int *fdjumpVal)
{
    FILE *fin;
    char *profileDir = new char[MAX_STRLEN];
    char *tt = new char[2*MAX_STRLEN+16];
    int nread=0,nread2,nObs=0,i;
    char *firstWord = new char[1000];
    char *line = new char[1000];
    char *dummy = new char[1000];

    strncpy(profileDir,"",1000);
    strncpy(line,"",1000);

    char param1[100];//,param2[100],param3[100],param4[100],param5[100];
    //char param6[100],param7[100],param8[100];
    //  int val1;
    int format,endit=0;
    int valid;
    int skip=0;

    // NOTE: These were static -- now they are not!
    double efac=1.0;
    double eset=-1.0;
    double emax=-1.0;  /* Maximum error used */
    double efloor=-1.0;
    double emin=-1.0;
    double equad = 0.0;
    double fmax=-1.0;
    double fmin=-1.0;
    double sigma=0.0;  /* Forced setting of TOA error */
    double time=0.0;
    int    add=0;
    static int infoNum=-1;
    char oldLine[1000];

    /* Attempt to open .tim file, exit if not possible */
    if (!(fin = fopen(timname,"r")))
    {
        printf("ERROR [FILE4]: Unable to open timfile >%s<\n",timname);
        exit(1);
    }
    /* Have successfully opened the file */

    /* Attempt to determine the file format */
    if (psr->fixedFormat > 0)
    {
        int k;
        format=1;
        for (k=0;k<psr->fixedFormat;k++)
            fgets(oldLine,1000,fin);
    }
    else
    {
        format=1;
        while (!feof(fin))
        {
            fscanf(fin,"%999s",firstWord);
            if (strcmp(firstWord,"FORMAT")==0) /* Tempo2 format */
                format = 0;
            else if (strcmp(firstWord,"HEAD")==0) /* .tpo format */
                format = 2;
        }
        fclose(fin);
        fin = fopen(timname,"r");  /* Could use rewind instead */
    }
    if (format==2)
    {
        do {
            fgets(line,1000,fin); removeCR2(line);
        } while (strcasecmp(line,"TOAS")!=0 && feof(fin)==0);
        format=1;

    }
    while (!feof(fin) && endit==0)
    {      
        valid=0;
        if (format==0) /* New Tempo2 (pseudo-)free format file format */
        {
            if (nread<5 && strlen(line)>1)
            {
                line[strlen(line)-1]='\0';
                /* printf("Warning: >%s< Not a complete line in .tim file. Ignoring line\n",line); */
            }
            nObs = psr->nobs;
            if (fgets(line,1000,fin)==NULL)             /* Read line from .tim file */
                valid=-2;
            //	  printf("Have read %s\n",line);

            /* Remove whitespace (or eol) characters at end of line */
            while (strlen(line) > 0 && isspace(line[strlen(line)-1])) line[strlen(line)-1]='\0';
            if (line[0]!='C' && line[0]!='#') /* Ignore comment lines     */
            {
                /* Read default columns */
                if (strlen(line)>0)
                {
                    char* sat_str = new char[1024];
                    psr->obsn[nObs].deleted=0;

                    if (line[0]=='I')
                    {
                        nread = sscanf(line,"%999s %1024s %lf %1023s %lf %99s",dummy,psr->obsn[nObs].fname,
                                &(psr->obsn[nObs].freq),sat_str,
                                &(psr->obsn[nObs].toaErr),psr->obsn[nObs].telID);	      
                        psr->obsn[nObs].deleted=1;
                    }		  
                    else
                    {
                        nread = sscanf(line,"%1024s %lf %1023s %lf %99s",psr->obsn[nObs].fname,
                                &(psr->obsn[nObs].freq),sat_str,
                                &(psr->obsn[nObs].toaErr),psr->obsn[nObs].telID);
                        if (strlen(profileDir)>0)
                        {
                            sprintf(tt,"%s/%s",profileDir,psr->obsn[nObs].fname);
                            strncpy(psr->obsn[nObs].fname,tt,500);
                        }
                    }

                    if (nread >=5) {
                        // try to parse SAT etc only if there are enough columns... otherwise it is junk

                        char sat_day_str[8]; // This must be >= 6 because  it needs a null terminator
                        char sat_sec_str[1024];

                        sat_sec_str[0] = '0';


                        /*
                         * This string copy is problematic as it doesn't terminate the strings properly
                         * Better to just use strncpy
                          for(int sindex = 0; sindex < 1024; sindex++){
                            if(sindex < 5){
                                sat_day_str[sindex] = sat_str[sindex];
                            }
                            if(sindex>=5){
                                sat_sec_str[sindex-4] = sat_str[sindex];
                            }
                        }
                        */
                        strncpy(sat_day_str,sat_str,5);
                        sat_day_str[5]='\0'; // we expect sat_str to be longer than 5, so we have to null terminate
                        // see definition of strncpy
                        strncpy(sat_sec_str+1,sat_str+5,1000);

                        psr->obsn[nObs].sat_day = parse_longdouble(sat_day_str);
                        psr->obsn[nObs].sat_sec = parse_longdouble(sat_sec_str);

                        psr->obsn[nObs].sat = parse_longdouble(sat_str);
                        psr->obsn[nObs].phaseOffset = 0.0;
                        /* Read the rest of the line */		  
                        psr->obsn[nObs].nFlags = 0;
                        psr->obsn[nObs].toaDMErr = 0;

                        /*		  strcpy(psr->obsn[nObs].flagID[0],"FLAGID");
                                  strcpy(psr->obsn[nObs].flagVal[0],"FLAGVAL"); */

                        for (i=0;i<(int)strlen(line)-1;i++)
                        {
                            if (line[i]=='-' && (line[i+1] < 48 || line[i+1] > 57))
                            {
                                strcpy(oldLine,line);
                                if (strchr(line+i,' ')!=NULL)
                                {
                                    strcpy(strchr(line+i,' '),"");
                                    strcpy(psr->obsn[nObs].flagID[psr->obsn[nObs].nFlags],line+i);
                                    i+=strlen(line+i)+1;
                                    strcpy(line,oldLine);
                                    if (strchr(line+i,' ')!=NULL)
                                    {
                                        int j;
                                        for (j=i;j<(int)strlen(line);j++)
                                        {
                                            if (line[j]!=' ')
                                            {
                                                i+=(j-i);
                                                break;
                                            }
                                        }
                                        strcpy(strchr(line+i,' '),"");
                                    }
                                    strcpy(psr->obsn[nObs].flagVal[psr->obsn[nObs].nFlags],line+i);
                                    // Check for offset to site arrival time
                                    if (strcmp(psr->obsn[nObs].flagID[psr->obsn[nObs].nFlags],"-addsat")==0)
                                    {
                                        double offVal;
                                        sscanf(psr->obsn[nObs].flagVal[psr->obsn[nObs].nFlags],"%lf",&offVal);
                                        psr->obsn[nObs].sat += (offVal/longdouble(86400.0));
                                    }

                                    // Check for DM changes
                                    if (strcmp(psr->obsn[nObs].flagID[psr->obsn[nObs].nFlags],"-dme")==0)
                                    {
                                        double dme,dm,freq,toaCorr;
                                        sscanf(psr->obsn[nObs].flagVal[psr->obsn[nObs].nFlags],"%lf",&dme);
                                        sscanf(psr->obsn[nObs].flagVal[psr->obsn[nObs].nFlags-1],"%lf",&dm);
                                        freq = psr->obsn[nObs].freq*1e6;

                                        toaCorr = dme/DM_CONST/1.0e-12/freq/freq;
                                        //				  printf("Have DM error %g %g %.5g %.5g %g\n",dm,dme,toaErr,freq,toaCorr);
                                        psr->obsn[nObs].toaDMErr = toaCorr*1e6;
                                        //				  psr->obsn[nObs].toaErr = sqrt(pow(psr->obsn[nObs].toaErr,2)+pow(toaCorr*1e6,2));
                                    }
                                    i+=strlen(line+i);
                                    strcpy(line,oldLine); 
                                    psr->obsn[nObs].nFlags++;
                                    if (psr->obsn[nObs].nFlags >= MAX_FLAGS)
                                    {
                                        printf("Number of different flags in the .tim file > MAX_FLAGS (%d)\n",MAX_FLAGS);
                                        exit(1);
                                    }
                                }
                            }
                        }
                        delete[] sat_str;
                    }
                }
                else
                    nread=0;

                /* Now check whether this is a valid observation */
                if (nread>=5 && valid==0) valid=1;
            }
            else
                nread=5;	      
        }
        else if (format==1) /* Parkes, princeton or ITOA file format */
        {
            if (fgets(line,1000,fin)==NULL)             /* Read line from .tim file */
                valid=-2;
            else
            {	    
                /* Remove \n character at end */
                if (line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0';
                nObs = psr->nobs;
                psr->obsn[nObs].nFlags = 0;
                nread = sscanf(line,"%999s",param1);
                add=0;

                if (strcmp(param1,"C")==0 || strcmp(param1,"c")==0) /* Comment line */
                {
                    valid=0;
                    /*		  add=1; */ /* NOW IGNORE COMMENTED LINES COMPLETELY */
                    /*		  psr->obsn[nObs].deleted=1; */
                }
                else
                {
                    psr->obsn[nObs].deleted=0;

                    if (strcasecmp(param1,"INCLUDE")!=0)
                    {
                        if (line[0+add]==' ')      /* Parkes format */
                        {
                            valid=1;
                            if (strlen(line)+add < 79) valid=-2;
                            else
                            {
                                strncpy(psr->obsn[nObs].fname,line+1+add,500); psr->obsn[nObs].fname[25]='\0';
                                strncpy(param1,line+25+add,99); param1[9]='\0'; if (strlen(param1)<2) valid=-2; 

                                if (sscanf(param1,"%lf",&(psr->obsn[nObs].freq))!=1) valid=-2;
                                strncpy(param1,line+34+add,99); param1[21]='\0';
                                psr->obsn[nObs].sat = parse_longdouble(param1);
                                strncpy(param1,line+55+add,99); param1[8]='\0';
                                if (sscanf(param1,"%lf",&(psr->obsn[nObs].phaseOffset))!=1) valid=-2;
                                strncpy(param1,line+63+add,99); param1[8]='\0';
                                if (sscanf(param1,"%lf",&(psr->obsn[nObs].toaErr))!=1) valid=-2;
                                sscanf(line+79+add,"%99s",psr->obsn[nObs].telID); psr->obsn[nObs].telID[1]='\0';
                            }
                        }
                        else if (line[1+add]==' ') /* Princeton format */
                        {
                            double dmoffset;
                            valid=1;
                            sscanf(line+add,"%99s",psr->obsn[nObs].telID); psr->obsn[nObs].telID[1]='\0';	
                            strcpy(psr->obsn[nObs].fname,"NOT SET");
                            strncpy(param1,line+15+add,99); param1[9]='\0';
                            if (sscanf(param1,"%lf",&(psr->obsn[nObs].freq))!=1) valid=-2;
                            strncpy(param1,line+25+add,99); param1[20]='\0';
                            psr->obsn[nObs].sat = parse_longdouble(param1);
                            strncpy(param1,line+45+add,99); param1[9]='\0';
                            if (sscanf(param1,"%lf",&(psr->obsn[nObs].toaErr))!=1) valid=-2;
                            strncpy(param1,line+68+add,99); param1[10]='\0'; /* SHOULD BE DM OFFSET */
                            if (sscanf(param1,"%lf",&dmoffset)==1)
                            {
                                strcpy(psr->obsn[nObs].flagID[psr->obsn[nObs].nFlags],"-dmo");
                                sprintf(psr->obsn[nObs].flagVal[psr->obsn[nObs].nFlags],"%.10g",(double)(dmoffset)); 
                                psr->obsn[nObs].nFlags++;

                            }
                        }
                        else if (strlen(line)>14 && line[14+add]=='.')  /* ITOA format */
                        {
                            valid=1;
                            strcpy(psr->obsn[nObs].fname,"NOT SET");		      
                            strcpy(param1,line+9+add); param1[19]='\0';
                            psr->obsn[nObs].sat = parse_longdouble(param1);
                            strcpy(param1,line+28+add); param1[6]='\0';
                            if (sscanf(param1,"%lf",&(psr->obsn[nObs].toaErr))!=1) valid=-2;
                            strcpy(param1,line+34+add); param1[11]='\0';
                            if (sscanf(param1,"%lf",&(psr->obsn[nObs].freq))!=1) valid=-2;
                            strcpy(param1,line+45+add); param1[10]='\0'; /* SHOULD BE DM OFFSET */
                            if (sscanf(param1,"%lf",&(psr->obsn[nObs].phaseOffset))!=1) valid=-2;
                            sscanf(line+57+add,"%99s",psr->obsn[nObs].telID); psr->obsn[nObs].telID[2]='\0';
                        }
                    }
                }
            }
        } // end of "else if (format == 1)"
        if (valid==1)  /* Check for validity of this TOA */
        {
            /* The above comment is misleading. "valid==1" implies that the TOA _is_ valid. 
               However, below, the observation undergoes a series of actions that 
               - Should be performed on all observations 
               - Are totally independent of the file format.

               Just so you know. 
               */
            if (psr->obsn[nObs].telID[0]=='@'            /* Have barycentric arrival time */
                    || strcasecmp(psr->obsn[nObs].telID,"bat")==0) 
            {
                psr->obsn[nObs].clockCorr=0;  /* therefore don't do clock corrections */
                psr->obsn[nObs].delayCorr=0;
            }
            else if (strcmp(psr->obsn[nObs].telID,"STL")==0)
            {
                psr->obsn[nObs].clockCorr=0;  /* don't do clock corrections */
                psr->obsn[nObs].delayCorr=1;
            }
            else if (strcmp(psr->obsn[nObs].telID,"STL_FBAT")==0)
            {
                psr->obsn[nObs].clockCorr=0;  /* don't do clock corrections */
                psr->obsn[nObs].delayCorr=1;
            }
            else
            {
                psr->obsn[nObs].clockCorr=1;
                psr->obsn[nObs].delayCorr=1;
            }	  

            /* Check for conditionals */


            /* Note ordering: EQUAD BEFORE EFAC */
            /* Joris is going to re-write the equad business now. So stay tuned. 
               In case everything goes wrong, you can put the following two lines back again.
               if (equad!=0.0)
               psr->obsn[nObs].toaErr = sqrt(pow(psr->obsn[nObs].toaErr,2)+equad*equad); */

            psr->obsn[nObs].toaErr = sqrt(pow(psr->obsn[nObs].toaErr,2)+equad*equad); 
            psr->obsn[nObs].equad = equad;

            //printf("%lg  ",psr->obsn[nObs].equad);

            // Furthermore, I suspect that the "EQUAD on certain flags" should be put before the
            // efac adding, but well, who am I anyway?
            if (sigma!=0.0)
                psr->obsn[nObs].toaErr = sigma;
            psr->obsn[nObs].toaErr *= (efac);
            psr->obsn[nObs].efac = (efac);
            /*	  for (i=0;i<nequadFlag;i++)
                  {
                  for (k=0;k<psr->obsn[nObs].nFlags;k++)
                  {
                  if (strcmp(psr->obsn[nObs].flagID[k],equadFlagID[i])==0)
                  {
                  if (strcmp(psr->obsn[nObs].flagVal[k],equadFlag[i])==0)
                  psr->obsn[nObs].toaErr = sqrt(pow(psr->obsn[nObs].toaErr,2)+pow(equadFlagVal[i],2));
                  }
                  }
                  }*/
            if (eset>-1)
                psr->obsn[nObs].toaErr = eset;

            psr->obsn[nObs].jump[0] = *jumpVal; // GH: Added jump[0]. Mar 2011
            psr->obsn[nObs].obsNjump = 1;


            psr->obsn[nObs].fdjump[0] = *jumpVal; // RMS: Added fdjump[0]. Apr 2020
            psr->obsn[nObs].obsNfdjump = 1;

            if (time!=0.0) 
            {
                psr->obsn[nObs].sat += time/60.0/60.0/24.0;
                logdbg("psr %s adding time >%g<",psr->name,time);
            }
            if (efloor!=-1 && psr->obsn[nObs].toaErr < efloor) psr->obsn[nObs].toaErr = efloor;
            if (emax!=-1 && psr->obsn[nObs].toaErr > emax) psr->obsn[nObs].deleted = 1;
            if (emin!=-1 && psr->obsn[nObs].toaErr < emin) psr->obsn[nObs].deleted = 1;
            if (fmax!=-1 && psr->obsn[nObs].freq > fmax) valid=0;
            if (fmin!=-1 && psr->obsn[nObs].freq < fmin) valid=0;
            if (infoNum!=-1)
            {
                strcpy(psr->obsn[nObs].flagID[psr->obsn[nObs].nFlags],"-i");
                sprintf(psr->obsn[nObs].flagVal[psr->obsn[nObs].nFlags],"%d",infoNum);
                psr->obsn[nObs].nFlags++;
                if (psr->obsn[nObs].nFlags >= MAX_FLAGS)
                {
                    printf("Number of different flags in the .tim file > MAX_FLAGS (%d)\n",MAX_FLAGS);
                    exit(1);
                }

            }
            if (skip==1) valid=0;
            psr->obsn[nObs].origErr = psr->obsn[nObs].toaErr;

            if (valid==1)(psr->nobs)++;
            if (psr->nobs > MAX_OBSN-2)
            {
                fprintf(stderr, "Too many TOAs! Use the -npsr and -nobs command line arguments!\n");
                exit(1);
            }
        }
        else if (valid!=-2) // This means: if the line does contain information, but not an observation
        {
            nread2 = sscanf(line,"%99s",param1);
            if (skip==0 && nread2 > 0)
            {
                if (strcasecmp(param1,"END")==0)
                    endit=1;	
                /* Global error multiplying factor */
                else if (strcasecmp(param1,"PROFILE_DIR")==0)
                    sscanf(line,"%99s %999s",param1,profileDir);
                else if (strcasecmp(param1,"GLOBAL_EFAC")==0)
                    sscanf(line,"%99s %lf",param1,&psr->T2globalEfac);
                else if (strcasecmp(param1,"EFAC")==0)  /* Error multiplying factor */    
                    sscanf(line,"%99s %lf",param1,&efac);
                else if (strcasecmp(param1,"EFLOOR")==0)  /* Minimum error                        */    
                    sscanf(line,"%99s %lf",param1,&efloor);
                else if (strcasecmp(param1,"T2EFAC")==0) /* EFAC for given flag                    */
                {
                    int nefacFlag = psr->nT2efac;
                    sscanf(line,"%99s %31s %31s %lf",param1,psr->T2efacFlagID[nefacFlag],psr->T2efacFlagVal[nefacFlag],
                            &psr->T2efacVal[nefacFlag]);
                    (psr->nT2efac)++;
                }
                else if (strcasecmp(param1,"T2EQUAD")==0) /* EQUAD for given flag                  */
                {
                    int nequadFlag = psr->nT2equad;
                    sscanf(line,"%99s %31s %31s %lf",param1,psr->T2equadFlagID[nequadFlag],psr->T2equadFlagVal[nequadFlag],
                            &psr->T2equadVal[nequadFlag]);
                    (psr->nT2equad)++;
                }
                else if (strcasecmp(param1,"EMAX")==0)  /* Maximum error                           */
                    sscanf(line,"%99s %lf",param1,&emax);
                else if (strcasecmp(param1,"EMIN")==0)  /* Minimum error                           */
                    sscanf(line,"%99s %lf",param1,&emin);
                else if (strcasecmp(param1,"ESET")==0)  /* Set all errors to given amount          */
                    sscanf(line,"%99s %lf",param1,&eset);
                else if (strcasecmp(param1,"FMAX")==0)  /* Maximum observing frequency             */
                    sscanf(line,"%99s %lf",param1,&fmax);
                else if (strcasecmp(param1,"FMIN")==0)  /* Minimum observing frequency             */
                    sscanf(line,"%99s %lf",param1,&fmin);
                else if (strcasecmp(param1,"INFO")==0)  /* Highlighting flag                       */
                {
                    sscanf(line,"%99s %d",param1,&infoNum);
                }
                else if (strcasecmp(param1,"PHASE")==0) /* Add phase jump                          */
                {
                    psr->phaseJump[psr->nPhaseJump] = psr->obsn[nObs-1].sat+1.0/SECDAY;		  
                    psr->phaseJumpID[psr->nPhaseJump] = nObs-1;
                    sscanf(line,"%99s %d",param1,&psr->phaseJumpDir[psr->nPhaseJump]);
                    psr->nPhaseJump++;		  
                }
                else if (strcasecmp(param1,"EQUAD")==0) /* Error to add in quadrature               */
                    sscanf(line,"%99s %lf",param1,&equad);
                else if (strcasecmp(param1,"SIGMA")==0) /* Set all errors to constant value         */
                    sscanf(line,"%99s %lf",param1,&sigma);
                else if (strcasecmp(param1,"TIME")==0)  /* Add a constant time to all arrival times */
                {
                    double dtime;
                    sscanf(line,"%99s %lf",param1,&dtime);
                    time+=dtime;
                    logdbg("Updating time: %g %g",dtime,time);

                }
                else if (strcasecmp(param1,"MODE")==0) /* Fit with errors */
                {
                    sscanf(line,"%99s %d",param1,&(psr->fitMode));
                    displayMsg(1,"TIM1","Please place MODE flags in the parameter file","",psr->noWarnings);
                }
                else if (strcasecmp(param1,"INCLUDE")==0) /* Include another .tim file */
                {
                    char *newtim = new char[MAX_FILELEN];
#ifdef PATH_MAX
                    char *relPath = new char[PATH_MAX];
#else
                    char *relPath;
#endif
                    if (sscanf(line,"%99s %1023s",param1,newtim)==2)
                    {
                        int ii;
                        // Relative file path
                        printf("Current filename = %s\n",timname);
#ifdef PATH_MAX
                        realpath(timname,relPath);
#else
                        relPath = realpath(timname,NULL);
#endif
                        // Remove filename
                        for (ii=strlen(relPath);ii>0;ii--)
                        {
                            if (relPath[ii]=='/')
                            {
                                relPath[ii+1]='\0';
                                break;
                            }
                        }
                        strcat(relPath,newtim);
                        printf("Rel path = %s\n",relPath);
                        //		      readTim(newtim,psr,jumpVal);
                        readTim(relPath,psr,jumpVal,fdjumpVal);
#ifndef PATH_MAX
                        free(relPath);
#endif
                    }
                    else
                        printf("Unable to parse INCLUDE line >%s<\n",line);
                    strcpy(param1,"");
#ifdef PATH_MAX
                    delete[] relPath;
#endif
                    delete[] newtim;
                }
                else if (strcasecmp(param1,"SKIP")==0) /* Skip data */
                    skip=1;
                else if (strcasecmp(param1,"JUMP")==0) /* JUMP */
                    (*jumpVal)++;
                else if (strcasecmp(param1,"NOSKIP")==0) /* Stop skipping data */
                    skip=0;
            }
            else
            {
                if (strcasecmp(param1,"NOSKIP")==0) /* Stop skipping data */
                    skip=0;
            }
        } // end of "else if (valid != -2)"
    } // end of "while (!feof(fin) && endit == 0)"
    fclose(fin);

    delete[] profileDir;
    delete[] tt;
    delete[] firstWord;
    delete[] line;
    delete[] dummy;
}


/* Write an arrival time file in the tempo2 format */
void writeTim(const char *timname,pulsar *psr,const char *fileFormat)
{
    FILE *fout;
    int i,j;
    char name[1000];
    double current_efac;
    double current_equad;
    double interim_error;
    longdouble oldsat;

    if (!(fout = fopen(timname,"w")))
    {
        printf("Sorry, unable to write to file %s\n",timname);
        return;
    }
    if (strcmp(fileFormat,"tempo2")==0) fprintf(fout,"FORMAT 1\n");
    if (psr->fitMode==1) fprintf(fout,"MODE 1\n");    

    // Print out the Tempo2 EFACs and EQUADs and GLOBAL_EFAC
    if (psr[0].T2globalEfac != 1)
        fprintf(fout,"GLOBAL_EFAC %.5g\n",psr[0].T2globalEfac);
    // Now print out the EQUADS
    for (i=0;i<psr[0].nT2equad;i++)
        fprintf(fout,"T2EQUAD %s %s %.5g\n",psr[0].T2equadFlagID[i],psr[0].T2equadFlagVal[i],psr[0].T2equadVal[i]);
    // Now print out the EFACs
    for (i=0;i<psr[0].nT2efac;i++)
        fprintf(fout,"T2EFAC %s %s %.5g\n",psr[0].T2efacFlagID[i],psr[0].T2efacFlagVal[i],psr[0].T2efacVal[i]);

    current_efac = psr[0].obsn[0].efac;
    current_equad = psr[0].obsn[0].equad;
    if(current_equad!=0.0)
        fprintf(fout,"EQUAD %lf\n",current_equad);
    if(current_efac!=1.0)
        fprintf(fout,"EFAC %lf\n",current_efac);

    for (i=0;i<psr->nobs;i++)
    {
        if((double)psr[0].obsn[i].equad !=  current_equad){
            fprintf(fout,"EQUAD %lf\n",psr[0].obsn[i].equad);
            current_equad = psr[0].obsn[i].equad;
        }
        if(psr[0].obsn[i].efac!=current_efac){
            fprintf(fout,"EFAC %lf\n",psr[0].obsn[i].efac);
            current_efac = psr[0].obsn[i].efac;
        }
        if (psr[0].obsn[i].deleted==1)
            fprintf(fout,"C");

        if (strcmp(fileFormat,"tempo2")==0) 
        {
            ///////////////////////////////////////////////////////////////////
            ///////////// Undo what was done in preProcessSimple.C ////////////
            ///////////////////////////////////////////////////////////////////
            oldsat = psr->obsn[i].sat;

            // If we don't do this then createRealisation fails
            if (psr->param[param_f].val[0] > 0)
                oldsat -= ((psr->obsn[i].phaseOffset/psr->param[param_f].val[0])/SECDAY);

            // Correctly taking note of the phase offset
            for (int k=0;k<psr->nToffset;k++) /* Calculate time offsets */
            {
                char offsetSite[256], obsSite[256];
                lookup_observatory_alias(psr->tOffsetSite[k], offsetSite);
                lookup_observatory_alias(psr->obsn[i].telID, obsSite);
                if ((psr->tOffset_f1[k]==0.0 || (psr->obsn[i].freq > psr->tOffset_f1[k])) &&
                        (psr->tOffset_f2[k]==0.0 || (psr->obsn[i].freq < psr->tOffset_f2[k])) &&
                        (psr->tOffset_t1[k]==0.0 || (psr->obsn[i].sat  > psr->tOffset_t1[k])) &&
                        (psr->tOffset_t2[k]==0.0 || (psr->obsn[i].sat  < psr->tOffset_t2[k])) && 
                        (strcmp(offsetSite,"0")==0 || 
                         strcmp(offsetSite,obsSite)==0))
                {
                    int use=1;
                    /* Check for flags */
                    if (strlen(psr->tOffsetFlags[k])>0)
                    {
                        char *myStr,str1[1000],flagID[100];
                        use=0;
                        strcpy(str1,psr->tOffsetFlags[k]);
                        myStr = strtok(str1," ");
                        do {
                            if (myStr[0]=='-' && (myStr[1]<48 || myStr[1]>57)) /* Look for flag */
                            {
                                strcpy(flagID,myStr);
                                myStr = strtok(NULL," ");
                                for (int l=0;l<psr->obsn[i].nFlags;l++)
                                {
                                    if (strcmp(flagID,psr->obsn[i].flagID[l])==0
                                            && strcmp(myStr,psr->obsn[i].flagVal[l])==0)
                                    {
                                        use=1;
                                        break;
                                    }
                                }
                            }
                        } while ((myStr = strtok(NULL," "))!=NULL);
                    }
                    if (use==1) {oldsat -= psr->tOffset[k]/SECDAY;}
                }  
            }
            for (int k=0;k<psr->obsn[i].nFlags;k++)
            {
                if (strcmp(psr->obsn[i].flagID[k],"-to")==0)
                {
                    longdouble v;
                    v = parse_longdouble(psr->obsn[i].flagVal[k]);
                    oldsat -= v/SECDAY;
                }
            }

            ///////////////////////////////////////////////////////////////////
            //////// Done undoing what was done in preProcessSimple.C /////////
            ///////////////////////////////////////////////////////////////////


            sscanf(psr->obsn[i].fname,"%999s",name);
            interim_error = psr->obsn[i].origErr/current_efac;
            interim_error = sqrt(pow(interim_error,2.0)-(pow(current_equad,2.0)));

            ld_fprintf(fout," %s %.8lf %.17Lf %.5lf %s ", name,psr->obsn[i].freq,
                    oldsat, interim_error,psr->obsn[i].telID);
            if(interim_error<longdouble(0.0)){
                logerr("TOAerror < 0!!\n");
                exit(1);
            }

            /* Now add flags */
            for (j=0;j<psr->obsn[i].nFlags;j++)
            {
                if (strcasecmp(psr->obsn[i].flagID[j],"-selectAll")!=0)
                    fprintf(fout,"%s %s ",psr->obsn[i].flagID[j],psr->obsn[i].flagVal[j]);
            }
            fprintf(fout,"\n");
        }
        else
        {
            if (strlen(psr[0].obsn[i].fname) > 22)
                psr[0].obsn[i].fname[22]='\0';

            ld_fprintf(fout," %-25.25s%8.3f  %.13Lf    0.00 %7.2f        %s",
                    psr[0].obsn[i].fname,(double)psr[0].obsn[i].freq,
                    psr[0].obsn[i].sat,(double)psr[0].obsn[i].origErr,psr[0].obsn[i].telID);
            fprintf(fout,"\n");
        }

        //      for (j=0;j<psr->nPhaseJump;j++)
        //	{
        //	  if (i<psr->nobs && psr->obsn[i].sat < psr->phaseJump[j] && 
        //	      psr->obsn[i+1].sat >= psr->phaseJump[j])
        //	    fprintf(fout,"PHASE %d\n",psr->phaseJumpDir[j]);
        //	}

    }
    fclose(fout);
}

/* Removes newline at end of string */
void removeCR2(char *str)
{
    if (str[strlen(str)-1]=='\n') str[strlen(str)-1] = '\0';
}
