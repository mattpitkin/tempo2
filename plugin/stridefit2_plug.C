
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
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "t2fit.h"
#include "constraints.h"

using namespace std;

void help() /* Display help */
{
    printf("Options:\n");
    printf("-start mjd      : MJD to start\n");
    printf("-end mjd        : MJD to end\n");
    printf("-dt days        : step size\n");
    printf("-width days     : window width\n");
    printf("-param label    : stride fit for this parameter\n");
    printf("-delta          : output offset from orig fit.\n");
    printf("-delta_epoch    : output offset from orig fit, corrected for the change in epoch.\n");
    printf("-pars           : write .par files with start/finish flags.\n");
    printf("                  Note - this will remove F1, PM in .par file etc from output.\n");
    printf("                  \n");
    printf("Use e.g. -nofit to prevent initial fit.\n");
    exit(1);
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
    char parFile[MAX_PSR][MAX_FILELEN];
    char timFile[MAX_PSR][MAX_FILELEN];

#define MAX_P 30
#define MAX_E 1000
    double globalParameter=0.0;
    double dt=100;
    double width=100;
    double start=0;
    double end=999999;
    int nparam=0;
    param_label params[MAX_P];
    int param_k[MAX_P];
    double values[MAX_E][MAX_P];
    double errors[MAX_E][MAX_P];
    double starts[MAX_E];
    double ends[MAX_E];

    double ovalues[MAX_P];
    bool delta=false;
    bool delta_correct_epoch=false;
    char pars=0;

    const char *CVS_verNum = "$Revision: 1.1 $";

    if (displayCVSversion == 1) CVSdisplayVersion("autoDM.C","plugin",CVS_verNum);

    printf("Graphical Interface: autoDM\n");
    printf("Author:              M. Keith\n");
    printf("CVS Version:         $Revision: 1.1 $\n");

    /* Obtain all parameters from the command line */
    for (int i=2;i<argc;i++)
    {
        if (strcmp(argv[i],"-f")==0) {
            strcpy(parFile[0],argv[++i]); 
            strcpy(timFile[0],argv[++i]);
        } else if(strcmp(argv[i],"-dt")==0){
            dt=atof(argv[++i]);
        } else if(strcmp(argv[i],"-width")==0){
            width=atof(argv[++i]);
        } else if(strcmp(argv[i],"-start")==0){
            start=atof(argv[++i]);
        } else if(strcmp(argv[i],"-end")==0){
            end=atof(argv[++i]);
        } else if(strcmp(argv[i],"-delta")==0){
            delta=true;
        } else if(strcmp(argv[i],"-pars")==0){
            pars=1;
        } else if(strcmp(argv[i],"-delta_epoch")==0){
            delta=true;
            delta_correct_epoch=true;
        } else if(strcmp(argv[i],"-help")==0){
            help();
        } else if(strcmp(argv[i],"-param")==0){
            ++i;
            bool found=false;
            logmsg("Enabling parameter '%s'",argv[i]);
            for (param_label ip = 0; ip < param_LAST; ++ip){
                for (int k=0; k < psr[0].param[ip].aSize; ++k){
                    if (strcasecmp(psr[0].param[ip].shortlabel[k],argv[i])==0){
                        params[nparam]=ip;
                        param_k[nparam]=k;
                        ++nparam;
                        found=true;
                        break;
                    }
                }
                if (found)break;
            }
            if (!found){
                logerr("Cannot find parameter '%s'",argv[i]);
                exit(1);
            }
        }
    }
    if (npsr==0){
        help();
    }

    readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
    readTimfile(psr,timFile,*npsr); /* Load the arrival times    */

    psr->param[param_start].val[0] = start;
    psr->param[param_start].paramSet[0] = 1;
    psr->param[param_start].fitFlag[0] = 1;
    psr->param[param_finish].val[0] = end;
    psr->param[param_finish].paramSet[0] = 1;
    psr->param[param_finish].fitFlag[0] = 1;

    preProcess(psr,*npsr,argc,argv);
    formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
    formResiduals(psr,*npsr,1);    /* Form the residuals                 */

    t2Fit(psr,*npsr,covarFuncFile);
    formBatsAll(psr,*npsr);                /* Form Barycentric arrival times */
    formResiduals(psr,*npsr,1);       /* Form residuals */
    textOutput(psr,*npsr,globalParameter,0,0,0,""); 

    for (int ip=0; ip < nparam; ++ip){
        ovalues[ip] = psr->param[params[ip]].val[param_k[ip]];
    }
    double t=start+width/2.0;
    unsigned int iepoch=0;
    char outpar[80];
    strcpy(outpar,"");

    logmsg("%lf %lf",t,end);
    while (t < end) {

        for (param_label ip = 0; ip < param_LAST; ++ip){
            for (int k=0; k < psr[0].param[ip].aSize; ++k){
                psr->param[ip].fitFlag[k]=0;
            }
        }
        for (int ij=0; ij <= psr->nJumps; ++ij) {
            psr->fitJump[ij]=0;
        }
        for (int ip=0; ip < nparam; ++ip){
            psr->param[params[ip]].fitFlag[param_k[ip]] = 1;
        }

        psr->param[param_start].val[0] = t-width/2.0;
        psr->param[param_start].paramSet[0] = 1;
        psr->param[param_start].fitFlag[0] = 1;
        psr->param[param_finish].val[0] = t+width/2.0;
        psr->param[param_finish].paramSet[0] = 1;
        psr->param[param_finish].fitFlag[0] = 1;


        updateEpoch(psr,0,t);

        if (delta_correct_epoch){
        for (int ip=0; ip < nparam; ++ip){
            ovalues[ip] = psr->param[params[ip]].val[param_k[ip]];
        }
        }



        formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
        formResiduals(psr,*npsr,1);    /* Form the residuals                 */

        int ndata=0;
        for (int i=0; i < psr->nobs; ++i) {
            if( psr->obsn[i].bat > psr->param[param_start].val[0] 
                    && psr->obsn[i].bat < psr->param[param_finish].val[0]
                    ) ndata++;
        }

        if(ndata < (nparam+2)){
            t += dt;
            continue;
        }


        starts[iepoch] = t-width/2.0;
        ends[iepoch] = t+width/2.0;

        t2Fit(psr,*npsr,covarFuncFile);
        formBatsAll(psr,*npsr);                /* Form Barycentric arrival times */
        formResiduals(psr,*npsr,1);       /* Form residuals */
        if(pars){
            sprintf(outpar,"stride.%lf",t);
        }
        textOutput(psr,*npsr,globalParameter,0,0,pars,outpar);

        for (int ip=0; ip < nparam; ++ip){
            values[iepoch][ip] = psr->param[params[ip]].val[param_k[ip]];
            if (delta){
                values[iepoch][ip] -= ovalues[ip];
            }
            errors[iepoch][ip] = psr->param[params[ip]].err[param_k[ip]];
        }

        // reset the parameters to prefit values
        for (param_label ip = 0; ip < param_LAST; ++ip){
            for (int k=0; k < psr[0].param[ip].aSize; ++k){
                psr->param[ip].val[k] = psr->param[ip].prefit[k];
            }
        }
        t += dt;
        ++iepoch;
    }

    unsigned nepoch=iepoch;

    t=start+width/2.0;
    FILE* out = fopen("stridefit.out","w");

    fprintf(out,"#MJD    \t");
    for (int ip=0; ip < nparam; ++ip){
        fprintf(out,"%s\t(err)\t",psr->param[params[ip]].shortlabel[param_k[ip]]);
    }

    fprintf(out,"\n");
    for(iepoch =0; iepoch < nepoch; ++iepoch){
        fprintf(out,"%.1f",t);
        for (int ip=0; ip < nparam; ++ip){
            fprintf(out,"\t%.20g\t%g",values[iepoch][ip], errors[iepoch][ip]);
        }

        fprintf(out,"\t%.20g\t%g",starts[iepoch], ends[iepoch]);
        t+=dt;
        fprintf(out,"\n");
    }
    fclose(out);

    return 0;
}

const char * plugVersionCheck = TEMPO2_h_VER;
