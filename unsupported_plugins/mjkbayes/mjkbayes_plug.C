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

/* Template for a tempo2 plugin */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include "tempo2.h"
#include "t2fit.h"
#include "mjknest.h"
#include "multinest.h"
#include "../../constraints_nestlike.h"
#include "../../t2fit_nestlike.h"
#include "../../TKmatrix.h"


using namespace std;





void help() /* Display help */
{
    /* This function should contain usage information about the plugin which should (in general) be accessed */
    /* by the user pressing 'h'                                                                              */
}

/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
    char parFile[MAX_PSR][MAX_FILELEN];
    char timFile[MAX_PSR][MAX_FILELEN];
    char cfg[MAX_FILELEN];
    strcpy(cfg,"mjkbayes.cfg");
    int i;
    *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */


    printf("Graphical Interface: mjkbayes\n");
    printf("Author:              Michael Keith\n");
    printf("Version:             v1.0\n");
    printf(" --- type 'h' for help information\n");

    /* Obtain all parameters from the command line */
    for (i=2;i<argc;i++)
    {
        if ((strcmp(argv[i],"-dcf")==0) || (strcmp(argv[i],"-chol")==0)){
            strcpy(covarFuncFile,argv[++i]);
        }

        if (strcmp(argv[i],"-f")==0)
        {
            strcpy(parFile[0],argv[++i]);
            strcpy(timFile[0],argv[++i]);
        }
        if (strcmp(argv[i],"-cfg")==0) {
            strcpy(cfg,argv[++i]);
        }

    }

    readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
    readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
    preProcess(psr,*npsr,argc,argv);

    formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
    formResiduals(psr,*npsr,1);    /* Form the residuals                 */

    if (false){
    t2Fit(psr,*npsr,covarFuncFile);

    formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
    formResiduals(psr,*npsr,1);    /* Form the residuals                 */

    textOutput(psr,*npsr,0,0,0,0,0);
    }

    mjkcontext* data = static_cast<mjkcontext*>(calloc(sizeof(mjkcontext),1));

    loadmjkbayescfg(cfg,psr,data);

    /*
    data->fittype[0]  = FITTYPE_PARAM;
    data->fitlabel[0]  = param_pmra;
    data->fitoffset[0] = 4;
    data->fitscale[0]  = 6.0;
    data->fitk[0]  = 0;

    data->fittype[1]  = FITTYPE_PARAM;
    data->fitlabel[1]  = param_pmdec;
    data->fitoffset[1] = 4;
    data->fitscale[1]  = 6;
    data->fitk[1]  = 0;


    data->fittype[2]  = FITTYPE_CHOL;
    data->fitoffset[2] = -26;
    data->fitscale[2]  = 8;
    data->fitk[2]  = FITTYPE_CHOL_K_AMP;

    data->fittype[3]  = FITTYPE_CHOL;
    data->fitoffset[3] = 4;
    data->fitscale[3]  = 2;
    data->fitk[3]  = FITTYPE_CHOL_K_ALPHA;

    data->fittype[2]  = FITTYPE_CHOL;
    data->fitoffset[2] = 0.01;
    data->fitscale[2]  = 0.04;
    data->fitk[2]  = FITTYPE_CHOL_K_FC;
*/


    data->psr = psr;
    data->debugfile = fopen("mjkres/debug","w");
    if(!data->debugfile){
        logerr("ERROR opening debug file");
    }


    int IS = 1;					// do Nested Importance Sampling?

    int mmodal = 0;					// do mode separation?

    int ceff = 0;					// run in constant efficiency mode?

    int nlive = 500;				// number of live points

    double efr = 0.2;				// set the required efficiency

    double tol = 0.5;				// tol, defines the stopping criteria

    int ndims = data->nfit;					// dimensionality (no. of free parameters)

    int nPar = ndims;					// total no. of parameters including free & derived parameters

    int nClsPar = ndims;				// no. of parameters to do mode separation on

    int updInt = 250;				// after how many iterations feedback is required & the output files should be updated
    // note: posterior files are updated & dumper routine is called after every updInt*10 iterations

    double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored

    int maxModes = 100;				// expected max no. of modes (used only for memory allocation)

    int pWrap[ndims];				// which parameters to have periodic boundary conditions?
    for(int i = 0; i < ndims; i++) pWrap[i] = 0;

    char root[100];

    int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock

    int fb = 1;					// need feedback on standard output?

    int resume = 0;					// resume from a previous job?

    int outfile = 1;				// write output files?

    int initMPI = 1;				// initialize MPI routines?, relevant only if compiling with MPI
    // set it to F if you want your main program to handle MPI initialization

    double logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest

    int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
    // has done max no. of iterations or convergence criterion (defined through tol) has been satisfied

    void *context = reinterpret_cast<void*>(data);				// not required by MultiNest, any additional information user wants to pass


    strcpy(root,"mjkres/mm");
    // shut tempo2 up whilst we do the sampling!
    quietFlag=1;

    // calling MultiNest

    nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LogLike, dumper, context);

    quietFlag=0;
    fclose(data->debugfile);
    free(data);

    return 0;
}



void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
    mjkcontext* data = reinterpret_cast<mjkcontext*>(context);
    pulsar* psr = data->psr;


    bool useTNPARAM=false;
    bool remakeResiduals=false;
    bool errorsChanged=false;


    assert((data->nfit)==ndim);

    for (int i=0; i < data->nfit; ++i){
        Cube[i] = Cube[i]*2.0 - 1.0;
        Cube[i] *= data->fitscale[i];
        Cube[i] += data->fitoffset[i];
        if (data->fittype[i] == FITTYPE_PARAM) {
            remakeResiduals=true;
            psr->param[data->fitlabel[i]].val[data->fitk[i]] = Cube[i];
        } else if (data->fittype[i] == FITTYPE_CHOL){
            useTNPARAM = true;
            switch(data->fitk[i]){
                case FITTYPE_CHOL_K_ALPHA:
                    psr->TNRedGam = Cube[i];
                    break;
                case FITTYPE_CHOL_K_AMP:
                    psr->TNRedAmp = pow(10.0,Cube[i]);
                    break;
                case FITTYPE_CHOL_K_FC:
                    psr->TNRedCorner = Cube[i];
                    break;
                default:
                    logwarn("BAD CHOL_K!!");
                    break;
            }
        } else if(data->fittype[i] == FITTYPE_EFAC){
            errorsChanged=true;
            for (int iobs=0; iobs < psr->nobs; ++iobs){
                if(data->flagmask[data->fitk[i]][iobs]){
                    psr->obsn[iobs].toaErr = Cube[i] * psr->obsn[iobs].toaErr;
                }
            }
        } else if(data->fittype[i] == FITTYPE_EQUAD){
            errorsChanged=true;
            for (int iobs=0; iobs < psr->nobs; ++iobs){
                if(data->flagmask[data->fitk[i]][iobs]){
                    psr->obsn[iobs].toaErr = (sqrt(pow(psr->obsn[iobs].toaErr,2)+pow(Cube[i],2)));
                }
            }
        }



    }

    if(remakeResiduals) {
        formBatsAll(psr,1);         /* Form the barycentric arrival times */
        formResiduals(psr,1,1);    /* Form the residuals                 */
    }

    if (useTNPARAM){
        t2Fit(psr,1,"TNPARAM");
    } else {
        t2Fit(psr,1,NULL);
    }



    lnew = -0.5 * (psr->fitChisq) + psr->detUinv;

//    fprintf(data->debugfile,"%lg %lg %lg %lg %lg\n",psr->TNRedAmp,psr->TNRedGam,psr->fitChisq,psr->detUinv,lnew);

//    textOutput(psr,1,0,0,0,1,"zz.par"); /* Output results to the screen */

//    exit(1);



    // copy pre-fit value back over so we always start `fresh'
    for (int iparam =0; iparam < MAX_PARAMS; ++iparam) {
        for (int k=0; k < psr->param[iparam].aSize; ++k) {
        psr->param[iparam].val[k] = psr->param[iparam].prefit[k];
        }
    }
    if (errorsChanged) {
    // copy the errors back also
    for (int iobs=0; iobs < psr->nobs; ++iobs){
        psr->obsn[iobs].toaErr = psr->obsn[iobs].origErr;
    }

    }
}


void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &o, void *context)
{
    // convert the 2D Fortran arrays to C++ arrays


    // the posterior distribution
    // postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns

    int i, j;

    double postdist[nSamples][nPar + 2];
    for( i = 0; i < nPar + 2; i++ )
        for( j = 0; j < nSamples; j++ )
            postdist[j][i] = posterior[0][i * nSamples + j];
    // last set of live points
    // pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column

    double pLivePts[nlive][nPar + 1];
    for( i = 0; i < nPar + 1; i++ )
        for( j = 0; j < nlive; j++ )
            pLivePts[j][i] = physLive[0][i * nlive + j];
}




void loadmjkbayescfg(const char* cfg, pulsar* psr, mjkcontext *data) {
    const unsigned maxlen=1024;
    FILE * infile = fopen(cfg,"r");
    char line[maxlen];
    char keyword[maxlen];
    char keyword2[maxlen];
    char keyword3[maxlen];
    char flag[maxlen];
    char flagval[maxlen];
    double centre=0;
    double halfrange=1;
    int k=0;
    int iflag=0;
    int ifit=0;
    while (!feof(infile)){
        int read = fscanf(infile,"%s",keyword);
        if (read<=0)continue;
        if (strncasecmp(keyword,"fit",3) == 0){
            // toggle a fit parameter. Need to know second word
            int read = fscanf(infile,"%s",keyword2);
            if (read==0){
                logerr("Error - need a second keyword for fit keyword");
                exit(1);
            }
            if (strncasecmp(keyword2,"param",5) == 0){
                // parameter fit
                fscanf(infile,"%s %d %lg %lg",keyword3, &k,&centre,&halfrange);
                int thelab=-1;
                for(int param=0; param < param_LAST; ++param){
                    if(psr->param[param].aSize > k && strcasecmp(keyword3,psr->param[param].shortlabel[k])==0){
                        thelab=param;
                        break;
                    }
                }
                if (thelab==-1){
                    logerr("no such parameter '%s'",keyword3);
                }
                data->fittype[ifit]  = FITTYPE_PARAM;
                data->fitlabel[ifit]  = (label)thelab;
                data->fitoffset[ifit] = centre;
                data->fitscale[ifit]  = halfrange;
                data->fitk[ifit]  = k;
                ++ifit;
                logmsg("Got: Fit for '%s' uniform prior over %lg to %lg",psr->param[thelab].shortlabel[k],centre-halfrange,centre+halfrange);
            }
            if (strncasecmp(keyword2,"chol",4) == 0){
                // cholesky fit
                fscanf(infile,"%s %lg %lg",keyword3, &centre,&halfrange);
                data->fittype[ifit]  = FITTYPE_CHOL;
                data->fitoffset[ifit] = centre;
                data->fitscale[ifit]  = halfrange;
                if(strcasecmp(keyword3,"amp")==0) {
                    data->fitk[ifit]  = FITTYPE_CHOL_K_AMP;
                    logmsg("Got: Fit for DCF[log(amp)], uniform prior over %lg to %lg",centre-halfrange,centre+halfrange);

                } else if(strcasecmp(keyword3,"fc")==0) {
                    data->fitk[ifit]  = FITTYPE_CHOL_K_FC;
                    logmsg("Got: Fit for DCF[fc], uniform prior over %lg to %lg",centre-halfrange,centre+halfrange);

                } else if(strcasecmp(keyword3,"alpha")==0) {
                    data->fitk[ifit]  = FITTYPE_CHOL_K_ALPHA;
                    logmsg("Got: Fit for DCF[alpha] uniform prior over %lg to %lg",centre-halfrange,centre+halfrange);

                } else {
                    logerr("invalid chol keyword '%s'",keyword3);
                    exit(1);
                }
                ++ifit;
            }
            if (strncasecmp(keyword2,"white",5) == 0){
                // white noise
                fscanf(infile,"%s",keyword3);
                if (strncasecmp(keyword3,"efac",8) == 0){
                    fscanf(infile,"%s %s %lg %lg",flag,flagval,&centre,&halfrange);
                    data->fittype[ifit]   = FITTYPE_EFAC;
                    data->fitoffset[ifit] = centre;
                    data->fitscale[ifit]  = halfrange;
                    data->fitk[ifit]      = iflag;
                    data->flagmask[iflag] = mjkbayesflagmask(psr,flag,flagval);
                    iflag++;
                    ifit++;
                    logmsg("Got: Fit for EFAC (%s %s), uniform prior over %lg to %lg",flag,flagval,centre-halfrange,centre+halfrange);
                } else if (strncasecmp(keyword3,"equad",8) == 0){
                    fscanf(infile,"%s %s %lg %lg",flag,flagval,&centre,&halfrange);
                    data->fittype[ifit]   = FITTYPE_EQUAD;
                    data->fitoffset[ifit] = centre;
                    data->fitscale[ifit]  = halfrange;
                    data->fitk[ifit]      = iflag;
                    data->flagmask[iflag] = mjkbayesflagmask(psr,flag,flagval);
                    iflag++;
                    ifit++;
                    logmsg("Got: Fit for EQUAD (%s %s), uniform prior over %lg to %lg",flag,flagval,centre-halfrange,centre+halfrange);
                }
            }
        }
    }
    data->nfit = ifit;

    if(iflag > 0){
        for (int iobs=0; iobs < psr->nobs; ++iobs){
            psr->obsn[iobs].origErr = psr->obsn[iobs].toaErr; // we cheat here!
        }
    }

}



char* mjkbayesflagmask(pulsar* psr, const char* flag, const char* flagval){
    char *mask= static_cast<char*>(calloc(sizeof(char),psr->nobs));
    int count=0;
    for(int iobs = 0; iobs < psr->nobs; ++iobs){
        for(int j = 0; j < psr->obsn[iobs].nFlags; ++j){
            if (strcmp(psr->obsn[iobs].flagID[j],flag)==0) {
                if (strcmp(psr->obsn[iobs].flagVal[j],flagval)==0) {
                    mask[iobs]=1;
                    count++;
                    break;
                }
            }
        }
    }
    logmsg("matched %d obsevations",count);
    return mask;
}

