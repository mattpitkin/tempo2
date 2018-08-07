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
#include <sstream>
#include <string>
#include <assert.h>
#include "../../tempo2.h"
#include "t2fit.h"
#include "mjkbayes.h"
#include "multinest.h"
#include "../../constraints_nestlike.h"
#include "../../t2fit_nestlike.h"
#include "../../TKmatrix.h"
#include "../../T2toolkit.h"
#include "../../constraints.h"

#include <vector>
#include <algorithm>
#include <mpi.h>


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
    bool actuallydoit = true;
    *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

    char quiet=1;


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
        if (strcmp(argv[i],"-v")==0) {
            quiet=0;
        }

        if (strcmp(argv[i],"-anaonly")==0) {
            actuallydoit=false;
        }

    }

    readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
    readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
    preProcess(psr,*npsr,argc,argv);

    formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
    formResiduals(psr,*npsr,1);    /* Form the residuals                 */

    //psr[0].dmoffs_fills_TN = 0;

    /// Slightly strange, but lets make sure every parameter has a sensible prefit value
    for (int iparam =0; iparam < MAX_PARAMS; ++iparam) {
        for (int k=0; k < psr->param[iparam].aSize; ++k) {
            psr->param[iparam].prefit[k] = psr->param[iparam].val[k];
        }
    }


    int myid=0;
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    mjkcontext* data = static_cast<mjkcontext*>(calloc(sizeof(mjkcontext),1));

    loadmjkbayescfg(cfg,psr,data);


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

    int ndims = data->params.size();					// dimensionality (no. of free parameters)

    int nPar = ndims+data->xtra.size();					// total no. of parameters including free & derived parameters

    int nClsPar = ndims;				// no. of parameters to do mode separation on

    int updInt = 250;				// after how many iterations feedback is required & the output files should be updated
    // note: posterior files are updated & dumper routine is called after every updInt*10 iterations

    double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored

    int maxModes = 100;				// expected max no. of modes (used only for memory allocation)

    int pWrap[ndims];				// which parameters to have periodic boundary conditions?
    for(int i = 0; i < ndims; i++) pWrap[i] = 0;


    int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock

    int fb = 1;					// need feedback on standard output?

    int resume = 0;					// resume from a previous job?

    int outfile = 1;				// write output files?

    int initMPI = 0;				// initialize MPI routines?, relevant only if compiling with MPI
    // set it to F if you want your main program to handle MPI initialization

    double logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest

    int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
    // has done max no. of iterations or convergence criterion (defined through tol) has been satisfied

    void *context = reinterpret_cast<void*>(data);				// not required by MultiNest, any additional information user wants to pass


    strcpy(data->root,"mjkres/mm-");

    // shut tempo2 up whilst we do the sampling!
    quietFlag=quiet;
    // calling MultiNest

    if (actuallydoit) {
    nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, data->root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LogLike, dumper, context);


    }

    quietFlag=0;
    logmsg("End %d",myid);
    if(myid==0){
        mjkbayes_analyse(psr, data);
    }

    quietFlag=0;
    fclose(data->debugfile);
    free(data);

    MPI_Finalize();

    return 0;
}



double computeLogLike(double *Cube, mjkcontext* data, const char* outpar);
double computeLogLike(double *Cube, mjkcontext* data){
    return computeLogLike(Cube,data,0);
}

void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
    mjkcontext* data = reinterpret_cast<mjkcontext*>(context);
    assert(data->params.size() == ndim);
    for (int i=0; i < data->params.size(); ++i){
        Cube[i] = Cube[i]*2.0 - 1.0;
        Cube[i] *= data->params[i].fitscale;
        Cube[i] += data->params[i].fitoffset;
        //printf("%lg %lg %lg %s\n",Cube[i],data->params[i].fitscale,data->params[i].fitoffset,data->params[i].shortlabel(data->psr).c_str());
    }
    lnew = computeLogLike(Cube, data);
}
void mass2dd(double am,double am2,double x,double ecc,double an,double *arr,double *ar,
        double *xk,double *si,double *gamma,double *pbdot);


//#define MJKBAYESDEBUG

double computeLogLike(double *Cube, mjkcontext* data, const char* outpar){
    pulsar* psr = data->psr;



    bool haveModel=false;
    bool remakeResiduals=false;
    bool errorsChanged=false;
    bool binarymodel=false;
    double asini=0,m2=0,ecc=0,pb=0,incl=0;
    double m1 = 0;
    double prior=0;

    double chol_alpha, chol_amp, chol_fc;

    double *prefitjump;
    std::stringstream model;

    if (psr->nJumps > 0){
        prefitjump = new double[psr->nJumps+1];
        std::copy(psr->jumpVal,psr->jumpVal+psr->nJumps+1,prefitjump);
    }



    for (int i=0; i < data->params.size(); ++i){
        const mjkparam* p = &(data->params[i]);

        if (p->fittype == FITTYPE_PARAM) {
            remakeResiduals=true;
            if(p->exp){
                psr->param[p->fitlabel].val[p->fitk] = pow(10.0,Cube[i]);
            } else {
                psr->param[p->fitlabel].val[p->fitk] = Cube[i];
            }
        } else if (p->fittype == FITTYPE_BIN){
            binarymodel=true;
            remakeResiduals=true;
            switch(p->fitk){
                case FITTYPE_BIN_K_ASINI:
                    asini=Cube[i];
                    break;
                case FITTYPE_BIN_K_M2:
                    m2=Cube[i];
                    break;
                case FITTYPE_BIN_K_ECC:
                    ecc=Cube[i];
                    break;
                case FITTYPE_BIN_K_PB:
                    pb=Cube[i];
                    break;
                case FITTYPE_BIN_K_INC:
                    incl=(M_PI*Cube[i]/180.);
                    break;
                case FITTYPE_BIN_K_M1:
                    logerr("Can't fit for M1 right now, use PB/ASINI instead");
                    exit(1);
                    break;
            }
        } else if (p->fittype == FITTYPE_CVM){
            haveModel = true;
            model << "MODEL " << p->txt;
            for (int k=0; k < p->fitk; ++k){
                p = &(data->params[i]);
                if (p->exp){
                    model << " " << pow(10.0,Cube[i]);
                } else {
                    model << " " << Cube[i];
                }
                ++i;
            }
            --i;
            model << "\n";
        } else if (p->fittype == FITTYPE_CONCVM){
            std::stringstream cmodel;
            cmodel << "MODEL " << p->txt;
            for (int k=0; k < p->fitk; ++k){
                p = &(data->params[i]);
                if (p->exp){
                    cmodel << " " << pow(10.0,Cube[i]);
                } else {
                    cmodel << " " << Cube[i];
                }
                ++i;
            }
            --i;
            cmodel << "\n";
            //printf("%s %p\n",cmodel.str().c_str(),psr->constraint_special[(int)p->fitlabel]);
            strncpy(psr->constraint_special[(int)p->fitlabel],cmodel.str().c_str(),1024);

        } else if(p->fittype == FITTYPE_EFAC){
            errorsChanged=true;
            for (int iobs=0; iobs < psr->nobs; ++iobs){
                if(p->flagmask[iobs]){
                    psr->obsn[iobs].toaErr = Cube[i] * psr->obsn[iobs].toaErr;
                }
            }
        }
    }


    // second round to get EQUAD that has to come after EFAC
    for (int i=0; i < data->params.size(); ++i){
        const mjkparam* p = &(data->params[i]);
        if(p->fittype == FITTYPE_EQUAD){
            errorsChanged=true;
            for (int iobs=0; iobs < psr->nobs; ++iobs){
                if(p->flagmask[iobs]){
                    psr->obsn[iobs].toaErr = (sqrt(pow(psr->obsn[iobs].toaErr,2)+pow(Cube[i],2)));
                }
            }
        }
    }


    if (binarymodel) {
        double a = asini/ sin(incl);
        double massfunc = 0.618 * a*a*a/pb/pb/24.0/24.0;

        m1 = sqrt(m2*m2*m2 / massfunc) - m2;

        double mtot= m1+m2;

        /*
        double omdot = 39.73 * pow(pb*24.0,-5./3.) * pow(mtot,2./3.) / (1-ecc*ecc);
        double To = 4.925490947e-6; // Pb in seconds, masses in solar masses
        double f_e = (1.0 + (73./24.)*ecc*ecc + (37./96.)*ecc*ecc*ecc*ecc)/pow(1-ecc*ecc,3.5);
//        double xdot = -(64.0/5.0)*SPEED_LIGHT*To*To*pow(2*M_PI/(pb*86400.0),2)*f_e*m1*m2*m2/mtot;
        double pbdot = -(192*M_PI/5.0)*pow(To,5.0/3.0)*pow(pb*86400.0/M_PI/2.0,-5.0/3.0)*f_e*m1*m2/pow(mtot,1.0/3.0);
        double gamma = ecc * pow((pb*86400.0)/2.0/M_PI,1.0/3.0) * pow(To,2.0/3.0) * pow(mtot,-1.0/3.0) * m2 * (1. + (m2/mtot));
*/


//void mass2dd(double am,double am2,double x,double ecc,double an,double *arr,double *ar,
//        double *xk,double *si,double *gamma,double *pbdot)
//
double arr,ar,xk,si,gamma,pbdot;
 mass2dd(m1+m2,m2,asini, ecc,2.*M_PI/(pb*SECDAY),&arr,&ar,&xk,&si,&gamma,&pbdot);

double omdot=360.0*365.25*xk/(pb);
incl=asin(si);



        prior += log(sin(incl));

        psr->param[param_pb].val[0] = pb;
        psr->param[param_ecc].val[0] = ecc;
        psr->param[param_a1].val[0] = asini;

        psr->param[param_m2].val[0] = m2;

        psr->param[param_omdot].val[0] = omdot;
        //psr->param[param_sini].val[0] = sin(incl);
        psr->param[param_sini].val[0] = si;
//        psr->param[param_a1dot].val[0] = 0;//xdot;
        psr->param[param_gamma].val[0] = gamma;;
        psr->param[param_pbdot].val[0] = pbdot;

#ifdef MJKBAYESDEBUG
           printf("PB     %lg\n",pb);
           printf("A:     %lg\n",a);
           printf("mf:    %lg\n",massfunc);
           printf("mfs3:  %lg\n",massfunc*sin(incl)*sin(incl)*sin(incl));
           printf("m1:    %lg\n",m1);
           printf("m2:    %lg\n",m2);
           printf("ecc:   %lg\n",ecc);
           printf("omdot: %lg\n",omdot);
           printf("pbdot: %lg\n",pbdot);
           printf("gamma: %lg\n",gamma);
           //printf("xdot:  %lg\n",xdot);
           printf("asini: %lg\n",asini);
           printf("incl:  %lg\n",180.0*incl/M_PI);

           printf("arr:   %lg\n",arr);
           printf("ar:    %lg\n",ar);
           printf("xk:    %lg\n",xk);
           printf("si:    %lg\n",si);
           printf("i:     %lg\n",asin(si)*180./M_PI);
//           printf("gam:   %lg\n",gam);
//           printf("pbd:   %lg\n",pbd);
//           printf("omd:   %lg\n",360.0*365.25*xk/(pb));
#endif
    }



    if(remakeResiduals) {
        formBatsAll(psr,1);         /* Form the barycentric arrival times */
        formResiduals(psr,1,1);    /* Form the residuals                 */
    }

    if (haveModel){
        char* str = (char*) malloc(model.str().length()+1);
        strcpy(str,model.str().c_str());
        t2Fit(psr,1,str);
        free(str);
    } else {
        t2Fit(psr,1,NULL);
    }



//   printf("% 10.4lg %10.4lg\n",psr->detUinv,psr->detBinv);

    double lnew = -0.5 * (psr->fitChisq) + psr->detUinv + prior;// + psr->detBinv;


#ifdef MJKBAYESDEBUG
        textOutput(psr,1,0,0,0,1,"zz.par"); /* Output results to the screen */

        exit(1);
#endif

    for (int ixtra=0; ixtra < data->xtra.size(); ++ixtra) {
        const mjkparam* xtra = &(data->xtra[ixtra]);
        const int idx = data->params.size() + ixtra;
        if (xtra->fittype == FITTYPE_BIN){
            if (xtra->fitk == FITTYPE_BIN_K_M1){
                Cube[idx] = m1;
            }
        }
        if (xtra->fittype == FITTYPE_PARAM){
            Cube[idx] = psr->param[xtra->fitlabel].val[xtra->fitk];
        }
    }
    if (errorsChanged) {
        // copy the errors back also
        for (int iobs=0; iobs < psr->nobs; ++iobs){
            psr->obsn[iobs].toaErr = psr->obsn[iobs].origErr;
        }


    }

    if(outpar != 0) {
        const int orig_nEF=psr->nTNEF;
        const int orig_nEQ=psr->nTNEQ;
        // disable the temporary used TNred parameters

        // temporarily add the white parameters
        for (int i=0; i < data->params.size(); ++i){
            const mjkparam* p = &(data->params[i]);
            if(p->fittype == FITTYPE_EFAC){
                strcpy(psr->TNEFFlagID[psr->nTNEF],p->flagid.c_str());
                strcpy(psr->TNEFFlagVal[psr->nTNEF],p->flagval.c_str());
                psr->TNEFVal[psr->nTNEF]=Cube[i];
                psr->nTNEF++;
            }
            if(p->fittype == FITTYPE_EQUAD){
                strcpy(psr->TNEQFlagID[psr->nTNEQ],p->flagid.c_str());
                strcpy(psr->TNEQFlagVal[psr->nTNEQ],p->flagval.c_str());
                psr->TNEQVal[psr->nTNEQ]=log10(Cube[i])-6.0;
                psr->nTNEQ++;
            }

        }
        textOutput(psr,1,0,0,0,1,outpar);
        psr->nTNEF = orig_nEF;
        psr->nTNEQ = orig_nEQ;

        if (haveModel){
            char modelout[1024];
            snprintf(modelout,1024,"%s.model",outpar);
            FILE* ff = fopen(modelout,"w");
            printf("Write model to '%s'\n",modelout);
            fprintf(ff,"%s",model.str().c_str());
            fclose(ff);
        }

    }
    // copy pre-fit value back over so we always start `fresh'
    for (int iparam =0; iparam < MAX_PARAMS; ++iparam) {
        for (int k=0; k < psr->param[iparam].aSize; ++k) {
            psr->param[iparam].val[k] = psr->param[iparam].prefit[k];
        }
    }
    if (psr->nJumps > 0){
        std::copy(prefitjump,prefitjump+psr->nJumps+1,psr->jumpVal);
        delete[] prefitjump;
    }
    if (psr->dmoffsDMnum > 0){
        memset(psr->dmoffsDM,0,sizeof(double)*psr->dmoffsDMnum);
    }
    if (psr->dmoffsCMnum > 0){
        memset(psr->dmoffsCM,0,sizeof(double)*psr->dmoffsCMnum);
    }




    return lnew;
}


void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &o, void *context)
{
    // convert the 2D Fortran arrays to C++ arrays


    // the posterior distribution
    // postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
    /*
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
    */
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
    bool havewhite=false;
    int k=0;
    while (!feof(infile)){
        int read = fscanf(infile,"%s",keyword);
        if (read<=0)continue;
        if (keyword[0] =='#'){
            fgets(line,1024,infile);
            continue;
        }
        logmsg("Read keyword '%s'",keyword);
        if (strncasecmp(keyword,"xtra",4) == 0){
            fscanf(infile,"%s",keyword2);
            if (strncasecmp(keyword2,"binary",6) == 0){
                fscanf(infile,"%s",keyword3);
                if(strcasecmp(keyword3,"M1")==0){
                    mjkparam p(FITTYPE_BIN,param_LAST,FITTYPE_BIN_K_M1);
                    data->xtra.push_back(p);
                    logmsg("Got: eXport 'bin_M1' no prior");
                }
            }
            if (strncasecmp(keyword2,"param",5) == 0){
                int thelab=-1;
                fscanf(infile,"%s",keyword3);
                for(int param=0; param < param_LAST; ++param){
                    for (int ik=0; ik < psr->param[param].aSize; ++ik)
                        if(strcasecmp(keyword3,psr->param[param].shortlabel[ik])==0){
                            k=ik;
                            thelab=param;
                            break;
                        }
                }
                if (thelab==-1){
                    logerr("no such parameter '%s'",keyword3);
                    exit(1);
                }
                data->xtra.emplace_back(FITTYPE_PARAM,(label)thelab,k);
                logmsg("Got: eXport '%s' no prior",psr->param[thelab].shortlabel[k]);
            }
        }
        if (strncasecmp(keyword,"fit",3) == 0){
            // toggle a fit parameter. Need to know second word
            int read = fscanf(infile,"%s",keyword2);
            if (read==0){
                logerr("Error - need a second keyword for fit keyword");
                exit(1);
            }
            if (strncasecmp(keyword2,"param",5) == 0){
                // parameter fit
                fscanf(infile,"%s",keyword3);
                int thelab=-1;
                for(int param=0; param < param_LAST; ++param){
                    for (int ik=0; ik < psr->param[param].aSize; ++ik)
                        if(strcasecmp(keyword3,psr->param[param].shortlabel[ik])==0){
                            k=ik;
                            thelab=param;
                            break;
                        }
                }
                if (thelab==-1){
                    logerr("no such parameter '%s'",keyword3);
                    exit(1);
                }

                char dat[1024];
                char* txt = fgets(dat,1024,infile);

                mjkparam p(FITTYPE_PARAM,(label)thelab,k);
                txt = p.parseScaleoffset(txt);

                // we are going to search over this parameter, so disable fitting.
                psr->param[thelab].paramSet[k]=1;
                psr->param[thelab].val[k] = p.fitoffset;
                psr->param[thelab].fitFlag[k]=0;

                // add the parameter to the params list.
                data->params.push_back(p);

                logmsg("%s",p.fitdesc(psr).c_str());
            }
            if (strncasecmp(keyword2,"cov",3) == 0){

                fscanf(infile,"%s",keyword3);
                char dat[1024];
                char* txt = fgets(dat,1024,infile);
                int ii=0;
                while(txt!=NULL){
                    mjkparam p(FITTYPE_CVM,param_LAST,ii);
                    p.txt = std::string(keyword3);
                    txt = p.parseScaleoffset(txt);
                    if (txt!=NULL){
                        data->params.push_back(p);
                        logmsg("%s",p.fitdesc(psr).c_str());
                        ++ii;
                    }
                }

                logmsg("Got %d cov parameters for %s",ii,keyword3);

                for (int i=0; i < ii; ++i){
                    (data->params.end()-i-1)->fitk=ii;
                }

            }
            if (strncasecmp(keyword2,"constrain",9) == 0){

                fscanf(infile,"%s",keyword3);
                constraint c =constraint_LAST;

                logmsg("Got constrain: %s",keyword3);

                //if (strcasecmp(keyword3,"DMMODEL_CMCOV")==0) c=constraint_dmmodel_cmcov;
                //if (strcasecmp(keyword3,"DMMODEL_DMCOV")==0) c=constraint_dmmodel_dmcov;

                if (c==constraint_LAST) {
                    logerr("Not sure what to do with '%s'",keyword3);
                    exit(1);
                }

                int ic;
                for (ic=0; ic < psr->nconstraints; ++ic){
                    if(psr->constraints[ic] == c)break;
                }
                psr->constraints[ic] = c;
                psr->constraint_special[ic] = (char*)malloc(1025);
                if(ic >= psr->nconstraints)psr->nconstraints=ic+1;

                fscanf(infile,"%s",keyword3);
                char dat[1024];
                char* txt = fgets(dat,1024,infile);
                int ii=0;
                while(txt!=NULL){
                    mjkparam p(FITTYPE_CONCVM,(label)ic,ii);
                    p.txt = std::string(keyword3);
                    txt = p.parseScaleoffset(txt);
                    if (txt!=NULL){
                        data->params.push_back(p);
                        logmsg("%s",p.fitdesc(psr).c_str());
                        ++ii;
                    }
                }
                logmsg("Got %d constrain cov parameters for %s (nconstraints=%d) %s",ii,keyword3,psr->nconstraints,get_constraint_name(c).c_str());
                for (int i=0; i < ii; ++i){
                    (data->params.end()-i-1)->fitk=ii;
                }

            }

            if (strncasecmp(keyword2,"binary",6) == 0){
                fscanf(infile,"%s",keyword3);

                psr->param[param_pb].paramSet[0]=1;
                psr->param[param_pb].fitFlag[0]=0;

                psr->param[param_a1].paramSet[0]=1;
                psr->param[param_a1].fitFlag[0]=0;

                psr->param[param_ecc].paramSet[0]=1;
                psr->param[param_ecc].fitFlag[0]=0;


                psr->param[param_omdot].paramSet[0]=1;
                psr->param[param_omdot].fitFlag[0]=0;

                psr->param[param_m2].paramSet[0]=1;
                psr->param[param_m2].fitFlag[0]=0;

//                psr->param[param_a1dot].paramSet[0]=1;
//                psr->param[param_a1dot].fitFlag[0]=0;

                psr->param[param_pbdot].paramSet[0]=1;
                psr->param[param_pbdot].fitFlag[0]=1;

                psr->param[param_gamma].paramSet[0]=1;
                psr->param[param_gamma].fitFlag[0]=0;

                psr->param[param_sini].paramSet[0]=1;
                psr->param[param_sini].fitFlag[0]=0;


                char dat[1024];
                char* txt = fgets(dat,1024,infile);

                mjkparam p(FITTYPE_BIN,param_LAST,0);
                txt = p.parseScaleoffset(txt);

                if(strcasecmp(keyword3,"asini")==0) {
                    p.fitk=FITTYPE_BIN_K_ASINI;
                    //logmsg("Got: Fit for binary asin(i), uniform prior over %lg to %lg",centre-halfrange,centre+halfrange);
                }else if(strcasecmp(keyword3,"m2")==0) {
                    p.fitk=FITTYPE_BIN_K_M2;
                    //logmsg("Got: Fit for binary M2, uniform prior over %lg to %lg",centre-halfrange,centre+halfrange);
                } else if(strcasecmp(keyword3,"ecc")==0) {
                    p.fitk=FITTYPE_BIN_K_ECC;
                    //logmsg("Got: Fit for binary ECC, uniform prior over %lg to %lg",centre-halfrange,centre+halfrange);
                } else if(strcasecmp(keyword3,"pb")==0) {
                    p.fitk=FITTYPE_BIN_K_PB;
                    //logmsg("Got: Fit for binary pb, uniform prior over %lg to %lg",centre-halfrange,centre+halfrange);
                } else if(strcasecmp(keyword3,"inc")==0) {
                    p.fitk = FITTYPE_BIN_K_INC;
                    //logmsg("Got: Fit for binary inclination, sin(i) prior over %lg to %lg deg",centre-halfrange,centre+halfrange);
                }
                // add the parameter to the params list.
                data->params.push_back(p);

                logmsg("%s",p.fitdesc(psr).c_str());


            }
            if (strncasecmp(keyword2,"white",5) == 0){
                // white noise
                havewhite=true;
                fscanf(infile,"%s",keyword3);
                if (strncasecmp(keyword3,"efac",8) == 0){
                    fscanf(infile,"%s %s",flag,flagval);
                    mjkparam p(FITTYPE_EFAC,param_LAST,0);

                    char dat[1024];
                    char* txt = fgets(dat,1024,infile);

                    txt = p.parseScaleoffset(txt);

                    p.flagmask = mjkbayesflagmask(psr,flag,flagval);
                    p.flagid = std::string(flag);
                    p.flagval = std::string(flagval);
                    data->params.push_back(p);
                    //logmsg("Got: Fit for EFAC (%s %s), uniform prior over %lg to %lg",flag,flagval,centre-halfrange,centre+halfrange);
                    logmsg("%s",p.fitdesc(psr).c_str());

                } else if (strncasecmp(keyword3,"equad",8) == 0){
                    fscanf(infile,"%s %s",flag,flagval);
                    mjkparam p(FITTYPE_EQUAD,param_LAST,0);
                    char dat[1024];
                    char* txt = fgets(dat,1024,infile);

                    txt = p.parseScaleoffset(txt);

                    p.flagmask = mjkbayesflagmask(psr,flag,flagval);
                    p.flagid = std::string(flag);
                    p.flagval = std::string(flagval);
                    data->params.push_back(p);
                    //logmsg("Got: Fit for EQUAD (%s %s), uniform prior over %lg to %lg",flag,flagval,centre-halfrange,centre+halfrange);
                    logmsg("%s",p.fitdesc(psr).c_str());
                }
            }
        }
    }

    if(havewhite){
        for (int iobs=0; iobs < psr->nobs; ++iobs){
            psr->obsn[iobs].origErr = psr->obsn[iobs].toaErr; // we cheat here!
            psr->nTNEF = 0;
            psr->nTNEQ = 0;

            psr->nT2efac = 0;
            psr->nT2equad = 0;
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


void mjkbayes_analyse(pulsar* psr, struct mjkcontext *context){
    logmsg("Analyse Results");
    char fname[1024];
    snprintf(fname,1024,"%spost_equal_weights.dat",context->root);
    FILE* postequal_f = fopen(fname,"r");

    snprintf(fname,1024,"%s.paramnames",context->root);
    FILE* labels_f = fopen(fname,"w");

    snprintf(fname,1024,"%s.results",context->root);
    FILE* results_f = fopen(fname,"w");
    std::vector<double> loglike;
    std::vector<std::vector<double>> params(context->params.size()+context->xtra.size());
    while (!feof(postequal_f)) {
        double val;
        for (int i=0; i < context->params.size()+context->xtra.size(); ++i){
            int n = fscanf(postequal_f,"%lg",&val);
            if (n != 1) break;
            params[i].push_back(val);
        }

        int n = fscanf(postequal_f,"%lg",&val);
        if (n != 1) break;
        loglike.push_back(val);
    }


    int imax = 0;
    for (int isamp=0; isamp < loglike.size(); ++isamp) {
        if (loglike[isamp] > loglike[imax]) {
            imax = isamp;
        }
    }
    double maxlike= loglike[imax];

    std::vector<double> maxcube;
    std::vector<double> meancube;
    std::vector<double> sigmacube;


    for (int iparam=0; iparam < context->params.size(); ++iparam) {
        double mean = TKmean_d(&params[iparam][0],params[iparam].size());
        double max  = params[iparam][imax];
        double sigma = sqrt(TKvariance_d(&params[iparam][0],params[iparam].size()));
        maxcube.push_back(max);
        meancube.push_back(mean);
        sigmacube.push_back(sigma);
    }
    for (int ixtra=0; ixtra < context->xtra.size(); ++ixtra) {
        int iparam=context->params.size()+ixtra;
        double mean = TKmean_d(&params[iparam][0],params[iparam].size());
        double max  = params[iparam][imax];
        double sigma = sqrt(TKvariance_d(&params[iparam][0],params[iparam].size()));
        maxcube.push_back(max);
        meancube.push_back(mean);
        sigmacube.push_back(sigma);
    }


    snprintf(fname,1024,"%smaxlike.par",context->root);
    double freshlike = computeLogLike(&maxcube[0],context,fname);
    snprintf(fname,1024,"%smean.par",context->root);
    double meanlike = computeLogLike(&meancube[0],context,fname);
    logmsg("LL %lg %lg %lg",maxlike,freshlike,meanlike);

    for (int iparam=0; iparam < context->params.size(); ++iparam) {
        double mean = TKmean_d(&params[iparam][0],params[iparam].size());
        double max  = params[iparam][imax];
        double sigma = sqrt(TKvariance_d(&params[iparam][0],params[iparam].size()));
        std::string lab = context->params[iparam].shortlabel(psr);
        printf("% 16s % 9.7lg % 9.7lg % 9.7lg\n",lab.c_str(),mean,max,sigma);
        fprintf(results_f,"% 16s % 9.7lg % 9.7lg % 9.7lg\n",lab.c_str(),mean,max,sigma);
        fprintf(labels_f,"%d %s\n",iparam,lab.c_str());
    }
    for (int ixtra=0; ixtra < context->xtra.size(); ++ixtra) {
        int iparam=context->params.size()+ixtra;
        double mean = TKmean_d(&params[iparam][0],params[iparam].size());
        double max  = params[iparam][imax];
        double sigma = sqrt(TKvariance_d(&params[iparam][0],params[iparam].size()));
        std::string lab = context->xtra[ixtra].shortlabel(psr);
        printf("% 16s % 9.7lg % 9.7lg % 9.7lg\n",lab.c_str(),mean,max,sigma);
        fprintf(results_f,"% 16s % 9.7lg % 9.7lg % 9.7lg\n",lab.c_str(),mean,max,sigma);
        fprintf(labels_f,"%d %s\n",iparam,lab.c_str());
    }

    fclose(labels_f);
    fclose(results_f);

}

