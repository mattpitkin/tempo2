#include "tempo2.h"

#define FITTYPE_PARAM 0
#define FITTYPE_CHOL 1

#define FITTYPE_CHOL_K_AMP 0
#define FITTYPE_CHOL_K_ALPHA 1
#define FITTYPE_CHOL_K_FC 2


struct mjkcontext {
    int nfit;
    int fittype[1024];
    label fitlabel[1024];
    int fitk[1024];
    double fitscale[1024];
    double fitoffset[1024];
    pulsar* psr;
    FILE* debugfile;
};



void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &o, void *context);
