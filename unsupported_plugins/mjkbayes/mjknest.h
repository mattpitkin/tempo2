#include "tempo2.h"

#define FITTYPE_PARAM 0
#define FITTYPE_CHOL 1
#define FITTYPE_EFAC 2
#define FITTYPE_EQUAD 3

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
    char* flagmask[1024];
    char root[1024];
};

void mjkbayes_analyse(pulsar* psr, struct mjkcontext *context);


void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &o, void *context);

void loadmjkbayescfg(const char* cfg, pulsar* psr, mjkcontext *context);

char* mjkbayesflagmask(pulsar* psr, const char* flag, const char* flagval);

