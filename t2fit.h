#include <tempo2.h>
#include <TKfit.h>

void t2Fit(pulsar *psr,unsigned int npsr, char *covarFuncFile);

unsigned int t2Fit_getFitData(const pulsar *psr, double* &x, double* &y,
        double* &e, int* &ip);

void t2Fit_fillGlobalFitInfo(pulsar* psr, unsigned int npsr,FitInfo &OUT);

void t2Fit_fillFitInfo(pulsar* psr, FitInfo &OUT);

void t2Fit_buildDesignMatrix(double x,double afunc[],int ma,pulsar *psr,int ipos,int ipsr);

void t2Fit_updateParameters(pulsar *psr,int ipsr,double *val,double *error);
