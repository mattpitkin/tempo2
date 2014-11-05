#include "tempo2.h"


void T2writeCovarFuncModel(double alpha,double fc,double val,double white,char *fname);
void T2get_covFunc_automatic(pulsar *psr, double expSmooth, char *outname, double *fc_w,double *fc_r, double *modelAlpha_out, double *modelVal, double *whiteNoiseLevel, int realflag, int dcmflag);
// Functions borrowed from spectralModel plugin
void T2cubicFit(double *resx,double *resy,double *rese,int nres,double *cubicVal,double *cubicErr);

void T2findSmoothCurve(double *resx,double *resy,double *rese,int nres,double *cubicVal,double *smoothModel,double expSmooth);
void T2interpolate(double *resx,double *resy,double *rese,
		 int nres,double *cubicVal,double *interpX,
		   double *interpY,int *nInterp,int interpTime,double expSmooth);
void T2getHighFreqRes(double *resy,double *smoothModel,int nres,double *highFreqRes);
int T2calculateSpectra(double *x,double *y,double *e,int n,int useErr,int preWhite,
		     int specType,double *specX,double *specY);

//int fitSpectra(double *preWhiteSpecX,double *preWhiteSpecY,int nPreWhiteSpec,double *modelAlpha,double *modelFc,int *modelNfit,double *modelScale,double *fitVar,int aval,int ipw,double ifc,double iexp,int inpt);

int T2fitSpectra(double *preWhiteSpecX,double *preWhiteSpecY,int nPreWhiteSpec,double *modelAlpha,double *modelFc,int *modelNfit,double *modelScale,double *fitVar,int aval,int ipw,double ifc,double iexp,int inpt,double amp, int useBeta,double *betaVal,double cutoff=0.0);

void T2calculateCholesky(double modelAlpha,double modelFc,double modelScale,double fitVar,double **uinv,double *covFunc, double *resx,double *resy,double *rese,int np,double *highFreqRes,double *errorScaleFactor, int dcmflag,int useBeta,double betaVal);
int T2calculateCovarFunc(double modelAlpha,double modelFc,double modelA,int useBeta,double betaVal,double *covFunc,double *resx,double *resy,double *rese,int np);

void T2getWhiteRes(double *resx,double *resy,double *rese,int nres,double **uinv,double *cholWhiteY);
void T2calculateDailyCovariance(double *x,double *y,double *e,int n,double *cv,int *in,double *zl,int usewt);
int T2obtainTimingResiduals(pulsar *psr,double *resx,double *resy,double *rese);


int T2guess_vals(double *x, double *y, int n, double *alpha, double *amp,   double *fc,  int *nfit, double wn, double *fc_white, int prewhite);
void T2getWhiteNoiseLevel(int n, double *y, int nlast, double *av);
void T2cholDecomposition(double **a, int n, double *p);

// some global variables that Ryan is still using for diagnostic purposes
extern double FCALPHA, WNLEVEL, EXPSMOOTH, UPW, NFIT, FCFINAL;
