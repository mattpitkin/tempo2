#include "tempo2.h"


// Functions borrowed from spectralModel plugin
void cubicFit(double *resx,double *resy,double *rese,int nres,double *cubicVal,double *cubicErr);

void findSmoothCurve(double *resx,double *resy,double *rese,int nres,double *cubicVal,double *smoothModel,double expSmooth);
  void interpolate(double *resx,double *resy,double *rese,
		 int nres,double *cubicVal,double *interpX,
		   double *interpY,int *nInterp,int interpTime,double expSmooth);
void getHighFreqRes(double *resy,double *smoothModel,int nres,double *highFreqRes);
int calculateSpectra(double *x,double *y,double *e,int n,int useErr,int preWhite,
		     int specType,double *specX,double *specY);

//int fitSpectra(double *preWhiteSpecX,double *preWhiteSpecY,int nPreWhiteSpec,double *modelAlpha,double *modelFc,int *modelNfit,double *modelScale,double *fitVar,int aval,int ipw,double ifc,double iexp,int inpt);

int fitSpectra(double *preWhiteSpecX,double *preWhiteSpecY,int nPreWhiteSpec,double *modelAlpha,double *modelFc,int *modelNfit,double *modelScale,double *fitVar,int aval,int ipw,double ifc,double iexp,int inpt,double amp);

void calculateCholesky(double modelAlpha,double modelFc,double modelScale,double fitVar,double **uinv,double *covFunc, double *resx,double *resy,double *rese,int np,double *highFreqRes,double *errorScaleFactor, int dcmflag);

void getWhiteRes(double *resx,double *resy,double *rese,int nres,double **uinv,double *cholWhiteY);
void calculateDailyCovariance(double *x,double *y,double *e,int n,double *cv,int *in,double *zl,int usewt);
void formCholeskyMatrix_pl(double *c,double *resx,double *resy,double *rese,int np,double **uinv);
int obtainTimingResiduals(pulsar *psr,double *resx,double *resy,double *rese);


int guess_vals(double *x, double *y, int n, double *alpha, double *amp,   double *fc,  int *nfit, double wn, double *fc_white, int prewhite);
void getWhiteNoiseLevel(int n, double *y, int nlast, double *av);

// some global variables that Ryan is still using for diagnostic purposes
double FCALPHA, WNLEVEL, EXPSMOOTH, UPW, NFIT, FCFINAL;
