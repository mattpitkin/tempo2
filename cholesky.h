#ifdef __cplusplus
extern "C" {
#endif


void cholesky_readFromCovarianceFunction(double **m, char* fname,double *resx,double *resy,double *rese,int np, int nc);
void cholesky_covarFunc2matrix(double** m, double* covarFunc, int ndays,double *resx,double *resy,double *rese,int np, int nc);
void cholesky_powerlawModel(double **m, double modelAlpha, double modelFc, double modelA,double *resx,double *resy,double *rese,int np, int nc);
  void cholesky_powerlawModel_withBeta(double **m, double modelAlpha, double beta, double modelFc, double modelA,double *resx,double *resy,double *rese,int np, int nc);
int cholesky_formUinv(double **uinv,double** m,int np);
void cholesky_dmModel(double **m, double D, double d, double ref_freq,double *resx,double *resy,double *rese,int np, int nc);
  void cholesky_ecm(double **m, char* fileName,double *resx,double *resy,double *rese,int np, int nc);
  void cholesky_dmModelCovarParam(double **m, double alpha, double a, double b,double *resx,double *resy,double *rese,int np, int nc);

#ifdef __Tempo2_h
void getCholeskyDiagonals(double **uinv, pulsar *psr, double *resx,double *resy,double *rese, int np, int nc,int* ip);
void cholesky_readT2Model1(double **m, FILE* file,double *resx,double *resy,double *rese,int np, int nc,int *ip, pulsar *psr);
void cholesky_readT2Model2(double **m, FILE* file,double *resx,double *resy,double *rese,int np, int nc,int *ip, pulsar *psr);
void cholesky_readT2CholModel(double **m, char* fname,double *resx,double *resy,double *rese,int np, int nc,int *ip, pulsar *psr);
void addCovar(double **m,double **mm,double *resx,double *resy,double *rese,int np, int nc,int *ip, pulsar *psr,double mjd_start,double mjd_end);
#endif

#ifdef __cplusplus
}
#endif


