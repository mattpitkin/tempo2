void cholesky_readFromCovarianceFunction(double **m, char* fname,double *resx,double *resy,double *rese,int np, int nc);
void cholesky_readT2CholModel(double **m, char* fname,double *resx,double *resy,double *rese,int np, int nc,int *ip, pulsar *psr);
void cholesky_covarFunc2matrix(double** m, double* covarFunc, int ndays,double *resx,double *resy,double *rese,int np, int nc);
void cholesky_powerlawModel(double **m, double modelAlpha, double modelFc, double modelA,double *resx,double *resy,double *rese,int np, int nc);
void cholesky_readT2Model1(double **m, FILE* file,double *resx,double *resy,double *rese,int np, int nc,int *ip, pulsar *psr);
void cholesky_formUinv(double **uinv,double** m,int np);
void cholesky_dmModel(double **m, double D, double d, double ref_freq,double *resx,double *resy,double *rese,int np, int nc);
void getCholeskyDiagonals(double **uinv, pulsar *psr, double *resx,double *resy,double *rese, int np, int nc,int* ip);

void addCovar(double **m,double **mm,double *resx,double *resy,double *rese,int np, int nc,int *ip, pulsar *psr,double mjd_start,double mjd_end);


