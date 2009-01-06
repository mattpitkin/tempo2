void TKleastSquares_svd_psr(double *x,double *y,double *sig,int n,double *p,double *e,int nf,double **cvm, double *chisq, void (*fitFuncs)(double, double [], int,pulsar *,int),int weight,pulsar *psr,double tol, int *ip);
void TKleastSquares_svd(double *x,double *y,double *sig,int n,double *p,double *e,int nf,double **cvm, double *chisq, void (*fitFuncs)(double, double [], int),int weight);
void TKremovePoly_f(float *px,float *py,int n,int m);
void TKremovePoly_d(double *px,double *py,int n,int m);
void TKleastSquares_svd_noErr(double *x,double *y,int n,double *p,int nf, void (*fitFuncs)(double, double [], int));     
