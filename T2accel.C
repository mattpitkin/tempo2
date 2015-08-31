#include "T2accel.h"
#include "tempo2.h"
#include "TKmatrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

char useT2accel=1;

#ifdef HAVE_LAPACK
#define F77_dpotf2 F77_FUNC (dpotf2, DPOTF2)
#define F77_dtptri F77_FUNC (dtptri, DTPTRI)
#define F77_dgels F77_FUNC (dgels, DGELS)
#define F77_dtrmm F77_FUNC (dtrmm, DTRMM)
#define F77_dgeequ F77_FUNC (dgeequb, DGEEQUB)
#define F77_dgesvj F77_FUNC (dgesvj, DGESVJ)


#define F77_dgemm F77_FUNC(dgemm,DGEMM)
extern "C" {
    extern void F77_dgemm(const char* ta, const char* tb, int* m, int* n, int* k, double* alpha, 
            double* a, int* lda, double* b, int* ldb, double* beta, double* c, int* ldc);

#define F77_dgemv F77_FUNC(dgemv,DGEMV)
    extern void F77_dgemv(const char* trans, int* m, int* n, double* alpha, 
            double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);
}


extern "C" {
    extern void F77_dpotf2(const char* uplo, int* n, double* a, int* lda, int* info);
    extern void F77_dtptri(const char* uplo,const char* diag, int* n, double* a, int* info);
    extern void F77_dgels(const char *trans, int *m, int *n, int *nhrs, double* A, int *lda, double* B, int *ldb, double* work, int *lwork, int *info);
    extern void F77_dtrmm(const char* lr,const char* uplo, const char* tr, const char* diag, int* n, int*m, double* alp, double* a, int* lda, double* b, int* ldb);
    extern void F77_dgeequ(int* m, int* n, double* A, int* lda, double* R, double* C, double* rowcnd, double* colcnd, double* amax, int* info);
    extern void F77_dgesvj(const char* joba, const char* jobu, const char* jobv, int* m, int* n, double* A, int* lda, double* SVA, int* mv, double* V, int* ldv, double* work, int* lwork, int* info);
}


/**
 * An accelerated cholesky decomposion to form uinv in plac.
 * uinv is a lower triangular, row-major, matrix.
 */
int accel_uinv(double* _m, int n){
    int i,j;

    double* _u=_m;
    double* _uinv=_m;


    // LAPACK Cholesky factorisation.
    // _u is a symetric matrix.
    F77_dpotf2("L",&n,_u,&n,&i);
    if(i!=0){
        logerr("Error in Cholesky Decomp i=%d",i);
        return i;
    }

    double* _t=(double*)malloc(sizeof(double)*(n*(n+1))/2);
    // This code taken from the LAPACK documentation
    // to pack a triangular matrix.
    int jc=0;
    for (j=0;j<n;j++){ // cols
        for (i=j;i<n;i++){ //rows
            // note that at this point _u is ordered column major
            // because the LAPACK routine is a FORTRAN code.
            _t[jc+i-j] = _u[j*n+i];
        }
        jc=jc+n-j;
    }

    logdbg("Done CholDecomp... Inverting...");
    F77_dtptri("L","N",&n,_t,&i);
    if(i!=0){
        logerr("Error in Invert i=%d",i);
        return i;
    }

    // Unpack the triangular matrix using reverse
    // of code above, but unpacking into a row-major matrix
    // for C compatibility.
    jc=0;
    for (j=0;j<n;j++){ // cols
        for (i=0; i < j; i++){
            // when unpacking we need to zero out the strict upper triangle.
            _u[i*n+j]=0;
        }
        for (i=j;i<n;i++){ //rows
            // here we arange _u in row-major order
            // to be C compatible on return.
            _u[i*n+j]=_t[jc+i-j];
        }
        jc=jc+n-j;
    }
    free(_t);

    logdbg("Done Invert.");
    return 0;
}

double accel_lsq_svd(double** DM, double* data, double* oparam, int ndata, int nparam, double** Ocvm){
    int dmy=nparam;
    int nwork = -1;
    int info=0;
    int i,j;
    logmsg("%p %p %p %p",DM,data,oparam,Ocvm);
    assert(ndata > 0);
    assert(ndata > nparam);
    assert(DM!=NULL);
    assert(data!=NULL);
    assert(getRows_TKmatrix_d(DM)==(size_t)ndata);
    assert(getCols_TKmatrix_d(DM)==(size_t)nparam);
    if(Ocvm!=NULL){
        logmsg("%p %p",DM[0],Ocvm[0]);
        assert(getRows_TKmatrix_d(Ocvm)==(size_t)nparam);
        assert(getCols_TKmatrix_d(Ocvm)==(size_t)nparam);
    }

    // transpose A = DM_T
    double** A = malloc_matrix_d(nparam,ndata);
    for (i=0; i < nparam; ++i){
        for(j=0; j < ndata; ++j){
            A[i][j] = DM[j][i];
        }
    }

    double rowcnd,colcnd,amax;
    double RowS[ndata];
    double ColS[nparam];
    F77_dgeequ(&ndata,&nparam,*A,&ndata, RowS,ColS,&rowcnd,&colcnd,&amax,&info);
    if(rowcnd < 0.0002){
        logwarn("Design matrix is ill-conditioned. Better to use SVD fitter.");
        logmsg("dgeequ reports: rowcnd = %lg colcnd = %lg amax = %lg",rowcnd,colcnd,amax);
    }

    for (j=0; j < nparam; ++j){
        for (i=0; i < ndata; ++i){
            A[j][i] *= ColS[j];
        }
    }


    double** V = malloc_matrix_sq_d(nparam);
    nwork=(ndata+nparam);
    if(nwork < 6)nwork=6;
    double* work = (double*)calloc(nwork,sizeof(double));
    double* SVA = (double*)calloc(nparam,sizeof(double));

    //extern void F77_dgesvj(const char* joba, const char* jobu, const char* jobv,
    //int* m, int n, double* A, int* lda, double* SVA, int* mv,
    //double* V, int* ldv, double* work, int* lwork, int* info);
    F77_dgesvj("G","U","V",&ndata,&nparam,*A,&ndata,SVA,&dmy,*V,&nparam,work,&nwork,&info);
    if (info!=0){
        logerr("SVD failed");
        return -1;
    }
    assert(info==0);

    // the number of non-zero singular values
    longdouble scale = work[0];
    int rankA = (int)round(work[1]);
    int ngood = (int)round(work[2]);
    free(work);

    logmsg("SVD, scale=%lg, nparam=%d rankA=%d, ngood=%d",(double)scale,nparam,rankA,ngood);
    if (nparam==1)ngood=1;
    assert(rankA <= nparam);
    if(ngood==0){
        logerr("SVD had no good values");
        return -1;
    }

    // UPT is U transpose, ignoring rows with zero SVs.
    double** UPT=malloc_matrix_d(rankA,ndata);
    double** VIS=malloc_matrix_d(nparam,rankA);
    double** VIST=malloc_matrix_d(rankA,nparam);
    for(j=0;j<rankA; j++){
        for(i=0; i < ndata;i++){
            UPT[j][i] = A[j][i];
            //        logmsg("UPT=%d %d %lg",j,i,UPT[j][i]);
        }
    }
    const double RCOND = 1e-27;//1e-30;

    longdouble PseudoInvSVA;
    for(i=0; i < rankA; i++){
        if(i >= ngood || (SVA[i]/SVA[0]) < RCOND){
            PseudoInvSVA=0;
            SVA[i]=0;
        }
        else {
            PseudoInvSVA = longdouble(1.0)/(longdouble)SVA[i]/scale;
        }
        //logmsg("PI=%lg",(double)PseudoInvSVA);
        for(j=0; j < nparam; j++){
            VIS[j][i] = (double)(V[i][j] * PseudoInvSVA); // undo implicit transpose due to fortran
            VIST[i][j] = V[i][j] * PseudoInvSVA;
            //  logmsg("VIS = %d %d %lg",i,j,VIS[i][j]);
        }
    }

    //oparam = VIS*UPT*y
    /*
       double vec[nparam];
       accel_multMatrixVec(*UPT,data,rankA,ndata,vec);
       for(j=0; j < rankA; j++){
       logmsg("vec=%lg",vec[i]);
       }*/

    if(oparam!=NULL){
        double** M = malloc_matrix_d(nparam,ndata);
        accel_multMatrix(*VIS,*UPT,nparam,rankA,ndata,*M);
        accel_multMatrixVec(*M,data,nparam,ndata,oparam);
        for (i=0; i < nparam; i++){
            oparam[i]*=ColS[i];
            //logmsg("P[%d] = %lg",i,oparam[i]);
        }
        free_matrix_d(M);
    }
    if(Ocvm != NULL){
        longdouble sum;
        longdouble wt[nparam];
        int k;
        for (i=0;i<nparam;i++)
        {
            if (SVA[i]!=0) wt[i] = longdouble(1.0)/(longdouble)SVA[i]/(longdouble)SVA[i];
            else wt[i] = 0.0;
        }
        for (i=0;i<nparam;i++)
        {
            for (j=0;j<=i;j++)
            {
                sum=0.0;
                for (k=0;k<nparam;k++)
                    sum+=(longdouble)V[k][i]*(longdouble)V[k][j]*wt[k];
                Ocvm[i][j] = Ocvm[j][i] = (double)sum;
            }
        }


        double** OcvmT = malloc_matrix_sq_d(nparam);
        accel_multMatrix(*VIS,*VIST,nparam,rankA,nparam,*OcvmT);
        for (i=0; i < nparam; i++){
            for (j=0; j < nparam; j++){
                Ocvm[i][j]*=ColS[i]*ColS[j];
                OcvmT[i][j]*=ColS[i]*ColS[j];
                //logmsg("OCVM %lg %lg",Ocvm[i][j],OcvmT[i][j]);
            }
        }
    }

    free_matrix_d(A);
    free_matrix_d(V);
    return ndata-nparam;
}

/**
 * Do the least squares using QR decomposition
 */
double accel_lsq_qr(double** DM, double* data, double* oparam, int ndata, int nparam, double** Ocvm){
    int nhrs=1;
    int nwork = -1;
    int info=0;
    int i,j;
    double iwork;

    // transpose A = DM_T
    double** A = malloc_matrix_d(nparam,ndata);
    for (i=0; i < nparam; ++i){
        for(j=0; j < ndata; ++j){
            A[i][j] = DM[j][i];
        }
    }

    //F77_dgequ(int* m, int* n, double* A, int* lda, double* R, double* C, double* rowcnd, double* colcnd, double* amax, int* info);
    double rowcnd,colcnd,amax;
    double RowS[ndata];
    double ColS[nparam];
    F77_dgeequ(&ndata,&nparam,*A,&ndata, RowS,ColS,&rowcnd,&colcnd,&amax,&info);
    if(rowcnd < 0.0002){
        logwarn("Design matrix is ill-conditioned. Better to use SVD fitter.");
        logmsg("dgeequ reports: rowcnd = %lg colcnd = %lg amax = %lg",rowcnd,colcnd,amax);
    }

    for (j=0; j < nparam; ++j){
        for (i=0; i < ndata; ++i){
            A[j][i] *= ColS[j];
        }
    }

    //  SUBROUTINE DGEEQU( M,   N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO )
    //

    // workspace query - works out optimal size for work array
    F77_dgels("N", &ndata, &nparam, &nhrs, *A, &ndata, data, &ndata, &iwork, &nwork, &info);
    nwork=(int)iwork;
    logdbg("nwork = %d %lg (%lf)",nwork,nwork/(double)(2*nparam),iwork);
    double* work = static_cast<double*>(malloc(sizeof(double)*nwork));
    logdbg("accel_lsq_qr ndata=%d nparam=%d",ndata,nparam);

    // as usual we prentend that the matrix is transposed to deal with C->fortran conversion
    // degls does a least-squares using QR decomposition. Fast and robust.
    F77_dgels("N", &ndata, &nparam, &nhrs, *A, &ndata, data, &ndata, work, &nwork, &info);
    double* params=data; // the params are stored in teh data vector

    free(work);
    // if info is not zero then the fit failed.
    if(info!=0){
        logerr("Error in lapack DEGLS. INFO=%d See full logs for explanation.",info);
        logmsg("");
        logmsg("From: http://www.netlib.org/lapack/explore-html/d8/dde/dgels_8f.html");
        logmsg("INFO is INTEGER");
        logmsg("  = 0:  successful exit");
        logmsg("  < 0:  if INFO = -i, the i-th argument had an illegal value");
        logmsg("  > 0:  if INFO =  i, the i-th diagonal element of the");
        logmsg("        triangular factor of A is zero, so that A does not have");
        logmsg("        full rank; the least squares solution could not be");
        logmsg("        computed.");
        logmsg("");
        if(info > 0){
            logerr("It appears that you are fitting for a 'bad' parameter - E.g A jump on a non-existant flag.");
            logmsg(" TEMPO2 will NOT attempt to deal with this!");
            logmsg("Cannot continue. Abort fit.");
            return -1;
        } else {
            logmsg("Cannot continue. Abort fit.");
            return -1;
        }
    }
    assert(info==0);

    if(oparam!=NULL){
        // copy out the output parameters, which are written into the "data" array.
        // scale by column scale factor
        for (j=0; j < nparam; ++j){
            oparam[j] = params[j]*ColS[j];
        }
    }
    double chisq=0;
    for( i = nparam; i < ndata; i++ ) chisq += data[i] * data[i];
    logdbg("chisq=%lf",chisq);
    if (Ocvm != NULL){

        int n=nparam;


        // packed triangular matrix.
        double* _t=(double*)malloc(sizeof(double)*(n*(n+1))/2);

        double **R = malloc_matrix_sq_d(n);
        // This code taken from the LAPACK documentation
        // to pack a triangular matrix.

        // we want the upper triangular matrix part of A.
        //
        // pack upper triangle like (r,c)
        // (1,1) (1,2) (2,2)
        //  i+(2n-j)(j-1)/2
        for (j=0;j<nparam;j++){ // cols
            for (i=0; i <=j; i++) { // rows
                // i=0 j=1 => row0, col1
                int fi = i+1;
                int fj = j+1;
                _t[fi + (fj-1)*fj/2 - 1] = A[j][i];
                //logmsg("%d %d %lg",i,j,_t[jc]);
            }
        }


        logdbg("Inverting...");
        F77_dtptri("U","N",&n,_t,&i);
        if(i!=0){
            logerr("Error in lapack DTPTRI. INFO=%d",i);
            logmsg("From: http://www.netlib.org/lapack/explore-html/d8/d05/dtptri_8f.html");
            logmsg("INFO is INTEGER");
            logmsg("  = 0:  successful exit");
            logmsg("  < 0:  if INFO = -i, the i-th argument had an illegal value");
            logmsg("  > 0:  if INFO = i, A(i,i) is exactly zero.  The triangular");
            logmsg("  matrix is singular and its inverse can not be computed.");
            logmsg("Cannot continue - abort fit");
            return -1;
        }

        double **Rinv = malloc_matrix_sq_d(n);
        double **RinvT= malloc_matrix_sq_d(n);

        // Unpack the triangular matrix using reverse of above
        // We will put it in fortran, so continue to use [col][row] order
        for (j=0;j<n;j++){ // cols
            for (i=0; i <=j; i++) { // rows
                int fi = i+1;
                int fj = j+1;
                Rinv[j][i] = _t[fi + (fj-1)*fj/2 - 1];
                RinvT[i][j] = Rinv[j][i];
            }
        }

        free(_t);

        double a=chisq/(double)(ndata-nparam);
        a=1;
        // (X^T X)^-1 = Rinv.Rinv^T gives parameter covariance matrix
        // Note that Ocvm is input and output
        // and that covar matrix will be transposed, but it is
        // symetric so it doesn't matter!
        F77_dtrmm(  "L",  "U",    "N",  "N", &n, &n,    &a, *Rinv,  &n, *RinvT, &n);
        // DTRMM ( SIDE, UPLO, TRANSA, DIAG,  M,  N, ALPHA,  A   , LDA,     B, LDB )

        double** cvm = RinvT;

        if(debugFlag){
            for(i=0;i<n;i++){
                logmsg("ROW");
                for(j=0;j<n;j++){
                    double corr = cvm[i][j] / sqrt(cvm[i][i]*cvm[j][j]);
                    logdbg("COVAR %02d\t%02d\t%08.4le %.2f",i,j,cvm[i][j],corr);
                }
            }
        }

        for(i=0;i<n;i++){
            for(j=0;j<n;j++){
                Ocvm[i][j]=cvm[i][j]*(ColS[i]*ColS[j]);// scale CVM
            }
        }
        free_matrix_d(Rinv);
        free_matrix_d(RinvT);
    }

    free_matrix_d(A);

    return chisq;
}


#endif



#ifdef HAVE_BLAS


void accel_multMatrixVec(double* m1,double* v, int ndata,int npol, double* out){

    int m,n,k;
    double alpha=1.0,beta=0;
    int inc=1;

    m=npol;
    n=ndata;
    /*
     * An important note about this. FORTRAN effectively transposes all matricies becaue
     * the memory is ordered differently. Therefore we specify "T" to do a transpose
     *
     *
     * M.Keith 2013.
     */
    F77_dgemv("T",&m,&n,&alpha ,m1,&m,v,&inc,&beta,out,&inc);

    /*
       m=ndata;
       n=npol;
       k=ndata;

       F77_dgemm("N","T",&m,&n,&k,&alpha,m1,&m,m2,&n,&beta,out,&m);
       */
}



void accel_multMatrix(double* m1,double* m2, int ndata,int ndata2,int npol, double* out){

    int m,n,k;
    double alpha=1.0,beta=0;

    m=npol;
    n=ndata;
    k=ndata2;
    /*
     * An important note about this. FORTRAN effectively transposes all matricies becaue
     * the memory is ordered differently.
     *
     * Therefore, to compute C=A.B we do C=B.A where there is an implicit transpose of all three
     * matricies.
     *
     * M.Keith 2013.
     */
    F77_dgemm("N","N",&m,&n,&k,&alpha,m2,&m,m1,&k,&beta,out,&m);


    /*
       m=ndata;
       n=npol;
       k=ndata;

       F77_dgemm("N","T",&m,&n,&k,&alpha,m1,&m,m2,&n,&beta,out,&m);
       */
}


#endif
