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
}


/**
 * An accelerated cholesky decomposion to form uinv in plac.
 * uinv is a lower triangular, row-major, matrix.
 */
int accel_uinv(double* _m, int n){
    int i,j;

    double* _u=_m;
// UNUSED VARIABLE //     double* _uinv=_m;


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

/**
 * Do the least squares using QR decomposition
 */
double accel_lsq_qr(double** A, double* data, double* oparam, int ndata, int nparam, double** Ocvm){
    int nhrs=1;
    int nwork = -1;
    int info=0;
    int i,j;
    double iwork;

    // workspace query - works out optimal size for work array
    F77_dgels("T", &nparam, &ndata, &nhrs, A[0], &nparam, data, &ndata, &iwork, &nwork, &info);
    nwork=(int)iwork;
    logdbg("nwork = %d (%lf)",nwork,iwork);
    double* work = static_cast<double*>(malloc(sizeof(double)*nwork));
    logdbg("accel_lsq_qr ndata=%d nparam=%d",ndata,nparam);

    // as usual we prentend that the matrix is transposed to deal with C->fortran conversion
    // degls does a least-squares using QR decomposition. Fast and robust.
    F77_dgels("T", &nparam, &ndata, &nhrs, A[0], &nparam, data, &ndata, work, &nwork, &info);

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
        memcpy(oparam,data,sizeof(double)*nparam);
    }
    double chisq=0;
    for( i = nparam; i < ndata; i++ ) chisq += data[i] * data[i];
    if (Ocvm != NULL){

        int n=nparam;


        // packed triangular matrix.
        double* _t=(double*)malloc(sizeof(double)*(n*(n+1))/2);

        // This code taken from the LAPACK documentation
        // to pack a triangular matrix.

        // we want the upper triangular matrix part of A.
        //
        // pack upper triangle like (r,c)
        // (1,1) (1,2) (2,2)
        //  i+(2n-j)(j-1)/2
        int jc=0;
        for (j=0;j<n;j++){ // cols
            for (i=0; i <=j; i++) { // rows
                _t[jc] = A[i][j]; // A came from fortran, so is in [col][row] ordering
                // BUT - we have transposed A, so we have to un-transpose it
                ++jc;
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

        double **Rinv = malloc_uinv(n);

        for (j=0;j<n;j++){ // cols
            for (i=0;i<n;i++){ //rows
                Ocvm[i][j]=0;
            }
        }
        // Unpack the triangular matrix using reverse of above
        // We will put it in fortran, so continue to use [col][row] order
        jc=0;
        for (j=0;j<n;j++){ // cols
            for (i=0; i <=j; i++) { // rows
                Rinv[j][i] = _t[jc];
                Ocvm[j][i] = _t[jc];
                ++jc;
            }
        }

        free(_t);

        double a=chisq/(double)(ndata-nparam);
        // (X^T X)^-1 = Rinv.Rinv^T gives parameter covariance matrix
        // Note that Ocvm is input and output
        // and that covar matrix will be transposed, but it is
        // symetric so it doesn't matter!
        F77_dtrmm(  "R",  "U",    "T",  "N", &n, &n,    &a, *Rinv,  &n, *Ocvm, &n);
        // DTRMM ( SIDE, UPLO, TRANSA, DIAG,  M,  N, ALPHA,  A   , LDA,     B, LDB )

        if(debugFlag){
            for(i=0;i<n;i++){
                for(j=0;j<n;j++){
                    logdbg("COVAR %d %d %lg",i,j,Ocvm[i][j]);
                }
            }
        }

        free_uinv(Rinv);
    }


    return chisq;
}


#endif



#ifdef HAVE_BLAS


void accel_multMatrixVec(double* m1,double* v, int ndata,int npol, double* out){

    int m,n;
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
