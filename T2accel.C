#include "T2accel.h"
#include "tempo2.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

char useT2accel=1;

#ifdef HAVE_LAPACK
#define F77_dpotf2 F77_FUNC (dpotf2, DPOTF2)
#define F77_dtptri F77_FUNC (dtptri, DTPTRI)

extern "C" {
   extern void F77_dpotf2(char* uplo, int* n, double* a, int* lda, int* info);
   extern void F77_dtptri(char* uplo,char* diag, int* n, double* a, int* info);
}


/**
 * An accelerated cholesky decomposion to form uinv in plac.
 * uinv is a lower triangular, row-major, matrix.
 */
int accel_uinv(double* _m, int n){
   int i,j;

   double* _u=_m;
   double* _uinv=_m;
   double* _t=(double*)malloc(sizeof(double)*(n*(n+1))/2);


   // LAPACK Cholesky factorisation.
   // _u is a symetric matrix.
   F77_dpotf2("L",&n,_u,&n,&i);
   if(i!=0){
	  logerr("Error in Cholesky Decomp i=%d",i);
      return i;
   }

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

   logdbg("Done CholDecomp... Inverting...",i);
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

   logdbg("Done Invert.",i);
   return 0;
}


#endif



#ifdef HAVE_BLAS

#define F77_dgemm F77_FUNC(dgemm,DGEMM)
extern "C" {
   extern void F77_dgemm(char* ta, char* tb, int* m, int* n, int* k, double* alpha, 
		 double* a, int* lda, double* b, int* ldb, double* beta, double* c, int* ldc);

#define F77_dgemv F77_FUNC(dgemv,DGEMV)
   extern void F77_dgemv(char* trans, int* m, int* n, double* alpha, 
		 double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);
}


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
