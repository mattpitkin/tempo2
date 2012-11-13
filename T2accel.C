#include "T2accel.h"
#include <tempo2.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#ifdef HAVE_LAPACK
extern "C" {
   extern void dpotf2_(char* uplo, int* n, double* a, int* lda, int* info);
   extern void dtptri_(char* uplo,char* diag, int* n, double* a, int* info);
}

void accel_uinv(double* _m, int n){
   int i,j;

   double* _u=_m;
   double* _uinv=_m;
   double* _t=(double*)malloc(sizeof(double)*(n*(n+1))/2);

   dpotf2_("L",&n,_u,&n,&i);
   if(i!=0){
	  logerr("Error in Cholesky Decomp i=%d",i);
	  exit(1);
   }

   double* pp=_t;
   for (i=0;i<n;i++){
	  for (j=0;j<n;j++){
		 if(j>i){
			_u[j*n+i]=0;
		 } else{
			*pp=_u[j*n+i];
			pp++;
		 }
	  }
   }

   logmsg("Done CholDecomp... Inverting...",i);
dtptri_("U","N",&n,_t,&i);
   if(i!=0){
	  logerr("Error in Invert i=%d",i);
	  exit(1);
   }



   pp=_t;
   for (i=0;i<n;i++){
	  for (j=0;j<n;j++){
		 if(j>i)_uinv[j*n+i]=0;
		 else {
			_uinv[j*n+i]= *pp;
			pp++;
		 }
	  }
   }
   free(_t);

   logmsg("Done Invert.",i);
}


#endif

   

#ifdef HAVE_BLAS

extern "C" {
   extern void dgemm_(char* ta, char* tb, int* m, int* n, int* k, double* alpha, 
		 double* a, int* lda, double* b, int* ldb, double* beta, double* c, int* ldc);
}

void accel_multMatrix(double* m1,double* m2, int ndata,int npol, double* out){

   int m,n,k;
   double alpha=1.0,beta=0;

  m=npol;
   n=ndata;
   k=ndata;

   dgemm_("N","T",&m,&n,&k,&alpha,m2,&m,m1,&n,&beta,out,&m);


/*
   m=ndata;
   n=npol;
   k=ndata;

   dgemm_("N","T",&m,&n,&k,&alpha,m1,&m,m2,&n,&beta,out,&m);
*/
}


#endif
