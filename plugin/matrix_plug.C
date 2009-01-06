#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"

void matinv(double **array,int norder,double *det);

extern "C" int tempoOutput(int argc,char *argv[],pulsar *psr,int npsr) 
{  
  int i,j,npol=0;
  int ip[1000],term[1000];
  //  double **error,**cpyError,**invErr,*col;
  double cvm[MAX_PARAMS][MAX_PARAMS];
  double d,err,det;
  int k;

  printf("\n\n Correlation matrix ... \n\n");

  printf("\t");

  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[0].param[i].aSize;k++)
	{
	  if (psr[0].param[i].fitFlag[k]==1)
	    {
	      if (i!=param_start && i!=param_finish)
		{
		  printf("%-9.9s\t",psr[0].param[i].shortlabel[k]);
		  ip[npol] = i;
		  term[npol] = k;
		  npol++;
		}
	    }
	}
    }

  for (i=1;i<=npol;i++)
    {
      for (j=1;j<=npol;j++)
	cvm[i-1][j-1]=psr[0].covar[i][j];
    }

  printf("\n");
  for (i=0;i<npol;i++){
      printf("%s\t",psr[0].param[ip[i]].shortlabel[term[i]]);
      for (j=0;j<npol;j++){
	  if (j<=i){
	    printf("%+.8f\t",((cvm[i][j])/sqrt((cvm[j][j])*(cvm[i][i]))));  
	  }	
	}
      printf("\n");
    }

  double **error;
  double **cpyError;
  double **invErr;
  int indx[npol];
  double col[npol];


  error = (double **)malloc(npol*sizeof(double *));
  cpyError = (double **)malloc(npol*sizeof(double *));
  invErr = (double **)malloc(npol*sizeof(double *));
  for (i=0;i<npol;i++)
    {
      error[i] = (double *)malloc(npol*sizeof(double));
      cpyError[i] = (double *)malloc(npol*sizeof(double));
      invErr[i] = (double *)malloc(npol*sizeof(double));
    }

   /* Determine the inverse matrix */
  for (i=0;i<npol;i++){
    for (j=0;j<npol;j++){
      invErr[i][j]=(cvm[i][j]/sqrt((cvm[j][j])*(cvm[i][i])));
    }
  }
  
   matinv(invErr,npol,&det);



   /*   ludcmp(cpyError,npol,indx,&d);

   for (j=1;j<=npol;j++)
     {
       for (i=1;i<=npol;i++) col[i]=0.0;
       col[j]=1.0;
       lubksb(cpyError,npol,indx,col);
       for (i=1;i<=npol;i++){
	 invErr[i][j]=col[i];
       }
       }	  */
   printf("\n"); 
   printf("gcor\t") ;
   for (i=0;i<npol;i++) 
     {
       //              err = (cvm[i][i])/sqrt((cvm[i][i])*(cvm[i][i]));
       //              printf("%+.8f\t",sqrt(fabs(1.0-1.0/(err*invErr[i][i]))));


       //       err = (cvm[i][i])/sqrt((cvm[i][i])*(cvm[i][i]));
       printf("%+.8f\t",sqrt(fabs(1.0-1.0/invErr[i][i])));
     }    
   printf("\n");
   printf("dp\t"); /* See 1991ApJ 371 739,Ryba & Taylor */
   for (i=0;i<npol;i++)
     {
       //       err = (cvm[i][i])/sqrt((cvm[i][i])*(cvm[i][i]));
       //       printf("%+8.1f\t",-log10(1-fabs(sqrt(1.0-1.0/(err*invErr[i][i])))));

       //       err = (cvm[i][i])/sqrt((cvm[i][i])*(cvm[i][i]));
              printf("%+8.1f\t",-log10(1-fabs(sqrt(1.0-1.0/(invErr[i][i])))));
     }    
   printf("\n");
  for (i=0;i<npol;i++)
    {
      free(error[i]);
      free(cpyError[i]);
      free(invErr[i]);
    }
  free(error);
  free(cpyError);
  free(invErr);

 }
 /* Inverse square matrix routine copied from the original matinv.f in tempo1 */
 void matinv(double **array,int norder,double *det)
 {
   int ik[norder],jk[norder];
   int i,j,k;
   double ss;
   double amax;
   int l;

   *det = 1.0;
   for (k=0;k<norder;k++)
     {
       amax=0.0;
     pos21:
       for (i=k;i<norder;i++)
	 {
	   for (j=k;j<norder;j++)
	     {
	       if (fabs(amax) <= fabs(array[i][j])) 
		 {
		   amax = array[i][j];
		   ik[k] = i;
		   jk[k] = j;
		 }	      
	     } // 30
	 } // 30

       // Interchange rows and columns to put amax in array(k,k)
       if (amax == 0)
	 {
	   *det = 0;
	   return;
	 }
       i = ik[k];
       if (i-k < 0) goto pos21;
       else if (i-k > 0)
	 {
	   for (j=0;j<norder;j++)
	     {
	       ss = array[k][j];
	       array[k][j] = array[i][j];
	       array[i][j] = -ss;
	     }
	 }
       j=jk[k];
       if (j-k < 0) goto pos21;
       else if (j-k > 0)
	 {
	   for (i=0;i<norder;i++)
	     {
	       ss = array[i][k];
	       array[i][k] = array[i][j];
	       array[i][j] = -ss;
	     }
	 }
       // Accumulate elements of inverse matrix
       for (i=0;i<norder;i++)
	 {
	   if (i!=k) array[i][k] = -array[i][k]/amax;
	 }
       for (i=0;i<norder;i++)
	 {
	   if (i!=k)
	     {
	       for (j=0;j<norder;j++)
		 {
		   if (j!=k)
		     array[i][j] += array[i][k]*array[k][j];
		 }
	     }
	 }
       for (j=0;j<norder;j++)
	 {
	   if (j!=k) array[k][j] = array[k][j]/amax;
	 }
       array[k][k]=1.0/amax;

       *det=(*det)*amax;
     } //  100

   // Restore normal ordering of matrix
   for (l=0;l<norder;l++)
     {
       k = norder-(l+1);
       j=ik[k];
       if (j >  k)
	 {
	   for (i=0;i<norder;i++)
	     {
	       ss = array[i][k];
	       array[i][k] = -array[i][j];
	       array[i][j] = ss;
	     }
	 }
       i = jk[k];
       if (i > k)
	 {
	   for (j=0;j<norder;j++)
	     {
	       ss = array[k][j];
	       array[k][j] = -array[i][j];
	       array[i][j] = ss;
	     }
	 }
     }
}

