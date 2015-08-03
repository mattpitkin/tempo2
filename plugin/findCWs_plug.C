//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

/*
*    This file is part of TEMPO2. 
* 
*    TEMPO2 is free software: you can redistribute it and/or modify 
*    it under the terms of the GNU General Public License as published by 
*    the Free Software Foundation, either version 3 of the License, or 
*    (at your option) any later version. 
*    TEMPO2 is distributed in the hope that it will be useful, 
*    but WITHOUT ANY WARRANTY; without even the implied warranty of 
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
*    GNU General Public License for more details. 
*    You should have received a copy of the GNU General Public License 
*    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
*/

/*
*    If you use TEMPO2 then please acknowledge it by citing 
*    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
*    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
*    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
*    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
*    timing model.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "TKspectrum.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

using namespace std;


void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  char covarFuncFile[128];
  int i, j;
  double globalParameter;
  const char *CVS_verNum = "$Revision: 1.1 $";
  FILE *fout;

  strcpy(covarFuncFile,"NULL");

  if (displayCVSversion == 1) CVSdisplayVersion((char *)"findCW.C",(char *)"plugin",CVS_verNum);

  *npsr = 0;

  printf("Graphical Interface: findCW\n");
  printf("Author:              X. Zhu, G. Hobbs\n");
  printf("CVS Version:         $Revision: 1.1 $\n");
  printf(" --- type 'h' for help information\n");

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[*npsr],argv[++i]); 
	  strcpy(timFile[*npsr],argv[++i]);
	  (*npsr)++;
	}
      else if (strcmp(argv[i],"-dcf")==0)
	strcpy(covarFuncFile,argv[++i]);
    }

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

    for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
  //  i=0;
    {
      logdbg("Calling formBatsAll");
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      logdbg("Calling formResiduals");
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      logdbg("Calling doFit");
      if (i==0) doFitAll(psr,*npsr,covarFuncFile);   /* Do the fitting     */
      else textOutput(psr,*npsr,globalParameter,0,0,0,(char *)"");  /* Display the output */
    }

  // Print A+ and Ax into a file
  fout = fopen("aplus_across.dat","w");
  for (i=0;i<psr[0].quad_ifuncN_p;i++)
    {
      fprintf(fout,"%.2lf %g %g %g %g\n",psr[0].quad_ifuncT_p[i],psr[0].quad_ifuncV_p[i],psr[0].quad_ifuncE_p[i],psr[0].quad_ifuncV_c[i],psr[0].quad_ifuncE_c[i]);
    }
  fclose(fout);
  //  return 0;
  // Calculate the Detection Statistics as a function of Fourier frequencies
  {
    double freq[1024];
    double DS[1024];
    double dt, Tspan, fmin, fmax;
    int    nSpec;    /* number of independent freq channels */
    int    nSpecOS4; /* number of freq channels with an OverSampling factor of 4 */
    const int    lp = psr[0].quad_ifuncN_p;     /* length of A+ or Ax */
    const int    lpc = 2*psr[0].quad_ifuncN_p;  /* length of the stacked data A+,x assuming A+&Ax have the same length*/
    const int    n_cst = 18; /* number of constraints set on A+&Ax */
    const int    mlen = lpc-n_cst; /* rank of the noise covariance matrix */
    // Assuming that the ifuncs are regularly sampled
    //
    dt = psr[0].quad_ifuncT_p[1]-psr[0].quad_ifuncT_p[0];
    Tspan = psr[0].quad_ifuncT_p[lp-1]-psr[0].quad_ifuncT_p[0];
    fmin = 1.0/Tspan/86400.0;
    fmax = 1.0/dt/86400.0/2;
    nSpec = floor (fmax/fmin);
    nSpecOS4 = 4*nSpec; /* somehow gsl can't handle f=fmax=36.5*fin, so f only goes up to 36.25*fmin */
    // check
    /*
    printf("%d\t%d\t%d\n", lp, lpc, mlen);
    printf("%g\t%g\t%d\n", fmin, fmax, nSpec);
    */

    gsl_vector *time = gsl_vector_alloc (lp);
    gsl_vector *Apn = gsl_vector_alloc (lp);
    gsl_vector *Acn = gsl_vector_alloc (lp);
    for (i = 0; i < lp; i++)
    {
      gsl_vector_set (time, i, psr[0].quad_ifuncT_p[i]);
      gsl_vector_set (Apn, i, psr[0].quad_ifuncV_p[i]);
      gsl_vector_set (Acn, i, psr[0].quad_ifuncV_c[i]);
    }

    gsl_matrix *Sigma_n = gsl_matrix_alloc (lpc, lpc); /* noise convariance matrix */
    gsl_matrix *Sigma_n1 = gsl_matrix_alloc (lpc, lpc); /* inverse noise convariance matrix */
    for (i = 0; i < psr[0].globalNfit; i++)
      {
        for (j = 0; j < psr[0].globalNfit; j++)
          {
            if ((psr[0].fitParamI[i] == param_quad_ifunc_p || psr[0].fitParamI[i] == param_quad_ifunc_c) && 
				(psr[0].fitParamI[j] == param_quad_ifunc_p || psr[0].fitParamI[j] == param_quad_ifunc_c))
              {
                gsl_matrix_set (Sigma_n, i, j, psr[0].covar[i][j]);
              }
          }
      }

    // print the covariance matrix to screen
    /* 
    for (j = 0; j < lpc; j++)
      printf ("Sigma_n(%d,%d) = %g\n", 0, j, gsl_matrix_get (Sigma_n, 0, j));
    */

    // eign-decomposition
    gsl_vector *eval = gsl_vector_alloc (lpc);
    gsl_matrix *evec = gsl_matrix_alloc (lpc, lpc);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (lpc);
    gsl_eigen_symmv (Sigma_n, eval, evec, w);
    gsl_eigen_symmv_free (w);
    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    gsl_matrix *D1 = gsl_matrix_alloc (mlen, mlen);
    gsl_matrix *A1 = gsl_matrix_alloc (lpc, mlen);
    gsl_matrix_set_zero (D1);

    /* print the eign-values, checked that they agree with results from Matlab
    for (j = 0; j < lpc; j++)  
      printf ("%d\t%g\n", j, gsl_vector_get (eval, j));
    */
    for (i = 0; i < mlen; i++)
      {
        double eval_i = gsl_vector_get (eval, i+n_cst);
        gsl_matrix_set (D1, i, i, 1.0/eval_i);
        gsl_vector *evec_i = gsl_vector_alloc (lpc);
        gsl_matrix_get_col (evec_i, evec, i+n_cst);
        gsl_matrix_set_col (A1, i, evec_i);
        gsl_vector_free (evec_i);
      }
    gsl_vector_free (eval);
    gsl_matrix_free (evec);
    gsl_matrix *AD1 = gsl_matrix_alloc (lpc, mlen);
    gsl_matrix_set_zero (AD1);
    gsl_matrix_set_zero (Sigma_n1);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, A1, D1, 0.0, AD1);
    gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, AD1, A1, 0.0, Sigma_n1); /* Sigma_n1=A1*D1*A1' */

    // print the inverse covariance matrix to screen
    /* checked and agree with Matlab results
    for (j = 0; j < lpc; j++)
      printf ("Sigma_n1(%d,%d) = %g\n", 1, j+1, gsl_matrix_get (Sigma_n1, 0, j));
    */
    gsl_matrix_free (Sigma_n);
    gsl_matrix_free (D1);
    gsl_matrix_free (A1);
    gsl_matrix_free (AD1);
    /* something (unknown) wrong with the following approach of getting submatrix
    gsl_matrix_view Sb1 = gsl_matrix_submatrix (Sigma_n1, 0, 0, lp, lp);
    gsl_matrix_view Sb2 = gsl_matrix_submatrix (Sigma_n1, 0, lp, lp, lp);
    gsl_matrix_view Sb3 = gsl_matrix_submatrix (Sigma_n1, lp, 0, lp, lp);
    gsl_matrix_view Sb4 = gsl_matrix_submatrix (Sigma_n1, lp, lp, lp, lp);
    gsl_matrix *S11 = &Sb1.matrix;
    gsl_matrix *S12 = &Sb2.matrix;
    gsl_matrix *S21 = &Sb3.matrix;
    gsl_matrix *S22 = &Sb4.matrix;
    */
    gsl_matrix *S11 = gsl_matrix_alloc (lp, lp);
    gsl_matrix *S12 = gsl_matrix_alloc (lp, lp);
    gsl_matrix *S21 = gsl_matrix_alloc (lp, lp);
    gsl_matrix *S22 = gsl_matrix_alloc (lp, lp);
    for (i = 0; i < lp; i++)
    {    
      for (j = 0; j < lp; j++)
      {
        gsl_matrix_set (S11, i, j, gsl_matrix_get (Sigma_n1, i, j));
        gsl_matrix_set (S12, i, j, gsl_matrix_get (Sigma_n1, i, (lp+j)));
        gsl_matrix_set (S21, i, j, gsl_matrix_get (Sigma_n1, (lp+i), j));
        gsl_matrix_set (S22, i, j, gsl_matrix_get (Sigma_n1, (lp+i), (lp+j)));
      }
    }
    gsl_matrix_free (Sigma_n1);
    /* print the submatrix
    for (j = 0; j < lp; j++)
      printf ("S11(%d,%d) = %g\n", 1, j+1, gsl_matrix_get (S11, 0, j));
    */
    // calculate the Detection Statistics for a set of frequencies
    for (i = 0; i < nSpecOS4; i++)
    {
      double fIndx = 0.5+0.25*i;
      double f = fmin*fIndx;
      double y11, y12, y21, y22, y31, y32, y41, y42;
      double r11, r12, r13, r14, r22, r23, r24, r33, r34, r44;
      double d1, d2, d3, d4, d5, d6, d7, d8;
      double y1, y2, y3, y4, r21, r31, r32, r41, r42, r43;
      gsl_vector *Sc = gsl_vector_alloc (lp);
      gsl_vector *Ss = gsl_vector_alloc (lp);
      for (j = 0; j < lp; j++)
      {
        gsl_vector_set (Sc, j, cos(2.0*86400.0*M_PI*f*psr[0].quad_ifuncT_p[j]));
        gsl_vector_set (Ss, j, sin(2.0*86400.0*M_PI*f*psr[0].quad_ifuncT_p[j]));
      }
      gsl_vector *tempy1 = gsl_vector_alloc (lp);
      gsl_vector *tempy2 = gsl_vector_alloc (lp);
      gsl_vector *tempy3 = gsl_vector_alloc (lp);
      gsl_vector *tempy4 = gsl_vector_alloc (lp);
      gsl_vector *tempy5 = gsl_vector_alloc (lp);
      gsl_vector *tempy6 = gsl_vector_alloc (lp);
      gsl_vector *tempy7 = gsl_vector_alloc (lp);
      gsl_vector *tempy8 = gsl_vector_alloc (lp);
      gsl_vector *tempr1 = gsl_vector_alloc (lp);
      gsl_vector *tempr2 = gsl_vector_alloc (lp);
      gsl_vector *tempr3 = gsl_vector_alloc (lp);
      gsl_vector *tempr4 = gsl_vector_alloc (lp);
      gsl_vector *tempr5 = gsl_vector_alloc (lp);
      gsl_vector *tempr6 = gsl_vector_alloc (lp);
      gsl_vector *tempr7 = gsl_vector_alloc (lp);
      gsl_vector *temp5 = gsl_vector_alloc (lp);
      gsl_vector *temp6 = gsl_vector_alloc (lp);
      gsl_vector *temp7 = gsl_vector_alloc (lp);
      gsl_vector *temp8 = gsl_vector_alloc (lp);

      gsl_blas_dgemv (CblasTrans, 1.0, S11, Apn, 0.0, tempy1);
      gsl_blas_dgemv (CblasTrans, 1.0, S21, Acn, 0.0, tempy2);
      gsl_blas_ddot (tempy1, Sc, &y11);
      gsl_blas_ddot (tempy2, Sc, &y12);
      y1 = y11 + y12;
      gsl_blas_ddot (tempy1, Ss, &y21);
      gsl_blas_ddot (tempy2, Ss, &y22);
      y2 = y21 + y22;
      gsl_blas_dgemv (CblasTrans, 1.0, S12, Apn, 0.0, tempy3);
      gsl_blas_dgemv (CblasTrans, 1.0, S22, Acn, 0.0, tempy4);
      gsl_blas_ddot (tempy3, Sc, &y31);
      gsl_blas_ddot (tempy4, Sc, &y32);
      y3 = y31 + y32;
      gsl_blas_ddot (tempy3, Ss, &y41);
      gsl_blas_ddot (tempy4, Ss, &y42);
      y4 = y41 + y42;

      gsl_blas_dgemv (CblasTrans, 1.0, S11, Sc, 0.0, tempr1);
      gsl_blas_dgemv (CblasTrans, 1.0, S12, Sc, 0.0, tempr2);
      gsl_blas_ddot (tempr1, Sc, &r11);
      gsl_blas_ddot (tempr1, Ss, &r12);
      gsl_blas_ddot (tempr2, Sc, &r13);
      gsl_blas_ddot (tempr2, Ss, &r14);
      gsl_blas_dgemv (CblasTrans, 1.0, S11, Ss, 0.0, tempr3);
      gsl_blas_dgemv (CblasTrans, 1.0, S21, Sc, 0.0, tempr4);
      gsl_blas_dgemv (CblasTrans, 1.0, S21, Ss, 0.0, tempr5);
      gsl_blas_ddot (tempr3, Ss, &r22);
      gsl_blas_ddot (tempr4, Ss, &r23);
      gsl_blas_ddot (tempr5, Ss, &r24);
      gsl_blas_dgemv (CblasTrans, 1.0, S22, Sc, 0.0, tempr6);
      gsl_blas_dgemv (CblasTrans, 1.0, S22, Ss, 0.0, tempr7);
      gsl_blas_ddot (tempr6, Sc, &r33);
      gsl_blas_ddot (tempr6, Ss, &r34);
      gsl_blas_ddot (tempr7, Ss, &r44);
      r21 = r12;
      r31 = r13;
      r32 = r23;
      r41 = r14;
      r42 = r24;
      r43 = r34;

      double a_data[] = { r11, r12, r13, r14,
                          r21, r22, r23, r24,
                          r31, r32, r33, r34,
                          r41, r42, r43, r44 };

      double b_data[] = { y1, y2, y3, y4 };
      gsl_matrix_view m = gsl_matrix_view_array (a_data, 4, 4);
      gsl_vector_view b = gsl_vector_view_array (b_data, 4);
      gsl_vector *Cc = gsl_vector_alloc (4);
      int s;
      gsl_permutation * p = gsl_permutation_alloc (4);
      gsl_linalg_LU_decomp (&m.matrix, p, &s);
      gsl_linalg_LU_solve (&m.matrix, p, &b.vector, Cc);
      gsl_permutation_free (p);
      /*
      for (j = 0; j < 4; j++)
        printf ("Cc(%d) = %g\n", j+1, gsl_vector_get (Cc, j));
      */

      double Cc1 = gsl_vector_get (Cc, 0);
      double Cc2 = gsl_vector_get (Cc, 1);
      double Cc3 = gsl_vector_get (Cc, 2);
      double Cc4 = gsl_vector_get (Cc, 3);
      gsl_vector_free (Cc);

      gsl_vector *Ap = gsl_vector_alloc (lp);
      gsl_vector *Ac = gsl_vector_alloc (lp);
      for (j = 0; j < lp; j++)
      {
        gsl_vector_set (Ap, j, Cc1*cos(2.0*86400.0*M_PI*f*psr[0].quad_ifuncT_p[j]) + Cc2*sin(2.0*86400.0*M_PI*f*psr[0].quad_ifuncT_p[j]));
        gsl_vector_set (Ac, j, Cc3*cos(2.0*86400.0*M_PI*f*psr[0].quad_ifuncT_p[j]) + Cc4*sin(2.0*86400.0*M_PI*f*psr[0].quad_ifuncT_p[j]));
      }

      gsl_blas_dgemv (CblasTrans, 1.0, S11, Ap, 0.0, tempy5);
      gsl_blas_dgemv (CblasTrans, 1.0, S21, Ac, 0.0, tempy6);
      gsl_blas_dgemv (CblasTrans, 1.0, S12, Ap, 0.0, tempy7);
      gsl_blas_dgemv (CblasTrans, 1.0, S22, Ac, 0.0, tempy8);
      gsl_blas_dgemv (CblasTrans, 1.0, S11, Apn, 0.0, temp5);
      gsl_blas_dgemv (CblasTrans, 1.0, S21, Acn, 0.0, temp6);
      gsl_blas_dgemv (CblasTrans, 1.0, S12, Apn, 0.0, temp7);
      gsl_blas_dgemv (CblasTrans, 1.0, S22, Acn, 0.0, temp8);

      gsl_blas_ddot (tempy5, Ap, &d1);
      gsl_blas_ddot (tempy6, Ap, &d2);
      gsl_blas_ddot (tempy7, Ac, &d3);
      gsl_blas_ddot (tempy8, Ac, &d4);
      gsl_blas_ddot (temp5, Ap, &d5);
      gsl_blas_ddot (temp6, Ap, &d6);
      gsl_blas_ddot (temp7, Ac, &d7);
      gsl_blas_ddot (temp8, Ac, &d8);

      freq[i] = f;
      DS[i] = 2.0*(d5+d6+d7+d8)-(d1+d2+d3+d4);
      gsl_vector_free (tempy1);
      gsl_vector_free (tempy2);
      gsl_vector_free (tempy3);
      gsl_vector_free (tempy4);
      gsl_vector_free (tempy5);
      gsl_vector_free (tempy6);
      gsl_vector_free (tempy7);
      gsl_vector_free (tempy8);
      gsl_vector_free (temp5);
      gsl_vector_free (temp6);
      gsl_vector_free (temp7);
      gsl_vector_free (temp8);
      gsl_vector_free (tempr1);
      gsl_vector_free (tempr2);
      gsl_vector_free (tempr3);
      gsl_vector_free (tempr4);
      gsl_vector_free (tempr5);
      gsl_vector_free (tempr6);
      gsl_vector_free (tempr7);
      gsl_vector_free (Sc);
      gsl_vector_free (Ss);
      gsl_vector_free (Ac);
      gsl_vector_free (Ap);
    }
    gsl_matrix_free (S11);
    gsl_matrix_free (S12);
    gsl_matrix_free (S21);
    gsl_matrix_free (S22);
    gsl_vector_free (time);
    gsl_vector_free (Apn);
    gsl_vector_free (Acn);

    fout = fopen("DectionSts.dat","w");
    for (i=0;i<nSpecOS4;i++)
      fprintf(fout,"%d %g %g\n",i+1,freq[i],DS[i]);
    fclose(fout);
  }
  return 0;
}

// char * plugVersionCheck = TEMPO2_h_VER;
