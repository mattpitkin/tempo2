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
#include <cpgplot.h>
#include "T2toolkit.h"
#include "TKspectrum.h"
#include "TKfit.h"
#include "fftw3.h"

using namespace std;

typedef struct sample {
  double x;
  double y;
  double e;
  double pred;
  int    actual;
} sample;



void plotResiduals(pulsar *psr,sample *samples,int nSample);
void plotModel(pulsar *psr,double startSample,double endSample,double spacingSample,sample *samples,int nSamples);
void getPowerSpectra(pulsar *psr,double modelA,double modelFc,double modelAlpha,double startSample,double endSample,double *covFunc,int *nCovFunc,sample *samples,int nSampleTimes);
void sortSamples(sample *s,int n);
void choldc(double **a, int n,double *p);
void lubksb(double **a, int n, int *indx, double b[]);
void ludcmp(double **a, int n, int *indx, double *d);
void matrixMult(double **m1,double **m2,int n);


void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,j,k;
  double globalParameter;
  double startSample,endSample,spacingSample;
  double modelA,modelFc,modelAlpha,x;
  int    setModelA,setModelFc,setModelAlpha;
  double covFunc[15000];
  sample samples[15000];
  int nSampleTimes=0;
  int nCovFunc;
  double **covMatrix,t,det;
  double **cn,**ucninv,**cninv,**mat,**imat;
  double *cholp,*vec;
  int t1,t2;
  double cint,sum;
  double rms;

  startSample = 1;
  endSample = 0;
  spacingSample = -1;
  setModelA=setModelFc=setModelAlpha=-1;

  const char *CVS_verNum = "$Revision$";


  if (displayCVSversion == 1) CVSdisplayVersion("interpolate.C","plugin",CVS_verNum);

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: interpolate\n");
  printf("Author:              Xinping Deng, G. Hobbs\n");
  printf("CVS Version:         $Revision$\n");
  printf(" --- type 'h' for help information\n");

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[0],argv[i+1]); 
	  strcpy(timFile[0],argv[i+2]);
	}
      else if (strcmp(argv[i],"-x1")==0)
	sscanf(argv[++i],"%lf",&startSample);
      else if (strcmp(argv[i],"-x2")==0)
	sscanf(argv[++i],"%lf",&endSample);
      else if (strcmp(argv[i],"-dx")==0)
	sscanf(argv[++i],"%lf",&spacingSample);
      else if (strcmp(argv[i],"-a")==0)
	{
	  sscanf(argv[++i],"%lf",&modelA);
	  setModelA = 1;
	}
      else if (strcmp(argv[i],"-fc")==0)
	{
	  sscanf(argv[++i],"%lf",&modelFc);
	  setModelFc = 1;
	}
      else if (strcmp(argv[i],"-alpha")==0)
	{
	  sscanf(argv[++i],"%lf",&modelAlpha);
	  setModelAlpha = 1;
	}
    }

  if (setModelA == 0)
    {
      printf("Must set A with -a option\n");
      exit(1);
    }
  if (setModelFc == 0)
    {
      printf("Must set Fc with -fc option\n");
      exit(1);
    }
  if (setModelAlpha == 0)
    {
      printf("Must set Alpha with -alpha option\n");
      exit(1);
    }

  if (endSample < startSample)
    {
      printf("x2 must be greater than x1 (use -x1 and -x2 options)\n");
      exit(1);
    }
  if (spacingSample < 1)
    {
      printf("Must set a sample spacing: use -dx option\n");
      exit(1);      
    }
  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      if (i==0) doFit(psr,*npsr,0);   /* Do the fitting     */
      else textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
    }
  rms = 0.0;

  for (i=0;i<psr[0].nobs;i++)
    {
      samples[i].x = (double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
      samples[i].y = (double)(psr[0].obsn[i].residual);
      samples[i].e = (double)(psr[0].obsn[i].toaErr*1.0e-6);
      samples[i].actual = 1;
      rms += 1.0/samples[i].e/samples[i].e;
    }
  rms /= (double)psr[0].nobs;
  rms = sqrt(1.0/rms);
  printf("rms = %g %d\n",rms,psr[0].nobs);

  nSampleTimes = psr[0].nobs;
  // Add sample times for requested data
  for (x=startSample;x<=endSample;x+=spacingSample)
    {
      samples[nSampleTimes].x = (double)(x-psr[0].param[param_pepoch].val[0]);
      samples[nSampleTimes].y = 0.0; // Set the requested residuals to zero
      samples[nSampleTimes].e = 1e6*rms;  // Set the requested errors to a large number -- SET PROPERLY
      samples[nSampleTimes].actual = 0;
      nSampleTimes++;
    }
  sortSamples(samples,nSampleTimes);
  //    for (i=0;i<nSampleTimes;i++)
  //      printf("%g %g\n",samples[i].x,samples[i].y);

  // Create the power spectrum
  getPowerSpectra(psr,modelA,modelFc,modelAlpha,startSample,endSample,covFunc,&nCovFunc,samples,nSampleTimes);

  // Create the covariance matrix
  covMatrix = (double **)malloc(sizeof(double *)*nSampleTimes);
  cn = (double **)malloc(sizeof(double *)*nSampleTimes);
  ucninv = (double **)malloc(sizeof(double *)*nSampleTimes);
  cninv = (double **)malloc(sizeof(double *)*nSampleTimes);
  mat = (double **)malloc(sizeof(double *)*nSampleTimes);
  imat = (double **)malloc(sizeof(double *)*nSampleTimes);
  cholp = (double *)malloc(sizeof(double)*nSampleTimes);
  vec = (double *)malloc(sizeof(double)*nSampleTimes);
  for (i=0;i<nSampleTimes;i++)
    {
      covMatrix[i] = (double *)malloc(sizeof(double)*nSampleTimes);
      cn[i] = (double *)malloc(sizeof(double)*nSampleTimes);
      ucninv[i] = (double *)malloc(sizeof(double)*nSampleTimes);
      cninv[i] = (double *)malloc(sizeof(double)*nSampleTimes);
      mat[i] = (double *)malloc(sizeof(double)*nSampleTimes);
      imat[i] = (double *)malloc(sizeof(double)*nSampleTimes);
    }
  for (i=0;i<nSampleTimes;i++)
    {
      for (j=0;j<nSampleTimes;j++)
	{
	  covMatrix[i][j] = fabs(samples[i].x-samples[j].x);
	  if (i==j)
	    cn[i][j] = pow(samples[i].e,2);
	  else
	    cn[i][j]=0.0;
	}
    }

   printf("\n\ncn ... \n\n");
  for (i=0;i<5;i++)
    {
      for (j=0;j<5;j++)
	{
	  printf("%g ",cn[i][j]);
	}
      printf("\n");
    }


  for (i=0;i<nSampleTimes;i++)
    {
      for (j=0;j<nSampleTimes;j++)
	{
	  t = covMatrix[i][j];
	  t1 = floor(t);
	  t2 = t1+1;
	  det = (double)(t-t1);
	  cint = covFunc[t1]*(1-det) + covFunc[t2]*det;  // Check starting from 0 and not 1
	  covMatrix[i][j] = cint;
	}
    }
  printf("\n\nCovariance matrix ... %d\n\n",nSampleTimes);
  for (i=0;i<5;i++)
    {
      for (j=0;j<5;j++)
	{
	  printf("%g ",covMatrix[i][j]);
	}
      printf("\n");
    }

  TKcholDecomposition(cn,nSampleTimes,cholp);
  // Now calculate inverse
  for (i=0;i<nSampleTimes;i++)
    {
      cn[i][i] = 1.0/cholp[i];
      ucninv[i][i] = cn[i][i];
      for (j=0;j<i;j++)
      	ucninv[i][j] = 0.0;
      for (j=i+1;j<nSampleTimes;j++)
	{
	  sum=0.0;
	  for (k=i;k<j;k++) sum-=cn[j][k]*cn[k][i];
	  cn[j][i]=sum/cholp[j];
	  ucninv[i][j] = cn[j][i];
	  cninv[i][j] = cn[j][i];
	}
    } 


   printf("\n\nUcninv ... \n\n");
  for (i=133;i<143;i++)
    {
      for (j=133;j<143;j++)
	{
	  printf("%g ",ucninv[i][j]);
	}
      printf("\n");
    }

  for (i=0;i<nSampleTimes;i++)
    {
      for (j=0;j<nSampleTimes;j++)
	cninv[i][j] = ucninv[j][i]*ucninv[i][j];
    }
   printf("\n\n cninv ... (%d)\n\n",nSampleTimes);
  for (i=133;i<143;i++)
    {
      for (j=133;j<143;j++)
	{
	  printf("%g ",cninv[i][j]);
	}
      printf("\n");
    }
  // Calculate final predicted residual
  for (i=0;i<nSampleTimes;i++)
    {
      for (j=0;j<nSampleTimes;j++)
	{
	  mat[i][j]=0.0;
	  for (k=0;k<nSampleTimes;k++)
	    mat[i][j]+=covMatrix[k][j]*cninv[i][k];
	}
    }
  printf("mat \n\n");
  for (i=0;i<5;i++)
    {
      for (j=0;j<5;j++)
	printf("%g ",mat[i][j]);
      printf("\n");
    }
  for (i=0;i<nSampleTimes;i++)
    {
      vec[i] = 0.0;
      for (j=0;j<nSampleTimes;j++)
	vec[i] += mat[j][i]*samples[j].y;
      mat[i][i] += 1.0;
      printf("vec: %d %g\n",i,vec[i]);
    }
  printf("mat2 \n\n");
  for (i=0;i<5;i++)
    {
      for (j=0;j<5;j++)
	printf("%g ",mat[i][j]);
      printf("\n");
    }
  // Calculate inverse of mat
  {
    double d,col[nSampleTimes];
    int indx[nSampleTimes];

    printf("Step 1\n");
    ludcmp(mat,nSampleTimes,indx,&d);
    printf("Step 2\n");
    // CHECK -1 IN NUMERICAL RECIPES
    for (j=0;j<nSampleTimes;j++)
      {
	for (i=0;i<nSampleTimes;i++) col[i] = 0.0;
	col[j]=1.0;
	lubksb(mat,nSampleTimes,indx,col);
	for (i=0;i<nSampleTimes;i++)
	  imat[i][j] = col[i];
      }
    printf("Step 3\n");
  }
  // Residual = mat^-1 * vec
  for (i=0;i<nSampleTimes;i++)
    {
      samples[i].pred = 0.0;
      for (j=0;j<nSampleTimes;j++)
	samples[i].pred += imat[j][i]*vec[j]; 
      printf("Answer = %g %g %g\n",samples[i].x,samples[i].y,samples[i].pred);
    }
  

  // Plot the residuals
  cpgbeg(0,"3/xs",1,1);
  plotResiduals(psr,samples,nSampleTimes);
  // Plot the model
  plotModel(psr,startSample,endSample,spacingSample,samples,nSampleTimes);
  cpgend();

  // Write out new par file with the IFUNC commands
  psr[0].param[param_ifunc].paramSet[0] = 1;
  psr[0].param[param_ifunc].val[0] = 2;
  psr[0].param[param_ifunc].fitFlag[0] = 0;
  psr[0].ifuncN = 0;
  for (i=0;i<nSampleTimes;i++)
    {
      if (samples[i].actual==0)
	{
	  psr[0].ifuncT[psr[0].ifuncN] = samples[i].x+(double)psr[0].param[param_pepoch].val[0];
	  psr[0].ifuncV[psr[0].ifuncN] = -samples[i].pred;
	  psr[0].ifuncE[psr[0].ifuncN] = 0;
	  psr[0].ifuncN++;
	}
    }
  textOutput(psr,1,0,0,0,1,"");

  // Free the memory
  for (i=0;i<nSampleTimes;i++)
    {
      free(covMatrix[i]);
      free(cn[i]);
      free(ucninv[i]);
      free(cninv[i]);
      free(mat[i]);
      free(imat[i]);
    }
  free(covMatrix);
  free(cn);
  free(ucninv);
  free(cninv);
  free(cholp);
  free(mat);
  free(imat);
  free(vec);
  return 0;
}

void getPowerSpectra(pulsar *psr,double modelA,double modelFc,double modelAlpha,double startSample,double endSample,double *covFunc,int *nCovFunc,sample *samples,int nSampleTimes)
{
  int i,j;
  int ndays;
  double pwr;
  float *fx,*fy,*freq;
  double *opf;
  double *pf;

  ndays = (int)(samples[nSampleTimes-1].x - samples[0].x + 0.5);
  printf("ndays = %d\n",ndays);
  // 6220
  if (!(fx = (float *)malloc(sizeof(float)*(ndays*2+2))))
    {
      printf("Unable to allocate memory\n");
      exit(1);
    }
  fy = (float *)malloc(sizeof(float)*(ndays*2+2));
  freq = (float *)malloc(sizeof(float)*(ndays*2+2));
  opf = (double *)malloc(sizeof(double)*(ndays*2+2));
  if (!(pf = (double *)malloc(sizeof(double)*(ndays*2+2))))
    {
      printf("Unable to allocate memory\n");
      exit(1);
    }

  *nCovFunc = ndays;

  for (i=0;i<ndays+1;i++)
    {
      freq[i] = (i)/(2*ndays/365.25); // remove +1
      pwr  = modelA/pow(1+pow(freq[i]/modelFc,modelAlpha/2.0),2);
      pf[i] = pwr;
      fx[i] = (float)(i+1);
      fy[i] = (float)log10(pwr);
    }
  j = ndays+1;
  for (i=ndays;i>0;i--)
    {
      // NOTE: NOT SETTING TO NEGATIVE FREQUENCIES ....

      freq[j] = freq[i];
      pwr  = modelA/pow(1+pow(freq[j]/modelFc,modelAlpha/2.0),2);
      printf("Here with %d %d %g %g\n",i,j,freq[j],pwr);
      pf[j] = pwr;
      fx[j] = (float)(j);
      fy[j] = (float)log10(pwr);
      j++;
    } 
  ndays=j;
  for (i=0;i<ndays;i++)
    printf("freq %d %g %g %d\n",i,freq[i],pf[i],ndays);
  printf("Got here\n");
  cpgbeg(0,"1/xs",1,2);
  cpgsch(1.4);
  cpgenv(TKfindMin_f(fx,ndays),TKfindMax_f(fx,ndays),TKfindMin_f(fy,ndays),TKfindMax_f(fy,ndays),0,20);
  cpglab("Frequency channel","PSD","");
  cpgline(ndays,fx,fy);

  // Get covariance function
  {
    fftw_complex* output;
    fftw_plan transform_plan;
    
    output = (fftw_complex*)opf;
    transform_plan = fftw_plan_dft_r2c_1d(ndays, pf, output, FFTW_ESTIMATE);
    fftw_execute(transform_plan);    
    fftw_destroy_plan(transform_plan);  
    for (i=0;i<ndays;i++) 
      {
	covFunc[i] = opf[2*i]/ndays*pow(86400.0*365.25,2)*365.25; // Note: ndays**2 scaling is because of IFFT --- CHECK THIS
	//	printf("covFunc: %d %g %g %d\n",i,opf[2*i],opf[2*i+1],j);
	fx[i] = (float)i;
	fy[i] = (float)covFunc[i];
	printf("Output %g %g\n",fx[i],fy[i]);
      }
    ndays/=2;
    printf("ndays = %d\n",ndays);
    cpgenv(TKfindMin_f(fx,ndays),TKfindMax_f(fx,ndays),TKfindMin_f(fy,ndays),TKfindMax_f(fy,ndays),0,1);
    cpglab("Lag","Covariance function","");
    cpgline(ndays,fx,fy);
  }

  cpgend();

  free(fx); free(fy); free(freq); free(opf); free(pf);


}

void plotModel(pulsar *psr,double startSample,double endSample,double spacingSample,sample *samples,int nSamples)
{
  int i;
  double x;
  float fx[1024],fy[1024];
  int n=0;

  for (i=0;i<nSamples;i++)
    {
      if (samples[i].actual == 0)
	{
	  fx[n] = (float)samples[i].x;
	  fy[n] = (float)samples[i].pred;
	  n++;
	}
    }
  cpgsci(2);  cpgline(n,fx,fy); cpgsci(1);
  //cpgpt(n,fx,fy,5);

}


void plotResiduals(pulsar *psr,sample *samples,int nSample)
{
  float fx[psr[0].nobs],fy[psr[0].nobs];
  float minx,maxx;
  int n,i;
  
  n = 0;

  for (i=0;i<nSample;i++)
    {
      if (i==0)
	{
	  minx = maxx = samples[i].x;
	}
      else
	{
	  if (minx > samples[i].x)  minx = samples[i].x;
	  if (maxx < samples[i].x)  maxx = samples[i].x;
	}
      if (samples[i].actual==1)
	{
	  fx[n] = (float)(samples[i].x);
	  fy[n] = (float)(samples[i].y);
	  n++;
	}
    }


  cpgenv(minx,maxx,TKfindMin_f(fy,n),TKfindMax_f(fy,n),0,1);
  cpglab("Day","Residual (s)","");
  cpgpt(n,fx,fy,9);

}

void sortSamples(sample *s,int n)
{
  sample t;
  int swapped=0;
  int i;

  do {
    swapped=0;
    for (i=0;i<n-1;i++)
      {
	if (s[i].x > s[i+1].x)
	  {
	    t = s[i+1];
	    s[i+1] = s[i];
	    s[i] = t;
	    swapped=1;
	  }
      }
  } while (swapped==1);
  
  
}

char * plugVersionCheck = TEMPO2_h_VER;

#define NRANSI
#define TINY 1.0e-20;
#define NR_END 1
#define FREE_ARG char*

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}
double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}



void ludcmp(double **a, int n, int *indx, double *d)
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv;
  
  vv=vector(0,n-1);
  *d=1.0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
    vv[i]=1.0/big;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n-1) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<n;i++) a[i][j] *= dum;
    }
  }
  free_vector(vv,0,n-1);
}
#undef TINY
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software #p21E6W)1.1&iE10(9p#. */


// MUST CHECK ALL STARTS FROM 0 CAREFULLY
void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=-1,ip,j;
	double sum;

	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii!=-1)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software #p21E6W)1.1&iE10(9p#. */
