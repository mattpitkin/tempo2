#include <stdlib.h>
#include <stdio.h>
#include <tempo2.h>
#include <math.h>
#include "TKfit.h"
#include "T2toolkit.h"
#include "choleskyRoutines.h"
#include "constraints.h"

void globalFITfuncs(double x,double afunc[],int ma,pulsar *psr,int ipos);
int gnpsr;
void multMatrix2(double **idcm,double **u,int ndata,int npol,double **uout);
void multMatrixVec2(double **idcm,double *b,int ndata,double *bout);
void TKsingularValueDecomposition_lsq2(double **designMatrix,int n,int nf,double **v,double *w,double **u);
void TKbacksubstitution_svd2(double **V, double *w,double **U,double *b,double *x,int n,int nf);
double TKpythag2(double a,double b);
void readUinv(int p,double **uinv,pulsar *psr,double *x,double *y,double *sig,int count,int nconstraints,int *ip);
void TKbidiagonal2(double **a,double *anorm,int ndata,int nfit,double **v,double *w,double **u,double *rv1);
void formCholeskyMatrix2(double *c,double *resx,double *resy,double *rese,int np,int nconstraints,double **uinv);

extern "C" int pluginFitFunc(pulsar *psr,int npsr,int writeModel) 
{
  int i,j,k,p,c0;
  double tol = 1.0e-27;  /* Tolerence for singular value decomposition routine */
  int npol=1,n,nf;
  int ip[MAX_OBSN];
  int nFitP[MAX_PSR];
  double *val,*error;
  double *x,*y,*sig,*sigOrig,**covar,***uinv;
  double whiteres[MAX_OBSN];
  long double toffset;
  double chisq;
  int weight=0;
  int count=0;
  FILE *fin,*fout;
  char fname[1000];
  int weightfit=0;
  double **v,**u,***uout;
  double sum,wmax;

  uinv = (double ***)malloc(sizeof(double **)*npsr);
  uout = (double ***)malloc(sizeof(double **)*npsr);

  //  printf("Setting weighting off\n");
  for (p=0;p<npsr;p++)
    psr[p].fitMode=1;  // Cholesky is a weighted fit

  gnpsr = npsr;

  x = (double *)malloc(MAX_OBSN*sizeof(double));
  y = (double *)malloc(MAX_OBSN*sizeof(double));
  sig = (double *)malloc(MAX_OBSN*sizeof(double));
  sigOrig = (double *)malloc(MAX_OBSN*sizeof(double));
  printf("HELLO\n");
  printf("About to undertake a global fit with whitening, number of pulsars = %d\n",npsr);


  // Whiten all the data sets
  int maxN=0;
  for (p=0;p<npsr;p++)
    {
      uinv[p] = malloc_uinv(psr[p].nobs+psr[p].nconstraints);      
      uout[p] = malloc_uinv(psr[p].nobs+psr[p].nconstraints);
      if (psr[p].nobs > maxN) maxN = psr[p].nobs+psr[p].nconstraints;
      printf("At this point: %d %d %d\n",p,psr[p].nobs,psr[p].nconstraints);
    }

  printf("Setting up matrices n(psr0)=%d\n",psr[0].nobs);
  p=0;


  //  for (i=0;i<maxN;i++)
  //    {
  //      uout[i] = (double *)malloc(sizeof(double)*maxN);
  //    }
  //      exit(1);

  printf("Done setting up matrices\n");
  // 
  /*  printf("Whitening the data\n");
  // Whiten the data
  fout = fopen("whiteres.dat","w");
  for (i=0;i<psr[0].nobs;i++)
    {
      sum=0.0;
      for (j=0;j<psr[0].nobs;j++)
	sum+=uinv[j][i]*psr[0].obsn[j].residual;
      whiteres[i] = sum;
      fprintf(fout,"%g %g\n",(double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]),whiteres[i]);
    }
    fclose(fout);*/
  
  // Form pre-fit residuals
  for (p=0;p<npsr;p++)
    {
      //      if (psr[p].fitMode==1) weightfit=1;
      //      if (weightfit==1 && psr[p].fitMode==0)  
      //	printf("WARNING: A weighted fit is being carried out, but PSR %s does not have MODE 1 in the parameter file\n",psr[p].name);
      for (i=0;i<psr[p].nobs;i++)
	{
	  if (psr[p].obsn[i].deleted!=0)
	    {
	      printf("ERROR: Please remove all deleted files in your TOA files before doing a global fit.  There is a problem with PSR %s\n",psr[p].name);
	      exit(1);	       
	    }
	  psr[p].obsn[i].prefitResidual = psr[p].obsn[i].residual;
	  x[count] = (double)(psr[p].obsn[i].bbat-psr[p].param[param_pepoch].val[0]);
	  y[count] = (double)(psr[p].obsn[i].prefitResidual);
	  ip[count] = i;
	  //	  if (psr[p].fitMode==0) sig[count] = 1.0;
	  //else sig[count] = 
	  sigOrig[count] = psr[p].obsn[i].toaErr*1.0e-6;
	  sig[count] =  1.0;
	  //	  printf("psr = %d: Reading in %g\n",p,sig[count]);
	  count++;
	}
      c0 = count;
           // add constraints as extra pseudo observations
     // These point to non-existant observations after the nobs array
     // These are later caught by getParamDeriv.
     for (i=0; i < psr[p].nconstraints; i++){
       ip[count] = psr[p].nobs+i; // psr->nobs+i;
	x[count]=x[c0-1];
	y[count]=0;
	sigOrig[count]=1e-12;
	sig[count]=1.0;    // Remember for the Cholesky that we switch off weighted fitting

	count++;
      }
     printf("count = %d\n",count);

    }
  for (p=0;p<npsr;p++)
    psr[p].nFit=count;
  printf("Total number of points = %d\n",count);

  // Determine number of fit parameters
  npol=0;
  int nGlobal=0;
  // Add global parameters
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[0].param[i].aSize;k++)
	{
	  if (psr[0].param[i].fitFlag[k]==2) {
	    {
	      npol++;
	      nGlobal++;
	    }
	  }
	}
    }
  printf("SHOULD ADD EXTRA PARAMETERS FOR JUMPS\n");
  /* Add extra parameters for sinusoidal whitening */
  if (psr[0].param[param_wave_om].fitFlag[0]==2)
    {
      printf("ADDING EXTRA PARAMETERS FOR WAVES: %d\n",psr[0].nWhite*2-1);
      npol+=psr[0].nWhite*2-1;
      nGlobal+=psr[0].nWhite*2-1;
    }

  if (psr[0].param[param_ifunc].fitFlag[0]==2){
    npol+=psr[0].ifuncN-1;
    nGlobal+=psr[0].ifuncN-1;
  }

  if (psr[0].param[param_quad_om].fitFlag[0]==2){
    npol+=psr[0].nQuad*4-1;
    nGlobal+=psr[0].nQuad*4-1;
  }

  if (psr[0].param[param_tel_dx].fitFlag[0]==2)	  
    {
      npol+=(psr[0].nTelDX-1);
      nGlobal+=(psr[0].nTelDX-1);
    }
  if (psr[0].param[param_tel_dy].fitFlag[0]==2)	  
    {
      npol+=(psr[0].nTelDY-1);      
      nGlobal+=(psr[0].nTelDY-1);
    }
  if (psr[0].param[param_tel_dz].fitFlag[0]==2)	  
    {
      npol+=(psr[0].nTelDZ-1);
      nGlobal+=(psr[0].nTelDZ-1);
    }
 if (psr->param[param_quad_ifunc_p].fitFlag[0]==2)
    {
      npol+=(psr->quad_ifuncN_p-1);
      nGlobal+=(psr->quad_ifuncN_p-1);
    }
  if (psr->param[param_quad_ifunc_c].fitFlag[0]==2)
    {
      npol+=(psr->quad_ifuncN_c-1);
      nGlobal+=(psr->quad_ifuncN_c-1);
    }
  if (psr[0].param[param_gwsingle].fitFlag[0]==2)
    {
      npol+=(4-1); 
      nGlobal+=(4-1);
    }
  printf("Number of global parameters = %d\n",nGlobal);



  // Add non-global parameters
  for (p=0;p<npsr;p++)
    {
      nFitP[p] = 1;
      npol++; // For the offset
      for (i=0;i<MAX_PARAMS;i++)
	{
	  for (k=0;k<psr[p].param[i].aSize;k++)
	    {
	      if (psr[p].param[i].fitFlag[k]==1) {
		if (i!=param_start && i!=param_finish)
		  {
		    nFitP[p]++;
		    npol++;
		  }
	      }
	    }
	}
      // Should do the same for jumps
      /* Add extra parameters for jumps */
      for (i=1;i<=psr[p].nJumps;i++)
	{
	  if (psr[p].fitJump[i]==1)
	    {
	      nFitP[p]++;
	      npol++;
	    }
	}
      /* Add extra parameters for sinusoidal whitening */
      if (psr[p].param[param_wave_om].fitFlag[0]==1)
	{
	  nFitP[p]+=psr[p].nWhite*2-1;
	  npol+=psr[p].nWhite*2-1;
	}
      if (psr[p].param[param_ifunc].fitFlag[0]==1){
	      npol+=psr[p].ifuncN-1;
	      nFitP[p]+=psr[p].ifuncN-1;
      }
      if (psr[p].param[param_quad_om].fitFlag[0]==1){
	      npol+=psr[p].nQuad*4-1;
	      nFitP[p]+=psr[p].nQuad*4-1;
      }
       if (psr[p].param[param_tel_dx].fitFlag[0]==1)
	{
	  npol+=(psr[p].nTelDX-1);
	  nFitP[p]+=(psr[p].nTelDX-1);
	}
       if (psr[p].param[param_tel_dy].fitFlag[0]==1)
	{
	  npol+=(psr[p].nTelDY-1);
	  nFitP[p]+=(psr[p].nTelDY-1);
	}
       if (psr[p].param[param_tel_dz].fitFlag[0]==1)
	{
	  npol+=(psr[p].nTelDZ-1);
	  nFitP[p]+=(psr[p].nTelDZ-1);
	}
       if (psr->param[param_quad_ifunc_p].fitFlag[0]==1)
	 {
	   npol+=(psr->quad_ifuncN_p-1);
	   nFitP[p]+=(psr->quad_ifuncN_p-1);
	 }
       if (psr->param[param_quad_ifunc_c].fitFlag[0]==1)
	 {
	   npol+=(psr->quad_ifuncN_c-1);
	   nFitP[p]+=(psr->quad_ifuncN_c-1);
	 }
       if (psr[p].param[param_gwsingle].fitFlag[0]==1)
	{
	  npol+=(4-1);
	  nFitP[p]+=(4-1);
	}

      printf("Number of non-global parameters for pulsar %d = %d\n",p,nFitP[p]);
    }
  val   = (double *)malloc(npol*sizeof(double));
  error = (double *)malloc(npol*sizeof(double));
  covar = (double **)malloc(npol*sizeof(double *));
  for (i=0;i<npol;i++)
    covar[i] = (double *)malloc(npol*sizeof(double));

  // Create design matrix
  n = count;
  nf = npol;
  printf("npol = %d\n",npol);

  double **designMatrix = (double **)malloc(n*sizeof(double *));
  double **designMatrixT = (double **)malloc(n*sizeof(double *));
  double w[nf],wt[nf];
  v = (double **)malloc(nf*sizeof(double *));
  u = (double **)malloc(n*sizeof(double *));
  for (i=0;i<n;i++) 
    {
      if ((designMatrix[i] = (double *)malloc(nf*sizeof(double))) == NULL) {printf("OUT OF MEMORY\n"); exit(1);}
      if ((designMatrixT[i] = (double *)malloc(nf*sizeof(double))) == NULL) {printf("OUT OF MEMORY\n"); exit(1);}
      if ((u[i] = (double *)malloc(nf*sizeof(double))) == NULL)  {printf("OUT OF MEMORY\n"); exit(1);}
    }
  for (i=0;i<nf;i++) v[i] = (double *)malloc(nf*sizeof(double));

  printf("Total number of fit parameters = %d\n",npol);
  
  // Now form individual design matrices
  int nfit=0;
  int offset=0;
  int joff=0,koff=0;
  double basisFunc[nf],b[n];
  double bout[n];

  for (i=0;i<n;i++)
    {
      for (j=0;j<nf;j++)
	designMatrix[i][j]=0.0;
    } 
  
  for (p=0;p<npsr;p++)
    {
      int nobs_and_constraints = psr[p].nobs+psr[p].nconstraints;
      printf("Pulsar %d, %d and %d offset = %d, nGlobal = %d, nFitP = %d\n",p,psr[p].nobs,psr[p].nconstraints,offset,nGlobal,nFitP[p]);
      for (i=0;i<nobs_and_constraints;i++)
	{
	  //	  printf("Here with %d %d\n",p,i);
	  // i was ip[i]
	  globalFITfuncs(x[i+offset],basisFunc,nf,psr,i+offset);
	  //	  	  for (j=0;j<nf;j++)
	  //	  	    printf("Have %d %d %g [%d %d %d]\n",p,j,basisFunc[j],koff,nGlobal,nFitP[p]);

	  for (j=0;j<nGlobal;j++)
	    designMatrixT[i][j]=basisFunc[j]/sig[i+offset];
	  for (j=0;j<nFitP[p];j++)
	    {
	      designMatrixT[i][j+nGlobal] = basisFunc[j+koff+nGlobal]/sig[i+offset];	      
	      //	      printf("p = %d, %g\n",p,basisFunc[j+koff+nGlobal]);
	    }
	  //	  printf("set up\n");
	  b[i+offset] = y[i+offset]/sig[i+offset];
	}
      koff+=nFitP[p];

      printf("Mult matrix\n");
      // Pre-whiten the model using the data covariance matrix
      printf("Attempting uinv1, count = %d\n",count);
      printf("ERROR: not checking for deleted points\n");
      readUinv(p,uinv[p],psr,x+offset,y+offset,sigOrig+offset,psr[p].nobs+psr[p].nconstraints,psr[p].nconstraints,ip+offset);
      printf("Read uinv1\n");
      multMatrix2(uinv[p],designMatrixT,nobs_and_constraints,nFitP[p]+nGlobal,uout[p]);
      printf("Finished multMatrix\n");
      for (i=0;i<nobs_and_constraints;i++)
	{
	  for (j=0;j<nGlobal;j++) designMatrix[i+offset][j] = uout[p][i][j];
	  for (j=0;j<nFitP[p];j++)
	    {
	      designMatrix[i+offset][j+joff+nGlobal] = uout[p][i][j+nGlobal];
	      //	      printf("design: %g\n",uout[p][i][j+nGlobal]);
	    }
	}
      joff+=nFitP[p];
      offset+=nobs_and_constraints;
    }
  //  exit(1);
  /*  for (i=0;i<n;i++)
    {
      for (j=0;j<nf;j++)
	printf("design %g\n",designMatrix[i][j]);
	}*/

  printf("Done calculating the design matrix\n");

  // Now pre-whiten the input data
  offset=0;
  for (p=0;p<npsr;p++)
    {
      int nobs_and_constraints = psr[p].nobs+psr[p].nconstraints;
      printf("Mult vec\n");
      readUinv(p,uinv[p],psr,x+offset,y+offset,sigOrig+offset,psr[p].nobs+psr[p].nconstraints,psr[p].nconstraints,ip+offset);
      printf("Read uinv2\n");
      multMatrixVec2(uinv[p],b+offset,nobs_and_constraints,bout);
      printf("Done mult vec\n");
      for (i=0;i<nobs_and_constraints;i++)
	{
	  b[i+offset] = bout[i];
	  //printf("Have finished with %g\n",b[i+offset]);
	}
      offset+=nobs_and_constraints;
    }
  /* Now carry out the singular value decomposition */
  printf("Doing singular value decomp\n");
  /*  for (i=0;i<n;i++)
    {
      for (j=0;j<nf;j++)
	printf("design %g\n",designMatrix[i][j]);
	}*/
  TKsingularValueDecomposition_lsq2(designMatrix,n,nf,v,w,u);
  printf("Done singular value decomp\n");
  wmax = TKfindMax_d(w,nf);
  for (i=0;i<nf;i++)
    {
      if (w[i] < tol*wmax) w[i]=0.0;
    }

  /* Back substitution */
  TKbacksubstitution_svd2(v, w, designMatrix, b, val, n, nf);
  for (i=0;i<nf;i++)
    {
      //      error[i]= 0.0;
      printf("val = %g\n",val[i]);
    }
  printf("Complete backsubstitution\n");

  // Now form the covariance matrix and determine uncertainties
  for (i=0;i<nf;i++)
    {
      if (w[i]!=0) wt[i] = 1.0/w[i]/w[i];
      else wt[i] = 0.0;     
    }
  for (i=0;i<nf;i++)
    {
      for (j=0;j<=i;j++)
	{
	  sum=0.0;
	  for (k=0;k<nf;k++)
	    sum+=v[i][k]*v[j][k]*wt[k];
	  covar[i][j] = covar[j][i] = sum;
	}
      error[i] = sqrt(covar[i][i]);
    }


  chisq = 0.0;


  // Reform the designMatrix
  for (i=0;i<n;i++)
    {
      for (j=0;j<nf;j++)
	designMatrix[i][j]=0.0;
    } 
  offset=0;
  koff=0, joff=0;
  for (p=0;p<npsr;p++)
    {
      int nobs_and_constraints = psr[p].nobs+psr[p].nconstraints;
      for (i=0;i<nobs_and_constraints;i++)
	{
	  // ip -> i
	  globalFITfuncs(x[i+offset],basisFunc,nf,psr,i+offset);
	  for (j=0;j<nGlobal;j++)
	    designMatrixT[i][j]=basisFunc[j];
	  for (j=0;j<nFitP[p];j++)
	    designMatrixT[i][j+nGlobal] = basisFunc[j+koff+nGlobal]; 
	  //	  b[i+offset] = y[i+offset];
	}
      koff+=nFitP[p];
      // Pre-whiten the model using the data covariance matrix
      printf("Trying uinv3\n");
      readUinv(p,uinv[p],psr,x+offset,y+offset,sigOrig+offset,psr[p].nobs+psr[p].nconstraints,psr[p].nconstraints,ip+offset);
      printf("Read uinv3\n");
      multMatrix2(uinv[p],designMatrixT,nobs_and_constraints,nFitP[p]+nGlobal,uout[p]);
      for (i=0;i<psr[p].nobs;i++)
	{
	  for (j=0;j<nGlobal;j++) designMatrix[i+offset][j] = uout[p][i][j];
	  for (j=0;j<nFitP[p];j++)
	    designMatrix[i+offset][j+joff+nGlobal] = uout[p][i][j+nGlobal];
	}
      joff+=nFitP[p];
      offset+=nobs_and_constraints;
    }


  double tmp;
  offset=0;
  //  for (p=0;p<npsr;p++)
  //    {
  for (j=0;j<count;j++)
    {
      sum=0.0;
      for (k=0;k<nf;k++)
	sum+=val[k]*designMatrix[j][k];
      chisq+=pow((b[j]-sum)/sig[j],2);
    }
      //      offset+=psr[p].nobs;
      //    }

  // MUST CALCULATE THIS WITH THE PRE-WHITENED DATA AND THE PRE-WHITENED MODEL
  /*  for (p=0;p<npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)
	{
	  globalFITfuncs(x[i+offset],basisFunc,nf,psr,ip[i+offset]);
	  sum=0.0;
	  for (k=0;k<nf;k++)
	    sum += val[k]*basisFunc[k];
	  //	  printf("func: %g %g %g %g\n",x[i+offset],y[i+offset],sig[i+offset],sum);
	  chisq += pow((y[i+offset]-sum)/sig[i+offset],2);
	}
      offset+=psr[p].nobs;
      } */
  //  for (i=0;i<n;i++)
  //    {
      //      globalFITfuncs(x[i],basisFunc,nf,psr,ip[i]);
      //      for (sum=0.0,j=0;j<nf;j++) sum += val[j]*designMatrix[i][j];
      //      printf("res1 = %g %g %g\n",x[i],b[i],sum);
  //      chisq += (tmp=(b[i]-sum),tmp*tmp);
  //    } 
  printf("chisq = %g, count = %d, nf = %d, red chisq = %g\n",chisq,count,nf,chisq/(count-nf));
  for (p=0;p<npsr;p++)
    {
      for (i=0;i<nf;i++)
	{
	  for (j=0;j<nf;j++)
	    psr[p].covar[i][j] = covar[i][j]*(chisq/(count-nf));
	}
    }
  printf("chisq = %g\n",chisq); 
  for (j=0;j<nf;j++)
    error[j]*=sqrt(chisq/(count-nf));
  // now update the parameters
  printf("Updating parameters\n");
  offset=0;
  // update global parameters
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[0].param[i].aSize;k++)
	{
	  if (psr[0].param[i].fitFlag[k] == 2)
	    {
	      printf("Have global parameter %d %d %d\n",i,param_wave_om,param_ifunc);
	      if (i==param_wave_om)
		{
		  int kk;
		  for (kk=0;kk<psr[0].nWhite;kk++)
		    {
		      for (p=0;p<npsr;p++)
			{
			  psr[p].wave_cos[kk]  -= val[offset]; 
			  psr[p].wave_cos_err[kk] = error[offset]; 
			  if (p==0) printf("Have wave %d %d %g %g\n",offset,kk,val[offset],error[offset]);
			}
		      offset++;
		      for (p=0;p<npsr;p++)
			{
			  psr[p].wave_sine[kk] -= val[offset]; 
			  psr[p].wave_sine_err[kk] = error[offset]; 
			}
		      offset++;
		    }
		  offset--;
		}
	      else if(i==param_ifunc) {
		printf("Updating %d point\n",psr[0].ifuncN);
		for (j=0;j<psr[0].ifuncN;j++)
		  {
		    printf("Updating %d %g\n",offset,val[offset]);
		    for (p=0;p<npsr;p++)
		      {
			psr[p].ifuncV[j]-=val[offset];
			psr[p].ifuncE[j]=error[offset];
		      }
		    offset++;
		  }
		offset--;
	      }
	      else if (i==param_quad_ifunc_p)
		  {
		  for (j=0;j<psr[0].quad_ifuncN_p;j++)
		    {
		      printf("Updating %g\n",val[offset]);
		      for (p=0;p<npsr;p++)
			{
			  psr[p].quad_ifuncV_p[j]-=val[offset];
			  psr[p].quad_ifuncE_p[j]=error[offset];
			}
		      offset++;
		    }
		  offset--;
		  }
	      else if (i==param_quad_ifunc_c)
		  {
		  for (j=0;j<psr[0].quad_ifuncN_c;j++)
		    {
		      printf("Updating %g\n",val[offset]);
		      for (p=0;p<npsr;p++)
			{
			  psr[p].quad_ifuncV_c[j]-=val[offset];
			  psr[p].quad_ifuncE_c[j]=error[offset];
			}
		      offset++;
		    }
		  offset--;
		  }
	      else if(i==param_quad_om) {
		printf("Updating %d point\n",psr[0].ifuncN);
		for (j=0;j<psr[0].nQuad;j++)
		  {
		    printf("Updating %d %g\n",offset,val[offset]);
		    for (p=0;p<npsr;p++)
		      {
			psr[p].quad_aplus_r[j]    -= val[offset];		      
			psr[p].quad_aplus_i[j]    -= val[offset+1];		      
			psr[p].quad_across_r[j]   -= val[offset+2];		      
			psr[p].quad_across_i[j]   -= val[offset+3];		      
			psr[p].quad_aplus_r_e[j]   = error[offset];		      
			psr[p].quad_aplus_i_e[j]   = error[offset+1];		      
			psr[p].quad_across_r_e[j]  = error[offset+2];		      
			psr[p].quad_across_i_e[j]  = error[offset+3];		      
		      }
		    offset+=4;
		  }
		offset--;
	      }
	      else if (i==param_tel_dx)
		{
		  for (j=0;j<psr[0].nTelDX;j++)
		    {
		      printf("Setting: %d %g\n",j,val[offset]);
		      for (p=0;p<npsr;p++)
			{
			  psr[p].telDX_v[j]-=val[offset];
			  psr[p].telDX_e[j]=error[offset];
			}
		      offset++;
		    }
		  offset--;
		}
	      else if (i==param_tel_dy)
		{
		  for (j=0;j<psr[0].nTelDY;j++)
		    {
		      printf("Setting: %d %g\n",j,val[offset]);
		      for (p=0;p<npsr;p++)
			{
			  psr[p].telDY_v[j]-=val[offset];
			  psr[p].telDY_e[j]=error[offset];
			}
		      offset++;
		    }
		  offset--;
		}
	      else if (i==param_tel_dz)
		{
		  for (j=0;j<psr[0].nTelDZ;j++)
		    {
		      printf("Setting: %d %g\n",j,val[offset]);
		      for (p=0;p<npsr;p++)
			{
			  psr[p].telDZ_v[j]-=val[offset];
			  psr[p].telDZ_e[j]=error[offset];
			}
		      offset++;
		    }
		  offset--;
		}
	      else if (i==param_gwsingle)
		{
		  for (p=0;p<npsr;p++)
		    {
		      psr[p].gwsrc_aplus_r -= val[offset];		      
		      psr[p].gwsrc_across_r -= val[offset+1];		      
		      psr[p].gwsrc_aplus_r_e = error[offset];		      
		      psr[p].gwsrc_across_r_e = error[offset+1];		      
		      psr[p].gwsrc_aplus_i -= val[offset+2];		      
		      psr[p].gwsrc_across_i -= val[offset+3];		      
		      psr[p].gwsrc_aplus_i_e = error[offset+2];		      
		      psr[p].gwsrc_across_i_e = error[offset+3];		      
		    }
		  offset+=3;
		}
	      else
		{
		  for (p=0;p<npsr;p++)
		    {
		      if (i==param_telx || i==param_tely || i==param_telz || i==param_gwm_amp)
			psr[p].param[i].val[k] -= val[offset];
		      else
			psr[p].param[i].val[k] += val[offset];
		      psr[p].param[i].err[k] = error[offset];
		    }
		}
	      offset++;
	    }
	}
    }
  for (p=0;p<npsr;p++)
    {
      updateParameters(psr,p,val+offset,error+offset);
      for (i=0;i<MAX_PARAMS;i++)
	{
	  for (k=0;k<psr[p].param[i].aSize;k++)
	    {
	      if (psr[p].param[i].fitFlag[k]==1) {
		if (i!=param_start && i!=param_finish)
		  offset++;
	      }
	    }
	}
      /* Add extra parameters for jumps */
      for (i=1;i<=psr[p].nJumps;i++)
	{
	  if (psr[p].fitJump[i]==1)
	    offset++;
	}
      /* Add extra parameters for sinusoidal whitening */
      if (psr[p].param[param_wave_om].fitFlag[0]==1)
	      offset+=psr[p].nWhite*2-1;
      if (psr[p].param[param_ifunc].fitFlag[0]==1)
	offset+=psr[p].ifuncN-1;
      if (psr[p].param[param_quad_om].fitFlag[0]==1)
	offset+=psr[p].nQuad*4-1;
      if (psr[p].param[param_quad_ifunc_p].fitFlag[0]==1)
	offset+=psr[p].quad_ifuncN_p-1;
      if (psr[p].param[param_quad_ifunc_c].fitFlag[0]==1)
	offset+=psr[p].quad_ifuncN_c-1;
      if (psr[p].param[param_tel_dx].fitFlag[0]==1)
	offset+=psr[p].nTelDX-1;
      if (psr[p].param[param_tel_dy].fitFlag[0]==1)
	offset+=psr[p].nTelDY-1;
      if (psr[p].param[param_tel_dz].fitFlag[0]==1)
	offset+=psr[p].nTelDZ-1;
      if (psr[p].param[param_gwsingle].fitFlag[0]==1)
	offset+=(4-1);

      offset++; // For arbitrary phase
    }

  // Now clean up
  free(error);
  free(val);
  free(sig);
  free(sigOrig);
  free(y);      
  free(x);
  //  free(bout);
  for (i=0;i<npsr;i++)
    {    
      free_uinv(uinv[i]);
      free_uinv(uout[i]);
    }
  //  for (i=0;i<maxN;i++)
  //    {
  //      free(uout[i]);
  //    }
  for (i=0;i<npol;i++)
    free(covar[i]);
  free(covar);
  //  free(uout);
  for (i=0;i<n;i++) 
    {
      free(designMatrix[i]);
      free(designMatrixT[i]);
      free(u[i]);
    }
  printf("Here\n");

  free(designMatrix);
  free(designMatrixT);
  free(u);
  printf("Got here\n");



  /*  exit(1);



  TKleastSquares_svd_psr(x,y,sig,count,val,error,npol,covar,&chisq,globalFITfuncs,weightfit,psr,tol,ip);

  // now update the parameters
  offset=0;
  // update global parameters
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[0].param[i].aSize;k++)
	{
	  if (psr[0].param[i].fitFlag[k] == 2)
	    {
	      for (p=0;p<npsr;p++)
		{
		  psr[p].param[i].val[k] += val[offset];
		  psr[p].param[i].err[k] = error[offset];
		}
	      offset++;
	    }
	}
    }
  for (p=0;p<npsr;p++)
    {
      updateParameters(psr,p,val+offset,error+offset);
      for (i=0;i<MAX_PARAMS;i++)
	{
	  for (k=0;k<psr[p].param[i].aSize;k++)
	    {
	      if (psr[p].param[i].fitFlag[k]==1) {
		if (i!=param_start && i!=param_finish)
		  offset++;
	      }
	    }
	}
      for (i=1;i<=psr[p].nJumps;i++)
	{
	  if (psr[p].fitJump[i]==1)
	    offset++;
	}
      if (psr[p].param[param_wave_om].fitFlag[0]==1)
	offset+=psr[p].nWhite*2-1;
      offset++; // For arbitrary phase
      }*/

  /*  free(x);
  free(y);
  free(sig);
  free(val);
  free(error);
  for (i=0;i<npol;i++)
    free(covar[i]);
    free(covar);*/
}

void readUinv(int p,double **uinv,pulsar *psr,double *x,double *y,double *sig,int count,int nconstraints,int *ip)
{
  int i,j;
  FILE *fin,*fout;
  char fname[1000];
  int ndays = (int)(x[count-nconstraints-1]-x[0])+2;
  double covarFunc[ndays];
  double whiteres[count];
  double sum;
  double escaleFactor = 1.0;
  char temp[1000];
  //  double **uinv2;

  //  uinv2=malloc_uinv(psr[p].nobs+psr[p].nconstraints);     
  //  printf("Using uinv2\n");
  // Read and obtain the covariance function
  printf("Getting covariance function for pulsar %d\n",p);
  
  //  strcpy(fname,covarFuncFile);
  strcpy(fname,"covarFunc.dat");
  sprintf(temp,"%s_%s",fname,psr[p].name);
  strcpy(fname,temp);
  //  for (i=0;i<count;i++)
  //    printf("Have %g %g %g %d\n",x[i],y[i],sig[i],ip[i]);
  printf("Using: %d %d %d\n",count,count-nconstraints,nconstraints);
  getCholeskyMatrix(uinv,fname,psr+p,x,y,sig,count,nconstraints,ip);
  //  if (p==1) exit(0);
  //  free_uinv(uinv2);
  /*  printf("Opening >%s<\n",fname);
  if (!(fin = fopen(fname,"r")))
    {
      printf("Unable to open covariance function file: %s\n",fname);
      exit(1);
    }
  fscanf(fin,"%lf",&escaleFactor);
  for (i=0;i<ndays;i++)
    fscanf(fin,"%lf",&covarFunc[i]);
  fclose(fin);
  printf("Read covariance function\n");
  printf("WARNING: scaling all errors by: %g\n",escaleFactor);
  for (i=0;i<count;i++)
    {
      sig[i]*=escaleFactor;
      //      printf("sig after correction = %g\n",sig[i]);
    }
  // Form the data covariance matrix
  formCholeskyMatrix2(covarFunc,x,y,sig,count,nconstraints,uinv); */

  sprintf(fname,"whitedata_%d.dat",p+1);
  //  if (p==1)
    {
      fout = fopen(fname,"w");
      for (i=0;i<count;i++)
	{
	  sum=0.0;
	  for (j=0;j<count;j++)
	    sum+=uinv[i][j]*psr[p].obsn[j].residual;
	  whiteres[i] = sum;
	  //            fprintf(fout,"%g %g %g\n",(double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]),(double)psr[0].obsn[i].residual,(double)psr[0].obsn[i].toaErr);
	  fprintf(fout,"%g %g %g %g\n",x[i],
		  whiteres[i],y[i],sig[i]);
	}
      fclose(fout);
      
      //      exit(1);
    }
    //    exit(1);
  /*
  strcpy(fname,"/DATA/BRAHE_1/hob044/idcm.dat_");
  strcat(fname,psr[p].name);
  printf("Looking for DCM file >%s<\n",fname);

  printf("Attempting psr %d\n",p);

  if (!(fin = fopen(fname,"r")))
    {
      printf("Unable to open inverse cholesky matrix: %s\n",fname);
      exit(1);
    }

  i=0;
  j=0;

  while (!feof(fin))
    {
      //      printf("Have %d %d %d %d\n",i,j,psr[p].nobs,p);
      if (fscanf(fin,"%lf",&uinv[j][i])==1)
	{
	  i++;
	  if (i==psr[p].nobs)
	    {
	      i=0;
	      j++;
	      if (j==psr[p].nobs+1)
		{
		  printf("The matrix file is the wrong size - too large\n");
		  printf("N_obs = %d\n",psr[p].nobs);
		  exit(1);
		}
	    }
	}      
    }
  if (j!=psr[p].nobs && i!=0)
    {
      printf("The matrix file is the wrong size - too short\n");
      printf("j = %d, i = %d, nobs = %d\n",j,i,psr[p].nobs);
      exit(1);
    }
    printf("Complete reading uinv\n"); */
}

void globalFITfuncs(double x,double afunc[],int ma,pulsar *psr,int counter)
{
  int i;
  int n=0;
  int p,pp;
  int new_ma;
  int ipos;
  int j,k;
  int tot=0;
  int nglobal=0;

  for (i=0;i<ma;i++) afunc[i]=0.0;

  for (p=0;p<gnpsr;p++)
    {
      if (counter < tot+psr[p].nobs + psr[p].nconstraints) break;
      tot+=psr[p].nobs+psr[p].nconstraints;
    }
  //  printf("In globalFit p = %d\n",p);
  ipos = counter-tot;
  //  printf("In glboalFit p = %d, ipos = %d, counter = %d, tot = %d\n",p,ipos,counter,tot);
  new_ma=1;
  // Add global parameters
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[0].param[i].aSize;k++)
	{
	  if (psr[0].param[i].fitFlag[k] == 2)
	    {
	      if (i==param_wave_om)
		nglobal+=psr[0].nWhite*2-1;
	      if (i==param_ifunc)
		 nglobal+=psr[p].ifuncN-1;
	      if (i==param_quad_om)
		 nglobal+=psr[p].nQuad*4-1;
	      if (i==param_quad_ifunc_p)
		nglobal+=psr[p].quad_ifuncN_p-1;
	      if (i==param_quad_ifunc_c)
		nglobal+=psr[p].quad_ifuncN_c-1;
	      if (i==param_tel_dx)
		 nglobal+=psr[p].nTelDX-1;
	      if (i==param_tel_dy)
		 nglobal+=psr[p].nTelDY-1;
	      if (i==param_tel_dz)
		 nglobal+=psr[p].nTelDZ-1;
	      if (i==param_gwsingle)
		nglobal+=(4-1);

	      nglobal++;
	    }
	}
    }
  new_ma+=nglobal;
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[p].param[i].aSize;k++)
	{
	  if (psr[p].param[i].fitFlag[k]==1) {
	    if (i!=param_start && i!=param_finish)
	      new_ma++;
	  }
	}
    }
  /* Add extra parameters for jumps */
  for (i=1;i<=psr[p].nJumps;i++)
    {
      if (psr[p].fitJump[i]==1)
	new_ma++;
    }
  /* Add extra parameters for sinusoidal whitening */
  if (psr[p].param[param_wave_om].fitFlag[0]==1)
    new_ma+=psr[p].nWhite*2-1;
  if (psr[p].param[param_ifunc].fitFlag[0]==1)
    new_ma+=psr[p].ifuncN-1;
  if (psr[p].param[param_quad_om].fitFlag[0]==1)
    new_ma+=psr[p].nQuad*4-1;
  if (psr[p].param[param_quad_ifunc_p].fitFlag[0]==1)
    new_ma+=psr[p].quad_ifuncN_p-1;
  if (psr[p].param[param_quad_ifunc_c].fitFlag[0]==1)
    new_ma+=psr[p].quad_ifuncN_c-1;
  if (psr[p].param[param_tel_dx].fitFlag[0]==1)
    new_ma+=psr[p].nTelDX-1;
  if (psr[p].param[param_tel_dy].fitFlag[0]==1)
    new_ma+=psr[p].nTelDY-1;
  if (psr[p].param[param_tel_dz].fitFlag[0]==1)
    new_ma+=psr[p].nTelDZ-1;
  if (psr[p].param[param_gwsingle].fitFlag[0]==1)
    new_ma+=(4-1);

  // Now calculate position in afunc array
  n=0;
  // Add global parameters
  n+=nglobal;
  for (pp=0;pp<p;pp++)
    {
      n++; // For the offset
      for (i=0;i<MAX_PARAMS;i++)
	{
	  for (k=0;k<psr[pp].param[i].aSize;k++)
	    {
	      if (psr[pp].param[i].fitFlag[k]==1) {
		if (i!=param_start && i!=param_finish)
		  n++;
	      }
	    }
	}
      // Should do the same for jumps
      /* Add extra parameters for jumps */
      for (i=1;i<=psr[pp].nJumps;i++)
	{
	  if (psr[pp].fitJump[i]==1)
	    n++;
	}
      /* Add extra parameters for sinusoidal whitening */
      if (psr[pp].param[param_wave_om].fitFlag[0]==1)
	n+=psr[pp].nWhite*2-1;
      if (psr[pp].param[param_ifunc].fitFlag[0]==1)
	n+=psr[pp].ifuncN-1;
      if (psr[pp].param[param_quad_om].fitFlag[0]==1)
	n+=psr[pp].nQuad*4-1;
      if (psr[pp].param[param_quad_ifunc_p].fitFlag[0]==1)
	n+=psr[pp].quad_ifuncN_p-1;
      if (psr[pp].param[param_quad_ifunc_c].fitFlag[0]==1)
	n+=psr[pp].quad_ifuncN_c-1;
      if (psr[pp].param[param_tel_dx].fitFlag[0]==1)
	n+=psr[pp].nTelDX-1;
      if (psr[pp].param[param_tel_dy].fitFlag[0]==1)
	n+=psr[pp].nTelDY-1;
      if (psr[pp].param[param_tel_dz].fitFlag[0]==1)
	n+=psr[pp].nTelDZ-1;
      if (psr[pp].param[param_gwsingle].fitFlag[0]==1)
	n+=(4-1);
    }

  // Global fit
  int c=0;
  int kk;
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[p].param[i].aSize;k++)
	{
	  if (psr[p].param[i].fitFlag[k] == 2)  
	    {
	      //	      afunc[c] = dotproduct(psr[p].posPulsar,psr[p].obsn[ipos].planet_ssb[4]);
	      if (i==param_wave_om)
		{
		  for (kk=0;kk<2*psr[0].nWhite;kk++)
		    {
		      afunc[c]= getParamDeriv(&psr[p],ipos,x+(double)psr[p].param[param_pepoch].val[0] - (double)psr[0].param[param_waveepoch].val[0],i,kk);
		      c++;
		    }
		}
	      else if(i==param_ifunc){
		for (j=0;j<psr[p].ifuncN;j++)
		  {
		    afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
		    //			      printf("ifc=%d %d %g\n",counter,c,afunc[c]);
		    c++;
		  }
		
	      }
	      else if(i==param_quad_ifunc_p){
		for (j=0;j<psr[p].quad_ifuncN_p;j++)
		  {
		    afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
		    c++;
		  }
		
	      }
	      else if(i==param_quad_ifunc_c){
		for (j=0;j<psr[p].quad_ifuncN_c;j++)
		  {
		    afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
		    c++;
		  }
		
	      }
	      else if(i==param_quad_om){
		for (j=0;j<psr[p].nQuad*4;j++)
		  {
		    afunc[c] = getParamDeriv(&psr[p],ipos,x+(double)psr[p].param[param_pepoch].val[0],i,j);
		    //			      printf("ifc=%d %d %g\n",counter,c,afunc[c]);
		    c++;
		  }
		
	      }
	      else if(i==param_tel_dx){
		for (j=0;j<psr[p].nTelDX;j++)
		  {
		    afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
		    c++;
		  }		
	      }
	      else if(i==param_tel_dy){
		for (j=0;j<psr[p].nTelDY;j++)
		  {
		    afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
		    c++;
		  }		
	      }
	      else if(i==param_tel_dz){
		for (j=0;j<psr[p].nTelDZ;j++)
		  {
		    afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
		    c++;
		  }		
	      }
	      else if (i==param_gwsingle)
		{
		  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,0);
		  c++;
		  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,1);
		  c++;
		  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,2);
		  c++;
		  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,3);
		  c++;
		}
	      else
		{
		  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,k);
		  c++;
		}
	    }
	} 
    }
  //  printf("Global fit = %g\n",afunc[0]);
  FITfuncs(x,afunc+n,new_ma-nglobal,&psr[p],ipos);
  //  printf("-----------------------------\n");
  //  for (i=0;i<ma;i++)
  //    printf("have %g\n",afunc[i]);
  //  n+=new_ma;
  //  exit(1);
}


void multMatrix2(double **idcm,double **u,int ndata,int npol,double **uout)
{
  int i,j,k;
  printf("HELLO MIKE\n");
  for (i=0;i<ndata;i++)
    {
      for (j=0;j<npol;j++)
	{
	  uout[i][j]=0.0;
	  for (k=0;k<ndata;k++)
	    uout[i][j]+=idcm[i][k]*u[k][j];
	}
    }
}

void multMatrixVec2(double **idcm,double *b,int ndata,double *bout)
{
  int i,j;
  for (i=0;i<ndata;i++)
    {
      bout[i] = 0.0;
      for (j=0;j<ndata;j++)
	bout[i]+=idcm[i][j]*b[j];
    }
}


/* Calculates SVD by following technique given in wikipedia */
void TKsingularValueDecomposition_lsq2(double **designMatrix,int n,int nf,double **v,double *w,double **u)
{
  double an;
  int i,j,k,its,l,nm,jj;
  int max_its = 40,pos1;
  double c,s,f,g,h,y,z,x;
  double rv1[nf];
  /* For A = U.W.V^T - obtain U, W and V */  
  
  /* Step 1: Reduce the matrix to a bidiagonal matrix */
  TKbidiagonal2(designMatrix,&an,n,nf,v,w,u,rv1);
  /* Step 2: Compute the SVD of the bidiagonal matrix */

  // Diagonalisation of the bidiagonal form: Loop over singular values
  // Code based on the numerical recipes routine
  for (k=nf-1;k>=0;k--)
    {
      for (its=1;its<=max_its;its++)
	{
	  pos1=0;
	  for (l=k;l>=0;l--)
	    {
	      nm=l-1;
	      if ((fabs(rv1[l])+an)==an) {pos1=2; break;}
	      if ((fabs(w[nm])+an)==an)  {pos1=1; break;}
	    }
	  if (pos1!=2)
	    {
	      c=0.0;
	      s=1.0;
	      for (i=l;i<=k;i++) // Check <= sign
		{
		  f=s*rv1[i];		  
		  rv1[i]=c*rv1[i];
		  if ((fabs(f)+an)==an) break;
		  g=w[i];
		  h=TKpythag2(f,g);
		  w[i]=h;
		  h=1.0/h;
		  c= (g*h);
		  s= -(f*h);
		  for (j=0;j<n;j++)
		    {
		      y = designMatrix[j][nm];
		      z = designMatrix[j][i];
		      designMatrix[j][nm] = (y*c)+(z*s);
		      designMatrix[j][i] = -(y*s)+(z*c);
		    }
		}
	    }	      
	  z = w[k];
	  if (l==k)
	    {
	      if (z < 0.0)
		{
		  w[k] = -z;
		  for (j=0;j<nf;j++)
		    v[j][k] = -v[j][k];
		}
	      break;
	    }
	  if (its==30)
	    {
	      printf("No convergence in singular value decomposition after 30 iterations\n");
	      exit(1);
	    }
	  x = w[l];
	  nm = k-1;
	  y = w[nm];
	  g =rv1[nm];
	  h= rv1[k];
	  f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	  g = TKpythag2(f,1.0);
	  f = ((x-z)*(x+z)+h*((y/(f+TKsign_d(g,f)))-h))/x;
	  c=1.0;
	  s=1.0;
	  for (j=l;j<=nm;j++)
	    {
	      i = j+1;
	      g = rv1[i];
	      y = w[i];
	      h = s*g;
	      g = c*g;
	      z = TKpythag2(f,h);
	      rv1[j] = z;
	      c = f/z;
	      s = h/z;
	      f = (x*c)+(g*s);
	      g = -(x*s)+(g*c);
	      h = y*s;
	      y = y*c;
	      for (jj=0;jj<nf;jj++)
		{
		  x = v[jj][j];
		  z = v[jj][i];
		  v[jj][j] = (x*c)+(z*s);
		  v[jj][i] = -(x*s)+(z*c);		  
		}
	      z = TKpythag2(f,h);
	      w[j] = z;
	      if (z != 0.0)
		{
		  z = 1.0/z;
		  c=f*z;
		  s=h*z;
		}
	      f = (c*g)+(s*y);
	      x = -(s*g)+(c*y);
	      for (jj=0;jj<n;jj++)
		{
		  y = designMatrix[jj][j];
		  z = designMatrix[jj][i];
		  designMatrix[jj][j] = (y*c)+(z*s);
		  designMatrix[jj][i] = -(y*s)+(z*c);		  
		}
	    }
	  rv1[l] = 0.0;
	  rv1[k] = f;
	  w[k]=x;
	}
    }

}


/* Use Householder reflections to do reduce the matrix to bidiagonal form       */
/* This is converted from the Fortran EISPACK code which is very similar        */
/* to the more recent numerical recipes code - originally this came from        */
/* the algol procedure svd, num. math. 14, 403-420 (1970) by Golub and Reinsch  */
/* The Fortran code was developed by Burton S. Garbow from the Argonne National */
/* laboratory (1983)                                                            */

void TKbidiagonal2(double **a,double *an,int ndata,int nfit,double **v,double *w,double **u,double *rv1)
{
  int i,j,k,l;
  double g=0.0;
  double scale=0.0,s,f,h;

  *an=0.0;
  for (i=0;i<nfit;i++)
    {
      l=i+1;
      rv1[i] = scale*g;
      g=0.0; s=0.0; scale = 0.0;
      if (i <= ndata)
	{
	  for (k=i;k<ndata;k++)
	    scale+=fabs(a[k][i]);
	  if (scale != 0.0)
	    {
	      for (k=i;k<ndata;k++)
		{
		  a[k][i]/=scale;
		  s+=a[k][i]*a[k][i];
		}
	      f=a[i][i];
	      g=-TKsign_d(sqrt(s),f);
	      h = f*g-s;
	      a[i][i]=f-g;
	      for (j=l;j<nfit;j++)
		{
		  s=0.0;
		  for (k=i;k<ndata;k++)
		    s+=a[k][i]*a[k][j];
		  f = s/h;
		  for (k=i;k<ndata;k++)
		      a[k][j]+=f*a[k][i];
		}
	      for (k=i;k<ndata;k++)
		  a[k][i]*=scale;	      

	    }
	}

      w[i]=scale*g;
      g=0.0;
      s=0.0;
      scale=0.0;
      if (i<=ndata && i != nfit)
	{
	  for (k=l;k<nfit;k++)
	    scale+=fabs(a[i][k]);
	  if (scale != 0.0)
	    {
	      for (k=l;k<nfit;k++)
		{
		  a[i][k] = a[i][k]/scale;
		  s+=a[i][k]*a[i][k];
		}
	      f=a[i][l];
	      g=-TKsign_d(sqrt(s),f);
	      h=f*g-s;
	      a[i][l] = f-g;
	      for (k=l;k<nfit;k++)
		rv1[k]=a[i][k]/h;
	      for (j=l;j<ndata;j++)
		{
		  s=0.0;
		  for (k=l;k<nfit;k++)
		    s+=a[j][k]*a[i][k];
		  for (k=l;k<nfit;k++)
		    a[j][k]+=s*rv1[k];
		}
	      for (k=l;k<nfit;k++)
		a[i][k]=scale*a[i][k];
	    }
	}
      *an=TKretMax_d(*an,(fabs(w[i])+fabs(rv1[i])));
    }

  // Accumulation of right-hand transformations
  for (i=nfit-1;i>=0;i--)
    {
      if (i < nfit-1)
	{
	  if (g != 0.0)
	    {
	      for (j=l;j<nfit;j++) v[j][i]=(a[i][j]/a[i][l])/g;
	      for (j=l;j<nfit;j++)
		{
		  s=0.0;
		  for (k=l;k<nfit;k++) s+=a[i][k]*v[k][j];
		  for (k=l;k<nfit;k++) v[k][j]+=s*v[k][i];
		}			 
	    }
	  for (j=l;j<nfit;j++)
	    v[i][j] = v[j][i] = 0.0;
	}
      v[i][i]=1.0;
      g = rv1[i];
      l=i;
    }

  // Accumulation of left-hand transformations
  for (i=TKretMin_i(ndata,nfit)-1;i>=0;i--)
    {
      l=i+1;
      g=w[i];
      for (j=l;j<nfit;j++) a[i][j]=0.0;
      if (g != 0.0)
	{
	  g=1.0/g;
	  for (j=l;j<nfit;j++)
	    {
	      s=0.0;
	      for (k=l;k<ndata;k++) s+=a[k][i]*a[k][j];
	      f = (s/a[i][i])*g;
	      for (k=i;k<ndata;k++) a[k][j]=a[k][j]+f*a[k][i];
	    }
	  for (j=i;j<ndata;j++) a[j][i]*=g;
	}
      else
	{
	  for (j=i;j<ndata;j++) a[j][i]=0.0;
	}
      a[i][i]++;
    }

}

/* Solves A.X = B for vector X using equation 2.6.7 in numerical recipes */
/* equation: x = V . [diag(1/w_j)] . (U^T.b)                             */ 
/*                                                                       */
/* Returns 'x'                                                           */
void TKbacksubstitution_svd2(double **V, double *w,double **U,double *b,double *x,int n,int nf)
{
  int i,j;
  double uTb[nf];

  /* Calculate [diag(1/w_j)] . U^T.b */
  for (i=0;i<nf;i++)
    {
      uTb[i]=0.0;
      if (w[i]!=0.0)
	{
	  for (j=0;j<n;j++)
	    uTb[i]+=U[j][i]*b[j];
	  uTb[i]/=w[i];
	}
    }
  /* Now multiply by V as in equation 2.6.7 */
  for (i=0;i<nf;i++)
    {
      x[i]=0.0;
      for (j=0;j<nf;j++)
	x[i]+=V[i][j]*uTb[j];
    }
}

/* Computes (a^2 + b^2)^1/2 */
double TKpythag2(double a,double b)
{
  double ret=0.0;
  double absa,absb;

  absa = fabs(a);
  absb = fabs(b);
  if (absa > absb)
    ret = absa*sqrt(1.0+pow(absb/absa,2));
  else
    {
      if (absb==0) ret = 0.0;
      else ret = absb*sqrt(1.0+pow(absa/absb,2));
    }
  return ret;
}

void formCholeskyMatrix2(double *c,double *resx,double *resy,double *rese,int np, int nc,double **uinv)
{
  double **m,**u,sum;
  double *cholp;
  int i,j,k,ix,iy;
  double t0,cint,t;
  int t1,t2;
  int debug=0;

  printf("Getting the covariance matrix in doFit, np = %d \n",np);
  m = (double **)malloc(sizeof(double *)*(np+1));
  u= (double **)malloc(sizeof(double *)*(np+1));
  cholp  = (double *)malloc(sizeof(double)*(np+1));  // Was ndays

  for (i=0;i<np+1;i++)
    {
      m[i] = (double *)malloc(sizeof(double)*(np+1));
      u[i] = (double *)malloc(sizeof(double)*(np+1));
    }
  


  for (ix=0;ix<np;ix++)
    {
      for (iy=0;iy<np;iy++)
	m[ix][iy] = fabs(resx[ix]-resx[iy]);
    }
  if (debug==1)
    {
      printf("First m = \n");
      for (i=np-5;i<np;i++)
	{ 
	  for (j=np-5;j<np;j++) printf("%10g ",m[i][j]); 
	  printf("\n");
	}

    }
  
  // Insert the covariance which depends only on the time difference.
  // Linearly interpolate between elements on the covariance function because
  // valid covariance matrix must have decreasing off diagonal elements.
  printf("Inserting into the covariance matrix\n");
  for (ix=0;ix<np;ix++)
    {
      for (iy=0;iy<np;iy++)
	{
	  t0 = m[ix][iy];
	  t1 = (int)floor(t0);
	  t2 = t1+1;
	  t  = t0-t1;
	  cint = c[t1]*(1-t)+c[t2]*t; // Linear interpolation
	  m[ix][iy] = cint;
	}
    }


  // add the values for the constraints
  // Constraints are not covariant with anything so it's all zero!
  for (i=np-nc; i < np; i++){
          for (j=0; j < np; j++){
		m[i][j]=0;
		m[j][i]=0;
          }
  }


  printf("Multiplying by errors\n");
  for (ix=0;ix<np;ix++)
    m[ix][ix]+=rese[ix]*rese[ix];
  printf("Complete the multipication\n");
  if (debug==1)
    {
      printf("m = \n\n");
      for (i=np-5;i<np;i++)
	{ 
	  for (j=np-5;j<np;j++) printf("% 10g ",m[i][j]); 
	  printf("\n");
	}
    }

  // Do the Cholesky
  printf("Doing the Cholesky decomposition\n");
  T2cholDecomposition(m,np,cholp);
  printf("Done the Cholesky decompositon\n");
  // Now calculate uinv
  for (i=0;i<np;i++)
    {
      m[i][i] = 1.0/cholp[i];
      uinv[i][i] = m[i][i];
      for (j=0;j<i;j++)
      	uinv[i][j] = 0.0;
      for (j=i+1;j<np;j++)
	{
	  sum=0.0;
	  for (k=i;k<j;k++) sum-=m[j][k]*m[k][i];
	  m[j][i]=sum/cholp[j];
	  uinv[i][j] = m[j][i];
	}
    } 
  printf("Created univ\n");
  if (debug==1)
    {
      printf("uinv = \n\n");
      for (i=0;i<5;i++)
	{ 
	  for (j=0;j<5;j++) printf("%10g ",uinv[i][j]); 
	  printf("\n");
	}
    }

  printf("Completed inverting the matrix\n");

  // Should free memory not required
  // (note: not freeing uinv)

  for (i=0;i<np+1;i++)
    {
      free(m[i]);
      free(u[i]);
    }
  free(m);
  free(u);
  free(cholp);
}
char * plugVersionCheck = TEMPO2_h_VER;
