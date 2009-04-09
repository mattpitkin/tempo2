//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russel Edwards

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
#include <tempo2.h>
#include "TKfit.h"

void globalFITfuncs(double x,double afunc[],int ma,pulsar *psr,int ipos);
int gnpsr;

extern "C" int pluginFitFunc(pulsar *psr,int npsr,int writeModel) 
{
  int i,j,k,p;
  double tol = 1.0e-40;  /* Tolerence for singular value decomposition routine */
  int npol=1;
  int ip[MAX_OBSN];
  double *val,*error;
  double *x,*y,*sig,**covar;
  long double toffset;
  double chisq;
  int offset;
  int count=0;
  int weightfit=0;

  gnpsr = npsr;

  x = (double *)malloc(MAX_OBSN*sizeof(double));
  y = (double *)malloc(MAX_OBSN*sizeof(double));
  sig = (double *)malloc(MAX_OBSN*sizeof(double));
  printf("HELLO\n");
  printf("About to undertake a global fit, number of pulsars = %d\n",npsr);

  // Form pre-fit residuals
  for (p=0;p<npsr;p++)
    {
      if (psr[p].fitMode==1) weightfit=1;
      if (weightfit==1 && psr[p].fitMode==0)  
	printf("WARNING: A weighted fit is being carried out, but PSR %s does not have MODE 1 in the parameter file\n",psr[p].name);
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
	  if (psr[p].fitMode==0) sig[count] = 1.0;
	  else sig[count] = psr[p].obsn[i].toaErr*1.0e-6;
	  ip[count]=count;
	  count++;
	}
    }
  for (p=0;p<npsr;p++)
    psr[p].nFit=count;
  printf("Total number of points = %d\n",count);

  // Determine number of fit parameters
  npol=0;
  // Add global parameters
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[0].param[i].aSize;k++)
	{
	  if (psr[0].param[i].fitFlag[k]==2) {
	    npol++;
	  }
	}
    }
  for (p=0;p<npsr;p++)
    {
      npol++; // For the offset
      for (i=0;i<MAX_PARAMS;i++)
	{
	  for (k=0;k<psr[p].param[i].aSize;k++)
	    {
	      if (psr[p].param[i].fitFlag[k]==1) {
		if (i!=param_start && i!=param_finish)
		  npol++;
	      }
	    }
	}
      // Should do the same for jumps
      /* Add extra parameters for jumps */
      for (i=1;i<=psr[p].nJumps;i++)
	{
	  if (psr[p].fitJump[i]==1)
	    npol++;
	}
      /* Add extra parameters for sinusoidal whitening */
      if (psr[p].param[param_wave_om].fitFlag[0]==1)
	npol+=psr[p].nWhite*2-1;
    }
  val   = (double *)malloc(npol*sizeof(double));
  error = (double *)malloc(npol*sizeof(double));
  covar = (double **)malloc(npol*sizeof(double *));
  for (i=0;i<npol;i++)
    covar[i] = (double *)malloc(npol*sizeof(double));

  printf("Number of fit parameters = %d\n",npol);
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
      /* Add extra parameters for jumps */
      for (i=1;i<=psr[p].nJumps;i++)
	{
	  if (psr[p].fitJump[i]==1)
	    offset++;
	}
      /* Add extra parameters for sinusoidal whitening */
      if (psr[p].param[param_wave_om].fitFlag[0]==1)
	offset+=psr[p].nWhite*2-1;
      offset++; // For arbitrary phase
    }

  free(x);
  free(y);
  free(sig);
  free(val);
  free(error);
  for (i=0;i<npol;i++)
    free(covar[i]);
  free(covar);
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
      if (counter < tot+psr[p].nobs) break;
      tot+=psr[p].nobs;
    }
  ipos = counter-tot;

  new_ma=1;
  // Add global parameters
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[0].param[i].aSize;k++)
	{
	  if (psr[0].param[i].fitFlag[k] == 2)
	    nglobal++;
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
    }

  // Global fit
  int c=0;
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[p].param[i].aSize;k++)
	{
	  if (psr[p].param[i].fitFlag[k] == 2)  
	    {
	      //	      afunc[c] = dotproduct(psr[p].posPulsar,psr[p].obsn[ipos].planet_ssb[4]);
	      afunc[c] = getParamDeriv(&psr[p],ipos,x,i,k);
	      c++;
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
