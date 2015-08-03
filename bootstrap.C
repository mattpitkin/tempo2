#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
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
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"

/* Routines to calculate the fitted parameter uncertainties using a Monte-Carlo bootstrap 
 * method.  These routines are based on the bootmc.f tempo1 algorithms
 */

#define MAX_ITER 4096

double random(long *idum);

/* Needs to be called without any commented out lines */

int bootstrap(pulsar *psr,int p,int npsr)
{
  longdouble param[MAX_PARAMS],err[MAX_PARAMS],psq[MAX_PARAMS],xmean[MAX_PARAMS];
  longdouble result[MAX_PARAMS][MAX_ITER];
  longdouble fac,x1,x2,xmid,sum,sumwt,wgt,x,dt,mean,meansq,sdev;
  int nFit=0,nFit2,npts,okay;
  int i,j,k,ii,nboot,iter,l;
  int il1,il2,ih1,ih2;
  double globalParam;
  long idum = -999;              /* Should be set be clock, or user */
  const char *CVS_verNum = "$Revision: 1.6 $";

  if (displayCVSversion == 1) CVSdisplayVersion("bootstrap.C","bootstrap()",CVS_verNum);

  printf("Bootstrap1 = %d\n",psr[0].bootStrap);
  copyPSR(psr,p,npsr);           /* Have a copy of the pulsar */
  for (i=0;i<MAX_PARAMS;i++)
    copyParam(psr[0].param[i],&(psr[npsr].param[i]));

  printf("Bootstrap = %d %d\n",psr[0].bootStrap,psr[1].bootStrap);
  nboot = (int)pow(2,psr[p].bootStrap);

  for (i=0;i<psr[p].nobs;i++)
    psr[p].obsn[i].residual = psr[p].obsn[i].prefitResidual;

  /* Store the current post-fit parameters and their errors */
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[p].param[i].aSize;k++)
	{
	  if (psr[p].param[i].fitFlag[k] == 1)
	    {
	      param[nFit] = psr[p].param[i].val[k]; /* - psr[p].param[i].prefit[k]; */
	      printf("Initial param = %s %Lf %Lf\n",psr[p].param[i].label[0],
		     psr[p].param[i].val[k], psr[p].param[i].prefit[k]);
	      err[nFit]   = psr[p].param[i].err[k];
	      psq[nFit]   = 0.0;
	      nFit++;
	      psr[p].param[i].val[k] = psr[p].param[i].prefit[k];
	    }
	}
    }

  /* Determine number of TOAs */
  npts=0;
  okay=0;
  for (i=0;i<psr[p].nobs;i++)
    {
      if (psr[p].obsn[i].deleted==0)
	okay=1;
      /* Check for START and FINISH flags */
      if (psr[p].param[param_start].paramSet[0]==1 && psr[p].param[param_start].fitFlag[0]==1 &&
	  psr[p].param[param_start].val[0] > psr[p].obsn[i].bat)
	okay=0;
      if (psr[p].param[param_finish].paramSet[0]==1 && psr[p].param[param_finish].fitFlag[0]==1 &&
	  psr[p].param[param_finish].val[0] < psr[p].obsn[i].bat)
	okay=0;
      
      if (okay==1)
	npts++;
    }

  
  /* Do the bootstrap monte-carlo */
  fac  = sqrt((double)npts);
  x1   = 0.342*nboot;
  x2   = 0.477*nboot;
  xmid = 0.5*(nboot+1); 
  il1  = (int)((xmid-x1)+0.5);
  il2  = (int)((xmid-x2)+0.5);
  ih1  = (int)((xmid+x1)+0.5);
  ih2  = (int)((xmid+x2)+0.5);

  for (iter=0;iter<nboot;iter++)
    {
      sum   = 0.0;
      sumwt = 0.0;
      for (j=0;j<nFit;j++)
	xmean[j] = 0.0;

      for (i=0;i<npts;i++)
	{
	  if (psr[npsr].fitMode==1) 
	    wgt = 1.0 /
	      (1.0e-6*psr[npsr].obsn[i].toaErr*psr[npsr].param[param_f].val[0]*
	       1.0e-6*psr[npsr].obsn[i].toaErr*psr[npsr].param[param_f].val[0]);
	  else wgt=1.0/(1.0e-6*psr[npsr].param[param_f].val[0]*1.0e-6*psr[npsr].param[param_f].val[0]);

	  dt = psr[npsr].obsn[i].residual;
	  
	  ii = (int)(npts*random(&idum));  /* Randomise the data index */
	  for (j=0;j<nFit;j++)
	    xmean[j]+=wgt;  /* *fctn[j]; --- NEEDS TO BE IN -- WHAT IS THIS FOR ANYWAY?? */
	  sum+=wgt*dt;
	  sumwt+=wgt;	 	  

	  /*                                              */
	  /* Randomly change around the observation order */
	  /*                                              */
	  psr[p].obsn[i].prefitResidual  = psr[npsr].obsn[ii].prefitResidual;
	  psr[p].obsn[i].residual        = psr[npsr].obsn[ii].residual;
	  psr[p].obsn[i].sat             = psr[npsr].obsn[ii].sat;
	  psr[p].obsn[i].bat             = psr[npsr].obsn[ii].bat;
	  psr[p].obsn[i].deleted         = psr[npsr].obsn[ii].deleted;
	  psr[p].obsn[i].freq            = psr[npsr].obsn[ii].freq;
	  psr[p].obsn[i].freqSSB         = psr[npsr].obsn[ii].freqSSB; 
	  psr[p].obsn[i].toaErr          = psr[npsr].obsn[ii].toaErr;
	  strcpy(psr[p].obsn[i].fname,psr[npsr].obsn[ii].fname);
	  strcpy(psr[p].obsn[i].telID,psr[npsr].obsn[ii].telID);	  
	  for (l=0;l<3;l++)
	    {
	      psr[p].obsn[i].earth_ssb[l] = psr[npsr].obsn[ii].earth_ssb[l];
	      psr[p].obsn[i].observatory_earth[l] = psr[npsr].obsn[ii].observatory_earth[l];
	    }
	}
      writeTim("testout.tim",psr,"fred");
      psr[p].bootStrap = 0;
      doFit(&psr[p],1,0);
      /*   textOutput(psr,npsr,globalParam,0,0,0,""); */ /* Output results to the screen */

      j=0;
      
      for (i=0;i<MAX_PARAMS;i++)
	{
	  for (k=0;k<psr[p].param[i].aSize;k++)
	    {
	      if (psr[p].param[i].fitFlag[k] == 1)
		{
		  /*		  x = fac*((psr[p].param[i].val[k] - psr[p].param[i].prefit[k])-param[j])/err[j]; */
		  /* WHY IS FACTOR USED HERE? */
		  x = ((psr[p].param[i].val[k] - psr[p].param[i].prefit[k])-param[j]);
		  result[j][iter] = psr[p].param[i].val[k]-param[j];
		  psq[j]+=x*x;
		  j++;
		}
	    }
	}
      /* Store the current post-fit parameters and their errors */
      for (i=0;i<psr[p].nobs;i++)
	psr[p].obsn[i].residual = psr[npsr].obsn[i].prefitResidual;

      for (i=0;i<MAX_PARAMS;i++)
	{
	  for (k=0;k<psr[p].param[i].aSize;k++)
	    {
	      if (psr[p].param[i].fitFlag[k] == 1)
		{
		  psr[p].param[i].prefit[k]   = psr[npsr].param[i].prefit[k];
		  psr[p].param[i].val[k]      = psr[npsr].param[i].prefit[k]; 
		}
	    }
	}
      printf("Finished iteration %d of %d, so %g percent done.\n", (int)(iter+1.5),(int)(nboot+0.5),(double)((double)(iter+1.5)*100/nboot));
    }
  
  /* Restore the post-fit parameters, but use the Monte-carlo error estimates */
  nFit2=0;
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[p].param[i].aSize;k++)
	{
	  if (psr[p].param[i].fitFlag[k] == 1)
	    {
	      mean=0.0L;
	      meansq=0.0L;
	      for (l=0;l<nboot;l++)
		{
		  mean  += (result[nFit2][l])/nboot;
		  meansq+= (result[nFit2][l]*result[nFit2][l])/nboot;
		  printf("bootstrap parameters [%s] = %.14Lg\n",psr[p].param[i].shortlabel[k],
			 result[nFit2][l]+param[nFit2]);
		}
	      /*	      mean/=(longdouble)nboot;
			      meansq/=(longdouble)nboot; */
	      sdev = (longdouble)sqrt(meansq-mean*mean);
	      printf("Bootstrap mean difference: %Lf mean squared: %.14Lf rms: %.14Lf; sigma: %.14Lf, mean value: %.14Lf\n",
		     mean,(mean*mean),meansq,sdev,(mean+param[nFit2]));	      
	      /* psr[p].param[i].val[k] = psr[npsr].param[i].val[k]; */
	      psr[p].param[i].val[k] = mean+param[nFit2];
	      psr[p].param[i].err[k] = sdev;

	      /* psr[p].param[i].err[k] = err[nFit2]*sqrt(psq[nFit2]/nboot); */
	      /* psr[p].param[i].err[k] = sqrt(psq[nFit2]/(nboot-1));        */
	      /*	      sort(nboot,a[1][j]);
			      fl1[nFit2] = a[il1][j];
			      fl2[nFit2] = a[il2][j];
			      fh1[nFit2] = a[ih1][j];
			      fh2[nFit2] = a[ih2][j]; */
	      nFit2++;
	    }
	}
    }
  

  return 0;
}

/* Based on ran1.f in original fortran */

double random(long *idum)
{
  int j;
  static longdouble r[100],result;
  long m1=259100;
  long ia1=7141;
  long ic1=54773;
  long m2=134456;
  long ia2=8121;
  long ic2=28411;
  long m3=243000;
  long ia3=4561;
  long ic3=51349;
  longdouble rm1,rm2;

  static int iff=0;
  static int ix1=0;
  static int ix2=0;
  static int ix3=0;

  rm1=1./m1;
  rm2=1./m2;
    
  if(*idum < 0 || iff == 0) 
    {
      iff=1;
      ix1=(int)fortran_mod(ic1-(*idum),m1);
      ix1=(int)fortran_mod(ia1*ix1+ic1,m1);
      ix2=(int)fortran_mod(ix1,m2);
      ix1=(int)fortran_mod(ia1*ix1+ic1,m1);
      ix3=(int)fortran_mod(ix1,m3);

      for (j=0;j<97;j++)
	{
	  ix1=(int)fortran_mod(ia1*ix1+ic1,m1);
	  ix2=(int)fortran_mod(ia2*ix2+ic2,m2);
	  r[j]=(ix1+ix2*rm2)*rm1;
	}
      *idum=1;
    }
  
  ix1=(int)fortran_mod(ia1*ix1+ic1,m1);
  ix2=(int)fortran_mod(ia2*ix2+ic2,m2);
  ix3=(int)fortran_mod(ia3*ix3+ic3,m3);
  j=(97*ix3)/m3;
  if(j>96 || j<0) {printf("Problem in bootstrap.C (%d)\n",j); exit(1);}

  result = r[j];
  r[j]=(ix1+ix2*rm2)*rm1;
  return result;
}
