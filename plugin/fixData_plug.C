#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russel Edwards, David Champion

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
#include <cpgplot.h>
#include "T2toolkit.h"
#include "TKspectrum.h"
#include "TKfit.h"
#include "tempo2.h"
#include "GWsim.h"

using namespace std;

void doPlugin1(pulsar *psr,char *flag,int removeQuad);
void doPlugin3(pulsar *psr,char *flag,char parFile[MAX_PSR_VAL][MAX_FILELEN],char timFile[MAX_PSR_VAL][MAX_FILELEN],int argc,char *argv[],float dstep);
int determine1dStructureFunction(float *x,float *y,float *ye,int nn,double *errfac1,
				 double *vsf,double *mverr);
void doPlugin2(pulsar *psr,char parFile[MAX_PSR_VAL][MAX_FILELEN],char timFile[MAX_PSR_VAL][MAX_FILELEN],int argc,char *argv[]);
float plotHistogram(float *x,int count,int *flagCol,int nFlag,char flagV[100][16]);
void doSummary(pulsar *psr,float errStep);

int nit = 1;
long double gwamp = 0;
long double alpha = 1;
char plotout[20]="/xs";
int plotoutSet = 0;
int script = 0;
int dayGap=1;

void help() /* Display help */
{
  printf("\n\n");
  printf("This plugin divides up a data set based on the different flags to indicate different\nbackend systems\n\n");
  printf("usage:\n");
  printf("\t-flag <flag>      defines which flag identifies different backends\n");
  printf("\t-plot <plottype>  = 1 to caclulate EFACs, = 2 to check whiteness\n");
  printf("For use with plot type 1:\n");
  printf("\t-removeQuad       remove quadratic from all the residuals with a given flag\n");
  printf("\t-daygap <int>     number of days to determine white level (default = 1d)\n");
  printf("For use with plot type 2:\n");
  printf("\t-numits <int>     number if realisations for noise simulation\n");
  printf("\t-gwamp <ldouble>  amplitude of GW noise source\n");
  printf("\t-alpha <ldouble>  power spectrum of noise (1=white)\n");
  printf("\t-plotout <grdev>  pgplot device for final plot\n");
  printf("\t-script           useful when running noninteractively - does not show any plots\n");
}


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  char flag[10]="-f";
  int i,summary;
  int plot=1;
  int removeQuad=0;
  double globalParameter;
  float errStep = 1e-6;
  float dstep = 100;
  //  int nit = 1;
  //  long double gwamp = 0;
  //  long double alpha = 1;

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */
  summary=0;

  printf("Graphical Interface: fixData\n");
  printf("Author:              G. Hobbs, D. Champion\n");
  printf("CVS Version:         $Revision $\n");

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[0],argv[i+1]); 
	  strcpy(timFile[0],argv[i+2]);
	}
      else if (strcmp(argv[i],"-flag")==0)
	strcpy(flag,argv[++i]);
      else if (strcmp(argv[i],"-removeQuad")==0)
	removeQuad=1;
      else if (strcmp(argv[i],"-plot")==0)
	sscanf(argv[++i],"%d",&plot);
      else if (strcmp(argv[i],"-h")==0)
	{
	  help();
	  exit(0);
	}
      else if (strcmp(argv[i],"-numits")==0)
	sscanf(argv[++i],"%d",&nit);
      else if (strcmp(argv[i],"-gwamp")==0)
	sscanf(argv[++i],"%Lf",&gwamp);
      else if (strcmp(argv[i],"-alpha")==0)
	sscanf(argv[++i],"%Lf",&alpha);
      else if (strcmp(argv[i],"-errstep")==0)
	sscanf(argv[++i],"%f",&errStep);
      else if (strcmp(argv[i],"-plotout")==0)
      {
	plotoutSet = 1;
	strcpy(plotout,argv[++i]);
      }
      else if (strcmp(argv[i],"-daygap")==0)
	sscanf(argv[++i],"%d",&dayGap);
      else if (strcmp(argv[i],"-dstep")==0)
		sscanf(argv[++i],"%f",&dstep);
      else if (strcmp(argv[i],"-script")==0)
	script = 1;
      else if (strcmp(argv[i],"-summary")==0)
	summary=1;
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

  if (summary==1)
    {
      doSummary(psr,errStep);
      exit(1);
    }
  if (plot==1)
    doPlugin1(psr,flag,removeQuad);
  else if (plot==2)
    doPlugin2(psr,parFile,timFile,argc,argv);
  else if (plot==3)
    doPlugin3(psr,flag,parFile,timFile,argc,argv,dstep);
  else
    printf("Unknown plot required\n");

  return 0;
}

void doPlugin3(pulsar *psr,char *flag,char parFile[MAX_PSR_VAL][MAX_FILELEN],char timFile[MAX_PSR_VAL][MAX_FILELEN],int argc,char *argv[],float dstep)
{
  int n = psr[0].nobs;
  int i,j,k;
  float resx[n],resy[n],rese[n],ye1[n],ye2[n];
  float wresx[n],wresy[n],wrese[n];
  float histX[100],histY[100];
  int nhist;
  float minResX,maxResX,minResY,maxResY;
  float actRMS;
  float tx[2],ty[2];
  float normY[n],wnormY[n];
  double xpos;
  double globalParameter=0;
  double reducedChisq;
  double equadv=0,equad=0;
  double efac=1,efacv=0;
  float highest;
  int loop=1;
  int nFlag=0,found;
  int col;
  char flagV[100][16];
  int flagCol[n],flagID[n];
  int nCount;
  int use_equad=1;
  double storeOrigErr[n];
  double px[n],py[n],specX[2*n],specY[2*n],pe[n],temp[2*n];
  float specX_f[2*n],specY_f[2*n];
  int np;
  double dx;
  int nSpec;

  printf("Use equad (1), efac (2) or both (3) ");
  scanf("%d",&use_equad);

  // b) Must produce table giving the variance of each backend normalised residuals

  for (i=0;i<psr[0].nobs;i++)
    {
      storeOrigErr[i]=psr[0].obsn[i].toaErr;
      for (j=0;j<psr[0].obsn[i].nFlags;j++)
	{
	  if (strcmp(psr[0].obsn[i].flagID[j],flag)==0)
	    {
	      found=0;
	      flagCol[i] = 1;
	      flagID[i] = 1;
	      for (k=0;k<nFlag;k++)
		{
		  if (strcmp(psr[0].obsn[i].flagVal[j],flagV[k])==0)
		    {
		      if (k+2 > 14)
			flagCol[i] = 14;
		      else
			flagCol[i] = k+2;
		      flagID[i] = k;
		      found=1;
		      break;
		    }
		}
	      if (found==0)
		{
		  strcpy(flagV[nFlag++],psr[0].obsn[i].flagVal[j]);
		  if (nFlag+1 > 14)
		    flagCol[i]=14;
		  else
		    flagCol[i] = nFlag+1;
		  flagID[i] = nFlag;
		}
	    }
	}
    }
  printf("Number of unique flags = %d\n",nFlag);
  for (k=0;k<nFlag;k++)
    printf(" -- %s\n",flagV[k]);


  for (i=0;i<n;i++)
    {
      resx[i] = (float)(psr[0].obsn[i].sat - psr[0].param[param_pepoch].val[0]);
      resy[i] = (float)(psr[0].obsn[i].residual);
      rese[i] = (float)(psr[0].obsn[i].toaErr*1.0e-6);
      ye1[i] = resy[i]-rese[i];
      ye2[i] = resy[i]+rese[i];
      normY[i] = resy[i]/rese[i];
    }
  minResX = TKfindMin_f(resx,n);
  maxResX = TKfindMax_f(resx,n);
  minResY = TKfindMin_f(resy,n);
  maxResY = TKfindMax_f(resy,n);

  // Plot original residuals
  cpgbeg(0,"23/xs",1,2);
  cpgsch(1.3);
  cpgenv(minResX,maxResX,minResY,maxResY,0,1);
  cpglab("Day","Residual","");
  for (i=0;i<n;i++)
    {
      cpgsci(flagCol[i]);
      cpgpt(1,resx+i,resy+i,4);
      cpgerry(1,resx+i,ye1+i,ye2+i,1);
      cpgsci(1);
    }
  minResY = TKfindMin_f(normY,n);
  maxResY = TKfindMax_f(normY,n);
  cpgenv(minResX,maxResX,minResY,maxResY,0,1);
  cpglab("Day","Normalised Residual","");
  for (i=0;i<n;i++)
    {
      cpgsci(flagCol[i]);
      cpgpt(1,resx+i,normY+i,4);
      cpgsci(1);
    }
  // Draw 1 sigma and 3 sigma lines
  tx[0] = minResX; tx[1] = maxResX;
  ty[0] = ty[1] = 1.0;  cpgsci(7); cpgsls(4); cpgline(2,tx,ty); cpgsls(1); cpgsci(1);
  ty[0] = ty[1] = -1.0;  cpgsci(7); cpgsls(4); cpgline(2,tx,ty); cpgsls(1); cpgsci(1);
  ty[0] = ty[1] = 3.0;  cpgsci(3); cpgsls(4); cpgline(2,tx,ty); cpgsls(1); cpgsci(1);
  ty[0] = ty[1] = -3.0;  cpgsci(3); cpgsls(4); cpgline(2,tx,ty); cpgsls(1); cpgsci(1);

  // Draw actual rms
  actRMS = TKfindRMS_f(normY,n);
  ty[0] = ty[1] = actRMS;  cpgsci(2);  cpgline(2,tx,ty); cpgsci(1);
  ty[0] = ty[1] = -actRMS;  cpgsci(2);  cpgline(2,tx,ty); cpgsci(1);
  printf("RMS of normalised residuals = %g\n",actRMS);
  cpgend();

  // Re-do the fit with constrained fit
  // 1. Turn off all fitting
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (j=0;j<psr[0].param[i].aSize;j++)
	{
	  if (psr[0].param[i].fitFlag[j]==1)
	    psr[0].param[i].fitFlag[j]=0;
	}
    }
  for (i=1;i<=psr[0].nJumps;i++)
    {
      if (psr[0].fitJump[i]==1)
	psr[0].fitJump[i]=0;
    }

  // 2. Turn on fitting for F0 and F1
  psr[0].param[param_f].fitFlag[0] = 1;
  psr[0].param[param_f].fitFlag[1] = 1;
  
  // 3. Add in constrained IFUNC values
  psr[0].param[param_ifunc].paramSet[0]=1;
  psr[0].param[param_ifunc].fitFlag[0]=1;
  psr[0].param[param_ifunc].val[0]=2;


  do {
    xpos = minResX + (double)psr[0].param[param_pepoch].val[0] - dstep/2;

    psr[0].param[param_f].val[0] = psr[0].param[param_f].prefit[0];
    psr[0].param[param_f].val[1] = psr[0].param[param_f].prefit[1];

    for (i=0;i<psr[0].nobs;i++)
      {
	// Check for T2efacs
	psr[0].obsn[i].toaErr = storeOrigErr[i];
	for (j=0;j<psr[0].obsn[i].nFlags;j++)
	  {	 
	    for (k=0;k<psr[0].nT2efac;k++)
	      {
		if (strcmp(psr[0].obsn[i].flagID[j],psr[0].T2efacFlagID[k])==0)
		  {
		    if (strcmp(psr[0].obsn[i].flagVal[j],psr[0].T2efacFlagVal[k])==0)
		      psr[0].obsn[i].toaErr *= psr[0].T2efacVal[k];
		  }
	      }
	  }
	if (use_equad==1 || use_equad==3)
	  {
	    psr[0].obsn[i].equad = equad;
	    psr[0].obsn[i].toaErr = sqrt(pow(psr[0].obsn[i].toaErr,2)+equad*equad);
	  }
	if (use_equad==2 || use_equad==3)
	  {
	    psr[0].obsn[i].efac = efac;
	    psr[0].obsn[i].toaErr *= efac;
	  }
	// Must do this because the DM correction error bar uses the original error
	psr[0].obsn[i].origErr = psr[0].obsn[i].toaErr;
      }

    psr[0].ifuncN = 0;
    
    do {
      if (psr[0].ifuncN > 0)
	{
	  nCount=0;
	  for (k=0;k<psr[0].nobs;k++)
	    {
	      if (psr[0].obsn[k].sat > psr[0].ifuncT[psr[0].ifuncN-1] &&
		  psr[0].obsn[k].sat < xpos)
		nCount++;
	    }
	  if (nCount>2)
	    {
	      psr[0].ifuncT[psr[0].ifuncN] = xpos;
	      psr[0].ifuncV[psr[0].ifuncN] = 0.0;
	      psr[0].ifuncE[psr[0].ifuncN] = 0.0;
	      psr[0].ifuncN++;
	    }
	}
      else
	{
	  psr[0].ifuncT[psr[0].ifuncN] = xpos;
	  psr[0].ifuncV[psr[0].ifuncN] = 0.0;
	  psr[0].ifuncE[psr[0].ifuncN] = 0.0;
	  psr[0].ifuncN++;
	}
      xpos += dstep;
    } while (xpos < maxResX + (double)psr[0].param[param_pepoch].val[0]);
    psr[0].constraints[psr[0].nconstraints++]= constraint_ifunc_0;
    psr[0].constraints[psr[0].nconstraints++]= constraint_ifunc_1;
    psr[0].constraints[psr[0].nconstraints++]= constraint_ifunc_2;
    
    formBatsAll(psr,1);         /* Form the barycentric arrival times */
    formResiduals(psr,1,1);    /* Form the residuals                 */
    doFit(psr,1,0);   /* Do the fitting     */
    formBatsAll(psr,1);         /* Form the barycentric arrival times */
    formResiduals(psr,1,1);    /* Form the residuals                 */
    textOutput(psr,1,globalParameter,0,0,0,"");  /* Display the output */
    reducedChisq = psr[0].fitChisq/(double)psr[0].fitNfree;

    cpgbeg(0,"24/xs",1,2);
    cpgsch(1.3);
    for (i=0;i<psr[0].nobs;i++)
      {
	wresx[i] = (float)(psr[0].obsn[i].sat - psr[0].param[param_pepoch].val[0]);
	wresy[i] = (float)(psr[0].obsn[i].residual);
	wrese[i] = (float)(psr[0].obsn[i].toaErr*1.0e-6);
	ye1[i] = wresy[i]-wrese[i];
	ye2[i] = wresy[i]+wrese[i];
	wnormY[i] = wresy[i]/wrese[i];	
	px[i] = (double)wresx[i];
	py[i] = (double)wnormY[i];
	pe[i] = 1.0;
	if (fabs(wnormY[i]) > 4) printf("> 4 sigma: %s %g %g %g\n",psr[0].obsn[i].fname,wresy[i],wrese[i],wnormY[i]);
      }

    minResY = TKfindMin_f(wresy,n);
    maxResY = TKfindMax_f(wresy,n);

    cpgenv(minResX,maxResX,minResY,maxResY,0,1);
    cpglab("Day","Residual","");

    for (i=0;i<n;i++)
      {
	cpgsci(flagCol[i]);
	cpgpt(1,wresx+i,wresy+i,4);
	cpgerry(1,wresx+i,ye1+i,ye2+i,1);
	cpgsci(1);
      }
    minResY = TKfindMin_f(wnormY,n);
    maxResY = TKfindMax_f(wnormY,n);
    cpgenv(minResX,maxResX,minResY,maxResY,0,1);
    cpglab("Day","Normalised Residual","");
    for (i=0;i<n;i++)
      {
	cpgsci(flagCol[i]);
	cpgpt(1,wresx+i,wnormY+i,4);
	cpgsci(1);
      }
    // Draw 1 sigma and 3 sigma lines
    tx[0] = minResX; tx[1] = maxResX;
    ty[0] = ty[1] = 1.0;  cpgsci(7); cpgsls(4); cpgline(2,tx,ty); cpgsls(1); cpgsci(1);
    ty[0] = ty[1] = -1.0;  cpgsci(7); cpgsls(4); cpgline(2,tx,ty); cpgsls(1); cpgsci(1);
    ty[0] = ty[1] = 3.0;  cpgsci(3); cpgsls(4); cpgline(2,tx,ty); cpgsls(1); cpgsci(1);
    ty[0] = ty[1] = -3.0;  cpgsci(3); cpgsls(4); cpgline(2,tx,ty); cpgsls(1); cpgsci(1);
    
    // Draw actual rms
    actRMS = TKfindRMS_f(wnormY,n);
    ty[0] = ty[1] = actRMS;  cpgsci(2);  cpgline(2,tx,ty); cpgsci(1);
    ty[0] = ty[1] = -actRMS;  cpgsci(2);  cpgline(2,tx,ty); cpgsci(1);
    printf("RMS of normalised residuals = %g\n",actRMS);

    // Provide table of normalised residuals for different backends
    printf("------------------------------------------------------------------------\n");
    printf(" -- %-16.16s %4.4s %9.9s %9.9s\n","Flag","Npts","Mean","RMS");
    printf("------------------------------------------------------------------------\n");
    for (i=0;i<nFlag;i++)
      {
	np=0;
	for (j=0;j<n;j++)
	  {
	    if (flagID[j] == i)
	      {
		py[np] = wnormY[j];
		np++;
	      }
	  }
	printf(" -- %-16.16s %4d %9.3g %9.3g\n",flagV[i],np,TKmean_d(py,np),TKfindRMS_d(py,np));
      }
    printf("------------------------------------------------------------------------\n");
    cpgend();

    // Plot spectrum of whitened, normalised residuals
    cpgbeg(0,"26/xs",1,1);
    TKspectrum(px,py,pe,n,0,0,0,0,0,2,1,1,1,specX,specY,&nSpec,0,0,temp,temp);
    for (i=0;i<nSpec;i++)
      {
	specX_f[i] = (float)specX[i]*365.25;
	specY_f[i] = (float)log10(specY[i]*pow(365.25*86400.0,2));
      }
    cpgenv(0,TKfindMax_f(specX_f,nSpec),TKfindMin_f(specY_f,nSpec),TKfindMax_f(specY_f,nSpec),0,20);
    cpglab("Frequency (yr\\u-1\\d)","PSD (yr)",psr[0].name);
    cpgline(nSpec,specX_f,specY_f);
    // Draw confidence levels
    tx[0] = 0; tx[1] = TKfindMax_f(specX_f,nSpec);
    ty[0] = ty[1] = log10(1.0/specX_f[nSpec-1]);  cpgsci(7); cpgsls(2); cpgline(2,tx,ty); cpgsls(1); cpgsci(1);
    ty[0] = ty[1] = log10(3.0/specX_f[nSpec-1]);  cpgsci(3); cpgsls(2); cpgline(2,tx,ty); cpgsls(1); cpgsci(1);
    ty[0] = ty[1] = log10(0.05/specX_f[nSpec-1]);  cpgsci(3); cpgsls(2); cpgline(2,tx,ty); cpgsls(1); cpgsci(1);

    tx[0] = tx[1] = 1.0/(dstep/365.25);
    ty[0] = TKfindMin_f(specY_f,nSpec);
    ty[1] = TKfindMax_f(specY_f,nSpec);
    cpgsci(2); cpgsls(1); cpgline(2,tx,ty); cpgsls(1); cpgsci(1);
    // Now plot expected curve
    for (i=0;i<nSpec;i++)
      {
	dx = specX_f[i]*dstep/365.25;
	specY_f[i] = log10((1-pow(sin(dx)/dx,4))*1.0/specX_f[nSpec-1]);
      }
    cpgsci(5); cpgsls(4); cpgline(nSpec,specX_f,specY_f); cpgsci(1); cpgsls(1);
    cpgend();


    cpgbeg(0,"25/xs",1,1);
    highest = plotHistogram(wnormY,n,flagCol,nFlag,flagV);
    // Overplot Gaussian with sdev = 1
    for (i=0;i<100;i++)
      {

	histX[i] = (i-50)/8.0;
	histY[i] = highest*exp(-pow(histX[i],2)/2);
      }
    cpgsci(7); cpgline(100,histX,histY); cpgsci(1);
    cpgend();
    

    printf("Reduced chisq = %g\n",reducedChisq);
    if (use_equad==1 || use_equad==3)
      {
	printf("New EQUAD: (%g; -1 = quit) ",equad);
	scanf("%lf",&equadv);
      
	if (equadv==-1)
	  loop=0;
	else
	  equad = equadv;
      }
    if (use_equad==2 || use_equad==3)
      {
	printf("New EFAC: (%g; -1 = quit) ",efac);
	scanf("%lf",&efacv);
	if (efacv==-1)
	  loop=0;
	else
	  efac=efacv;
      }

    psr[0].nconstraints-=3;
  } while (loop==1); 
  if (use_equad==1) printf("Final choice of EQUAD = %g us\n",equad);
  else if (use_equad==2) printf("Final choice of EFAC = %g\n",efac);
}

void doSummary(pulsar *psr,float errStep)
{
  int i;
  float x[psr[0].nobs];
  float y[psr[0].nobs];
  float e[psr[0].nobs];
  float plot_x[psr[0].nobs];
  float plot_y[psr[0].nobs];
  float fx[2],fy[2];
  int n,nplot=0;
  float minErr=psr[0].obsn[0].toaErr*1e-6;
  float maxErr=psr[0].obsn[0].toaErr*1e-6;
  float minx,maxx,miny,maxy;
  float err;
  double errfac,vsf,mverr;
  long double sat1[psr[0].nobs];
  int it,nit;
  int npts_sf;
  long seed = TKsetSeed();

  // Find minimum error value
  for (i=0;i<psr[0].nobs;i++)
    {
      if (minErr > psr[0].obsn[i].toaErr*1.0e-6)
	minErr = psr[0].obsn[i].toaErr*1.0e-6;
      if (maxErr < psr[0].obsn[i].toaErr*1.0e-6)
	maxErr = psr[0].obsn[i].toaErr*1.0e-6;
    }
  printf("Minimum error bar = %g seconds\n",minErr);
  printf("Maximum error bar = %g seconds\n",maxErr);
  for (err = minErr;err < maxErr;err+=errStep)
    {
      n=0;
      for (i=0;i<psr[0].nobs;i++)
	{
	  if (psr[0].obsn[i].toaErr*1.0e-6 >= err &&
	      psr[0].obsn[i].toaErr*1.0e-6 < err+errStep)
	    {
	      x[n] = (float)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
	      y[n] = (float)(psr[0].obsn[i].residual);
	      e[n] = (float)(psr[0].obsn[i].toaErr*1.0e-6);
	      n++;
	    }
	}
	      printf("%g < Err < %g has %d points\n",err,err+errStep,n);

      if (n > 2)
	{
	  npts_sf = determine1dStructureFunction(x,y,e,n,&errfac,&vsf,&mverr);
	  if (npts_sf > 1)
	    {
	      printf("Have %d %g %g %g\n",npts_sf,errfac,vsf,mverr);
	  

	      plot_x[nplot] = log10(pow(err+errStep/2.0,2));
	      //	      plot_y[nplot] = log10(errfac);
	      plot_y[nplot] = log10((vsf));
	      if (nplot==0)
		{
		  minx = maxx = plot_x[nplot];
		  miny = maxy = plot_y[nplot];
		}
	      else
		{
		  if (minx > plot_x[nplot]) minx = plot_x[nplot];
		  if (maxx < plot_x[nplot]) maxx = plot_x[nplot];
		  if (miny > plot_y[nplot]) miny = plot_y[nplot];
		  if (maxy < plot_y[nplot]) maxy = plot_y[nplot];
		}
	      nplot++;
	    }
	}
    }
  cpgbeg(0,"/xs",1,1);
  cpgsch(1.4);
  minx -= 0.1;
  maxx += 0.1;
  miny -= 2;
  maxy += 0.1;

  cpgenv(minx,maxx,miny,maxy,0,0);
  cpglab("log10(Error bar size\\u2\\d)","log10(var)","");
  cpgpt(nplot,plot_x,plot_y,5);
  fx[0] = minx; fx[1] = maxx;
  fy[0] = fx[0];
  fy[1] = fx[1]; 
  cpgsls(4); cpgsci(7); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
  fy[0] = log10(pow(10,fx[0])*pow(2,2)); // EFAC = 2
  fy[1] = log10(pow(10,fx[1])*pow(2,2)); 
  cpgsls(4); cpgsci(2); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
  fy[0] = log10(pow(10,fx[0])*pow(4,2)); // EFAC = 4
  fy[1] = log10(pow(10,fx[1])*pow(4,2)); 
  cpgsls(4); cpgsci(3); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
  printf("errfac = %g\n",errfac);

  // Calculate errors
  for (i=0;i<psr[0].nobs;i++)
    sat1[i] = psr[0].obsn[i].sat - psr[0].obsn[i].residual/86400.0L;
  nit = 100;
  for (it=0;it<nit;it++)
    {
      nplot=0;
      for (i=0;i<psr[0].nobs;i++)
	psr[0].obsn[i].sat = sat1[i] + psr[0].obsn[i].toaErr*1.0e-6*TKgaussDev(&seed)/86400.0L;
      formBatsAll(psr,1);         /* Form the barycentric arrival times */
      formResiduals(psr,1,1);    /* Form the residuals                 */

      // Now do the process again
      for (err = minErr;err < maxErr;err+=errStep)
	{
	  n=0;
	  for (i=0;i<psr[0].nobs;i++)
	    {
	      if (psr[0].obsn[i].toaErr*1.0e-6 >= err &&
		  psr[0].obsn[i].toaErr*1.0e-6 < err+errStep)
		{
		  x[n] = (float)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
		  y[n] = (float)(psr[0].obsn[i].residual);
		  e[n] = (float)(psr[0].obsn[i].toaErr*1.0e-6);
		  n++;
		}
	    }
	  //	  printf("%g < Err < %g has %d points (part 2)\n",err,err+errStep,n);
	  
	  if (n > 2)
	    {
	      npts_sf = determine1dStructureFunction(x,y,e,n,&errfac,&vsf,&mverr);
	      if (npts_sf > 1)
		{
		  //		  printf("Have part 2 %d %g %g %g\n",npts_sf,errfac,vsf,mverr);
		  
		  plot_x[nplot] = log10(pow(err+errStep/2.0,2));
		  //	      plot_y[nplot] = log10(errfac);
		  plot_y[nplot] = log10((vsf));
		  nplot++;
		}
	    }
	  cpgpt(nplot,plot_x,plot_y,1);
	}     
    }
  

  cpgend();
}

// Check whether the data are white
void doPlugin2(pulsar *psr,char parFile[MAX_PSR_VAL][MAX_FILELEN],char timFile[MAX_PSR_VAL][MAX_FILELEN],int argc,char *argv[])
{
  long double sat0[MAX_OBSN];
  int i,j,p;
  long seed = -23;
  int it;
  int count=0;
  double measuredRMS;
  float rms[nit];
  double actSpecX[MAX_OBSN],actSpecY[MAX_OBSN];
  double specX[MAX_OBSN];
  double ix[MAX_OBSN],iy[MAX_OBSN],ie[MAX_OBSN],outY_re[MAX_OBSN],outY_im[MAX_OBSN];
  double **specY;
  float fx[MAX_OBSN],fy[MAX_OBSN],t95[MAX_OBSN],t5[MAX_OBSN];
  int specN;
  // GW
  long double a;
  long double toffset;
  long double kp[3];
  long double flo,fhi;
  long double res[MAX_OBSN],mean;
  double dist;
  int addGW=0,ngw,k;
  gwSrc *gw;
  //  long double gwamp = 18e-14;

  //alpha = 0.8; //2.5/2.0;
  a = (long double)gwamp*pow(86400.0*365.25,alpha);
  dist =  3.08568025e19; // 1 kpc in m
  setupPulsar_GWsim(psr[0].param[param_raj].val[0],
		    psr[0].param[param_decj].val[0],kp);
  flo = 1.0L/(30*365.25*86400.0L);
  fhi = 1.0L/(2.0*86400.0L);
  ngw = 1000;
  if((gw = (gwSrc *)malloc(sizeof(gwSrc)*ngw))==NULL){
    printf("Unable to allocate memory for gwSrc.\n");
    exit(1);
  }
  toffset = psr[0].param[param_pepoch].val[0];

  specY = (double **)malloc(sizeof(double *)*nit);
  for (i=0;i<nit;++i)
    specY[i] = (double *)malloc(sizeof(double)*MAX_OBSN);
  // Obtain spectrum
  for (i=0;i<psr[0].nobs;i++)
    {
      ix[i] = (double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
      iy[i] = (double)(psr[0].obsn[i].residual);
      ie[i] = (double)psr[0].obsn[i].toaErr*1.0e-6;
    }
  printf("Producing spectrum\n");
  TKspectrum(ix,iy,ie,psr[0].nobs,0,0,0,0,0,2,1,1,1,actSpecX,actSpecY,&specN,0,0,outY_re,outY_im);
  printf("Complete producing spectrum\n");

  if (script == 0) 
    cpgbeg(0,"/xs",1,1);
  else 
    cpgbeg(0,"/null",1,1);
  cpgask(0);
  // Obtain idealised TOAs
  for (j=0;j<5;j++)
    {
      psr[0].nJumps = 0;
      for(i=0;i<MAX_PARAMS;i++){
	psr[0].param[i].nLinkTo = 0;
	psr[0].param[i].nLinkFrom = 0;
      }
      readParfile(psr,parFile,timFile,1); /* Load the parameters       */
      formBatsAll(psr,1);         /* Form the barycentric arrival times */
      formResiduals(psr,1,0);    /* Form the residuals                 */
      for (i=0;i<psr[0].nobs;i++)
	psr[0].obsn[i].sat -= (long double)psr[0].obsn[i].residual/86400.0L;
    }
  for (i=0;i<psr[0].nobs;i++)
    sat0[i] = psr[0].obsn[i].sat;


  measuredRMS = (double)psr[0].rmsPost;

  for (it=0;it<nit;it++)
    {
      psr[0].nJumps = 0;
      for(i=0;i<MAX_PARAMS;i++){
	psr[0].param[i].nLinkTo = 0;
	psr[0].param[i].nLinkFrom = 0;
      }
      for (i=0;i<psr[0].nobs;i++)
	psr[0].obsn[i].sat = sat0[i]+(TKgaussDev(&seed)*psr[0].obsn[i].toaErr*1.0e-6)/86400.0;

      // Add GW
      GWbackground(gw,ngw,&seed,flo,fhi,a,alpha,1);
      mean=0.0L;
      //printf("Calc residuals 1\n");
      for (j=0;j<psr[0].nobs;j++)
	{
	  res[j]=0.0L;
	  for (k=0;k<ngw;k++)
	    res[j]+=calculateResidualGW(kp,&gw[k],
					(psr[0].obsn[j].sat-toffset)*86400.0L,
					dist);	  
	  mean+=res[j];
	}
      for (j=0;j<psr[0].nobs;j++)
	{
	  psr[0].obsn[j].sat+=(res[j]-mean/psr[0].nobs)/86400.0L;
	  if (it==0)
	    printf("gwres %g %g\n",(double)psr[0].obsn[j].sat,(double)res[j]);
	}
      if (it==0)
	writeTim("sim.tim",psr,"tempo2");


      readParfile(psr,parFile,timFile,1); /* Load the parameters                */
      preProcess(psr,1,argc,argv);
      formBatsAll(psr,1);                 /* Form the barycentric arrival times */
      formResiduals(psr,1,0);             /* Form the residuals                 */
      doFit(psr,1,0);                     /* Do the fitting                     */ 
      formBatsAll(psr,1);                 /* Form the barycentric arrival times */
      formResiduals(psr,1,0);             /* Form the residuals                 */
      textOutput(psr,1,0,0,0,0,"");
      rms[it] = (float)psr[0].rmsPost;
      if (rms[it] > measuredRMS) count++;
      //      plotHistogram(rms,it+1);
      printf("[%d/%d] rms = %g, measured rms = %g. percent rms > measured value = %g\n",it+1,nit,rms[it],measuredRMS,count/(double)(it+1)*100.0);
      for (i=0;i<psr[0].nobs;i++)
	{
	  ix[i] = (double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
	  iy[i] = (double)(psr[0].obsn[i].residual);
	  ie[i] = (double)psr[0].obsn[i].toaErr*1.0e-6;
	}
      printf("Making spectrum\n");
      TKspectrum(ix,iy,ie,psr[0].nobs,0,0,0,0,0,2,1,1,1,specX,specY[it],&specN,0,0,outY_re,outY_im);
      printf("%d Complete making spectrum\n",it);
    }
  cpgend();
 
  if (script == 0 || plotoutSet == 1) 
    cpgbeg(0,plotout,1,2);
  else
    cpgbeg(0,"/null",1,2);
  cpgsch(1.4);
  //  plotHistogram(rms,nit);
  fx[0] = fx[1] = (float)measuredRMS;
  fy[0] = 0; fy[1] = 1000;
  cpgsls(4); cpgline(2,fx,fy); cpgsls(1);
  // Plot spectrum
  {
    float avY[MAX_OBSN];
    double temp[MAX_OBSN];
    for (i=0;i<specN;i++)
      {
	fx[i] = log10(actSpecX[i]);
	fy[i] = log10(actSpecY[i]);
      }
    cpgenv(TKfindMin_f(fx,specN),TKfindMax_f(fx,specN),TKfindMin_f(fy,specN)-2,TKfindMax_f(fy,specN)+2,0,30);
    cpglab("Frequency (d\\u-1\\d)","PSD","");
    cpgline(specN,fx,fy);
    for (i=0;i<specN;i++)
      {
	avY[i] = 0.0;
	for (j=0;j<nit;j++)
	  avY[i] += specY[j][i];
	avY[i] =  log10(avY[i]/(double)nit);
      }
    cpgsci(2); cpgline(specN,fx,avY); cpgsci(1);
    // Plot 95% confidence levels
    printf("plotting 95 percent levels\n");
    for (i=0;i<specN;i++)
      {
	for (j=0;j<nit;j++)
	  temp[j] = specY[j][i];
	TKsort_d(temp,nit);
	t95[i] = (float)log10(temp[(int)(95.0/100.0*nit+0.5)]);
	t5[i] = (float)log10(temp[(int)(5.0/100.0*nit+0.5)]);
      }
    cpgsls(4); cpgsci(2); cpgline(specN,fx,t95);
    cpgline(specN,fx,t5); cpgsci(1); cpgsls(1);
  }
  cpgend();

  for (i=0;i<nit;++i)
    free(specY[i]);
  free(specY);
}

void doPlugin1(pulsar *psr,char *flag,int removeQuad)
{
  int nflag=0;
  char flagV[MAX_JUMPS][100];
  double errfacs[MAX_JUMPS];
  float fx[MAX_OBSN],fy[MAX_OBSN],fe1[MAX_OBSN],fe2[MAX_OBSN],ye[MAX_OBSN],y[MAX_OBSN];
  int id[MAX_OBSN];
  double vsf,mverr;
  float yrange;
  float minx,miny;
  float maxx,maxy;
  int i,j,k,n;
  int found;
  float valx,valy;
  float tt;
  double errfac;
  float span;
  char xlabel[100];
  char text[100];
  int nefacFlag=0;
  int npts_sf;

  // Get range
  for (i=0;i<psr[0].nobs;i++)
    {
      valx = (float)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
      valy = (float)(psr[0].obsn[i].residual);
      if (i==0)
	{
	  minx = maxx = valx;	  
	  maxy = fabs(valy);
	}
      else
	{
	  if (minx > valx) minx = valx;
	  if (maxx < valx) maxx = valx;
	  if (maxy < fabs(valy)) maxy = fabs(valy);
	}
    }
  yrange = 2*maxy;
  tt = maxx-minx;
  minx -= tt*0.05;
  maxx += tt*0.05;
  // Count flags
  for (i=0;i<psr[0].nobs;++i)
    {
      if (psr[0].obsn[i].deleted==0)
	{
	  for (j=0;j<psr[0].obsn[i].nFlags;++j)
	    {
	      if (strcmp(psr[0].obsn[i].flagID[j],flag)==0)
		{
		  found=0;
		  for (k=0;k<nflag;k++) 
		    {
		      if (strcmp(flagV[k],psr[0].obsn[i].flagVal[j])==0)
			{
			  found=1;
			  break;
			}
		    }
		  if (found==0)
		    strcpy(flagV[nflag++],psr[0].obsn[i].flagVal[j]);
		}
	    }
	}
    }
  if (script == 0) 
    cpgbeg(0,"/xs",1,1);
  else 
    cpgbeg(0,"/null",1,1);
  cpgsch(1.4);
  cpgask(0);
  cpgenv(minx,maxx,0,nflag,0,-1);
  cpgbox("ATNSBC",0.0,0,"",0.0,0);
  sprintf(xlabel,"MJD - %.2f",(float)psr[0].param[param_pepoch].val[0]);
  cpglab(xlabel,"","");
  printf("------------------------------------------------------------------------------\n");
  printf("%-15s %5s %6s %7s %8s %8s %8s %8s %5s\n","Data set","Npts","Span","Npts sf","W. Var","M. err^2","W. RMS","M. err","EFAC");
  printf("%-15s %5s %6s %7s %8s %8s %8s %8s %5s\n","","","(d)","","(s^2)","(s^2)","(s)","(s)","");
  printf("------------------------------------------------------------------------------\n");

  for (j=0;j<nflag;j++)
    {
      //      cpgenv(-2000,2000,-1e-6,1e-6,0,1);
      fx[0] = minx;
      fx[1] = maxx;
      fy[0] = fy[1] = j+0.5;
      cpgsci(1);
      cpgline(2,fx,fy);

      n=0;
      for (i=0;i<psr[0].nobs;i++)
	{
	  found=0;
	  for (k=0;k<psr[0].obsn[i].nFlags;k++)
	    {
	      if (strcmp(psr[0].obsn[i].flagVal[k],flagV[j])==0)
		{
		  found=1;
		  break;
		}
	    }
	  if (found==1)
	    {
	      fx[n] = (float)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
	      y[n] = (float)psr[0].obsn[i].residual;
	      id[n] = i;
	      n++;
	    }
	}
      cpgsci((j%2)+1);
      for (i=0;i<n;i++)
	fy[i] = (float)psr[0].obsn[id[i]].residual/yrange;

      if (removeQuad==1)
	{
	  TKremovePoly_f(fx,fy,n,3);
	  TKremovePoly_f(fx,y,n,3);
	}
      for (i=0;i<n;i++)
	{
	  fy[i] += (0.5+j);
	  fe1[i] = (float)(fy[i]-psr[0].obsn[id[i]].toaErr*1.0e-6/yrange);
	  fe2[i] = (float)(fy[i]+psr[0].obsn[id[i]].toaErr*1.0e-6/yrange);
	  ye[i] = psr[0].obsn[id[i]].toaErr*1.0e-6;

	}
      cpgpt(n,fx,fy,1);
      cpgerry(n,fx,fe1,fe2,1);
      cpgsch(0.8); cpgtext(minx+(maxx-minx)*0.05,j+0.8,flagV[j]); cpgsch(1.4);
      //      for (i=0;i<nflag;++i)
      printf("%-15s %5d ",flagV[j],n);
      span = TKrange_f(fx,n);
      printf("%6d ",(int)span);
      npts_sf = determine1dStructureFunction(fx,y,ye,n,&errfac,&vsf,&mverr);
      errfacs[j] = errfac;
      printf("%7d %8.2e %8.2e %8.2e %8.2e %5.2f ",npts_sf,vsf,mverr,sqrt(vsf),sqrt(mverr),errfac);
      //      strcpy(efacFlagID[nefacFlag],flag);
      //      strcpy(efacFlag[nefacFlag],flagV[j]);
      //      efacFlagVal[nefacFlag]=errfac;
      printf("\n");
    }
  printf("------------------------------------------------------------------------------\n");
  writeTim("fixData.tim",psr,"tempo2");

  cpgend();

  for (j=0;j<nflag;j++)
    {
      printf("T2EFAC %s %s %5.2f\n",flag, flagV[j],errfacs[j]);
    }
  printf("------------------------------------------------------------------------------\n");

}

int determine1dStructureFunction(float *x,float *y,float *ye,int nn,double *errfac1,
				 double *mvsf,double *mverr)
{
  int i,p,j,n=0;
  double sf,verr,vsf;
  float tt1,tt2,tt3;
  //  int dayGap=7;
  int changed=0;
  // CHECK SORT 
  do {
    changed=0;
    for (i=0;i<nn-1;i++)
      {
	if (x[i+1] < x[i])
	  {
	    tt1 = x[i+1];
	    tt2 = y[i+1];
	    tt3 = ye[i+1];
	    x[i+1] = x[i];
	    y[i+1] = y[i];
	    ye[i+1] = ye[i];
	    x[i] = tt1;
	    y[i] = tt2;
	    ye[i] = tt3;
	    changed=1;
	    //	  printf("The data set for pulsar is not time sorted - this may, or may not, be a problem\n");
	    //	      exit(1);
	  }
      }
  } while (changed==1);

  sf=0;
  vsf=0;
  verr=0.0;
  n=0;
  for (i=0;i<nn;i++)
    verr += pow(ye[i],2);

  for (i=0;i<nn;i++)
    {
      //           verr += 1.0;
      for (j=i+1;j<nn;j++)
	{
	  if (fabs(x[i] - x[j]) < dayGap)
	    {
	      sf += pow(y[i]-y[j],2);
	      n++;
	    }
	  else
	    {
	      i+=(j-i-1);
	      break;
	    }
	}
    }

  *errfac1 = sqrt(0.5*sf/(double)n/(verr/(double)nn));
  *mvsf = 0.5*sf/(double)n;
  *mverr = verr/(double)nn;
  //      printf("errfac1 = %d %d %g\n",n,nn,errfac1);
  return n;
}

float plotHistogram(float *x,int count,int *flagCol,int nFlag,char flagV[100][16])
{
  int i,j,k;
  float binvalX[100],binvalY[100];
  float binvalY_flag[nFlag][100];
  float area;
  float highest,meanx;
  float mousex,mousey,fontSize;
  float tbin;
  char key;
  int nbin;
  float minx;
  float maxx,ominx;
  int normalise;
  float maxy;
  int colour;
  int log;
  float offset;
  int endit=0;
  int b;
  float fx[2],fy[2];

  nbin = 21;
  normalise = 0;
  colour = -1;
  log = 0;
  //  offset = 0.5;
  offset = 0;

  maxy = -1;
  
  meanx = 0.0;
  for (i=0;i<count;i++)
    meanx+=x[i];
  meanx/=count;
  meanx = 0.0; // Don't remove a mean!

    minx = x[0]-meanx;
      maxx = x[0]-meanx; 
  for (i=0;i<count;i++)
    {
      if (minx > x[i]-meanx) minx = x[i]-meanx;
      if (maxx < x[i]-meanx) maxx = x[i]-meanx;
      } 
  if (fabs(minx) > maxx) maxx = fabs(minx);
  else minx = -maxx;
  //  tbin = maxx*2/(double)nbin;

  printf("Histogram: have minx = %g, maxx = %g\n",minx,maxx);
  ominx = minx;
  minx -= 2.0*(maxx-ominx)/nbin;
  maxx += 2.0*(maxx-ominx)/nbin;

  for (i=0;i<nbin;i++)
    {
      binvalX[i]=i*(maxx-minx)/nbin+minx+(maxx-minx)/nbin/2.0*offset;
      binvalY[i]=0.0;
      for (j=0;j<nFlag;j++)
	binvalY_flag[j][i]=0.0;
    }
  for (i=0;i<count;i++)
    {
      b = (int)((x[i]-meanx-minx)/((maxx-minx)/nbin));
      if (b>=0 && (int)(b<nbin))
	{
	  binvalY[b]++;
	  binvalY_flag[flagCol[i]-2][b]++;
	}
    }

  if (normalise>0)
    {
      area=0.0;
      for (i=0;i<nbin;i++)
	area+=binvalY[i];
      for (i=0;i<nbin;i++)
	{
	  binvalY[i]=binvalY[i]/area*normalise;
	  for (j=0;j<nFlag;j++)
	    binvalY_flag[j][i]=binvalY_flag[j][i]/area*normalise;
	}
    }
  /* Normalise so that highest point lies at 1 */
  highest=0.0;
  for (i=0;i<nbin;i++)
    {
      if (highest<binvalY[i])
	highest=binvalY[i];
    }
  if (maxy<0) maxy = (float)highest+0.1*highest;
  if (normalise<0)
    {
      if (highest!=0.0)
	{
	  for (i=0;i<nbin;i++)
	    {
	      binvalY[i]/=highest;
	      for (j=0;j<nFlag;j++)
		binvalY_flag[j][i]/=highest;
	    }
	}
    }
  if (colour==-1)
    {
      if (log==1)
	cpgenv(minx, maxx, 0, maxy, 0, 10); 
      else
	cpgenv(minx, maxx, 0, maxy, 0, 0); 
      cpglab("Normalised residual","Number","");
      cpgqch(&fontSize); cpgsch(0.5); cpgmtxt("B",8,0.92,0.0,"fixData v.1.0 (G. Hobbs)"); cpgsch(fontSize); 
      
      cpgsls(1);
      cpgsci(1);
    }
  else if (colour==-2)
    {
      cpgswin(minx,maxx,0,maxy);
      cpgsci(1);
      cpgsls(1);
      cpgslw(1);
      
    }
  else
    {
      cpgsci(colour);
      /*      cpgsci(1); */
      cpgsls(1);
      cpgslw(1);
      if (colour==2) cpgslw(5);
      if (colour==3) cpgsls(3);
      if (colour==4) cpgsls(5);
      
    }
  
  cpgbin(nbin,binvalX,binvalY,0);
  cpgsfs(1);

  // Make solid histogram
  for (i=nFlag-1;i>=0;i--)
    {
      for (j=0;j<i;j++)
	{
	  for (k=0;k<nbin;k++)
	    binvalY_flag[i][k]+=binvalY_flag[j][k];
	}
    }
  for (i=nFlag-1;i>=0;i--)
    {
      cpgsci(i+2);
      //      cpgbin(nbin,binvalX,binvalY_flag[i],0);
      for (j=0;j<nbin-1;j++)
	{
	  cpgrect(binvalX[j],binvalX[j+1],0,binvalY_flag[i][j]);
	}
    }
  cpgsls(1);
  cpgsci(1);
  cpgbin(nbin,binvalX,binvalY,0);
  cpgslw(1);

  // Plot legend
  fx[0] = minx+(maxx-minx)*0.05;
  fx[1] = minx+(maxx-minx)*0.1;

  for (i=0;i<nFlag;i++)
    {
      fy[0] = fy[1] = highest-i*highest/20.0;
      if (i+2<14)
	cpgsci(i+2);
      else
	cpgsci(14);
      cpgline(2,fx,fy);
      cpgsci(1);
      cpgtext(fx[1],fy[1],flagV[i]);
    }
  cpgsci(1);
  /*  {
    // Plot Gaussian
    float fx[1000],fy[1000];
    float x;
    int n=0;
    for (x=-10;x<10;x+=0.1)
      {
	fx[n] = x;
	fy[n] = highest*exp(-x*x/2.0);
	n++;
      }
    cpgsci(2); cpgslw(4); cpgline(n,fx,fy); cpgslw(1); cpgsci(1);
    }*/
  return highest;
}
char * plugVersionCheck = TEMPO2_h_VER;
