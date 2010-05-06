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
#include "tempo2.h"
#include "GWsim.h"

using namespace std;

void doPlugin1(pulsar *psr,char *flag);
void determine1dStructureFunction(float *x,float *y,float *ye,int nn,double *errfac1);
void doPlugin2(pulsar *psr,char parFile[MAX_PSR_VAL][MAX_FILELEN],char timFile[MAX_PSR_VAL][MAX_FILELEN],int argc,char *argv[]);
void plotHistogram(float *x,int count);
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
  printf("\t-daygap <int>     number of days to determine white level\n");
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
  int i;
  int plot=1;
  double globalParameter;
  //  int nit = 1;
  //  long double gwamp = 0;
  //  long double alpha = 1;

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: fixData\n");
  printf("Author:              G. Hobbs, D. Champion\n");
  printf("Version:             1.0\n");

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
      else if (strcmp(argv[i],"-plotout")==0)
      {
	plotoutSet = 1;
	strcpy(plotout,argv[++i]);
      }
      else if (strcmp(argv[i],"-daygap")==0)
	sscanf(argv[++i],"%d",&dayGap);
      else if (strcmp(argv[i],"-script")==0)
      {
	script = 1;
      }
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

  if (plot==1)
    doPlugin1(psr,flag);
  else if (plot==2)
    doPlugin2(psr,parFile,timFile,argc,argv);
  else
    printf("Unknown plot required\n");

  return 0;
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
      plotHistogram(rms,it+1);
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
  plotHistogram(rms,nit);
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

void doPlugin1(pulsar *psr,char *flag)
{
  int nflag=0;
  char flagV[MAX_JUMPS][100];
  double errfacs[MAX_JUMPS];
  float fx[MAX_OBSN],fy[MAX_OBSN],fe1[MAX_OBSN],fe2[MAX_OBSN],ye[MAX_OBSN],y[MAX_OBSN];
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
  printf("-----------------------------------------------------------------------\n");
  printf("%-15s %5s %6s %5s\n","Data set","Npts","Span","EFAC");
  printf("%-15s %5s %6s %5s\n","","","(d)","");
  printf("-----------------------------------------------------------------------\n");

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
	      fy[n] = (float)psr[0].obsn[i].residual/yrange+0.5+j;
	      fe1[n] = (float)(fy[n]-psr[0].obsn[i].toaErr*1.0e-6/yrange);
	      fe2[n] = (float)(fy[n]+psr[0].obsn[i].toaErr*1.0e-6/yrange);
	      ye[n] = psr[0].obsn[i].toaErr*1.0e-6;
	      n++;
	    }
	}
      cpgsci((j%2)+1);
      cpgpt(n,fx,fy,1);
      cpgerry(n,fx,fe1,fe2,1);
      cpgsch(0.8); cpgtext(minx+(maxx-minx)*0.05,j+0.8,flagV[j]); cpgsch(1.4);
      //      for (i=0;i<nflag;++i)
      printf("%-15s %5d ",flagV[j],n);
      span = TKrange_f(fx,n);
      printf("%6d ",(int)span);
      determine1dStructureFunction(fx,y,ye,n,&errfac);
      errfacs[j] = errfac;
      printf("%5.2f ",errfac);
      //      strcpy(efacFlagID[nefacFlag],flag);
      //      strcpy(efacFlag[nefacFlag],flagV[j]);
      //      efacFlagVal[nefacFlag]=errfac;
      printf("\n");
    }
  printf("-----------------------------------------------------------------------\n");
  writeTim("fixData.tim",psr,"tempo2");

  cpgend();

  for (j=0;j<nflag;j++)
    {
      printf("T2EFAC %s %s %5.2f\n",flag, flagV[j],errfacs[j]);
    }
  printf("-----------------------------------------------------------------------\n");

}

void determine1dStructureFunction(float *x,float *y,float *ye,int nn,double *errfac1)
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
      //      printf("errfac1 = %d %d %g\n",n,nn,errfac1);
}

void plotHistogram(float *x,int count)
{
  int i;
  float binvalX[100],binvalY[100];
  float area;
  float highest,meanx;
  float mousex,mousey,fontSize;
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


  nbin = 20;
  normalise = 1;
  colour = -1;
  log = 0;
  offset = 0.5;

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
  ominx = minx;
  minx -= 2.0*(maxx-ominx)/nbin;
  maxx += 2.0*(maxx-ominx)/nbin;

  for (i=0;i<100;i++)
    {
      binvalX[i]=i*(maxx-minx)/nbin+minx+(maxx-minx)/nbin/2.0*offset;
      binvalY[i]=0.0;
    }
  for (i=0;i<count;i++)
    {
      if ((int)((x[i]-meanx-minx)/((maxx-minx)/nbin))>=0 
        && (int)((x[i]-meanx-minx)/((maxx-minx)/nbin)) <100)
      binvalY[(int)((x[i]-meanx-minx)/((maxx-minx)/nbin))]++;
    }

  if (normalise>0)
    {
      area=0.0;
      for (i=0;i<nbin;i++)
      area+=binvalY[i];
      for (i=0;i<nbin;i++)
      binvalY[i]=binvalY[i]/area*normalise;
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
	    binvalY[i]/=highest;
	}
    }
  if (colour==-1)
    {
      if (log==1)
      cpgenv(minx, maxx, 0, maxy, 0, 10); 
      else
      cpgenv(minx, maxx, 0, maxy, 0, 0); 
      cpglab("RMS (\\gms)","Number","");
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
  cpgsls(1);
  cpgsci(1);
  cpgslw(1);

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
}
