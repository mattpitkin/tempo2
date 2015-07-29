#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <cpgplot.h>
#include "tempo2.h"
#include "T2toolkit.h"
#include "TKspectrum.h"
#include "GWsim.h"

using namespace std;

#define MAX_POLY 30
void lombScargle(pulsar *psr);
void plotResiduals(pulsar *psr);
void shufflePoints(pulsar *psr,long idum);

void shuffle(double *R, double *err, double *R2, double *shuffledE, int N,long *idum);
double calcStat(double *x,double *y,double *e,int n,int type);
void plotHistogram(float *x,int count);
void average(pulsar *psr);
void corr2pt(pulsar *psr,long idum);

void help() /* Display help */
{
  printf("The following tests are available:\n\n");
  printf("-res   Plot residuals (default)\n");
  printf("-ls    Lomb-Scargle periodogram\n");
}


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i;
  int analysis=1;
  long idum=-432;

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: checkWhite\n");
  printf("Author:              G. Hobbs\n");
  printf("Version:             1.0\n");
  printf(" --- type 'h' for help information\n");
  printf(" --- type \"tempo2 -gr checkWhite -h\" for usage instructions\n");


  /* Obtain the .par and the .tim file from the command line */
  if (argc==4) /* Only provided .tim name */
    {
      strcpy(timFile[0],argv[3]);
      strcpy(parFile[0],argv[3]);
      parFile[0][strlen(parFile[0])-3] = '\0';
      strcat(parFile[0],"par");
    }

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[0],argv[i+1]); 
	  strcpy(timFile[0],argv[i+2]);
	}
      else if (strcmp(argv[i],"-res")==0) // Timing residuals
	analysis=1;
      else if (strcmp(argv[i],"-ls")==0) // Lomb-Scargle periodogram
	analysis=2;
      else if (strcmp(argv[i],"-average")==0) // Average data
	analysis=3;
      else if (strcmp(argv[i],"-shuffle")==0) // Shuffle data
	analysis=5;
      else if (strcmp(argv[i],"-corr")==0) // 2-pt correlation
	analysis=6;
      else if (strcmp(argv[i],"-idum")==0)
	sscanf(argv[++i],"%d",&idum);
      else if (strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-help")==0){
	help();
	exit(0);
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
      else textOutput(psr,*npsr,0,0,0,0,"");  /* Display the output */
    }

  if (analysis==1)
    plotResiduals(psr);
  else if (analysis==2)
    lombScargle(psr);
  else if (analysis==3)
    average(psr);
  else if (analysis==5)
    shufflePoints(psr,idum);
  else if (analysis==6)
    corr2pt(psr,idum);
  return 0;
}

void corr2pt(pulsar *psr,long idum)
{
  double x[MAX_OBSN],y[MAX_OBSN];
  double corr=0.0,y2=0.0;
  double err[MAX_OBSN];
  double sY[MAX_OBSN];
  double sE[MAX_OBSN];
  double simCorr[10000];
  int nit=10000;
  int percent=0;
  int it;
  int n;
  int i;

  for (i=0;i<psr[0].nobs;i++)
    {
      x[i] = (double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
      y[i] = (double)psr[0].obsn[i].residual;
      err[i] = (double)psr[0].obsn[i].toaErr*1.0e-6;
    }
  n = psr[0].nobs;
  corr=0.0;
  y2=0.0;
  for (i=0;i<n-1;i++)
    {
      corr+=y[i]*y[i+1];
      y2+=y[i]*y[i];
    }
  corr/=(double)(n-1);
  y2/=(double)(n-1);
  corr/=y2;
  printf("Actual covariance = %g\n",corr);
  for (it=0;it<nit;it++)
    {
      shuffle(y, err, sY, sE, n,&idum);
      simCorr[it]=0;
      y2=0.0;
      for (i=0;i<n-1;i++)
	{
	  simCorr[it]+=sY[i]*sY[i+1];
	  y2+=sY[i]*sY[i];
	}
      y2/=(double)(n-1);
      simCorr[it]/=(double)(n-1);
      simCorr[it]/=y2;
      if (fabs(simCorr[it]) > fabs(corr)) percent++;
      printf("corr: %g\n",simCorr[it]);
    }
  printf("percent above = %g percent: corr = %g\n",percent/(double)(nit)*100.0,corr);
}

void average(pulsar *psr)
{
  double x[MAX_OBSN],y[MAX_OBSN];
  double x2[MAX_OBSN],y2[MAX_OBSN];
  float px[MAX_OBSN],py[MAX_OBSN];
  int n,n2;
  double rms;
  int i,j;
  int n3=0;
  int divide=1;
  double mean=0.0;
  float maxy,miny;

  for (i=0;i<psr[0].nobs;i++)
    {
      x[i] = (double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
      y[i] = (double)psr[0].obsn[i].residual;
    }
  n = psr[0].nobs;
  for (divide=1;divide < 100;divide++)
    {
      n2=0;
      for (i=0;i<n-(divide-1);i+=divide)
	{
	  mean=0.0;
	  for (j=i;j<i+divide;j++)
	    mean+=y[j];
	  y2[n2]=mean/divide;
	  n2++;
	}
      rms = TKfindRMS_d(y2,n2);
      px[n3] = (float)log10(divide);
      py[n3] = (float)log10(rms);
	n3++;
      printf("%d %d %g\n",divide,n2,rms);
    }
  cpgbeg(0,"?",1,1);
  cpgsch(1.4);
  cpgsfs(2);
  miny = TKfindMin_f(py,n3);
  maxy = TKfindMax_f(py,n3);

  cpgenv(TKfindMin_f(px,n3),TKfindMax_f(px,n3),miny-(maxy-miny)*0.5,maxy+(maxy-miny)*0.1,0,30);
  cpglab("Number of points averaged","rms","");
  cpgpt(n3,px,py,9);
  cpgline(n3,px,py);

  px[1] = px[n3-1];
  py[1] = (py[0]+0.5*px[0])-0.5*px[1];
  cpgsci(2); cpgline(2,px,py); cpgsci(1);
  cpgend();

}

void shufflePoints(pulsar *psr,long idum)
{
  int i,it;
  double x[MAX_OBSN],y[MAX_OBSN],e[MAX_OBSN],sy[MAX_OBSN],se[MAX_OBSN];;
  double statReal;
  float statSim[10000];
  float fx[MAX_OBSN],fy[MAX_OBSN];
  int n;
  int nits=1000;

  n=psr[0].nobs;
  for (i=0;i<n;i++)
    {
      x[i] = (double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
      y[i] = (double)(psr[0].obsn[i].residual);
      e[i] = (double)(psr[0].obsn[i].toaErr*1e-6);
    }
  statReal = calcStat(x,y,e,n,2);
  for (it=0;it<nits;it++)
    {
      shuffle(y,e,sy,se,n,&idum);
      statSim[it] = (float)calcStat(x,sy,se,n,2);
      printf("Compare %g %g\n",statReal,statSim[it]);
    }
  cpgbeg(0,"?",1,1);
  cpgsch(1.4);
  cpgsfs(2);
  plotHistogram(statSim,nits);
  cpgend();
}

double calcStat(double *x,double *y,double *e,int n,int type)
{
  int i;
  double ret;
  if (type==1) // Unweighted rms
    ret = TKfindRMS_d(y,n);
  else if (type==2) // Summed power spectrum
    {
      double px[MAX_OBSN],py[MAX_OBSN];
      int nout,jmax;
      double prob,var;
      TKlomb_d(x,y,n,1,1,px,py,&nout,&var);
      ret = 0.0;
      for (i=0;i<nout;i++)
	ret+=py[i];
    }
  return ret;
}

void plotResiduals(pulsar *psr)
{
  int i,n;
  float fx[MAX_OBSN],fy[MAX_OBSN],fe1[MAX_OBSN],fe2[MAX_OBSN],minx,maxx,miny,maxy;

  n = psr[0].nobs;
  for (i=0;i<n;i++)
    {
      fx[i] = (float)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
      fy[i] = (float)(psr[0].obsn[i].residual);
      fe1[i] = (float)(fy[i]-psr[0].obsn[i].toaErr*1e-6);
      fe2[i] = (float)(fy[i]+psr[0].obsn[i].toaErr*1e-6);
    }
  minx = TKfindMin_f(fx,n);
  maxx = TKfindMax_f(fx,n);
  miny = TKfindMin_f(fy,n);
  maxy = TKfindMax_f(fy,n);
  cpgbeg(0,"?",1,1);
  cpgsch(1.4);
  cpgsfs(2);
  cpgenv(minx,maxx,miny,maxy,0,0);
  cpglab("Day","Residual (s)","");
  cpgpt(n,fx,fy,5);
  cpgerry(n,fx,fe1,fe2,1);
  cpgend();
}

void lombScargle(pulsar *psr)
{
  int i;
  double x[MAX_OBSN],y[MAX_OBSN];
  double px[MAX_OBSN],py[MAX_OBSN],var,prob;
  float fx[MAX_OBSN],fy[MAX_OBSN],minx,maxx,miny,maxy;
  int nout,jmax;
  double effm;

  for (i=0;i<psr[0].nobs;i++)
    {
      if (psr[0].obsn[i].deleted==1)
	{
	  printf("Please remove all deleted points from the .tim file and rerun\n");
	  exit(1);
	}
      x[i] = (double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
      y[i] = (double)(psr[0].obsn[i].residual);
    }
  //  TKperiod_m(x,y,psr[0].nobs,4,1,px,py,MAX_OBSN,&nout,&jmax,&prob,&var,&effm);  
  printf("Effective number of points m = %g\n",effm);
  TKconvertFloat1(px,fx,nout);
  TKconvertFloat1(py,fy,nout);
  minx = TKfindMin_f(fx,nout);
  maxx = TKfindMax_f(fx,nout);
  miny = TKfindMin_f(fy,nout);
  maxy = TKfindMax_f(fy,nout);
  cpgbeg(0,"?",1,1);
  cpgsch(1.4);
  cpgsfs(2);
  cpgenv(minx,maxx,0,maxy+0.1*maxy,0,1);
  cpglab("Frequency (d\\u-1\\d)","Power","");
  printf("Maximum point at (%g,%g) with probability of %g (time scale = %f days)\n",px[jmax],py[jmax],prob,1.0/px[jmax]);
  cpgline(nout,fx,fy);
  fx[0] = minx;
  fx[1] = maxx;
  fy[0] = fy[1] = -log(0.05/effm);
  printf("95 percent confidence level is approximately %g\n",fy[0]);
  cpgslw(4); cpgsls(4); cpgline(2,fx,fy); cpgsls(1); cpgslw(1);
  cpgend();
}

void shuffle(double *R, double *err, double *R2, double *shuffledE, int N,long *idum)
{
  int j,k,i;
  double temp;
  double temp2;

  for (i=0; i<N; i++)
    {
      R2[i] = R[i];
      shuffledE[i] = err[i];
    }

  for (i=0; i<N-1; i++)
    {
      j= (int)(TKranDev(idum)*((N-1)-i)); 
      temp = R2[i+j];
      temp2 = err[i+j];
      
      R2[i+j]=R2[i];
      R2[i]=temp;

      shuffledE[i+j]=shuffledE[i];
      shuffledE[i]=temp2;     
    }
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
  normalise = 0;
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
      cpglab("Statistic","Number","");
      cpgqch(&fontSize); cpgsch(0.5); cpgmtxt("B",8,0.92,0.0,"checkWhite v.1.0 (G. Hobbs)"); cpgsch(fontSize); 

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
}
char * plugVersionCheck = TEMPO2_h_VER;
