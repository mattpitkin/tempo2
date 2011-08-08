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

/* Plugin to determine the GW background limit described by Jenet et al. (2006) */
/* The user provides data-sets that are white. (Note: other plugins can be used */
/* to check whether a data-set is white)                                        */
/* **************************************************************************** */
/* Typical usage:
 *
 * get threshold values
 * > tempo2 -gr GWlimit -f mypar1.par mytim1.tim -f mypar2.par mytim2.tim -npsr 2 
 *    -threshold -fast -idum -543       
 *
 * now get limit
 *
 * > tempo2 -gr GWlimit -f mypar1.par mytim1.tim -f mypar2.par mytim2.tim -npsr 2
 *    -fast -ngw 1000 -idum -43 -dist 1 -dist 1 -alpha -0.6666666 -GSbackground -maxamp 1e-20
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "GWsim.h"
#include "T2toolkit.h"
#include "TKspectrum.h"
#include <time.h>

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

using namespace std;

#define MAX_POLY 8
#define MAX_FLAG 10

#define MAX_FREQ 1000
#define MAX_ITERATION 25000

void getThreshold(pulsar *psr,double *freqVal,int *nFreq,double *threshold,long *idum,int fast,int npsr,long double **Func);
void cumulativeHistogram(double val[MAX_ITERATION][MAX_FREQ],int nval,double *freqVal,int nFreq,double limit,double *threshold);
void cumulativeHistogram2(double *val,int nval,double limit,double *threshold);
void sortit(int n, double array[], double rasort[]);
void shuffle(long double *R, double *toaE, long double *R2, double *toaE2, int N,long *idum);
void checkReal(pulsar psr,double *freqVal,int *nFreq,double *threshold,double alpha);
void getLimits(pulsar *psr,double *freqVal,int *nFreq,double *threshold,long *idum,int checkBackground,double alpha,double *dist,int distNum,double maxAmp,int fast,int npsr,int numberGW,long double alpha,long double **Func);
void setupPulsar(long double ra_p,long double dec_p,long double *kp);
//void GramSchmidt(long double x[], long double y[], long double err[], int ObsAmt, int Npoly, long double CoeffArray[],int wtyn,long double **Func);
void GramSchmidt(long double x[], long double y[], long double err[], int ObsAmt, int Npoly, long double CoeffArray[],int wtyn);
void writeCommands(int argc, char *argv[]);

double storeVal[MAX_ITERATION][MAX_FREQ];

void help() /* Display help */
{
  /* This function should contain usage information about the plugin which should (in general) be accessed */
  /* by the user pressing 'h'                                                                              */
}


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{

  writeCommands(argc,argv);

  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,j,k;
  double globalParameter;
  double threshold[MAX_FREQ],freqVal[MAX_FREQ],threshold2[MAX_FREQ];
  int freqKnown=0;
  int nFreq;
  int fast=0;
  FILE *fin;
  double alpha=1;
  long idum=-847;
  int checkRealFlag=0;
  int getThresholdFlag=0;
  int checkBackground=0;
  double dist[MAX_PSR];  /* Distance to pulsar */
  int distNum=0;
  double maxAmp=-1;
  long double index=-2.0L/3.0L;
  double scale;
  int p;
  int numberGW=1000;
  /* The value of function x at observation y */
  long double **Func;
  Func = (long double **)malloc(MAX_POLY*sizeof(long double *));
  for (i=0;i<MAX_POLY;i++)
    Func[i] = (long double *)malloc(MAX_OBSN_VAL*sizeof(long double));

  /* Set all distances to zero */
  for (i=0;i<MAX_PSR;i++) dist[i]=0.0;

  *npsr = 0;  

  printf("Graphical Interface: GWlimit\n");
  printf("Author:              G. Hobbs, F. Jenet, C. Shettigara\n");
  printf("Version:             v1.0 \n");
  printf(" --- type 'h' for help information\n");

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[*npsr],argv[i+1]); 
	  strcpy(timFile[*npsr],argv[i+2]);
	  (*npsr)++;
	}
      else if (strcmp(argv[i],"-threshold")==0)
	getThresholdFlag=1;
      else if (strcasecmp(argv[i],"-freqknown")==0)
	freqKnown=1;
      else if (strcmp(argv[i],"-real")==0)
	checkRealFlag=1;
      else if (strcmp(argv[i],"-idum")==0)
	sscanf(argv[++i],"%d",&idum);
      else if (strcmp(argv[i],"-ngw")==0)
	sscanf(argv[++i],"%d",&numberGW);
      else if (strcasecmp(argv[i],"-maxamp")==0)
	sscanf(argv[++i],"%lf",&maxAmp);
      else if (strcmp(argv[i],"-fast")==0)
	fast=1;
      else if (strcmp(argv[i],"-dist")==0) /* In kpc */
	{
	  sscanf(argv[++i],"%lf",&dist[distNum]);
	  dist[distNum++] *= 3.08568025e19; /* in m */
	}
      else if (strcmp(argv[i],"-mult")==0)
	sscanf(argv[++i],"%lf",&alpha);
      else if (strcmp(argv[i],"-alpha")==0)
	sscanf(argv[++i],"%Lf",&index);
      else if (strcmp(argv[i],"-background")==0)
	checkBackground=1;
      else if (strcmp(argv[i],"-GSbackground")==0)
	checkBackground=2;
    }

  scale = pow(86400.0*365.25,index);
  maxAmp *= scale;
  printf("GW amplitude = %g\n",maxAmp);

  if (dist[0]==0)
    {
      printf("---------------------------------\n");
      printf("WARNING: pulsar term not included\n");
      printf("---------------------------------\n");
    }

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr);         /* Load the arrival times    */

  preProcess(psr,*npsr,argc,argv);
  /* Find threshold values */
  if (getThresholdFlag==1)
    {
      getThreshold(psr,freqVal,&nFreq,threshold,&idum,fast,*npsr,Func); 
      printf("Completed getting the threshold\n");
      exit(1);    
    }
  /* Read threshold values from a file */
  if (checkBackground==1)
    {
      printf("Reading threshold3.dat\n");
      fin = fopen("threshold3.dat","r");
      fscanf(fin,"%lf",&threshold[0]);
      fclose(fin);
      nFreq=11;
    }
  else if (checkBackground==2)
    {
      printf("Reading threshold4.dat\n");
      fin = fopen("threshold4.dat","r");
      fscanf(fin,"%lf",&threshold[0]);
      fclose(fin);
      nFreq=11;
    }
  else
    {
      fin = fopen("threshold.dat","r");
      while (!feof(fin))
	{
	  if (fscanf(fin,"%lf %lf",&freqVal[nFreq],&threshold[nFreq])==2)
	    nFreq++;
	}
      fclose(fin);
    }

  /*  if (freqKnown==0 && getThresholdFlag==1)
    {
      getThreshold3(psr[0],freqVal,&nFreq,threshold,threshold2,&idum); 
      exit(1);
      } */
  /* Check with real data */

  /* For safety - reload the data */
  if (checkRealFlag==1)
    {
      readParfile(psr,parFile,timFile,*npsr); /* Load the parameters                */
      readTimfile(psr,timFile,*npsr);         /* Load the arrival times             */
      preProcess(psr,*npsr,argc,argv);
      formBatsAll(psr,1);                     /* Form the barycentric arrival times */
      formResiduals(psr,1,0);               /* Form the residuals                 */
      /* Fit the timing model */
      doFit(psr,1,0); 
      formBatsAll(psr,1);                     /* Form the barycentric arrival times */
      formResiduals(psr,1,0);               /* Form the residuals                 */
      
      checkReal(psr[0],freqVal,&nFreq,threshold,alpha);
      exit(1);
    }
  for (p=0;p<*npsr;p++)
    {
      psr[p].nJumps = 0;
      for(i=0;i<MAX_PARAMS;i++){
	psr[p].param[i].nLinkTo = 0;
	psr[p].param[i].nLinkFrom = 0;
      }
    }

  /* Obtain limits on the existence of a single GW source */
  /* For safety - reload the data */
  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters                */
  readTimfile(psr,timFile,*npsr);         /* Load the arrival times             */
  preProcess(psr,*npsr,argc,argv);
  formBatsAll(psr,*npsr);                     /* Form the barycentric arrival times */
  formResiduals(psr,*npsr,0);               /* Form the residuals                 */
  printf("Calling getLimits\n");
  getLimits(psr,freqVal,&nFreq,threshold,&idum,checkBackground,alpha,dist,distNum,maxAmp,fast,*npsr,numberGW,index,Func);
  
  /*  for (i=0;i<psr[0].nobs;i++)
      printf("Residual = %f %g\n",(double)psr[0].obsn[i].bat,(double)psr[0].obsn[i].residual); */

  return 0;
}

/* This function is used to obtain limits on the existence of a single GW source in the data */
/* It works in a similar fashion to the function to obtain the threshold values              */

void getLimits(pulsar *psr,double *freqVal,int *nFreq,double *threshold,long *idum,int checkBackground,double alphaMult,double *dist,int distNum,double maxAmp,int fast,int npsr,int numberGW,long double alpha,long double **Func)
{
  printf("Starting getLimit\n");

  long double **sat1;
  long double **residual,***shuffled;
  long double ***res;
  double ***err;
  long idumStore;
  int nFlags[MAX_PSR],foundFlag;
  int nPts[MAX_PSR][10];
  char flag[MAX_PSR][10][100];
  double xval[MAX_OBSN],yval[MAX_OBSN],py[MAX_OBSN],var[MAX_OBSN];
  long double GSmean,variance,meanX,meanY,rangeX,minXval;
  double globalParameter;
  double gwAmp,minAmp=0.0,amp1,amp2,amp3;
  double limit[MAX_FREQ];
  double **toaE,***shuffledToaE;
  long double kp[MAX_PSR][3];            /* Vector pointing to pulsar           */
  long double gx[MAX_OBSN],gy[MAX_OBSN],ge[MAX_OBSN],coeffArray[10];
  double gsVal,gsValp[MAX_POLY][MAX_PSR],GSspec[MAX_POLY];
  int  nPoly=8;
  FILE *fout,*fout2;
  gwSrc *gw;
  double mean;
  int it,numIt,itAmp,k;
  int detect=0;
  int i,j,p,iFreq,sig, stopIt;
  int jmax;
  double prob,val2,val1,val3;
  long double fhi,flo,max,min;
  double storeJumpVal[MAX_PSR][MAX_JUMPS];
  int wtyn;

  /* Parameters for root finding */
  double eps=3e-8,a,b,c,d,e,min1,min2;
  double fa,fb,fc,pr,q,r,s,tol1,xm,tol=0.1e-22,tol2=0.15;    
  int endit=0;
  
  idumStore = *idum;

  FILE *fout3;

  /* Memory allocation */
  printf("Attempting to allocate memory\n");
  fout3 = fopen("results.dat","a");

  sat1     = (long double **)malloc(MAX_PSR*sizeof(long double*));
  residual = (long double **)malloc(MAX_PSR*sizeof(long double*));
  toaE     = (double **)malloc(MAX_PSR*sizeof(double*));
  shuffled = (long double ***)malloc(MAX_PSR*sizeof(long double **));
  res      = (long double ***)malloc(MAX_PSR*sizeof(long double **));
  err      = (double ***)malloc(MAX_PSR*sizeof(double **));
  shuffledToaE = (double ***)malloc(MAX_PSR*sizeof(double **));

  for (i=0;i<MAX_PSR;i++)
    {
      sat1[i]     = (long double *)malloc(MAX_OBSN*sizeof(long double));
      residual[i] = (long double *)malloc(MAX_OBSN*sizeof(long double));
      toaE[i]     = (double *)malloc(MAX_OBSN*sizeof(double));
      shuffled[i] = (long double **)malloc(MAX_FLAG*sizeof(long double *));
      res[i]      = (long double **)malloc(MAX_FLAG*sizeof(long double *));
      err[i]      = (double **)malloc(MAX_FLAG*sizeof(double *));
      shuffledToaE[i] = (double **)malloc(MAX_FLAG*sizeof(double *));
      for (j=0;j<MAX_FLAG;j++)
	{
	  shuffled[i][j] = (long double *)malloc(MAX_OBSN*sizeof(long double));
	  res[i][j]      = (long double *)malloc(MAX_OBSN*sizeof(long double));
	  err[i][j]      = (double *)malloc(MAX_OBSN*sizeof(double));
	  shuffledToaE[i][j] = (double *)malloc(MAX_OBSN*sizeof(double));	  
	}
    }
  printf("Completed allocating memory\n");



  printf("alpha= %Lg\n",alpha);

  for (p=0;p<npsr;p++)
    {
      for (i=0;i<psr[p].nJumps;i++)
	storeJumpVal[p][i] = psr[p].jumpVal[i];
    }

  if ((gw = (gwSrc *)malloc(sizeof(gwSrc)*numberGW))==NULL)
    {
      printf("Unable to allocate memory for %d GW sources\n",numberGW);
      exit(1);
    }

  if (checkBackground==0)
    fout2 = fopen("limit.dat","w");
  else
    fout2 = fopen("limit_background.dat","w");

  a  = 0;
  fa = 0.1-95.0; /* Percentage of detection given by the threshold value */  
  b  = maxAmp;
  c  = maxAmp;
  
  for (p=0;p<npsr;p++)
    {
      nPts[p][0] = 0;
      nFlags[p] = 0;
      for (i=0;i<psr[p].nobs;i++)
	{
	  toaE[p][i] = psr[p].obsn[i].toaErr;
	  residual[p][i] = 0.0;
	  foundFlag=0;
	  for (j=0;j<nFlags[p];j++)
	    {
	      if (strcmp(flag[p][j],psr[p].obsn[i].flagVal[0])==0)
		{
		  nPts[p][j]++;
		  foundFlag=1;
		  break;
		}
	    }
	  if (foundFlag==0)
	    {
	      strcpy(flag[p][nFlags[p]],psr[p].obsn[i].flagVal[0]);
	      nPts[p][nFlags[p]] ++;
	      nFlags[p]++;
	      nPts[p][nFlags[p]] = 0;
	    }
	}
      printf("Flags are (pulsar number, flag, number of points):\n");
      for (j=0;j<nFlags[p];j++)
	printf("%d ... %s (%d)\n",j,flag[p][j],nPts[p][j]);
    }
  /* Obtain residuals from the observed TOAs and the given timing model */
  printf("Forming bats\n");
  formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
  printf("Forming residuals\n");
  formResiduals(psr,npsr,0);   /* Form the residuals                 */
  
  /* Now obtain TOA2 - a set of perfect TOAs - by removing the residuals from the original TOAs */
  for (j=0;j<2;j++) /* Should actually keep iterating until the rms residual is below some level ... */
    {
      for (p=0;p<npsr;p++)
	{
	  for (i=0;i<psr[p].nobs;i++)
	    {
	      if (psr[p].obsn[i].deleted == 1)
		{
		  printf("PROBLEM: REMOVE COMMENT LINES IN .TIM FILE AND RERUN\n");
		  exit(1);
		}
	      residual[p][i] += psr[p].obsn[i].residual;
	      psr[p].obsn[i].sat -= psr[p].obsn[i].residual/SECDAY;
	      sat1[p][i] = psr[p].obsn[i].sat;
	    }
	  formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
	  formResiduals(psr,npsr,0);   /* Form the residuals                 */
	}
    }
  for (p=0;p<npsr;p++)
    {
      for (j=0;j<nFlags[p];j++)
	nPts[p][j]=0;
      for (i=0;i<psr[p].nobs;i++)
	{
	  for (j=0;j<nFlags[p];j++)
	    {
	      if (strcmp(psr[p].obsn[i].flagVal[0],flag[p][j])==0)
		{
		  res[p][j][nPts[p][j]]=residual[p][i];
		  err[p][j][nPts[p][j]]=toaE[p][i];
		  nPts[p][j]++;
		}
	    }
	}
    }


  numIt = 1000;
  sig = numIt-(int)(0.90*numIt);
  printf("sig = %i \n", sig);

  /*  flo = 3.0e-10;
      fhi = 1.0e-6; */

  /* Setup flo and fhi based on the first pulsar data-set */
  /* Number of points -> power of 2 -> Largest amount of time  */

  max = psr[0].obsn[0].sat;
  min = psr[0].obsn[0].sat;
  for (j=0;j<npsr;j++)
    {
      for (i=0;i<psr[j].nobs;i++)
	{
	  if (max < psr[j].obsn[i].sat) max = psr[j].obsn[i].sat;
	  if (min > psr[j].obsn[i].sat) min = psr[j].obsn[i].sat;
	}
    }

  flo = 0.01L/((max-min)*86400.0L);  

  /* Rick: half a day */
  /*  fhi = 1.0L/(86400.0L);  */
  fhi = 2.0L/(86400.0L); 

  /*  alpha = -2.0L/3.0L;*/

  for (p=0;p<npsr;p++)
    setupPulsar(psr[p].param[param_raj].val[0],psr[p].param[param_decj].val[0],kp[p]);

  //  for (iFreq = 10; iFreq < *nFreq;iFreq++) 
  for (iFreq = 0; iFreq < 1;iFreq++) 
    {
      endit=0;
      if (checkBackground==0)
	{
	  minAmp = 0.0;
	  maxAmp = 1e-6;
	}
      
      itAmp=0;

      amp1 = maxAmp;
      gwAmp = amp1;
      amp2 = -1;

      do
	{
	  detect = 0;

	  for (it=0;it<numIt;it++)
	    {	      
	      if (it%100==0) printf("Iteration GW sim %d %g detect = %d, percent = %g\n",it+1,gwAmp,detect,detect/(double)(it)*100.0); 
	      if (checkBackground==1 || checkBackground==2)  /* Create a new stochastic background with amplitude gwAmp */
		GWbackground(gw,numberGW,idum,flo,fhi,gwAmp,alpha,1);
	      

	      for (p=0;p<npsr;p++)
		{
		  for (j=0;j<nFlags[p];j++)
		    shuffle(res[p][j],err[p][j],shuffled[p][j],shuffledToaE[p][j],nPts[p][j],idum);
		  
		  for (j=0;j<nFlags[p];j++)
		    nPts[p][j]=0;
		  
		  for (i=0;i<psr[p].nobs;i++)
		    {
		      for (j=0;j<nFlags[p];j++)
			{
			  if (strcmp(psr[p].obsn[i].flagVal[0],flag[p][j])==0)
			    {
			      psr[p].obsn[i].sat = sat1[p][i]+shuffled[p][j][nPts[p][j]]/SECDAY;
			      psr[p].obsn[i].toaErr = shuffledToaE[p][j][nPts[p][j]];
			      nPts[p][j]++;
			      break;
			    }
			}
		      if (checkBackground==0)
			psr[p].obsn[i].sat += gwAmp*sin(freqVal[iFreq]*2.0*M_PI*(psr[p].obsn[i].sat-psr[p].param[param_pepoch].val[0]))/SECDAY;
		      else  /* Simulate a stochastic background */
			{
			  for (k=0;k<numberGW;k++)
			    {
			      psr[p].obsn[i].sat += calculateResidualGW(kp[p],&gw[k],(psr[p].obsn[i].sat-psr[0].param[param_pepoch].val[0])*86400.0L,dist[p])/86400.0L; 
			    } 
			}
		    }	
		}
	      /*	      writeTim("test.tim",&psr,"tempo2"); */

	      if (fast==1)
		{
		  vectorPulsar(psr,npsr);  
		  calculate_bclt(psr,npsr);
		  formBats(psr,npsr);      
		  formResiduals(psr,npsr,0);
		}
	      else
		{
		  formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
		  formResiduals(psr,npsr,0);   /* Form the residuals                 */
		}
	      /* SHOULD FIT */
	      doFit(psr,npsr,0);  
	      if (fast==1)
		{
		  vectorPulsar(psr,npsr);   /* 1. Form a vector pointing at the pulsar */
		  calculate_bclt(psr,npsr);           /* 3. Calculate bclt (WHAT IS THIS?) */
		  formBats(psr,npsr);                   /* Form Barycentric arrival times */
		  formResiduals(psr,npsr,0);   /* Form the residuals                 */  
		}
	      else
		{
		  formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
		  formResiduals(psr,npsr,0);   /* Form the residuals                 */
		}
	      /* Restore pulsar parameters */
	      for (i=0;i<MAX_PARAMS;i++)
		{
		  for (p=0;p<npsr;p++)
		    {
		      for (j=0;j<psr[p].param[i].aSize;j++)
			{
			  if (psr[p].param[i].fitFlag[j]==1)
			    psr[p].param[i].val[j] = psr[p].param[i].prefit[j];
			}
		    }
		}
	      for (p=0;p<npsr;p++)
		{
		  for (i=0;i<psr[p].nJumps;i++)
		    psr[p].jumpVal[i] = storeJumpVal[p][i];
		}

	      /* Now obtain the Lomb periodogram of these data */
	      for (p=0;p<npsr;p++)
		{
		  mean=0.0;
		  meanX = 0.0L;
		  meanY = 0.0L;

		  for (i=0;i<psr[p].nobs;i++)
		    {
		      xval[i] = (double)(psr[p].obsn[i].bat - psr[p].param[param_pepoch].val[0]);
		      yval[i] = (double)(psr[p].obsn[i].residual); 
		      gx[i] = (psr[p].obsn[i].bat - psr[p].param[param_pepoch].val[0]);
		      gy[i] = psr[p].obsn[i].residual;
		      ge[i] = psr[p].obsn[i].toaErr*1.0e-6;
		      meanX += gx[i];
		      meanY += gy[i];
		      mean  += yval[i];
		    }
		  minXval = gx[0];
		  rangeX = gx[psr[p].nobs-1]-gx[0];
		  /* Remove mean from residuals */
		  mean/=(double)psr[p].nobs;
		  meanY /=(long double)psr[p].nobs;
		  meanX /=(long double)psr[p].nobs;
		  
		  
		  for (i=0;i<psr[p].nobs;i++) 
		    {
		      yval[i]-=mean; 
		      gy[i]-=meanY;  
		      gx[i]=(gx[i]-minXval)/(rangeX)*2.0-1.0L; 
		    }


		  if (it<100 && p==0) 
		    {
		      char fname[100];
		      sprintf(fname,"residual_gw%d.dat",it);
		      
		      fout = fopen(fname,"w");
		      for (i=0;i<psr[p].nobs;i++)
			fprintf(fout,"%f %g %s\n",(double)gx[i],(double)gy[i],psr[p].obsn[i].flagVal);  
		      fclose(fout);
		    }  
		  
		  if (checkBackground==2)
		    {
		      wtyn = 1;
		      //		      GramSchmidt(gx, gy, ge, psr[p].nobs, nPoly, coeffArray,wtyn,Func);		      
		      GramSchmidt(gx, gy, ge, psr[p].nobs, nPoly, coeffArray,wtyn);
		      /* Get variance of data set */
		      variance=0.0;
		      GSmean = 0.0;
		      for (i=0;i<psr[p].nobs;i++)
			GSmean += gy[i];
		      GSmean/=(long double)psr[p].nobs;
		      //		      for (i=0;i<psr[p].nobs;i++)
		      //			variance+=((gy[i]-GSmean)*(gy[i]-GSmean));
		      //		      variance/=(psr[p].nobs-1.0);
		      for (i=0;i<psr[p].nobs;i++)
			variance+=(pow(gy[i]-GSmean,2)/powl(ge[i],2));
		      variance/=(long double)(psr[p].nobs);

		      /* Scale coefficients */
		      for (i=0;i<nPoly;i++)	
			{
			  coeffArray[i]/=sqrt(variance); 
			  gsValp[i][p] = (double)(coeffArray[i]*coeffArray[i]);
			}
		    }
		  else
		    {
		      double var;
		      //		      TKperiod(xval-1,yval-1,psr[p].nobs,1,1,freqVal-1,py-1,MAX_OBSN,nFreq,&jmax,&prob,&var);  
		      TKlomb_d(xval,yval,psr[p].nobs,1,1,freqVal,py,nFreq,&var);  
		    }
		  /*      TKlomb_d(xval,yval,psr.nobs,1,1,freqVal,py,nFreq,var);   */
		  /*      TKlomb_slow_d(xval,yval,psr.nobs,freqVal,py,nFreq,var);*/
		  if (it<100 && p==0) 
		    {
		      char fname[100];
		      sprintf(fname,"lomb_gw%d.dat",it);
		      fout = fopen(fname,"w");
		      for (i=0;i<*nFreq;i++)
			fprintf(fout,"%g %g %d\n",freqVal[i],py[i],i);     
		      fclose(fout);
		    }  
		}
	      /*	      printf("Checking %g %g %g\n",py[iFreq],threshold[iFreq],threshold[22]); */
	      if (checkBackground==0)
		{
		  if (py[iFreq] > threshold[iFreq]*alphaMult)
		    detect++;
		}
	      else if (checkBackground==1)
		{
		  double sum1,sum2;
		  sum1 = 0.0;
		  sum2 = threshold[0];
		  for (k=0;k<8;k++)
		    sum1+=py[k];
		  if (sum1 > sum2)
		    detect++;
		}
	      else
		{
		  gsVal = 0.0;
		  for (i=3;i<nPoly;i++)
		    {
		      GSspec[i] = 0.0;
		      for (p=0;p<npsr;p++)
			GSspec[i] += gsValp[i][p];
		      
		      gsVal += GSspec[i];
		    }

		  if (gsVal > threshold[0])
		    detect++;
		}
	    }
	  printf("(%g) Percentage of detection = %g (%d) [%g]\n",gwAmp,(double)detect/(double)numIt*100.0,detect,threshold[0]);
	  fprintf(fout3,"%d %g %g\n",idumStore,gwAmp,(double)detect/(double)numIt*100.0);
	  endit=1;
	  

	  if (checkBackground==0 || checkBackground==1 || checkBackground==2)
	    {
	      /*		  if ((double)detect/(double)numIt*100.0 > 95)
				  maxAmp = gwAmp;
				  else
				  minAmp = gwAmp; */
	      /* Use Brent's method for root finding */
	      fb = (double)detect/(double)numIt*100.0-95.0;
	      if (itAmp==0) /* First time through */
		{
		  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		    {
		      printf("Root not bracketed\n");
		      gwAmp*=10.0;
		      itAmp=-1;
		    }		      
		  fc = fb;
		}
	      
	      if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
		{
		  c=a;
		  fc=fa;
		  e=d=b-a;
		  printf("Setting d at 1 %g\n",d);
		}
	      if (fabs(fc) < fabs(fb))
		{
		  a=b;
		  b=c;
		  printf("Setting b = %g\n",b);
		  c=a;
		  fa=fb;
		  fb=fc;
		  fc=fa;
		}
	      printf("c b (2) = %g %g\n",c,b);
	      tol1 = 2.0*eps*fabs(b)+0.5*tol;
	      xm=0.5*(c-b);
	      printf("Have %g %g %g\n",fabs(xm),tol1,fb);
	      /* if (fabs(xm/1.0e-20) <= tol1 || fabs(fb) < tol2) */
	      if (fabs(fb) < tol2)
		{
		  printf("GOT A SOLUTION %g %g %g\n",fabs(xm/1.0e-20),tol1,fb);
		  endit=1;
		}
	      else
		{
		  if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
		    {
		      s = fb/fa;
		      if (a==c)
			{
			  pr = 2.0*xm*s;
			  q = 1.0-s;
			}
		      else
			{
			  q = fa/fc;
			  r = fb/fc;
			  pr = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
			  q=(q-1.0)*(r-1.0)*(s-1.0);
			}
		      printf("b (3) = %g\n",b);
		      if (pr > 0.0) q = -q;
		      pr=fabs(pr);
		      min1 = 3.0*xm*q-fabs(tol1*q);
		      min2 = fabs(e*q);
		      if (2.0*pr < (min1 < min2 ? min1 : min2))
			{
			  e = d;
			  d = pr/q;
			  printf("Setting d at 2 %g\n",d);
			}
		      else
			{
			  d=xm;
			  e=d;
			  printf("Setting d at 3 %g\n",d);
			}			  
		      printf("b (4) = %g\n",b);
		    }
		  else
		    {
		      d=xm;
		      e=d;
		      printf("Setting d at 4 %g\n",d);
		    }
		  printf("b (5) = %g %g\n",b,d);
		  a=b;
		  printf("b (5.1) = %g %g\n",b,d);
		  fa=fb;
		  printf("b (5.2) = %g %g\n",b,d);
		  if (fabs(d) > tol1)
		    {
		      printf("b (5.3) = %g %g\n",b,d);
		      printf("Setting b in 1 %g %g\n",b,d);
		      b += d;
		    }
		  else			
		    {
		      printf("Setting b in 2 %g %g\n",b,SIGN(tol1,xm));
		      b += SIGN(tol1*1e-20,xm);
		    }
		  printf("b (6) = %g\n",b);
		
		  if (itAmp!=-1)
		    gwAmp = b;
		  printf("b (7) = %g\n",b);
		  itAmp++;     
		}
	    }	
	}
      while (endit==0); 
      
      if (itAmp == 100)
	limit[iFreq] = -1;
      else
	limit[iFreq] = gwAmp;
      printf("Complete after %d iterations with amplitude of %g\n",itAmp,gwAmp);
      if (checkBackground==0)
	fprintf(fout2,"%g %g\n",freqVal[iFreq],gwAmp);
      else
	fprintf(fout2,"%g\n",gwAmp);
      fflush(fout2);
      if (checkBackground==1 || checkBackground==2)
	break;
    }

if (checkBackground==0)
    {
      for (i=0;i<*nFreq;i++)
	printf("%g %g\n",freqVal[i],limit[i]);
    }
  fclose(fout2);
  fclose(fout3);

}

void checkReal(pulsar psr,double *freqVal,int *nFreq,double *threshold,double alpha)
{
  double xval[MAX_OBSN],yval[MAX_OBSN],py[MAX_FREQ],var[MAX_OBSN];
  int i;
  FILE *fout;
  int jmax;
  double prob;

  for (i=0;i<psr.nobs;i++)
    {
      xval[i] = (double)(psr.obsn[i].bat - psr.param[param_pepoch].val[0]);
      yval[i] = (double)(psr.obsn[i].residual);
    }
  {
    double var;
    //  TKperiod(xval-1,yval-1,psr.nobs,1,1,freqVal-1,py-1,MAX_OBSN,nFreq,&jmax,&prob,&var);  
  TKlomb_d(xval,yval,psr.nobs,1,1,freqVal,py,nFreq,&var);  
  }
  /*  TKlomb_d(xval,yval,psr.nobs,1,1,freqVal,py,nFreq,var);   */
  fout = fopen("checkReal.dat","w");
  for (i=0;i<*nFreq;i++)
    {
      if (py[i] > threshold[i]*alpha)
	printf("Significant detection in real data at %g %g (threshold = %g) [period = %g (d) = %g (yr)]\n",freqVal[i],py[i],threshold[i],1.0/freqVal[i],1.0/freqVal[i]/365.25);
      fprintf(fout,"%g %g %g\n",freqVal[i],py[i],threshold[i]);
    }
  fclose(fout);
}

/* Obtain a threshold power for each frequency bin in a Lomb periodogram */
void getThreshold(pulsar *psr,double *freqVal,int *nFreq,double *threshold,long *idum,int fast,int npsr,long double **Func)
{
  int i,j,k,it=0;
  long double **sat1; 
  long double **residual;
  /*  long double shuffled[MAX_PSR][MAX_FLAG][MAX_OBSN]; */
  long double ***shuffled;
  long double ***res,meanX,meanY,rangeX,minXval;
  double ***err;
  double mean;
  char flag[MAX_PSR][MAX_FLAG][20];
  int nFlags[MAX_PSR],foundFlag;
  double **toaE,***shuffledToaE;
  int    nPts[MAX_PSR][MAX_FLAG];
  long double gx[MAX_OBSN],gy[MAX_OBSN],ge[MAX_OBSN],coeffArray[100];
  double gsVal[MAX_ITERATION],***gsValp,GSthreshold,GSspec[MAX_POLY][MAX_ITERATION];
  int    nPoly = 8;
  double xval[MAX_OBSN],yval[MAX_OBSN],py[MAX_OBSN],var[MAX_OBSN];
  long double variance,GSmean;
  double globalParameter;
  double alpha,alphaMax,alphaMin;
  int printRes = 10;
  int printLomb = 10;
  int nobs=0;
  FILE *fout;
  parameter param[MAX_PARAMS];
  int jmax,p;
  double prob;
  int detect;
  double sum1,sum2,sval[MAX_OBSN],BGthreshold;
  int wtyn=1;
  int ngreater;
  long double datastat = 0.0;

  /* Memory allocation */
  printf("Attempting to allocate memory\n");
  
  sat1     = (long double **)malloc(MAX_PSR*sizeof(long double*));
  residual = (long double **)malloc(MAX_PSR*sizeof(long double*));
  toaE     = (double **)malloc(MAX_PSR*sizeof(double*));
  shuffled = (long double ***)malloc(MAX_PSR*sizeof(long double **));
  res      = (long double ***)malloc(MAX_PSR*sizeof(long double **));
  err      = (double ***)malloc(MAX_PSR*sizeof(double **));
  shuffledToaE = (double ***)malloc(MAX_PSR*sizeof(double **));
  gsValp = (double ***)malloc(MAX_POLY*sizeof(double **));


  for (i=0;i<MAX_POLY;i++)
    {
      gsValp[i] = (double **)malloc(MAX_PSR*sizeof(double *));
      for (j=0;j<MAX_PSR;j++)
	gsValp[i][j] = (double *)malloc(MAX_ITERATION*sizeof(double));
    }

  for (i=0;i<MAX_PSR;i++)
    {
      sat1[i]     = (long double *)malloc(MAX_OBSN*sizeof(long double));
      residual[i] = (long double *)malloc(MAX_OBSN*sizeof(long double));
      toaE[i]     = (double *)malloc(MAX_OBSN*sizeof(double));
      shuffled[i] = (long double **)malloc(MAX_FLAG*sizeof(long double *));
      res[i]      = (long double **)malloc(MAX_FLAG*sizeof(long double *));
      err[i]      = (double **)malloc(MAX_FLAG*sizeof(double *));
      shuffledToaE[i] = (double **)malloc(MAX_FLAG*sizeof(double *));
      for (j=0;j<MAX_FLAG;j++)
	{
	  shuffled[i][j] = (long double *)malloc(MAX_OBSN*sizeof(long double));
	  res[i][j]      = (long double *)malloc(MAX_OBSN*sizeof(long double));
	  err[i][j]      = (double *)malloc(MAX_OBSN*sizeof(double));
	  shuffledToaE[i][j] = (double *)malloc(MAX_OBSN*sizeof(double));	  
	}
    }
  printf("Completed allocating memory\n");

  for (p=0;p<npsr;p++)
    {
      nFlags[p]=0;
      nPts[p][0] = 0;
      for (i=0;i<psr[p].nobs;i++)  /* FIRST FLAG MUST BE IDENTIFIER */
	{
	  residual[p][i] = 0.0;
	  toaE[p][i] = psr[p].obsn[i].toaErr;
	  foundFlag=0;
	  for (j=0;j<nFlags[p];j++)
	    {
	      if (strcmp(flag[p][j],psr[p].obsn[i].flagVal[0])==0)
		{
		  nPts[p][j]++;
		  foundFlag=1;
		  break;
		}
	    }
	  if (foundFlag==0)
	    {
	      strcpy(flag[p][nFlags[p]],psr[p].obsn[i].flagVal[0]);
	      nPts[p][nFlags[p]] ++;
	      (nFlags[p])++;
	      nPts[p][nFlags[p]] = 0;
	    }
	}
    }

  printf("Flags are:\n");
  for (p=0;p<npsr;p++)
    {
      printf("psr = %d %d\n",p,npsr);
      for (j=0;j<nFlags[p];j++)
	printf("%d ... %s (%d)\n",p,flag[p][j],nPts[p][j]);
    }
  /* Obtain residuals from the observed TOAs and the given timing model */
  formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
  formResiduals(psr,npsr,0);   /* Form the residuals                 */


  /* Get the GramSchmidt statistic from the original data. -ANL 10/1/07 */
  for (p=0;p<npsr;p++)
     {
       mean=0.0;
       meanY=0.0L;
       meanX=0.0L;
       for (i=0;i<psr[p].nobs;i++)
	 {
	   gx[i] = (psr[p].obsn[i].bat - psr[p].param[param_pepoch].val[0]);
	   gy[i] = psr[p].obsn[i].residual;
	   ge[i] = psr[p].obsn[i].toaErr*1e-6; 
	   
	   meanX += gx[i];
	   meanY += (gy[i]);
	   mean+=(double)(psr[p].obsn[i].residual);  
	 }
       minXval = gx[0];
       rangeX = gx[psr[p].nobs-1]-gx[0];
       /* Remove mean from residuals */
       mean/=(double)psr[p].nobs; 
       meanY /=(long double)psr[p].nobs;
       meanX /=(long double)psr[p].nobs;
       
       for (i=0;i<psr[p].nobs;i++) 
	 {
	   gy[i]-=meanY;  
	   gx[i]=(gx[i]-minXval)/(rangeX)*2.0-1.0L; 
	 }
       GramSchmidt(gx, gy, ge, psr[p].nobs, nPoly, coeffArray,wtyn);
       /* Get variance of data set */
       variance=0.0L;
       GSmean = 0.0L;
       for (i=0;i<psr[p].nobs;i++)
	 GSmean += gy[i];
       GSmean/=(long double)psr[p].nobs;
       
       for (i=0;i<psr[p].nobs;i++)
	 variance+=(pow(gy[i]-GSmean,2)/powl(ge[i],2));
       variance/=(long double)(psr[p].nobs);
       /* Scale coefficients */
       for (i=0;i<nPoly;i++)	
	 {
	   coeffArray[i]/=sqrt(variance); 
	   gsValp[i][p][it] = (double)(coeffArray[i]*coeffArray[i]);
	 }
     }      
  datastat = 0.0;
  for (i=3;i<nPoly;i++)
    {
      GSspec[i][it] = 0.0;
      for (p=0;p<npsr;p++)
	GSspec[i][it] += gsValp[i][p][it];
      
      datastat += GSspec[i][it];
    }
  for (i=0;i<nPoly;i++)
    {
      printf("GSspec %d %g\n", i, (double)GSspec[i][it]);
    }
  printf("Statistic in the data is %Lg\n", datastat);
  
  /* Now obtain TOA2 - a set of perfect TOAs - by removing the residuals from the original TOAs */
  for (j=0;j<2;j++) /* Should acutally keep iterating until the rms residual is below some level ... */
    {
      for (p=0;p<npsr;p++)
	{	 
	  for (i=0;i<psr[p].nobs;i++)
	    {
	      if (psr[p].obsn[i].deleted == 1)
		{
		  printf("PROBLEM: REMOVE COMMENT LINES IN .TIM FILE AND RERUN\n");
		  exit(1);
		}
	      residual[p][i] += psr[p].obsn[i].residual;
	      psr[p].obsn[i].sat -= psr[p].obsn[i].residual/SECDAY;
	      sat1[p][i] = psr[p].obsn[i].sat;
	    }
	}
      formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,npsr,0);   /* Form the residuals                 */
    }

  for (p=0;p<npsr;p++)  /* Reset flags */
    {
      for (j=0;j<nFlags[p];j++) nPts[p][j]=0;
      for (i=0;i<psr[p].nobs;i++)
	{
	  for (j=0;j<nFlags[p];j++)
	    {
	      if (strcmp(psr[p].obsn[i].flagVal[0],flag[p][j])==0)
		{
		  res[p][j][nPts[p][j]]=residual[p][i];
		  err[p][j][nPts[p][j]]=toaE[p][i];
		  if (err[p][j][nPts[p][j]]==0)
		    {
		      printf("Error = 0, %f %d\n",toaE[p][i],i);
		      exit(1);
		    }
		  nPts[p][j]++;
		}
	    }
	}
    }
  
  /* Now add to the perfect TOAs either a shuffled set of residuals or Gaussian noise of given amplitude */
  for (it=0;it<10000;it++)
    {
      if ((it+1)%100 == 0) {printf("Iteration %d                \r",it+1); fflush(stdout);}
     
      for (p=0;p<npsr;p++)
	{
	  for (j=0;j<nFlags[p];j++)
	    shuffle(res[p][j],err[p][j],shuffled[p][j],shuffledToaE[p][j],nPts[p][j],idum);
	  for (j=0;j<nFlags[p];j++)
	    nPts[p][j]=0;
	  for (i=0;i<psr[p].nobs;i++)
	    {
	      for (j=0;j<nFlags[p];j++)
		{
		  if (strcmp(psr[p].obsn[i].flagVal[0],flag[p][j])==0)
		    {
		      psr[p].obsn[i].sat    = sat1[p][i]+shuffled[p][j][nPts[p][j]]/SECDAY;
		      psr[p].obsn[i].toaErr = shuffledToaE[p][j][nPts[p][j]];
		      if (psr[p].obsn[i].toaErr == 0.0)
			{
			  printf("Error on TOA (%f) %d %d %d = 0.0\n",(double)psr[p].obsn[i].sat,i,j,nPts[p][j]);
			  exit(1);
			}
		      nPts[p][j]++;
		      break;
		    }
		}
	    }
	}
      if (fast==1) /* Skip re-doing clock corrections and reading the ephemeris */
	{
	  vectorPulsar(psr,npsr);   
	  calculate_bclt(psr,npsr); 
	  formBats(psr,npsr);       
	  formResiduals(psr,npsr,0); 
	}
      else  /* Do full procedure */
	{
	  formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
	  formResiduals(psr,npsr,0);   /* Form the residuals                 */
	}
      /* SHOULD FIT */
      doFit(psr,npsr,0); 
      /* Form postfit residuals */
      if (fast==1)
	{
	  vectorPulsar(psr,npsr);   /* 1. Form a vector pointing at the pulsar */
	  calculate_bclt(psr,npsr);           /* 3. Calculate bclt (WHAT IS THIS?) */
	  formBats(psr,npsr);                   /* Form Barycentric arrival times */
	  formResiduals(psr,npsr,0);   /* Form the residuals                 */	  
	}
      else
	{
	  formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
	  formResiduals(psr,npsr,0);   /* Form the residuals                 */
	}
      /* Restore pulsar parameters */
      for (i=0;i<MAX_PARAMS;i++)
	{
	  for (p=0;p<npsr;p++)
	    {
	      for (j=0;j<psr[p].param[i].aSize;j++)
		{
		  if (psr[p].param[i].fitFlag[j]==1)
		    psr[p].param[i].val[j] = psr[p].param[i].prefit[j];
		}
	    }
	} 
      /* Now obtain the Lomb periodogram and Gram-Schmidt polynomials of these data */
      for (p=0;p<npsr;p++)
	{
	  mean=0.0;
	  meanY=0.0L;
	  meanX=0.0L;
	  for (i=0;i<psr[p].nobs;i++)
	    {
	      gx[i] = (psr[p].obsn[i].bat - psr[p].param[param_pepoch].val[0]);
	      gy[i] = psr[p].obsn[i].residual;
	      ge[i] = psr[p].obsn[i].toaErr*1e-6; 
	      	
	      meanX += gx[i];
	      meanY += (gy[i]);
	      xval[i] = (double)(psr[p].obsn[i].bat - psr[p].param[param_pepoch].val[0]);
	      yval[i] = (double)(psr[p].obsn[i].residual); 
	      mean+=(double)(psr[p].obsn[i].residual);  
	    }
	  minXval = gx[0];
	  rangeX = gx[psr[p].nobs-1]-gx[0];
	  /* Remove mean from residuals */
	  mean/=(double)psr[p].nobs; 
	  meanY /=(long double)psr[p].nobs;
	  meanX /=(long double)psr[p].nobs;

	  for (i=0;i<psr[p].nobs;i++) 
	    {
	      yval[i]-=mean; 
	      gy[i]-=meanY;  
	      gx[i]=(gx[i]-minXval)/(rangeX)*2.0-1.0L; 
	    }
	  
	  if (p==0)
	    {
	      /*	      TKperiod(xval-1,yval-1,psr[p].nobs,1,1,freqVal-1,py-1,MAX_OBSN,nFreq,&jmax,&prob); */
	    }
	  /*      TKlomb_d(xval,yval,psr.nobs,1,1,freqVal,py,nFreq,var); */

	  /*	  if (it<100 && p==0) 
	    {
	      char fname[100];
	      sprintf(fname,"lomb%d.dat",it);
	      fout = fopen(fname,"w");
	      for (i=0;i<*nFreq;i++)
		fprintf(fout,"%g %g %d\n",freqVal[i],py[i],i);     
	      fclose(fout);
	      }*/
	  for (i=0;i<*nFreq;i++)
	    storeVal[it][i] = py[i];  
	  if (it<100 && p==0) 
	    {
	      char fname[100];
	      sprintf(fname,"residual%d.dat",it);
	      
	      fout = fopen(fname,"w");
	      for (i=0;i<psr[p].nobs;i++)
		fprintf(fout,"%Lg %Lg %g\n",gx[i],gy[i],psr[p].obsn[i].toaErr/1e6);  
	      fclose(fout);
	    }
	  /* Now look at orthogonal polynomials */
	  //	  GramSchmidt(gx, gy, ge, psr[p].nobs, nPoly, coeffArray,wtyn,Func);
	  GramSchmidt(gx, gy, ge, psr[p].nobs, nPoly, coeffArray,wtyn);
	  /* Get variance of data set */
	  variance=0.0L;
	  GSmean = 0.0L;
	  for (i=0;i<psr[p].nobs;i++)
	    GSmean += gy[i];
	  GSmean/=(long double)psr[p].nobs;
	  //	  for (i=0;i<psr[p].nobs;i++)
	  //	    variance+=((gy[i]-GSmean)*(gy[i]-GSmean));
	  //	    variance/=(psr[p].nobs-1.0);

	  for (i=0;i<psr[p].nobs;i++)
	    variance+=(pow(gy[i]-GSmean,2)/powl(ge[i],2));
	  variance/=(long double)(psr[p].nobs);
	  /* Scale coefficients */
	  /*	  printf("Variance = %g %g\n",(double)variance,(double)(sqrt(variance)*1.0e9));*/
	  for (i=0;i<nPoly;i++)	
	    {
	      coeffArray[i]/=sqrt(variance); 
	      gsValp[i][p][it] = (double)(coeffArray[i]*coeffArray[i]);
	    }
	  //	  for (i=0; i<nPoly;i++)
	  //  gsValp[p][it] += (double)(coeffArray[i]*coeffArray[i]);
	  /*	  for (i=0;i<nPoly;i++)
	    printf("Have %d %g %g %g\n",i,(double)coeffArray[i],(double)(coeffArray[i]*coeffArray[i]),
	    (double)(coeffArray[i]*sqrt(variance)*sqrt(psr[p].nobs)));  */
	}      
      gsVal[it] = 0.0;
      for (i=3;i<nPoly;i++)
	{
	  GSspec[i][it] = 0.0;
	  for (p=0;p<npsr;p++)
	    GSspec[i][it] += gsValp[i][p][it];

	  gsVal[it] += GSspec[i][it];
	}

       if (it < 100) 
       {
	      char fname[100];
	      sprintf(fname,"gs%d.dat",it);
	      fout = fopen(fname,"w");
              for (i=0;i<nPoly;i++) fprintf(fout, "%d %g\n", i, (double)GSspec[i][it]);
	      fclose(fout);
        }
    }

  /* Now find threshold */
  cumulativeHistogram(storeVal,it,freqVal,*nFreq,0.999,threshold);

  /* Now find threshold value for any point being above the level */
  printf("Searching for second threshold value\n");

  alphaMax = 10;
  alphaMin = 1;
  
  printf("Final alpha = %g\n",alpha);
  fout = fopen("threshold2.dat","w");
  for (i=0;i<*nFreq;i++)
    fprintf(fout,"%g %g\n",freqVal[i],threshold[i]*alpha); 
  fclose(fout);
  /*  for (i=0;i<psr.nobs;i++)
      printf("Residual0 = %f %g\n",xval[i],yval[i]);  */
  /*      for (i=0;i<*nFreq;i++)
	printf("Freq %g %g %d\n",freqVal[i],py[i],i);     
      exit(1);
      printf("Storing %d %g\n",it,freqVal[123]); */

  /* Now get threshold for background */
  /* Now find threshold value for any point being above the level */
  /*  printf("Searching for the background threshold value \n");
  for (i=0;i<it;i++)
    {
      sval[i] = 0.0;
      for (j=0;j<8;j++)
	sval[i]+=storeVal[i][j];
    }
  cumulativeHistogram2(sval,it,0.999,&BGthreshold);  
  fout = fopen("threshold3.dat","w");
  fprintf(fout,"%g\n",BGthreshold);
  printf("BGthreshold = %g\n",BGthreshold);
  fclose(fout); */


  /* Now get threshold for background using Gram-Schmidt*/
  double gsmean=0.0;
  for (i=0;i<it;i++)
    gsmean+=(double)gsVal[i];
  printf("mean GSthreshold = %g\n",(double)gsmean/(double)it);

  cumulativeHistogram2(gsVal,it,0.999,&GSthreshold);  
  fout = fopen("threshold4.dat","w");
  fprintf(fout,"%g\n",GSthreshold);
  printf("GSthreshold = %g\n",GSthreshold);
  fclose(fout);

  /* Andrea's addition to print out what fraction of the shuffled */
  /* statistics are larger than the statistic in the data */
  ngreater=0;
  for (i=0; i< it; i++)
     {
	if (gsVal[i] > datastat) ++ngreater;
     }
  printf("%d out of %d\n", ngreater, it);
  printf("of the shuffled statistics are larger than the original statistic\n");
  printf("Original statistic (from data) is %Lg\n", datastat);

  /* Now must deallocate the memory */
  

}

void cumulativeHistogram(double val[MAX_ITERATION][MAX_FREQ],int nval,double *freqVal,int nFreq,double limit,double *threshold)
{
  int i,j;
  int bin[MAX_FREQ];
  double array[MAX_ITERATION],sortArray[MAX_ITERATION];
  double binwidth;
  double max,min;
  int nbin,f;
  FILE *fout,*fout2;
  char str[1000];

  fout = fopen("threshold.dat","w");
  for (f=0;f<nFreq;f++)
    {
      for (i=0;i<nval;i++)
	array[i] = val[i][f];
      /*pass array into function sortit. receive sorted array*/
      sortit(nval,array,sortArray);
      threshold[f] = sortArray[(int)(limit*nval)];

      printf("Threshold = %d %g (%d)\n",f,threshold[f],nFreq);
      fprintf(fout,"%g %g\n",freqVal[f],threshold[f]); 
    }
  fclose(fout);
}

void cumulativeHistogram2(double *val,int nval,double limit,double *threshold)
{
  int i,j;
  int bin[MAX_FREQ];
  double array[MAX_ITERATION],sortArray[MAX_ITERATION];
  double binwidth;
  double max,min;
  int nbin,f;
  
  for (i=0;i<nval;i++)
    array[i] = val[i];
  /*pass array into function sortit. receive sorted array*/
  sortit(nval,array,sortArray);
  *threshold = sortArray[(int)(limit*nval)];
}



void shuffle(long double *R, double *toaE, long double *R2, double *toaE2, int N,long *idum)
{
  int j,k,i;
  long double temp;
  double temp2;

  for (i=0; i<N; i++)
    {
      R2[i] = R[i];
      toaE2[i] = toaE[i];
    }

  for (i=0; i<N-1; i++)
    {
      j= (int)(TKranDev(idum)*((N-1)-i));
      temp = R2[i+j];
      temp2 = toaE2[i+j];

      R2[i+j]=R2[i];
      toaE2[i+j] = toaE2[i];
      R2[i]=temp;
      toaE2[i] = temp2;
    }
}

void sortit(int n, double array[], double rasort[])
{
  int l, j, ir, i;
  double rra;
  
  for(i=0; i<n; i++)  /*copy array into rasort for sorting*/
    {
      rasort[i]= array[i];
    }

  l=(n >> 1)+1;
  printf("\n n %i l %i \n \n", n, l);
  ir=n;

  for(;;)
    {
      if (l > 1)
	rra=rasort[(--l)-1];
      else
	{
	  rra=rasort[ir-1];
	  rasort[ir-1]=rasort[0];
	  if(--ir == 1)
	    {
	      rasort[0]=rra;
	      return;
	    }
	}
      i=l;
      j=l<<1;
      while(j<= ir)
	{
	  if (j< ir && rasort[j-1] < rasort[j])
	    ++j;
	  if (rra < rasort[j-1])
	    {
	      rasort[i-1]=rasort[j-1];
	      j += (i=j);
	    }
	  else j = ir+1;
	}
      rasort[i-1]=rra;
    }

}

/* Set up pulsar: note: in KJ this is based in Gwave.cpp: Load_Pulsar_Data */
/* Sets up the vector pointing at the pulsar                               */
void setupPulsar(long double ra_p,long double dec_p,long double *kp)
{
  long double deg2rad = M_PI/180.0;

  kp[0] = cos(dec_p)*cos(ra_p);
  kp[1] = cos(dec_p)*sin(ra_p);
  kp[2] = sin(dec_p);
}


//void GramSchmidt(long double *x, long double *y, long double *err, int ObsAmt, int Npoly, long double *CoeffArray, int wtyn,long double **Func){
void GramSchmidt(long double *x, long double *y, long double *err, int ObsAmt, int Npoly, long double *CoeffArray, int wtyn) {
  /* x is an array with the TOA's (any unit will do.)
     y is an array with residuals (units irrelevant, but same units as err)
     err is an array with 1sigma standard deviations on y
     ObsAmt is how many measurement points there are (stored in x/y/err[0] to [Npts-1])
     Npoly is the total number of polynomials (the highest order polynomial will be of order (Npoly-1))
     Coeffarray is an array that will contain the C values of the Npoly polynomials.*/

  /* Counting variables for the polynomial number, the observation number and 
     the order of the polynomial term under consideration */
  int poly, obs, odr,i;

  long double a[Npoly][Npoly]; 

  /* The value of function x at observation y */

long double **Func;
Func = (long double **)malloc(MAX_POLY*sizeof(long double *));
for (i=0;i<MAX_POLY;i++)
    Func[i] = (long double *)malloc(MAX_OBSN_VAL*sizeof(long double));
  long double N1[Npoly], D1[Npoly], D2[Npoly];
  long double sw;

  /*  for (i=0;i<ObsAmt;i++)
      err[i] = 1.0; */

  /* Initialising values (putting to zero or to the known values) */
  /* ============================================================ */
  for (poly=0;poly<Npoly;poly++){
      for (odr=0;odr<Npoly;odr++){
	  a[poly][odr] = 0.0L;
	}
      for (obs=0;obs<ObsAmt;obs++){	   
	Func[poly][obs] = 0.0L;				
      }							
      N1[poly] = 0.0L;					
      D1[poly] = 0.0L;					
      D2[poly] = 0.0L;					
  }							

  a[0][0] = 1.0L;					
  a[1][1] = 1.0L;
  sw = 0.0;

  if (wtyn ==1){
    for (obs=0;obs<ObsAmt;obs++){ 			    
      a[1][0] -= (x[obs]/pow(err[obs],2));
      sw      += 1.0/pow(err[obs],2);
    }				
    a[1][0]=a[1][0]/sw;
  }
  else if (wtyn == 0) {
    for (obs=0;obs<ObsAmt;obs++)
      a[1][0] -= x[obs]/ObsAmt;
    sw = ObsAmt;
  }

  // We can already determine the polynomial values for the 0th and 1st order polynomials.
  for (poly = 0;poly<=1;poly++){	
    for (obs=0;obs<ObsAmt;obs++){
      for (odr = 0; odr<=poly; odr++){	    
	Func[poly][obs] += a[poly][odr]*powl(x[obs],odr);
      }
    }
  }	  				

  /* Calculating the polynomial coefficients */				
  /* ======================================= */ 
  for (poly=2;poly<Npoly;poly++){
    
    /* N1, D1 and D2 are intermediate values in the calculation */	
    for (obs=0;obs<ObsAmt;obs++){		
      if(wtyn==1){
	N1[poly-1] += powl(Func[poly-1][obs],2) * x[obs]/powl(err[obs],2);
	D1[poly-1] += powl(Func[poly-1][obs],2)/powl(err[obs],2);		
      }		
      else if(wtyn ==0){
	N1[poly-1] += powl(Func[poly-1][obs],2) * x[obs];
	D1[poly-1] += powl(Func[poly-1][obs],2);
      }
    }							
    
    if (poly == 2){
      D2[poly] = sw; 
    }
    else if (poly > 2) D2[poly] = D1[poly-2];
    
    /* Calculating the actual coefficients */					      
    a[poly][0] = -a[poly-1][0]*N1[poly-1]/D1[poly-1] - a[poly-2][0]*D1[poly-1]/D2[poly];  
    for (odr = 1;odr<Npoly;odr++){					         
      a[poly][odr] = a[poly-1][odr-1]-a[poly-1][odr]*N1[poly-1]/D1[poly-1]-a[poly-2][odr]*D1[poly-1]/D2[poly];  
    }										    
    
    for (obs=0;obs<ObsAmt;obs++){					    
      for (odr=0;odr<=poly;odr++){					    
	Func[poly][obs] += a[poly][odr]*pow(x[obs],odr);		    
      }									    
    }								
  } // Finished calculating the polynomial coefficients (and hence the polynomials)
  
  /* Normalising the polynomials 					
     =========================== */					
  D1[0] = sw;								
  D1[Npoly-1] = 0.0L;							

  if (wtyn == 1){
    for (obs = 0; obs < ObsAmt; obs++)
      D1[Npoly-1] += powl(Func[Npoly-1][obs],2)/powl(err[obs],2);
  }
  else if (wtyn == 0){
    for (obs = 0; obs < ObsAmt; obs++)
      D1[Npoly-1] += powl(Func[Npoly-1][obs],2);
  }
  
  for (poly=0;poly<Npoly;poly++){					
    for (odr=0;odr<Npoly;odr++){					
      a[poly][odr] = a[poly][odr]/sqrt(D1[poly]);      
    }									
    for (obs=0;obs<ObsAmt;obs++){					
      Func[poly][obs] = Func[poly][obs]/sqrt(D1[poly]);
    }  									
  } 

  if(wtyn==1){
    for (poly=0;poly<Npoly;poly++){
      CoeffArray[poly]=0.0L;
      for (obs = 0; obs < ObsAmt; obs++){
	CoeffArray[poly] += Func[poly][obs]*y[obs]/powl(err[obs],2);	
	/*      printf("err = %g \n",(double)err[obs]); */
      }
      /*    printf("CoeffArray = %g\n",(double)CoeffArray[poly]); */
    }
  }
  else if (wtyn == 0){
    for (poly = 0; poly<Npoly;poly++){
      CoeffArray[poly]=0.0L;
      for (obs = 0; obs <ObsAmt; obs++){
	CoeffArray[poly] += Func[poly][obs]*y[obs];
      }
    }
  }

  // Now Func is multiplied to result in the function values of the fitted polynomials.
  //  for (poly = 0; poly<Npoly;poly++){
  //  for (obs = 0; obs <ObsAmt; obs++){
  //    Func[poly][obs] = CoeffArray[poly]*Func[poly][obs];
  //  }
  //  }




  /* Check orthonormality */
  /*  {
    long double sum=0.0L;
    int i,j,k;
    
    for (i=0;i<Npoly;i++)
      {
	for (j=0;j<Npoly;j++)
	  {
	    sum = 0.0L;
	    for (k=0;k<ObsAmt;k++)
	      sum+=Func[i][k]*Func[j][k]/powl(err[k],2);
	    printf("%d %d %g\n",i,j,(double)sum);
	  }	
      }  
      }  
      exit(1);  */

  for (i=0;i<MAX_POLY;i++)
    free(Func[i]);
  free(Func);
  //  exit(1);
}

void writeCommands(int argc, char *argv[]){
  int ii;

  char commandfile[200] = "T2command.input";
  FILE *fout;
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  char timeval[200];
  strcpy(timeval,asctime (timeinfo));
  strcpy(&timeval[(int)strlen(timeval)-1],"");

  fout = fopen(commandfile,"a");
  fprintf(fout,"[%s]>> ",timeval);
  for(ii = 0;ii<argc;ii++){
    fprintf(fout," %s ",argv[ii]);
  }
  fprintf(fout,"\n");
  fclose(fout);
}
char * plugVersionCheck = TEMPO2_h_VER;
