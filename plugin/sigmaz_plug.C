/* Direct copy of the Matsakis, Taylor & Eubanks (1997) sigmaz software                          */
/* The code has been directly copied from the fortran - and looks rather ugly in its 'c' form
 * - at some stage I intend to write this using more standard 'c'
 *
 * Use "tempo2 -gr sigmaz -h" to get some help information
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <cpgplot.h>
#include "tempo2.h"
#include "T2toolkit.h"
#include "TKfit.h"
#include "GWsim.h"

#define MAX_GWS 10000

using namespace std;


/* Global variables used in the fortran code */
int npt,nusewt,nxunits,ntunits,nformat,nwriteres,nbintype;
int npt1last,npt2last,ncubic,ncubics,ntau,linfile,indx[90000],ndim;
double data[90000],utjd[90000],taumin,sigmai[90000],permax,root2;
double utjd1,utjd2,tmin,tmax,xmin,xmax,utjdlast,tausec,taumax,tauday;
double prtl[5],utmean,secyear,taulog,addvar,tauyear,tauensure,tdiffmin;
double utfirst,utlast;

void readin(pulsar psr);
void getprtj(int n);
void indexx8(int n,double *arrin,int *indx);
void getweights(int n, double *wt);
void fit4(int *nfit,double *p4,double *cov4,int ndostats,double *chidf,double *avewt);
void mat20(double sam[21][21],double a[21][21],int n,double *determ,int *nbad);
void simWhiteFunc(pulsar *psr,long *idum,char parFile[MAX_PSR_VAL][MAX_FILELEN],char timFile[MAX_PSR_VAL][MAX_FILELEN],int weights,double mintau);

void calcSigmaz(pulsar psr,int weights,double *tau,double *szbias, double *e1,
		double *e2,int *nval, double mintau);
void doplot(pulsar *psr,int npsr,char *grDev,float mint,float maxt,float minsz,float maxsz,
	    int style,int average,double tau[MAX_PSR_VAL][100],double szbias[MAX_PSR_VAL][100],
	    double e1[MAX_PSR_VAL][100],double e2[MAX_PSR_VAL][100],int nval[MAX_PSR_VAL],
	    int nWhite,float *white,char *cline,int slopes,int bound);


//void calcSigmaz(pulsar *psr,int npsr,int time,int psrfile,double mintau,
//		int nSel,int *selPsr,int weights,int style,int average,
//		float *avX,float *avY,int *nav);
void sortTimes(pulsar psr,int *nobs,double *times,double *resid,double *error);
void fitv(double x,double afunc[],int ma,pulsar *psr,int ipos);
void plotOmega_g(double omega,float *px,float *py);
void plotA_g(double a,double alpha,float *px,float *py);
void shufflePts(long double *R, double *toaE, long double *R2, double *toaE2, int N,long *idum);
void convert_gravWaveBackground_noFit(pulsar *psr,int npsr,double convertGW,long *idum,int sameBackground);
void convert_gravWaveBackground_fit(pulsar *psr,int npsr,double convertGW,long *idum,int sameBackground,  
				    char parFile[MAX_PSR_VAL][MAX_FILELEN],char timFile[MAX_PSR_VAL][MAX_FILELEN]);
void calcSpline(float *px,float *py,int count);
float SplineBlend(int k,int t,int *u,float v);
void calculateGWlim(pulsar *psr,long *idum,double obs,char parFile[MAX_PSR_VAL][MAX_FILELEN],char timFile[MAX_PSR_VAL][MAX_FILELEN],int weights,double mintau,long double gwamp,float mint,float maxt,float minsz,float maxsz,double *szbias,double *e1obs,double *e2obs,int nit,int ngw,long double lowAmp,long double hiAmp);
 
typedef struct XY
{
  double x;
  double y;
}XY;

void help() /* Display help */
{
  printf("sigmaz:\n\n");
  printf("Common usage: tempo2 -gr sigmaz -f mypar.par mytim.tim -mintau 1 -maxt 50\n");
  printf("\n\n");
  printf("Command line arguments:\n\n");
  printf("-allpartim       select all the .par and .tim files in the current directory instead of using -f option\n");
  printf("-average         calculate average sigmaz values\n");
  printf("-bound           draw line at sigmaz values that provides bounds of current sensitivity \n");
  printf("-convertWhite x  convert the pulsar residuals to white noise with an rms of x\n");
  printf("-g x             set pgplot graphics device\n");
  printf("-image           create plot suitable for presentation\n");
  printf("-list            give output listing of tau, szbias values\n");
  printf("-minsz x         set minimum log(sz) value in plot\n");
  printf("-maxsz x         set maximum log(sz) value in plot\n");
  printf("-mint  x         set minimum tau value in years for plot\n");
  printf("-maxt  x         set maximum tau value in years for plot\n");
  printf("-noweights       do not use TOA uncertainties when determining sigmaz\n");
  printf("-p n             only plot the sigmaz plots for pulsar 'n' (can be used multiple times, e.g. -p 3 -p 5)\n");
  printf("-ran x           set random number seed\n");
  printf("-slopes          draw 'wheel' showing different power law slopes\n");


  exit(1);
}


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  char grDev[100]="?";
  int i,n,j;
  int res=0;
  char resFile[100][MAX_FILELEN],dummy[100];
  double globalParameter;
  double dpx[MAX_OBSN],dpy[MAX_OBSN],de[MAX_OBSN];
  double cVal[6],eVal[6];
  char name[10][100];
  int nname=0;
  FILE *fin,*fout;
  int it;
  int sim=0;
  int publish=0;
  int nlabel=0;
  int autoScale=0;
  char label[100][100];
  float lx[100],ly[100];
  double mintau=-1;
  int style=0;
  long idum = 0;         
  int nit=100;
  int shuffle=0;
  int nres=0;
  int ngw=10000;
  float maxt = 20;  // Maximum tau in years
  float mint = 0.001; // Minimum tau in years
  float minsz = -17;
  float maxsz = -9;
  char cline[1000]="";
  int allpartim=0;
  int nSel=0;
  int selPsr[MAX_PSR];
  int nWhite=0;
  float white[100];
  float convertWhite=-1.0;
  int weights=-1;
  int average=0;
  int p;
  int slopes=0;
  int listing=0;
  int bound=0;
  double convertGW=0.0;
  long double gwamp=0.0;
  long double lowAmp=0.0;
  long double hiAmp=0.0;
  int simWhite=0;

  double tau[MAX_PSR_VAL][100];
  double szbias[MAX_PSR_VAL][100];
  double e1[MAX_PSR_VAL][100];
  double e2[MAX_PSR_VAL][100];
  int nval[MAX_PSR_VAL];

  int calculateGWlimit=0;

  FILE *pin;
  
  *npsr = 0;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: sigmaz\n");
  printf("Author:              George Hobbs\n");
  printf("Version:             1.0\n");
  printf(" --- type 'h' for help information\n");


  /* Obtain the .par and the .tim file from the command line */
  if (argc==4) /* Only provided .tim name */
    {
      strcpy(timFile[0],argv[3]);
      strcpy(parFile[0],argv[3]);
      parFile[0][strlen(parFile[0])-3] = '\0';
      strcat(parFile[0],"par");
    }

  /* Obtain all parameters from the command line */
  for (i=0;i<argc;i++)
    {
      strcat(cline,argv[i]);
      strcat(cline," ");
      if (strcmp(argv[i],"-h")==0)
	{
	  help();
	  exit(0);
	}
    }
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[*npsr],argv[i+1]); 
	  strcpy(timFile[*npsr],argv[i+2]);
	  (*npsr)++;
	}
      else if (strcmp(argv[i],"-image")==0) // For talks
	style=1;
      else if (strcmp(argv[i],"-publish")==0) // For publicatinos
	style=2;
      else if (strcmp(argv[i],"-image2")==0) // For talks
	style=3;
      else if (strcmp(argv[i],"-image3")==0) 
	style=4;
      else if (strcmp(argv[i],"-res")==0)
	{
	  strcpy(resFile[nres],argv[i+1]);
	  (nres)++;
	}
      else if (strcmp(argv[i],"-name")==0)
	{
	  strcpy(name[nname++],argv[i+1]);
	}
      else if (strcmp(argv[i],"-average")==0)
	sscanf(argv[++i],"%d",&average);
      else if (strcmp(argv[i],"-simWhite")==0)
	simWhite=1;
      else if (strcmp(argv[i],"-g")==0)
	strcpy(grDev,argv[++i]);
      else if (strcmp(argv[i],"-mintau")==0) // Minimum value of tau for stopping in days
	sscanf(argv[++i],"%lf",&mintau);
      else if (strcmp(argv[i],"-maxt")==0) // Maximum tau in years
	sscanf(argv[++i],"%f",&maxt);
      else if (strcmp(argv[i],"-mint")==0) // Minimum tau in years for plotting
	sscanf(argv[++i],"%f",&mint);
      else if (strcmp(argv[i],"-maxsz")==0) // Maximum log(sz)
	sscanf(argv[++i],"%f",&maxsz);
      else if (strcmp(argv[i],"-minsz")==0) // Minimum log(sz)
	sscanf(argv[++i],"%f",&minsz);
      else if (strcmp(argv[i],"-allpartim")==0) // Select all par and tim files in directory
	allpartim=1;
      else if (strcmp(argv[i],"-p")==0) // Select pulsars to plot
	sscanf(argv[++i],"%d",&selPsr[nSel++]);
      else if (strcmp(argv[i],"-w")==0) // Choose white noise levels to plot (s)
	sscanf(argv[++i],"%g",&white[nWhite++]);
      else if (strcmp(argv[i],"-ran")==0 || strcmp(argv[i],"-idum")==0) 
	sscanf(argv[++i],"%d",&idum);
      else if (strcasecmp(argv[i],"-convertWhite")==0)
	sscanf(argv[++i],"%g",&convertWhite);
      else if (strcasecmp(argv[i],"-convertGW")==0)
	sscanf(argv[++i],"%lf",&convertGW);
      else if (strcasecmp(argv[i],"-noweights")==0)
	weights=0;
      else if (strcasecmp(argv[i],"-gwamp")==0)
	sscanf(argv[++i],"%Lf",&gwamp);
      else if (strcasecmp(argv[i],"-nit")==0)
	sscanf(argv[++i],"%d",&nit);
      else if (strcasecmp(argv[i],"-ngw")==0)
	sscanf(argv[++i],"%d",&ngw);
      else if (strcmp(argv[i],"-list")==0)
	listing=1;
      else if (strcmp(argv[i],"-slopes")==0)
	slopes=1;
      else if (strcmp(argv[i],"-gwlimit")==0)
	calculateGWlimit=1;
      else if (strcasecmp(argv[i],"-gwlow")==0)
	sscanf(argv[++i],"%Lf",&lowAmp);
      else if (strcasecmp(argv[i],"-gwhigh")==0)
	sscanf(argv[++i],"%Lf",&hiAmp);
      else if (strcasecmp(argv[i],"-gwclose")==0)
	{
	  sscanf(argv[++i],"%Lf",&gwamp);
	  hiAmp = gwamp+0.1*gwamp;
	  lowAmp = gwamp*0.9;	  
	  printf("searching between %g and %g\n",(double)lowAmp,(double)hiAmp);
	}
      else if (strcmp(argv[i],"-bound")==0)
	bound=1;
    }

  if (idum==0)
    idum = TKsetSeed();

  if (allpartim==1)
    {
      n=0;
      pin = popen("ls *.par","r");
      while (!feof(pin))
	{
	  if (fscanf(pin,"%s",parFile[n])==1)
	    n++;
	}
      pclose(pin);
      n=0;
      pin = popen("ls *.tim","r");
      while (!feof(pin))
	{
	  if (fscanf(pin,"%s",timFile[n])==1)
	    n++;
	}
      pclose(pin);
      *npsr = n;
    }
  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */  
  preProcess(psr,*npsr,argc,argv);
  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,0);    /* Form the residuals                 */
      if (i==0) doFit(psr,*npsr,0);   /* Do the fitting     */
      //	  else textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
    }
  printf("Read pulsars\n");
  if (nres>0)
    {
      for (p=0;p<nres;p++)
	{
	  fin = fopen(resFile[p],"r");
	  psr[*npsr].nobs=0;
	  while(!feof(fin))
	    {
	      if (fscanf(fin,"%Lf %Lf %lf",&(psr[*npsr].obsn[psr[*npsr].nobs].sat),
			 &(psr[*npsr].obsn[psr[*npsr].nobs].residual),
			 &(psr[*npsr].obsn[psr[*npsr].nobs].toaErr))==3)
		{
		  psr[*npsr].nobs++;
		}
	    }
	  fclose(fin);
	  (*npsr)++;
	  //	  calcSigmaz(psr[0],weights,tau[p],szbias[p],e1[p],e2[p],&nval[p],mintau);	  
	}      
    }
  if (nname>0)
    {
      for (p=0;p<*npsr;p++)
	strcpy(psr[p].name,name[p]);
    }
  if (convertWhite>0 && convertGW > 0.0)
    {
      printf("ERROR: Cannot use -convertWhite and -convertGW\n");
      exit(1);
    }
  if (convertWhite > 0)
    {
      for (p=0;p<*npsr;p++)
	{
	  for (i=0;i<psr[p].nobs;i++)
	    {
	      psr[p].obsn[i].residual = TKgaussDev(&idum)*convertWhite;
	      psr[p].obsn[i].toaErr = convertWhite/1.0e-6;
	    }
	}
    }
  if (convertGW > 0.0)
    {
      convert_gravWaveBackground_fit(psr,*npsr,convertGW,&idum,1,parFile,timFile);
      //      convert_gravWaveBackground_noFit(psr,*npsr,convertGW,&idum,1);
    }
  printf("npsr=%d\n",*npsr);
  if (*npsr==0)
    {
      printf("You have not read in any pulsars or data files\n");
      exit(1);
    }
  
  
  for (p=0;p<*npsr;p++)
    {
      //      if (p==0) weights=1;
      //      if (p==1) weights=0;
      for (i=0;i<psr[p].nobs;i++)
	printf("res: %g %g\n",(double)(psr[p].obsn[i].sat-psr[p].param[param_pepoch].val[0]),(double)psr[p].obsn[i].residual);
      printf("Calculating sigmaz\n");
      calcSigmaz(psr[p],weights,tau[p],szbias[p],e1[p],e2[p],&nval[p],mintau);
      printf("Completed calculating sigmaz\n");
    }
  if (listing==1)
    {
      for (p=0;p<*npsr;p++)
	{
	  printf("Displaying p = %d, nval = %d\n",p,nval[p]);
	  for (i=0;i<nval[p];i++)
	    printf("%d %d %g %g\n",p,i,tau[p][i],szbias[p][i]);
	}
    }
  cpgbeg(0,grDev,1,1);
  //  cpgbeg(0,"/xs",1,1);
  cpgscf(2);
  cpgslw(2);
  cpgsch(1.4);

  doplot(psr,*npsr,grDev,mint,maxt,minsz,maxsz,style,average,tau,szbias,e1,e2,nval,nWhite,white,cline,
	 slopes,bound);

  if (calculateGWlimit==1)
    {
      calculateGWlim(psr,&idum,szbias[0][0],parFile,timFile,weights,mintau,gwamp,mint,maxt,minsz,maxsz,szbias[0],e1[0],e2[0],nit,ngw,lowAmp,hiAmp);
    }
  if (simWhite==1)
    simWhiteFunc(psr,&idum,parFile,timFile,weights,mintau);


  cpgend();
  return 0;
}

void simWhiteFunc(pulsar *psr,long *idum,char parFile[MAX_PSR_VAL][MAX_FILELEN],char timFile[MAX_PSR_VAL][MAX_FILELEN],int weights,double mintau)
{
  long double sat0[MAX_OBSN];
  int i,j,it;
  double tau[1000],szbias[1000],e1[1000],e2[1000];
  double avszbias[1000];
  float fx[1000],fy[1000];
  int nit=100;
  int nval;
  long double min,max;

  for (i=0;i<1000;i++)
    avszbias[i]=0.0;

  // Hack it
  min = psr[0].obsn[0].sat;
  max = psr[0].obsn[psr[0].nobs-1].sat;
  for (i=0;i<psr[0].nobs;i++)
    {
      psr[0].obsn[i].sat = min+(max-min)*(long double)i/(long double)psr[0].nobs/10.0;
    }
  

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

  for (it=0;it<nit;it++)
    {
      printf("iteration %d/%d\n",it,nit);
      psr[0].nJumps = 0;
      for(i=0;i<MAX_PARAMS;i++){
	psr[0].param[i].nLinkTo = 0;
	psr[0].param[i].nLinkFrom = 0;
      }
      readParfile(psr,parFile,timFile,1); /* Load the parameters                */
      for (i=0;i<psr[0].nobs;i++)
	psr[0].obsn[i].sat = sat0[i]+(TKgaussDev(idum)*166e-9)/86400.0;
	//	psr[0].obsn[i].sat = sat0[i]+(TKgaussDev(idum)*psr[0].obsn[i].toaErr*1.0e-6)/86400.0;
      formBatsAll(psr,1);                 /* Form the barycentric arrival times */
      formResiduals(psr,1,0);             /* Form the residuals                 */
      doFit(psr,1,0);                     /* Do the fitting                     */ 
      formBatsAll(psr,1);                 /* Form the barycentric arrival times */
      formResiduals(psr,1,0);             /* Form the residuals                 */
      calcSigmaz(psr[0],weights,tau,szbias,e1,e2,&nval,mintau);
      for (i=0;i<nval;i++)
	avszbias[i]+=szbias[i];
    }  
  for (i=0;i<nval;i++)
    {
      fx[i] = (float)tau[i];
      fy[i] = (float)avszbias[i]/(float)nit;
    }
  cpgline(nval,fx,fy);

}

void calculateGWlim(pulsar *psr,long *idum,double obs,char parFile[MAX_PSR_VAL][MAX_FILELEN],char timFile[MAX_PSR_VAL][MAX_FILELEN],int weights,double mintau,long double gwamp,float mint,float maxt,float minsz,float maxsz,double *szbiasObs,double *eObs1,double *eObs2,int nit,int ngw,long double lowAmp,long double hiAmp)
{
  int i,j,k,p;
  long double a;
  long double alpha = -2.0/3.0,toffset;
  long double kp[3];
  long double flo,fhi;
  long double res[MAX_OBSN],mean;
  double dist[MAX_PSR];
  double weight[MAX_OBSN];
  int addWhite=1;
  int fast=0;
  long storeSeed = *idum;
  long double sat0[MAX_OBSN];
  float my[MAX_OBSN];
  gwSrc *gw;
  double tau[100],**szbias,e1[100],e2[100];
  float fx[100],fy[100],fy95[100],fy5[100],fe1[100],fe2[100],fy2[100];
  float yval[nit];
  int nval;
  int it=0;
  int above=0;
  int doSearch=0;

  if (hiAmp!=0.0) doSearch=1;

  cpgask(0);
  szbias = (double **)malloc(sizeof(double *)*nit);
  for (i=0;i<nit;i++)
    szbias[i] = (double *)malloc(sizeof(double)*100);

  printf("Calculating limit. Observation = %g\n",obs);
  printf("Forming ideal sats ...\n");
  // Form ideal SATs
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

  writeTim("ideal.tim",&psr[0],"tempo2");  

  psr[0].nJumps = 0;
  for(i=0;i<MAX_PARAMS;i++){
    psr[0].param[i].nLinkTo = 0;
      psr[0].param[i].nLinkFrom = 0;
  }
  readParfile(psr,parFile,timFile,1); /* Load the parameters       */
  
  printf("GW background calculation ...\n");
  a = (long double)gwamp*pow(86400.0*365.25,alpha);
  hiAmp *=  pow(86400.0*365.25,alpha);
  lowAmp *= pow(86400.0*365.25,alpha);

  if((gw = (gwSrc *)malloc(sizeof(gwSrc)*ngw))==NULL){
    printf("Unable to allocate memory for gwSrc.\n");
    exit(1);
  }
  // Setup vector pointing at each pulsar
  dist[0] =  3.08568025e19; // 1 kpc in m
  setupPulsar_GWsim(psr[0].param[param_raj].val[0],
		    psr[0].param[param_decj].val[0],kp);
  flo = 1.0L/(30*365.25*86400.0L);
  fhi = 1.0L/(2.0*86400.0L);

  toffset = psr[0].param[param_pepoch].val[0];

  printf("A = %g, ngw = %d, flo = %g, fhi = %g, alpha = %g\n",(double)a,ngw,(double)flo,(double)fhi,
	 (double)alpha);
  do {
    if (doSearch==1)
      {
	a = (hiAmp+lowAmp)/2.0;
	above=0;
	printf("Searching with a = %g\n",(double)a);
      }
    for (it=0;it<nit;it++)
      {
	//      printf("Calc backgroud\n");
	GWbackground(gw,ngw,idum,flo,fhi,a,alpha,1);
	//printf("done Calc backgroud\n");
	mean=0.0L;
	//printf("Calc residuals 1\n");
	for (j=0;j<psr[0].nobs;j++)
	  {
	    res[j]=0.0L;
	    for (k=0;k<ngw;k++)
	      res[j]+=calculateResidualGW(kp,&gw[k],
					  (psr[0].obsn[j].sat-toffset)*86400.0L,
					dist[0]);	  
	  mean+=res[j];
	}
      //printf("Calc residuals 2\n");
      for (j=0;j<psr[0].nobs;j++)
	{
	  psr[0].obsn[j].sat = (sat0[j]+(res[j]-mean/psr[0].nobs)/86400.0);
	  if (addWhite>0)
	    {
	      psr[0].obsn[j].sat += (long double)((addWhite*TKgaussDev(idum)*psr[0].obsn[j].toaErr*1e-6)/86400.0);
	    }
	}
      //printf("Done calc residuals\n");
      psr[0].nJumps = 0;
      for(i=0;i<MAX_PARAMS;i++){
	psr[0].param[i].nLinkTo = 0;
	psr[0].param[i].nLinkFrom = 0;
      }
      readParfile(psr,parFile,timFile,1); /* Load the parameters                */
      //printf("Getting bats\n");
      if (it==0 || fast==0)
	{
	  formBatsAll(psr,1);                 /* Form the barycentric arrival times */
	  formResiduals(psr,1,0);             /* Form the residuals                 */
	}
      else
	{
	  vectorPulsar(psr,1);   
	  calculate_bclt(psr,1); 
	  formBats(psr,1);       
	  formResiduals(psr,1,0); 
	}
      //printf("Doing fit\n");
      doFit(psr,1,0);                     /* Do the fitting                     */
      //printf("Getting postfit\n");
      if (it==0 || fast==0)
	{
	  formBatsAll(psr,1);                 /* Form the barycentric arrival times */
	  formResiduals(psr,1,0);             /* Form the residuals                 */
	}
      else
	{
	  vectorPulsar(psr,1);   /* 1. Form a vector pointing at the pulsar */
	  calculate_bclt(psr,1);           /* 3. Calculate bclt (WHAT IS THIS?) */
	  formBats(psr,1);                   /* Form Barycentric arrival times */
	  formResiduals(psr,1,0);   /* Form the residuals                 */	  
	}
      //printf("calculating sigmaz\n");
      calcSigmaz(psr[0],weights,tau,szbias[it],e1,e2,&nval,mintau);
      //printf("Doing plots\n");
      TKconvertFloat1(tau,fx,nval);
      if (it==0)
	{
	  double pv[4],v[4];
	  int ma=3;

	  writeTim("sigmaz_gwsim.tim",&psr[0],"tempo2");  
	  // Model the actual observations
	  TKleastSquares_svd_noErr(tau,szbiasObs,nval,pv,ma,TKfitPoly);
	  for (i=0;i<nval;i++)	    
	    {
	      my[i] = 0.0;
	      TKfitPoly(tau[i],v,ma);
	      for (j=0;j<ma;j++)
		my[i]+=v[j]*pv[j];
	    }
	}
      // Calculate mean szbias
      if (it%100 == 0 || it == nit-1)
	{
	  for (i=0;i<nval;i++)
	    {
	      fy[i] = 0.0;
	      for (j=0;j<it+1;j++)
		{
		  fy[i]  += (float)szbias[j][i];
		  yval[j] = (float)szbias[j][i];
		  //	      printf("Have %d %d %g\n",i,j,szbias[j][i]);
		}
	      fy[i] /= (float)(it+1);
	      TKsort_f(yval,it+1);
	      fy95[i] = yval[(int)(95.0/100.0*(it+1)+0.5)-1];
	      fy5[i] = yval[(int)(5.0/100.0*(it+1)+0.5)-1];
	    }
	  //      TKconvertFloat1(szbias,fy,nval);
	  cpgenv(log10(mint),log10(maxt),minsz,maxsz,0,30);
	  cpglab("\\gt (yr)","\\gs\\dz\\u","");
	  cpgline(nval-1,fx,fy);
	  cpgsci(7); cpgline(nval-1,fx,my); cpgsci(1);
	  cpgsls(4); cpgline(nval-1,fx,fy95);
	  cpgline(nval-1,fx,fy5); cpgsls(1);
	  TKconvertFloat1(szbiasObs,fy2,nval);
	  cpgpt(nval-1,fx,fy2,5);
	  TKconvertFloat1(eObs1,fe1,nval);
	  TKconvertFloat1(eObs2,fe2,nval);
	  cpgpt(nval-1,fx,fy2,5);
	  for (i=0;i<nval-1;i++)
	    {
	      cpgerr1(2,fx[i],fy2[i],fe1[i],1); 
	      cpgerr1(2,fx[i],fy2[i],fe2[i],1); 
	    }
	}
      if (szbias[it][0] > obs) above++;
      printf("unweighted percentage above observed statistic = %d/%d = %g\% (total number of iterations = %d, A = %g)\n",above,(it+1),above/(double)(it+1)*100.0,nit,(double)(a/pow(86400.0*365.25,alpha)));
      //      printf("Looping\n");
    }

  // Now calculate weighting function
  for (i=0;i<nval-1;i++)
    {
      weight[i] = pow(10,fy[i])/pow(10,my[i])/pow(10,my[i]);
      printf("weights: %g %g %g %g\n",fx[i],weight[i],fy[i],my[i]);
    }
  
  // Now recalculate the statistics
  {
    double measuredStat=0.0;
    double simStat=0.0;
    int nuse;
    int nmax;

    if (doSearch==0)
      {
	nmax = nval;
	for (nuse=1;nuse<nmax;nuse++)
	  {
	    above=0;
	    measuredStat=0.0;
	    for (i=0;i<nmax;i++)
	      measuredStat += (((szbiasObs[i]))*weight[i]);
	    for (it=0;it<nit;it++)
	      {
		simStat=0.0;
		for (i=0;i<nmax;i++) simStat+=(((szbias[it][i]))*weight[i]);
		//	    printf("simStat = %g %g\n",simStat,measuredStat);
		if (simStat > measuredStat) above++;
	      }  
	    printf("nuse = %d, Detected %d out of %d = %g percent.  Measured stat = %g\n",nuse,above,nit,(double)above/(double)nit*100.0,measuredStat);
	    //	exit(1);
	  } 
      }
  }
  if (doSearch==1)
    {
      if ((double)above/(double)nit*100.0 > 95.2)
	hiAmp = a;
      if ((double)above/(double)nit*100.0 < 94.8)
	lowAmp = a;
    }
  } while (doSearch!=0 && ((double)above/(double)nit*100.0 > 95.2 || (double)above/(double)nit*100.0 < 94.8));
  // Now output results
  {
    FILE *fout;
    fout = fopen("GWlimit_results.dat","a");
    fprintf(fout,"%g %g %d %d %g %s %s %d\n",(double)a/pow(86400.0*365.25,alpha),(double)alpha,ngw,nit,(double)above/(double)nit*100.0,parFile[0],timFile[0],storeSeed);
    fclose(fout);
  }

  // Now make a postscript plot
  cpgend();
  cpgbeg(0,"sigmaz.ps/ps",1,1);
  cpgslw(4);
  cpgsch(1.4);
  cpgsfs(2);
  cpgenv(log10(mint),log10(maxt),minsz,maxsz,0,30);
  cpglab("\\gt (yr)","\\gs\\dz\\u","");
  cpgline(nval-1,fx,fy);
  cpgsls(4); cpgline(nval-1,fx,fy95);
  cpgline(nval-1,fx,fy5); cpgsls(1);
  TKconvertFloat1(szbiasObs,fy,nval);
  cpgpt(nval-1,fx,fy,5);
  TKconvertFloat1(eObs1,fe1,nval-1);
  TKconvertFloat1(eObs2,fe2,nval-1);
  cpgpt(nval-1,fx,fy2,5);
  for (i=0;i<nval-1;i++)
    {
      cpgerr1(2,fx[i],fy2[i],fe1[i],1); 
      cpgerr1(2,fx[i],fy2[i],fe2[i],1); 
    }
  cpgend();
}


// If sameBackground = 1 then the same background is used for all pulsar residuals
void convert_gravWaveBackground_fit(pulsar *psr,int npsr,double convertGW,long *idum,int sameBackground,  
				    char parFile[MAX_PSR_VAL][MAX_FILELEN],char timFile[MAX_PSR_VAL][MAX_FILELEN])
{
  int i,j,k,p;
  long double a;
  long double alpha = -2.0/3.0,toffset;
  int ngw = 10000;
  long double kp[MAX_PSR][3];
  long double flo,fhi;
  long double res[MAX_OBSN],mean;
  double dist[MAX_PSR];
  gwSrc *gw;
  int addWhite=0;


  printf("Forming ideal sats ...\n");
  // Form ideal SATs
  for (j=0;j<5;j++)
    {
      for(p=0;p<npsr;p++){
	psr[p].nJumps = 0;
	for(i=0;i<MAX_PARAMS;i++){
	  psr[p].param[i].nLinkTo = 0;
	  psr[p].param[i].nLinkFrom = 0;
	}
      }
      readParfile(psr,parFile,timFile,npsr); /* Load the parameters       */
      formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,npsr,0);    /* Form the residuals                 */
      for (p=0;p<npsr;p++)
	{
	  for (i=0;i<psr[p].nobs;i++)
	    psr[p].obsn[i].sat -= (long double)psr[p].obsn[i].residual/86400.0L;
	}
    }
  writeTim("ideal.tim",&psr[0],"tempo2");  

  if (addWhite==1)
    {
      for (p=0;p<npsr;p++)
	{
	  for (i=0;i<psr[p].nobs;i++)
	    psr[p].obsn[i].sat += (long double)(psr[p].obsn[i].toaErr*1.0e-6*TKgaussDev(idum))/86400.0;
	}
      writeTim("idealPlusWhite.tim",&psr[0],"tempo2");  
    }
  for(p=0;p<npsr;p++){
    psr[p].nJumps = 0;
    for(i=0;i<MAX_PARAMS;i++){
      psr[p].param[i].nLinkTo = 0;
      psr[p].param[i].nLinkFrom = 0;
    }
  }
  readParfile(psr,parFile,timFile,npsr); /* Load the parameters       */
  
  printf("GW background calculation ...\n");
  a = (long double)convertGW*pow(86400.0*365.25,alpha);

  if((gw = (gwSrc *)malloc(sizeof(gwSrc)*ngw))==NULL){
    printf("Unable to allocate memory for gwSrc.\n");
    exit(1);
  }
  // Setup vector pointing at each pulsar
  for (p=0;p<npsr;p++)
    {
      dist[p] =  3.08568025e19; // 1 kpc in m
      setupPulsar_GWsim(psr[p].param[param_raj].val[0],
			psr[p].param[param_decj].val[0],kp[p]);
    }
  flo = 1.0L/(30*365.25*86400.0L);
  fhi = 1.0L/(2.0*86400.0L);

  toffset = psr[0].param[param_pepoch].val[0];
  if (sameBackground==1)
    {
      printf("A = %g, ngw = %d, flo = %g, fhi = %g, alpha = %g\n",(double)a,ngw,(double)flo,(double)fhi,
	     (double)alpha);
      GWbackground(gw,ngw,idum,flo,fhi,a,alpha,1);
      for (p=0;p<npsr;p++)
	{
	  mean=0.0L;
	  for (j=0;j<psr[p].nobs;j++)
	    {
	      res[j]=0.0L;
	      for (k=0;k<ngw;k++)
		res[j]+=calculateResidualGW(kp[p],&gw[k],
					 (psr[p].obsn[j].sat-toffset)*86400.0L,
					    dist[p]);	  
	      mean+=res[j];
	    }
	  for (j=0;j<psr[p].nobs;j++)
	    {
	      psr[p].obsn[j].sat += (res[j]-mean/psr[p].nobs)/86400.0;
	    }
	}
      formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,npsr,0);    /* Form the residuals                 */
      doFit(psr,npsr,0);   /* Do the fitting     */
      formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,npsr,0);    /* Form the residuals                 */
      for (i=0;i<psr[0].nobs;i++)
	printf("res: %g %g\n",(double)psr[0].obsn[i].sat,(double)psr[0].obsn[i].residual);
      textOutput(psr,1,0,0,0,0,"");
    }
  writeTim("gw.tim",&psr[0],"tempo2");
}

// If sameBackground = 1 then the same background is used for all pulsar residuals
void convert_gravWaveBackground_noFit(pulsar *psr,int npsr,double convertGW,long *idum,int sameBackground)
{
  int i,j,k,p;
  long double a;
  long double alpha = -2.0/3.0,toffset;
  int ngw = 10000;
  long double kp[MAX_PSR][3];
  long double flo,fhi;
  double res,mean;
  double dist[MAX_PSR];
  gwSrc *gw;

  printf("GW background calculation ...\n");
  a = (long double)convertGW; //*pow(86400.0*365.25,alpha);

  if((gw = (gwSrc *)malloc(sizeof(gwSrc)*ngw))==NULL){
    printf("Unable to allocate memory for gwSrc.\n");
    exit(1);
  }
  // Setup vector pointing at each pulsar
  for (p=0;p<npsr;p++)
    {
      dist[p] =  3085.68025e19; // 1000 kpc in m
      setupPulsar_GWsim(psr[p].param[param_raj].val[0],
			psr[p].param[param_decj].val[0],kp[p]);
    }
  flo = 1.0L/(100*365.25*86400.0L);
  fhi = 2.0L/(86400.0L);

  toffset = psr[0].param[param_pepoch].val[0];
  printf("In here with a = %g %d\n",(double)a,sameBackground);
  if (sameBackground==1)
    {
      GWbackground(gw,ngw,idum,flo,fhi,a,alpha,1);
      for (p=0;p<npsr;p++)
	{
	  mean=0.0;
	  for (j=0;j<psr[p].nobs;j++)
	    {
	      res=0.0;
	      for (k=0;k<ngw;k++)
		res+=calculateResidualGW(kp[p],&gw[k],
					 (psr[p].obsn[j].sat-toffset)*86400.0L,
					 dist[p])/86400.0L;	  
	      psr[p].obsn[j].residual = res;
	      mean+=res;
	    }
	}

	      printf("Haveres %Lg %g\n",psr[p].obsn[j].sat,res);

    }
}

void doplot(pulsar *psr,int npsr,char *grDev,float mint,float maxt,float minsz,float maxsz,
	    int style,int average,double tau[MAX_PSR_VAL][100],double szbias[MAX_PSR_VAL][100],
	    double e1[MAX_PSR_VAL][100],double e2[MAX_PSR_VAL][100],int nval[MAX_PSR_VAL],
	    int nWhite,float *white,char *cline,int slopes,int bound)
{
  float avX[100],avY[100];
  int   nav[100];
  int   i,j,p;
  float plotX[MAX_OBSN],plotY[MAX_OBSN],yerr1[MAX_OBSN],yerr2[MAX_OBSN];
  float px[5],py[5];

  for (i=0;i<average;i++)
    {
      avX[i]=log10(mint)+(log10(maxt)-log10(mint))*(float)i/(float)average;
      avY[i]=0.0;
      nav[i]=0;
    }


    cpgenv(log10(mint),log10(maxt),minsz,maxsz,0,30);
    /*  cpgsvp(0.1,0.95,0.15,0.90);
  cpgswin(log10(mint*365.25),log10(maxt*365.25),minsz,maxsz);
  cpgbox("LCMS",0.0,0,"",0.0,0);
  cpgswin(log10(mint),log10(maxt),minsz,maxsz);
  cpgbox("BLNS",0.0,0,"BCNSV",0.0,0);*/
    //  cpglab("\\gt (yr)","log\\d10\\u[\\gs\\dz\\u]","\\gt (day)");
  cpglab("\\gt (yr)","\\gs\\dz\\u","");
  if (style==0){
  cpgsch(0.5); cpglab("","",cline); cpgsch(1.4);
  }

  if (bound==1)
    {
      // Plot typical clock stability (from nist-ptb)
      px[0]=-3; py[0] =-14.6;
      px[1]=2;  py[1] =-14.6;
      px[2]=px[1]; py[2] = -20;
      px[3]=px[0]; py[3] = -20;
      cpgsci(1); cpgline(2,px,py);
      cpgsci(7); cpgsfs(4); cpgpoly(4,px,py); cpgsfs(1); cpgsci(1);

      /* Expected GW background: Omega_g h^2 = 10^-9 */
      px[0] = -3;
      py[0] = -15.2253+0.91397*px[0];
      px[1] = 2;
      py[1] = -15.2253+0.91397*px[1];
      px[2] = 2;
      py[2] = -20;
      px[3] = -3;
      py[3] = -20;
      cpgsls(1); cpgline(2,px,py); cpgsls(1);
      cpgsci(7); cpgsfs(4); cpgpoly(4,px,py); cpgsfs(1); cpgsci(1);

      /* Plot 50ns white noise for 0437 sampling*/
      px[0] = -3;
      py[0] = log10((pow(pow(10,px[0]),-3.0/2.0))*50e-15/100);
      px[1] = 2;
      py[1] = log10((pow(pow(10,px[1]),-3.0/2.0))*50e-15/100);
      px[2] = -3;
      py[2] = py[1];
      cpgsls(1); cpgline(2,px,py); cpgsls(1);
      cpgsci(7); cpgsfs(4); cpgpoly(3,px,py); cpgsfs(1); cpgsci(1);
    }




  for (p=0;p<npsr;p++)
    {
      if (average > 0)
	{
	  for (i=0;i<nval[p];i++)
	    {
	      for (j=0;j<average-1;j++)
		{
		  if (tau[p][i] >= avX[j] && tau[p][i] < avX[j+1])
		    {
		      avY[j]+=szbias[p][i];
		      nav[j]++;
		    }
		}
	    }      
	}
    }
  // Now do the plot
  for (p=0;p<npsr;p++)
    {
      for (i=0;i<nval[p];i++)
	{
	  plotX[i] = (float)tau[p][i];
	  plotY[i] = (float)szbias[p][i];
	  yerr1[i] = (float)e1[p][i];
	  yerr2[i] = (float)e2[p][i];
	}
      if (style==0)
	{ 
	  cpgpt(nval[p],plotX,plotY,5);
	  for (i=0;i<nval[p];i++)
	    {
	      cpgerr1(2,plotX[i],plotY[i],yerr1[i],1); 
	      cpgerr1(2,plotX[i],plotY[i],yerr2[i],1); 
	    }
	}
      else if (style==1)
	{
	  if (p>2)
	    cpgsci(2+p);
	  else
	    cpgsci(1+p);
	  if (p<2)
	    cpgpt(nval[p],plotX,plotY,5);
	  else
	    cpgpt(nval[p],plotX,plotY,16);
	  for (i=0;i<nval[p];i++)
	    {
	      cpgerr1(2,plotX[i],plotY[i],yerr1[i],1); 
	      cpgerr1(2,plotX[i],plotY[i],yerr2[i],1); 
	    }
	  if (p>1)
	    cpgline(nval[p],plotX,plotY);
	  cpgsci(1);
	}
      else if (style==4)
	{
	  cpgsci(1+p);
	  cpgpt(nval[p],plotX,plotY,16);
	  cpgline(nval[p],plotX,plotY);
	  for (i=0;i<nval[p];i++)
	    {
	      cpgerr1(2,plotX[i],plotY[i],yerr1[i],1); 
	      cpgerr1(2,plotX[i],plotY[i],yerr2[i],1); 
	    }
	}
      else if (style==2)
	{
	  if (p==0)
	    cpgpt(nval[p],plotX,plotY,18);
	  else if (p==1)
	    cpgpt(nval[p],plotX,plotY,17);
	  else if (p==2)
	    cpgpt(nval[p],plotX,plotY,21);
	  else
	    cpgpt(nval[p],plotX,plotY,5+p);
	  for (i=0;i<nval[p];i++)
	    {
	      if (fabs(yerr2[i]-yerr1[i]) > 0.1)
		{		  
		  cpgerr1(2,plotX[i],plotY[i],yerr1[i],1); 
		  cpgerr1(2,plotX[i],plotY[i],yerr2[i],1); 
		}
	    }
	  calcSpline(plotX,plotY,nval[p]);
	}
      else if (style==3)
	{
	  //	  if (p>2)
	  //	    cpgsci(2+p);
	  //	  else
	  //	    cpgsci(1+p);
	  if (p<2)
	    {
	      if (p==0)
		cpgpt(nval[p],plotX,plotY,5);
	      else if (p==1)
		cpgpt(nval[p],plotX,plotY,6);
	    }
	  else
	    cpgpt(nval[p],plotX,plotY,16);
	  for (i=0;i<nval[p];i++)
	    {
	      cpgerr1(2,plotX[i],plotY[i],yerr1[i],1); 
	      cpgerr1(2,plotX[i],plotY[i],yerr2[i],1); 
	    }
	  if (p>1)
	    cpgline(nval[p],plotX,plotY);
	  cpgsci(1);
	}
    }
  if (style==1)
    {
      for (p=0;p<npsr;p++)
	{
	  if (p>2)
	    cpgsci(2+p);
	  else
	    cpgsci(1+p);


	  cpgsch(0.8); cpgtext(log10(maxt)-0.3*(log10(maxt)-log10(mint)),maxsz-(p+1)*0.05*(maxsz-minsz),psr[p].name); cpgsch(1.4);
	}
      cpgsci(1);
    }
  if (style==0)
    {
      char allnames[1000]="";
      for (p=0;p<npsr;p++)
	{
	  strcat(allnames,psr[p].name);
	  if (p!=npsr-1) strcat(allnames,", ");
	}
      cpgsch(1); cpgsci(1); cpgtext(-2.5,-8,allnames); cpgsch(1.4);
    }
  //  cpgsch(1); cpgptxt(1.1,-14.52,24,0.0,"\\gW\\dg\\uh\\u2\\d = 10\\u-9\\d"); cpgsch(1.4);
  
  /* plot white noise levels */
  /* The scalings have been done from simulation -- should do analytically */
  for (i=0;i<nWhite;i++)
    {
      px[0] = -3;
      py[0] = log10((pow(pow(10,px[0]),-3.0/2.0))*6e-15/100e-9*white[i]);
      px[1] = 2;
      py[1] = log10((pow(pow(10,px[1]),-3.0/2.0))*6e-15/100e-9*white[i]);
      cpgsls(4); cpgline(2,px,py); cpgsls(1);
    }
  
  /* Plot average */
  if (average > 0)
    {
      float fx[1000],fy[1000];
      int n=0;

      for (i=0;i<average;i++)
	{
	  if (nav[i]>1)
	    {
	      fx[n] = avX[i]+(avX[1]-avX[0])/2.0;
	      fy[n] = avY[i]/(float)nav[i];
	      n++;
	    }
	}
      cpgline(n,fx,fy);
      cpgpt(n,fx,fy,9);
    }


  // plot slopes
  if (slopes==1)
    {
      int a;
      /*      for (a=-6;a<=0;a+=2)
	{
	  px[0] = -3; py[0] = log10((pow(pow(10,px[0]),-(a+3.0)/2.0))*1e-10);
	  px[1] = 2;  py[1] = log10((pow(pow(10,px[1]),-(a+3.0)/2.0))*1e-10);
	  cpgsci(14); cpgsls(3); cpgline(2,px,py); cpgsls(1); cpgsci(1);
	  }          */
      cpgslw(3);
      /*      px[0] = -3; py[0] = log10(0.26*(pow(pow(10,px[0]),-(-5+3.0)/2.0))*1e-10);
      px[1] = 2;  py[1] = log10(0.26*(pow(pow(10,px[1]),-(-5+3.0)/2.0))*1e-10);
      cpgsci(14); cpgsls(3); cpgline(2,px,py); cpgsls(1); cpgsci(1);

      px[0] = -3; py[0] = log10(2*(pow(pow(10,px[0]),-(-1+3.0)/2.0))*1e-10);
      px[1] = 2;  py[1] = log10(2*(pow(pow(10,px[1]),-(-1+3.0)/2.0))*1e-10);
      cpgsci(14); cpgsls(2); cpgline(2,px,py); cpgsls(1); cpgsci(1);*/

      px[0] = -3; py[0] = log10(2*(pow(pow(10,px[0]),-(0+3.0)/2.0))*5e-15);
      px[1] = 2;  py[1] = log10(2*(pow(pow(10,px[1]),-(0+3.0)/2.0))*5e-15);
      cpgsci(14); cpgsls(4); cpgline(2,px,py); cpgsls(1); cpgsci(1);
      cpgslw(1);
    }
}

/* Direct copy of the original Fortran */
/* using identical fitting routines etc. */
/* Note that all arrays start from 1 */
//void calcSigmaz(pulsar *psr,int npsr,int time,int psrfile,double mintau,
//		int nSel,int *selPsr,int weights,int style,int average,
//		float *avX,float *avY,int *nav)
void calcSigmaz(pulsar psr,int weights,double *ret_tau,double *ret_szbias,double *ret_e1,
		double *ret_e2,int *ret_nval,double mintau)
{
  int n,nbad,nfit,nread,idiag,niter,ntable,ios,nsigs,nuncor,ntot;
  int ndostats,nptot,nuncorbins,ntries;
  double x,szlog,sz,cov4,rootsz,p4,sumwt,sumwtuncor,error,scale,bias;
  double utjd2nuncor,wt,dplus,dminus,errlog,chidf,avewt,scalelog;
  double taulog1,addvsq,szbias,dp,dm,tabledata[101];
  int endit;
  int   j;

  *ret_nval=0;
  nxunits = 0;
  nusewt = weights; 
  tmin =-1e50;
  tmax = 1.0e50;
  xmin =-1.0e50;
  xmax = 1.0e50;
  root2 = sqrt(2.0);
  secyear = 86400.0*365.25;
  utfirst = 1.0e50;
  utlast = -1.0e50;
  npt1last = 0;
  npt2last = 0;
  utjd1 = 0.0;
  ndim = 0;
  if (nxunits >= 0) ndim=4;
  if (nxunits <= -1 && nxunits >= -3) ndim=3;
  if (nxunits <= -4 && nxunits >= -6) ndim=1;

  int nobs = psr.nobs;
  double t2times[MAX_OBSN],t2resid[MAX_OBSN],t2error[MAX_OBSN],t0;
  int i;
  scale = 2.0*sqrt(5.0);
  scalelog = log10(scale);
  for (n=1;n<=100;n++)
    tabledata[n] = 0.0;
  
  readin(psr);
  tausec = taumax;
  tauday = tausec/86400.0;
  ntau = 0;
  ncubics=0;
  if (nusewt >= 0)
    {
      tauyear = tauday/365.25;
      utjd1 = utjd[indx[1]]-1.0e-2;
      utjd2 = utjd[indx[npt]]+1.0e-2;
      utmean = (utjd1+utjd2)/2.0;
      for (niter=1;niter<=nusewt+1;niter++)
	{
	  fit4(&nfit,&p4,&cov4,ndostats,&chidf,&avewt);
	  addvsq=addvar*addvar+(chidf-1.0)/avewt;
	  if (addvsq > 0.0)
	    addvar=sqrt(addvsq);
	  else
	    {
	      if (nusewt != 0 && nusewt != -3) addvar = 0.0;
	    }
	}
    }
  // Should check nbintype
  ndostats=0;

  // 110 continue statement
  do {
    ntau++;
    ncubic = 0;
    tauday = tausec/86400.0;
    tauyear = tauday/365.25;
    taulog = log10(tausec);
    sumwt = 0.0;
    sz=0.0;
    utjd2nuncor = 0.0;
    sumwtuncor = 0.0;
    nuncorbins = 0;
    ntot = 0;
    ntries = 0;
    nptot = 0;
    npt1last = 0;
    npt2last = 0;
    
    if (nbintype==0) utjd1=utjd[indx[1]]-tausec;
    if (nbintype==1) utjd1=utjd[indx[1]]-taumin;
    if (nbintype==2) utjd1=utjd[indx[1]]-tausec;
    endit=0;

    do {
      //      printf("utjd1 = %g\n",utjd1);
      if (nbintype==0) utjd1+=tausec;
      if (nbintype==1) utjd1+=taumin;
      if (nbintype==2) utjd1+=tausec;
      if (utjd1 < utjd[indx[npt]])
	{
	  utjd2=utjd1+tausec;
	  utmean=(utjd1+utjd2)/2.0;
	  if (utjd1 >= utjd2nuncor)
	    {
	      nuncor=1;
	      utjd2nuncor=utjd2;
	    }
	  else
	    {
	      // Should check overlap
	      nuncor=0;
	    }
	  ncubics++;		  
	  fit4(&nfit,&p4,&cov4,ndostats,&chidf,&avewt);
	  ntries++;
	  if (nfit==1)
	    {
	      wt=1.0/cov4;
	      sz=sz+wt*p4*p4;
	      sumwt=sumwt+wt;
	      ntot++;
	      nptot=nptot+npt2last-npt1last+1;
	      if (nuncor==1)
		{
		  sumwtuncor=sumwtuncor+wt;
		  nuncorbins=nuncorbins+1;
		}	     
	      //	      printf("Here %g %g %d %d\n",(double)sumwtuncor,(double)wt,(int)nuncorbins,nuncor);
	    }
	}
      else
	endit=1;
    } while (endit==0); // 130 loop
    //    printf("Add or not (Exit): %d %g %g\n",ntot,sz,sumwt);
    if (ntot==0 || sz <=0.0 || sumwt<=0.0)
      {
      }
    else
      {
	rootsz = tausec*tausec*sqrt(sz/sumwt)/scale;
	szlog = log10(rootsz);
	errlog=0.0;
	bias=0.0;
	dp=0.0;
	dm=0.0;
	if (sumwtuncor > 1.0e-50)
	  {
	    error=nuncorbins/sumwtuncor;
	    if (error > 0.0) errlog=2.0*taulog+0.5*log10(error)-scalelog;
	    x=nuncorbins;
	    bias=0.17/x;
	    dm = 0.31/sqrt(x);
	    if (nuncorbins > 1)
	      dp=0.31/sqrt(x-1.0);
	    else
	      dp=0.52;
	    
	  }
	szbias = szlog-bias;
	dplus = szbias + dp;
	dminus = szbias - dm;

	ret_szbias[*ret_nval] = szbias;
	ret_tau[*ret_nval]    = log10(tauyear);
	ret_e1[*ret_nval]     = dp;
	ret_e2[*ret_nval]     = -dm;
	//	printf("Have %g %g %g %g\n",szbias,log10(tauyear),dp,-dm);

	//	exit(1);
	(*ret_nval)++;
      }
    if (nbintype==0) tausec=tausec/2.0;
    if (nbintype==1) tausec=tausec-taumin;
    if (nbintype==2) tausec=tausec+taumin;
    //    exit(1);
    //    printf("Exited (testing) with %d %g %g %g %g %d\n",ntot,tausec,tauensure,tausec/86400.0/365.25,mintau,*ret_nval);
    //	fflush(stdout);
  } while (ntot > 0 || (tausec >= tauensure && tauensure >= 0.0) || (tausec/86400.0/365.25 > mintau));
  //  printf("Exited with %d %g %g %g %g %d\n",ntot,tausec,tauensure,tausec/86400.0/365.25,mintau,*ret_nval);
  //	fflush(stdout);
}



/* WHAT AM I RETURNING IN AND OUT OF THIS FUNCTION? */
void fit4(int *nfit,double *p4,double *cov4,int ndostats,double *chidf,double *avewt)
{
  int n,i,j,nbad,nok,nf,npt1,npt2,m,k;
  static int npass=0;
  static int ntaulast=-1;
  double am[21][21],aminv[21][21],fvec[21],determ,y,ylast,x,ave,rms,alv,yf;
  double res,sumwt,par[5],wt,chisqu,degfree,utave,err4;
  static double scalefact=1.0;
  static double scalecon=1.0;

  npass++;
  ncubic++;
  ncubics++;
  if (ndim > 1) scalefact = pow(86400.0,ndim-1);
  alv = 0.0;
  for (i=1;i<=ndim;i++)
    {
      fvec[i] = 0.0;
      par[i] = 0.0;
      for (j=1;j<=ndim;j++)
	am[i][j]=0.0;
    }
  sumwt=0.0;
  npt1=0;
  npt2=0;
  *nfit=0;
  if (ndim != 4)
    {
      if (ntunits >= -2)
	scalecon=3.0;
      else
	{
	  if (ndim==3) scalecon = 3.0*tdiffmin;
	  if (ndim==2) scalecon = 6.0*tdiffmin*tdiffmin;
	  if (ndim==1) scalecon = 6.0*tdiffmin*tdiffmin*tdiffmin;
	}
    }
  for (n=1;n<=npt;n++)
    {
      if (utjd[indx[n]] >= utjd1-1.0e-5)
	{
	  if (utjd[indx[n]]+(4-ndim)*tdiffmin > utjd2+1.0e-5) goto pos90;
	  if (npt2 == 0) npt1=n;
	  npt2=n;
	}
    }
 pos90:
  if (npt1 == 0 || npt2 == 0 || npt2-npt1+1 < ndim || npt2 == 0 || (npt1 == npt1last && npt2 == npt2last) 
      || (utjd[indx[npt2]]+(4.0-ndim)*tdiffmin-utjd[indx[npt1]] < tausec/root2))
    {
      goto pos1000;
    }
  npt1last = npt1;
  npt2last = npt2;
  for (n=npt1;n<=npt2;n++)
    {
      getweights(n,&wt);
      sumwt=sumwt+wt;
      getprtj(n);
      y = data[indx[n]];
      for (i=1;i<=ndim;i++)
	{
	  fvec[i] = fvec[i] + y*prtl[i]*wt;
	  for (j=1;j<=ndim;j++)
	    {
	      am[i][j]=am[i][j]+prtl[i]*prtl[j]*wt;
	    }
	}
    }
  mat20(am,aminv,ndim,&determ,&nbad);
  if (nbad != 0.0)
    {
      goto pos1000;     
    }
  for (i=1;i<=ndim;i++)
    {
      for (j=1;j<=ndim;j++)
	par[i]=par[i]+aminv[i][j]*fvec[j];
    }
  *nfit=1;
  *p4=par[ndim]/scalefact;
  *cov4=aminv[ndim][ndim]/scalefact/scalefact;
  if (ndim!=4)
    {
      *p4=(*p4)/scalecon;
      *cov4=(*cov4)/scalecon/scalecon;
    }
  utave=(utjd1/86400.0+utjd2/86400.0)/2.0;
  err4=sqrt(*cov4);
  if (ntau!=ntaulast) ntaulast=ntau;
  ave=0.0;
  nok=0;
  sumwt=0.0;
  chisqu=0.0;
  for (n=npt1;n<=npt2;n++)
    {
      yf = 0.0;
      getweights(n,&wt);
      getprtj(n);
      for (i=1;i<=ndim;i++)
	yf=yf+prtl[i]*par[i];
      res=data[indx[n]]-yf;
      ave=ave+wt*res;
      chisqu=chisqu+wt*res*res;
      sumwt=sumwt+wt;
      nok++;
      if (nok>1) alv=alv+(res-ylast)*(res-ylast);
      ylast=res;
    }
  // 401 Continue
  if (nok > 1 && sumwt > 1.0e-50)
    {
      ave = ave/sumwt;
      rms = chisqu/sumwt-ave*ave;
      if (rms > 0.0) rms=sqrt(rms);
      if (alv > 0.0) alv=sqrt(alv/2.0/(nok-1.0));
      *avewt=sumwt/nok;
      degfree=nok-ndim;
      if (nok > 4)
	*chidf = chisqu/degfree;
      else
	*chidf=1.0;
    }
 pos1000:
  return;
}

void mat20(double sam[21][21],double a[21][21],int n,double *determ,int *nbad)
{
  int index[21][3],ipivot[21];
  double pivot[21],amax,swap,x,t;
  int i,j,k,l,l1,irow,icolum,jrow,jcolum;

  *nbad=0;
  for (i=1;i<=n;i++)
    {
      for (j=1;j<=n;j++)
	a[i][j]=sam[i][j];
    }
  *determ=0.0;
  for (j=1;j<=n;j++)
    ipivot[j]=0;
  for (i=1;i<=n;i++) // End 550
    {
      amax=0.0;
      for (j=1;j<=n;j++) // End 105
	{
	  if (ipivot[j]-1 != 0) // CHECK THIS <<<
	    {
	      for (k=1;k<=n;k++)
		{
		  if (ipivot[k]-1 < 0) goto pos_mat80;
		  if (ipivot[k]-1 == 0) goto pos_mat100;
		  if (ipivot[k]-1 > 0) goto pos_mat750;
		pos_mat80:
		  if (fabs(amax)-fabs(a[j][k]) < 0)
		    {
		      irow=j;
		      icolum=k;
		      amax=a[j][k];
		    }
		pos_mat100:
		  {
		  }
		} // 100
	    }
	} // 105
      ipivot[icolum]=ipivot[icolum]+1;
      if (irow-icolum != 0)
	{
	  for (l=1;l<=n;l++)
	    {
	      swap = a[irow][l];
	      a[irow][l]=a[icolum][l];
	      a[icolum][l]=swap;
	    }
	} // 260
      index[i][1]=irow;
      index[i][2]=icolum;
      pivot[i]=a[icolum][icolum];
      x=fabs(pivot[i]);
      if (x<=1.0e-70 || x>=1e70)
	{
	  *nbad=1;
	  return;
	}
      *determ=*determ+log10(fabs(pivot[i]));
      a[icolum][icolum]=1.0;
      for (l=1;l<=n;l++)
	a[icolum][l]=a[icolum][l]/pivot[i];
      for (l1=1;l1<=n;l1++)
	{
	  if (l1-icolum!=0)
	    {
	      t=a[l1][icolum];
	      a[l1][icolum]=0.0;
	      for (l=1;l<=n;l++)
		a[l1][l]=a[l1][l]-a[icolum][l]*t;
	    }
	}
    }// End 550

  for (i=1;i<=n;i++)
    {
      l=n+1-i;
      if ((index[l][1]-index[l][2])!=0.0)
	{
	  jrow=index[l][1];
	  jcolum=index[l][2];
	  for (k=1;k<=n;k++)
	    {
	      swap=a[k][jrow];
	      a[k][jrow]=a[k][jcolum];
	      a[k][jcolum]=swap;
	    }
	}
    }
 pos_mat750:
  return;
}

void getprtj(int n)
{
  double x;
  x = (utjd[indx[n]]-utmean)/86400.0;
  prtl[1] = 1.0;
  prtl[2] = x;
  prtl[3] = x*x;
  prtl[4] = x*x*x;
}

void getweights(int n, double *wt)
{
  if (nusewt == 0) *wt = 1.0/addvar/addvar;
  if (nusewt == -1) *wt = 1.0/pow(sigmai[indx[n]],2);
  if (nusewt == -2 || nusewt > 0) *wt = 1.0/(pow(sigmai[indx[n]],2)+pow(addvar,2));
  if (nusewt == -3) *wt = 1.0/addvar/addvar;
}

void readin(pulsar psr)
{
  int nalv,nread,n;
  static int idiag = 0;
  double tlast,ave,rms,var,sumwt,alv,t,x,sig,wt,t1,xr,xi,tautest;
  double xold,rootvar,permin;
  static double unitfact = 1.0;
  double mean=0.0;
  int i;
  int loop=0;

  for (i=0;i<psr.nobs;i++)
    mean+=psr.obsn[i].residual;
  mean/=(double)psr.nobs;

  npt = 0;
  utjdlast=-1.0; 
  ave=0.0;
  rms=0.0;
  var=0.0;
  sumwt=0.0;
  alv=0.0;
  nalv=0;
  tlast=0.0;
  sig=1.0;
  wt=1.0;

  // Should check ntunits and nxunits - assume both are zero
  ntunits=0;

  tauensure = -1.0;
  tauensure = tauensure*86400.0;
  unitfact = 1.0e-6;

  // 10 continue statement
  int ic=0;
  int breakit=0;

  do {
    do {
      do {
	// Assuming that nformat = 0
	
	if (nusewt == 0 || nusewt == -3)
	  {
	    t = (double)psr.obsn[ic].sat;
	    x = (double)(psr.obsn[ic].residual-mean)/1e-6;
	    if (ic>0) {
	      if (psr.obsn[ic-1].sat > psr.obsn[ic].sat)
		{
		  printf("1: ERROR: The TOAs need to be time sorted, sorry! - sats %g and %g\n",(double)psr.obsn[ic-1].sat,(double)psr.obsn[ic].sat);
		  //		  exit(1);
		}
	    }
	    ic++;
	    if (ic == psr.nobs) breakit=1;
	  }
	else
	  {
	    t = (double)psr.obsn[ic].sat;
	    x = (double)(psr.obsn[ic].residual-mean)/1e-6;
	    sig = (double)(psr.obsn[ic].toaErr);
	    if (ic>0) {
	      if (psr.obsn[ic-1].sat > psr.obsn[ic].sat)
		{
		  printf("2: ERROR: The TOAs need to be time sorted, sorry! sats: %.15g and %.15g\n",(double)psr.obsn[ic-1].sat,(double)psr.obsn[ic].sat);
		  //		  exit(1);
		}
	    }
	    ic++;
	    if (ic == psr.nobs) breakit=1;
	  }
	nread++;
      } while ((t < tmin || t > tmax) && breakit==0);
      if (npt == 0) t1=t;
      if (t < utfirst) utfirst=t;
      if (t > utlast) utlast=t;
      t=t-t1; // to avoid losing precision
      if (ntunits==0) t=t*86400.0;
    } while ((x < xmin || x > xmax) && breakit==0);
    x = unitfact*x;
    if (nusewt != 0 && nusewt != -3) sig = unitfact*sig;
    npt++;
    indx[npt] = npt;
    utjd[indx[npt]]=t;
    data[indx[npt]]=x;
    sigmai[indx[npt]]=sig;
    utjdlast = t;
    wt = 1.0/sig/sig;
    ave=ave+wt*x;
    var=var+wt*x*x;
    sumwt=sumwt+wt;
    if (npt > 1)
      {
	nalv++;
	alv = alv+pow(x-xold,2);
      }
  } while (breakit==0);
  // Ave different because I have a better mean removal or something similar
  if (npt <= 0) exit(1);
  indexx8(npt,utjd,indx);
  tdiffmin = utjd[indx[2]]-utjd[indx[1]];
  for (n=1;n<=npt-1;n++)
    {
      if (tdiffmin > utjd[indx[n+1]]-utjd[indx[n]])
	tdiffmin = utjd[indx[n+1]]-utjd[indx[n]];
    }
  ave=ave/sumwt;
  var=var/sumwt;
  rootvar=sqrt(var);
  rms=var-ave*ave;
  rms = sqrt(rms);
  if (nalv > 0 && alv > 0.0) alv=sqrt(alv/2.0/nalv);
  permax = utjd[indx[npt]]-utjd[indx[1]];
  permin = permax/(npt-1.0);
  taumin = permin;
  taumax = permax;
  if (tauensure > 0.0) {
    taumax = tauensure;
    //105 CONTINUE
    loop=0;
    do
      {
	if (taumax < permax)
	  {
	    loop=1;
	    taumax=2.0*taumax;
	  }
      }while (loop==1);
  }
  // Should check nbintype
  tautest = taumax;
  do {
    tautest=tautest/2.0;
  } while (tautest > taumin);
  taumin=tautest;

  if (nusewt >= 0.0) addvar = rms;
}


// Check conversion from Fortran carefully
void indexx8(int n,double *arrin,int *indx)
{
  //     from numerical recipes
  double q;
  int j,l,ir,i=0,indxt;

  for (j=1;j<=n;j++)
    indx[j] = j;
  l = n/2+1;
  ir = n;

 pos10:
  if (l > 1){
    l--;
    indxt = indx[l];
    q=arrin[indxt];
  } else {
    indxt=indx[ir];
    q = arrin[indxt];
    indx[ir]=indx[1];
    ir--;
    if (ir == 1){
      indx[1] = indxt;
      return;
    }
  }
  i=l;
  j=l+l;
 pos20:
 if (j <= ir) {
    if (j < ir) {
      if (arrin[indx[j]] < arrin[indx[j+1]])j++;
    }
    if (q < arrin[indx[j]]){
      indx[i] = indx[j];
      i=j;
      j=j+j;
    } else {
      j=ir+1;
    }
    goto pos20; // Sorry about this - copying the fortran!
 }
 indx[i]=indxt;
 goto pos10;
}
    




void sortTimes(pulsar psr,int *nobs,double *times,double *resid,double *error)
{
  int i,changed,tot=0;
  double t1,t2,t3;
  double tt2[MAX_OBSN],r2[MAX_OBSN],e2[MAX_OBSN];

  for (i=0;i<*nobs;i++)
    {
      times[i] = (double)(psr.obsn[i].sat - psr.param[param_pepoch].val[0]);
      resid[i] = (double)psr.obsn[i].residual;
      error[i] = (double)psr.obsn[i].toaErr;     
    }
  /* Now sort */
  do {
    changed=0;

    for (i=0;i<*nobs-1;i++)
      {
	if (times[i] > times[i+1])
	  {
	    t1 = times[i]; t2 = resid[i]; t3 = error[i];
	    times[i] = times[i+1];
	    resid[i] = resid[i+1];
	    error[i] = error[i+1];
	    times[i+1] = t1;
	    resid[i+1] = t2;
	    error[i+1] = t3;
	    changed=1;
	  }
      }
  } while (changed==1);
  /* Remove simultaneous observations */
  /* for (i=0;i<*nobs;i++)
    {
      tt2[tot] = times[i];
      r2[tot] = resid[i];
      e2[tot] = error[i];
      printf("REMOVE %g\n",tt2[tot]);
      if (tot==0 || (tt2[tot]-tt2[tot-1])>1.0e-5)
	tot++;
      else
	printf("REMOVE\n");
    }
  for (i=0;i<tot;i++)
    {
      times[i] = tt2[i];
      resid[i] = r2[i];
      error[i] = e2[i];
      } 
      *nobs = tot; */
}

void fitv(double x,double afunc[],int ma,pulsar *psr,int ipos)
{
  int i;
  for (i=0;i<ma;i++)
    afunc[i+1] = pow(x,i);
}




/*
   This returns the point "output" on the spline curve.
   The parameter "v" indicates the position, it ranges from 0 to n-t+2
   
*/
void SplinePoint(int *u,int n,int t,float v,XY *control,XY *output)
{
   int k;
   float b;

   output->x = 0;
   output->y = 0;

   for (k=0;k<=n;k++) {
      b = SplineBlend(k,t,u,v);
      output->x += control[k].x * b;
      output->y += control[k].y * b;
   }
}

/*
   Calculate the blending value, this is done recursively.
   
   If the numerator and denominator are 0 the expression is 0.
   If the deonimator is 0 the expression is 0
*/
float SplineBlend(int k,int t,int *u,float v)
{
   float value;

   if (t == 1) {
      if ((u[k] <= v) && (v < u[k+1]))
         value = 1;
      else
         value = 0;
   } else {
      if ((u[k+t-1] == u[k]) && (u[k+t] == u[k+1]))
         value = 0;
      else if (u[k+t-1] == u[k]) 
         value = (u[k+t] - v) / (u[k+t] - u[k+1]) * SplineBlend(k+1,t-1,u,v);
      else if (u[k+t] == u[k+1])
         value = (v - u[k]) / (u[k+t-1] - u[k]) * SplineBlend(k,t-1,u,v);
     else
         value = (v - u[k]) / (u[k+t-1] - u[k]) * SplineBlend(k,t-1,u,v) + 
                 (u[k+t] - v) / (u[k+t] - u[k+1]) * SplineBlend(k+1,t-1,u,v);
   }
   return(value);
}

/*
   The positions of the subintervals of v and breakpoints, the position
   on the curve are called knots. Breakpoints can be uniformly defined
   by setting u[j] = j, a more useful series of breakpoints are defined
   by the function below. This set of breakpoints localises changes to
   the vicinity of the control point being modified.
*/
void SplineKnots(int *u,int n,int t)
{
   int j;

   for (j=0;j<=n+t;j++) {
     if (j < t)
       u[j] = 0;
     else if (j <= n)
       u[j] = j - t + 1;
     else if (j > n)
       u[j] = n - t + 2;	
   }
}

/*-------------------------------------------------------------------------
   Create all the points along a spline curve
   Control points "inp", "n" of them.
   Knots "knots", degree "t".
   Ouput curve "outp", "res" of them.
*/
void SplineCurve(XY *inp,int n,int *knots,int t,XY *outp,int res)
{
   int i;
   float interval,increment;

   interval = 0;
   increment = (n - t + 2) / (float)(res - 1);
   for (i=0;i<res-1;i++) {
      SplinePoint(knots,n,t,interval,inp,&(outp[i]));
      interval += increment;
   }
   //   printf("Setting %d %d %g\n",res-1,n,inp[n]);
   outp[res-1] = inp[n];
}

/*
   Example of how to call the spline functions
	Basically one needs to create the control points, then compute
   the knot positions, then calculate points along the curve.
*/

void calcSpline(float *px,float *py,int count)
{
   int i;
   int T = 4;
   int RESOLUTION=20;
   float fx[RESOLUTION],fy[RESOLUTION];
   XY inp[count+1];
   int knots[count+T+1];
   XY outp[RESOLUTION+1];

   for (i=0;i<count;i++)
     {
       inp[i].x = px[i];
       inp[i].y = py[i];	 
     }

   SplineKnots(knots,count-1,T);
   SplineCurve(inp,count-1,knots,T,outp,RESOLUTION);
   for (i=0;i<RESOLUTION;i++)
     {
       fx[i] = outp[i].x;
       fy[i] = outp[i].y;
       //       printf("OUTPUT = %g %g\n",fx[i],fy[i]);
     }
   // PUT BACK IN
   //   cpgline(RESOLUTION,fx,fy);
   /* Display the curve, in this case in OOGL format for GeomView */
   //   for (i=0;i<RESOLUTION;i++)
   //      printf("out %g %g\n",outp[i].x,outp[i].y);


}
char * plugVersionCheck = TEMPO2_h_VER;
