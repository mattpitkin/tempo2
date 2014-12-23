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

/* Template for a tempo2 plugin */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "T2toolkit.h"
#include "TKfit.h"
#include <cpgplot.h>
#include <float.h>

using namespace std;

typedef struct glitchS {
  double glep;
  double glph;
  double glf0;
  double glf1;
  double glf0d;
  double gltd; 
  int fitph;
  int fitf0;
  int fitf1;
  int fitf0d;
  int fittd;
} glitchS;

void defineGlitchVal(glitchS *glitch,int nglt);
void doPlot(double *epoch,double *f0,double *f0e,double *f1,double *f1e,int fitf1,int *nFit,int *id,int n,float *gt,int ngt,int *plotType,int nplot,double plotOffset,double *plotResX,double *plotResY,double *plotResE,int nplotVal,int combine,float fontSize,char *title,float *yscale_min,float *yscale_max,int *yscale_set,int interactive);
void plot1(double *epoch,double *f0,float *yerr1,float *yerr2,
	   int n,double plotOffset,int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx);          
void plot8(double *epoch,double *f0,float *yerr1,float *yerr2,
	   int n,double plotOffset,int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx);          
void plot2(double *epoch,double *f0,double *f0e,int *nFit,int *id,int n,double plotOffset,int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx);  
void interactivePlot(double *epoch,double *f0,double *f0e,int *nFit,int *id,int n);
void fitFuncs(double x,double *p,int m);
void changeFit(glitchS *glitch,int nglt);
void drawMenu(float minx,float maxx,float miny,float maxy,glitchS *glitch, int nglt,
	      int fitf0,int fitf1);
void checkMenu(float minx,float maxx,float miny,float maxy,glitchS *glitch, int nglt,
	       float mx,float my,int *fitf0,int *fitf1,char key);
double nonlinearFunc( double t, const double *par,int obsNum );
void plot3(double *epoch,double *f1,double *f1e,int *nFit,int *id,int n,double plotOffset,int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx);
void plot6(double *epoch,double *f1,double *f1e,int *nFit,int *id,int n,double plotOffset,int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx);
void plot4(double *epoch,double *f0,double *f0e,int *nFit,int *id,int n,double plotOffset,float *gt,int ngt,int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx);
void plot5(double *epoch,double *f0,double *f0e,int *nFit,int *id,int n,double plotOffset,float *gt,int ngt,int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx);
void plot7(double *plotResX,double *plotResY,double *plotResE,int nplotVal,double plotOffset,int combine,int pos,int nplot,double start,double end,double psrF0,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx);
void plot9(double *epoch,double *f0,double *f0e,int *nFit,int *id,int n,double plotOffset,
	   int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx);


// Global parameters for the fit
glitchS *global_glitch;
int global_nglt;
int global_fitf0;
int global_fitf1;
double global_valf0;
double global_valf1;
double global_footer;
double global_header;


// Functions for nonlinear fit
// Various setup routines for the LM-fit procedure
/* Collection of control (input) parameters. */
typedef struct {
    double ftol;      /* relative error desired in the sum of squares. */
    double xtol;      /* relative error between last two approximations. */
    double gtol;      /* orthogonality desired between fvec and its derivs. */
    double epsilon;   /* step used to calculate the jacobian. */
    double stepbound; /* initial bound to steps in the outer loop. */
    int maxcall;      /* maximum number of iterations. */
    int scale_diag;   /* UNDOCUMENTED, TESTWISE automatical diag rescaling? */
    int printflags;   /* OR'ed to produce more noise */
} lm_control_struct;

/* Collection of status (output) parameters. */
typedef struct {
    double fnorm;     /* norm of the residue vector fvec. */
    int nfev;	      /* actual number of iterations. */
    int info;	      /* status of minimization. */
} lm_status_struct;

/* Recommended control parameter settings. */
extern const lm_control_struct lm_control_double;
extern const lm_control_struct lm_control_float;



/* Standard monitoring routine. */
void lm_printout_std( int n_par, const double *par, int m_dat,
                      const void *data, const double *fvec,
                      int printflags, int iflag, int iter, int nfev);

/* Refined calculation of Eucledian norm, typically used in printout routine. */
double lm_enorm( int, const double * );

/* The actual minimization. */
void lmmin( int n_par, double *par, int m_dat, const void *data, 
            void (*evaluate) (const double *par, int m_dat, const void *data,
                              double *fvec, int *info),
            const lm_control_struct *control, lm_status_struct *status,
            void (*printout) (int n_par, const double *par, int m_dat,
                              const void *data, const double *fvec,
                              int printflags, int iflag, int iter, int nfev) );


/** Legacy low-level interface. **/

/* Alternative to lm_minimize, allowing full control, and read-out
   of auxiliary arrays. For usage, see implementation of lm_minimize. */
void lm_lmdif( int m, int n, double *x, double *fvec, double ftol,
	       double xtol, double gtol, int maxfev, double epsfcn,
	       double *diag, int mode, double factor, int *info, int *nfev,
	       double *fjac, int *ipvt, double *qtf, double *wa1,
	       double *wa2, double *wa3, double *wa4,
               void (*evaluate) (const double *par, int m_dat, const void *data,
                                 double *fvec, int *info),
               void (*printout) (int n_par, const double *par, int m_dat,
                                 const void *data, const double *fvec,
                                 int printflags, int iflag, int iter, int nfev),
               int printflags, const void *data );

extern const char *lm_infmsg[];
extern const char *lm_shortmsg[];

typedef struct {
    const double *t;
    const double *y;
    double (*f) (double t, const double *par,int obsNum);
} lmcurve_data_struct;

void lmcurve_fit( int n_par, double *par, int m_dat,
                  const double *t, const double *y,
                  double (*f)( double t, const double *par, int obsNum ),
                  const lm_control_struct *control, lm_status_struct *status );


#define LM_MACHEP     DBL_EPSILON   /* resolution of arithmetic */
#define LM_DWARF      DBL_MIN       /* smallest nonzero number */
#define LM_SQRT_DWARF sqrt(DBL_MIN) /* square should not underflow */
#define LM_SQRT_GIANT sqrt(DBL_MAX) /* square should not overflow */
//#define LM_USERTOL    30*LM_MACHEP  /* users are recommended to require this */
#define LM_USERTOL 1e-9 // IMPORTANT TO SET TO SOMETHING SENSIBLE <<<<<

const lm_control_struct lm_control_double = {
    LM_USERTOL, LM_USERTOL, LM_USERTOL, LM_USERTOL, 100., 100, 1, 0 };
const lm_control_struct lm_control_float = {
    1.e-7, 1.e-7, 1.e-7, 1.e-7, 100., 100, 0, 0 };

void lm_lmpar( int n, double *r, int ldr, int *ipvt, double *diag,
	       double *qtb, double delta, double *par, double *x,
	       double *sdiag, double *aux, double *xdi );
void lm_qrfac( int m, int n, double *a, int pivot, int *ipvt,
	       double *rdiag, double *acnorm, double *wa );
void lm_qrsolv( int n, double *r, int ldr, int *ipvt, double *diag,
	        double *qtb, double *x, double *sdiag, double *wa );


#define MIN(a,b) (((a)<=(b)) ? (a) : (b))
#define MAX(a,b) (((a)>=(b)) ? (a) : (b))
#define SQR(x)   (x)*(x)
// END SETUP FOR THE LM ROUTINE





void help() /* Display help */
{
  /* This function should contain usage information about the plugin which should (in general) be accessed */
  /* by the user pressing 'h'                                                                              */
  printf("-combine  combine multiple plots\n");
  printf("-fitf1    fit for F0 and F1 (default just F0)\n");  
  printf("-font x   set font size\n");
  printf("-foot x   set fraction of page as footer\n");
  printf("-gt x     set glitch at MJD x (can use -gr multiple times)\n");
  printf("-h        this help file\n");
  printf("-head x   set fraction of page as header\n");
  printf("-loadResults x Load a result.dat file from 'x'\n");
  printf("-offset x set time offset\n");
  printf("-p x      make plot of type x (can be used multiple times)\n");
  printf("-removeF2 remove the pre-fit F2 value from the resulting plots\n");
  printf("-t        data file containing MJD ranges and .par files\n");
  printf("-title x  set title\n");

  printf("\n\n");
  printf("plot: 1 - F0 versus day\n");
  printf("      2 - F0 with gradient removed and data information provided\n");
  printf("      3 - F1\n");
  printf("      4 - F0 with the pre-glitch gradient removed\n");
  printf("      5 - F0 with the mean post-glitch values removed\n");
  printf("      6 - F1 with mean removed\n");
  printf("      7 - post-fit timing residuals\n");
  printf("      8 - F0 with mean removed\n");
  printf("      9 - F0 with gradient removed but no data information\n");
  printf("\n\n");
  printf("Example usage: tempo2 -gr glitch -f mypar.par mytim.tim -t jianping.dat -fitf1 -offset 52600 -gt 53444 -combine -font 1.2 -foot 0.1 -head 0.1 -title \"PSR B1800-21\" -p 6 -p 4\n");
  exit(1);
}

#define MAX_TIMES 2000

/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  char origArgV[100][MAX_FILELEN];
  int i,j,k;
  double globalParameter;
  char timeList[MAX_STRLEN];
  long double *mjd1,*mjd2;
  long double centreMJD;
  char parFileName[MAX_TIMES][MAX_STRLEN];
  char timFileName[MAX_TIMES][MAX_STRLEN];
  double *epoch;
  double *f0;
  double *f0e;
  double *f1;
  double *f1e;
  double plotOffset=0;
  char plotPar[100]="NULL";
  int fitf1=0;
  int nFit[MAX_TIMES];
  int id[MAX_TIMES];
  int nStride=0;
  int n=0;
  int argn=0;
  float gt[MAX_TIMES];
  double plotResX[MAX_OBSN],plotResY[MAX_OBSN],plotResE[MAX_OBSN];
  int  nplotVal;
  int plotType[100];
  int nplot=0;
  int ngt=0;
  int okay=0;
  int interactive=0;
  int combine=0;
  float fontSize=1.4;
  char title[100]="";
  char tname[100]="";
  int nread;
  char resultFname[1024];
  int loadResults=0;
  float yscale_min[100];
  float yscale_max[100];
  int   yscale_set[100];
  int removeF2=0;
  long double storeEpoch;
  long double storeExpectedF0,storeExpectedF1,storeExpectedF2;
  long double storeOrigF0,storeOrigF1,storeOrigF2;
  int removeExpected=0;

  // Allocate memory
  printf("Allocating memory\n");
  if (!(mjd1 = (long double *)malloc(sizeof(long double)*MAX_TIMES)))
    {printf("Sorry: out of memory\n"); exit(1);}
  if (!(mjd2 = (long double *)malloc(sizeof(long double)*MAX_TIMES)))
    {printf("Sorry: out of memory\n"); exit(1);}
  if (!(epoch = (double *)malloc(sizeof(double)*MAX_TIMES)))
    {printf("Sorry: out of memory\n"); exit(1);}
  if (!(f0 = (double *)malloc(sizeof(double)*MAX_TIMES)))
    {printf("Sorry: out of memory\n"); exit(1);}
  if (!(f0e = (double *)malloc(sizeof(double)*MAX_TIMES)))
    {printf("Sorry: out of memory\n"); exit(1);}
  if (!(f1 = (double *)malloc(sizeof(double)*MAX_TIMES)))
    {printf("Sorry: out of memory\n"); exit(1);}
  if (!(f1e = (double *)malloc(sizeof(double)*MAX_TIMES)))
    {printf("Sorry: out of memory\n"); exit(1);}


  for (i=0;i<100;i++)
    yscale_set[i]=0;

  global_footer = 0.15;
  global_header = 0.1;

  FILE *fin;
  FILE *fout;

  for (i=0;i<argc;i++)
    strcpy(origArgV[i],argv[i]);

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: glitch\n");
  printf("Author:              George Hobbs\n");
  printf("CVS Version:         $Revision $\n");
  printf(" --- type 'h' for help information\n");


  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  argn=i;
	  strcpy(parFile[0],argv[++i]); 
	  strcpy(timFile[0],argv[++i]);
	}
      else if (strcmp(argv[i],"-h")==0)
	help();
      else if (strcmp(argv[i],"-title")==0)
	strcpy(title,argv[++i]);
      else if (strcmp(argv[i],"-t")==0)
	strcpy(timeList,argv[++i]);
      else if (strcmp(argv[i],"-gt")==0)
	sscanf(argv[++i],"%f",&gt[ngt++]);
      else if (strcmp(argv[i],"-p")==0)
	sscanf(argv[++i],"%d",&plotType[nplot++]);
      else if (strcmp(argv[i],"-offset")==0)
	sscanf(argv[++i],"%lf",&plotOffset);
      else if (strcmp(argv[i],"-i")==0)
	interactive=1;
      else if (strcmp(argv[i],"-fitf1")==0)
	fitf1=1;
      else if (strcmp(argv[i],"-font")==0)
	sscanf(argv[++i],"%f",&fontSize);
      else if (strcmp(argv[i],"-head")==0)
	sscanf(argv[++i],"%lf",&global_header);
      else if (strcmp(argv[i],"-foot")==0)
	sscanf(argv[++i],"%lf",&global_footer);
      else if (strcmp(argv[i],"-combine")==0)
	combine=1;
      else if (strcmp(argv[i],"-plotPar")==0)
	strcpy(plotPar,argv[++i]);
      else if (strcmp(argv[i],"-removeF2")==0)
	removeF2=1;
      else if (strcmp(argv[i],"-removeExpected")==0)
	removeExpected=1;
      else if (strcmp(argv[i],"-loadResults")==0)
	{
	  loadResults=1;
	  strcpy(resultFname,argv[++i]);
	}
      else if (strcmp(argv[i],"-yscale")==0)
	{
	  int n;
	  sscanf(argv[++i],"%d",&n);
	  sscanf(argv[++i],"%f",&yscale_min[n]);
	  sscanf(argv[++i],"%f",&yscale_max[n]);
	  yscale_set[n] = 1;
	}
    }
  if (loadResults==0)
    {
      readParfile(psr,parFile,timFile,1); /* Load the parameters       */
      readTimfile(psr,timFile,1); /* Load the arrival times    */
      storeEpoch = psr[0].param[param_pepoch].val[0];
      storeOrigF0 = psr[0].param[param_f].val[0];
      storeOrigF1 = psr[0].param[param_f].val[1];
      storeOrigF2 = psr[0].param[param_f].val[2];
    
      // Read the list of strides
      fin = fopen(timeList,"r");
      while (!feof(fin))
	{
	  nread = fscanf(fin,"%Lf %Lf %s %s",&mjd1[nStride],&mjd2[nStride],parFileName[nStride],tname);
	  if (nread==3 || nread==4)
	    {
	      if (nread==4)
		strcpy(timFileName[nStride],tname);
	      else
		strcpy(timFileName[nStride],timFile[0]);
	      nStride++;
	      if (nStride > MAX_TIMES)
		{
		  printf("Sorry the maximum number of regions that can be loaded is %d\n",MAX_TIMES);
		  exit(1);
		}
	    }
	}
      fclose(fin);
      printf("Have read %d strides\n",nStride);
      if (nStride == 0)
	{
	  printf("ERROR: require at least one stride - use -t option\n");
	  exit(1);
	}
    }

  // Now do the stride fitting
  if (loadResults==0)
    {
      if (!(fout = fopen("result.dat","w")))
	{
	  printf("Unable to open output file result.dat\n");
	  exit(1);
	}
      for (i=0;i<nStride;i++)
	{
	  centreMJD = (mjd1[i]+mjd2[i])/2.0;
	  printf("Analysing %g %g %g\n",(double)mjd1[i],(double)mjd2[i],(double)centreMJD);
	  
	  strcpy(parFile[0],parFileName[i]);
	  strcpy(timFile[0],timFileName[i]);
	  
	  psr[0].nJumps=0;
	  psr[0].nobs=0;
	  for(j=0;j<MAX_PARAMS;j++){
	    psr[0].param[j].nLinkTo = 0;
	    psr[0].param[j].nLinkFrom = 0;
	    for (k=0;k<psr[0].param[j].aSize;k++)
	      {
		psr[0].param[j].paramSet[k] = 0;
		psr[0].param[j].fitFlag[k] = 0;
		psr[0].param[j].prefit[k] = 0.0;
		psr[0].param[j].val[k] = 0.0;
	      }
	  }
	  readParfile(psr,parFile,timFile,1); /* Load the parameters       */
	  readTimfile(psr,timFile,1); /* Load the arrival times    */
	  printf("ntoas = %d\n",psr[0].nobs);
	  // Update the epoch in the par file for centreMJD
	  strcpy(argv[argn],"-epoch");
	  sprintf(argv[argn+1],"%.5f",(double)centreMJD);
	  preProcess(psr,1,argc,argv);      
	  
	  storeExpectedF0 = storeOrigF0 + storeOrigF1*(centreMJD-storeEpoch)*86400.0L + 0.5*storeOrigF2*pow((centreMJD-storeEpoch)*86400.0L,2);
	  storeExpectedF1 = storeOrigF1 + storeOrigF2*(centreMJD-storeEpoch)*86400.0L;
	  
	  // Turn off all fitting
	  for (j=0;j<MAX_PARAMS;j++)
	    {
	      for (k=0;k<psr[0].param[j].aSize;k++)
		psr[0].param[j].fitFlag[k] = 0;
	    }
	  
	  // Update the start and finish flags
	  psr[0].param[param_start].val[0] = mjd1[i];
	  psr[0].param[param_start].fitFlag[0] = 1;
	  psr[0].param[param_start].paramSet[0] = 1;
	  
	  psr[0].param[param_finish].val[0] = mjd2[i];
	  psr[0].param[param_finish].fitFlag[0] = 1;
	  psr[0].param[param_finish].paramSet[0] = 1;
	  
	  // Turn on required fitting
	  psr[0].param[param_f].fitFlag[0] = 1;
	  if (fitf1==1)
	    psr[0].param[param_f].fitFlag[1] = 1;
	  
	  // Do the fitting
	  for (j=0;j<2;j++)                   /* Do two iterations for pre- and post-fit residuals*/
	    {
	      formBatsAll(psr,1);         /* Form the barycentric arrival times */
	      formResiduals(psr,1,1);    /* Form the residuals                 */
	      if (j==0) doFit(psr,1,0);   /* Do the fitting     */
	      else textOutput(psr,1,globalParameter,0,0,0,(char *)"");  /* Display the output */
	    }
	  if ((psr[0].nFit>1 && fitf1==0) || (fitf1==1 && psr[0].nFit>2))
	    {
	      epoch[n] = (double)centreMJD;
	      f0[n]    = (double)psr[0].param[param_f].val[0];
	      if ((psr[0].nFit==2 && fitf1==0))
		f0e[n]=-1; // Set a negative error
	      else
		f0e[n]   = (double)psr[0].param[param_f].err[0];
	      if (removeF2==1)
		f0[n] -= (double)(0.5*psr[0].param[param_f].val[2]*pow((centreMJD-storeEpoch)*86400.0,2));
	      if (removeExpected==1)
		f0[n] -= storeExpectedF0;
	      
	      if (fitf1==1)
		{
		  f1[n] = (double)psr[0].param[param_f].val[1];
		  if (psr[0].nFit==3)
		    {
		      f1e[n] = -1;
		      f0e[n] = -1;
		    }
		  else
		    f1e[n] = (double)psr[0].param[param_f].err[1];
		  if (removeF2==1)
		    {
		      //		  printf("Updating by: %g\n",(double)(psr[0].param[param_f].val[2]*(centreMJD-psr[0].param[param_pepoch].val[0])*86400.0));
		      f1[n] -= (double)(psr[0].param[param_f].val[2]*(centreMJD-storeEpoch)*86400.0);
		    }
		  if (removeExpected==1)
		    f1[n] -= storeExpectedF1;
		}
	      nFit[n]  = psr[0].nFit;
	      id[n]  = i;
	      if (fitf1==0)
		fprintf(fout,"%g %.10f %g %d %d\n",epoch[n],f0[n],f0e[n],nFit[n],id[n]+1);
	      else
		fprintf(fout,"%g %.10f %g %g %g %d %d\n",epoch[n],f0[n],f0e[n],f1[n],f1e[n],nFit[n],id[n]+1);
	      n++;    
	    }
	  else 
	    fprintf(fout,"# %g %d failed: nfit = %d \n",(double)centreMJD,i+1,psr[0].nFit);
	}
      fclose(fout);
    }
  else
    {
      FILE *fin;
      printf("Loading results from %s\n",resultFname);
      fin = fopen(resultFname,"r");
      n=0;
      while (!feof(fin))
	{
	  if (fitf1==0)
	    {
	      if (fscanf(fin,"%lf %lf %lf %d %d",&epoch[n],&f0[n],&f0e[n],&nFit[n],&id)==5)
		{
		  id[n]--;
		  n++;
		}
	    }
	  else
	    {
	      char str[1024];
	      printf("At this point %d\n",n);
	      fgets(str,1024,fin);
	      if (str[0]=='#')
		printf("Warning line: %s\n",str);
	      else if (sscanf(str,"%lf %lf %lf %lf %lf %d %d",&epoch[n],&f0[n],&f0e[n],&f1[n],&f1e[n],&nFit[n],&id[n])==7)
		{
		  id[n]--;
		  n++;
		}
	      printf("Read %g %d\n",epoch[n-1],id[n-1]);
	    }
	}
      fclose(fin);
      printf("Complete\n");
    }
  // Get timing residuals for plotting
  int needRes=0;
  for (i=0;i<nplot;i++)
    {
      if (plotType[i] == 7)
	needRes=1;
    }
  if (needRes==1)
    {
      strcpy(argv[argn],origArgV[argn]);
      strcpy(argv[argn+1],origArgV[argn+1]);
      strcpy(argv[argn+2],origArgV[argn+2]);
      strcpy(timFile[0],origArgV[argn+2]);
      strcpy(parFile[0],origArgV[argn+1]);
      if (strcmp(plotPar,"NULL")==0)
	strcpy(parFile[0],parFileName[0]);
      else
	strcpy(parFile[0],plotPar);
      
      psr[0].nJumps=0;
      for(j=0;j<MAX_PARAMS;j++){
	psr[0].param[j].nLinkTo = 0;
	psr[0].param[j].nLinkFrom = 0;
      }
      printf("\n\n------------------------------------------\n");
      printf("Reading the par file after the stridefits\n");
      for (i=0;i<argc;i++)
	strcpy(argv[i],origArgV[i]);
      
      for (i=0;i<argc;i++)
	printf("%s\n",argv[i]);
      readParfile(psr,parFile,timFile,1); /* Load the parameters       */
      readTimfile(psr,timFile,1); /* Load the arrival times    */
      //  psr[0].param[param_start].fitFlag[0] = 0;
      //  psr[0].param[param_start].paramSet[0] = 0;
      //  psr[0].param[param_finish].fitFlag[0] = 0;
      //  psr[0].param[param_finish].paramSet[0] = 0;
      
      preProcess(psr,1,argc,argv);        
      // Do the fitting
      for (j=0;j<2;j++)                   /* Do two iterations for pre- and post-fit residuals*/
	{
	  formBatsAll(psr,1);         /* Form the barycentric arrival times */
	  formResiduals(psr,1,1);    /* Form the residuals                 */
	  if (j==0) doFit(psr,1,0);   /* Do the fitting     */
	  else textOutput(psr,1,globalParameter,0,0,0,(char *)"");  /* Display the output */
	}
    
      nplotVal=0;
      for (i=0;i<psr[0].nobs;i++)
	{
	  okay=1;
	  if (psr[0].obsn[i].deleted!=0)
	    okay=0;
	  if (psr[0].param[param_start].paramSet[0]==1 && psr[0].param[param_start].fitFlag[0]==1 &&
	      (psr[0].param[param_start].val[0] > psr[0].obsn[i].sat))
	    okay=0;
	  if (psr[0].param[param_finish].paramSet[0]==1 && psr[0].param[param_finish].fitFlag[0]==1 &&
	      psr[0].param[param_finish].val[0] < psr[0].obsn[i].sat)
	    okay=0;
	  if (okay==1)
	    {
	      plotResX[nplotVal] = (double)(psr[0].obsn[i].sat);
	      plotResY[nplotVal] = (double)(psr[0].obsn[i].residual);
	      plotResE[nplotVal] = (double)(psr[0].obsn[i].toaErr*1.0e-6);
	      nplotVal++;
	    }
	}
    }
      // Do the plotting
  printf("Doing the plot\n");
  //  if (interactive==0)
    doPlot(epoch,f0,f0e,f1,f1e,fitf1,nFit,id,n,gt,ngt,plotType,nplot,plotOffset,plotResX,plotResY,plotResE,nplotVal,combine,fontSize,title,yscale_min,yscale_max,yscale_set,interactive);
  //  else
  //    interactivePlot(epoch,f0,f0e,nFit,id,n);
  printf("Done the plot\n");

  // De-allocate memory
  free(mjd1); free(mjd2); free(epoch); free(f0); free(f0e); free(f1); free(f1e);


  return 0;
}

void doPlot(double *epoch,double *f0,double *f0e,double *f1,double *f1e,int fitf1,int *nFit,int *id,int n,float *gt,int ngt,int *plotType,int nplot,double plotOffset,double *plotResX,double *plotResY,double *plotResE,int nplotVal,int combine,float fontSize,char *title,float *yscale_min,float *yscale_max,int *yscale_set,int interactive)
{
  float yerr1[n],yerr2[n];
  float fx[2],fy[2];
  float mx,my;
  char key;
  int i,k;
  int plot=2;
  int last;
  double minx=-1,maxx=-1;

  minx = epoch[0];
  maxx = epoch[n-1];


  for (i=0;i<n;i++)
    {
      yerr1[i] = (float)(f0[i]-f0e[i]);
      yerr2[i] = (float)(f0[i]+f0e[i]);
    }

  cpgbeg(0,"?",1,1);
  
  /*  if (combine==1)
    {
      cpgsubp(1,nplot);
      }*/
  cpgsch(fontSize);
  printf("Here setting font\n");
  cpgscf(2);
  cpgslw(2);

  do {
    for (k=0;k<nplot;k++)
      {
	last=0;
	if (plotType[k]==1)
	  plot1(epoch,f0,yerr1,yerr2,n,plotOffset,combine,k,nplot,yscale_min[1],yscale_max[1],yscale_set[1],minx,maxx);           
	else if (plotType[k]==2)
	  plot2(epoch,f0,f0e,nFit,id,n,plotOffset,combine,k,nplot,yscale_min[2],yscale_max[2],yscale_set[2],minx,maxx);           
	else if (plotType[k]==3)
	  plot3(epoch,f1,f1e,nFit,id,n,plotOffset,combine,k,nplot,yscale_min[3],yscale_max[3],yscale_set[3],minx,maxx);
	else if (plotType[k]==4)
	  plot4(epoch,f0,f0e,nFit,id,n,plotOffset,gt,ngt,combine,k,nplot,yscale_min[4],yscale_max[4],yscale_set[4],minx,maxx);           
	else if (plotType[k]==5)
	  plot5(epoch,f0,f0e,nFit,id,n,plotOffset,gt,ngt,combine,k,nplot,yscale_min[5],yscale_max[5],yscale_set[5],minx,maxx);           
	else if (plotType[k]==6)
	  plot6(epoch,f1,f1e,nFit,id,n,plotOffset,combine,k,nplot,yscale_min[6],yscale_max[6],yscale_set[6],minx,maxx);
	else if (plotType[k]==7)
	  plot7(plotResX,plotResY,plotResE,nplotVal,plotOffset,combine,k,nplot,epoch[0],epoch[n-1],f0[0],yscale_min[7],yscale_max[7],yscale_set[7],minx,maxx);
	else if (plotType[k]==8)
	  plot8(epoch,f0,yerr1,yerr2,n,plotOffset,combine,k,nplot,yscale_min[8],yscale_max[8],yscale_set[8],minx,maxx);           
	else if (plotType[k]==9)
	  plot9(epoch,f0,f0e,nFit,id,n,plotOffset,combine,k,nplot,yscale_min[9],yscale_max[9],yscale_set[9],minx,maxx);           
	
	for (i=0;i<ngt;i++)
	  {
	    fx[0] = fx[1] = gt[i]-plotOffset;
	    fy[0] = -1000; fy[1] = 1000;
	    cpgsls(4); cpgline(2,fx,fy); cpgsls(1);
	    
	  }
      }
    if (interactive==1)
      {
	float mx2,my2;
	cpgcurs(&mx,&my,&key);
	if (key=='z')
	  {
	    cpgband(4,0,mx,my,&mx2,&my2,&key);
	    if (mx2 > mx) { minx = mx; maxx=mx2;}
	    else {minx = mx2; maxx=mx;}
	    cpgeras();
	  }
	else if (key=='u')
	  {
	    cpgeras();
	    minx = maxx = -1;
	  }
      }
      } while (interactive==1 && key!='q');
  cpglab("","",title);
  cpgend();
}

void plot7(double *plotResX,double *plotResY,double *plotResE,int n,double plotOffset,int combine,int pos,int nplot,double start,double end,double psrF0,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx)
{
  float miny,maxy;
  float borderx,bordery;
  float fx[n],fy[n],yerr1[n],yerr2[n];
  int i;
  char xlabel[100]="MJD";
  int n2=0;

  if (plotOffset!=0)
    sprintf(xlabel,"Days after MJD %.0f",plotOffset);

  for (i=0;i<n;i++)
    {
      fx[n2] = (float)(plotResX[i]-plotOffset);
      fy[n2] = (float)(plotResY[i]*1000.0);
      yerr1[n2] = fy[i] - plotResE[i]*1000.0;
      yerr2[n2] = fy[i] + plotResE[i]*1000.0;
      if (minx==maxx || (fx[n2] > minx && fx[n2] <= maxx))
	n2++;
    }
  if (minx == maxx)
    {
      minx = start-plotOffset;
      maxx = end-plotOffset;
    }
  if (yscale_set==0)
    {
      miny = TKfindMin_f(fy,n2);
      maxy = TKfindMax_f(fy,n2);
      bordery = (maxy-miny)*0.1;      
    }
  else
    {
      miny = yscale_min;
      maxy = yscale_max;
      bordery = 0;
    }
      borderx = (maxx-minx)*0.1;

  if (combine==0)
    {
      cpgenv(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery,0,1);
      cpglab(xlabel,"Residual (ms)","");
    }
  else
    {
      cpgsvp(0.1,0.95,global_footer+pos*(1.0-(global_footer+global_header))/(double)nplot,global_footer+(pos+1)*(1.0-(global_footer+global_header))/(double)nplot);
      cpgswin(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery);
    
      if (pos!=0)
	{
	  cpgbox("BCTS",0,0,"BCNTS",0,0);
	  cpglab("","Residual (ms)","");
	}
      else
	{
	  cpgbox("BCNTS",0,0,"BCNTS",0,0);
	  cpglab(xlabel,"Residual (ms)","");
	  } 
    }
  cpgpt(n2,fx,fy,16);
  cpgerry(n2,fx,yerr1,yerr2,1);
  // Draw line representing one phase
  /*  fx[0] = fx[1] = minx-borderx/2.0;
  fy[0] = -1.0/psrF0/2.0*1000.0;
  fy[1] = 1.0/psrF0/2.0*1000.0;
  cpgsci(2); cpgsls(4); cpgarro(fx[0],fy[0],fx[1],fy[1]); cpgsci(1); cpgsls(1); */
}


void plot1(double *epoch,double *f0,float *yerr1,float *yerr2,
	   int n,double plotOffset,int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx)           
{
  float miny,maxy;
  float borderx,bordery;
  float fx[n],fy[n];
  int i;
  char xlabel[100]="MJD";
  int n2=0;

  if (plotOffset!=0)
    sprintf(xlabel,"Days after MJD %.0f",plotOffset);

  for (i=0;i<n;i++)
    {
      fx[n2] = (float)epoch[i]-plotOffset;
      fy[n2] = (float)f0[i];
      if (minx==maxx || (fx[n2] > minx && fx[n2] <= maxx))
	n2++;
    }
  
  if (minx == maxx)
    {
      minx = TKfindMin_f(fx,n2);
      maxx = TKfindMax_f(fx,n2);
    }
 if (yscale_set==0)
    {
      miny = TKfindMin_f(fy,n2);
      maxy = TKfindMax_f(fy,n2);
      bordery = (maxy-miny)*0.1;
    }
  else
    {
      miny = yscale_min;
      maxy = yscale_max;
      bordery = 0;
    }
 borderx = (maxx-minx)*0.1;
  

  if (combine==0)
    {
      cpgenv(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery,0,1);
      cpglab(xlabel,"F0 (Hz)","");
    }
  else
    {
      cpgsvp(0.1,0.95,global_footer+pos*(1.0-(global_footer+global_header))/(double)nplot,global_footer+(pos+1)*(1.0-(global_footer+global_header))/(double)nplot);
      cpgswin(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery);
      if (pos!=0)
	{
	  cpgbox("BCTS",0,0,"BCNTS",0,0);
	  cpglab("","F0 (Hz)","");
	}
      else
	{
	  cpgbox("BCNTS",0,0,"BCNTS",0,0);
	  cpglab(xlabel,"F0 (Hz)","");
	}
    }
  cpgpt(n2,fx,fy,16);
  cpgerry(n2,fx,yerr1,yerr2,1);

}

void plot8(double *epoch,double *f0,float *yerr1,float *yerr2,
	   int n,double plotOffset,int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx)           
{
  float miny,maxy;
  float borderx,bordery;
  float fx[n],fy[n];
  int i;
  char xlabel[100]="MJD";
  double mean=0;
  int n2=0;

  if (plotOffset!=0)
    sprintf(xlabel,"Days after MJD %.0f",plotOffset);

  for (i=0;i<n;i++)
    {
      fx[n2] = (float)epoch[i]-plotOffset;
      if (minx==maxx || (fx[n2] > minx && fx[n2] <= maxx))
	{
	  mean+=(double)f0[n2];
	  n2++;
	}
    }
  n2=0;
  for (i=0;i<n;i++)
    {
      fy[n2] = (float)(f0[i]-mean/(double)n);
      if (minx==maxx || (epoch[i]-plotOffset > minx && epoch[i]-plotOffset <= maxx))
	n2++;
    }
  if (minx == maxx)
    {  
      minx = TKfindMin_f(fx,n2);
      maxx = TKfindMax_f(fx,n2);
    }

  if (yscale_set==0)
    {
      miny = TKfindMin_f(fy,n2);
      maxy = TKfindMax_f(fy,n2);
      bordery = (maxy-miny)*0.1;
    }
  else
    {
      miny = yscale_min;
      maxy = yscale_max;
      bordery = 0;
    }
  borderx = (maxx-minx)*0.1;
  

  if (combine==0)
    {
      cpgenv(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery,0,1);
      cpglab(xlabel,"F0 (Hz)","");
    }
  else
    {
      cpgsvp(0.1,0.95,global_footer+pos*(1.0-(global_footer+global_header))/(double)nplot,global_footer+(pos+1)*(1.0-(global_footer+global_header))/(double)nplot);
      cpgswin(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery);
      if (pos!=0)
	{
	  cpgbox("BCTS",0,0,"BCNTS",0,0);
	  cpglab("","F0 (Hz)","");
	}
      else
	{
	  cpgbox("BCNTS",0,0,"BCNTS",0,0);
	  cpglab(xlabel,"F0 (Hz)","");
	}
    }
  cpgpt(n2,fx,fy,16);
  cpgerry(n2,fx,yerr1,yerr2,1);

}

void plot2(double *epoch,double *f0,double *f0e,int *nFit,int *id,int n,double plotOffset,
	   int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx)           
{
  float miny,maxy;
  float borderx,bordery;
  double ty[n];
  float fx[n],fy[n],yerr1[n],yerr2[n];
  char str[100];
  int i;
  int n2=0;
  char xlabel[100]="MJD";

  if (plotOffset!=0)
    sprintf(xlabel,"Days after MJD %.0f",plotOffset);

  for (i=0;i<n;i++)
    ty[i] = f0[i];
  
  
  TKremovePoly_d(epoch,ty,n,2);
  for (i=0;i<n;i++)
    {
      fx[n2] = (float)epoch[i]-plotOffset;
      fy[n2] = (float)ty[i];
      yerr1[n2] = ty[i]-f0e[i];
      yerr2[n2] = ty[i]+f0e[i];
      if (minx==maxx || (fx[n2] > minx && fx[n2] <= maxx))
	n2++;
    }

  if (minx == maxx)
    {
      minx = TKfindMin_f(fx,n2);
      maxx = TKfindMax_f(fx,n2);
    }
  if (yscale_set==0)
    {
      miny = TKfindMin_f(fy,n2);
      maxy = TKfindMax_f(fy,n2);
      bordery = (maxy-miny)*0.1;
    }
  else
    {
      miny = yscale_min;
      maxy = yscale_max;
      bordery = 0;
    }
  borderx = (maxx-minx)*0.1;
  


  
  if (combine==0)
    {
      cpgenv(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery,0,1);
      cpglab(xlabel,"F0 - slope and offset (Hz)","");
    }
  else
    {
      cpgsvp(0.1,0.95,global_footer+pos*(1.0-(global_footer+global_header))/(double)nplot,global_footer+(pos+1)*(1.0-(global_footer+global_header))/(double)nplot);
      cpgswin(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery);
      if (pos!=0)
	{
	  cpgbox("BCTS",0,0,"BCNTS",0,0);
	  cpglab("","F0 - slope and offset (Hz)","");
	}
      else
	{
	  cpgbox("BCNTS",0,0,"BCNTS",0,0);
	  cpglab(xlabel,"F0 - slope and offset (Hz)","");
	  }
    }





  cpgpt(n2,fx,fy,16);
  cpgerry(n2,fx,yerr1,yerr2,1);
  cpgsci(2);
  cpgsch(0.7);
    for (i=0;i<n2;i++)
    {
      cpgsci(2);
      sprintf(str,"%d",nFit[i]);
      cpgtext(fx[i]+borderx/5,fy[i],str);
      cpgsci(3);
      sprintf(str,"%d",id[i]+1);
      cpgtext(fx[i]-borderx/4,fy[i],str);
      } 
  cpgsci(1);
  //  cpgsch(1.4);

}

void plot9(double *epoch,double *f0,double *f0e,int *nFit,int *id,int n,double plotOffset,
	   int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx)           
{
  float miny,maxy;
  float borderx,bordery;
  double ty[n];
  float fx[n],fy[n],yerr1[n],yerr2[n];
  char str[100];
  int i;
  int n2=0;
  char xlabel[100]="MJD";

  if (plotOffset!=0)
    sprintf(xlabel,"Days after MJD %.0f",plotOffset);

  for (i=0;i<n;i++)
    ty[i] = f0[i];
  
  
  TKremovePoly_d(epoch,ty,n,2);
  for (i=0;i<n;i++)
    {
      fx[n2] = (float)epoch[i]-plotOffset;
      fy[n2] = (float)ty[i]*1e9;
      yerr1[n2] = ty[i]*1e9-f0e[i]*1e9;
      yerr2[n2] = ty[i]*1e9+f0e[i]*1e9;
      if (minx==maxx || (fx[n2] > minx && fx[n2] <= maxx))
	n2++;
    }

  if (minx == maxx)
    {
      minx = TKfindMin_f(fx,n2);
      maxx = TKfindMax_f(fx,n2);
    }
  if (yscale_set==0)
    {
      miny = TKfindMin_f(fy,n2);
      maxy = TKfindMax_f(fy,n2);
      bordery = (maxy-miny)*0.1;
    }
  else
    {
      miny = yscale_min;
      maxy = yscale_max;
      bordery = 0;
    }
  borderx = (maxx-minx)*0.1;
      


  
  if (combine==0)
    {
      cpgenv(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery,0,1);
      cpglab(xlabel,"F0 - slope and offset (Hz)","");
    }
  else
    {
      cpgsvp(0.1,0.95,global_footer+pos*(1.0-(global_footer+global_header))/(double)nplot,global_footer+(pos+1)*(1.0-(global_footer+global_header))/(double)nplot);
      cpgswin(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery);
      if (pos!=0)
	{
	  cpgbox("BCTS",0,0,"BCNTS",0,0);
	  cpglab("","\\gn (10\\u-9\\d Hz)","");
	}
      else
	{
	  cpgbox("BCNTS",0,0,"BCNTS",0,0);
	  cpglab(xlabel,"\\gn (10\\u-9\\d Hz)","");
	  }
    }





  cpgpt(n2,fx,fy,16);
  cpgerry(n2,fx,yerr1,yerr2,1);
  cpgsci(1);
  //  cpgsch(1.4);

}

void plot4(double *epoch,double *f0,double *f0e,int *nFit,int *id,int n,double plotOffset,float *gt,int ngt,
	   int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx)
{
  float miny,maxy;
  float borderx,bordery;
  float mingt;
  double ty[n];
  float fx[n],fy[n],yerr1[n],yerr2[n];
  char str[100];
  double par[3];
  int i;
  char xlabel[100]="MJD";
  double vx[n],vy[n],ve[n];
  int nfit=0;
  int n2=0;

  mingt = TKfindMin_f(gt,ngt);

  if (plotOffset!=0)
    sprintf(xlabel,"Days after MJD %.0f",plotOffset);

  for (i=0;i<n;i++)
    {
      if (epoch[i] < mingt)
	{
	  vx[nfit] = epoch[i];
	  vy[nfit] = f0[i];
	  ve[nfit] = f0e[i];
	  nfit++;
	}
    }

               
  TKleastSquares_svd_noErr(vx,vy,nfit, par, 2, TKfitPoly);  
  for (i=0;i<n;i++)
    ty[i] = f0[i]-par[0]-par[1]*epoch[i];

  for (i=0;i<n;i++)
    {
      fx[n2] = (float)epoch[i]-plotOffset;
      fy[n2] = (float)ty[i]/1.0e-6;
      yerr1[n2] = ty[i]/1.0e-6-f0e[i];
      yerr2[n2] = ty[i]/1.0e-6+f0e[i];
      if (minx==maxx || (fx[n2] > minx && fx[n2] <= maxx))
	n2++;
    }


  if (minx == maxx)
    {
      minx = TKfindMin_f(fx,n2);
      maxx = TKfindMax_f(fx,n2);
    }
 if (yscale_set==0)
    {
      miny = TKfindMin_f(fy,n2);
      maxy = TKfindMax_f(fy,n2);
      bordery = (maxy-miny)*0.1;
    }
  else
    {
      miny = yscale_min;
      maxy = yscale_max;
      bordery = 0;
    }
 borderx = (maxx-minx)*0.1;
      
  

  if (combine==0)
    {
      cpgenv(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery,0,1);
      cpglab(xlabel,"\\gD\\gn(\\gmHz)","");
    }
  else
    {
      cpgsvp(0.1,0.95,global_footer+pos*(1.0-(global_footer+global_header))/(double)nplot,global_footer+(pos+1)*(1.0-(global_footer+global_header))/(double)nplot);
      cpgswin(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery);
      if (pos!=0)
	{
	  cpgbox("BCTS",0,0,"BCNTS",0,0);
	  cpglab("","\\gD\\gn(\\gmHz)","");
	}
      else
	{
	  cpgbox("BCNTS",0,0,"BCNTS",0,0);
	  cpglab(xlabel,"\\gD\\gn(\\gmHz)","");
	}
    }

  for (i=0;i<n2;i++)
    {
      if (f0e[i]>0)
	{
	  cpgpt(1,fx+i,fy+i,16);
	  cpgerry(1,fx+i,yerr1+i,yerr2+i,1);
	}
      else
	{
	  cpgsci(2);
	  cpgpt(1,fx+i,fy+i,16);
	  cpgsci(1);
	}
    }
}

void plot5(double *epoch,double *f0,double *f0e,int *nFit,int *id,int n,double plotOffset,float *gt,int ngt,
	   int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx)
{
  float miny,maxy;
  float borderx,bordery;
  float mingt;
  double ty[n];
  float fx[n],fy[n],yerr1[n],yerr2[n];
  char str[100];
  double par[3];
  int i;
  char xlabel[100]="MJD";
  double vx[n],vy[n],ve[n];
  double start = 0;
  double end;
  double maxTime;
  int nfit=0;
  int k;
  int n2=0;

  maxTime = TKfindMax_d(epoch,n);

  mingt = TKfindMin_f(gt,ngt);

  for (i=0;i<n;i++)
    {
      if (epoch[i] < mingt)
	{
	  vx[nfit] = epoch[i];
	  vy[nfit] = f0[i];
	  ve[nfit] = f0e[i];
	  nfit++;
	}
    }

               
  TKleastSquares_svd_noErr(vx,vy,nfit, par, 2, TKfitPoly);  
  for (i=0;i<n;i++)
    ty[i] = f0[i]-par[0]-par[1]*epoch[i];


  for (k=0;k<=ngt;k++)
    {
      if (k==ngt)
	end = maxTime+1;
      else
	end = gt[k];

      if (k==0)
	start = 0;
      else
	start = gt[k-1];

      nfit=0;
      for (i=0;i<n;i++)
	{
	  if (epoch[i] > start && epoch[i] <= end)
	    {
	      vx[nfit] = epoch[i];
	      vy[nfit] = ty[i];
	      ve[nfit] = f0e[i];
	      nfit++;
	    }
	}
      TKleastSquares_svd_noErr(vx,vy,nfit, par, 1, TKfitPoly);  
      for (i=0;i<n;i++)
	{
	  if (epoch[i] > start && epoch[i] <= end)
	    {
	      fx[i] = (float)epoch[i]-plotOffset;
	      //	      fy[i] = (float)(f0[i]-par[0]-par[1]*epoch[i]);
	      fy[i] = (float)(ty[i]-par[0])/1.0e-6;
	      yerr1[i] = fy[i]/1.0e-6-f0e[i];
	      yerr2[i] = fy[i]/1.0e-6+f0e[i];
	    }
	}
      
    }


  if (plotOffset!=0)
    sprintf(xlabel,"Days after MJD %.0f",plotOffset);


               


  if (minx == maxx)
    {
      minx = TKfindMin_f(fx,n);
      maxx = TKfindMax_f(fx,n);
    }
  if (yscale_set==0)
    {
      if (minx==maxx)
	{
	  miny = TKfindMin_f(fy,n);
	  maxy = TKfindMax_f(fy,n);	  
	}
      else
	{
	  int t=1;
	  for (i=0;i<n;i++)
	    {
	      if (fx[i] > minx && fx[i] <= maxx)
		{
		  if (t==1){miny = fy[i]; maxy = fy[i]; t=2;}
		  if (miny > fy[i]) miny = fy[i];
		  if (maxy < fy[i]) maxy = fy[i];
		}
	    }
	  printf("In here with %g %g\n",miny,maxy);
	}

      bordery = (maxy-miny)*0.1;
    }
  else
    {
      miny = yscale_min;
      maxy = yscale_max;
      bordery = 0;
    }
 borderx = (maxx-minx)*0.1;
  
  

  if (combine==0)
    {
      cpgenv(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery,0,1);
      cpglab(xlabel,"\\gD\\gn(\\gmHz)","");
    }
  else
    {
      cpgsvp(0.1,0.95,global_footer+pos*(1.0-(global_footer+global_header))/(double)nplot,global_footer+(pos+1)*(1.0-(global_footer+global_header))/(double)nplot);
      cpgswin(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery);
      if (pos!=0)
	{
	  cpgbox("BCTS",0,0,"BCNTS",0,0);
	  cpglab("","\\gD\\gn(\\gmHz)","");
	}
      else
	{
	  cpgbox("BCNTS",0,0,"BCNTS",0,0);
	  cpglab(xlabel,"\\gD\\gn(\\gmHz)","");
	}
    }


  cpgpt(n,fx,fy,16);
  cpgerry(n,fx,yerr1,yerr2,1);
}


void plot3(double *epoch,double *f1,double *f1e,int *nFit,int *id,int n,double plotOffset,int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx)           
{
  float miny,maxy;
  float borderx,bordery;
  double ty[n];
  float fx[n],fy[n],yerr1[n],yerr2[n];
  char str[100];
  int i;
  int n2=0;
  char xlabel[100]="MJD";

  if (plotOffset!=0)
    sprintf(xlabel,"Days after MJD %.0f",plotOffset);

  for (i=0;i<n;i++)
    ty[i] = f1[i];
    
  //  TKremovePoly_d(epoch,ty,n,2);
  for (i=0;i<n;i++)
    {
      fx[n2] = (float)epoch[i]-plotOffset;
      fy[n2] = (float)ty[i];
      yerr1[n2] = ty[i]-f1e[i];
      yerr2[n2] = ty[i]+f1e[i];
      if (minx==maxx || (fx[n2] > minx && fx[n2] <= maxx))
	n2++;
    }


  if (minx == maxx)
    {
      minx = TKfindMin_f(fx,n2);
      maxx = TKfindMax_f(fx,n2);
    }

  if (yscale_set==0)
    {
      miny = TKfindMin_f(fy,n2);
      maxy = TKfindMax_f(fy,n2);
      bordery = (maxy-miny)*0.1;
    }
  else
    {
      miny = yscale_min;
      maxy = yscale_max;
      bordery = 0;
    }
 borderx = (maxx-minx)*0.1;
  
  
  if (combine==0)
    {
      cpgenv(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery,0,1);
      cpglab(xlabel,"F1","");
    }
  else
    {
      cpgsvp(0.1,0.95,global_footer+pos*(1.0-(global_footer+global_header))/(double)nplot,global_footer+(pos+1)*(1.0-(global_footer+global_header))/(double)nplot);
      cpgswin(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery);
      if (pos!=0)
	{
	  cpgbox("BCTS",0,0,"BCNTS",0,0);
	  cpglab("","F1","");
	}
      else
	{
	  cpgbox("BCNTS",0,0,"BCNTS",0,0);
	  cpglab(xlabel,"F1","");
	}
    }
  cpgpt(n2,fx,fy,16);
  cpgerry(n2,fx,yerr1,yerr2,1);
}

void plot6(double *epoch,double *f1,double *f1e,int *nFit,int *id,int n,double plotOffset,int combine,int pos,int nplot,float yscale_min,float yscale_max,int yscale_set,float minx,float maxx)           
{
  float miny,maxy;
  float borderx,bordery;
  double ty[n];
  float fx[n],fy[n],yerr1[n],yerr2[n];
  char str[100];
  int i;
  char xlabel[100]="MJD";
  char ylabel[100];
  double mean;
  int n2=0;

  if (plotOffset!=0)
    sprintf(xlabel,"Days after MJD %.0f",plotOffset);

  for (i=0;i<n;i++)
    ty[i] = f1[i]/1.0e-15;
    
  mean = TKmean_d(ty,n);
  for (i=0;i<n;i++)
    {
      fx[n2] = (float)epoch[i]-plotOffset;
      fy[n2] = (float)(ty[i]-mean);
      yerr1[n2] = (ty[i]-mean)-f1e[i]/1.0e-15;
      yerr2[n2] = (ty[i]-mean)+f1e[i]/1.0e-15;
      if (minx==maxx || (fx[n2] > minx && fx[n2] <= maxx))
	n2++;
    }


  if (minx == maxx)
    {
      minx = TKfindMin_f(fx,n2);
      maxx = TKfindMax_f(fx,n2);
    }
if (yscale_set==0)
    {
      miny = TKfindMin_f(fy,n2);
      maxy = TKfindMax_f(fy,n2);
      bordery = (maxy-miny)*0.1;
    }
  else
    {
      miny = yscale_min;
      maxy = yscale_max;
      bordery = 0;
    }
 borderx = (maxx-minx)*0.1;
  
  
  if (combine==0)
    {
      cpgenv(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery,0,1);
      sprintf(ylabel,"\\gn\\u\\u\\b\\b . \\d\\d%+.0f(10\\u-15\\d s\\u-2\\d)",-mean);
      cpglab(xlabel,ylabel,"");
    }
  else
    {
      cpgsvp(0.1,0.95,global_footer+pos*(1.0-(global_footer+global_header))/(double)nplot,global_footer+(pos+1)*(1.0-(global_footer+global_header))/(double)nplot);
      cpgswin(minx-borderx,maxx+borderx,miny-bordery,maxy+bordery);
      if (pos!=0)
	{
	  cpgbox("BCTS",0,0,"BCNTS",0,0);
	  sprintf(ylabel,"\\gn\\u\\u\\b\\b . \\d\\d%+.0f(10\\u-15\\d s\\u-2\\d)",-mean);
	  cpglab("",ylabel,"");
	}
      else
	{
	  cpgbox("BCNTS",0,0,"BCNTS",0,0);
	  sprintf(ylabel,"\\gn\\u\\u\\b\\b . \\d\\d%+.0f(10\\u-15\\d s\\u-2\\d)",-mean);
	  cpglab(xlabel,ylabel,"");

	}
    }

  for (i=0;i<n2;i++)
    {
      if (f1e[i] > 0)
	{
	  cpgpt(1,fx+i,fy+i,16);
	  cpgerry(1,fx+i,yerr1+i,yerr2+i,1);
	}
      else
	{
	  cpgsci(2);
	  cpgpt(1,fx+i,fy+i,16);
	  cpgsci(1);
	}
    }
}

void interactivePlot(double *epoch,double *f0,double *f0e,int *nFit,int *id,int n)
{
  float fx[n],fy[n],yerr1[n],yerr2[n];
  double dx[n],dy[n],de[n];
  double dx2[n],dy2[n],de2[n];
  double olda,oldb;
  float minx,maxx,miny,maxy;
  int i,j,k;
  float mx,my,border;
  float start=-1;
  float end=-1;
  int istart=0;
  int iend=0;
  long double removeA=0.0,removeB=0.0,removeC=0.0; // slope removal
  double f0val,f1val,f2val,forig;
  double toffset=epoch[0];
  double dt;
  double model;
  char key='a';
  glitchS glitch[10];
  double **cvm,chisq;
  int nglt=0;
  int fitf0=0;
  int fitf1=0;
  
  cvm = (double **)malloc(sizeof(double *)*10);
  for (i=0;i<10;i++)
    cvm[i] = (double *)malloc(sizeof(double));

  f0val = 0.0; //forig = f0[0];
  f1val = 0.0;
  f2val = 0.0;

  for (i=0;i<n;i++)
    {
      dx[i] = epoch[i];
      dy[i] = f0[i];
      de[i] = f0e[i];
    }

  cpgbeg(0,"/xs",1,1);
  cpgask(0);
  while (key!='q')
    {
      printf("q             - quit\n");
      printf("middle button - add glitch\n");
      printf("F             - non-linear fit\n");
      printf("-----------------------------------\n");
      printf("F0  %.10f\n",f0val-forig);
      printf("F1  %g\n",f1val);
      printf("F2  %g\n",f2val);
      printf("\n");
      for (i=0;i<nglt;i++)
	{
	  printf("GLEP_%d %.5f\n",i+1,glitch[i].glep);
	  printf("GLPH_%d 0\n",i+1);
	  printf("GLF0_%d %.8g %d\n",i+1,glitch[i].glf0,glitch[i].fitf0);
	  printf("GLF1_%d %.8g %d\n",i+1,glitch[i].glf1,glitch[i].fitf1);
	  printf("GLF0D_%d %.8g %d\n",i+1,glitch[i].glf0d,glitch[i].fitf0d);
	  printf("GLTD_%d %.8g %d\n",i+1,glitch[i].gltd,glitch[i].fittd);
	}
      printf("-----------------------------------\n");

      for (i=0;i<n;i++)
	{
	  fx[i] = (float)(dx[i]-toffset);
	  dx2[i] = dx[i];
	  dy2[i] = dy[i];

	  if (removeA!=0 || removeB!=0)
	    {
	      dy2[i] = (float)(dy2[i] - (removeA + removeB*(dx[i]-toffset) + 
					 removeC*pow(dx[i]-toffset,2)));
	    }
	}
      // Do glitch model
      for (i=0;i<n;i++)
	{
	  model = f0val+f1val*dx[i]*86400.0;
	  for (k=0;k<nglt;k++)
	    {
	      
	      if (dx[i] > glitch[k].glep)
		{
		  dt = dx[i]-glitch[k].glep;
		  model += glitch[k].glf0+glitch[k].glf1*(dt*86400.0)+
		    glitch[k].glf0d*exp(-dt/glitch[k].gltd);
		}
	    }
	  dy2[i] -= model;
		  
	}
      for (i=0;i<n;i++)
	{
	  yerr1[i] = (float)(dy2[i]-de[i]);
	  yerr2[i] = (float)(dy2[i]+de[i]);	  
	  fy[i] = (float)dy2[i];	  
	}
      minx = TKfindMin_f(fx,n);
      maxx = TKfindMax_f(fx,n);
      miny = TKfindMin_f(fy,n);
      maxy = TKfindMax_f(fy,n);      
      border = (maxy-miny)*0.1; maxy += border; miny -= border;
      border = (maxx-minx)*0.1; maxx += border; minx -= border;

      if (istart==0)
	start = minx;
      if (iend==0)
	end = maxx;

      cpgenv(minx,maxx,miny,maxy,0,0);
      cpglab("MJD","F0","");
      cpgpt(n,fx,fy,16);
      cpgerry(n,fx,yerr1,yerr2,1);

      for (i=0;i<nglt;i++)
	{
	  fx[0] = fx[1] = (float)(glitch[i].glep-toffset);
	  fy[0] = miny; fy[1] = maxy;
	  cpgsci(7); cpgsls(4); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);	  
	}

      if (istart==1)
	{
	  fx[0] = fx[1] = start;
	  fy[0] = miny; fy[1] = maxy;
	  cpgsci(3); cpgsls(4); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
	}
      if (iend==1)
	{
	  fx[0] = fx[1] = end;
	  fy[0] = miny; fy[1] = maxy;
	  cpgsci(3); cpgsls(4); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
	}
      /*      if (removeA!=0 || removeB!=0)
	{
	  for (i=0;i<4;i++)
	    {
	      fx[i] = minx+(maxx-minx)*(float)i/4.0;
	      fy[i] =(removeA + removeB*fx[i]);
	    }
	  cpgline(4,fx,fy);
	  }*/

      drawMenu(minx,maxx,miny,maxy,glitch,nglt,fitf0,fitf1);
      cpgcurs(&mx,&my,&key);
      if (my > maxy)
	{
	  if (key=='A' || key=='X')
	    checkMenu(minx,maxx,miny,maxy,glitch,nglt,mx,my,&fitf0,&fitf1,key);
	}
      else
	{
	  if (key=='1' || key=='2' || key=='3')
	    {
	      double vx[n],vy[n],ve[n];
	      double par[3];
	      int m;
	      int npt=0;
	      
	      if (key=='1') m=1;
	      if (key=='2') m=2;
	      if (key=='3') m=3;
	      
	      for (i=0;i<n;i++)
		{
		  if (dx[i]-toffset > start && dx[i]-toffset < end)
		    {
		      vx[npt] = dx[i]-toffset;
		      vy[npt] = dy2[i];
		      ve[npt] = de[i];
		      npt++;
		    }
		}
	      printf("Fitting with %d points and polynomial order %d\n",npt,m);
	      TKleastSquares_svd_noErr(vx,vy,npt, par, m, TKfitPoly);  
	      printf("Result: %g %g %g\n",par[0],par[1],par[2]);
	      f0val += par[0];
	      if (m>1) f1val += (par[1]/86400.0);
	      if (m>2) f2val += (par[2]/86400.0/86400.0);
	      //	      removeA += par[0];
	      //	      if (m==2 || m==3) removeB += par[1];
	      //	      if (m==3) removeC += par[2];
	    }
	  else if (key=='d')
	    defineGlitchVal(glitch,nglt);
	  else if (key=='c')
	    changeFit(glitch,nglt);
	  else if (key=='s')
	    {
	      start=mx; istart=1;
	    }
	  else if (key=='e')
	    {
	      iend=1;
	      end=mx;
	    }
	      else if (key=='o') // Overlay
		{
		  float tx[100],ty[100];
		  double dt;
		  double func=0.0;
		  int endit=0;
		  int yesno=0;
		  do {
		    printf("Enter parameters: ");
		    scanf("%lf %lf %lf %lf %lf",&glitch[0].glep,&glitch[0].glf0,&glitch[0].glf1,&glitch[0].glf0d,&glitch[0].gltd);
		    func=0.0;
		    for (i=0;i<100;i++)
		      {
			tx[i] = minx+(maxx-minx)/100.0*i;
			dt = tx[i]+toffset - glitch[0].glep; 
			if (dt > 0)
			  {
			    func = glitch[0].glf0+glitch[0].glf1*(dt*86400.0)+
			      glitch[0].glf0d*exp(-dt/glitch[0].gltd);
			  }
			ty[i] = func;
		      }
		    cpgsls(4); cpgsci(3); cpgline(100,tx,ty); cpgsci(1); cpgsls(1);
		    printf("0 to try again, 1 to stop\n");
		    scanf("%d",&endit);
		  } while (endit==0);
		}
	  else if (key=='f' || key=='F') // Fit!
	    {
	      double pval[10],eval[10];
	      int nfit=0;
	      
	      global_glitch = glitch;
	      global_nglt = nglt;
	      global_fitf0 = fitf0;
	      global_fitf1 = fitf1;
	      global_valf0 = f0val;
	      global_valf1 = f1val;

	      // Count parameters to fit
	      printf("Setting %g %g\n",f0val,f1val);
	      if (fitf0 == 1) {pval[nfit] = f0val; nfit++;}
	      if (fitf1 == 1) {pval[nfit] = f1val; nfit++;}
	      for (i=0;i<nglt;i++)
		{
		  if (glitch[i].fitf0==1) {pval[nfit] = glitch[i].glf0;nfit++;}
		  if (glitch[i].fitf1==1) {pval[nfit] = glitch[i].glf1; nfit++;}
		  if (glitch[i].fitf0d==1) {pval[nfit] = glitch[i].glf0d; nfit++;}
		  if (glitch[i].fittd==1) {pval[nfit] = glitch[i].gltd; nfit++;}
		}
	      printf("Fitting for %d parameters\n",nfit);
	      if (key=='f')
		{
		  TKleastSquares_svd(dx,dy2,de,n,pval,eval,nfit,cvm,&chisq,fitFuncs,1);
		}
	      else if (key=='F')
		{
		  lm_status_struct status;
		  lm_control_struct control = lm_control_double;
		  control.printflags = 3; // monitor status (+1) and parameters (+2)

		  //		  lmcurve_fit( nfit, pval, n, dx, dy2, nonlinearFunc, &control, &status );
		  lmcurve_fit( nfit, pval, n, dx, dy, nonlinearFunc, &control, &status );
		  
		  for (i=0;i<nfit;i++)
		    printf("Parameter %d = %g\n",i,pval[i]);

		}
	      for (i=0;i<nfit;i++)
		printf("%d %g %g\n",i+1,pval[i],eval[i]);
	      printf("red chisq = %g\n",chisq/(double)(n-nfit));
	      // Update parameters
	      nfit=0;
	      if (fitf0==1)
		{
		  if (key=='f')
		    f0val+=pval[nfit];
		  else
		    f0val=pval[nfit];
		  nfit++;
		}
	      if (fitf1==1)
		{
		  if (key=='f')
		    f1val+=pval[nfit]/86400.0;
		  else
		    f1val=pval[nfit]/86400.0;
		  nfit++;
		}
	      for (i=0;i<nglt;i++)
		{
		  if (glitch[i].fitf0==1)
		    {
		      if (key=='f') glitch[i].glf0 += pval[nfit];
		      else glitch[i].glf0 = pval[nfit];
		      nfit++;
		    }
		  if (glitch[i].fitf1==1) 
		    {
		      if (key=='f')
			glitch[i].glf1 += pval[nfit]/86400.0;
		      else
			glitch[i].glf1 = pval[nfit]/86400.0;
		      nfit++;
		    }
		  if (glitch[i].fitf0d==1) 
		    {
		      if (key=='f')
			glitch[i].glf0d += pval[nfit];
		      else
			glitch[i].glf0d = pval[nfit];
		      nfit++;
		    }
		  if (glitch[i].fittd==1) 
		    {
		      if (key=='f')
			glitch[i].gltd += pval[nfit];
		      else
			glitch[i].gltd = pval[nfit]; 
		      nfit++;
		    }
		}
	    }
	  else if (key=='D') // Middle button -- add glitch
	    {
	      glitch[nglt].glep = (double)(mx+toffset);
	      glitch[nglt].glf0 = 0.0;
	      glitch[nglt].glph = 0.0;
	      glitch[nglt].glf1 = 0.0;
	      glitch[nglt].glf0d = 0.0;
	      glitch[nglt].gltd = 0.0;
	      
	      glitch[nglt].fitf0 = 0;
	      glitch[nglt].fitph = 0;
	      glitch[nglt].fitf1 = 0;
	      glitch[nglt].fitf0d =0;
	      glitch[nglt].fittd = 0;
	      
	      nglt++;
	    }
	  else if (key=='q')
	    printf("Goodbye\n");
	  else
	    printf("Unknown key >%d< >%c<\n",(int)key,key);
	}
    }
  cpgend();

}

void checkMenu(float minx,float maxx,float miny,float maxy,glitchS *glitch, int nglt,
	       float mx,float my,int *fitf0,int *fitf1,char key)
{
  int xpos,ypos;
  float cheight,cwidth;
  double newval;

  cpgqcs(4,&cwidth,&cheight);

  xpos = (int)((mx-minx)/((maxx-minx)/6.0));
  ypos = (int)((maxy+(maxy-miny)*0.1+cheight-my)/((maxy-miny)*0.05));

  if (key=='X')
    {
      printf("New value ");
      scanf("%lf",&newval);
      if (xpos==2 && ypos==0)
	glitch[0].glf0 = newval;
      if (xpos==3 && ypos==0)
	glitch[0].glf1 = newval;
      if (xpos==4 && ypos==0)
	glitch[0].glf0d = newval;
      if (xpos==5 && ypos==0)
	glitch[0].gltd = newval;
    }
  else
    {
      if (xpos==0 && ypos==0)
	{
	  if (*fitf0 == 1)
	    *fitf0 = 0;
	  else
	    *fitf0 = 1;
	}
      if (xpos==1 && ypos==0)
	{
	  if (*fitf1 == 1)
	    *fitf1 = 0;
	  else
	    *fitf1 = 1;
	}
      
      if (xpos==2 && ypos==0)
	{
	  if (glitch[0].fitf0==1)
	    glitch[0].fitf0=0;
	  else
	    glitch[0].fitf0=1;
	}
      else if (xpos==3 && ypos==0)
	{
	  if (glitch[0].fitf1==1)
	    glitch[0].fitf1=0;
	  else
	    glitch[0].fitf1=1;
	}
      else if (xpos==4 && ypos==0)
	{
	  if (glitch[0].fitf0d==1)
	    glitch[0].fitf0d=0;
	  else
	    glitch[0].fitf0d=1;
	}
      else if (xpos==5 && ypos==0)
	{
	  if (glitch[0].fittd==1)
	    glitch[0].fittd=0;
	  else
	    glitch[0].fittd=1;
	}
    } 
}

void drawMenu(float minx,float maxx,float miny,float maxy,glitchS *glitch, int nglt,
	      int fitf0,int fitf1)
{
  int i;
  float xpos,ypos;
  char str[100];
  xpos = minx;
  ypos = maxy+(maxy-miny)*0.1;
  if (fitf0==1) cpgsci(2); else cpgsci(1);
  cpgtext(xpos,ypos,"F0"); xpos+=(maxx-minx)/6.0;
  if (fitf1==1) cpgsci(2); else cpgsci(1);
  cpgtext(xpos,ypos,"F1"); xpos+=(maxx-minx)/6.0;
  for (i=0;i<nglt;i++)
    {
      if (glitch[i].fitf0==1) cpgsci(2); else cpgsci(1);
      sprintf(str,"GLF0_%d",i+1); cpgtext(xpos,ypos,str); xpos+=(maxx-minx)/6.0;
      if (glitch[i].fitf1==1) cpgsci(2); else cpgsci(1);
      sprintf(str,"GLF1_%d",i+1);cpgtext(xpos,ypos,str); xpos+=(maxx-minx)/6.0;
      if (glitch[i].fitf0d==1) cpgsci(2); else cpgsci(1);
      sprintf(str,"GLF0D_%d",i+1);cpgtext(xpos,ypos,str); xpos+=(maxx-minx)/6.0;
      if (glitch[i].fittd==1) cpgsci(2); else cpgsci(1);
      sprintf(str,"GLTD_%d",i+1);cpgtext(xpos,ypos,str); xpos+=(maxx-minx)/6.0;
      xpos=minx;
      ypos-=(maxy-miny)*0.05;
    }
  cpgsci(1);
}

void defineGlitchVal(glitchS *glitch,int nglt)
{
  int g,k;
  double newval;

  printf("Which glitch to modify (1->%d)? ",nglt);
  scanf("%d",&g);
  g--;
  printf("Press number 1-5 to modify the following values:\n");
  printf("1 GLEP %.5f\n",glitch[g].glep);
  printf("2 GLF0 %10g\n",glitch[g].glf0);
  printf("3 GLF1 %10g\n",glitch[g].glf1);
  printf("4 GLF0D %10g\n",glitch[g].glf0d);
  printf("5 GLTD %10g\n",glitch[g].gltd);
  printf("6 all\n");
  scanf("%d",&k);
  if (k!=6)
    {
      printf("Enter new value: ");
      scanf("%lf",&newval);
    }
  if (k==1) glitch[g].glep=newval;
  else if (k==2) glitch[g].glf0=newval;
  else if (k==3) glitch[g].glf1=newval;
  else if (k==4) glitch[g].glf0d=newval;
  else if (k==5) glitch[g].gltd=newval;
  else if (k==6)
    {
      printf("Enter GLEP GLF0 GLF1 GLF0D GLTD\n");
      scanf("%lf %lf %lf %lf %lf",&glitch[g].glep,&glitch[g].glf0,&glitch[g].glf1,&glitch[g].glf0d,&glitch[g].gltd);
    }
}


void changeFit(glitchS *glitch,int nglt)
{
  int g,k;
  double newval;

  printf("Which glitch to modify (1->%d)? ",nglt);
  scanf("%d",&g);
  g--;
  printf("Press number 1-5 to modify whether fitting is on or off:\n");
  if (glitch[g].fitf0==1)
    printf("1 GLF0 ON\n");
  else
    printf("1 GLF0 OFF\n");

  if (glitch[g].fitf1==1)
    printf("2 GLF1 ON\n");
  else
    printf("2 GLF1 OFF\n");

  if (glitch[g].fitf0d==1)
    printf("3 GLF0D ON\n");
  else
    printf("3 GLF0D OFF\n");

  if (glitch[g].fittd==1)
    printf("4 GLTD ON\n");
  else
    printf("4 GLTD OFF\n");
  scanf("%d",&k);
  if (k==1) 
    {
      if (glitch[g].fitf0==1)
	glitch[g].fitf0=0;
      else
	glitch[g].fitf0=1;
    }
  else if (k==2) 
    {
      if (glitch[g].fitf1==1)
	glitch[g].fitf1=0;
      else
	glitch[g].fitf1=1;
    }
  else if (k==3) 
    {
      if (glitch[g].fitf0d==1)
	glitch[g].fitf0d=0;
      else
	glitch[g].fitf0d=1;
    }
  else if (k==4) 
    {
      if (glitch[g].fittd==1)
	glitch[g].fittd=0;
      else
	glitch[g].fittd=1;
    }
}

void fitFuncs(double x,double *p,int m)
{
  int i;
  double dt,retF0=0.0,retF1=0.0,retF0D=0.0,retTD=0.0;
  int nfit=0;

  if (global_fitf0==1)
    {
      p[nfit] = 1;
      nfit++;
    }
  if (global_fitf1==1)
    {
      p[nfit] = x;
      nfit++;
    }

  dt = x - global_glitch[0].glep;
  if (dt > 0)
    {
      retF0 = 1;
      retF1 = dt;
      retF0D = exp(-dt/global_glitch[0].gltd);
      retTD = global_glitch[0].glf0d*exp(-dt/global_glitch[0].gltd)*dt/global_glitch[0].gltd/global_glitch[0].gltd;
    }
  else
    {
      retF0 = 0;
      retF1 = 0;
      retF0D = 0;
      retTD = 0;
    }
  // Find the parameters
  for (i=0;i<global_nglt;i++)
    {
      if (global_glitch[i].fitf0==1) 
	{
	  p[nfit] = retF0;
	  nfit++;
	}
      if (global_glitch[i].fitf1==1) 
	{
	  p[nfit] = retF1;
	  nfit++;
	}
      if (global_glitch[i].fitf0d==1) 
	{
	  p[nfit] = retF0D;
	  nfit++;
	}
      if (global_glitch[i].fittd==1)
	{
	  p[nfit] = retTD;
	  nfit++;
	}
    }
}


// -------------
// LM functions
// -------------
void lmcurve_evaluate( const double *par, int m_dat, const void *data,
                       double *fvec, int *info )
{
    int i;
    for ( i = 0; i < m_dat; i++ )
	fvec[i] =
            ((lmcurve_data_struct*)data)->y[i] -
            ((lmcurve_data_struct*)data)->f(
                ((lmcurve_data_struct*)data)->t[i], par ,i);
    // *info = *info; /* to prevent a 'unused variable' warning */
}


void lmcurve_fit( int n_par, double *par, int m_dat, 
                  const double *t, const double *y,
                  double (*f)( double t, const double *par, int obsNum ),
                  const lm_control_struct *control, lm_status_struct *status )
{
    lmcurve_data_struct data = { t, y, f };

    lmmin( n_par, par, m_dat, (const void*) &data,
           lmcurve_evaluate, control, status, lm_printout_std );
}



void lm_printout_std( int n_par, const double *par, int m_dat,
                      const void *data, const double *fvec,
                      int printflags, int iflag, int iter, int nfev)
/*
 *       data  : for soft control of printout behaviour, add control
 *                 variables to the data struct
 *       iflag : 0 (init) 1 (outer loop) 2(inner loop) -1(terminated)
 *       iter  : outer loop counter
 *       nfev  : number of calls to *evaluate
 */
{
    if( !printflags )
        return;

    int i;

    if( printflags & 1 ){
        /* location of printout call within lmdif */
        if (iflag == 2) {
            printf("trying step in gradient direction  ");
        } else if (iflag == 1) {
            printf("determining gradient (iteration %2d)", iter);
        } else if (iflag == 0) {
            printf("starting minimization              ");
        } else if (iflag == -1) {
            printf("terminated after %3d evaluations   ", nfev);
        }
    }

    if( printflags & 2 ){
        printf("  par: ");
        for (i = 0; i < n_par; ++i)
            printf(" %18.11g", par[i]);
        printf(" => norm: %18.11g", lm_enorm(m_dat, fvec));
    }

    if( printflags & 3 )
        printf( "\n" );

    if ( (printflags & 8) || ((printflags & 4) && iflag == -1) ) {
	printf("  residuals:\n");
	for (i = 0; i < m_dat; ++i)
	    printf("    fvec[%2d]=%12g\n", i, fvec[i] );
    }
}


/*****************************************************************************/
/*  lm_minimize (intermediate-level interface)                               */
/*****************************************************************************/

void lmmin( int n_par, double *par, int m_dat, const void *data, 
            void (*evaluate) (const double *par, int m_dat, const void *data,
                              double *fvec, int *info),
            const lm_control_struct *control, lm_status_struct *status,
            void (*printout) (int n_par, const double *par, int m_dat,
                              const void *data, const double *fvec,
                              int printflags, int iflag, int iter, int nfev) )
{

/*** allocate work space. ***/

    double *fvec, *diag, *fjac, *qtf, *wa1, *wa2, *wa3, *wa4;
    int *ipvt;

    int n = n_par;
    int m = m_dat;

    if ( (fvec = (double *) malloc(m * sizeof(double))) == NULL ||
	 (diag = (double *) malloc(n * sizeof(double))) == NULL ||
	 (qtf  = (double *) malloc(n * sizeof(double))) == NULL ||
	 (fjac = (double *) malloc(n*m*sizeof(double))) == NULL ||
	 (wa1  = (double *) malloc(n * sizeof(double))) == NULL ||
	 (wa2  = (double *) malloc(n * sizeof(double))) == NULL ||
	 (wa3  = (double *) malloc(n * sizeof(double))) == NULL ||
	 (wa4  = (double *) malloc(m * sizeof(double))) == NULL ||
	 (ipvt = (int *)    malloc(n * sizeof(int)   )) == NULL    ) {
	status->info = 9;
	return;
    }

    int j;
    if( ! control->scale_diag )
        for( j=0; j<n_par; ++j )
            diag[j] = 1;

/*** perform fit. ***/

    status->info = 0;

    /* this goes through the modified legacy interface: */
    lm_lmdif( m, n, par, fvec, control->ftol, control->xtol, control->gtol,
              control->maxcall * (n + 1), control->epsilon, diag,
              ( control->scale_diag ? 1 : 2 ),
              control->stepbound, &(status->info),
              &(status->nfev), fjac, ipvt, qtf, wa1, wa2, wa3, wa4,
              evaluate, printout, control->printflags, data );

    if ( printout )
        (*printout)( n, par, m, data, fvec,
                     control->printflags, -1, 0, status->nfev );
    status->fnorm = lm_enorm(m, fvec);
    if ( status->info < 0 )
	status->info = 11;

/*** clean up. ***/

    free(fvec);
    free(diag);
    free(qtf);
    free(fjac);
    free(wa1);
    free(wa2);
    free(wa3);
    free(wa4);
    free(ipvt);
} /*** lm_minimize. ***/

/*****************************************************************************/
/*  lm_enorm (Euclidean norm)                                                */
/*****************************************************************************/

double lm_enorm(int n, const double *x)
{
/*     Given an n-vector x, this function calculates the
 *     euclidean norm of x.
 *
 *     The euclidean norm is computed by accumulating the sum of
 *     squares in three different sums. The sums of squares for the
 *     small and large components are scaled so that no overflows
 *     occur. Non-destructive underflows are permitted. Underflows
 *     and overflows do not occur in the computation of the unscaled
 *     sum of squares for the intermediate components.
 *     The definitions of small, intermediate and large components
 *     depend on two constants, LM_SQRT_DWARF and LM_SQRT_GIANT. The main
 *     restrictions on these constants are that LM_SQRT_DWARF**2 not
 *     underflow and LM_SQRT_GIANT**2 not overflow.
 *
 *     Parameters
 *
 *	n is a positive integer input variable.
 *
 *	x is an input array of length n.
 */
    int i;
    double agiant, s1, s2, s3, xabs, x1max, x3max, temp;

    s1 = 0;
    s2 = 0;
    s3 = 0;
    x1max = 0;
    x3max = 0;
    agiant = LM_SQRT_GIANT / n;

    /** sum squares. **/

    for (i = 0; i < n; i++) {
	xabs = fabs(x[i]);
	if (xabs > LM_SQRT_DWARF) {
            if ( xabs < agiant ) {
                s2 += xabs * xabs;
            } else if ( xabs > x1max ) {
		temp = x1max / xabs;
		s1 = 1 + s1 * SQR(temp);
		x1max = xabs;
	    } else {
		temp = xabs / x1max;
		s1 += SQR(temp);
	    }
	} else if ( xabs > x3max ) {
	    temp = x3max / xabs;
	    s3 = 1 + s3 * SQR(temp);
	    x3max = xabs;
	} else if (xabs != 0.) {
            temp = xabs / x3max;
            s3 += SQR(temp);
	}
    }

    /** calculation of norm. **/

    if (s1 != 0)
	return x1max * sqrt(s1 + (s2 / x1max) / x1max);
    else if (s2 != 0)
        if (s2 >= x3max)
            return sqrt(s2 * (1 + (x3max / s2) * (x3max * s3)));
        else
            return sqrt(x3max * ((s2 / x3max) + (x3max * s3)));
    else
        return x3max * sqrt(s3);

} /*** lm_enorm. ***/

void lm_lmdif( int m, int n, double *x, double *fvec, double ftol,
	       double xtol, double gtol, int maxfev, double epsfcn,
	       double *diag, int mode, double factor, int *info, int *nfev,
	       double *fjac, int *ipvt, double *qtf, double *wa1,
	       double *wa2, double *wa3, double *wa4,
               void (*evaluate) (const double *par, int m_dat, const void *data,
                                 double *fvec, int *info),
               void (*printout) (int n_par, const double *par, int m_dat,
                                 const void *data, const double *fvec,
                                 int printflags, int iflag, int iter, int nfev),
	       int printflags, const void *data )
{
/*
 *   The purpose of lmdif is to minimize the sum of the squares of
 *   m nonlinear functions in n variables by a modification of
 *   the levenberg-marquardt algorithm. The user must provide a
 *   subroutine evaluate which calculates the functions. The jacobian
 *   is then calculated by a forward-difference approximation.
 *
 *   The multi-parameter interface lm_lmdif is for users who want
 *   full control and flexibility. Most users will be better off using
 *   the simpler interface lmmin provided above.
 *
 *   Parameters:
 *
 *	m is a positive integer input variable set to the number
 *	  of functions.
 *
 *	n is a positive integer input variable set to the number
 *	  of variables; n must not exceed m.
 *
 *	x is an array of length n. On input x must contain an initial
 *        estimate of the solution vector. On OUTPUT x contains the
 *        final estimate of the solution vector.
 *
 *	fvec is an OUTPUT array of length m which contains
 *	  the functions evaluated at the output x.
 *
 *	ftol is a nonnegative input variable. Termination occurs when
 *        both the actual and predicted relative reductions in the sum
 *        of squares are at most ftol. Therefore, ftol measures the
 *        relative error desired in the sum of squares.
 *
 *	xtol is a nonnegative input variable. Termination occurs when
 *        the relative error between two consecutive iterates is at
 *        most xtol. Therefore, xtol measures the relative error desired
 *        in the approximate solution.
 *
 *	gtol is a nonnegative input variable. Termination occurs when
 *        the cosine of the angle between fvec and any column of the
 *        jacobian is at most gtol in absolute value. Therefore, gtol
 *        measures the orthogonality desired between the function vector
 *        and the columns of the jacobian.
 *
 *	maxfev is a positive integer input variable. Termination
 *	  occurs when the number of calls to lm_fcn is at least
 *	  maxfev by the end of an iteration.
 *
 *	epsfcn is an input variable used in choosing a step length for
 *        the forward-difference approximation. The relative errors in
 *        the functions are assumed to be of the order of epsfcn.
 *
 *	diag is an array of length n. If mode = 1 (see below), diag is
 *        internally set. If mode = 2, diag must contain positive entries
 *        that serve as multiplicative scale factors for the variables.
 *
 *	mode is an integer input variable. If mode = 1, the
 *	  variables will be scaled internally. If mode = 2,
 *	  the scaling is specified by the input diag.
 *
 *	factor is a positive input variable used in determining the
 *	  initial step bound. This bound is set to the product of
 *	  factor and the euclidean norm of diag*x if nonzero, or else
 *	  to factor itself. In most cases factor should lie in the
 *	  interval (0.1,100.0). Generally, the value 100.0 is recommended.
 *
 *	info is an integer OUTPUT variable that indicates the termination
 *        status of lm_lmdif as follows:
 *
 *        info < 0  termination requested by user-supplied routine *evaluate;
 *
 *	  info = 0  fnorm almost vanishing;
 *
 *	  info = 1  both actual and predicted relative reductions
 *		    in the sum of squares are at most ftol;
 *
 *	  info = 2  relative error between two consecutive iterates
 *		    is at most xtol;
 *
 *	  info = 3  conditions for info = 1 and info = 2 both hold;
 *
 *	  info = 4  the cosine of the angle between fvec and any
 *		    column of the jacobian is at most gtol in
 *		    absolute value;
 *
 *	  info = 5  number of calls to lm_fcn has reached or
 *		    exceeded maxfev;
 *
 *	  info = 6  ftol is too small: no further reduction in
 *		    the sum of squares is possible;
 *
 *	  info = 7  xtol is too small: no further improvement in
 *		    the approximate solution x is possible;
 *
 *	  info = 8  gtol is too small: fvec is orthogonal to the
 *		    columns of the jacobian to machine precision;
 *
 *	  info =10  improper input parameters;
 *
 *	nfev is an OUTPUT variable set to the number of calls to the
 *        user-supplied routine *evaluate.
 *
 *	fjac is an OUTPUT m by n array. The upper n by n submatrix
 *	  of fjac contains an upper triangular matrix r with
 *	  diagonal elements of nonincreasing magnitude such that
 *
 *		pT*(jacT*jac)*p = rT*r
 *
 *              (NOTE: T stands for matrix transposition),
 *
 *	  where p is a permutation matrix and jac is the final
 *	  calculated jacobian. Column j of p is column ipvt(j)
 *	  (see below) of the identity matrix. The lower trapezoidal
 *	  part of fjac contains information generated during
 *	  the computation of r.
 *
 *	ipvt is an integer OUTPUT array of length n. It defines a
 *        permutation matrix p such that jac*p = q*r, where jac is
 *        the final calculated jacobian, q is orthogonal (not stored),
 *        and r is upper triangular with diagonal elements of
 *        nonincreasing magnitude. Column j of p is column ipvt(j)
 *        of the identity matrix.
 *
 *	qtf is an OUTPUT array of length n which contains
 *	  the first n elements of the vector (q transpose)*fvec.
 *
 *	wa1, wa2, and wa3 are work arrays of length n.
 *
 *	wa4 is a work array of length m, used among others to hold
 *        residuals from evaluate.
 *
 *      evaluate points to the subroutine which calculates the
 *        m nonlinear functions. Implementations should be written as follows:
 *
 *        void evaluate( double* par, int m_dat, void *data,
 *                       double* fvec, int *info )
 *        {
 *           // for ( i=0; i<m_dat; ++i )
 *           //     calculate fvec[i] for given parameters par;
 *           // to stop the minimization, 
 *           //     set *info to a negative integer.
 *        }
 *
 *      printout points to the subroutine which informs about fit progress.
 *        Call with printout=0 if no printout is desired.
 *        Call with printout=lm_printout_std to use the default implementation.
 *
 *      printflags is passed to printout.
 *
 *      data is an input pointer to an arbitrary structure that is passed to
 *        evaluate. Typically, it contains experimental data to be fitted.
 *
 */
    int i, iter, j;
    double actred, delta, dirder, eps, fnorm, fnorm1, gnorm, par, pnorm,
	prered, ratio, step, sum, temp, temp1, temp2, temp3, xnorm;
    static double p1 = 0.1;
    static double p0001 = 1.0e-4;

    *nfev = 0;			/* function evaluation counter */
    iter = 0;			/* outer loop counter */
    par = 0;			/* levenberg-marquardt parameter */
    delta = 0;	 /* to prevent a warning (initialization within if-clause) */
    xnorm = 0;	 /* ditto */
    temp = MAX(epsfcn, LM_MACHEP);
    eps = sqrt(temp); /* for calculating the Jacobian by forward differences */

/*** lmdif: check input parameters for errors. ***/

    if ((n <= 0) || (m < n) || (ftol < 0.)
	|| (xtol < 0.) || (gtol < 0.) || (maxfev <= 0) || (factor <= 0.)) {
	*info = 10;		// invalid parameter
	return;
    }
    if (mode == 2) {		/* scaling by diag[] */
	for (j = 0; j < n; j++) {	/* check for nonpositive elements */
	    if (diag[j] <= 0.0) {
		*info = 10;	// invalid parameter
		return;
	    }
	}
    }
#ifdef LMFIT_DEBUG_MESSAGES
    printf("lmdif\n");
#endif

/*** lmdif: evaluate function at starting point and calculate norm. ***/

    *info = 0;
    (*evaluate) (x, m, data, fvec, info);
    ++(*nfev);
    if( printout )
        (*printout) (n, x, m, data, fvec, printflags, 0, 0, *nfev);
    if (*info < 0)
	return;
    fnorm = lm_enorm(m, fvec);
    if( fnorm <= LM_DWARF ){
        *info = 0;
        return;
    }

/*** lmdif: the outer loop. ***/

    do {
#ifdef LMFIT_DEBUG_MESSAGES
	printf("lmdif/ outer loop iter=%d nfev=%d fnorm=%.10e\n",
	       iter, *nfev, fnorm);
#endif

/*** outer: calculate the Jacobian. ***/

	for (j = 0; j < n; j++) {
	    temp = x[j];
	    step = MAX(eps*eps, eps * fabs(temp));
	    x[j] = temp + step; /* replace temporarily */
	    *info = 0;
	    (*evaluate) (x, m, data, wa4, info);
            ++(*nfev);
            if( printout )
                (*printout) (n, x, m, data, wa4, printflags, 1, iter, *nfev);
	    if (*info < 0)
		return;	/* user requested break */
	    for (i = 0; i < m; i++)
		fjac[j*m+i] = (wa4[i] - fvec[i]) / step;
	    x[j] = temp; /* restore */
	}
#ifdef LMFIT_DEBUG_MATRIX
	/* print the entire matrix */
	for (i = 0; i < m; i++) {
	    for (j = 0; j < n; j++)
		printf("%.5e ", fjac[j*m+i]);
	    printf("\n");
	}
#endif

/*** outer: compute the qr factorization of the Jacobian. ***/

	lm_qrfac(m, n, fjac, 1, ipvt, wa1, wa2, wa3);
        /* return values are ipvt, wa1=rdiag, wa2=acnorm */

	if (!iter) { 
            /* first iteration only */
	    if (mode != 2) {
                /* diag := norms of the columns of the initial Jacobian */
		for (j = 0; j < n; j++) {
		    diag[j] = wa2[j];
		    if (wa2[j] == 0.)
			diag[j] = 1.;
		}
	    }
            /* use diag to scale x, then calculate the norm */
	    for (j = 0; j < n; j++)
		wa3[j] = diag[j] * x[j];
	    xnorm = lm_enorm(n, wa3);
            /* initialize the step bound delta. */
	    delta = factor * xnorm;
	    if (delta == 0.)
		delta = factor;
	} else {
            if (mode != 2) {
                for (j = 0; j < n; j++)
                    diag[j] = MAX( diag[j], wa2[j] );
            }
        }

/*** outer: form (q transpose)*fvec and store first n components in qtf. ***/

	for (i = 0; i < m; i++)
	    wa4[i] = fvec[i];

	for (j = 0; j < n; j++) {
	    temp3 = fjac[j*m+j];
	    if (temp3 != 0.) {
		sum = 0;
		for (i = j; i < m; i++)
		    sum += fjac[j*m+i] * wa4[i];
		temp = -sum / temp3;
		for (i = j; i < m; i++)
		    wa4[i] += fjac[j*m+i] * temp;
	    }
	    fjac[j*m+j] = wa1[j];
	    qtf[j] = wa4[j];
	}

/*** outer: compute norm of scaled gradient and test for convergence. ***/

	gnorm = 0;
        for (j = 0; j < n; j++) {
            if (wa2[ipvt[j]] == 0)
                continue;
            sum = 0.;
            for (i = 0; i <= j; i++)
                sum += fjac[j*m+i] * qtf[i];
            gnorm = MAX( gnorm, fabs( sum / wa2[ipvt[j]] / fnorm ) );
        }

	if (gnorm <= gtol) {
	    *info = 4;
	    return;
	}

/*** the inner loop. ***/
	do {
#ifdef LMFIT_DEBUG_MESSAGES
	    printf("lmdif/ inner loop iter=%d nfev=%d\n", iter, *nfev);
#endif

/*** inner: determine the levenberg-marquardt parameter. ***/

	    lm_lmpar( n, fjac, m, ipvt, diag, qtf, delta, &par,
                      wa1, wa2, wa4, wa3 );
            /* used return values are fjac (partly), par, wa1=x, wa3=diag*x */

	    for (j = 0; j < n; j++)
		wa2[j] = x[j] - wa1[j]; /* new parameter vector ? */

	    pnorm = lm_enorm(n, wa3);

            /* at first call, adjust the initial step bound. */

	    if (*nfev <= 1+n)
		delta = MIN(delta, pnorm);

/*** inner: evaluate the function at x + p and calculate its norm. ***/

	    *info = 0;
	    (*evaluate) (wa2, m, data, wa4, info);
            ++(*nfev);
            if( printout )
                (*printout) (n, wa2, m, data, wa4, printflags, 2, iter, *nfev);
	    if (*info < 0)
		return; /* user requested break. */

	    fnorm1 = lm_enorm(m, wa4);
#ifdef LMFIT_DEBUG_MESSAGES
	    printf("lmdif/ pnorm %.10e  fnorm1 %.10e  fnorm %.10e"
		   " delta=%.10e par=%.10e\n",
		   pnorm, fnorm1, fnorm, delta, par);
#endif

/*** inner: compute the scaled actual reduction. ***/

	    if (p1 * fnorm1 < fnorm)
		actred = 1 - SQR(fnorm1 / fnorm);
	    else
		actred = -1;

/*** inner: compute the scaled predicted reduction and 
     the scaled directional derivative. ***/

	    for (j = 0; j < n; j++) {
		wa3[j] = 0;
		for (i = 0; i <= j; i++)
		    wa3[i] -= fjac[j*m+i] * wa1[ipvt[j]];
	    }
	    temp1 = lm_enorm(n, wa3) / fnorm;
	    temp2 = sqrt(par) * pnorm / fnorm;
	    prered = SQR(temp1) + 2 * SQR(temp2);
	    dirder = -(SQR(temp1) + SQR(temp2));

/*** inner: compute the ratio of the actual to the predicted reduction. ***/

	    ratio = prered != 0 ? actred / prered : 0;
#ifdef LMFIT_DEBUG_MESSAGES
	    printf("lmdif/ actred=%.10e prered=%.10e ratio=%.10e"
		   " sq(1)=%.10e sq(2)=%.10e dd=%.10e\n",
		   actred, prered, prered != 0 ? ratio : 0.,
		   SQR(temp1), SQR(temp2), dirder);
#endif

/*** inner: update the step bound. ***/

	    if (ratio <= 0.25) {
		if (actred >= 0.)
		    temp = 0.5;
		else
		    temp = 0.5 * dirder / (dirder + 0.55 * actred);
		if (p1 * fnorm1 >= fnorm || temp < p1)
		    temp = p1;
		delta = temp * MIN(delta, pnorm / p1);
		par /= temp;
	    } else if (par == 0. || ratio >= 0.75) {
		delta = pnorm / 0.5;
		par *= 0.5;
	    }

/*** inner: test for successful iteration. ***/

	    if (ratio >= p0001) {
                /* yes, success: update x, fvec, and their norms. */
		for (j = 0; j < n; j++) {
		    x[j] = wa2[j];
		    wa2[j] = diag[j] * x[j];
		}
		for (i = 0; i < m; i++)
		    fvec[i] = wa4[i];
		xnorm = lm_enorm(n, wa2);
		fnorm = fnorm1;
		iter++;
	    }
#ifdef LMFIT_DEBUG_MESSAGES
	    else {
		printf("ATTN: iteration considered unsuccessful\n");
	    }
#endif

/*** inner: test for convergence. ***/

            if( fnorm<=LM_DWARF ){
                *info = 0;
                return;
            }

	    *info = 0;
	    if (fabs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1)
		*info = 1;
	    if (delta <= xtol * xnorm)
		*info += 2;
	    if (*info != 0)
		return;

/*** inner: tests for termination and stringent tolerances. ***/

	    if (*nfev >= maxfev){
		*info = 5;
                return;
            }
	    if (fabs(actred) <= LM_MACHEP &&
		prered <= LM_MACHEP && 0.5 * ratio <= 1){
		*info = 6;
                return;
            }
	    if (delta <= LM_MACHEP * xnorm){
		*info = 7;
                return;
            }
	    if (gnorm <= LM_MACHEP){
		*info = 8;
		return;
            }

/*** inner: end of the loop. repeat if iteration unsuccessful. ***/

	} while (ratio < p0001);

/*** outer: end of the loop. ***/

    } while (1);

} /*** lm_lmdif. ***/


/*****************************************************************************/
/*  lm_qrfac (QR factorisation, from lapack)                                 */
/*****************************************************************************/

void lm_qrfac(int m, int n, double *a, int pivot, int *ipvt,
	      double *rdiag, double *acnorm, double *wa)
{
/*
 *     This subroutine uses householder transformations with column
 *     pivoting (optional) to compute a qr factorization of the
 *     m by n matrix a. That is, qrfac determines an orthogonal
 *     matrix q, a permutation matrix p, and an upper trapezoidal
 *     matrix r with diagonal elements of nonincreasing magnitude,
 *     such that a*p = q*r. The householder transformation for
 *     column k, k = 1,2,...,min(m,n), is of the form
 *
 *	    i - (1/u(k))*u*uT
 *
 *     where u has zeroes in the first k-1 positions. The form of
 *     this transformation and the method of pivoting first
 *     appeared in the corresponding linpack subroutine.
 *
 *     Parameters:
 *
 *	m is a positive integer input variable set to the number
 *	  of rows of a.
 *
 *	n is a positive integer input variable set to the number
 *	  of columns of a.
 *
 *	a is an m by n array. On input a contains the matrix for
 *	  which the qr factorization is to be computed. On OUTPUT
 *	  the strict upper trapezoidal part of a contains the strict
 *	  upper trapezoidal part of r, and the lower trapezoidal
 *	  part of a contains a factored form of q (the non-trivial
 *	  elements of the u vectors described above).
 *
 *	pivot is a logical input variable. If pivot is set true,
 *	  then column pivoting is enforced. If pivot is set false,
 *	  then no column pivoting is done.
 *
 *	ipvt is an integer OUTPUT array of length lipvt. This array
 *	  defines the permutation matrix p such that a*p = q*r.
 *	  Column j of p is column ipvt(j) of the identity matrix.
 *	  If pivot is false, ipvt is not referenced.
 *
 *	rdiag is an OUTPUT array of length n which contains the
 *	  diagonal elements of r.
 *
 *	acnorm is an OUTPUT array of length n which contains the
 *	  norms of the corresponding columns of the input matrix a.
 *	  If this information is not needed, then acnorm can coincide
 *	  with rdiag.
 *
 *	wa is a work array of length n. If pivot is false, then wa
 *	  can coincide with rdiag.
 *
 */
    int i, j, k, kmax, minmn;
    double ajnorm, sum, temp;

/*** qrfac: compute initial column norms and initialize several arrays. ***/

    for (j = 0; j < n; j++) {
	acnorm[j] = lm_enorm(m, &a[j*m]);
	rdiag[j] = acnorm[j];
	wa[j] = rdiag[j];
	if (pivot)
	    ipvt[j] = j;
    }
#ifdef LMFIT_DEBUG_MESSAGES
    printf("qrfac\n");
#endif

/*** qrfac: reduce a to r with householder transformations. ***/

    minmn = MIN(m, n);
    for (j = 0; j < minmn; j++) {
	if (!pivot)
	    goto pivot_ok;

        /** bring the column of largest norm into the pivot position. **/

	kmax = j;
	for (k = j + 1; k < n; k++)
	    if (rdiag[k] > rdiag[kmax])
		kmax = k;
	if (kmax == j)
	    goto pivot_ok;

	for (i = 0; i < m; i++) {
	    temp = a[j*m+i];
	    a[j*m+i] = a[kmax*m+i];
	    a[kmax*m+i] = temp;
	}
	rdiag[kmax] = rdiag[j];
	wa[kmax] = wa[j];
	k = ipvt[j];
	ipvt[j] = ipvt[kmax];
	ipvt[kmax] = k;

      pivot_ok:
        /** compute the Householder transformation to reduce the
            j-th column of a to a multiple of the j-th unit vector. **/

	ajnorm = lm_enorm(m-j, &a[j*m+j]);
	if (ajnorm == 0.) {
	    rdiag[j] = 0;
	    continue;
	}

	if (a[j*m+j] < 0.)
	    ajnorm = -ajnorm;
	for (i = j; i < m; i++)
	    a[j*m+i] /= ajnorm;
	a[j*m+j] += 1;

        /** apply the transformation to the remaining columns
            and update the norms. **/

	for (k = j + 1; k < n; k++) {
	    sum = 0;

	    for (i = j; i < m; i++)
		sum += a[j*m+i] * a[k*m+i];

	    temp = sum / a[j + m * j];

	    for (i = j; i < m; i++)
		a[k*m+i] -= temp * a[j*m+i];

	    if (pivot && rdiag[k] != 0.) {
		temp = a[m * k + j] / rdiag[k];
		temp = MAX(0., 1 - temp * temp);
		rdiag[k] *= sqrt(temp);
		temp = rdiag[k] / wa[k];
		if ( 0.05 * SQR(temp) <= LM_MACHEP ) {
		    rdiag[k] = lm_enorm(m-j-1, &a[m*k+j+1]);
		    wa[k] = rdiag[k];
		}
	    }
	}

	rdiag[j] = -ajnorm;
    }
}

/*****************************************************************************/
/*  lm_lmpar (determine Levenberg-Marquardt parameter)                       */
/*****************************************************************************/

void lm_lmpar(int n, double *r, int ldr, int *ipvt, double *diag,
	      double *qtb, double delta, double *par, double *x,
	      double *sdiag, double *aux, double *xdi)
{
/*     Given an m by n matrix a, an n by n nonsingular diagonal
 *     matrix d, an m-vector b, and a positive number delta,
 *     the problem is to determine a value for the parameter
 *     par such that if x solves the system
 *
 *	    a*x = b  and  sqrt(par)*d*x = 0
 *
 *     in the least squares sense, and dxnorm is the euclidean
 *     norm of d*x, then either par=0 and (dxnorm-delta) < 0.1*delta,
 *     or par>0 and abs(dxnorm-delta) < 0.1*delta.
 *
 *     Using lm_qrsolv, this subroutine completes the solution of the problem
 *     if it is provided with the necessary information from the
 *     qr factorization, with column pivoting, of a. That is, if
 *     a*p = q*r, where p is a permutation matrix, q has orthogonal
 *     columns, and r is an upper triangular matrix with diagonal
 *     elements of nonincreasing magnitude, then lmpar expects
 *     the full upper triangle of r, the permutation matrix p,
 *     and the first n components of qT*b. On output
 *     lmpar also provides an upper triangular matrix s such that
 *
 *	    pT*(aT*a + par*d*d)*p = sT*s.
 *
 *     s is employed within lmpar and may be of separate interest.
 *
 *     Only a few iterations are generally needed for convergence
 *     of the algorithm. If, however, the limit of 10 iterations
 *     is reached, then the output par will contain the best
 *     value obtained so far.
 *
 *     parameters:
 *
 *	n is a positive integer input variable set to the order of r.
 *
 *	r is an n by n array. on input the full upper triangle
 *	  must contain the full upper triangle of the matrix r.
 *	  on OUTPUT the full upper triangle is unaltered, and the
 *	  strict lower triangle contains the strict upper triangle
 *	  (transposed) of the upper triangular matrix s.
 *
 *	ldr is a positive integer input variable not less than n
 *	  which specifies the leading dimension of the array r.
 *
 *	ipvt is an integer input array of length n which defines the
 *	  permutation matrix p such that a*p = q*r. column j of p
 *	  is column ipvt(j) of the identity matrix.
 *
 *	diag is an input array of length n which must contain the
 *	  diagonal elements of the matrix d.
 *
 *	qtb is an input array of length n which must contain the first
 *	  n elements of the vector (q transpose)*b.
 *
 *	delta is a positive input variable which specifies an upper
 *	  bound on the euclidean norm of d*x.
 *
 *	par is a nonnegative variable. on input par contains an
 *	  initial estimate of the levenberg-marquardt parameter.
 *	  on OUTPUT par contains the final estimate.
 *
 *	x is an OUTPUT array of length n which contains the least
 *	  squares solution of the system a*x = b, sqrt(par)*d*x = 0,
 *	  for the output par.
 *
 *	sdiag is an array of length n which contains the
 *	  diagonal elements of the upper triangular matrix s.
 *
 *	aux is a multi-purpose work array of length n.
 *
 *	xdi is a work array of length n. On OUTPUT: diag[j] * x[j].
 *
 */
    int i, iter, j, nsing;
    double dxnorm, fp, fp_old, gnorm, parc, parl, paru;
    double sum, temp;
    static double p1 = 0.1;

#ifdef LMFIT_DEBUG_MESSAGES
    printf("lmpar\n");
#endif

/*** lmpar: compute and store in x the gauss-newton direction. if the
     jacobian is rank-deficient, obtain a least squares solution. ***/

    nsing = n;
    for (j = 0; j < n; j++) {
	aux[j] = qtb[j];
	if (r[j * ldr + j] == 0 && nsing == n)
	    nsing = j;
	if (nsing < n)
	    aux[j] = 0;
    }
#ifdef LMFIT_DEBUG_MESSAGES
    printf("nsing %d ", nsing);
#endif
    for (j = nsing - 1; j >= 0; j--) {
	aux[j] = aux[j] / r[j + ldr * j];
	temp = aux[j];
	for (i = 0; i < j; i++)
	    aux[i] -= r[j * ldr + i] * temp;
    }

    for (j = 0; j < n; j++)
	x[ipvt[j]] = aux[j];

/*** lmpar: initialize the iteration counter, evaluate the function at the
     origin, and test for acceptance of the gauss-newton direction. ***/

    iter = 0;
    for (j = 0; j < n; j++)
	xdi[j] = diag[j] * x[j];
    dxnorm = lm_enorm(n, xdi);
    fp = dxnorm - delta;
    if (fp <= p1 * delta) {
#ifdef LMFIT_DEBUG_MESSAGES
	printf("lmpar/ terminate (fp<p1*delta)\n");
#endif
	*par = 0;
	return;
    }

/*** lmpar: if the jacobian is not rank deficient, the newton
     step provides a lower bound, parl, for the 0. of
     the function. otherwise set this bound to 0.. ***/

    parl = 0;
    if (nsing >= n) {
	for (j = 0; j < n; j++)
	    aux[j] = diag[ipvt[j]] * xdi[ipvt[j]] / dxnorm;

	for (j = 0; j < n; j++) {
	    sum = 0.;
	    for (i = 0; i < j; i++)
		sum += r[j * ldr + i] * aux[i];
	    aux[j] = (aux[j] - sum) / r[j + ldr * j];
	}
	temp = lm_enorm(n, aux);
	parl = fp / delta / temp / temp;
    }

/*** lmpar: calculate an upper bound, paru, for the 0. of the function. ***/

    for (j = 0; j < n; j++) {
	sum = 0;
	for (i = 0; i <= j; i++)
	    sum += r[j * ldr + i] * qtb[i];
	aux[j] = sum / diag[ipvt[j]];
    }
    gnorm = lm_enorm(n, aux);
    paru = gnorm / delta;
    if (paru == 0.)
	paru = LM_DWARF / MIN(delta, p1);

/*** lmpar: if the input par lies outside of the interval (parl,paru),
     set par to the closer endpoint. ***/

    *par = MAX(*par, parl);
    *par = MIN(*par, paru);
    if (*par == 0.)
	*par = gnorm / dxnorm;
#ifdef LMFIT_DEBUG_MESSAGES
    printf("lmpar/ parl %.4e  par %.4e  paru %.4e\n", parl, *par, paru);
#endif

/*** lmpar: iterate. ***/

    for (;; iter++) {

        /** evaluate the function at the current value of par. **/

	if (*par == 0.)
	    *par = MAX(LM_DWARF, 0.001 * paru);
	temp = sqrt(*par);
	for (j = 0; j < n; j++)
	    aux[j] = temp * diag[j];

	lm_qrsolv( n, r, ldr, ipvt, aux, qtb, x, sdiag, xdi );
        /* return values are r, x, sdiag */

	for (j = 0; j < n; j++)
	    xdi[j] = diag[j] * x[j]; /* used as output */
	dxnorm = lm_enorm(n, xdi);
	fp_old = fp;
	fp = dxnorm - delta;
        
        /** if the function is small enough, accept the current value
            of par. Also test for the exceptional cases where parl
            is zero or the number of iterations has reached 10. **/

	if (fabs(fp) <= p1 * delta
	    || (parl == 0. && fp <= fp_old && fp_old < 0.)
	    || iter == 10)
	    break; /* the only exit from the iteration. */
        
        /** compute the Newton correction. **/

	for (j = 0; j < n; j++)
	    aux[j] = diag[ipvt[j]] * xdi[ipvt[j]] / dxnorm;

	for (j = 0; j < n; j++) {
	    aux[j] = aux[j] / sdiag[j];
	    for (i = j + 1; i < n; i++)
		aux[i] -= r[j * ldr + i] * aux[j];
	}
	temp = lm_enorm(n, aux);
	parc = fp / delta / temp / temp;

        /** depending on the sign of the function, update parl or paru. **/

	if (fp > 0)
	    parl = MAX(parl, *par);
	else if (fp < 0)
	    paru = MIN(paru, *par);
	/* the case fp==0 is precluded by the break condition  */
        
        /** compute an improved estimate for par. **/
        
	*par = MAX(parl, *par + parc);
        
    }

} /*** lm_lmpar. ***/

/*****************************************************************************/
/*  lm_qrsolv (linear least-squares)                                         */
/*****************************************************************************/

void lm_qrsolv(int n, double *r, int ldr, int *ipvt, double *diag,
	       double *qtb, double *x, double *sdiag, double *wa)
{
/*
 *     Given an m by n matrix a, an n by n diagonal matrix d,
 *     and an m-vector b, the problem is to determine an x which
 *     solves the system
 *
 *	    a*x = b  and  d*x = 0
 *
 *     in the least squares sense.
 *
 *     This subroutine completes the solution of the problem
 *     if it is provided with the necessary information from the
 *     qr factorization, with column pivoting, of a. That is, if
 *     a*p = q*r, where p is a permutation matrix, q has orthogonal
 *     columns, and r is an upper triangular matrix with diagonal
 *     elements of nonincreasing magnitude, then qrsolv expects
 *     the full upper triangle of r, the permutation matrix p,
 *     and the first n components of (q transpose)*b. The system
 *     a*x = b, d*x = 0, is then equivalent to
 *
 *	    r*z = qT*b,  pT*d*p*z = 0,
 *
 *     where x = p*z. If this system does not have full rank,
 *     then a least squares solution is obtained. On output qrsolv
 *     also provides an upper triangular matrix s such that
 *
 *	    pT *(aT *a + d*d)*p = sT *s.
 *
 *     s is computed within qrsolv and may be of separate interest.
 *
 *     Parameters
 *
 *	n is a positive integer input variable set to the order of r.
 *
 *	r is an n by n array. On input the full upper triangle
 *	  must contain the full upper triangle of the matrix r.
 *	  On OUTPUT the full upper triangle is unaltered, and the
 *	  strict lower triangle contains the strict upper triangle
 *	  (transposed) of the upper triangular matrix s.
 *
 *	ldr is a positive integer input variable not less than n
 *	  which specifies the leading dimension of the array r.
 *
 *	ipvt is an integer input array of length n which defines the
 *	  permutation matrix p such that a*p = q*r. Column j of p
 *	  is column ipvt(j) of the identity matrix.
 *
 *	diag is an input array of length n which must contain the
 *	  diagonal elements of the matrix d.
 *
 *	qtb is an input array of length n which must contain the first
 *	  n elements of the vector (q transpose)*b.
 *
 *	x is an OUTPUT array of length n which contains the least
 *	  squares solution of the system a*x = b, d*x = 0.
 *
 *	sdiag is an OUTPUT array of length n which contains the
 *	  diagonal elements of the upper triangular matrix s.
 *
 *	wa is a work array of length n.
 *
 */
    int i, kk, j, k, nsing;
    double qtbpj, sum, temp;
    double _sin, _cos, _tan, _cot; /* local variables, not functions */

/*** qrsolv: copy r and (q transpose)*b to preserve input and initialize s.
     in particular, save the diagonal elements of r in x. ***/

    for (j = 0; j < n; j++) {
	for (i = j; i < n; i++)
	    r[j * ldr + i] = r[i * ldr + j];
	x[j] = r[j * ldr + j];
	wa[j] = qtb[j];
    }
#ifdef LMFIT_DEBUG_MESSAGES
    printf("qrsolv\n");
#endif

/*** qrsolv: eliminate the diagonal matrix d using a Givens rotation. ***/

    for (j = 0; j < n; j++) {

/*** qrsolv: prepare the row of d to be eliminated, locating the
     diagonal element using p from the qr factorization. ***/

	if (diag[ipvt[j]] == 0.)
	    goto L90;
	for (k = j; k < n; k++)
	    sdiag[k] = 0.;
	sdiag[j] = diag[ipvt[j]];

/*** qrsolv: the transformations to eliminate the row of d modify only 
     a single element of qT*b beyond the first n, which is initially 0. ***/

	qtbpj = 0.;
	for (k = j; k < n; k++) {

            /** determine a Givens rotation which eliminates the
                appropriate element in the current row of d. **/

	    if (sdiag[k] == 0.)
		continue;
	    kk = k + ldr * k;
	    if (fabs(r[kk]) < fabs(sdiag[k])) {
		_cot = r[kk] / sdiag[k];
		_sin = 1 / sqrt(1 + SQR(_cot));
		_cos = _sin * _cot;
	    } else {
		_tan = sdiag[k] / r[kk];
		_cos = 1 / sqrt(1 + SQR(_tan));
		_sin = _cos * _tan;
	    }

            /** compute the modified diagonal element of r and
                the modified element of ((q transpose)*b,0). **/

	    r[kk] = _cos * r[kk] + _sin * sdiag[k];
	    temp = _cos * wa[k] + _sin * qtbpj;
	    qtbpj = -_sin * wa[k] + _cos * qtbpj;
	    wa[k] = temp;

            /** accumulate the tranformation in the row of s. **/

	    for (i = k + 1; i < n; i++) {
		temp = _cos * r[k * ldr + i] + _sin * sdiag[i];
		sdiag[i] = -_sin * r[k * ldr + i] + _cos * sdiag[i];
		r[k * ldr + i] = temp;
	    }
	}

      L90:
        /** store the diagonal element of s and restore
            the corresponding diagonal element of r. **/

	sdiag[j] = r[j * ldr + j];
	r[j * ldr + j] = x[j];
    }

/*** qrsolv: solve the triangular system for z. if the system is
     singular, then obtain a least squares solution. ***/

    nsing = n;
    for (j = 0; j < n; j++) {
	if (sdiag[j] == 0. && nsing == n)
	    nsing = j;
	if (nsing < n)
	    wa[j] = 0;
    }

    for (j = nsing - 1; j >= 0; j--) {
	sum = 0;
	for (i = j + 1; i < nsing; i++)
	    sum += r[j * ldr + i] * wa[i];
	wa[j] = (wa[j] - sum) / sdiag[j];
    }

/*** qrsolv: permute the components of z back to components of x. ***/

    for (j = 0; j < n; j++)
	x[ipvt[j]] = wa[j];

} /*** lm_qrsolv. ***/

double nonlinearFunc( double x, const double *par,int obsNum )
{
  int i;
  double dt;
  double fitfunc;
  int nfit=0;
  int f0dv,tdv;

  fitfunc = 0.0;

  if (global_fitf0==1)
    fitfunc+=par[nfit++];
  else
    fitfunc+=global_valf0;

  if (global_fitf1==1)
    fitfunc+=x*par[nfit++];
  else
    fitfunc+=x*global_valf1*86400.0;

  for (i=0;i<global_nglt;i++)
    {
      dt = x - global_glitch[i].glep;
      if (dt > 0)
	{
	  if (global_glitch[i].fitf0==1) fitfunc += par[nfit++]; 
	  if (global_glitch[i].fitf1==1) fitfunc += par[nfit++]*dt;
	  if (global_glitch[i].fitf0d==1) {f0dv = nfit; fitfunc += par[nfit++]*exp(-dt/global_glitch[0].gltd);}
	  if (global_glitch[i].fittd==1) {tdv = nfit; fitfunc += par[f0dv]*exp(-dt/par[tdv]); nfit++;}
	}
    }
  return fitfunc;
}
char * plugVersionCheck = (char *)TEMPO2_h_VER;
