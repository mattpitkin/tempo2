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

/* This plugin is a simple graphical interface for TEMPO2          */
/* It demonstrates methods to interface with the main TEMPO2 code  */
/* and provides a very quick and easy to use interface that allows */
/* pulsar timing residuals to be plotted.  The user can delete     */ 
/* points and redo the fit.                                        */
/* The plotting commands are the same as in the 'plk' package      */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"
#include "T2toolkit.h"
#include <time.h>
#include <cpgplot.h>

using namespace std;   /* Is this required for a plugin ? Yes, for linux */

char dcmFile[MAX_FILELEN];
char covarFuncFile[MAX_FILELEN];

void overPlotN(int overN,float overX[], float overY[],float overYe[]);
void doPlot(pulsar *psr,int npsr,char *gr,double unitFlag,char parFile[][MAX_FILELEN],
	    char timFile[][MAX_FILELEN],float lockx1,float lockx2,float locky1,float locky2,int xplot,int yplot,int publish,int argc,char *argv[],int menu,char *setupFile,
            int showChisq,int nohead,char* flagColour,char *bandsFile,int displayPP);
int setPlot(float *x,int count,pulsar *psr,int iobs,double unitFlag,int plotPhase,int plot,int *userValChange,
	    char *userCMD,char *userValStr,float *userX,longdouble centreEpoch,int log,char *flagStr);
void drawAxisSel(float x,float y,char *str,int sel1,int sel2);
float findMinY(float *y,float *x,int count,float xmin,float xmax);
float findMaxY(float *y,float *x,int count,float xmin,float xmax);
float findMean(float *x,pulsar *psr,int i1,int i2);
double findMeanD(float *x,pulsar *psr,int i1,int i2);
void callFit(pulsar *psr,int npsr);
void reFit(int fitFlag,int setZoomX1,int setZoomX2,float zoomX1,float zoomX2,longdouble origStart,
	   longdouble origFinish,longdouble centreEpoch,pulsar *psr,int npsr,int plotX,char *dcmFile,char *covarFuncFile,int zoom);
float deletePoint(pulsar *psr,float *x,float *y,int *id,int count,float mouseX,float mouseY);
void displayStatistics(float *x,float *y,int count,float plotx1,float plotx2,float ploty1,float ploty2);
int idPoint(pulsar *psr,float *x,float *y,int *id,int count,float mouseX,float mouseY);
double fortranMod(double a,double p);
void sort(float *x,float *y,float *yerr1,float *yerr2,float *freq,int *id,int count);
void changeParameters(pulsar *psr);
void changeFitParameters(pulsar *psr);
void averagePts(float *x,float *y,int n,int width,float *meanX,float *meanY,int *nMean);
void overPlotShapiro(pulsar *psr,float offset,longdouble centreEpoch);
void binResiduals(pulsar *psr,int npsr,float *x,float *y,int count,int *id, int *overN,
		  float overX[], float overY[], float overYe[],int xplot,int yplot,
		  float errBar[],double unitFlag,int plotPhase,double centreEpoch);
void drawMenu(pulsar *psr,float plotx1,float plotx2,float ploty1,float ploty2,int menu,int paramOffset);
void drawMenu3(pulsar *psr, float plotx1,float plotx2,float ploty1,float ploty2,int menu,int xplot,int yplot);
void slaClyd ( int iy, int im, int id, int *ny, int *nd, int *jstat );
void slaCalyd ( int iy, int im, int id, int *ny, int *nd, int *j );
void drawMenu3_2(pulsar *psr, float plotx1,float plotx2,float ploty1,float ploty2,int menu,int xplot,int yplot,int jumpOffset,int iFlagColour, int nFlags);
void checkMenu(pulsar *psr,float mx,float my,int button,int fitFlag,int setZoomX1,int setZoomX2,
	       float zoomX1,float zoomX2,longdouble origStart,longdouble origFinish,longdouble centreEpoch,
	       int menu,int plotx,char parFile[][MAX_FILELEN], char timFile[][MAX_FILELEN],int argc,char *argv[],int *xplot,int *yplot,int *graphics, char highlightID[100][100],char  highlightVal[100][100],int *highlightNum,float aspect,int fontType,int lineWidth,char *bkgrdColour,char *lineColour,int *jumpOffset,int zoom,int *paramOffset);
void checkMenu3(pulsar *psr,float mx,float my,int button,int fitFlag,int setZoomX1,int setZoomX2,
	       float zoomX1,float zoomX2,longdouble origStart,longdouble origFinish,longdouble centreEpoch,
	       int menu,int plotx,char parFile[][MAX_FILELEN], char timFile[][MAX_FILELEN],int argc,char *argv[],int *xplot,int *yplot,int *graphics,char highlightID[100][100],char  highlightVal[100][100],int *highlightNum,float aspect,int fontType,int lineWidth,char *bkgrdColour,char *lineColour,int *jumpOffset);
void setLabel(char *ystr,int yplot,int plotPhase,double unitFlag,longdouble centreEpoch,char *userValStr,char *flagStr);
void drawOption(float x,float y,char *str,int fit);
void swapFit(pulsar *psr,int par,int k,int button);
void newTim(pulsar *psr);
void plotFITWAVES_spec();
void viewModels(pulsar *psr,float x1,float x2,longdouble centreEpoch,int removeMean,double mean,int count,int *id,
		int fitFlag,float *x,float *y);
double lmst2(double mjd,double olong,double *tsid,double *tsid_der);

bool cholmode=false;

/* GLOBAL VARIABLES FOR FITWAVES */
double FITWAVES_omega;
int    FITWAVES_n;
int    FITWAVES_harmonicStep;
double FITWAVES_par[1000];
char flagStore[100][100];

void help() /* Display help */
{
  printf("\nFitting and Calculating Options\n"); /* Playing around with individual TOA's and global TOA trends */
  printf("===============================\n");
  printf("b          Bin TOAs within certain time bin\n");
  printf("c          Change fitting parameters\n");
  printf("d (or right mouse) delete point\n");
  printf("ctrl-d     delete highlighted points\n");
  printf("e          multiply all TOA errors by given amount\n");
  printf("F          run FITWAVES\n"); 
  printf("ctrl-f     remove FITWAVES curve from residuals\n");  
  printf("i (or left mouse) identify point\n");  
  printf("M          toggle removing mean from the residuals\n");
  printf("ctrl-n     Add white noise to site-arrival-times\n");
  printf("p          Change model parameter values\n");
  printf("ctrl-p     Toggle plotting versus pulse phase\n");
  printf("r          Reset (reload .par and .tim file)\n");
  printf("ctrl-r     Select regions in MJDs and write to file\n");
  printf("w          toggle fitting using weights\n");
  printf("x          redo fit using post-fit parameters\n");
  printf("+          add positive phase jump\n");
  printf("-          add negative phase jump\n");
  printf("BACKSPACE  remove all phase jumps\n");
  printf("ctrl-=     add period to residuals above cursor\n");
  printf("/          re-read .par file\n");

  printf("\nPlot Selection\n"); /* Determines WHAT (X and Y axes) will be displayed */
  printf("==============\n");
  printf("D (or middle mouse) view profile\n");
  printf("s          start of zoom section\n");
  printf("f          finish of zoom section\n"); 
  printf("Ctrl-u     Overplot Shapiro delay\n");
  printf("u          unzoom\n");
  printf("v          view profiles for highlighted points\n");
  printf("V          define the user parameter\n");
  printf("Ctrl-v     for pre-fit plotting, decompose the timing model fits\n");
  printf("           (i.e. overplot the fitted curves - only for prefit plots\n");
  printf("ctrl-X     select x-axis specifically\n");
  printf("y          Rescale y-axis only\n");
  printf("Y          set y-scale exactly\n");
  printf("ctrl-Y     select y-axis specifically\n");
  printf("z          Zoom using mouse\n");
  printf("<          in zoom mode include previous observation\n");  
  printf(">          in zoom mode include next observation\n");	     
  printf("1          plot pre-fit  residuals vs date\n");	     
  printf("2          plot post-fit residuals vs date\n");	     
  printf("3          plot pre-fit  residuals vs orbital phase\n");   
  printf("4          plot post-fit residuals vs orbital phase\n");   
  printf("5          plot pre-fit  residuals serially\n");	     
  printf("6          plot post-fit residuals serially\n");	     
  printf("7          plot pre-fit  residuals vs day of year\n");     
  printf("8          plot post-fit residuals vs day of year\n");     
  printf("9          plot pre-fit  residuals vs frequency\n");	     
  printf("a          plot post-fit residuals vs frequency\n");	     
  printf("!          plot pre-fit  residuals vs TOA error\n");	     
  printf("@          plot post-fit residuals vs TOA error\n");       
  printf("#          plot pre-fit  residuals vs user values\n");     
  printf("$          plot post-fit residuals vs user values\n");     
  printf("\%%          plot pre-fit  residuals vs year\n");	     
  printf("^          plot post-fit residuals vs year\n");  	     
  printf("&          plot pre-fit residuals vs elevation\n");	     
  printf("*          plot post-fit residuals vs elevation\n");       
  printf("(          plot pre-fit residuals vs rounded MJD\n");
  printf(")          plot post-fit residuals vs rounded MJD\n");
  printf("\n");
  printf("Options for selecting x and y axes individually\n");
  printf("Ctrl-X n   set x-axis\n");
  printf("Ctrl-Y n   set y-axis\n");
  printf("where n = \n\n");
  printf("1         plot pre-fit residuals\n");
  printf("2         plot post-fit residuals\n");
  printf("3         plot centred MJD\n");
  printf("4         plot orbital phase\n");
  printf("5         plot TOA number\n");
  printf("6         plot day of year\n");
  printf("7         plot frequency\n");
  printf("8         plot TOA error\n");
  printf("9         plot user value\n");
  printf("0         plot year\n");
  printf("-         plot elevation\n");
  printf("K         subtract/add TN red noise\n");
  printf("F         subtract/add TN DM Variations\n");

  printf("\nDisplay Options\n"); /* Determines HOW to display the things */
  printf("===============\n");
  printf("B          place periodic marks on the x-scale\n");
  printf("ctrl-c     Toggle between period epoch and centre for the reference epoch\n"); 
  printf("E          toggle plotting error bars\n");
  printf("g          change graphics device\n");
  printf("G          change gridding on graphics device\n");
  printf("ctrl-e     highlight points more than 3 sigma from the mean\n");
  printf("H          highlight points with specific flag using symbols\n");
  printf("ctrl-i     highlight points with specific flag using colours\n");
  printf("I          indicate individual observations\n");
  printf("j          draw line between points \n");
  printf("J          toggle plotting points\n");
  printf("L          add label to plot\n");
  printf("ctrl-l     add line to plot\n");
  printf("ctrl-m     toggle menu bar\n");
  printf("N          highlight point with a given filename\n");
  printf("o          obtain/highlight all points currently in plot\n");
  printf("ctrl-T     set text size\n");
  printf("U          unhighlight selected points\n");
  printf("[          toggle plotting x-axis on log scale\n");	     
  printf("]          toggle plotting y-axis on log scale\n");	     
  
  printf("\nOutput Options\n");
  printf("==============\n");
  printf("Ctrl-J     output listing of residuals in Jodrell format\n");
  printf("Ctrl-O     output listing of residuals in simple format\n");
  printf("l          list all data points in zoomed region\n");
  printf("m          measure distance between two points\n");
  printf("P          write new .par file\n");
  printf("Ctrl-w     over-write input .par file\n");
  printf("S          save a new .tim file\n");
  printf("Ctrl-S     overwrite input.tim file\n");
  printf("t          Toggle displaying statistics for zoomed region\n");
  printf("Ctrl-z     Listing of all highlighted points\n");

  printf("\nVarious Options\n");
  printf("===============\n");
  printf("C          run unix command with filenames for highlighted observations\n");
  printf("h          this help file\n");
  printf("q          quit\n");
}

/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  int i,gotTim=0;
  int newpar = 0;
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  char newParFile[MAX_FILELEN];
  char gr[100]="/xs";
  double unitFlag=1.0;  /* plot in seconds */
  float lockx1=0,lockx2=0;
  float locky1=0,locky2=0;
  int   xplot=3;
  int   yplot=1;
  int   publish=0;
  int   menu=3;
  int   nohead=0;
  char  setupFile[100]="";
  char  bandsFile[100]="";
  int   displayPP = 1;
  //display chisq?
  int showChisq = 0;
  char flagColour[100];
  const char *CVS_verNum = "$Revision: 1.61 $";

  if (displayCVSversion == 1) CVSdisplayVersion("plk_plug.C","plugin",CVS_verNum);

  strcpy(flagColour,"");

  *npsr = 1;  /* This graphical interface will only show results for one pulsar */

  strcpy(dcmFile,"NULL");
  strcpy(covarFuncFile,"NULL");

  printf("Graphical Interface: plk emulator\n");
  printf("Authors:             George Hobbs, J. Verbiest (v4. 3 Aug 2007)\n");
  printf("CVS Version:         $Revision: 1.61 $\n");
  printf(" --- type 'h' for help information\n");
  /* Obtain the .par and the .tim file from the command line */

  if (argc==4) /* Only provided .tim name */
    {
      strcpy(timFile[0],argv[3]);
      strcpy(parFile[0],argv[3]);
      parFile[0][strlen(parFile[0])-3] = '\0';
      strcat(parFile[0],"par");
    }
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-grdev")==0)
        strcpy(gr,argv[++i]);
      else if (strcmp(argv[i],"-menu")==0)
        sscanf(argv[i+1],"%d",&menu);
      else if (strcmp(argv[i],"-nophase")==0)
	displayPP=0;
      else if (strcmp(argv[i],"-locky")==0)
        {
          sscanf(argv[++i],"%f",&locky1);
          sscanf(argv[++i],"%f",&locky2);
        }
      else if (strcmp(argv[i],"-lockx")==0)
        {
          sscanf(argv[++i],"%f",&lockx1);
          sscanf(argv[++i],"%f",&lockx2);
        }
      else if (strcmp(argv[i],"-setup")==0)
        strcpy(setupFile,argv[++i]);
      else if (strcmp(argv[i],"-xplot")==0)
        sscanf(argv[++i],"%d",&xplot);
      else if (strcmp(argv[i],"-yplot")==0)
        sscanf(argv[++i],"%d",&yplot);
      else if (strcmp(argv[i],"-publish")==0)
        publish=1;
      else if (strcmp(argv[i],"-image")==0)
        publish=2;
      else if (strcmp(argv[i],"-ms")==0)
        unitFlag=1.0e-3;
      else if (strcmp(argv[i],"-us")==0)
        unitFlag=1.0e-6;
      else if (strcmp(argv[i],"-ns")==0)
        unitFlag=1.0e-9;
      else if (strcmp(argv[i],"-min")==0)
        unitFlag=60.0;
      else if (strcmp(argv[i],"-bands")==0)
	strcpy(bandsFile,argv[++i]);
      else if (strcmp(argv[i],"-nohead")==0)
        nohead=1;
      else if (strcmp(argv[i],"-dcm")==0)
        strcpy(dcmFile,argv[++i]);
      else if (strcmp(argv[i],"-dcf")==0){
		cholmode=true;
        strcpy(covarFuncFile,argv[++i]);
	  }
      else if (strcmp(argv[i],"-newparS")==0)//this is just for use with calcDMe
        {
          newpar = 1;
          strcpy(newParFile,argv[++i]);
        }
      else if (strcmp(argv[i],"-period")==0)
        unitFlag=-1; /* In units of pulse period */
      else if (strcmp(argv[i],"-f")==0)
        {
          strcpy(parFile[0],argv[++i]); 
          if (argv[i+1][0]!='-')
            {
              strcpy(timFile[0],argv[++i]);
              gotTim=1;
            }
        }
      else if (strcmp(argv[i],"-showchisq") == 0)
        {
          showChisq = 1;
        }
      else if (strcmp(argv[i],"-colour") == 0)
        {
          strcpy(flagColour,argv[++i]);
        }
      else if (strcmp(argv[i],"-h")==0||strcmp(argv[i],"--help")==0){
        printf("\n TEMPO2 plk plugin\n");
        printf("===================\n");
        printf("\nUSAGE: \n\t tempo2 -gr plk -f par.par tim.tim\n");
        printf("\n Command line options:\n");
        printf("\t -grdev, followed by the graphics device of choise (e.g. 3/xs; TOAs.ps/cps)\n");
        printf("\t -locky -1e-6 1e-6: lock the y axis to this range\n");
        printf("\t -showchisq: show the chisq of the fit\n");
        printf("\t -xplot 3 or -yplot 5: Determine the x and y axes\n");
        printf("\t\t The arguments for this command are as follows:\n");
        printf("\t\t  1\t prefit residuals\n");
        printf("\t\t  2\t postfit residuals\n");
        printf("\t\t  3\t Modified Julian Date\n");
        printf("\t\t  4\t Orbital phase\n");
        printf("\t\t  5\t TOA number\n");
        printf("\t\t  6\t Day of year\n");
        printf("\t\t  7\t Frequency\n");
        printf("\t\t  8\t TOA error\n");
        printf("\t\t  9\t DPA (derived parallactic angle)\n");
        printf("\t\t 10\t Year\n");
        printf("\t\t 11\t Elevation\n");
        printf("\t\t 12\t Round MJD values\n\n");
        printf("\t -publish: publication style graphics (larger fontsize, no menu,...)\n");
        printf("\t -image: graphics suitable for presentations (larger fontsize, no menu,...)\n");
        printf("\t -ms or -us or -ns: specifies residual units (default is seconds)\n");
        printf("\t -period: calculates residuals in pulse periods, instead of seconds.\n");
        printf("\t -h or --help: this help. More help is available by pressing 'h' while running plk.\n");
        printf("\t -epoch centre : automatically puts epochs such as:\n");
        printf("\t\t PEPOCH (period determination epoch), \n");
        printf("\t\t POSEPOCH (position determination epoch) and \n");
        printf("\t\t DMEPOCH (DM determination epoch) \n");
        printf("\t                 equal to the central MJD of the dataset.\n");
        printf("\n\n===============================================");
        printf("===============================================\n");
        exit(0);
      }
      else if (i==argc-1 && gotTim == 0)
        {
          strcpy(timFile[0],argv[i]);
          gotTim=1;
        }
    }
  if (debugFlag==1) printf("plk: calling readParfile\n");
  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
    
  if (debugFlag==1) printf("plk: calling readTimfile\n");
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  if (psr[0].nobs==0)
    {
      printf("Error: no arrival times loaded\n");
      exit(1);
    }
  if (unitFlag==-1) unitFlag = 1.0/psr[0].param[param_f].val[0]; /* Units of pulse period */
  if (debugFlag==1) printf("plk: calling preProcess %d\n",psr[0].nobs);
  
  preProcess(psr,*npsr,argc,argv);
  if (debugFlag==1) printf("plk: calling callFit %d\n",psr[0].nobs);
  callFit(psr,*npsr);             /* Do all the fitting routines */
  if (newpar==1)
    textOutput(psr,*npsr,0,0,0,1,newParFile);
  if (debugFlag==1) printf("plk: calling doPlot\n");
  doPlot(psr,*npsr,gr,unitFlag,parFile,timFile,lockx1,lockx2,locky1,locky2,xplot,yplot,
	 publish,argc,argv,menu,setupFile,showChisq,nohead,flagColour,bandsFile,displayPP);  /* Do plot */
  if (debugFlag==1) printf("plk: End\n");
  return 0;
}  

/* This function calls all of the fitting routines.             */
/* The function is looped twice, the first time for the pre-fit */
/* residuals and the second time for the post-fit residuals     */

void callFit(pulsar *psr,int npsr)
{
  int iteration,i,p,it,k;
  double globalParameter = 0.0;
  FILE *pin;

  for (it=0;it<psr[0].nits;it++)
    {
      if (it>0) /* Copy post-fit values to pre-fit values */
	{
	  for (i=0;i<MAX_PARAMS;i++)
	    {
	      for (p=0;p<npsr;p++)
		{
		  
		  for (k=0;k<psr[p].param[i].aSize;k++)
		    {
		      psr[p].param[i].prefit[k] = psr[p].param[i].val[k];
		      psr[p].param[i].prefitErr[k] = psr[p].param[i].err[k];
		    }
		}
	    }
	}
      for (iteration=0;iteration<2;iteration++)
        {
          formBatsAll(psr,npsr);
          
          /* Form residuals */
          formResiduals(psr,npsr,1); // iteration);
          
          /* Do the fitting */
          if (iteration==0) 
	    doFitAll(psr,npsr,covarFuncFile);
          else 
	    textOutput(psr,npsr,globalParameter,0,0,0,"");
        }
    }
  if (psr[0].nFit==0)
    {
      printf("No timing residuals to plot.  Please check your TOA file and filter commands\n");
      exit(1);
    }
}


void doPlot(pulsar *psr,int npsr,char *gr,double unitFlag, char parFile[][MAX_FILELEN],
	    char timFile[][MAX_FILELEN],float lockx1, float lockx2, float locky1, float locky2,int xplot,int yplot,
	    int publish,int argc,char *argv[],int menu,char *setupFile, int showChisq,int nohead,char* flagColour,char *bandsFile,int displayPP)
{
  int i,fitFlag=1,exitFlag=0,scale1=0,scale2=psr[0].nobs,count,ncount,j,k;
  longdouble centreEpoch;
  char xstr[1000],ystr[1000],title[1000];
  float lx1[100],lx2[100],ly1[100],ly2[100];
  int overN;
  int   nline=0;
  int id[MAX_OBSN];
  float minx,maxx,miny,maxy,plotx1,plotx2,ploty1,ploty2,mean;
  float mouseX,mouseY,mouseX2,mouseY2;
  char bkgrdColour[100],lineColour[100];
  float zoomX1=0.0,zoomX2=0.0,zoomY1=0.0,zoomY2=0.0;
  float delX1=0.0, delX2=0.0, delY1=0.0,delY2=0.0;
  float aspect=1.0;
  float viewport_x0 = 0.1,viewport_x1 = 0.95,viewport_y0 = 0.15, viewport_y1 = 0.85;
  int   fontType=1;
  int   lineWidth=1;
  char  userValStr[1000]="DPA";
  int   userValChange=1;
  int   setZoomX1=0,setZoomX2=0,setZoomY1=0,setZoomY2=0,graphics=0;
  int   zoom=-1;
  int   noreplot=0;
  int overPlotn = 0;
  int   overPlotS = -1;
  int   plotPhase=-1;
  int   centre=1;
  int   indicateObs=0;
  int   logx=-1,logy=-1;
  int   viewModel=-1;
  int   join=-1;
  int   plotPoints=1;
  int   placeMarks=-1;
  int   plotErr=1;
  int   iFlagColour=0;
  char  key;
  float plottingx=-1e10;
  int   plotFITWAVES=0;
  int   removeFITWAVES=-1;
  int   highlightNum=0;
  char  highlightID[100][100];
  char  highlightVal[100][100];
  float fontSize=1.0;
  float mark0,markT;
  int   label=0;
  char  labelStr[100];
  char  fitType[100];
  char  userCMD[100]="vap";
  float labelX,labelY,shapiroOffset;
  longdouble min1,max1;
  longdouble origStart=-1,origFinish=-1;
  double dmean;
  int   errtype=0;
  int   removeMean=1;
  int   statistics=-1;
  int   bad=0;
  int   bw=-1;
  int   out;
  int jumpOffset=0;
  int okay;

  int  freqColourNum=0;
  float minFreqCol[100], maxFreqCol[100];
  int   freqCol[100],freqStyle[100];

  int  flagColourNum=0;
  char flagIDstr[100][100],flagIDval[100][100];
  int  flagCol[100],flagStyle[100];

  float *userX,*errBar;
  userX = (float *)malloc(sizeof(float)*MAX_OBSN);
  errBar = (float *)malloc(sizeof(float)*MAX_OBSN);
  float *px,*py,*x,*y,*yerr1,*yerr2,*x2,*y2,*yerr1_2,*yerr2_2,*freq;
  float *overX,*overY,*overYe;
  float bandsX1[128],bandsX2[128];
  int nbands=0;
  int flagN=0;
  int paramOffset=0;

  char flagStrX[1024]="NULL";
  char flagStrY[1024]="NULL";

  for (i=0;i<100;i++)
    flagCol[i]= 1;

  if (strlen(bandsFile)>0)
    {
      FILE *fin;
      if (!(fin = fopen(bandsFile,"r")))
	{
	  printf("Unable to open file >%s<\n",bandsFile);
	  exit(1);
	}
      while (!feof(fin))
	{
	  if (fscanf(fin,"%f %f",&bandsX1[nbands],&bandsX2[nbands])==2)
	    nbands++;
	}
    }

  overX = (float *)malloc(sizeof(float)*MAX_OBSN);
  overY = (float *)malloc(sizeof(float)*MAX_OBSN);
  overYe = (float *)malloc(sizeof(float)*MAX_OBSN);

  px = (float *)malloc(sizeof(float)*MAX_OBSN);
  py = (float *)malloc(sizeof(float)*MAX_OBSN);
  x = (float *)malloc(sizeof(float)*MAX_OBSN);
  y = (float *)malloc(sizeof(float)*MAX_OBSN);
  yerr1 = (float *)malloc(sizeof(float)*MAX_OBSN);
  yerr2 = (float *)malloc(sizeof(float)*MAX_OBSN);
  x2 = (float *)malloc(sizeof(float)*MAX_OBSN);
  y2 = (float *)malloc(sizeof(float)*MAX_OBSN);
  yerr1_2 = (float *)malloc(sizeof(float)*MAX_OBSN);
  yerr2_2 = (float *)malloc(sizeof(float)*MAX_OBSN);
  freq = (float *)malloc(sizeof(float)*MAX_OBSN);

  // Setup colour scheme
  char fname[1000];
  char str[1000];
  FILE *fin;
  if (strlen(setupFile)>0)
    sprintf(fname,"%s",setupFile);
    //    sprintf(fname,"%s/plugin_data/%s",getenv(TEMPO2_ENVIRON),setupFile);
  else if (publish==1)
    sprintf(fname,"%s/plugin_data/plk_setup_publish.dat",getenv(TEMPO2_ENVIRON));
  else if (publish==2)
    sprintf(fname,"%s/plugin_data/plk_setup_image.dat",getenv(TEMPO2_ENVIRON));  
  else
    sprintf(fname,"%s/plugin_data/plk_setup.dat",getenv(TEMPO2_ENVIRON));	  
  if (!(fin = fopen(fname,"r")))
    {
      printf("WARNING: Unable to open file %s\n",fname);
      printf("Default plotting setup will be used\n");
      fontSize=1.0;
      aspect=0.8;
      viewport_x0 = 0.1;
      viewport_x1 = 0.95;
      viewport_y0 = 0.15;
      viewport_y1 = 0.85;
      fontType=1;
      lineWidth=1;
      minFreqCol[0] = 0;    maxFreqCol[0] = 500;   freqCol[0] = 1; freqStyle[0] = 16;
      minFreqCol[1] = 500;  maxFreqCol[1] = 1000;  freqCol[1] = 2; freqStyle[1] = 16;
      minFreqCol[2] = 1000; maxFreqCol[2] = 1700;  freqCol[2] = 3; freqStyle[2] = 16;
      minFreqCol[3] = 1700; maxFreqCol[3] = 3300;  freqCol[3] = 4; freqStyle[3] = 16;
      minFreqCol[4] = 3300; maxFreqCol[4] = 10000; freqCol[4] = 5; freqStyle[4] = 16;
      strcpy(bkgrdColour,"black");
      strcpy(lineColour,"white");
      freqColourNum=5;
    }
  else
    {
      while (!feof(fin))
	{
	  if (fscanf(fin,"%s",str)==1)
	    {
	      if (strcasecmp(str,"fontsize")==0)
		fscanf(fin,"%f",&fontSize);
	      else if (strcasecmp(str,"aspect")==0)
		fscanf(fin,"%f",&aspect);
	      else if (strcasecmp(str,"viewport_x0")==0)
		fscanf(fin,"%f",&viewport_x0);
	      else if (strcasecmp(str,"viewport_x1")==0)
		fscanf(fin,"%f",&viewport_x1);
	      else if (strcasecmp(str,"viewport_y0")==0)
		fscanf(fin,"%f",&viewport_y0);
	      else if (strcasecmp(str,"viewport_y1")==0)
		fscanf(fin,"%f",&viewport_y1);
	      else if (strcasecmp(str,"fonttype")==0)
		fscanf(fin,"%d",&fontType);
	      else if (strcasecmp(str,"linewidth")==0)
		fscanf(fin,"%d",&lineWidth);
	      else if (strcasecmp(str,"menu")==0)
		fscanf(fin,"%d",&menu);
	      else if (strcasecmp(str,"background")==0)
		fscanf(fin,"%s",bkgrdColour);
	      else if (strcasecmp(str,"line")==0)
		fscanf(fin,"%s",lineColour);
	      else if (strcasecmp(str,"freq")==0)
		{
		  fscanf(fin,"%f %f %d %d",&minFreqCol[freqColourNum],
			 &maxFreqCol[freqColourNum],&freqCol[freqColourNum],
			 &freqStyle[freqColourNum]);	     
		  freqColourNum++;
		}
	      else if (strcasecmp(str,"flag")==0)
		{
		  fscanf(fin,"%s %s %d %d",flagIDstr[flagColourNum],
			 &flagIDval[flagColourNum],&flagCol[flagColourNum],
			 &flagStyle[flagColourNum]);	     
		  flagColourNum++;
		}
	    }
	}
      fclose(fin);
    }


  if (psr[0].param[param_start].fitFlag[0]==1) 
    origStart = psr[0].param[param_start].val[0];
  if (psr[0].param[param_finish].fitFlag[0]==1) 
      origFinish = psr[0].param[param_finish].val[0];

  /* Obtain a graphical PGPLOT window */
  cpgbeg(0,gr,1,1);
  cpgpap(0,aspect);
  cpgscf(fontType);
  cpgslw(lineWidth);
  cpgscrn(0,bkgrdColour,&out);
  cpgscrn(1,lineColour,&out);
  if (strstr(gr,"/xs")==NULL)
    {
      exitFlag=1;
      fitFlag=2;
    }
  cpgask(0);


  sprintf(highlightID[highlightNum],"%s","-selectAll");
  sprintf(highlightVal[highlightNum],"%s","on");
  highlightNum++;
  min1 = psr[0].obsn[0].bat;
  max1 = psr[0].obsn[0].bat;
  for (i=0;i<psr[0].nobs;i++)
    {
      strcpy(psr[0].obsn[i].flagID[psr[0].obsn[i].nFlags],"-selectAll");
      strcpy(psr[0].obsn[i].flagVal[psr[0].obsn[i].nFlags],"off");
      psr[0].obsn[i].nFlags++;
      if (min1 > psr[0].obsn[i].bat) min1 = psr[0].obsn[i].bat;
      if (max1 < psr[0].obsn[i].bat) max1 = psr[0].obsn[i].bat;
    }



  if (strcmp(flagColour,"")!=0){
	  int found=0;
	  iFlagColour=1;
	  flagN=0;
	  for (i=0;i<psr[0].nobs;i++)
	  {
		  for (k=0;k<psr[0].obsn[i].nFlags;k++)
		  {
			  if (strcmp(psr[0].obsn[i].flagID[k],flagColour)==0)
			  {
				  found=0;
				  for (j=0;j<flagN;j++)
				  {
					  if (strcmp(psr[0].obsn[i].flagVal[k],flagStore[j])==0)
					  {
						  found=1;
						  break;
					  }
				  }
				  if (found==0)
				  {
					  strcpy(flagStore[flagN],psr[0].obsn[i].flagVal[k]);
					  flagN++;
				  }			
			  }
		  }	
	  }
  }


  do {
    if(debugFlag) 
      printf("Fitflag = %d\n",fitFlag);
    if (centre==-1)     centreEpoch = psr[0].param[param_pepoch].val[0];
    else if (centre==1)	centreEpoch = (min1+max1)/2.0;

    setLabel(xstr,xplot,plotPhase,unitFlag,centreEpoch,userValStr,flagStrX);
    setLabel(ystr,yplot,plotPhase,unitFlag,centreEpoch,userValStr,flagStrY);

    count=0;
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
	    freq[count]=(float)(psr[0].obsn[i].freq);
	    id[count] = i;
	    bad = setPlot(x,count,psr,i,unitFlag,plotPhase,xplot,&userValChange,userCMD,userValStr,userX,
			  centreEpoch,logx,flagStrX);
	    
	    bad = setPlot(y,count,psr,i,unitFlag,plotPhase,yplot,&userValChange,userCMD,userValStr,userX,
			  centreEpoch,logy,flagStrY);

	    if (yplot==1)                            /* Get pre-fit residual */
	      strcpy(fitType,"pre-fit");
	    else if (yplot==2)                       /* Post-fit residual    */
	      strcpy(fitType,"post-fit");
	    if (errtype==0)
		    errBar[count] = psr[0].obsn[i].toaErr;
	    else if (errtype==1)
		    errBar[count] = psr[0].obsn[i].origErr;
	    else 
		    errBar[count] = psr[0].obsn[i].toaDMErr;
	    if(yplot==16){errBar[count] = (float) psr[0].obsn[i].TNRedErr/1e-6;}
	    if(yplot==17){errBar[count] = (float) psr[0].obsn[i].TNDMErr/1e-6;}
	    if (bad==0) count++;	    
	  }
      }
    if (userValChange==2) userValChange=0;
    
    /* Remove mean from the residuals and calculate error bars */
    mean = findMean(y,psr,scale1,count);
    dmean = findMeanD(y,psr,scale1,count);
    for (i=0;i<count;i++)
      {
	if ((yplot==1 || yplot==2) && removeMean == 1) y[i]-=mean;

	if (plotPhase==-1)
	  {
	    yerr1[i] = y[i]-errBar[i]*1e-6/unitFlag;
	    yerr2[i] = y[i]+errBar[i]*1e-6/unitFlag;
	  }
	else
	  {
	    yerr1[i] = y[i]-errBar[i]*1e-6*(double)psr[0].param[param_f].val[0];
	    yerr2[i] = y[i]+errBar[i]*1e-6*(double)psr[0].param[param_f].val[0];
	  }
      }

    /* Sort into ascending x order */
    sort(x,y,yerr1,yerr2,freq,id,count);  
    if (origStart==-1)
      origStart = x[0] + centreEpoch - 1;
    if (origFinish==-1)
      origFinish = x[count-1] + 1 + centreEpoch;
  
    /* Get scaling for graph */
    minx = x[0];
    maxx = x[count-1];

    if (setZoomX1==0) plotx1 = minx-fabs(maxx-minx)*0.1;
    else plotx1 = zoomX1;

    if (setZoomX2==0) plotx2 = maxx+fabs(maxx-minx)*0.1;
    else plotx2 = zoomX2;

    miny = findMinY(y,x,count,plotx1,plotx2);
    maxy = findMaxY(y,x,count,plotx1,plotx2);
	
   if(yplot==16 || yplot == 17){
	miny = findMinY(yerr1,x,count,plotx1,plotx2);
    	maxy = findMaxY(yerr2,x,count,plotx1,plotx2);
	}
    
    if (setZoomY1==0) ploty1 = miny-fabs(maxy-miny)*0.1;
    else ploty1 = zoomY1;

    if (setZoomY2==0) ploty2 = maxy+fabs(maxy-miny)*0.1;
    else ploty2 = zoomY2;

    if (locky1 != locky2)
      {
	ploty1 = locky1;
	ploty2 = locky2;
      }

    if (lockx1 != lockx2)
    {
      plotx1 = lockx1;
      plotx2 = lockx2;
    }

    /* Plot the residuals */
    if (noreplot==0)
      {
	cpgsch(fontSize);
	if (menu>0 && graphics!=2)
	  {
	    cpgeras();
	    cpgsvp(0.0,1.0,0.8,1.0);
	    cpgswin(0,1,0,1);
	    /*	    cpgbox("",0.0,0,"",0.0,0); */
	    drawMenu(psr,plotx1,plotx2,ploty1,ploty2,menu,paramOffset);
	    if (menu==3)
	      {
		cpgsvp(0.0,0.3,0.3,0.8);
		cpgswin(0,1,0,1);
		    drawMenu3(psr,plotx1,plotx2,ploty1,ploty2,menu,xplot,yplot);
		    
		    cpgsvp(0.0,1.0,0.0,0.2);
		    cpgswin(0,1,0,1);
		    drawMenu3_2(psr,plotx1,plotx2,ploty1,ploty2,menu,xplot,yplot,jumpOffset,iFlagColour,flagN);
		    
		    cpgsvp(0.3,0.9,0.3,0.8);
		    cpgswin(plotx1,plotx2,ploty1,ploty2);
        //cpgbox("BCNST1",0.0,0,"BCNST1",0.0,0);
	      }
	    else
	      {
		cpgeras();
		cpgsvp(viewport_x0,viewport_x1,viewport_y0,viewport_y1);

		cpgswin(plotx1,plotx2,ploty1,ploty2);
		//cpgbox("BCNST1",0.0,0,"BCNST1",0.0,0);
	      }
	    // Print period axis (right-hand y-axis)
	    if( publish==0 && (yplot == 1 || yplot == 2 ) && displayPP==1){
        cpgbox( "BCNST", 0.0, 0, "BNST1", 0.0, 0 );
        cpgaxis( "N", plotx2, ploty1, plotx2, ploty2,
                 ploty1*psr[0].param[param_f].val[0]*unitFlag,
                 ploty2*psr[0].param[param_f].val[0]*unitFlag, 0.0, 0, 0.5, 0.0,
                 0.5, 0.3, 0.0 );
        cpgptxt( plotx2+(plotx2-plotx1)/15.0, (ploty2+ploty1)/2.0, 90.0, 
                 0.5, "Residual in pulse periods" );
	    }else
	      cpgbox( "BCNST1", 0.0, 0, "BCNST1", 0.0, 0 );
	  }
	else{
	  cpgeras();
	  cpgsvp(viewport_x0,viewport_x1,viewport_y0,viewport_y1);
	  cpgswin(plotx1,plotx2,ploty1,ploty2);
	  //cpgbox("BCNST1",0.0,0,"BCNST1",0.0,0);
	  // Print period axis (right-hand y-axis )
	  if(publish==0 &&( yplot == 1 || yplot == 2 ) && displayPP == 1){
	    cpgbox( "BCNST1", 0.0, 0, "BNST1", 0.0, 0 );
	    cpgaxis( "N", plotx2, ploty1, plotx2, ploty2,
		     ploty1*psr[0].param[param_f].val[0]*unitFlag,
		     ploty2*psr[0].param[param_f].val[0]*unitFlag, 0.0, 0, 0.5, 0.0,
		     0.5, 0.3, 0.0 );
	    cpgptxt( plotx2+(plotx2-plotx1)/15.0, (ploty2+ploty1)/2.0, 90.0, 
		     0.5, "Residual in pulse periods" );
	  }else
	    cpgbox( "BCNST1", 0.0, 0, "BCNST1", 0.0, 0 );
	}

  char rmsStr[10];
  if (cholmode) strcpy(rmsStr,"Crms");
  else {
	 if (psr[0].fitMode==1) strcpy(rmsStr,"Wrms");
	 else strcpy(rmsStr,"rms");
  }

  if (yplot==2) sprintf(title,"%s (%s = %.3f \\gms) %s",psr[0].name,rmsStr,psr[0].rmsPost,fitType);
  else sprintf(title,"%s (%s = %.3f \\gms) %s",psr[0].name,rmsStr,psr[0].rmsPre,fitType);

  if (showChisq == 1)
	 sprintf(title,"%s chisq=%.2f",title,psr[0].fitChisq/psr[0].fitNfree);
  if (nohead==0)
	 cpglab(xstr,ystr,title);
  else
	 cpglab(xstr,ystr,"");

  if (publish==0 && (strlen(setupFile)==0))
  {
	 cpgsch(0.5); 
	 cpgmtxt("B",6.5,0.92,0.0,"plk v.3.0 (G. Hobbs)"); 
	 cpgsch(fontSize); 
  }

  if (nbands > 0)
  {
	 int i;
	 float ffx[2],ffy[2];
	 printf("nbands = %d\n",nbands);
	 cpgsls(2);
	 for (i=0;i<nbands;i++)
	 {
		ffx[0] = bandsX1[i]- (double)centreEpoch;
		ffx[1] = bandsX2[i]- (double)centreEpoch;
		ffy[0] = ffy[1] = ploty2-(ploty2-ploty1)*0.05 - i*(ploty2-ploty1)*0.02;
		cpgsci((i%2)+1); cpgline(2,ffx,ffy);
	 }
	 cpgsls(1);
  }

  if (placeMarks==1)
  {
	 i=0;
	 cpgsls(3);
	 do
	 {
		px[0] = mark0 + i*markT;
		py[0] = ploty1;
		px[1] = px[0];
		py[1] = ploty2;
		cpgline(2,px,py);

		px[2] = mark0 - i*markT;
		py[2] = ploty1;
		px[3] = px[2];
		py[3] = ploty2;
		cpgline(2,px+2,py+2);
		i++;
	 } while (px[0] < plotx2 || px[2] > plotx1);
	 cpgsls(1);
  }
  if (psr[0].nPhaseJump > 0 && (xplot==3 || xplot==5)) // Display phase jumps
  {
	 char tstr[10];
	 float xch,ych;

	 cpgsls(4); cpgsci(7);
	 cpgsch(0.5); 
	 cpgqcs(4,&xch,&ych);
	 for (i=0;i<psr[0].nPhaseJump;i++)
	 {
		px[0] = (float)(psr[0].obsn[psr[0].phaseJumpID[i]].bat-centreEpoch);
		px[1] = px[0];
		py[0] = ploty1;
		py[1] = ploty2;
		if (psr[0].phaseJumpDir[i]==0)
		   cpgsci(8);
		else 
		   cpgsci(7);
		cpgline(2,px,py);
		sprintf(tstr,"%+d",psr[0].phaseJumpDir[i]);
		cpgtext(px[0]-xch/2.0,py[1]+ych/2.0,tstr);
	 }
	 cpgsls(1); cpgsci(1);
	 cpgsch(fontSize); 
  }
  if (nline > 0)
  {
	 //	    py[0]=miny;
	 //	    py[1]=maxy;
	 for (i=0;i<nline;i++)
	 {
		px[0] = lx1[i];
		px[1] = lx2[i];
		py[0] = ly1[i];
		py[1] = ly2[i];
		cpgline(2,px,py);
	 }			
  }
	  }
	noreplot=0;
	float x2[MAX_OBSN],y2[MAX_OBSN],yerr1_2[MAX_OBSN],yerr2_2[MAX_OBSN];
	if(plotPoints==1)
	{
	   for (j=0;j<freqColourNum;j++)
	   {
		  ncount=0;
		  for (i=0;i<count;i++)
		  {
			 if (freq[i] > minFreqCol[j] && freq[i] <= maxFreqCol[j])
			 {
				x2[ncount] = x[i];
				y2[ncount] = y[i];
				yerr1_2[ncount] = yerr1[i];
				yerr2_2[ncount] = yerr2[i];
				ncount++;
			 }
		  }
		  cpgsci(freqCol[j]);
		  //if (plotPoints==1)
		  //  {
		  cpgpt(ncount,x2,y2,freqStyle[j]);
		  if (plotErr==1 && (yplot==1 || yplot==2 || yplot==16 || yplot==17)) cpgerry(ncount,x2,yerr1_2,yerr2_2,1);
		  if (plotErr==2 && (yplot==1 || yplot==2)) cpgerry(ncount,x2,yerr1_2,yerr2_2,0);
		  //  }
		  if (join==1)
			 cpgline(ncount,x2,y2);
	   }
	   // Plot based on colours for flags 
	   if (flagColourNum>0)
	   {
		  for (i=0;i<count;i++)
		  {
			 for (k=0;k<psr[0].obsn[id[i]].nFlags;k++)
			 {
				for (j=0;j<flagColourNum;j++)
				{
				   if (strcasecmp(psr[0].obsn[id[i]].flagID[k],flagIDstr[j])==0
						 && strcasecmp(psr[0].obsn[id[i]].flagVal[k],flagIDval[j])==0)
				   {
					  cpgsci(flagCol[j]);
					  cpgpt(1,&x[i],&y[i],flagStyle[j]);
				   }
				}
			 }
		  }
	   }

	   if (iFlagColour==1)
	   {
		  int found=0,col;
		  for (i=0;i<count;i++)
		  {
			 for (k=0;k<psr[0].obsn[id[i]].nFlags;k++)
			 {
				if (strcmp(psr[0].obsn[id[i]].flagID[k],flagColour)==0)
				{
				   found=0;
				   for (j=0;j<flagN;j++)
				   {
					  if (strcmp(psr[0].obsn[id[i]].flagVal[k],flagStore[j])==0)
					  {
						 found=1;
						 col=j;
						 break;
					  }
				   }
				   if (found==0)
				   {
					  printf("WARNING: FLAG NOT FOUND\n");
					  col=1;
				   }
				   if (col+1 >= 15) col-=13;
				   if (publish==0)
				   {
					  cpgsci(col+1);
					  cpgpt(1,&x[i],&y[i],16);
				   }
				   else
					  cpgpt(1,&x[i],&y[i],col-12);
				   //		    cpgsci(7);
				   if (plotErr==1 && (yplot==1 || yplot==2)) cpgerry(1,&x[i],&yerr1[i],&yerr2[i],1);
				   if (plotErr==2 && (yplot==1 || yplot==2)) cpgerry(1,&x[i],&yerr1[i],&yerr2[i],0);

				}
			 }	
		  }
	   }
	}
	if (statistics==1)
	   displayStatistics(x,y,count,plotx1,plotx2,ploty1,ploty2);
	if (overPlotS==1)      
	   overPlotShapiro(psr,shapiroOffset,centreEpoch);
	if(plottingx==x[0] && overPlotn==0)overPlotn=1;
	if(overPlotn==1){
	   overPlotN(overN,overX,overY,overYe);
	   overPlotn=0;
	}
	if (indicateObs>0) /* Indicate individual observations */
	{
	   float px[2],py[2];
	   char str1[1000],str2[1000];
	   for (i=1;i<count;i++)
	   {
		  sscanf(psr[0].obsn[id[i-1]].fname,"%s",str1);
		  sscanf(psr[0].obsn[id[i]].fname,"%s",str2);
		  if (strcmp(str1,str2)!=0)
		  {
			 px[0] = 0.5*(x[i]+x[i-1]);
			 px[1] = px[0];
			 py[0] = ploty1;
			 py[1] = ploty2;
			 cpgline(2,px,py);
			 if (indicateObs==2 && px[0] > plotx1 && px[0] < plotx2)
			 {
				cpgsch(0.8);
				cpgptxt(px[0]-1,py[1],90,1.0,str1);
				cpgsch(fontSize);
			 }
		  }
	   }
	}
	if (plotFITWAVES==1 && removeFITWAVES==-1)
	{
	   float ymod;
	   float px[1000],py[1000];

	   cpgsls(3);
	   for (k=0;k<1000;k++)
	   {
		  px[k]=x[0]+((x[count-1]-x[0])/1000.0)*k; 
		  ymod=0.0;
		  for (j=0;j<FITWAVES_n/2;j++)
			 ymod+=FITWAVES_par[2*j]*cos(FITWAVES_omega/FITWAVES_harmonicStep*(j)*px[k])+
				FITWAVES_par[2*j+1]*sin(FITWAVES_omega/FITWAVES_harmonicStep*(j)*px[k]);
		  py[k]=ymod-mean;
	   }
	   cpgsls(3);
	   cpgline(1000,px,py);
	   cpgsls(1);
	}
	for (j=0;j<highlightNum;j++)
	{	
	   ncount=0;
	   for (i=0;i<count;i++)
	   {
		  for (k=0;k<psr[0].obsn[id[i]].nFlags;k++)
		  {
			 if (strcmp(psr[0].obsn[id[i]].flagID[k],highlightID[j])==0
				   && strcmp(psr[0].obsn[id[i]].flagVal[k],highlightVal[j])==0)
			 {
				x2[ncount] = x[i];
				y2[ncount] = y[i];
				yerr1_2[ncount] = yerr1[i];
				yerr2_2[ncount] = yerr2[i];
				ncount++;
			 }
		  }
	   }
	   cpgsci(1);	  
	   cpgsch(2); 	
	   cpgpt(ncount,x2,y2,4+j);
	   /*	if (plotErr==1) cpgerry(ncount,x2,yerr1_2,yerr2_2,1);	   */
	   cpgsch(fontSize); 
	}
	if (label==1)
	{
	   cpgtext(labelX,labelY,labelStr);
	}
	if (viewModel==1 && (fitFlag==1 || fitFlag==3 || fitFlag==7)){
	   viewModels(psr,minx,maxx,centreEpoch,removeMean,dmean,count,id,fitFlag,x,y);
	}

	cpgsci(1);
	if (strstr(gr,"/xs")!=NULL && graphics!=2)
	{
	   //	cpgcurs(&mouseX,&mouseY,&key);
	   cpgband(0,0,0,0,&mouseX,&mouseY,&key);
	   /* Check key press */
	   if (key=='q') exitFlag=1;
	   else if (key=='1') {if ((xplot!=3 && yplot!=1)){setZoomX1 = 0; setZoomX2 = 0;} xplot=3; yplot=1;fitFlag=1; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='2') {if ((xplot!=3 && yplot!=2)){setZoomX1 = 0; setZoomX2 = 0;} xplot=3; yplot=2;fitFlag=2;setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='3' && psr[0].param[param_pb].paramSet[0]==1) {xplot=4;yplot=1;fitFlag=3;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='4' && psr[0].param[param_pb].paramSet[0]==1) {xplot=4;yplot=2;fitFlag=4;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='5') {xplot=5; yplot=1;fitFlag=5;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='6') {xplot=5; yplot=2;fitFlag=6;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='7') {xplot=6; yplot=1;fitFlag=7;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='8') {xplot=6; yplot=2;fitFlag=8;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='9') {xplot=7; yplot=1;fitFlag=9;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='a') {xplot=7; yplot=2;fitFlag=10;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='!') {xplot=8; yplot=1;fitFlag=11;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='@') {xplot=8; yplot=2;fitFlag=12;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='#') {xplot=9; yplot=1;fitFlag=13;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='$') {xplot=9; yplot=2;fitFlag=14;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='%') {xplot=10; yplot=1;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='^') {xplot=10; yplot=2;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='&') {xplot=11; yplot=1;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='*') {xplot=11; yplot=2;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='(') {xplot=12; yplot=1;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key==')') {xplot=12; yplot=2;setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='.') {xplot=2; setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key==16) {plotPhase*=-1; setZoomY1 = 0; setZoomY2 =0;}
	   else if (key=='O') { // Select flag ID for display
	     printf("Enter flag ID for x-axis (NULL for nothing) ");
	     scanf("%s",flagStrX);
	     printf("Enter flag ID for y-axis (NULL for nothing) ");
	     scanf("%s",flagStrY);
	     if (strcmp(flagStrX,"NULL")!=0)
	       xplot=18;
	     if (strcmp(flagStrY,"NULL")!=0)
	       yplot=18;
	   }
	    
	   else if (key=='>') /* Select next point to right if available */
	   {
		  if (setZoomX2 != 0)
		  {
			 for (i=0;i<count;i++)
			 {
				if (x[i] > zoomX2)
				{
				   zoomX2 = x[i]+1.0;
				   break;
				}
			 }
		  }
		  setZoomY1 = 0; setZoomY2 = 0;
	   }
	   else if (key=='<') /* Select next point to left if available */
	   {
		  if (setZoomX1 != 0)
		  {
			 for (i=count-1;i>0;i--)
			 {
				if (x[i] < zoomX1)
				{
				   zoomX1 = x[i]-1.0;
				   break;
				}
			 }
		  }
		  setZoomY1 = 0; setZoomY2 = 0;
	   }
		else if(key=='F'){
			if(psr[0].TNsubtractDM==0){
				printf("will substract PL DM Variations on next Fit \n");
				psr[0].TNsubtractDM=1;
			}
			else if(psr[0].TNsubtractDM==1){
                                printf("will Re-add PL DM Variations on next Fit \n");
                                psr[0].TNsubtractDM=0;
                        }

		}
		else if(key=='K'){
			if(psr[0].TNsubtractRed==0){
	                        printf("will substract Red Noise on next Fit \n");
        	                psr[0].TNsubtractRed=1;
			}
			else if(psr[0].TNsubtractRed==1){
                                printf("will Re-add Red Noise on next Fit \n");
                                psr[0].TNsubtractRed=0;
                        }

                }

	   else if (key==3) /* Change central point */
		  centre*=-1;
	   else if (key=='E') /* Toggle plotting error bars */
	   {
		  plotErr+=1;
		  if (plotErr==3) plotErr=0;
	   }
	   else if (key=='\''){
		  errtype++;
		  switch (errtype){
			 case 1:
				printf("Orig Error\n");
				break;
			 case 2:
				printf("DM Error\n");
				break;
			 default:
				printf("Combined Error\n");
				errtype=0;
				break;
		  }
	   }
	   else if (key=='"'){
		  for(i=0; i < psr[0].dmoffsCMnum; i++){
			 cpgsci(7);
			 cpgsls(4);
			 cpgmove((float)(psr[0].dmoffsCM_mjd[i]-centreEpoch),-1e-2);
			 cpgdraw((float)(psr[0].dmoffsCM_mjd[i]-centreEpoch),1e-2);
			 cpgsls(1);
		  }
		  for(i=0; i < psr[0].dmoffsDMnum; i++){
			 cpgsci(6);
			 cpgsls(3);
			 cpgmove((float)(psr[0].dmoffsDM_mjd[i]-centreEpoch),-1e-2);
			 cpgdraw((float)(psr[0].dmoffsDM_mjd[i]-centreEpoch),1e-2);
			 cpgsls(1);
		  }


		  noreplot=1;
		  continue;
	   }

	   else if (key=='P') /* New parameter file */
	   {
		  textOutput(psr,npsr,0,0,0,1,(char *)"");
	   }
	   else if (key==23) /* over-write par file */
	   {
		  textOutput(psr,npsr,0,0,0,1,parFile[0]);
	   }
	   /*	else if (key==61) // Add phase jump above the cursor 
			{
			for (i=0;i<count;i++)
			{
			if (y[i]>mouseY)
			{
			psr[0].obsn[id[i]].sat-=(1.0/psr[0].param[param_f].val[0])/2.0;
			}
			}
			callFit(psr,npsr);
			}*/
			else if (key==2) // ctrl-b (run pdv to get backend jumps)
			{
			   char str[1024];
			   int i;
			   for (i=0;i<count;i++)
			   {
				  printf("Have %s\n",psr[0].obsn[id[i]].fname);
			   }
			}
			else if (key=='b'){ /* bin the residuals */
			   binResiduals(psr,npsr,x,y,count,id,&overN,overX,overY,overYe,xplot,yplot,
					 yerr1,unitFlag,plotPhase,(double)centreEpoch);
			   overPlotn = 1;
			   plottingx = x[0]; 
			}
			else if (key=='h') help();
			else if (key=='D') //view profile (same as middle button)
			{
			   char str[1000];
			   int closest;
			   closest = idPoint(psr,x,y,id,count,mouseX,mouseY); /* Identify closest point */
			   strcpy(psr[0].obsn[closest].flagVal[psr[0].obsn[closest].nFlags-1],"on");
			   if (psr[0].jboFormat==1)
			   {
				  sprintf(str,"viewJodrell psrav.prf %.14lf",(double)psr[0].obsn[closest].sat);
				  printf("Attempting to view Jodrell: %s\n",str);
			   }
			   else
			   {
				  int ichan=-1;
				  int isub=-1;
				  char options1[100]="";
				  char options2[100]="";
				  char options4[1000]="";
				  // Check flags to see if -chan, -sub are set
				  for (k=0;k<psr[0].obsn[closest].nFlags;k++)
				  {
					 if (strcmp(psr[0].obsn[closest].flagID[k],"-chan")==0)
						sscanf(psr[0].obsn[closest].flagVal[k],"%d",&ichan);
					 if (strcmp(psr[0].obsn[closest].flagID[k],"-sub")==0
						   || strcmp(psr[0].obsn[closest].flagID[k],"-subint")==0)
						sscanf(psr[0].obsn[closest].flagVal[k],"%d",&isub);
				  }
				  if (ichan!=-1) sprintf(options1," -H %d ",ichan);
				  if (isub!=-1) sprintf(options2," -I %d ",isub);
				  if (ichan==-1 && isub==-1)
					 sprintf(options4," -FTp");
				  else
					 sprintf(options4," -p ");
				  sprintf(str,"pav -D%s %s %s -g10/xs %s",options4,options1,options2,psr[0].obsn[closest].fname);
				  printf("RUNNING %s\n",str);
			   }
			   system(str);
			}
			else if (key==4) /* Delete selected points (ctrl-D) */
			{
			   for (i=0;i<psr[0].nobs;i++)
			   {
				  if (strcmp(psr[0].obsn[i].flagVal[psr[0].obsn[i].nFlags-1],"on")==0)
					 psr[0].obsn[i].deleted=1;
			   }
			   userValChange=1;
			}
			else if (key==26) /* ctrl-z list selected points (ctrl-D) */
			{
			   int len=0;
			   char tstr[1000],nstr[1000];

			   for (i=0;i<psr[0].nobs;i++)
			   {
				  if (strcmp(psr[0].obsn[i].flagVal[psr[0].obsn[i].nFlags-1],"on")==0)
				  {
					 if (len < strlen(psr[0].obsn[i].fname)) 
						len = strlen(psr[0].obsn[i].fname);
				  }
			   }
			   sprintf(tstr,"\%%-%d.%ds \%%9.3f \%%.5f \%%7.3f \%%10.3f",len+2,len+2);
			   for (i=0;i<psr[0].nobs;i++)
			   {
				  if (strcmp(psr[0].obsn[i].flagVal[psr[0].obsn[i].nFlags-1],"on")==0)
				  {
					 sprintf(nstr,tstr,psr[0].obsn[i].fname,(double)psr[0].obsn[i].freq,(double)psr[0].obsn[i].sat,(double)psr[0].obsn[i].toaErr,(double)psr[0].obsn[i].residual/1e-6);
					 printf("%s\n",nstr);
				  }
			   }
			   userValChange=1;
			}
			else if (key=='v') /* View all highlighted profiles */
			{
			   char str[1000],str2[1000];
			   int xn,yn;
			   printf("Enter Nx Ny "); scanf("%d %d",&xn,&yn);
			   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
			   strcpy(str,"pav -DFTCp -g 10/xs ");
			   sprintf(str2,"-N %d,%d ",xn,yn);
			   strcat(str,str2);
			   for (j=0;j<highlightNum;j++)
			   {	
				  ncount=0;
				  for (i=0;i<count;i++)
				  {
					 for (k=0;k<psr[0].obsn[id[i]].nFlags;k++)
					 {
						if (strcmp(psr[0].obsn[id[i]].flagID[k],highlightID[j])==0
							  && strcmp(psr[0].obsn[id[i]].flagVal[k],highlightVal[j])==0)
						{
						   strcat(str,psr[0].obsn[id[i]].fname);
						   strcat(str," ");
						}
					 }
				  }
			   }
			   system(str);
			}
			else if (key=='L') /* Add label to plot */
			{
			   char labelS[1024];
			   printf("Enter string:\n");
			   fgets(labelS,1024,stdin);
			   labelS[strlen(labelS)-1]='\0';
			   strcpy(labelStr,labelS);
			   printf("Enter x y: ");
			   scanf("%f %f",&labelX,&labelY);
			   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)

			   label=1;
			}
			else if (key=='C') /* Run UNIX command */
			{
			   char str[1000];
			   printf("Command: "); gets(str);
			   for (j=0;j<highlightNum;j++)
			   {	
				  ncount=0;
				  for (i=0;i<count;i++)
				  {
					 for (k=0;k<psr[0].obsn[id[i]].nFlags;k++)
					 {
						if (strcmp(psr[0].obsn[id[i]].flagID[k],highlightID[j])==0
							  && strcmp(psr[0].obsn[id[i]].flagVal[k],highlightVal[j])==0)
						{
						   strcat(str," ");
						   strcat(str,psr[0].obsn[id[i]].fname);
						}
					 }
				  }
			   }
			   printf("Running: %s\n",str);
			   system(str);
			}
			else if (key=='j') /* Draw line between points */
			   join*=-1;
			else if (key=='J') /* Toggle plotting points */
			   plotPoints*=-1;
			else if (key=='Q') /* Output deleted files for MySQL query */
			{
			   FILE *foutSQL;
			   foutSQL = fopen("sql_delete","w");
			   printf("Deleted\n");
			   for (i=0;i<psr[0].nobs;i++)
			   {
				  if (psr[0].obsn[i].deleted==1)
				  {
					 printf("%s\n",psr[0].obsn[i].fname);
					 fprintf(foutSQL,"%s\n",psr[0].obsn[i].fname);
				  }
			   }
			   fclose(foutSQL);
			}
			else if (key=='m') /* Measure */
			{
			   float mouseX2,mouseY2;
			   cpgband(1,0,mouseX,mouseY,&mouseX2,&mouseY2,&key);	
			   printf("Measure x = %.3g days\t y = %.3g\t gradient = %.3g\n",mouseX2-mouseX,mouseY2-mouseY,(mouseY2-mouseY)/(mouseX2-mouseX));
			}
			else if (key=='M') /* Toggle removing mean */
			{
			   removeMean*=-1;
			   if (removeMean==1)
				  printf("Removing mean from timing residuals\n");
			   else
				  printf("Not removing mean from timing residuals\n");
			}
			else if (key==9) /* ctrl-i - highlight points with colour */
			{
			   if (iFlagColour==1)
				  iFlagColour=0;
			   else
			   {
				  int found=0;
				  printf("Enter flag identifier ");
				  scanf("%s",flagColour);
				  getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
				  iFlagColour=1;
				  flagN=0;
				  for (i=0;i<psr[0].nobs;i++)
				  {
					 for (k=0;k<psr[0].obsn[i].nFlags;k++)
					 {
						if (strcmp(psr[0].obsn[i].flagID[k],flagColour)==0)
						{
						   found=0;
						   for (j=0;j<flagN;j++)
						   {
							  if (strcmp(psr[0].obsn[i].flagVal[k],flagStore[j])==0)
							  {
								 found=1;
								 break;
							  }
						   }
						   if (found==0)
						   {
							  strcpy(flagStore[flagN],psr[0].obsn[i].flagVal[k]);
							  flagN++;
						   }			
						}
					 }	
				  }
			   }

			}
			else if (key=='S') /* Save new .tim file */
			   newTim(psr);
			else if (key=='n') /* Save new .tim file with pulse numbers */
			{
			   logmsg("Writing pulse numbers to withpn.tim");
			   for (i=0;i<psr[0].nobs;i++)
			   {
				  int flagid=psr[0].obsn[i].nFlags;
				  for(int k=0; k < flagid ; k++){
					 if(strcmp(psr[0].obsn[i].flagID[k],"-pnadd")==0){
						printf("Removing -pnadd flag\n");
						for (int kk=k; kk < flagid-1; kk++){
						   strcpy(psr[0].obsn[i].flagID[kk],psr[0].obsn[i].flagID[kk+1]);
						   strcpy(psr[0].obsn[i].flagVal[kk],psr[0].obsn[i].flagVal[kk+1]);
						}
						flagid-=1;
						k-=1;
					 }
				  }
				  psr[0].obsn[i].nFlags = flagid;
				  for(int k=0; k < flagid ; k++){
					 if(strcmp(psr[0].obsn[i].flagID[k],"-pn")==0){
						flagid=k;
						break;
					 }
				  }
				  //printf("%g %lld\n",(double)psr[0].obsn[i].sat,(psr[0].obsn[i].pulseN - psr[0].obsn[0].pulseN));
				  strcpy(psr[0].obsn[i].flagID[flagid],"-pn");
				  sprintf(psr[0].obsn[i].flagVal[flagid],"%lld",psr[0].obsn[i].pulseN-psr[0].obsn[0].pulseN);
				  if (flagid==psr[0].obsn[i].nFlags)
					 psr[0].obsn[i].nFlags++;
			   }
			   writeTim("withpn.tim",psr,"tempo2");
			}
			else if (key==19) /* over-write .tim file */
			   writeTim(timFile[0],psr,"tempo2");
			else if (key=='t') /* Toggle displaying statistics */
			   statistics*=-1;
			else if (key=='H') /* Highlight */
			{
			   printf("FlagID = ");  scanf("%s",highlightID[highlightNum]);
			   printf("FlagVal = "); scanf("%s",highlightVal[highlightNum++]);
			   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
			}
			else if (key=='-' || key=='+') /*  phase jump */
			{
			   // Find closest point to the left
			   int i,iclosest;
			   float closest,x1,x2,x3,x4,y1,y2,y3,y4,xscale,yscale,xpos,ypos;

			   cpgqvp(3,&x1,&x2,&y1,&y2);
			   cpgqwin(&x3,&x4,&y3,&y4);
			   xscale = (x2-x1)/(x4-x3);
			   yscale = (y2-y1)/(y4-y3);
			   mouseX = (mouseX-x3)*xscale;
			   mouseY = (mouseY-y3)*yscale;
			   iclosest=-1;
			   for (i=0;i<count;i++)
			   {
				  xpos = (x[i]-x3)*xscale;
				  ypos = (y[i]-y3)*yscale;
				  if (iclosest==-1 && xpos < mouseX)
				  {
					 iclosest=id[i];
					 closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
				  }
				  else if (fabs(xpos-mouseX)<closest && xpos < mouseX)
				  {
					 iclosest=id[i];
					 closest = fabs(xpos-mouseX);
				  }
			   }
			   psr[0].phaseJump[psr[0].nPhaseJump] = psr[0].obsn[iclosest].sat+1.0/SECDAY;
			   psr[0].phaseJumpID[psr[0].nPhaseJump] = iclosest;
			   //	    printf("Have %s %d\n",psr[0].obsn[iclosest].fname,iclosest);
			   if (key=='-') psr[0].phaseJumpDir[psr[0].nPhaseJump] = -1;
			   else psr[0].phaseJumpDir[psr[0].nPhaseJump] = 1;
			   for (i=0;i<psr[0].nobs;i++)
			   {
				  k=psr[0].nPhaseJump;
				  if (psr[0].obsn[i].sat > psr[0].obsn[psr[0].phaseJumpID[k]].sat)
				  {
					 psr[0].obsn[i].residual += (double)psr[0].phaseJumpDir[k]/psr[0].param[param_f].val[0];
					 psr[0].obsn[i].prefitResidual += (double)psr[0].phaseJumpDir[k]/psr[0].param[param_f].val[0];
					 psr[0].obsn[i].pulseN -= psr[0].phaseJumpDir[k];
				  }
			   }
			   //if (key=='-') 
			   //  {
			   setZoomY1 = 0; setZoomY2 = 0;
			   //    /*      zoomY1 -= 1.0/psr[0].param[param_f0].val;
			   //            zoomY2 += 1.0/psr[0].param[param_f0].val; */
			   //  } 
			   //if (key=='+') 
			   //  {
			   //    setZoomY1 = 0; setZoomY2 = 0;
			   //    /*      zoomY1 -= 1.0/psr[0].param[param_f0].val;
			   //            zoomY2 += 1.0/psr[0].param[param_f0].val; */
			   //  } 
			   psr[0].nPhaseJump++;
			}
			else if (key==8) // Backspace - remove all phase jumps
			{
			   for (k=0;k<psr[0].nPhaseJump;k++)
			   {
				  for (i=0;i<psr[0].nobs;i++) /* Get rid of existing jumps */
				  {

					 if (psr[0].obsn[i].sat > psr[0].phaseJump[k])
					 {
						psr[0].obsn[i].residual -= (double)psr[0].phaseJumpDir[k]/psr[0].param[param_f].val[0];
						psr[0].obsn[i].prefitResidual -= (double)psr[0].phaseJumpDir[k]/psr[0].param[param_f].val[0];
					 }
				  }
			   }
			   setZoomY1 = 0;
			   setZoomY2 = 0;
			   psr[0].nPhaseJump = 0;

			}
			else if (key=='o') /* Highlight all */
			{
			   for (i=0;i<count;i++)
			   {
				  if ((setZoomX1==0 || x[i]>zoomX1) &&
						(setZoomX2==0 || x[i]<zoomX2) &&
						(setZoomY1==0 || y[i]>zoomY1) &&
						(setZoomY2==0 || y[i]<zoomY2))
					 strcpy(psr[0].obsn[id[i]].flagVal[psr[0].obsn[id[i]].nFlags-1],"on");
			   }
			}
			else if (key==5) /* Highlight points more than 3 sigma from the mean */
			{	    
			   for (i=0;i<count;i++)
			   {
				  if (fabs(y[i])/(psr[0].obsn[id[i]].toaErr*1e-6) > 3)
					 strcpy(psr[0].obsn[id[i]].flagVal[psr[0].obsn[id[i]].nFlags-1],"on");
				  else
					 strcpy(psr[0].obsn[id[i]].flagVal[psr[0].obsn[id[i]].nFlags-1],"off");
			   }
			}
			else if (key=='N') /* Highlight point with a given filename */
			{	    
			   char filen[100];
			   printf("Please enter filename (or part of filename) ");
			   scanf("%s",filen);
			   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
			   for (i=0;i<count;i++)
			   {
				  if (strstr(psr[0].obsn[id[i]].fname,filen)!=NULL)
					 strcpy(psr[0].obsn[id[i]].flagVal[psr[0].obsn[id[i]].nFlags-1],"on");
			   }
			}
			else if (key==21)  /* Overplot Shapiro delay */
			{     // Ctrl+u
			   overPlotS*=-1;
			   if (overPlotS==1)
			   {
				  printf("Enter constant offset ");
				  scanf("%f",&shapiroOffset);
				  getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
			   }
			}
			else if (key=='U') /* Unselect all */
			{
			   for (i=0;i<psr[0].nobs;i++)
			   {
				  if (psr[0].obsn[i].deleted == 0 &&
						(psr[0].param[param_start].paramSet[0]!=1 || psr[0].param[param_start].fitFlag[0]!=1 ||
						 psr[0].param[param_start].val[0] < psr[0].obsn[i].bat) &&
						(psr[0].param[param_finish].paramSet[0]!=1 || psr[0].param[param_finish].fitFlag[0]!=1 ||
						 psr[0].param[param_finish].val[0] > psr[0].obsn[i].bat))

				  {
					 strcpy(psr[0].obsn[i].flagVal[psr[0].obsn[i].nFlags-1],"off");
				  }
			   }
			}
			else if (key=='e') /* Multiply errors */
			{
			   double errMult;
			   printf("Enter error multiplication factor "); scanf("%lf",&errMult);
			   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
			   for (i=0;i<psr[0].nobs;i++){
				  psr[0].obsn[i].efac *= errMult;
				  psr[0].obsn[i].toaErr*=errMult;
			   }
			   for (i=0;i<count;i++)
				  errBar[i]*=errMult;
			}
			else if (key=='w')  /* Toggle fitting with weights */
			{
			   if (psr[0].fitMode==0) psr[0].fitMode=1;
			   else psr[0].fitMode=0;
			   if (psr[0].fitMode==1)
			   {
				  for (i=0;i<psr[0].nobs;i++)
				  {
					 if (psr[0].obsn[i].toaErr == 0)
					 {
						printf("ERROR: the error on the arrival time for MJD %f is zero - cannot fit using weights \n",(double)psr[0].obsn[i].sat);
						psr[0].fitMode=0;
					 }
				  }
			   }
			}
			else if (key=='[') /* Toggle x-axis log scale */
			   logx*=-1;
			else if (key==']')
			   logy*=-1;
			else if (key=='l') { /* List all */
			   for (i=0;i<count;i++)
			   {
				  if ((setZoomX1==0 || x[i]>zoomX1) &&
						(setZoomX2==0 || x[i]<zoomX2) &&
						(setZoomY1==0 || y[i]>zoomY1) &&
						(setZoomY2==0 || y[i]<zoomY2))
				  {
					 printf("%s %s\t%7.5g\t%.5g\t%.5g %.5g ",psr[0].obsn[id[i]].fname,
						print_longdouble(psr[0].obsn[id[i]].sat).c_str(),x[i],y[i],
						(double)psr[0].obsn[id[i]].toaErr*1e-6,
						(double)psr[0].obsn[id[i]].freq);
					 for (j=0;j<psr[0].obsn[id[i]].nFlags;j++)
					   printf("%s %s ",psr[0].obsn[id[i]].flagID[j],psr[0].obsn[id[i]].flagVal[j]);
					 printf("\n");
				  }
			   }
			   printf("Mean residual = %f\n",mean);
			}
			else if (key==15) {  /* Simple output listing */
			   char fname[1000];
			   FILE *outfile;
			   float recX[MAX_OBSN];
			   int smooth=0;
			   int  nrec=0,used;
			   float meanX[MAX_OBSN],meanY[MAX_OBSN];
			   int nMean=0;

			   printf("Enter filename ");
			   scanf("%s",fname);
			   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
			   printf("Enter smoothing (0 for no smoothing) ");
			   scanf("%d",&smooth);
			   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
			   outfile = fopen(fname,"w");
			   if (smooth>0)
			   {
				  averagePts(x,y,count,smooth,meanX,meanY,&nMean);

				  for (i=0;i<nMean;i++)
					 fprintf(outfile,"%12.5f %13.6g 1.0\n",meanX[i],meanY[i]);
			   }
			   else
			   {
				  for (i=0;i<psr[0].nobs;i++)
				  {
					 if (psr[0].obsn[i].deleted == 0 &&
						   (psr[0].param[param_start].paramSet[0]!=1 || psr[0].param[param_start].fitFlag[0]!=1 ||
							psr[0].param[param_start].val[0] < psr[0].obsn[i].bat) &&
						   (psr[0].param[param_finish].paramSet[0]!=1 || psr[0].param[param_finish].fitFlag[0]!=1 ||
							psr[0].param[param_finish].val[0] > psr[0].obsn[i].bat))
						fprintf(outfile,"%12.5f %13.6g %12.5f\n",(double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]),(double)psr[0].obsn[i].residual,(double)psr[0].obsn[i].toaErr);
				  }
				  /*	      for (i=0;i<count;i++)
							  {
							  if ((setZoomX1==0 || x[i]>zoomX1) &&
							  (setZoomX2==0 || x[i]<zoomX2) &&
							  (setZoomY1==0 || y[i]>zoomY1) &&
							  (setZoomY2==0 || y[i]<zoomY2))
							  {
							  used=0;
							  for (j=0;j<nrec;j++)
							  {
							  if (x[i]==recX[j])
							  used=1;
							  }
							  if (used==0)
							  {
							  fprintf(outfile,"%12.5f %13.6g %12.5f",x[i],y[i],
							  psr[0].obsn[id[i]].toaErr); 

							  recX[nrec]=x[i];
							  nrec++;
							  }
							  else
							  fprintf(outfile,"%12.5f %13.6g %12.5f **",x[i],y[i],
							  psr[0].obsn[id[i]].toaErr); 

							  fprintf(outfile,"\n");
							  }
							  }*/ 
			   }
			   fclose(outfile); 

			}
			else if (key==10) { /* List all in Jodrell format */
			   char fname[1000];
			   FILE *outfile;

			   printf("Enter filename ");
			   scanf("%s",fname);
			   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
			   outfile = fopen(fname,"w");
			   for (i=0;i<count;i++)
			   {
				  if ((setZoomX1==0 || x[i]>zoomX1) &&
						(setZoomX2==0 || x[i]<zoomX2) &&
						(setZoomY1==0 || y[i]>zoomY1) &&
						(setZoomY2==0 || y[i]<zoomY2))
				  {
					 printf("%12.5f %12.5f %12.5f %12.5f\n",x[i],y[i]*1000.0,psr[0].obsn[id[i]].toaErr/1000.0,
						   (double)psr[0].obsn[id[i]].bat);
					 fprintf(outfile,"%12.5f %13.6f %12.5f %20.15f %20.15f %12.5f %s",x[i],y[i]*1000.0,
						   psr[0].obsn[id[i]].toaErr/1000.0, 
						   (double)psr[0].obsn[id[i]].bat,(double)(psr[0].obsn[id[i]].bat-
							  psr[0].param[param_pepoch].val[0]),(double)psr[0].obsn[id[i]].freq,psr[0].obsn[id[i]].fname); 
					 /*		  fprintf(outfile,"%g %g ",y[i],(double)psr[0].obsn[id[i]].residual-mean); */
					 for (j=0;j<psr[0].obsn[id[i]].nFlags;j++)
					 {
						if (strcmp(psr[0].obsn[id[i]].flagID[j],"-selectAll")!=0)
						   fprintf(outfile," %s",psr[0].obsn[id[i]].flagVal[j]);
					 }
					 fprintf(outfile,"\n");
				  }
			   }
			   fclose(outfile);
			}
			else if (key=='y')
			{
			   setZoomY1 = 0;
			   setZoomY2 = 0;
			}
			else if (key=='Y')
			{
			   printf("Enter y-range: minimum maximum ");
			   scanf("%f %f",&zoomY1,&zoomY2);
			   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
			   setZoomY1 = 1;
			   setZoomY2 = 1;
			}
			else if (key=='B') { /* Place periodic marks on the x-axis */
			   placeMarks*=-1;
			   if (placeMarks==1)
			   {
				  mark0 = mouseX;
				  printf("Enter periodicity (d) ");
				  scanf("%f",&markT);
				  getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
			   }     
			}
			else if (key==24) /* cntr-X -- set X-scale */
			{      
			   cpgcurs(&mouseX,&mouseY,&key);
			   setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;
			   if (key=='1') xplot=1; // Prefit residuals	     
			   if (key=='2') xplot=2; // Postfit residuals      
			   if (key=='3') xplot=3; // MJD		     
			   if (key=='4') xplot=4; // Orbital Phase	     
			   if (key=='5') xplot=5; // TOA number	     
			   if (key=='6') xplot=6; // Day of Year	     
			   if (key=='7') xplot=7; // Frequency		     
			   if (key=='8') xplot=8; // TOA error		     
			   if (key=='9') xplot=9; // User-defined variable  
			   if (key=='0') xplot=10;// Year		     
			   if (key=='-') xplot=11;// Elevation		     
			   if (key=='=') xplot=12;// round MJD              
			   setLabel(xstr,xplot,plotPhase,unitFlag,centreEpoch,userValStr,flagStrX);
			   setLabel(ystr,yplot,plotPhase,unitFlag,centreEpoch,userValStr,flagStrY);
			}
			else if (key==25) /* cntr-Y -- set Y-scale */
			{      
			   cpgcurs(&mouseX,&mouseY,&key);
			   setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 =0;
			   if (key=='1') yplot=1;  // Prefit residuals	     
			   if (key=='2') yplot=2;  // Postfit residuals          
			   if (key=='3') yplot=3;  // MJD		     
			   if (key=='4') yplot=4;  // Orbital Phase	     
			   if (key=='5') yplot=5;  // TOA number	     
			   if (key=='6') yplot=6;  // Day of Year	     
			   if (key=='7') yplot=7;  // Frequency		     
			   if (key=='8') yplot=8;  // TOA error		     
			   if (key=='9') yplot=9;  // User-defined variable  
			   if (key=='0') yplot=10; // Year		          
			   if (key=='-') yplot=11; // Elevation		     
			   if (key=='=') yplot=12; // round MJD              
			   setLabel(xstr,xplot,plotPhase,unitFlag,centreEpoch,userValStr,flagStrX);
			   setLabel(ystr,yplot,plotPhase,unitFlag,centreEpoch,userValStr,flagStrY);

			}
			else if (key=='x') {  /* Do fit, but define start and finish by the zoom           */
			   /* Need to update doFit.C to take notice of START and FINISH */
			   /* for this to work properly                                 */

			   reFit(fitFlag,setZoomX1,setZoomX2,zoomX1,zoomX2,origStart,origFinish,centreEpoch,psr,npsr,xplot,dcmFile,covarFuncFile,zoom);
			   setZoomY1 = 0;
			   setZoomY2 = 0;
			}
			else if (key==14) { /* ctrl-n: add noise */
			   double noiseAdd=0.0;
			   long idum=TKsetSeed();
			   printf("Please enter amount of noise to add in micro-sec ");
			   scanf("%lf",&noiseAdd);
			   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
			   for (i=0;i<psr[0].nobs;i++)
				  psr[0].obsn[i].sat+=((noiseAdd*1e-6*TKgaussDev(&idum))/86400.0);
			   reFit(fitFlag,setZoomX1,setZoomX2,zoomX1,zoomX2,origStart,origFinish,centreEpoch,psr,npsr,xplot,dcmFile,covarFuncFile,zoom);

			   setZoomY1 = 0;
			   setZoomY2 = 0;
			}
			else if (key=='g') {  /* Change graphics device */
			   cpgend();
			   cpgbeg(0,"?",1,1);
			   cpgpap(0,aspect);
			   //	  	  cpgpap(0,0.4);
			   cpgscf(fontType);
			   cpgslw(lineWidth);
			   //	  cpgscrn(0,bkgrdColour,&out);
			   //	  cpgscrn(1,lineColour,&out);
			   graphics=1;
			   /*      cpgenv(plotx1,plotx2,ploty1,ploty2,0,0);
					   cpglab(xstr,ystr,psr[0].name);
					   cpgpt(count,x,y,16);
					   if (plotErr==1) cpgerry(count,x,yerr1,yerr2,1);
					   cpgend();
					   cpgbeg(0,"/xs",1,1);
					   cpgask(0); */
			}
			else if (key=='G') { /* Change gridding on graphics device */
			   int nx,ny;
			   printf("Enter number of columns and number of rows (e.g. 2 2) ");
			   scanf("%d %d",&nx,&ny);
			   printf("Enter fontsize (e.g. 1) ");
			   scanf("%f",&fontSize);
			   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
			   cpgend();
			   cpgbeg(0,"/xs",nx,ny);
			   cpgask(0);
			   graphics=-1;
			   if (nx == 1 && ny == 1) graphics=1;
			}
			else if (key==20) /* ctrl-T - set text size */
			{
			   printf("Enter fontsize (e.g. 1) ");
			   scanf("%f",&fontSize);
			   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
			}
			else if (key=='d' || key=='X') 
			{
			   // Look for updating phase jumps
			   if (mouseY > ploty2 && mouseX > plotx1 && mouseX < plotx2 
					 && (xplot==3 || xplot==5) && psr[0].nPhaseJump>0)
			   {
				  float xch,ych;
				  cpgsch(0.5);cpgqcs(4,&xch,&ych);cpgsch(fontSize);
				  if (mouseY < ploty2+1.5*ych)
				  {
					 // Find closest
					 int close=0;
					 double diff = fabs(mouseX-(psr[0].phaseJump[0]-centreEpoch));
					 for (i=1;i<psr[0].nPhaseJump;i++)
					 {
						if (fabs(mouseX-(psr[0].phaseJump[i]-centreEpoch))<diff)
						{
						   close = i;
						   diff = fabs(mouseX-(psr[0].phaseJump[i]-centreEpoch));
						}
					 }
					 if (diff < xch)
					 {
						//printf("close = %d\n",close);
						psr[0].phaseJumpDir[close]++;
						for (i=0;i<psr[0].nobs;i++)
						{
						   k=close;
						   if (psr[0].obsn[i].sat > psr[0].phaseJump[k])
						   {
							  psr[0].obsn[i].residual += (double)1.0/psr[0].param[param_f].val[0];
							  psr[0].obsn[i].prefitResidual += (double)1.0/psr[0].param[param_f].val[0];
						   }
						}
					 }
				  }
			   }
			   if (menu>0 && mouseY > ploty2)
				  checkMenu(psr,mouseX,mouseY,2,fitFlag,setZoomX1,setZoomX2,zoomX1,zoomX2,origStart,origFinish,centreEpoch,menu,xplot,parFile,timFile,argc,argv,&xplot,&yplot,&graphics, highlightID,highlightVal,&highlightNum,aspect,fontType,lineWidth,bkgrdColour,lineColour,&jumpOffset,zoom,&paramOffset);
			   else
				  deletePoint(psr,x,y,id,count,mouseX,mouseY); /* Delete closest point */
			}
			else if (key=='p') changeParameters(psr);   /* Change parameter values */
			else if (key==12) /* Add line to plot */
			{
			   printf("Enter x1 y2 x2 y2 ");
			   scanf("%f %f %f %f",&lx1[nline],&ly1[nline],&lx2[nline],&ly2[nline]);
			   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
			   nline++;
			}
			else if (key=='c') changeFitParameters(psr); /* Change fit parameters   */
			else if (key=='V') //define the user parameters
			{
			   printf("Enter program (e.g. vap or vip) ");
			   scanf("%s",userCMD);
			   printf("Enter parameter understood by '%s' ",userCMD);
			   scanf("%s",userValStr);
			   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
			   userValChange=1;
			}
			else if (key==22) /* Ctrl-v: View models */
			   viewModel*=-1;
			else if (key==18) /* Control-R - select regions */
			{
			   int maxRegions = 5000;
			   float rx1[maxRegions],rx2[maxRegions];
			   int   nregion=0;
			   char key2;
			   float mx2,my2,mx1,my1;
			   char temp[1000];
			   float fx[2],fy[2],fx2[2];
			   FILE *fout,*fin;

			   // Find already recorded regions
			   if (!(fin = fopen("regions.dat","r")))
			   {
				  printf("Creating new regions.dat file\n");
			   }
			   else
			   {
				  while (!feof(fin))
				  {
				    long double lx1,lx2;
				    if (fscanf(fin,"%Lf %Lf %s %s",&lx1,&lx2,temp,temp)==4)
				      {

					fx[0] = (float)(lx1 - centreEpoch);
					fx2[0] = (float)(lx2 - centreEpoch);
					fx[1] = fx[0];
					fx2[1] = fx2[0];
					printf("Loaded %.15Lf %.15Lf,%g %g %.15Lf\n",lx1,lx2,fx[0],fx2[0],centreEpoch);
					fy[0] = ploty1; fy[1] = ploty2;
					cpgsci(2); cpgline(2,fx,fy);
					cpgsls(4); cpgline(2,fx2,fy); cpgsls(1); cpgsci(1);
				      }
				    
				  }
				  fclose(fin);
			   }

			   printf("Use left button to select regions\n");
			   printf("Click right button to finish\n");
			   do {
				  cpgcurs(&mx1,&my1,&key2);
				  if (key2!='X')
				  {
					 cpgband(4,0,mx1,my1,&mx2,&my2,&key2);
					 rx1[nregion] = (float)(mx1+centreEpoch);
					 rx2[nregion] = (float)(mx2+centreEpoch);
					 fx[0] = mx1;
					 fx2[0] = mx2;
					 fx[1] = fx[0];
					 fx2[1] = fx2[0];
					 fy[0] = ploty1; fy[1] = ploty2;
					 cpgsci(2); cpgline(2,fx,fy);
					 cpgsls(4); cpgline(2,fx2,fy); cpgsls(1); cpgsci(1);

					 nregion++;
					 if (nregion > maxRegions)
					 {
						printf("ERROR: maximum number of regions at any given time = %d -- please save and start again\n",maxRegions);
						nregion--;

					 }
				  }
			   } while (key2 != 'X');
			   printf("Goodbye\n");
			   fout = fopen("regions.dat","a");
			   for (i=0;i<nregion;i++)
				  fprintf(fout,"%10.5f %10.5f %s %s\n",rx1[i],rx2[i],parFile[0],timFile[0]);
			   fclose(fout);
			}
			else if (key=='r') {  /* RESET */
			   readParfile(psr,parFile,timFile,1); /* Load the parameters       */
			   readTimfile(psr,timFile,1); /* Load the arrival times    */
			   preProcess(psr,1,argc,argv);
			   callFit(psr,npsr);
			}
			else if (key=='I') /* indicate individual observations */
			{
			   if (indicateObs==0)
				  indicateObs=1;
			   else if (indicateObs==1)
				  indicateObs=2;
			   else
				  indicateObs=0;
			}
			else if (key=='i' || key=='A') 
			{
			   // Look for updating phase jumps
			   if (mouseY > ploty2 && mouseX > plotx1 && mouseX < plotx2 
					 && (xplot==3 || xplot==5) && psr[0].nPhaseJump>0)
			   {
				  float xch,ych;
				  cpgsch(0.5);cpgqcs(4,&xch,&ych);cpgsch(fontSize);
				  if (mouseY < ploty2+1.5*ych)
				  {
					 // Find closest
					 int close=0;
					 double diff = fabs(mouseX-(psr[0].phaseJump[0]-centreEpoch));
					 for (i=1;i<psr[0].nPhaseJump;i++)
					 {
						if (fabs(mouseX-(psr[0].phaseJump[i]-centreEpoch))<diff)
						{
						   close = i;
						   diff = fabs(mouseX-(psr[0].phaseJump[i]-centreEpoch));
						}
					 }
					 if (diff < xch)
					 {
						//printf("close = %d\n",close);
						psr[0].phaseJumpDir[close]--;
						for (i=0;i<psr[0].nobs;i++)
						{
						   k=close;
						   if (psr[0].obsn[i].sat > psr[0].phaseJump[k])
						   {
							  psr[0].obsn[i].residual -= (double)1.0/psr[0].param[param_f].val[0];
							  psr[0].obsn[i].prefitResidual -= (double)1.0/psr[0].param[param_f].val[0];
						   }
						}
					 }
				  }
			   }
			   if ((menu==1 || menu==2) && mouseY > ploty2)
				  checkMenu(psr,mouseX,mouseY,1,fitFlag,setZoomX1,setZoomX2,zoomX1,zoomX2,origStart,origFinish,centreEpoch,menu,xplot,parFile,timFile,argc,argv,&xplot,&yplot,&graphics, highlightID,highlightVal,&highlightNum,aspect,fontType,lineWidth,bkgrdColour,lineColour,&jumpOffset,zoom,&paramOffset);
			   else if (menu==3 && (mouseY > ploty2 || mouseX < plotx1 || mouseY < ploty1))
			   {
				  checkMenu(psr,mouseX,mouseY,1,fitFlag,setZoomX1,setZoomX2,zoomX1,zoomX2,origStart,origFinish,centreEpoch,menu,xplot,parFile,timFile,argc,argv,&xplot,&yplot,&graphics, highlightID,highlightVal,&highlightNum,aspect,fontType,lineWidth,bkgrdColour,lineColour,&jumpOffset,zoom,&paramOffset);
			   }
			   else
				  idPoint(psr,x,y,id,count,mouseX,mouseY); /* Identify closest point */
			}
			else if (key=='s')  /* Start zoom */
			{
			   setZoomX1 = 1;
			   zoom=1;
			   zoomX1 = mouseX;
			}
			else if (key=='f')  /* Finish zoom */
			{
			   setZoomX2 = 1;
			   zoom=1;
			   zoomX2 = mouseX;
			}
			else if (key==6)   /* ctrl-F: remove FITWAVES curve */
			{
			   double px,py,ymod;

			   removeFITWAVES*=-1;
			   for (k=0;k<count;k++)
			   {
				  px=x[k]; 
				  ymod=0.0;
				  for (j=0;j<FITWAVES_n/2;j++)
					 ymod+=FITWAVES_par[2*j]*cos(FITWAVES_omega/FITWAVES_harmonicStep*(j)*px)+
						FITWAVES_par[2*j+1]*sin(FITWAVES_omega/FITWAVES_harmonicStep*(j)*px);
				  py=ymod;
				  psr[0].obsn[k].residual-=removeFITWAVES*py;
				  psr[0].obsn[k].sat-=removeFITWAVES*py/60.0/60.0/24.0;
				  psr[0].obsn[k].bat-=removeFITWAVES*py/60.0/60.0/24.0;
				  /*	    psr[0].obsn[k].prefitResidual-=py; */
			   }
			}
			else if (key=='/') 
			{
			   printf("Parameters reloaded from .par file\n");
			   readParfile(psr,parFile,timFile,1); /* Load the parameters       */
			   callFit(psr,1);
			}
			else if (key=='k')   /* translate timing residuals to dm*/
			{
			   for (i=0;i<count;i++)
			   {
				  if ((setZoomX1==0 || x[i]>zoomX1) &&
						(setZoomX2==0 || x[i]<zoomX2) &&
						(setZoomY1==0 || y[i]>zoomY1) &&
						(setZoomY2==0 || y[i]<zoomY2))
					 printf("%s %s\t%7.5g\t%.5g\t%.5g\t%.5g\n",psr[0].obsn[id[i]].fname,
						   print_longdouble(psr[0].obsn[id[i]].sat).c_str(),x[i],
						   y[i]*DM_CONST*psr[0].obsn[id[i]].freq*psr[0].obsn[id[i]].freq,
						   psr[0].obsn[id[i]].toaErr*DM_CONST*psr[0].obsn[id[i]].freq*psr[0].obsn[id[i]].freq*1.0e-6,
						   (double)psr[0].obsn[id[i]].freq);
			   }	    
			}
			else if (key=='z') /* Define zoom using mouse */
			{
			   cpgband(2,0,mouseX,mouseY,&mouseX2,&mouseY2,&key);
			   setZoomX1 = 1; setZoomX2 = 1; setZoomY1 = 1; setZoomY2 = 1;
			   zoom=1;
			   zoomX1 = TKretMin_f(mouseX,mouseX2);
			   zoomX2 = TKretMax_f(mouseX,mouseX2);
			   zoomY1 = TKretMin_f(mouseY,mouseY2);
			   zoomY2 = TKretMax_f(mouseY,mouseY2);
			}
			else if (key=='u') /* Unzoom */
			{
			   zoom=0;
			   setZoomX1 = 0; setZoomX2 = 0; setZoomY1 = 0; setZoomY2 = 0;
			   psr[0].param[param_start].val[0] = origStart;
			   psr[0].param[param_finish].val[0] = origFinish;
			   psr[0].param[param_start].fitFlag[0] = 0;
			   psr[0].param[param_finish].fitFlag[0] = 0;
			}
			else if (key==13) /* Toggle menu bar */
			{
			   if (menu==-1)     menu=1;
			   else if (menu==1) menu=2;
			   else if (menu==2) 
			   {
				  menu=3;
				  cpgend();
				  cpgbeg(0,"/xs",1,1);
				  cpgpap(0,aspect);
				  cpgscf(fontType);
				  cpgslw(lineWidth);
				  cpgscrn(0,bkgrdColour,&out);
				  cpgscrn(1,lineColour,&out);
				  cpgask(0);
				  graphics=-1;
			   }
			   else if (menu==3) menu=-1;
			}
			else if (key=='Z') //delete points in a region selected by mouse
			{
			   cpgband(2,0,mouseX,mouseY,&mouseX2,&mouseY2,&key);
			   delX1 = TKretMin_f(mouseX,mouseX2);
			   delX2 = TKretMax_f(mouseX,mouseX2);
			   delY1 = TKretMin_f(mouseY,mouseY2);
			   delY2 = TKretMax_f(mouseY,mouseY2);
			   printf("%f %f %f %f\n",delX1,delX2,delY1,delY2);

			   for (i=0;i<count;++i) {
				  //	    printf("Here with %d %f %f %d\n",i,x[i],y[i],id[i]);
				  if ((x[i]>delX1) &&
						(x[i]<delX2) &&
						(y[i]>delY1) &&
						(y[i]<delY2)) {
					 psr[0].obsn[id[i]].deleted = 1;
				  }
			   }
			}

			else printf("Unknown key press %c (%d)\n",key,(int)key);
	}
	if (graphics==1) graphics=2;
	else if (graphics==2)
	{
	   cpgend();
	   cpgbeg(0,"/xs",1,1);
	   cpgpap(0,aspect);
	   cpgscf(fontType);
	   cpgslw(lineWidth);
	   cpgscrn(0,bkgrdColour,&out);
	   cpgscrn(1,lineColour,&out);
	   cpgask(0);
	   graphics=0;
	}
  } while (exitFlag==0);

  cpgend();
}

void changeParameters(pulsar *psr)
{
   int i,k,j,deriv=0;
   char yesno[100];

   printf("Setting parameters\n\n");

   printf("Enter name of parameter to change, or hit enter for a full list.\n");
   gets(yesno);
   if(strcmp(yesno,"")!=0){
	  /* Determine the parameter */
	  for(i=0;i<MAX_PARAMS;i++){
		 if(psr[0].param[i].aSize>1){
			for(j=0;j<psr[0].param[i].aSize;j++){
			   if(strcmp(psr[0].param[i].shortlabel[j],yesno)==0){
				  k=i;
				  deriv=j;
				  break;
			   }
			}
		 }
		 else {
			if(strcmp(*psr[0].param[i].shortlabel,yesno)==0){
			   k=i;
			   break;
			}
		 }
	  }
	  /* Get value to change to */
	  printf("Enter desired value for parameter:\n");
	  printf("%s = %.14lg : change to (don't change: n)\t",yesno,
			(double)psr[0].param[k].val[deriv]);
	  gets(yesno);
	  if(!(yesno[0]=='n'))
		 psr[0].param[k].val[deriv] = parse_longdouble(yesno);
   }
   else if(strcmp(yesno,"")==0){
	  /* Display entire list */
	  for (i=0;i<MAX_PARAMS;i++)
	  {
		 for (k=0;k<psr[0].param[i].aSize;k++)
		 {
			if (psr[0].param[i].paramSet[k]==1) /* If parameter set */
			{
			   printf("%s = %.14f : change [(n = default)]\t",psr[0].param[i].label[k],(double)psr[0].param[i].val[k]);
			   gets(yesno);
			   if (!(yesno[0]=='n' || yesno[0]=='N' || strlen(yesno)==0))
				  psr[0].param[i].val[k] = parse_longdouble(yesno);

			}
		 }
	  }
   }
   callFit(psr,1);
}

void binResiduals(pulsar *psr,int npsr,float *x,float *y,int count,int *id,int *overN,
	  float overX[], float overY[], float overYe[],int xplot,int yplot,
	  float yerr1[],double unitFlag,int plotPhase,double centreEpoch){
   float binSize, x1,x2, pos, overrms=0.0, err=0.0;
   int   j,num,tot=0,noerr=0;
   float errBar[count];

   printf("\nEnter positive width (current units) of each bin.\n");
   printf("Or enter -1 for specification of the number of bins.\n");
   printf("\t(Might be off by one due to rounding.)\n");
   scanf("%f",&binSize);
   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
   if(binSize==-1){
	  printf("Enter number of bins: ");
	  scanf("%f",&binSize);
	  binSize = (x[count-1]-x[0])/binSize;
   }

   /* Weighting has been implemented: the Y-values are divided
	  by the absolute value of their errors (brought to the same units).

	  Since now, 
	  {1}  overY = [sum(y/err)] / [sum(1/err)], 
	  the errorbars of the binned values, will become:
	  overYe = sqrt{N*[1/sum(1/err)]^2}, which is the same as:
	  {2}  overYe = sqrt(N) / sum(1/err).

	  Simple illustration to show this error is not wrong: 
	  suppose all points have the same error - the case where weighting 
	  or not doesn't matter. We expect then:
	  overY = sum(y)/N;
	  overYe = err/sqrt(N);

	  While the formulae above ({1} and {2}), become:
	  overY = [1/err*sum(y)] / [N/err] = sum(y)/N;
	  overYe = sqrt(N)/[N/err] = err/sqrt(N);
	  */

   /* Note y[] and x[] are the unbinned values, overX[] and overY[] are the 
	  binned values. errBar[] are the unbinned error bars on y[], but they
	  are not yet in the same units. overYe[] are the errorbars on overY[].

	  I'm sorry for the confusing names. "over" comes from the fact that 
	  overX, overY and overYe are not the basic things to be plotted, but 
	  will be plotted \emph{over} the original data.
	  */
   for(j=0;j<count;j++){
	  if(plotPhase ==-1)
		 errBar[j] = (y[j]-yerr1[j])/1.0e-6*unitFlag;
	  else 
		 errBar[j] = (y[j]-yerr1[j])/1.0e-6*(double)psr[0].param[param_f].val[0];
	  if(errBar[j] == 0){
		 printf("ERROR: point %d has TOA error 0!!\n");
		 printf("Will not measure error bars for bins\n");
		 noerr = 1;
	  }
   }

   x1  = x[0];
   x2  = x[count-1];
   for (pos=x1;pos<x2;pos+=binSize){
	  num=0;
	  overX[tot] = 0.0;
	  overY[tot] = 0.0;
	  overYe[tot] = 0.0;
	  err = 0.0;
	  for (j=0;j<count;j++){
		 if (x[j] >= pos && x[j] < pos+binSize){
			overX[tot] += (float) x[j];

			/*if(noerr==0)
			  overYe[tot] += powf(errBar[j]*1.0e-6/(float)unitFlag,2.0);
			  else if(noerr == 1)
			  overYe[tot] += 0.0;*/

			if(noerr==0)
			   overY[tot] += (float)(y[j])/fabs(errBar[j]*1.0e-6/
					 (float)unitFlag);
			else if(noerr==1)
			   overY[tot] += (float)y[j];

			if(noerr==0) 
			   err += 1.0/fabs(errBar[j]*1.0e-6/(float)unitFlag);
			else if(noerr==1)
			   err += 1.0;
			num++;
		 }	    
	  }
	  if (num>0){
		 overX[tot] = overX[tot]/(float)num;
		 overY[tot] = overY[tot]/err;
		 overrms += powf(overY[tot],2.0);
		 //overYe[tot] = sqrt(overYe[tot])/(float)num;
		 overYe[tot] = sqrt((float)num)/err;
		 printf("bin %d %.4f %g %g %g\n",tot+1,overX[tot]+centreEpoch,overX[tot],overY[tot],overYe[tot]);
		 err = 0.0;
		 tot++;
	  }
   }  

   overrms = sqrt(overrms/tot);
   printf("Total number of bins = %d; rms: %lg\n",tot, overrms);
   *overN = tot;
}

void changeFitParameters(pulsar *psr)
{
   int i,k=0;
   char yesno[100];

   /* CLEAR STDIN */
   fflush(stdin);

   printf("Setting fitting parameters\n\n");

   printf("Enter name of parameter to switch fitting for, or hit enter for a full list.\n");
   gets(yesno);
   if(strcmp(yesno,"")!=0){
	  /* Determine the parameter */
	  for(i=0;i<MAX_PARAMS;i++){
		 if(strcmp(*psr[0].param[i].shortlabel,yesno)==0){
			k=i;
			break;
		 }
	  }

	  /* Turn fitting on/off */
	  if(psr[0].param[k].fitFlag[0]==1){
		 printf("Turning fitting off for %s.\n",yesno);
		 psr[0].param[k].fitFlag[0]=0;
	  }
	  else if(psr[0].param[k].fitFlag[0]==0){
		 printf("Turning fitting on for %s.\n",yesno);
		 psr[0].param[k].fitFlag[0]=1;
	  }
   }
   else if(strcmp(yesno,"")==0){
	  /* Display entire list */
	  for (i=0;i<MAX_PARAMS;i++)
	  {
		 for (k=0;k<psr[0].param[i].aSize;k++)
		 {
			if (psr[0].param[i].paramSet[k]==1) /* If parameter set */
			{
			   if (psr[0].param[i].fitFlag[k]==1)
			   {
				  printf("%s fitting (turn off fitting 'y', otherwise hit enter)\t",psr[0].param[i].label[k]);
				  gets(yesno);
				  if (yesno[0]=='y' || yesno[0]=='Y')
					 psr[0].param[i].fitFlag[k]=0;
			   }
			   else if (psr[0].param[i].fitFlag[k]==0)
			   {
				  printf("%s not fitting (turn on fitting 'y', otherwise hit enter)\t",psr[0].param[i].label[k]);
				  gets(yesno);
				  if (yesno[0]=='y' || yesno[0]=='Y')
					 psr[0].param[i].fitFlag[k]=1;
			   }
			}
		 }
	  }
   }
   callFit(psr,1);
}


void sort(float *x,float *y,float *yerr1,float *yerr2,float *freq,int *id,int count)
{
   int i,changed=0;
   float temp;
   int itemp;

   do {
	  changed=0;
	  for (i=0;i<count-1;i++)
	  {
		 if (x[i] > x[i+1])
		 {
			temp = x[i]; x[i] = x[i+1]; x[i+1] = temp;
			temp = y[i]; y[i] = y[i+1]; y[i+1] = temp;
			temp = yerr1[i]; yerr1[i] = yerr1[i+1]; yerr1[i+1] = temp;
			temp = yerr2[i]; yerr2[i] = yerr2[i+1]; yerr2[i+1] = temp;
			temp = freq[i]; freq[i] = freq[i+1]; freq[i+1] = temp;
			itemp = id[i]; id[i] = id[i+1]; id[i+1] = itemp;
			changed=1;
		 }
	  }

   }while (changed==1);

}

float deletePoint(pulsar *psr,float *x,float *y,int *id,int count,float mouseX,float mouseY)
{
   int i,iclosest;
   float closest,x1,x2,x3,x4,y1,y2,y3,y4,xscale,yscale,xpos,ypos;

   cpgqvp(3,&x1,&x2,&y1,&y2);
   cpgqwin(&x3,&x4,&y3,&y4);
   xscale = (x2-x1)/(x4-x3);
   yscale = (y2-y1)/(y4-y3);
   mouseX = (mouseX-x3)*xscale;
   mouseY = (mouseY-y3)*yscale;
   iclosest=-1;
   for (i=0;i<count;i++)
   {
	  xpos = (x[i]-x3)*xscale;
	  ypos = (y[i]-y3)*yscale;
	  if (iclosest==-1)
	  {
		 iclosest=id[i];
		 closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
	  }
	  else if (pow(xpos-mouseX,2)+pow(ypos-mouseY,2)<closest)
	  {
		 iclosest=id[i];
		 closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
	  }
   }
   psr[0].obsn[iclosest].deleted=1;
   return 0.0;
}



int idPoint(pulsar *psr,float *x,float *y,int *id,int count,float mouseX,float mouseY)
{
   int i,iclosest,l;
   float closest,x1,x2,x3,x4,y1,y2,y3,y4,xscale,yscale,xpos,ypos;

   cpgqvp(3,&x1,&x2,&y1,&y2);
   cpgqwin(&x3,&x4,&y3,&y4);
   xscale = (x2-x1)/(x4-x3);
   yscale = (y2-y1)/(y4-y3);
   mouseX = (mouseX-x3)*xscale;
   mouseY = (mouseY-y3)*yscale;
   iclosest=-1;
   for (i=0;i<count;i++)
   {
	  xpos = (x[i]-x3)*xscale;
	  ypos = (y[i]-y3)*yscale;
	  if (iclosest==-1)
	  {
		 iclosest=id[i];
		 closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
	  }
	  else if (pow(xpos-mouseX,2)+pow(ypos-mouseY,2)<closest)
	  {
		 iclosest=id[i];
		 closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
	  }
   }

   printf("---------------------------------------------------\n");
   printf("Closest point has TOA number %d (starting from 1)\n",iclosest+1);
   printf("SAT = %s\n",psr[0].obsn[iclosest].fname);
   /*  printf("SAT = %s\n",print_longdouble(psr[0].obsn[iclosest].sat).c_str());
	   printf("BAT = %s\n",print_longdouble(psr[0].obsn[iclosest].bat).c_str()); */
   printf("SAT = %.14Lf, TOA error = %.3f (us)\n",psr[0].obsn[iclosest].sat,(double)psr[0].obsn[iclosest].toaErr);
   printf("BAT = %.14Lf\n",psr[0].obsn[iclosest].bat);
   printf("Pre-fit residual = %lg\n",(double)psr[0].obsn[iclosest].prefitResidual);
   printf("Post-fit residual = %lg\n",(double)psr[0].obsn[iclosest].residual);
   printf("Observing frequency = %f\n",psr[0].obsn[iclosest].freq);
   printf("Flags = ");
   for (i=0;i<psr[0].obsn[iclosest].nFlags;i++)
	  printf(" \"%s\" %s ",psr[0].obsn[iclosest].flagID[i],psr[0].obsn[iclosest].flagVal[i]);
   printf("\n");

	  printf(psr[0].obsn[iclosest].telID);
   for (l=0;l<psr[0].obsn[iclosest].nclock_correction;l++){
	  printf(" -> ");
	  printf(psr[0].obsn[iclosest].correctionsTT[l].corrects_to);
   }	
   printf("\n");
   printf("---------------------------------------------------\n");

   /* Remove period from this point */
   /*  psr[0].obsn[iclosest].sat-=1.0/psr[0].param[param_f0].val;
	   psr[0].obsn[iclosest].bat-=1.0/psr[0].param[param_f0].val;
	   psr[0].obsn[iclosest].prefitResidual-=1.0/psr[0].param[param_f0].val;
	   psr[0].obsn[iclosest].residual-=1.0/psr[0].param[param_f0].val; */
   return iclosest;
}




float findMaxY(float *y,float *x,int count,float xmin,float xmax)
{
   int i,t=1;
   float max;

   for (i=0;i<count;i++)
   {
	  if (x[i] > xmin && x[i] < xmax)
	  {
		 if (t==1)
		 {
			t=0;
			max = y[i];
		 }
		 else if (max < y[i])
			max = y[i];	    
	  }

   }
   return max;
}

float findMinY(float *y,float *x,int count,float xmin,float xmax)
{
   int i,t=1;
   float min;

   for (i=0;i<count;i++)
   {
	  if (x[i] > xmin && x[i] < xmax)
	  {
		 if (t==1)
		 {
			t=0;
			min = y[i];
		 }
		 else if (min > y[i])
			min = y[i];	    
	  }

   }
   return min;
}

float findMean(float *x,pulsar *psr,int i1,int i2)
{
   int i,count;
   float mean;

   count=0;
   mean=0.0;
   for (i=i1;i<i2;i++)
   {
	  mean+=x[i];
	  count++;
   }
   mean/=count;
   return mean;
}

double findMeanD(float *x,pulsar *psr,int i1,int i2)
{
   int i,count;
   float mean;

   count=0;
   mean=0.0;
   for (i=i1;i<i2;i++)
   {
	  mean+=x[i];
	  count++;
   }
   mean/=count;
   return mean;
}

double fortranMod(double a,double p)
{
   double ret;

   ret = a - (int)(a/p)*p;
   return ret;
}

void overPlotN(int overN,float overX[], float overY[],float overYe[]){
   cpgsci(2);
   cpgpt(overN,overX,overY,3);
   cpgerrb(6,overN,overX,overY,overYe,1.0);
}

void overPlotShapiro(pulsar *psr,float offset,longdouble centreEpoch)
{
   float x[MAX_OBSN],y[MAX_OBSN];
   int i;
   char temp[100];

   for (i=0;i<psr[0].nobs;i++)
   {
	  x[i] = (float)(double)(psr[0].obsn[i].bat-centreEpoch);
	  y[i] = (float)(double)psr[0].obsn[i].shapiroDelaySun*1e6-offset; /* For mus plotting */
	  printf("Have %d %f %f\n",i,x[i],y[i]);
   }
   cpgsls(2);
   cpgline(psr[0].nobs,x,y); 
   cpgsls(1);
}

void checkMenu(pulsar *psr,float mx,float my,int button,int fitFlag,int setZoomX1,int setZoomX2,
	  float zoomX1,float zoomX2,longdouble origStart,longdouble origFinish,
	  longdouble centreEpoch,int menu,int plotx,char parFile[][MAX_FILELEN],
	  char timFile[][MAX_FILELEN],int argc,char *argv[],int *xplot,int *yplot,int *graphics,
	  char highlightID[100][100],char  highlightVal[100][100],int *highlightNum,float aspect,int fontType,int lineWidth,char *bkgrdColour,char *lineColour,int *jumpOffset,int zoom,int *paramOffset) 
{
   float x1,x2,y1,y2,x3,x4,y3,y4,x7,y7,xscale,yscale,xscale2,x0;
   float x5,x6,y5,y6;
   int mouseX,mouseY;

   cpgqvp(3,&x1,&x2,&y1,&y2);
   cpgqwin(&x3,&x4,&y3,&y4);
   cpgqvsz(3,&x5,&x6,&y5,&y6);
   xscale = (x2-x1)/(x4-x3);
   yscale = (y2-y1)/(y4-y3);

   x7 = (mx-x3)*xscale+x1; 
   y7 = (my-y3)*yscale+y1; 
   mouseX = (int)(x7/(x6-x5)*10.0);
   mouseY = (int)((1-(y7/(y6-y5)))/0.2*5.0); // 0.2 because the top menu bar has height of 0.2
   if (((mx-x3)*xscale+x1)/(x2+x1)*10 < mouseX+0.2 || mouseY==3)
   {
	  if (menu==1)
	  {
		 if (mouseX==0 && mouseY==0) swapFit(psr,param_raj,0,button);
		 if (mouseX==1 && mouseY==0) swapFit(psr,param_decj,0,button);
		 if (mouseX==2 && mouseY==0) swapFit(psr,param_f,0,button);
		 if (mouseX==3 && mouseY==0) swapFit(psr,param_f,1,button);
		 if (mouseX==4 && mouseY==0) swapFit(psr,param_f,2,button);
		 if (mouseX==5 && mouseY==0) swapFit(psr,param_dm,0,button);
		 if (mouseX==6 && mouseY==0) swapFit(psr,param_pmra,0,button);
		 if (mouseX==7 && mouseY==0) swapFit(psr,param_pmdec,0,button);
		 if (mouseX==8 && mouseY==0) swapFit(psr,param_px,0,button);

		 if (mouseX==0 && mouseY==1) swapFit(psr,param_pb,0,button);
		 if (mouseX==1 && mouseY==1) swapFit(psr,param_a1,0,button);
		 if (strcmp(psr->binaryModel,"ELL1")==0
			   || (strcmp(psr->binaryModel,"T2")==0 && psr->param[param_eps1].paramSet[0]==1))
		 {
			if (mouseX==2 && mouseY==1) swapFit(psr,param_tasc,0,button);
			if (mouseX==3 && mouseY==1) swapFit(psr,param_eps1,0,button);
			if (mouseX==4 && mouseY==1) swapFit(psr,param_eps2,0,button);
			if (mouseX==5 && mouseY==1) swapFit(psr,param_pbdot,0,button);
			if (mouseX==7 && mouseY==1) swapFit(psr,param_kom,0,button);
			if (mouseX==8 && mouseY==1) swapFit(psr,param_kin,0,button);
		 }
		 else if (strcmp(psr->binaryModel,"T2")==0)
		 {
			if (mouseX==2 && mouseY==1) swapFit(psr,param_t0,0,button);
			if (mouseX==3 && mouseY==1) swapFit(psr,param_ecc,0,button);
			if (mouseX==4 && mouseY==1) swapFit(psr,param_om,0,button);
			if (mouseX==5 && mouseY==1) swapFit(psr,param_omdot,0,button);
			if (mouseX==6 && mouseY==1) swapFit(psr,param_pbdot,0,button);
			if (mouseX==7 && mouseY==1) swapFit(psr,param_kom,0,button);
			if (mouseX==8 && mouseY==1) swapFit(psr,param_kin,0,button);
			if (mouseX==9 && mouseY==1) swapFit(psr,param_sini,0,button);
			if (mouseX==9 && mouseY==0) swapFit(psr,param_m2,0,button);
		 }
		 else
		 {
			if (mouseX==2 && mouseY==1) swapFit(psr,param_t0,0,button);
			if (mouseX==3 && mouseY==1) swapFit(psr,param_ecc,0,button);
			if (mouseX==4 && mouseY==1) swapFit(psr,param_om,0,button);
			if (mouseX==5 && mouseY==1) swapFit(psr,param_omdot,0,button);
			if (mouseX==6 && mouseY==1) swapFit(psr,param_pbdot,0,button);
		 }
	  }
	  else if (menu==2)
	  {
		 if (mouseX==0 && mouseY==0) swapFit(psr,param_f,0,button);
		 if (mouseX==1 && mouseY==0) swapFit(psr,param_f,1,button);
		 if (mouseX==2 && mouseY==0) swapFit(psr,param_f,2,button);
		 if (mouseX==3 && mouseY==0) swapFit(psr,param_f,3,button);
		 if (mouseX==4 && mouseY==0) swapFit(psr,param_f,4,button);
		 if (mouseX==5 && mouseY==0) swapFit(psr,param_f,5,button);
		 if (mouseX==6 && mouseY==0) swapFit(psr,param_f,6,button);
		 if (mouseX==7 && mouseY==0) swapFit(psr,param_f,7,button);
		 if (mouseX==8 && mouseY==0) swapFit(psr,param_f,8,button);
	  }
	  else if (menu==3)
	  {
		 int count=0;
		 int xpos,ypos;
		 int i,j;
		 xpos=0;
		 ypos=0;
		 //	  if (*paramOffset > 0)
		 {
			if (mouseX == 9 && mouseY == 0 && *paramOffset > 0)		
			   (*paramOffset)--;
			else if (mouseX == 9 && mouseY == 2)
			   (*paramOffset)++;
			printf("Setting offset in menu to be %d\n",*paramOffset);
		 }
		 for (i=0;i<MAX_PARAMS;i++)
		 {
			for (j=0;j<psr[0].param[i].aSize;j++)
			{
			   if (psr[0].param[i].paramSet[j]==1)
			   {
				  if (strcmp(psr[0].param[i].shortlabel[j],"DMEPOCH")!=0 &&
						strcmp(psr[0].param[i].shortlabel[j],"PEPOCH")!=0 &&
						strcmp(psr[0].param[i].shortlabel[j],"START")!=0 &&
						strcmp(psr[0].param[i].shortlabel[j],"FINISH")!=0 &&
						strcmp(psr[0].param[i].shortlabel[j],"TRACK")!=0 &&
						strcmp(psr[0].param[i].shortlabel[j],"TZRMJD")!=0 &&
						strcmp(psr[0].param[i].shortlabel[j],"TZRFRQ")!=0 &&
						strcmp(psr[0].param[i].shortlabel[j],"TRES")!=0 &&
						strcmp(psr[0].param[i].shortlabel[j],"EPHVER")!=0 &&
						strcmp(psr[0].param[i].shortlabel[j],"POSEPOCH")!=0)
				  {
					 if (mouseX==xpos && mouseY==ypos-(*paramOffset))
						swapFit(psr,i,j,button);
					 xpos++;
					 if (xpos > 8)
					 {
						xpos=0;
						ypos++;
					 }
				  }
			   }
			}
		 }

	  }
	  if (mouseX==0 && mouseY==3) 
	  {
		 // Doing the refit

		 reFit(fitFlag,setZoomX1,setZoomX2,zoomX1,zoomX2,origStart,origFinish,centreEpoch,psr,1,plotx,dcmFile,covarFuncFile,zoom);

		 /*     callFit(psr,1); */ /* MUST FIX FOR ZOOM MODES */
	  }
	  if (mouseX==1 && mouseY==3) textOutput(psr,1,0,0,0,1,""); 
	  if (mouseX==2 && mouseY==3) newTim(psr);
	  if (mouseX==3 && mouseY==3)  /* Reload */
	  {
		 int k,l;
		 psr[0].nJumps=0;
		 for (k=0;k<MAX_PARAMS;k++)
		 {
			psr[0].param[k].nLinkTo=0;
			psr[0].param[k].nLinkFrom=0;
		 }
		 initialise(psr, 1);
		 readParfile(psr,parFile,timFile,1); /* Load the parameters       */
		 readTimfile(psr,timFile,1); /* Load the arrival times    */
		 preProcess(psr,1,argc,argv);
		 callFit(psr,1);	  
	  }
	  /*      if (mouseX==4 && mouseY==3)
			  {
			  int i,k;
			  for (i=0;i<MAX_PARAMS;i++)
			  {
			  for (k=0;k<psr[0].param[i].aSize;k++)
			  {
			  if (psr[0].param[i].fitFlag[k]==1)
			  psr[0].param[i].val[k] = psr[0].param[i].prefit[k];
			  }
			  }
			  for (i=0;i<psr[0].nobs;i++)
			  psr[0].obsn[i].residual = psr[0].obsn[i].prefitResidual;
			  }*/
   }
   if (menu==3)
	  checkMenu3(psr,mx,my,1,fitFlag,setZoomX1,setZoomX2,zoomX1,zoomX2,origStart,origFinish,centreEpoch,menu,plotx,parFile,timFile,argc,argv,xplot,yplot,graphics, highlightID,highlightVal,highlightNum,aspect,fontType,lineWidth,bkgrdColour,lineColour,jumpOffset);

}

void checkMenu3(pulsar *psr,float mx,float my,int button,int fitFlag,int setZoomX1,int setZoomX2,
	  float zoomX1,float zoomX2,longdouble origStart,longdouble origFinish,
	  longdouble centreEpoch,int menu,int plotx,char parFile[][MAX_FILELEN],
	  char timFile[][MAX_FILELEN],int argc,char *argv[],int *xplot,int *yplot,int *graphics,
	  char highlightID[100][100],char  highlightVal[100][100],int *highlightNum,float aspect,int fontType,int lineWidth,char *bkgrdColour,char *lineColour,int *jumpOffset) 
{
   float x1,x2,y1,y2,x3,x4,y3,y4,x7,y7,xscale,yscale,xscale2,x0;
   float x5,x6,y5,y6;
   int mouseX,mouseY;

   cpgqvp(3,&x1,&x2,&y1,&y2);
   cpgqwin(&x3,&x4,&y3,&y4);
   cpgqvsz(3,&x5,&x6,&y5,&y6);
   xscale = (x2-x1)/(x4-x3);
   yscale = (y2-y1)/(y4-y3);

   x7 = (mx-x3)*xscale+x1; 
   y7 = (my-y3)*yscale+y1; 
   mouseX = (int)((x7/(x6-x5)/0.3-0.5)/0.1);
   mouseY = (int)(((1-(y7/(y6-y5)))-0.2)*2.0/0.6*10.0); // 0.2 because the top menu bar has height of 0.2
   //  printf("MENU3 mouseX = %d, mouseY = %d [%g]\n",mouseX,mouseY,((1-(y7/(y6-y5)))-0.2)*2.0/0.6*10.0);
   if (mouseX==0)
   {
	  if (mouseY==0) *xplot=1;
	  else if (mouseY==1)  *xplot=2;
	  else if (mouseY==2)  *xplot=3;
	  else if (mouseY==3)  *xplot=4;
	  else if (mouseY==4)  *xplot=5;
	  else if (mouseY==5)  *xplot=6;
	  else if (mouseY==6)  *xplot=7;
	  else if (mouseY==7)  *xplot=8;
	  else if (mouseY==8)  *xplot=10;
	  else if (mouseY==9)  *xplot=11;
	  else if (mouseY==10) *xplot=12;
	  else if (mouseY==11) *xplot=13;
	  else if (mouseY==12) *xplot=14;
	  else if (mouseY==13) *xplot=15;
	  else if (mouseY==14) *xplot=16;
	  else if (mouseY==15) *xplot=17;

   }
   else if (mouseX==1)
   {
	  if (mouseY==0) *yplot=1;
	  else if (mouseY==1)  *yplot=2;
	  else if (mouseY==2)  *yplot=3;
	  else if (mouseY==3)  *yplot=4;
	  else if (mouseY==4)  *yplot=5;
	  else if (mouseY==5)  *yplot=6;
	  else if (mouseY==6)  *yplot=7;
	  else if (mouseY==7)  *yplot=8;
	  else if (mouseY==8)  *yplot=10;
	  else if (mouseY==9)  *yplot=11;
	  else if (mouseY==10) *yplot=12;
	  else if (mouseY==11) *yplot=13;
	  else if (mouseY==12) *yplot=14;
	  else if (mouseY==13) *yplot=15;
	  else if (mouseY==14 &&  psr[0].TNRedAmp != 0 && psr[0].TNRedGam != 0) *yplot=16;
	  else if (mouseY==15 && psr[0].TNDMAmp != 0 && psr[0].TNDMGam != 0) *yplot=17;

   }

   // Now check the bottom menu
   mouseX = (int)(x7/(x6-x5)*10.0);
   mouseY = (int)(((1-(y7/(y6-y5)))-0.8)/0.2*5); // 0.2 because the top menu bar has height of 0.2
   //  printf("MENU3 mouseX = %d, mouseY = %d [%g]\n",mouseX,mouseY,((1-(y7/(y6-y5)))-0.2)*2.0/0.6*10.0);
   if (mouseX==0 && mouseY==0)
   {
	  printf("Goodbye\n");
	  exit(1);
   }
   else if (mouseX==1 && mouseY==0)
   {      
	  int i;
	  int k;
	  for (i=0;i<MAX_PARAMS;i++)
	  {
		 for (k=0;k<psr[0].param[i].aSize;k++)	    
			psr[0].param[i].fitFlag[k]=0;
	  }
	  for (i=1;i<psr[0].nJumps;i++)
		 psr[0].fitJump[i]=0;
   }
   else if (mouseX==2 && mouseY==0)
   {
	  float mouse_x1,mouse_x2,mouse_y1,mouse_y2;
	  char key;
	  cpgband(7,0,0,0,&mouse_x1,&mouse_y1,&key);
	  cpgband(1,0,mouse_x1,mouse_y1,&mouse_x2,&mouse_y2,&key);	
	  printf("Measure x = %.3g days\t y = %.3g\t gradient = %.3g\n",mouse_x2-mouse_x1,mouse_y2-mouse_y1,(mouse_y2-mouse_y1)/(mouse_x2-mouse_x1));
   }
   else if (mouseX==3 && mouseY==0)
	  help();
   else if (mouseX==0 && mouseY==1)
   {
	  int i;
	  for (i=0;i<psr[0].nobs;i++)
	  {
		 if (psr[0].obsn[i].deleted==0)
			printf("%12.5f %13.8g %12.5g\n",(double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]),(double)psr[0].obsn[i].residual,(double)psr[0].obsn[i].toaErr*1e-6);
	  }
   }
   else if (mouseX==1 && mouseY==1)
   {
	  int out;
	  cpgend();
	  cpgbeg(0,"/cps",1,1);
	  cpgpap(0,aspect);
	  cpgscf(fontType);
	  cpgslw(lineWidth);
	  //      cpgscrn(0,bkgrdColour,&out);
	  //      cpgscrn(1,lineColour,&out);

	  *graphics=1;
   }
   else if (mouseX==2 && mouseY==1)
   {
	  printf("FlagID = ");  scanf("%s",highlightID[*highlightNum]);
	  printf("FlagVal = "); scanf("%s",highlightVal[(*highlightNum)++]);
	  getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
   }

   /* Now check jumps */
   if (psr[0].nJumps > 0)
   {
	  int yj,xj;
	  int c=1;

	  if (psr[0].nJumps > 15)
	  {
		 printf("%g %g %d %d\n",mx,my,mouseX,mouseY);
		 if (mouseX == 9 && mouseY == 2)
		 {
			if (*jumpOffset>=5)
			   *jumpOffset-=5;
			c=0;
		 }
		 else if (mouseX == 9 && mouseY == 4)
		 {
			if (*jumpOffset+15 <= psr[0].nJumps)
			   *jumpOffset+=5;
			c=0;
		 }
	  }
	  if (c==1 && mouseX < 10)
	  {
		 xj = (int)(mouseX/2);
		 yj = (int)(mouseY-2);
		 if (psr[0].fitJump[xj+(yj*5)+1+(*jumpOffset)]==1) 
		 {
			psr[0].fitJump[xj+(yj*5)+1+(*jumpOffset)]=0;
			psr[0].jumpValErr[xj+(yj*5)+1+(*jumpOffset)] = 0.0;
		 }
		 else psr[0].fitJump[xj+(yj*5)+1+(*jumpOffset)]=1;
		 //      printf("%d %d %d %d\n",mouseX,mouseY,xj,xj+(yj*5)+1);
	  }
   }
}

void swapFit(pulsar *psr,int par,int k,int button)
{
   if (button==1) /* Left mouse button */
   {
	  if (psr->param[par].fitFlag[k]==0)
	  {
		 psr->param[par].fitFlag[k]=1;
		 if (psr->param[par].paramSet[k]==0)
		 {
			psr->param[par].paramSet[k]=1;
			psr->param[par].prefit[k]=0.0;
			psr->param[par].val[k] = 0.0;
			psr->param[par].err[k] = 0.0;
		 }
	  }
	  else
		 psr->param[par].fitFlag[k]=0;
   }
   else
   {
	  if (psr->param[par].paramSet[k]==0)
	  {
		 psr->param[par].paramSet[k]=1;
		 psr->param[par].val[k]=0.0;
		 psr->param[par].err[k]=0.0;
	  }

	  printf("Current value for %s = %.14lg\n",psr->param[par].label[k],(double)psr->param[par].val[k]);
	  printf("Please enter new value "); 
	  char inpstr[1024];
	  scanf("%s",inpstr);
	  getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
	  psr->param[par].val[k]  = parse_longdouble(inpstr);

   }
}

void drawMenu(pulsar *psr, float plotx1,float plotx2,float ploty1,float ploty2,int menu,int paramOffset)
{
   float x1,y1,x2,y2,xscale,yscale,x3,x4,y3,y4;
   float x0,y0,xscale2,yscale2;
   float mx,my;
   float xpos,ypos;
   char key;

   /*  cpgrect(plotx1,plotx1+2,ploty1,10); */
   /*  cpgqvp(3,&x1,&x2,&y1,&y2);
	   cpgqwin(&x3,&x4,&y3,&y4);
	   xscale = (x4-x3)/(x2-x1);
	   yscale = (y2-y1)/(y4-y3); 


	   cpgrect(plotx1,plotx2,ploty1,10);

	   printf("Have %f %f %f %f %f %f %f\n",x1,x2,x3,x4,xscale,plotx1,plotx2);
	   printf("Scales = %f %f\n",xscale2,yscale2);
	   cpgcurs(&mx,&my,&key);
	   printf("Have %f %f\n",mx,my); */
   if (menu==1)
   {
	  drawOption(0,0.9,"RAJ",psr->param[param_raj].fitFlag[0]);
	  drawOption(0.1,0.90,"DECJ",psr->param[param_decj].fitFlag[0]);
	  drawOption(0.2,0.90,"F0",psr->param[param_f].fitFlag[0]);
	  drawOption(0.3,0.90,"F1",psr->param[param_f].fitFlag[1]);
	  drawOption(0.4,0.90,"F2",psr->param[param_f].fitFlag[2]);
	  drawOption(0.5,0.90,"DM",psr->param[param_dm].fitFlag[0]);
	  drawOption(0.6,0.90,"PMRA",psr->param[param_pmra].fitFlag[0]);
	  drawOption(0.7,0.90,"PMDEC",psr->param[param_pmdec].fitFlag[0]);
	  drawOption(0.8,0.90,"PX",psr->param[param_px].fitFlag[0]);
	  if (strcmp(psr->binaryModel,"ELL1")==0 
			|| (strcmp(psr->binaryModel,"T2")==0 && psr->param[param_eps1].paramSet[0]==1))
	  {
		 drawOption(0.2,0.70,"TASC",psr->param[param_tasc].fitFlag[0]);
		 drawOption(0.3,0.70,"EPS1",psr->param[param_eps1].fitFlag[0]);  
		 drawOption(0.4,0.70,"EPS2",psr->param[param_eps2].fitFlag[0]);
		 drawOption(0.5,0.70,"PBDOT",psr->param[param_pbdot].fitFlag[0]);    
		 drawOption(0.7,0.70,"KOM",psr->param[param_kom].fitFlag[0]);
		 drawOption(0.8,0.70,"KIN",psr->param[param_kin].fitFlag[0]);
		 drawOption(0.0,0.70,"PB",psr->param[param_pb].fitFlag[0]);
		 drawOption(0.1,0.70,"A1",psr->param[param_a1].fitFlag[0]);
	  }
	  else if (strcmp(psr->binaryModel,"T2")==0)
	  {
		 drawOption(0.2,0.70,"T0",psr->param[param_t0].fitFlag[0]);
		 drawOption(0.3,0.70,"ECC",psr->param[param_ecc].fitFlag[0]);  
		 drawOption(0.4,0.70,"OM",psr->param[param_om].fitFlag[0]);
		 drawOption(0.5,0.70,"OMDOT",psr->param[param_omdot].fitFlag[0]);
		 drawOption(0.6,0.70,"PBDOT",psr->param[param_pbdot].fitFlag[0]);
		 drawOption(0.7,0.70,"KOM",psr->param[param_kom].fitFlag[0]);
		 drawOption(0.8,0.70,"KIN",psr->param[param_kin].fitFlag[0]);
		 drawOption(0.9,0.90,"M2",psr->param[param_m2].fitFlag[0]);
		 drawOption(0.9,0.70,"SINI",psr->param[param_sini].fitFlag[0]);
		 drawOption(0.0,0.70,"PB",psr->param[param_pb].fitFlag[0]);
		 drawOption(0.1,0.70,"A1",psr->param[param_a1].fitFlag[0]);
	  }
	  else if (psr->param[param_pb].paramSet[0] == 1)
	  {
		 drawOption(0.2,0.70,"T0",psr->param[param_t0].fitFlag[0]);
		 drawOption(0.3,0.70,"ECC",psr->param[param_ecc].fitFlag[0]);  
		 drawOption(0.4,0.70,"OM",psr->param[param_om].fitFlag[0]);
		 drawOption(0.5,0.70,"OMDOT",psr->param[param_omdot].fitFlag[0]);
		 drawOption(0.6,0.70,"PBDOT",psr->param[param_pbdot].fitFlag[0]);
		 drawOption(0.0,0.70,"PB",psr->param[param_pb].fitFlag[0]);
		 drawOption(0.1,0.70,"A1",psr->param[param_a1].fitFlag[0]);
	  }
	  else 
		 if(debugFlag==1) printf("No binary model\n");
   }
   else if (menu==2)
   {
	  drawOption(0.0,0.90,"F0",psr->param[param_f].fitFlag[0]);
	  drawOption(0.1,0.90,"F1",psr->param[param_f].fitFlag[1]);
	  drawOption(0.2,0.90,"F2",psr->param[param_f].fitFlag[2]);
	  drawOption(0.3,0.90,"F3",psr->param[param_f].fitFlag[3]);
	  drawOption(0.4,0.90,"F4",psr->param[param_f].fitFlag[4]);
	  drawOption(0.5,0.90,"F5",psr->param[param_f].fitFlag[5]);
	  drawOption(0.6,0.90,"F6",psr->param[param_f].fitFlag[6]);
	  drawOption(0.7,0.90,"F7",psr->param[param_f].fitFlag[7]);
	  drawOption(0.8,0.90,"F8",psr->param[param_f].fitFlag[8]);
   }
   else if (menu==3)
   {
	  int count=0;
	  float xpos,ypos;
	  float xposv,yposv;
	  int i,j;
	  xpos=0.0;
	  ypos=0.9;
	  for (i=0;i<MAX_PARAMS;i++)
	  {
		 for (j=0;j<psr[0].param[i].aSize;j++)
		 {
			if (psr[0].param[i].paramSet[j]==1)
			{
			   if (strcmp(psr[0].param[i].shortlabel[j],"DMEPOCH")!=0 &&
					 strcmp(psr[0].param[i].shortlabel[j],"PEPOCH")!=0 &&
					 strcmp(psr[0].param[i].shortlabel[j],"START")!=0 &&
					 strcmp(psr[0].param[i].shortlabel[j],"FINISH")!=0 &&
					 strcmp(psr[0].param[i].shortlabel[j],"TRACK")!=0 &&
					 strcmp(psr[0].param[i].shortlabel[j],"TZRMJD")!=0 &&
					 strcmp(psr[0].param[i].shortlabel[j],"TZRFRQ")!=0 &&
					 strcmp(psr[0].param[i].shortlabel[j],"TRES")!=0 &&
					 strcmp(psr[0].param[i].shortlabel[j],"EPHVER")!=0 &&
					 strcmp(psr[0].param[i].shortlabel[j],"POSEPOCH")!=0)
			   {
				  if (ypos+paramOffset*0.2 <= 0.9 && ypos+paramOffset*0.2 >= 0.5)
					 drawOption(xpos,ypos+paramOffset*0.2,psr[0].param[i].shortlabel[j],psr->param[i].fitFlag[j]);
				  xpos+=0.1;
				  if (xpos > 0.9)
				  {
					 xpos=0.0;
					 ypos-=0.2;
					 if (paramOffset > 0)
					 {
						xposv = 0.985; yposv = 0.9;
						cpgpt(1,&xposv,&yposv,30);			      
					 }

					 if (ypos+paramOffset*0.2 < 0.5) // Don't plot any more
					 {
						i = MAX_PARAMS;
						xposv = 0.985; yposv = 0.5;
						cpgpt(1,&xposv,&yposv,31);
					 }
				  }
			   }
			}
		 }
	  }
   }
   cpgrect(0,0.08,0.28,0.3+0.1); cpgsci(0); cpgtext(0,0.30,"RE-FIT"); cpgsci(1);
   cpgrect(0.1,0.18,0.28,0.3+0.1); cpgsci(0); cpgtext(0.1,0.30,"New par"); cpgsci(1);
   cpgrect(0.2,0.28,0.28,0.3+0.1); cpgsci(0); cpgtext(0.2,0.30,"New tim"); cpgsci(1);
   cpgrect(0.3,0.38,0.28,0.3+0.1); cpgsci(0); cpgtext(0.3,0.30,"Restart"); cpgsci(1);
   /*  cpgrect(0.4,0.48,0.28,0.3+0.1); cpgsci(0); cpgtext(0.4,0.50,"Undo"); cpgsci(1);*/

}

void drawMenu3(pulsar *psr, float plotx1,float plotx2,float ploty1,float ploty2,int menu,int xplot,int yplot)
{
   float x1,y1,x2,y2,xscale,yscale,x3,x4,y3,y4;
   float x0,y0,xscale2,yscale2;
   float mx,my;
   char key;

   cpgsci(1);
   cpgtext(0.5,1.03,"x");
   cpgtext(0.6,1.03,"y");
   drawAxisSel(0,1.0,"pre-fit",xplot==1,yplot==1);
   drawAxisSel(0,0.94,"post-fit",xplot==2,yplot==2);
   drawAxisSel(0,0.88,"date",xplot==3,yplot==3);
   drawAxisSel(0,0.82,"orbital phase",xplot==4,yplot==4);
   drawAxisSel(0,0.76,"serial",xplot==5,yplot==5);
   drawAxisSel(0,0.70,"day of year",xplot==6,yplot==6);
   drawAxisSel(0,0.64,"frequency",xplot==7,yplot==7);
   drawAxisSel(0,0.58,"TOA error",xplot==8,yplot==8);
   drawAxisSel(0,0.52,"year",xplot==10,yplot==10);
   drawAxisSel(0,0.46,"elevation",xplot==11,yplot==11);
   drawAxisSel(0,0.40,"rounded MJD",xplot==12,yplot==12);
   drawAxisSel(0,0.34,"sidereal time",xplot==13,yplot==13);
   drawAxisSel(0,0.28,"hour angle",xplot==14,yplot==14);
   drawAxisSel(0,0.22,"para. angle",xplot==15,yplot==15);
   if(psr[0].TNRedAmp != 0 && psr[0].TNRedGam != 0){
	drawAxisSel(0,0.16,"Red Noise",xplot==16,yplot==16);
   }
   if(psr[0].TNDMAmp != 0 && psr[0].TNDMGam != 0){
        drawAxisSel(0,0.10,"DM Var",xplot==17,yplot==17);
   }

}

void drawMenu3_2(pulsar *psr, float plotx1,float plotx2,float ploty1,float ploty2,int menu,int xplot,int yplot,int jumpOffset, int iFlagColour, int nFlags)
{
   float x1,y1,x2,y2,xscale,yscale,x3,x4,y3,y4;
   float x0,y0,xscale2,yscale2;
   float mx,my;
   char key;
   char jumps[100],dummy[100];
   float xpos,ypos;
   int i,j;

   cpgrect(0,0.08,1.0,0.90-0.03);   cpgsci(0); cpgtext(0,0.9,"Quit"); cpgsci(1);
   cpgrect(0.1,0.18,1.0,0.90-0.03); cpgsci(0); cpgtext(0.1,0.90,"Clear"); cpgsci(1);
   cpgrect(0.2,0.28,1.0,0.90-0.03); cpgsci(0); cpgtext(0.2,0.90,"Measure"); cpgsci(1);
   cpgrect(0.3,0.38,1.0,0.90-0.03); cpgsci(0); cpgtext(0.3,0.90,"Help"); cpgsci(1);

   cpgrect(0,0.08,0.8,0.7-0.03); cpgsci(0); cpgtext(0,0.7,"List"); cpgsci(1);
   cpgrect(0.1,0.18,0.80,0.7-0.03); cpgsci(0); cpgtext(0.1,0.7,"Print"); cpgsci(1);
   cpgrect(0.2,0.28,0.80,0.7-0.03); cpgsci(0); cpgtext(0.2,0.7,"Highlight"); cpgsci(1);

   if (psr[0].nJumps > 0)
   {
	  xpos=0;
	  ypos=0.45;
	  for (i=jumpOffset+1;i<=psr[0].nJumps;i++)
	  {
		 sscanf(psr[0].jumpStr[i],"%s %s",dummy,jumps);
		 if (iFlagColour == 1)
		 {
			for (j=0;j<nFlags; j++)
			   if (strcmp(jumps,flagStore[j])==0)
			   {
				  cpgsci(j+1);
				  if (j>=14)
					 cpgsci(j-12);
				  break;
			   }
			cpgmove(xpos,ypos-0.03);
			cpgdraw(xpos+0.18,ypos-0.03);
			cpgdraw(xpos+0.18,ypos+0.15);
			cpgdraw(xpos,ypos+0.15);
			cpgdraw(xpos,ypos-0.03);
		 }
		 if (psr[0].fitJump[i]==1) cpgsci(3);
		 else cpgsci(2);

		 cpgtext(xpos,ypos,jumps);
		 xpos+=0.2;
		 if (xpos==1)
		 {
			xpos=0.0;
			ypos-=0.2;
		 }
	  }
	  cpgsci(1);
	  //if ctrl-i was pressed, add the zeroth backend (the one we're using as reference for jumps)
	  //WARNING: if the psr[0].nJumps happens to be multiple of 15, this won't be displayed (FIX IT!)
	  /*      if (iFlagColour == 1 && psr[0].nJumps - jumpOffset < 15) {
			  cpgtext(xpos,ypos,flagStore[0]);
			  cpgmove(xpos,ypos-0.03);
			  cpgdraw(xpos+0.18,ypos-0.03);
			  cpgdraw(xpos+0.18,ypos+0.15);
			  cpgdraw(xpos,ypos+0.15);
			  cpgdraw(xpos,ypos-0.03);
			  } */
   }
   if (psr[0].nJumps > 15)
   {
	  if (jumpOffset >= 5) {
		 xpos = 0.985; ypos = 0.45;
		 cpgpt(1,&xpos,&ypos,30);
	  }
	  if (psr[0].nJumps - jumpOffset > 15) {
		 xpos = 0.985; ypos = 0.15;
		 cpgpt(1,&xpos,&ypos,31);
	  }
   }
}

void drawAxisSel(float x,float y,char *str,int sel1,int sel2)
{
   cpgsci(1);
   cpgtext(x,y-0.03,str);
   if (sel1==0) cpgsci(1);
   else cpgsci(2);
   cpgrect(x+0.5,x+0.55,y,y-0.04);
   if (sel2==0) cpgsci(1);
   else cpgsci(2);
   cpgrect(x+0.6,x+0.65,y,y-0.04); 
   cpgsci(1);
}
void drawOption(float x,float y,char *str,int fit)
{
   if (fit==0) cpgsci(1);
   else cpgsci(2);
   cpgrect(x,x+0.015,y,y+0.08);
   cpgsci(1);
   if (strlen(str) > 7)
	  cpgsch(0.8);
   else
	  cpgsch(1);

   cpgtext(x+0.017,y,str);
   cpgsch(1);
}

void newTim(pulsar *psr)
{
   int i,j;
   char newtim[1000],format[1000];
   FILE *fout;
   printf("Enter format (tempo2/parkes) ");
   scanf("%s",format);
   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
   if(strcmp(format,"tempo2")*strcmp(format,"parkes")!=0){
	  printf("Enter format (tempo2/parkes) ");
	  scanf("%s",format);
	  getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
   }

   printf("New .tim file name ");
   scanf("%s",newtim);
   getchar(); // Apparently gcc doesn't flush stdin with fflush(stdin)
   writeTim(newtim,psr,format);

   /*  if (!(fout = fopen(newtim,"w")))
	   printf("Unable to open file %s\n",newtim);
	   else
	   {
	   if (psr[0].fitMode==1) fprintf(fout,"MODE 1\n");
	   for (i=0;i<psr[0].nobs;i++)
	   {
	   if (psr[0].obsn[i].deleted==1)
	   fprintf(fout,"C");
	   fprintf(fout," %s %8.3f  %.13Lf    0.00 %7.2f        %s",psr[0].obsn[i].fname,(double)psr[0].obsn[i].freq,
	   psr[0].obsn[i].sat,(double)psr[0].obsn[i].toaErr,psr[0].obsn[i].telID);
	   fprintf(fout,"\n");
	   for (j=0;j<psr[0].nPhaseJump;j++)
	   {
	   if (i<psr[0].nobs && psr[0].obsn[i].bat < psr[0].phaseJump[j] && psr[0].obsn[i+1].bat > psr[0].phaseJump[j])
	   fprintf(fout,"PHASE %d\n",psr[0].phaseJumpDir[j]);
	   }
	   }	
	   fclose(fout);
	   }*/
}

/* Routine to overplot the timing model over the pre-fit values */
void viewModels(pulsar *psr,float x1,float x2,longdouble centreEpoch,int removeMean,double mean,int count,
	  int *id,int fitFlag,float *x,float *y)
{
   longdouble diff;
   //  float px[MAX_OBSN],py[MAX_OBSN],pyt[MAX_OBSN],pxt[MAX_OBSN];
   float *px,*py,*pyt,*pxt;
   double xval,yval;
   double afunc[MAX_OBSN];
   int i,j,k,n=1,npts,arr,arr2;
   printf("Starting view models\n");
   if ((px = (float *)malloc(sizeof(float)*MAX_OBSN))==NULL){
	  printf("Unable to allocate enough memory to view models\n");
	  exit(1);
   }
   if ((py = (float *)malloc(sizeof(float)*MAX_OBSN))==NULL){
	  printf("Unable to allocate enough memory to view models\n");
	  exit(1);
   }
   if ((pyt = (float *)malloc(sizeof(float)*MAX_OBSN))==NULL){
	  printf("Unable to allocate enough memory to view models\n");
	  exit(1);
   }
   if ((pxt = (float *)malloc(sizeof(float)*MAX_OBSN))==NULL){
	  printf("Unable to allocate enough memory to view models\n");
	  exit(1);
   }
   printf("Allocated memory\n");
   for (i=0;i<MAX_OBSN;i++)
	  pyt[i] = 0.0;

   if (MAX_PSR==1)
   {
	  printf("Sorry, to view models you require -npsr 2 (or larger) on the command line\n");
	  return;
   }

   for (j=0;j<MAX_PARAMS;j++)
   {
	  for (arr=0;arr<psr[0].param[j].aSize;arr++)
	  {
		 copyPSR(psr,0,1); /* Create a new pulsar */
		 for (i=0;i<MAX_PARAMS;i++)
			copyParam(psr[0].param[i],&(psr[1].param[i]));
		 strcpy(psr[1].JPL_EPHEMERIS,psr[0].JPL_EPHEMERIS);
		 strcpy(psr[1].clock,psr[0].clock);
		 strcpy(psr[1].clockFromOverride,psr[0].clockFromOverride);
		 for (i=0;i<MAX_PARAMS;i++)
		 {
			for (arr2=0;arr2<psr[0].param[i].aSize;arr2++)
			   psr[1].param[i].fitFlag[arr2] = 0;
		 }
		 if (psr[0].param[j].fitFlag[arr]==1)
		 {
			n++;
			psr[1].param[j].fitFlag[arr]=1;
			diff = (psr[0].param[j].val[arr] - psr[0].param[j].prefit[arr]);
			//	      printf("diff = %g\n",(double)diff);
			if (j==param_f && arr==1) diff*=-1.0/psr[0].param[param_f].prefit[0]*SECDAY*SECDAY;
			else if (j==param_f && arr==0) diff*=-1.0;
			else if (j==param_f && arr==2) diff/=-((psr[0].param[param_f].prefit[0]/pow(SECDAY,3)));
			else if (j==param_om) diff/=(180.0/M_PI); // degrees -> radians
			else if (j==param_omdot) diff/=((SECDAY*365.25)*180.0/M_PI); // deg/sec -> rad/yr
			else if (j==param_pmra) diff/=(180.0/M_PI*60.0*60.0*
				  1000.0*365.25*cos(psr[0].param[param_decj].val[0]));
			else if (j==param_pmdec) diff/=(180.0/M_PI*60.0*60.0*1000.0*365.25);
			//	      else if (j==param_px) diff/=(180.0/M_PI*60.0*60.0*1000.0);

			if (j!=param_wave_om)
			{
			   for (i=0;i<count;i++)
			   {
				  /* xval is the point in time we are looking at.
					 i is the number of the observation we are observing.
					 id[i] is the observation number (What is the difference between i and id[i]?)
					 afunc will hold the derivative of the parameter versus time.
					 afunc[1] = 1: arbitrary offset in the parameter
					 afunc[2] = first-order derivative, calculated in FITfuncs (doFit.C)
					 and, deeper, in getParamDeriv (doFit.C)
					 py is the y-value correction by the fit. It is calculated as follows:
					 py = afunc[2]*diff+offset.
					 I suspect offset is generally 0
					 diff is the change in the parameter prefit-postfit. 
					 psr[1] is used as a dummy to be 'corrupted' by FITfuncs
					 */
				  xval = (double)(psr[0].obsn[id[i]].bat - psr[0].param[param_pepoch].val[0]);
				  if (fitFlag==1){ // X-axis = MJD
					 px[i] = xval-(double)(centreEpoch - psr[0].param[param_pepoch].val[0]);
				  }
				  else if (fitFlag==3) // X-axis = Binary Phase
				  {
					 double pbdot=0.0;
					 double tpb;
					 if( psr[0].param[param_t0].paramSet[0] ){
						tpb = ( psr[0].obsn[id[i]].bat - psr[0].param[param_t0].val[0])
						   / ( psr[0].param[param_pb].val[0] );
					 }else if( psr[0].param[param_tasc].paramSet[0] ){
						tpb = ( psr[0].obsn[id[i]].bat - psr[0].param[param_tasc].val[0] )
						   / ( psr[0].param[param_pb].val[0] );
					 }else{
						printf( "ERROR: Neither Tasc not T0 set...\n");
						tpb = ( psr[0].obsn[id[i]].bat - psr[0].param[param_t0].val[0] )
						   / ( psr[0].param[param_pb].val[0] );
					 }
					 double phase;

					 if (psr[0].param[param_pb].paramSet[0]==0)
					 {
						printf("This is not a binary pulsar\n");
						px[i] = 0.0;
					 }
					 else
					 {
						if (psr[0].param[param_pbdot].paramSet[0] == 1)
						   pbdot = psr[0].param[param_pbdot].val[0];

						/* Add 1000000 to make sure that the number is positive?) */
						phase = fortranMod(tpb-0.5*pbdot*tpb*tpb+1000000.0,1.0);
						if (phase < 0.0) phase+=1.0; 
						px[i] = (float)phase;
					 }
				  }
				  else if(fitFlag==7){ // X-axis = Day of Year
					 px[i] = (float)fmod((float)(double)psr[0].obsn[id[i]].bat,365.25);
				  }
				  pxt[i] = px[i];
				  FITfuncs(xval,afunc,2,&psr[1],id[i],0); 
				  py[i] = (double)(afunc[2]*diff+(double)psr[0].offset);
				  pyt[i] += (double)(afunc[2]*diff);		  
				  if (removeMean==1) py[i]-=mean;
			   }
			   npts = count;
			}
			else  /* plot fitwaves curve */
			{
			   double om;
			   for (i=0;i<256;i++){
				  /*  xval = (double)(psr[0].obsn[id[i]].bat - psr[0].param[param_pepoch].val);*/
				  /*  if (fitFlag==1) px[i] = xval-(double)(centreEpoch - psr[0].param[param_pepoch].val); */

				  xval = x[0]+i*(x[count-1]-x[0])/256.0 + (centreEpoch - psr[0].param[param_pepoch].val[0]); /* -centreEpoch+psr[0].param[param_pepoch].val; */
				  yval = 0.0;
				  om = psr[0].param[param_wave_om].val[0];
				  for (k=0;k<psr[0].nWhite;k++)
					 yval+=-psr[0].wave_sine[k]*sin(om*(k+1)*xval)-psr[0].wave_cos[k]*cos(om*(k+1)*xval);
				  px[i] = xval-(double)(centreEpoch - (double)psr[0].param[param_pepoch].val[0]);
				  py[i] = yval;
				  if (removeMean==1) py[i]-=mean; 
				  py[i]+=psr[0].offset;
				  printf("%f %f wave\n",px[i],py[i]);
			   }
			   for (i=0;i<count;i++) /* Now obtain the curve at each observation */
			   {
				  xval = (double)(psr[0].obsn[id[i]].bat - psr[0].param[param_pepoch].val[0]);
				  pxt[i] = xval-(double)(centreEpoch - psr[0].param[param_pepoch].val[0]);
				  om = psr[0].param[param_wave_om].val[0];
				  yval=0.0;
				  for (k=0;k<psr[0].nWhite;k++)
					 yval+=-psr[0].wave_sine[k]*sin(om*(k+1)*xval)-psr[0].wave_cos[k]*cos(om*(k+1)*xval);

				  pyt[i] += yval;
			   }
			   npts = 256;	      
			}
			// Somehow, certain parameters are clearly DC offset. I therefore
			// recalculate (and subtract) the mean.
			// Joris V., 20 Dec 2006.
			int tempParam;
			float tPmean=0.0;
			for(tempParam =0;tempParam<npts;tempParam++){
			   tPmean+=py[tempParam]/npts;
			}
			for(tempParam=0;tempParam<npts;tempParam++)
			   py[tempParam]-=tPmean;
			// Here the old code starts again. If you don't trust me, comment out the prevous 5 lines.
			cpgsci(n); cpgline(npts,px,py); cpgsci(1);
			cpgtext(px[count-1],py[count-1],psr[0].param[j].shortlabel[arr]);
		 }
	  }
   }


   for (i=0;i<count;i++)
   {
	  if (removeMean==1)
		 pyt[i]-=mean;
	  pyt[i]+=psr[0].offset;
   }
   // Same as before
   int tempParam;
   float tPmean=0.0;
   for(tempParam=0;tempParam<count;tempParam++){
	  tPmean+=pyt[tempParam]/count;
   }
   for(tempParam=0;tempParam<count;tempParam++)
	  pyt[tempParam]-=tPmean;
   // Joris.
   cpgsci(1); cpgline(count,pxt,pyt);
}

void reFit(int fitFlag,int setZoomX1,int setZoomX2,float zoomX1,float zoomX2,longdouble origStart,
	  longdouble origFinish,longdouble centreEpoch,pulsar *psr,int npsr,int plotX,char *dcmFile,char *covarFuncFile,int zoom)
{
   int i,k;
   printf("Redoing the fit %d (%g %d %d)\n",zoom,(double)origStart,setZoomX1,fitFlag);
   if (fitFlag==1 || fitFlag==2)
   {
	  if (setZoomX1 == 0) {
		 if (origStart==-1)
		 {
			psr[0].param[param_start].fitFlag[0]=0;
		 }
		 else
		 {
			psr[0].param[param_start].paramSet[0]=1;
			psr[0].param[param_start].val[0] = origStart;
			psr[0].param[param_start].fitFlag[0]=1;
		 }
	  }

	  if (setZoomX2 == 0) {
		 if (origFinish==-1)
		 {
			psr[0].param[param_finish].fitFlag[0]=0;
		 }
		 else
		 {
			psr[0].param[param_finish].paramSet[0]=1;
			psr[0].param[param_finish].val[0] = origFinish;
			psr[0].param[param_finish].fitFlag[0]=1;
		 }
	  }
	  if (setZoomX1 == 1) {
		 psr[0].param[param_start].paramSet[0]=1;
		 psr[0].param[param_start].fitFlag[0]=1;
		 if (plotX==12)
			psr[0].param[param_start].val[0] = zoomX1;
		 else
			psr[0].param[param_start].val[0] = zoomX1+centreEpoch;
		 psr[0].param[param_start].prefit[0] = psr[0].param[param_start].val[0];
	  }
	  if (setZoomX2 == 1) {
		 psr[0].param[param_finish].paramSet[0]=1; 
		 psr[0].param[param_finish].fitFlag[0]=1;
		 if (plotX==12)
			psr[0].param[param_finish].val[0] = zoomX2;
		 else
			psr[0].param[param_finish].val[0] = zoomX2+centreEpoch;
		 psr[0].param[param_finish].prefit[0] = psr[0].param[param_finish].val[0];
	  }
	  if (zoom==0)
	  {
		 psr[0].param[param_start].fitFlag[0]=0;
		 psr[0].param[param_finish].fitFlag[0]=0;
	  }

   }
   /* Convert all prefit values to postfit values */
   for (i=0;i<MAX_PARAMS;i++)
   {
	  for (k=0;k<psr->param[i].aSize;k++)
	  {
		 if (psr->param[i].paramSet[k]==1)
		 {
			psr->param[i].prefit[k] = psr->param[i].val[k];
			psr->param[i].prefitErr[k] = psr->param[i].prefitErr[k];
		 }
	  }
   }
   doFitAll(psr,npsr,covarFuncFile);
   formBatsAll(psr,npsr);	  
   /* Form residuals */
   formResiduals(psr,npsr,1); // iteration);
   textOutput(psr,npsr,0,0,0,0,"");

   printf("Leaving with %g %d\n",(double)origFinish,psr[0].param[param_finish].fitFlag[0]);
}


void displayStatistics(float *x,float *y,int count,float plotx1,float plotx2,
	  float ploty1,float ploty2)
{
   int i;
   float px[MAX_OBSN],py[MAX_OBSN];
   double sum=0.0,sumsq=0.0,mean,sdev;
   int n=0;

   for (i=0;i<count;i++)
   {
	  if (x[i]>plotx1 && x[i]<plotx2 && y[i]>ploty1 && y[i]<ploty2)
	  {
		 px[n]  = x[i];
		 py[n]  = y[i];
		 sum   += y[i];
		 sumsq += y[i]*y[i];
		 n++;
	  }
   }
   mean = sum/(double)n;
   sdev = sqrt((sumsq/n-sum/n*sum/n));

   printf("Number of points   = %d\n",n);
   printf("Minimum x value    = %g\n",px[0]);  
   printf("Maximum x value    = %g\n",px[n-1]);
   printf("x range            = %g\n",px[n-1]-px[0]);
   printf("Minimum y value    = %g\n",py[0]);  
   printf("Maximum y value    = %g\n",py[n-1]);
   printf("y range            = %g\n",py[n-1]-py[0]);
   printf("Sum y              = %g\n",sum);
   printf("Sum y^2            = %g\n",sumsq);
   printf("Mean y             = %g\n",mean);
   printf("Standard deviation = %g\n",sdev);
   printf("\n\n");
}

int setPlot(float *x,int count,pulsar *psr,int iobs,double unitFlag,int plotPhase,int plot,int *userValChange,char *userCMD,char *userValStr,float *userX,longdouble centreEpoch,int log,char *flagStr)
{


   if (plot==1)
   {
	  if (plotPhase==-1)
		 x[count] = (float)(double)psr[0].obsn[iobs].prefitResidual/unitFlag;
	  else
		 x[count] = (float)(double)psr[0].obsn[iobs].prefitResidual/unitFlag*psr[0].param[param_f].val[0];
   }
   else if (plot==2)
   {
	  if (plotPhase==-1)
		 x[count] = (float)(double)psr[0].obsn[iobs].residual/unitFlag;
	  else
		 x[count] = (float)(double)psr[0].obsn[iobs].residual/unitFlag*psr[0].param[param_f].val[0];
   }
   else if (plot==3)       /* Get barycentric arrival time */
	  x[count] = (float)(double)(psr[0].obsn[iobs].bat-centreEpoch);
   else if (plot==4)       /* Orbital phase */
   {
	  double pbdot=0.0;
	  double tpb;
	  if( psr[0].param[param_t0].paramSet[0] ){
		 tpb = ( psr[0].obsn[iobs].bat - psr[0].param[param_t0].val[0] )
			/ ( psr[0].param[param_pb].val[0] );
	  }else if( psr[0].param[param_tasc].paramSet[0] ){
		 tpb = ( psr[0].obsn[iobs].bat - psr[0].param[param_tasc].val[0] )
			/ ( psr[0].param[param_pb].val[0] );
	  }else{
		 printf( "ERROR: Neither T0 nor tasc set...\n" );
		 tpb = ( psr[0].obsn[iobs].bat - psr[0].param[param_t0].val[0] )
			/ ( psr[0].param[param_pb].val[0] );
	  }
	  double phase;
	  if (psr[0].param[param_pb].paramSet[0]==0)
	  {
		 printf("WARNING: This is not a binary pulsar\n");
		 x[count]=0.0;
	  }
	  else
	  {
		 if (psr[0].param[param_pbdot].paramSet[0] == 1)
			pbdot = psr[0].param[param_pbdot].val[0];

		 /*		phase = 2.0*M_PI*fortranMod(tpb-0.5*pbdot*tpb*tpb,1.0); */

		 /* Add 1000000 to make sure that the number is positive?) */
		 phase = fortranMod(tpb+1000000.0,1.0);
		 if (phase < 0.0) phase+=1.0; 
		 x[count] = (float)phase;
	  }
   }
   else if (plot==5)  /* TOA number */
	  x[count] = (float)iobs;
   else if (plot==7)  /* Observing frequency */
	  x[count] = (float)psr[0].obsn[iobs].freq;
   else if (plot==8)  /* TOA error */
	  x[count] = (float)psr[0].obsn[iobs].toaErr;
   else if (plot==9)  /* User value */
   {
	  char str[1000];
	  FILE *pin;		
	  if (*userValChange!=0)
	  {
		 if (psr[0].obsn[iobs].deleted == 0)
		 {
			if (strcasecmp(userCMD,"vip")==0)
			   sprintf(str,"%s -c %s %s | awk '{print $3}' |tail -1",userCMD,userValStr,psr[0].obsn[iobs].fname);
			else
			{
			   if (strcasecmp(userValStr,"snr")==0)
			   {
				  //sprintf(str,"vap -r %s | awk '{print $2}' |tail -1",psr[0].obsn[iobs].fname); OUTDATED
				  sprintf(str,"psrstat -Q -c snr -j FTp  %s | awk '{print $2}' |tail -1",psr[0].obsn[iobs].fname);
			   }
			   else
				  sprintf(str,"%s -c %s %s | awk '{print $2}' |tail -1",userCMD,userValStr,psr[0].obsn[iobs].fname);
			}
			printf("Running %s\n",str);
			pin = popen(str,"r");
			fscanf(pin,"%f",&x[count]);
			userX[count] = x[count];
			pclose(pin);
			*userValChange=2;
		 }
	  }
	  else
		 x[count] = userX[count];
   }
   else if (plot==6 || plot==10)  /* Day of the year     */
   {
	  double jd,fjd,day;
	  int ijd,b,c,d,e,g,month,year;
	  int retYr,retDay,stat;

	  jd = psr[0].obsn[iobs].bat + 2400000.5;
	  ijd = (int)(jd+0.5);
	  fjd = (jd+0.5)-ijd;
	  if (ijd > 2299160)
	  {
		 int a;
		 a = (int)((ijd-1867216.25)/36524.25);
		 b = ijd + 1 + a - (int)(a/4.0);
	  }
	  else
		 b = ijd;

	  c = b + 1524;
	  d = (int)((c - 122.1)/365.25);
	  e = (int)(365.25*d);
	  g = (int)((c-e)/30.6001);
	  day = c-e+fjd-(int)(30.6001*g);
	  if (g<13.5)
		 month = g-1;
	  else
		 month = g-13;
	  if (month>2.5)
		 year = d-4716;
	  else
		 year = d-4715;
	  slaCalyd(year, month, (int)day, &retYr, &retDay, &stat);

	  //      double J,Y,C,M,Day,Month,Year;
	  //      J = psr[0].obsn[iobs].bat + 2400001 + 68569;
	  //      C = 4 * J / 146097;
	  //      J = J - (146097 * C + 3) / 4;
	  //      Y = 4000 * (J + 1) / 1461001;
	  //      J = J - 1461 * Y / 4 + 31;
	  //      M = 80 * J / 2447;
	  //      Day = J - 2447 * M / 80;
	  //      J = M / 11;
	  //      Month = M + 2 - (12 * J);
	  //      Year = 100 * (C - 49) + Y + J;
	  if (plot==10)
		 x[count] = retYr+(retDay+(day-(int)day))/365.25;
	  else if (plot==6)
		 x[count] = retDay+(day-(int)day);
   }
   else if (plot==11) /* Source elevation */
   {
	  observatory *obs;
	  double source_elevation;

	  obs = getObservatory(psr[0].obsn[iobs].telID);
	  // get source elevation neglecting proper motion
	  source_elevation = asin(dotproduct(psr[0].obsn[iobs].zenith,
			   psr[0].posPulsar)
			/ obs->height_grs80);
	  x[count] = source_elevation*180.0/M_PI;
   }
   else if (plot==12) /* Nice MJD units */
   {
	  x[count] = (float)psr[0].obsn[iobs].bat;
   }
   else if (plot==13 || plot==14 || plot==15) 
	  /* 13 = Sidereal time, 14 = hour angle, 15 = parallactic angle */
   {
	  double tsid,sdd,erad,hlt,alng,hrd;
	  double toblq,oblq,pc,ph;
	  double siteCoord[3],ha;
	  observatory *obs;
	  //      printf("In here\n");
	  obs = getObservatory(psr[0].obsn[iobs].telID);
	  erad = sqrt(obs->x*obs->x+obs->y*obs->y+obs->z*obs->z);//height(m)
	  hlt  = asin(obs->z/erad); // latitude
	  alng = atan2(-obs->y,obs->x); // longitude
	  hrd  = erad/(2.99792458e8*499.004786); // height (AU)
	  siteCoord[0] = hrd * cos(hlt) * 499.004786; // dist from axis (lt-sec)
	  siteCoord[1] = siteCoord[0]*tan(hlt); // z (lt-sec)
	  siteCoord[2] = alng; // longitude

	  toblq = (psr[0].obsn[iobs].sat+2400000.5-2451545.0)/36525.0;
	  oblq = (((1.813e-3*toblq-5.9e-4)*toblq-4.6815e1)*toblq +84381.448)/3600.0;

	  pc = cos(oblq*M_PI/180.0+psr[0].obsn[iobs].nutations[1])*psr[0].obsn[iobs].nutations[0];

	  lmst2(psr[0].obsn[iobs].sat+psr[0].obsn[iobs].correctionUT1/SECDAY,0.0,&tsid,&sdd);
	  tsid*=2.0*M_PI;
	  /* Compute the local, true sidereal time */
	  ph = tsid+pc-siteCoord[2];  
	  ha = (fmod(ph,2*M_PI)-psr[0].param[param_raj].val[0])/M_PI*12;
	  if (plot==13)
		 x[count] = (float)fmod(ph/M_PI*12,24.0);
	  else if (plot==14)
		 x[count] = (float)(fmod(ha,12));
	  else if (plot==15)
	  {
		 double cp,sqsz,cqsz,pa;
		 double phi,dec;
		 phi =  hlt;
		 dec =  psr[0].param[param_decj].val[0];
		 cp =   cos(phi);
		 sqsz = cp*sin(ha*M_PI/12.0);
		 cqsz = sin(phi)*cos(dec)-cp*sin(dec)*cos(ha*M_PI/12.0);
		 if (sqsz==0 && cqsz==0) cqsz=1.0;
		 pa=atan2(sqsz,cqsz);
		 x[count] = (float)(pa*180.0/M_PI);
	  }
	  //      printf("Local sidereal time = %s %g %g %g %g %g\n",psr[0].obsn[iobs].fname,ph,tsid,pc,siteCoord[2],x[count]);
   }
   else if(plot==16){
		x[count]=(float)psr[0].obsn[iobs].TNRedSignal;
	}
     else if(plot==17){
		double DMKappa=2.410*pow(10.0,-16);
		double freq=(double)psr[0].obsn[iobs].freqSSB;
                long double yrs = (psr[0].obsn[iobs].sat - psr[0].param[param_dmepoch].val[0])/365.25;
                long double arg = 1.0;
                double dmDot=0;
                double dmDotErr=0;
                for (int d=1;d<9;d++){
                        arg *= yrs;
                        if (psr[0].param[param_dm].paramSet[d]==1){
                                dmDot+=(double)(psr[0].param[param_dm].val[d]*arg);
                        }
                }
	
                x[count]=(float)(psr[0].obsn[iobs].TNDMSignal*(DMKappa*pow(freq,2)) + dmDot);
        }
     else if (plot==18) // Plot on flag value
       {
	 int i;
	 int found=0;
	 float val;
	 for (i=0;i<psr[0].obsn[iobs].nFlags;i++)
	   {
	     if (strcmp(psr[0].obsn[iobs].flagID[i],flagStr)==0)
	       {
		 found=1;
		 sscanf(psr[0].obsn[iobs].flagVal[i],"%f",&val);
	       }
	   }
	 if (found==1)
	   x[count] = val;
	 else
	   {
	     printf("ERROR: Cannot get value for observation %s -- setting to -1\n",psr[0].obsn[iobs].fname);
	     x[count] = -1;
	   }

       }
   if (log==1 && x[count]>0)
	  x[count] = log10(x[count]);
   else if (log==1 && x[count]<0)
	  return 1;

   return 0;
}

void setLabel(char *str,int plot,int plotPhase,double unitFlag,longdouble centreEpoch,char *userValStr,char *flagStr)
{
   if (plot==1)
   {
	  if (unitFlag==1)           sprintf(str,"Prefit Residual (sec)");
	  else if (unitFlag==1.0e-3) sprintf(str,"Prefit Residual (ms)");
	  else if (unitFlag==1.0e-6) sprintf(str,"Prefit Residual (\\gms)");
	  else if (unitFlag==1.0e-9) sprintf(str,"Prefit Residual (ns)");
	  else if (unitFlag==60) sprintf(str,"Prefit Residual (minutes)");
	  else sprintf(str,"Prefit Residual (rotational period)");
	  if (plotPhase==1)	sprintf(str,"Prefit Residual (phase)");
   }
   else if (plot==2)
   {
	  if (unitFlag==1)           sprintf(str,"Postfit Residual (sec)");
	  else if (unitFlag==1.0e-3) sprintf(str,"Postfit Residual (ms)");
	  else if (unitFlag==1.0e-6) sprintf(str,"Postfit Residual (\\gms)");
	  else if (unitFlag==1.0e-9) sprintf(str,"Postfit Residual (ns)");
	  else if (unitFlag==60) sprintf(str,"Postfit Residual (minutes)");
	  else sprintf(str,"Postfit Residual (rotational period)");
	  if (plotPhase==1)	sprintf(str,"Postfit Residual (phase)");
   }
   else if (plot==3)   sprintf(str,"MJD-%.1lf",(double)centreEpoch); 
   else if (plot==4)   sprintf(str,"Orbital phase");
   else if (plot==5)   sprintf(str,"TOA number");
   else if (plot==6)   sprintf(str,"Day of year");
   else if (plot==7)   sprintf(str,"Frequency (MHz)");
   else if (plot==8)   sprintf(str,"TOA error (\\gms)");
   else if (plot==9)   sprintf(str,userValStr);
   else if (plot==10)  sprintf(str,"Year");
   else if (plot==11)  sprintf(str,"Elevation (deg)");
   else if (plot==12)  sprintf(str,"MJD");
   else if (plot==13)  sprintf(str,"Local sidereal time (hour)");
   else if (plot==14)  sprintf(str,"Hour angle (hour)");
   else if (plot==15)  sprintf(str,"Parallactic angle (deg)");
   else if (plot==16) sprintf(str,"Red Noise (sec)");
   else if (plot==17) sprintf(str,"DM Variations (cm^-3 pc)");
   else if (plot==18) sprintf(str,flagStr);
}

void averagePts(float *x,float *y,int n,int width,float *meanX,float *meanY,int *nMean)
{
   float max,min;
   int i;
   printf("nMean = %d, n = %d\n",*nMean,n);

   meanX[*nMean] = 0.0;
   meanY[*nMean] = 0.0;

   max = x[0],min=x[0];
   for (i=0;i<n;i++)
   {
	  if (max < x[i]) max=x[i];
	  if (min > x[i]) min=x[i];
	  meanY[*nMean]+=y[i];
	  meanX[*nMean]+=x[i];
   }
   meanY[*nMean]/=(double)n;
   meanX[*nMean]/=(double)n;
   if (max-min < (float)width)
   {
	  printf("RETURNING %f %g %g\n",max-min,meanX[*nMean],meanY[*nMean]);
	  (*nMean)++;
	  return;
   }  
   else
   {
	  // Find largest gap in the data-points
	  //      printf("Finding largest gap\n");
	  float gap=-1;
	  int   igap=-1;

	  for (i=1;i<n;i++)
	  {
		 if (gap < fabs(x[i]-x[i-1]))
		 {
			gap = fabs(x[i]-x[i-1]);
			igap=i;
		 }
	  }
	  //      printf("LEFT SIDE %d %d %g [%g %g]\n",n,igap,gap,x[igap],x[igap-1]);
	  // Consider left hand side
	  averagePts(x,y,igap,width,meanX,meanY,nMean);
	  // Consider right hand side
	  //      printf("RIGHT SIDE %d %d %g [%g %g]\n",n,igap,gap,x[igap],x[igap-1]);
	  averagePts(x+igap,y+igap,(n-igap),width,meanX,meanY,nMean);            
   }
}



void slaCalyd ( int iy, int im, int id, int *ny, int *nd, int *j )
   /*
	**  - - - - - - - - -
	**   s l a C a l y d
	**  - - - - - - - - -
	**
	**  Gregorian calendar date to year and day in year (in a Julian
	**  calendar aligned to the 20th/21st century Gregorian calendar).
	**
	**  (Includes century default feature:  use slaClyd for years
	**   before 100AD.)
	**
	**  Given:
	**     iy,im,id   int    year, month, day in Gregorian calendar
	**                       (year may optionally omit the century)
	**  Returned:
	**     *ny        int    year (re-aligned Julian calendar)
	**     *nd        int    day in year (1 = January 1st)
	**     *j         int    status:
	**                         0 = OK
	**                         1 = bad year (before -4711)
	**                         2 = bad month
	**                         3 = bad day (but conversion performed)
	**
	**  Notes:
	**
	**  1  This routine exists to support the low-precision routines
	**     slaEarth, slaMoon and slaEcor.
	**
	**  2  Between 1900 March 1 and 2100 February 28 it returns answers
	**     which are consistent with the ordinary Gregorian calendar.
	**     Outside this range there will be a discrepancy which increases
	**     by one day for every non-leap century year.
	**
	**  3  Years in the range 50-99 are interpreted as 1950-1999, and
	**     years in the range 00-49 are interpreted as 2000-2049.
	**
	**  Called:  slaClyd
	**
	**  Last revision:   22 September 1995
	**
	**  Copyright P.T.Wallace.  All rights reserved.
	*/
{
   int i;

   /* Default century if appropriate */
   if ( ( iy >= 0 ) && ( iy <= 49 ) )
	  i = iy + 2000;
   else if ( ( iy >= 50 ) && ( iy <= 99 ) )
	  i = iy + 1900;
   else
	  i = iy;

   /* Perform the conversion */
   slaClyd ( i, im, id, ny, nd, j );
}

void slaClyd ( int iy, int im, int id, int *ny, int *nd, int *jstat )
   /*
	**
	**  Returned:
	**     ny          int    year (re-aligned Julian calendar)
	**     nd          int    day in year (1 = January 1st)
	**     jstat       int    status:
	**                          0 = OK
	**                          1 = bad year (before -4711)
	**                          2 = bad month
	**                          3 = bad day (but conversion performed)
	**
	**  Notes:
	**
	**  1  This routine exists to support the low-precision routines
	**     slaEarth, slaMoon and slaEcor.
	**
	**  2  Between 1900 March 1 and 2100 February 28 it returns answers
	**     which are consistent with the ordinary Gregorian calendar.
	**     Outside this range there will be a discrepancy which increases
	**     by one day for every non-leap century year.
	**
	**  3  The essence of the algorithm is first to express the Gregorian
	**     date as a Julian Day Number and then to convert this back to
	**     a Julian calendar date, with day-in-year instead of month and
	**     day.  See 12.92-1 and 12.95-1 in the reference.
	**
	**  Reference:  Explanatory Supplement to the Astronomical Almanac,
	**              ed P.K.Seidelmann, University Science Books (1992),
	**              p604-606.
	**
	**  Last revision:   26 November 1994
	**
	**  Copyright P.T.Wallace.  All rights reserved.
	*/
{
   long i, j, k, l, n, iyL, imL;

   /* Month lengths in days */
   static int mtab[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };



   /* Validate year */
   if ( iy < -4711 ) { *jstat = 1; return; }

   /* Validate month */
   if ( ( im < 1 ) || ( im > 12 ) ) { *jstat = 2; return; }

   /* Allow for (Gregorian) leap year */
   mtab[1] = ( ( ( iy % 4 ) == 0 ) &&
		 ( ( ( iy % 100 ) != 0 ) || ( ( iy % 400 ) == 0 ) ) ) ?
	  29 : 28;

   /* Validate day */
   *jstat = ( id < 1 || id > mtab[im-1] ) ? 3 : 0;

   /* Perform the conversion */
   iyL = (long) iy;
   imL = (long) im;
   i = ( 14 - imL ) /12L;
   k = iyL - i;
   j = ( 1461L * ( k + 4800L ) ) / 4L
	  + ( 367L * ( imL - 2L + 12L * i ) ) / 12L
	  - ( 3L * ( ( k + 4900L ) / 100L ) ) / 4L + (long) id - 30660L;
   k = ( j - 1L ) / 1461L;
   l = j - 1461L * k;
   n = ( l - 1L ) / 365L - l / 1461L;
   j = ( ( 80L * ( l - 365L * n + 30L ) ) / 2447L ) / 11L;
   i = n + j;
   *nd = 59 + (int) ( l -365L * i + ( ( 4L - n ) / 4L ) * ( 1L - j ) );
   *ny = (int) ( 4L * k + i ) - 4716;
}

// Get sidereal time
double lmst2(double mjd,double olong,double *tsid,double *tsid_der)
{
   double xlst,sdd;
   double gmst0;
   double a = 24110.54841;
   double b = 8640184.812866;
   double c = 0.093104;
   double d = -6.2e-6;
   double bprime,cprime,dprime;
   double tu0,fmjdu1,dtu,tu,seconds_per_jc,gst;
   int nmjdu1;

   nmjdu1 = (int)mjd;
   fmjdu1 = mjd - nmjdu1;

   tu0 = ((double)(nmjdu1-51545)+0.5)/3.6525e4;
   dtu  =fmjdu1/3.6525e4;
   tu = tu0+dtu;
   gmst0 = (a + tu0*(b+tu0*(c+tu0*d)))/86400.0;
   seconds_per_jc = 86400.0*36525.0;

   bprime = 1.0 + b/seconds_per_jc;
   cprime = 2.0 * c/seconds_per_jc;
   dprime = 3.0 * d/seconds_per_jc;

   sdd = bprime+tu*(cprime+tu*dprime);

   gst = gmst0 + dtu*(seconds_per_jc + b + c*(tu+tu0) + d*(tu*tu+tu*tu0+tu0*tu0))/86400;
   xlst = gst - olong/360.0;
   xlst = fortran_mod(xlst,1.0);

   if (xlst<0.0)xlst=xlst+1.0;

   *tsid = xlst;
   *tsid_der = sdd;
   return 0.0;
}


char * plugVersionCheck = TEMPO2_h_VER;
