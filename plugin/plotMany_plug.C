#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"
#include <cpgplot.h>

/* using namespace std; */  /* Is this required for a plugin ? */

void doPlot(pulsar *psr,int npsr);
float findMin(float *x,pulsar *psr,int p,int i1,int i2);
float findMax(float *x,pulsar *psr,int p,int i1,int i2);
float findMean(float *x,pulsar *psr,int p,int i1,int i2);
void callFit(pulsar *psr,int npsr);
float findMinVal(float *a, int n);
float findMaxVal(float *a, int n);
double fortranMod(double a,double p);

/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i;
  FILE *pin;
  char str[1000];

  *npsr = 0;  /* This graphical interface will only show results for one pulsar */

  printf("Graphical Interface: plk emulator\n");
  printf("Author:              George Hobbs (21 Dec 2003)\n");
  printf("Version:             1.0\n");

  /* Obtain the .par and the .tim file from the command line */
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{	  
	  strcpy(parFile[*npsr],argv[i+1]); 
	  strcpy(timFile[*npsr],argv[i+2]);	  
	  (*npsr)++;
	  if (*npsr > MAX_PSR)
	    {
	      printf("Error, npsr > MAX_PSR.  Must increase MAX_PSR\n");
	      exit(1);
	    } 
	}
    }

if (*npsr==0) /* Select all files */
    {
      printf("Using all available .par and .tim files\n");
      sprintf(str,"ls `ls *.par | sed s/par/tim/` | sed s/.tim/\"\"/");
      pin = popen(str,"r");
      while (!feof(pin))
	{
	  if (fscanf(pin,"%s",str)==1)
	    {
	      sprintf(parFile[*npsr],"%s.par",str);
	      sprintf(timFile[*npsr],"%s.tim",str);
	      (*npsr)++;
	    }
	}
      pclose(pin);
      printf("Obtained files for %d pulsars\n",*npsr);
    }

  initialise(psr,0);              /* Initialise the structures */
  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);
  callFit(psr,*npsr);             /* Do all the fitting routines */
  doPlot(psr,*npsr);              /* Do plot */
}

/* This function calls all of the fitting routines.             */
/* The function is looped twice, the first time for the pre-fit */
/* residuals and the second time for the post-fit residuals     */

void callFit(pulsar *psr,int npsr)
{
  int iteration;
  double globalParameter = 0.0;

  for (iteration=0;iteration<2;iteration++)
    {
      formBatsAll(psr,npsr);
      formResiduals(psr,npsr,1);
      /* Do the fitting */
      if (iteration==0) doFit(psr,npsr,0);
      else textOutput(psr,npsr,globalParameter,0,0,0,"");
    }

}


void doPlot(pulsar *psr,int npsr)
{
  int i,j,fitFlag=2,exitFlag=0,scale1=0,scale2,count[MAX_PSR],p,xautoscale=0,k,graphics=1;
  int yautoscale=0,plotpre=1;
  int time=0;
  char xstr[1000],ystr[1000];
  float px[2],py[2],pye1[2],pye2[2];
  float x[MAX_PSR][MAX_OBSN],y[MAX_PSR][MAX_OBSN],yerr1[MAX_PSR][MAX_OBSN],yerr2[MAX_PSR][MAX_OBSN],tmax,tmin,tmaxy1,tminy1,tmaxy2,tminy2;
  float minx[MAX_PSR],maxx[MAX_PSR],miny[MAX_PSR],maxy[MAX_PSR],plotx1,plotx2,ploty1,ploty2,mean;
  float mouseX,mouseY;
  float fontSize=1.4;
  char key;
  float widthPap=0.0,aspectPap=0.618;

  /* Obtain a graphical PGPLOT window */
  cpgbeg(0,"/ps",1,1);
  cpgpap(widthPap,aspectPap);
  cpgsch(fontSize);
  cpgask(0);

  for (p=0;p<npsr;p++)
    {
      scale2 = psr[p].nobs;
      
      /*      sprintf(xstr,"MJD-%.1Lf",psr[0].param[param_pepoch].val[0]); */
      sprintf(xstr,"MJD-%.1f",51000.0); 
      sprintf(ystr,"Residual (\\gmsec)");
      
      count[p]=0;
      printf("points = %d\n",psr[p].nobs);
      for (i=0;i<psr[p].nobs;i++)
	{	  
	  if (psr[p].obsn[i].deleted == 0)
	    {
	      /* x[p][count[p]] = (double)(psr[p].obsn[i].bat-psr[0].param[param_pepoch].val[0]);	     	       */
	      x[p][count[p]] = (double)(psr[p].obsn[i].bat-51000);	     	      
	      y[p][count[p]] = (double)psr[p].obsn[i].residual*1.0e6;
	      count[p]++;
	    }
	}
      /* Remove mean from the residuals and calculate error bars */
      mean = findMean(y[p],psr,p,scale1,count[p]);
      count[p]=0;
      for (i=0;i<psr[p].nobs;i++)
	{
	  if (psr[p].obsn[i].deleted==0   &&
	      (psr[p].param[param_start].paramSet[0]!=1 || psr[p].param[param_start].fitFlag[0]!=1 ||
	       psr[p].param[param_start].val[0] < psr[p].obsn[i].bat) &&
	      (psr[p].param[param_finish].paramSet[0]!=1 || psr[p].param[param_finish].fitFlag[0]!=1 ||
	       psr[p].param[param_finish].val[0] > psr[p].obsn[i].bat))
	    {
	      psr[p].obsn[i].residual-=mean/1.0e6;
	      y[p][count[p]]-=mean;
	      yerr1[p][count[p]] = y[p][count[p]]-(float)psr[p].obsn[i].toaErr;
	      yerr2[p][count[p]] = y[p][count[p]]+(float)psr[p].obsn[i].toaErr;
	      count[p]++;
	    }
	}
    	  
      /* Get scaling for graph */
      minx[p] = findMin(x[p],psr,p,scale1,count[p]);
      maxx[p] = findMax(x[p],psr,p,scale1,count[p]);

      miny[p] = findMin(y[p],psr,p,scale1,count[p]);
      maxy[p] = findMax(y[p],psr,p,scale1,count[p]);
    }

  tmin = findMinVal(minx,npsr);
  tmax = findMaxVal(maxx,npsr);

  tminy2 = findMinVal(miny,npsr);
  tmaxy2 = findMaxVal(maxy,npsr);

  plotx1 = tmin-(tmax-tmin)*0.1;
  plotx2 = tmax+(tmax-tmin)*0.3;
  
  ploty1 = tminy2-(tmaxy2-tminy2)*0.1;
  ploty2 = tmaxy2+(tmaxy2-tminy2)*0.1;
	
  for (p=0;p<npsr;p++)
    {
      for (i=0;i<count[p];i++)
	{
	  y[p][i]=(p+1)+y[p][i]/(ploty2-ploty1);
	  yerr1[p][i]=(p+1)+yerr1[p][i]/(ploty2-ploty1);
	  yerr2[p][i]=(p+1)+yerr2[p][i]/(ploty2-ploty1);
	}
    } 
  
  printf("ytick = %g\n",ploty2-ploty1);
      /*  cpgenv(plotx1,plotx2,ploty1,ploty2+(ploty2-ploty1)*(npsr-1),0,0); */
  cpgenv(plotx1,plotx2,0,npsr+1,0,-1);
  cpgbox("ATNSBC",0.0,0,"",0.0,0);
  cpglab(xstr,"","");	    

  for (p=0;p<npsr;p++)
    {
      float xx[MAX_OBSN],yy[MAX_OBSN],yyerr1[MAX_OBSN],yyerr2[MAX_OBSN];
      int num=0,colour;
      cpgsch(1);
      cpgtext(tmax+(tmax-tmin)*0.05,p+1,psr[p].name);
      cpgsch(1);
      
      cpgsls(4);
      px[0] = plotx1;
      px[1] = tmax+(tmax-tmin)*0.03;
      py[0] = p+1;
      py[1] = p+1;
      cpgline(2,px,py);
      px[0] = plotx1+0.01*(plotx2-plotx1);
      py[0] = p+1;
      pye1[0] = p+1 + 5/(ploty2-ploty1);
      pye2[0] = p+1 - 5/(ploty2-ploty1);
      cpgsls(1);
      cpgsch(3);
      cpgerry(1,px,pye1,pye2,1); 
      cpgsch(1);

      for (colour=0;colour<5;colour++)
	{
	  num=0;
	  for (i=0;i<count[p];i++)
	    {
	      if ((colour==0 && psr[p].obsn[i].freq<=500) ||
		  (colour==1 && psr[p].obsn[i].freq>500 && psr[p].obsn[i].freq<=1000) ||
		  (colour==2 && psr[p].obsn[i].freq>1000 && psr[p].obsn[i].freq<=1500) ||
		  (colour==3 && psr[p].obsn[i].freq>1500 && psr[p].obsn[i].freq<=3300) ||
		  (colour==4 && psr[p].obsn[i].freq>3300))
		{
		  xx[num]=x[p][i];
		  yy[num]=y[p][i];
		  yyerr1[num]=yerr1[p][i];
		  yyerr2[num]=yerr2[p][i];
		  num++;
		}
	    }
	  cpgsci(colour+1);
	  cpgpt(num,xx,yy,16);
	  cpgerry(num,xx,yyerr1,yyerr2,1);
	}
      cpgsci(1);
    }

  
  cpgend();
}


float findMax(float *x,pulsar *psr,int p,int i1,int i2)
{
  int i,count;
  float max;

  max = x[i1];
  count=0;
  for (i=i1;i<i2;i++)
    {
      if (x[count]>max) max=x[count];
      count++;
    }
  return max;
}

float findMin(float *x,pulsar *psr,int p,int i1,int i2)
{
  int i,count;
  float min;

  min = x[i1];

  count=0;
  for (i=i1;i<i2;i++)
    {
      if (x[count]<min) min=x[count];
      count++;
    }
  return min;
}

float findMean(float *x,pulsar *psr,int p,int i1,int i2)
{
  int i,count;
  float mean;

  count=0;
  mean=0.0;
  for (i=i1;i<i2;i++)
    {
      mean+=x[count];
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

float findMaxVal(float *a, int n)
{
  int i;
  float max=a[0];

  for (i=0;i<n;i++)
    {
      if (max < a[i]) max = a[i];
    }
  return max;
}

float findMinVal(float *a, int n)
{
  int i;
  float min=a[0];

  for (i=0;i<n;i++)
    {
      if (min > a[i]) min = a[i];
    }
  return min;
}
