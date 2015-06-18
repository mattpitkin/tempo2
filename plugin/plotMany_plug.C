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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"
#include <cpgplot.h>

// using namespace std;   /* Is this required for a plugin ? */
char covarFuncFile2[MAX_FILELEN];

void doPlot(pulsar *psr,int npsr,float *scale,int nScale,char *grDev,int plotUs,float fontsize,float centreMJD,int ptStyle,float ptSize,int error,float miny,float maxy,float minx,float maxx,int nOverlay,float labelsize,float fracX);
float findMin(float *x,pulsar *psr,int p,int i1,int i2);
float findMax(float *x,pulsar *psr,int p,int i1,int i2);
float findMean(float *x,pulsar *psr,int p,int i1,int i2);
void callFit(pulsar *psr,int npsr);
float findMinVal(float *a, int n);
float findMaxVal(float *a, int n);
double fortranMod(double a,double p);
double lmst2(double mjd,double olong,double *tsid,double *tsid_der);
void slaClyd ( int iy, int im, int id, int *ny, int *nd, int *jstat );
void slaCalyd ( int iy, int im, int id, int *ny, int *nd, int *j );
float calcYr(float mjd);


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */

void help()
{
  printf("plotMany: plugin to plot the timing residuals for multiple pulsars\n\n");
  
  printf("-centremjd val Set the central mjd value to 'val'.  Use -1 to use years\n");
  printf("-fontsize val  Set the fontsize to 'val'\n");
  printf("-h             This help\n");
  printf("-g val         Set the graphics device to 'val'\n");
  printf("-plotus        Scale in microseconds\n");
  printf("-ptsize val    Set the point size to 'val'\n");
  printf("-ptstyle val   Set the point style to 'val'\n");
  printf("-reverse       Reverse the ordering of the pulsars\n");

  exit(1);
}

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char oparFile[MAX_PSR][MAX_FILELEN]; // Overlay files
  char timFile[MAX_PSR][MAX_FILELEN];
  char otimFile[MAX_PSR][MAX_FILELEN];
  char tpar[MAX_PSR][MAX_FILELEN];
  char ttim[MAX_PSR][MAX_FILELEN];
  int i;
  int plotUs=0;
  float scale[1000];
  int nScale=0;
  char grDev[100]="/vps";
  //  char grDev[100]="/xs";
  FILE *pin;
  FILE *fin;
  char str[1000];
  float fontsize = 1.2;
  float centreMJD=54500;
  int ptStyle=16;
  float ptSize = 1;
  int reverse=0;
  int error=1;
  float miny=-1;
  float maxy=-1;
  float minx=-1;
  float maxx=-1;
  int nOverlay=0;
  float labelsize=-1;
  float fracX = 0.85;
  *npsr = 0; 
  strcpy(covarFuncFile2,"NULL");

  printf("Graphical Interface: plotMany\n");
  printf("Author:              George Hobbs\n");
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
      else if (strcmp(argv[i],"-overlay")==0)
	sscanf(argv[i+1],"%d",&nOverlay);
      else if (strcmp(argv[i],"-dcf")==0){
        strcpy(covarFuncFile2,argv[++i]);
	  }
      else if (strcmp(argv[i],"-g")==0)
	strcpy(grDev,argv[++i]);
      else if (strcmp(argv[i],"-h")==0)
	help();
      else if (strcmp(argv[i],"-scale")==0)
	{
	  sscanf(argv[i+1],"%f",&scale[nScale]);
	  nScale++;
	}
      else if (strcmp(argv[i],"-fontsize")==0)
	sscanf(argv[i+1],"%f",&fontsize);
      else if (strcasecmp(argv[i],"-fracX")==0)
	sscanf(argv[i+1],"%f",&fracX);
      else if (strcasecmp(argv[i],"-labelSize")==0)
	sscanf(argv[i+1],"%f",&labelsize);
      else if (strcmp(argv[i],"-ptstyle")==0)
	  sscanf(argv[i+1],"%d",&ptStyle);
      else if (strcmp(argv[i],"-ptsize")==0)
	  sscanf(argv[i+1],"%f",&ptSize);
      else if (strcmp(argv[i],"-miny")==0)
	  sscanf(argv[i+1],"%f",&miny);
      else if (strcmp(argv[i],"-maxy")==0)
	  sscanf(argv[i+1],"%f",&maxy);
      else if (strcmp(argv[i],"-minx")==0)
	  sscanf(argv[i+1],"%f",&minx);
      else if (strcmp(argv[i],"-maxx")==0)
	  sscanf(argv[i+1],"%f",&maxx);
      else if (strcmp(argv[i],"-reverse")==0)
	reverse=1;
      else if (strcmp(argv[i],"-centremjd")==0)
	sscanf(argv[i+1],"%f",&centreMJD);
      else if (strcmp(argv[i],"-plotus")==0)
	plotUs=1;
      else if (strcmp(argv[i],"-noerror")==0)
	error=0;
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
 if (reverse==1) // Reverse order of .par and .tim files
   {
     for (i=0;i<*npsr;i++)
       {
     	 strcpy(tpar[i],parFile[i]);
     	 strcpy(ttim[i],timFile[i]);	 
       }
     for (i=0;i<*npsr;i++)
       {
	 //	 printf("Have %d %d %d >%s< >%s<\n",i,*npsr,(*npsr)-i-1,tpar[(*npsr)-i-1],ttim[(*npsr)-i-1]);
	 strcpy(parFile[i],tpar[(*npsr)-i-1]);
	 strcpy(timFile[i],ttim[(*npsr)-i-1]);
	 //	 strcpy(timFile[i],"J0900-3144.new2.tim");
       }
   }

 printf("Starting plugin1\n");
 //  initialise(psr,0);              /* Initialise the structures */
 readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
 printf("Starting plugin2\n");
 readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
 printf("Starting plugin3\n");
  preProcess(psr,*npsr,argc,argv);
 printf("Starting plugin4\n");
  callFit(psr,*npsr);             /* Do all the fitting routines */
 printf("Starting plugin5\n");
 doPlot(psr,*npsr,scale,nScale,grDev,plotUs,fontsize,centreMJD,ptStyle,ptSize,error,miny,maxy,minx,maxx,nOverlay,labelsize,fracX);              /* Do plot */
 printf("Starting plugin6\n");
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
      if (iteration==0) doFitAll(psr,npsr,covarFuncFile2);
      else textOutput(psr,npsr,globalParameter,0,0,0,"");
    }

}


void doPlot(pulsar *psr,int npsr,float *scale,int nScale,char *grDev,int plotUs,float fontSize,float centreMJD,int ptStyle,float ptSize,int error,float minyv,float maxyv,float minxv,float maxxv,int nOverlay,float labelsize,float fracX)
{
  int i,j,fitFlag=2,exitFlag=0,scale1=0,scale2,count[MAX_PSR],p,xautoscale=0,k,graphics=1;
  int yautoscale=0,plotpre=1;
  int ps,pe,pi;
  int time=0;
  char xstr[1000],ystr[1000];
  float px[2],py[2],pye1[2],pye2[2];
  float x[MAX_PSR][MAX_OBSN],y[MAX_PSR][MAX_OBSN],yerr1[MAX_PSR][MAX_OBSN],yerr2[MAX_PSR][MAX_OBSN],tmax,tmin,tmaxy1,tminy1,tmaxy2,tminy2;
  float sminy[MAX_PSR],smaxy[MAX_PSR];
  float minx[MAX_PSR],maxx[MAX_PSR],miny[MAX_PSR],maxy[MAX_PSR],plotx1,plotx2,ploty1,ploty2,mean;
  float fx[2],fy[2];
  float mouseX,mouseY;
  char key;
  //  float widthPap=0.0,aspectPap=0.618;
  float widthPap=0.0,aspectPap=1;
  float xx[MAX_OBSN],yy[MAX_OBSN],yyerr1[MAX_OBSN],yyerr2[MAX_OBSN];
  int num=0,colour;

  /* Obtain a graphical PGPLOT window */
  cpgbeg(0,grDev,1,1);
  //    cpgpap(widthPap,aspectPap);
  cpgsch(fontSize);
  cpgscf(2);
  cpgslw(2);
  cpgask(0);

  for (p=0;p<npsr;p++)
    {
      scale2 = psr[p].nobs;
      
      /*      sprintf(xstr,"MJD-%.1Lf",psr[0].param[param_pepoch].val[0]); */
      if (centreMJD == -1)
	sprintf(xstr,"Year"); 
      else
	sprintf(xstr,"MJD-%.1f",centreMJD); 

      sprintf(ystr,"Residual (\\gmsec)");
      
      count[p]=0;
      printf("points = %d\n",psr[p].nobs);
      for (i=0;i<psr[p].nobs;i++)
	{	  
	  if (psr[p].obsn[i].deleted == 0 &&
	      (psr[p].param[param_start].paramSet[0]!=1 || psr[p].param[param_start].fitFlag[0]!=1 ||
	       psr[p].param[param_start].val[0] < psr[p].obsn[i].bat) &&
	      (psr[p].param[param_finish].paramSet[0]!=1 || psr[p].param[param_finish].fitFlag[0]!=1 ||
	       psr[p].param[param_finish].val[0] > psr[p].obsn[i].bat))
	    {
	      /* x[p][count[p]] = (double)(psr[p].obsn[i].bat-psr[0].param[param_pepoch].val[0]);	     	       */
	      if (centreMJD == -1)
		x[p][count[p]] = calcYr(psr[p].obsn[i].bat);
	      else
		x[p][count[p]] = (double)(psr[p].obsn[i].bat-centreMJD); 
	      y[p][count[p]] = (double)psr[p].obsn[i].residual*1.0e6;
	      if (nScale>0)
		y[p][count[p]] *= scale[p];
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
      if (minxv == maxxv) {
	minx[p] = findMin(x[p],psr,p,scale1,count[p]);
	maxx[p] = findMax(x[p],psr,p,scale1,count[p]);
      }
      else {
	minx[p] = minxv;
	maxx[p] = maxxv;
      }
      if (minyv == maxyv){
	miny[p] = findMin(y[p],psr,p,scale1,count[p]);
	maxy[p] = findMax(y[p],psr,p,scale1,count[p]);
      }
      else {
	miny[p] = minyv;
	maxy[p] = maxyv;
      }
      sminy[p] = miny[p]/1e6;
      smaxy[p] = maxy[p]/1e6;
    }
  for (p=0;p<npsr;p++)
    {
      for (i=0;i<count[p];i++)
	{
	  y[p][i] = (y[p][i]-miny[p])/(maxy[p]-miny[p]);
	  yerr1[p][i] = (yerr1[p][i]-miny[p])/(maxy[p]-miny[p]);
	  yerr2[p][i] = (yerr2[p][i]-miny[p])/(maxy[p]-miny[p]);
	}
      //      maxy[p] = 1.0;
      //      miny[p] = 0.0;
    }
  

  tmin = findMinVal(minx,npsr);
  tmax = findMaxVal(maxx,npsr);

  tminy2 = 0.0; //findMinVal(miny,npsr);
  tmaxy2 = 1.0; //findMaxVal(maxy,npsr);

  plotx1 = tmin-(tmax-tmin)*0.1;
  plotx2 = tmax+(tmax-tmin)*0.1;
  
  //  ploty1 = tminy2-(tmaxy2-tminy2)*0.1;
  //  ploty2 = tmaxy2+(tmaxy2-tminy2)*0.1;
	
  ploty1 = 0.1;
  ploty2 = 0.9;

  for (p=0;p<npsr;p++)
    {
      for (i=0;i<count[p];i++)
	{
	  y[p][i]=(p)+ploty1+y[p][i]*(ploty2-ploty1);
	  yerr1[p][i]=(p)+ploty1+yerr1[p][i]*(ploty2-ploty1);
	  yerr2[p][i]=(p)+ploty1+yerr2[p][i]*(ploty2-ploty1);
	}
    } 
  
  printf("ytick = %g\n",ploty2-ploty1);
      /*  cpgenv(plotx1,plotx2,ploty1,ploty2+(ploty2-ploty1)*(npsr-1),0,0); */
  //  cpgenv(plotx1,plotx2,0,npsr+1,0,-1);

  if (labelsize!=-1)
    cpgsch(labelsize);
  cpgsvp(fracX,1.0,0.1,1.0);
  cpgswin(0,1,0,npsr);
  cpgbox("ABC",0.0,0,"C",0.0,0);
  cpgsch(fontSize);
  char str[1000];
  for (p=0;p<npsr;p++)
    {
      cpgsch(fontSize);
      //      cpgtext(tmax+(tmax-tmin)*0.05,p+1.5-0.5,psr[p].name);
      cpgtext(0,p+0.6,psr[p].name);
      //      cpgsch(fontSize);
      if (plotUs==0)
	{
	  sprintf(str,"%.2f",(double)((smaxy[p]-sminy[p])*psr[p].param[param_f].val[0]));
	  cpgtext(0,p+0.4,str);
	  //	  cpgtext(tmax+(tmax-tmin)*0.05,p+1.1-0.5,str);
	}
      else
	{
	  sprintf(str,"%.2f\\gms",(double)((smaxy[p]-sminy[p])/1e-6));
	  //	  cpgtext(tmax+(tmax-tmin)*0.05,p+1.1-0.5,str);
	  cpgtext(0,p+0.1,str);
	}
      cpgsch(1);
      px[0] = 0;
      //      px[1] = tmax; //+(tmax-tmin)*0.03;
	px[1] = 1;
      py[0] = p;
      py[1] = p;
      cpgline(2,px,py);
      
    }
  if (labelsize!=-1)
    cpgsch(labelsize);

  cpgsvp(0.1,fracX,0.1,1.0);
  cpgswin(plotx1,plotx2,0,npsr);
  cpgbox("ATNSBC",0.0,0,"B",0.0,0);
  cpglab(xstr,"","");	    
  cpgsch(fontSize);

  for (p=0;p<npsr;p++)
    {
      cpgsls(1);
      px[0] = plotx1;
      //      px[1] = tmax; //+(tmax-tmin)*0.03;
      px[1] = plotx2;
      py[0] = p;
      py[1] = p;
      cpgline(2,px,py);
      cpgsls(4);
      px[0] = tmin;
      px[1] = tmax+(tmax-tmin)*0.03;

      py[0]=py[1] =(p)+ploty1+(-miny[p]/(maxy[p]-miny[p]))*(ploty2-ploty1);
      //      py[0]=py[1] = (p)+ploty1;
      //      py[0] = py[1] = (0-miny[p])/(maxy[p]-miny[p])/(ploty2-ploty1)+p;
      cpgline(2,px,py);

      px[0] = plotx1+0.005*(plotx2-plotx1);
      py[0] = p;
      pye1[0] = p + 5/(ploty2-ploty1);
      pye2[0] = p - 5/(ploty2-ploty1);
      cpgsls(1);
      cpgsch(3);
      //      cpgerry(1,px,pye1,pye2,1); 
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
		  //		  printf("plotting: %g\n",yy[num]);

		  num++;
		}
	    }
	  cpgsci(colour+1);
	  cpgsch(ptSize);
	  cpgpt(num,xx,yy,ptStyle);
	  if (error==1)
	    cpgerry(num,xx,yyerr1,yyerr2,1);
	  cpgsch(fontSize);
	  // Plot arrow giving one period
	  fx[0] = fx[1] = tmin-(tmax-tmin)*0.05;
	  //	  fy[0] = (p+1)+0.5-(float)(1.0/psr[p].param[param_f].val[0])/2.0/(ploty2-ploty1);
	  //	  fy[1] = (p+1)+0.5+(float)(1.0/psr[p].param[param_f].val[0])/2.0/(ploty2-ploty1);

	  //	  fy[0] = (-(float)(1.0/psr[p].param[param_f].val[0])/2.0/1.0e6 - miny[p])/(maxy[p]-miny[p])/(ploty2-ploty1) + (p+1)+0.5;
	  //	  fy[1] = ((float)(1.0/psr[p].param[param_f].val[0])/2.0/1.0e6 - miny[p])/(maxy[p]-miny[p])/(ploty2-ploty1) + (p+1)+0.5;
	  fy[0] = (p+1)+0.5+(float)(1.0/psr[p].param[param_f].val[0])/2.0/(maxy[p]-miny[p])*1e6;
	  fy[1] = (p+1)+0.5-(float)(1.0/psr[p].param[param_f].val[0])/2.0/(maxy[p]-miny[p])*1e6;
	  if (fy[0] > (p+1)+1) fy[0] = (p+1)+1;
	  if (fy[1] < (p+1)) fy[1] = (p+1);
	  
	  //	  cpgsls(1); cpgline(2,fx,fy); cpgsls(1);
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

float calcYr(float mjd)
{
  double jd,fjd,day;
  int ijd,b,c,d,e,g,month,year;
  int retYr,retDay,stat;

  jd = mjd + 2400000.5;
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
  
  return retYr+(retDay+(day-(int)day))/365.25;

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
