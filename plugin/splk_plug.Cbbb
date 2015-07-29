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

/* splk: SuperPLK    */
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
#include <cpgplot.h>

/* using namespace std; */  /* Is this required for a plugin ? */

void doPlot(pulsar *psr,int npsr,int overlay);
float findMin(float *x,pulsar *psr,int p,int i1,int i2);
float findMax(float *x,pulsar *psr,int p,int i1,int i2);
float findMean(float *x,pulsar *psr,int p,int i1,int i2);
void callFit(pulsar *psr,int npsr);
float deletePoint(pulsar *psr,int npsr,float *x,float *y,float mouseX,float mouseY);
float idPoint(pulsar *psr,int npsr,float *x,float *y,float mouseX,float mouseY);
double fortranMod(double a,double p);

/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */

char dcmFile[MAX_FILELEN];
char covarFuncFile[MAX_FILELEN];


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i;
  int overlay=0;
  FILE *pin;
  char str[1000];
  *npsr = 0;  /* This graphical interface will only show results for one pulsar */


  strcpy(dcmFile,"NULL");
  strcpy(covarFuncFile,"NULL");

  printf("Graphical Interface: plk emulator\n");
  printf("Author:              George Hobbs (21 Dec 2003)\n");
  printf("CVS Version:         $Revision $\n");

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
      else if (strcmp(argv[i],"-dcm")==0)
	strcpy(dcmFile,argv[++i]);
      else if (strcmp(argv[i],"-dcf")==0)
	strcpy(covarFuncFile,argv[++i]);
      else if (strcmp(argv[i],"-overlay")==0)
	overlay=1;
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

//  initialise(psr,0);              /* Initialise the structures */
  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);
  callFit(psr,*npsr);             /* Do all the fitting routines */
  doPlot(psr,*npsr,overlay);              /* Do plot */
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
      if (iteration==0) doFitAll(psr,npsr,covarFuncFile);
      else textOutput(psr,npsr,globalParameter,0,0,0,"");
    }

}


void doPlot(pulsar *psr,int npsr,int overlay)
{
  int i,j,fitFlag=1,exitFlag=0,scale1=0,scale2,count,p,xautoscale=1,k,graphics=1;
  int yautoscale=1,plotpre=1;
  int time=0;
  char xstr[1000],ystr[1000];
  float x[MAX_OBSN],y[MAX_OBSN],yerr1[MAX_OBSN],yerr2[MAX_OBSN],tmax,tmin,tmaxy1,tminy1,tmaxy2,tminy2;
  float minx,maxx,miny,maxy,plotx1,plotx2,ploty1,ploty2,mean;
  float mouseX,mouseY;
  float fontSize=1.8;
  char key;
  float widthPap=0.0,aspectPap=0.618;

  /* Obtain a graphical PGPLOT window */
  if (overlay==1)
    cpgbeg(0,"?",2,1);
  else
    cpgbeg(0,"?",2,npsr);
  cpgpap(widthPap,aspectPap);
  cpgsch(fontSize);
  cpgask(0);

  do {
    for (p=0;p<npsr;p++)
      {
	scale2 = psr[p].nobs;
	for (j=0;j<2;j++)
	  {
	    if (j==0) fitFlag=1;
	    else if (j==1) fitFlag=2;

	    sprintf(xstr,"MJD-%.1Lf",psr[0].param[param_pepoch].val[0]); 
	    sprintf(ystr,"Residual (\\gmsec)");

	    count=0;
	    for (i=0;i<psr[p].nobs;i++)
	      {
		if (psr[p].obsn[i].deleted == 0  &&
	    (psr[p].param[param_start].paramSet[0]!=1 || psr[p].param[param_start].fitFlag[0]!=1 ||
	      psr[p].param[param_start].val[0] < psr[p].obsn[i].bat) &&
	    (psr[p].param[param_finish].paramSet[0]!=1 || psr[p].param[param_finish].fitFlag[0]!=1 ||
	     psr[p].param[param_finish].val[0] > psr[p].obsn[i].bat))
		  {
		    if (xautoscale==1)
		      x[count] = (double)(psr[p].obsn[i].bat-psr[p].param[param_pepoch].val[0]);
		    else
		      x[count] = (double)(psr[p].obsn[i].bat-psr[0].param[param_pepoch].val[0]);
		    
		    if (fitFlag==1)  /* Get pre-fit residual */
		      y[count] = (double)psr[p].obsn[i].prefitResidual*1.0e6;
		    else if (fitFlag==2) /* Post-fit residual */
		      y[count] = (double)psr[p].obsn[i].residual*1.0e6;
		    count++;
		  }
	      }
	    /* Remove mean from the residuals and calculate error bars */
	    mean = findMean(y,psr,p,scale1,count);
	    count=0;
	    for (i=0;i<psr[p].nobs;i++)
	      {
		if (psr[p].obsn[i].deleted==0   &&
	    (psr[p].param[param_start].paramSet[0]!=1 || psr[p].param[param_start].fitFlag[0]!=1 ||
	      psr[p].param[param_start].val[0] < psr[p].obsn[i].bat) &&
	    (psr[p].param[param_finish].paramSet[0]!=1 || psr[p].param[param_finish].fitFlag[0]!=1 ||
	     psr[p].param[param_finish].val[0] > psr[p].obsn[i].bat))
		  {
		    psr[p].obsn[i].residual-=mean/1.0e6;
		    y[count]-=mean;
		    yerr1[count] = y[count]-(float)psr[p].obsn[i].toaErr;
		    yerr2[count] = y[count]+(float)psr[p].obsn[i].toaErr;
		    count++;
		  }
	      }
	    
	    /* Get scaling for graph */
	    minx = findMin(x,psr,p,scale1,count);
	    maxx = findMax(x,psr,p,scale1,count);
	    if (xautoscale==1)
	      {
		plotx1 = minx-(maxx-minx)*0.1;
		plotx2 = maxx+(maxx-minx)*0.1;
	      }
	    else
	      {
		plotx1 = tmin-(tmax-tmin)*0.1;
		plotx2 = tmax+(tmax-tmin)*0.1;
	      }
	    miny = findMin(y,psr,p,scale1,count);
	    maxy = findMax(y,psr,p,scale1,count);

	    if (yautoscale==1)
	      {
		ploty1 = miny-(maxy-miny)*0.1;
		ploty2 = maxy+(maxy-miny)*0.1;
	      }
	    else
	      {
		if (j==0)
		  {
		    ploty1 = tminy1-(tmaxy1-tminy1)*0.1;
		    ploty2 = tmaxy1+(tmaxy1-tminy1)*0.1;
		  }
		else
		  {
		    ploty1 = tminy2-(tmaxy2-tminy2)*0.1;
		    ploty2 = tmaxy2+(tmaxy2-tminy2)*0.1;
		  }
	      }

	    /* Plot the residuals */	    
	    if (plotpre==1 || j!=0)
	      {
		float xx[MAX_OBSN],yy[MAX_OBSN],yyerr1[MAX_OBSN],yyerr2[MAX_OBSN];
		int num=0,colour;
		if (overlay==0 || (overlay==1 && p==0))
		  {
		    cpgenv(plotx1,plotx2,ploty1,ploty2,0,0);
		    cpglab(xstr,ystr,psr[p].name);	    
		  }

		for (colour=0;colour<5;colour++)
		  {
		    num=0;
		    for (i=0;i<count;i++)
		      {
			if ((colour==0 && psr[p].obsn[i].freq<=500) ||
			    (colour==1 && psr[p].obsn[i].freq>500 && psr[p].obsn[i].freq<=1000) ||
			    (colour==2 && psr[p].obsn[i].freq>1000 && psr[p].obsn[i].freq<=1500) ||
			    (colour==3 && psr[p].obsn[i].freq>1500 && psr[p].obsn[i].freq<=3300) ||
			    (colour==4 && psr[p].obsn[i].freq>3300))
			  {
			    xx[num]=x[i];
			    yy[num]=y[i];
			    yyerr1[num]=yerr1[i];
			    yyerr2[num]=yerr2[i];
			    num++;
			  }
		      }
		    cpgsci(colour+1);
		    if (overlay==1)
		      cpgsci(p+1);
		    cpgpt(num,xx,yy,16);
		    cpgerry(num,xx,yyerr1,yyerr2,1);
		  }
		cpgsci(1);
	      }
	  }
      }
    printf("------------------------------\n");
    printf("`a'     set aspect ratio\n");
    printf("`f'     set font size\n");
    printf("`g'     set graphics device\n");
    printf("`q'     quit\n");
    printf("`x'     toggle autoscale x axis\n");
    printf("`y'     toggle autoscale y axis\n");
    printf("`p'     toggle prefit plotting\n");
    printf("`r'     output residuals to file\n");

    if (graphics==1)
      {
	cpgcurs(&mouseX,&mouseY,&key);
	
	/* Check key press */
	if (key=='q') exitFlag=1;
	if (key=='p') {
	  plotpre*=-1;
	  if (plotpre==-1)
	    {
	      cpgend();
	      if (overlay==1)
		cpgbeg(0,"/xs",1,1);
	      else
		cpgbeg(0,"/xs",1,npsr);
	      cpgpap(widthPap,aspectPap);
	      cpgsch(fontSize);
	      cpgask(0);	      
	    }
	  else
	    {
	      cpgend();
	      if (overlay==1)
		cpgbeg(0,"/xs",2,1);
	      else
		cpgbeg(0,"/xs",2,npsr);
	      cpgpap(widthPap,aspectPap);
	      cpgsch(fontSize);
	      cpgask(0);	      
	    }
	}
	else if (key=='a') /* Change aspect ratio */
	  {
	    printf("Please enter a new aspect ratio ");
	    scanf("%f",&aspectPap);
	    cpgend();
	    cpgbeg(0,"/xs",2,npsr);
	    cpgpap(widthPap,aspectPap);
	    cpgsch(fontSize);
	    cpgask(0);	      
	  }
	else if (key=='f') /* Change font size */
	  {
	    printf("Please enter a new font size ");
	    scanf("%f",&fontSize);
	    cpgend();
	    cpgbeg(0,"/xs",2,npsr);
	    cpgpap(widthPap,aspectPap);
	    cpgsch(fontSize);
	    cpgask(0);	      
	  }
	else if (key=='g')
	  {
	    graphics=0;
	    cpgend();
	    if (plotpre==-1)
	      {
		cpgend();
		if (overlay==1)
		  cpgbeg(0,"?",1,1);
		else
		  cpgbeg(0,"?",1,npsr);
		cpgpap(widthPap,aspectPap);
		cpgsch(fontSize);
		cpgask(0);	      
	      }
	    else
	      {
		cpgend();
		if (overlay==1)
		  cpgbeg(0,"?",1,1);
		else
		  cpgbeg(0,"?",2,npsr);
		cpgpap(widthPap,aspectPap);
		cpgsch(fontSize);
		cpgask(0);	      
	      }
	  }
	else if (key=='r') /* Output residuals to file */
	  {
	    FILE *fout;
	    char fname[1000];
	    int ii,jj;

	    for (ii=0;ii<npsr;ii++)
	      {
		sprintf(fname,"%s.res",psr[ii].name);
		fout = fopen(fname,"w");
		/* Print header */
		fprintf(fout,"#PSR %s\n",psr[ii].name);
		fprintf(fout,"#F0  %.14Lf\n",psr[ii].param[param_f].val[0]);
		fprintf(fout,"#RAJ %s\n",psr[ii].rajStrPre);
		fprintf(fout,"#DECJ %s\n",psr[ii].decjStrPre);
		for (jj=0;jj<psr[ii].nobs;jj++)
		  fprintf(fout,"%.5lf %.5lg %.5lg\n",
			  (double)(psr[ii].obsn[jj].bat-psr[0].param[param_pepoch].val[0]),
			  (double)(psr[ii].obsn[jj].residual),(double)(psr[ii].obsn[jj].toaErr)/1.0e6);
		fclose(fout);
	      }
	  }
	else if (key=='x') 
	  {
	    xautoscale*=-1;
	    if (xautoscale==-1)
	      {
		for (k=0;k<npsr;k++)
		  {
		    count=0;
		    for (i=0;i<psr[k].nobs;i++)
		      {
			if (psr[k].obsn[i].deleted==0   &&
			    (psr[k].param[param_start].paramSet[0]!=1 || psr[k].param[param_start].fitFlag[0]!=1 ||
			     psr[k].param[param_start].val[0] < psr[k].obsn[i].bat) &&
			    (psr[k].param[param_finish].paramSet[0]!=1 || psr[k].param[param_finish].fitFlag[0]!=1 ||
			     psr[k].param[param_finish].val[0] > psr[k].obsn[i].bat))

			  {x[count] = (double)(psr[k].obsn[i].bat-psr[0].param[param_pepoch].val[0]); count++;}
		      }
		    minx = findMin(x,psr,k,scale1,count);
		    maxx = findMax(x,psr,k,scale1,count);
		    if (k==0)
		      {
			tmin = minx;
			tmax = maxx;
			printf("Have1 tmin = %f, tmax = %f\n",tmin,tmax);
		      }
		    else
		      {
			if (tmin > minx) tmin = minx;
			if (tmax < maxx) tmax = maxx;
			printf("Have2 tmin = %f, tmax = %f\n",tmin,tmax);

		      }
		  }
	      }
	  }       
	else if (key=='y') 
	  {
	    yautoscale*=-1;
	    if (yautoscale==-1)
	      {
		for (k=0;k<npsr;k++)
		  {
		    count=0;
		    for (i=0;i<psr[k].nobs;i++)
		      {
			if (psr[k].obsn[i].deleted==0   &&
			    (psr[k].param[param_start].paramSet[0]!=1 || psr[k].param[param_start].fitFlag[0]!=1 ||
			     psr[k].param[param_start].val[0] < psr[k].obsn[i].bat) &&
			    (psr[k].param[param_finish].paramSet[0]!=1 || psr[k].param[param_finish].fitFlag[0]!=1 ||
			     psr[k].param[param_finish].val[0] > psr[k].obsn[i].bat))			  
			  {y[count] = (double)psr[k].obsn[i].prefitResidual*1e6; count++;}
		      }
		    miny = findMin(y,psr,k,scale1,count);
		    maxy = findMax(y,psr,k,scale1,count);
		    if (k==0)
		      {
			tminy1 = miny;
			tmaxy1 = maxy;
		      }
		    else
		      {
			if (tminy1 > miny) tminy1 = miny;
			if (tmaxy1 < maxy) tmaxy1 = maxy;
		      }

		    count=0;
		    for (i=0;i<psr[k].nobs;i++)
		      {
			if (psr[k].obsn[i].deleted==0   &&
			    (psr[k].param[param_start].paramSet[0]!=1 || psr[k].param[param_start].fitFlag[0]!=1 ||
			     psr[k].param[param_start].val[0] < psr[k].obsn[i].bat) &&
			    (psr[k].param[param_finish].paramSet[0]!=1 || psr[k].param[param_finish].fitFlag[0]!=1 ||
			     psr[k].param[param_finish].val[0] > psr[k].obsn[i].bat))
			  {y[count] = (double)psr[k].obsn[i].residual*1e6; count++;}
		      }
		    miny = findMin(y,psr,k,scale1,count);
		    maxy = findMax(y,psr,k,scale1,count);
		    if (k==0)
		      {
			tminy2 = miny;
			tmaxy2 = maxy;
		      }
		    else
		      {
			if (tminy2 > miny) tminy2 = miny;
			if (tmaxy2 < maxy) tmaxy2 = maxy;
		      }
		  }
		printf("Have tminy2 = %g %g\n",tminy2,tmaxy2);
	      }
	  }       
	else printf("Unknown key press %c\n",key);
      }
    else
      {
	graphics=1;

	cpgend();
	if (plotpre==-1)
	  {
	    cpgend();
	    if (overlay==1)
	      cpgbeg(0,"/xs",1,1);
	    else
	      cpgbeg(0,"/xs",1,npsr);
	    cpgpap(widthPap,aspectPap);
	    cpgsch(fontSize);
	    cpgask(0);	      
	  }
	else
	  {
	    cpgend();
	    if (overlay==1)
	      cpgbeg(0,"/xs",2,1);
	    else
	      cpgbeg(0,"/xs",2,npsr);
	    cpgpap(widthPap,aspectPap);
	    cpgsch(fontSize);
	    cpgask(0);	      
	  }
      }
  } while (exitFlag==0);
  
  cpgend();
}

float deletePoint(pulsar *psr,int npsr,float *x,float *y,float mouseX,float mouseY)
{
  int i,iclosest,count;
  float closest,x1,x2,x3,x4,y1,y2,y3,y4,xscale,yscale,xpos,ypos;

  cpgqvp(3,&x1,&x2,&y1,&y2);
  cpgqwin(&x3,&x4,&y3,&y4);
  xscale = (x2-x1)/(x4-x3);
  yscale = (y2-y1)/(y4-y3);
  mouseX = (mouseX-x3)*xscale;
  mouseY = (mouseY-y3)*yscale;
  iclosest=-1;
  count=0;
  for (i=0;i<psr[0].nobs;i++)
    {
      if (psr[0].obsn[i].deleted==0)
	{
	  xpos = (x[count]-x3)*xscale;
	  ypos = (y[count]-y3)*yscale;
	  count++;
	  if (iclosest==-1)
	    {
	      iclosest=i;
	      closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
	    }
	  else if (pow(xpos-mouseX,2)+pow(ypos-mouseY,2)<closest)
	    {
	      iclosest=i;
	      closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
	    }
	}
    }
  psr[0].obsn[iclosest].deleted=1;
}

float idPoint(pulsar *psr,int npsr,float *x,float *y,float mouseX,float mouseY)
{
  int i,iclosest,count;
  float closest,x1,x2,x3,x4,y1,y2,y3,y4,xscale,yscale,xpos,ypos;

  cpgqvp(3,&x1,&x2,&y1,&y2);
  cpgqwin(&x3,&x4,&y3,&y4);
  xscale = (x2-x1)/(x4-x3);
  yscale = (y2-y1)/(y4-y3);
  mouseX = (mouseX-x3)*xscale;
  mouseY = (mouseY-y3)*yscale;
  iclosest=-1;
  count=0;
  for (i=0;i<psr[0].nobs;i++)
    {
      if (psr[0].obsn[i].deleted==0)
	{
	  xpos = (x[count]-x3)*xscale;
	  ypos = (y[count]-y3)*yscale;
	  count++;
	  if (iclosest==-1)
	    {
	      iclosest=i;
	      closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
	    }
	  else if (pow(xpos-mouseX,2)+pow(ypos-mouseY,2)<closest)
	    {
	      iclosest=i;
	      closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
	    }
	}
    }
  printf("---------------------------------------------------\n");
  printf("Closest point has TOA number %d (starting from 1)\n",iclosest+1);
  printf("SAT = %Lf\n",psr[0].obsn[iclosest].sat);
  printf("BAT = %Lf\n",psr[0].obsn[iclosest].bat);
  printf("Pre-fit residual = %Lf\n",psr[0].obsn[iclosest].prefitResidual);
  printf("Post-fit residual = %Lf\n",psr[0].obsn[iclosest].residual);
  printf("Observing frequency = %f\n",psr[0].obsn[iclosest].freq);
  printf("---------------------------------------------------\n");
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

char * plugVersionCheck = TEMPO2_h_VER;
