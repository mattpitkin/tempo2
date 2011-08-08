//  Copyright (C) 2004,2006,2007,2008,2009, George Hobbs, Russel Edwards

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

/* On linux: g++-3.2.3 -fPIC -c delays_plug.C -I/usr/local/pgplot -I/usr/local/include ; g77 -shared -o $TEMPO2PLUG/delays_linux_plug.so delays_plug.o -L/packages/linux/pgplot/ -lcpgplot -lpgplot -ltk8.0 -ltcl8.0 -L/usr/X11R6/lib -lX11 -lpng -lz -lg2c

/* On sungcc: g++-3.1.1 -fPIC -c delays_plug.C ; g++-3.1.1 -shared -o ../delays_sungcc_plug.so delays_plug.o -lcpgplot -lpgplot -lM77 -lF77 -lsunmath */
/* CC -c -I. -I../../include -I/psr/packages/include -I/psr/packages/sun4sol/include -I/usr/local/include -I/psr/include -I/psr/packages/sun4sol/include -I/usr/local/lib -o delays_plug.o delays_plug.C ; CC -G -o ../delays_plug.so delays_plug.o -L/usr/openwin/lib -lcpgplot -lpgplot -lM77 -lF77 -lsunmath -lXext -lX11 -lresolv -lsocket -lnsl */

/* This plug-in is used to plot different delays and clock errors versus TOA number etc. */
/* Author:  George Hobbs                                                                 */
/* Date:    13 Jan 2004                                                                  */
/* Version: 1.0                                                                          */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"
#include <cpgplot.h>

#define MAX_HIGHLIGHT 100 /* Maximum number of points that can be highlighted */

void doPlot(pulsar *psr,int npsr);
float findMin(float *x,pulsar *psr,int i1,int i2);
float findMax(float *x,pulsar *psr,int i1,int i2);
float findMean(float *x,pulsar *psr,int i1,int i2);
void callFit(pulsar *psr,int npsr);
float deletePoint(pulsar *psr,int npsr,float *x,float *y,float mouseX,float mouseY);
float idPoint(pulsar *psr,int npsr,float *x,float *y,float mouseX,float mouseY,int *ihighlight,int *nhighlight);
double fortranMod(double a,double p);
void createNewArrivalTimes(pulsar *psr,int npsr);


void help()
{
  printf("\n\n");
  printf("Delays plugin for tempo2\n");
  printf("\n\n");
  printf("q                          quit\n");
  printf("y-1 ('y' followed by '1')  plot first clock correction\n");
  printf("y-2 ('y' followed by '2')  plot second clock correction\n");
  printf("y-3                        plot third clock correction\n");
  printf("y-4                        plot fourth clock correction\n");
  printf("y-5                        plot fifth clock correction\n");
  printf("y-6                        plot UT1\n");
  printf("y-7                        plot Shapiro delay due to Sun\n");
  printf("y-8                        plot dispersion delay in solar system\n");
  printf("y-9                        plot dispersion delay in ISM\n");
  printf("y-0                        plot Roemer delay\n");
  printf("y-a                        plot pre-fit residuals\n");
  printf("y-b                        plot post-fit residuals\n");
  printf("y-c                        plot Shapiro delay due to Jupiter\n");
  printf("y-d                        plot Shapiro delay due to Saturn\n");
  printf("y-e                        plot Shapiro delay due to Uranus\n");
  printf("y-f                        plot Shapiro delay due to Neptune\n");
  printf("y-g                        plot total planetary Shapiro delay\n");
  printf("y-h                        plot tropospheric propagation delay\n");
  printf("1       plot TOA number on x-axis\n");
  printf("2       plot MJD day on x-axis\n");
  printf("3       plot observing frequency on x-axis\n");
  printf("4       plot day-of-year on x-axis\n");
  printf("5       plot binary phase on x-axis\n");
  printf("6       plot year on x-axis\n");
  printf("x       redo fit\n");
  printf("h       display this help\n");
  printf("f       change font size\n");
  printf("g       change the graphics terminal\n");
  printf("right mouse button    delete TOA\n");
  printf("left mouse button     identify TOA\n");
  printf("c                     clear highlights\n");
  printf("s                     start zoom region\n");
  printf("f                     finish zoom region\n");
  printf("t                     simulate new arrival times\n");
  printf("u                     unzoom\n");
}

/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int ii; // generic counter
  *npsr = 1;  /* This graphical interface will only show results for one pulsar */

  printf("Graphical Interface: delays \n");
  printf("Author:              George Hobbs (13 Jan 2004)\n");
  printf("Version:             1.0\n");

  /* Obtain the .par and the .tim file from the command line */
  if(argc==5){
      // that means: tempo2 -gr delays par.par tim.tim
      strcpy(parFile[0],argv[4]); 
      strcpy(timFile[0],argv[5]);
  }
  else{
      for(ii=0;ii<argc;ii++){
	  if(strcmp(argv[ii],"-f")==0){
	      strcpy(parFile[0],argv[++ii]);
	      strcpy(timFile[0],argv[++ii]);
	  }
	  else if(strcmp(argv[ii],"-h")==0){
	      printf("Use: \"tempo2 -gr delays par.par tim.tim\".\n");
	      printf("Then hit \"h\" in the pgplot window for more help.\n");
	  }
      }
  }

  initialise(psr,0);              /* Initialise the structures */
  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);
  callFit(psr,*npsr);             /* Do all the fitting routines */
  doPlot(psr,*npsr);              /* Do plot */
  return 0;
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
      /* Clock corrections */
      toa2utc(psr,npsr);
      //      utc2tai(psr,npsr);
      // tai2tt(psr,npsr);
      tai2ut1(psr,npsr);
      tt2tb(psr,npsr);
      
      /* Ephemeris routines */
      vectorPulsar(psr,npsr);
      readEphemeris(psr,npsr,0);
      get_obsCoord(psr,npsr);
      tt2tb(psr,npsr);          /* Observatory/time-dependent part of TT-TB */
      readEphemeris(psr,npsr,0);  /* Re-evaluate ephemeris with correct TB */ 

      /* extra delays */
      calculate_bclt(psr,npsr);
      
      /* Form barycentric arrival times */
      formBats(psr,npsr);
      
      /* Form residuals */
      formResiduals(psr,npsr,1);    
      
      /* Do the fitting */
      if (iteration==0) doFit(psr,npsr,0);
      else textOutput(psr,npsr,globalParameter,0,0,0,"");
    }

}


void doPlot(pulsar *psr,int npsr)
{
  int i,xfitFlag=1,yfitFlag=1,exitFlag=0,scale1=0,scale2=psr[0].nobs,count,ihighlight[MAX_HIGHLIGHT];
  float fontsize=1.0;
  int nhighlight=0,hcount=0;
  char xstr[1000],ystr[1000];
  float x[MAX_OBSN],y[MAX_OBSN],yerr1[MAX_OBSN],yerr2[MAX_OBSN],hx[MAX_OBSN],hy[MAX_OBSN];
  float minx,maxx,miny,maxy,plotx1,plotx2,ploty1,ploty2,mean;
  float mouseX,mouseY;
  char key;
  int changeplot=0;
  printf("Here in doplot\n");
  /* Obtain a graphical PGPLOT window */
  cpgbeg(0,"/xs",1,1);
  cpgask(0);

  scale2=psr[0].nobs;

  do {

    /*    if (fitFlag==1 || fitFlag==2)
	  sprintf(xstr,"MJD-%.1Lf",psr[0].param[param_pepoch].val); 
	  else if (fitFlag==3 || fitFlag==4)
	  sprintf(xstr,"Orbital phase");
	  else if (fitFlag==5 || fitFlag==6)
	  sprintf(xstr,"TOA number");
	  else if (fitFlag==7 || fitFlag==8)
	  sprintf(xstr,"Day of year");
	  else if (fitFlag==9 || fitFlag==10)
	  sprintf(xstr,"Frequency (MHz)");
	  sprintf(ystr,"Residual (sec)"); */
    if (xfitFlag==1) sprintf(xstr,"TOA number");
    else if (xfitFlag==2) sprintf(xstr,"MJD-%.1Lf",psr[0].param[param_pepoch].val[0]);
    else if (xfitFlag==3) sprintf(xstr,"Observing Frequency (MHz)");
    else if (xfitFlag==4) sprintf(xstr,"Day of year");
    else if (xfitFlag==5) sprintf(xstr,"Orbital phase");
    else if (xfitFlag==6) sprintf(xstr,"Year");

    if (yfitFlag==1) sprintf(ystr,psr[0].obsn[0].correctionsTT[0].corrects_to);
    else if (yfitFlag==2) sprintf(ystr,psr[0].obsn[0].correctionsTT[1].corrects_to);
    else if (yfitFlag==3) sprintf(ystr,psr[0].obsn[0].correctionsTT[2].corrects_to);
    else if (yfitFlag==4) sprintf(ystr,psr[0].obsn[0].correctionsTT[3].corrects_to);
    else if (yfitFlag==5) sprintf(ystr,psr[0].obsn[0].correctionsTT[4].corrects_to);
    else if (yfitFlag==6) sprintf(ystr,"    -> UT1 (s)");
    else if (yfitFlag==7) sprintf(ystr,"Shapiro delay due to Sun (s)");
    else if (yfitFlag==8) sprintf(ystr,"Dispersion delay SS (us)");
    else if (yfitFlag==9) sprintf(ystr,"Dispersion delay ISM (s)");
    else if (yfitFlag==10) sprintf(ystr,"Roemer delay (s)");
    else if (yfitFlag==11) sprintf(ystr,"Pre-fit residual (s)");
    else if (yfitFlag==12) sprintf(ystr,"Post-fit residual (s)");
    else if (yfitFlag==13) sprintf(ystr,"Shapiro delay due to Jupiter (ns)");
    else if (yfitFlag==14) sprintf(ystr,"Shapiro delay due to Saturn (ns)");
    else if (yfitFlag==15) sprintf(ystr,"Shapiro delay due to Uranus (ns)");
    else if (yfitFlag==16) sprintf(ystr,"Shapiro delay due to Neptune (ns)");
    else if (yfitFlag==17) sprintf(ystr,"Total planetary Shapiro delay (ns)");
    else if (yfitFlag==18) sprintf(ystr,"Tropospheric delay (ns)");
    count=0;
    printf("Number of observations = %d\n",psr[0].nobs);
    for (i=0;i<psr[0].nobs;i++)
      {
	if (psr[0].obsn[i].deleted == 0)
	  {
	    /* Plotting on x-axis */
	    if (xfitFlag==1) x[count] = (float)i;  /* TOA Number */
	    else if (xfitFlag==2)                  /* Get barycentric arrival time */
	      x[count] = (double)(psr[0].obsn[i].bat-psr[0].param[param_pepoch].val[0]);
	    else if (xfitFlag==3)                  /* Observing frequency */
	      x[count] = (double)psr[0].obsn[i].freq;
	    else if (xfitFlag==4)                  /* Day of the year */
	      x[count] = fmod((double)psr[0].obsn[i].bat,365.25);
	    else if (xfitFlag==5)                  /* Orbital phase */
	      {
		double pbdot=0.0;
		double tpb = (psr[0].obsn[i].bat-psr[0].param[param_t0].val[0])*86400.0/(psr[0].param[param_pb].val[0]*SECDAY);
		double phase;

		if (psr[0].param[param_pbdot].paramSet[0] == 1)
		  pbdot = psr[0].param[param_pbdot].val[0];

		/*		phase = 2.0*M_PI*fortranMod(tpb-0.5*pbdot*tpb*tpb,1.0); */

		/* Add 1000000 to make sure that the number is positive?) */
		phase = fortranMod(tpb+1000000.0,1.0);
		if (phase < 0.0) phase+=1.0; 
		x[count] = (float)phase;
	      }
	    else if (xfitFlag==6)
	      {
		double J,Y,C,M,Day,Month,Year;
		J = psr[0].obsn[i].bat + 2400001 + 68569;
		C = 4 * J / 146097;
		J = J - (146097 * C + 3) / 4;
		Y = 4000 * (J + 1) / 1461001;
		J = J - 1461 * Y / 4 + 31;
		M = 80 * J / 2447;
		Day = J - 2447 * M / 80;
		J = M / 11;
		Month = M + 2 - (12 * J);
		Year = 100 * (C - 49) + Y + J;
		x[count] = Year;
	      }
	    
	    /* Plotting on y-axis */
	    if (yfitFlag==1)                
	      y[count] = (double)psr[0].obsn[i].correctionsTT[0].correction;
	    else if (yfitFlag==2)
	      y[count] = (double)psr[0].obsn[i].correctionsTT[1].correction;
	    else if (yfitFlag==3)                
	      y[count] = (double)psr[0].obsn[i].correctionsTT[2].correction;
	    else if (yfitFlag==4)                
	      y[count] = (double)psr[0].obsn[i].correctionsTT[3].correction;
	    else if (yfitFlag==5)                
	      y[count] = (double)psr[0].obsn[i].correctionsTT[4].correction;
	    else if (yfitFlag==6)                
	      y[count] = (double)psr[0].obsn[i].correctionUT1;
	    else if (yfitFlag==7)
	      y[count] = (double)psr[0].obsn[i].shapiroDelaySun;
	    else if (yfitFlag==8)
	      y[count] = (double)psr[0].obsn[i].tdis2*1e6;
	    else if (yfitFlag==9)
	      {
		y[count] = (double)psr[0].obsn[i].tdis1;
		printf("Have %d %Lf %g %g %g\n",i,psr[0].obsn[i].sat,x[count],y[count],psr[0].obsn[i].tdis1);
	      }
	    else if (yfitFlag==10)
	      y[count] = (double)psr[0].obsn[i].roemer;
	    else if (yfitFlag==11)         /* Get pre-fit residual */
	      y[count] = (double)psr[0].obsn[i].prefitResidual;
	    else if (yfitFlag==12)           /* Post-fit residual */
	      y[count] = (double)psr[0].obsn[i].residual;
	    else if (yfitFlag==13)
	      y[count] = (double)psr[0].obsn[i].shapiroDelayJupiter*1e9;
	    else if (yfitFlag==14)
	      y[count] = (double)psr[0].obsn[i].shapiroDelaySaturn*1e9;
	    else if (yfitFlag==15)
	      y[count] = (double)psr[0].obsn[i].shapiroDelayUranus*1e9;
	    else if (yfitFlag==16)
	      y[count] = (double)psr[0].obsn[i].shapiroDelayNeptune*1e9;
	    else if (yfitFlag==17)
	      y[count] = (double)(psr[0].obsn[i].shapiroDelayNeptune+psr[0].obsn[i].shapiroDelayUranus
				 +psr[0].obsn[i].shapiroDelaySaturn+psr[0].obsn[i].shapiroDelayJupiter)*1e9;
	    else if (yfitFlag==18)
	      y[count] = (double)psr[0].obsn[i].troposphericDelay*1e9;

	    count++;
	  }
      }

    /*    if (xfitFlag==2)
      {
	for (i=0;i<count;i++)
	  printf("Have %g %g\n",x[i],y[i]);
	printf("Count = %d\n",count);
	} */

    /* Remove mean from the residuals and calculate error bars */
    if (yfitFlag==11 || yfitFlag==12)
      {
	mean = findMean(y,psr,scale1,scale2);
	count=0;
	for (i=0;i<psr[0].nobs;i++)
	  {
	    if (psr[0].obsn[i].deleted==0)
	      {
		y[count]-=mean;
		yerr1[count] = y[count]-(float)psr[0].obsn[i].toaErr*1e-6;
		yerr2[count] = y[count]+(float)psr[0].obsn[i].toaErr*1e-6;
		count++;
	      }
	  }
      }
    hcount=0;
    for (i=0;i<nhighlight;i++)
      {
	hx[hcount] = x[ihighlight[i]];
	hy[hcount] = y[ihighlight[i]];
	hcount++;
      }
    /* Get scaling for graph */
    minx = findMin(x,psr,0,count);
    maxx = findMax(x,psr,0,count);
    printf("minx = %f, maxx = %f (%d %d)\n",minx,maxx,scale1,scale2);
    plotx1 = minx-(maxx-minx)*0.1; 
    plotx2 = maxx+(maxx-minx)*0.1; 
    
    miny = findMin(y,psr,0,count);
    maxy = findMax(y,psr,0,count);
    ploty1 = miny-(maxy-miny)*0.1;
    ploty2 = maxy+(maxy-miny)*0.1;
    printf("Minimum = %f %f\n",ploty1,ploty2);
    /* Plot the residuals */
    cpgsch(fontsize);
    cpgenv(plotx1,plotx2,ploty1,ploty2,0,0);
    cpglab(xstr,ystr,psr[0].name);
    cpgpt(count,x,y,16);
    cpgsci(2); cpgsch(3); cpgpt(hcount,hx,hy,-4); cpgsci(1); cpgsch(fontsize);
    if (yfitFlag==10 || yfitFlag==11) cpgerry(count,x,yerr1,yerr2,1);

    if (changeplot==1)
      {
	cpgbeg(0,"/xs",1,1);
	cpgask(0);
	changeplot=0;
      }
    else
      {
	cpgcurs(&mouseX,&mouseY,&key);
	
	/* Check key press */
	if (key=='q') exitFlag=1;
	else if (key=='y')
	  {
	    cpgcurs(&mouseX,&mouseY,&key);
	    if (key=='1') yfitFlag=1;
	    else if (key=='2') yfitFlag=2;
	    else if (key=='3') yfitFlag=3;
	    else if (key=='4') yfitFlag=4;
	    else if (key=='5') yfitFlag=5;
	    else if (key=='6') yfitFlag=6;
	    else if (key=='7') yfitFlag=7;
	    else if (key=='8') yfitFlag=8;
	    else if (key=='9') yfitFlag=9;
	    else if (key=='0') yfitFlag=10;
	    else if (key=='a') yfitFlag=11;
	    else if (key=='b') yfitFlag=12;
	    else if (key=='c') yfitFlag=13;
	    else if (key=='d') yfitFlag=14;
	    else if (key=='e') yfitFlag=15;
	    else if (key=='f') yfitFlag=16;
	    else if (key=='g') yfitFlag=17;
	    else if (key=='h') yfitFlag=18;
	  }
	else if (key=='1') xfitFlag=1;
	else if (key=='2') xfitFlag=2;
	else if (key=='3') xfitFlag=3;
	else if (key=='4') xfitFlag=4;
	else if (key=='5' && psr[0].param[param_t0].paramSet[0]==1) xfitFlag=5;
	else if (key=='6') xfitFlag=6;
	else if (key=='x') callFit(psr,npsr);
	else if (key=='h') help();
	else if (key=='g')
	  {
	    cpgbeg(0,"?",1,1);
	    cpgask(0);
	    changeplot=1;
	  }
	else if (key=='d' || key=='X') deletePoint(psr,npsr,x,y,mouseX,mouseY); /* Delete closest point */
	else if (key=='i' || key=='A') idPoint(psr,npsr,x,y,mouseX,mouseY,ihighlight,&nhighlight); /* Identify closest point */
	else if (key=='c') nhighlight = 0; /* Clear highlights */
	else if (key=='f')
	  {
	    printf("Enter fontsize ");
	    scanf("%f",&fontsize);
	  }
	else if (key=='s')  /* Start zoom */
	  {
	    int icount=0;
	    for (i=0;i<psr[0].nobs;i++)
	      {
		if (psr[0].obsn[i].deleted==0)
		  {
		    if (x[icount]>mouseX) {scale1=i; break;}
		    icount++;
		  }
	      }
	  }
	else if (key=='t')
	  {
	    createNewArrivalTimes(psr,npsr);
	    scale2=psr[0].nobs;
	  }
	else if (key=='f')  /* Finish zoom */
	  {
	    int icount=0;
	    for (i=0;i<psr[0].nobs;i++)
	      {
		if (psr[0].obsn[i].deleted==0)
		  {
		    if (x[icount]>mouseX) {scale2=i-1; break;}
		    icount++;
		  }
	      }
	  }
	else if (key=='u') /* Unzoom */
	  {
	    scale1=0; scale2=count;
	  }
	else printf("Unknown key press %c\n",key);
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
  return 0.0;
}

float idPoint(pulsar *psr,int npsr,float *x,float *y,float mouseX,float mouseY,int *ihighlight,int *nhighlight)
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
	  if (iclosest==-1)
	    {
	      iclosest=count;
	      closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
	    }
	  else if (pow(xpos-mouseX,2)+pow(ypos-mouseY,2)<closest)
	    {
	      iclosest=count;
	      closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
	    }
	  count++;

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
  ihighlight[*nhighlight] = iclosest;
  (*nhighlight)++;
  return 0.0;
}


float findMax(float *x,pulsar *psr,int i1,int i2)
{
  int i,count;
  float max;

  max = x[i1];
  for (i=i1;i<i2;i++)
    {
      if (x[i]>max) max=x[i];
    }
  return max;
}

float findMin(float *x,pulsar *psr,int i1,int i2)
{
  int i,count;
  float min;

  min = x[i1];
  for (i=i1;i<i2;i++)
    {
      if (x[i]<min) min=x[i];
    }
  return min;
}

float findMean(float *x,pulsar *psr,int i1,int i2)
{
  int i,count;
  float mean;

  count=0;
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

void createNewArrivalTimes(pulsar *psr,int npsr)
{
  int i;
  double start,step;
  int nobs;

  printf("Enter start (MJD) step (d) number of observations ");
  scanf("%lf %lf %d",&start,&step,&nobs);
  psr[0].nobs = nobs;
  for (i=0;i<nobs;i++)
    {
      psr[0].obsn[i].sat = start+i*step;
      psr[0].obsn[i].clockCorr = 1;
      psr[0].obsn[i].delayCorr = 1;
      psr[0].obsn[i].deleted = 0;
      psr[0].obsn[i].freq = 1400.0;
      psr[0].obsn[i].toaErr = 1e-6;
      sprintf(psr[0].obsn[i].fname,"SIMULATED_%d\n",i);
      strcpy(psr[0].obsn[i].telID,"7");
    }
  callFit(psr,npsr);             /* Do all the fitting routines */
}
char * plugVersionCheck = TEMPO2_h_VER;
