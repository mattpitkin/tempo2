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
#include <cpgplot.h>
#include "T2toolkit.h"
#include "TKfit.h"

using namespace std;


void help() /* Display help */
{
  /* This function should contain usage information about the plugin which should (in general) be accessed */
  /* by the user pressing 'h'                                                                              */
}
#define MAX_TIMES 1000

void doPlot(double *epoch,double *dmVal,double *dmE,int *id,int n,pulsar *psr);
double mjd2year(double mjd);
void selectData(pulsar *psr,float *rx,float *ry,double f1,double f2,float *plotX,float *plotY,float *plotE,int *nplot);

void slaClyd ( int iy, int im, int id, int *ny, int *nd, int *jstat );
void slaCalyd ( int iy, int im, int id, int *ny, int *nd, int *j );


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  FILE *fin,*fout;
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  char timeList[100];
  int i,j,k,nread;
  long double mjd1[MAX_TIMES],mjd2[MAX_TIMES];
  long double centreMJD;
  double epoch[MAX_TIMES],dmVal[MAX_TIMES],dmE[MAX_TIMES];
  char parFileName[MAX_TIMES][MAX_STRLEN];
  char tname[100]="";
  int nStride=0;
  double globalParameter;
  int argn;
  int fitf0=0;
  int n=0;
  int id[MAX_TIMES];

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: dm\n");
  printf("Author:              George Hobbs\n");
  printf("Version:             v1.0\n");
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
      else if (strcmp(argv[i],"-t")==0)
	strcpy(timeList,argv[++i]);
      else if (strcmp(argv[i],"-fitf0")==0)
	fitf0=1;
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

  // Now calculate the DMs

  // Read the list of strides
  fin = fopen(timeList,"r");
  while (!feof(fin))
    {
      //      nread = fscanf(fin,"%Lf %Lf %s %s",&mjd1[nStride],&mjd2[nStride],parFileName[nStride],tname);
      nread = fscanf(fin,"%Lf %Lf %s %s",&mjd1[nStride],&mjd2[nStride],tname,tname);
      if (nread==3 || nread==4)
	{
	  strcpy(parFileName[nStride],parFile[0]);
	  nStride++;
	}
    }
  fclose(fin);

  // Now do the fitting
  if (!(fout = fopen("dmvals.dat","w")))
    {
      printf("ERROR: Unable to open dmvals.dat for writing\n");
      exit(1);
    }
  for (i=0;i<nStride;i++)
    {
      centreMJD = (mjd1[i]+mjd2[i])/2.0;
      printf("Analysing %g\n",(double)centreMJD);
      strcpy(parFile[0],parFileName[i]);
      psr[0].nJumps=0;
      readParfile(psr,parFile,timFile,1); /* Load the parameters       */
      readTimfile(psr,timFile,1); /* Load the arrival times    */
      printf("Number of jumps = %d %g\n",psr[0].nJumps,(double)centreMJD);
      // Update the epoch in the par file for centreMJD
      strcpy(argv[argn],"-epoch");
      sprintf(argv[argn+1],"%.5f",(double)centreMJD);
      printf("Do the preprocessing\n");
      psr[0].dmOffset=0;
      preProcess(psr,1,argc,argv);      
      printf("Done the preprocessing\n");

      // Turn off all fitting
      for (j=0;j<MAX_PARAMS;j++)
	{
	  for (k=0;k<psr[0].param[j].aSize;k++)
	    psr[0].param[j].fitFlag[k] = 0;
	}
      for (j=0;j<psr[0].nJumps;j++)
	psr[0].fitJump[j]=0;

      // Update the start and finish flags
      psr[0].param[param_start].val[0] = mjd1[i];
      psr[0].param[param_start].fitFlag[0] = 1;
      psr[0].param[param_start].paramSet[0] = 1;

      psr[0].param[param_finish].val[0] = mjd2[i];
      psr[0].param[param_finish].fitFlag[0] = 1;
      psr[0].param[param_finish].paramSet[0] = 1;

      // Turn on required fitting
      psr[0].param[param_dm].fitFlag[0] = 1;
      if (fitf0==1)
	psr[0].param[param_f].fitFlag[0] = 1;
      // Do the fitting
      printf("Do the fitting %g %g\n",(double)psr[0].param[param_start].val[0],
	     (double)psr[0].param[param_finish].val[0]);
      if (psr[0].param[param_finish].val[0] <= psr[0].param[param_start].val[0])
	{
	  printf("ERROR: In %s file, the finish position %g <= the end position %g\n",timeList,(double)psr[0].param[param_finish].val[0],(double)psr[0].param[param_start].val[0]);
	  exit(1);
	}
      for (j=0;j<2;j++)                   /* Do two iterations for pre- and post-fit residuals*/
	{
	  formBatsAll(psr,1);         /* Form the barycentric arrival times */
	  formResiduals(psr,1,1);    /* Form the residuals                 */
	  if (j==0) doFit(psr,1,0);   /* Do the fitting     */
	  else textOutput(psr,1,globalParameter,0,0,0,"");  /* Display the output */
	}
      printf("Done the fitting: nfit = %d\n",psr[0].nFit);
      if (psr[0].nFit>2)
	{
	  epoch[n]=(double)centreMJD;
	  dmVal[n]=(double)psr[0].param[param_dm].val[0];
	  dmE[n]=(double)psr[0].param[param_dm].err[0];	  
	  id[n] = i;
	  printf("%g %g %g %d\n",epoch[n],dmVal[n],dmE[n],id[n]);
	  fprintf(fout,"%g %g %g %d\n",epoch[n],dmVal[n],dmE[n],id[n]);
	  n++;
	}
    }  
  fclose(fout);
  printf("Completed getting DM values\n");
  // Redo the initial fit
  psr[0].nJumps=0;
  readParfile(psr,parFile,timFile,1); /* Load the parameters       */
  readTimfile(psr,timFile,1); /* Load the arrival times    */
      for (j=0;j<MAX_PARAMS;j++)
	{
	  for (k=0;k<psr[0].param[j].aSize;k++)
	    psr[0].param[j].fitFlag[k] = 0;
	}
      for (j=0;j<psr[0].nJumps;j++)
	psr[0].fitJump[j]=0;

  preProcess(psr,1,argc,argv);      
  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      if (i==0) doFit(psr,*npsr,0);   /* Do the fitting     */
      else textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
    }
  printf("About to plot\n");

  // Now plot the results
  doPlot(epoch,dmVal,dmE,id,n,psr);

  return 0;
}

void doPlot(double *epoch,double *dmVal,double *dmE,int *id,int n,pulsar *psr)
{
  int i,j,k;
  int removePoly=0;
  float earliest,latest;
  float px[n],py[n],pe1[n],pe2[n];
  float ix[MAX_OBSN],iy[MAX_OBSN],iix[MAX_OBSN];
  float plotX[MAX_OBSN],plotY[MAX_OBSN],plotY1[MAX_OBSN],plotY2[MAX_OBSN],plotE[MAX_OBSN];
  float dmx[2],dmy[2];
  double dm0;
  int nplot=0;
  float minv,maxv;
int nInterp=0;
  float minx,maxx,miny,maxy;
  float minxr,maxxr,minyr,maxyr;
  float rx[MAX_OBSN],ry[MAX_OBSN],ry1[MAX_OBSN],ry2[MAX_OBSN];
  int plotType=1;
  float mx,my;
  float first20,first50,first10;
  char key='a';
  
  int dmCurvePoly=1;


  dm0 = (double)psr[0].param[param_dm].val[0];
  // Find earliest and latest point in the .tim file
  earliest = latest = (float)psr[0].obsn[0].sat;

  for (i=0;i<psr[0].nobs;i++)
    {
      rx[i] = mjd2year((double)psr[0].obsn[i].sat);
      ry[i] = (float)psr[0].obsn[i].residual*1.0e6;
      ry1[i] = ry[i]-(float)psr[0].obsn[i].toaErr;
      ry2[i] = ry[i]+(float)psr[0].obsn[i].toaErr;

      if (earliest > (float)psr[0].obsn[i].sat) earliest = (float)psr[0].obsn[i].sat;
      if (latest < (float)psr[0].obsn[i].sat) latest = (float)psr[0].obsn[i].sat;
    }
  earliest-=0.5;
  latest+=0.5;
  //  earliest = epoch[0];
  //  latest = epoch[n-1];

  cpgbeg(0,"/xs",1,1);
  cpgsch(1.4); cpgsfs(2);
  cpgask(0);
  do 
    {
      cpgeras();
      for (i=0;i<n;i++)
	{
	  px[i] = (float)mjd2year(epoch[i]);
	  py[i] = (float)dmVal[i];
	  pe1[i] = (float)(dmVal[i]-dmE[i]);
	  pe2[i] = (float)(dmVal[i]+dmE[i]);
	}
      minx = TKfindMin_f(px,n);
      maxx = TKfindMax_f(px,n);
      miny = TKfindMin_f(py,n);
      maxy = TKfindMax_f(py,n);

      //      minyr = TKfindMin_f(ry,psr[0].nobs);
      //      maxyr = TKfindMax_f(ry,psr[0].nobs);

      // Plot residuals
      // Plot the different frequencies           
      selectData(psr,rx,ry,0,800,plotX,plotY,plotE,&nplot);
      if (removePoly!=0 && nplot>0) TKremovePoly_f(plotX,plotY,nplot,removePoly);
      minv = TKfindMin_f(plotY,nplot);
      maxv = TKfindMax_f(plotY,nplot);
      selectData(psr,rx,ry,800,1600,plotX,plotY,plotE,&nplot);
      if (removePoly!=0 && nplot>0) TKremovePoly_f(plotX,plotY,nplot,removePoly);
      if (minv > TKfindMin_f(plotY,nplot)) minv = TKfindMin_f(plotY,nplot);
      if (maxv < TKfindMax_f(plotY,nplot)) maxv = TKfindMax_f(plotY,nplot);
      selectData(psr,rx,ry,1600,5000,plotX,plotY,plotE,&nplot);
      if (removePoly!=0 && nplot>0) TKremovePoly_f(plotX,plotY,nplot,removePoly);
      if (minv > TKfindMin_f(plotY,nplot)) minv = TKfindMin_f(plotY,nplot);
      if (maxv < TKfindMax_f(plotY,nplot)) maxv = TKfindMax_f(plotY,nplot);

      cpgsvp(0.1,0.95,0.70,0.95);
      minyr = minv;
      maxyr = maxv;
      cpgswin(minx-(maxx-minx)*0.1,maxx+(maxx-minx)*0.1,
	      minyr-(maxyr-minyr)*0.05,maxyr+(maxyr-minyr)*0.05);      
      cpgbox("BCTS",0,0,"BCNTS",0,0);
      cpglab("","Residual (\\gms)","");


      selectData(psr,rx,ry,0,800,plotX,plotY,plotE,&nplot);
      if (removePoly!=0 && nplot > 0) TKremovePoly_f(plotX,plotY,nplot,removePoly);
      for (i=0;i<nplot;i++)
	{
	  plotY1[i] = plotY[i]-plotE[i];
	  plotY2[i] = plotY[i]+plotE[i];
	}
      first50 = plotY[0];
      cpgsci(2); cpgsch(0.5); cpgpt(nplot,plotX,plotY,4); cpgsch(1.4); cpgerry(nplot,plotX,plotY1,plotY2,1); cpgsci(1);

      selectData(psr,rx,ry,800,1600,plotX,plotY,plotE,&nplot);
      if (removePoly!=0 && nplot > 0) TKremovePoly_f(plotX,plotY,nplot,removePoly);
      for (i=0;i<nplot;i++)
	{
	  plotY1[i] = plotY[i]-plotE[i];
	  plotY2[i] = plotY[i]+plotE[i];
	}
      first20 = plotY[0];
      cpgsci(3); cpgsch(0.5); cpgpt(nplot,plotX,plotY,4); cpgsch(1.4); cpgerry(nplot,plotX,plotY1,plotY2,1); cpgsci(1);

      selectData(psr,rx,ry,1600,5000,plotX,plotY,plotE,&nplot);
      if (removePoly!=0 && nplot > 0) TKremovePoly_f(plotX,plotY,nplot,removePoly);
      for (i=0;i<nplot;i++)
	{
	  plotY1[i] = plotY[i]-plotE[i];
	  plotY2[i] = plotY[i]+plotE[i];
	}
      first10 = plotY[0];
      cpgsci(4); cpgsch(0.5); cpgpt(nplot,plotX,plotY,4); cpgsch(1.4); cpgerry(nplot,plotX,plotY1,plotY2,1); cpgsci(1);

      // If available: overplot the DM curve
      if (nInterp > 0)
	{
	  float tt;
	  float ty[nInterp];
	  double freqf;
	  double mean;

	  freqf=1400e6;
	  for (i=0;i<nInterp;i++)
	    ty[i] = (iy[i]-dm0)/DM_CONST/1.0e-12/freqf/freqf*1.0e6;
	  TKremovePoly_f(ix,ty,nInterp,dmCurvePoly);
	  tt = ty[0];
	  for (i=0;i<nInterp;i++)
	    ty[i]=ty[i]-tt+first20;
	  cpgsci(3); cpgline(nInterp,ix,ty); cpgsci(1);

	  freqf=650e6; 
	  for (i=0;i<nInterp;i++)
	    ty[i] = (iy[i]-dm0)/DM_CONST/1.0e-12/freqf/freqf*1.0e6;
	  TKremovePoly_f(ix,ty,nInterp,dmCurvePoly);
	  tt = ty[0];
	  for (i=0;i<nInterp;i++)
	    ty[i]=ty[i]-tt+first50;
	  cpgsci(2); cpgline(nInterp,ix,ty); cpgsci(1);
	  
	}


      // Plot DM
      cpgsvp(0.1,0.95,0.15,0.70);
      cpgswin(minx-(maxx-minx)*0.1,maxx+(maxx-minx)*0.1,
      	     miny-(maxy-miny)*0.1,maxy+(maxy-miny)*0.1);
      cpgbox("BCNTS",0,0,"BCNTS",0,0);
      cpglab("Year","DM (cm\\u-3\\dpc)","");
      dmx[0] = minx-(maxx-minx)*0.1;
      dmx[1] = maxx+(maxx-minx)*0.1;
      dmy[0] = dmy[1] = (float)psr[0].param[param_dm].val[0];
      cpgsls(4); cpgline(2,dmx,dmy); cpgsls(1);

      cpgpt(n,px,py,4);
      cpgerry(n,px,pe1,pe2,1);
      if (nInterp>0)
	cpgline(nInterp,ix,iy);
      printf("--------------------------\n");
      printf("0 = don't remove polynomial from residuals\n");
      printf("1 = remove mean from residuals\n");
      printf("2 = remove linear term from residuals\n");
      printf("3 = remove quadratic term from residuals\n");
      printf("! = remove mean from DM curve\n");
      printf("@ = remove linear term from DM curve\n");
      printf("# = remove quadratic term from DM curve\n");
      
      printf("f = interpolate and smooth\n");
      printf("q = quit\n");
      printf("s = save DM file\n");
      printf("--------------------------\n");
      cpgcurs(&mx,&my,&key);
      if (key=='0')
	removePoly=0;
      else if (key=='1')
	removePoly=1;
      else if (key=='2')
	removePoly=2;
      else if (key=='3')
	removePoly=3;
      else if (key=='!')
	dmCurvePoly=1;
      else if (key=='@')
	dmCurvePoly=2;
      else if (key=='#')
	dmCurvePoly=3;
	else if (key=='f')
	{
	  double sampleTime;
	  double smoothTime;
	  double dt;
	  double tl,bl;
	  int smooth=2;
	  char type;
	  int nread;
	  //	  printf("Enter: sample time (e.g. 5 days) ");
	  //	  scanf("%lf",&sampleTime);
	  sampleTime = 1;
	  printf("Enter: smooth time (days) type\n");
	  printf("type = 'e' for exponential smoother, 'g' for Gaussian smoother \n");
	  printf("Example: 30 g\n");
	  printf("Value: ");
	  nread = scanf("%lf %c",&smoothTime,&type);
	  if (nread==2)
	    {
	      if (type=='e') smooth=1;
	      else if (type=='c') smooth=2;
	      else {
		printf("Unknown smooth type: using Gaussian smoother\n");
		smooth=2;
	      }
	    }
	  nInterp = (int)((latest-earliest)/sampleTime+0.5);
	  printf("nInterp = %d\n",nInterp);
	  for (i=0;i<nInterp;i++)
	    {
	      iix[i] = earliest + sampleTime*i;
	      tl = 0.0;
	      bl = 0.0;
	      for (j=0;j<n;j++)
		{
		  dt = iix[i]-epoch[j];
		  if (smooth==1) // Exponential smoother
		    {
		      tl += dmVal[j]/dmE[j]/dmE[j]*exp(-fabs(dt/smoothTime));
		      bl += 1.0/dmE[j]/dmE[j]*exp(-fabs(dt/smoothTime));
		    }
		  else if (smooth==2) // Gaussian smoother
		    {
		      tl += dmVal[j]/dmE[j]/dmE[j]*exp(-0.5*pow(dt/smoothTime,2));
		      bl += 1.0/dmE[j]/dmE[j]*exp(-0.5*pow(dt/smoothTime,2));
		    }
		}
	      ix[i] = mjd2year(iix[i]);
	      iy[i] = tl/bl;	      
	    }	 	  
	}
      else if (key=='s')
	{
	  if (nInterp==0)
	    {
	      printf("Must use 'f' option to fit a smooth line first\n");
	    }
	  else
	    {
	      FILE *fout;
	      char fname[100];
	      sprintf(fname,"%s.dm",psr[0].name);
	      fout = fopen(fname,"w");
	      for (i=0;i<nInterp;i++)
		fprintf(fout,"%.4f %.8f\n",iix[i],iy[i]);
	      printf("Output written to %s\n",fname);
	      fclose(fout);
	    }
	}
    } while (key!='q');
  cpgend();
}

double mjd2year(double mjd)
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


      return  retYr+(retDay+(day-(int)day))/365.25;

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

void selectData(pulsar *psr,float *rx,float *ry,double f1,double f2,float *plotX,float *plotY,float *plotE,int *nplot)
{
  int i;

  *nplot=0;
  for (i=0;i<psr[0].nobs;i++)
    {
      if (psr[0].obsn[i].deleted==0 && psr[0].obsn[i].freq > f1 && psr[0].obsn[i].freq <= f2)
	{
	  plotX[*nplot] = rx[i];
	  plotY[*nplot] = ry[i];
	  plotE[*nplot] = psr[0].obsn[i].toaErr;
	  (*nplot)++;
	}
    }     

}
char * plugVersionCheck = TEMPO2_h_VER;
