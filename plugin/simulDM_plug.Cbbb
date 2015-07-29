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

using namespace std;

#define MAX_DM 1000 // Maximum number of DM measurements possible

void help() /* Display help */
{
  /* This function should contain usage information about the plugin which should (in general) be accessed */
  /* by the user pressing 'h'                                                                              */
}

void doplugin(pulsar *psr);
void doPlot(pulsar *psr,int *highFreq_id,int *lowFreq_id,double *dmVal,
	    double *timeVal,int nDM,double *dmValE);

/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i;
  double globalParameter;

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: simulDM\n");
  printf("Author:              G. Hobbs\n");
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
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[0],argv[i+1]); 
	  strcpy(timFile[0],argv[i+2]);
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
      else textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
    }

  doplugin(psr);
  
  return 0;
}

void doplugin(pulsar *psr)
{
  int i,j;
  double f10_l = 2000; // Low bound of the 10cm band  (MHz)
  double f10_h = 4000; // High bound of the 10cm band (MHz)
  double f50_l = 500;  // Low bound of the 50cm band (MHz)
  double f50_h = 900;  // High bound of the 50cm band (MHz)
  double simulTime = 120; // Time between points to indicate a simultaneous point (sec)
  int highFreq_id[MAX_DM];
  int lowFreq_id[MAX_DM];
  double dmVal[MAX_DM];
  double dmValE[MAX_DM];
  double timeVal[MAX_DM];
  int nDM = 0;


  // Find simultaneous 10 and 50CM points
  for (i=0;i<psr[0].nobs;i++)
    {
      if (psr[0].obsn[i].freq > f10_l && psr[0].obsn[i].freq < f10_h)
	{
	  // Found a 10cm point
	  // Now see if there is a simultaneous 50cm point
	  for (j=0;j<psr[0].nobs;j++)
	    {
	      if (psr[0].obsn[j].freq > f50_l && psr[0].obsn[j].freq < f50_h)
		{
		  // Have a 50cm point ... is it simultaneous?
		  if (fabs(psr[0].obsn[i].sat - psr[0].obsn[j].sat) < simulTime/SECDAY)
		    {
		      // Have a simultaneous 10cm and 50cm point
		      highFreq_id[nDM] = i;
		      lowFreq_id[nDM] = j;
		      dmVal[nDM] = 
			DM_CONST*(psr[0].obsn[j].residual-psr[0].obsn[i].residual)/
			(pow(psr[0].obsn[j].freq,-2)-pow(psr[0].obsn[i].freq,-2));
		      dmValE[nDM] = DM_CONST*(sqrt(pow(psr[0].obsn[i].toaErr*1.0e-6,2) + pow(psr[0].obsn[j].toaErr*1.0e-6,2)))/(pow(psr[0].obsn[j].freq,-2)-pow(psr[0].obsn[i].freq,-2));
		      timeVal[nDM] = (double)(psr[0].obsn[i].sat+psr[0].obsn[j].sat)*0.5;
		      nDM++;
		    }
		}
	    }
	}
    }
  printf("Number of dispersion measures calculated = %d\n",nDM);
  for (i=0;i<nDM;i++)
    printf("%d %d %g %g %g\n",highFreq_id[i],lowFreq_id[i],timeVal[i],dmVal[i],dmValE[i]);

  // Now make a plot of the DM values
  doPlot(psr,highFreq_id,lowFreq_id,dmVal,timeVal,nDM,dmValE);

}

void doPlot(pulsar *psr,int *highFreq_id,int *lowFreq_id,double *dmVal,
	    double *timeVal,int nDM,double *dmValE)
{
  float minx,maxx,miny,maxy;
  float xval[MAX_OBSN],yval[MAX_OBSN],yerr1[MAX_OBSN],yerr2[MAX_OBSN];
  float yval2[MAX_OBSN];
  double freq;
  int i;

  cpgbeg(0,"?",1,1);
  cpgslw(2);
  cpgsch(1.4);
  cpgsfs(2);

  for (i=0;i<nDM;i++)
    {
      xval[i] = (float)timeVal[i];
      yval[i] = (float)dmVal[i];
      yerr1[i] = (float)(dmVal[i]-dmValE[i]);
      yerr2[i] = (float)(dmVal[i]+dmValE[i]);
    }
  
  miny = TKfindMin_f(yval,nDM);
  maxy = TKfindMax_f(yval,nDM);
  minx = TKfindMin_f(xval,nDM);
  maxx = TKfindMax_f(xval,nDM);

  //  cpgenv(minx-0.1*(maxx-minx),maxx+0.1*(maxx-minx),
  //	 miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny),0,1);
  cpgsvp(0.1,0.95,0.5,0.9);
  cpgswin(minx-0.1*(maxx-minx),maxx+0.1*(maxx-minx),
	  miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny));
  cpgsch(1);
  cpgbox("BCTS",0,0,"BCNTS",0,0);
  cpglab("","\\gDDM",psr[0].name);
  cpgsch(1.4);
  cpgpt(nDM,xval,yval,9);
  cpgerry(nDM,xval,yerr1,yerr2,1);

  // Should also plot dDM caused by the Solar Wind
  for (i=0;i<psr[0].nobs;i++)
    {
      xval[i] = psr[0].obsn[i].sat;
      freq = psr[0].obsn[i].freq;
      yval[i] = (miny+maxy)*0.5+psr[0].obsn[i].tdis2*freq*freq*DM_CONST;
    }
  cpgsci(2); cpgsls(4); cpgline(psr[0].nobs,xval,yval); cpgsci(1); cpgsls(1);

  // Should also plot induced residuals at 10,20,50cm
  for (i=0;i<nDM;i++)
    {
      xval[i] = (float)timeVal[i];
      freq = 1400.0; // MHz
      yval[i] = (float)dmVal[i]/DM_CONST/freq/freq/1e-6;      
    }
  for (i=0;i<nDM;i++)
    {
      freq = 3000.0; // MHz
      yval2[i] = (float)dmVal[i]/DM_CONST/freq/freq/1e-6;      
    }
  miny = TKfindMin_f(yval,nDM);
  maxy = TKfindMax_f(yval,nDM);
  if (miny > TKfindMin_f(yval2,nDM)) miny = TKfindMin_f(yval2,nDM);
  if (maxy < TKfindMax_f(yval2,nDM)) maxy = TKfindMax_f(yval2,nDM);
  minx = TKfindMin_f(xval,nDM);
  maxx = TKfindMax_f(xval,nDM);

  cpgsvp(0.1,0.95,0.15,0.5);
  cpgswin(minx-0.1*(maxx-minx),maxx+0.1*(maxx-minx),
	  miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny));
  cpgsch(1);
  cpgbox("BNCTS",0,0,"BCNTS",0,0);
  cpglab("MJD","Delay (\\gms)","");
  cpgsch(1.4);
  cpgsci(3); cpgline(nDM,xval,yval);
  cpgsci(4); cpgline(nDM,xval,yval2);

  cpgend();  
}
char * plugVersionCheck = TEMPO2_h_VER;
