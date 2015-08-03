#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"

#define USE_NEW_CLK_CORR


double toa2utc_nist(double toa,char *clockFile);
double linearInterpolate(double x1,double y1,double x2,double y2,double x3);
double convertTOA(double mjd,char *clks);

/* ******************************************** */
/* toa2utc                                      */
/* Author:  G. Hobbs (05 May 2003)              */
/* Purpose: Converts TOAs which have been       */
/*          obtained using an observatory time  */
/*          standard and convert to a uniform   */
/*          time standard.                      */
/* Inputs:  obsn - structure of observations    */
/*          nObs number of lines in .tim file   */
/* Outputs: Calculates                          */
/*          psr[p].obsn[i].correctionUTC        */
/*                                              */
/* Notes:                                       */
/*                                              */
/* Changes:                                     */
/* ******************************************** */

#ifdef USE_NEW_CLK_CORR
void toa2utc(pulsar *psr,int npsr)
{
  int i,p;
  const char *CVS_verNum = "$Revision: 1.7 $";

  if (displayCVSversion == 1) CVSdisplayVersion("toa2utc.C","toa2utc__1()",CVS_verNum);

  for (p=0;p<npsr;p++)
      for (i=0;i<psr[p].nobs;i++)
	{
	  if (psr[p].obsn[i].clockCorr!=0 && psr[p].obsn[i].clockCorr!=2)
	  {
	    //	    logdbg("Getting clock corrections for %d",i);
	    getClockCorrections(psr[p].obsn+i, psr[p].clockFromOverride,
				psr[p].clock, psr[p].noWarnings);
	  }
	  // else clock is presumed to be TT/TDB/TCG already
	}
}
#else

void toa2utc(pulsar *psr,int npsr)
{
  int i,j,p,found=0;
  const char *CVS_verNum = "$Revision: 1.7 $";

  if (displayCVSversion == 1) CVSdisplayVersion("toa2utc.C","toa2utc__2()",CVS_verNum);

  for (p=0;p<npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)
	{
	  found=0;
	  /* Find correction between UTC(PKS) and UTC(NIST) (output in seconds) */
	  if (psr[p].obsn[i].clockCorr==0 || psr[p].obsn[i].clockCorr==2) /* If no clock 
									     corrections being applied 
									     or already UTC */
	    psr[p].obsn[i].correctionUTC = 0.0;   /* then set correction to zero */
	  else
	    {
	      /* Determine clock correction file list */
	      for (j=0;j<psr[p].obsn[i].nFlags;j++)
		{
		  if (strcmp(psr[p].obsn[i].flagID[j],"-c")==0)
		    {
		      psr[p].obsn[i].correctionUTC = convertTOA((double)psr[p].obsn[i].sat,
								psr[p].obsn[i].flagVal[j]);
		      found=1;
		      /* printf("Correction = %g\n",(float)psr[p].obsn[i].correctionUTC); */
		      break;
		    }
		}
	      if (found==0) /* Using old technique for determining clock correction */
		psr[p].obsn[i].correctionUTC = toa2utc_nist((double)psr[p].obsn[i].sat,
							    psr[p].OBSERVATORY_CLOCK_2_UTC_NIST); 
	    }
	  printf("Correction %.15g\n",  psr[p].obsn[i].correctionUTC);
	}
    }
}

/* Converts a TOA to UTC */
double convertTOA(double mjd,char *clks)
{
  char *token;
  int count=0,i,j,endit,clk,lti,gti;
  char fname[100][100];
  char first[100],store2[100];
  char store[100];
  static char readFile[10][100];  /* List of clock files already loaded into memory */
  static double readMJD[10][MAX_CLKCORR];
  static double readClock[10][MAX_CLKCORR];
  static int clockNum[10];
  static int fileNum=0;
  static int warned=0;
  FILE *fin;
  double correction=0.0;
  double bestclkcorr;
  char line[1000];
  int found=0;

  strcpy(store,clks);

  /* Split up the clock correction file list to determine the files required */
  token=strtok(clks,"_"); strcpy(first,token);
  while ( (token=strtok(NULL, "_")) != NULL)
    {            
      strcpy(store2,getenv(TEMPO2_ENVIRON));
      sprintf(fname[count],"%s/clock/%s_2_%s.clk",store2,first,token);
      strcpy(first,token);
      count++;
    }
  /* Test whether all these files have already been opened, if not open them */
  for (i=0;i<count;i++)
    {
      for (j=0;j<fileNum;j++)
	{
	  if (strcmp(readFile[j],fname[i])==0)
	    {
	      found=1;
	      break;
	    }
	}
      if (found==0) /* Must read the file */
	{
	  if (!(fin = fopen(fname[i],"r")))
	    {
	      printf("Unable to open clock correction file >%s<\n",fname[i]);
	      exit(1);
	    }
	  strcpy(readFile[fileNum],fname[i]);

	  clockNum[fileNum] = 0;
	  /* Read header line */
	  fgets(line,1000,fin);
	  /* Now read the file */
	  while (!feof(fin))
	    {
	      if (fscanf(fin,"%lf %lf",&readMJD[fileNum][clockNum[fileNum]],
			 &readClock[fileNum][clockNum[fileNum]])==2)
		clockNum[fileNum]++;		
	    }
	  fclose(fin);
	  fileNum++;
	}
    }
  

  /* Now calculate the clock correction */
  correction=0.0;
  for (i=0;i<count;i++)
    {
      clk=-1;
      for (j=0;j<fileNum;j++)
	{
	  if (strcmp(readFile[j],fname[i])==0)
	    {
	      clk=j;
	      break;
	    }
	}
      if (clk==-1) {printf("Big problem in toa2utc.C\n"); exit(1);}
      lti = 0;
      gti = 0;
      /* Now have the array of the clock correction required */
      for (j=0;j<clockNum[clk];j++)
	{
	  if (mjd < readMJD[clk][j] &&
	      fabs(readMJD[clk][j]-mjd) < fabs(readMJD[clk][gti]-mjd))
	    gti = j;
	  else if (mjd > readMJD[clk][j] &&
		   fabs(mjd - readMJD[clk][j]) < fabs(mjd - readMJD[clk][lti]))
	    lti = j; 
	}
      /* Give warning if TOA is greater than last MJD in clock file */
      if (readMJD[clk][gti] == 0.0 && warned == 0)
	{
	  warned=1;
	  printf("Warning [CLK1]: MJD later than last entry in time.dat.  Imprecise clock corrections will be applied\n");
	}
      /* Do a linear interpolation between these two points */
      bestclkcorr = linearInterpolate(readMJD[clk][lti],readClock[clk][lti],readMJD[clk][gti],
				      readClock[clk][gti],mjd);
      /* Convert from microseconds into seconds */
      bestclkcorr /= 1.0e6;
      correction+=bestclkcorr;
    }

  strcpy(clks,store);

  return correction;  
}

/* ********************************************* */
/* Converts a TOA to UTC(NIST)                   */
/* Author: G. Hobbs (ATNF)                       */
/* Version 1: June 05 2003                       */
/*                                               */
/* The subroutine reads the time.dat file        */
/* containing clock corrections between          */
/* UTC(PKS) and UTC(NIST)                        */
/* The clock correction in seconds is            */
/* returned                                      */
/* A linear interpolation is used to             */
/* interpolate between points in the time.dat    */
/* file.                                         */
/*                                               */
/* The first time this function is called        */
/* the time.dat file is read and arrays are      */
/* filled.                                       */
/*                                               */
/* Important variables:                          */ 
/*                                               */
/* clkcorrMJD -  MJD of lines in clock correction*/
/*               file                            */
/* clkcorrVal -  Clock correction in microsec    */
/* clkcorrSite - Observatory site for clock      */
/*               correction                      */
/* lti - array element with MJD just less than   */
/*       the TOA                                 */
/* gti - array element with MJD just greater than*/
/*       the TOA                                 */
/*                                               */
/* Precision issues: TO BE COMPLETED             */
/* ********************************************* */

double toa2utc_nist(double toa,char *clockFile)
{
  static double clkcorrMJD[MAX_CLKCORR],clkcorrVal[MAX_CLKCORR];
  static int clkcorrSite[MAX_CLKCORR],clkcorrN,time=0;
  static int warned = 0;
  double bestclkcorr;
  int nread,i,lti,gti;
  char dummy[100],line[1000];
  FILE *fin;

  /* Read the clock corrections between UTC(PARKES) and UTC(NIST) if not already done*/
  if (time==0)
    {
      time=1; /* Don't do see on any subsequent calls to toa2utc_nist */
      if (!(fin = fopen(clockFile,"r")))
	{
	  printf("ERROR [FILE5]: Unable to open file %s\n",clockFile);
	  exit(1);
	}
      while (!feof(fin))
	{
	  /* Read information from the clock correction file */
	  /* Note, the fourth column contains a comment which is ignored here */
	  fgets(line,1000,fin);
	  nread = sscanf(line,"%lf %lf %d %s\n",&clkcorrMJD[clkcorrN],
			 &clkcorrVal[clkcorrN],&clkcorrSite[clkcorrN],dummy);
	  /* Attempting to fix crazy fixed format in time.dat (see newsrc.f, line 290) */
	  /* MUST FIND A BETTER WAY OF DOING THIS */
	  if (line[16]==' ')  clkcorrVal[clkcorrN] *= -1; 
	  if (nread==4) clkcorrN++;
	  if (clkcorrN > MAX_CLKCORR)
	    {
	      printf("ERROR [CLK2] Must increase MAX_CLKCORR\n");
	      exit(1);
	    }
	}
      fclose(fin);
    }
  /* Now have TOA */
  /* Must find the closest times to the TOA in the clock file and interpolate */
  /* to find the clock correction                                             */
  /* Have SITE hardcoded to 7 -- PARKES !!! MUST CHANGE */
  lti = 0;
  gti = 0;
  for (i=0;i<clkcorrN;i++)
    {
      if (clkcorrSite[i] == 7 && toa < clkcorrMJD[i] && 
	  fabs(clkcorrMJD[i]-toa) < fabs(clkcorrMJD[gti]-toa))
	gti  = i;
      else if (clkcorrSite[i] == 7 && toa > clkcorrMJD[i] && 
	       fabs(toa - clkcorrMJD[i]) < fabs(toa - clkcorrMJD[lti]))
	lti = i;
    }
  /* Give warning if TOA is greater than last MJD in time.dat file */
  if (clkcorrVal[gti] == 0.0 && warned==0)
    {
      warned=1;
      printf("Warning [CLK1]: MJD later than last entry in time.dat.  Imprecise clock corrections will be applied\n");
    }
  /* Do a linear interpolation between these two points */
  bestclkcorr = linearInterpolate(clkcorrMJD[lti],clkcorrVal[lti],clkcorrMJD[gti],
				  clkcorrVal[gti],toa);
  printf("Hmmm %.10lg\n", clkcorrVal[lti]*1.0e-6);
  /* Convert from microseconds into seconds */
  bestclkcorr /= 1.0e6;

  return -bestclkcorr;

}
/* ****************************************************************** */
/* G. Hobbs (ATNF)                                                    */
/* Does a linear interpolation between (x1,y1) and (x2,y2) to predict */
/* the value of y3 at x3.  The value y3 is returned                   */
/*                                                                    */
/* This routine is based on the following                             */
/* y1 = m*x1 + c                                                      */
/* y2 = m*x2 + c                                                      */
/* m = (y2-y1)/(x2-x1); c = y1 - m*x1                                 */
/* therefore y3 = (y2-y1)/(x2-x1)*(x3-x1) + y1                        */
/*                                                                    */
/* Precision comments:                                                */
/* ****************************************************************** */

double linearInterpolate(double x1,double y1,double x2,double y2,double x3)
{
  double y3;
  y3 = ((y2-y1)/(x2-x1))*(x3-x1) + y1;
  return y3;
}

#endif
