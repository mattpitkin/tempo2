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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"


#if 0 // not needed any more -- done generically in clkcorr.C


/* ******************************************** */
/* utc2tai                                      */
/* Author:  G. Hobbs (05 May 2003)              */
/* Purpose: Converts TOAs in UTC format         */
/*          to atomic time (TAI)                */
/* Inputs:  obsn - structure of observations    */
/*          nObs number of lines in .tim file   */
/* Outputs: Fills correction_utc_tai            */
/*                                              */
/* Notes:   Difference between UTC(NIST) and    */
/*          TAI is an integral number of leap   */
/*          seconds.                            */
/*                                              */
/*          This routine is very different to   */
/*          the original TEMPO. Must understand */
/*          why ...                             */
/*                                              */
/* Changes:                                     */
/* ******************************************** */
void utc2tai(pulsar *psr,int npsr)
{
  int i,j,p;
  int leap_mjd[MAX_LEAPSEC],nleap=0,nread;
  double timediff;
  /*  double a1utcf; */
  FILE *fin;
  char fname[MAX_FILELEN];
  const char *CVS_verNum = "$Revision: 1.5 $";

  if (displayCVSversion == 1) CVSdisplayVersion("utc2tai.C","utc2tai()",CVS_verNum);
  
  /* Read the leap second file and store in arrays */
  strcpy(fname,getenv(TEMPO2_ENVIRON));
  strcat(fname,LEAPSECOND_FILE);
  if (!(fin = fopen(fname,"r")))
    {
      printf("ERROR [FILE6]: Unable to open the leap second file: %s\n",LEAPSECOND_FILE);
      exit(1);
    }
  while (!feof(fin))
    {
      nread = fscanf(fin,"%d",&leap_mjd[nleap]);
      if (nread == 1) /* Read successfully */
	{
	  nleap++;
	  if (nleap > MAX_LEAPSEC)
	    {printf("ERROR [CLK3], must increase size of MAX_LEAPSEC\n"); exit(1);}
	}
    }
  fclose(fin);
  for (p=0;p<npsr;p++)
    {
      /* In 1972, when the leap second procedure was introduced, TAI = UTC+10[sec] */            

      for (i=0;i<psr[p].nobs;i++)
	{
	  if (psr[p].obsn[i].clockCorr==0)
	    {
	      psr[p].obsn[i].correctionUTC_TAI = 0.0;
	    }
	  else
	    {
	      /*	  timediff=10.0343817;  */
	      timediff = 10.0; 
	      if (psr[p].obsn[i].sat < 41317.0)
		{
		  printf("Warning: cannot calculated UTC -> TAI prior to 1972\n");
		  psr[p].obsn[i].correctionUTC_TAI = 0.0;
		  /*	      psr[p].obsn[i].a1utcf = 0.0; */
		}
	      else
		{
		  for (j=0;j<nleap;j++)
		    {
		      if (psr[p].obsn[i].sat>leap_mjd[j])
			timediff+=1.0; /* Add one second */
		    }
		  psr[p].obsn[i].correctionUTC_TAI = timediff;	   
		  /* Where is this 32.15 coming from? */
		  /* psr[p].obsn[i].correctionUTC_TAI = timediff+32.15; */
		  /*		  psr[p].obsn[i].a1utcf = timediff; */
		}
	    }
	}
    }
}

#endif
