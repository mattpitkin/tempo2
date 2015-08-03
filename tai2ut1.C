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

double ut1red(double mjd,int warnings);


/* ******************************************** */
/* utc2tai                                      */
/* Author:  G. Hobbs (05 May 2003)              */
/* Purpose: Converts TAI to UT1                 */
/*                                              */
/* Inputs:  obsn - structure of observations    */
/*          nObs number of lines in .tim file   */
/*                                              */
/*                                              */
/*                                              */
/* Changes:                                     */
/* ******************************************** */
void tai2ut1(pulsar *psr,int npsr)
{
  int p,i;
  const char *CVS_verNum = "$Revision: 1.6 $";

  if (displayCVSversion == 1) CVSdisplayVersion("tai2ut1.C","tai2ut1()",CVS_verNum);

  for (p=0;p<npsr;p++)
    {
      if (psr[p].t2cMethod==T2C_TEMPO)
       for (i=0;i<psr[p].nobs;i++)
	{
	  if (psr[p].obsn[i].clockCorr!=0)
	    {
	      /* Original TEMPO calculation */
	      /*	      psr[p].obsn[i].correctionUT1 = ut1red((double)(psr[p].obsn[i].sat+
							   psr[p].obsn[i].correctionUTC/SECDAY+
							   psr[p].obsn[i].correctionUTC_TAI/SECDAY),
							   psr[p].obsn[i].a1utcf);*/
	      
	      /* redwards note, the frame of the argument to the IERS does
		 not seem to be specified. We will simply use the whatever
		 version of TT was specified, which could be out from UTC by up
		 to a couple of minutes within our lifetimes. TAI lays between
		 the two. 
		 The maximum rate of change of the EOP are as follows:
                         x, y   : 5 mas/day
                         UT1-UTC:  4 ms/day  (excluding leap seconds!!)
                 For x & y, 1 ns is 9 mas, so no problem.
		 For UT1-UTC, 1 ns is 330 microsec, so we must have the
		 argument right to within 7000 seconds ... no problem! */

	      /* Where does the 0.0343817 come from?? See end of ut1red ... */
	      /* redwards -- it is due to tempo's internal use of A.1,
		 which differs from TAI by this many seconds */

	      /* redwards, changed ut1red to return UT1-TAI. So now we have
		 UT1-TOA = UT1-TAI + (TAI-TOA) */


	      psr[p].obsn[i].correctionUT1 =  
		ut1red((double)(psr[p].obsn[i].sat),
		       psr[p].noWarnings)
			+ getCorrection(psr[p].obsn+i, 
					psr[p].clockFromOverride,
					"TAI", psr[p].noWarnings);
		      
	    }
	  else
	    psr[p].obsn[i].correctionUT1 = 0.0;
	}
    }
}


/* ******************************************** */
/* ut1red                                       */
/* Author:  G. Hobbs (01 Sept 2003)             */
/* Purpose: Reads UT1 values from file (ut1.dat)*/
/*                                              */
/*                                              */
/* Notes: based on ut1red.f                     */
/*        Not checked very well                 */
/* Changes:                                     */
/*  Redwards eliminated use of UTC & A.1, now it*/
/* returns UT1-TAI as provided (-ve) in the table     */
/* (previously it took AT1-UTC, & knew AT1-TAI, */
/*  and returned UT1-UTC, uggggh!               */
/* ******************************************** */
double ut1red(double mjd, int warnings)
{
  static int read=0;
  static double mjd1,mjd2;
  static int mjdVal[1000];
  static int entry[10000];
  static int count2=0,count=0;
  double t,s,f2,y2[2],y1[2],ut1;
  double units = 1.0e-4;
  double tab[4];
  int nread,i,iint=5,it,itsav=-1,j,nr;
  char line[1000],fname[MAX_FILELEN];
  FILE *fin;

  if (read==0) /* Read data file */
    {
      strcpy(fname,getenv(TEMPO2_ENVIRON));
      strcat(fname,UT1_FILE);
	
      fin = fopen(fname,"r");
      if (!fin)
      {
	fprintf(stderr, "Could not open UT1 file: %s\n", fname);
	exit(1);
      }
      fgets(line,1000,fin); /* Read two header lines (I hope that these don't change!!!!) */
      fgets(line,1000,fin); /* ... no checks are done ... */
      nread=1;
      while ( (nread==1) && (fgets(line,1000,fin)!=NULL) )
	{
          // Chop off everything past char 57 to avoid the extra 'item
          // count' value being read in as either a UT1 value or an MJD.
          line[57]='\0';
	  nread = sscanf(line,"%d %d %d %d %d %d %d",&mjdVal[count],
                  &entry[count2+0], &entry[count2+1], &entry[count2+2],
                  &entry[count2+3], &entry[count2+4], &entry[count2+5]);
	  if (nread>1) 
	    {
	      count++;
	      for (i=0;i<nread-1;i++)
		{
		  if (entry[count2]!=0)
		    count2++;
		  else
		    nread=-1;
		}	  
              if (nread>0) { nread=1; }
	    }
          else
            {
              nread=-1;
            }
	  if (count>999) /* Too many lines? */
	    {
	      printf("ERROR: Too many lines in ut1 file -- increase array size in tai2tdb\n");
	      exit(1);
	    }
	}
      fclose(fin);
      mjd1 = mjdVal[0]+iint; /* iint = interval in days between tabular values */
      mjd2 = mjdVal[count-1]-iint; /* Not sure about this + and - iint  *******************/ 
    }

    read=1;

    /* Is MJD within range of table */
    if (mjd>mjd1 && mjd<mjd2) 
      {
	/* Calculate interpolation times and value of tabular points */
	t = (mjd-mjd1)/iint;   
	it = (int)t;
	t = t - it;
	s=1.0-t;
	if (it!=itsav)
	  {
	    for (i=0;i<4;i++)
	      {
		j=it+i;
		tab[i]=entry[j];
	      }
	    for (i=0;i<2;i++)
	      {
		nr = i+1;
		f2 = 1.0/6.0 * (tab[nr+1]+tab[nr-1]);
		y1[i] = 4.0/3.0 * tab[nr] -f2;
		y2[i] = -1.0/3.0 * tab[nr]+f2;		
	      }
	    itsav = it;
	  }
	/* second difference interpolation */
	ut1 = (t*(y1[1]+t*t*y2[1]) + s*(y1[0]+s*s*y2[0]))*units;
	/*	printf("UT1 for %f = %f\n",mjd,ut1); */
      }
    else
      {
	char msg[1000],msg2[1000];
	sprintf(msg,"MJD is outside of UT1 table range ");
	sprintf(msg2,"(MJD %f)",mjd);
	displayMsg(1,"CLK10",msg,msg2,warnings);

	/* Extrapolate using last value in table */
	if (mjd<=mjd1)
	  ut1 = entry[0]*units;
	else
	  ut1 = entry[count2-1]*units;	
      }

//     ut1 = atut-0.0343817 - ut1; /* Table TAI-UT1 */
    return -ut1;
    
}
