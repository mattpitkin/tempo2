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

// redwards code to get the Earth orientation parameters.
// These are meant to be taken from a file holding the EOP C04 series
// of the IERS. This is somewhat important since not all the series
// can be linearly interpolated

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include "dynarr.h"
#include "tempo2.h"

typedef struct
{
  double mjd;
  double xp;
  double yp;
  double dut1;
} EOPSample;


void
load_EOP(DynamicArray *EOPsamples,char *eopcFile)
{
  char fname[1024], line[1024], mjd_s[1024], xp_s[1024], yp_s[1024], dut1_s[1024];
  int year;
  EOPSample newSample;
  FILE *f;
  int iline,idummy;
  int format=0;
  const char *CVS_verNum = "$Revision: 1.10 $";

  if (displayCVSversion == 1) CVSdisplayVersion("eop.C","load_EOP()",CVS_verNum);

  // init array
  DynamicArray_init(EOPsamples, sizeof(EOPSample));

  // open file
  //  sprintf(fname, "%s/%s", getenv("TEMPO2"), EOPC04_FILE);
  sprintf(fname, "%s/%s", getenv("TEMPO2"), eopcFile);
  f = fopen(fname, "r");
  if (!f)
  {
    fprintf(stderr, "Fatal Error: Unable to open file %s for reading: %s\n",
	      fname, strerror(errno));
    exit(1);
  }

  // parse file
  
  for (iline=1; fgets(line, 1024, f)!=NULL; iline++)
  {
    // ignore header lines, ancient history (pre-PSR), etc
    if (sscanf(line, "%d", &year)==1 && year > 1966 && year < 3000)
    {
      if (format==0)
	{
	  // break it up, ugly fixed format fortran I/O
	  strncpy(mjd_s, line+17, 5); mjd_s[5]='\0';
	  strncpy(xp_s, line+22, 9); xp_s[9]='\0';
	  strncpy(yp_s, line+31, 9); yp_s[9]='\0';
	  strncpy(dut1_s, line+40, 10); dut1_s[10]='\0';
	  // put into sample, in radians and seconds
	  if (sscanf(mjd_s, "%lf", &newSample.mjd)!=1 ||
	      sscanf(xp_s, "%lf", &newSample.xp)!=1 ||
	      sscanf(yp_s, "%lf", &newSample.yp)!=1 ||
	      sscanf(dut1_s, "%lf", &newSample.dut1)!=1)
	    {
	      fprintf(stderr, "Error parsing line %d of %s", iline, fname);
	      exit(1);
	    }
	}
      else if (format==1)
	{
	  if (sscanf(line,"%d %d %d %lf %lf %lf %lf",&year,&idummy,&idummy,&newSample.mjd,&newSample.xp,&newSample.yp,&newSample.dut1)!=7)
	    {
	      fprintf(stderr, "Error parsing line %d of %s", iline, fname);
	      exit(1);
	    }
	}
      newSample.xp *= M_PI/(180.0*60.0*60.0);
      newSample.yp *= M_PI/(180.0*60.0*60.0);
      DynamicArray_push_back(EOPsamples, &newSample);
    }
    else {  // Check which format we are using
      if (strstr(line,"FORMAT(2X,I4,2X,A4,I3,2X,I5,2F9.6,F10.7,2X,F10.7,2X,2F9.6)")!=NULL)
	{
	  format=0;
	  printf("Warning: using old format for $TEMPO2/eopc04_IAU2000.62-now - should update\n");
	}
      else if (strstr(line,"FORMAT(3(I4),I7,2(F11.6),2(F12.7),2(F11.6),2(F11.6),2(F11.7),2F12.6)")!=NULL)
	{
	  format=1;
	  logdbg("Using new format for $TEMPO2/eopc04_IAU2000.62-now");
	}
    }
  }
  fclose(f); 
}

void
get_EOP(double mjd, double *xp, double *yp, double *dut1, 
	      double *dut1dot,
	int dispWarnings,char *eopcFile)
{
  static int first = 1;
  static DynamicArray EOPsamples;
  EOPSample *samp;
  int proxy = -1, isamp;
  double f;
  double leap=0.0;
  // load EOP series first
  if (first)
  {
    load_EOP(&EOPsamples,eopcFile);
    first = 0;
  }
  samp =  (EOPSample *)EOPsamples.data;

  // check for out of bounds
  if (samp[0].mjd > mjd)
    proxy = 0;
  else if (samp[EOPsamples.nelem-1].mjd <= mjd)
    proxy = EOPsamples.nelem-1;
  if (proxy >=0)
  {
    *xp = samp[proxy].xp;
    *yp = samp[proxy].yp;
    *dut1 = samp[proxy].dut1;
    {
      char msg2[1000];
      sprintf(msg2,"%.2f. Using %.1lf instead",mjd,samp[proxy].mjd);
      if( dispWarnings != 2 )
        displayMsg(1,"EPH1","No Earth orientation parameters for MJD ",msg2,dispWarnings);      
    }
    // GHOBBS: Added next line (and removed the 'return' after) otherwise dut1dot never gets set!
    mjd = samp[proxy].mjd;
    //    return;
  }

    /* find first sample to fall after requested time */
  for (isamp = 0; 
       isamp < (int)EOPsamples.nelem  && samp[isamp].mjd <= mjd; 
       isamp++)
    ;

  // cope with leap second ... take it off second point as jump happens
  // right AT second point
  if (samp[isamp].dut1 - samp[isamp-1].dut1 > 0.5)
    leap = 1.0;
  else if (samp[isamp].dut1 - samp[isamp-1].dut1 < -0.5)
    leap = -1.0;
  else
    leap = 0.0;

  *dut1dot = (samp[isamp].dut1 - leap - samp[isamp-1].dut1) / 86400.0;

  /* interpolate */
  f = (mjd - samp[isamp-1].mjd) / (samp[isamp].mjd - samp[isamp-1].mjd);
  *xp = samp[isamp-1].xp + f*(samp[isamp].xp - samp[isamp-1].xp);
  *yp = samp[isamp-1].yp + f*(samp[isamp].yp - samp[isamp-1].yp);
  *dut1 = samp[isamp-1].dut1 + f*(samp[isamp].dut1 - leap - samp[isamp-1].dut1);
 
 
  return;
}

