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
#include <string.h>
#include "tempo2.h"
//#include "tempo2Util.h"
//#include "tempo2pred.h"
//#include "tempo2pred_int.h"

char TEMPO2_ENVIRON[MAX_STRLEN]="TEMPO2";
char TEMPO2_ERROR[MAX_STRLEN]="";

int MAX_PSR   = MAX_PSR_VAL;    /* Maximum number of pulsars to fit simultaneously  */
int MAX_OBSN  = MAX_OBSN_VAL;
double ECLIPTIC_OBLIQUITY = ECLIPTIC_OBLIQUITY_VAL;
int forceGlobalFit = 0;
int veryFast = 0;
int displayCVSversion = 0;
char tempo2MachineType[MAX_FILELEN] = "";

char dcmFile[MAX_FILELEN]="NULL";
char covarFuncFile[MAX_FILELEN]="NULL";
char tempo2_plug_path[32][MAX_STRLEN];
int tempo2_plug_path_len=0;

// From choleskyRoutines.h ... "some global variables that Ryan is still using for diagnostic purposes"
double FCALPHA, WNLEVEL, EXPSMOOTH, UPW, NFIT, FCFINAL;




#define MAX_FUNCTIONS 1024 /* Maximum functions in tempo2 */

void extra_delays(pulsar *psr,int npsr)
{  
  const char *CVS_verNum = "$Revision: 1.23 $";
  if (displayCVSversion == 1) CVSdisplayVersion((char *)"global.C",(char *)"extra_delays()",CVS_verNum);

  calculate_bclt(psr,npsr);/* 3. Calculate bclt  */
  /*  shapiro_delay(psr,npsr); */ /* 1. Calculate the Shapiro delay */
  /* dm_delays(psr,npsr); */    /* 2. Extra dispersion measure delays */  
}

void clock_corrections(pulsar *psr,int npsr)
{  
  const char *CVS_verNum = "$Revision: 1.23 $";
  if (displayCVSversion == 1) CVSdisplayVersion((char *)"global.C",(char *)"clock_corrections()",CVS_verNum);

  logdbg("Calling toa2utc");
  toa2utc(psr,npsr);        /* 1. UTC(Observatory) -> UTC(NIST) */
  logdbg("Calling tai2ut1");
//   utc2tai(psr,npsr);     /* 2. UTC(NIST) -> TAI              */
  tai2ut1(psr,npsr);        /* 3. TAI -> UT1                    */
//   tai2tt(psr,npsr);      /* 4. TAI -> TT                     */
  logdbg("Calling tt2tb");
  tt2tb(psr,npsr);          /* 5. Rough estimate of TT-TB (+-2.2 microsec) */
  logdbg("Done clock corrections");
}

void ephemeris_routines(pulsar *psr,int npsr)
{ 
  const char *CVS_verNum = "$Revision: 1.23 $";
  if (displayCVSversion == 1) CVSdisplayVersion((char *)"global.C",(char *)"ephemeris_routines()",CVS_verNum);

  logtchk("call vectorPulsar()");
  vectorPulsar(psr,npsr);   /* 1. Form a vector pointing at the pulsar */
  logtchk("call readEphemeris()");
  if (psr[0].useCalceph == 0)
    readEphemeris(psr,npsr,0);/* 2. Read the ephemeris */
  else
    readEphemeris_calceph(psr,npsr);
  logtchk("call get_obsCoord()");
  get_obsCoord(psr,npsr);   /* 3. Get Coordinate of observatory relative to Earth's centre */
  logtchk("call tt2tb()");
  tt2tb(psr,npsr);          /* Observatory/time-dependent part of TT-TB */
  logtchk("call readEphemeris()");
  if (psr[0].useCalceph == 0)
    readEphemeris(psr,npsr,0);  /* Re-evaluate ephemeris with correct TB */ 
  else
    readEphemeris_calceph(psr,npsr);
}

void formBatsAll(pulsar *psr,int npsr)
{
  const char *CVS_verNum = "$Revision: 1.23 $";
  if (displayCVSversion == 1) CVSdisplayVersion((char *)"global.C",(char *)"formBatsAll()",CVS_verNum);

  logtchk("enter formBatsAll()");
  logdbg("Calling clock corrections");
  logtchk("call clock_corrections()");
  clock_corrections(psr,npsr);          /* Clock corrections  ... */  
  logdbg("Reading ephemeris routines");
  logtchk("call ephemeris_routines()");
  ephemeris_routines(psr,npsr);         /* Ephemeris routines ... */
  logdbg("Reading extra delays");
  logtchk("call extra_delays()");
  extra_delays(psr,npsr);               /* Other time delays  ... */
  logtchk("call formBats()");
  formBats(psr,npsr);                   /* Form Barycentric arrival times */
  logtchk("call secularMotion()");
  secularMotion(psr,npsr); 
  logtchk("exit formBatsAll()");
}

// Only recalculate that which is likely
// to change if psr position has been altered.
void updateBatsAll(pulsar *psr, int npsr)
{
  const char *CVS_verNum = "$Revision: 1.23 $";
  if (displayCVSversion == 1) CVSdisplayVersion((char *)"global.C",(char *)"updateBatsAll()",CVS_verNum);

  vectorPulsar(psr, npsr);
  calculate_bclt(psr, npsr);
  formBats(psr, npsr);
  secularMotion(psr,npsr); 
}


// Display the version number if it hasn't already been displayed
void CVSdisplayVersion(char *file,char *func,const char *verNum)
{
  static char alreadyFunc[MAX_FUNCTIONS][64],alreadyFile[MAX_FUNCTIONS][64];
  static int counter=0;
  int i,have=0;
  for (i=0;i<counter;i++)
    {
      if (strcmp(alreadyFunc[i],func)==0 &&
	  strcmp(alreadyFile[i],file)==0)
	{
	  have=1;
	  break;
	}
    }
  if (have==0)
    {
      if (counter==0) // First go display tempo2.h version
	printf("[TEMPO2 VERSION:] [%-20.20s] [%-20.20s] [%-20.20s]\n","tempo2.h","",TEMPO2_h_VER);
      printf("[TEMPO2 VERSION:] [%-20.20s] [%-20.20s] [%-20.20s]\n",file,func,verNum);
      strcpy(alreadyFunc[counter],func);
      strcpy(alreadyFile[counter],file);
      counter++;
    }
}


