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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tempo2.h"
#include "read_fortran.h"

#include "ifteph.h"

/* TDB is at mean sea level -- what about height of observatory? Backer & Hellings (1986) */

/* ******************************************** */
/* tt2tb                                       */
/* Author:  G. Hobbs (05 May 2003)              */
/* Purpose: Converts TOAs in TAI format         */
/*          to barycentric dynamical time       */
/* Inputs:                                      */
/*                                              */
/* Outputs: Fills correctionTT_TDB             */
/*                                              */
/* Notes: based on tdbread.f and interp.f       */
/*        calculates ctatv                      */
/* Changes:                                     */
/*    RTE: -now computes TDB or TCB depending on */
/*         options. Name changed to reflect this.*/
/*         -can use IF or FB time ephemeris     */
/*         -computes Teph also                  */
/* ******************************************** */


void init_ifte();
double IF_deltaT(longdouble mjd_tt);
double FB_deltaT(longdouble mjd_tt);


void tt2tb(pulsar *psr,int npsr)
{
  int i, p,k;
  longdouble mjd_tt, mjd_teph;
  double deltaT=0.0, obsTerm, earthVel[3];
  double deltaTDot, obsTermDot, earthVelDot[3];
  int first=1;
  const char *CVS_verNum = "$Revision: 1.12 $";

  if (displayCVSversion == 1) CVSdisplayVersion("tt2tdb.C","tt2tdb()",CVS_verNum);

  init_ifte();
  for (p=0;p<npsr;p++)
  {
    for (i=0;i<psr[p].nobs;i++)
    {
      if (psr[p].obsn[i].clockCorr==0 && strcmp(psr[p].obsn[i].telID,"STL")!=0)
	psr[p].obsn[i].correctionTT_TB = 0.0;
      else
      {
	if (strcmp(psr[p].obsn[i].telID,"STL")==0)
	  mjd_tt = psr[p].obsn[i].sat;
	else
	  mjd_tt = psr[p].obsn[i].sat + getCorrectionTT(psr[p].obsn+i)/SECDAY;
	/* Evaluate the time ephemeris and its derivative */
	if (psr[p].timeEphemeris == IF99_TIMEEPH)
	  {
	    deltaT = IF_deltaT(mjd_tt);
	  }
	if (psr[p].timeEphemeris == FB90_TIMEEPH)
	  {
	    deltaT = FB_deltaT(mjd_tt) - IFTE_TEPH0; /* FB code returns Teph-TDB not TT-TDB */
	  }
	/* Get term involving observatory position (if known) */

	/*	for (k=0;k<3;k++)
		varray[k] = -psr[p].obsn[i].earthMoonBary_earth[k+3];
		vectorsum(earthVel, 
		psr[p].obsn[i].earthMoonBary_ssb+3,
		varray);   */
	for (k=0;k<3;k++)
	  earthVel[k] = psr[p].obsn[i].earth_ssb[k+3];
	//logdbg("i %d Vearth = %g %g %g",i,earthVel[0],earthVel[1],earthVel[2]);
	obsTerm = dotproduct(earthVel, psr[p].obsn[i].observatory_earth)
	  /  (1.0-IFTE_LC) ;
	/* both vectors were meant to be defined in Ephemeris units ...
	 This correction is neglible but worth doing so we can test with K=2*/
	if (psr[p].units == SI_UNITS)
	  obsTerm /= (double)(IFTE_K*IFTE_K); // both were in SI 
	else
	  obsTerm /= (double)IFTE_K;  // obs_earth was in SI
	//	printf("obsTerm = %g\n",(double)obsTerm);
	/* Compute Teph : Irwin & Fukushima (1999) eq 13*/	
	/* Note, DeltaT(Teph0) ~ -2x10^-14 s so is neglected */
	psr[p].obsn[i].correctionTT_Teph = IFTE_TEPH0 + obsTerm
	  + deltaT / (1.0-IFTE_LC);
	/* Compute TCB or TDB as requested */
	if (psr[p].units == TDB_UNITS)
	  {
	    psr[p].obsn[i].correctionTT_TB = 
	      psr[p].obsn[i].correctionTT_Teph - (1? 0.0 : IFTE_TEPH0); 
	  }
	else
	  {
	    psr[p].obsn[i].correctionTT_TB = 
	      IFTE_KM1 * (double)(mjd_tt-IFTE_MJD0)*86400.0/* linear drift term */
	      + IFTE_K * (psr[p].obsn[i].correctionTT_Teph-(longdouble)IFTE_TEPH0);
	  }
	/* Compute dTB/dTT for correct frequency transfer if needed */
	if (psr[p].dilateFreq==1)
	{
	  if (first==1)
	    {
	      //	      init_ifte();
	      first=0;

	    }
	  mjd_teph = mjd_tt + psr[p].obsn[i].correctionTT_Teph/86400.0;

	  deltaTDot = IFTE_DeltaTDot(2400000.0+(int)mjd_teph, 
				     0.5+(mjd_teph-(int)mjd_teph));
	  // XXX seems deltaTDot is coming out too big
	  IFTE_get_vEDot(2400000.0+(int)mjd_teph, 
			 0.5+(mjd_teph-(int)mjd_teph), earthVelDot);
	  vectorscale(earthVelDot, 1.0/86400.0);
	  obsTermDot = 
	    (dotproduct(earthVelDot, psr[p].obsn[i].observatory_earth)
	     + dotproduct(earthVel, psr[p].obsn[i].siteVel))
	    /  (1.0-IFTE_LC) ;
	  if (psr[p].units == SI_UNITS)
	    obsTerm /= (double) (IFTE_K*IFTE_K); // both were in SI 
	  else
	    obsTerm /= (double)IFTE_K;  // obs_earth was in SI
	  
	  if (psr[p].units == TDB_UNITS)	    
	    psr[p].obsn[i].einsteinRate = 
	      1.0 + obsTermDot + deltaTDot/(1.0-IFTE_LC);
	  else
	    {
	      psr[p].obsn[i].einsteinRate = 
		IFTE_K * (1.0 + obsTermDot + deltaTDot/(1.0-IFTE_LC));
	      //	  obsTermDot = 
	      //	    (dotproduct(earthVelDot, psr[p].obsn[i].observatory_earth)
	      //	     + dotproduct(earthVel, psr[p].obsn[i].siteVel))
	      //	    /  (1.0-IFTE_LC) ;

							   
	    }

	logdbg("corrn TT_Teph %llg",psr[p].obsn[i].correctionTT_Teph);
	logdbg("corrn TT_TB %llg",psr[p].obsn[i].correctionTT_TB);
	logdbg("corrn einstin %llg",psr[p].obsn[i].einsteinRate-1);
// 	  if (!first) printf("%g %g %g %g %g Einstein\n", (double)mjd_tt,
// 		 deltaT, deltaTDot, obsTerm, obsTermDot);
//	  if (!first)
//	    printf("ER-1 %lg\n", psr[p].obsn[i].einsteinRate-1.0);
	}	  
      }
    }
  }
  IFTE_close_file();
}


void init_ifte()
{
  //  static int first=1;
  //  if (first)
  //  {
  //    first = 0;
    char fname[1024];
    strcpy(fname,getenv(TEMPO2_ENVIRON));
    strcat(fname,IFTEPH_FILE);
    IFTE_init(fname);
    //  } 
}

double IF_deltaT(longdouble mjd_tt)
{
  //  init_ifte();
  return IFTE_DeltaT(2400000.0+(int)mjd_tt, 0.5+(mjd_tt-(int)mjd_tt))*86400.0;
  /* Note, DeltaT(Teph0) ~ -2x10^-14 s so is neglected */
}

double FB_deltaT(longdouble mjd_tt)
{
  double ctatv;    /* output TDB-TDT */
  longdouble tdt; /* jd1 + jd2  (Julian date) */
  static int tdbnrl=-1;
  int nr,nrecl,j,k,np,nv;
  double jda,jdb,tdbd1,tdbd2,t[2];
  int tdbdt,tdbncf;
  char fname[MAX_FILELEN];
  double dna,dt1,temp,pc[18],tc,twot,dummy;
  static double buf[16];
  int l;

  /* Initialise the TDB-TDT file (tdbinit.f) */
  /* Set up the TDB-TDT ephemeris file for reading */
  nrecl = 4; /* If recl is in bytes (for Sun -- what about LINUX ????) */

  strcpy(fname,getenv(TEMPO2_ENVIRON));
  strcat(fname,TDBTDT_FILE);

  /* Now do the calculations */
  open_file(fname); /* Open a Fortran made file for reading in C */
  tdbd1 = read_double();
  tdbd2 = read_double();
  tdbdt = read_int();
  tdbncf = read_int();
  dummy = read_double();
  dummy = read_double();
  dummy = read_double();
  dummy = read_double();
  dummy = read_double();
	      
  /* Use the corrected TT time and convert to Julian date */
  tdt = mjd_tt + 2400000.5; 
  if (tdt - (int)tdt >= 0.5)
  {
    jda = (int)tdt + 0.5;
    jdb = tdt - (int)tdt - 0.5;
  }
  else
  {
    jda = (int)tdt - 0.5;
    jdb = tdt - (int)tdt + 0.5;
  }
  nr = (int)((jda-tdbd1)/tdbdt)+2; 
  if (nr < 1 || tdt > tdbd2) 
  {
    printf("ERROR [CLK4]: Date %.10f out of range of TDB-TDT table\n",(double)tdt);
    exit(1);
  }
  if (nr!=tdbnrl)
  {
    tdbnrl = nr;
    /* MUST JUMP TO RECORD NR IN THE FILE */
    for (j=0;j<nr-1;j++)
    {
      for (k=0;k<tdbncf;k++)
	buf[k] = read_double();
    }
  }
  t[0] = ((jda-((nr-2)*tdbdt+tdbd1))+jdb)/tdbdt; /* Fraction within record */
  t[1] = 1.0; /* Unused */
	      
	      /* Interpolation: call interp(buf,t,tdbncf,1,  1,   1,   ctatv) */
  np = 2;
  nv = 3;
  twot = 0.0; 
	      
  pc[0] = 1.0; pc[1]=0.0;
	      
  dna = 1.0;
  dt1 = (int)(t[0]);
  temp = dna * t[0];
  l = (int)(temp - dt1)+1;
  tc = 2.0*(fortran_mod(temp,1.0)+dt1)-1.0;
	      
  if (tc != pc[1])
  {
    np = 2;
    nv = 3; 
    pc[1] = tc;
    twot = tc+tc;	  
  } 
  if (np<tdbncf)
  {
    for (k=np+1;k<=tdbncf;k++)
      pc[k-1] = twot*pc[k-2]-pc[k-3];
    np = tdbncf;
  }
  ctatv = 0.0;
  for (j=tdbncf;j>=1;j--)
    ctatv = ctatv +pc[j-1]*buf[j-1];
  close_file();

  return ctatv;
}


