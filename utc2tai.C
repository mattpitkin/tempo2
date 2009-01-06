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
