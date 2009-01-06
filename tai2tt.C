#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"

#if 0 // not needed any more, done generically in clkcorr.C
/* ******************************************** */
/* utc2tai                                      */
/* Author:  G. Hobbs (Dec 23, 2003)             */
/* Purpose: Converts atomic time (TAI)  to      */
/*          terrestrial time (TT) [also called  */
/*          TDT]                                */
/*                                              */
/* Inputs:  obsn - structure of observations    */
/*          nObs number of lines in .tim file   */
/* Outputs: Fills correction_tai_tt             */
/*                                              */
/* Notes: Conversion is a constant offset       */
/*        of 32.184 seconds                     */
/*                                              */
/* Changes:                                     */
/* ******************************************** */
void tai2tt(pulsar *psr,int npsr)
{
  int i,p;
  for (p=0;p<npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)
	{
	  if (psr[p].obsn[i].clockCorr==0)
	    psr[p].obsn[i].correctionTAI_TT = 0.0;
	  else
	    psr[p].obsn[i].correctionTAI_TT = 32.184;
	}
    }
}
#endif
