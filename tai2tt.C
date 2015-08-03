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
  const char *CVS_verNum = "$Revision: 1.5 $";

  if (displayCVSversion == 1) CVSdisplayVersion("tai2tt.C","tai2tt()",CVS_verNum);

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
