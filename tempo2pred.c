//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russel Edwards

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


#include "tempo2pred.h"
#include "tempo2pred_int.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

unsigned tempo2_verbose = 1;

int
T2Predictor_Read(T2Predictor *t2p, char *fname)
{
  FILE *f = fopen(fname, "r");
  int ret;
  if (!f)
    return -1;
  ret = T2Predictor_FRead(t2p, f);
  fclose(f);
  return ret;
}
  
int
T2Predictor_FRead(T2Predictor *t2p, FILE *f)
{

  // determine the kind of file we're dealing with by trial and error
  long fpos = ftell(f);
  if (ChebyModelSet_Read(&t2p->modelset.cheby, f)==0)
  {
    t2p->kind = Cheby;
    return 0;
  }
  else
    fseek(f, fpos, SEEK_SET); // go back to where we started
 
  if (T1PolycoSet_Read(&t2p->modelset.t1, f)==0)
  {
    t2p->kind = T1;
    return 0;
  }
  else
    fseek(f, fpos, SEEK_SET); // go back to where we started

    return -1;
}
  
void
T2Predictor_Write(const T2Predictor *t2p, char *fname)
{
  FILE *f = fopen(fname, "w");
  if (!f)
  {
    fprintf(stderr, "Fatal error: T2Predictor_Write unable to open file %s\n",
	    fname);
    exit(1);
  }
  T2Predictor_FWrite(t2p, f);
  fclose(f);
}

void 
T2Predictor_FWrite(const T2Predictor *t2p, FILE *f)
{
  switch (t2p->kind)
  {
  case Cheby:
    ChebyModelSet_Write(&t2p->modelset.cheby, f);
    break;
  case T1:
    T1PolycoSet_Write(&t2p->modelset.t1, f);
    break;
  }
}

void 
T2Predictor_Init(T2Predictor *t2p)
{
  t2p->kind = None;
  memset (&t2p->modelset, 0, sizeof(t2p->modelset));
}

void
T2Predictor_Copy(T2Predictor *into_t2p, const T2Predictor *from_t2p)
{
  T2Predictor_Destroy(into_t2p);
  T2Predictor_Insert(into_t2p, from_t2p);
}

int
T2Predictor_Insert(T2Predictor *into_t2p, const T2Predictor *from_t2p)
{
  if (into_t2p->kind == None)
    into_t2p->kind = from_t2p->kind;

  if (into_t2p->kind != from_t2p->kind) {
    fprintf (stderr, "T2Predictor_Insert: not of the same kind!\n");
    return -1;
  }

  switch (into_t2p->kind)
  {
  case Cheby:
    ChebyModelSet_Insert(&into_t2p->modelset.cheby,
			 &from_t2p->modelset.cheby);
    break;
  case T1:
    fprintf (stderr, "T2Predictor_Insert: not implemented for T1 polyco!\n");
    exit (-1);
    break;
  }
}

void T2Predictor_Keep(T2Predictor* t2p, unsigned nmjd, const long double* mjd)
{
  switch (t2p->kind)
  {
  case Cheby:
    ChebyModelSet_Keep(&t2p->modelset.cheby, nmjd, mjd);
    break;
  case T1:
    fprintf (stderr, "T2Predictor_Keep: not implemented for T1 polyco!\n");
    exit (-1);
    break;
  }
}

void
T2Predictor_Destroy(T2Predictor *t2p)
{
  switch (t2p->kind)
  {
  case Cheby:
    ChebyModelSet_Destroy(&t2p->modelset.cheby);
    break;
  case T1:
    T1PolycoSet_Destroy(&t2p->modelset.t1);
    break;
  }
  t2p->kind = None;
}

  /* Information */
char * 
T2Predictor_GetPSRName(T2Predictor *t2p)
{
  char *ret;
  switch (t2p->kind)
  {
  case Cheby:
    ret = t2p->modelset.cheby.segments[0].psrname;
    break;
  case T1:
    ret = t2p->modelset.t1.segments[0].psrname;
    break;
  }
  return ret;
}
char * T2Predictor_GetSiteName(T2Predictor *t2p)
{
  char *ret; 
  switch (t2p->kind)
  {
  case Cheby:
    ret = t2p->modelset.cheby.segments[0].sitename;
    break;
  case T1:
    ret = t2p->modelset.t1.segments[0].sitename;
    break;
  }
  return ret;
}
long double T2Predictor_GetStartMJD(T2Predictor *t2p)
{
  long double ret;
  switch (t2p->kind)
  {
  case Cheby:
    ret = t2p->modelset.cheby.segments[0].mjd_start;
    break;
  case T1:
    ret = t2p->modelset.t1.segments[0].mjd_mid - 
      t2p->modelset.t1.segments[0].span/2880.0L;
    break;
  }
  return ret;
}
long double T2Predictor_GetEndMJD(T2Predictor *t2p)
{
  long double ret;
  switch (t2p->kind)
  {
  case Cheby:
    ret = t2p->modelset.cheby.segments
      [t2p->modelset.cheby.nsegments-1].mjd_end;
    break;
  case T1:
    ret = t2p->modelset.t1.segments
      [t2p->modelset.t1.nsegments-1].mjd_mid + 
      t2p->modelset.t1.segments[t2p->modelset.t1.nsegments-1].span/2880.0L;
    break;
  }
  return ret;
} 
long double T2Predictor_GetStartFreq(T2Predictor *t2p) // MHz
{
  long double ret;
  switch (t2p->kind)
  {
  case Cheby:
    ret = t2p->modelset.cheby.segments[0].freq_start;    
    break;
  case T1:
    ret = t2p->modelset.t1.segments[0].frequency_obs;    
    break;
  }
  return ret;
}
long double T2Predictor_GetEndFreq(T2Predictor *t2p)  // MHz
{
  long double ret;
  switch (t2p->kind)
  {
  case Cheby:
    ret = t2p->modelset.cheby.segments[0].freq_end;    
    break;
  case T1:
    ret = t2p->modelset.t1.segments[0].frequency_obs;    
    break;
  }
  return ret;
}
T2PredictorKind T2Predictor_Kind(T2Predictor *t2p)
{
  return t2p->kind;
}

/* Prediction */
long double T2Predictor_GetPhase(const T2Predictor *t2p, long double mjd,
				 long double freq) // freq in MHz
{
  long double ret;
  switch (t2p->kind)
  {
  case Cheby:
    ret = ChebyModelSet_GetPhase(&t2p->modelset.cheby, mjd, freq);    
    break;
  case T1:
    ret = T1PolycoSet_GetPhase(&t2p->modelset.t1, mjd, freq); 
    break;
  }
  return ret;
}

long double T2Predictor_GetFrequency(const T2Predictor *t2p, long double mjd,
				     long double freq) // freq in MHz
{
  long double ret;
  switch (t2p->kind)
  {
  case Cheby:
    ret = ChebyModelSet_GetFrequency(&t2p->modelset.cheby, mjd, freq);    
    break;
  case T1:
    ret = T1PolycoSet_GetFrequency(&t2p->modelset.t1, mjd, freq);    
    break;
  }
  return ret;
}


int T2Predictor_GetPlan(char *filename,
			long double mjd_start,
			long double mjd_end,
			long double step, // seconds
			long double freq, // MHz
			// remaining arguments are returned
			long double *phase0,
			int *nsegments,
			long double *pulse_frequencies)
{
  T2Predictor_GetPlan_Ext(filename, mjd_start, mjd_end, step, freq,
			  NULL, NULL, phase0, nsegments, pulse_frequencies);
}

int T2Predictor_GetPlan_Ext(char *filename,
			long double mjd_start,
			long double mjd_end,
			long double step, // seconds
			long double freq, // MHz
			// remaining arguments are returned
			char *psrname, char *sitename,
			long double *phase0,
			int *nsegments,
			long double *pulse_frequencies)
{
  int iseg, ret;
  long double lastphase, phase, mjd, lastmjd;
  T2Predictor pred;

  if ((ret=T2Predictor_Read(&pred, filename)))
    return ret;
  if (psrname)
    strcpy(psrname, T2Predictor_GetPSRName(&pred));
  if (sitename)
    strcpy(sitename, T2Predictor_GetSiteName(&pred));

  *nsegments = (int)ceill((mjd_end-mjd_start)*86400.0L/step);
  *phase0 = lastphase = T2Predictor_GetPhase(&pred, mjd_start, freq);
  lastmjd = mjd_start;

  for (iseg=0; iseg < *nsegments; iseg++)
  {
    mjd = mjd_start +(iseg+1)*step/86400.0L;
    if (mjd > mjd_end)
      mjd = mjd_end;
    phase = T2Predictor_GetPhase(&pred, mjd, freq);
    pulse_frequencies[iseg] = (phase-lastphase)/(mjd-lastmjd)/86400.0L;
    lastphase = phase;
    lastmjd = mjd;
  }

  T2Predictor_Destroy(&pred);

  return 0;
}

