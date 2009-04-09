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
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
 
/* Tempo1 polyco functions destined for tempo2pred library (C linkage) */

long double T1Polyco_GetPhase(T1Polyco *t1p, long double mjd, long double freq)
{
  long double dt, arg, phase, spin_freq;
  int ic;

  /* Evaluate polynomial in phase and frequency */
  dt = (mjd-t1p->mjd_mid)*1440;
  printf("WARNIGN WARNGING: Silly update in polyco read\n");
  printf("Have %.20Lg %.20Lg\n",mjd,t1p->mjd_mid);
  printf("Diff = %.15Lg\n",mjd-t1p->mjd_mid);

  phase = t1p->reference_phase + dt*60*t1p->frequency_psr_0 + t1p->coeff[0];
  spin_freq = t1p->frequency_psr_0 + t1p->coeff[1]/60;
  arg = dt;
  for (ic=1; ic<t1p->ncoeff; ic++)
  {
    phase += t1p->coeff[ic]*arg;
    if (ic < t1p->ncoeff-1)
      spin_freq += (ic+1)*t1p->coeff[ic+1]*arg/60;
    arg *= dt;
  }
  /* compute dispersion phase delay in observatory frame */
  phase -= spin_freq * t1p->dm/2.41e-4L / (1.0L+t1p->doppler*1.0e-4L) * 
    ( 1.0L/(freq*freq) - 1.0L/(t1p->frequency_obs*t1p->frequency_obs));
  return phase;
}

long double T1Polyco_GetFrequency(T1Polyco *t1p, long double mjd, 
				  long double freq)
{
  long double dt, arg, spin_freq;
  int ic;

  dt = (mjd-t1p->mjd_mid)*1440;
  spin_freq = t1p->frequency_psr_0 + t1p->coeff[1]/60;
  arg = dt;

  for (ic=2; ic<t1p->ncoeff; ic++)
  {
    spin_freq += ic*t1p->coeff[ic]*arg/60;
    arg *= dt;
  }

  return spin_freq;
}

void T1Polyco_Write(T1Polyco *t1p, FILE *fout)
{
  int i;

  if (t1p->psrname[0]=='J') 
    fprintf(fout,"%-10.10s ",t1p->psrname+1);
  else 
    fprintf(fout,"%-10.10s ",t1p->psrname);
  fprintf(fout,"%9.9s",t1p->date_string);
  fprintf(fout,"%11.11s",t1p->utc_string);
  fprintf(fout,"%20.11Lf",t1p->mjd_mid);
  fprintf(fout,"%21.6lf ",t1p->dm);
  fprintf(fout,"%6.3lf",t1p->doppler);
  fprintf(fout,"%7.3lf",t1p->log10rms);
  fprintf(fout,"\n");
  fprintf(fout,"%20.6Lf",t1p->reference_phase);
  fprintf(fout,"%18.12Lf",t1p->frequency_psr_0);
  fprintf(fout,"%5s",t1p->sitename);
  fprintf(fout,"%5d",t1p->span);
  fprintf(fout,"%5d",t1p->ncoeff);
  fprintf(fout,"%10.3lf",t1p->frequency_obs);
  fprintf(fout,"%16.14lf",t1p->binary_phase);
  fprintf(fout,"\n");
  for (i=0;i<t1p->ncoeff;i++)
  {
    fprintf(fout,"%25.17Le", t1p->coeff[i]);
    if ((i+1)%3==0) fprintf(fout,"\n");
  }
}

/* helper functions for T1Polyco_Read */
void T1P_grabString(char *str, int istart, int nchar, char *out)
{
  char grabbed[64];

  strncpy(grabbed, str+istart, nchar);
  grabbed[nchar] = '\0';
  sscanf(grabbed, "%s", out);// remove leading/trailing whitespace
}

long double T1P_grabLongDouble(char *str, int istart, int nchar)
{
  char grabbed[64];
  long double res;
  char *c;

  T1P_grabString(str, istart, nchar, grabbed);
 
  /* change D to e for exponential notation*/
  for (c=grabbed; *c; c++)
    if (*c == 'D')
      *c = 'e';

  sscanf(grabbed, "%Lf", &res);

  return res;
}

long double T1P_grabInt(char *str, int istart, int nchar)
{
  char grabbed[64];
  int res;

  T1P_grabString(str, istart, nchar, grabbed);
  sscanf(grabbed, "%d", &res);

  return res;
}

int T1Polyco_Read_NewFormat(T1Polyco *t1p, FILE *f)
{
  int ic;
  char junk[1024];

  if (fscanf(f, "%s %s %s %Lf %lf %lf %lf %Lf %Lf %s %d %d %lf %lf %lf",
	     t1p->psrname, t1p->date_string, t1p->utc_string,
	     &t1p->mjd_mid, &t1p->dm, &t1p->doppler, &t1p->log10rms,
	     &t1p->reference_phase, &t1p->frequency_psr_0,
	     t1p->sitename, &t1p->span, &t1p->ncoeff,
	     &t1p->frequency_obs, &t1p->binary_phase, &t1p->binary_frequency)
      != 15)
    return -1;

  for (ic=0; ic < t1p->ncoeff; ic++)
    if (fscanf(f, "%Lf", &t1p->coeff[ic])!=1)
      return -1;
  fgets(junk, 1024, f); /* eat final newline */
  return 0;
}

int T1Polyco_Read(T1Polyco *t1p, FILE *f)
{
  char line[128];
  int ic;

  if (fgets(line, 128, f)!=line)
    return -1;

  if (!strncasecmp(line, "TEMPO2:", 6))
    return T1Polyco_Read_NewFormat(t1p, f);

  T1P_grabString(line, 0, 11, t1p->psrname);
  T1P_grabString(line, 11, 9, t1p->date_string);
  T1P_grabString(line, 20, 11, t1p->utc_string);
  t1p->mjd_mid = T1P_grabLongDouble(line, 31, 20);
  t1p->dm = T1P_grabLongDouble(line, 51, 21);
  t1p->doppler = T1P_grabLongDouble(line, 73, 6);
  t1p->log10rms = T1P_grabLongDouble(line, 79, 7);
  
  if (fgets(line, 128, f)!=line)
    return -1;
  t1p->reference_phase = T1P_grabLongDouble(line, 0, 20);
  t1p->frequency_psr_0 = T1P_grabLongDouble(line, 21, 18);
  T1P_grabString(line, 38, 5, t1p->sitename);
  t1p->span = T1P_grabInt(line, 43, 5);
  t1p->ncoeff = T1P_grabInt(line, 48, 5);
  t1p->frequency_obs = T1P_grabLongDouble(line, 53, 10);
  t1p->binary_phase = T1P_grabLongDouble(line, 63, 7);
  t1p->binary_frequency = T1P_grabLongDouble(line, 70, 9);
  
  for (ic=0; ic < t1p->ncoeff; ic+=3)
  {
    if (fgets(line, 128, f)!=line)
      return -1;
    t1p->coeff[ic] = T1P_grabLongDouble(line, 0, 25);
    if (ic+1 < t1p->ncoeff)
      t1p->coeff[ic+1] = T1P_grabLongDouble(line, 25, 25);
    if (ic+2 < t1p->ncoeff)
      t1p->coeff[ic+2] = T1P_grabLongDouble(line, 50, 25);
  }

  return 0;
}
 
T1Polyco *T1PolycoSet_GetNearest(T1PolycoSet *t1ps, long double mjd)
{
  int inearest=0;
  long double best_offset=1e6, offset;
  int iseg;
  for (iseg=0; iseg < t1ps->nsegments ; iseg++)
  {
    offset = fabs(t1ps->segments[iseg].mjd_mid - mjd);
    if (offset < best_offset)
    {
      inearest = iseg;
      best_offset = offset;
    }
  }

  return &t1ps->segments[inearest];
}

long double T1PolycoSet_GetPhase(T1PolycoSet *t1ps, long double mjd, long double freq)
{
  return T1Polyco_GetPhase(T1PolycoSet_GetNearest(t1ps, mjd), mjd, freq);
}

long double T1PolycoSet_GetFrequency(T1PolycoSet *t1ps, long double mjd, long double freq)
{
  return T1Polyco_GetFrequency(T1PolycoSet_GetNearest(t1ps, mjd), mjd, freq);
}

void T1PolycoSet_Write(T1PolycoSet *t1ps, FILE *f)
{
  int iseg;
  for (iseg=0; iseg < t1ps->nsegments ; iseg++)
    T1Polyco_Write(&t1ps->segments[iseg], f);
}

int T1PolycoSet_Read(T1PolycoSet *t1ps, FILE *f)
{
  int nallocated = 5;
  int iseg=0;
  
  t1ps->segments = (T1Polyco *)malloc(nallocated*sizeof(T1Polyco));

  while (T1Polyco_Read(&t1ps->segments[iseg], f)==0)
  {
    iseg++;
    if (iseg >= nallocated)
    {
      nallocated += 5;
      t1ps->segments = (T1Polyco *)realloc(t1ps->segments,
					   nallocated*sizeof(T1Polyco));
    }
  }
  t1ps->nsegments = iseg;

  if (t1ps->nsegments > 0)
    return 0;
  else
    return -1;
}

void T1PolycoSet_Destroy(T1PolycoSet *t1ps)
{
  free(t1ps->segments);
}


