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

/* This file contains useful functions that are commonly used throughout the TEMPO2 */
/* software */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "tempo2.h"

/* Computes the dot product of two vectors */
double dotproduct(double *v1,double *v2)
{
  double dot;
  dot = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
  return dot;
}
/* Computes the sum of two vectors */
void vectorsum(double *res, double *v1,double *v2)
{
  res[0] = v1[0] + v2[0];
  res[1] = v1[1] + v2[1];
  res[2] = v1[2] + v2[2];
}

void vectorscale(double *v, double k)
{
  v[0] *= k;
  v[1] *= k;
  v[2] *= k;
}

int turn_hms(double turn, char *hms){
 
  /* Converts double turn to string " hh:mm:ss.ssss" */
  
  int hh, mm, isec;
  double sec;

  hh = (int)(turn*24.);
  mm = (int)((turn*24.-hh)*60.);
  sec = ((turn*24.-hh)*60.-mm)*60.;
  isec = (int)((sec*10000. +0.5)/10000);
    if(isec==60){
      sec=0.;
      mm=mm+1;
      if(mm==60){
        mm=0;
        hh=hh+1;
        if(hh==24){
          hh=0;
        }
      }
    }

  sprintf(hms," %02d:%02d:%010.7f",hh,mm,sec);
  return 0;
}

int turn_dms(double turn, char *dms){
  
  /* Converts double turn to string "sddd:mm:ss.sss" */
  
  int dd, mm, isec;
  double trn, sec;
  char sign;
  
  sign=' ';
  if (turn < 0.){
    sign = '-';
    trn = -turn;
  }
  else{
    sign = '+';
    trn = turn;
  }
  dd = (int)(trn*360.);
  mm = (int)((trn*360.-dd)*60.);
  sec = ((trn*360.-dd)*60.-mm)*60.;
  isec = (int)((sec*1000. +0.5)/1000);
    if(isec==60){
      sec=0.;
      mm=mm+1;
      if(mm==60){
        mm=0;
        dd=dd+1;
      }
    }
  sprintf(dms,"%c%02d:%02d:%08.5f",sign,dd,mm,sec);
  return 0;
}

double turn_deg(double turn){
 
  /* Converts double turn to string "sddd.ddd" */
  return turn*360.0;
}

double hms_turn(char *line) {

  /* Converts string " hh:mm:ss.ss" or " hh mm ss.ss" to double turn */
  
  int i;
  double hr, min, sec, turn=0.0;

  /* Get rid of ":" */
  for(i=0; *(line+i) != '\0'; i++)if(*(line+i) == ':')*(line+i) = ' ';

  i = sscanf(line,"%lf %lf %lf", &hr, &min, &sec);

  if(i > 0){
    turn = hr/24.;
    if(i > 1)turn += min/1440.;
    if(i > 2)turn += sec/86400.;
  }
  if(i == 0 || i > 3)turn = 1.0;

  return turn;
}

double dms_turn(char *line){

  /* Converts string "-dd:mm:ss.ss" or " -dd mm ss.ss" to double turn */
  
  int i;
  char *ic, ln[40];
  double deg, min, sec, sign, turn=0.0;

  /* Copy line to internal string */
  strcpy(ln,line);

  /* Get rid of ":" */
  for(i=0; *(ln+i) != '\0'; i++)if(*(ln+i) == ':')*(ln+i) = ' ';

  /* Get sign */
  if((ic = strchr(ln,'-')) == NULL)
     sign = 1.;
  else {
     *ic = ' ';
     sign = -1.;
  }

  /* Get value */
  i = sscanf(ln,"%lf %lf %lf", &deg, &min, &sec);
  if(i > 0){
    turn = deg/360.;
    if(i > 1)turn += min/21600.;
    if(i > 2)turn += sec/1296000.;
    if(turn >= 1.0)turn = turn - 1.0;
    turn *= sign;
  }
  if(i == 0 || i > 3)turn =1.0;

  return turn;
}

/* Emulate the Fortran "mod" operator */
longdouble fortran_mod(longdouble a,longdouble p)
{
  longdouble ret;
  ret = a - (int)(a/p)*p;
  return ret;
}
double fortran_mod(double a,double p)
{
  longdouble ret;
  ret = a - (int)(a/p)*p;
  return ret;
}

int fortran_nint(double x)
{
  int i;
  if(x>0.){
    i=(int)(x+0.5);
  }
  else{
    i=(int)(x-0.5);
  }
  return(i);
}

long fortran_nlong(longdouble x)
{
  long i;
  if(x>0.){
    i=(long)(x+0.5);
  }
  else{
    i=(long)(x-0.5);
  }
  return(i);
}

/* print out a long double to a std::string */
std::string print_longdouble(const longdouble &ld)
{
  char buf[1024];
#ifdef USE_BUILTIN_LONGDOUBLE
  sprintf(buf, "%Lg", ld);
#else
  ld.write(buf);
#endif
  return std::string(buf);

}

longdouble parse_longdouble(const char *str)
{
  longdouble ld;
#ifdef USE_BUILTIN_LONGDOUBLE
  sscanf(str, "%Lf", &ld);
#else
  ld = str;
#endif
  return ld;
}

/* Rotates the vector x(0,1,2) from equatorial coordinates to ecliptic coordinates 
 * i.e. rotate about "x" axis by angle epsilon 
 *
 * original tempo: epsilon = 84381.412 arc sec at epoch of J2000.0 (IERS tech note 21, p. 19)
 * tempo2 using ECLIPTIC_OBLIQUITY
 * 
 * ce = cos(epsilon)
 * se = sin(epsilon)
 *
 */

void equ2ecl(double *x)
{
  longdouble tmpy,tmpz;
  longdouble ce,se;
  longdouble arcsec2rad = 1.0/60.0/60.0*M_PI/180.0;

  /*  ce = 0.91748213149438;
      se = 0.39777699580108;  */
  ce = cosl(ECLIPTIC_OBLIQUITY*arcsec2rad);
  se = sinl(ECLIPTIC_OBLIQUITY*arcsec2rad); 

  /* x[0] remains unchanged */

  tmpy = x[1];
  tmpz = x[2];

  x[1] = ce*tmpy+se*tmpz;
  x[2] = -se*tmpy+ce*tmpz;
}

/* Long double support routines */
#ifndef USE_BUILTIN_LONGDOUBLE
dd_real pow(const dd_real &a, const dd_real &b)
{ return exp(b*log(a)); }
// operator float(const dd_real &a) 
// {return (float)(double)a;}
#endif

/* Function to copy one pulsar to another */
void copyPSR(pulsar *p,int p1,int p2)
{
  int i;

  strcpy(p[p2].name,p[p1].name);
  strcpy(p[p2].binaryModel,p[p1].binaryModel);
  p[p2].fitMode = p[p1].fitMode;
  p[p2].nobs = p[p1].nobs;
  p[p2].nits = p[p1].nits;
  for (i=0;i<p[p1].nobs;i++)
    p[p2].obsn[i] = p[p1].obsn[i]; 
  for (i=0;i<3;i++)
    {
      p[p2].posPulsar[i]=p[p1].posPulsar[i];
      p[p2].velPulsar[i]=p[p1].velPulsar[i];     
    }
  p[p2].eclCoord = p[p1].eclCoord;
  p[p2].nJumps = p[p1].nJumps;
  for (i=0; i<=p[p2].nJumps; i++) 
    {
      p[p2].jumpVal[i] = p[p1].jumpVal[i];
      p[p2].fitJump[i] = p[p1].fitJump[i];
      p[p2].jumpValErr[i] = p[p1].jumpValErr[i];
      strcpy(p[p2].jumpStr[i], p[p1].jumpStr[i]);
    }
  p[p2].nFit = p[p1].nFit;
  p[p2].fitMode = p[p1].fitMode;
  p[p2].offset = p[p1].offset;
  p[p2].calcShapiro = p[p1].calcShapiro;
  p[p2].units = p[p1].units;
  p[p2].tempo1 = p[p1].tempo1;
  p[p2].dilateFreq = p[p1].dilateFreq;
  p[p2].timeEphemeris = p[p1].timeEphemeris;
  p[p2].t2cMethod = p[p1].t2cMethod;
  p[p2].correctTroposphere = p[p1].correctTroposphere;
  strcpy(p[p2].clock,p[p1].clock);
  strcpy(p[p2].clockFromOverride,p[p1].clockFromOverride);
  strcpy(p[p2].JPL_EPHEMERIS,p[p1].JPL_EPHEMERIS);
  strcpy(p[p2].ephemeris,p[p1].ephemeris);

  // Note, this causes memory leaks.  Tempo2 always allocates memory for
  // MAX_PSR pulsars to start.
  //allocateMemory(p+p2,1);
}


/* Function to copy one parameter into another */
void copyParam(parameter p1,parameter *p2)
{
  int i;
  for (i=0;i<p1.aSize;i++)
    {
      strcpy(p2->label[i],p1.label[i]);
      strcpy(p2->shortlabel[i],p1.shortlabel[i]);
      p2->val[i]       = p1.val[i];
      p2->err[i]       = p1.err[i];
      p2->fitFlag[i]   = p1.fitFlag[i];
      p2->paramSet[i]  = p1.paramSet[i];
      p2->prefit[i]    = p1.prefit[i];
      p2->prefitErr[i] = p1.prefitErr[i];      
    }
  p2->aSize = p1.aSize;
}

/* Function to display warning and error messages */
/* type = 1 for warning, = 2 for error            */
/* key = e.g. CLK1                                */

void displayMsg(int type,char *key,char *searchStr,char *variableStr,int noWarnings)
{
  static char msg[MAX_MSG][1000];
  static char keyRec[MAX_MSG][10];
  static int warningNum=0;
  static int warnedAboutMultiple=0;
  int i,got;

  if (type==1)
    {
  if (noWarnings!=2)
    {

      if (noWarnings>=0)
	{
	  got=0;
	  for (i=0;i<warningNum;i++)
	    {
	      if (strcasecmp(keyRec[i],key)==0 &&
		  strcasecmp(msg[i],searchStr)==0)
		{
		  got=1;
		  break;
		}
	    }

	  if (got==0)
	    {
	      if (warningNum+1 > MAX_MSG)
		{
		  printf("Too many warnings: bailing out\n");
		  exit(1);
		}
	      strcpy(keyRec[warningNum],key);
	      strcpy(msg[warningNum],searchStr);
	      warningNum++;
	      
	      printf("WARNING [%s]: %s %s\n",key,searchStr,variableStr);
	    }
	  else if (!warnedAboutMultiple)
	  {
	    printf("WARNING: duplicated warnings have been suppressed.\n");
	    warnedAboutMultiple=1;
	  }
	}
      else
	printf("WARNING [%s]: %s %s\n",key,searchStr,variableStr);
    }
    }
  else
    {
      printf("\n\n");
      printf("-------------------------------------------------\n");
      printf("ERROR [%s]: %s %s\n",key,searchStr,variableStr);
      printf("-------------------------------------------------\n");
      printf("\n\n");
      exit(1);
    }
}

/* Gets a parameter value and checks whether this parameter is linked to */
/* another parameter                                                     */

longdouble getParameterValue(pulsar *psr,int param,int arr)
{
  if (psr->param[param].nLinkTo > 0)
    {
      //      printf("Here checking param %d %s\n",param,psr[0].param[param].shortlabel[0]);
      if (psr->param[param].nLinkTo > 1)
	{
	  printf("Sorry, can currently only link two parameters together\n");
	  exit(1);
	}
      if (param==param_daop){
	if (psr->param[param].linkTo[0] == param_pbdot)
	  return (psr->param[param_pbdot].val[0]*SPEED_LIGHT
		  /(psr->param[param_pb].val[0]*SECDAY*
		    (pow(psr->param[param_pmra].val[0]*MASYR2RADS,2)
		     +pow(psr->param[param_pmdec].val[0]*MASYR2RADS,2))))/PCM;
      }
      if (param==param_kin)
	{
	  if (psr->param[param].linkTo[0] == param_stig)
	    {
	      long double stig = psr->param[param_stig].val[0];
	      //	      printf("kin = %g\n",(double)((2.0*stig)/(1+stig*stig)));
	      return asin((2.0*stig)/(1+stig*stig))*180.0/M_PI;
	    }
	}
      if (param==param_sini)
	{
	  if (psr->param[param].linkTo[0] == param_kin)
	    return sin(psr->param[param_kin].val[0]/180*M_PI);
	}
      if (param==param_dshk)
	{
	  if (psr->param[param].linkTo[0] == param_px)
	    return 1.0/psr->param[param_px].val[0];
	}
      if (param==param_dmmodel){
	      if (psr->param[param].linkTo[0] == param_dm)
		      return psr->param[param_dm].val[0];

      }
    }
  return psr->param[param].val[arr];
}
