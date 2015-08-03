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

/*-*-C-*- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <float.h>
#include <glob.h>

#include "tempo2.h"
#include "dynarr.h"
#include "tabulatedfunction.h"

/* a sampled function connecting two clocks */
typedef struct
{
  TabulatedFunction table;
  char clockFrom[16];
  char clockTo[16];
  float badness;
  //  DynamicArray samples; 
} ClockCorrectionFunction;


/* Local static variables */

static int clockCorrections_initialized = 0;
DynamicArray clockCorrectionFunctions; /* All loaded functions */
DynamicArray clockCorrectionSequences; /* All auto-computed sequences.
					 NOTE: elements are Dynamic arrays of
					 POINTERS
				         into clockCorrectionFunctions
				         array */

/************************************************************************/

void
ClockCorrectionFunction_load(ClockCorrectionFunction *func,
			     char *fileName)
{
  int narg;

  TabulatedFunction_load(&func->table, fileName);

  narg = sscanf(func->table.header_line+1, // skip # 
		"%s %s %f", func->clockFrom, func->clockTo, &func->badness);
  if (narg < 2)
  {
    fprintf(stderr, 
	    "Error parsing clock file %s: first line must be of form # clock_from clock_to\n",
	    fileName);
    exit(1);
  }
  else if (narg != 3)
    func->badness = 1.0;
}

double
ClockCorrectionFunction_getCorrection(ClockCorrectionFunction *func,
				      double mjd)
{
  return TabulatedFunction_getValue(&func->table, mjd);
}

double
ClockCorrectionFunction_getStartMJD(ClockCorrectionFunction *func)
{
  return TabulatedFunction_getStartX(&func->table);
}
double
ClockCorrectionFunction_getEndMJD(ClockCorrectionFunction *func)
{
  return TabulatedFunction_getEndX(&func->table);
}


double
ClockCorrectionSequence_getStartMJD(DynamicArray *sequence)
{
  size_t ifunc;
  ClockCorrectionFunction *func;
  double mjd = 0.0, thismjd;

  for (ifunc=0; ifunc < sequence->nelem; ifunc++)
  {
    func = ((ClockCorrectionFunction **)sequence->data)[ifunc];
    thismjd = ClockCorrectionFunction_getStartMJD(func);
    if (thismjd > mjd)
      mjd = thismjd;
  }

  return mjd;
}
double
ClockCorrectionSequence_getEndMJD(DynamicArray *sequence)
{
  size_t ifunc;
  ClockCorrectionFunction *func;
  double mjd = 1e6, thismjd;

  for (ifunc=0; ifunc < sequence->nelem; ifunc++)
  {
    func = ((ClockCorrectionFunction **)sequence->data)[ifunc];
    thismjd =  ClockCorrectionFunction_getEndMJD(func);
    if (thismjd < mjd)
      mjd = thismjd;
  }

  return mjd;
}

void
initialize_ClockCorrections(int dispWarnings)
{
  glob_t g;
  char pattern[1024];
  char **pfname;
  int globRet;

  ClockCorrectionFunction func;
  DynamicArray_init(&clockCorrectionFunctions, sizeof(ClockCorrectionFunction));
  DynamicArray_init(&clockCorrectionSequences, sizeof(DynamicArray));

  /* load all .clk files in $TEMPO2/clock */
  sprintf(pattern, "%s/clock/*.clk", getenv(TEMPO2_ENVIRON));
  globRet = glob(pattern, 0, NULL, &g);
  if (globRet == GLOB_NOSPACE)
    { printf("Out of memory in clkcorr.C\n"); exit(1);}
#ifdef GLOB_ABORTED
  else if (globRet == GLOB_ABORTED)
    { printf("Read error in clkcorr.C\n"); exit(1); }
#endif
#ifdef GLOB_NOMATCH
  else if (globRet == GLOB_NOMATCH)
    { printf("No clock correction files in $TEMPO2/clock\n"); exit(1); }
#endif

  for (pfname = g.gl_pathv; *pfname != NULL; pfname++)
  {
    ClockCorrectionFunction_load(&func, *pfname);
    DynamicArray_push_back(&clockCorrectionFunctions, &func);
    if (dispWarnings==0)
      printf("Loaded clock correction %s : %s -> %s badness %.3g\n",
	     func.table.fileName, func.clockFrom, func.clockTo, func.badness);
  }
  clockCorrections_initialized = 1;
  globfree(&g);
}

void
defineClockCorrectionSequence(char *fileList_in,int dispWarnings)
{
  ClockCorrectionFunction *func;
  size_t ifunc;
  DynamicArray seq;
  char fileList[2048], *c, *white=" \t\r\n";
  
  if  ( clockCorrections_initialized != 1 )
    initialize_ClockCorrections(dispWarnings);

  DynamicArray_init(&seq, sizeof(ClockCorrectionFunction *));

  /* for each file name, find the function from the list of loaded ones,
     and add a pointer into that list to the
     array defining the present sequence */
  strcpy(fileList, fileList_in);
  for (c=strtok(fileList, white); c!=NULL; c=strtok(NULL, white))
  {
    for (ifunc=0; ifunc < clockCorrectionFunctions.nelem; ifunc++)
    {
      func = ((ClockCorrectionFunction *)clockCorrectionFunctions.data)+ifunc;
      if (!strcmp(func->table.fileName, c))
	  break;
    }
    if (ifunc == clockCorrectionFunctions.nelem)
    {
      printf( "Requested clock correction file %s not found!\n",c);
      exit(1);
    }
    DynamicArray_push_back(&seq, &func);
  }

  /* save the sequence in the list of sequences*/
  DynamicArray_push_back(&clockCorrectionSequences, &seq);
}

/* Function to make a clock correction sequence using Dijkstra's
   shortest path algorithm , see 
  http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm */
DynamicArray *makeClockCorrectionSequence(char *clockFrom, char *clockTo,
					  double mjd,int warnings)
{
  size_t slen = 16;
  DynamicArray names; /* Array of clock names. Other values index this arr */
  int *S; /* Clocks for which shortest path are known */
  int n_known;
  int *Q; /* Remaining clocks */
  float *d; /* Best distance to a clock */
  int *edge_to_previous; /* edge to use to link to previous node */ 
  int s, t; /* indices of clockFrom, clockTo */
  size_t ifunc;
  int iclock;
  ClockCorrectionFunction *func;
  char *clock;
  int to_present, from_present;
  int ibest;
  float best;
  /* make a list of edges */
  typedef struct
  {
    int from;
    int to;
    ClockCorrectionFunction *func;
  } Edge;
  Edge *edges;
  int iedge, v, nedges;
  int istep, nsteps;
  DynamicArray seq;

  //  printf("HERE %s %s %lf %d\n",clockFrom,clockTo,mjd,warnings);

  DynamicArray_init(&names, slen);

  /* initialize list of all known clocks, and list of edges
     for functions that cover the requested mjd */
  edges = (Edge *)malloc(sizeof(Edge)*clockCorrectionFunctions.nelem);
  iedge = 0;

  for (ifunc=0; ifunc < clockCorrectionFunctions.nelem; ifunc++)
  {
    func = ((ClockCorrectionFunction *)clockCorrectionFunctions.data)+ifunc;
    if (ClockCorrectionFunction_getStartMJD(func) <= mjd &&
	ClockCorrectionFunction_getEndMJD(func) >= mjd)
    {
      /* see if to and from clocks there already */
      to_present = from_present = 0;
      for (iclock=0;  
	   iclock < (int)names.nelem && !(to_present && from_present); iclock++)
      {
	clock = ((char *)names.data)+iclock*slen;
	if (!strcasecmp(clock, func->clockTo))
	{
	  to_present = 1;
	  edges[iedge].to = iclock;
	}
	if (!strcasecmp(clock, func->clockFrom))
	{
	  from_present = 1;
	  edges[iedge].from = iclock;
	  // 	  printf("%s present already\n", func->clockFrom); 
	}
/* 	else */
/* 	  printf("%s != %s!\n", clock, func->clockFrom); */
      }
      /* add to list if not already present */
      if (!to_present)
      {
	DynamicArray_push_back(&names, func->clockTo);
	edges[iedge].to = names.nelem-1;
      }
      if (!from_present)
      {
	DynamicArray_push_back(&names, func->clockFrom);
	edges[iedge].from = names.nelem-1;
      }
      //        printf("Add edge %d <%s> <%s> %d %d\n", iedge, func->clockTo, func->clockFrom,  
      //  	     edges[iedge].to, edges[iedge].from);  
      edges[iedge].func = func;
      iedge++;
    }
  }
  nedges = iedge;

  /* find requested start (s) and end (t) clock in array of names */
  s = t = -1;
  for (iclock=0; iclock < (int)names.nelem && (t < 0 || s < 0); iclock++)
  {
    clock = ((char *)names.data)+iclock*slen;
    if (!strcasecmp(clock, clockFrom))
      s = iclock;
    if (!strcasecmp(clock, clockTo))
      t = iclock;
  }
  if (s < 0)
  {
    char msg[1000],msg2[1000];
    sprintf(msg,"no clock corrections available for clock %s for MJD",
	    clockFrom);
    sprintf(msg2,"%.1f",mjd);
    displayMsg(1,"CLK3",msg,msg2,warnings);
	//logerr("no clock corrections available for clock %s -> %s for MJD %.1f",clockFrom,clockTo,mjd);
    //    printf("EXITING 1\n");
    free(edges); edges = NULL;
    DynamicArray_free(&names);
    return NULL;
  }
  if (t < 0)
  {
    char msg[1000],msg2[1000];
    sprintf(msg,"no clock corrections available for clock %s for MJD",
	    clockTo);
    sprintf(msg2,"%.1f",mjd);
    displayMsg(1,"CLK3",msg,msg2,warnings);
	//logerr("no clock corrections available for clock %s -> %s for MJD %.1f",clockFrom,clockTo,mjd);

    //    printf("EXITING 2\n");
    free(edges); edges = NULL;
  DynamicArray_free(&names);
    return NULL;
  }
 
  /* initialize S, Q, d, p */
  n_known = 0;
  S = (int *)malloc(sizeof(int)*names.nelem);
  Q = (int *)malloc(sizeof(int)*names.nelem);
  d = (float *)malloc(sizeof(float)*names.nelem);
  edge_to_previous = (int *)malloc(sizeof(int)*names.nelem);
  for (iclock=0; iclock < (int)names.nelem; iclock++)
  {
    S[iclock] = 0;
    Q[iclock] = 1;
    d[iclock] = FLT_MAX;
    edge_to_previous[iclock] = -1;
  }
  d[s] = 0.0;
  /* main algorithm */
  while (1)
  {
    /* find the member of Q with the lowest d */
    best = FLT_MAX;
    ibest = -1;
    for (iclock=0; iclock < (int)names.nelem; iclock++)
      if (Q[iclock] && d[iclock] < best)
      {
	best = d[iclock];
	ibest = iclock;
      }
    if (ibest == -1)
      break;
    Q[ibest] = 0; /* remove it from Q */
    S[ibest] = 1; /* add it to S */
    if (ibest == t)
      break; /* shortest path found */
    /* for every edge connected to ibest */
/*     printf("Looking for edges connected to %s d[%d] = %f\n",   */
/*  	   ((char *)names.data)+ibest*slen, ibest, best);  */
    for (iedge=0; iedge < nedges; iedge++)
      if (edges[iedge].from==ibest || edges[iedge].to==ibest)
      {
	/* set v to where the edge connects to */
	if (edges[iedge].from==ibest)
	  v = edges[iedge].to;
	else
	  v = edges[iedge].from;
	/* use this edge  if it is better */
	// 	printf("Try %s %s badness %f\n", 
	// 	       edges[iedge].func->clockFrom, edges[iedge].func->clockTo,  
	// 	       edges[iedge].func->badness); 
	if (d[v] > d[ibest] + edges[iedge].func->badness)
	{
	  d[v] = d[ibest] + edges[iedge].func->badness;
	  edge_to_previous[v] = iedge;
/*  	  printf("  edge %d --> %d %f\n", iedge, v, d[v]); */
	}
	/*	else */
	  /*  printf("Reject %s %s\n", edges[iedge].func->clockFrom, edges[iedge].func->clockTo); */
      }
  }
  
  //    printf("Badness total = %f\n", best);
  if (ibest == -1)
  {
    char msg[1000],msg2[1000];
    sprintf(msg,"no clock corrections available from %s to %s for MJD",
	    clockFrom,clockTo);
    sprintf(msg2,"%.1f",mjd);
	//logerr("no clock corrections available for clock %s -> %s for MJD %.1f",clockFrom,clockTo,mjd);
    displayMsg(1,"CLK7",msg,msg2,warnings);
    //    printf("EXITING 3\n");
    free(edges); edges = NULL;
  DynamicArray_free(&names);
    return NULL;
  }
  else
  {

    /* Write out the best path in reverse order, as an edge list to
       "S" */
    nsteps = 0;
    v = t;
    while (v != s)
    {
      S[nsteps] = edge_to_previous[v];
      v = (edges[edge_to_previous[v]].to==v ? 
	   edges[edge_to_previous[v]].from : edges[edge_to_previous[v]].to);
      nsteps++;
    }

    /* set up the actual sequence now in the correct order */
    DynamicArray_init(&seq, sizeof(ClockCorrectionFunction *));
    if (warnings==0)
      printf("Using the following chain of clock corrections for %s -> %s\n",
	   clockFrom, clockTo);
    for (istep=nsteps-1; istep >=0; istep--)
    {
      /* not any more since some functions are excluded (wrong mjd range) */
      /*     func = ((ClockCorrectionFunction *)clockCorrectionFunctions.data)  */
      /*       + S[istep]; */
      func = edges[S[istep]].func;
      DynamicArray_push_back(&seq, &func);
      if (warnings==0) printf("%s : %s <-> %s\n", func->table.fileName, func->clockFrom, func->clockTo);
    }
  }
  /* free memory */
  free(edges); edges = NULL;
  free(edge_to_previous); edge_to_previous = NULL;
  free(d); d = NULL;
  free(Q); Q = NULL;
  free(S); S = NULL;
  DynamicArray_free(&names);
  //  printf("END HERE\n");
  return (DynamicArray *)DynamicArray_push_back(&clockCorrectionSequences, &seq);
}

DynamicArray *getClockCorrectionSequence(char *clockFrom, char *clockTo,
					 double mjd,int warnings)
{
  size_t iseq;
  DynamicArray *seq;
  ClockCorrectionFunction *firstFunc, *lastFunc;

  //  return makeClockCorrectionSequence(clockFrom, clockTo, mjd,warnings);

  if  ( clockCorrections_initialized != 1 )
    initialize_ClockCorrections(warnings);

  /* search for the given clocks in the list of already defined sequences */
  for (iseq=0; iseq < clockCorrectionSequences.nelem; iseq++)
  {
    seq = ((DynamicArray *)clockCorrectionSequences.data) + iseq;
    firstFunc = ((ClockCorrectionFunction **)seq->data)[0];
    lastFunc  = ((ClockCorrectionFunction **)seq->data)[seq->nelem - 1];

    if ( 
	// from clock matches first function (backwards or forwards)
	(!strcasecmp(firstFunc->clockFrom, clockFrom)
	 || !strcasecmp(firstFunc->clockTo, clockFrom))
	&&
	// to clock matches last function (backwards or forwards)
	(!strcasecmp(lastFunc->clockFrom, clockTo)
	 || !strcasecmp(lastFunc->clockTo, clockTo))
	&&
	// date is within possible range
	ClockCorrectionSequence_getStartMJD(seq) <= mjd
	&& ClockCorrectionSequence_getEndMJD(seq) >= mjd
	)
      return seq;
  }

  /* no pre-defined sequence found. Search for one using Dijkstra's
     shortest-path algorithm. */
  return makeClockCorrectionSequence(clockFrom, clockTo, mjd,warnings);

}
   

#if 0 /* no longer used */
double 
getClockCorrection(double mjd, char *clockFrom, char *clockTo)
{
  DynamicArray *sequence = getClockCorrectionSequence(clockFrom, clockTo);
  double correction = 0;
  size_t ifunc;
  char *currClock = clockFrom;
  ClockCorrectionFunction *func;

  for (ifunc=0; ifunc < sequence->nelem; ifunc++)
  {
    func = ((ClockCorrectionFunction **)sequence->data)[ifunc];

    /* add correction using sign based on direction of correction */
    correction += ClockCorrectionFunction_getCorrection(func, mjd)
      * (!strcasecmp(currClock, func->clockFrom) ? 1.0 : -1.0);

    currClock = func->clockTo;
  }

  return correction;
}
#endif

void 
getClockCorrections(observation *obs, char *clockFrom, 
		    char *clockTo,int warnings)
{
  DynamicArray *sequence;
  size_t ifunc;
  char *currClock;
  ClockCorrectionFunction *func;
  double correction = 0.0;
  obs->nclock_correction = 0;

  char *clockFromOrig; 
  //  logdbg("In getClockCorrections");

  //logdbg("Getting clockFrom >%s<",obs->telID);
  if (clockFrom[0]=='\0') // get from observatory instead
    clockFrom = getObservatory(obs->telID)->clock_name;
  //  logdbg("Got clockFrom");
  clockFromOrig = clockFrom;

  if (!strcasecmp(clockTo, clockFrom))
  {
    obs->nclock_correction = 0;
    return;
  }
  //  logdbg("In getClockCorrections calling sequence 1");
  // Memory leak here
  sequence = getClockCorrectionSequence(clockFrom, clockTo, obs->sat,warnings);

  if (sequence == NULL)
  {
    char msg[1000],msg2[1000];
    sprintf(msg,"Trying assuming UTC =");
    sprintf(msg2,"%s",clockFrom);
    displayMsg(1,"CLK4",msg,msg2,warnings);
    clockFrom="UTC";
    if (!strcasecmp(clockTo, clockFrom))
    {
      obs->nclock_correction = 0;
      return;
    }
    logdbg("In getClockCorrections calling sequence 2");
    sequence = getClockCorrectionSequence(clockFrom, clockTo, obs->sat,
					  warnings); 
    if (sequence == NULL)
    {
      char msg[1000],msg2[1000];
      sprintf(msg,"Trying TT(TAI) instead of ");
      sprintf(msg2,"%s",clockTo);
      displayMsg(1,"CLK5",msg,msg2,warnings);
      clockFrom = clockFromOrig;
      clockTo="TT(TAI)";
      if (!strcasecmp(clockTo, clockFrom))
      {
	obs->nclock_correction = 0;
	return;
      }
      logdbg("In getClockCorrections calling sequence 3");
      sequence = getClockCorrectionSequence(clockFrom, clockTo, obs->sat,
					    warnings); 
     if (sequence == NULL)
     {
       displayMsg(1,"CLK8","Trying both","",warnings);       
        clockFrom="UTC";
	if (!strcasecmp(clockTo, clockFrom))
	{
	  obs->nclock_correction = 0;
	  return;
	}
	logdbg("In getClockCorrections calling sequence 4");
	sequence = getClockCorrectionSequence(clockFrom, clockTo, obs->sat,
					      warnings); 

      if (sequence==NULL)
      {
       obs->nclock_correction = 0;
       if (warnings==0)
	 printf( "Warning [CLK:7], no clock correction available for TOA @ MJD %.4lf!\n",
		 (double)obs->sat);
       return;
     }
     }
    }
    displayMsg(1,"CLK9","... ok, using stated approximation","",warnings);       
  }
  currClock = clockFrom;
  // Already bad
  for (ifunc=0; ifunc < sequence->nelem 
	 && strcasecmp(currClock, clockTo); ifunc++)
  {
    func = ((ClockCorrectionFunction **)sequence->data)[ifunc];
    bool backwards = strcasecmp(currClock, func->clockFrom);
    if (backwards && strcasecmp(currClock, func->clockTo))
    {
      printf("Programming error! Broken clock correction chain!!\n");
      exit(1);
    }
    /* add correction using sign based on direction of correction */
    obs->correctionsTT[obs->nclock_correction].correction = 
      ClockCorrectionFunction_getCorrection(func, obs->sat+correction/SECDAY)
      * (backwards ? -1.0 : 1.0);
    correction += obs->correctionsTT[obs->nclock_correction].correction;
    strcpy(obs->correctionsTT[obs->nclock_correction].corrects_to, 
	   (backwards ? func->clockFrom : func->clockTo)); 
    obs->nclock_correction++;
    currClock = (backwards ? func->clockFrom : func->clockTo);
  }
  //  logdbg("leaving getClockCorrections");
/*   printf("Correction: %lg\n", correction); */
}

double 
getCorrectionTT(observation *obs)
{
  double correction = 0.0;
  int ic;
  for (ic=0; ic < obs->nclock_correction; ic++)
    correction += obs->correctionsTT[ic].correction;
  return correction;
}

double
getCorrection(observation *obs, char *clockFrom, char *clockTo, int warnings)
{
    observatory *site;
  DynamicArray *sequence;
  size_t ifunc;
  char *currClock;
  ClockCorrectionFunction *func;
  double correction = 0.0;
  const char *CVS_verNum = "$Revision: 1.10 $";

  if (clockFrom[0]=='\0')
    site = getObservatory(obs->telID);

  if (displayCVSversion == 1) CVSdisplayVersion("clkcorr.C","getCorrection()",CVS_verNum);

  if (clockFrom[0]=='\0')
    clockFrom = site->clock_name;
  currClock = clockFrom;

  if (!strcasecmp(clockTo, clockFrom))
    return 0.0;

/*   printf("Getting %s<->%s\n", site->clock_name, clockTo); */
  sequence = getClockCorrectionSequence(clockFrom, clockTo, obs->sat,
					warnings);

  if (sequence == NULL)
  {
    char msg[1000],msg2[1000];
    sprintf(msg,"Proceeding assuming UTC = ");
    sprintf(msg2,"%s",currClock);
    displayMsg(1,"CLK6",msg,msg2,warnings);

    if (!strcasecmp(clockTo, "UTC"))
      return 0.0;
    sequence = getClockCorrectionSequence("UTC", clockTo, obs->sat,
					  warnings); 
    currClock = "UTC";
    if (sequence == NULL)
    {
      if (warnings==0) 
	printf( "!Warning [CLK:9], no clock correction available for TOA @ MJD %.4lf!\n",
		(double)obs->sat);
      return 0.0;
    }
  }

  for (ifunc=0; ifunc < sequence->nelem
	 && strcasecmp(currClock, clockTo); ifunc++)
  {
    func = ((ClockCorrectionFunction **)sequence->data)[ifunc];
    bool backwards = strcasecmp(currClock, func->clockFrom);
    /* add correction using sign based on direction of correction */
    correction += 
      ClockCorrectionFunction_getCorrection(func, obs->sat+correction/SECDAY)
      * (backwards ? -1.0 : 1.0);
    currClock = (backwards ? func->clockFrom : func->clockTo);
/*     printf("--> %s\n", currClock); */
  }


  return correction;
}
