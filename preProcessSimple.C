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

#include <stdio.h>
#include <stdlib.h>
#include "tempo2.h"
#include <math.h>
#include <string.h>

void preProcessSimple (pulsar *psr)
{
  const char *CVS_verNum = "$Revision: 1.9 $";

  if (displayCVSversion == 1) CVSdisplayVersion("preProcessSimple.C","preProcessSimple()",CVS_verNum);

  preProcessSimple1 (psr, 0, -1);
  preProcessSimple2 (psr, 0, 0, 0, 0);
  preProcessSimple3 (psr);
}

void preProcessSimple1 (pulsar *psr, int tempo1, double thelast)
{  
  logdbg("In preProcessSimple1");

  ECLIPTIC_OBLIQUITY = ECLIPTIC_OBLIQUITY_VAL;

  int i = 0;

  /* Check whitening */
  if (psr->param[param_wave_om].paramSet[0] == 1 && psr->param[param_wave_om].val[0] == 0.0) /* Set fundamental frequency */
    {  
      longdouble first=-1,last=-1;
      for (i=0;i<psr->nobs;i++)
	{
	  if (psr->obsn[i].deleted==0 && first==-1) first=psr->obsn[i].sat;
	  if (psr->obsn[i].deleted==0 && last==-1)  last=psr->obsn[i].sat;
	  if (psr->obsn[i].deleted==0 && psr->obsn[i].sat < first) first = psr->obsn[i].sat;
	  if (psr->obsn[i].deleted==0 && psr->obsn[i].sat > last)   last = psr->obsn[i].sat;
	}
      psr->param[param_wave_om].val[0] = 2.0*M_PI/(last-first)/
	(1.0+4.0/(double)(psr->nWhite));
      fprintf(stderr, "%.3e\n",   psr->param[param_wave_om].val[0]);
     
    }


   /* Check whitening */
  if (psr->param[param_wave_dm].paramSet[0] == 1 && psr->param[param_wave_dm].val[0] == 0.0) /* Set fundamental frequency */
    {  
      longdouble first=-1,last=-1;
      for (i=0;i<psr->nobs;i++)
	{
	  if (psr->obsn[i].deleted==0 && first==-1) first=psr->obsn[i].sat;
	  if (psr->obsn[i].deleted==0 && last==-1)  last=psr->obsn[i].sat;
	  if (psr->obsn[i].deleted==0 && psr->obsn[i].sat < first) first = psr->obsn[i].sat;
	  if (psr->obsn[i].deleted==0 && psr->obsn[i].sat > last)   last = psr->obsn[i].sat;
	}
      psr->param[param_wave_dm].val[0] = 2.0*M_PI/(last-first)/
	(1.0+4.0/(double)(psr->nWhite));
    }



  /* Set tempo emulation mode */
  if (psr->param[param_ephver].paramSet[0]==1)
    {
      if (psr->param[param_ephver].val[0] < 5)
	{
	  printf("************************************************* \n");
	  printf("Warning: you are running in tempo1 emulation mode \n");
	  printf("************************************************* \n");
	  tempo1 = 1;
	}
      else
	tempo1 = 0;
    }
  else if (tempo1 == 1)
    {
      psr->param[param_ephver].paramSet[0]=1;
      psr->param[param_ephver].val[0]=2;
    }
  else 
    {
      psr->param[param_ephver].paramSet[0]=1;
      psr->param[param_ephver].val[0]=5;
    }
  if (thelast>0)
    /* Define the start parameter based upon the last observation */
    {
      psr->param[param_start].val[0] = psr->obsn[psr->nobs-1].sat-thelast;
      psr->param[param_start].fitFlag[0] = 1;
      psr->param[param_start].paramSet[0] = 1;
    }

  psr->tempo1 = 0;
  if (tempo1 == 1)
    {
      psr->tempo1 = 1;
      psr->units = TDB_UNITS;
      psr->timeEphemeris = FB90_TIMEEPH;
      psr->dilateFreq = 0;
      psr->planetShapiro=0;
      psr->t2cMethod = T2C_TEMPO;
      psr->correctTroposphere = 0;
      psr->ne_sw = 9.961;
      ECLIPTIC_OBLIQUITY = 84381.412;
    }
  /* XXXX Hack!! -- removed - should use transform plugin*/
  /* Problem, if a function is not called then the compiler does not include
   * it in the library - so the plugin cannot find it */
  if (psr->units==100) /* SI_UNITS)  */
    transform_units(psr, TDB_UNITS, SI_UNITS);
}

void preProcessSimple2 (pulsar *psr,
			float startdmmjd, int ndm, float* dmvals, int trimonly)
{  
  logdbg("In preProcessSimple2");

  int i=0, j=0, k=0;

  for (i=0;i<psr->nobs;i++)
    {
      //psr->obsn[i].efac = 1.0;
      //psr->obsn[i].equad = 0.0;
      /* Check whether any error = 0 */
      if (psr->fitMode==1 && psr->obsn[i].toaErr==0)
	{
	  printf("Error: the uncertainty for TOA %d (MJD %f) is equal to zero\n ",i+1,(double)psr->obsn[i].sat);
	  exit(1);
	}

      for (k=0; k<6; k++) /* zero these vectors for first call to tt2tb */
	{
	  psr->obsn[i].earthMoonBary_ssb[k] = 0.0;
	  psr->obsn[i].earthMoonBary_earth[k] = 0.0;
	  psr->obsn[i].observatory_earth[k] = 0.0;
	  psr->obsn[i].earth_ssb[k] = 0.0;
	}

      // Correctly taking note of the phase offset
      psr->obsn[i].sat += ((psr->obsn[i].phaseOffset/psr->param[param_f].val[0])/SECDAY); 
      for (k=0;k<psr->nToffset;k++) /* Calculate time offsets */
	{
	  char offsetSite[256], obsSite[256];
	  lookup_observatory_alias(psr->tOffsetSite[k], offsetSite);
	  lookup_observatory_alias(psr->obsn[i].telID, obsSite);
	  if ((psr->tOffset_f1[k]==0.0 || (psr->obsn[i].freq > psr->tOffset_f1[k])) &&
	      (psr->tOffset_f2[k]==0.0 || (psr->obsn[i].freq < psr->tOffset_f2[k])) &&
	      (psr->tOffset_t1[k]==0.0 || (psr->obsn[i].sat  > psr->tOffset_t1[k])) &&
	      (psr->tOffset_t2[k]==0.0 || (psr->obsn[i].sat  < psr->tOffset_t2[k])) && 
	      (strcmp(offsetSite,"0")==0 || 
	       strcmp(offsetSite,obsSite)==0))
	    {
	      int use=1;
	      /* Check for flags */
	      if (strlen(psr->tOffsetFlags[k])>0)
		{
		  char *myStr,str1[1000],flagID[100];
		  use=0;
		  strcpy(str1,psr->tOffsetFlags[k]);
		  myStr = strtok(str1," ");
		  do {
		    if (myStr[0]=='-' && (myStr[1]<48 || myStr[1]>57)) /* Look for flag */
		      {
			strcpy(flagID,myStr);
			myStr = strtok(NULL," ");
			for (int l=0;l<psr->obsn[i].nFlags;l++)
			  {
			    if (strcmp(flagID,psr->obsn[i].flagID[l])==0
				&& strcmp(myStr,psr->obsn[i].flagVal[l])==0)
			      {
				use=1;
				break;
			      }
			  }
		      }
		  } while ((myStr = strtok(NULL," "))!=NULL);
		}
	      if (use==1) {psr->obsn[i].sat += psr->tOffset[k]/SECDAY;}
	    }  
	}
      /* Check filtering */
      if (strlen(psr->filterStr)>0)
	{
	  char flag[100],filtS[1000],*filt;
	  strcpy(filtS,psr->filterStr);
	  filt = strtok(filtS," ");
	  strcpy(flag,filt);
	  while (filt != NULL)
	    {
	      if (filt[0]=='-') strcpy(flag,filt);
	      else
		{
		  for (j=0;j<psr->obsn[i].nFlags;j++)
		    {
		      if (strcmp(psr->obsn[i].flagID[j],flag)==0 &&
			  strcmp(psr->obsn[i].flagVal[j],filt)==0)
			psr->obsn[i].deleted=1;
		    }
		}
	      filt = strtok(NULL," ");
	    }
	}

      /* Check filtering for 'pass' */
      if (strlen(psr->passStr)>0)
	{
	  char flag[100],filtS[1000],*filt;
	  int found;
	  strcpy(filtS,psr->passStr);
	  filt = strtok(filtS," ");
	  found=0;
	  while ( filt != NULL && found==0)
	    {
	      if (filt[0]=='-') // Have a flag
		strcpy(flag,filt);
	      else
		{
		  for (j=0;j<psr->obsn[i].nFlags;j++)
		    {
		      if (strcmp(psr->obsn[i].flagID[j],flag)==0 &&
			  strcmp(psr->obsn[i].flagVal[j],filt)==0)
			{
			  found=1;
			  break;
			}
		    }
		}
	      filt = strtok(NULL," ");
	    }

	  if (found==0)
	    psr->obsn[i].deleted=1;		 
	}

      /* CHECK JUMPS */	  
      for (k=1;k<=psr->nJumps;k++)
	{
	  int yes=0;

	  char str1[100],str2[100],str3[100],str4[100],str5[100];
	  double val1, val2;

	  /* Must parse jump string to determine whether this observation should jump or not */	      
	  int nread = sscanf(psr->jumpStr[k],"%s %s %s",str1,str2,str3);
	  if (strcasecmp(str1,"MJD")==0 || strcasecmp(str1,"FREQ")==0)
	    {
	      sscanf(str2,"%lf",&val1);
	      sscanf(str3,"%lf",&val2);

	      if (strcasecmp(str1,"MJD")==0) {
		if (psr->obsn[i].sat >= val1 && psr->obsn[i].sat < val2) {
		  yes=1;
		}
	      }
	      else if (strcasecmp(str1,"FREQ")==0) {
		if (psr->obsn[i].freq >= val1 && psr->obsn[i].freq < val2) {
		  yes=1;
		}
	      }
	    }
	  if (strcasecmp(str1,"NAME")==0 && strstr(psr->obsn[i].fname,str2)!=NULL)
	    yes=1;	
	  if (strcasecmp(str1,"TEL")==0)
	    {
	      char selectedSite[256], obsSite[256];

	      lookup_observatory_alias(str2, selectedSite);
	      lookup_observatory_alias(psr->obsn[i].telID, obsSite);
	      if (strcasecmp(selectedSite, obsSite)!=0)
		yes=1;	    
	    }
	  else if (str1[0]=='-')
	    {
	      int jj;
	      for (jj=0;jj<psr->obsn[i].nFlags;jj++)
		{
		  if (strcmp(psr->obsn[i].flagID[jj],str1)==0 &&
		      strcmp(psr->obsn[i].flagVal[jj],str2)==0)
		    yes=1;
		}
	    }
	    
	  if (yes==1) {psr->obsn[i].jump[psr->obsn[i].obsNjump]=k; (psr->obsn[i].obsNjump)++;}
	}
      // Check for time offset flags
      for (k=0;k<psr->obsn[i].nFlags;k++)
	{
	  if (strcmp(psr->obsn[i].flagID[k],"-to")==0)
	    {
	      long double v;
	      sscanf(psr->obsn[i].flagVal[k],"%Lf",&v);
	      psr->obsn[i].sat += v/SECDAY;
	    }
	}
      // Check for dm updates
      if (startdmmjd && ndm && dmvals)
	{
	  //float dm; //not needed any more
	  if ((int)(psr->obsn[i].sat-startdmmjd+0.5) > 0)
	    {
	      if ((int)(psr->obsn[i].sat-startdmmjd+0.5) >= ndm)
		psr->obsn[i].deleted=1;
	      else if (trimonly == 0) // set flags only if really correcting for dm, not only trimming the data
		{
		  strcpy(psr->obsn[i].flagID[psr->obsn[i].nFlags],"-dm");
		  sprintf(psr->obsn[i].flagVal[psr->obsn[i].nFlags],"%g",dmvals[(int)(psr->obsn[i].sat-startdmmjd+0.5)]);
		  psr->obsn[i].nFlags++;
		  if (psr->obsn[i].nFlags >= MAX_FLAGS)
		    {
		      printf("Number of different flags in the .tim file > MAX_FLAGS (%d)\n",MAX_FLAGS);
		      exit(1);
		    }
		  
		  //dm = dmvals[(int)(psr->obsn[i].sat-startdmmjd+0.5)];             
		}
	    }
	  else
	    {
	      psr->obsn[i].deleted=1;
	    }
	}
    }
  logdbg("Complete preProcessSimple2");
}

void preProcessSimple3 (pulsar *psr)
{  
  int  i,j,k,l;
  logdbg("Start preProcessSimple3");      
  //
  // Check fjump
  //
  if (strlen(psr->fjumpID)>0)
    {
      char val[MAX_JUMPS][16];
      int  nval[MAX_JUMPS];
      int  nf=0;
      int found;
      // Find flag corresponding to longest data set
      for (i=0;i<psr->nobs;i++)
	{
	  if (psr->obsn[i].deleted==0)
	    {
	      for (k=0;k<psr->obsn[i].nFlags;k++)
		{
		  if (strcmp(psr->obsn[i].flagID[k],psr->fjumpID)==0)
		    {
		      found=0;
		      for (l=0;l<nf;l++)
			{
			  if (strcmp(psr->obsn[i].flagVal[k],val[l])==0)
			    {
			      found=1;
			      break;
			    }
			}
		      if (found==0)
			{
			  strcpy(val[nf],psr->obsn[i].flagVal[k]);
			  nf++;
			  if (nf >= MAX_JUMPS)
			    {
			      printf("preProcessSimple3: Number of different flags in the .tim file > MAX_JUMPS (%d)\n",MAX_JUMPS);
			      exit(1);
			    }
			  
			}
		    }
		}
	    }
	}
      // Now create jumps
      for (l=1;l<nf;l++)
	{
	  psr->nJumps++;
	  sprintf(psr->jumpStr[psr->nJumps],"%s %s",psr->fjumpID,val[l]);
	  for (i=0;i<psr->nobs;i++)
	    {
	      for (k=0;k<psr->obsn[i].nFlags;k++)
		{
		  if (strcmp(psr->obsn[i].flagVal[k],val[l])==0)
		    psr->obsn[i].jump[(psr->obsn[i].obsNjump)++]=psr->nJumps;
		}
	    }
	  psr->fitJump[psr->nJumps]=1;
	}
    }
      
  // Now link phase jumps to a particular site-arrival-time
  for (i=0;i<psr->nPhaseJump;i++)
    {
      if (psr->phaseJumpID[i] = -1)
	{
	  // Find closest TOA
	  if (i==0) printf("WARNING: Use of phase jumps => .tim file must be sorted in time order\n");
	  for (k=1;k<psr->nobs;k++)
	    {
	      if ((double)psr->obsn[k].sat > (double)psr->phaseJump[i])
		{
		  psr->phaseJumpID[i]=k-1;
		  break;
		}		  
	    }
	}
    }
  logdbg("Complete preProcessSimple3");
}
