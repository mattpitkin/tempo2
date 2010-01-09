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

void preProcess(pulsar *psr,int npsr,int argc,char *argv[])
{
  int p,i,k,l,nread,yes,fitN=0,setN=0,j;
  char fitStr[10][100];
  char setStr[10][100];
  float dmvals[10000];
  float startdmmjd;
  int ndm;
  longdouble setVal[10];
  FILE *fdmin;
  char newEpoch[100]="NONE";
  char selectFname[1000]="";
  char globalFname[1000]="";
  char select1[100][100];
  char select2[100][100];
  char line[MAX_STRLEN];
  char hashcheck;
  double select3[100];
  double select4[100];
  int nSelect=0;
  char str1[100],str2[100],str3[100],str4[100],str5[100];
  int v5;
  char name[100];
  char dmfile[100]="";
  int setName=0;
  double val1,val2;
  double last=-1;
  int tempo1=0;
  int nojump=0;
  int nofit=0;
  int modify=0;
  //trim data to match dm correction, but don't correct for dm
  int trimonly = 0;
  char modifyFname[100];
  double simulate=0;
  
  if (debugFlag==1) printf("In preProcess\n");

  //MAX_PSR   = MAX_PSR_VAL;    /* Maximum number of pulsars to fit simultaneously  */
  //MAX_OBSN  = MAX_OBSN_VAL;
  ECLIPTIC_OBLIQUITY = ECLIPTIC_OBLIQUITY_VAL;
  //  debugFlag = 0;

  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-epoch")==0)
	strcpy(newEpoch,argv[++i]);
      else if (strcmp(argv[i],"-last")==0)
	sscanf(argv[++i],"%lf",&last);
      else if (strcmp(argv[i],"-setdm")==0) 
	sscanf(argv[++i],"%s",dmfile); // Should deal with multiple pulsars
      else if (strcmp(argv[i],"-trimonly")==0)
	trimonly = 1;
      else if (strcmp(argv[i],"-tempo1")==0)
	tempo1=1;
      else if (strcmp(argv[i],"-nojump")==0)
	nojump=1;
      else if (strcmp(argv[i],"-select")==0)
	sscanf(argv[++i],"%s",selectFname);
      else if (strcmp(argv[i],"-global")==0)
	sscanf(argv[++i],"%s",globalFname);
      else if (strcmp(argv[i],"-modify")==0)
	{
	  modify=1;
	  sscanf(argv[++i],"%s",modifyFname);
	}
       else if (strcmp(argv[i],"-name")==0)
	{
	  setName=1;
	  strcpy(name,argv[i+1]);
	}
      else if (strcmp(argv[i],"-fit")==0)
	strcpy(fitStr[fitN++],argv[i+1]);
      else if (strcmp(argv[i],"-set")==0)
	{
	  strcpy(setStr[setN],argv[i+1]);
	  sscanf(argv[i+2],"%Lf",&setVal[setN]);
	  setN++;
	}
      else if (strcmp(argv[i],"-simulate")==0)
	{
	  sscanf(argv[i+1],"%lf",&simulate);
	}
      else if (strcmp(argv[i],"-nofit")==0)
	nofit = 1;
      else if (strcmp(argv[i],"-clock")==0)
	{
	  for (p=0;p<npsr;p++)
	    strcpy(psr[p].clock,argv[i+1]);
	}
    }
  for (p=0;p<npsr;p++)
    {
      for (i=0;i<MAX_PARAMS;i++){
	if(psr[p].param[i].nLinkTo>0){
	  psr[p].param[i].val[0] = getParameterValue(&psr[p],i,0);
	  psr[p].param[i].prefit[0] = getParameterValue(&psr[p],i,0);
	}
      }

      if (setName==1)
	strcpy(psr[p].name,name);
      if (nojump==1)
	psr[p].nJumps=0;

      if (nofit==1)
	{
	  for (i=0;i<MAX_PARAMS;i++)
	    {
	      for (k=0;k<psr[p].param[i].aSize;k++)
		psr[p].param[i].fitFlag[k] = 0;
	    }
	  // Turn off fitting for jumps
	  for (i=0;i<=psr[p].nJumps;i++)
	    psr[p].fitJump[i]=0;
	}
      /* Select command line fitting */
      for (i=0;i<fitN;i++)
	{
	  for (j=0;j<MAX_PARAMS;j++)
	    {
	      for (k=0;k<psr[p].param[j].aSize;k++)
		{
		  if (strcasecmp(fitStr[i],psr[p].param[j].shortlabel[k])==0)
		    {
		      if (psr[p].param[j].paramSet[k]!=1)
			{
			  psr[p].param[j].paramSet[k]=1;
			  psr[p].param[j].val[k]=0.0;
			  psr[p].param[j].prefit[k]=0.0;		      
			}
		      psr[p].param[j].fitFlag[k]=1;
		    }
		}
	    }
	}
      /* Simulate global parameter */
      if (simulate!=0.0)
	{
	  for (i=0;i<psr[p].nobs;i++)
	    {
	      psr[p].obsn[i].sat += simulate/SECDAY*sin(0.003*psr[p].obsn[i].sat);
	    }
	}

      /* Set jump values if already set */
      /* MOVED FOLLOWING INTO READPARFILE.C */
      /*      for (k=1;k<=psr[p].nJumps;k++)
	{	
	  v5 = -1;
	  nread = sscanf(psr[p].jumpStr[k],"%s %s %s %s %s",str1,str2,str3,str4,str5);

	  if (strcasecmp(str1,"MJD")==0 || strcasecmp(str1,"FREQ")==0)
	    {
	      if (nread>3)
		{
		  sscanf(str4,"%lf",&(psr[p].jumpVal[k]));
		  if (sscanf(str5,"%d",&v5)==1)
		    {
		      if (v5!=1) psr[p].fitJump[k]=0;
		    }
		  else
		    psr[p].fitJump[k]=0;
		}
	    }
	  else if (strcasecmp(str1,"NAME")==0 || strcasecmp(str1,"TEL")==0 || str1[0]=='-')
	    {
	      if (nread>2)
		{
		  sscanf(str3,"%lf",&(psr[p].jumpVal[k]));
		  if (sscanf(str4,"%d",&v5)==1)
		    {
		      if (v5!=1) psr[p].fitJump[k]=0;
		    }
		  else
		    psr[p].fitJump[k]=0;
		}
	    }
	    } */

      /* Select command line parameter setting */
      for (i=0;i<setN;i++)
	{
	  if (strcasecmp(setStr[i],"NITS")==0)
	    psr[p].nits = (int)setVal[i];
	  
	  for (j=0;j<MAX_PARAMS;j++)
	    {
	      for (k=0;k<psr[p].param[j].aSize;k++)
		{
		  if (strcasecmp(setStr[i],psr[p].param[j].shortlabel[k])==0)
		    {
		      psr[p].param[j].val[k]=setVal[i];
		      psr[p].param[j].prefit[k]=setVal[i];
		      psr[p].param[j].paramSet[k]=1;
		    }
		}
	    }
	}
      /* Check whitening */
      if (psr[p].param[param_wave_om].paramSet[0] == 1 && psr[p].param[param_wave_om].val[0] == 0.0) /* Set fundamental frequency */
	{  
	  longdouble first=-1,last=-1;
	  for (i=0;i<psr[p].nobs;i++)
	    {
	      if (psr[p].obsn[i].deleted==0 && first==-1) first=psr[p].obsn[i].sat;
	      if (psr[p].obsn[i].deleted==0 && last==-1)  last=psr[p].obsn[i].sat;
	      if (psr[p].obsn[i].deleted==0 && psr[p].obsn[i].sat < first) first = psr[p].obsn[i].sat;
	      if (psr[p].obsn[i].deleted==0 && psr[p].obsn[i].sat > last)   last = psr[p].obsn[i].sat;
	    }
	  //	  printf("WHITE: %Lg %d\n",(last-first)/365.25,psr[p].nWhite);
	  //	  psr[p].param[param_wave_om].val[0] = 2.0*M_PI/(last-first)/
	  //	    (1.0+4.0/(double)((psr[p].nWhite+1)*2.0));
	  psr[p].param[param_wave_om].val[0] = 2.0*M_PI/(last-first)/
	    (1.0+4.0/(double)(psr[p].nWhite));
	}

      /* Set tempo emulation mode */
      if (psr[p].param[param_ephver].paramSet[0]==1)
	{
	  if (psr[p].param[param_ephver].val[0] < 5)
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
	  psr[p].param[param_ephver].paramSet[0]=1;
	  psr[p].param[param_ephver].val[0]=2;
	}
      else 
	{
	  psr[p].param[param_ephver].paramSet[0]=1;
	  psr[p].param[param_ephver].val[0]=5;
	}
      if (last>0) /* Define the start parameter based upon the last observation */
	{
	  psr[p].param[param_start].val[0] = psr[p].obsn[psr[p].nobs-1].sat-last;
	  psr[p].param[param_start].fitFlag[0] = 1;
	  psr[p].param[param_start].paramSet[0] = 1;
	}

      psr[p].tempo1 = 0;
      if (tempo1 == 1)
      {
	psr[p].tempo1 = 1;
	psr[p].units = TDB_UNITS;
	psr[p].timeEphemeris = FB90_TIMEEPH;
	psr[p].dilateFreq = 0;
	psr[p].planetShapiro=0;
	psr[p].t2cMethod = T2C_TEMPO;
	psr[p].correctTroposphere = 0;
	psr[p].ne_sw = 9.961;
	ECLIPTIC_OBLIQUITY = 84381.412;
      }
      /* XXXX Hack!! -- removed - should use transform plugin*/
      /* Problem, if a function is not called then the compiler does not include
       * it in the library - so the plugin cannot find it */
      if (psr[p].units==100) /* SI_UNITS)  */
	transform_units(&psr[p], TDB_UNITS, SI_UNITS); 

      /* Update period epoch if necessary */
      if (strcmp(newEpoch,"NONE")!=0)
	{	  
	  longdouble nMJD,dt;
	  longdouble earliest=-1;
	  longdouble latest=-1;
	  int okay=1;

	  //=psr[p].obsn[0].sat,latest=psr[p].obsn[0].sat;
	  for (i=0;i<psr[p].nobs;i++)
	    {
	      okay=1;
	      if (psr[p].obsn[i].deleted==1) okay=0;
	      if (psr[p].param[param_start].paramSet[0]==1 && psr[p].param[param_start].fitFlag[0]==1 &&
		  (psr[p].param[param_start].val[0] > psr[p].obsn[i].sat))
		okay=0;
	      if (psr[p].param[param_finish].paramSet[0]==1 && psr[p].param[param_finish].fitFlag[0]==1 &&
		  psr[p].param[param_finish].val[0] < psr[p].obsn[i].sat)
		okay=0;
	
	      if (okay==1)
		{
		  if (earliest==-1)
		    {
		      earliest = psr[p].obsn[i].sat;
		      latest   = psr[p].obsn[i].sat;
		    }
		  if (earliest > psr[p].obsn[i].sat) earliest = psr[p].obsn[i].sat;
		  if (latest < psr[p].obsn[i].sat) latest = psr[p].obsn[i].sat;
		}
	    }
	

	  if (strcasecmp(newEpoch,"CENTRE")==0 ||
	      strcasecmp(newEpoch,"CENTER")==0) /* Find centre of data */
	    nMJD = (int)((earliest+latest)/2.0);
	  else if (strcasecmp(newEpoch,"LEFT")==0)
	    nMJD = (int)(earliest);
	  else if (strcasecmp(newEpoch,"RIGHT")==0)
	    nMJD = (int)(latest);
	  else
	    nMJD = parse_longdouble(newEpoch);

	  dt = (nMJD - psr[p].param[param_pepoch].val[0])*86400.0;
	  psr[p].param[param_f].val[0] = psr[p].param[param_f].val[0]+
	  psr[p].param[param_f].val[1]*dt + 0.5*psr[p].param[param_f].val[2]*dt*dt;
	  psr[p].param[param_f].val[1] += psr[p].param[param_f].val[2]*dt;
	  if (psr[p].param[param_f].paramSet[3]==1)
	    {
	      psr[p].param[param_f].val[0] += 1.0/6.0*psr[p].param[param_f].val[3]*dt*dt*dt;
	      psr[p].param[param_f].val[1] += 0.5*psr[p].param[param_f].val[3]*dt*dt;
	      psr[p].param[param_f].val[2] += psr[p].param[param_f].val[3]*dt;
	    }
	      
	  psr[p].param[param_f].prefit[0] = psr[p].param[param_f].val[0];
	  psr[p].param[param_f].prefit[1] = psr[p].param[param_f].val[1];
	  psr[p].param[param_pepoch].val[0] = nMJD;
	  psr[p].param[param_pepoch].prefit[0] = nMJD;

	  /* Update position epoch */
	  dt = (nMJD - psr[p].param[param_posepoch].val[0])/365.25;

	  if (psr[p].param[param_pmra].paramSet[0]==1)
	    {
	      char retstr[1000];
	      printf("Updating RAJ\n");
	      psr[p].param[param_raj].val[0] = psr[p].param[param_raj].val[0]+psr[p].param[param_pmra].val[0]
		/cos(psr[p].param[param_decj].val[0])/1000.0*(M_PI/180.0)/60.0/60.0*dt; 
	      psr[p].param[param_raj].prefit[0] = psr[p].param[param_raj].val[0];
	      /* Must obtain this in hms form */
	      turn_hms(psr[p].param[param_raj].val[0]/(2.0*M_PI), retstr);
	      strcpy(psr[p].rajStrPost,retstr);
	      strcpy(psr[p].rajStrPre,retstr);
	    }

	  if (psr[p].param[param_pmdec].paramSet[0]==1)
	    {
	      char retstr[1000];
	      psr[p].param[param_decj].val[0] = psr[p].param[param_decj].val[0]+psr[p].param[param_pmdec].val[0]/1000.0*(M_PI/180.0)/60.0/60.0*dt; 
	      psr[p].param[param_decj].prefit[0] = psr[p].param[param_decj].val[0];
	      turn_dms(psr[p].param[param_decj].val[k]/(2.0*M_PI), retstr);
	      strcpy(psr[p].decjStrPost,retstr);
	      strcpy(psr[p].decjStrPre,retstr);
	    }

	  psr[p].param[param_posepoch].val[0] = nMJD;
	  psr[p].param[param_posepoch].prefit[0] = nMJD;

	  /* Update dmepoch */
	  dt = (nMJD - psr[p].param[param_dmepoch].val[0])/365.25;
	  psr[p].param[param_dm].val[0] = psr[p].param[param_dm].val[0]+
	  psr[p].param[param_dm].val[1]*dt + 0.5*psr[p].param[param_dm].val[2]*dt*dt;
	  psr[p].param[param_dmepoch].val[0] = nMJD;
	  psr[p].param[param_dmepoch].prefit[0] = nMJD;


	  /* Update binary parameters if necessary */
	  if (psr[p].param[param_pb].paramSet[0]==1)  /* Binary pulsar */
	    {
	      longdouble orbits,pb,tt0,pbdot,xpbdot,t0p,t0m=0.0;
	      int        norbits;

	      if (psr[p].param[param_pbdot].paramSet[0]==1) pbdot = psr[p].param[param_pbdot].val[0];
	      else pbdot = 0.0;

	      if (psr[p].param[param_xpbdot].paramSet[0] == 1) xpbdot = psr[p].param[param_xpbdot].val[0];
	      else xpbdot = 0.0;

	      pb = psr[p].param[param_pb].val[0]*SECDAY;
	      
	      if (psr[p].param[param_tasc].paramSet[0]==1)
		tt0 = (nMJD-psr[p].param[param_tasc].val[0])*SECDAY;
	      else
		tt0 = (nMJD-psr[p].param[param_t0].val[0])*SECDAY;
	      printf("tt0 = %.14Lf %.14Lf %.14Lf %.14Lf\n",tt0,nMJD,psr[p].param[param_t0].val[0],
		     nMJD-psr[p].param[param_t0].val[0]);

	      orbits = tt0/pb - 0.5*(pbdot+xpbdot)*(tt0/pb)*(tt0/pb);
	      norbits = (int)(orbits+0.5);
	      printf("Orbits = %.5Lf  (%d)\n",orbits,norbits);

	      if (xpbdot > 0 || pbdot > 0)
		{
		  /*		  t0p = (1.0/pb + sqrtl(1.0/(pb*pb)-2.0*(pbdot+xpbdot)*norbits/pb/pb))
		    /((pbdot+xpbdot)/pb/pb);
		    t0m = (1.0/pb - sqrtl(1.0/(pb*pb)-2.0*(pbdot+xpbdot)*norbits/pb/pb))/((pbdot+xpbdot)/pb/pb); */
		  t0p = pb/(pbdot+xpbdot)*(1.0+sqrtl(1.0-2.0*(pbdot+xpbdot)*(longdouble)norbits));
		  t0m = pb/(pbdot+xpbdot)*(1.0-sqrtl(1.0-2.0*(pbdot+xpbdot)*(longdouble)norbits));

		  if (psr[p].param[param_tasc].paramSet[0]==1)
		    {
		      t0p = psr[p].param[param_tasc].val[0]+t0p/SECDAY;
		      t0m = psr[p].param[param_tasc].val[0]+t0m/SECDAY;     
		      
		      if (fabs(t0p-psr[p].param[param_tasc].val[0]) > fabs(t0m-psr[p].param[param_tasc].val[0]))
			t0p = t0m;
		    }
		  else
		    {
		      t0p = psr[p].param[param_t0].val[0]+t0p/SECDAY;
		      t0m = psr[p].param[param_t0].val[0]+t0m/SECDAY;     
		    
		      if (fabs(t0p-psr[p].param[param_t0].val[0]) > fabs(t0m-psr[p].param[param_t0].val[0]))
			t0p = t0m;
		    }
		}
	      else
		{
		  t0p = norbits*pb;
		  if (psr[p].param[param_tasc].paramSet[0]==1)		  
		    t0p = psr[p].param[param_tasc].val[0]+t0p/SECDAY;
		  else
		    t0p = psr[p].param[param_t0].val[0]+t0p/SECDAY;
		}
	       
	      if (psr[p].param[param_tasc].paramSet[0]==1)		  
		{
		  psr[p].param[param_tasc].val[0] = t0p;
		  psr[p].param[param_tasc].prefit[0] = psr[p].param[param_tasc].val[0];
		}
	      else
		{
		  psr[p].param[param_t0].val[0] = t0p;
		  psr[p].param[param_t0].prefit[0] = psr[p].param[param_t0].val[0];
		}
	      printf("Result = %.14Lf %.14Lf\n",t0p,t0m);
	      
	      psr[p].param[param_pb].val[0] += psr[p].param[param_pbdot].val[0]*tt0/SECDAY;
	      psr[p].param[param_pb].prefit[0] = psr[p].param[param_pb].val[0];

	      psr[p].param[param_om].val[0] += psr[p].param[param_omdot].val[0]*tt0/SECDAY/365.25;
	      psr[p].param[param_om].prefit[0] = psr[p].param[param_om].val[0];

	      if (psr[p].param[param_a1dot].paramSet[0]==1)
		{
		  psr[p].param[param_a1].val[0] += psr[p].param[param_a1dot].val[0]*1.0e-12*tt0;
		  psr[p].param[param_a1].prefit[0] = psr[p].param[param_a1].val[0];
		}



	    }
	}    

      if (psr[p].param[param_pepoch].paramSet[0]==1 && psr[p].param[param_pepoch].fitFlag[0]==1)
	{printf("Warning: Cannot fit for pepoch\n"); psr[p].param[param_pepoch].fitFlag[0]=0;}
      if (psr[p].param[param_posepoch].paramSet[0]==1 && psr[p].param[param_posepoch].fitFlag[0]==1)
	{printf("Warning: Cannot fit for posepoch\n"); psr[p].param[param_posepoch].fitFlag[0]=0;}
      if (psr[p].param[param_dmepoch].paramSet[0]==1 && psr[p].param[param_dmepoch].fitFlag[0]==1)
	{printf("Warning: Cannot fit for dmepoch\n"); psr[p].param[param_dmepoch].fitFlag[0]=0;}
      if (psr[p].param[param_track].paramSet[0]==1 && psr[p].param[param_track].fitFlag[0]==1)
	{printf("Warning: Cannot fit for track\n"); psr[p].param[param_track].fitFlag[0]=0;}
      if (strlen(dmfile)>0)
	{
	  float tt;
	  fdmin = fopen(dmfile,"r");
	  ndm=0;
	  while (!feof(fdmin))
	    {
	      //check wheter a line is a comment (= starts with a hash)
	      fgets(line,MAX_STRLEN,fdmin);
	      sscanf(line,"%c",&hashcheck);
	      if (hashcheck == '#') {
		//do nothing, perhaps give a debug message
		if (debugFlag)
		  printf("preProces():skipping line in dmfile\n");
	      }
	      else if (sscanf(line,"%f %f",&tt,&dmvals[ndm])==2)
		{
		  if (ndm==0)
		    startdmmjd = tt;
		  ndm++;
		}
	    }
	  fclose(fdmin);
	}
      
      // Modify TOA flags if required
      if (modify==1)
	{
	  FILE *fin;
	  double mjd1,mjd2;
	  char flag1[MAX_STRLEN],flag2[MAX_STRLEN],flag3[MAX_STRLEN];
	  if (!(fin = fopen(modifyFname,"r")))
	    {
	      printf("Unable to open >%s< to modify the flags\n",modifyFname);
	      exit(1);
	    }
	  while (!feof(fin))
	    {
	      if (fscanf(fin,"%s %s %lf %lf %s",flag1,flag2,&mjd1,&mjd2,flag3)==5)
		{
		  for (i=0;i<psr[p].nobs;i++)
		    {
		      for (j=0;j<psr[p].obsn[i].nFlags;j++)
			{
			  if (strcmp(psr[p].obsn[i].flagID[j],flag1)==0 &&
			      strcmp(psr[p].obsn[i].flagVal[j],flag2)==0 &&
			      (double)psr[p].obsn[i].sat > mjd1 &&
			      (double)psr[p].obsn[i].sat < mjd2)
			    strcpy(psr[p].obsn[i].flagVal[j],flag3);
			}
		    }
		}
	    }
	  fclose(fin);
	}
      
      for (i=0;i<psr[p].nobs;i++)
	{
	  //psr[p].obsn[i].efac = 1.0;
	  //psr[p].obsn[i].equad = 0.0;
	  /* Check whether any error = 0 */
	  if (psr[p].fitMode==1 && psr[p].obsn[i].toaErr==0)
	    {
	      printf("Error: the uncertainty for TOA %d (MJD %f) is equal to zero\n ",i+1,(double)psr[p].obsn[i].sat);
	      exit(1);
	    }

	  for (k=0; k<6; k++) /* zero these vectors for first call to tt2tb */
	  {
	    psr[p].obsn[i].earthMoonBary_ssb[k] = 0.0;
	    psr[p].obsn[i].earthMoonBary_earth[k] = 0.0;
	    psr[p].obsn[i].observatory_earth[k] = 0.0;
	    psr[p].obsn[i].earth_ssb[k] = 0.0;
	  }


	  /*	  psr[p].obsn[i].sat += psr[p].obsn[i].phaseOffset/psr[p].param[param_f0].val[0]; */
	  for (k=0;k<psr[p].nToffset;k++) /* Calculate time offsets */
	    {
	      char offsetSite[256], obsSite[256];
	      lookup_observatory_alias(psr[p].tOffsetSite[k], offsetSite);
	      lookup_observatory_alias(psr[p].obsn[i].telID, obsSite);
	      if ((psr[p].tOffset_f1[k]==0.0 || (psr[p].obsn[i].freq > psr[p].tOffset_f1[k])) &&
		  (psr[p].tOffset_f2[k]==0.0 || (psr[p].obsn[i].freq < psr[p].tOffset_f2[k])) &&
		  (psr[p].tOffset_t1[k]==0.0 || (psr[p].obsn[i].sat  > psr[p].tOffset_t1[k])) &&
		  (psr[p].tOffset_t2[k]==0.0 || (psr[p].obsn[i].sat  < psr[p].tOffset_t2[k])) && 
		  (strcmp(offsetSite,"0")==0 || 
 		   strcmp(offsetSite,obsSite)==0))
		 {
		   int use=1;
		   /* Check for flags */
		   if (strlen(psr[p].tOffsetFlags[k])>0)
		     {
		       char *myStr,str1[1000],flagID[100];
		       use=0;
		       strcpy(str1,psr[p].tOffsetFlags[k]);
		       myStr = strtok(str1," ");
		       do {
			 if (myStr[0]=='-' && (myStr[1]<48 || myStr[1]>57)) /* Look for flag */
			   {
			     strcpy(flagID,myStr);
			     myStr = strtok(NULL," ");
			     for (l=0;l<psr[p].obsn[i].nFlags;l++)
			       {
				 if (strcmp(flagID,psr[p].obsn[i].flagID[l])==0
				     && strcmp(myStr,psr[p].obsn[i].flagVal[l])==0)
				   {
				     use=1;
				     break;
				   }
			       }
			   }
		       } while ((myStr = strtok(NULL," "))!=NULL);
		     }
		   if (use==1) {psr[p].obsn[i].sat += psr[p].tOffset[k]/SECDAY;}
		 }  
	    }
	  /* Check filtering */
	  if (strlen(psr[p].filterStr)>0)
	    {
	      char flag[100],filtS[1000],*filt;
	      strcpy(filtS,psr[p].filterStr);
	      filt = strtok(filtS," ");
	      strcpy(flag,filt);
	      while (filt != NULL)
		{
		  if (filt[0]=='-') strcpy(flag,filt);
		  else
		    {
		      for (j=0;j<psr[p].obsn[i].nFlags;j++)
			{
			  if (strcmp(psr[p].obsn[i].flagID[j],flag)==0 &&
			      strcmp(psr[p].obsn[i].flagVal[j],filt)==0)
			    psr[p].obsn[i].deleted=1;
			}
		    }
		  filt = strtok(NULL," ");
		}
	    }

	  /* Check filtering for 'pass' */
	  if (strlen(psr[p].passStr)>0)
	    {
	      char flag[100],filtS[1000],*filt;
	      int found;

       	      strcpy(filtS,psr[p].passStr);
	      filt = strtok(filtS," ");
	      found=0;
	      while ( filt != NULL && found==0)
		{
		  if (filt[0]=='-') // Have a flag
		    strcpy(flag,filt);
		  else
		    {
		      for (j=0;j<psr[p].obsn[i].nFlags;j++)
			{
			  if (strcmp(psr[p].obsn[i].flagID[j],flag)==0 &&
			      strcmp(psr[p].obsn[i].flagVal[j],filt)==0)
			    {
			      found=1;
			      break;
			    }
			}
		    }
		  filt = strtok(NULL," ");
		}

	      if (found==0)
		psr[p].obsn[i].deleted=1;		 
	    }

	  /* CHECK JUMPS */	  
	  for (k=1;k<=psr[p].nJumps;k++)
	    {
	      yes=0;

	      /* Must parse jump string to determine whether this observation should jump or not */	      
	      nread = sscanf(psr[p].jumpStr[k],"%s %s %s",str1,str2,str3);
	      if (strcasecmp(str1,"MJD")==0 || strcasecmp(str1,"FREQ")==0)
		{
		  sscanf(str2,"%lf",&val1);
		  sscanf(str3,"%lf",&val2);

		  if (strcasecmp(str1,"MJD")==0) {
		    if (psr[p].obsn[i].sat >= val1 && psr[p].obsn[i].sat < val2) {
		      yes=1;
		    }
		  }
		  else if (strcasecmp(str1,"FREQ")==0) {
		    if (psr[p].obsn[i].freq >= val1 && psr[p].obsn[i].freq < val2) {
		      yes=1;
		    }
		  }
		}
	      if (strcasecmp(str1,"NAME")==0 && strstr(psr[p].obsn[i].fname,str2)!=NULL)
		yes=1;	
	      if (strcasecmp(str1,"TEL")==0)
	      {
		char selectedSite[256], obsSite[256];

		lookup_observatory_alias(str2, selectedSite);
		lookup_observatory_alias(psr[p].obsn[i].telID, obsSite);
		if (strcasecmp(selectedSite, obsSite)!=0)
		  yes=1;	    
	      }
	      else if (str1[0]=='-')
		{
		  int jj;
		  for (jj=0;jj<psr[p].obsn[i].nFlags;jj++)
		    {
		      if (strcmp(psr[p].obsn[i].flagID[jj],str1)==0 &&
			  strcmp(psr[p].obsn[i].flagVal[jj],str2)==0)
			yes=1;
		    }
		}
	    
	      if (yes==1) psr[p].obsn[i].jump=k;
	    }
	  // Check for time offset flags
	  for (k=0;k<psr[p].obsn[i].nFlags;k++)
	    {
	      if (strcmp(psr[p].obsn[i].flagID[k],"-to")==0)
		{
		  long double v;
		  sscanf(psr[p].obsn[i].flagVal[k],"%Lf",&v);
		  psr[p].obsn[i].sat += v/SECDAY;
		}
	    }
	  // Check for dm updates
	  if (strlen(dmfile)>0)
	    {
	      //float dm; //not needed any more
	      if ((int)(psr[p].obsn[i].sat-startdmmjd+0.5) > 0)
		{
		  if ((int)(psr[p].obsn[i].sat-startdmmjd+0.5) >= ndm)
		    psr[p].obsn[i].deleted=1;
		  else if (trimonly == 0) // set flags only if really correcting for dm, not only trimming the data
		  {
		    strcpy(psr[p].obsn[i].flagID[psr[p].obsn[i].nFlags],"-dm");
		    sprintf(psr[p].obsn[i].flagVal[psr[p].obsn[i].nFlags],"%g",dmvals[(int)(psr[p].obsn[i].sat-startdmmjd+0.5)]);
		    psr[p].obsn[i].nFlags++;
		    //dm = dmvals[(int)(psr[p].obsn[i].sat-startdmmjd+0.5)];             
		  }
		}
	      else
		{
		  psr[p].obsn[i].deleted=1;
		}
	    }
	}
  // Check for select file
      if (strlen(selectFname) > 0)
	{
	  int k,l;
	  int okay=0,found=0;
	  double low,high,tdiff;
	  char str[100],str1[100],str2[100],str3[100],str4[100],str5[100];
	  
	  FILE *fin;
	  if (!(fin = fopen(selectFname,"r")))
	    {
	      printf("Unable to open select file: >%s<\n",selectFname);
	      exit(1);
	    }
	  while (!feof(fin))
	    {
	      //		  if (fscanf(fin,"%s %s %lf %lf",select1[nSelect],select2[nSelect],
	      //			     &select3[nSelect],&select4[nSelect])==4)
	      //		    nSelect++;
	      if (fscanf(fin,"%s",str)==1)
		{
		  if (strcasecmp(str,"ONLY")==0)
		    {
		      fscanf(fin,"%s",str1);
		      if (strcasecmp(str1,"SIMUL")==0)
			{
			  nread=fscanf(fin,"%s %s %s %lf",str2,str3,str4,&tdiff);
			  if (nread==3) tdiff=60; // Default to 60 seconds
			  for (i=0;i<psr[p].nobs;i++)
			    {
			      found=0;
			      for (l=0;l<psr[p].obsn[i].nFlags;l++)		      
				{
				  if (strcmp(psr[p].obsn[i].flagID[l],str2)==0 &&
				      strcmp(psr[p].obsn[i].flagVal[l],str3)==0)
				    {
				      for (j=0;j<psr[p].nobs;j++)
					{
					  if (i!=j && strcmp(psr[p].obsn[j].flagID[l],str2)==0 &&
					      strcmp(psr[p].obsn[j].flagVal[l],str4)==0 &&
					      fabs(psr[p].obsn[i].sat - psr[p].obsn[j].sat)<=tdiff/SECDAY)
					    {
					      psr[p].obsn[j].deleted=-2;
					      found=1;
					      break;
					    }
					}
				    }
				}
			      if (found==0 && psr[p].obsn[i].deleted!=-2)
				psr[p].obsn[i].deleted=1;
			      else if (psr[p].obsn[i].deleted==-2)
				psr[p].obsn[i].deleted=0;
			    }
			  for (i=0;i<psr[p].nobs;i++)
			    {
			      if (psr[p].obsn[i].deleted==-2)
				psr[p].obsn[i].deleted=0;
			    }
			}
		    }
		  else if (str[0]=='-') // Filter flags
		    {
		      fscanf(fin,"%s %lf %lf",str1,&low,&high);
		      for (i=0;i<psr[p].nobs;i++)
			{
			  found=0;
			  for (l=0;l<psr[p].obsn[i].nFlags;l++)		      
			    {
			      if (strcmp(psr[p].obsn[i].flagID[l],str)==0 &&
				  strcmp(psr[p].obsn[i].flagVal[l],str1)==0)
				{
				  found=-1;
				  if (psr[p].obsn[i].sat >= low && 
				      psr[p].obsn[i].sat < high)
				    {
				      found=1;
				      break;
				    }		
				}	     			      
			    }
			  if (found==-1)
			    psr[p].obsn[i].deleted=1;
			}
		    }
		  else if (strcasecmp(str,"FILTER")==0)
		    {
		      fscanf(fin,"%s",str);
		      if (strcasecmp(str,"FREQ")==0)
			{
			  if (fscanf(fin,"%lf %lf",&low,&high)==2)
			    {
			      for (i=0;i<psr[p].nobs;i++)
				{
				  if (psr[p].obsn[i].freq >= low && psr[p].obsn[i].freq < high)
				    psr[p].obsn[i].deleted=1;
				}
			    }
			}
		      else if (strcasecmp(str,"MJD")==0)
			{
			  if (fscanf(fin,"%lf %lf",&low,&high)==2)
			    {
			      for (i=0;i<psr[p].nobs;i++)
				{
				  if ((double)psr[p].obsn[i].sat >= low && (double)psr[p].obsn[i].sat < high)
				    psr[p].obsn[i].deleted=1;
				}
			    }
			}
		      else if (strcasecmp(str,"TOAERR")==0)
			{
			  if (fscanf(fin,"%lf %lf",&low,&high)==2)
			    {
			      for (i=0;i<psr[p].nobs;i++)
				{
				  if ((double)psr[p].obsn[i].toaErr >= low && (double)psr[p].obsn[i].toaErr < high)
				    psr[p].obsn[i].deleted=1;
				}
			    }
			}
		    }
		  if (strcasecmp(str,"PASS")==0)
		    {
		      fscanf(fin,"%s",str1);
		      if (strcasecmp(str1,"SIMUL")==0)
			{
			  nread=fscanf(fin,"%s %s %s %lf",str2,str3,str4,&tdiff);
			  if (nread==3) tdiff=60; // Default to 60 seconds
			  for (i=0;i<psr[p].nobs;i++)
			    {
			      found=0;
			      for (l=0;l<psr[p].obsn[i].nFlags;l++)		      
				{
				  if (strcmp(psr[p].obsn[i].flagID[l],str2)==0 &&
				      strcmp(psr[p].obsn[i].flagVal[l],str3)==0)
				    {
				      for (j=0;j<psr[p].nobs;j++)
					{
					  if (i!=j && strcmp(psr[p].obsn[j].flagID[l],str2)==0 &&
					      strcmp(psr[p].obsn[j].flagVal[l],str4)==0 &&
					      fabs(psr[p].obsn[i].sat - psr[p].obsn[j].sat)<=tdiff/SECDAY)
					    {
					      if (psr[p].obsn[i].toaErr > psr[p].obsn[j].toaErr)
						psr[p].obsn[i].deleted=1;
					      else
						psr[p].obsn[j].deleted=1;
					      break;
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	  fclose(fin);
	}
      
      //
      // Check fjump
      //
      if (strlen(psr[p].fjumpID)>0)
	{
	  char val[MAX_FLAGS][16];
	  int  nval[MAX_FLAGS];
	  int  nf=0;
	  int  k,l;
	  int found;
	  // Find flag corresponding to longest data set
	  for (i=0;i<psr[p].nobs;i++)
	    {
	      if (psr[p].obsn[i].deleted==0)
		{
		  for (k=0;k<psr[p].obsn[i].nFlags;k++)
		    {
		      if (strcmp(psr[p].obsn[i].flagID[k],psr[p].fjumpID)==0)
			{
			  found=0;
			  for (l=0;l<nf;l++)
			    {
			      if (strcmp(psr[p].obsn[i].flagVal[k],val[l])==0)
				{
				  found=1;
				  break;
				}
			    }
			  if (found==0)
			    {
			      strcpy(val[nf],psr[p].obsn[i].flagVal[k]);
			      nf++;
			    }
			}
		    }
		}
	    }
	  // Now create jumps
	  for (l=1;l<nf;l++)
	    {
	      psr[p].nJumps++;
	      sprintf(psr[p].jumpStr[psr[p].nJumps],"%s %s",psr[p].fjumpID,val[l]);
	      for (i=0;i<psr[p].nobs;i++)
		{
		  for (k=0;k<psr[p].obsn[i].nFlags;k++)
		    {
		      if (strcmp(psr[p].obsn[i].flagVal[k],val[l])==0)
			  psr[p].obsn[i].jump=psr[p].nJumps;
		    }
		}
	      psr[p].fitJump[psr[p].nJumps]=1;
	    }
	}
    }





  // Now check for global parameters
  if (strlen(globalFname) > 0)
    {
      FILE *fin;
      char tpar[MAX_STRLEN][MAX_FILELEN];
      char ttim[MAX_STRLEN][MAX_FILELEN];

      //      fin = fopen(globalFname,"r");      
      printf("Setting global parameters\n");
      strcpy(tpar[0],globalFname);
      for (p=0;p<npsr;p++)
	{
	  readParfile(psr+p,tpar,ttim,1);
	}
      //      fclose(fin);
    }

} 
