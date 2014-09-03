//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards
//
// Plugin to read a set of arrival time files and produce a list of corrections based on state changes


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
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "T2toolkit.h"
#include "tempo2.h"
#include "toasim.h"
#include "TKfit.h"

#define MAX_STEP 100

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  double epochs[MAX_OBSN];
  int i,nit,j,p;
  char fname[MAX_FILELEN];
  double globalParameter;
  long double result;
  long seed = TKsetSeed();
  
  //
  // For the output file
  //
  toasim_header_t* header;
  toasim_header_t* read_header;
  FILE* file;
  double offsets[MAX_OBSN]; // Will change to doubles - should use malloc
  double mean=0.0;
  // Create a set of corrections.
  toasim_corrections_t* corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));
  
  double f1_1;
  double f1_2;
  double t0 = 0;
  double t;
  double phaseStart=0;
  double phase0=0;
  double f0;
  double phase;
  double nudot;
  double nu;
  int nStep;
  double stepTime[MAX_STEP];
  int sw;
  int nSwitch;
  double nu0,initialT,nudot_now;
  double nuSwitch[MAX_STEP+1];
  double fa,fb;
  double deltaT;
  int ii;

  nStep = 0;


  corr->offsets=offsets;
  corr->params=""; // Normally leave as NULL. Can store this along with each realisation. 
  // Same length string in every iteration - defined in r_param_length see below
  corr->a0=0; // constant
  corr->a1=0; // a1*x
  corr->a2=0; // a2*x*X
  
  *npsr = 0;
  nit = 1;

  printf("Graphical Interface: add2state\n");
  printf("Author:              G. Hobbs\n");
  printf("Version:             1.0\n");
  
  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      
      if (strcmp(argv[i],"-nreal")==0){
	nit=atoi(argv[++i]);
      }
      if (strcmp(argv[i],"-step")==0){
	sscanf(argv[++i],"%lf",&stepTime[nStep]);
	nStep++;
      }
      if (strcmp(argv[i],"-fa")==0){
	sscanf(argv[++i],"%lf",&f1_1);
      }
      if (strcmp(argv[i],"-fb")==0){
	sscanf(argv[++i],"%lf",&f1_2);
      }
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[*npsr],argv[++i]); 
	  strcpy(timFile[*npsr],argv[++i]);
	  (*npsr)++;
	}
      if (strcmp(argv[i],"-seed")==0){
	sscanf(argv[++i],"%d",&seed);
      }
    }
  
  if (seed > 0)seed=-seed;
  
  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  // Now read in all the .tim files
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  
  preProcess(psr,*npsr,argc,argv);
  
  for (p=0;p<*npsr;p++)
    {
      f0 = (double)psr[p].param[param_f].val[0];


      printf("NTOA = %d\n",psr[p].nobs);
      header = toasim_init_header();
      strcpy(header->short_desc,"add2state");
      strcpy(header->invocation,argv[0]);
      strcpy(header->timfile_name,timFile[p]);
      strcpy(header->parfile_name,"Unknown");
      header->idealised_toas="NotSet"; // What should this be
      header->orig_parfile="NA";
      header->gparam_desc=""; // Global parameters
      header->gparam_vals="";
      header->rparam_desc=""; // Desciprtion of the parameters
      header->rparam_len=0; // Size of the string
      header->seed = seed;
		
      header->ntoa = psr[p].nobs;
      header->nrealisations = nit;
      
      // First we write the header...
      sprintf(fname,"%s.add2state",timFile[p]);
      file = toasim_write_header(header,fname);
      
      
      for (i=0;i<nit;i++)
	{
	  mean=0;
	  if(i%10 == 0){
	    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	    printf("Iteration %d/%d",i+1,nit);
	    fflush(stdout);
	  }
	  
	  t0 = (double)psr[p].obsn[0].sat;

	  for (j=0;j<psr[p].nobs;j++)
	    {
	      t= (psr[p].obsn[j].sat-t0)*86400.0;
	      // Calculate nudot at time t
	      sw=1;
	      for (ii=0;ii<nStep;ii++)
		{
		  if (psr[p].obsn[j].sat > stepTime[ii])
		    sw*=-1;
		}
	      if (sw == 1)
		nudot = f1_1;
	      else 
		nudot = f1_2;
	      
	      // Calculate nu at time t
	      // Integrate from time t0 to the current time
	      // Find necessary number of switches
	      nSwitch=0;
	      for (ii=0;ii<nStep;ii++)
		{
		  if (stepTime[ii] < psr[p].obsn[j].sat) {nSwitch++;}
		}
	      nu0 = f0;
	      initialT = t0;
	      nudot_now = f1_1;
	      sw=1;
	      nuSwitch[0] = nu0;
	      for (ii=0;ii<nSwitch;ii++)
		{
		  nu0 = nu0 + nudot_now*(stepTime[ii]-initialT)*86400.0;
		  nuSwitch[ii+1] = nu0;
		  sw*=-1;
		  if (sw==1)
		    nudot_now = f1_1;
		  else
		    nudot_now = f1_2;
		  initialT = stepTime[ii];
		}
	      // Now do the last bit
	      nu = nu0 + nudot_now*(psr[p].obsn[j].sat - initialT)*86400.0;
	      
	      // Now calculate phase at time t
	      // Dealing with all the switches
	      phase = phaseStart;
	      fb = f0;
	      initialT = t0;
	      nudot_now = f1_1;
	      sw=1;
	      
	      for (ii=0;ii<nSwitch;ii++)
		{
		  fa = nuSwitch[ii];
		  fb = nuSwitch[ii+1];
		  if (ii>0)
		    deltaT = stepTime[ii]-stepTime[ii-1];
		  else
		    deltaT = stepTime[ii] - t0;
		  sw*=-1;
		  if (sw==1)
		    nudot_now = f1_1;
		  else
		    nudot_now = f1_2;
		  
		  phase += (deltaT*86400.0*fb)+0.5*deltaT*86400.0*fabs(fa-fb);
		  initialT = stepTime[ii];
		}
	      // Now do last bit
	      deltaT = (psr[p].obsn[j].sat-initialT)*86400.0;
	      fa = fb;
	      fb = fa+deltaT*nudot_now;
	      phase += fb*deltaT+0.5*deltaT*fabs(fa-fb);
	      //	      phase = phase - floor(phase);
	      
	      printf("initial %g %g %d %.15f %.15f %g\n",(double)psr[p].obsn[j].sat,nudot,nSwitch,nu,phase,f0);
	      epochs[j] = (double)(psr[p].obsn[j].sat-psr[p].param[param_pepoch].val[0]);
	      offsets[j] = (double)(phase);
	      mean+=offsets[j];
	    }
	  for (j=0;j<psr[p].nobs;j++)
	    offsets[j]-=mean/(double)psr[p].nobs;
	  TKremovePoly_d(epochs,offsets,psr[p].nobs,3);

	  for (j=0;j<psr[p].nobs;j++)
	    printf("value = %g %.15f\n",epochs[j],offsets[j]);



  //			  //			  printf("Error = %g\n",(double)psr[p].obsn[j].toaErr);
  //			}
	    toasim_write_corrections(corr,header,file);

	  printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	  printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	  printf("Iteration %d/%d\n",i,nit);
	}	  
      printf("Close file\n");
      fclose(file);
    }
  return 0;
}


