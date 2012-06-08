//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards
//
// Plugin to read a set of arrival time files and produce the offsets required to simulated a gravitational wave background 

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
#include "GWsim.h"

using namespace std;

long double getTspan(pulsar *psr,int npsr);

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,nit,j,p,k;
  char fname[MAX_FILELEN];
  double globalParameter;
  long double result;
  long seed = TKsetSeed();
  // For the simulation
  gwSrc *gw;
  long double timeOffset; 
  long double scale;
  long double alpha = -2.0/3.0;
  long double gwAmp = 1e-14;
  long double ra_p,dec_p;
  long double flo=0.0,fhi=0.0;
  long double kp[3];            /* Vector pointing to pulsar           */
  long double tspan;
  long double time;
  long double gwRes[MAX_OBSN];
  long double dist[MAX_PSR];
  long double mean;
  int distNum=0;
  int logspacing=1;
  int ngw=1000;
  char readGW=0;
  char writeGW=0;
  char gwFileName[MAX_FILELEN];
  FILE *gwFile;


  fname[0]='\0';


  //
  // For the output file
  //
  toasim_header_t* header;
  toasim_header_t* read_header;
  FILE* file;
  double offsets[MAX_OBSN]; // Will change to doubles - should use malloc
  double epochs[MAX_OBSN]; // Will change to doubles - should use malloc
  // Create a set of corrections.
  toasim_corrections_t* corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));

  corr->offsets=offsets;
  corr->params=""; // Normally leave as NULL. Can store this along with each realisation. 
                   // Same length string in every iteration - defined in r_param_length see below
  corr->a0=0; // constant
  corr->a1=0; // a1*x
  corr->a2=0; // a2*x*X
  
  *npsr = 0;
  nit = 1;

  printf("Graphical Interface: addGWB\n");
  printf("Author:              G. Hobbs, M. Keith\n");
  printf("Version:             1.0\n");

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
	    if (strcmp(argv[i],"-nreal")==0){
		    nit=atoi(argv[++i]);
	    }

	    else if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[*npsr],argv[++i]); 
	  strcpy(timFile[*npsr],argv[++i]);
	  (*npsr)++;
	}
      else if (strcmp(argv[i],"-dist")==0) // Distance in kpc
	{
	  sscanf(argv[++i],"%Lf",&dist[distNum]);
	  dist[distNum]*=3.086e19;
	  distNum++;
	}
      else if (strcmp(argv[i],"-gwamp")==0)
	{sscanf(argv[++i],"%Lf",&gwAmp);}
      else if (strcmp(argv[i],"-alpha")==0)
	{sscanf(argv[++i],"%Lf",&alpha);}
      else if (strcmp(argv[i],"-ngw")==0)
	{sscanf(argv[++i],"%d",&ngw);}
      else if (strcmp(argv[i],"-flo")==0)
	sscanf(argv[++i],"%Lf",&flo);
      else if (strcmp(argv[i],"-fhi")==0)
	sscanf(argv[++i],"%Lf",&fhi);
      else if (strcmp(argv[i],"-seed")==0)
	sscanf(argv[++i],"%d",&seed);
      else if (strcmp(argv[i],"-readGW")==0){
	sscanf(argv[++i],"%s",&gwFileName);
	readGW=1;
      } else if (strcmp(argv[i],"-outf")==0){
	      sscanf(argv[++i],"%s",&fname);
      } else if (strcmp(argv[i],"-writeGW")==0){
	sscanf(argv[++i],"%s",&gwFileName);
	writeGW=1;
      }

    }

  if (seed > 0)seed=-seed;
  scale = pow(86400.0*365.25,alpha);
  gwAmp *= scale;
  if (distNum!=*npsr)
    {
      printf("ERROR: Distances not provided for all the pulsars: Npsr = %d, Ndist = %d\nUse -dist to provide distances (in kpc) on the command line.\n",*npsr,distNum);
      exit(1);
    }

  if ((gw = (gwSrc *)malloc(sizeof(gwSrc)*ngw))==NULL)
    {
      printf("Unable to allocate memory for %d GW sources\n",ngw);
      exit(1);
    }

  if(readGW){
	  gwFile = fopen(gwFileName,"r");
  }
  if(writeGW){
	  gwFile=fopen(gwFileName,"w");
  }

  // Now read in all the .tim files
  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

  // Set range of frequencies for GW simulation
  tspan=getTspan(psr,*npsr)*SECDAY;
  if (flo==0)
    {
      flo=0.01/tspan;
      printf("flo = %.5Lg\n",flo);
    }
  if (fhi==0)
    {
      fhi = 1.0/(long double)SECDAY;
      printf("fhi = %.5Lg\n",fhi);
    }
//  timeOffset = psr[0].param[param_pepoch].val[0];
    timeOffset=56000; // this needs to be the same for all pulsars!


  for (p=0;p<*npsr;p++)
    {
      ra_p   = psr[p].param[param_raj].val[0]; // Get position of the pulsar
      dec_p  = psr[p].param[param_decj].val[0];
      setupPulsar_GWsim(ra_p,dec_p,kp);

      header = toasim_init_header();
      strcpy(header->short_desc,"addGWB");
      strcpy(header->invocation,argv[0]);
      strcpy(header->timfile_name,timFile[p]);
      header->idealised_toas="NA"; // What should this be
      header->gparam_desc=""; // Global parameters
      header->gparam_vals="";
      header->rparam_desc=""; // Desciprtion of the parameters
      header->rparam_len=0; // Size of the string
      header->seed = seed;

      header->ntoa = psr[p].nobs;
      header->nrealisations = nit;

      // First we write the header...
      if(fname[0]=='\0')
	      sprintf(fname,"%s.addGWB",timFile[p]);
      file = toasim_write_header(header,fname);

      for (j=0;j<nit;j++)
	{
		if(readGW){
			ngw=GWbackground_read(gw,gwFile,j);
			if(ngw<1)exit(1);
		} else {
			GWbackground(gw,ngw,&seed,flo,fhi,gwAmp,alpha,logspacing);
			if(writeGW){
				GWbackground_write(gw,gwFile,ngw,j);
			}
		}
	  for (i=0;i<ngw;i++)
	    setupGW(&gw[i]);
	  
	  mean = 0.0;

	  printf("Iteration %d/%d\n",j+1,nit);
	  for (i=0;i<psr[p].nobs;i++)
	    {
	      time = (psr[p].obsn[i].sat - timeOffset)*SECDAY;
	      gwRes[i] = 0.0;
	      for (k=0;k<ngw;k++)
		gwRes[i]+=calculateResidualGW(kp,&gw[k],time,dist[p]);
	      mean += gwRes[i];
	    }
	  mean /= (double)psr[p].nobs;

	  for (i=0;i<psr[p].nobs;i++)
	    {
              epochs[i]=(double)psr[p].obsn[i].sat;
	      offsets[i] = (double)((gwRes[i]-mean));
	    }


	  // remove quadratic to make the total variation smaller.
	  TKremovePoly_d(epochs,offsets,psr[p].nobs,2);

	  printf("Write '%s'\n",fname);
	  toasim_write_corrections(corr,header,file);
	}
      fclose(file);
    }
  if (writeGW || readGW){
	  fclose(gwFile);
  }
  return 0;
}


long double getTspan(pulsar *psr,int npsr)
{
  long double first,last;
  int i,p;
    
  
  first = psr[0].obsn[0].sat;
  last = psr[0].obsn[0].sat;

  for (p=0;p<npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)
	{
	  if (first > psr[p].obsn[i].sat)
	    first = psr[p].obsn[i].sat;
	  if (last < psr[p].obsn[i].sat)
	    last = psr[p].obsn[i].sat;
	}
    }

  return last-first;
}
char * plugVersionCheck = TEMPO2_h_VER;
