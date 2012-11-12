//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards
//
// Plugin to read a set of arrival time files and produce a list of Gaussian random numbers based on the TOA uncertainties 


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
#include "TKfit.h"
#include "tempo2.h"
#include "toasim.h"
#include "makeRedNoise.h"

using namespace std;

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
   char parFile[MAX_PSR][MAX_FILELEN];
   char timFile[MAX_PSR][MAX_FILELEN];
   int i,nit,j,p;
   char fname[MAX_FILELEN];
   double globalParameter;
   long double result;
   long seed = TKsetSeed();

   double secperyear=365*86400.0;
   // my parameters
   double alpha= 1e99;
   int npts=1024;
   float p_1yr=-1; // s^2 yr
   char writeTextFiles=0;
   float cnr_flat=0;
   float cnr_cut=0;
   float old_fc=-1;


   //
   // For the output file
   //
   toasim_header_t* header;
   toasim_header_t* read_header;
   FILE* file;
   double offsets[MAX_OBSN]; // should use malloc
   double mjds[MAX_OBSN]; //  should use malloc
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

   printf("Graphical Interface: addRedNoise\n");
   printf("Author:              M. Keith\n");
   printf("Version:             1.0\n");

   /* Obtain all parameters from the command line */
   for (i=2;i<argc;i++)
   {

	  if (strcmp(argv[i],"-nreal")==0){
		 nit=atoi(argv[++i]);
	  }
	  if (strcmp(argv[i],"-f")==0)
	  {
		 strcpy(parFile[*npsr],argv[++i]); 
		 strcpy(timFile[*npsr],argv[++i]);
		 (*npsr)++;
	  }
	  if (strcmp(argv[i],"-npts")==0){
		 npts=atoi(argv[++i]);
	  }

	  if (strcmp(argv[i],"-Pus2yr")==0){
		 p_1yr=atof(argv[++i])*1e-12; // convert to s^2yr
	  }
	  if (strcmp(argv[i],"-Pyr3")==0){
		 p_1yr=atof(argv[++i]) * secperyear*secperyear; // convert to s^2yr
	  }

	  if (strcmp(argv[i],"-fc")==0){
		 old_fc=atof(argv[++i]);
		 cnr_flat=old_fc;
	  }

	  if (strcmp(argv[i],"-cnr_flat")==0){
		 cnr_flat=atof(argv[++i]);
	  }

	  if (strcmp(argv[i],"-cnr_cut")==0){
		 cnr_cut=atof(argv[++i]);
	  }

	  if (strcmp(argv[i],"-a")==0){
		 alpha=atof(argv[++i]);
		 if(alpha > 0){
			logmsg("Warning: alpha should normally be negative!!");
		 }
	  }

	  if (strcmp(argv[i],"-txt")==0){
		 writeTextFiles=1;
	  }


	  else if (strcmp(argv[i],"-seed")==0){
		 sscanf(argv[++i],"%d",&seed);
		 if (seed > 0)seed=-seed;
	  }

   }


   if(p_1yr < 0 || alpha>1e90){
	  printf("Invalid parameters\n");
	  exit(1);
   }


   readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
   // Now read in all the .tim files
   readTimfile(psr,timFile,*npsr); /* Load the arrival times    */

   preProcess(psr,*npsr,argc,argv);
   formBatsAll(psr,*npsr);

   for (p=0;p<*npsr;p++)
   {
	  printf("NTOA = %d\n",psr[p].nobs);
	  header = toasim_init_header();
	  strcpy(header->short_desc,"addRedNoise");
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
	  sprintf(fname,"%s.addRedNoise",timFile[p]);
	  file = toasim_write_header(header,fname);

	  double mjd_start=1000000.0;
	  double mjd_end=-10000000.0;
	  for (j=0;j<psr[p].nobs;j++){
		 // find the start and end times
		 if(psr[p].obsn[j].bat < mjd_start)mjd_start=(double)psr[p].obsn[j].bat;
		 if(psr[p].obsn[j].bat > mjd_end)mjd_end=(double)psr[p].obsn[j].bat;
	  }


	  printf("start    = %f (mjd)\n",mjd_start);
	  printf("end      = %f (mjd)\n",mjd_end  );
	  printf("npts     = %d (days)\n",npts     );
	  if(old_fc > 0){
		 printf("Mode     : T2Chol\n");
		 printf("P        = %g (yr^3)\n",p_1yr/secperyear/secperyear);
		 printf("P        = %g (us^2yr)\n",p_1yr*1e12);
		 printf("fc       = %f (yr^-1)\n",cnr_flat);
	  }else{
		 printf("Mode     : Simple\n");
		 printf("P(1yr)   = %g (yr^3)\n",p_1yr/secperyear/secperyear);
		 printf("P(1yr)   = %g (us^2yr)\n",p_1yr*1e12);
		 printf("Flatten  @ %f (yr^-1)\n",cnr_flat);
		 printf("Cutoff   @ %f (yr^-1)\n",cnr_cut);
	  }
	  printf("index    = %f\n",alpha);

	  printf("seed     = %d\n",seed);

	  printf("\n");
	  printf("Generating red noise...\n");

	  rednoisemodel_t* model = setupRedNoiseModel(mjd_start,mjd_end,npts,nit,p_1yr,alpha);

	  model->cutoff=cnr_cut;
	  model->flatten=cnr_flat;
	  if(old_fc>0)
		 model->mode=MODE_T2CHOL;

	  populateRedNoiseModel(model,seed);

	  if (writeTextFiles){
		 FILE *log_spec = fopen("red.spec","w");
		 float* pwr_spec=getPowerSpectrum(model);
		 float ps_fres=1.0/((model->end-model->start)/365.25); //yr^-1
		 for (j=0;j<(model->npt/2+1);j++){
			fprintf(log_spec,"%10.10g %10.10g\n",ps_fres*j,pwr_spec[j]);
		 }

		 fclose(log_spec);
		 FILE *log_ts2 = fopen("red.ts2","w");
		 for (j=0;j<(model->npt);j++){
			fprintf(log_ts2,"%10.10g %10.10g\n",model->start+model->tres*j,model->data[j]);
		 }

		 fclose(log_ts2);
	  }

	  int itjmp=nit/50;
	  if (itjmp<1)itjmp=1;
	  int dots=0;
	  printf("v");
	  for (i=0;i<nit/itjmp;i++){
		 printf("_");
	  }
	  printf("v\n");
	  printf("[");
	  fflush(stdout);
	  for (i=0;i<nit;i++)
	  {
		 if (i%itjmp==0){
			int v = i/itjmp;
			v-=dots;
			while (v > 0){
			   printf(".");
			   fflush(stdout);
			   v--;
			   dots++;
			}
		 }

		 for (j=0;j<psr[p].nobs;j++){
			offsets[j]=getRedNoiseValue(model,psr[p].obsn[j].bat,i);
		 }
		 FILE *log_ts;
		 if (writeTextFiles)
			log_ts = fopen("red.ts","w");
		 double sum=0;
		 for (j=0;j<psr[p].nobs;j++){
			sum+=offsets[j];
		 }
		 sum/=psr[p].nobs;
		 for (j=0;j<psr[p].nobs;j++){
			mjds[j]=(double)psr[p].obsn[j].bat;
			offsets[j]-=sum;
			if (writeTextFiles)
			   fprintf(log_ts,"%lg %lg %lg %lg\n",(double)psr[p].obsn[j].bat,offsets[j]);
		 }
		 TKremovePoly_d(mjds,offsets,psr[p].nobs,2); // remove a quadratic to reduce the chances of phase wraps
		 // The above is ok because it's linear with F0/F1
		 toasim_write_corrections(corr,header,file);
		 if (writeTextFiles)
			fclose(log_ts);
	  }
	  int v = i/itjmp;
	  v-=dots;
	  while (v > 0){
		 printf(".");
		 v--;
		 dots++;
	  }

	  printf("]\n");

	  printf("Close file\n");
	  fclose(file);
   }
   return 0;
}

char * plugVersionCheck = TEMPO2_h_VER;




