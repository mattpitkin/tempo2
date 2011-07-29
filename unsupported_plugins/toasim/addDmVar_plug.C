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
#include "tempo2.h"
#include "toasim.h"

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

	// my parameters
	double alpha=-11.0/3.0;
	double cal_scale=1.0/100;
	double cal_val=1e-4;
	int resolution=1024;
	double is = 10; // days
	double os = -10; // 10 times data length


	//
	// For the output file
	//
	toasim_header_t* header;
	toasim_header_t* read_header;
	FILE* file;
	double offsets[MAX_OBSN]; // Will change to doubles - should use malloc
	double dms[MAX_OBSN]; // Will change to doubles - should use malloc
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

	printf("Graphical Interface: addDmVar\n");
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
		if (strcmp(argv[i],"-strength")==0){
			cal_val=atof(argv[++i]);
		}
		if (strcmp(argv[i],"-alpha")==0){
			alpha=atof(argv[++i]);
		}
		if (strcmp(argv[i],"-os")==0){
			os=atof(argv[++i]);
		}

	}

	alpha/=2.0; // we are computing amplitudes

	readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
	// Now read in all the .tim files
	readTimfile(psr,timFile,*npsr); /* Load the arrival times    */

	preProcess(psr,*npsr,argc,argv);
	formBatsAll(psr,*npsr);

	for (p=0;p<*npsr;p++)
	{
		printf("NTOA = %d\n",psr[p].nobs);
		header = toasim_init_header();
		strcpy(header->short_desc,"addDmVar");
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
		sprintf(fname,"%s.addDmVar",timFile[p]);
		file = toasim_write_header(header,fname);

		double mjd_start=1000000.0;
		double mjd_end=-10000000.0;
		for (j=0;j<psr[p].nobs;j++){
			// find the start and end times
			if(psr[p].obsn[j].bat < mjd_start)mjd_start=(double)psr[p].obsn[j].bat;
			if(psr[p].obsn[j].bat > mjd_end)mjd_end=(double)psr[p].obsn[j].bat;
		}
		if(os < 0)os=-os*(mjd_end-mjd_start);

		double fft_binsize = 1.0/(mjd_end-mjd_start);
		int nwav=(int)(os/is);
		double fstart=1.0/os; // per day
		double fend=1.0/is;
		double fstep=(fend-fstart)/double(nwav);

		double fft_sf = 1e-6*fstep/fft_binsize;


		printf("start    = %f (mjd)\n",mjd_start);
		printf("end      = %f (mjd)\n",mjd_end  );
		printf("OS       = %f (days)\n",os       );
		printf("IS       = %f (days)\n",is       );
		printf("nwav     = %d\n",nwav     );
		printf("fmin     = %f (1/day)\n",fstart);
		printf("fmax     = %f (1/day)\n",fend);
		printf("cal_scale= %f (1/day)\n",cal_scale);
		printf("cal_value= %f (?)\n",cal_val);
		printf("fft_sf   = %f (?)\n",fft_sf);

		

		for (i=0;i<nit;i++)
		{
			FILE *log_spec = fopen("dmvar.spec","w");
			if(i%10 == 0){
				printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
				printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
				printf("Iteration %d/%d",i+1,nit);
				fflush(stdout);
			}
			for (j=0;j<psr[p].nobs;j++){
				offsets[j] =0;
				dms[j] =0;
			}
			double a=pow(cal_scale,alpha);
			double scalefactor = sqrt(fft_sf*cal_val/(2*a*a));
			for (int iwav = 0; iwav < nwav; iwav++){
				double f=fstart+(double)iwav*fstep;
				a = scalefactor*pow(f,alpha);
				double a2 =a*TKgaussDev(&seed);
				a*=TKgaussDev(&seed);
				fprintf(log_spec,"%g %g\n",f,a*a+a2*a2);
				double phase=TKranDev(&seed) * 2*M_PI;
				for (j=0;j<psr[p].nobs;j++){
					double t=psr[p].obsn[j].bat - mjd_start;
					double dmv=a*sin(2*M_PI*t*f + phase) + a2*cos(2*M_PI*t*f + phase);
					dms[j]+=dmv;
				}
			}
			FILE *log_ts = fopen("dmvar.ts","w");
			double sum=0;
			for (j=0;j<psr[p].nobs;j++){
				sum+=dms[j];
			}
			sum/=psr[p].nobs;
			for (j=0;j<psr[p].nobs;j++){
				dms[j]-=sum;
				double ofreq=psr[p].obsn[j].freqSSB;
				offsets[j] = (double)(dms[j]/DM_CONST/ofreq/ofreq)*1e12;
				fprintf(log_ts,"%lg %lg %lg %lg\n",(double)psr[p].obsn[j].bat,dms[j],offsets[j],(double)ofreq);
			}
			toasim_write_corrections(corr,header,file);
			fclose(log_spec);
			fclose(log_ts);
		}
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("Iteration %d/%d\n",i,nit);

		printf("Close file\n");
		fclose(file);
	}
	return 0;
}

