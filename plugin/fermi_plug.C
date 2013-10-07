#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include "../tempo2.h"
#include <cpgplot.h>
#include <fitsio.h>
#include <time.h>

#define SECDAY 86400.0

using namespace std;

void extra_delays_fermi(pulsar *psr,int npsr);
void clock_corrections_fermi(pulsar *psr,int npsr);
void ephemeris_routines_fermi(pulsar *psr,int npsr);
void formBatsAll_fermi(pulsar *psr,int npsr);

static char random_letter(int is_cap);
static char random_number();
static void random_string(int length, char *str);
float HTest(int Nphotons, float phases[]);
void cpgpt(int n, const float *xpts, const longdouble *ypts, int symbol);

// met2mjd: conversion of Mission Elapsed Time in TT to MJD
longdouble met2mjd(double met)
{
	longdouble mjd_ref = 51910.0007428703703703703;
	longdouble mjd = mjd_ref + met / 86400.;
	
	return mjd;
}

// mjd2met : conversion of MJD in TT to Mission Elapsed Time
double mjd2met(longdouble mjd)
{
	longdouble mjd_ref = 51910.0007428703703703703;
	double met = (double)(86400. * (mjd - mjd_ref));
	
	return met;
}

double inner_product(double vect_x[], double vect_y[])
{
  	return vect_x[0]*vect_y[0] + vect_x[1]*vect_y[1] + vect_x[2]*vect_y[2];
}

void outer_product(double vect_x[], double vect_y[], double vect_z[])
{
  	vect_z[0] = vect_x[1]*vect_y[2] - vect_x[2]*vect_y[1];
  	vect_z[1] = vect_x[2]*vect_y[0] - vect_x[0]*vect_y[2];
  	vect_z[2] = vect_x[0]*vect_y[1] - vect_x[1]*vect_y[0];
}

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
	// Provide basic information
	printf("\n");
	printf("------------------------------------------\n");
	printf("Output interface:    fermi\n");
	printf("Author:              Lucas Guillemot\n");
	printf("Updated:             4 October 2013\n");
	printf("Version:             5.10\n");
	printf("------------------------------------------\n");
	printf("\n");

	int i,j,k,l,event = 0,event2;

	char parFile[1][MAX_FILELEN];
	char timFile[1][MAX_FILELEN];
	char tempchar[9];
	char tempfits[14];
	char gr[256]="/xs";					// default graphical output is x-window

	srand((unsigned)time(NULL));
	tempchar[8] = '\0';
	random_string(8,tempchar);
	strcpy(timFile[0],tempchar);
	strcat(timFile[0],".tim");
	strcpy(tempfits,tempchar);
	strcat(tempfits,".fits");
	
	int nbins  = 25;					// default number of bins for the output phase histogram
	int Hbins  = 25;					// default number of bins for the H-test vs time plot		

	char FT1[MAX_FILELEN];
	char FT2[MAX_FILELEN];
	char output[MAX_FILELEN];
	char output_pos_file[MAX_FILELEN];
	char output_prof_file[MAX_FILELEN];
	char phasecol[32];
	strcpy(phasecol,"PULSE_PHASE");
    char barycol[32];
    strcpy(barycol,"BARY_TIME");
	char command[128];
	
    char error_buffer[128];
    char history[128+MAX_FILELEN];
	int  par_file        = 0;
	int  FT1_file	     = 0;
	int  FT2_file	     = 0;
    int  using_tt        = 1;
	
	int  ophase          = 0;
	int  graph	         = 1;
	int  output_file     = 0;
	int  output_pos      = 0;
	int  output_prof     = 0;
	int  phase_replace   = 0;
    int  phasecol_set    = 0;
    int  phasecol_exists = 0;
    int  bary_replace    = 0;
    int  barycol_set     = 0;
    int  barycol_exists  = 0;
    int  mint_set        = 0;
    int  maxt_set        = 0;


	double intpart;

	FILE *outputf;
	FILE *output_posf;
	FILE *output_proff;
	FILE *temp_tim;

	/* ------------------------------------------------- //
	// fits file definitions
	// ------------------------------------------------- */

	char colname[FLEN_VALUE];
	int open_status = 0, status = 0, status2 = 0, hdupos, ncols_FT1,ncols_FT2,anynul;
	long nrows_FT1,nrows_FT2;
	double nulval = 0.;

	int FT1_time_col,FT1_phase_col,FT1_bary_col,FT2_time_col1,FT2_time_col2,FT2_pos_col;
	int nrows2, nrows3, rows_status, rows_left;
	int max_rows = MAX_OBSN_VAL / 2;
	
	/* ------------------------------------------------- //
	// Time and satellite position definitions
	// ------------------------------------------------- */

	double temptime;
	double minFT1time = 999999999., maxFT1time = 0.;
	double minFT2time = 999999999., maxFT2time = 0.;
	
	double mint = 0., maxt = 999999999.;

	double time_MET_TT[max_rows];
	longdouble time_MJD_TT;
	double obs_earth[max_rows][3];

	double sctime1 = 0., sctime2 = 0.;
	double scposn1[3] = {0., 0., 0.};
	double scposn2[3] = {0., 0., 0.};
	double intposn[3] = {0., 0., 0.};

	double fract = 0.;
	double length1, length2, length12, intlength;
    double vector12[3], vectprod_out[3];
	double inttheta, factor_cos, factor_sin;

	longdouble lasttime, tzrmjd_bary;
	double lastpos[3];
	
	double tt0, tpb, fac;

	/* ------------------------------------------------- //
	// cpgplot definitions
	// ------------------------------------------------- */

	float fontsize = 0.75;
	int linewidth = 1;

	/* ------------------------------------------------- //
	// Additional definitions
	// ------------------------------------------------- */

	psr[0].correctTroposphere = 0;

	/* ------------------------------------------------- //
	// Command-line arguments
	// ------------------------------------------------- */

	for (i=2;i<argc;i++)
	{
		if (strcmp(argv[i],"-f")==0)
		{
			par_file = 1;
			strcpy(parFile[0],argv[i+1]); 
			i++;
		}
		else if (strcmp(argv[i],"-ft1")==0)
		{
			FT1_file = 1;
			strcpy(FT1,argv[i+1]);
			i++;
		}
		else if (strcmp(argv[i],"-ft2")==0)
		{
			FT2_file = 1;
			strcpy(FT2,argv[i+1]);
			i++;
		}
		else if (strcmp(argv[i],"-graph")==0)
		{
			sscanf(argv[i+1],"%d",&graph);
			i++;
		}
      	else if (strcmp(argv[i],"-grdev")==0)
		{
			strcpy(gr,argv[i+1]);
			i++;
		}
		else if (strcmp(argv[i],"-bins")==0)
		{
			sscanf(argv[i+1],"%d",&nbins);
			i++;
		}
		else if (strcmp(argv[i],"-Hbins")==0)
		{
			sscanf(argv[i+1],"%d",&Hbins);
			i++;
		}
		else if (strcmp(argv[i],"-output")==0)
		{
			output_file = 1;
			strcpy(output,argv[i+1]);
			i++;
		}
		else if (strcmp(argv[i],"-pos")==0)
		{
			output_pos = 1;
			strcpy(output_pos_file,argv[i+1]);
			i++;
		}
		else if (strcmp(argv[i],"-profile")==0)
		{
			output_prof = 1;
			strcpy(output_prof_file,argv[i+1]);
			i++;
		}
		else if (strcmp(argv[i],"-phase")==0)
		{
			phase_replace = 1;
		}
		else if (strcmp(argv[i],"-ophase")==0)
		{
			ophase = 1;
		}
		else if (strcmp(argv[i],"-col")==0)
		{
			strcpy(phasecol,argv[i+1]);
			phasecol_set = 1;
			i++;
		}
		else if (strcmp(argv[i],"-tmin")==0)
		{
			sscanf(argv[i+1],"%lf",&mint);
			mint_set = 1;
			i++;
		}
		else if (strcmp(argv[i],"-tmax")==0)
		{
			sscanf(argv[i+1],"%lf",&maxt);
			maxt_set = 1;
			i++;
		}
		else if (strcmp(argv[i],"-bary")==0)
		{
			bary_replace = 1;
		}
		else if (strcmp(argv[i],"-barycol")==0)
		{
			strcpy(barycol,argv[i+1]);
			barycol_set = 1;
			i++;
		}
		else if (strcmp(argv[i],"-h")==0)
		{
			printf("\n TEMPO2 fermi plugin\n");
			printf("======================\n");
			printf("\n USAGE: \n\t tempo2 -gr fermi -ft1 FT1.fits -ft2 FT2.fits -f par.par\n");
			printf("\n Command line options:\n");
			printf("\t -grdev XXX/YYY, where XXX/YYY is the graphical device of your choice (e.g. a.ps/cps). Default mode is XW\n");
			printf("\t -graph 0: no output graph is drawn\n");
			printf("\t -bins N: number of bins of the phase histogram\n");
			printf("\t -Hbins N: number of time bins for the H-test vs time plot\n");
			printf("\t -output XXX: writes event times and phases in file XXX\n");
			printf("\t -pos XXX: puts the satellite (X,Y,Z) positions as a function of time in file XXX\n");
			printf("\t -profile XXX: stores the binned light curve in file XXX (plotting mode only)\n");
			printf("\t -phase: stores phases in the FT1 by the ones calculated by TEMPO2\n");
			printf("\t -ophase: will calculate orbital phases instead of pulse phases. Changes default column name to ORBITAL_PHASE.\n");
			printf("\t -col XXX: phases will be stored in column XXX. Default is PULSE_PHASE\n");
			printf("\t -tmin XXX: event times smaller than MJD XXX will be skipped.\n");
			printf("\t -tmax XXX: event times larger than MJD XXX will be skipped.\n");
			printf("\t -bary: write out barycenter arrival times.\n");
			printf("\t -barycol XXX: bats will be stored in column XXX. Default is BARY_TIME.\n");
			printf("\t -h: this help.\n");
			printf("===============================================");
			printf("===============================================\n");
			exit(0);			
		}
    }
    
	if (FT1_file == 0)
	{
		printf("No input FT1 file !\n");
		printf("At least it should look like :\n");
		printf("\t tempo2 -gr fermi -ft1 FT1.fits -ft2 FT2.fits -f par.par\n");
		printf("More options are available. If you need help, please type:\n");
		printf("\t tempo2 -gr fermi -h\n");
		exit(1);
	}

	if (FT2_file == 0)
	{
		printf("No spacecraft file !\n");
		printf("At least it should look like :\n");
		printf("\t tempo2 -gr fermi -ft1 FT1.fits -ft2 FT2.fits -f par.par\n");
		printf("More options are available. If you need help, please type:\n");
		printf("\t tempo2 -gr fermi -h\n");		
		exit(2);
	}

	if (par_file == 0)
	{
		printf("No ephemeris file !\n");
		printf("At least it should look like :\n");
		printf("\t tempo2 -gr fermi -ft1 FT1.fits -ft2 FT2.fits -f par.par\n");
		printf("More options are available. If you need help, please type:\n");
		printf("\t tempo2 -gr fermi -h\n");
		exit(3);
	}

	if ((graph == 0) && (output_prof))
	{
		printf("Can only save the profile in plotting mode (graph set to 1)!\n");
		printf("If you need help, please type:\n");
		printf("\t tempo2 -gr fermi -h\n");
		exit(4);
	}

	if (output_file)
	{
		outputf = fopen(output,"w+");
	}

	if (output_pos)
	{
		output_posf = fopen(output_pos_file,"w+");
	}

	if (ophase)
	{
		if (phasecol_set == 0) strcpy(phasecol,"ORBITAL_PHASE");
	}

	if (mint_set)
	{
		mint = mjd2met(mint);
	}

	if (maxt_set)
	{
		maxt = mjd2met(maxt);
	}

	/* ------------------------------------------------- //
	// FT1 file
	// ------------------------------------------------- */

    fitsfile *ft1;
	
  	if (!fits_open_file(&ft1,FT1, READWRITE, &open_status))
    {
        // Check that the input file has not been barycentered.
        int kw_status;
        char value[80];
        char comment[80];
        
        kw_status = 0;
        fits_read_key(ft1, TSTRING, "TIMESYS",(void*)value, comment, &kw_status);
        
        if (kw_status <= 0)
        {
            if (strcmp(value,"TDB")==0) {
                fprintf(stdout,"FT1 file has already been barycentered; will disable clock and delay corrections in phase computation.\n");
                using_tt = 0;
            }

            else if (strcmp(value,"TT"))
            {
                fprintf(stderr,"Error: TIMESYS is %s rather than TT or TDB; unable to proceed.\n", value);
                exit(2);
            }
        }
        else
        {
            fprintf(stderr, "Warning: unable to determine whether data are in a compatible time format; assuming TT. (Plugin will give wrong results if this is not the case!).\n");
        }
        
        // Load photon arrival times.
		fits_movabs_hdu(ft1,2,NULL,&status);

		fits_get_num_rows(ft1, &nrows_FT1, &status);
		fits_get_num_cols(ft1, &ncols_FT1, &status);

		fits_get_colname(ft1,CASESEN,(char *)"TIME",colname,&FT1_time_col,&status);
		
		//
		
		for (i=1;i<=nrows_FT1;i++)
		{
			fits_read_col_dbl(ft1,FT1_time_col,i,1,1,nulval,&temptime,&anynul,&status);
			if (temptime < minFT1time) minFT1time = temptime;
			if (temptime > maxFT1time) maxFT1time = temptime;
		}
		
		//

		if (phase_replace)
		{
		    fits_get_colname(ft1,CASESEN,phasecol,colname,&FT1_phase_col,&status);
            
            if (status > 0)
            {
                phasecol_exists = 0;
            }
            else
            {
                phasecol_exists = 1;
            }
    	}

		if (bary_replace)
		{
            status = 0; // reset status in case added phase column above
		    fits_get_colname(ft1,CASESEN,barycol,colname,&FT1_bary_col,&status);
            barycol_exists = (status == 0);
            
    	}
	}
	
	if (open_status != 0)
	{
		printf("Can't find %s !\n",FT1);
		exit(-1);
	}

    fits_close_file(ft1, &status);

    printf("First photon date in FT1: %lf MET (s)\n",minFT1time);
	printf(" Last photon date in FT1: %lf MET (s)\n\n",maxFT1time);
	
	
	
	
	// This part was adapted from A. Archibald's photons plugin
	fitsfile *ft_out;
    
	if (phase_replace)
    {
        char FT_in_filter[MAX_FILELEN+20];
        status = 0; 
        int ii = 1;
        
        if (phasecol_exists)
        {
            printf("Replacing existing column %s.\n\n", phasecol);
            strcpy(FT_in_filter, FT1);
        } 
        else 
        {
            printf("Adding new column %s.\n\n", phasecol);
            sprintf(FT_in_filter, "%s[%d][col *;%s = -1.0]", FT1, 1, phasecol);
        }
        
        if (!fits_open_file(&ft1, FT_in_filter, READONLY, &open_status))
        {
            // Record the column number for later
            if (fits_movabs_hdu(ft1,2,NULL,&status) || fits_get_colname(ft1,CASEINSEN,phasecol,colname,&FT1_phase_col,&status))
            {
                fprintf(stderr, "Unable to find column %s even after trying to create it.\n", phasecol);
                exit(-1);
            }
            
            status = 0;
            
            if (!fits_create_file(&ft_out, tempfits,  &status))
            {
                // Copy all HDUs one-by-one. In the process, make room for more header comments
                status = 0;
                i      = 1;
                
                while (!fits_movabs_hdu(ft1, i++, NULL, &status)) fits_copy_hdu(ft1, ft_out, 5, &status);
        
                if (status==END_OF_FILE) status=0;
                else printf("Error happened while transferring\n");
        
                fits_close_file(ft_out, &status);
            }
            
            fits_close_file(ft1, &status);
        }
        
        if (status != 0)
        {
            fits_report_error(stderr, status);
            exit(-1);
        }
        
        if (rename(tempfits,FT1))
        {
            if (errno != ENOENT)
            {
                perror("Problem renaming temporary fits file");
                exit(-1);
            }
        }
    }

	if (bary_replace)
    {
        char FT_in_filter[MAX_FILELEN+20];
        status = 0; 
        int ii = 1;
        
        if (barycol_exists)
        {
            printf("Replacing existing column %s.\n\n", barycol);
            strcpy(FT_in_filter, FT1);
        } 
        else 
        {
            printf("Adding new column %s.\n\n", barycol);
            sprintf(FT_in_filter, "%s[%d][col *;%s = -1.0]", FT1, 1, barycol);
        }
        
        if (!fits_open_file(&ft1, FT_in_filter, READONLY, &open_status))
        {
            // Record the column number for later
            if (fits_movabs_hdu(ft1,2,NULL,&status) || fits_get_colname(ft1,CASEINSEN,barycol,colname,&FT1_bary_col,&status))
            {
                fprintf(stderr, "Unable to find column %s even after trying to create it.\n", barycol);
                exit(-1);
            }
            
            status = 0;
            
            if (!fits_create_file(&ft_out, tempfits,  &status))
            {
                // Copy all HDUs one-by-one. In the process, make room for more header comments
                status = 0;
                i      = 1;
                
                while (!fits_movabs_hdu(ft1, i++, NULL, &status)) fits_copy_hdu(ft1, ft_out, 5, &status);
        
                if (status==END_OF_FILE) status=0;
                else printf("Error happened while transferring\n");
        
                fits_close_file(ft_out, &status);
            }
            
            fits_close_file(ft1, &status);
        }
        
        if (status != 0)
        {
            fits_report_error(stderr, status);
            exit(-1);
        }
        
        if (rename(tempfits,FT1))
        {
            if (errno != ENOENT)
            {
                perror("Problem renaming temporary fits file");
                exit(-1);
            }
        }
    }
    
    
	/* ------------------------------------------------- //
	// FT2 file
	// ------------------------------------------------- */

    fitsfile *ft2;
    open_status = 0, status = 0;

	if (!fits_open_file(&ft2,FT2, READONLY, &open_status))
	{
		fits_movabs_hdu(ft2,2,NULL,&status);

		fits_get_num_rows(ft2, &nrows_FT2, &status);
		fits_get_num_cols(ft2, &ncols_FT2, &status);

		fits_get_colname(ft2,CASESEN,(char *)"START",colname,&FT2_time_col1,&status);
		fits_get_colname(ft2,CASESEN,(char *)"STOP",colname,&FT2_time_col2,&status);
		fits_get_colname(ft2,CASESEN,(char *)"SC_POSITION",colname,&FT2_pos_col,&status);
		
		for (i=1;i<=nrows_FT2;i++)
		{
			fits_read_col_dbl(ft2,FT2_time_col1,i,1,1,nulval,&temptime,&anynul,&status);
			if (temptime < minFT2time) minFT2time = temptime;
			if (temptime > maxFT2time) maxFT2time = temptime;
		}
	}

	if (open_status != 0)
	{
		printf("Can't find %s !\n",FT2);
		exit(-1);
	}

	fits_close_file(ft2, &status);
	
	printf("First START date in FT2:  %lf MET (s)\n",minFT2time);
	printf(" Last START date in FT2:  %lf MET (s)\n\n",maxFT2time);
	
	if (minFT1time < minFT2time)
	{
		printf("Warning: no spacecraft data for date %lf MET. Check the FT2 file.\n",minFT1time);
		printf("All photon phases before %lf MET will be set to -1.\n\n",minFT2time);
	}
	
	if (maxFT1time > maxFT2time)
	{
		printf("Warning: no spacecraft data for date %lf MET. Check the FT2 file.\n",maxFT1time);
		printf("All photon phases after %lf MET will be set to -1.\n\n",maxFT2time);
	}

	if (mint_set)
	{
		printf("Will skip event times smaller than %lf MET.\n\n",mint);

		if (mint > maxFT1time)
		{
			printf("Error: requested tmin larger than all event times in the FT1 file! Exitting...\n");
			exit(5);
		}
	}

	if (maxt_set)
	{
		printf("Will skip event times larger than %lf MET.\n\n",maxt);

		if (maxt < minFT1time)
		{
			printf("Error: requested tmax smaller than all event times in the FT1 file! Exitting...\n");
			exit(6);
		}
	}
	

	// The FT1 file is cut into small tim files with # rows <= (MAX_OBSN_VAL / 2)
	rows_status = 1;
	rows_left	= nrows_FT1 - rows_status + 1;
	nrows2 		= min(rows_left,max_rows);

	float *phase;
	longdouble *times;
    int arr_size = max_rows;

    if (graph == 1)
    {
        if (nrows_FT1 > 300000)         // To avoid core dumps
        {
            printf("WARNING: large FT1 file, turning off graphical output.\n");
            graph = 0;
        }
        else
        {
            arr_size = nrows_FT1;
        }
    }
         
    phase  = (float *)calloc(arr_size,sizeof(float));
    times  = (longdouble *)calloc(arr_size,sizeof(longdouble));
	float tmin   = 100000., tmax   = -100000.;
	
	
	/* ------------------------------------------------- //
	// Barycentric TZRMJD
	// ------------------------------------------------- */
	
	// A temporary file is created. It is first used to get the barycentered TZR
	temp_tim = fopen(timFile[0],"w+");
	fprintf(temp_tim,"FORMAT 1\n");
	fprintf(temp_tim," fermi 0.0 %.12lf 0.00000 BAT\n",1.);
	fclose(temp_tim);

	// Load the arrival times
	readParfile(psr,parFile,timFile,*npsr);
	psr[0].nobs = 0;
	readTimfile(psr,timFile,*npsr);
	
	if (ophase && (strcmp(psr[0].binaryModel,"NONE")==0))
	{
		printf("Error: no binary parameters found in %s !\n",parFile[0]);
		sprintf(command,"rm -f %s",timFile[0]);
		system(command);
		exit(-1);
	}
	
	if (psr->tzrsite[0] == '@')
	{
		tzrmjd_bary = psr->param[param_tzrmjd].val[0];
	}
	else
	{
		psr->obsn[0].sat = psr->param[param_tzrmjd].val[0];
		psr->obsn[0].freq = psr->param[param_tzrfrq].val[0];
		strcpy(psr->obsn[0].telID, psr->tzrsite); 
		psr->obsn[0].deleted = 0;
		psr->obsn[0].nFlags = 0;
		psr->obsn[0].delayCorr = 1;
		psr->obsn[0].clockCorr = 1;
		
		// Form barycentric TZR
		formBatsAll(psr,*npsr);
		tzrmjd_bary = psr->obsn[0].bat;
	}
	
	

	/* ------------------------------------------------- //
	// Main loop
	// ------------------------------------------------- */	

	while (rows_left > 0)
	{
		printf("Treating events # %d to %d... \n",rows_status,rows_status + nrows2 - 1);

		// Acquisition of the non-barycentered TOAs in MET TT
		nrows3 = nrows2;
		j      = 0;

		if (!fits_open_file(&ft1,FT1, READONLY, &status))
		{
			fits_movabs_hdu(ft1,2,NULL,&status);

			for (i=rows_status;i<rows_status + nrows2;i++)
			{
				fits_read_col_dbl(ft1,FT1_time_col,i,1,1,nulval,&time_MET_TT[j],&anynul,&status);

				if (mint_set && (time_MET_TT[j] < mint)) nrows3--;
				if (maxt_set && (time_MET_TT[j] > maxt)) nrows3--;

				j++;
			}
		}

		fits_close_file(ft1, &status);

		// If all event times are outside the requested range, continue to the next set of photons TBD
		if (nrows3 == 0)
		{
			printf("All events are outside the requested time range. Continuing...\n");

			if (graph == 0)
			{
				event = 0;
			}
			else
			{
				for (i=0;i<nrows2;i++)
				{
					phase[event] = -1.;
					event++;
				}
			}

			rows_status	+= nrows2;
			rows_left	= nrows_FT1 - rows_status + 1;
			nrows2 		= min(rows_left,max_rows);
			continue;
		}

		// A temporary file is created. TOAs in MJD TT are stored is this file
		temp_tim = fopen(timFile[0],"w+");
		fprintf(temp_tim,"FORMAT 1\n");

		for (i=0;i<nrows2;i++)
		{
			time_MJD_TT = met2mjd(time_MET_TT[i]);		
			fprintf(temp_tim," fermi 0.0 %.12Lf 0.00000 BAT\n",time_MJD_TT);
		}

		fclose(temp_tim);

		// Load the arrival times
		readParfile(psr,parFile,timFile,*npsr); 
		psr[0].nobs = 0;
		readTimfile(psr,timFile,*npsr);
		
		/* ------------------------------------------------- //
		// Satellite position as a function of time
		// ------------------------------------------------- */

		j = 1;

		for (i=0;i<nrows2;i++)
		{
			if ((time_MET_TT[i] < minFT2time) || (time_MET_TT[i] > maxFT2time))
			{
				for (k=0;k<3;k++)
				{
					obs_earth[i][k] = 0.;
				}
			}
			else if (mint_set && (time_MET_TT[i] < mint))
			{
				for (k=0;k<3;k++)
				{
					obs_earth[i][k] = 0.;
				}
			}
			else if (maxt_set && (time_MET_TT[i] > maxt))
			{
				for (k=0;k<3;k++)
				{
					obs_earth[i][k] = 0.;
				}
			}
			else
			{
				if (!fits_open_file(&ft2,FT2, READONLY, &status))
				{
					fits_movabs_hdu(ft2,2,NULL,&status);
	
					if (time_MET_TT[i] < sctime1) j = 1;		// Important, if times in FT1 are not monotonically increasing
	
					for (;j<=nrows_FT2;j++)
					{
						fits_read_col_dbl(ft2,FT2_time_col1,j,1,1,nulval,&sctime2,&anynul,&status);
	
						if (time_MET_TT[i] < sctime2) break;
					}
	
					fits_read_col_dbl(ft2,FT2_time_col1,j-1,1,1,nulval,&sctime1,&anynul,&status);
	
					fits_read_col_dbl(ft2,FT2_pos_col,j-1,1,3,nulval,scposn1,&anynul,&status);
					fits_read_col_dbl(ft2,FT2_pos_col,j,1,3,nulval,scposn2,&anynul,&status);
				}
	
				fits_close_file(ft2, &status);
				
				// Interpolation
				fract = (time_MET_TT[i] - sctime1) / (sctime2 - sctime1);
	
				// linear interpolation for vector length
				length1 = sqrt(inner_product(scposn1, scposn1));
				length2 = sqrt(inner_product(scposn2, scposn2));
				intlength = length1 + fract*(length2 - length1);
	
				// Compute a base vector on the orbital plane (vector12)
				outer_product(scposn1, scposn2, vectprod_out);
				outer_product(vectprod_out, scposn1, vector12);
				length12 = sqrt(inner_product(vector12, vector12));
		    
				// Check vectors scposn1 and scposn2
				if ((length1 == 0.0) && (length2 == 0.0))
				{
					// both vectors are null 
					for (k=0;k<3;k++)
					{
						intposn[k] = 0.0;
					}
				} 
				else if (length1 == 0.0)
				{
					// scposn1 is null, but scposn2 is not
					for (k=0;k<3;k++)
					{
						intposn[k] = scposn2[k] / length2 * intlength;
					}
				} 
				else if ((length2 == 0.0) || (length12 == 0.0))
				{
					// left:  scposn2 is null, but scposn1 is not 
					// right: either vector is not null, but they are parallel
					for (k=0;k<3; k++)
					{
						intposn[k] = scposn1[k] / length1 * intlength;
					}
				} 
				else
				{ 
					// both has a non-zero length, and they are not parallel
					
					// linear interpolation for orbital phase 
					inttheta = fract * acos(inner_product(scposn1, scposn2) / length1 / length2);
					factor_cos = cos(inttheta);
					factor_sin = sin(inttheta);
	
					for (k=0;k<3;k++)
					{
						intposn[k] = intlength * (scposn1[k] / length1 * factor_cos + vector12[k] / length12 * factor_sin);
					}
				}
	
				for (k=0;k<3;k++)
				{
					obs_earth[i][k] = intposn[k] / SPEED_LIGHT;
				}
				
				// Output TT time and observatory position
				if (output_pos == 1)
				{
					fprintf(output_posf,"%.6lf\t",time_MET_TT[i]);
					
					for (k=0;k<3;k++)
					{
						fprintf(output_posf,"%12.3f ",intposn[k]);
					}
					
					fprintf(output_posf,"\n");
				}
			}
		}


		/* ------------------------------------------------- //
		// Delays, corrections and positions
		// ------------------------------------------------- */

		for (i=0;i<nrows2;i++)
		{
			psr->obsn[i].deleted = 0;
			psr->obsn[i].nFlags = 0;
			psr->obsn[i].delayCorr = using_tt;	// Correct delays for event i
			psr->obsn[i].clockCorr = using_tt;	// Also make clock correction TT -> TDB

			// Position replacements
			for (k=0;k<3;k++) 
			{
				psr[0].obsn[i].observatory_earth[k] = obs_earth[i][k];
			}
		}


		/* ------------------------------------------------- //
		// Calculation of the event phases - step 1
		// ------------------------------------------------- */

		// keep track of the last TOA before shifting the others
		lasttime = psr->obsn[psr[0].nobs-1].sat;
		
		for (k=0;k<3;k++) lastpos[k] = psr->obsn[psr[0].nobs-1].observatory_earth[k];
		
		for (i=psr[0].nobs-1;i>0;i--)
		{
			psr->obsn[i].sat = psr->obsn[i-1].sat;
			
			for (k=0;k<3;k++) psr[0].obsn[i].observatory_earth[k] = psr[0].obsn[i-1].observatory_earth[k];
		}
		
		/* ------------------------------------------------- //
		// Stick in a fake obs to get the reference phase
		// ------------------------------------------------- */

		psr->obsn[0].sat = tzrmjd_bary;
		psr->obsn[0].freq = 0.;
		psr->obsn[0].deleted = 0;
		psr->obsn[0].nFlags = 0;
		psr->obsn[0].delayCorr = 0;
		psr->obsn[0].clockCorr = 0;
		
		// ------------------------------------------------- //
		// Form barycentric arrival times - step 1
		// ------------------------------------------------- //

		formBatsAll_fermi(psr,*npsr);

		// ------------------------------------------------- //
		// Calculate event phases - step 1
		// ------------------------------------------------- //		

		formResiduals(psr,*npsr,0.0);

		for (i=1;i<nrows2;i++)
		{
			if ((time_MET_TT[i-1] < minFT2time) || (time_MET_TT[i-1] > maxFT2time))
			{
				phase[event] = -1.;
			}
			else if (mint_set && (time_MET_TT[i-1] < mint))
			{
				phase[event] = -1.;
			}
			else if (maxt_set && (time_MET_TT[i-1] > maxt))
			{
				phase[event] = -1.;
			}
			else
			{
				if (ophase == 0)
				{
					phase[event] = modf(psr[0].obsn[i].phase,&intpart);			
					if (phase[event] < 0.) phase[event]++;
				}
				else
				{
					if (strcmp(psr[0].binaryModel,"BTX")==0)
					{
						tt0 = (psr[0].obsn[i].bat-psr[0].param[param_t0].val[0])*SECDAY;
						tpb = psr[0].param[param_fb].val[0]*tt0;
						fac = 1.0;
						
						for (j=1;j<psr[0].param[param_fb].aSize;j++)
						{
							fac = fac/(double)(j+1);
							if (psr[0].param[param_fb].paramSet[j]==1) tpb += fac*psr[0].param[param_fb].val[j]*powl(tt0,j+1);
						}
					}
					else if (psr[0].param[param_pb].paramSet[0] && psr[0].param[param_tasc].paramSet[0])
					{
						tpb = (psr[0].obsn[i].bat-psr[0].param[param_tasc].val[0])/(psr[0].param[param_pb].val[0]);
						if (psr[0].param[param_pbdot].paramSet[0]) tpb -= 0.5*psr[0].param[param_pbdot].val[0]*pow((psr[0].obsn[i].bat-psr[0].param[param_tasc].val[0])/(psr[0].param[param_pb].val[0]),2.);
					}
					else if (psr[0].param[param_pb].paramSet[0] && psr[0].param[param_t0].paramSet[0])
					{
						tpb = (psr[0].obsn[i].bat-psr[0].param[param_t0].val[0])/(psr[0].param[param_pb].val[0]);
						if (psr[0].param[param_pbdot].paramSet[0]) tpb -= 0.5*psr[0].param[param_pbdot].val[0]*pow((psr[0].obsn[i].bat-psr[0].param[param_t0].val[0])/(psr[0].param[param_pb].val[0]),2.);
					}
					else
					{
						printf("Warning: PB, or T0, or TASC not set.\n");
					}
					
					phase[event] = modf(tpb+1000000.0,&intpart);
					if (phase[event] < 0.) phase[event]++;
				}
	
				times[event] = psr[0].obsn[i].bat;
				if (times[event] > tmax) tmax = times[event];
				if (times[event] < tmin) tmin = times[event];
			}
				
			if (output_file)
			{
				if (graph == 0)	fprintf(outputf,"%d\t",event + rows_status);
				else fprintf(outputf,"%d\t",event + 1);
			
				fprintf(outputf,"%20.12Lf\t%20.12Lf\t%12.10le\t%.0lf\n",psr[0].obsn[i].sat,psr[0].obsn[i].bat,phase[event],intpart);
			}
	
			event++;
		}		

		// ------------------------------------------------- //
		// Calculation of the event phases - step 2
		// ------------------------------------------------- //

		psr[0].obsn[1].sat = lasttime;
		for (k=0;k<3;k++) psr[0].obsn[1].observatory_earth[k] = lastpos[k];
	
		// ------------------------------------------------- //
		// Form barycentric arrival times - step 2
		// ------------------------------------------------- //

		formBatsAll_fermi(psr,*npsr);

		// ------------------------------------------------- //
		// Calculate event phases - step 2
		// ------------------------------------------------- //		

		formResiduals(psr,*npsr,0.0);
		
		if ((time_MET_TT[psr[0].nobs-1] < minFT2time) || (time_MET_TT[psr[0].nobs-1] > maxFT2time))
		{
			phase[event] = -1.;
		}
		else if (mint_set && (time_MET_TT[psr[0].nobs-1] < mint))
		{
			phase[event] = -1.;
		}
		else if (maxt_set && (time_MET_TT[psr[0].nobs-1] > maxt))
		{
			phase[event] = -1.;
		}
		else
		{
			if (ophase == 0)
			{
				phase[event] = modf(psr[0].obsn[1].phase,&intpart);
				if (phase[event] < 0.) phase[event]++;
			}
			else
			{
				if (strcmp(psr[0].binaryModel,"BTX")==0)
				{
					tt0 = (psr[0].obsn[1].bat-psr[0].param[param_t0].val[0])*SECDAY;
					tpb = psr[0].param[param_fb].val[0]*tt0;
					fac = 1.0;
					
					for (j=1;j<psr[0].param[param_fb].aSize;j++)
					{
						fac = fac/(double)(j+1);
						if (psr[0].param[param_fb].paramSet[j]==1) tpb += fac*psr[0].param[param_fb].val[j]*powl(tt0,j+1);
					}
				}
				else if (psr[0].param[param_pb].paramSet[0] && psr[0].param[param_tasc].paramSet[0])
				{
					tpb = (psr[0].obsn[1].bat-psr[0].param[param_tasc].val[0])/(psr[0].param[param_pb].val[0]);
					if (psr[0].param[param_pbdot].paramSet[0]) tpb -= 0.5*psr[0].param[param_pbdot].val[0]*pow((psr[0].obsn[1].bat-psr[0].param[param_tasc].val[0])/(psr[0].param[param_pb].val[0]),2.);
				}
				else if (psr[0].param[param_pb].paramSet[0] && psr[0].param[param_t0].paramSet[0])
				{
					tpb = (psr[0].obsn[1].bat-psr[0].param[param_t0].val[0])/(psr[0].param[param_pb].val[0]); 
					if (psr[0].param[param_pbdot].paramSet[0]) tpb -= 0.5*psr[0].param[param_pbdot].val[0]*pow((psr[0].obsn[i].bat-psr[0].param[param_t0].val[0])/(psr[0].param[param_pb].val[0]),2.);
				}
				else
				{
					printf("Warning: PB, or T0, or TASC not set.\n");
				}
				
				phase[event] = modf(tpb+1000000.0,&intpart);
				if (phase[event] < 0.) phase[event]++;
			}
		
			times[event] = psr[0].obsn[1].bat;		
			if (times[event] > tmax) 	tmax = times[event];
			if (times[event] < tmin)	tmin = times[event];
		}

		if (output_file)
		{
			if (graph == 0) fprintf(outputf,"%d\t",event + rows_status);
			else fprintf(outputf,"%d\t",event + 1);
			
			fprintf(outputf,"%20.12Lf\t%20.12Lf\t%12.10le\t%.0lf\n",psr[0].obsn[1].sat,psr[0].obsn[1].bat,phase[event],intpart);
		}

		// ------------------------------------------------- //
		// Put new phases into the input FT1 file
		// ------------------------------------------------- //

		if (phase_replace)
		{
  			if (!fits_open_file(&ft1,FT1, READWRITE, &open_status))
    		{
				fits_movabs_hdu(ft1,2,NULL,&status);
				
				for (event2=rows_status;event2<rows_status+nrows2;event2++)
				{
					if (graph == 0) i = event2 - rows_status;
					else i = event2 - 1;

					if (phase[i] != -1.) fits_write_col_flt(ft1,FT1_phase_col,event2,1,1,&phase[i],&status);
				}
			}

			fits_close_file(ft1, &status);		
		}

		if (bary_replace)
		{
  			if (!fits_open_file(&ft1,FT1, READWRITE, &open_status))
    		{
				fits_movabs_hdu(ft1,2,NULL,&status);
				
				for (event2=rows_status;event2<rows_status+nrows2;event2++)
				{
					if (graph == 0) i = event2 - rows_status;
					else i = event2 - 1;

                    // write out -1 for out-of-bounds events
                    double met_TDB = phase[i]<0?-1:mjd2met(times[i]);
					fits_write_col_dbl(ft1,FT1_bary_col,event2,1,1,&met_TDB,&status);
				}
			}

			fits_close_file(ft1, &status);		
		}

		/* ------------------------------------------------- //
		// End of the loop
		// ------------------------------------------------- */	

		rows_status	+= nrows2;
		rows_left	= nrows_FT1 - rows_status + 1;
		nrows2 		= min(rows_left,max_rows);
		
		if (graph == 0)
		{
			event = 0;
		}
		else
		{
			event++;
		}
	}
        
    // ------------------------------------------------- //
	// Add a little bit of history to the header
	// ------------------------------------------------- //

	if (phase_replace)
	{
		sprintf(history,"Pulse phases calculated with the TEMPO2 Fermi plugin using ephemeris %s",parFile[0]);

		if (!fits_open_file(&ft1,FT1, READWRITE, &open_status))
		{
			status = 0;
			fits_movabs_hdu(ft1,2,NULL,&status);
			fits_write_history(ft1,history,&status);

			if (status != 0)
			{
				fits_get_errstatus(status,error_buffer);
				printf( "fits_insert_col: %s\n", error_buffer);
				exit(-1);
			}
		}

		fits_close_file(ft1, &status);
	}

    if (bary_replace) 
    {
		sprintf(history,"Barycenter arrival times calculated with the TEMPO2 Fermi plugin using ephemeris %s",parFile[0]);
		if (!fits_open_file(&ft1,FT1, READWRITE, &open_status))
        {
            // change TIMESYS and TIMEREF if we have clobbered the TIME column
            if (strcasecmp(barycol,"TIME")==0)
            {
                char value[80];
                char comment[80];
                for(int hdr_ct=1;hdr_ct<3;hdr_ct++) {
                    fits_movabs_hdu(ft1,hdr_ct,NULL,&status);
                    status = 0;
                    strcpy(value,"TDB");
                    strcpy(comment,"type of time system that is used");
                    fits_update_key(ft1,TSTRING,"TIMESYS",value,comment,&status);
                    strcpy(value,"SOLARSYSTEM");
                    strcpy(comment,"reference frame used for times");
                    fits_update_key(ft1,TSTRING,"TIMEREF",value,comment,&status);
                }
            }

            // now write history
			fits_movabs_hdu(ft1,2,NULL,&status);
			fits_write_history(ft1,history,&status);

			if (status != 0)
			{
				fits_get_errstatus(status,error_buffer);
				printf( "fits_insert_col: %s\n", error_buffer);
				exit(-1);
			}
        }

        fits_close_file(ft1, &status);
    }
        

	if (output_file) fclose(outputf);
	if (output_pos)  fclose(output_posf);

	// ------------------------------------------------- //
	// Graphical output
	// ------------------------------------------------- //	

	if (graph)
	{	
		// ------------------------------------------------- //
		// 1D histogram of phase
		// ------------------------------------------------- //							

		float bins[nbins*2];
		float lc[nbins*2];
		float errX[nbins*2],errY[nbins*2];
		float binmax = -1000;
		
		for (i=0;i<nbins*2;i++)
		{
			bins[i] = ((float)i) / ((float)nbins);
			lc[i] = 0.;
		}
		
		for (i=0;i<nrows_FT1;i++)
		{
			for (j=0;j<nbins;j++)
			{
				if ((phase[i] >= (((float)j)/((float)nbins))) && (phase[i] < (((float)j+1.)/((float)nbins))))
				{
					lc[j]++;
					lc[j+nbins]++;
					break;
				}
			}
		}
		
		for (i=0;i<nbins;i++)
		{
			errX[i] = bins[i] + 0.5*(bins[1]-bins[0]);
			errX[i+nbins] = errX[i] + 1.;
			
			errY[i] = sqrt(lc[i]);
			errY[i+nbins] = errY[i];
			
			if (lc[i] > binmax)	binmax = lc[i];
		}
		
		binmax += sqrt(binmax);
		
		// Output the profile in a file
		if (output_prof)
		{
			output_proff = fopen(output_prof_file,"w+");

			for (i=0;i<nbins;i++)
			{
				fprintf(output_proff,"%f\t%f\t%f\n",bins[i],lc[i],errY[i]);
			}

			fclose(output_proff);
		}

		// ------------------------------------------------- //
		// H-test TS calculation
		// ------------------------------------------------- //	

		float Hbin[Hbins+1], HTS[Hbins+1], hmin = 1.e6, hmax = -1.e6;
		int Nphot[Hbins+1];
		
		Hbin[0] = tmin, HTS[0] = 0;
		
		for (i=1;i<Hbins+1;i++)
		{
			Hbin[i]  = tmin + ((float)i) * (tmax - tmin) / ((float)Hbins);
			Nphot[i] = 0;
			
			for (j=0;j<nrows_FT1;j++) 
			{
				if ((times[j] <= Hbin[i]) && (phase[j] >= 0.)) Nphot[i]++;
			}
			
			float phases[Nphot[i]];
			l = 0;
			
			for (j=0;j<nrows_FT1;j++) 
			{
				if ((times[j] <= Hbin[i]) && (phase[j] >= 0.))
				{
					phases[l] = phase[j];
					l++;
				}
			}
			
			HTS[i] = HTest(Nphot[i],phases);
			
			if (HTS[i] < hmin) hmin = HTS[i];
			if (HTS[i] > hmax) hmax = HTS[i];
		}
			
		float xmoy[2],ymoy[2];
		xmoy[0] = 0., xmoy[1] = 2.;
		ymoy[0] = (float)l/(float)nbins;
		ymoy[1] = ymoy[0];
		

		// ------------------------------------------------- //
		// Graphical output initialization
		// ------------------------------------------------- //							

		cpgbeg(0,gr,1,1);
		cpgask(0);

		// Invert foreground and background: black on white
		cpgscr(0,1,1,1); 
		cpgscr(1,0,0,0);
	
		if (strstr(gr,"/ps") != NULL)  linewidth = 2;
		if (strstr(gr,"/cps") != NULL) linewidth = 2;

		// ------------------------------------------------- //
		// First panel: phase histogram
		// ------------------------------------------------- //

		cpgsvp(0.05,0.5,0.55,0.975);
		cpgswin(0.,2.,0.,1.1*binmax);

		cpgsci(1);
		cpgslw(linewidth);
		cpgsch(fontsize);
		cpgbox("ABCNTS",0.5,5,"ABCNTS",pow(10,int(log10(binmax))),5);

		if (ophase == 0) cpglab("Pulse phase","Number of events","");
		else			 cpglab("Orbital phase","Number of events","");

		cpgsci(2);
		cpgbin(nbins*2,bins,lc,0);
		cpgerrb(6,nbins*2,errX,lc,errY,1.0);

		cpgsci(1);
		cpgsls(2);
		cpgline(2,xmoy,ymoy);
		cpgsls(1);

		// ------------------------------------------------- //
		// Second panel: time vs phase
		// ------------------------------------------------- //

		cpgsvp(0.6,0.95,0.1,0.975);
		cpgswin(0.,1.,tmin,tmax);

		cpgsci(1);

		cpgbox("ABCNTS",0.5,5,"ABCNTS1",0.,10);

		if (ophase == 0) cpglab("Pulse phase","Event time (MJD)","");
		else			 cpglab("Orbital phase","Event time (MJD)","");

		cpgsci(2);
		cpgslw(5);
		cpgpt(nrows_FT1,phase,times,-1);

		// ------------------------------------------------- //
		// Third panel: H-test TS vs time
		// ------------------------------------------------- //

		cpgsvp(0.05,0.5,0.1,0.45);
		cpgswin(tmin,tmax,0,hmax + 0.1*(hmax-hmin));

		cpgsci(1);
		cpgslw(linewidth);
		cpgbox("ABCNTS1",0.,10,"ABCNTS",0.,5);
		cpglab("Event time (MJD)","H-test TS","");
		cpgsci(2);
		cpgslw(5);
		cpgpt(Hbins+1,Hbin,HTS,3);
		cpgslw(1);
		cpgline(Hbins+1,Hbin,HTS);
			
		cpgend();
	}

	printf("Done with %s\n",psr[0].name);

	// ------------------------------------------------- //
	// Clean temporary tim file
	// ------------------------------------------------- //	

	sprintf(command,"rm -f %s",timFile[0]);
	system(command);
	
	return 0;
}





/* ------------------------------------------------- //
// Barycentering tools
// ------------------------------------------------- */

void clock_corrections_fermi(pulsar *psr,int npsr)
{
	// toa2utc(psr,npsr);                      // UTC(Observatory) -> UTC(NIST) 
	// tai2ut1(psr,npsr);                      // TAI -> UT1			 
	tt2tb(psr,npsr);                        // Rough estimate of TT-TB (+-2.2 microsec) 
}

void ephemeris_routines_fermi(pulsar *psr,int npsr)
{
	vectorPulsar(psr,npsr);                 // Form a vector pointing at the pulsar 
	readEphemeris(psr,npsr,0);              // Read the ephemeris 
	// get_obsCoord(psr,npsr);                 // Get Coordinate of observatory relative to Earth's centre 
	tt2tb(psr,npsr);                        // Observatory/time-dependent part of TT-TB 
	readEphemeris(psr,npsr,0);              // Re-evaluate ephemeris with correct TB 
}

void extra_delays_fermi(pulsar *psr,int npsr)
{
	calculate_bclt(psr,npsr);               // Calculate bclt
}

void formBatsAll_fermi(pulsar *psr,int npsr)
{
	clock_corrections_fermi(psr,npsr);            	// Clock corrections  ... 
	
	//	printf("Reading ephemeris routines...\n");
	ephemeris_routines_fermi(psr,npsr);           	// Ephemeris routines ... 
	
	//	printf("Reading extra delays...\n");
	extra_delays_fermi(psr,npsr);                 	// Other time delays  ... 
	
	//	printf("Forming barycentric arrival times...\n");
	formBats(psr,npsr);                     	// Form Barycentric arrival times 
	
	//	printf("Evaluating secular motion...\n");
	secularMotion(psr,npsr);
}

//


static char random_letter(int is_cap)
{
   int letter = (int)(26 * (rand() / (RAND_MAX + 1.0)));
   return((char)((is_cap == 1) ? (letter + 65) : (letter + 97)));
}

static char random_number()
{
   int number = (int)(10 * (rand() / (RAND_MAX + 1.0)));
   return((char)(number + 48));
}

static void random_string(int length, char *str)
{
	int i;
	int char_type;
	
	for(i = 0; i < length; i++)
	{
		char_type = (int)(3 * (rand() / (RAND_MAX + 1.0)));
	
		switch(char_type)
		{
			case 0:
				str[i] = random_letter(0);
				break;
			case 1:
				str[i] = random_letter(1);
				break;
			case 2:
				str[i] = random_number();
				break;
			default:
				str[i] = random_number();
				break;
		}
	}  
}

float HTest(int Nphotons, float phases[])
{
        int i_,j_;
        const int Nharm = 20;

        float alpha_i[Nharm+1];
        float beta_i[Nharm+1];
        float z2[Nharm+1];

        for (i_=0;i_<=Nharm;i_++)
        {
                alpha_i[i_]  = 0.;
                beta_i[i_]   = 0.;
                z2[i_]       = 0.;
        }

        for (i_=1;i_<=Nharm;i_++)
        {
                for (j_=0;j_<Nphotons;j_++)
                {
                        alpha_i[i_]  += cos(float(i_)*phases[j_]*2.*3.14159265359);
                        beta_i[i_]   += sin(float(i_)*phases[j_]*2.*3.14159265359);
                }

                alpha_i[i_]  = alpha_i[i_] / float(Nphotons);
                beta_i[i_]   = beta_i[i_]  / float(Nphotons);
        }

        for (i_=1;i_<=Nharm;i_++)
        {
                for (j_=1;j_<=i_;j_++)
                {
                        z2[i_]  += alpha_i[j_]*alpha_i[j_]+beta_i[j_]*beta_i[j_];
                }

                z2[i_] = 2.*Nphotons*z2[i_];
        }

        float h      = 0.;
        float max_h  = 0.;
        int harm_max = 0;

        for (i_=1;i_<=Nharm;i_++)
        {
                h =z2[i_] -(4.*i_)+4.;

                if (h > max_h)
                {
                        max_h    = h;
                        harm_max = i_;
                }
        }

        return max_h;
}

// version of cpgpt that handles longdouble arrays via copy
void cpgpt(int n, const float *xpts, const longdouble *ypts, int symbol)
{
    float farray[n];
    for (int cpy_ctr=0; cpy_ctr < n; cpy_ctr++)
    {
        farray[cpy_ctr] = float(ypts[cpy_ctr]);
    }
    cpgpt(n,xpts,farray,symbol);
}

char * plugVersionCheck = TEMPO2_h_VER;
