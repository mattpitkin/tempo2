#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include "../tempo2.h"
#include "../ifteph.h"
#include <fitsio.h>
#include <time.h>

using namespace std;        /* Is this required for a plugin ? Yes, for Linux */

static char random_letter(int is_cap);
static char random_number();
static void random_string(int length, char *str);

#define IAU_TEPH0 (-6.55e-5/86400) /* s */
#define IAU_K 1.550519768e-8
#define IAU_KINV 1.55051974395888e-8 /* 1 - 1/(1-IAU_K) */
longdouble tcb2tdb(longdouble mjd) {
    return mjd - IAU_K*(mjd-IFTE_MJD0) + IAU_TEPH0;
    // same as (1-IAU_K)*mjd + IAU_K*IFTE_MJD0 + IFTE_TEPH0
}
longdouble tdb2tcb(longdouble mjd) {
    // Few hundred nanosec roundtrip error
    mjd -= IAU_TEPH0;
    return mjd + IAU_KINV*(mjd-IFTE_MJD0);
}

int find_event_hdu(fitsfile*ft_in) {
    int event_hdu;
    int status = 0;
    if (fits_movnam_hdu(ft_in, ANY_HDU, (char*)"EVENTS", 0, &status)) {
        fits_report_error(stderr, status);
        fprintf(stderr, "Unable to find HDU called 'EVENTS'; try using -hdu to specify the number\n");
        exit(13);
    }
    fits_get_hdu_num(ft_in, &event_hdu);
    return event_hdu;
}
longdouble get_mjdref(fitsfile*ft_in) {
    int status = 0;
    int mjdreff_set = 0;
    int mjdrefi_set = 0;
    double mjdref;
    double mjdrefi;
    double mjdreff;
    char comment[81];
    if(fits_movabs_hdu(ft_in,1,NULL,&status)) {
        fits_report_error(stderr,status);
        exit(17);
    }
    while (1) {

        status = 0;
        if (!fits_read_key(ft_in, TDOUBLE, (char*)"MJDREFI", 
                    (void*)&mjdrefi, comment, &status)) {
            mjdrefi_set = 1;
            if (mjdreff_set)
                return ((longdouble)mjdrefi)+mjdreff;
        }

        status = 0;
        if (!fits_read_key(ft_in, TDOUBLE, (char*)"MJDREFF", 
                    (void*)&mjdreff, comment, &status)) {
            mjdreff_set = 1;
            if (mjdrefi_set)
                return ((longdouble)mjdrefi)+mjdreff;
        }

        status = 0;
        if (!fits_read_key(ft_in, TDOUBLE, (char*)"MJDREF", 
                    (void*)&mjdref, comment, &status)) {
            return mjdref;
        }

        status=0;
        if(fits_movrel_hdu(ft_in,1,NULL,&status)) {
            if (status==END_OF_FILE)
                break;
            fits_report_error(stderr,status);
            exit(18);
        }
    }

    fprintf(stderr,"Error: MJDREF or MJDREFI/MJDREFF not found; try specifying -mjdref\n");
    exit(3);
}
void check_barycentered(fitsfile*ft_in, int event_hdu) {
    char comment[81];
    char value[81];
    int status = 0;
    if (fits_movabs_hdu(ft_in,event_hdu,NULL,&status)) {
        fits_report_error(stderr,status);
        exit(14);
    }
    status = 0;
    if (fits_read_key(ft_in, TSTRING, (char*)"TIMESYS", 
            (void*)value, comment, &status)) {
        fprintf(stderr,"Assuming TIMESYS is TDB\n");
    } else {
        if (strcmp(value,"TDB")) {
            fprintf(stderr,"Error: TIMESYS is %s rather than TDB; input file has not been barycentered\n", value);
            exit(2);
        }
    } 
    status = 0;
    if( fits_read_key(ft_in, TSTRING, (char*)"TIMEREF", 
            (void*)value, comment, &status)) {
        fprintf(stderr,"Warning: assuming TIMEREF is SOLARSYSTEM\n");
    } else {
        if (strcmp(value,"SOLARSYSTEM")) {
            fprintf(stderr,"Error: TIMEREF is %s rather than SOLARSYSTEM\n", value);
            exit(2);
        }
    }
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
    // Provide basic information
    printf("\n");
    printf("------------------------------------------\n");
    printf("Output interface:    photons\n");
    printf("Author:              Anne Archibald\n");
    printf("                     based on the fermi plugin v5.2\n");
    printf("                     by Lucas Guillemot\n");
    printf("Updated:             30 November 2010\n");
    printf("Version:             0.3\n");
    printf("------------------------------------------\n");
    printf("\n");

    int i,j,k,l,event = 0,event2;

    char parFile[1][MAX_FILELEN];
    char timFile[1][MAX_FILELEN];
    char temptim[9];

    srand((unsigned)time(NULL));
    temptim[8] = '\0';
    random_string(8,temptim);
    strcpy(timFile[0],temptim);
    strcat(timFile[0],".tim");

    char FT_in[MAX_FILELEN];
    char FT_out[MAX_FILELEN];
    char output[MAX_FILELEN];
    char phasecol[32];
    strcpy(phasecol,"PULSE_PHASE");
    char timecol[32];
    strcpy(timecol,"TIME");
    char command[128];
    
    char error_buffer[128];
    char history[128];
    int  par_file      = 0;
    int  in_file      = 0;
    
    int  ophase        = 0;
    int  graph         = 0;
    int  output_file   = 0;
    int  phase_replace = 0;
    int  phase_col_set = 0;

    double intpart;

    FILE *outputf;
    FILE *temp_tim;

    /* ------------------------------------------------- //
    // fits file definitions
    // ------------------------------------------------- */

    fitsfile *ft_in;
    fitsfile *ft_out;

    char colname[FLEN_VALUE];
    int open_status = 0, status = 0, status2 = 0, hdupos, ncols_FT1,anynul;
    long nrows_FT1;
    double nulval = 0.;

    int FT1_time_col,FT1_phase_col;
    int nrows2, rows_status, rows_left;
    int max_rows = MAX_OBSN_VAL - 100;
    int event_hdu;
    int event_hdu_set = 0;

    /* ------------------------------------------------- //
    // Time and satellite position definitions
    // ------------------------------------------------- */

    longdouble mjd_ref = 51910.0007428703703703703; 
    int mjd_ref_set = 0;
    int timesys_found = 0;
    int timeref_found = 0;
    int phasecol_exists = 0;
    double temptime;
    double minFT1time = 999999999., maxFT1time = 0.;
    
    double time_MET_TDB[max_rows];
    longdouble time_MJD_TDB;
    double obs_earth[max_rows][3];

    longdouble lasttime, tzrmjd_bary;
    double lastpos[3];
    
    double tpb;

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
            if (++i>=argc) {
                fprintf(stderr, "Option -f requires an argument\n");
                exit(5);
            }
            strcpy(parFile[0],argv[i]); 
        }
        else if (strcmp(argv[i],"-in")==0)
        {
            in_file = 1;
            if (++i>=argc) {
                fprintf(stderr, "Option -in requires an argument\n");
                exit(5);
            }
            strcpy(FT_in,argv[i]);
        }
        else if (strcmp(argv[i],"-textoutput")==0)
        {
            output_file = 1;
            if (++i>=argc) {
                fprintf(stderr, "Option -textoutput requires an argument\n");
                exit(5);
            }
            strcpy(output,argv[i]);
        }
        else if (strcmp(argv[i],"-fitsoutput")==0)
        {
            phase_replace = 1;
            if (++i>=argc) {
                fprintf(stderr, "Option -fitsoutput requires an argument\n");
                exit(5);
            }
            strcpy(FT_out,argv[i]);
        }
        else if (strcmp(argv[i],"-ophase")==0)
        {
            ophase = 1;
        }
        else if (strcmp(argv[i],"-col")==0)
        {
            phase_col_set = 1;
            if (++i>=argc) {
                fprintf(stderr, "Option -col requires an argument\n");
                exit(5);
            }
            strcpy(phasecol,argv[i]);
        }   
        else if (strcmp(argv[i],"-hdu")==0)
        {
            event_hdu_set = 1;
            if (++i>=argc) {
                fprintf(stderr, "Option -hdu requires an argument\n");
                exit(5);
            }
            event_hdu = atoi(argv[i]);
        }   
        else if (strcmp(argv[i],"-timecol")==0)
        {
            if (++i>=argc) {
                fprintf(stderr, "Option -timecol requires an argument\n");
                exit(5);
            }
            strcpy(timecol,argv[i]);
        }   
        else if (strcmp(argv[i],"-mjdref")==0)
        {
            if (++i>=argc) {
                fprintf(stderr, "Option -mjdref requires an argument\n");
                exit(5);
            }
            sscanf(argv[i],"%Lf",&mjd_ref);
            mjd_ref_set = 1;
        }   
        else if (strcmp(argv[i],"-h")==0)
        {
            // FIXME: option to accept a time column in MJDs
            printf("\n TEMPO2 photons plugin\n");
            printf("========================\n");
            printf("\n USAGE: \n\t tempo2 -gr photons -in FT1.fits -f par.par\n");
            printf("\n Command line options:\n");
            printf("\t -in XXX: read events from fits file XXX\n");
            printf("\t -textoutput XXX: writes event times and phases in text file XXX\n");
            printf("\t -fitsoutput XXX: stores calculated phases in a column of XXX\n");
            printf("\t -ophase: will calculate orbital phases instead of pulse phases. Changes default column name to ORBITAL_PHASE.\n");
            printf("\t -col XXX: phases will be stored in column XXX. Default is PULSE_PHASE\n");
            printf("\t -timecol XXX: times will be read from column XXX. Default is TIME\n");
            printf("\t -mjdref XXX: use XXX as MJDREF instead of looking in header\n");
            printf("\t -h: this help.\n");
            printf("===============================================\n");
            exit(0);            
        }
    }
    
    if (!in_file)
    {
        printf("No input FITS file !\n");
        printf("At least it should look like :\n");
        printf("\t tempo2 -gr photons -in FT1.fits -f par.par\n");
        printf("More options are available. If you need help, please type:\n");
        printf("\t tempo2 -gr photons -h\n");
        exit(1);
    }
    if (!par_file)
    {
        printf("No ephemeris file !\n");
        printf("At least it should look like :\n");
        printf("\t tempo2 -gr photons -in FT1.fits -f par.par\n");
        printf("More options are available. If you need help, please type:\n");
        printf("\t tempo2 -gr photons -h\n");
        exit(1);
    }
    if (!phase_replace && !output_file) {
        printf("Warning: neither -fitsoutput nor -textoutput were selected\n");
        printf("so while phases will be calculated they will not be\n");
        printf("recorded anywhere.\n");
        // FIXME: should this be an error?
    }
    if (output_file)
    {
        outputf = fopen(output,"w+");
    }
    if (ophase)
    {
        if (!phase_col_set) strcpy(phasecol,"ORBITAL_PHASE");
    }

    /* ------------------------------------------------- //
    // FT1 file
    // ------------------------------------------------- */

    if (!fits_open_file(&ft_in,FT_in, READONLY, &open_status))
    {

        if (!event_hdu_set) {
            event_hdu = find_event_hdu(ft_in);
        }
        check_barycentered(ft_in, event_hdu);
        mjd_ref = get_mjdref(ft_in);
        
        fits_movabs_hdu(ft_in,event_hdu,NULL,&status);

        fits_get_num_rows(ft_in, &nrows_FT1, &status);
        fits_get_num_cols(ft_in, &ncols_FT1, &status);

        fits_get_colname(ft_in,CASEINSEN,timecol,colname,&FT1_time_col,&status);
        if (status>0) {
            int c, c_status;
            char colname[80];
            char *templt=(char*)"*";
            fprintf(stderr, "No column called %s in %s; try using -timecol\n",
                    timecol, FT_in);
            fprintf(stderr, "Column names are:\n");
            for (c=1;c<=ncols_FT1;c++) {
                fits_get_colname(ft_in, CASESEN, templt, colname, &c, &c_status);
                fprintf(stderr,"\t%s\n",colname);
            }
            exit(4);
        }
        if (phase_replace) {
            int c, c_status=0;
            char colname[80];
            fits_get_colname(ft_in,CASEINSEN,phasecol,colname,&c,&c_status);
            if (c_status>0) {
                phasecol_exists = 0;
            } else {
                phasecol_exists = 1;
            }
        }
        
        //
        
        for (i=1;i<=nrows_FT1;i++)
        {
            fits_read_col_dbl(ft_in,FT1_time_col,i,1,1,nulval,&temptime,&anynul,&status);
            if (temptime < minFT1time) minFT1time = temptime;
            if (temptime > maxFT1time) maxFT1time = temptime;
        }
    }
    
    if (open_status != 0)
    {
        fprintf(stderr,"Can't read input file '%s'!\n",FT_in);
        fits_report_error(stderr,status);
        exit(-1);
    }

    fits_close_file(ft_in, &status);

    printf("MJDREF: %Lf\n",mjd_ref);
    printf("First photon date in input: %lf MET (s), MJD %Lf\n",minFT1time,minFT1time/86400.+mjd_ref);
    printf(" Last photon date in input: %lf MET (s), MJD %Lf\n\n",maxFT1time,maxFT1time/86400.+mjd_ref);
    
    // Copy input file to output file, possibly adding room for a new column
    if (phase_replace) {
        char FT_in_filter[MAX_FILELEN+20];
        int status=0, ii=1;

        if (phasecol_exists) {
            printf("Replacing existing column %s\n", phasecol);
            strcpy(FT_in_filter, FT_in);
        } else {
            printf("Adding new column %s\n", phasecol);
            sprintf(FT_in_filter, "%s[%d][col *;%s = 0.0]", FT_in, event_hdu-1, phasecol);
        }

        printf("Copying photons from %s to %s\n", FT_in_filter, FT_out);
        status = 0;
        if (!fits_open_file(&ft_in, FT_in_filter, READONLY, &status)) {
            // Record the column number for later
            if (fits_movabs_hdu(ft_in,event_hdu,NULL,&status) || fits_get_colname(ft_in,CASEINSEN,phasecol,colname,&FT1_phase_col,&status)) {
                fprintf(stderr, "Bizarre: unable to find phase column %s even after trying to create it.\n", phasecol);
                fits_report_error(stderr, status);
                exit(11);
            }
            if (unlink(FT_out)) {
                if (errno!=ENOENT) {
                    perror("Problem deleting old output file");
                    exit(15);
                }
            } else {
                printf("Deleted old version of %s\n", FT_out);
            }
            if (!fits_create_file(&ft_out, FT_out,  &status)) {
                // Copy all HDUs one-by-one
                // In the process, make room for more header comments
                status = 0;
                i = 1;
                while (!fits_movabs_hdu(ft_in, i++, NULL, &status))
                    fits_copy_hdu(ft_in, ft_out, 5, &status);

                if (status==END_OF_FILE) status=0;
                else printf("Error happened while transferring\n");

                fits_close_file(ft_out, &status);
            }
            fits_close_file(ft_in, &status);
        }
        if (status) {
            fits_report_error(stderr, status);
            exit(12);
        }
        printf("Copied.\n");
    }
    // The FT1 file is cut into small tim files with # rows <= 10000

    rows_status = 1;
    rows_left   = nrows_FT1 - rows_status + 1;
    nrows2      = min(rows_left,max_rows);

    float *phase;
    float *times;
        
        phase  = (float *)calloc(max_rows,sizeof(float));
        times  = (float *)calloc(max_rows,sizeof(float));
    
    float tmin   = 100000., tmax   = -100000.;
    
    
    /* ------------------------------------------------- //
    // Barycentric TZRMJD
    // ------------------------------------------------- */
    
    // A temporary file is created. It is first used to get the barycentered TZR
    temp_tim = fopen(timFile[0],"w+");
    fprintf(temp_tim,"FORMAT 1\n");
    fprintf(temp_tim," photons 0.0 %.12lf 0.00000 BAT\n",1.);
    fclose(temp_tim);

    // Load the arrival times
    readParfile(psr,parFile,timFile,*npsr);
    readTimfile(psr,timFile,*npsr);
    
    if (ophase && (psr[0].param[param_pb].paramSet[0] == 0))
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
    if (fits_open_file(&ft_in,FT_in, READONLY, &status)) {
        fprintf(stderr, "Problem reading from %s\n", FT_in);
        exit(6);
    }
    fits_movabs_hdu(ft_in,event_hdu,NULL,&status);

    if (phase_replace) {
        if (fits_open_file(&ft_out,FT_out, READWRITE, &open_status)) {
            fprintf(stderr, "Problem writing to %s\n", FT_out);
            fits_report_error(stderr, open_status);
            exit(7);
        }
        fits_movabs_hdu(ft_out,event_hdu,NULL,&status);
    }

    while (rows_left > 0)
    {
        printf("Treating events # %d to %d... \n",rows_status,rows_status + nrows2 - 1);

        // Acquisition of the TOAs in MET TDB
        j = 0;

        for (i=rows_status;i<rows_status + nrows2;i++)
        {
            fits_read_col_dbl(ft_in,FT1_time_col,i,1,1,nulval,&time_MET_TDB[j],&anynul,&status);
            j++;
        }



        // A temporary file is created. TOAs in MJD TT are stored is this file
        temp_tim = fopen(timFile[0],"w+");
        fprintf(temp_tim,"FORMAT 1\n");

        for (i=0;i<nrows2;i++)
        {
            time_MJD_TDB = time_MET_TDB[i]/86400.+mjd_ref;      
            if (psr[0].units == TDB_UNITS) {
                fprintf(temp_tim," photons 0.0 %.12Lf 0.00000 @\n",time_MJD_TDB);
            } else {
                fprintf(temp_tim," photons 0.0 %.12Lf 0.00000 @\n",tdb2tcb(time_MJD_TDB));
                //printf("TDB->TCB->TDB roundtrip error: %Lg s\n",
                //       (tcb2tdb(tdb2tcb(time_MJD_TDB))-time_MJD_TDB)*86400);
            }
        }

        fclose(temp_tim);

        // Load the arrival times
        readParfile(psr,parFile,timFile,*npsr); 
        readTimfile(psr,timFile,*npsr);
        
        


        /* ------------------------------------------------- //
        // Calculation of the event phases - step 1
        // Step 1 is for all but the last photon
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
        formBatsAll(psr,*npsr);
        //printf("%Lf\t%Lf\n",psr->obsn[0].sat,psr->obsn[0].bat);
        //printf("%Lf\t%Lf\n",psr->obsn[1].sat,psr->obsn[1].bat);

        // ------------------------------------------------- //
        // Calculate event phases - step 1
        // ------------------------------------------------- //     
        formResiduals(psr,*npsr,0.0);

        for (i=1;i<nrows2;i++)
        {
                        if (ophase == 0)
                        {
                                phase[event] = modf(psr[0].obsn[i].phase,&intpart);         
                                if (phase[event] < 0.) phase[event]++;
                        }
                        else
                        {
                                if(psr[0].param[param_tasc].paramSet[0])
                                {
                                        tpb = (psr[0].obsn[i].bat-psr[0].param[param_tasc].val[0])/(psr[0].param[param_pb].val[0]);
                                }
                                else if(psr[0].param[param_t0].paramSet[0])
                                {
                                        tpb = (psr[0].obsn[i].bat-psr[0].param[param_t0].val[0])/(psr[0].param[param_pb].val[0]);
                                }
                                
                                phase[event] = modf(tpb+1000000.0,&intpart);;
                                if (phase[event] < 0.) phase[event]++;
                        }

                        times[event] = psr[0].obsn[i].bat;
                        if (times[event] > tmax)    tmax = times[event];
                        if (times[event] < tmin)    tmin = times[event];
                
            if (output_file)
            {
                fprintf(outputf,"%d\t",event + rows_status);
                fprintf(outputf,"%20.15Lf %12.10le\n",psr[0].obsn[i].bat,phase[event]);
            }
    
            event++;
        }       

        // ------------------------------------------------- //
        // Calculation of the event phases - step 2
        // Step 2 repeats the above code for the last photon
        // ------------------------------------------------- //
        psr[0].obsn[1].sat = lasttime;
        for (k=0;k<3;k++) psr[0].obsn[1].observatory_earth[k] = lastpos[k];
    
        // ------------------------------------------------- //
        // Form barycentric arrival times - step 2
        // ------------------------------------------------- //
        formBatsAll(psr,*npsr);

        // ------------------------------------------------- //
        // Calculate event phases - step 2
        // ------------------------------------------------- //     
        formResiduals(psr,*npsr,0.0);
    
                if (ophase == 0)
                {
                        phase[event] = modf(psr[0].obsn[1].phase,&intpart);
                        if (phase[event] < 0.) phase[event]++;
                }
                else
                {
                        if(psr[0].param[param_tasc].paramSet[0])
                        {
                                tpb = (psr[0].obsn[1].bat-psr[0].param[param_tasc].val[0])/(psr[0].param[param_pb].val[0]);
                        }
                        else if(psr[0].param[param_t0].paramSet[0])
                        {
                                tpb = (psr[0].obsn[1].bat-psr[0].param[param_t0].val[0])/(psr[0].param[param_pb].val[0]);
                        }
                        
                        phase[event] = modf(tpb+1000000.0,&intpart);;
                        if (phase[event] < 0.) phase[event]++;
                }
        
                times[event] = psr[0].obsn[1].bat;      
                if (times[event] > tmax)    tmax = times[event];
                if (times[event] < tmin)    tmin = times[event];

        if (output_file)
        {
            fprintf(outputf,"%d\t",event + rows_status);
            fprintf(outputf,"%20.15Lf %12.10le\n",psr[0].obsn[1].bat,phase[event]);
        }

        // ------------------------------------------------- //
        // Put new phases into the input FT1 file
        // ------------------------------------------------- //
        if (phase_replace)
        {
            for (event2=rows_status;event2<rows_status+nrows2;event2++)
            {
                if (graph == 0) i = event2 - rows_status;
                else i = event2 - 1;
            
                fits_write_col_flt(ft_out,FT1_phase_col,event2,1,1,&phase[i],&status);
            }

        }

        /* ------------------------------------------------- //
        // End of the loop
        // ------------------------------------------------- */ 
        rows_status += nrows2;
        rows_left   = nrows_FT1 - rows_status + 1;
        nrows2      = min(rows_left,max_rows);
        
                event = 0;
    }
    if (phase_replace)
        fits_close_file(ft_out, &status);      
    fits_close_file(ft_in, &status);
        
    // ------------------------------------------------- //
    // Add a little bit of history to the header
    // ------------------------------------------------- //
    if (phase_replace)
    {
            sprintf(history,"Pulse phases calculated with the TEMPO2 photons plugin using ephemeris %s",parFile[0]);
            
            if (!fits_open_file(&ft_out,FT_out, READWRITE, &open_status))
            {
                    status = 0;
                    fits_movabs_hdu(ft_out,event_hdu,NULL,&status);
                    fits_write_history(ft_out,history,&status);
                    
                    if (status != 0)
                    {
                            fits_get_errstatus(status,error_buffer);
                            printf( "fits_insert_col: %s\n", error_buffer);
                            exit(-1);
                    }
            }
            
            fits_close_file(ft_out, &status);
    }

    if (output_file) fclose(outputf);

    printf("Done with %s\n",psr[0].name);

    // ------------------------------------------------- //
    // Clean temporary tim file
    // ------------------------------------------------- // 

    sprintf(command,"rm -f %s",timFile[0]);
    system(command);
    
    return 0;
}








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

