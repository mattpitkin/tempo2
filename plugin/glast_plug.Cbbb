/* GLAST tempo2 plugin */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include <cpgplot.h>
using namespace std;

void sla_CALDJ(int IY, int IM, int ID, double *DJM, int *J);
void sla_CLDJ (int IY, int IM, int ID, double *DJM, int *J);
void slaClyd ( int iy, int im, int id, int *ny, int *nd, int *jstat );
void slaCalyd ( int iy, int im, int id, int *ny, int *nd, int *j );
int getParameter(pulsar psr,char *param, double *value);
float fitwave_function(pulsar *psr, float x, float fitwaves_omega, float fitwaves_epoch);
void indexx_patrick(unsigned long n, float arr[], unsigned long indx[]);

void help() /* Display help */
{
  /* This function should contain usage information about the plugin which should (in general) be accessed */
  /* by the user pressing 'h'                                                                              */
    printf("-D        Specify ploting device.\n");
    printf("-w        Specify window width and height (in pixels).\n");
    printf("-p        Also show this alternative par file.\n");
    printf("-P        Write out post fit par file (newpar.par).\n");
    printf("-c        Color code the frequencies.\n");
    printf("-t        Write text in image.\n");
    printf("-fw       Overplot fitwaves fit.\n");
    printf("-bat      Use bat instead of sat to calculate date.\n");
    printf("-dewrap   Try to solve wraps.\n");
    printf("-pgscript Write out a pgplot script instead of calling pgplot.\n");
    printf("-v        verbose.\n");
}


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  char months[12][10];
  int i, j, k, n, bat, pgscript;
  double globalParameter;
  float *resids[2], *resids_err[2], *resids_freq[2], *sats[2], *sats2[2], mean_resids;
  float x, y, y_10, y_20, y_50,y_40, ymin, ymax, ymin2, ymax2, xmin, xmax, rms[2], ticklength;
  int deviceID, windowwidth, windowheight, colorcode, plot_fitwaves, verbose;
  int solvewraps, plotmonths;
  char PlotDevice[100], txt[1000], parfile2[MAX_PSR][MAX_FILELEN], comment_str[100], *resids_backend[2];
  int yr, mt, day, hour, min, sec, npar, outputNewPar, n_10, n_20, n_50, n_40; 
  double djm, glitch_epoch, fitwaves_omega, fitwaves_epoch, mean2, period, pdot;
  unsigned long *indx[2];
  FILE *fin, *fout_pgplot;

  printf("Graphical Interface: glast\n");
  printf("Author:              P. Weltevrede\n");
  printf("Version:             1.0\n");

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */
  outputNewPar = 0;
  strcpy(PlotDevice, "?");
  windowwidth = windowheight = -1;
  comment_str[0] = 0;
  colorcode = 0;
  bat = 0;
  solvewraps = 0;
  plot_fitwaves = 0;
  sprintf(months[0], "Jan");
  sprintf(months[1], "Feb");
  sprintf(months[2], "Mar");
  sprintf(months[3], "Apr");
  sprintf(months[4], "May");
  sprintf(months[5], "Jun");
  sprintf(months[6], "Jul");
  sprintf(months[7], "Aug");
  sprintf(months[8], "Sep");
  sprintf(months[9], "Oct");
  sprintf(months[10], "Nov");
  sprintf(months[11], "Dec");
  verbose = 0;
  pgscript = 0;

  /* Obtain the .par and the .tim file from the command line */
  if (argc==4) /* Only provided .tim name */
    {
      strcpy(timFile[0],argv[3]);
      strcpy(parFile[0],argv[3]);
      parFile[0][strlen(parFile[0])-3] = '\0';
      strcat(parFile[0],"par");
    }
  /* Obtain all parameters from the command line */
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0) {
	strcpy(parFile[0],argv[i+1]); 
	strcpy(timFile[0],argv[i+2]);
	i += 2;
      }else if(strcmp(argv[i], "-D") == 0) {
	strcpy(PlotDevice, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-h") == 0) {
	help();
	return 0;
      }else if(strcmp(argv[i], "-v") == 0) {
	verbose = 1; 
      }else if(strcmp(argv[i], "-dewrap") == 0) {
	solvewraps = 1; 
      }else if(strcmp(argv[i], "-gr") == 0) {
	/* Ignore the -gr option */
	i += 1; 
      }else if(strcmp(argv[i], "-epoch") == 0) {
	/* Ignore the -epoch option */
	i += 1; 
      }else if(strcmp(argv[i], "-p") == 0) {
	strcpy(parfile2[0], argv[i+1]);
	fin = fopen(parfile2[0], "r");
	if(fin == NULL) {
	  fprintf(stderr, "WARNING: CANNOT OPEN %s\n", parfile2[0]);
	  parfile2[0][0] = 0;
	}else {
	  fclose(fin);
	}
	i++;
      }else if(strcmp(argv[i], "-t") == 0) {
	strcpy(comment_str, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-P") == 0) {
	outputNewPar = 1;
      }else if(strcmp(argv[i], "-c") == 0) {
	colorcode = 1;
      }else if(strcmp(argv[i], "-fw") == 0) {
	plot_fitwaves = 1;
      }else if(strcmp(argv[i], "-bat") == 0) {
	bat = 1;
      }else if(strcmp(argv[i], "-pgscript") == 0) {
	pgscript = 1;
      }else if(strcmp(argv[i], "-w") == 0) {
	j = sscanf(argv[i+1], "%d %d", &windowwidth, &windowheight);
	if(j != 2) {
	  fprintf(stderr, "Error parsing -w option\n");
	  return 0;
	}
	i++;
      }else {
	fprintf(stderr, "Didn't recognize option %s\n", argv[i]);
	return 0;
      }
    }
  if(pgscript) {
    fout_pgplot = fopen("residuals.pgplot", "w");
    if(fout_pgplot == NULL) {
      fprintf(stderr, "Cannot open residuals.pgplot\n");
      return 0;
    }
  }
  for(npar = 0; npar < 2; npar++) {
    if(npar == 0) /* Load the parameters       */
      readParfile(psr,parFile,timFile,*npsr); 
    else
      readParfile(psr,parfile2,timFile,*npsr); 
    readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
    preProcess(psr,*npsr,argc,argv);
    for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
      {
        if(verbose) printf("There are %d timing points.\n", psr[0].nobs);
	formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
	formResiduals(psr,*npsr,1);    /* Form the residuals                 */
	if (i==0) doFit(psr,*npsr,0);   /* Do the fitting     */
	else textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
      }
    if(npar == 0) {
      period = 1.0/(double)psr[0].param[param_f].val[0];
      pdot = -(double)psr[0].param[param_f].val[1]*period*period;
      printf("Sensible parameters: name=%s P=%f Pdot=%e age=%e edot=%e\n",psr[0].name, period, pdot, (0.5*period/pdot)/(365.25*24.0*3600.0), 4.0*M_PI*M_PI*1e45*pdot/(period*period*period)); 
      /*
      printf("RAJ = %s (pre) %s (post)\n", psr[0].rajStrPre, psr[0].rajStrPost);
      printf("pos = %f %f %f\n", psr[0].posPulsar[0], psr[0].posPulsar[1], psr[0].posPulsar[2]);
      */
      if(verbose) printf("There are %d timing points.\n", psr[0].nobs);
      for (i=0;i<psr[0].nobs;i++)
	{
	  if(bat == 1) {
	    if(verbose) printf("bat = %lf residual = %lf phase = %lf toaErr = %lf filename = %s\n",(double)psr[0].obsn[i].bat,(double)psr[0].obsn[i].residual,(double)psr[0].obsn[i].residual*(double)psr[0].param[param_f].val[0],(double)psr[0].obsn[i].toaErr, psr[0].obsn[i].fname);
	  }else {
	    if(verbose) printf("sat = %lf, residual = %lf, toaErr = %lf, filename = %s\n",(double)psr[0].obsn[i].sat,(double)psr[0].obsn[i].residual,(double)psr[0].obsn[i].toaErr, psr[0].obsn[i].fname);
	    
	    strcpy(txt,psr[0].obsn[i].fname+1);
	    txt[2] = 0;
	    sscanf(txt, "%d", &yr);
	    strcpy(txt,psr[0].obsn[i].fname+3);
	    txt[2] = 0;
	    sscanf(txt, "%d", &mt);
	    strcpy(txt,psr[0].obsn[i].fname+5);
	    txt[2] = 0;
	    sscanf(txt, "%d", &day);
	    yr += 2000;
	    sla_CALDJ(yr, mt, day, &djm, &j);
	    if(j != 0) {
	      printf("SLA ERROR %d (%d-%d-%d)\n", j, yr,mt,day);
	      return 0;
	    }
	    if(verbose) printf("mjd = %lf (%d-%d-%d)\n", djm, yr, mt, day);
	    /* It can be out by more as one day if observation starts a few minutes before midnight. This is because difference in mjd of middle of observation and the start of the observation. */
	    if(fabs(psr[0].obsn[i].sat - djm) > 2.01) {
	      printf("SLA ERROR: More than one day difference between sat and observing day.\n");
	      return 0;
	    }
	  }
	}
    }


    resids[npar] = (float *)malloc(sizeof(float)*(psr[0].nobs));
    resids_err[npar] = (float *)malloc(sizeof(float)*(psr[0].nobs));
    resids_freq[npar] = (float *)malloc(sizeof(float)*(psr[0].nobs));
    resids_backend[npar] = (char *)malloc(sizeof(float)*(psr[0].nobs));
    sats[npar] = (float *)malloc(sizeof(float)*(psr[0].nobs));
    sats2[npar] = (float *)malloc(sizeof(float)*(psr[0].nobs));
    indx[npar] = (unsigned long *)malloc(sizeof(unsigned long)*(psr[0].nobs));
    if(resids[npar] == NULL || resids_err[npar] == NULL || resids_freq[npar] == NULL || resids_backend[npar] == NULL || sats[npar] == NULL || sats2[npar] == NULL || indx[npar] == NULL) {
      fprintf(stderr, "Cannot allocate enough memory.\n");
      return 0;
    }

    for (i = 0; i < psr[0].nobs; i++) {
      if(bat == 1)
	sats2[npar][i] = (float)psr[0].obsn[i].bat;
      else
	sats2[npar][i] = (float)psr[0].obsn[i].sat;
    }
    indexx_patrick(psr[0].nobs, &sats2[npar][0]-1, &indx[npar][0] - 1);
    for (i = 0; i < psr[0].nobs; i++) {
      sats[npar][i] = sats2[npar][indx[npar][i]-1];
      resids[npar][i] = (float)psr[0].obsn[indx[npar][i]-1].residual;
      resids_err[npar][i] = (float)psr[0].obsn[indx[npar][i]-1].toaErr*0.000001;
      resids_freq[npar][i] = (float)psr[0].obsn[indx[npar][i]-1].freq;
      resids_backend[npar][i] = psr[0].obsn[indx[npar][i]-1].fname[0];
    }
    if(solvewraps) {
      float dy, dx, dy2, dx2, dy3, totaljump, residold;
      int dojump;
      totaljump = 0;
      dy2 = 0;
      dx2 = 1;
      for (i = 0; i < psr[0].nobs; i++) {
	residold = resids[npar][i];
	resids[npar][i] += totaljump/psr[0].param[param_f].val[0];
	if(i > 1 && psr[0].nPhaseJump < 1e10) {
	  dy = resids[npar][i] - resids[npar][i-1];
	  dx = sats[npar][i] - sats[npar][i-1];
	  dy *= psr[0].param[param_f].val[0];
	  if(verbose) printf("dy=%f at %f (%f - %f)\n", dy, (float)sats[npar][i], (float)(resids[npar][i]*psr[0].param[param_f].val[0]), (float)(resids[npar][i-1]*psr[0].param[param_f].val[0]));
	  dojump = 0;
	  if(dx < 10) {
	    if(fabs(dy) > 0.4) {
	      if(dy > 0)
		dojump = -1;
	      else
		dojump = 1;
	    }
	  }else if(i > 2) {
	    dy3 = dy-(dy2*dx/dx2);
	    if(verbose) printf("  dy3 = %f\n", dy3);
	    if(fabs(dy3) > 0.4) {
	      if(dy3 > 0)
		dojump = -1;
	      else
		dojump = 1;
	    }
	    dy += dojump;
	  }
	  if(dojump != 0) {
	    psr[0].phaseJump[psr[0].nPhaseJump] = 0.5*(sats[npar][i] - sats[npar][i-1])+sats[npar][i-1];
	    psr[0].phaseJumpDir[psr[0].nPhaseJump] = dojump;
	    totaljump += (float)psr[0].phaseJumpDir[psr[0].nPhaseJump];
	    resids[npar][i] = residold + totaljump/psr[0].param[param_f].val[0];
	    if(verbose)
	      printf("*** Add jump at %f (%d -> %.0f %f)\n", psr[0].phaseJump[psr[0].nPhaseJump], psr[0].phaseJumpDir[psr[0].nPhaseJump], totaljump, (float)(resids[npar][i]*psr[0].param[param_f].val[0]));
	    psr[0].nPhaseJump++;
	  }
	}
	if(dx > 10) {
	  dx2 = dx;
	  dy2 = dy;
	}
      }
    }

    y_10 = 0;
    y_20 = 0;
    y_40 =0;
    y_50 = 0;
    n_10 = 0;
    n_20 = 0;
    n_40 =0;
    n_50 = 0;
    mean_resids = 0;
    for (i = 0; i < psr[0].nobs; i++) {
      mean_resids += psr[0].obsn[i].residual;
      if(resids_freq[npar][i] > 600 && resids_freq[npar][i] < 700) {
	y_50 += resids[npar][i];
	n_50++;
      }
      if(resids_freq[npar][i] > 700 && resids_freq[npar][i] < 800) {
	y_40 += resids[npar][i];
	n_40++;
      }
      else if(resids_freq[npar][i] > 1300 && resids_freq[npar][i] < 1500) {
	y_20 += resids[npar][i];
	n_20++;
      }else if(resids_freq[npar][i] > 3000 && resids_freq[npar][i] < 3300) {
	y_10 += resids[npar][i];
	n_10++;
      }
    }
    if(n_10 > 0)
      y_10 /= (float)n_10;
    if(n_20 > 0)
      y_20 /= (float)n_20;
    if(n_50 > 0)
      y_50 /= (float)n_50;
    if(n_40 > 0)
      y_40 /= (float)n_40;
    for (i = 0; i < psr[0].nobs; i++) {
      resids[npar][i] -= y_20;
    }
    rms[npar] = psr[0].rmsPost;
    printf("10cm: offset = %f msec = %f phase\n", 1000*(y_10-y_20), (y_10-y_20)*(double)psr[0].param[param_f].val[0]);
    printf("50cm: offset = %f msec = %f phase\n", 1000*(y_50-y_20), (y_50-y_20)*(double)psr[0].param[param_f].val[0]);
    printf("40cm: offset = %f msec = %f phase\n", 1000*(y_40-y_20), (y_40-y_20)*(double)psr[0].param[param_f].val[0]);
    if(strlen(parfile2[0]) == 0)  /* Only one par file, so quit loop */
      break;
  }
  
  mean_resids /= (float)psr[0].nobs;
  if(strlen(parfile2[0]) == 0)  /* Only one par file */
    npar = 1;
  else
    npar = 2;

  ymin = ymax = resids[0][0];
  xmin = xmax = sats[0][0];
  for(n = 0; n < npar; n++) {
    for (i = 0; i < psr[0].nobs; i++) {
      y = resids[n][i] + resids_err[n][i];
      if(y > ymax)
	ymax = y;
      y = resids[n][i] - resids_err[n][i];
      if(y < ymin)
	ymin = y;
      if(sats[n][i] > xmax)
	xmax = sats[n][i];
      if(sats[n][i] < xmin)
	xmin = sats[n][i];
    }
  }

  if(pgscript) 
    fprintf(fout_pgplot, "pgopen %s\n", PlotDevice);
  else
    deviceID = cpgopen(PlotDevice);
  
  if(deviceID < 0) {
    fprintf(stderr, "Cannot open PGPLOT device\n");
  }
  if(verbose) printf("Opened PGPLOT device '%s'\n", PlotDevice);
  /*  cpgslct(deviceID);  */
  if(windowwidth > 0) {
    x = windowwidth*0.01175548589341692789994673739445429916373;
    y = (windowheight-1)/(float)windowwidth;
    if(pgscript) 
      fprintf(fout_pgplot, "pgpap %e %e\n", x, y);
    else
      cpgpap(x,y); 
  }
  xmin -= 10;
  xmax += 10;
  y = ymax - ymin;
  ymax += 0.05*y;
  ymin -= 0.05*y;
  ymin2 = ymin*(float)psr[0].param[param_f].val[0];
  ymax2 = ymax*(float)psr[0].param[param_f].val[0];
  if(pgscript) 
    fprintf(fout_pgplot, "pgask 0\n");
  else
    cpgask(0);
  if(pgscript) 
    fprintf(fout_pgplot, "pgslw 1\n");
  else
    cpgslw(1);
  if(pgscript) 
    fprintf(fout_pgplot, "pgpage\n");
  else
    cpgpage(); 
  if(pgscript) 
    fprintf(fout_pgplot, "pgsvp 0.1 0.9 0.15 0.9\n");
  else
    cpgsvp(0.1, 0.9, 0.15, 0.9); 
  if(pgscript) 
    fprintf(fout_pgplot, "pgswin %e %e %e %e\n", xmin,xmax,ymin,ymax);
  else
    cpgswin(xmin,xmax,ymin,ymax);
  if(pgscript) 
    fprintf(fout_pgplot, "pgsch 1.5\n");
  else
    cpgsch(1.5); 
  /*  cpgbox("bcnst",0.0,0,"bcnt",0.0,0); */
  if(npar == 2)
    sprintf(txt, "%s: RMS = %.0f us (%.1f times better)", psr[0].name, rms[1], rms[0]/rms[1]);
  else
    sprintf(txt, "%s: RMS = %.0f us", psr[0].name, rms[0]);
  printf("Post RMS parfile 1 = %lf us = %.1lf milli-phase\n", rms[0], 0.001*rms[0]*(float)psr[0].param[param_f].val[0]);
  if(npar == 2)
    printf("Post RMS parfile 2 = %lf us = %.1lf milli-phase\n", rms[1], 0.001*rms[1]*(float)psr[0].param[param_f].val[0]);
  if(pgscript) 
    fprintf(fout_pgplot, "pglab 'MJD - 50000' 'Residual(s)' '%s'\n", txt);
  else
    cpglab("MJD - 50000", "Residual(s)", txt); 
  if(pgscript) 
    fprintf(fout_pgplot, "pgmtxt T 2.0 0.0 0.3 '%s'\n", comment_str);
  else
    cpgmtxt("T", 2.0, 0.0, 0.3, comment_str); 
  /* Left */
  if(pgscript) 
    fprintf(fout_pgplot, "pgaxis 'n' %e %e %e %e %e %e 0.0 0 0.0 0.7 0.3 1.0 180.0\n", xmin, ymax, xmin, ymin,ymax, ymin);
  else
    cpgaxis("n", xmin, ymax, xmin, ymin,ymax, ymin, 0.0, 0, 0.0, 0.7, 0.3, 1.0, 180.0);
  /* Right */
  if(pgscript) 
    fprintf(fout_pgplot, "pgaxis 'n' %e %e %e %e %e %e 0.0 0 0.0 0.7 0.3 1.0 0.0\n", xmax, ymin, xmax, ymax,ymin2, ymax2);
  else
    cpgaxis("n", xmax, ymin, xmax, ymax,ymin2, ymax2, 0.0, 0, 0.0, 0.7, 0.3, 1.0, 0.0);
  /* Bottom */
  if(pgscript) 
    fprintf(fout_pgplot, "pgaxis 'n1' %e %e %e %e %e %e 0.0 0 0.7 0.0 0.3 1.0 0.0\n", xmin, ymin, xmax, ymin,xmin, xmax);
  else
    cpgaxis("n1", xmin, ymin, xmax, ymin,xmin, xmax, 0.0, 0, 0.7, 0.0, 0.3, 1.0, 0.0);
  /* Top */
  /*  cpgaxis("", xmin, ymax, xmax, ymax,xmin, xmax, 0.0, 0, 0, 0.7, 0.3, 1.0, 0.0); */
  if(pgscript) 
    fprintf(fout_pgplot, "pgaxis '' %e %e %e %e %e %e 100000.0 0 0 0.7 0.3 1.0 0.0\n", xmin, ymax, xmax, ymax,xmin, xmax);
  else
    cpgaxis("", xmin, ymax, xmax, ymax,xmin, xmax, 100000.0, 0, 0, 0.7, 0.3, 1.0, 0.0); 

  if(xmax-xmin < 365*0.5)
    plotmonths = 1;
  else if(xmax-xmin < 365)
    plotmonths = 4;
  else if(xmax-xmin < 2*365)
    plotmonths = 6;
  else
    plotmonths = 12;
  day = 1;
  for(yr = 1967; yr < 2050; yr++) {
    for(mt = 1; mt <= 12; mt++) {
      sla_CALDJ(yr, mt, day, &djm, &j);
      if(j != 0) {
	printf("SLA ERROR %d (%d-%d-%d)\n", j, yr,mt,day);
	return 0;
      }
      if(djm >= xmin && djm <= xmax) {
	if(mt%plotmonths == 1) {
	  sprintf(txt, "%s %d", months[mt-1], yr);
	  ticklength = 0.7;
	}else {
	  txt[0] = 0;
	  ticklength = 0.3;
	}
	if(pgscript) 
	  fprintf(fout_pgplot, "pgtick %e %e %e %e 0.5 0.0 %f -0.2 0 '%s'\n", djm-0.5, ymax, djm+0.5, ymax, ticklength, txt);
	else
	  cpgtick(djm-0.5, ymax, djm+0.5, ymax, 0.5, 0.0, ticklength, -0.2, 0, txt);
      }
    }
  }

  if(pgscript) 
    fprintf(fout_pgplot, "pgscr 16 0.3 0.3 0.3\n");
  else
    cpgscr(16,0.3,0.3,0.3);
  if(pgscript) 
    fprintf(fout_pgplot, "pgscr 17 1 1 1\n");
  else
    cpgscr(17,1,1,1);
  if(pgscript) 
    fprintf(fout_pgplot, "pgscr 18 0.3 0.0 0.0\n");
  else
    cpgscr(18,0.3,0.0,0.0);
  if(pgscript) 
    fprintf(fout_pgplot, "pgscr 19 1 0 0\n");
  else
    cpgscr(19,1,0,0);
  if(pgscript) 
    fprintf(fout_pgplot, "pgscr 20 0.0 0.3 0.0\n");
  else
    cpgscr(20,0.0,0.3,0.0);
  if(pgscript) 
    fprintf(fout_pgplot, "pgscr 21 0 1 0\n");
  else
    cpgscr(21,0,1,0);
  if(pgscript) 
    fprintf(fout_pgplot, "pgscr 22 0.0 0.0 0.3\n");
  else
    cpgscr(22,0.0,0.0,0.3);
  if(pgscript) 
    fprintf(fout_pgplot, "pgscr 23 0.5 0.5 1\n");
  else
    cpgscr(23,0.5,0.5,1);
  if(pgscript) 
    fprintf(fout_pgplot, "pgsch 3\n");
  else
    cpgsch(3);
  if(pgscript) 
    fprintf(fout_pgplot, "pgslw 4\n");
  else
    cpgslw(4);
  for(n = 0; n < npar; n++) {
    /* Plot glitch epochs */
    if(n == 0) {   /* Somehow getParamer gives segmentation fault if n!=0 */
      for(j = 0; j < psr[n].param[param_glep].aSize; j++) {
	/*	sprintf(txt, "GLEP_%d", j+1); */
	glitch_epoch = psr[n].param[param_glep].val[j];
	if(glitch_epoch > 1) {
	  printf("Found glitch at MJD %lf\n", glitch_epoch);
	  if(pgscript) 
	    fprintf(fout_pgplot, "pgsci 6\n");
	  else
	    cpgsci(6);
	  if(pgscript) 
	    fprintf(fout_pgplot, "pgmove %e %e\n", glitch_epoch, ymin);
	  else
	    cpgmove(glitch_epoch, ymin);
	  if(pgscript) 
	    fprintf(fout_pgplot, "pgdraw %e %e\n", glitch_epoch, ymax);
	  else
	    cpgdraw(glitch_epoch, ymax);
	  if(pgscript) 
	    fprintf(fout_pgplot, "pgsci 1\n");
	  else
	    cpgsci(1);
	}else {
	  if(verbose) printf("No GLEP_%d in par file\n", j+1);
	  break;
	}
      }
      /* Plot phase jumps */
      for(j = 0; j < psr[n].nPhaseJump; j++) {
	if(pgscript) 
	  fprintf(fout_pgplot, "pgsci 7\n");
	else
	  cpgsci(7);
	if(pgscript) 
	  fprintf(fout_pgplot, "pgmove %e %e\n", psr[n].phaseJump[j], ymin);
	else
	  cpgmove(psr[n].phaseJump[j], ymin);
	if(pgscript) 
	  fprintf(fout_pgplot, "pgdraw %e %e\n", psr[n].phaseJump[j], ymax);
	else
	  cpgdraw(psr[n].phaseJump[j], ymax);
	if(pgscript) 
	  fprintf(fout_pgplot, "pgsci 1\n");
	else
	  cpgsci(1);
      }
    }
    /* Plot fitwaves */
    if(n == 0) {   /* Somehow getParamer gives segmentation fault if n!=0 */
      if(plot_fitwaves) {
	if(getParameter(psr[n], "WAVE_OM", &fitwaves_omega) == 1) {
	  if(fitwaves_omega > 1e-8) {
	    printf("WAVE_OM = %lf\n", fitwaves_omega);
	    for(j = 0; j < psr[n].nWhite; j++) {
	      printf("WAVE%d amp sin=%lf amp cos=%lf\n", j+1, psr[n].wave_sine[j], psr[n].wave_cos[j]);
	    }	
	    if(pgscript) 
	      fprintf(fout_pgplot, "pgsci 6\n");
	    else
	      cpgsci(6);
	    getParameter(psr[n], "PEPOCH", &fitwaves_epoch);
	    fitwaves_omega *= 1;
	    mean2 = 0;
	    for (i = 0; i < psr[n].nobs; i++) {
	      mean2 += fitwave_function(&psr[n], sats[n][i], fitwaves_omega, fitwaves_epoch);
	    }
	    mean2 /= (double)psr[n].nobs;
	    for(k = 0; k <= 1000; k++) {
	      x = xmin+(xmax-xmin)*k/1000.0;
	      y = fitwave_function(&psr[n], x, fitwaves_omega, fitwaves_epoch) + mean_resids - mean2;
	      if(k == 0)
		if(pgscript) 
		  fprintf(fout_pgplot, "pgmove %e %e\n", x, y);
		else
		  cpgmove(x, y);
	      else	
		if(pgscript) 
		  fprintf(fout_pgplot, "pgdraw %e %e\n", x, y);
		else
		  cpgdraw(x, y);
	    }
	  }
	  if(pgscript) 
	    fprintf(fout_pgplot, "pgsci 1\n");
	  else
	    cpgsci(1);
	}
      }
    }

    for (i = 0; i < psr[0].nobs; i++) {
      /*            printf("%f\n", resids_freq[n][i]);  */
      if(colorcode == 0) {
        j = 16+n;
      }else {
	if(resids_freq[n][i] > 600 && resids_freq[n][i] < 700) {
	  j = 18+n;
	}
	else if(resids_freq[n][i] > 700 && resids_freq[n][i] < 800)
	  {
	    // RMS make 40 cm data this colour
	    j= 19+n;
	  }
	else if(resids_freq[n][i] > 1300 && resids_freq[n][i] < 1500) {
	  j = 16+n;
	}else if(resids_freq[n][i] > 3000 && resids_freq[n][i] < 3300) {
	  j = 22+n;
	}else{
          j = 20+n;
	}
      }
      if(npar == 1)
	j += 1;
      if(pgscript) 
	fprintf(fout_pgplot, "pgsci %d\n", j);
      else
	cpgsci(j);
      if(resids_backend[n][i] == 'r') {
	if(pgscript) 
	  fprintf(fout_pgplot, "pgpt1 %e %e 3\n", sats[n][i], resids[n][i]);
	else
	  cpgpt1(sats[n][i], resids[n][i], 3);
      }else if(resids_backend[n][i] == 's') {
	if(pgscript) 
	  fprintf(fout_pgplot, "pgpt1 %e %e 4\n", sats[n][i], resids[n][i]);
	else
	  cpgpt1(sats[n][i], resids[n][i], 4);
      }
      if(pgscript) 
	fprintf(fout_pgplot, "pgerr1 6 %e %e %e 1\n", sats[n][i], resids[n][i], resids_err[n][i]);
      else
	cpgerr1(6, sats[n][i], resids[n][i], resids_err[n][i], 1);
    }
    if(outputNewPar == 1) {
      if(n == 1 || (n == 0 && strlen(parfile2[0]) == 0)) {
	printf("Written out newpar.par");
	textOutput(psr,*npsr,0,0,0,1,"newpar.par");
      }
    }
  }

  /*
  if(pgscript) 
    fprintf(fout_pgplot, "pg\n");
  else

blaat
  */

  if(pgscript) 
    fprintf(fout_pgplot, "pgsch 1\n");
  else
    cpgsch(1);
  free(sats[0]);
  free(resids[0]);
  free(resids_err[0]);
  free(resids_freq[0]);
  free(resids_backend[0]);
  if(npar > 1) {
    free(sats[1]);
    free(resids[1]);
    free(resids_err[1]);
    free(resids_freq[1]);
    free(resids_backend[1]);
  }
  if(pgscript) 
    fprintf(fout_pgplot, "pgend\n");
  else
    cpgend();
  if(pgscript) {
    fclose(fout_pgplot);
  }
  if(verbose) printf("I did properly end.\n");
  return 0;
}

/*
     - - - - - -
      C A L D J
     - - - - - -

  Gregorian Calendar to Modified Julian Date

  (Includes century default feature:  use sla_CLDJ for years
   before 100AD.)

  Given:
     IY,IM,ID     int    year, month, day in Gregorian calendar

  Returned:
     DJM          dp     modified Julian Date (JD-2400000.5) for 0 hrs
     J            int    status:
                           0 = OK
                           1 = bad year   (MJD not computed)
                           2 = bad month  (MJD not computed)
                           3 = bad day    (MJD computed)

  Acceptable years are 00-49, interpreted as 2000-2049,
                       50-99,     "       "  1950-1999,
                       100 upwards, interpreted literally.

  Called:  sla_CLDJ

  P.T.Wallace   Starlink   November 1985
*/

void sla_CALDJ(int IY, int IM, int ID, double *DJM, int *J)
{
  int NY;
  /*  Default century if appropriate */
  if(IY >= 0 && IY < 49) {
    NY=IY+2000;
  }else if(IY>=50 && IY <99) {
    NY=IY+1900;
  }else {
    NY=IY;
  }
  /*  Modified Julian Date */
  sla_CLDJ(NY,IM,ID,DJM,J);
}


/*
     - - - - -
      C L D J
     - - - - -

  Gregorian Calendar to Modified Julian Date

  Given:
     IY,IM,ID     int    year, month, day in Gregorian calendar

  Returned:
     DJM          dp     modified Julian Date (JD-2400000.5) for 0 hrs
     J            int    status:
                           0 = OK
                           1 = bad year   (MJD not computed)
                           2 = bad month  (MJD not computed)
                           3 = bad day    (MJD computed)

  The year must be -4699 (i.e. 4700BC) or later.

  The algorithm is derived from that of Hatcher 1984
  (QJRAS 25, 53-55).

  P.T.Wallace   Starlink   December 1985
-
*/
void sla_CLDJ (int IY, int IM, int ID, double *DJM, int *J) 
{

  /*  Month lengths in days */
  int MTAB[12] = {31,28,31,30,31,30,31,31,30,31,30,31};


  /*  Preset status */
  *J=0;

  /*  Validate year */
  if (IY<-4699)
    *J=1;
  else {
    /*     Validate month */
    if (IM>=1 && IM<=12) {
      /*        Allow for leap year */
      if ((IY % 4) == 0) 
	MTAB[1]=29;
      else
	MTAB[1]=28;
      if ((IY % 100) == 0 && (IY % 400) != 0)
	MTAB[1]=28;

      /*        Validate day */
      if (ID < 1 || ID > MTAB[IM-1]) 
	*J=3;

      /*        Modified Julian Date */
      *DJM=((1461*(IY-(12-IM)/10+4712))/4
	    +(306*((IM+9)%12)+5)/10
	    -(3*((IY-(12-IM)/10+4900)/100))/4
	    +ID-2399904);

      /*        Bad month */
    }else {
      *J=2;
    }
  }  /* End else year */
}

int getParameter(pulsar psr,char *param, double *value)
{
  int i;
  int gotit=-1;
  for (i=0;i<MAX_PARAMS;i++) {
    if (strcasecmp(psr.param[i].shortlabel[0],param)==0)
      gotit=i;
  }
  if (gotit != -1) {
    *value = *(psr.param[gotit].val);
    return 1;
  }else
    return 0;
}

float fitwave_function(pulsar *psr, float x, float fitwaves_omega, float fitwaves_epoch)
{
  int j;
  float y;
  y = 0;
  for(j = 0; j < psr[0].nWhite; j++) {
    y += psr[0].wave_sine[j]*sin(fitwaves_omega*(double)(j+1)*(x-fitwaves_epoch));
    y += psr[0].wave_cos[j]*cos(fitwaves_omega*(double)(j+1)*(x-fitwaves_epoch));
  }
  y *= -1;
  return y;
}

/* Patrick: taken from NR */
#define NRANSI
/*#include "glast_nrutil.h"*/
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50
#define NR_END 1
#define FREE_ARG char*
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}
int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}
void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
void indexx_patrick(unsigned long n, float arr[], unsigned long indx[])
{
	unsigned long i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	float a;

	istack=ivector(1,NSTACK);
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=l;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[l+1]]) {
				SWAP(indx[l],indx[l+1])
			}
			i=l+1;
			j=ir;
			indxt=indx[l+1];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l+1]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in indexx.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI
