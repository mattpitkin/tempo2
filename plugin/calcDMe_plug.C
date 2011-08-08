//  Copyright (C) 2006,2007,2008,2009, Stefan Oslowski, George Hobbs, Russel Edwards

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

/* This plugin was written by Stefan Oslowski (CAS, Swinburne) 
 * and George Hobbs (ATNF, CSIRO)
 * Feel free to use and modify the code but do not remove this information
 * (that includes all the comments until first #include statement)
 * and give credit when you use the results obtained with this code
 * in a publication. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"
#include "T2toolkit.h"
#include <time.h>
#include <cpgplot.h>
#include "TKspectrum.h"

/* TODO:
 * identifying jumps when there's no FJUMP -f in the parfile
 */

char parFile[23][MAX_FILELEN], timFile[23][MAX_FILELEN];
char *outFileName;
//xObs: fit - freqency used for fit, dm - frequency used to get the dm variations, imp -frequency to be improved
int fitObs[MAX_OBSN_VAL], dmObs[MAX_OBSN_VAL], impObs[MAX_OBSN_VAL];
//xCount: counts of observations of given type (see above)
int fitCount = 0, dmCount = 0, impCount = 0;
//auxilary array holding frequencies: (used only in init and setFreq*)
//frequencies are: 0: fitFreq 1: dmFreq 2: fitFreqWidth 3:dmFreqWidth
//4: fitFreqMin 5: fitFreqMax 6: dmFreqMin 7: dmFreqMax
double freqArray[8];
//control variable for setting frequencies above:
int freqOffset[2]; //accepts values 0, 1, 2 and 3
//if freqOffset[i] == 3, then its a flag identifying observations at given freq:
char freq1f[MAX_STRLEN], freq2f[MAX_STRLEN];
//DM after first callFit (which fits after choosing just 10cm data) and it's error
long double dm0 = -1.0, dm0_err = -1.0;
//same for F0:
long double f0_0 = -1.0, f0_0_err = -1.0;
//size of the bin:
long double binSizeDays = 14;
//variables needed to use only appropriate jumps
int nf = 0;
int valID[MAX_FLAGS];
//should f0 be fitted? yes by default
int f0fit = 1;
//array holding ids of fitObs and dmObs in given 10day bi
int binObs[MAX_OBSN_VAL][2];
//counts of obs in bins and the last change
//IMPORTANT: bin_fitCount must not equal bin_fitCount_inc at the beginning!!
//(otherwise the fitting won't start)
int bin_fitCount = 0, bin_dmCount = 0;
int bin_fitCount_inc = -1, bin_dmCount_inc = -1;
//starting point of bin
long double binStart = 1e10;
//arrays holding the calculated delta dm, its error and corresponding mjd
long double ddm[MAX_OBSN_VAL], ddmMJD[MAX_OBSN_VAL], ddmErr[MAX_OBSN_VAL];
int ddmCount = 0;
//arrays for the daily ddm variations (both interpolated and smoothed)
//we obtain them by interpolating ddm* and resampling daily and smoothing
double outX[2 * MAX_OBSN_VAL], outY[2 * MAX_OBSN_VAL];
int outInterpCount = -1, outSmoothCount = -1;
//smoothing width:
int smoothWidth = 100;
//output as ASCII?
int ascii = 1;
//auxilary variable for output file name
int gotOut = 0;
//output header?
int header = 0;
//output full DM instead of its variations?
int outDM = 0;
//remove means from output?
int mean = 0, meanMJD = 0;
double meanMJDval = 0, meanVal = 0;
//display results on screen? x/y labels, title
int doDisplay = 0;
char *xlab, *ylab, *title;
//create a postscript?
int hardcopy = 0;
char *gr;
//output spline or raw?
int splineOut = 0, rawOut = 0;
//how many days separate an observing sessions?
double sessionSeparation = 3;
//sessions starts and ends:
double *start_sessions;
double *finish_sessions;
//which session was last used? if none or we don't use them directly, then should equal -1
int lastUsedSession = -1;
//how many session were found?
int nSessions = 0;
//write all the par and tim files for each fit?
int allParTim = 0;

using namespace std; /* Is this required for a plugin ? Yes, for linux */

char dcmFile[MAX_FILELEN];

void callFit(pulsar *psr, int npsr);
void get_binObs(pulsar *psr);
void setFitParams(pulsar *psr);
void setAllDeleted(pulsar *psr);
void resetDMandF0(pulsar *psr);
void findFirst(pulsar *psr);
void interpolateSplineSmooth();
void interpolateWeightedSmooth();
void init(int argc, char* argv[]);
void handleFreqPoints(pulsar *psr);
void findSessions(pulsar *psr);
void describe();
double findMean(double *x, int count);
void display(char *gr, int publish, double *xx, double *yy, long double *ddmMJD, long double *ddm, long double *ddmErr, int outInterpCount, int outSmoothCount, int ddmCount, char *xlab, char *ylab, char *title, double meanMJDval, double meanVal);
void output(char *outFileName, int ascii, double dm0, int header, int outDM, double *outX, double *outY, int outInterpCount, int outSmoothCount, int mean, int meanMJD, double *meanMJDval, double *meanVal, int splineOut, int rawOut, long double *ddmMJD, long double *ddm, long double *ddmErr, int ddmCount);

void help() /* Display help */ {
  printf("\n\npress ENTER to continue\n");
  scanf("%*c");
  printf("\nOptions:\n");
  printf("===========\n");
  printf("\n-f par_file tim_file\n\n   provide par and tim file\n   (tim file can also be given as the last argument)\n");
  printf("FREQUENCY CHOICE:\n===========\n");
  printf("\n-freq1 central_freq freq_width\n\n   provide first frequency used to calculate dm variations\n");
  printf("   central_freq is the frequency one wants to use (MHz)\n   freq_width is the acceptable offset (MHz)\n");
  printf("   anything in range (central_freq - freq_width; central_freq + freq_width)\n   will be included as freq1\n");
  printf("   ALTERNATIVE SYNTAX: -freq1 alias\n       where alias is either 10cm, 20cm or 50cm\n");
  printf("   DEFAULT: central_freq = 3100MHz (10cm) freq_width = 100MHz\n");
  printf("\n-freq1r freq_min freq_max\n\n   provide the range of first frequency used to calculate dm variations\n");
  printf("   anything within range (freq_min (MHz):freq_max (MHz))\n   will be included as freq1\n");
  printf("\n-freq1f flag\n\n   provide a flag identifying observations at given frequency\n");
  printf("   flag could for example be \"-F 10cm\". Note that the flag includes the flag identifier as well\n");
  printf("   if you use this flag, you have to use -freq2f as well\n");
  printf("\n-freq2 central_freq freq_width\n\n   provide second frequency used to calculate dm variations\n");
  printf("   central_freq is the frequency one wants to use (MHz)\n   freq_width is the acceptable offset (MHz)\n");
  printf("   anything in range (central_freq - freq_width; central_freq + freq_width)\n   will be included as freq2\n");
  printf("   ALTERNATIVE SYNTAX: -freq2 alias\n       where alias is either 10cm, 20cm or 50cm\n");
  printf("   DEFAULT: central_freq = 685MHz (50cm) freq_width = 100MHz\n");
  printf("\n-freq2r freq_min freq_max\n\n   provide the range of second frequency used to calculate dm variations\n");
  printf("   anything within range (freq_min (MHz):freq_max (MHz))\n   will be included as freq2\n");
  printf("\n-freq2f flag\n\n   provide a flag identifying observations at given frequency\n");
  printf("   flag could for example be \"-F 50cm\". Note that the flag includes the flag identifier as well\n");
  printf("   if you use this flag, you have to use -freq1f as well\n");
  printf("\nFITTING:\n===========\n");
  printf("\n-bs bin_size\n\n   defines width of bin in which \n   dm variations are measured (in days)\n");
  printf("   if bin_size < 0 then single observing sessions (see -ssep) will be used\n");
  printf("   DEFAULT: bin_size = 14days\n");
  printf("\n-ssep session_separation\n\n   defines the time-span of observing sessions (in days). This will affect the binning of data, if \n");
  printf("   session_separation > 0 then the observations from given sessions will not be split into separate binsn \n");
  printf("   DEFAULT: session_separation = 3days\n");
  printf("\n-sw smoothingWidth\n\n   defines smoothing width (days > 1)\n");
  printf("   DEFAULT: smoothingWidth=100days\n");
  printf("\n-noffit\n\n   don't fit for pulsar frequency\n");
  printf("\n-details\n\n   if this flag is provided, the plugin will give details\n   of fit in each of the bins\n");
  printf("\nOUTPUT\n===========\n");
  printf("\n-display\n\n   show a plot on screen of data points (green),\n   interpolated (red, dashed) and smoothed (white or black, solid) deltaDMs / DMs\n");
  printf("\n-hardcopy\n\n   produce a ps version of the plot described in '-display'\n");
  printf("   it will be stored in a file dm_PSR.ps, where PSR is the name of the pulsar\n");
  printf("   postscript version of the plot doesn't contain copyrights note\n");
  printf("\n-binary\n\n   output files are binary files\n");
  printf("\n-o out_file\n\n   output file name \n   DEFAULT: dm_smoothPSR.out, where PSR is the pulsar id)\n");
  printf("   FORMAT:\n   next nSmooth lines: MJD deltaDM\n");
  printf("   where all the deltaDM are given with respect to dm0\n   and MJDs are at the start of bin\n");
  printf("   WARNING: this will overwrite files without asking\n");
  printf("\n-header\n\n   output a header in the output file: one line containing dm0 and number of output lines\n");
  printf("\n-outDM\n\n   output DM instead of deltaDM\n");
  printf("\n-mean\n\n   remove mean from deltaDM (or DM)\n");
  printf("\n-meanMJD\n\n   remove mean from MJD\n");
  printf("\n-splineout\n\n   this flag turns on output of dm variations before they are smoothed\n   the output file is same as set by -o but with '.int' extension\n   this file is always ascii\n");
  printf("\n-rawout\n\n   this flag turns on output of dm variations before they are interpolated\n   the output file is same as set by -o but with '.raw' extension\n   this file is always plain ascii\n");
  printf("   it will contain three columns: MJD deltaDM deltaDM_error, where deltaDM can be also DM (see -outDM)\n");
  printf("\n-allout\n\n   output par and tim file for each fit (one per bin) + after the initial fit to higher frequency\n");
  printf("   the file names will be n.par and n.tim for n-th session\n"); 
  printf("\n Example of usage:  tempo2 -gr calcDMe -bs 50 -sw 100 -f J0437-4715.par J0437-4715.tim -outDM -header\n");
  printf("   this will calculate dm variations based on 10 and 50cm data and output files suitable for use with \n");
  printf("   tempo2 -setdm option\n");
  printf("\n KNOWN PROBLEMS: \n   if for any reason the plugin fails to bin the data properly, it is likely \n");
  printf("   it will fall into an infinite loop. Possible solution: check the start of the current bin and if it\n");
  printf("   exceeds the last TOA just quit - this would require the data to be sorted though, perhaps there's a better way?\n");
  exit(0);
} //help()

void describe() /* display description*/ {
  printf("                        calcDMe\n");
  printf("This plugin was written by Stefan Oslowski (CAS, Swinburne)\n");
  printf("and George Hobbs (ATNF, CSIRO).\n");
  printf("Feel free to use and modify the code but do not remove this information\n");
  printf("and give credit when you use the results obtained with this code in a \n");
  printf("publication.\n\n");
  printf("calcDMe calculates DM variations using pulsar TOAs at two different\n");
  printf("frequencies (see -freq1 and -freq2 in help). This DM variations are\n");
  printf("interpolated, resampled at 1 day period and smoothed at provided\n");
  printf("smoothing period (see -sw in help). The results are written into output\n");
  printf("file (see -o in help) and can be used with tempo2 -setdm option)\n");
  printf("to improve timing residuals or timing parameters.\n");
  printf("Any suggestions which flags should be on by default are most welcome.\n\n");
  printf("calcDMe plugin for tempo2 Copyright (C) 2009 Stefan Oslowski, George Hobbs\n");
} //describe()

/* The main function called from the TEMPO2 package is 'graphicalInterface' */

/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc, char *argv[], pulsar *psr, int *npsr) {
  int i, k;

  //variables for displaying / hardcopying results:
  outFileName = (char*) malloc(500 * sizeof(char));
  gr = (char*) malloc(100 * sizeof(char));
  xlab = (char*) malloc(100 * sizeof(char));
  ylab = (char*) malloc(100 * sizeof(char));
  title = (char*) malloc(100 * sizeof(char));

  *npsr = 1; /* This graphical interface will only show results for one pulsar */

  //initilize freqArray:
  for (i = 0; i < 8; ++i) {
    freqArray[i] = -1;
  }

  //handle command line arguments:
  if (debugFlag == 1) printf("calcDMe: calling init()\n");
  init(argc, argv);
  if (debugFlag == 1) printf("calcDMe: calling readParfile\n");
  readParfile(psr, parFile, timFile, *npsr); /* Load the parameters       */
  if (debugFlag == 1) printf("calcDMe: calling readTimfile\n");
  readTimfile(psr, timFile, *npsr); /* Load the arrival times    */
  if (debugFlag == 1) printf("calcDMe: calling preProcess %d\n", psr[0].nobs);
  preProcess(psr, *npsr, argc, argv);

  //set the output file name to default, if user didn't provide:
  if (debugFlag == 1) printf("checking if outFileName was set...\n");
  if (gotOut == 0) {
    outFileName = strcpy(outFileName, "dm_smooth");
    outFileName = strcat(outFileName, psr[0].name);
    outFileName = strcat(outFileName, ".out");
  }

  if (debugFlag == 1) printf("outFileName = %s\n", outFileName);

  //find freq1 and freq2 frequency data to use for calculating deltaDMs
  handleFreqPoints(psr);
  if (debugFlag >= 1) {
    printf("found %d fitObs and %d dmObs\n", fitCount, dmCount);
    printf("found jumps:\n");
    for (i = 0; i <= psr[0].nJumps; ++i) {
      printf("%s\n", psr[0].jumpStr[i]);
    }
  }

  //do the initial fit
  //make sure we are not fitting for dm:
  // (otherwise jumps become all wrong
  psr[0].param[param_dm].fitFlag[0] = 0;
  if (debugFlag == 1) printf("calcDMe: calling callFit %d\n", psr[0].nobs);
  callFit(psr, *npsr);
  if (allParTim == 1) {
    textOutput(psr,*npsr,0,0,0,1,"afterFit.par");
    writeTim("afterFit.tim",psr,"tempo2");
  }
  /* Convert all prefit values to postfit values */
  for (i = 0; i < MAX_PARAMS; i++) {
    for (k = 0; k < psr->param[i].aSize; k++) {
      if (psr->param[i].paramSet[k] == 1) {
	psr->param[i].prefit[k] = psr->param[i].val[k];
	psr->param[i].prefitErr[k] = psr->param[i].prefitErr[k];
      }
    }
  }

  //exit(0);
  //store initial dm0 and f0:
  dm0 = psr[0].param[param_dm].val[0];
  dm0_err = psr[0].param[param_dm].err[0];
  f0_0 = psr[0].param[param_f].val[0];
  f0_0_err = psr[0].param[param_f].err[0];

  //locate the first observation from the fitObs - that's the beginning of our first point
  findFirst(psr);
  //mark all observations as deleted
  setAllDeleted(psr);
  //now we have to fit for DM (and F0 possibly), also turn off jump-fit
  setFitParams(psr);
  //determine observing sessions
  if (sessionSeparation > 0 ) 
  {
    start_sessions = (double*) malloc((fitCount + dmCount)*sizeof(double));
    finish_sessions = (double*) malloc((fitCount + dmCount)*sizeof(double));
    if (debugFlag==1) printf("allocated memory for start and finish\n");
    findSessions(psr);
    if (debugFlag==1) printf("found sessions: %d points in %d sessions\n\n",fitCount+dmCount,nSessions);
  }

  //start the actual binning and fitting:
    while (fitCount != bin_fitCount && lastUsedSession < nSessions) { 
    //find all the points from fit- and dmObs within this bin
    if (debugFlag==1) printf("trying to find a bin...\n");
    get_binObs(psr);
    //ensure we've got enough points for fitting:
    if (bin_dmCount_inc > 0 && bin_fitCount_inc > 0 && (bin_fitCount_inc + bin_dmCount_inc) > 1 + f0fit) {
      //do the fit
      if (debugFlag==1) printf("found enough points, trying to fit!\n");
      if (debugFlag==1) printf("using session %d\n",lastUsedSession-1);
      if (debugFlag==1) printf("Calling fit ------\n");
      callFit(psr, *npsr);
      if (allParTim == 1) {
	char *tmpstring;
	tmpstring = (char*) malloc(100*sizeof(char));
	sprintf(tmpstring,"%d.par",lastUsedSession);
	textOutput(psr,*npsr,0,0,0,1,tmpstring);
	free(tmpstring);
	tmpstring = (char*) malloc(100*sizeof(char));
	sprintf(tmpstring,"%d.tim",lastUsedSession);
	writeTim(tmpstring,psr,"tempo2");
	free(tmpstring);
      }
      ddm[ddmCount] = psr[0].param[param_dm].val[0] - dm0;
      ddmErr[ddmCount] = psr[0].param[param_dm].err[0];
      ddmMJD[ddmCount] = binStart;
      ddmCount++;
      //reset dm and f0: not necessary, actually...
      //            resetDMandF0(psr);
    } else if (debugFlag==1) {
      printf("==============================\nNOT ENOUGH POINTS FOR FITTING IN THE BIN\n");
      printf("STARTING AT %Lf (%d freq1 and %d freq2 points)   \n==============================\n", binStart, bin_fitCount_inc, bin_dmCount_inc);
      lastUsedSession++;
    }
    if (debugFlag==1) printf("moving to next bin\n");
    binStart += binSizeDays;
  }
  //interpolate the delta DMs and smooth them afterwards
  /* if (debugFlag == 1) printf("calling interpolateSplineSmooth()\n");
     interpolateSplineSmooth();*/
  if (debugFlag == 1) printf("calling interpolateWeightedSmooth()\n");
  interpolateWeightedSmooth();
//just a check, remove it:
//TODO
k=0;
for (i=0;i<fitCount;i++)
if (psr[0].obsn[fitObs[i]].sat < start_sessions[0] || psr[0].obsn[fitObs[i]].sat > finish_sessions[89])
k++;
  //write output to a file:
  if (debugFlag == 1) printf("calling output()\n");
  output(outFileName, ascii, (double) dm0, header, outDM, outX, outY, outInterpCount, outSmoothCount, mean, meanMJD, &meanMJDval, &meanVal, splineOut, rawOut, ddmMJD, ddm, ddmErr, ddmCount);

  //display result on screen, when needed
  if (debugFlag == 1) printf("calling display()\n");
  if (doDisplay == 1) {
    if (hardcopy == 1) {
      gr = strcpy(gr, "dm_");
      gr = strcat(gr, psr[0].name);
      gr = strcat(gr, ".ps/CPS");
    } else
      gr = strcpy(gr, "/XS");
    if (outDM == 1) {
      ylab = strcpy(ylab, "DM (cm\\u-3\\dpc)");
    } else
      ylab = strcpy(ylab, "\\gDDM(cm\\u-3\\dpc)");
    xlab = strcpy(xlab, "MJD");
    if (meanMJD == 1) {
      xlab = strcat(xlab, " - ");
      char *tmp;
      tmp = (char*) malloc(20 * sizeof(char));
      sprintf(tmp, "%.2lf", meanMJDval);
      xlab = strcat(xlab, tmp);
      free(tmp);
    }

    display(gr, 0, outX, outY, ddmMJD, ddm, ddmErr, outInterpCount, outSmoothCount, ddmCount, xlab, ylab, psr[0].name, meanMJDval, meanVal);
  }
  free(outFileName);
  free(gr);
  free(xlab);
  free(ylab);
  free(title);

  return 0;
} //graphicalInterface()

//handle plugin initilization

void init(int argc, char* argv[]) {
  int i, k, gotTim = 0;
  strcpy(dcmFile, "NULL");

  printf("calcDMe plugin\n");
  printf("Author(s):           Stefan Oslowski, George Hobbs\n");
  printf("Contact: soslowski@astro.swin.edu.au\n");
  printf("Version:             0.4.1\n");
  printf(" --- type '-h' or '--help' for help information\n\n");

  describe();


  //initialize freqOffset - this variables control in what way we identify points at given frequency, related to options -freq:
  freqOffset[0] = freqOffset[1] = 0;
  /* Obtain the .par and the .tim file from the command line */
  if (argc < 4) {
    printf("\n not enough arguments provided\n you have to provide at least .tim file\n displaying help:\n\n");
    help();
  } else if (argc == 4 && strcmp(argv[3], "-h") != 0) /* Only provided .tim name */ {
    strcpy(timFile[0], argv[3]);
    strcpy(parFile[0], argv[3]);
    parFile[0][strlen(parFile[0]) - 3] = '\0';
    strcat(parFile[0], "par");
  } else if (argc == 4 && strcmp(argv[3], "-h") == 0) {
    help();
  }    
  //handle the command line arguments:
  for (i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-bs") == 0) {
      sscanf(argv[++i], "%Lf", &binSizeDays);
    } else if (strcmp(argv[i], "-sw") == 0)
      sscanf(argv[++i], "%d", &smoothWidth);
    else if (strcmp(argv[i], "-ssep") == 0)
      sscanf(argv[++i], "%lf",&sessionSeparation);
    else if (strcmp(argv[i], "-noffit") == 0)
      f0fit = 0;
    else if (strcmp(argv[i], "-binary") == 0)
      ascii = 0;
    else if (strcmp(argv[i], "-details") == 0)
      if (debugFlag == 0) {
	debugFlag = 2;
      } else
	debugFlag = 1;
      else if (strcmp(argv[i], "-o") == 0) {
	if (argv[++i][0] != '-') {
	  strcpy(outFileName, argv[i]);
	  gotOut = 1;
	} else {
	  printf("wrong usage, displaying help:\n\n\n");
	  help();
	}
      } else if (strcmp(argv[i], "-display") == 0) {
	doDisplay = 1;
      } else if (strcmp(argv[i], "-hardcopy") == 0) {
	hardcopy = 1;
	doDisplay = 1;
      } else if (strcmp(argv[i], "-header") == 0) {
	header = 1;
      } else if (strcmp(argv[i], "-outDM") == 0) {
	outDM = 1;
      } else if (strcmp(argv[i], "-mean") == 0) {
	mean = 1;
      } else if (strcmp(argv[i], "-meanMJD") == 0) {
	meanMJD = 1;
      } else if (strcmp(argv[i], "-f") == 0) {
	strcpy(parFile[0], argv[++i]);
	if (argv[i + 1][0] != '-') {
	  strcpy(timFile[0], argv[++i]);
	  gotTim = 1;
	}
      } else if (strcmp(argv[i], "-splineout") == 0) {
	splineOut = 1;
      } else if (strcmp(argv[i], "-rawout") == 0) {
	rawOut = 1;
      } else if (strcmp(argv[i], "-allout") == 0) {
	allParTim = 1;
      } else if (strcmp(argv[i], "-freq1f") == 0) {
	strncpy(freq1f, argv[++i], MAX_STRLEN);
	strcat(freq1f, " ");
	freqOffset[0] = 3;
      } else if (strcmp(argv[i], "-freq2f") == 0) {
	strncpy(freq2f, argv[++i], MAX_STRLEN);
	strcat(freq2f, " ");
	freqOffset[1] = 3;
      } else if (strcmp(argv[i], "-freq1") == 0) {
	k = i + 1;
	if (strcmp(argv[k], "10cm") == 0) {
	  freqArray[0] = 3100.0;
	  freqArray[2] = 100.0;
	  i++;
	} else if (strcmp(argv[k], "20cm") == 0) {
	  freqArray[0] = 1410.0;
	  freqArray[2] = 100.0;
	  i++;
	} else if (strcmp(argv[k], "50cm") == 0) {
	  freqArray[0] = 685.0;
	  freqArray[2] = 100.0;
	  i++;
	} else {
	  sscanf(argv[++i], "%lf", &freqArray[0]);
	  sscanf(argv[++i], "%lf", &freqArray[2]);
	}
	freqOffset[0] = 1;
      } else if (strcmp(argv[i], "-freq2") == 0) {
	k = i + 1;
	if (strcmp(argv[k], "10cm") == 0) {
	  freqArray[1] = 3100.0;
	  freqArray[3] = 100.0;
	  i++;
	} else if (strcmp(argv[k], "20cm") == 0) {
	  freqArray[1] = 1410.0;
	  freqArray[3] = 100.0;
	  i++;
	} else if (strcmp(argv[k], "50cm") == 0) {
	  freqArray[1] = 685.0;
	  freqArray[3] = 100.0;
	  i++;
	} else {
	  sscanf(argv[++i], "%lf", &freqArray[1]);
	  sscanf(argv[++i], "%lf", &freqArray[3]);
	}
	freqOffset[1] = 1;
      } else if (strcmp(argv[i], "-freq1r") == 0) {
	sscanf(argv[++i], "%lf", &freqArray[4]);
	sscanf(argv[++i], "%lf", &freqArray[5]);
	freqOffset[0] = 2;
      } else if (strcmp(argv[i], "-freq2r") == 0) {
	sscanf(argv[++i], "%lf", &freqArray[6]);
	sscanf(argv[++i], "%lf", &freqArray[7]);
	freqOffset[1] = 2;
      } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0 || strcmp(argv[i],"-help") == 0) {
	help();
      } else if (i == argc - 1 && gotTim == 0) {
	strcpy(timFile[0], argv[i]);
	gotTim = 1;
      }
  }
  //if -freq1f was provided, then -freq2f has to be provided as well
  //(or vice versa)
  if (freqOffset[0] != 3 && freqOffset[1] == 3) {
    printf("\n\n   ERROR: only -freq1f was provided, you must use also -freq2f!!!  \n\n");
    exit(-1);
  }
  if (freqOffset[1] != 3 && freqOffset[0] == 3) {
    printf("\n\n   ERROR: only -freq2f was provided, you must use also -freq1f!!!  \n\n");
    exit(-1);
  }
} //init

void handleFreqPoints(pulsar *psr) {
  int i, j, k, l, found, foundFlag;
  //which frequencies are to be used as 'fit' and 'dm' data
  double fitFreq, dmFreq;
  //alternatively, their accepted range:
  long double fitFreqMin = 0, fitFreqMax = 0;
  long double dmFreqMin = 0, dmFreqMax = 0;
  char _tmp[100];

  //ensure freq1 > freq2 (i.e. use 10cm as fit and 50cm as dm)
  if (freqArray[0] < freqArray[1]) {
    printf("swapping freqs...\n");
    std::swap(freqArray[0], freqArray[1]);
    std::swap(freqArray[2], freqArray[3]);
    std::swap(freqArray[4], freqArray[6]);
    std::swap(freqArray[5], freqArray[7]);
    std::swap(freqOffset[0], freqOffset[1]);
  }

  //check what flags where provided and set the freqs aproprietly:
  //adjust the range for first frequency:
  if (freqOffset[0] == 0) {
    fitFreq = 3100.0;
    fitFreqMin = (long double) (fitFreq - 100.0);
    fitFreqMax = (long double) (fitFreq + 100.0);
  } else if (freqOffset[0] == 1) {
    fitFreq = freqArray[0];
    fitFreqMin = (long double) (fitFreq - freqArray[2]);
    fitFreqMax = (long double) (fitFreq + freqArray[2]);
  } else if (freqOffset[0] == 2) {
    fitFreqMin = (long double) (freqArray[4]);
    fitFreqMax = (long double) (freqArray[5]);
    fitFreq = (fitFreqMax - fitFreqMin) / 2.0;
  } else {
    dmFreq = -1.0;
    fitFreq = -1.0;
    fitFreqMin = fitFreqMax = dmFreqMin = dmFreqMax = -1.0;
  }
  //adjust the range for second frequency:
  if (freqOffset[1] == 0) {
    dmFreq = 685.0;
    dmFreqMin = (long double) (dmFreq - 100.0);
    dmFreqMax = (long double) (dmFreq + 100.0);
  } else if (freqOffset[1] == 1) {
    dmFreq = freqArray[1];
    dmFreqMin = (long double) (dmFreq - freqArray[3]);
    dmFreqMax = (long double) (dmFreq + freqArray[3]);
  } else if (freqOffset[1] == 2) {
    dmFreqMin = (long double) (freqArray[6]);
    dmFreqMax = (long double) (freqArray[7]);
    dmFreq = (dmFreqMax - dmFreqMin) / 2.0;
  } else {
    dmFreq = -1.0;
    fitFreq = -1.0;
    fitFreqMin = fitFreqMax = dmFreqMin = dmFreqMax = -1.0;
  }
  nf = 0;
  //identify observations at given freqs:
  if (freqOffset[0] != 3 && freqOffset[1] != 3) {
    for (i = 0; i < psr[0].nobs; ++i) {
      if (!(psr[0].obsn[i].freq > fitFreqMin && psr[0].obsn[i].freq < fitFreqMax) && psr[0].obsn[i].deleted == 0) {
	if (psr[0].obsn[i].freq > dmFreqMin && psr[0].obsn[i].freq < dmFreqMax) {
	  psr[0].obsn[i].deleted = 1;
	  dmObs[dmCount] = i;
	  (dmCount)++;
	} else {
	  psr[0].obsn[i].deleted = 1;
	  impObs[impCount] = i;
	  (impCount)++;
	}
      } else if (psr[0].obsn[i].deleted == 0) {
	fitObs[fitCount] = i;
	(fitCount)++;
	//that's needed for the jumps:
	//it goes through all observations and turns on jumps only at
	// used frequencies
	// TODO:
	// it works only for the FJUMP -f, doesn't work with JUMP in the par file
	if (strlen(psr[0].fjumpID)>0) { // this for the FJUMP case
	  for (k = 0; k < psr[0].obsn[i].nFlags; ++k) {
	    found = 0;
	    if (strcmp(psr[0].obsn[i].flagID[k], psr[0].fjumpID) == 0) {
	      for (l = 0; l<nf; l++) {
		strncpy(_tmp, "", sizeof (_tmp));
		strcat(_tmp, psr[0].fjumpID);
		strcat(_tmp, " ");
		strcat(_tmp, psr[0].obsn[i].flagVal[k]);
		if (strcmp(_tmp, psr[0].jumpStr[l]) == 0) {
		  psr[0].obsn[i].jump = l;
		  found = 1;
		  break;
		}
	      }
	      if (found == 0) {
		sprintf(psr[0].jumpStr[nf], "%s %s", psr[0].fjumpID, psr[0].obsn[i].flagVal[k]);
		psr[0].obsn[i].jump = nf;
		(nf)++;
	      }
	    }
	  }
	} else // this if for the case were we have JUMP flag, not FJUMP -f
	{
	  // TODO:
	  // implement it :) base it on readParFile()
	  // psr->nJumps should be set already by readParfile, need to set nf:
	  // and preProcess around 669 (CHECK JUMPS)
	  nf = psr[0].nJumps+1;
	  for (k=0; k < psr[0].obsn[i].nFlags; ++k)
	  {
	    strncpy(_tmp, "", sizeof (_tmp));
	    strcat(_tmp, psr[0].obsn[i].flagID[k]);
	    strcat(_tmp, " ");
	    strcat(_tmp, psr[0].obsn[i].flagVal[k]);
	    for (l = 0; l < nf; l++)
	    {
	      if (strcmp(_tmp, psr[0].jumpStr[l]) == 0) 
	      {
		psr[0].obsn[i].jump = l;
		break;
	      }
	    }
	  }

	}
      }
    }
  }
  else {
    char flag1[100], filtS1[MAX_STRLEN], *filt1;
    char flag2[100], filtS2[MAX_STRLEN], *filt2;
    strcpy(filtS1, freq1f);
    strcpy(filtS2, freq2f);
    filt1 = strtok(filtS1, " ");
    if (filt1[0] == '-')
      strcpy(flag1, filt1);
    else {
      printf("ERROR: wrong flag provided for -freq1f\n");
      exit(-1);
    }
    filt1 = strtok(NULL, " ");

    filt2 = strtok(filtS2, " ");
    if (filt2[0] == '-')
      strcpy(flag2, filt2);
    else {
      printf("ERROR: wrong flag provided for -freq2f\n");
      exit(-1);
    }
    filt2 = strtok(NULL, " ");
    for (i = 0; i < psr[0].nobs; ++i) {
      foundFlag = 0;
      if (psr[0].obsn[i].deleted == 0) {
	for (j = 0; j < psr[0].obsn[i].nFlags; ++j) {
	  if (strcmp(psr[0].obsn[i].flagID[j], flag1) == 0 && strcmp(psr[0].obsn[i].flagVal[j], filt1) == 0) {
	    fitObs[fitCount] = i;
	    (fitCount)++;
	    //that's needed for the jumps:
	    //it goes through all observations and turns on jumps only at
	    // used frequencies
	    found = 0;
	    for (k = 0; k < psr[0].obsn[i].nFlags; ++k) {
	      if (strcmp(psr[0].obsn[i].flagID[k], psr[0].fjumpID) == 0) {
		for (l = 0; l<nf; ++l) {
		  strncpy(_tmp, "", sizeof (_tmp));
		  strcat(_tmp, psr[0].fjumpID);
		  strcat(_tmp, " ");
		  strcat(_tmp, psr[0].obsn[i].flagVal[k]);
		  if (strcmp(_tmp, psr[0].jumpStr[l]) == 0) {
		    found = 1;
		    psr[0].obsn[i].jump = l;
		    break;
		  }
		}
		if (found == 0) {
		  sprintf(psr[0].jumpStr[nf], "%s %s", psr[0].fjumpID, psr[0].obsn[i].flagVal[k]);
		  psr[0].fitJump[nf] = 1;
		  psr[0].obsn[i].jump = nf;
		  if (1==1 || debugFlag >= 1) {
		    printf("found jump %d: %s in obs %d\n", nf, psr[0].jumpStr[nf], i);
		    printf("data file: %s flag(%d):%s\n", psr[0].obsn[i].fname, k, psr[0].obsn[i].flagVal[k]);
		  }
		  (nf)++;
		}
	      }
	    }
	    foundFlag = 1;
	  } else if (strcmp(psr[0].obsn[i].flagID[j], flag2) == 0 && strcmp(psr[0].obsn[i].flagVal[j], filt2) == 0) {
	    psr[0].obsn[i].deleted = 1;
	    dmObs[dmCount] = i;
	    (dmCount)++;
	    foundFlag = 1;
	  }
	}
      }
      if (foundFlag == 0) {
	psr[0].obsn[i].deleted = 1;
	impObs[impCount] = i;
	(impCount)++;
      }
    }
    if (debugFlag >= 1) {
      printf("flag1=%s flag2=%s filt1=%s filt2=%s\n", flag1, flag2, filt1, filt2);
      printf("fitCount=%d dmCount=%d impCount=%d nf=%d\n", fitCount, dmCount, impCount, nf);
    }
  }
  if (fitCount == 0 || dmCount == 0) {
    fprintf(stderr,"\nERROR!\nFound %d points at first specified frequency and %d at the second\n ", fitCount, dmCount);
    fprintf(stderr,"More than zero points required at both frequencies!!\n");
    fprintf(stderr,"Check how you defined the frequencies to be used\n");
    fprintf(stderr,"(see flags freq1, freq1r, freq1f, freq2, freq2r and freq2f)\n");
    fprintf(stderr,"calcDMe will exit now\n");
    exit(-1);
  }
  psr[0].nJumps = (nf) - 1;
} //handleFreqPoints()

//interpolation (spline) and smoothing (Hann):
//this function interpolates calculated deltaDMs using constrained spline, smooths it and resamples at 1/day frequency

void interpolateSplineSmooth() {
  //array needed by TKcmonot
  double yd[MAX_OBSN][4];
  //arrays to cast input onto doubles
  double _ddm[MAX_OBSN], _ddmMJD[MAX_OBSN];
  //output arrays with interpolation result
  double interpX[MAX_OBSN], interpY[MAX_OBSN];
  //auxilary 'i' and count of interpolated points
  int i, nInterp;
  //smoothed output and count of smoothed points
  double smoothX[MAX_OBSN], smoothY[MAX_OBSN];
  int nSmooth;

  //cast input to doubles
  for (i = 0; i < ddmCount; ++i) {
    _ddm[i] = (double) ddm[i];
    _ddmMJD[i] = (double) ddmMJD[i];
  }

  //interpolation (from TKspecturm.h)
  TKcmonot(ddmCount, _ddmMJD, _ddm, yd);
  nInterp = 0;
  do {
    interpX[nInterp] = _ddmMJD[0] + nInterp;
    nInterp++;
  } while (interpX[nInterp - 1] < _ddmMJD[ddmCount - 1]);
  nInterp--;
  TKspline_interpolate(ddmCount, _ddmMJD, _ddm, yd, interpX, interpY, nInterp);
  //smooth the output
  if (smoothWidth < 1) {
    smoothWidth = 1;
    printf("smoothWidth can't be smaller than 1 => resetting to one..\n");
  }
  TKhann(interpX, interpY, nInterp, smoothX, smoothY, &nSmooth, smoothWidth);

  for (i = 0; i < nInterp; ++i) {
    outX[i] = interpX[i];
    outY[i] = interpY[i];
    if (i < nSmooth) {
      outX[i + nInterp] = smoothX[i];
      outY[i + nInterp] = smoothY[i];
    }
  }
  outInterpCount = nInterp;
  outSmoothCount = nSmooth;

} //interpolateSplineSmooth

//this function follows Daniel Yardley's weighted interpolation/smoothing scheme
void interpolateWeightedSmooth() {
  //int separation = 1; //number of days between post-interpolated data points (I guess we don't really need that - in our case it will be 1 always?)
  //number of post-interpolation residuals:
  int nInterp = (int) ceill(floorl(ddmMJD[ddmCount-1])-ceill(ddmMJD[0]));
  int currentDay = -1;
  double sum;
  double weight;
  //output arrays with interpolation result
  double interpX[MAX_OBSN], interpY[MAX_OBSN], interpYErr[MAX_OBSN];
  //smoothed output and count of smoothed points
  double smoothX[MAX_OBSN], smoothY[MAX_OBSN], smoothYErr[MAX_OBSN];
  outSmoothCount = outInterpCount = nInterp;
  for ( int ii = 0; ii<nInterp ; ++ii) {
    sum=0;
    currentDay =  (int) (ddmMJD[0] + ii);
    interpX[ii] = (double) currentDay;
    interpYErr[ii] = interpY[ii] = 0;

    for (int kk = 0; kk < ddmCount; ++kk) {
      weight = exp(-pow(ddmMJD[kk] - currentDay,2) / 2.0 / pow(smoothWidth,2)) / pow((double)ddmErr[kk],2); 
      interpY[ii] += weight * ((double)ddm[kk]);
      interpYErr[ii] += pow(weight*((double)ddmErr[kk]),2);
      sum += weight;
    }
    interpY[ii] /= sum;
    interpYErr[ii] = sqrt(interpYErr[ii] / pow(sum,2));
  }
  //actually store it in proper variables:
  for (int ii = 0; ii < nInterp; ++ii) {
    outX[ii+nInterp] = outX[ii] = interpX[ii];
    outY[ii+nInterp] = outY[ii] = interpY[ii];
  }
  /*  free(&nInterp);
  free(interpX);
  free(interpY);
  free(interpYErr);*/
}//interpolateWeightedSmooth()

//This function writes output
//TODO: I think there's sth wrong with -mean -outDM or sth

void output(char *outFileName, int ascii, double dm0, int header, int outDM, double *outX, double *outY, int outInterpCount, int outSmoothCount, int mean, int meanMJD, double *meanMJDval, double *meanVal, int splineOut, int rawOut, long double *ddmMJD, long double *ddm, long double *ddmErr, int ddmCount) {
  int i;
  double smoothX[MAX_OBSN], smoothY[MAX_OBSN];
  FILE *outFile;
  //auxilary char array
  char *tmp;

  for (i = 0; i < outSmoothCount; ++i) {
    smoothX[i] = outX[outInterpCount + i];
    smoothY[i] = outY[outInterpCount + i];
  }
  // removing means, when necessary
  // output ddm or dm?
  if (meanMJD == 1) {
    *meanMJDval = findMean(smoothX, outSmoothCount);
  } else
    *meanMJDval = 0;
  if (mean == 1) {
    if (outDM == 1) {
      for (i = 0; i < outSmoothCount; ++i) {
	smoothY[i] += dm0;
      }
      *meanVal = findMean(smoothY, outSmoothCount);
    } else
      *meanVal = findMean(smoothY, outSmoothCount);
  } else
    *meanVal = 0;

  //proper output
  outFile = fopen(outFileName, "w");
  if (outFile != NULL) {
    if (debugFlag==1) printf("opened output file\n");
    if (ascii == 0) {
      if (header == 1) {
	fwrite(&dm0, sizeof (double), 1, outFile);
	fwrite(&outSmoothCount, sizeof (int), 1, outFile);
      }
      for (i = 0; i < outSmoothCount; ++i) {
	smoothX[i] -= *meanMJDval;
	fwrite(&smoothX[i], sizeof (double), 1, outFile);
	smoothY[i] -= *meanVal + (double) (outDM - mean * outDM) * dm0;
	fwrite(&smoothY[i], sizeof (double), 1, outFile);
      }
    } else {
      if (header == 1) fprintf(outFile, "#%lf %d %lf\n", dm0, outSmoothCount, findMean(smoothY, outSmoothCount));
      for (i = 0; i < outSmoothCount; ++i) {
	fprintf(outFile, "%lf %lf\n", smoothX[i] - *meanMJDval, smoothY[i] + (double) (outDM - mean * outDM) * dm0 - *meanVal);
      }
    }
    fclose(outFile);
  } else {
    printf("=========WARNING=========\nFAILED TO OPEN THE OUTPUT FILE\n");
    printf("attempting to print the dm0,numberOfOutPoints followed by mjd,smoothed_dm variations into file /tmp/smoothedTMP::\n\n");
    printf("to extract this data do:\n");
    printf("awk '{if (NR==1) print $2; else print $2,$3}' /tmp/smoothedTMP\n");
    outFile = fopen("/tmp/smoothedTMP", "w");
    if (outFile != NULL) {
      fprintf(outFile, "%lf %d\n", dm0, outSmoothCount);
      for (i = 0; i < outSmoothCount; ++i) {
	fprintf(outFile, "SMOOTH: %lf %lf\n", smoothX[i], smoothY[i]);
      }
    }
  }
  //output data before smoothing?
  if (splineOut == 1) {
    tmp = (char*) malloc(100 * sizeof(char));
    tmp = strcpy(tmp, outFileName);
    tmp = strcat(tmp, ".int");

    outFile = fopen(tmp, "w");
    free(tmp);

    if (header == 1) fprintf(outFile, "#%lf %d %lf\n", dm0, outInterpCount, findMean(outY, outInterpCount));
    for (i = 0; i < outInterpCount; ++i) {
      fprintf(outFile, "%lf %lf\n", outX[i] - (double) (meanMJD * (*meanMJDval)), outY[i] + (double) ((outDM - mean * outDM) * dm0) - *meanVal);
    }

    fclose(outFile);
  }
  //output the observed data?
  if (rawOut == 1) {
    tmp = (char*) malloc(100 * sizeof(char));
    tmp = strcpy(tmp, outFileName);
    tmp = strcat(tmp, ".raw");

    outFile = fopen(tmp, "w");
    free(tmp);

    if (header == 1) fprintf(outFile, "#%lf %d\n", dm0, ddmCount);
    for (i = 0; i < ddmCount; ++i) {
      fprintf(outFile, "%Lf %Lf %Lf\n", ddmMJD[i] - (long double) (meanMJD * (*meanMJDval)), ddm[i] + (long double) ((outDM - mean * outDM) * dm0) - *meanVal, ddmErr[i]);
    }

    fclose(outFile);
  }
} //output

//locate the first observation from the fitObs - that's the beginning of our first point

void findFirst(pulsar *psr) {
  int l, i;
  l = fitObs[0];
  binStart = psr[0].obsn[l].sat;
  for (i = 1; i < fitCount; ++i) {
    l = fitObs[i];
    if (psr[0].obsn[l].sat < binStart) {
      binStart = psr[0].obsn[l].sat;
    }
  }
} //findFirst()

//This function searches for data points in a given bin, coming from observations at freq1 and freq2

void get_binObs(pulsar *psr) {
  int i, l;
  bin_fitCount_inc = 0;
  bin_dmCount_inc = 0;
  //auxilary variable for checking whether a point was used before
  int wasUsed;
  //effective bin size (it may be bigger than binSizeDays if sessionSeparation > 0)
  long double effectiveBinEnd;

  //set the bin size:
  if (sessionSeparation < 0) // fix the < and or > (next else if)
  {
    effectiveBinEnd = binStart + binSizeDays;
    //find the session during which the bin would normally end 
    //and adjust the end of bin to include the whole session
    for (i = 0; i < nSessions; ++i) {
      if ( effectiveBinEnd > start_sessions[i] && effectiveBinEnd < finish_sessions[i]) 
      {
	effectiveBinEnd = finish_sessions[i];
	break;
      }
    }
  } else if (sessionSeparation > 0) 
  {
    binStart = start_sessions[lastUsedSession+1];
    effectiveBinEnd = finish_sessions[lastUsedSession+1];
    lastUsedSession++;
  } else 
  {
    effectiveBinEnd = binStart + binSizeDays;
  }

  //first the fitObs:
  for (i = 0; i<fitCount; ++i) {
    l = fitObs[i];
    psr[0].obsn[l].deleted = 1;
    if (psr[0].obsn[l].sat >= binStart && psr[0].obsn[l].sat < effectiveBinEnd) {
      //check if it wasn't used before
      wasUsed = 0;
      //if wasn't used, use it in the fit:
      if (wasUsed == 1) {
	psr[0].obsn[l].deleted = 1;
      } else {
	psr[0].obsn[l].deleted = 0;
	binObs[bin_fitCount][0] = l;
	(bin_fitCount)++;
	(bin_fitCount_inc)++;
      }
    }
  }
  //and now the dmObs:
  for (i = 0; i<dmCount; ++i) {
    l = dmObs[i];
    psr[0].obsn[l].deleted = 1;
    if (psr[0].obsn[l].sat >= binStart && psr[0].obsn[l].sat < effectiveBinEnd) {
      //check if it wasn't used before
      wasUsed = 0;
      //if wasn't used, use it in the fit:
      if (wasUsed == 1) {
	psr[0].obsn[l].deleted = -1;
      } else {
	psr[0].obsn[l].deleted = 0;
	binObs[bin_dmCount][1] = l;
	(bin_dmCount)++;
	(bin_dmCount_inc)++;
      }
    }
  }


} //get_binObs()

// this function turns off fitting of everything
// (including jumps and dm derivatives)
// except for dm and f0

void setFitParams(pulsar *psr) {
  int i, k;
  //turn all off:
  for (i = 0; i < MAX_PARAMS; ++i) {
    for (k = 0; k < psr[0].param[i].aSize; ++k) {
      psr[0].param[i].fitFlag[k] = 0;
    }
  }
  //reenable dm and f
  psr[0].param[param_dm].fitFlag[0] = 1;
  if (f0fit == 1)
    psr[0].param[param_f].fitFlag[0] = 1;
  //turn off jumps
  for (i = 0; i <= psr[0].nJumps; ++i) { // for some reason there has to be <= not < (the jumps start at 1st not 0th entry)
    psr[0].fitJump[i] = 0;
  }
  //make sure no dm derivatives are used:
  for (i = 1; i < psr[0].param[param_dm].aSize; ++i) {
    psr[0].param[param_dm].paramSet[i] = 0;
  }
} //setFitParams();

// this function marks all observations as deleted:

void setAllDeleted(pulsar *psr) {
  int i;
  for (i = 0; i < psr[0].nobs; i++) {
    psr[0].obsn[i].deleted = 1;
  }
} //setAllDeleted()

// this functions sets the dm and f0 values to the ones obtained from the first fit

void resetDMandF0(pulsar *psr) {
  psr[0].param[param_dm].val[0] = dm0;
  psr[0].param[param_dm].err[0] = dm0_err;
  psr[0].param[param_f].val[0] = f0_0;
  psr[0].param[param_f].err[0] = f0_0_err;
} //resetDMandF0()

/* This function calls all of the fitting routines.             */
/* The function is looped twice, the first time for the pre-fit */

/* residuals and the second time for the post-fit residuals     */
void callFit(pulsar *psr, int npsr) {
  int iteration, i, p, it, k;
  double globalParameter = 0.0;
  static int call=0;
  char str[100];

  for (it = 0; it < psr[0].nits; it++) {
    if (it > 0) /* Copy post-fit values to pre-fit values */ {
      for (i = 0; i < MAX_PARAMS; i++) {
	for (p = 0; p < npsr; p++) {

	  for (k = 0; k < psr[p].param[i].aSize; k++) {
	    psr[p].param[i].prefit[k] = psr[p].param[i].val[k];
	    psr[p].param[i].prefitErr[k] = psr[p].param[i].err[k];
	  }
	}
      }
    }
    for (iteration = 0; iteration < 2; iteration++) {
      formBatsAll(psr, npsr);
      /* Form residuals */
      formResiduals(psr, npsr, iteration);
      /* Do the fitting */
      if (iteration == 0) {
	if (strcmp(dcmFile, "NULL") == 0)
	  doFit(psr, npsr, 0);
	else
	  doFitDCM(psr, dcmFile, npsr, 0);
      } else if (debugFlag >= 1 || 1==1) 
	{
	  sprintf(str,"fit%d.par",call);
	  textOutput(psr, npsr, globalParameter, 0, 0, 1, str);
	}
      /*      if (iteration==1)
	{
	  for (i=0;i<psr[0].nobs;i++)
	    {
	      if (psr[0].obsn[i].deleted==0)
		printf("Out: %g\n",(double)psr[0].obsn[i].sat);
	    }
	    }*/
    }
  }
  sprintf(str,"test%d",call);
  writeTim(str,psr,"tempo2");
  /*  exit(1);*/
  if (psr[0].nFit == 0) {
    printf("No timing residuals to plot.  Please check your TOA file and filter commands\n");
    exit(1);
  }
  call++;
} //callFit()

//This function calculates the mean of data

double findMean(double *x, int count) {
  int i = 0;
  double mean = 0;

  for (i = 0; i < count; ++i)
    mean += x[i];

  return mean / (double) count;
} // findMean()

//this function generates plot of raw, interpolated and smoothed results
//displays it either on screen (using /XS device) or saves it to a postscript file

void display(char *gr, int publish, double *xx, double *yy, long double *ddmMJD, long double *ddm, long double *ddmErr, int outInterpCount, int outSmoothCount, int ddmCount, char *xlab, char *ylab, char *title, double meanMJDval, double meanVal) {
  int i;
  float fontSize = 1.0;
  int fontType = 1;
  int lineWidth = 1;

  float x[2 * MAX_OBSN], y[2 * MAX_OBSN];
  float yerr1[MAX_OBSN], yerr2[MAX_OBSN];

  float minx, maxx, miny, maxy, plotx1, plotx2, ploty1, ploty2;

  cpgopen(gr);
  cpgscf(fontType);
  cpgslw(lineWidth);

  cpgask(0);
  cpgsch(fontSize);
  if (strcmp(gr, "/XS") == 0) {
    cpgeras();
  }

  for (i = 0; i < outInterpCount; ++i) {
    x[i] = (double) (xx[i]) - meanMJDval;
    y[i] = (double) (yy[i]) - meanVal;
  }
  minx = x[0];
  maxx = x[outInterpCount - 1];

  plotx1 = minx - fabs(maxx - minx)*0.1;
  plotx2 = maxx + fabs(maxx - minx)*0.1;

  miny = y[0];
  maxy = y[0];
  for (i = 1; i < outInterpCount; ++i) {
    if (y[i] < miny)
      miny = y[i];
    if (y[i] > maxy)
      maxy = y[i];
  }
  ploty1 = miny - fabs(maxy - miny)*0.1;
  ploty2 = maxy + fabs(maxy - miny)*0.1;

  cpgenv(plotx1, plotx2, ploty1, ploty2, 0, -1);
  cpgbox("BCNST1", 0.0, 0, "BCNST1", 0.0, 0);

  cpglab(xlab, ylab, title);
  cpgsch(0.6);
  if (strcmp(gr, "/XS") == 0) {
    cpgmtxt("B", 6.5, 0.9, 0.5, "calcDMe, Copyright (C) 2009 Stefan Oslowski");
  }
  cpgsch(fontSize);

  //plot interpolated data
  cpgsci(2);
  cpgsls(2);
  cpgline(outInterpCount, x, y);

  //plot smoothed data
  cpgsci(1);
  cpgsls(1);
  for (i = 0; i < outSmoothCount; ++i) {
    x[i] = (double) (xx[i + outInterpCount]) - meanMJDval;
    y[i] = (double) (yy[i + outInterpCount]) - meanVal;
  }
  cpgline(outSmoothCount, x, y);

  //plot actual data points with errors
  cpgsci(3);
  for (i = 0; i < ddmCount; ++i) {
    x[i] = (double) (ddmMJD[i]) - meanMJDval;
    y[i] = (double) (ddm[i] - meanVal);
    yerr1[i] = (double) (ddm[i] + ddmErr[i] - meanVal);
    yerr2[i] = (double) (ddm[i] - ddmErr[i] - meanVal);
  }
  cpgpt(ddmCount, x, y, 0);
  cpgerry(ddmCount, x, yerr1, yerr2, 1);

  cpgend();
} // display()

void findSessions(pulsar *psr)
{
  double obsTimes[fitCount + dmCount];
  double gapSize[fitCount + dmCount];
  int pointsInSession[fitCount + dmCount];
  int i,j;

  //setting the obsTimes from fitObs
  for (i=0;i<fitCount;++i)
  {
    obsTimes[i] = (double)psr[0].obsn[fitObs[i]].sat;
    pointsInSession[i] = 0;
  }
  //setting the obsTimes from dmObs
  for (i=0;i<dmCount;++i) 
  {
    obsTimes[i+fitCount] = (double)psr[0].obsn[dmObs[i]].sat;
    pointsInSession[i+fitCount] = 0;
  }
  // Sort the observing times
  TKsort_d(obsTimes,fitCount+dmCount);

  // Find gap sizes
  for (i = 0; i < dmCount + fitCount - 1; ++i)
  {
    gapSize[i] = obsTimes[i+1]-obsTimes[i];
  }
  // Now determine the sessions

  start_sessions[0]=obsTimes[0]-0.5;
  for (i=0;i<dmCount + fitCount - 2;i++)
  {
      pointsInSession[nSessions]++;
      if (gapSize[i] > sessionSeparation)
      {
	//	pointsInSession[nSessions]--;
	finish_sessions[nSessions] = obsTimes[i]+0.5;
	nSessions++;
	start_sessions[nSessions] = obsTimes[i+1]-0.5;
      }
  }
  finish_sessions[nSessions++] = obsTimes[dmCount+fitCount-1]+0.5;
}//findSessions

//END
char * plugVersionCheck = TEMPO2_h_VER;
