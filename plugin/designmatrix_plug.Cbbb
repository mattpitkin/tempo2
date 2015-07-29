/* designmatrix_plug.C -- This plugin outputs the design matrix for the timing
 * model parameters, given that the parameters are marked with '1' in the .par
 * file.
 *
 * Rutger van Haasteren August 10th 2011 vhaasteren@gmail.com
 *
 * Copyright (C) 2011 Rutger van Haasteren.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"

using namespace std;


/* The number of fitting parameters for a pulsar */
int tempo2_GetNumberOfParameters(pulsar *pPsr);

/* Construct the linearised timing model, and form timing residuals */
void ProcessTempo2Objects(int argc, char *argv[], pulsar *pPsr, int *pnPsr, char parFile[][MAX_FILELEN], char timFile[][MAX_FILELEN]);

/* Write the design matrix to an ascii file */
void WriteDesignMatrix(const char *strFileName, pulsar *pPsr, int nTargetPulsar);



/* Description
 * This function returns the number of parameters of a single pulsar
 *
 * Input
 * pPsr:	Pointer to the pulsar structure
 *
 * Output
 *
 * Returns an integer with the number of parameters
 * */
int tempo2_GetNumberOfParameters(pulsar *pPsr) {
  int nPol;
  int i, k;

  nPol = 1;
  for(i=0; i<MAX_PARAMS; i++) {
    for(k=0; k<pPsr->param[i].aSize; k++) {
      if(pPsr->param[i].fitFlag[k]==1) {
	if(i!=param_start && i!=param_finish) {
	  nPol++;
	}
      } // if fitFlag
    } // for k
  } // for i
  
  /* Add extra parameters for jumps */
  for(i=1; i<=pPsr->nJumps; i++) {
    if(pPsr->fitJump[i]==1)
      nPol++;
  }

  /* Add extra parameters for sinusoidal whitening */
  if(pPsr->param[param_wave_om].fitFlag[0]==1)
    nPol += pPsr->nWhite*2-1;
  
  return nPol;
} /* tempo2_GetNumberOfParameters */



/* Description
 * This function constructs the linearised timing model so that we can extract
 * the design matrix, and the timing residuals are formed.
 *
 * Input
 * argc:	Number of arguments on command line
 * argv:	The commend line arguments
 * pPsr:	Structures corresponding to the par-files
 * pnPsr:	Number of pulsars in pPsr
 * parFile:	Filenames of parfiles
 * timFile:	Filenames of timfiles
 *
 * Output
 *
 * Return
 * */
void ProcessTempo2Objects(int argc, char *argv[], pulsar *pPsr, int *pnPsr, char parFile[][MAX_FILELEN], char timFile[][MAX_FILELEN]) {
  if(*pnPsr > 0) {
    readParfile(pPsr, parFile, timFile, *pnPsr); /* Load the parameters       */
    readTimfile(pPsr, timFile, *pnPsr);         /* Load the arrival times    */
    preProcess(pPsr, *pnPsr, argc, argv);
    formBatsAll(pPsr, *pnPsr);                 /* Form the barycentric arrival times */
    formResiduals(pPsr, *pnPsr, 0.0);           /* Form the residuals                 */
  } /* if pnPsr */
} /* ProcessTempo2Objects */



/* Description
 * This function writes the design matrix of pulsar pPsr[nTargetPulsar] to the
 * ascii file strFileName. The file strFileName is truncated without warning if
 * it already exists.
 *
 * Makeup of the file strFileName:
 * 
 * nObservations  nParameters
 * - empty line -
 * TOA[i]  residual[i]  errorbar[i]  (nObservations lines)
 * - empty line -
 * PrefitParameter[p]  # Parameter description  (nParameters lines)
 * - empty line -
 * Designmatrix  (nObservations lines, nParameters columns)
 *
 * Input
 * strFileName:		String to file name of designmatrix ascii file
 * pPsr:		Array of pulsars as given by tempo2
 * nTargetPulsar:	Index of pulsar we want to know design matrix of
 *
 * Output
 *
 * Return
 *
 * */
void WriteDesignMatrix(const char *strFileName, pulsar *pPsr, int nTargetPulsar) {
  int i, j, k;
  int nObservations, nParameters, nTempo2ParameterIndex;
  double pdParamDeriv[MAX_PARAMS], dMultiplication;

  /* Now we need to set the matrix with all the marginalization and fitting
   * parameter derivatives. We get them from the tempo2 routines. */
  FILE *pFile;

  nObservations = pPsr[nTargetPulsar].nobs;
  nParameters = tempo2_GetNumberOfParameters(&pPsr[nTargetPulsar]);

  if(! (pFile = fopen(strFileName, "w+")) ) {
    /* unable to open file */
    printf("Unable to open file: %s\n", strFileName);
    return;
  } /* if pFile */

  /* First write the number of observations, and the number of timing model
   * parameters */
  fprintf(pFile, "%i \t %i\n", nObservations, nParameters);
  fprintf(pFile, "\n");

  /* Write the TOAs, timing residuals, and the error bars */
  for(i=0; i<nObservations; i++) {
    fprintf(pFile, "%22.20e \t %22.20e \t %22.20e\n", ((double)(pPsr[nTargetPulsar].obsn[i].bat)-(double)(pPsr[nTargetPulsar].param[param_pepoch].val[0])) * 86400.0, (double)(pPsr[nTargetPulsar].obsn[i].residual), (double)(1.0E-6*pPsr[nTargetPulsar].obsn[i].toaErr));
  } /* for i */
  fprintf(pFile, "\n");

  /* Write the best-guess parameters, around which we linearise the timing
   * model. At each line we also print the parameter description */
  fprintf(pFile, "%22.20e \t # %s\n", 0.0, "Offset [sec]");
  nTempo2ParameterIndex = 1;

  /* Loop over all Tempo2 parameters to find the right ones */
  for(j=0; j<MAX_PARAMS; j++) {
    for(k=0; k<pPsr[nTargetPulsar].param[j].aSize; k++) {
      if(pPsr[nTargetPulsar].param[j].fitFlag[k]==1) {
	/* This parameter is fitted for */
	if(j!=param_start && j!=param_finish) {
	  /* Add this parameter to the descriptions */

	  /* Units are as in textOutput.C of tempo2: correct it */
	  if(j == param_pmra) {
	    dMultiplication = 1.0;
	  } else if(j == param_pmdec) {
	    dMultiplication = 1.0;
	  } else if((j == param_raj || j == param_decj) &&
	             pPsr[nTargetPulsar].eclCoord == 1) {
	    dMultiplication = 180.0 / M_PI;
	  } else if(j == param_sini &&
	            pPsr[nTargetPulsar].param[param_sini].nLinkTo != 0) {
	    dMultiplication = M_PI/180.0;
	  } else {
	    dMultiplication = 1.0;
	  } // if j

	  fprintf(pFile, "%22.20e \t # %s\n",
	      (double)pPsr[nTargetPulsar].param[j].prefit[k],
	      pPsr[nTargetPulsar].param[j].label[k]);
	    nTempo2ParameterIndex++;
	} // if this is not the start/finish parameter
      } // if fitFlag
    } // for k
  } // for j
  fprintf(pFile, "\n");


  /* Write the parameter derivatives (design matrix row) for each observation */
  for(i=0; i<nObservations; i++) {
    FITfuncs(pPsr[nTargetPulsar].obsn[i].bat - pPsr[nTargetPulsar].param[param_pepoch].val[0], pdParamDeriv, nParameters, pPsr+nTargetPulsar, i,0);
    for(j=0; j<nParameters; j++) {
      fprintf(pFile, "%22.20e", pdParamDeriv[j]);
      if(j!=nParameters-1) {
	fprintf(pFile, " \t ");
      } /* if j */
    } /* for j */
    fprintf(pFile, "\n");
  } /* for i */

  fclose(pFile);
} /* WriteDesignMatrix */




/* This function should contain usage information about the plugin which
 * should (in general) be accessed by the user pressing 'h'
 * */
void help() {
  printf("--------------------------------------\n");
  printf("Command line inputs:\n\n");
  printf("-f     parfile.par timfile.tim:  input par and tim files\n");
  printf("-h     this help file\n");
  printf("-p     Select the pulsar we need to write the design matrix of\n");
  printf("\n\n");
  printf("Example usage: tempo2 -gr designmatrix -f tt.par 0437.2048.tim\n\n");
  printf("               tempo2 -gr designmatrix -f tt.par 0437.2048.tim -p 0\n\n");
  printf("--------------------------------------\n");
} /* help */


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc, char *argv[], pulsar *pPsr, int *pnPsr) {
  int i;
  int nTargetPulsar=0;
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  const char strDesignMatrixFile[80] = "designmatrix.txt";

  printf("Graphical Interface: designmatrix\n");
  printf("Author:              R. van Haasteren\n");
  printf("Version:             1.0\n");
  printf(" --- type 'h' for help information\n");

  /* Obtain the number of pulsars ourselves */
  *pnPsr = 0;

  for(i=2; i<argc; i++) {
    if(strcmp(argv[i], "-f")==0 && i+2 < argc) {
      strcpy(parFile[*pnPsr], argv[i+1]); 
      strcpy(timFile[*pnPsr], argv[i+2]);	  
      (*pnPsr)++;
    } /* if strcmp */
  } /* for i */

  /* Obtain all parameters from the command line */
  for(i=2; i<argc; i++) {
    if (strcmp(argv[i],"-f") == 0 && i+2 < argc) {
      strcpy(parFile[*pnPsr], argv[i+1]); 
      strcpy(timFile[*pnPsr], argv[i+2]);	  
      (*pnPsr)++;
    } else if(strcmp(argv[i],"-h") == 0) {
      help();
      return 0;
    } else if(strcmp(argv[i],"-p") == 0 && i+1 < argc) { 
      /* Set the target pulsar */
      nTargetPulsar = atoi(argv[i+1]);
    } /* if strcmp */
  } /* for i */

  /* Obtain the number of pulsars, and process the par/tim files */
  ProcessTempo2Objects(argc, argv, pPsr, pnPsr, parFile, timFile);
  if(*pnPsr < 1) {
    printf("No pulsars to process\n");
    return 0;
  } /* if pnPsr */

  if(nTargetPulsar >= 0 && nTargetPulsar < *pnPsr) {
    /* Work for pulsar nTargetPulsar */
    WriteDesignMatrix(strDesignMatrixFile, pPsr, nTargetPulsar);
  } else {
    printf("No valid pulsar selected\n");
  }/* if nTargetPulsar */

  return 0;
} /* graphicalInterface */
