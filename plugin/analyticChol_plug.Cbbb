//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

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
#include "tempo2.h"
#include "TKspectrum.h"
#include "TKfit.h"
#include "T2toolkit.h"
#include "choleskyRoutines.h"

using namespace std;


void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i;
  double globalParameter;
  const char *CVS_verNum = "$Revision: 1.9 $";
  double modelAlpha,modelFc,modelA;
  int np,ndays;
  double resx[MAX_OBSN],resy[MAX_OBSN],rese[MAX_OBSN];
  double covFunc[MAX_OBSN];
  double whiteNoise;
  char fname[128];
  double rms;
  FILE *fout;
  double gwamp=1.0e-15;
  modelAlpha = 13.0/3.0;

  modelFc = -1;
  modelA = -1; // 1.819e-26;
  whiteNoise = -1;
  rms = 100e-9;

  if (displayCVSversion == 1) CVSdisplayVersion("analyticChol.C","plugin",CVS_verNum);

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: analyticChol\n");
  printf("Author:              G. Hobbs\n");
  printf("CVS Version:         $Revision: 1.9 $\n");
  printf(" --- type 'h' for help information\n");

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[0],argv[i+1]); 
	  strcpy(timFile[0],argv[i+2]);
	}
      else if (strcmp(argv[i],"-fc")==0)
	sscanf(argv[++i],"%lf",&modelFc);
      else if (strcmp(argv[i],"-gwamp")==0)
	sscanf(argv[++i],"%lf",&gwamp);
      else if (strcmp(argv[i],"-amp")==0)
	sscanf(argv[++i],"%lf",&modelA);
      else if (strcmp(argv[i],"-rms")==0)
	sscanf(argv[++i],"%lf",&rms);
      else if (strcmp(argv[i],"-alpha_res")==0)
	{
	  sscanf(argv[++i],"%lf",&modelAlpha);
	  modelAlpha = -modelAlpha;
	}
    }

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);


  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      if (i==0) doFit(psr,*npsr,0);   /* Do the fitting     */
      else textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
    }
  for (i=0;i<psr[0].nobs;i++)
    {
      resx[i] = (double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
      resy[i] = (double)(psr[0].obsn[i].residual);
      rese[i] = (double)(psr[0].obsn[i].toaErr*1.0e-6);
    }
  np=psr[0].nobs;

  if (modelFc == -1)
    {
      int ndays = ceil(resx[np-1]-resx[0]);
      modelFc = 0.5/(ndays/365.25);
    }
  if (modelA == -1)
    {
	//modelA = gwamp*gwamp/12.0/M_PI/M_PI*pow(1.0+pow(1.0/modelFc,2),modelAlpha/2.0);
      modelA = gwamp*gwamp/12.0/M_PI/M_PI*pow(modelFc,-modelAlpha);
    }
  if (whiteNoise == -1)
    {
      double tspan;
      int npts;
      
      tspan = (resx[np-1]-resx[0])/365.25;
      npts = np;
      whiteNoise = pow(rms/86400.0/365.25,2)*2*tspan/npts;
    }

  ndays = ceil((resx[np-1])-(resx[0])); 
  T2calculateCovarFunc(modelAlpha,modelFc,modelA,0,0,covFunc,resx,resy,rese,np);

  // Write it out
  sprintf(fname,"covarFunc.dat_%s",psr[0].name);
  fout = fopen(fname,"w");
  fprintf(fout,"1\n");
  for (i=0;i<=ndays;i++)
    fprintf(fout,"%.15g\n",covFunc[i]);
  fclose(fout);

  sprintf(fname,"%s.model",psr[0].name);
  fout = fopen(fname,"w");
  fprintf(fout,"MODEL 1\n");
  fprintf(fout,"ALPHA %g\n",modelAlpha);
  fprintf(fout,"FC %g\n",modelFc);
  fprintf(fout,"AMP %g\n",modelA);
  fprintf(fout,"WHITENOISE %g\n",whiteNoise);
  return 0;
}




char * plugVersionCheck = TEMPO2_h_VER;
