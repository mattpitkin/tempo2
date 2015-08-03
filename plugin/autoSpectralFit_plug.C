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
#include <TKspectrum.h>
#include "T2toolkit.h"
#include "TKfit.h"
#include "choleskyRoutines.h"

using namespace std;

#define MAX_FREQ 10000


void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i;
  double globalParameter;
  double tsmooth;
  double expSmooth,fc_w,modelAlpha_out,modelVal,whiteNoiseLevel,fc_r;
  char covarFuncName[128];
  char modelName[128];

  const char *CVS_verNum = "$Revision: 1.4 $";

  if (displayCVSversion == 1) CVSdisplayVersion("autoSpectralFit_plug.C","plugin",CVS_verNum);

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: autoSpectralFit\n");
  printf("Author:              R. Shannon, G. Hobbs\n");
  printf("CVS Version:         $Revision: 1.4 $\n");
  printf(" --- type 'h' for help information\n");


  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[0],argv[i+1]); 
	  strcpy(timFile[0],argv[i+2]);
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

  sprintf(covarFuncName,"covarFunc.dat_%s",psr[0].name);
  T2get_covFunc_automatic(psr,0,covarFuncName,&fc_w,&fc_r,&modelAlpha_out,&modelVal,&whiteNoiseLevel,0,0);
  sprintf(modelName,"%s.model",psr[0].name);
  T2writeCovarFuncModel(modelAlpha_out,fc_r,modelVal,whiteNoiseLevel,modelName);
  printf("Complete with\n");
  printf("fc              = %g\n",fc_r);
  printf("modelAlpha      = %g\n",modelAlpha_out);
  printf("modelVal        = %g\n",modelVal);
  printf("whiteNoiseLevel = %g\n",whiteNoiseLevel);
  return 0;
}
