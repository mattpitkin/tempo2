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

using namespace std;


void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  char covarFuncFile[128];
  int i;
  double globalParameter;
  const char *CVS_verNum = "$Revision: 1.2 $";
  FILE *fout;

  strcpy(covarFuncFile,"NULL");

  if (displayCVSversion == 1) CVSdisplayVersion((char *)"findCW.C",(char *)"plugin",CVS_verNum);

  *npsr = 0;

  printf("Graphical Interface: findCW\n");
  printf("Author:              X. Zhu, G. Hobbs\n");
  printf("CVS Version:         $Revision: 1.2 $\n");
  printf(" --- type 'h' for help information\n");

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[*npsr],argv[++i]); 
	  strcpy(timFile[*npsr],argv[++i]);
	  (*npsr)++;
	}
      else if (strcmp(argv[i],"-dcf")==0)
	strcpy(covarFuncFile,argv[++i]);
    }

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

    for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
  //  i=0;
    {
      logdbg("Calling formBatsAll");
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      logdbg("Calling formResiduals");
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      logdbg("Calling doFit");
      if (i==0) doFitAll(psr,*npsr,covarFuncFile);   /* Do the fitting     */
      else textOutput(psr,*npsr,globalParameter,0,0,0,(char *)"");  /* Display the output */
    }

  // Print A+ and Ax into a file
  fout = fopen("aplus_across.dat","w");
  for (i=0;i<psr[0].quad_ifuncN_p;i++)
    {
      fprintf(fout,"%.2lf %g %g %g %g\n",psr[0].quad_ifuncT_p[i],psr[0].quad_ifuncV_p[i],psr[0].quad_ifuncE_p[i],psr[0].quad_ifuncV_c[i],psr[0].quad_ifuncE_c[i]);
    }
  fclose(fout);
  //  return 0;
  // Calculate a power spectrum
  {
    double specX[1024];
    double specY_ap[1024];
    double specY_re_ap[1024];
    double specY_im_ap[1024];
    double specY_ac[1024];
    double specY_re_ac[1024];
    double specY_im_ac[1024];
    double specY_pfitr[1024];
    double specY_r[1024];
    double dt;
    int    nSpec;
    double x[MAX_OBSN],y[MAX_OBSN],e[MAX_OBSN],y2[MAX_OBSN];
    int    n;
    double specY[1024];

    n = psr[0].nobs;
    for (i=0;i<n;i++)
      {
	x[i] = (double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
	y[i] = (double)(psr[0].obsn[i].prefitResidual);
	y2[i] = (double)(psr[0].obsn[i].residual);
	e[i] = (double)psr[0].obsn[i].toaErr*1e-6;
      }

    // Note: unweighted DFT
    TKspectrum(psr[0].quad_ifuncT_p,psr[0].quad_ifuncV_p,psr[0].quad_ifuncE_p,psr[0].quad_ifuncN_p,0,0,0,0,0,1,1,1,1,specX,specY_ap,&nSpec,0,0,specY_re_ap,specY_im_ap);
    TKspectrum(psr[0].quad_ifuncT_c,psr[0].quad_ifuncV_c,psr[0].quad_ifuncE_c,psr[0].quad_ifuncN_c,0,0,0,0,0,1,1,1,1,specX,specY_ac,&nSpec,0,0,specY_re_ac,specY_im_ac);
    //
    // Assuming that the ifuncs are regularly sampled
    //
    dt = psr[0].quad_ifuncT_p[1]-psr[0].quad_ifuncT_p[0];
    for (i=0;i<nSpec;i++)
      {
	specY_ap[i] = specY_ap[i]*pow(sin(M_PI*dt*specX[i])/(M_PI*dt*specX[i]),4);
	specY_ac[i] = specY_ac[i]*pow(sin(M_PI*dt*specX[i])/(M_PI*dt*specX[i]),4);
      }

    fout = fopen("spectrum.dat","w");
    for (i=0;i<nSpec;i++)
      specY[i] = specY_ap[i]+specY_ac[i];

    for (i=0;i<nSpec;i++)
      fprintf(fout,"%g %g %g %g\n",specX[i],specY_ap[i],specY_ac[i],specY[i]);
    fclose(fout);
  }
  return 0;
}

char * plugVersionCheck = TEMPO2_h_VER;
