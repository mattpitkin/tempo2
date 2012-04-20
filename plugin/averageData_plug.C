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

using namespace std;

#define MAX_TIMES 1000

void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,j,k,argn=0;
  double globalParameter;
  long double centreMJD;
  long double newTOA;
  const char *CVS_verNum = "$Revision$";
  char timeList[MAX_STRLEN];
  char parFileName[MAX_TIMES][MAX_STRLEN];
  char timFileName[MAX_TIMES][MAX_STRLEN];
  int nStride=0;
  FILE *fin;
  int nread;
  long double mjd1[MAX_TIMES],mjd2[MAX_TIMES];
  char tname[100];
  FILE *fout;
  FILE *fout2;
  FILE *fout3;
  char str[1024];
  int closeID,t;
  double distance=0;

  if (displayCVSversion == 1) CVSdisplayVersion("averageData.C","plugin",CVS_verNum);

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: averageData\n");
  printf("Author:              G. Hobbs, R. Manchester\n");
  printf("CVS Version:         $Revision$\n");
  printf(" --- type 'h' for help information\n");

  printf("Plugin to average data to produce a new .tim file\n");
  printf("Algorithm: \n");
  printf(" 1. Update PEPOCH to be at the centre of the data set\n");
  printf(" 2. Fit only for F0 and phase\n");

  fout = fopen("avpts.tim","w");
  fprintf(fout,"FORMAT 1\n");
  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  argn=i;
	  strcpy(parFile[0],argv[++i]); 
	  strcpy(timFile[0],argv[++i]);
	}
      else if (strcmp(argv[i],"-t")==0)
	strcpy(timeList,argv[++i]);
    }
  // Read the list of strides
  fin = fopen(timeList,"r");
  while (!feof(fin))
    {
      nread = fscanf(fin,"%Lf %Lf %s %s",&mjd1[nStride],&mjd2[nStride],parFileName[nStride],tname);
      if (nread==3 || nread==4)
	{
	  if (nread==4)
	    strcpy(timFileName[nStride],tname);
	  else
	    strcpy(timFileName[nStride],timFile[0]);
	  nStride++;
	}
    }
  fclose(fin);
  printf("Have read %d strides\n",nStride);
  if (nStride == 0)
    {
      printf("ERROR: require at least one stride - use -t option\n");
      exit(1);
    }
  for (i=0;i<nStride;i++)
    {
      printf("***** Stride %d *****\n",i);
      sprintf(str,"stride%d",i);
      // Find the centre of the data span
      centreMJD = 0.5*(mjd1[i]+mjd2[i]);
      strcpy(parFile[0],parFileName[i]);
      strcpy(timFile[0],timFileName[i]);

      psr[0].nJumps=0;
      psr[0].nconstraints = 0;
      for(j=0;j<MAX_PARAMS;j++){
	psr[0].param[j].nLinkTo = 0;
	psr[0].param[j].nLinkFrom = 0;
	for (k=0;k<psr[0].param[j].aSize;k++)
	  {
	    psr[0].param[j].paramSet[k] = 0;
	    psr[0].param[j].fitFlag[k] = 0;
	    psr[0].param[j].prefit[k] = 0.0;
	    psr[0].param[j].val[k] = 0.0;
	  }
      }

      readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
      readTimfile(psr,timFile,*npsr); /* Load the arrival times    */

      // Make a new timFile that only contains the TOAs in the region
      fout3 = fopen("short.tim","w");
      fprintf(fout3,"FORMAT 1\n");
      for (k=0;k<psr[0].nobs;k++)
	{
	  if (psr[0].obsn[k].sat >= mjd1[i] && psr[0].obsn[k].sat < mjd2[i])
	    {
	      fprintf(fout3,"%s %.8f %.17Lf %.5f %s\n",psr[0].obsn[k].fname,
		      psr[0].obsn[k].freq,psr[0].obsn[k].sat,
		      psr[0].obsn[k].toaErr,psr[0].obsn[k].telID);
	    }
	}	  
      fclose(fout3);
      strcpy(timFile[0],"short.tim");
      readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
      
      
      // Update the epoch in the par file for centreMJD
      strcpy(argv[argn],"-epoch");
      sprintf(argv[argn+1],"%.5f",(double)centreMJD);
      printf("Before preProcess\n");
      preProcess(psr,1,argc,argv);      
      printf("After preProcess\n");
      // Turn off all fitting
      for (j=0;j<MAX_PARAMS;j++)
	{
	  for (k=0;k<psr[0].param[j].aSize;k++)
	    psr[0].param[j].fitFlag[k] = 0;
	}
      
      // Update the start and finish flags
      psr[0].param[param_start].val[0] = mjd1[i];
      psr[0].param[param_start].fitFlag[0] = 1;
      psr[0].param[param_start].paramSet[0] = 1;
      
      psr[0].param[param_finish].val[0] = mjd2[i];
      psr[0].param[param_finish].fitFlag[0] = 1;
      psr[0].param[param_finish].paramSet[0] = 1;
    
      // Find closest point
      t=1;
      distance=-1;
      for (k=0;k<psr[0].nobs;k++)
	{
	  if (psr[0].obsn[k].sat > psr[0].param[param_start].val[0]
	      && psr[0].obsn[k].sat < psr[0].param[param_finish].val[0])
	    {
	      if (t==1)
		{
		  closeID = k;
		  distance = (double)fabs(psr[0].obsn[k].sat - centreMJD);
		  t=2;
		}
	      else
		{
		  if ((double)fabs(psr[0].obsn[k].sat - centreMJD) < distance)
		    {
		      closeID=k;
		      distance = (double)fabs(psr[0].obsn[i].sat - centreMJD);
		    }
		}
	    }
	  
	}
      if (distance == -1)
	{
	  printf("Cannot identify any points within the current region: mjd = %g\n",(double)centreMJD);
	  exit(1);
	}
      
      // Turn on required fitting
      //      psr[0].param[param_f].fitFlag[0] = 1;
      
      for (k=0;k<2;k++)                   /* Do two iterations for pre- and post-fit residuals*/
	{
	  formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
	  formResiduals(psr,*npsr,0);    /* Form the residuals                 */
	  if (k==0) 
	    {
	      char str2[1024];
	      sprintf(str2,"%s.res",str);
	      fout2 = fopen(str2,"w");
	      for (j=0;j<psr[0].nobs;j++)
		fprintf(fout2,"%g %g %g\n",(double)(psr[0].obsn[j].sat-centreMJD),(double)psr[0].obsn[j].residual,psr[0].obsn[j].toaErr*1.0e-6);
	      fclose(fout2);
	      doFit(psr,*npsr,0);   /* Do the fitting     */
	    }
	  else textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
	}
      // Now fit once more to the post-fit residuals
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,0);    /* Form the residuals                 */
      doFit(psr,*npsr,0);   /* Do the fitting     */

      printf("Offset = %g, offset_e = %g\n",psr[0].offset,psr[0].offset_e);

      //      for (i=0;i<psr[0].nobs;i++)
      //	printf("res1: %g %g %g\n",(double)(psr[0].obsn[i].sat-centreMJD),(double)psr[0].obsn[i].residual-psr[0].offset,(double)psr[0].obsn[i].toaErr*1.0e-6);
      
      // Now add in a pseudo point at centreMJD
      psr[0].obsn[psr[0].nobs].sat  = centreMJD;
      psr[0].obsn[psr[0].nobs].freq = psr[0].obsn[closeID].freq;
      strcpy(psr[0].obsn[psr[0].nobs].fname,"avpt");
      strcpy(psr[0].obsn[psr[0].nobs].telID,psr[0].obsn[closeID].telID);
      psr[0].obsn[psr[0].nobs].phaseOffset=0.0;
      psr[0].obsn[psr[0].nobs].deleted=0.0;
      psr[0].obsn[psr[0].nobs].toaErr=psr[0].offset_e/1.0e-6;
      psr[0].obsn[psr[0].nobs].clockCorr=1;
      psr[0].obsn[psr[0].nobs].delayCorr=1;
      psr[0].obsn[psr[0].nobs].efac=1;
      

      psr[0].nobs++;
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,0);    /* Form the residuals                 */
      //      for (i=0;i<psr[0].nobs;i++)
      //      	printf("res2: %g %g %g\n",(double)(psr[0].obsn[i].sat-centreMJD),(double)psr[0].obsn[i].residual-psr[0].offset,(double)psr[0].obsn[i].toaErr*1.0e-6);
      
      for (j=0;j<3;j++) // Iterate to converge
	{
	  psr[0].obsn[psr[0].nobs-1].sat  -= (long double)(psr[0].obsn[psr[0].nobs-1].residual-psr[0].offset)/SECDAY;
	  formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
	  formResiduals(psr,*npsr,0);    /* Form the residuals                 */
	}
      //      psr[0].obsn[psr[0].nobs-1].sat += (long double)psr[0].offset/SECDAY;
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,0);    /* Form the residuals                 */
      if (psr[0].nobs==2) // Only dealing with one point
	{ 
	  fprintf(fout,"%s %.8f %.17Lf %.5f %s -av 1\n",psr[0].obsn[psr[0].nobs-2].fname,
		  psr[0].obsn[psr[0].nobs-2].freq,psr[0].obsn[psr[0].nobs-2].sat,
		  psr[0].obsn[psr[0].nobs-2].toaErr,psr[0].obsn[psr[0].nobs-2].telID);
	}
      else
	{
	  fprintf(fout,"%s %.8f %.17Lf %.5f %s -av 1 -nfit %d\n",psr[0].obsn[psr[0].nobs-1].fname,
		  psr[0].obsn[psr[0].nobs-1].freq,psr[0].obsn[psr[0].nobs-1].sat,
		  psr[0].obsn[psr[0].nobs-1].toaErr,psr[0].obsn[psr[0].nobs-1].telID,psr[0].nobs);
	}
      fout2 = fopen(str,"w");
      for (k=0;k<psr[0].nobs;k++)
	//	fprintf(fout2,"res3: %g %g %g\n",(double)(psr[0].obsn[k].sat-centreMJD),(double)psr[0].obsn[k].residual-psr[0].offset,(double)psr[0].obsn[k].toaErr*1.0e-6);
      fprintf(fout2,"res3: %g %g %g\n",(double)(psr[0].obsn[k].sat-centreMJD),(double)psr[0].obsn[k].residual,(double)psr[0].obsn[k].toaErr*1.0e-6);
      fflush(fout);
      fclose(fout2);
    }
  fclose(fout);
  return 0;
}

char * plugVersionCheck = TEMPO2_h_VER;
