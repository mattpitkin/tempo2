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
#include "T2toolkit.h"

using namespace std;


void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,p,j,it,nit;
  double globalParameter;
  long seed=TKsetSeed();
  double mean;
  double ra,dec,raStep;
  int ndec,nra;
  int addWhite=0;
  const char *CVS_verNum = "$Revision: 1.1 $";
  long double **sat0;
  double gwm_amp=0.0;
  double gwm_phi;
  int count=0;
  int setPhi=0;
  char fname[128];

  FILE *fout;

  nit = 1;

  if (displayCVSversion == 1) CVSdisplayVersion((char *)"gwm.C",(char *)"plugin",CVS_verNum);

  *npsr = 0; 

  printf("Graphical Interface: gwm\n");
  printf("Author:              J. Wang, G. Hobbs\n");
  printf("CVS Version:         $Revision: 1.1 $\n");
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
      else if (strcmp(argv[i],"-addWhite")==0)
	addWhite=1;
      else if (strcmp(argv[i],"-amp")==0)
	sscanf(argv[++i],"%lf",&gwm_amp);
      else if (strcmp(argv[i],"-nit")==0)
	sscanf(argv[++i],"%d",&nit);
      else if (strcmp(argv[i],"-phi")==0)
	{
	  sscanf(argv[++i],"%lf",&gwm_phi);
	  setPhi=1;
	}
    }


  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

  // Allocate memory for ideal sats
  if (addWhite==1)
    {
      sat0 = (long double **)malloc(sizeof(long double*)*(*npsr));
      for (p=0;p<*npsr;p++)
	{
	  sat0[p] = (long double *)malloc(sizeof(long double)*(psr[p].nobs));
	  for (i=0;i<psr[p].nobs;i++)
	    sat0[p][i] = psr[p].obsn[i].sat;
	}
    }
  ndec = 40;
  nra  = 40;

  fout = fopen("output.dat","w");
  // List the pulsars
  fprintf(fout,"%d\n",*npsr);
  for (i=0;i<*npsr;i++)
    fprintf(fout,"%g %g %s\n",(double)psr[i].param[param_raj].val[0],(double)psr[i].param[param_decj].val[0],psr[i].name);
  for (dec = -M_PI/2.0;dec < M_PI/2;dec+=M_PI/ndec)
    {
      if (cos(dec)!=0.0)
	raStep = 2*M_PI/nra/cos(dec);
      else
	raStep = 2*M_PI;
      for (ra = 0;ra < 2*M_PI;ra += raStep)
	{
	  mean = 0.0;
	  //	  ra = 2.6915;
	  //	  dec = -1.33518; // 1.1362
	  for (it = 0;it <nit ;it++)
	    {
	      if (setPhi==0)
		gwm_phi = TKranDev(&seed)*2.0*M_PI;
	      //	      gwm_phi = 0.2;
	       for (p=0;p<*npsr;p++)
		 {
		   psr[p].param[param_gwm_amp].val[0] = gwm_amp;
		   psr[p].param[param_gwm_amp].fitFlag[0] = 2;
		   psr[p].param[param_gwm_amp].paramSet[0] = 1;
		   psr[p].gwm_raj = ra;
		   psr[p].gwm_decj = dec;
		   psr[p].gwm_epoch = 52187;
		   psr[p].gwm_phi = gwm_phi;
		 }
	       // Simulate ideal arrival times		   
	       if (addWhite==1)
		 {
		   int wnit=2,wit;
		   for (p=0;p<*npsr;p++)
		     {
		       for (i=0;i<psr[p].nobs;i++)
			 psr[p].obsn[i].sat = sat0[p][i];// -psr[p].obsn[i].residual/86400.0; // + (TKgaussDev(&seed)*(psr[p].obsn[i].toaErr*1e-6)/86400.0L);
		     }		 
		   for (wit = 0;wit<wnit;wit++)
		     {
		       formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
		       formResiduals(psr,*npsr,1);     /* Form the residuals                 */
		       for (p=0;p<*npsr;p++)
			 {
			   for (i=0;i<psr[p].nobs;i++)
			     psr[p].obsn[i].sat -= psr[p].obsn[i].residual/86400.0; // + (TKgaussDev(&seed)*(psr[p].obsn[i].toaErr*1e-6)/86400.0L);
			   
			 }
		     }
		   for (p=0;p<*npsr;p++)
		     {
		       for (i=0;i<psr[p].nobs;i++)
			 psr[p].obsn[i].sat += (TKgaussDev(&seed)*(psr[p].obsn[i].toaErr*1e-6)/86400.0L);
		       psr[p].param[param_gwm_amp].val[0] = 0.0;		       
		     }
		 }
	       
	       //	       sprintf(fname,"try%d.tim",count);
	       //	       count++;
	       //	       writeTim(fname,psr,"tempo2");
	       //	       exit(1);
	       formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
	       formResiduals(psr,*npsr,1);     /* Form the residuals                 */
	       doFit(psr,*npsr,0);   /* Do the fitting     */
	       //	       	       textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
	       //	       exit(1);
	       printf("Significance = %g %g %g\n",(double)psr[0].param[param_gwm_amp].val[0],
		      (double)psr[0].param[param_gwm_amp].err[0],(double)(psr[0].param[param_gwm_amp].val[0]/psr[0].param[param_gwm_amp].err[0]));
	       mean += (double)(psr[0].param[param_gwm_amp].val[0]/psr[0].param[param_gwm_amp].err[0]);
	       for (p=0;p<*npsr;p++)
		 {
		   for (i=0;i<MAX_PARAMS;i++)
		     {
		       for (j=0;j<psr[p].param[i].aSize;j++)
			 {
			   if (psr[p].param[i].paramSet[j]==1)
			     {
			       if (psr[p].param[i].fitFlag[j]==1)
				 psr[p].param[i].val[j] = psr[p].param[i].prefit[j];
			     }
			 }
		     }
		   psr[p].param[param_gwm_amp].paramSet[0]=0;
		 }
	    }
	  mean /= (double)nit;
	  printf("Mean significance = %g\n",mean);
	  fprintf(fout,"%g %g %g 0\n",ra,dec,mean);
	  fflush(fout);
	}
    }
  fclose(fout);
  if (addWhite==1)
    {
      for (p=0;p<*npsr;p++)
	free(sat0[p]);
      free(sat0);
    }
  return 0;
}

char * plugVersionCheck = (char *)TEMPO2_h_VER;
