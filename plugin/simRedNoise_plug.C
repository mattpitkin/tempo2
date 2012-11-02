// Plugin to simulate red noise data sets with real sampling and TOA uncertainties
//
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

/* Template for a tempo2 plugin */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "T2toolkit.h"
#include "TKfit.h"
#include "TKspectrum.h"

using namespace std;

void doPlugin(pulsar *psr,int npsr,double amp,double alpha,double fc,int removeQuad);
void getRedNoiseRealisation(pulsar psr,double amp,double alpha,double fc,long *seed,double *redNoise,int *nRedNoise,double *minx,double *delta);

void help() /* Display help */
{
  /* This function should contain usage information about the plugin which should (in general) be accessed */
  /* by the user pressing 'h'                                                                              */
}


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i;
  double globalParameter;
  double amp=0,alpha=0,fc = 0;
  int removeQuad=0;

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: simRedNoise\n");
  printf("Author:              G. Hobbs\n");
  printf("Version:             1.0\n");
  printf(" --- type 'h' for help information\n");


  /* Obtain the .par and the .tim file from the command line */
  if (argc==4) /* Only provided .tim name */
    {
      strcpy(timFile[0],argv[3]);
      strcpy(parFile[0],argv[3]);
      parFile[0][strlen(parFile[0])-3] = '\0';
      strcat(parFile[0],"par");
    }

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[0],argv[i+1]); 
	  strcpy(timFile[0],argv[i+2]);
	}
      else if (strcmp(argv[i],"-alpha")==0)
	sscanf(argv[++i],"%lf",&alpha);
      else if (strcmp(argv[i],"-fc")==0)
	sscanf(argv[++i],"%lf",&fc);
      else if (strcmp(argv[i],"-a")==0)
	sscanf(argv[++i],"%lf",&amp);
      else if (strcmp(argv[i],"-removeQuad")==0)
	removeQuad=1;
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

  doPlugin(psr,*npsr,amp,alpha,fc,removeQuad);

  return 0;
}

void doPlugin(pulsar *psr,int npsr,double amp,double alpha,double fc,int removeQuad)
{
  int i,p,j;
  int nit=2;
  long double sat0[MAX_OBSN];
  long seed = TKsetSeed();
  char fname[100];
  int addToaError=1;
  long origSeed = seed;
  //  double amp=5e-25,alpha=3.5,fc = 0.09;

  //   double amp=5e-31,alpha=2.0,fc = 0.5;

  // double amp=0,alpha=3.5,fc = 0.09;
  double redNoise[MAX_OBSN],delta;
  int nRedNoise;
  int iclosest;
  double minx;

  printf("Seed = %d\n",origSeed);

  // Form idealised site arrival times
  for (j=0;j<nit;j++)
    {
      for (p=0;p<npsr;p++)
	{
	  for (i=0;i<psr[p].nobs;i++)
	    psr[p].obsn[i].sat -= (long double)psr[p].obsn[i].residual/SECDAY;
	}
      formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,npsr,1);     /* Form the residuals                 */
    }
  
  // Make new data sets
  for (p=0;p<npsr;p++)
    {
      // 
      if (amp != 0.0) // Add in red noise based on spectral model
	{
	  getRedNoiseRealisation(psr[p],amp,alpha,fc,&seed,redNoise,&nRedNoise,&minx,&delta);
	  if (removeQuad==1)
	    {
	      double xv[nRedNoise];
	      for (i=0;i<nRedNoise;i++)
		xv[i] = i;
	      TKremovePoly_d(xv,redNoise,nRedNoise,3);
	    }
	}
      for (i=0;i<psr[p].nobs;i++)
	{
	  // Add on TOA error bar
	  if (addToaError==1)
	    psr[p].obsn[i].sat += TKgaussDev(&seed)*psr[p].obsn[i].toaErr*1.0e-6/SECDAY;
	  if (amp!=0)
	    {

	      iclosest = (int)((double)(psr[p].obsn[i].sat-minx)/delta);
	      if (iclosest == nRedNoise) iclosest = nRedNoise-1;
	      	      printf("Have %g %g %g %d %g\n",(double)psr[p].obsn[i].sat,minx,delta,iclosest,redNoise[iclosest]);
	      psr[p].obsn[i].sat += redNoise[iclosest]/SECDAY;
	    }
	}
      // Write out new .tim file
      sprintf(fname,"simulate_%s.tim",psr[p].name);
      writeTim(fname,psr+p,"tempo2");
    }


}

void getRedNoiseRealisation(pulsar psr,double amp,double alpha,double fc,long *seed,double *redNoise,int *nRedNoise,double *minx,double *delta)
{
  double dspan;
  double re[MAX_OBSN],im[MAX_OBSN];
  double freq,psd,theta;
  int close2n;
  int i;
  double maxx;
  
  *minx=(double)psr.obsn[0].sat;
  maxx=*minx;

  for (i=1;i<psr.nobs;i++)
    {
      if (maxx < (double)psr.obsn[i].sat) maxx = (double)psr.obsn[i].sat;
      if (*minx > (double)psr.obsn[i].sat) *minx = (double)psr.obsn[i].sat;
    }
  dspan = maxx-(*minx);
  printf("Data span = %g (days)\n",dspan);
  close2n = (int)(log(dspan)/log(2.0)+1);
  printf("Closest 2**N number = %d giving %g\n",close2n,pow(2,close2n));
  // Maybe consider simulating e.g. 10x longer and then truncating so don't have to mess around with sub harmonics

  *delta = dspan/pow(2,close2n);
  for (i=0;i<(int)pow(2,close2n);i++)
    {
      //      freq = i/2.0/(*delta); // Units of d^-1
      //      freq = i/(*delta); // Units of d^-1 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      freq = i/dspan;
      //     psd  = amp*pow(365.25*86400.0,3)/pow(1.0+pow(freq*365.25/fc,2),alpha/2.0);
      //      psd = amp*pow(365.25*86400.0,2)/2.0/(dspan/365.25)*pow(2,close2n)*pow(2,close2n)*pow(2,close2n)*M_PI/pow(1.0+pow(freq*365.25/fc,2),alpha/2.0);
      //      psd = amp*pow(365.25*86400.0,2)/(dspan/(365.25*86400.0))/pow(1.0+pow(freq*365.25/fc,2),alpha/2.0);
      //      psd = amp/pow(1.0+pow(freq*365.25/fc,2),alpha/2.0);

      if (i==0)
	psd=0.0;
      else
	psd = amp/pow(1.0+pow(freq*365.25/fc,2),alpha/2.0); // model in yr^3
      printf("%d amp = %g, psd = %g, freq = %g, alpha = %g, fc = %g\n",i,amp,psd,freq,alpha,fc);

      if (i<=((int)pow(2,close2n)/2))
	{
	  	  re[i] = 0.5*sqrt(psd/(dspan/365.25)*pow(365.25*86400.0,2))*TKgaussDev(seed);
	  	  im[i] = 0.5*sqrt(psd/(dspan/365.25)*pow(365.25*86400.0,2))*TKgaussDev(seed);
	  //	  re[i] = sqrt(psd/2.0/(dspan/365.25)*pow(365.25*86400.0,2))*TKgaussDev(seed);
	  //	  im[i] = sqrt(psd/2.0/(dspan/365.25)*pow(365.25*86400.0,2))*TKgaussDev(seed);
	  printf("Real val = %g %g %g  %g %d %g %g %g\n",psd,*delta,(double)dspan,freq*365.25,(int)pow(2,close2n),1.0/pow(1.0+pow(freq*365.25/fc,2),alpha/2.0),re[i],im[i]);
	}
      else
	{
	  re[i] = re[(int)pow(2,close2n)-i];
	  im[i] = -1.0*im[(int)pow(2,close2n)-i];
	}
      printf("spectrum = %g %lg %lg %lg %g\n",freq,psd,re[i],im[i],re[i]*re[i]+im[i]*im[i]);
    }
  TK_fft(-1,(int)pow(2,close2n),re,im);
  *nRedNoise = (int)pow(2,close2n);
  for (i=0;i<(int)pow(2,close2n);i++)
    {
      //      redNoise[i] = re[i]*pow(86400.0*365.25,2)*2*365.25/dspan; ///1e25;      
      //      redNoise[i] = re[i]/pow(86400.0*365.25,2)/2/365.25*dspan; ///1e25;      
      redNoise[i] = re[i]; //*86400.0*365.25; //*86400.0*365.25/pow(2,close2n);
                  printf("output %d %g %g %g\n",i,re[i],im[i],redNoise[i]);
    }

}
char * plugVersionCheck = TEMPO2_h_VER;
