/* Plugin to add red noise to a data set */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "GWsim.h"

using namespace std;


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
  long double a;
  long double alpha,toffset;
  long double kp[3];
  long double flo,fhi;
  long double res[MAX_OBSN],mean;
  double dist;
  int addGW=0,ngw,k;
  int j;
  long seed = TKsetSeed();
  double addwLevel=0;
  double varywLevel=0;
  double whitet,whitett;
  gwSrc *gw;
  long double amp;
  


  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: addRed\n");
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
      else if (strcmp(argv[i],"-amp")==0)
	sscanf(argv[++i],"%Lf",&amp);
      else if (strcmp(argv[i],"-alpha")==0)
	sscanf(argv[++i],"%Lf",&alpha);
      else if (strcmp(argv[i],"-addw")==0) // Add white
	sscanf(argv[++i],"%lf",&addwLevel);
      else if (strcmp(argv[i],"-varyw")==0) // Let white noise vary
	sscanf(argv[++i],"%lf",&varywLevel);
    }

  alpha=(alpha+3)/2.0L; // Convert to a gravitational wave alpha

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

  a = (long double)amp*pow(86400.0*365.25,alpha);
  dist =  3.08568025e19; // 1 kpc in m
  setupPulsar_GWsim(psr[0].param[param_raj].val[0],
		    psr[0].param[param_decj].val[0],kp);
  flo = 1.0L/(30*365.25*86400.0L);
  fhi = 1.0L/(2.0*86400.0L);
  ngw = 1000;
  if((gw = (gwSrc *)malloc(sizeof(gwSrc)*ngw))==NULL){
    printf("Unable to allocate memory for gwSrc.\n");
    exit(1);
  }
  toffset = psr[0].param[param_pepoch].val[0];
  GWbackground(gw,ngw,&seed,flo,fhi,a,alpha,1);
  mean=0.0L;
  //printf("Calc residuals 1\n");
  for (j=0;j<psr[0].nobs;j++)
    {
      res[j]=0.0L;
      for (k=0;k<ngw;k++)
	res[j]+=calculateResidualGW(kp,&gw[k],
				    (psr[0].obsn[j].sat-toffset)*86400.0L,
				    dist);	  
      mean+=res[j];
    }
  for (j=0;j<psr[0].nobs;j++)
    {
      psr[0].obsn[j].sat+=(res[j]-mean/psr[0].nobs)/86400.0L;
      if (addwLevel>0)
	{
	  whitett = addwLevel;
	  if (varywLevel>0)
	    {
	      whitett = sqrt(pow(addwLevel+TKgaussDev(&seed)*varywLevel,2)+pow(100e-9,2));
	    }
	  whitet  =(TKgaussDev(&seed)*whitett);
	  psr[0].obsn[j].sat+=(whitet)/86400.0L;
	  psr[0].obsn[j].toaErr = whitett/1.0e-6;
	}
    }
  readParfile(psr,parFile,timFile,1); /* Load the parameters                */
  preProcess(psr,1,argc,argv);
  formBatsAll(psr,1);                 /* Form the barycentric arrival times */
  formResiduals(psr,1,0);             /* Form the residuals                 */
  doFit(psr,1,0);                     /* Do the fitting                     */ 
  formBatsAll(psr,1);                 /* Form the barycentric arrival times */
  formResiduals(psr,1,0);             /* Form the residuals                 */
  textOutput(psr,1,0,0,0,1,"new.par");
  writeTim("new.tim",psr,"tempo2");
  
  return 0;
}

char * plugVersionCheck = TEMPO2_h_VER;
