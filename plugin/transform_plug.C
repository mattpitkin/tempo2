/* Template for a tempo2 plugin */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"

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
  char outputFile[1000];
  int back=0;

  *npsr = 1; 

  printf("Graphical Interface: transform\n");
  printf("Author:              George Hobbs\n");
  printf("Version:             1\n");
  printf("Transforms parameter files from old format to new format\n");
  
  if (argc!=5 && argc!=6)
    {
      printf("Usage tempo2 -gr transform parFile outputFile [back]\n");
      exit(1);
    }

  if (argc==6) back=1;

  strcpy(parFile[0],argv[3]);
  strcpy(outputFile,argv[4]);

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  psr[0].nobs = 0;
  preProcess(psr,*npsr,argc,argv);

  /* Now do the conversion */
  if (back==0) 
    {
      transform_units(&psr[0],TDB_UNITS,SI_UNITS);
      psr[0].param[param_ephver].val[0]=5;
      psr[0].units = SI_UNITS;
      psr[0].timeEphemeris = IF99_TIMEEPH;
      psr[0].t2cMethod = T2C_IAU2000B;
      psr[0].dilateFreq = 1;
      psr[0].correctTroposphere = 1;
      psr[0].planetShapiro = 1;
    }
  else 
    {
      transform_units(&psr[0],SI_UNITS,TDB_UNITS);   
      psr[0].param[param_ephver].val[0]=2;
    }

  /* Convert position to strings */
  turn_hms(psr[0].param[param_raj].val[0]/(2.0*M_PI), psr[0].rajStrPost);
  turn_dms(psr[0].param[param_decj].val[0]/(2.0*M_PI), psr[0].decjStrPost);

  /* Now write out the new parameter file */
  textOutput(psr,*npsr,0,0,0,1,outputFile);

  return 0;
}

