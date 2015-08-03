#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
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
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include "tempo2.h"


void printplugs(bool full);

/* ******************************************** */
/* getInputs                                    */
/* Author:  G. Hobbs (02 May 2003)              */
/* Purpose: Parse the command line              */
/* Inputs:  argc, argv - command line arguments */
/* Outputs: File path for .tim and .par files   */
/*                                              */
/* Notes:                                       */
/*                                              */
/* Changes:                                     */
/* ******************************************** */

void getInputs(pulsar *psr,int argc, char *argv[],char timFile[][MAX_FILELEN],
	       char parFile[][MAX_FILELEN],int *list,int *npsr,
	       int *nGlobal,int *outRes,int *writeModel,char *outputSO,
	       int *polyco, char *polyco_args, char *polyco_file,
	       int *newpar,int *onlypre,char *dcmFile,char *covarFuncFile,char* newparname)
{
  int i,p;
  int gr=0;
  int parfile_num = 0;  /* Have we got a parfile? */
  int timfile_num = 0;  /* Have we got a timfile? */
  int gotTim=0;
  *list = 0;  /* Don't list parameters */
  const char *CVS_verNum = "$Revision: 1.17 $";

  if (displayCVSversion == 1) CVSdisplayVersion("getInputs.C","getInputs()",CVS_verNum);
  //  *nGlobal=0; /* How many global parameters are we fitting? */


  if (argc==2 && strcasecmp(argv[1],"-allParTim")!=0) /* Just have .tim file name */
    {
      if (strcmp(argv[1],"-h")==0) /* Some help */
	{
	  printf("\n");
	  printf("%s\n\n",PACKAGE_STRING);
	  printf("examples: %s mytim.tim\n",argv[0]);
	  printf("          %s -f mypar.par mytim.tim\n",argv[0]);
	  printf("          %s -gr plk -f mypar.par mytim.tim\n",argv[0]);
	  printf("\n");
	  printf("Options: \n\n");

	  printf("-epoch centre     Centres the PEPOCH in the fit\n");
	  printf("-f parFile        Selects parameter file\n");
	  printf("-dcm dcmFile      Data covariance matrix file\n");
	  printf("-gr name          Uses 'name' plugin for graphical interface\n");
	  printf("-h                This help\n");
	  printf("-list             Provides listing of clock corrections and residuals etc.\n");
	  printf("-output name      Uses 'name' plugin for output format\n");
	  printf("-pred \"args\"    Creates a predictive 2D Chebyshev polynomial.\n");
	  printf("      args = \"sitename mjd1 mjd2 freq1 freq2 ntimecoeff nfreqcoeff seg_length (s)\"\n"); 
	  printf("-polyco \"args\"  Creates a TEMPO1-style polyco file.\n");
	  printf("                   args = \"mjd1 mjd2 nspan ncoeff maxha sitename freq\"\n");
      printf("-polyco_file      Specify a leading string for file outputs.\n");
	  printf("-residuals        Outputs the residuals\n");
	  printf("-allInfo          Prints out clock, Earth orientation and similar information\n");
	  printf("-reminder         Saves the command line to T2command.input for future reference.\n");
	  printf("-norescale        Do not rescale parameter uncertainties by the sqrt(red. chisq)\n");
	  printf("-displayVersion   Display detailed CVS version number of every file used.\n");
	  printf("-v                Print verson number.\n");
	  printf("\n\n");
	  printf("Available plugins\n");
	  printplugs(false);
//	  system("ls $TEMPO2/plugins/ | grep plug | sed s/\"_\"/\" \"/ | awk '{print \"  - \" $1}' | sort | uniq");
	  printf("-----------------\n");
	  exit(1);
	}	    
      strcpy(timFile[timfile_num],argv[1]);
      strcpy(parFile[parfile_num],argv[1]);
      strcpy(parFile[parfile_num]+strlen(parFile[parfile_num])-3,"par");
      timfile_num++; parfile_num++;
      (*npsr)++;
    }
  else /* Have multiple command line arguments */
    {
      for (i=1;i<argc;i++)
	{
	  if (strcmp(argv[i],"-f")==0) /* Have .par file */
	    {
	      strcpy(parFile[parfile_num],argv[++i]);
	      if (argv[i+1][0]!='-') /* Must be tim file */
		{
		  strcpy(timFile[timfile_num],argv[++i]);
		  timfile_num++;
		  gotTim=1;
		}
	      /*	      i+=2; */
	      parfile_num++;
	    }
	  else if (strcasecmp(argv[i],"-allParTim")==0)
	    {
	      FILE *pin;
	      char str[1000];

	      // Load in all available par and tim files
	      printf("Using all available .par and .tim files\n");
	      sprintf(str,"ls `ls *.par | sed s/par/tim/` | sed s/.tim/\"\"/");
	      pin = popen(str,"r");
	      while (!feof(pin))
		{
		  if (fscanf(pin,"%s",str)==1)
		    {
		      sprintf(parFile[*npsr],"%s.par",str);
		      sprintf(timFile[*npsr],"%s.tim",str);
		      (*npsr)++;
		      parfile_num++;
		      timfile_num++;
		      gotTim=1;

		    }
		}
	      pclose(pin);
	      printf("Obtained files for %d pulsars\n",*npsr);

	    }
	  else if (strcmp(argv[i],"-fitfunc")==0)
	    strcpy(psr[0].fitFunc,argv[++i]);
	  else if (strcasecmp(argv[i],"-norescale")==0)
	    {
	      for (p=0;p<MAX_PSR;p++)
		psr[p].rescaleErrChisq=0;
	    }
	  else if (strcmp(argv[i],"-gr")==0)
	    gr=1;
	  else if (strcmp(argv[i],"-dcm")==0)
	    strcpy(dcmFile,argv[++i]);
	  else if ((strcmp(argv[i],"-dcf")==0) || (strcmp(argv[i],"-chol")==0))
	    strcpy(covarFuncFile,argv[++i]);
	  else if (strcmp(argv[i],"-filter")==0)
	    {
	      strcat(psr[0].filterStr,argv[++i]);
	      strcat(psr[0].filterStr," ");
	    }
	  else if (strcmp(argv[i],"-pass")==0)
	    {
	      strcat(psr[0].passStr,argv[++i]);
	      strcat(psr[0].passStr," ");
	    }
	  else if (strcmp(argv[i],"-list")==0)
	    *list = 1;	  
	  else if (strcmp(argv[i],"-output")==0)
	      strcpy(outputSO,argv[++i]);
	  else if (strcmp(argv[i],"-residuals")==0)
	    *outRes = 1;
	  else if (strcmp(argv[i],"-del")==0)
	    strcpy(psr[0].deleteFileName,argv[++i]);	  
	  else if (strcmp(argv[i],"-model")==0)
	    *writeModel = 1;
	  else if (strcmp(argv[i],"-newpar")==0){
		 strcpy(newparname,"new.par");
	    *newpar=1;
	  }
	  else if (strcmp(argv[i],"-outpar")==0){
	    *newpar=1;
		strcpy(newparname,argv[++i]);
	  }
	  else if (strcmp(argv[i],"-pre")==0) /* Don't iterate to form post-fit residuals */
	    *onlypre = 1;
	  else if (strcmp(argv[i],"-pred")==0)
	    {
	      *polyco=2;
	      strcpy(polyco_args, argv[++i]);
	    }
	  else if (strcmp(argv[i],"-polyco")==0)
	    {
	      *polyco=1;
	      strcpy(polyco_args, argv[++i]);
	    }
	  else if (strcmp(argv[i],"-polyco_file")==0)
	    {
	      strcpy(polyco_file, argv[++i]);
	    }
	  else if (i==argc-1 && gotTim==0) /* Must be .tim file name */
	    {
	      strcpy(timFile[timfile_num],argv[i]);	    
	      //printf("Have %s\n",timFile[timfile_num]);
	      timfile_num++;
	    }
	}
      if (timfile_num==0 && *polyco==0 && gr==0)
	{
	  printf("ERROR [FILE1]: No .tim file given on command line\n"); 
	  exit(1);
	}      
      if (parfile_num==0 && gr==0)
	{
	  printf("ERROR [FILE2]: No .par file given on command line\n"); 
	  exit(1);
	}
      *npsr = parfile_num;
    }
}










void printplugs(bool full){
	char pname[MAX_STRLEN];
	char matchstr[MAX_STRLEN];
	char pname_list[MAX_STRLEN][MAX_STRLEN];
	int np=0;

	sprintf(matchstr,"%s_plug.t2",tempo2MachineType);
	for (int i=0; i < tempo2_plug_path_len; i++){
		printf("in '%s'\n",tempo2_plug_path[i]);
		struct dirent *pent = NULL;
		DIR* d = NULL;
		d=opendir(tempo2_plug_path[i]);
		if(d==NULL){
			printf("(dir not readable)\n");
			continue;
		}
		int count=0;
		while(pent = readdir(d)){
			char* name=pent->d_name;
			strcpy(pname,"");
			int j=0;
			char flag=' ';
			int len=strlen(name);
			while (j < len){
				if (name[j]=='_'){
					memcpy(pname,name,j);
					pname[j]='\0';
					break;
				}
				j++;
			}
			while (j < len){
				if (strcmp(name+j,matchstr)==0){
					for (int k=0; k < np; k++){
						if(strcmp(pname_list[k],pname)==0){
							flag='#';
							break;
						}
					}
					if(full || flag==' '){
						printf(" - %s%c\n",pname,flag);
						count++;
					}
					strcpy(pname_list[np++],pname);
					break;
				}
				j++;
			}
		}
		if (!count){
			printf("(none or all hidden)\n");
		}
	}
}


void setPlugPath(){
   int i;
  if (getenv("TEMPO2_PLUG_PATH")!=NULL){
    char *p_path = (char*)malloc(MAX_STRLEN*32);
    strcpy(p_path,getenv("TEMPO2_PLUG_PATH"));
    int len= strlen(p_path);
    for (i=0; i < len; i++){
      if (p_path[i] == ':')p_path[i]='\0';
    }
    i=0;
    while(i < len){
      strcpy(tempo2_plug_path[tempo2_plug_path_len++],p_path+i);
      i+=strlen(p_path+i)+1;
    }
    free(p_path);
  }
#ifdef TEMPO2_CONFIGURE_PLUG  
  // If tempo2 was compiled with a non-standard plugin path, add it to the default search path
  strcpy(tempo2_plug_path[tempo2_plug_path_len++],TEMPO2_CONFIGURE_PLUG);
#endif
  sprintf(tempo2_plug_path[tempo2_plug_path_len++],"%s/plugins/",getenv(TEMPO2_ENVIRON));

}
