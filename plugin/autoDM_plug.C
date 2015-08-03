
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
#include "constraints.h"

using namespace std;

void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
   char parFile[MAX_PSR][MAX_FILELEN];
   char timFile[MAX_PSR][MAX_FILELEN];
   char newparname[MAX_FILELEN];
   double dm_dt=100;
   double cm_dt=100;
   int mode = 1;
   int p,i;


   const char *CVS_verNum = "$Revision: 1.1 $";

   if (displayCVSversion == 1) CVSdisplayVersion("autoDM.C","plugin",CVS_verNum);

   printf("Graphical Interface: autoDM\n");
   printf("Author:              M. Keith\n");
   printf("CVS Version:         $Revision: 1.1 $\n");

   /* Obtain all parameters from the command line */
   for (i=2;i<argc;i++)
   {
	  if (strcmp(argv[i],"-f")==0) {
		 strcpy(parFile[0],argv[++i]); 
		 strcpy(timFile[0],argv[++i]);
	  } else if(strcmp(argv[i],"-dt")==0){
		 dm_dt=atof(argv[++i]);
		 cm_dt=dm_dt;
	  } else if(strcmp(argv[i],"-DMdt")==0){
		 dm_dt=atof(argv[++i]);
	  } else if(strcmp(argv[i],"-CMdt")==0){
		 cm_dt=atof(argv[++i]);
	  } else if(strcmp(argv[i],"-mode")==0){
		 mode=atoi(argv[++i]);
	  }
   }

   readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
   readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
   preProcess(psr,*npsr,argc,argv);

   formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
   formResiduals(psr,*npsr,1);    /* Form the residuals                 */



   logmsg("TEST %d",*npsr);
   for (p=0; p < *npsr; p++){
	  longdouble c=0;
	  for(i=0; i < psr[p].nobs; i++){
		 if (psr[p].obsn[i].sat < c){
			logerr("Data need to be sorted");
			exit(1);

		 }
		 c=psr[p].obsn[i].sat;
	  }

	  sprintf(newparname,"%s_dm.par",psr[p].name);
	  double start=(double)(psr[p].obsn[0].sat);
	  double end=(double)(psr[p].obsn[psr[p].nobs-1].sat);
	  double len = end-start;
	  if (mode==1){
		 double centre=start+len/2.0;
		 double frac = len/dm_dt - (int)(len/dm_dt);
		 if (frac > 0.5){
			len+=1-frac;
		 } else{
			len-=frac;
		 }
		 start = centre-len/2.0;
		 end=centre+len/2.0;
	  } else if (mode==2){
		 double ndm=floor(len/dm_dt+0.5);
		 double ncm = floor(len/cm_dt+0.5);
		 dm_dt=len/ndm;
		 cm_dt=len/ncm;
	  }
	  logmsg("%s dmdt=%lf cmdt=%lf start=%lf end=%lf",psr[p].name,dm_dt,cm_dt,start,end);
	  psr[p].param[param_dmmodel].fitFlag[0] = 1;
	  psr[p].param[param_dmmodel].paramSet[0] = 1;
	  psr[p].param[param_dmmodel].val[0] = 2;
	  psr[p].param[param_dmmodel].nLinkTo=0;
	  psr[p].param[param_dm].nLinkFrom=0;
	  psr[p].param[param_dmmodel].linkTo[(psr[p].param[param_dmmodel].nLinkTo)++] = param_dm;
	  psr[p].param[param_dm].linkFrom[(psr[p].param[param_dm].nLinkFrom)++]=param_dmmodel;

	  autosetDMCM(psr+p,dm_dt,cm_dt,start,end,false);
	  psr[p].nconstraints = 0;
	  psr[p].fitMode=1;

	  psr[p].constraints[psr[p].nconstraints++] = constraint_dmmodel_mean;
	  psr[p].constraints[psr[p].nconstraints++] = constraint_dmmodel_cw_0;
	  psr[p].constraints[psr[p].nconstraints++] = constraint_dmmodel_cw_1;
	  psr[p].constraints[psr[p].nconstraints++] = constraint_dmmodel_cw_2;


	  textOutput(psr+p,1,0.0,0,0,1,newparname);  /* Display the output */
   }
   return 0;
}

char * plugVersionCheck = TEMPO2_h_VER;
