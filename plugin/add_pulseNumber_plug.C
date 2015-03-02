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

/* This plugin allows the user to define the output format for listing observational based
 * values
 *
 * On the command line: use -s "....." option
 * In a file: use -file filename
 *
 * options {resPre} - prefit residuals
 *         {resPost} - postfit residuals
 * 
 * 
 * 
 *
 * Other commands
 *
 * {ERRMULT x} multiple all errors by x
 * {FORMAT c}  sets the format for displaying values and errors (printf-type format)
 * {TAB x}     move the cursor to position x
 * {NULL c}    set the null string to display if a parameter is not set
 *
 * All other characters are printed 
 * \n newline
 * \{ print a { character
 * \} print a } character
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"

extern "C" int tempoOutput(int argc,char *argv[],pulsar *psr,int npsr) 
{  
  int i;
  int flagid;

  for (i=0;i<psr[0].nobs;i++)
    {
      flagid=psr[0].obsn[i].nFlags;
      for(int k=0; k < flagid ; k++){
	if(strcmp(psr[0].obsn[i].flagID[k],"-pnadd")==0){
	  printf("Removing -pnadd flag\n");
	  for (int kk=k; kk < flagid-1; kk++){
	    strcpy(psr[0].obsn[i].flagID[kk],psr[0].obsn[i].flagID[kk+1]);
	    strcpy(psr[0].obsn[i].flagVal[kk],psr[0].obsn[i].flagVal[kk+1]);
	  }
	  flagid-=1;
	  k-=1;
	}
      }
      psr[0].obsn[i].nFlags = flagid;
      for(int k=0; k < flagid ; k++){
	if(strcmp(psr[0].obsn[i].flagID[k],"-pn")==0){
	  flagid=k;
	  break;
	}
      }

       strcpy(psr[0].obsn[i].flagID[flagid],"-pn");
       sprintf(psr[0].obsn[i].flagVal[flagid],"%lld",psr[0].obsn[i].pulseN-psr[0].obsn[0].pulseN);
       if (flagid==psr[0].obsn[i].nFlags)
	 psr[0].obsn[i].nFlags++;
       
    }
  writeTim("withpn.tim",psr,"tempo2");
  printf("In here\n");
}
