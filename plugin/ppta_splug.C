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
#include <string.h>

using namespace std;


extern "C" int selectInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  int i;
  int j,k,p,nread;
  int found1,found2;
  int have_mcpsr2,have_ncpsr2;
  double tdiff,frange;
  char selectFname[1000]="";
  int simulID[MAX_OBSN];
  int nsimul;

  printf("Plugin to select points for the PPTA project\n");
  printf("Step 1) uses the 'select' file is present\n");
  printf("Step 2) identify simultaneous points and then choose which to keep\n");

  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-select")==0)
	sscanf(argv[++i],"%s",selectFname);
    }
  if (strlen(selectFname) > 0) // Use a select file
    {
      printf("Using select file >%s<\n",selectFname);
      useSelectFile(selectFname,psr,*npsr);
    }
  else
    printf("No select file present\n");

  // Search for simultaneous points
  tdiff = 600; // Seconds
  frange = 400; // MHz
  for (p=0;p<*npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)
	{
	  nsimul=0;
	  for (j=i+1;j<psr[p].nobs;j++)
	    {
	      if (fabs(psr[p].obsn[i].sat-psr[p].obsn[j].sat)*SECDAY<tdiff && 
		  fabs(psr[p].obsn[i].freq - psr[p].obsn[j].freq) < frange && psr[p].obsn[j].deleted==0 && psr[p].obsn[i].deleted==0)
		{
		  if (nsimul==0)
		    simulID[nsimul++] = i;
		  simulID[nsimul++] = j;
		}
	    }
	  if (nsimul>0)
	    {
	      printf("Found %d simultaneous points: ",nsimul);
	      have_mcpsr2=-1;
	      have_ncpsr2=-1;
	      for (j=0;j<nsimul;j++)
		{
		  for (k=0;k<psr[p].obsn[simulID[j]].nFlags;k++)
		    {
		      printf(" %s %s ",psr[p].obsn[simulID[j]].flagID[k],psr[p].obsn[simulID[j]].flagVal[k]);
		      if (strstr(psr[p].obsn[simulID[j]].flagVal[k],"CPSR2n")!=NULL)
			have_ncpsr2=j;
		      if (strstr(psr[p].obsn[simulID[j]].flagVal[k],"CPSR2m")!=NULL)
			have_mcpsr2=j;
		    }
		}	      
	      printf("\n");
	      if (have_ncpsr2!=-1 && have_mcpsr2!=-1)
		{
		  double refErr;
		  double smallestV;
		  int    smallestI;
		  printf("Do have CPSR2 n and m points\n");
		  // Obtain typical error bar size / (sqrt(2)) for the CPSR2 points
		  refErr = (psr[p].obsn[simulID[have_ncpsr2]].toaErr+psr[p].obsn[simulID[have_mcpsr2]].toaErr)/2.0/sqrt(2);
		  // Now identify the smallest toaErr
		  smallestI = -1;
		  smallestV = refErr;
		  for (j=0;j<nsimul;j++)
		    {
		      if (j!=have_ncpsr2 && j!=have_mcpsr2 && psr[p].obsn[simulID[j]].toaErr < smallestV)
			{
			  smallestV = psr[p].obsn[simulID[j]].toaErr;
			  smallestI = j;
			}
		    }
		  if (smallestI == -1) // Keep CPSR2 - delete the rest
		    {
		      for (j=0;j<nsimul;j++)
			{
			  if (j!=have_ncpsr2 && j!=have_mcpsr2)
			      psr[p].obsn[simulID[j]].deleted=1;
			}
		    }
		  else
		    {
		      for (j=0;j<nsimul;j++)
			{
			  if (j!=smallestI)
			      psr[p].obsn[simulID[j]].deleted=1;
			}
		    }
		}
	      else
		{
		  double smallestV;
		  int    smallestI;
		  printf("Don't have CPSR2 points: selecting smallest TOA uncertainty\n");
		  
		  smallestI = 0;
		  smallestV = psr[p].obsn[simulID[0]].toaErr;
		  for (j=1;j<nsimul;j++)
		    {
		      if (psr[p].obsn[simulID[j]].toaErr < smallestV)
			{
			  smallestI = j;
			  smallestV = psr[p].obsn[simulID[j]].toaErr;
			}
		    }
		  // Now delete the other points
		  for (j=0;j<nsimul;j++)
		    {
		      if (smallestI != j)
			psr[p].obsn[simulID[j]].deleted=1;
		    }
		}
	    }

	}
	     
    }
  
  return 0;
}

