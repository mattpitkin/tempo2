//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russel Edwards

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
#include <math.h>
#include "tempo2.h"

void getLabel(pulsar *psr,char *lab,int i);

extern "C" int tempoOutput(int argc,char *argv[],pulsar *psr,int npsr) 
{  
  int p,i,j;
  char lab1[128],lab2[128];
  int resultType=1;
  double tol = 0.9;
  char useParam[128];

  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-result")==0)
	sscanf(argv[++i],"%d",&resultType);
      else if (strcmp(argv[i],"-tol")==0)
	sscanf(argv[++i],"%lf",&tol);
      else if (strcmp(argv[i],"-param")==0)
	strcpy(useParam,argv[++i]);
    }

  printf("\n\n");
  for (p=0;p<npsr;p++)
    {
      printf("Pulsar: %s\n",psr[p].name);
      printf("\n");
      printf("Number of parameters in the fit: %d\n",psr[p].nParam);
      printf("Number of data points in the fit: %d\n",psr[p].nFit);
      printf("chisq of fit: %g\n",psr[p].fitChisq);
      printf("reduced chisq of fit: %g\n",psr[p].fitChisq/psr[p].fitNfree);

      // Print top headers
      if (resultType==1)
	{
	  printf("% 11s\t","Param");
	  for (j=0;j<psr[p].nParam;j++)
	    {
	      getLabel(&psr[p],lab2,j);
	      printf("% 11s\t",lab2);
	    }
	  printf("\n");
	  
	  for (i=0;i<psr[p].nParam;i++)
	    {
	      getLabel(&psr[p],lab1,i);
	      printf("% 11s\t",lab1);
	      for (j=0;j<=i;j++)
		{
		  printf("%+10.8f\t",psr[p].covar[i][j]/sqrt(psr[p].covar[i][i]*psr[p].covar[j][j]));
		}
	      printf("\n");
	    }
	  
	}

      if (resultType==2)
	{
	  for (i=0;i<psr[p].nParam;i++)
	    {
	      for (j=0;j<i;j++)
		{
		  if (psr[p].covar[i][j]/sqrt(psr[p].covar[i][i]*psr[p].covar[j][j]) > tol ||
		      psr[p].covar[i][j]/sqrt(psr[p].covar[i][i]*psr[p].covar[j][j]) < -tol)
		    {
		      getLabel(&psr[p],lab1,i);
		      getLabel(&psr[p],lab2,j);
		  printf("%s %s %+.8f\n",lab1,lab2,psr[p].covar[i][j]/sqrt(psr[p].covar[i][i]*psr[p].covar[j][j]));
		    }
		}
	    }
	}
      if (resultType==3)
	{
	  int iuse=-1;
	  for (i=0;i<psr[p].nParam;i++)
	    {
	      getLabel(&psr[p],lab1,i);
	      if (strcasecmp(lab1,useParam)==0)
		iuse = i;	      
	    }
	  if (iuse!=-1)
	    {
	      for (j=0;j<psr[p].nParam;j++)
		{
		  getLabel(&psr[p],lab2,j);
		  printf("%s-%s: %+.8f\n",useParam,lab2,psr[p].covar[iuse][j]/sqrt(psr[p].covar[iuse][iuse]*psr[p].covar[j][j]));
		}

	    }
	}
    }
  printf("\n");
}

void getLabel(pulsar *psr,char *lab,int i)
{
  if (psr->fitParamI[i]==-1) strcpy(lab,"JUMP");
  else if (psr->fitParamI[i]==-2) sprintf(lab,"_DM_%.1f",(float)(psr->dmoffsDM_mjd[psr->fitParamK[i]]));
  else if (psr->fitParamI[i]==-3) sprintf(lab,"_CM_%.1f",(float)(psr->dmoffsCM_mjd[psr->fitParamK[i]]));
  else if (i==0)sprintf(lab,"%s","YOFF");
  else sprintf(lab,"%s",psr->param[psr->fitParamI[i]].shortlabel[psr->fitParamK[i]]);
  
}

char * plugVersionCheck = TEMPO2_h_VER;
