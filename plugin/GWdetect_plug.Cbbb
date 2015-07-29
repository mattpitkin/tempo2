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
#include "GWsim.h"
#include "T2toolkit.h"
#include <cpgplot.h>
#include "tempo2.h"
#include "TKspectrum.h"

using namespace std;

void searchGridPos(double dlat,double dlong,int gridPos,pulsar *psr,int npsr,char *addname);

void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,j,p,k;
  double globalParameter;
  int nlong = 60, nlat = 60;
  double dlat,dlong,dlongStep;
  int gridPos=0;
  char fname[128];
  FILE *fout;
  int posFileN=0;
  char posFile[128];
  char app[128]="";
  int randPos=0;
  FILE *fin;
  long idum = TKsetSeed();

  *npsr = 0; 

  printf("Graphical Interface: GWdetect\n");
  printf("Author:              G. Hobbs, R. Shannon\n");
  printf("Version:             2.0\n");
  printf(" --- type 'h' for help information\n");


  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[*npsr],argv[i+1]); 
	  strcpy(timFile[*npsr],argv[i+2]);
	  (*npsr)++;
	}
      else if (strcmp(argv[i],"-posfile")==0)
	{
	  strcpy(posFile,argv[++i]);
	  posFileN++;
	}
      else if (strcmp(argv[i],"-app")==0)
	strcpy(app,argv[++i]);
      else if (strcmp(argv[i],"-nlat")==0)
	sscanf(argv[++i],"%d",&nlat);
      else if (strcmp(argv[i],"-nlong")==0)
	sscanf(argv[++i],"%d",&nlong);
      else if (strcmp(argv[i],"-randpos")==0)
	randPos=1;
    }

  gridPos=0;

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

  // If requested randomise the pulsar positions for use in the GW angles
  if (randPos==1)
    {
      for (p=0;p<*npsr;p++)
	{
	  psr[p].param[param_raj].val[1] = TKranDev(&idum)*2*M_PI;   
	  psr[p].param[param_decj].val[1] = acos((TKranDev(&idum)-0.5)*2);
	  psr[p].param[param_raj].fitFlag[1] = 0;
	  psr[p].param[param_decj].fitFlag[1] = 0;
	  psr[p].param[param_raj].paramSet[1] = 1;
	  psr[p].param[param_decj].paramSet[1] = 1;
	}
    }
  
  if (posFileN==1)
    {
      gridPos=0;
      fin = fopen(posFile,"r");
      while (!feof(fin))
	{
	  if (fscanf(fin,"%lf %lf",&dlong,&dlat)==2)
	    {
	      searchGridPos(dlat,dlong,gridPos,psr,*npsr,app);
	      gridPos++;
	    }
	}
    }
  else
    {
      for (dlat = -M_PI/2.0;dlat < M_PI/2.0;dlat+=(M_PI/nlat))    
	{      
	  if (cos(dlat)!=0.0)
	    dlongStep = 2.0*M_PI/nlong/cos(dlat);
	  else
	    dlongStep = 2.0*M_PI;
	  
	  for (dlong=0;dlong<2.0*M_PI;dlong+=dlongStep)
	    {
	      searchGridPos(dlat,dlong,gridPos,psr,*npsr,app);
	      gridPos++;
	    }
	}
    }

  return 0;
}

void searchGridPos(double dlat,double dlong,int gridPos,pulsar *psr,int npsr,char *addname)
{
  FILE *fout;
  int i,j,k,p;
  char fname[128];
  double omega;
  int useQuad;

  if (psr[0].param[param_quad_om].paramSet[0]==1)
    useQuad=1;
  else
    useQuad=0;

  // Change the position of the GW source
  printf("ABOUT TO TRY: %.15f %.15f\n",dlong,dlat);
  if (gridPos!=0) // Reset the parameters to pre-fit values
    {
      for (p=0;p<npsr;p++)
	{
	  for (i=0;i<MAX_PARAMS;i++)
	    {
	      for (k=0;k<psr[p].param[i].aSize;k++)
		{
		  if (psr[p].param[i].fitFlag[k] == 1)
		    psr[p].param[i].val[k] = psr[p].param[i].prefit[k];
		}
	    }
	}
    }
  if (useQuad==1)
    {    
      for (p=0;p<npsr;p++)
	{
	  psr[p].quadRA  = dlong;
	  psr[p].quadDEC = dlat;
	  for (i=0;i<psr[p].nQuad;i++)
	    {
	      psr[p].quad_aplus_r[i] = 0;
	      psr[p].quad_aplus_i[i] = 0;
	      psr[p].quad_across_r[i] = 0;
	      psr[p].quad_across_i[i] = 0;
	    }
	}
    }  
  else
    {
      for (p=0;p<npsr;p++)
	{
	  psr[p].quad_ifunc_p_RA  = psr[p].quad_ifunc_c_RA = dlong;
	  psr[p].quad_ifunc_p_DEC = psr[p].quad_ifunc_c_DEC = dlat;
	  for (i=0;i<psr[p].quad_ifuncN_p;i++)
	    psr[p].quad_ifuncV_p[i] = 0;
	  for (i=0;i<psr[p].quad_ifuncN_c;i++)
	    psr[p].quad_ifuncV_c[i] = 0;
	  psr[p].quad_ifunc_geom_p = 0;
	  psr[p].quad_ifunc_geom_c = 0;
	}
    }
  if (gridPos==0 || useQuad==0) // if useQuad = 0 then need to recalculate the geometrical factors
    {
      formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,npsr,1);     /* Form the residuals                 */
    }
  // If fitting for pulsar parameters then must re-set the parameters each time
  doFit(psr,npsr,0);             /* Do the fitting                     */
  sprintf(fname,"gwdetect_grid%d.dat%s",gridPos,addname);
  if (!(fout = fopen(fname,"w")))
    {
      printf("Unable to open output file: %s\n",fname);
      exit(1);
    }
  if (useQuad==1)
    {
      fprintf(fout,"# QUAD\n");
      fprintf(fout,"# %.15g %.15g\n",dlong,dlat);
      fprintf(fout,"# %.15g %.15g %.15g %d\n",(double)psr[0].quadEpoch,(double)psr[0].obsn[0].sat,
	      (double)psr[0].obsn[psr[0].nobs-1].sat,psr[0].nQuad);
      fprintf(fout,"# %d\n",npsr);
      for (i=0;i<npsr;i++)
	fprintf(fout,"# %s %.15g %.15g\n",psr[i].name,(double)psr[i].param[param_raj].val[0],(double)psr[i].param[param_decj].val[0]);
      for (i=0;i<psr[0].nQuad;i++)
	{
	  omega = (double)(psr[0].param[param_quad_om].val[0]*(i+1));
	  fprintf(fout,"%g %g %g %g %g %g %g %g %g\n",omega,
		  psr[0].quad_aplus_r[i],psr[0].quad_aplus_i[i],
		  psr[0].quad_across_r[i],psr[0].quad_across_i[i],
		  psr[0].quad_aplus_r_e[i],psr[0].quad_aplus_i_e[i],
		  psr[0].quad_across_r_e[i],psr[0].quad_across_i_e[i]);
	}
    }
  else
    {
      fprintf(fout,"# QUAD_IFUNC\n");
      fprintf(fout,"# %.15g %.15g\n",dlong,dlat);
      fprintf(fout,"# %d\n",npsr);
      for (i=0;i<npsr;i++)
	fprintf(fout,"# %s %.15g %.15g\n",psr[i].name,(double)psr[i].param[param_raj].val[0],(double)psr[i].param[param_decj].val[0]);
      for (i=0;i<psr[0].quad_ifuncN_p;i++)
	fprintf(fout,"plus %g %g %g\n",psr[0].quad_ifuncT_p[i],psr[0].quad_ifuncV_p[i],psr[0].quad_ifuncE_p[i]);
      for (i=0;i<psr[0].quad_ifuncN_c;i++)
	fprintf(fout,"cross %g %g %g\n",psr[0].quad_ifuncT_c[i],psr[0].quad_ifuncV_c[i],psr[0].quad_ifuncE_c[i]);
    }
  
  fclose(fout); 
}
