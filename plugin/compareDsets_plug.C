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
#include <cpgplot.h>

using namespace std;

void compareDatasets(pulsar *psr,int *npsr,char parFile[MAX_PSR_VAL][MAX_FILELEN],
		     char timFile[MAX_PSR_VAL][MAX_FILELEN]);
int findOverlap(pulsar *psr,int *npsr,int *overlap1,int *overlap2);
int idPoint(pulsar *psr,int np,float *x,float *y,int *id,int count,float mouseX,float mouseY);
int idPoint2(pulsar *psr,int np,float *x,float *y,int *id1,int count,float mouseX,float mouseY,int *overlap1,int *overlap2);

void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i;
  double globalParameter;
  const char *CVS_verNum = "$Revision$";

  if (displayCVSversion == 1) CVSdisplayVersion("compareDsets.C","plugin",CVS_verNum);

  *npsr = 0;

  printf("Graphical Interface: compareDsets\n");
  printf("Author:              G. Hobbs, X. You\n");
  printf("CVS Version:         $Revision$\n");
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

  compareDatasets(psr,npsr,parFile,timFile);

  return 0;
}

void compareDatasets(pulsar *psr,int *npsr,char parFile[MAX_PSR_VAL][MAX_FILELEN],
		     char timFile[MAX_PSR_VAL][MAX_FILELEN])
{
  float x1[MAX_OBSN],y1[MAX_OBSN],e1[MAX_OBSN],y1e_u[MAX_OBSN],y1e_l[MAX_OBSN];
  float x2[MAX_OBSN],y2[MAX_OBSN],e2[MAX_OBSN],y2e_u[MAX_OBSN],y2e_l[MAX_OBSN];
  int id1[MAX_OBSN];
  int overlap1[MAX_OBSN];
  int overlap2[MAX_OBSN];
  int nOverlap=0;
  float miny,maxy,minx,maxx;
  float mx,my;
  char key;
  int n1=0,n2=0;
  int i;
  int plot=1;
  int zoom=0;

  do {

    if (plot==1)
      {
	n1=0,n2=0;
	for (i=0;i<psr[0].nobs;i++)
	  {
	    if (psr[0].obsn[i].deleted==0)
	      {
		id1[n1] = i;
		x1[n1] = psr[0].obsn[i].sat;
		y1[n1] = psr[0].obsn[i].residual;
		e1[n1] = psr[0].obsn[i].toaErr*1e-6;
		y1e_u[n1] = y1[n1]+e1[n1];
		y1e_l[n1] = y1[n1]-e1[n1];
		if (zoom==0)
		  {
		    if (n1==0)
		      {
			minx = maxx = x1[n1];
			miny = maxy = y1[n1];	  
		      }
		    else 
		      {
			if (minx > x1[n1]) minx = x1[n1];
			if (maxx < x1[n1]) maxx = x1[n1];
			if (miny > y1[n1]) miny = y1[n1];
			if (maxy < y1[n1]) maxy = y1[n1];
		      }
		  }
		n1++;
	      }
	  }
	
	for (i=0;i<psr[1].nobs;i++)
	  {
	    if (psr[1].obsn[i].deleted==0)
	      {
		x2[n2] = psr[1].obsn[i].sat;
		y2[n2] = psr[1].obsn[i].residual;
		e2[n2] = psr[1].obsn[i].toaErr*1e-6;
		y2e_u[n2] = y2[n2]+e2[n2];
		y2e_l[n2] = y2[n2]-e2[n2];

		if (zoom==0)
		  {
		    if (minx > x2[n2]) minx = x2[n2];
		    if (maxx < x2[n2]) maxx = x2[n2];
		    if (miny > y2[n2]) miny = y2[n2];
		    if (maxy < y2[n2]) maxy = y2[n2];
		  }
		n2++;
	      }
	  }
      }
    printf("n1 = %d, n2 = %d, nOverlap = %d\n",n1,n2,nOverlap);
    if (plot==2)
      {
	n1=0;
	for (i=0;i<nOverlap;i++)
	  {
	    id1[n1] = overlap1[i];
	    x1[n1] = psr[0].obsn[overlap1[i]].sat;
	    y1[n1] = psr[0].obsn[overlap1[i]].residual - psr[1].obsn[overlap2[i]].residual;
	    e1[n1] = sqrt(pow(psr[0].obsn[overlap1[i]].toaErr*1e-6,2)+pow(psr[1].obsn[overlap2[i]].toaErr*1e-6,2));
	    y1e_u[n1] = y1[n1]+e1[n1];
	    y1e_l[n1] = y1[n1]-e1[n1];

	    if (zoom==0)
	      {
		if (n1==0)
		  {
		    minx = maxx = x1[n1];
		    miny = maxy = y1[n1];	  
		  }
		else 
		  {
		    if (minx > x1[n1]) minx = x1[n1];
		    if (maxx < x1[n1]) maxx = x1[n1];
		    if (miny > y1[n1]) miny = y1[n1];
		    if (maxy < y1[n1]) maxy = y1[n1];
		  }
	      }
	    
	    n1++;
	  }
	printf("n1 = %d\n",n1);
      }

    if (plot==3)
      {
	n1=0;
	for (i=0;i<nOverlap;i++)
	  {
	    id1[n1] = overlap1[i];
	    x1[n1] = psr[0].obsn[overlap1[i]].toaErr*1e-6;
	    y1[n1] = psr[1].obsn[overlap2[i]].toaErr*1e-6;

	    if (zoom==0)
	      {
		if (n1==0)
		  {
		    minx = maxx = x1[n1];
		    miny = maxy = y1[n1];	  
		  }
		else 
		  {
		    if (minx > x1[n1]) minx = x1[n1];
		    if (maxx < x1[n1]) maxx = x1[n1];
		    if (miny > y1[n1]) miny = y1[n1];
		    if (maxy < y1[n1]) maxy = y1[n1];
		  }
	      }
	    
	    n1++;
	  }
	printf("n1 = %d\n",n1);
      }

    cpgbeg(0,"/xs",1,1);
    cpgsch(1.4);
    cpgscf(1);
    cpgslw(2);
    cpgenv(minx,maxx,miny,maxy,0,1);
    if (plot==2)
      cpglab("MJD","Difference","");
    if (plot==1)
      cpglab("MJD","Residual (s)","");
    if (plot==3)
      {
	char str1[128],str2[128];
	sprintf(str1,"%s, %s, error bar (us)",parFile[0],timFile[0]);
	sprintf(str2,"%s, %s, error bar (us)",parFile[1],timFile[1]);
	cpglab(str1,str2,"");
      }
    cpgsci(2); cpgpt(n1,x1,y1,16);
    if (plot==1 || plot==2)
      cpgerry(n1,x1,y1e_u,y1e_l,1);
    if (plot==1) {
      cpgsci(3); cpgpt(n2,x2,y2,4);
      cpgerry(n2,x2,y2e_u,y2e_l,1);

      printf("Red = %s %s\n",parFile[0],timFile[0]);
      printf("Green = %s %s\n",parFile[1],timFile[1]);
    }
    if (plot==3)
      {
	float fx[2],fy[2];
	fx[0] = minx; fx[1] = maxx;
	fy[0] = minx; fy[1] = maxx;
	cpgsci(3); cpgline(2,fx,fy); cpgsci(1);
      }

    cpgsci(1);
    cpgcurs(&mx,&my,&key);
    if (key=='o') // find overlap points
      {	
	nOverlap = findOverlap(psr,npsr,overlap1,overlap2);
	printf("Number of overlapping points = %d\n",nOverlap);
      }
    else if (key=='t') // Output overlap points to file
      {
	FILE *fout;
	fout = fopen("overlap.dat","w");
	for (i=0;i<nOverlap;i++)
	  {
	    fprintf(fout,"%.5Lf %g %g\n",psr[0].obsn[overlap1[i]].sat,
		    (double)(psr[0].obsn[overlap1[i]].residual - psr[1].obsn[overlap2[i]].residual), sqrt(pow(psr[0].obsn[overlap1[i]].toaErr*1e-6,2)+pow(psr[1].obsn[overlap2[i]].toaErr*1e-6,2)));
	  }
	fclose(fout);

      }
    else if (key=='A')
      {
	if (plot==1)
	  idPoint(psr,0,x1,y1,id1,n1,mx,my);
	else if (plot==2 || plot==3)
	  idPoint2(psr,0,x1,y1,id1,n1,mx,my,overlap1,overlap2);
      }
    else if (key=='z')
      {
	float mx2,my2;
	cpgband(2,0,mx,my,&mx2,&my2,&key);
	zoom=1;
	minx = TKretMin_f(mx,mx2);
	maxx = TKretMax_f(mx,mx2);
	miny = TKretMin_f(my,my2);
	maxy = TKretMax_f(my,my2);
      }
    else if (key=='w')
      {
	writeTim((char *)"compare_dset1.tim",&psr[0],(char *)"tempo2");
	writeTim((char *)"compare_dset2.tim",&psr[1],(char *)"tempo2");
      }
    else if (key=='u')
      zoom=0;
    else if (key=='1') 
      plot=1;
    else if (key=='2') // Plot difference between residuals
      plot=2;
    else if (key=='3') // Plot difference between error bar sizes
      plot=3;
  } while (key!='q');
  cpgend();
}

int findOverlap(pulsar *psr,int *npsr,int *overlap1,int *overlap2)
{
  int i,j,found;
  long double tdiff;
  double maxDiff = 120; // Number of seconds corresponding to overlap
  int nOverlap=0;

  for (i=0;i<psr[0].nobs;i++)
    {
      if (psr[0].obsn[i].deleted==0)
	{
	  found=0;
	  for (j=0;j<psr[1].nobs;j++)
	    {
	      tdiff = fabs(psr[0].obsn[i].sat - psr[1].obsn[j].sat)*86400.0;
	      if (tdiff < maxDiff)
		{
		  overlap1[nOverlap] = i;
		  overlap2[nOverlap] = j;
		  nOverlap++;
		  found=1;
		  break;
		}
	    }
	  if (found==0)
	    psr[0].obsn[i].deleted=1;
	}
    }

  for (i=0;i<psr[1].nobs;i++)
    {
      if (psr[1].obsn[i].deleted==0)
	{
	  found=0;
	  for (j=0;j<psr[0].nobs;j++)
	    {
	      tdiff = fabs(psr[1].obsn[i].sat - psr[0].obsn[j].sat)*86400.0;
	      if (tdiff < maxDiff)
		{
		  found=1;
		  break;
		}
	    }
	  if (found==0)
	    psr[1].obsn[i].deleted=1;
	}
    }  
  return nOverlap;
}

int idPoint(pulsar *psr,int np,float *x,float *y,int *id1,int count,float mouseX,float mouseY)
{
   int i,iclosest,l;
   float closest,x1,x2,x3,x4,y1,y2,y3,y4,xscale,yscale,xpos,ypos;

   cpgqvp(3,&x1,&x2,&y1,&y2);
   cpgqwin(&x3,&x4,&y3,&y4);
   xscale = (x2-x1)/(x4-x3);
   yscale = (y2-y1)/(y4-y3);
   mouseX = (mouseX-x3)*xscale;
   mouseY = (mouseY-y3)*yscale;
   iclosest=-1;
   for (i=0;i<count;i++)
   {
	  xpos = (x[i]-x3)*xscale;
	  ypos = (y[i]-y3)*yscale;
	  if (iclosest==-1)
	  {
		 iclosest=id1[i];
		 closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
	  }
	  else if (pow(xpos-mouseX,2)+pow(ypos-mouseY,2)<closest)
	  {
		 iclosest=id1[i];
		 closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
	  }
   }

   printf("---------------------------------------------------\n");
   printf("Closest point has TOA number %d (starting from 1)\n",iclosest+1);
   printf("ToA value = %f\n",(double)psr[0].obsn[iclosest].sat);
   return iclosest;
}


int idPoint2(pulsar *psr,int np,float *x,float *y,int *id1,int count,float mouseX,float mouseY,int *overlap1,int *overlap2)
{
   int i,iclosest,l;
   float closest,x1,x2,x3,x4,y1,y2,y3,y4,xscale,yscale,xpos,ypos;

   cpgqvp(3,&x1,&x2,&y1,&y2);
   cpgqwin(&x3,&x4,&y3,&y4);
   xscale = (x2-x1)/(x4-x3);
   yscale = (y2-y1)/(y4-y3);
   mouseX = (mouseX-x3)*xscale;
   mouseY = (mouseY-y3)*yscale;
   iclosest=-1;
   for (i=0;i<count;i++)
   {
	  xpos = (x[i]-x3)*xscale;
	  ypos = (y[i]-y3)*yscale;
	  if (iclosest==-1)
	  {
		 iclosest=i;
		 closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
	  }
	  else if (pow(xpos-mouseX,2)+pow(ypos-mouseY,2)<closest)
	  {
		 iclosest=i;
		 closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
	  }
   }

   printf("---------------------------------------------------\n");
   printf("Closest point represents: %s %s\n",psr[0].obsn[overlap1[iclosest]].fname,
	  psr[1].obsn[overlap2[iclosest]].fname);
   printf("ToA value = %f\n",(double)psr[0].obsn[iclosest].sat);
   return iclosest;
}


char * plugVersionCheck = TEMPO2_h_VER;

