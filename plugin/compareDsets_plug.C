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
		     char timFile[MAX_PSR_VAL][MAX_FILELEN],double maxDiff,char *compare,char *compare2);
int checkSecondComparison(pulsar *psr,int i,int j,char *compare2);
int findOverlap(pulsar *psr,int *npsr,int *overlap1,int *overlap2,double maxDiff,char *compare,char *compare2);
int idPoint(pulsar *psr,int np,float *x_1,float *y_1,int *id_1,int count_1,float *x_2,float *y_2,int *id_2,int count_2,float mouseX,float mouseY,char parFile[MAX_PSR_VAL][MAX_FILELEN], char timFile[MAX_PSR_VAL][MAX_FILELEN],int view);
int idPoint2(pulsar *psr,int np,float *x,float *y,int *id1,int count,float mouseX,float mouseY,int *overlap1,int *overlap2,int view);

void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  char compare[128]="time";
  char compare2[128]="null";
  int i;
  int autoRun=0;
  double globalParameter;
  double maxDiff = 120; // Number of seconds corresponding to overlap
  const char *CVS_verNum = "$Revision: 1.5 $";

  if (displayCVSversion == 1) CVSdisplayVersion("compareDsets.C","plugin",CVS_verNum);

  *npsr = 0;

  printf("Graphical Interface: compareDsets\n");
  printf("Author:              G. Hobbs, X. You\n");
  printf("CVS Version:         $Revision: 1.5 $\n");
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
      else if (strcmp(argv[i],"-auto")==0)
	autoRun=1;
      else if (strcmp(argv[i],"-compare")==0)
	strcpy(compare,argv[++i]);
      else if (strcmp(argv[i],"-compare2")==0)
	strcpy(compare2,argv[++i]);
      else if (strcmp(argv[i],"-delta")==0)
	sscanf(argv[++i],"%lf",&maxDiff);
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

  if (autoRun==1)
    {
      int nOverlap;
      int overlap1[MAX_OBSN];
      int overlap2[MAX_OBSN];
      
      nOverlap = findOverlap(psr,npsr,overlap1,overlap2,maxDiff,compare,compare2);
      printf("Have automatically found %d overlapping points\n",nOverlap);
      writeTim((char *)"compare_dset1.tim",&psr[0],(char *)"tempo2");
      writeTim((char *)"compare_dset2.tim",&psr[1],(char *)"tempo2");
      printf("Written out compare_dset1.tim and compare_dset2.tim\n");
    }
  else
    compareDatasets(psr,npsr,parFile,timFile,maxDiff,compare,compare2);

  return 0;
}

void compareDatasets(pulsar *psr,int *npsr,char parFile[MAX_PSR_VAL][MAX_FILELEN],
		     char timFile[MAX_PSR_VAL][MAX_FILELEN],double maxDiff,char *compare,char *compare2)
{
  float x1[MAX_OBSN],y1[MAX_OBSN],e1[MAX_OBSN],y1e_u[MAX_OBSN],y1e_l[MAX_OBSN];
  float x2[MAX_OBSN],y2[MAX_OBSN],e2[MAX_OBSN],y2e_u[MAX_OBSN],y2e_l[MAX_OBSN];
  int id1[MAX_OBSN],id2[MAX_OBSN];
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
  char grDev[128];
  char xstr[128];
  int changePlot=0;
  char flagStr[128]="-length";
  int found1,found2;

  strcpy(grDev,"/xs");

  do {

    if (plot==1)
      {
	n1=0,n2=0;
	for (i=0;i<psr[0].nobs;i++)
	  {
	    if (psr[0].obsn[i].deleted==0)
	      {
		id1[n1] = i;
		x1[n1] = psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0];
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
		id2[n2] = i;
		x2[n2] = psr[1].obsn[i].sat-psr[0].param[param_pepoch].val[0];
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
    if (plot==2)
      {
	n1=0;
	for (i=0;i<nOverlap;i++)
	  {
	    id1[n1] = overlap1[i];
	    x1[n1] = psr[0].obsn[overlap1[i]].sat-psr[0].param[param_pepoch].val[0];
	    y1[n1] = (psr[0].obsn[overlap1[i]].residual - psr[1].obsn[overlap2[i]].residual)*1e6;
	    e1[n1] = (sqrt(pow(psr[0].obsn[overlap1[i]].toaErr*1e-6,2)+pow(psr[1].obsn[overlap2[i]].toaErr*1e-6,2)))*1e6;
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

    if (plot==3 || plot==4 || plot==5)
      {
	n1=0;
	for (i=0;i<nOverlap;i++)
	  {
	    id1[n1] = overlap1[i];
	    if (plot==3)
	      {
		x1[n1] = psr[0].obsn[overlap1[i]].toaErr;
		y1[n1] = psr[1].obsn[overlap2[i]].toaErr;
	      }
	    else if (plot==4)
	      {
		x1[n1] = psr[0].obsn[overlap1[i]].freq;
		y1[n1] = psr[1].obsn[overlap2[i]].freq;
	      }
	    else if (plot==5)
	      {
		int kk;
		found1=found2=0;
 		for (kk=0;kk<psr[0].obsn[overlap1[i]].nFlags;kk++)
		  {
		    if (strcmp(psr[0].obsn[overlap1[i]].flagID[kk],flagStr)==0)
		      {
			sscanf(psr[0].obsn[overlap1[i]].flagVal[kk],"%f",&x1[n1]);
			found1=1;
			break;
		      }
		  }

 		for (kk=0;kk<psr[1].obsn[overlap2[i]].nFlags;kk++)
		  {
		    if (strcmp(psr[1].obsn[overlap2[i]].flagID[kk],flagStr)==0)
		      {
			sscanf(psr[1].obsn[overlap2[i]].flagVal[kk],"%f",&y1[n1]);
			found2=1;
			break;
		      }
		  }
	      }
	    //	    printf("found = %d %d %g %g\n",found1,found2,y1[n1],y2[n1]);
	    if (plot!=5 || (found1==1 && found2 == 1))	    
	      {
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
	printf("n1 = %d\n",n1);
	if (zoom==0)
	  {
	    if (minx > miny) minx = miny;
	    else miny = minx;

	    if (maxx < maxy) maxx = maxy;
	    else maxy = maxx;
	  }
      }

    cpgbeg(0,grDev,1,1);
    cpgsch(1.4);
    cpgscf(1);
    cpgslw(2);
    cpgenv(minx,maxx,miny,maxy,0,1);
    sprintf(xstr,"MJD - %.1f",(double)psr[0].param[param_pepoch].val[0]);
    if (plot==2)
      cpglab(xstr,"Difference (\\gms)","");
    if (plot==1)
      cpglab(xstr,"Residual (s)","");
    if (plot==3)
      {
	char str1[128],str2[128];
	sprintf(str1,"%s, %s, error bar (us)",parFile[0],timFile[0]);
	sprintf(str2,"%s, %s, error bar (us)",parFile[1],timFile[1]);
	cpglab(str1,str2,"");
      }
    if (plot==4)
      {
	char str1[128],str2[128];
	sprintf(str1,"%s, %s, frequency (MHz)",parFile[0],timFile[0]);
	sprintf(str2,"%s, %s, frequency (MHz)",parFile[1],timFile[1]);
	cpglab(str1,str2,"");
      }
    if (plot==5)
      {
	char str1[128],str2[128];
	sprintf(str1,"%s, %s, %s",flagStr,parFile[0],timFile[0]);
	sprintf(str2,"%s, %s, %s",flagStr,parFile[1],timFile[1]);
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
    if (plot==3 || plot==4 || plot==5) 
      {
	float fx[2],fy[2];
	fx[0] = minx; fx[1] = maxx;
	fy[0] = minx; fy[1] = maxx;
	cpgsci(3); cpgline(2,fx,fy); cpgsci(1);
      }

    cpgsci(1);
    if (changePlot==0)
      {
	cpgcurs(&mx,&my,&key);
	if (key=='o') // find overlap points
	  {	
	    nOverlap = findOverlap(psr,npsr,overlap1,overlap2,maxDiff,compare,compare2);
	    printf("Number of overlapping points = %d\n",nOverlap);
	  }
	else if (key=='t') // Output overlap points to file
	  {
	    FILE *fout;
	    fout = fopen("overlap.dat","w");
	    for (i=0;i<nOverlap;i++)
	      {
		fprintf(fout,"%s %s %.5Lf %g %g\n",psr[0].obsn[overlap1[i]].fname,psr[1].obsn[overlap2[i]].fname,psr[0].obsn[overlap1[i]].sat,
			(double)(psr[0].obsn[overlap1[i]].residual - psr[1].obsn[overlap2[i]].residual), sqrt(pow(psr[0].obsn[overlap1[i]].toaErr*1e-6,2)+pow(psr[1].obsn[overlap2[i]].toaErr*1e-6,2)));
	      }
	    fclose(fout);
	    
	  }
	else if (key=='A')
	  {
	    if (plot==1)
	      idPoint(psr,0,x1,y1,id1,n1,x2,y2,id2,n2,mx,my,parFile,timFile,0);
	    else if (plot==2 || plot==3 || plot==4 || plot==5)
	      idPoint2(psr,0,x1,y1,id1,n1,mx,my,overlap1,overlap2,0);
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
	else if (key=='4') // Plot difference between frequencies
	  plot=4;
	else if (key=='5') // Plot difference based on flag
	  plot=5;
	else if (key=='g') 
	  {
	    changePlot=1;
	    printf("Enter graphical device: ");
	    scanf("%s",grDev);
	  }
	else if (key=='D') // View point
	  {
	    if (plot==1)
	      idPoint(psr,0,x1,y1,id1,n1,x2,y2,id2,n2,mx,my,parFile,timFile,1);
	    else if (plot==2 || plot==3 || plot==4 || plot==5)
	      idPoint2(psr,0,x1,y1,id1,n1,mx,my,overlap1,overlap2,1);
	  }
	else if (key=='c') // Command point
	  {
	    if (plot==1)
	      idPoint(psr,0,x1,y1,id1,n1,x2,y2,id2,n2,mx,my,parFile,timFile,2);
	    else if (plot==2 || plot==3 || plot==4 || plot==5)
	      idPoint2(psr,0,x1,y1,id1,n1,mx,my,overlap1,overlap2,2);
	  }
	else  if (key!='q')
	  {
	    printf("Unknown key stroke %c\n",key);
	  }
      } 
    else {
      strcpy(grDev,"/xs");
      changePlot=0;
    }
  } while (key!='q');
  cpgend();
}

int findOverlap(pulsar *psr,int *npsr,int *overlap1,int *overlap2,double maxDiff,char *compare,char *compare2)
{
  int i,j,found;
  long double tdiff;
  int nOverlap=0;
  int good;

  // Check first comparison
  for (i=0;i<psr[0].nobs;i++)
    {
      if (psr[0].obsn[i].deleted==0)
	{
	  found=0;
	  for (j=0;j<psr[1].nobs;j++)
	    {
	      good = checkSecondComparison(psr,i,j,compare2);
	      if (good != -1)
		{
		  if (strcmp(compare,"time")==0){
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
		  else if (strcmp(compare,"file2")==0){
		    if (strstr(psr[0].obsn[i].fname,psr[1].obsn[j].fname)!=NULL)
		      {
			overlap1[nOverlap] = i;
			overlap2[nOverlap] = j;
			nOverlap++;
			found=1;
			break;
		      }
		  }
		  else if (strcmp(compare,"file1")==0){
		    if (strstr(psr[1].obsn[j].fname,psr[0].obsn[i].fname)!=NULL)
		      {
			overlap1[nOverlap] = i;
			overlap2[nOverlap] = j;
			nOverlap++;
			found=1;
			break;
		      }
		  }
		  else {
		    printf("Unknown comparison\n");
		  }
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

int checkSecondComparison(pulsar *psr,int i,int j,char *compare2)
{
  int good=-1;
  int k;
  double v1=0,v2=0;

  if (strcmp(compare2,"length")==0)
    {
      for (k=0;k<psr[0].obsn[i].nFlags;k++)
	{
	  if (strcmp(psr[0].obsn[i].flagID[k],"-length")==0)
	    sscanf(psr[0].obsn[i].flagVal[k],"%lf",&v1);
	}
      for (k=0;k<psr[1].obsn[j].nFlags;k++)
	{
	  if (strcmp(psr[1].obsn[j].flagID[k],"-length")==0)
	    sscanf(psr[1].obsn[j].flagVal[k],"%lf",&v2);
	}
      printf("Checking %g %g %g\n",v2,v1,v2/v1);
      if (fabs(v2/v1) > 0.8 && fabs(v2/v1) < 1.2)
	return 1;
      else
	return -1;
    }
  else if (strcmp(compare2,"null")==0)
    {
    }
  else {
    printf("Unknown second check\n");
    exit(1);
  }
  
  
  return 0;
}

int idPoint(pulsar *psr,int np,float *x_1,float *y_1,int *id_1,int count_1,float *x_2,float *y_2,int *id_2,int count_2,float mouseX,float mouseY,char parFile[MAX_PSR_VAL][MAX_FILELEN], char timFile[MAX_PSR_VAL][MAX_FILELEN],int view)
{
  int i,iclosest,l,d;
  float closest,x1,x2,x3,x4,y1,y2,y3,y4,xscale,yscale,xpos,ypos;

   cpgqvp(3,&x1,&x2,&y1,&y2);
   cpgqwin(&x3,&x4,&y3,&y4);
   xscale = (x2-x1)/(x4-x3);
   yscale = (y2-y1)/(y4-y3);
   mouseX = (mouseX-x3)*xscale;
   mouseY = (mouseY-y3)*yscale;
   iclosest=-1;
   for (i=0;i<count_1;i++)
   {
	  xpos = (x_1[i]-x3)*xscale;
	  ypos = (y_1[i]-y3)*yscale;
	  if (iclosest==-1)
	  {
		 iclosest=id_1[i];
		 closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
		 d=1;
	  }
	  else if (pow(xpos-mouseX,2)+pow(ypos-mouseY,2)<closest)
	  {
		 iclosest=id_1[i];
		 closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
		 d=0;
	  }
   }

   for (i=0;i<count_2;i++)
   {
	  xpos = (x_2[i]-x3)*xscale;
	  ypos = (y_2[i]-y3)*yscale;
	  if (pow(xpos-mouseX,2)+pow(ypos-mouseY,2)<closest)
	  {
		 iclosest=id_2[i];
		 closest = pow(xpos-mouseX,2)+pow(ypos-mouseY,2);
		 d=1;
	  }
   }

   printf("---------------------------------------------------\n");
   printf("Closest point is from %s %s\n",parFile[d],timFile[d]);
   printf("Filename = %s\n",psr[d].obsn[iclosest].fname);
   printf("ToA value = %f\n",(double)psr[d].obsn[iclosest].sat);
   printf("Error = %f (us)\n",(double)psr[d].obsn[iclosest].toaErr);
   printf("Frequency = %f (MHz)\n",(double)psr[d].obsn[iclosest].freq);
   printf("Flags = ");
   for (i=0;i<psr[d].obsn[iclosest].nFlags;i++)
	  printf(" \"%s\" %s ",psr[d].obsn[iclosest].flagID[i],psr[d].obsn[iclosest].flagVal[i]);
   printf("\n");

   printf("---------------------------------------------------\n");
   if (view==1)
     {
       char str[1024];
       sprintf(str,"pav -CDFTp %s -g10/xs\n",psr[d].obsn[iclosest].fname);
       system(str);
     }
   else if (view==2)
     {
       char str[1024];
       char cmd[1024];
       printf("Enter command (e.g. vap -c nchan) ");
       //       scanf("%s",cmd);
       fgets(cmd,1024,stdin);
       cmd[strlen(cmd)-1]='\0';
       sprintf(str,"%s %s",cmd,psr[d].obsn[iclosest].fname);
       printf("Running >%s<\n",str);
       system(str);
     }
   return iclosest;
}


int idPoint2(pulsar *psr,int np,float *x,float *y,int *id1,int count,float mouseX,float mouseY,int *overlap1,int *overlap2,int view)
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
   printf("ToA value 1 = %f\n",(double)psr[0].obsn[overlap1[iclosest]].sat);
   printf("ToA value 2 = %f\n",(double)psr[1].obsn[overlap2[iclosest]].sat);

   printf("Error 1 = %f (us)\n",(double)psr[0].obsn[overlap1[iclosest]].toaErr);
   printf("Error 2 = %f (us)\n",(double)psr[1].obsn[overlap2[iclosest]].toaErr);

   printf("Frequency 1 = %f (MHz)\n",(double)psr[0].obsn[overlap1[iclosest]].freq);
   printf("Frequency 2 = %f (MHz)\n",(double)psr[1].obsn[overlap2[iclosest]].freq);

   printf("Flags1 = ");
   for (i=0;i<psr[0].obsn[overlap1[iclosest]].nFlags;i++)
	  printf(" \"%s\" %s ",psr[0].obsn[overlap1[iclosest]].flagID[i],psr[0].obsn[overlap1[iclosest]].flagVal[i]);
   printf("\n");

   printf("Flags2 = ");
   for (i=0;i<psr[1].obsn[overlap2[iclosest]].nFlags;i++)
	  printf(" \"%s\" %s ",psr[1].obsn[overlap2[iclosest]].flagID[i],psr[1].obsn[overlap2[iclosest]].flagVal[i]);
   printf("\n");
   printf("---------------------------------------------------\n");
   if (view==1)
     {
       char str[1024];
       sprintf(str,"pav -CDFTp -N1,2 %s %s -g10/xs\n",psr[0].obsn[overlap1[iclosest]].fname,psr[1].obsn[overlap2[iclosest]].fname);
       system(str);
     }
   else if (view==2)
     {
       char str[1024];
       char cmd[1024];
       printf("Enter command (e.g. vap -c nchan) ");
       //       scanf("%s",cmd);
       fgets(cmd,1024,stdin);
       cmd[strlen(cmd)-1]='\0';
       sprintf(str,"%s %s %s",cmd,psr[0].obsn[overlap1[iclosest]].fname,psr[1].obsn[overlap2[iclosest]].fname);
       printf("Running >%s<\n",str);
       system(str);
     }

   return iclosest;
}


char * plugVersionCheck = TEMPO2_h_VER;

