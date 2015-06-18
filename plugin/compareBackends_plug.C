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
#include <cpgplot.h>
#include "T2toolkit.h"

using namespace std;


void runPlugin(pulsar *psr,int npsr,char *flagID1,char *flagID2,char *flagVal1,char *flagVal2,char *grDev,double maxTimeDiff);

void help() /* Display help */
{
  printf("-dt set threshold for time difference between ToAs (minutes)\n");
  printf("-grDev set graphical device\n");
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i;
  double globalParameter;
  char flagID1[128],flagID2[128];
  char flagVal1[128],flagVal2[128];
  char grDev[128]="/xs";
  double maxTimeDiff = 5; // Minutes

  strcpy(flagID1,"-f");
  strcpy(flagID2,"-f");
  strcpy(flagVal1,"MULTI_PDFB2");
  strcpy(flagVal2,"MULTI_CPSR2n");

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: compareBackends\n");
  printf("Author:              G. Hobbs\n");
  printf("Version:             v1.0\n");
  printf(" --- type 'h' for help information\n");

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
 	  strcpy(parFile[0],argv[++i]); 
	  strcpy(timFile[0],argv[++i]);
	}
      else if (strcasecmp(argv[i],"-grDev")==0)
	strcpy(grDev,argv[++i]);
      else if (strcasecmp(argv[i],"-dt")==0)
	sprintf(argv[++i],"%lf",&maxTimeDiff);
      else if (strcmp(argv[i],"-h")==0)
	help();
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
  printf("\n\n");
  printf("Please type in the flag ID and value for the two backends that should be compared\n");
  printf("e.g., -f MULTI_PDFB4\n");
  printf("\n\n");
  printf("Backend 1: enter flagID flagVal ");
  scanf("%s %s",flagID1,flagVal1);
  printf("Backend 2: enter flagID flagVal ");
  scanf("%s %s",flagID2,flagVal2);

  runPlugin(psr,*npsr,flagID1,flagID2,flagVal1,flagVal2,grDev,maxTimeDiff);


  return 0;
}

void runPlugin(pulsar *psr,int npsr,char *flagID1,char *flagID2,char *flagVal1,char *flagVal2,char *grDev,double maxTimeDiff)
{
  int i;
  int j,k;
  int found1=-1;
  int found2=-1;
  int ploterr=1;
  //  double maxTimeDiff = 5.0; // Minutes
  double dres,ebar;
  float xval[10][MAX_OBSN],yval[10][MAX_OBSN],e1[10][MAX_OBSN],e2[10][MAX_OBSN];
  int id1[10][MAX_OBSN],id2[10][MAX_OBSN];
  int n[10];
  char key;
  float mx,my;
  FILE *fout;
  int drawAxis=1;
  int recalc=1;
  float minx,maxx;
  float miny,maxy;
  int it=0;
  int xaxis=1;
  int zoom=0;
  float zoomX1,zoomX2,zoomY1,zoomY2;
  float mouseX,mouseY,mouseX2,mouseY2;

  cpgbeg(0,grDev,1,1);
  cpgask(0);
  printf("\n\nYou should now see a pgplot display that shows the difference between the two backends\n");
  printf("Also a file: output.dat has been written to disk giving:\n");
  printf(" - filenames, residual difference, error bar on ToA 1, error bar on ToA 2\n");
  printf("\n\n");
  printf("Use the following keys in the pgplot window\n\n");
  printf("q = quit\n");
  printf("z = zoom in a region using the mouse cursor\n");
  printf("u = unzoom\n");
  printf("s = provide statistics on the data points within the zoomed region\n");
  printf("r = define a new set of backends to compare\n");
  printf("o = produce a text file with the difference values\n\n"); 
  do {
    if (recalc==1)
      {
	fout = fopen("output.dat","w");
	n[it]=0;
	for (i=0;i<psr[0].nobs;i++)
	  {
	    //	    printf("Processing observation: %d out of %d\n",i+1,psr[0].nobs);
	    found1=found2=-1;
	    for (k=0;k<psr[0].obsn[i].nFlags;k++)
	      {
		if (strcmp(flagID1,psr[0].obsn[i].flagID[k])==0 &&
		    strcmp(flagVal1,psr[0].obsn[i].flagVal[k])==0)
		  {
		    found1=i;
		    break;
		  }
	      }
	    if (found1!=-1)
	      {
		for (j=0;j<psr[0].nobs;j++)
		  {
		    for (k=0;k<psr[0].obsn[j].nFlags;k++)
		      {
			if (strcmp(flagID2,psr[0].obsn[j].flagID[k])==0 &&
			    strcmp(flagVal2,psr[0].obsn[j].flagVal[k])==0)
			  {
			    found2=j;
			    if (fabs(psr[0].obsn[i].sat - psr[0].obsn[j].sat) < maxTimeDiff/60.0/24.0)
			      {
				dres = psr[0].obsn[i].residual - psr[0].obsn[j].residual;
				if (xaxis==1)
				  xval[it][n[it]] = (float)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
				yval[it][n[it]] = (float)dres;
				ebar = sqrt(pow(psr[0].obsn[i].toaErr*1.0e-6,2)+pow(psr[0].obsn[j].toaErr*1.0e-6,2));
				e1[it][n[it]] = yval[it][n[it]] - ebar;
				e2[it][n[it]] = yval[it][n[it]] + ebar;
				id1[it][n[it]]  = i;
				id2[it][n[it]]  = j;
				n[it]++;
				fprintf(fout,"Have match: %s %s %g %g %g %g\n",psr[0].obsn[i].fname,psr[0].obsn[j].fname,(double)psr[0].obsn[i].sat,dres,psr[0].obsn[i].toaErr*1.0e-6,psr[0].obsn[j].toaErr*1.0e-6);
			      }
			    break;
			  }
		      }
		  }
		
	      }
	    
	  }
	fclose(fout);
      }
    
    if (drawAxis==1)
      {
	for (i=0;i<=it;i++)
	  {
	    for (j=0;j<n[i];j++)
	      {
		if (i==0 && j==0)
		  {
		    minx = maxx = xval[i][j];
		    miny = maxy = yval[i][j];
		  }
		else
		  {
		    if (minx > xval[i][j]) minx = xval[i][j];
		    if (miny > yval[i][j]) miny = yval[i][j];
		    if (maxx < xval[i][j]) maxx = xval[i][j];
		    if (maxy < yval[i][j]) maxy = yval[i][j];
		    
		  }
	      }
	  }    
	cpgsci(1);
	if (zoom==1)
	  {
	    minx = zoomX1;
	    miny = zoomY1;
	    maxx = zoomX2;
	    maxy = zoomY2;
	  }
	cpgenv(minx,maxx,miny,maxy,0,1);
	cpglab("days since PEPOCH","Difference (sec)","");
	//	drawAxis=0;
      }


    for (j=0;j<=it;j++)
      {
	cpgsci(j+2);
	cpgpt(n[j],xval[j],yval[j],4);
	if (ploterr==1) cpgerry(n[j],xval[j],e1[j],e2[j],1);
    
      }
    cpgcurs(&mx,&my,&key);
    if (key=='r')
      {
	printf("Backend 1: enter flagID flagVal ");
	scanf("%s %s",flagID1,flagVal1);
	printf("Backend 2: enter flagID flagVal ");
	scanf("%s %s",flagID2,flagVal2);
	it++;
      }
    else if (key=='z')
      {
	cpgband(2,0,mx,my,&mouseX2,&mouseY2,&key);
	zoom=1;
	zoomX1 = TKretMin_f(mx,mouseX2);
	zoomX2 = TKretMax_f(mx,mouseX2);
	zoomY1 = TKretMin_f(my,mouseY2);
	zoomY2 = TKretMax_f(my,mouseY2);	
      }
    else if (key=='o') // Output
      {
	int k;
	FILE *fout;
	char fname[128];
	fout = fopen("result.dat","w");
	for (k=0;k<=it;k++)
	  {
	    for (j=0;j<=n[k];j++)
	      fprintf(fout,"%d %d %.4f %g\n",k,j,xval[k][j],yval[k][j]);
	  }
	fclose(fout);
      }
    else if (key=='s') // Statistics
      {
	double meanX=0.0,meanY=0.0;
	double y2=0.0;
	int nc=0,k;
	for (k=0;k<=it;k++)
	  {
	    meanX = meanY = y2 = 0.0;
	    nc  = 0;
	    for (j=0;j<=n[k];j++)
	      {
		if (xval[k][j] > minx && xval[k][j] < maxx && yval[k][j] > miny && yval[k][j] < maxy)
		  {
		    meanX += xval[k][j];
		    meanY += yval[k][j];
		    y2 += pow(yval[k][j],2);
		    nc++;
		  }
	      }
	    printf("[%d] Range processed: minx = %g, maxx = %g, miny = %g, maxy = %g\n",k+1,minx,maxx,miny,maxy); 
	    printf("[%d] MeanX = %g\n",k+1,meanX/(double)nc);
	    printf("[%d] MeanY = %g (s) = %g (us) = %g (ns)\n",k+1,meanY/(double)nc,meanY/(double)nc*1e6,meanY/(double)nc*1e9);
	    printf("[%d] Standard deviation in Y = %g (s) = %g (ns)\n",k+1,sqrt(y2/nc-pow(meanY/nc,2)),sqrt(y2/nc-pow(meanY/nc,2))*1e9);
	  }
	printf("-----------------------\n");
      }
    else if (key=='u')
      zoom=0;
  } while (key!='q');
  cpgend();

}
char * plugVersionCheck = TEMPO2_h_VER;
