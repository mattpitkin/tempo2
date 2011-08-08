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

/* Plugin to produce a power spectrum of the input timing residuals */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "T2toolkit.h"
#include "TKspectrum.h"
#include <cpgplot.h>

using namespace std;
char dcmFile[MAX_FILELEN];
char covarFuncFile[MAX_FILELEN];

#define MAX_ID 50

void drawOption(int on,float x,float y,char *str);
void doPlugin(pulsar *psr,int npsr, char parFile[MAX_PSR_VAL][MAX_FILELEN], char timFile[MAX_PSR_VAL][MAX_FILELEN],int white,int filter);
void drawMenu(int specType,int xaxis,int logv,int specOut);
void checkMenu(float mx,float my,int *change,int *xaxis,int *logv,int *specType,int *specOut);
void identify(float mx,float my,float px[MAX_PSR_VAL][MAX_OBSN_VAL],float py[MAX_PSR_VAL][MAX_OBSN_VAL],
	      int *sn, int idV[MAX_ID],int idP[MAX_ID], int *iN,int npsr);
void model(pulsar *psr, char parFile[MAX_PSR_VAL][MAX_FILELEN], char timFile[MAX_PSR_VAL][MAX_FILELEN]);


void help() /* Display help */
{
  printf("-----------------------------------------\n");
  printf("left click:\t identify point\n");
  printf("h\t\t this help\n");
  printf("l\t\t text output of spectrum\n");
  printf("L\t\t text output of input residuals\n");
  printf("o\t\t set oversampling factor\n");
  printf("p\t\t toggle plotting points\n");
  printf("P\t\t toggle plotting line\n");
  printf("q\t\t quit\n");
  printf("u\t\t unzoom\n");
  printf("z\t\t zoom\n");
  printf("-----------------------------------------\n");
}


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,n=0,white=0,filter=0,res=0;
  char resFile[MAX_FILELEN];
  double px[MAX_OBSN],py[MAX_OBSN],pe[MAX_OBSN];

  *npsr = 0;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: spectrum\n");
  printf("Author:              G. Hobbs\n");
  printf("Version:             1.0\n");
  printf(" --- type 'h' for help information\n");

  strcpy(dcmFile,"NULL");
  strcpy(covarFuncFile,"NULL");

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[*npsr],argv[i+1]); 
	  strcpy(timFile[*npsr],argv[i+2]);
	  (*npsr)++;
	}
      else if (strcmp(argv[i],"-dcm")==0)
	strcpy(dcmFile,argv[++i]);
      else if (strcmp(argv[i],"-dcf")==0)
	strcpy(covarFuncFile,argv[++i]);
      else if (strcmp(argv[i],"-white")==0)
	sscanf(argv[++i],"%d",&white);
      else if (strcmp(argv[i],"-fil")==0)
	filter=1;
      else if (strcmp(argv[i],"-res")==0)
	{
	  res=1;
	  strcpy(resFile,argv[++i]);
	}
    }

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      if (i==0) 
	{
	  if (strcmp(dcmFile,"NULL")==0 && strcmp(covarFuncFile,"NULL")==0)
	    doFit(psr,*npsr,0);
	  else
	    doFitDCM(psr,dcmFile,covarFuncFile,*npsr,0);
	}
      else textOutput(psr,*npsr,0,0,0,0,"");  /* Display the output */
    }

  if (res==1) // Replace residuals with residual file
    {
      FILE *fin;
      psr[0].nobs=0;
      fin = fopen(resFile,"r");
      while (!feof(fin))
	{
	  if (fscanf(fin,"%Lf %Lf %lf",&psr[0].obsn[psr[0].nobs].sat,
		     &psr[0].obsn[psr[0].nobs].residual,
		     &psr[0].obsn[psr[0].nobs].toaErr)==3)
	    psr[0].nobs++;
	}
      fclose(fin);
    }
 

  doPlugin(psr,*npsr,parFile,timFile,white,filter);

  return 0;
}

void doPlugin(pulsar *psr,int npsr, char parFile[MAX_PSR_VAL][MAX_FILELEN], char timFile[MAX_PSR_VAL][MAX_FILELEN],int white,int filter)
{
  float px[MAX_PSR_VAL][MAX_OBSN_VAL],py[MAX_PSR_VAL][MAX_OBSN_VAL];
  double resx[MAX_OBSN],resy[MAX_OBSN],rese[MAX_OBSN];
  int n;
  float minx,maxx,miny,maxy;
  double specX[MAX_OBSN],specY[MAX_OBSN];
  float mouseX,mouseY;
  int p;
  char key;
  int specType=2;
  int specN[MAX_PSR];
  int xaxis=3;
  int endit=0;
  int change=1;
  int logv=0;
  int plotPoints=-1;
  int plotLine=1;
  int i;
  int idV[MAX_ID],idP[MAX_ID];
  int iN=0;
  int fontSize=1;
  char str[100];
  float chx,chy;
  int zoom=0;
  int ofac=1;
  int specOut=1; // 1 = PSD, 2 = Amplitude, 3 = Power
  char ylabel[100];
  char gr[100]="/xs";
  int changeDevice=0;
  double outY_re[MAX_OBSN],outY_im[MAX_OBSN];

  do {
    cpgbeg(0,gr,1,1);
    cpgsch(fontSize);
    cpgask(0);
    cpgeras();
    cpgslw(1);
    cpgsvp(0.0,1.0,0.8,1.0);
    cpgswin(0,1,0,1);
    drawMenu(specType,xaxis,logv,specOut);
    if (changeDevice==1)
      {
	strcpy(gr,"/xs");
	cpgslw(4);
	changeDevice=0;
      }
    cpgsvp(0.1,0.9,0.1,0.8);
    
    if (change==1)
      {
	for (p=0;p<npsr;p++)
	  {
	    n=0;
	    for (i=0;i<psr[p].nobs;i++)
	      {
		if (psr[p].obsn[i].deleted==0)
		  {
		    resx[n] = (double)(psr[p].obsn[i].sat-psr[p].param[param_pepoch].val[0]);
		    resy[n] = (double)psr[p].obsn[i].residual;
		    rese[n] = (double)psr[p].obsn[i].toaErr*1.0e-6;
		    //		    rese[n] = 1.0;
		    //	  py[n] = 3*sin(2.0*M_PI/2041.34587909003*10*px[n]);
		    n++;
		  }
	      }
	    if (filter==1)
	      {
		double smx[MAX_OBSN],smy[MAX_OBSN],sme[MAX_OBSN],top,bot,t;
		double wx[MAX_OBSN],wy[MAX_OBSN],we[MAX_OBSN];
		double ix[MAX_OBSN],iy[MAX_OBSN],ie[MAX_OBSN];
		int intp=21;
		int s=30;
		int n2=0,j;
		double yd[MAX_OBSN][4],h;
		
		// Smooth and interpolate
		//		    for (i=0;i<n;i++)
		//		      printf("orig %g %g\n",resx[i],resy[i]);
		//		    exit(1);
		/*		    do
				    {
				    sme[n2]=1.0;
				    smx[n2]=resx[0]+n2*intp;
				    top=0.0; bot=0.0;
				    for (i=0;i<n;i++)
				    {
				    t = 1.0/rese[i]/rese[i]*exp(-fabs(smx[n2]-resx[i])/s);
				    top+=t*resy[i];
				    bot+=t;
				    }
				    smy[n2]=top/bot;
				    n2++;
				    } while (smx[n2-1] < resx[n-1]);*/
		for (i=0;i<n;i++)
		  {
		    smx[i] = resx[i];
		    sme[i] = 1.0;
		    top=0.0; bot=0.0;
		    for (j=0;j<n;j++)
		      {
			t = 1.0/rese[j]/rese[j]*exp(-fabs(smx[i]-resx[j])/s);
			top+=t*resy[j];
			bot+=t;
		      }
		    smy[i]=top/bot;
		    wx[i] = resx[i];
		    we[i] = rese[i];
		    wy[i] = resy[i]-smy[i];
		  }
		for (i=0;i<n;i++) printf("smx %g %g %g %g\n",smx[i],smy[i],wy[i],we[i]);
		// Now interpolate onto a regular grid
		TKcmonot(n,smx,smy,yd);
		
		n2=0;
		do {
		  ie[n2]=1.0;
		  ix[n2]=smx[0]+n2*intp;		      
		  n2++;
		} while (ix[n2-1] < resx[n-1]);
		TKspline_interpolate(n,smx,smy,yd,ix,iy,n2);
		for (i=0;i<n2;i++)
		  printf("ix %g %g\n",ix[i],iy[i]);
		
		TKspectrum(ix,iy,ie,n2,0,0,0,0,white,specType,ofac,1,specOut,specX,specY,&specN[p],0,0,outY_re,outY_im);

	      }

	    if (specType==3)
	      {
		TKspectrum(resx,resy,rese,n,1,0,0,1,0,specType,ofac,1,specOut,specX,specY,&specN[p],0,0,outY_re,outY_im);
	      }
	    else
	      {
		printf("In here with %d\n",white);
		if (white==0)
		  TKspectrum(resx,resy,rese,n,0,0,0,0,0,specType,ofac,1,specOut,specX,specY,&specN[p],0,0,outY_re,outY_im);
		else
		  TKspectrum(resx,resy,rese,n,1,0,0,1,white,specType,ofac,1,specOut,specX,specY,&specN[p],0,0,outY_re,outY_im);
	      }
	    TKconvertFloat2(specX,specY,px[p],py[p],specN[p]);
	    if (xaxis==1) // Convert x-axis to s^-1
	      {
		for (i=0;i<specN[p];i++)
		  px[p][i]/=SECDAY;
	      }
	    else if (xaxis==3) // Convert x-axis to yr^-1
	      {
		for (i=0;i<specN[p];i++)
		  px[p][i]*=365.25;
	      }
	    else if (xaxis==4) // Convert x-axis to s
	      {
		for (i=0;i<specN[p];i++)
		  px[p][i]=SECDAY/px[p][i];
	      }
	    else if (xaxis==5) // Convert x-axis to d
	      {
		for (i=0;i<specN[p];i++)
		  px[p][i]=1.0/px[p][i];
	      }
	    else if (xaxis==6) // Convert x-axis to d
	      {
		for (i=0;i<specN[p];i++)
		  px[p][i]=1.0/px[p][i]/365.25;
	      }
	  
	    if (logv > 0 && change==1)
	      {
		for (i=0;i<specN[p];i++)
		  {
		    if (logv==1 || logv==3)
		      px[p][i] = log10(px[p][i]);
		    if (logv==2 || logv==3)
		      py[p][i] = log10(py[p][i]);
		  }
	      }
	  }
      }

    for (p=0;p<npsr;p++)
      {
	if (p==0)
	  {
	    minx = maxx = px[0][0];
	    miny = maxy = py[0][0];
	  }
	if (zoom==0 || change==1)
	  {
	    for (i=0;i<specN[p];i++)
	      {
		if (logv==0)
		  {
		    minx = 0.0;
		    miny = 0.0;
		  }
		else
		  {
		    if (minx > px[p][i]) minx = px[p][i];
		    if (miny > py[p][i]) miny = py[p][i];
		  }
		if (maxx < px[p][i]) maxx = px[p][i];
		if (maxy < py[p][i]) maxy = py[p][i];
	      }
	  }
      }
  
  
        cpgswin(minx,maxx,miny,maxy+0.1*(maxy-miny));
    //    cpgswin(minx,maxx,-35,-28);
    if (logv==0)
      cpgbox("BCNST1G",0.0,0,"BCNST1G",0.0,0);
    else if (logv==1)
      cpgbox("BCNST1LG",0.0,0,"BCNST1G",0.0,0);
    else if (logv==2)
      cpgbox("BCNST1G",0.0,0,"BCNST1LG",0.0,0);
    else if (logv==3)
      cpgbox("BCNST1LG",0.0,0,"BCNST1LG",0.0,0);
    
    if (specOut==1) strcpy(ylabel,"PSD");
    else if (specOut==2) strcpy(ylabel,"Amplitude");
    else if (specOut==3) strcpy(ylabel,"Power");
    
    if (xaxis==1) cpglab("Frequency (s\\u-1\\d)",ylabel,"");
    else if (xaxis==2) cpglab("Frequency (d\\u-1\\d)",ylabel,"");
    else if (xaxis==3) cpglab("Frequency (yr\\u-1\\d)",ylabel,"");
    else if (xaxis==4) cpglab("Sec",ylabel,"");
    else if (xaxis==5) cpglab("Day",ylabel,"");
    else if (xaxis==6) cpglab("Yr",ylabel,"");
    for (p=0;p<npsr;p++)
      {
	cpgsci(p+1);
	if (plotLine==1)
	  cpgline(specN[p],px[p],py[p]);
	if (plotPoints==1)
	  cpgpt(specN[p],px[p],py[p],9);
      }
    // Label identified points
    cpgsch(fontSize/1.8);
    cpgqcs(4,&chx,&chy);
    cpgsci(3);
    for (i=0;i<iN;i++)
      {
	sprintf(str,"id = %d",i+1); cpgtext(px[idP[i]][idV[i]]+chx,py[idP[i]][idV[i]]+chy,str);	
	sprintf(str,"index = %d",idV[i]+1); cpgtext(px[idP[i]][idV[i]]+chx,py[idP[i]][idV[i]],str);	
	sprintf(str,"(%.2g,%-.2g)",px[idP[i]][idV[i]],py[idP[i]][idV[i]]); 
	cpgtext(px[idP[i]][idV[i]]+chx,py[idP[i]][idV[i]]-chy,str);	
      }
    cpgsci(1);
    cpgsch(fontSize);

    change=0;

    cpgcurs(&mouseX,&mouseY,&key);
    checkMenu(mouseX,mouseY,&change,&xaxis,&logv,&specType,&specOut);
    if (key=='q') endit=1;
    else if (key=='A' && change==0) // Identify point 
      identify(mouseX,mouseY,px,py,specN,idV,idP,&iN,npsr);
    else if (key=='u') // Unzoom
      zoom=0;
    else if (key=='h') // Help
      help();
    else if (key=='g') // Change graphical device
      {
	printf("Graphical device: ");
	scanf("%s",gr);
	changeDevice=1;
      }
    else if (key=='o') // New oversampling factor
      {
	printf("Enter oversampling factor: ");
	scanf("%d",&ofac);
	change=1;
      }
    else if (key=='p') // Toggle plotting points
      plotPoints*=-1;
    else if (key=='P') // Toggle plotting line
      plotLine*=-1;
    else if (key=='m') // Model data
      model(psr,parFile,timFile);
    else if (key=='l') // List points
      {
	for (p=0;p<npsr;p++)
	  {
	    for (i=0;i<specN[p];i++)
	      printf("%d %d %g %g\n",p+1,i+1,px[p][i],py[p][i]);
	  }
      }
    else if (key=='z') // Zoom
      {
	float my,mx;
	zoom=1;
	cpgband(2,0,mouseX,mouseY,&mx,&my,&key);
	minx = TKretMin_f(mx,mouseX);
	maxx = TKretMax_f(mx,mouseX);
	miny = TKretMin_f(my,mouseY);
	maxy = TKretMax_f(my,mouseY);
      }
    else if (key=='L') // List input residuals
      {
	FILE *fout;
	fout = fopen("spectrum_res","w");
	for (i=0;i<n;i++)
	  {
	    printf("res %d %g %g\n",i+1,resx[i],resy[i]);
	    fprintf(fout,"%g %g\n",resx[i],resy[i]);
	  }
	fclose(fout);
      }
  } while (endit==0);
  printf("Goodbye\n");
  cpgend();
}

void identify(float mx,float my,float px[MAX_PSR_VAL][MAX_OBSN_VAL],float py[MAX_PSR_VAL][MAX_OBSN_VAL],
	      int *sn, int idV[MAX_ID],int idP[MAX_ID], int *iN,int npsr)
{
  int i,iclosest,p,pclosest;
  float closest,x1,x2,x3,x4,y1,y2,y3,y4,xscale,yscale,xpos,ypos;

  cpgqvp(3,&x1,&x2,&y1,&y2);
  cpgqwin(&x3,&x4,&y3,&y4);
  xscale = (x2-x1)/(x4-x3);
  yscale = (y2-y1)/(y4-y3);
  mx = (mx-x3)*xscale;
  my = (my-y3)*yscale;
  iclosest=-1;
  pclosest=-1;
  for (p=0;p<npsr;p++)
    {
      for (i=0;i<sn[p];i++)
	{
	  xpos = (px[p][i]-x3)*xscale;
	  ypos = (py[p][i]-y3)*yscale;
	  if (iclosest==-1)
	    {
	      iclosest=i;
	      pclosest=p;
	      closest = pow(xpos-mx,2)+pow(ypos-my,2);
	    }
	  else if (pow(xpos-mx,2)+pow(ypos-my,2)<closest)
	    {
	      iclosest=i;
	      pclosest=p;
	      closest = pow(xpos-mx,2)+pow(ypos-my,2);
	    }
	}
    }
  printf("Closest point has index %d (starting from 1)\n",iclosest+1);
  printf("This has coordinates (%g,%g)\n",px[p][iclosest],py[p][iclosest]);
  idV[*iN] = iclosest;
  idP[*iN] = pclosest;
  (*iN)++;
}

void drawMenu(int specType,int xaxis,int logv,int specOut)
{
  drawOption(specType,0.0,0.8,"DFT");
  drawOption(specType-1,0.1,0.8,"Lomb");
  drawOption(specType-2,0.2,0.8,"FFT");
  drawOption(specType-3,0.3,0.8,"WLS");
  drawOption(specType-4,0.4,0.8,"Sine");
  drawOption(xaxis,0.0,0.6,"s\\u-1\\d");
  drawOption(xaxis-1,0.1,0.6,"d\\u-1\\d");
  drawOption(xaxis-2,0.2,0.6,"yr\\u-1\\d");  
  drawOption(xaxis-3,0.3,0.6,"s");
  drawOption(xaxis-4,0.4,0.6,"d");
  drawOption(xaxis-5,0.5,0.6,"yr");
  drawOption(logv+1,0.0,0.4,"lin xy");
  drawOption(logv,0.1,0.4,"log x");
  drawOption(logv-1,0.2,0.4,"log y");
  drawOption(logv-2,0.3,0.4,"log xy");
  drawOption(specOut,0.0,0.2,"PSD");
  drawOption(specOut-1,0.1,0.2,"Amp");
  drawOption(specOut-2,0.2,0.2,"Pow");
}

void drawOption(int on,float x,float y,char *str)
{
  if (on==1) cpgsci(2);
  else cpgsci(1);

  cpgrect(x,x+0.015,y,y+0.08);
  cpgsci(1); cpgtext(x+0.017,y,str);
}

void checkMenu(float mx,float my,int *change,int *xaxis,int *logv,int *specType,int *specOut)
{
  float x1,x2,y1,y2,x3,x4,y3,y4,x7,y7,xscale,yscale,xscale2,x0;
  float x5,x6,y5,y6;
  int mouseX,mouseY;

  *change=0;
  cpgqvp(3,&x1,&x2,&y1,&y2);
  cpgqwin(&x3,&x4,&y3,&y4);
  cpgqvsz(3,&x5,&x6,&y5,&y6);
  xscale = (x2-x1)/(x4-x3);
  yscale = (y2-y1)/(y4-y3);

  x7 = (mx-x3)*xscale+x1; 
  y7 = (my-y3)*yscale+y1; 
  mouseX = (int)(x7/(x6-x5)*10.0);
  mouseY = (int)((1-(y7/(y6-y5)))/0.2*5.0); // 0.2 because the top menu bar has height of 0.2
  //  if (((mx-x3)*xscale+x1)/(x2+x1)*10 < mouseX+0.2 || mouseY==3)
  if (mouseY < 4)
    {
      *change=1;
      if (mouseY==0)
	{
	  if (mouseX==0) *specType=1;
	  else if (mouseX==1) *specType=2;
	  else if (mouseX==2) *specType=3;
	  else if (mouseX==3) *specType=4;
	  else if (mouseX==4) *specType=5;
	}
      if (mouseY==1)
	{
	  if (mouseX==0) *xaxis=1;
	  else if (mouseX==1) *xaxis=2;
	  else if (mouseX==2) *xaxis=3;
	  else if (mouseX==3) *xaxis=4;
	  else if (mouseX==4) *xaxis=5;
	  else if (mouseX==5) *xaxis=6;
	}
      else if (mouseY==2)
	{
	  if (mouseX==0) *logv=0;
	  else if (mouseX==1) *logv=1;
	  else if (mouseX==2) *logv=2;
	  else if (mouseX==3) *logv=3;
	}
      else if (mouseY==3)
	{
	  if (mouseX==0) *specOut=1;
	  else if (mouseX==1) *specOut=2;
	  else if (mouseX==2) *specOut=3;
	}
    }
}

void model(pulsar *psr, char parFile[MAX_PSR_VAL][MAX_FILELEN], char timFile[MAX_PSR_VAL][MAX_FILELEN])
{
  int i,j,k,p;
  int n;
  double resx[MAX_OBSN],resy[MAX_OBSN],rese[MAX_OBSN];
  double specX[MAX_OBSN],specY[MAX_OBSN];
  double storeRes[MAX_OBSN];
  float px[MAX_OBSN],py[MAX_OBSN],fy[MAX_OBSN],fx[MAX_OBSN];
  float minx,maxx,miny,maxy;
  double tau=120;
  double whiteLevel;
  double meanActualWhiteLevel,meanActual;
  double area1,area2;
  int ofac=1,nv;
  int specType=2;
  int specOut=3;
  int specN;
  long idum=TKsetSeed();
  int it,nit;
  double outY_re[MAX_OBSN],outY_im[MAX_OBSN];
  nit = 100;

  n=0;
  for (i=0;i<psr[0].nobs;i++)
    {
      if (psr[0].obsn[i].deleted==0)
	{
	  resx[n] = (double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
	  resy[n] = (double)psr[0].obsn[i].residual;
	  rese[n] = (double)psr[0].obsn[i].toaErr*1.0e-6;
	  storeRes[n] = (double)psr[0].obsn[i].residual;
	  n++;
	}
    }
  
  TKspectrum(resx,resy,rese,n,0,0,0,0,0,specType,ofac,1,specOut,specX,specY,&specN,0,0,outY_re,outY_im);
  TKconvertFloat2(specX,specY,px,py,specN);
  meanActual=0;
  nv=0;
  for (i=0;i<specN;i++)
    {
      if (specX[i] > 1.0/tau)
	{
	  nv++;
	  meanActual+=specY[i];
	}
      px[i]=log10(px[i]);
      py[i]=log10(py[i]);
    }
  printf("Mean spectrum for f > 1/tau = %g\n",meanActual/(double)nv);

  minx = TKfindMin_f(px,specN);
  maxx = TKfindMax_f(px,specN);
  miny = TKfindMin_f(py,specN);
  miny = -12;

  maxy = TKfindMax_f(py,specN);
  cpgend();
  cpgbeg(0,"/cps",1,1);
  cpgslw(4);
  //  cpgeras();
  cpgswin(minx,maxx,miny,maxy+0.1*(maxy-miny));  
  cpgbox("BCNST1L",0.0,0,"BCNST1L",0.0,0);
  cpglab("Frequency (d\\u-1\\d)","Power","");
  cpgline(specN,px,py);

  fx[0] = fx[1] = log10(1.0/365.25);
  fy[0] = miny; fy[1] = maxy+0.1*(maxy-miny);
  cpgsci(8); cpgsls(3); cpgline(2,fx,fy); cpgsls(1); cpgsci(1);

  // Draw slope of -2, -4 and -6
  fx[0] = minx; fx[1] = maxx;
  fy[0] = maxy; fy[1] = -2*fx[1]+(fy[0]-(-2*fx[0]));
  cpgsls(4); cpgline(2,fx,fy); cpgsls(1);
  fy[0] = maxy; fy[1] = -4*fx[1]+(fy[0]-(-4*fx[0]));
  cpgsls(4); cpgline(2,fx,fy); cpgsls(1);
  fy[0] = maxy; fy[1] = -6*fx[1]+(fy[0]-(-6*fx[0]));
  cpgsls(4); cpgline(2,fx,fy); cpgsls(1);

  printf("Colours\n\n");
  printf("white = power spectrum of post-fit timing residuals\n");

  // Interpolate to split red and white noise
  {
     double interpX[MAX_OBSN],interpY[MAX_OBSN],interpE[MAX_OBSN],tx;
    double specWX[MAX_OBSN],specWY[MAX_OBSN];
    double white[MAX_OBSN];
    double w,sw;
    int specWN;
    int nAwhite;

    printf("Obtaining a good model of the white noise\n");

    //    tau = 30.0;
    // First interpolate onto the same grid as the actual data
    for (i=0;i<psr[0].nobs;i++)
      interpX[i] = (psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
    for (j=0;j<psr[0].nobs;j++)
      {
	sw = 0.0;
	interpY[j] = 0.0;
	for (i=0;i<psr[0].nobs;i++)
	  {
	    tx = (double)(psr[0].obsn[i].sat - psr[0].param[param_pepoch].val[0]);
	    w = 1.0/pow(psr[0].obsn[i].toaErr*1.0e-6,2)*exp(-fabs(tx-interpX[j])/tau);
	    sw += w;
	    interpY[j] += psr[0].obsn[i].residual*w;
	  }
	interpY[j]/=sw;
	interpE[j] = psr[0].obsn[j].toaErr*1.0e-6; // Not using weights in the spectral analysis
	white[j] = psr[0].obsn[j].residual-interpY[j];
	//	printf("interp %g %g %g %g\n",interpX[j],interpY[j],(double)psr[0].obsn[j].residual,white[j]);
      }

    //    exit(1);
    // Get a spectrum of this white noise
    TKspectrum(interpX,white,interpE,psr[0].nobs,0,0,0,0,0,specType,ofac,1,specOut,specWX,specWY,&specWN,0,0,outY_re,outY_im);    
    TKconvertFloat2(specWX,specWY,px,py,specWN);
    for (i=0;i<specWN;i++)
      {
	px[i]=log10(px[i]);
	py[i]=log10(py[i]);
	printf("Have %g %g\n",px[i],py[i]);
      }
    cpgsci(2); cpgline(specWN,px,py); cpgsci(1);
    fx[0]=fx[1]=log10(1.0/tau);
    fy[0]=miny;fy[1]=maxy+0.1*(maxy-miny);
    cpgsls(4); cpgsci(2); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);


    // Calculate mean actual white level
    meanActualWhiteLevel = 0.0;
    nAwhite=0;
    for (i=0;i<specWN;i++)
      {
	//	printf("have %g %g\n",specWX[i],1.0/tau);

	if (specWX[i] > 1.0/tau)
	  {
	    meanActualWhiteLevel+=specWY[i];
	    nAwhite++;
	  }
      }
    meanActualWhiteLevel/=(double)nAwhite;
    area1 = meanActualWhiteLevel;
    printf("mean measured white level = %g %d\n",meanActualWhiteLevel,nAwhite);
  }
  
  // Obtain white level
  /*  {
    long double sat0[MAX_OBSN];
    double avSpecY[MAX_OBSN];
    int count;

    for (i=0;i<psr[0].nobs;i++)
      avSpecY[i]=0.0;
      // Form ideal SATs

    printf("Starting white\n");
    for (j=0;j<5;j++)
      {
	psr[0].nJumps = 0;
	for(i=0;i<MAX_PARAMS;i++){
	  psr[0].param[i].nLinkTo = 0;
	  psr[0].param[i].nLinkFrom = 0;
	}
	readParfile(psr,parFile,timFile,1); 
	formBatsAll(psr,1);         
	formResiduals(psr,1,0);   
	for (i=0;i<psr[0].nobs;i++)
	  psr[0].obsn[i].sat -= (long double)psr[0].obsn[i].residual/86400.0L;
      }
    printf("Getting sat0\n");
    for (i=0;i<psr[0].nobs;i++)
      sat0[i] = psr[0].obsn[i].sat;

    for (it=0;it<nit;it++)
      {
	printf("iteration %d/%d\n",it+1,nit);
	for (i=0;i<psr[0].nobs;i++)
	  psr[0].obsn[i].sat = sat0[i] + (TKgaussDev(&idum)*psr[0].obsn[i].toaErr*1.0e-6)/SECDAY;
	psr[0].nJumps = 0;
	for(i=0;i<MAX_PARAMS;i++){
	  psr[0].param[i].nLinkTo = 0;
	  psr[0].param[i].nLinkFrom = 0;
	}
	readParfile(psr,parFile,timFile,1); 
	formBatsAll(psr,1);                 
	formResiduals(psr,1,0);             
	doFit(psr,1,0);                     
	formBatsAll(psr,1);                 
	formResiduals(psr,1,0);             
	n=0;
	for (i=0;i<psr[0].nobs;i++)
	  {
	    if (psr[0].obsn[i].deleted==0)
	      {
		resx[n] = (double)(psr[0].obsn[i].sat-psr[0].param[param_pepoch].val[0]);
		resy[n] = (double)psr[0].obsn[i].residual;
		rese[n] = (double)psr[0].obsn[i].toaErr*1.0e-6;
		n++;
	      }
	  }
	
	TKspectrum(resx,resy,rese,n,0,0,0,0,0,specType,ofac,1,specOut,specX,specY,&specN,0,0,outY_re,outY_im);
	for (i=0;i<specN;i++)
	  avSpecY[i] += specY[i];
      }
    for (i=0;i<specN;i++)
      avSpecY[i]/=(double)nit;
    TKconvertFloat2(specX,avSpecY,px,py,specN);

    for (i=0;i<specN;i++)
      {
	px[i]=log10(px[i]);
	py[i]=log10(py[i]);
      }
    cpgsci(3);
    cpgline(specN,px,py);

    // Now form white model
    count=0;
    whiteLevel=0.0;
    for (i=0;i<specN;i++)
      {
	if (specX[i] > 1.0/tau)
	  {
	    count++;
	    whiteLevel+=(specY[i]);
	  }
      }
    whiteLevel/=(double)(count);
    area2 = whiteLevel;
    printf("whiteLevel = %g\n",whiteLevel);
    printf("EFAC = %g\n",sqrt(area1/area2));
    for (i=0;i<specN;i++)
      py[i] = log10(whiteLevel);
    cpgsls(4);
    cpgline(specN,px,py);
    for (i=0;i<specN;i++) fy[i] = log10(pow(10,py[i])*3.67);
    cpgline(specN,px,fy);
    for (i=0;i<specN;i++) fy[i] = log10(pow(10,py[i])*0.025);
    cpgline(specN,px,fy);
    cpgsls(1);
  }*/

  // Now process the red component

  // Resample onto daily grid
  double start = (double)(psr[0].obsn[0].sat-psr[0].param[param_pepoch].val[0]);
  int interpN=(int)((psr[0].obsn[psr[0].nobs-1].sat-psr[0].obsn[0].sat)+0.5);
    double interpX[MAX_OBSN],interpY[MAX_OBSN],interpE[MAX_OBSN],tx;
    double specWX[MAX_OBSN],specWY[MAX_OBSN];
    double white[MAX_OBSN];
    double w,sw;
    int specWN;
    int nAwhite;

  for (i=0;i<interpN;i++)
    interpX[i] = start+i;
  for (j=0;j<interpN;j++)
    {
      sw = 0.0;
      interpY[j] = 0.0;
      for (i=0;i<psr[0].nobs;i++)
	{
	  tx = (double)(psr[0].obsn[i].sat - psr[0].param[param_pepoch].val[0]);
	  w = 1.0/pow(psr[0].obsn[i].toaErr*1.0e-6,2)*exp(-fabs(tx-interpX[j])/tau);
	  sw += w;
	  interpY[j] += psr[0].obsn[i].residual*w;
	}
      interpY[j]/=sw;
      interpE[j] = 0.0;
      //      printf("interp %g %g %g %g\n",interpX[j],interpY[j]);
    }
    TKspectrum(interpX,interpY,interpE,interpN,0,0,0,0,0,specType,ofac,1,specOut,specWX,specWY,&specWN,0,0,outY_re,outY_im);    
    TKconvertFloat2(specWX,specWY,px,py,specWN);
    for (i=0;i<specWN;i++)
      {
	px[i]=log10(px[i]);
	py[i]=log10(py[i]);
      }
    cpgsci(3); cpgline(specWN,px,py); cpgsci(1);

    TKspectrum(interpX,interpY,interpE,interpN,0,0,0,0,1,specType,ofac,1,specOut,specWX,specWY,&specWN,0,0,outY_re,outY_im);    
    TKconvertFloat2(specWX,specWY,px,py,specWN);
    for (i=0;i<specWN;i++)
      {
	px[i]=log10(px[i]);
	py[i]=log10(py[i]);
      }
    cpgsci(4); cpgline(specWN,px,py); cpgsci(1);

    TKspectrum(interpX,interpY,interpE,interpN,0,0,0,0,2,specType,ofac,1,specOut,specWX,specWY,&specWN,0,0,outY_re,outY_im);    
    TKconvertFloat2(specWX,specWY,px,py,specWN);
    for (i=0;i<specWN;i++)
      {
	px[i]=log10(px[i]);
	py[i]=log10(py[i]);
      }
        cpgsci(5); cpgline(specWN,px,py); cpgsci(1);
  

  cpgend();
  exit(1);
}
char * plugVersionCheck = TEMPO2_h_VER;
