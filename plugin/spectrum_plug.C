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

#define MAX_ID 50

void drawOption(int on,float x,float y,char *str);
void doPlugin(pulsar *psr,int npsr);
void drawMenu(int specType,int xaxis,int logv,int specOut);
void checkMenu(float mx,float my,int *change,int *xaxis,int *logv,int *specType,int *specOut);
void identify(float mx,float my,float px[MAX_PSR_VAL][MAX_OBSN_VAL],float py[MAX_PSR_VAL][MAX_OBSN_VAL],
	      int *sn, int idV[MAX_ID],int idP[MAX_ID], int *iN,int npsr);

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
  int i,n=0;
  double px[MAX_OBSN],py[MAX_OBSN],pe[MAX_OBSN];

  *npsr = 0;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: spectrum\n");
  printf("Author:              G. Hobbs\n");
  printf("Version:             1.0\n");
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
    }

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      if (i==0) doFit(psr,*npsr,0);   /* Do the fitting     */
      else textOutput(psr,*npsr,0,0,0,0,"");  /* Display the output */
    }

  doPlugin(psr,*npsr);

  return 0;
}

void doPlugin(pulsar *psr,int npsr)
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

  do {
    cpgbeg(0,"/xs",1,1);
    cpgsch(fontSize);
    cpgask(0);
    cpgeras();
    
    cpgsvp(0.0,1.0,0.8,1.0);
    cpgswin(0,1,0,1);
    drawMenu(specType,xaxis,logv,specOut);
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
		    //	  py[n] = 3*sin(2.0*M_PI/2041.34587909003*10*px[n]);
		    n++;
		  }
	      }

	    TKspectrum(resx,resy,rese,n,0,0,0,0,0,specType,ofac,specOut,specX,specY,&specN[p],0,0);
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
    if (logv==0)
      cpgbox("BCNST1",0.0,0,"BCNST1",0.0,0);
    else if (logv==1)
      cpgbox("BCNST1L",0.0,0,"BCNST1",0.0,0);
    else if (logv==2)
      cpgbox("BCNST1",0.0,0,"BCNST1L",0.0,0);
    else if (logv==3)
      cpgbox("BCNST1L",0.0,0,"BCNST1L",0.0,0);
    
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
