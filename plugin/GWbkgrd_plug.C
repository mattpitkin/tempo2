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

/* This plugin allows the user to simulate a gravitational wave background */
/*
 * Inputs: 
 *  - ag:    amplitude of GW background
 *  - alpha: spectral exponent of the characteristic strain spectrum
 *
 * 
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "T2toolkit.h"
#include "GWsim.h"
#include "TKfit.h"
#include <cpgplot.h>

using namespace std;

void doPlot(pulsar *psr,int npsr,gwSrc *gw,long double **gwRes,long double timeOffset,int ngw,long double tspan,long double alpha,long double gwAmp);
long double getTspan(pulsar *psr,int npsr);
void plotResiduals(pulsar *psr,long double **gwRes,int p,long double timeOffset,int plotType);
void plotSpectrum(gwSrc *gw,int ngw,long double tspan,long double alpha,long double gwAmp);
void plotPosn(pulsar *psr,int npsr,gwSrc *gw,int ngw);
void draw_grid(double start_gl,double end_gl,double start_gb,double end_gb,double gstep,double bstep,int celestialCoords);
void convertXY_celestial(double raj,double decj,double *retx,double *rety);

void help() /* Display help */
{
  printf("--------------------------------------\n");
  printf("Command line inputs:\n\n");
  printf("-addWhite add white noise based on the TOA error bar size\n");
  printf("-alpha GW spectral exponent (usually -0.666) \n");
  printf("-clock simulate clock errors instead of a GWB\n");
  printf("-dist  pulsar distance in kpc\n");
  printf("-f     parfile.par timfile.tim:  input par and tim files\n");
  printf("-flo   lowest GW frequency to simulate (default = 0.01/Tspan) (Hz) \n");
  printf("-fhi   highest GW frequency to simulate (default = d^-1) (Hz)\n");
  printf("-gwamp GWamplitude in dimensionless units\n");
  printf("-h     this help file\n");
  printf("-linear Use linear spacing for GW frequencies (default = logarithmic spacing)\n");
  printf("-ngw   Number of gravitational waves to simulate\n");
  printf("-plot  Plot the GW positions and induced residuals\n");
  printf("-seed  Set random number seed (choose negative integer)\n");
  printf("-zero  Set all the post-fit residuals to zero before continuing\n");
  printf("\n\n");
  printf("Example usage: tempo2 -gr GWbkgrd -f tt.par 0437.2048.tim -dist 1 -gwamp 1e-15 -alpha -0.666666 -ngw 1000 -plot\n\n");
  printf("Plot options\n\n");
  printf("1 Plot sky plot with pulsar position indicated and all GW source positions\n");
  printf("2 Plot GW power versus GW frequency\n");
  printf("3 Plot the induced residuals caused by the GW background\n");
  printf("4 Plot the induced residuals caused by the GW backgroudn after quadratic removal\n");
  printf("5 Plot the pre-fit timing residuals\n");
  printf("6 Plot the post-fit timing residuals\n");
  printf("h This help file\n");
  printf("q Quit\n");
  printf("s Save residuals and site arrival times\n");
  printf("p Select which pulsar number to plot\n");
  printf("--------------------------------------\n");
  
}


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  FILE *fout,*fout2;
  int i,j,k,p;
  double globalParameter;
  int setgwAmp=0,setAlpha=0;
  int plotIt=0;
  long double gwAmp,alpha;
  long double timeOffset;
  long double ra_p,dec_p;
  long double flo=0.0,fhi=0.0;
  long double kp[3];            /* Vector pointing to pulsar           */
  long double tspan;
  long double time;
  long double **gwRes;
  long double dist[MAX_PSR];
  long double mean;
  int clock=0;
  int distNum=0;
  int logspacing=1;
  int ngw=0;
  int addWhite=0;
  long seed=TKsetSeed();
  int zeroResiduals=0;
  double scale;
  char fname[100];
  int outGW=0;
  int identicalTimes=0;
  gwSrc *gw;

  *npsr = 0;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: GWbkgrd\n");
  printf("Author:              G. Hobbs\n");
  printf("Version:             1.0\n");
  printf(" --- type 'h' for help information\n");


  /* Obtain the .par and the .tim file from the command line */
  if (argc==4) /* Only provided .tim name */
    {
      strcpy(timFile[0],argv[3]);
      strcpy(parFile[0],argv[3]);
      parFile[0][strlen(parFile[0])-3] = '\0';
      strcat(parFile[0],"par");
    }

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)  // Read parameter file and arrival time files
	{
 	  strcpy(parFile[*npsr],argv[i+1]); 
	  strcpy(timFile[*npsr],argv[i+2]);	  
	  (*npsr)++;
	}
      else if (strcmp(argv[i],"-identicalTimes")==0)
	identicalTimes=1;
      else if (strcmp(argv[i],"-outGW")==0)
	outGW=1;
      else if (strcmp(argv[i],"-clock")==0)
	clock=1;
      else if (strcmp(argv[i],"-addWhite")==0)
	{
	  printf("Using TOA error bars to add white noise\n");
	  addWhite=1;
	}
      else if (strcmp(argv[i],"-h")==0)
	{
	  help();
	  exit(1);
	}
      else if (strcmp(argv[i],"-dist")==0) // Distance in kpc
	{
	  sscanf(argv[++i],"%Lf",&dist[distNum]);
	  dist[distNum]*=3.086e19;
	  distNum++;
	}
      else if (strcmp(argv[i],"-gwamp")==0)
	{sscanf(argv[++i],"%Lf",&gwAmp); setgwAmp=1;}
      else if (strcmp(argv[i],"-alpha")==0)
	{sscanf(argv[++i],"%Lf",&alpha); setAlpha=1;}
      else if (strcmp(argv[i],"-ngw")==0)
	{sscanf(argv[++i],"%d",&ngw);}
      else if (strcmp(argv[i],"-plot")==0)
	plotIt=1;
      else if (strcmp(argv[i],"-flo")==0)
	sscanf(argv[++i],"%Lf",&flo);
      else if (strcmp(argv[i],"-zero")==0)
	zeroResiduals=1;
      else if (strcmp(argv[i],"-fhi")==0)
	sscanf(argv[++i],"%Lf",&fhi);
      else if (strcmp(argv[i],"-seed")==0)
	sscanf(argv[++i],"%d",&seed);
      else if (strcmp(argv[i],"-linear")==0)
	logspacing=0;
    }
  

  // Check that all the parameters are set
  if (setgwAmp==0)
    {
      gwAmp = 1e-20;
      printf("WARNING: automatically setting the GW amplitude to be 1e-20\nUse -gwamp option to set manually.\n");
    }

  if (setAlpha==0)
    {
      alpha=-2.0/3.0;
      printf("WARNING: automatically setting the characteristic strain exponent to -2/3\nUse -alpha option to set manually\n");
    }
  scale = pow(86400.0*365.25,alpha);
  gwAmp *= scale;


  if (ngw==0)
    {
      printf("WARNING: automatically setting the number of gravitational waves to 1000\nUse -ngw option to set manually\n");
      ngw=1000;
    }
  if (*npsr==0)
    {
      printf("ERROR: Must use -f option to provide a pulsar timing model\n");
      exit(1);
    }
  if (distNum!=*npsr)
    {
      printf("ERROR: Distances not provided for all the pulsars: Npsr = %d, Ndist = %d\nUse -dist to provide distances (in kpc) on the command line.\n",*npsr,distNum);
      exit(1);
    }

  if ((gw = (gwSrc *)malloc(sizeof(gwSrc)*ngw))==NULL)
    {
      printf("Unable to allocate memory for %d GW sources\n",ngw);
      exit(1);
    }
  gwRes = (long double **)malloc(MAX_PSR*sizeof(long double*));
  for (i=0;i<MAX_PSR;i++)
    gwRes[i] = (long double *)malloc(MAX_OBSN*sizeof(long double));

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr);         /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);
  formBatsAll(psr,*npsr);                 /* Form the barycentric arrival times */
  formResiduals(psr,*npsr,1);           /* Form the residuals                 */

  if (zeroResiduals==1)
    {
      for (j=0;j<5;j++)
	{
	  for (p=0;p<*npsr;p++)
	    {
	      for (i=0;i<psr[p].nobs;i++)
		psr[p].obsn[i].sat -= psr[p].obsn[i].residual/SECDAY;
	    }
	  formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
	  formResiduals(psr,*npsr,1);   /* Form the residuals                 */
	}
    }


  // Set range of frequencies for GW simulation
  tspan=getTspan(psr,*npsr)*SECDAY;
  if (flo==0)
    {
      flo=0.01/tspan;
      printf("flo = %.5Lg, tspan = %.5Lg\n",flo,tspan);
    }
  if (fhi==0)
    {
      fhi = 1.0/(long double)SECDAY;
      printf("fhi = %.5Lg\n",fhi);
    }

  timeOffset = psr[0].param[param_pepoch].val[0];
  GWbackground(gw,ngw,&seed,flo,fhi,gwAmp,alpha,logspacing);

  for (i=0;i<ngw;i++)
    {
      setupGW(&gw[i]);
    }
  if (clock==1)
    fout = fopen("signal.dat","w");
  for (p=0;p<*npsr;p++)
    {
      if (outGW==1)
	{
	  char str[128];
	  sprintf(str,"%s.gw",psr[p].name);
	  fout2 = fopen(str,"w");
	}

      if (clock==0)
	{
	  ra_p   = psr[p].param[param_raj].val[0];
	  dec_p  = psr[p].param[param_decj].val[0];
	}
      else
	{
	  ra_p   = psr[0].param[param_raj].val[0];
	  dec_p  = psr[0].param[param_decj].val[0];
	}
      setupPulsar_GWsim(ra_p,dec_p,kp);
      mean=0.0;
      for (i=0;i<psr[p].nobs;i++) 
	{
	  if (identicalTimes==1)
	    time = (psr[0].obsn[i].sat - timeOffset)*SECDAY;
	  else
	    time = (psr[p].obsn[i].sat - timeOffset)*SECDAY;
	  gwRes[p][i] = 0.0;
	  for (k=0;k<ngw;k++)
	    gwRes[p][i]+=calculateResidualGW(kp,&gw[k],time,dist[p]);
	  mean += gwRes[p][i];
	}
      mean /= (long double)psr[p].nobs;
      for (i=0;i<psr[p].nobs;i++)
	{
	  if (clock==1)
	    fprintf(fout,"%.10f %.10g\n",(double)psr[p].obsn[i].sat,(double)gwRes[p][i]);
	  if (outGW==1)
	    {
	      if (identicalTimes==1)
		fprintf(fout2,"%.10f %.10g\n",(double)psr[0].obsn[i].sat,(double)gwRes[p][i]);
	      else
		fprintf(fout2,"%.10f %.10g\n",(double)psr[p].obsn[i].sat,(double)gwRes[p][i]);
	    }
	  gwRes[p][i]-=mean;
	  psr[p].obsn[i].sat += (gwRes[p][i]/(long double)SECDAY);
	  if (addWhite==1)
	    psr[p].obsn[i].sat += (long double)TKgaussDev(&seed)*(psr[p].obsn[i].toaErr*1.0e-6)/(long double)SECDAY;
	}
    }
  if (clock==1)
    fclose(fout);
  if (outGW==2)
    fclose(fout2);
  for (p=0;p<*npsr;p++)
    {
      sprintf(fname,"%s.gwsim.tim",psr[p].name);
      writeTim(fname,psr+p,"tempo2");
    }
  if (plotIt==1)
    {
      formBatsAll(psr,*npsr);                 /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);           /* Form the residuals                 */
      doFit(psr,*npsr,0);
      formBatsAll(psr,*npsr);                 /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);           /* Form the residuals                 */
      
      doPlot(psr,*npsr,gw,gwRes,timeOffset,ngw,tspan,alpha,gwAmp);
    }
  return 0;
}

void doPlot(pulsar *psr,int npsr,gwSrc *gw,long double **gwRes,long double timeOffset,int ngw,long double tspan,long double alpha,long double gwAmp)
{
  int plot=1;
  int pulsar=0;
  int endit=0;
  int i,j;
  float mx,my;
  char key;
  char fname[100];
  FILE *fout;

  cpgbeg(0,"/xs",1,1);
  cpgpap(0,0.5);
  cpgask(0);
  do {
  if (plot==3 || plot==4 || plot==5 || plot==6)
    plotResiduals(psr,gwRes,pulsar,timeOffset,plot-2);
  else if (plot==2)
    plotSpectrum(gw,ngw,tspan,alpha,gwAmp);
  else if (plot==1)
    plotPosn(psr,npsr,gw,ngw);

  cpgcurs(&mx,&my,&key);
  if (key=='q') endit=1;
  else if (key=='1') plot=1;
  else if (key=='2') plot=2;
  else if (key=='5') plot=5;
  else if (key=='6') plot=6;
  else if (key=='3')
    plot=3;
  else if (key=='4')
    plot=4;
  else if (key=='h')
    help();
  else if (key=='s') // Save residuals and site arrival times
    {
      for (j=0;j<npsr;j++)
	{
	  printf("Writing new site arrival times\n");
	  sprintf(fname,"%s.gwsim.tim",psr[j].name);
	  writeTim(fname,psr+j,"tempo2");
	  printf("Writing new residual file\n");
	  sprintf(fname,"%s.res",psr[j].name);
	  fout = fopen(fname,"w");
	  for (i=0;i<psr[j].nobs;i++)
	    {
	      fprintf(fout,"%.15Lf %g %g\n",psr[j].obsn[i].sat,(double)psr[j].obsn[i].residual,(double)psr[j].obsn[i].toaErr);
	    }
	  fclose(fout);
	}
    }
  else if (key=='p') // Select pulsar
    {
      printf("Enter pulsar number (from 0 to %d)\n",npsr-1);
      scanf("%d",&pulsar);
    }

  } while (endit==0);
  cpgend();
}

void plotPosn(pulsar *psr,int npsr,gwSrc *gw,int ngw)
{
  int i;
  double px[ngw],py[ngw];
  double rad2deg = 180.0/M_PI;
  float fx[ngw],fy[ngw];

  draw_grid(-180,180,-90,90,60,30,1);
  cpglab("","","Gravitational wave source and pulsar positions");
  // Plot the GW source positions
  for (i=0;i<ngw;i++)
    {
      convertXY_celestial((double)(gw[i].phi_g*rad2deg)-180,
			  (double)(gw[i].theta_g*rad2deg)-90,&px[0],&py[0]);
      fx[0] = (float)px[0];
      fy[0] = (float)py[0];
      cpgpt(1,fx,fy,1); 
    }

  // Plot the pulsar positions  
  for (i=0;i<npsr;i++)
    {
      convertXY_celestial((double)(psr[i].param[param_raj].val[0]*rad2deg)-180,
			  (double)psr[i].param[param_decj].val[0]*rad2deg,&px[0],&py[0]);
      fx[0] = (float)px[0];
      fy[0] = (float)py[0];
      cpgsch(2); cpgsci(2); cpgpt(1,fx,fy,12); cpgsci(1); cpgsch(1);
    }
}

void plotSpectrum(gwSrc *gw,int ngw,long double tspan,long double alpha,long double gwAmp)
{
  float maxx,minx,maxy,miny;
  float fx[ngw],fy[ngw];
  int i;
  for (i=0;i<ngw;i++)
    {
      fx[i] = (float)(log10((gw[i].omega_g)/(2.0*M_PI)));
      fy[i] = (float)(log10(pow(gw[i].aplus_g,2)+pow(gw[i].across_g,2)));
    }
  minx = TKfindMin_f(fx,ngw);
  maxx = TKfindMax_f(fx,ngw);
  miny = TKfindMin_f(fy,ngw);
  maxy = TKfindMax_f(fy,ngw);
  cpgenv(minx-1,maxx+1,miny,maxy,0,1);
  cpglab("log\\d10\\u[f\\dg\\u] (Hz)","log\\d10\\u[A\\d+\\u(f)\\u2\\d+A\\dx\\u(f)\\u2\\d]","");
  cpgpt(ngw,fx,fy,1);
  fx[0] = log10(1.0/86400.0); fx[1] = fx[0];
  fy[0] = miny; fy[1] = maxy;
  cpgsci(2); cpgsls(4); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
  fx[0] = log10(1.0/tspan); fx[1] = fx[0];
  fy[0] = miny; fy[1] = maxy;
  cpgsci(2); cpgsls(4); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
  
  // Draw line with slope of 2*"alpha"
  //  fx[0] = minx-1;
  //  fy[0] = maxy;
  //  fx[1] = maxx+1;
  //  fy[1] = 2.0*alpha*fx[1]+(fy[0]-2.0*alpha*fx[0]);
  
  fx[0] = minx-1;
  fx[1] = maxx+1;
  fy[0] = log10(gwAmp*gwAmp/12.0/M_PI/M_PI*pow(pow(10,fx[0]),2*alpha));
  fy[1] = log10(gwAmp*gwAmp/12.0/M_PI/M_PI*pow(pow(10,fx[1]),2*alpha));
  printf("line from %g %g to %g %g\n",fx[0],fy[0],fx[1],fy[1]);
  cpgsci(3); cpgsls(3); cpgline(2,fx,fy); cpgsci(1); cpgsls(1);
}

void plotResiduals(pulsar *psr,long double **gwRes,int p,long double timeOffset,int plotType)
{
  float px[MAX_OBSN],py[MAX_OBSN];
  float minx,maxx,miny,maxy;
  float yerr1[MAX_OBSN],yerr2[MAX_OBSN];
  int i;
  for (i=0;i<psr[p].nobs;i++)
    {
      px[i] = (float)(psr[p].obsn[i].sat - timeOffset);
      if (plotType==1 || plotType==2)
	py[i] = (float)gwRes[p][i];
      else if (plotType==3)
	py[i] = (float)psr[p].obsn[i].prefitResidual;
      else if (plotType==4)
	py[i] = (float)psr[p].obsn[i].residual;
    }
  TKremovePoly_f(px,py,psr[p].nobs,1);
  for (i=0;i<psr[p].nobs;i++)
    {
      yerr1[i] = py[i] - psr[p].obsn[i].toaErr*1e-6;
      yerr2[i] = py[i] + psr[p].obsn[i].toaErr*1e-6;
    }
  if (plotType==2)
    TKremovePoly_f(px,py,psr[p].nobs,3);
  minx = TKfindMin_f(px,psr[p].nobs);
  maxx = TKfindMax_f(px,psr[p].nobs);
  miny = TKfindMin_f(py,psr[p].nobs);
  maxy = TKfindMax_f(py,psr[p].nobs);
  cpgenv(minx,maxx,miny,maxy,0,1);
  if (plotType==1)
    cpglab("Day","GW residual (s)","Induced residuals due to GW background");
  else if (plotType==2)
    cpglab("Day","GW residual (s)","Induced residuals due to GW background after quadratic removed");
  else if (plotType==3)
    cpglab("Day","GW residual (s)","Pre-fit timing residuals");
  else if (plotType==4)
    cpglab("Day","GW residual (s)","Post-fit timing residuals");
    
  cpgpt(psr[p].nobs,px,py,9);
  if (plotType==3 || plotType==4)
    cpgerry(psr[p].nobs,px,yerr1,yerr2,1);
}

long double getTspan(pulsar *psr,int npsr)
{
  long double first,last;
  int i,p;
    
  
  first = psr[0].obsn[0].sat;
  last = psr[0].obsn[0].sat;

  for (p=0;p<npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)
	{
	  if (first > psr[p].obsn[i].sat)
	    first = psr[p].obsn[i].sat;
	  if (last < psr[p].obsn[i].sat)
	    last = psr[p].obsn[i].sat;
	}
    }

  return last-first;
}

void draw_grid(double start_gl,double end_gl,double start_gb,double end_gb,double gstep,double bstep,int celestialCoords)
{
  double l,b,x,y;
  float plotx[1000],ploty[1000];
  int count=0;
  char str[100];
  cpgenv(start_gl,end_gl,start_gb,end_gb,0,-2);
  cpgsls(4);

  /* Plot lines of latitude */
  for (b=start_gb;b<=end_gb;b+=bstep)
    {
      count=0;
      for (l=start_gl;l<=end_gl;l=l+1.0)
	{
	  if (celestialCoords==1) convertXY_celestial(l,b,&x,&y);
	  /*get_xy(l,b,&x,&y); */
	  plotx[count] = (float)x;
	  ploty[count] = (float)y;
	  /*	  printf("%d %f %f\n",count,plotx[count],ploty[count]); */
	  count++;
	}
      cpgline(count,plotx,ploty);
    }

  /* Plot lines of longitude */
  for (l=start_gl;l<=end_gl;l+=gstep)
    {
      count=0;
      for (b=start_gb;b<=end_gb;b=b+1.0)
	{
	  if (celestialCoords==1) convertXY_celestial(l,b,&x,&y);
	  /*	  get_xy(l,b,&x,&y); */
	  plotx[count] = (float)x;
	  ploty[count] = (float)y;
	  count++;
	}
      if (l==-180 || l==180)
	cpgsls(1);
      else
	cpgsls(4);
      cpgline(count,plotx,ploty);
    }
  

  /* Label axes */
  cpgsci(3);
  for (l=0;l<360;l+=gstep)
    {
      if (celestialCoords==1) convertXY_celestial(l-180,-45,&x,&y);

      /*      if (celestialCoords==1)
	get_xy(l+180.0,-45,&x,&y);
      else
      get_xy(l,-45,&x,&y); */
      if (l!=180.0 || celestialCoords==1)
	{
	  if (celestialCoords==0 || l!=0)
	    {
	      if (celestialCoords==0) sprintf(str,"%.0f\\uo\\d",l);
	      else sprintf(str,"%.0f\\uh\\d",l/360.0*24.0);
	      cpgptxt((float)x,(float)y,0,0.5,str);
	    }
	}
    }
  for (b=-60;b<=60;b+=bstep)
    {
      if (celestialCoords==1) convertXY_celestial(-180,b,&x,&y);
      /*      get_xy(180,b,&x,&y); */
      if (b>0)
	{
	  sprintf(str,"+%.0f\\uo\\d",b);
	  cpgptxt((float)x,(float)y,0,1.0,str);
	}
      else if (b==0)
	{
	  sprintf(str,"%.0f\\uo\\d",b);
	  cpgptxt((float)x-2,(float)y,0,1.0,str);
	}
      else
	{
	  sprintf(str,"%.0f\\uo\\d",b);
	  cpgptxt((float)x,(float)y-7,0,1.0,str);
	}
    }
  cpgsci(1);
  cpgsls(1);
}

/* Convert from RAJ, DECJ to x,y using Aitoff projection */
void convertXY_celestial(double raj,double decj,double *retx,double *rety)
{
  double sa;
  double r2deg = 180.0/M_PI;
  double alpha2,delta;
  double r2,f,cgb,denom;
  double x_ret,y_ret;

  sa = raj;
  alpha2 = sa/(2*r2deg);
  delta = decj/r2deg;   

  r2 = sqrt(2.);    
  f = 2*r2/M_PI;    

  cgb = cos(delta);    
  denom =sqrt(1. + cgb*cos(alpha2));

  x_ret = cgb*sin(alpha2)*2.*r2/denom;
  y_ret = sin(delta)*r2/denom;

  x_ret = x_ret*r2deg/f;
  y_ret = y_ret*r2deg/f;

  *retx = x_ret;
  *rety = y_ret;
}
char * plugVersionCheck = TEMPO2_h_VER;
