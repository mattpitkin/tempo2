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

/* This plugin allows the user to simulate a single GW source*/
/*
 * Inputs: 
 *
 * 
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "GWsim.h"
#include "T2toolkit.h"
#include "TKfit.h"
#include "TKspectrum.h"
#include <cpgplot.h>

using namespace std;

void doPlot(pulsar *psr,int npsr,gwSrc gw,long double **gwRes,long double timeOffset,long double tspan);
long double getTspan(pulsar *psr,int npsr);
void plotResiduals(pulsar *psr,long double **gwRes,int p,long double timeOffset,int plotType);
void plotSpectrum(pulsar *psr,int p,long double timeOffset);
void plotPosn(pulsar *psr,int npsr,gwSrc gw);
void draw_grid(double start_gl,double end_gl,double start_gb,double end_gb,double gstep,double bstep,int celestialCoords);
void convertXY_celestial(double raj,double decj,double *retx,double *rety);

void help() /* Display help */
{
  printf("----------------------------------\n");
  printf("Command line options\n\n");
  printf("-dist   Distance in kpc\n");
  printf("-f parfile.par timfile.tim Reads in .par and .tim file\n");
  printf("-h      This help\n");
  printf("-plot   Interactive plotting\n");
  printf("-across GW Ax\n");
  printf("-aplus  GW A+\n");
  printf("-omega  GW frequency\n");
  printf("-gra    GW right ascension (hours)\n");

  printf("----------------------------------\n");
}


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,j,k,p;
  double globalParameter;
  int plotIt=0;
  int toas=0;
  long double timeOffset;
  long double ra_p,dec_p;
  long double flo=0.0,fhi=0.0;
  long double kp[3];            /* Vector pointing to pulsar           */
  long double tspan;
  long double time;
  long double **gwRes;
  long double dist[MAX_PSR];
  long double mean;
  long double aplus,across;
  long double omega;
  long double gra,gdec;
  int setAplus=0,setAcross=0,setOmega=0,setGra=0,setGdec=0;
  int distNum=0;
  int seed=-431;

  gwSrc gw;

  *npsr = 0;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: GWsingle\n");
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
      else if (strcmp(argv[i],"-plot")==0)
	plotIt=1;
      else if (strcmp(argv[i],"-toas")==0)
	toas=1;
      else if (strcmp(argv[i],"-aplus")==0)
	{ sscanf(argv[++i],"%Lf",&aplus); setAplus=1;}
      else if (strcmp(argv[i],"-across")==0)
	{ sscanf(argv[++i],"%Lf",&across); setAcross=1;}
      else if (strcmp(argv[i],"-omega")==0)
	{ sscanf(argv[++i],"%Lf",&omega); setOmega=1;}
      else if (strcmp(argv[i],"-gra")==0) // hr
	{ sscanf(argv[++i],"%Lf",&gra); setGra=1;}
      else if (strcmp(argv[i],"-gdec")==0) // Degrees
	{ sscanf(argv[++i],"%Lf",&gdec); setGdec=1;}
      else if (strcmp(argv[i],"-seed")==0)
	sscanf(argv[++i],"%d",&seed);
    }

  // Check that all the parameters are set
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

  gwRes = (long double **)malloc(MAX_PSR*sizeof(long double*));
  for (i=0;i<MAX_PSR;i++)
    gwRes[i] = (long double *)malloc(MAX_OBSN*sizeof(long double));
  
  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr);         /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);
  formBatsAll(psr,*npsr);                 /* Form the barycentric arrival times */
  formResiduals(psr,*npsr,0);           /* Form the residuals                 */

  timeOffset = psr[0].param[param_pepoch].val[0];

  gw.phi_g       = gra*M_PI/12.0;
  gw.theta_g     = M_PI/2.0-(gdec*M_PI/180.0);
  gw.omega_g     = omega;
  gw.phase_g     = 0;
  gw.aplus_g     = aplus;
  gw.across_g    = across;
  gw.aplus_im_g  = 0.0;
  gw.across_im_g = 0.0;
  gw.phi_polar_g = 0;

  setupGW(&gw);
  for (p=0;p<*npsr;p++)
    {
      ra_p   = psr[p].param[param_raj].val[0];
      dec_p  = psr[p].param[param_decj].val[0];
      setupPulsar_GWsim(ra_p,dec_p,kp);
      for (i=0;i<psr[p].nobs;i++) 
	{
	  time = (psr[p].obsn[i].sat - timeOffset)*SECDAY;
	  gwRes[p][i] = 0.0;
	  gwRes[p][i]+=calculateResidualGW(kp,&gw,time,dist[p]);
	}
      for (i=0;i<psr[p].nobs;i++)	
	psr[p].obsn[i].sat += (gwRes[p][i]/(long double)SECDAY);
    }
  formBatsAll(psr,*npsr);                 /* Form the barycentric arrival times */
  formResiduals(psr,*npsr,0);           /* Form the residuals                 */
  doFit(psr,*npsr,0);
  formBatsAll(psr,*npsr);                 /* Form the barycentric arrival times */
  formResiduals(psr,*npsr,0);           /* Form the residuals                 */
  if (toas==1)
    {
      char fname[1000];
      for (p=0;p<*npsr;p++)
	{
	  sprintf(fname,"%s.gw.tim",psr[p].name);
	  writeTim(fname,psr+p,"tempo2");
	}
    }
  if (plotIt==1)
    doPlot(psr,*npsr,gw,gwRes,timeOffset,tspan);
  return 0;
}

void doPlot(pulsar *psr,int npsr,gwSrc gw,long double **gwRes,long double timeOffset,long double tspan)
{
  int plot=1;
  int pulsar=0;
  int endit=0;
  float mx,my;
  char key;

  cpgbeg(0,"/xs",1,1);

  cpgask(0);
  do {
  if (plot==2 || plot==3 || plot==4 || plot==5)
    plotResiduals(psr,gwRes,pulsar,timeOffset,plot-1);
  else if (plot==6)
    plotSpectrum(psr,pulsar,timeOffset);
  else if (plot==1)
    plotPosn(psr,npsr,gw);

  cpgcurs(&mx,&my,&key);
  if (key=='q') endit=1;
  else if (key=='1') plot=1;
  else if (key=='6') plot=6;
  else if (key=='4') plot=4;
  else if (key=='5') plot=5;
  else if (key=='2')
    plot=2;
  else if (key=='3')
    plot=3;
  else if (key=='p') // Select pulsar
    {
      printf("Enter pulsar number (from 0 to %d)\n",npsr-1);
      scanf("%d",&pulsar);
    }

  } while (endit==0);
  cpgend();
}

void plotPosn(pulsar *psr,int npsr,gwSrc gw)
{
  int i;
  double px[MAX_PSR],py[MAX_PSR];
  double rad2deg = 180.0/M_PI;
  float fx[MAX_PSR],fy[MAX_PSR];

  draw_grid(-180,180,-90,90,60,30,1);
  cpglab("","","Gravitational wave source and pulsar positions");
  // Plot the GW source positions
  convertXY_celestial((double)(gw.phi_g*rad2deg)-180,
		      (double)(gw.theta_g*rad2deg)-90,&px[0],&py[0]);
  fx[0] = (float)px[0];
  fy[0] = (float)py[0];
  cpgpt(1,fx,fy,8); 

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

void plotSpectrum(pulsar *psr,int p,long double timeOffset)
{
  float maxx,minx,maxy,miny;
  float fx[MAX_OBSN],fy[MAX_OBSN];
  double px[MAX_OBSN],py[MAX_OBSN];
  double ox[MAX_OBSN],oy[MAX_OBSN];
  int outN;
  double pwr;
  int i;

  for (i=0;i<psr[p].nobs;i++)
    {
      px[i] = (double)(psr[p].obsn[i].sat - timeOffset);
      py[i] = (double)psr[p].obsn[i].residual;
    }
  TKremovePoly_d(px,py,psr[p].nobs,1);  
  TKlomb_d(px,py,psr[p].nobs,1,1,ox,oy,&outN,&pwr);
  for (i=0;i<outN;i++)
    {
      fx[i] = (float)ox[i];
      fy[i] = (float)oy[i];
   }
  minx = TKfindMin_f(fx,outN);
  maxx = TKfindMax_f(fx,outN);
  miny = TKfindMin_f(fy,outN);
  maxy = TKfindMax_f(fy,outN);
  cpgenv(minx,maxx,miny,maxy,0,1);
  cpglab("Frequency (d\\u-1\\d)","DFT","");
  cpgpt(outN,fx,fy,9);
  cpgline(outN,fx,fy);
}

void plotResiduals(pulsar *psr,long double **gwRes,int p,long double timeOffset,int plotType)
{
  float px[MAX_OBSN],py[MAX_OBSN];
  float minx,maxx,miny,maxy;
  float yerr1[MAX_OBSN],yerr2[MAX_OBSN];
  int i;
  printf("MAXOBS = %d %d\n",MAX_OBSN,psr[p].nobs);
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
    cpglab("Day","GW residual (s)","Induced residuals due to GW source");
  else if (plotType==2)
    cpglab("Day","GW residual (s)","Induced residuals due to GW source after quadratic removed");
  else if (plotType==3)
    cpglab("Day","GW residual (s)","Pre-fit timing residuals");
  else if (plotType==4)
    cpglab("Day","GW residual (s)","Post-fit timing residuals");
    
  cpgpt(psr[p].nobs,px,py,9);
  if (plotType==1 || plotType==2)
    {
      TKsort_2f(px,py,psr[p].nobs);
      cpgline(psr[p].nobs,px,py);
    }
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
