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

/* Plugin for TEMPO2 that plots the pulsar's position on a P-Pdot diagram */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include <cpgplot.h>

using namespace std; 

void callFit(pulsar *psr,int npsr);
void plot_ppdot(pulsar *psr);

void help() /* Display help */
{
  printf("c         displays catalogued ephemeris for pulsar\n");
  printf("d         set declination range to select pulsars within given declination range\n");
  printf("g         set graphics device\n");
  printf("h         this help\n");
  printf("l         list pulsars in selected region\n");
  printf("q         quit\n");
  printf("u         unzoom\n");
  printf("z         zoom using mouse\n");
}

/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  char gr[100]="/xs";
  double unitFlag=1.0;  /* plot in seconds */
  double p0,p1,age;
  int i;

  *npsr = 1;  /* This graphical interface will only show results for one pulsar */

  printf("Graphical Interface: basic\n");
  printf("Author:              George Hobbs (May 2004)\n");
  printf("Version:             1.0\n");
  printf(" --- type 'H' for help information\n");
  /* Obtain the .par and the .tim file from the command line */

  if (argc==4) /* Only provided .tim name */
    {
      strcpy(timFile[0],argv[3]);
      strcpy(parFile[0],argv[3]);
      parFile[0][strlen(parFile[0])-3] = '\0';
      strcat(parFile[0],"par");
    }
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-term")==0)
	strcpy(gr,argv[++i]);
      else if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[0],argv[i+1]); 
	  strcpy(timFile[0],argv[i+2]);
	}
    }

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */

  preProcess(psr,*npsr,argc,argv);

  callFit(psr,*npsr);             /* Do all the fitting routines */


  printf("\n");
  p0 = (double)(1.0/psr[0].param[param_f].val[0]);
  p1 = (double)(-1.0/(psr[0].param[param_f].val[0]*psr[0].param[param_f].val[0])*psr[0].param[param_f].val[1]);
    
  printf("Spin-period = %g sec\n",p0);
  printf("P1 = %g\n",p1);

  age = p0/2.0/p1/60.0/60.0/24.0/365.0/1.0e6;

  printf("Characteristic age = %f Myr\n",age);
  plot_ppdot(psr);
  cpgend();
  


  return 0;
}

/* This function calls all of the fitting routines.             */
/* The function is looped twice, the first time for the pre-fit */
/* residuals and the second time for the post-fit residuals     */

void callFit(pulsar *psr,int npsr)
{
  int iteration,i;
  double globalParameter = 0.0;
  FILE *pin;

  for (iteration=0;iteration<2;iteration++)
    {
      /* Clock corrections */
      toa2utc(psr,npsr);
      // utc2tai(psr,npsr);
      //      tai2tt(psr,npsr);
      tai2ut1(psr,npsr);
      tt2tb(psr,npsr);
      /* Ephemeris routines */
      vectorPulsar(psr,npsr);
      readEphemeris(psr,npsr,0);
      get_obsCoord(psr,npsr);
      tt2tb(psr,npsr);
      readEphemeris(psr,npsr,0);
      /* extra delays */
      /*      shapiro_delay(psr,npsr);
	      dm_delays(psr,npsr); */
      /* Used to remove effect of DM model */
      /*      for (i=0;i<psr[0].nobs;i++)
	      psr[0].obsn[i].tdis = psr[0].param[param_dm].val/2.41e-16/psr[0].obsn[i].freqSSB/psr[0].obsn[i].freqSSB; */
      calculate_bclt(psr,npsr);
      
      /* Form barycentric arrival times */
      formBats(psr,npsr);
      
      /* Form residuals */
      formResiduals(psr,npsr,1);
      
      /* Do the fitting */
      if (iteration==0) doFit(psr,npsr,0);
      else textOutput(psr,npsr,globalParameter,0,0,0,"");
    }

}


/* psrcat -c "p0 p1" -o short -nonumber -nohead -l "exist(p0) && exist(p1)" > ppdot.lis */
void plot_ppdot(pulsar *psr)
{
  char fname[1024];
  float np[5000],npdot[5000],bs,pval,myp,period[5000],pdot[5000];
  float plotP[5000],plotPd[5000];
  float dec[5000],dec1=-100,dec2;
  int ip[5000];
  char str[100];
  char psrname[5000][1000];
  int i,count=0,nread,endit=0,npsr;
  float minx,maxx,miny,maxy;
  float mouseX,mouseY;
  char key;
  int setTerm;
  FILE *fin;
  char gr[100]="/xs";

  cpgbeg(0,gr,1,1);
  cpgask(0);

  cpgsci(1);
  cpgpap(0.0,1.0);
  minx = -3.0; maxx = 1.1;
  miny = -22.0; maxy = -8;

  sprintf(fname, "%s/plugin_data/ppdot.lis", getenv("TEMPO2"));
  fin = fopen(fname,"r");
  npsr=0;
  while (!feof(fin))
    {
      nread = fscanf(fin,"%f %f %s %f",&period[npsr],&pdot[npsr],psrname[npsr],&dec[npsr]);
      if (nread==4) npsr++;
    }
  fclose(fin);
  do
    {
      if (setTerm==1) setTerm=2;
      else setTerm=0;

      cpgenv(minx,maxx,miny,maxy,10,10);
      cpglab("Period (s)","log\\d10\\u(Period derivative)","" );
      count=0;
      for (i=0;i<npsr;i++)
	{
	  if (period[i]>0 && pdot[i]>0) 
	    {
	      if (dec1<-91 || (dec1 < dec[i] && dec2 > dec[i]))
		{
		  plotP[count] = log10(period[i]);
		  plotPd[count] = log10(pdot[i]);	  
		  ip[count] = i;
		  count++;
		}
	    }
	}

      cpgpt(count,plotP,plotPd,4);
      /* Draw death line */
      cpgsch(2);
      cpgslw(1);
      cpgsci(14);  
      for (i=0;i<1000;i++)
	{
	  np[i] = log10(0.0001+5*log10(i+1.0));
	  npdot[i] = log10(pow(pow(10,11.14285)/3.2e19,2)*pow(pow(10,np[i]),19.0/7.0));
	}
      cpgsls(2);  
      cpgline(i-1,np,npdot);
      cpgsls(4);
      /* Draw surface magnetic field lines */
      for (bs=8;bs<16;bs++)
	{
	  i=0;
	  for (pval=-3;pval<=1.1;pval+=0.1)
	    {
	      np[i] = pval;
	      myp = pow(10,pval);
	      npdot[i] = log10(pow(pow(10,bs)/3.2e19,2)/myp);
	      i++;
	    }
	  cpgline(i,np,npdot);
	  if (bs > 9 && bs < 16)
	    {
	      cpgsch(0.8);
	      sprintf(str,"10\\u%.0f\\dG",bs);
	      cpgtext(1.15,npdot[i-2]-0.3,str);
	      cpgsch(2);
	    }
	}
      cpgsls(1);
      cpgsci(2);
      if (psr[0].param[param_f].paramSet[0]==1)
	{
	  np[0] = (float)log10(1.0/(double)psr[0].param[param_f].val[0]);
	  np[1] = (float)log10(1.0/(double)psr[0].param[param_f].val[0]);
	  npdot[0] = -22;
	  npdot[1] = -8;
	  cpgline(2,np,npdot);
	}
      cpgsch(1);
      if (psr[0].param[param_f].paramSet[1]==1 && psr[0].param[param_f].paramSet[0]==1)
	{
	  npdot[0] = (float)(double)log10(-1.0/(psr[0].param[param_f].val[0]*psr[0].param[param_f].val[0])*psr[0].param[param_f].val[1]);
	  npdot[1] = (float)(double)log10(-1.0/(psr[0].param[param_f].val[0]*psr[0].param[param_f].val[0])*psr[0].param[param_f].val[1]);
	  np[0] = -3;
	  np[1] = 1.1;
	  cpgline(2,np,npdot);
	  np[0] = (float)(double)log10(1.0/psr[0].param[param_f].val[0]);
	  cpgpt(1,np,npdot,8);
	}
      cpgsci(1);
      cpgsch(1);
      cpgcurs(&mouseX,&mouseY,&key);
      if (key=='q')
	endit=1;
      else if (key=='g')
	{
	  cpgend();
	  cpgbeg(0,"?",1,1);
	  cpgpap(0.0,1.0);
	  cpgask(0);
	  setTerm=1;
	}
      else if (key=='h')
	help();
      else if (key=='l') /* List points in region */
	{
	  printf("\n\n");
	  for (i=0;i<count;i++)
	    {
	      if (plotP[i] > minx && plotP[i] < maxx && plotPd[i]>miny && plotPd[i]<maxy)
		printf("%s\t%f\t%g\n",psrname[ip[i]],period[ip[i]],pdot[ip[i]]);
	    }
	  printf("\n\n");
	}
      else if (key=='d') /* Set declination range */
	{
	  printf("Enter declination limit from X to Y ");
	  scanf("%f %f",&dec1,&dec2);
	}
      else if (key=='u') /* unzoom */
	{
	  minx = -3.0; maxx = 1.1;
	  miny = -22.0; maxy = -8;
	}
      else if (key=='z') /* Zoom */
	{
	  float mouseX2,mouseY2;
	  cpgband(2,0,mouseX,mouseY,&mouseX2,&mouseY2,&key);
	  minx = mouseX; maxx = mouseX2;
	  miny = mouseY; maxy = mouseY2;
	}
      else if (key=='c') /* Catalogue */
	{
	  char str[1000];
	  sprintf(str,"psrcat -all -e2 %s",psr[0].name);
	  system(str);
	}
      if (setTerm==2)
	{
	  cpgend();
	  cpgbeg(0,"/xs",1,1);
	  cpgpap(0.0,1.0);
	  cpgask(0);
	}
    } while (endit==0);

}
char * plugVersionCheck = TEMPO2_h_VER;
