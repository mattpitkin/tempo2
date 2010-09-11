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

/* Template for a tempo2 plugin */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "TKfit.h"
#include "T2toolkit.h"
#include <cpgplot.h>

using namespace std;


void help() /* Display help */
{
  /* This function should contain usage information about the plugin which should (in general) be accessed */
  /* by the user pressing 'h'                                                                              */
}
double mjd2year(double mjd);
void slaClyd ( int iy, int im, int id, int *ny, int *nd, int *jstat );
void slaCalyd ( int iy, int im, int id, int *ny, int *nd, int *j );


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  char overlay[100]="NULL";
  int i;
  double globalParameter;
  double maxy,miny;
  double wmean,wmeanE[MAX_OBSN];
  int setmaxy=0,setminy=0;
  int overlayInvert=1;
  int nstep=14;
  int simplePlot=0;

  *npsr = 0;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: name\n");
  printf("Author:              author\n");
  printf("Version:             version number\n");
  printf(" --- type 'h' for help information\n");


  printf("SHOULD BE USING ORTHONORMAL POLYNOMIALS OR BILL'S TECHNIQUE OF INDIVIDUAL POINTS?? \n");

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
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[*npsr],argv[++i]); 
	  strcpy(timFile[*npsr],argv[++i]);
	  (*npsr)++;
	}
       else if (strcmp(argv[i],"-maxy")==0)
	{
	  sscanf(argv[++i],"%lf",&maxy);
	  setmaxy=1;
	}
       else if (strcmp(argv[i],"-miny")==0)
	{
	  sscanf(argv[++i],"%lf",&miny);
	  setminy=1;
	}
      else if (strcmp(argv[i],"-invert")==0)
	overlayInvert=-1;
      else if (strcmp(argv[i],"-overlay")==0)
	strcpy(overlay,argv[++i]);
      else if (strcmp(argv[i],"-simple")==0)
	simplePlot=1;
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

  // Write new par files
  {
    char ttName[100];
    int p;
    for (p=0;p<*npsr;p++)
      {
	sprintf(ttName,"%s_clock.par",psr[p].name);
	textOutput(psr+p,1,globalParameter,0,0,1,ttName);
      }
  }


  // Show clock function
  double sx[MAX_OBSN],sy[MAX_OBSN],sy2[MAX_OBSN],sye[MAX_OBSN];
  double xval;
  double earliestTime,latestTime;
  float fx[MAX_OBSN],fy[MAX_OBSN],fye1[MAX_OBSN],fye2[MAX_OBSN];
  float ex[MAX_OBSN],ey[MAX_OBSN],ey1[MAX_OBSN],ey2[MAX_OBSN],ey0[MAX_OBSN];
  float ey2_1[MAX_OBSN],ey2_2[MAX_OBSN];
  float px[MAX_OBSN];
  int npt,j,p,k;
  int nptsClk=psr[0].nWhite*2;
  float cvX[nptsClk],cvY[nptsClk],cvY1[nptsClk],cvY2[nptsClk];
  FILE *fout_clkpts;
  FILE *fout_clkcurve;

  fout_clkpts = fopen("clock_pts.dat","w");
  fout_clkcurve = fopen("clock_curve.dat","w");

  earliestTime = (double)psr[0].obsn[0].sat;
  latestTime = (double)psr[0].obsn[psr[0].nobs-1].sat;
  for (p=1;p<*npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)
	{
	  if ((double)psr[p].obsn[i].sat < earliestTime) earliestTime = (double)psr[p].obsn[i].sat;
	  if ((double)psr[p].obsn[i].sat > latestTime) latestTime = (double)psr[p].obsn[i].sat;
	}
    }
  npt = (int)((latestTime-earliestTime)/14.0);
  printf("Earliest and latest times = %g %g\n",(double)earliestTime,(double)latestTime);
  
  for (i=0;i<npt;i++)
    {
      if (psr[0].param[param_wave_om].paramSet[0]==1)
	{
	  //	  printf("In here withasdw %g\n",(double)psr[0].param[param_waveepoch].val[0]);
	  sx[i] = earliestTime+(i*14.0)-(double)psr[0].param[param_waveepoch].val[0];
	}
      else
	sx[i] = earliestTime+(i*14.0);

      sy[i] = 0.0;
      xval = sx[i];
      sye[i] = 0.0;
      if (psr[0].param[param_wave_om].paramSet[0]==1)
	{
	  for (j=0;j<psr[0].nWhite;j++)
	    {
	      sy[i]+=(psr[0].wave_cos[j]*cos((j+1)*psr[0].param[param_wave_om].val[0]*xval)+psr[0].wave_sine[j]*sin((j+1)*psr[0].param[param_wave_om].val[0]*xval));
	      sye[i] += pow(psr[0].wave_cos_err[j]*cos((j+1)*psr[0].param[param_wave_om].val[0]*xval),2);
	      sye[i] += pow(psr[0].wave_sine_err[j]*sin((j+1)*psr[0].param[param_wave_om].val[0]*xval),2);
	    }
	  sy2[i] = sy[i];
	  sye[i] = sqrt(sye[i]);
	}
      else
	{
	  long double speriod,tt,dt;
	  long double t1=0.0L;
	  speriod = (long double)(psr[0].ifuncT[1]-psr[0].ifuncT[0]); 
	  xval = sx[i];
	  for (k=0;k<psr[0].ifuncN;k++)
	    {
	      dt = xval-(long double)psr[0].ifuncT[k];
	      tt = M_PI/speriod*dt;
	      t1+=(long double)psr[0].ifuncV[k]*sinl(tt)/tt;
	    }
	  sy[i] = t1;
	  sy2[i] = sy[i];
	  //	  printf("Here with %g %g %Lg %g\n",sx[i],sy[i],speriod,xval);
	}
    }
  TKremovePoly_d(sx,sy2,npt,3);
  for (i=0;i<npt;i++)
    {
      if (psr[0].param[param_wave_om].paramSet[0]==1)
	px[i] = (float)mjd2year(sx[i]+(double)psr[0].param[param_waveepoch].val[0]);
      else
	px[i] = (float)mjd2year(sx[i]);

      fx[i] = (float)sx[i];
      fy[i] = (float)sy2[i];
      fye1[i] = (float)(sy2[i] - sye[i]);
      fye2[i] = (float)(sy2[i] + sye[i]);
      //      printf("fye = %g %g\n",fye1[i],fye2[i]);
      fprintf(fout_clkcurve,"%g %g %g %g\n",px[i],sx[i],sy[i],sy2[i]);
    }

  // Calculate variances
  int ne=0;

  if (psr[0].param[param_wave_om].paramSet[0]==1)
    {
      double var;
      int nvar;
      double time;
      double value[*npsr][MAX_OBSN];
      double evalue[*npsr][MAX_OBSN];
      
      double wmean,mean;
      int nvalue[*npsr];
      double wi,s1,s2;
      long double ls1,ls2;
      int np;
      
      for (i=0;i<npt-1;i+=nstep)
	{
	  var = 0.0;
	  nvar = 0;
	  
	  for (p=0;p<*npsr;p++)
	    {
	      nvalue[p]=0;
	      
	      for (j=0;j<psr[p].nobs;j++)
		{
		  if ((double)(psr[p].obsn[j].sat-psr[0].param[param_waveepoch].val[0]) > fx[i] &&
		      (double)(psr[p].obsn[j].sat-psr[0].param[param_waveepoch].val[0]) < fx[i+nstep])
		    {
		      value[p][nvalue[p]]=(double)psr[p].obsn[j].residual;
		      evalue[p][nvalue[p]]=(double)psr[p].obsn[j].toaErr*1.0e-6;
		      nvalue[p]++;
		      //		  var += pow(psr[p].obsn[j].residual,2);
		      //		  nvar++;
		    }
		}
	    }
	  // Remove means
	  for (p=0;p<*npsr;p++)
	    {
	      mean = TKmean_d(value[p],nvalue[p]);
	      for (j=0;j<nvalue[p];j++)
		value[p][j]-=mean;
	    }
	  // Find weighted mean
	  np=0;
	  wmean=0;
	  s1=0;
	  s2=0;
	  for (p=0;p<*npsr;p++)
	    {
	      for (j=0;j<nvalue[p];j++)
		{
		  wi = 1.0/pow(evalue[p][j],2);
		  s1 += wi*value[p][j];
		  s2 += wi;
		  np++;
		}
	    }
	  if (np>0)
	    {
	      wmean = s1/s2;
	      printf("weighted mean = %g\n",wmean);
	    }
	  // Find variance
	  ls1=0.0; ls2=0.0;
	  nvar=0;
	  for (p=0;p<*npsr;p++)
	    {
	      for (j=0;j<nvalue[p];j++)
		{
		  wi = 1.0/pow(evalue[p][j],2);
		  
		  ls1 += (wi*pow(value[p][j]-wmean,2));
		  ls2 += (wi);

		  nvar++;
		}
	    }
		  
	  if (nvar > 2)
	    {
	      var = ls1/(ls2)/(((double)nvar-1.0)/(double)nvar); 
	      //	  var = ls1/(double)nvar;
	      printf("var = %g %Lg %Lg\n",var,ls1,ls2);
	      ex[ne]=0.5*(px[i]+px[i+nstep]);
	      ey[ne] = sqrt(var);
	      ey0[ne] = fy[i+(int)(nstep/2.0+0.5)];
	      ey1[ne] = fy[i+(int)(nstep/2.0+0.5)]-ey[ne];
	      ey2[ne] = fy[i+(int)(nstep/2.0+0.5)]+ey[ne];

	      // Calculate error in weighted mean
	      int c=0;
	      ls1=0.0;
	      ls2=0.0;
	      for (p=0;p<*npsr;p++)
		{
		  for (j=0;j<nvalue[p];j++)
		    {
		      wi = 1.0/pow(evalue[p][j],2);
		      ls1 += wi;
		      c++;
		    }
		}
	      wmeanE[ne] = sqrt(1.0/ls1);
	      printf("This point %d %g %d %d %d %d\n",ne,wmeanE[ne],c,nvalue[0],nvalue[1],nvalue[2]);
	      ey2_1[ne] = fy[i+(int)(nstep/2.0+0.5)]-wmeanE[ne];
	      ey2_2[ne] = fy[i+(int)(nstep/2.0+0.5)]+wmeanE[ne];
	      ne++;
	    }
	}
      //
      // Determine errors using covariance matrix
      //
      int nfit=psr[0].nWhite*2;
      double covp[nfit][nfit],t;
      double matM[nptsClk][nfit];
      double matMT[nfit][nptsClk];

      
      // Must set properly knowing whether the covariances are
      printf("Covariance matrix, nptsClk = %d\n",nptsClk);
      for (i=0;i<nfit;i++)
	{
	  for (j=0;j<nfit;j++)
	    {
	      covp[i][j] = psr[0].covar[i][j];
	      printf("%g ",covp[i][j]);	      
	    }
	  printf("\n");
	}
      // Set matrix of cosine and sine terms
      for (i=0;i<nptsClk;i++)
	{
	  t = earliestTime+i*(latestTime-earliestTime)/(double)nptsClk-psr[0].param[param_waveepoch].val[0]+(latestTime-earliestTime)/(double)nptsClk*0.5;
	  cvX[i] = (float)mjd2year(t+psr[0].param[param_waveepoch].val[0]);
	  cvY[i] = 0.0;
	  for (j=0;j<psr[0].nWhite;j++)
	    {
	      cvY[i]+=(psr[0].wave_cos[j]*cos((j+1)*psr[0].param[param_wave_om].val[0]*t)+psr[0].wave_sine[j]*sin((j+1)*psr[0].param[param_wave_om].val[0]*t));
	    }
	  for (j=0;j<nfit/2;j++)
	    {
	      matM[i][2*j] = cos((j+1)*psr[0].param[param_wave_om].val[0]*t);
	      matM[i][2*j+1] = sin((j+1)*psr[0].param[param_wave_om].val[0]*t);
	    }
	}
      TKremovePoly_f(cvX,cvY,nptsClk,3);
      printf("Matrix of sine and cosine terms\n");
      for (i=0;i<nptsClk;i++)
	{
	  for (j=0;j<nfit;j++)
	    printf("%g ",matM[i][j]);
	  printf("\n");
	}
      // Get transpose of the matrix
      for (i=0;i<nptsClk;i++)
	{
	  for (j=0;j<nfit;j++)
	    matMT[j][i] = matM[i][j];
	}
      printf("Transpose of matrix\n");
      for (j=0;j<nfit;j++)
	{
	  for (i=0;i<nptsClk;i++)
	    printf("%g ",matMT[j][i]);
	  printf("\n");
	}
      // Multiply the matrices together
      // covP.M^T
      double add;
      double m1[nfit][nptsClk];
      //double m1[nptsClk][nfit];
      double m2[nptsClk][nptsClk];

      
	for (i=0;i<nptsClk;i++)
	{
	  for (j=0;j<nfit;j++)
	    {
	      add=0.0;
	      for (k=0;k<nfit;k++)
		{
		  if (i==0 && j==0)
		    {
		      printf("%g %g ",covp[k][j],matMT[i][k]);
		      printf("\n");
		    }
		  add+=covp[j][k]*matMT[k][i];
		}
	      m1[j][i]=add;
	    }
	}
      printf("First multiplication\n");
      for (j=0;j<nfit;j++)
	{
	  for (i=0;i<nptsClk;i++)
	    printf("%g ",m1[j][i]);
	  printf("\n");
	}

      for (i=0;i<nptsClk;i++)
	{
	  for (j=0;j<nptsClk;j++)
	    {
	      add=0.0;
	      for (k=0;k<nfit;k++)
		add+=matM[i][k]*m1[k][j];
	      m2[i][j] = add;
	    }	  
	}
      printf("Second multiplication\n");
      for (i=0;i<nptsClk;i++)
	{
	  for (j=0;j<nptsClk;j++)
	    printf("%g ",m2[i][j]);
	  printf("\n");
	}
      


      for (i=0;i<nptsClk;i++)
	{
	  printf("diag = %g, error = %g, y = %g\n",m2[i][i],sqrt(m2[i][i]),cvY[i]);
	  cvY1[i] = cvY[i] - sqrt(m2[i][i]);
	  cvY2[i] = cvY[i] + sqrt(m2[i][i]);
	  fprintf(fout_clkpts,"%g %g %g\n",cvX[i],cvY[i],sqrt(m2[i][i]));
	}
    }
  else
    {
      for (k=0;k<psr[0].ifuncN;k++)
	{
	  ex[k] = (float)mjd2year(psr[0].ifuncT[k]);
	  ey0[k] = (float)psr[0].ifuncV[k];
	}
      ne = psr[0].ifuncN;
      TKremovePoly_f(ex,ey0,ne,3);
      for (k=0;k<psr[0].ifuncN;k++)
	{
	  ey1[k] = (float)(ey0[k]-psr[0].ifuncE[k]);
	  ey2[k] = (float)(ey0[k]+psr[0].ifuncE[k]);
	  printf("Error has size %g\n",psr[0].ifuncE[k]);
	}

    }
  fclose(fout_clkpts);
  fclose(fout_clkcurve);

  cpgbeg(0,"?",1,1);
  cpgsfs(2);
  cpgslw(3);
  cpgsch(1.4);
  if (setmaxy==0)
    maxy = TKfindMax_f(fy,npt);
  if (setminy==0)
  miny = TKfindMin_f(fy,npt);
  if (simplePlot==0)
    {
      cpgsvp(0.1,0.95,0.50,0.95);
      cpgswin(px[0],px[npt-1],0,1);
      cpgbox("BCST",0,0,"BCS",0,0);
      // Plot data spans
      {
	float dspanX[MAX_OBSN],dspanY[MAX_OBSN];
	for (p=0;p<*npsr;p++)
	  {
	    for (i=0;i<psr[p].nobs;i++)
	      {
		dspanX[i] = mjd2year((double)psr[p].obsn[i].sat);
		dspanY[i] = 0.9-(0.8)*(double)p/(double)(*npsr);
	      }
	    cpgsci((p%3)+1); 
	    cpgsch(0.6);
	    cpgtext(px[0]-(px[npt-1]-px[0])*0.10,dspanY[0]-0.8/(*npsr)*0.25,psr[p].name);
	    cpgpt(psr[p].nobs,dspanX,dspanY,16);
	    cpgsch(1.4);
	  }
	cpgsci(6);
	if (psr[0].param[param_wave_om].paramSet[0]==1)
	  {
	    for (i=0;i<npt-1;i+=nstep)
	      {
		dspanX[0] = mjd2year(fx[i]+psr[0].param[param_waveepoch].val[0]);
		dspanX[1] = mjd2year(fx[i+nstep]+psr[0].param[param_waveepoch].val[0]);
		dspanY[0] = dspanY[1] = 0.05;
		cpgsah(2,45,0.3);
		cpgsch(0.5);cpgarro(dspanX[0],0.05,dspanX[1],0.05); cpgsch(1.4);
		cpgsch(0.5);cpgarro(dspanX[1],0.05,dspanX[0],0.05); cpgsch(1.4);
	      }
	  }
	cpgsci(1);
      }
      cpgsci(1);


      cpgsvp(0.1,0.95,0.13,0.50);
  cpgswin(px[0],px[npt-1],miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny));
  cpgsch(1);
  cpgbox("BNCTS",0,0,"BCNTS",0,0);
  cpglab("Year","Clock difference (sec)","");
    }
  else
    {
      cpgenv(px[0],px[npt-1],miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny),0,1);
  cpglab("Year","Clock difference (sec)","");
    }

  cpgsch(1.4);
  cpgsci(2);
  cpgsls(4);cpgline(npt,px,fy);cpgsls(1);
  cpgsci(1);
  //  cpgsls(4);
  //  cpgline(npt,px,fye1);
  //  cpgline(npt,px,fye2);
  //  cpgsls(1);
  //  cpgerry(ne,ex,ey1,ey2,1);
  //  cpgpt(ne,ex,ey0,4);
  //  cpgslw(2);  cpgsci(3); cpgerry(ne,ex,ey2_1,ey2_2,1);cpgsci(1); cpgslw(1);
  
  // Covariance
  cpgsci(7); cpgpt(nptsClk,cvX,cvY,4); cpgerry(nptsClk,cvX,cvY1,cvY2,1); cpgsci(1);
  
  if (strcmp(overlay,"NULL")!=0)
    {
      FILE *fin;
      FILE *fout;
      float frx[MAX_OBSN],fry[MAX_OBSN];
      double rx0[MAX_OBSN];
      double rx[MAX_OBSN],ry[MAX_OBSN],ry2[MAX_OBSN];

      int nr=0;

      fout = fopen("expected.dat","w");

      if (!(fin = fopen(overlay,"r")))
	{
	  printf("Unable to open file: >%s<\n",overlay);
	  exit(1);
	}
      while (!feof(fin))
	{
	  if (fscanf(fin,"%lf %lf",&rx[nr],&ry[nr]))
	    {
	      if (rx[nr] > earliestTime && rx[nr] < latestTime)
		{
		  rx0[nr] = (double)(rx[nr]-psr[0].param[param_pepoch].val[0]);
		  rx[nr] = mjd2year(rx[nr]);
		  ry[nr] *= overlayInvert;
		  nr++;
		}
	    }
	}
      fclose(fin);
      printf("Have read %d lines from the clock file\n",nr);
      // Remove quadratic
      TKremovePoly_d(rx,ry,nr,3);
      for (i=0;i<nr;i++)
	{
	  frx[i]=(float)rx[i];
	  fry[i]=(float)ry[i]; // Check not including a MINUS SIGN!!!
	  //	  printf("Expected: %g %g %g\n",frx[i],rx0[i],fry[i]);
	  fprintf(fout,"%g %g\n",rx[i],ry[i]);
	}
      fclose(fout);
      cpgsci(3);
      //      cpgsls(2); 
      cpgslw(4); cpgline(nr,frx,fry);  cpgslw(2);
      //      cpgsls(1);
      cpgsci(1);
      {
	double dist=0;
	// Now only obtain the overlay file at the points from the covariance matrix
	for (i=0;i<nptsClk;i++)
	  {
	    dist = (fabs)(rx[0] - cvX[i]); // Simply take closest point -- SHOULD DO AN INTERPOLATION
	    ry2[i] = ry[0];
	    for (j=0;j<nr;j++)
	      {
		if ((fabs)(rx[j] - cvX[i]) < dist)
		  {
		    dist = fabs(rx[j]-cvX[i]);
		    ry2[i] = ry[j];
		  }
	      }
	  }
	for (i=0;i<nptsClk;i++)
	  {
	    frx[i] = (float)cvX[i];
	    fry[i] = (float)ry2[i];
	  }
	TKremovePoly_f(frx,fry,nptsClk,3);	
	//	cpgpt(nptsClk,frx,fry,16);
	//	cpgline(nptsClk,frx,fry);
      }
    }
  fx[0] = px[0];
  fx[1] = px[npt-1];
  fy[0] = 0.0;
  fy[1] = 0.0;
  cpgsls(4); cpgline(2,fx,fy); cpgsls(1);
  cpgend();
  return 0;
}

double mjd2year(double mjd)
{
      double jd,fjd,day;
      int ijd,b,c,d,e,g,month,year;
      int retYr,retDay,stat;

      jd = mjd + 2400000.5;
      ijd = (int)(jd+0.5);
      fjd = (jd+0.5)-ijd;
      if (ijd > 2299160)
	{
	  int a;
	  a = (int)((ijd-1867216.25)/36524.25);
	  b = ijd + 1 + a - (int)(a/4.0);
	}
      else
	b = ijd;

      c = b + 1524;
      d = (int)((c - 122.1)/365.25);
      e = (int)(365.25*d);
      g = (int)((c-e)/30.6001);
      day = c-e+fjd-(int)(30.6001*g);
      if (g<13.5)
	month = g-1;
      else
	month = g-13;
      if (month>2.5)
	year = d-4716;
      else
	year = d-4715;
      slaCalyd(year, month, (int)day, &retYr, &retDay, &stat);


      return  retYr+(retDay+(day-(int)day))/365.25;

}
void slaCalyd ( int iy, int im, int id, int *ny, int *nd, int *j )
/*
**  - - - - - - - - -
**   s l a C a l y d
**  - - - - - - - - -
**
**  Gregorian calendar date to year and day in year (in a Julian
**  calendar aligned to the 20th/21st century Gregorian calendar).
**
**  (Includes century default feature:  use slaClyd for years
**   before 100AD.)
**
**  Given:
**     iy,im,id   int    year, month, day in Gregorian calendar
**                       (year may optionally omit the century)
**  Returned:
**     *ny        int    year (re-aligned Julian calendar)
**     *nd        int    day in year (1 = January 1st)
**     *j         int    status:
**                         0 = OK
**                         1 = bad year (before -4711)
**                         2 = bad month
**                         3 = bad day (but conversion performed)
**
**  Notes:
**
**  1  This routine exists to support the low-precision routines
**     slaEarth, slaMoon and slaEcor.
**
**  2  Between 1900 March 1 and 2100 February 28 it returns answers
**     which are consistent with the ordinary Gregorian calendar.
**     Outside this range there will be a discrepancy which increases
**     by one day for every non-leap century year.
**
**  3  Years in the range 50-99 are interpreted as 1950-1999, and
**     years in the range 00-49 are interpreted as 2000-2049.
**
**  Called:  slaClyd
**
**  Last revision:   22 September 1995
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   int i;

/* Default century if appropriate */
   if ( ( iy >= 0 ) && ( iy <= 49 ) )
      i = iy + 2000;
   else if ( ( iy >= 50 ) && ( iy <= 99 ) )
      i = iy + 1900;
   else
      i = iy;

/* Perform the conversion */
   slaClyd ( i, im, id, ny, nd, j );
}

void slaClyd ( int iy, int im, int id, int *ny, int *nd, int *jstat )
/*
**
**  Returned:
**     ny          int    year (re-aligned Julian calendar)
**     nd          int    day in year (1 = January 1st)
**     jstat       int    status:
**                          0 = OK
**                          1 = bad year (before -4711)
**                          2 = bad month
**                          3 = bad day (but conversion performed)
**
**  Notes:
**
**  1  This routine exists to support the low-precision routines
**     slaEarth, slaMoon and slaEcor.
**
**  2  Between 1900 March 1 and 2100 February 28 it returns answers
**     which are consistent with the ordinary Gregorian calendar.
**     Outside this range there will be a discrepancy which increases
**     by one day for every non-leap century year.
**
**  3  The essence of the algorithm is first to express the Gregorian
**     date as a Julian Day Number and then to convert this back to
**     a Julian calendar date, with day-in-year instead of month and
**     day.  See 12.92-1 and 12.95-1 in the reference.
**
**  Reference:  Explanatory Supplement to the Astronomical Almanac,
**              ed P.K.Seidelmann, University Science Books (1992),
**              p604-606.
**
**  Last revision:   26 November 1994
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   long i, j, k, l, n, iyL, imL;

/* Month lengths in days */
   static int mtab[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };



/* Validate year */
   if ( iy < -4711 ) { *jstat = 1; return; }

/* Validate month */
   if ( ( im < 1 ) || ( im > 12 ) ) { *jstat = 2; return; }

/* Allow for (Gregorian) leap year */
   mtab[1] = ( ( ( iy % 4 ) == 0 ) &&
             ( ( ( iy % 100 ) != 0 ) || ( ( iy % 400 ) == 0 ) ) ) ?
             29 : 28;

/* Validate day */
   *jstat = ( id < 1 || id > mtab[im-1] ) ? 3 : 0;

/* Perform the conversion */
   iyL = (long) iy;
   imL = (long) im;
   i = ( 14 - imL ) /12L;
   k = iyL - i;
   j = ( 1461L * ( k + 4800L ) ) / 4L
     + ( 367L * ( imL - 2L + 12L * i ) ) / 12L
     - ( 3L * ( ( k + 4900L ) / 100L ) ) / 4L + (long) id - 30660L;
   k = ( j - 1L ) / 1461L;
   l = j - 1461L * k;
   n = ( l - 1L ) / 365L - l / 1461L;
   j = ( ( 80L * ( l - 365L * n + 30L ) ) / 2447L ) / 11L;
   i = n + j;
   *nd = 59 + (int) ( l -365L * i + ( ( 4L - n ) / 4L ) * ( 1L - j ) );
   *ny = (int) ( 4L * k + i ) - 4716;
}
