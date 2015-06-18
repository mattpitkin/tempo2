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

bool cholmode=false;
char covarFuncFile[MAX_FILELEN];

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
  int i,j,k,p;
  double globalParameter;
  double maxy,miny;
  double wmean,wmeanE[MAX_OBSN];
  int setmaxy=0,setminy=0;
  int overlayInvert=1;
  int invertRes=0;
  int nstep=14;
  int simplePlot=0;
  int harmTaper=0;
  int covarError=0;
  int shuffle=0;
  double whiteErr[MAX_IFUNC];
  int showName=1;
  long seed = TKsetSeed();
  char grDev[1024]="?";

  strcpy(covarFuncFile,"NULL");

  *npsr = 0;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: clock\n");
  printf("Author:              G. Hobbs\n");
  printf("Version:             version number\n");
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
      else if (strcmp(argv[i],"-noname")==0)
	showName=0;
       else if (strcmp(argv[i],"-maxy")==0)
	{
	  sscanf(argv[++i],"%lf",&maxy);
	  setmaxy=1;
	}
      else if (strcmp(argv[i],"-harmTaper")==0)
	sscanf(argv[++i],"%d",&harmTaper);
      else if (strcmp(argv[i],"-covarError")==0)
	covarError=1;
       else if (strcmp(argv[i],"-miny")==0)
	{
	  sscanf(argv[++i],"%lf",&miny);
	  setminy=1;
	}
      else if (strcmp(argv[i],"-invert")==0)
	overlayInvert=-1;
      else if (strcmp(argv[i],"-invertRes")==0)
	invertRes=1;
      else if (strcmp(argv[i],"-overlay")==0)
	strcpy(overlay,argv[++i]);
      else if (strcmp(argv[i],"-simple")==0)
	simplePlot=1;
      else if (strcmp(argv[i],"-shuffle")==0)
	shuffle=1;
      else if (strcmp(argv[i],"-g")==0)
	strcpy(grDev,argv[++i]);
      else if (strcmp(argv[i],"-dcf")==0){
		cholmode=true;
		strcpy(covarFuncFile,argv[++i]);
      }

    }
  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);
  if (shuffle==1) // Shuffle all the ToAs by their error bars
    {
      printf("WARNING: shuffling ToAs\n");
      for (p=0;p<*npsr;p++)
	{
	  for (i=0;i<psr[p].nobs;i++)
	    {
	      psr[p].obsn[i].sat += (TKgaussDev(&seed)*psr[p].obsn[i].toaErr*1e-6)/86400.0L;
	    }
	}

    }


  if (psr[0].param[param_wave_om].paramSet[0]==1 && harmTaper==0)
    {
      printf("Must set a harmonic taper number using -harmTaper (make smaller than the number of waves being fit)\n");
      exit(1);
    }

  // If using IFUNCS get the error bars for the white component by turning off all
  // fits except for the global parameters
  if (psr[0].ifuncN > 0 && covarError==1)
    {
      FILE *fout_nofit;
      //      char tpar[MAX_STRLEN][MAX_FILELEN];
      //      char ttim[MAX_STRLEN][MAX_FILELEN];

      //      sprintf(tpar[0],"global_clock1.par");
      fout_nofit = fopen("clk_nofit.dat","w");
      for (p=0;p<*npsr;p++)
	{
	  for (i=0;i<MAX_PARAMS;i++)
	    {
	      for (j=0;j<psr[p].param[i].aSize;j++)
		{
		  if (psr[p].param[i].fitFlag[j] == 1)
		    psr[p].param[i].fitFlag[j] = 0;
		}
	    }

	  // Turn off jump fitting
	  for (i=0;i<psr[p].nJumps;i++)
	    {
	      if (psr[p].fitJump[i] == 1)
		psr[p].fitJump[i]=0;
	    }
	}
      // Now do the fit
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      doFitAll(psr,*npsr,covarFuncFile);   /* Do the fitting     */
      //      textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
      // Store error bars
      for (i=0;i<psr[0].ifuncN;i++)
	{
	  printf("Values without fitting: %g %g %g\n",psr[0].ifuncT[i],psr[0].ifuncV[i],psr[0].ifuncE[i]);
	  fprintf(fout_nofit,"Values without fitting: %g %g %g\n",psr[0].ifuncT[i],psr[0].ifuncV[i],psr[0].ifuncE[i]);
	  whiteErr[i] = psr[0].ifuncE[i];
	}
      fclose(fout_nofit);
      // Re-read the par file
      for (p=0;p<*npsr;p++)
	{
	  psr[p].nconstraints = 0;
	  psr[p].nJumps = 0;
	  psr[p].nT2efac = 0;
	  psr[p].nT2equad = 0;
	  psr[p].T2globalEfac = 1;
	  for(i=0;i<MAX_JUMPS;i++){
	    psr[p].jumpVal[i] = 0.0;
	    psr[p].jumpValErr[i] = 0.0;
	  }
	  for(i=0;i<MAX_PARAMS;i++){
	    psr[p].param[i].nLinkTo = 0;
	    psr[p].param[i].nLinkFrom = 0;
	  }
	}
      
      readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
      readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
      preProcess(psr,*npsr,argc,argv);
      // Now re-read the global par file
      //      readParfileGlobal(psr,*npsr,tpar,ttim);
      // Continue as normal
    }


  printf("Step 1\n");
  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      if (i==0) doFitAll(psr,*npsr,covarFuncFile);   /* Do the fitting     */
      else textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
    }
  printf("Step 2:\n");
  for (k=0;k<psr[0].ifuncN;k++)
    {
      psr[0].ifuncE[k] *= 1e6; // Convert to microseconds
      psr[0].ifuncV[k] *= 1e6;
    }

  // Remove last ifunc if using no interpolation
  if (psr[0].param[param_ifunc].val[0]==0)
    psr[0].ifuncN--;
    

  //  for (p=0;p<*npsr;p++)
  //    {
  //      for (i=0;i<psr[p].nobs;i++)
  //	printf("psr%d %g %g\n",p+1,(double)psr[p].obsn[i].sat,(double)psr[p].obsn[i].residual);
  //    }



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
  double taperY1[MAX_OBSN];
  double xval;
  double earliestTime,latestTime;
  float fx[MAX_OBSN],fy[MAX_OBSN],fye1[MAX_OBSN],fye2[MAX_OBSN];
  float fyTaper[MAX_OBSN];
  double taperVal;
  float ex[MAX_OBSN],ey[MAX_OBSN],ey1[MAX_OBSN],ey2[MAX_OBSN],ey0[MAX_OBSN];
  float ey2_1[MAX_OBSN],ey2_2[MAX_OBSN];
  float px[MAX_OBSN];
  int npt,nfit;
  int nptsClk=psr[0].nWhite*2;
  int nr=0;
  float cvX[nptsClk],cvY[nptsClk],cvY1[nptsClk],cvY2[nptsClk],cvTime[nptsClk];
  FILE *fout_clkpts;
  FILE *fout_clkcurve;
  FILE *fout_newclk;
  FILE *fout_covar;
  char str[128];
  float frx[MAX_OBSN],fry[MAX_OBSN];
      
  // Look at the covariance matrix of the fitted parameters
  if (psr[0].ifuncN > 0)
    {
      nfit = psr[0].ifuncN;
      for (i=0;i<nfit;i++)
	{
	  sprintf(str,"covarpt_%d",i+1);
	  fout_covar = fopen(str,"w");
	  for (j=0;j<nfit;j++)
	    {
	      printf("%g ",psr[0].covar[i][j]);
	      fprintf(fout_covar,"%d %g\n",j-i,psr[0].covar[i][j]);
	    }
	  printf("\n");
	  fclose(fout_covar);
	}
    }


  fout_clkpts = fopen("clock_pts.dat","w");
  fout_clkcurve = fopen("clock_curve.dat","w");
  fout_newclk = fopen("tai2tt_ppta2010.clk","w");
  fprintf(fout_newclk,"# TAI TT(PPTA2010)\n");
  fprintf(fout_newclk,"# Obtained using the clock plugin\n");
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
      taperY1[i] = 0.0;
      xval = sx[i];
      sye[i] = 0.0;
      if (psr[0].param[param_wave_om].paramSet[0]==1)
	{
	  for (j=0;j<psr[0].nWhite;j++)
	    {
	      sy[i]+=(psr[0].wave_cos[j]*cos((j+1)*psr[0].param[param_wave_om].val[0]*xval)+psr[0].wave_sine[j]*sin((j+1)*psr[0].param[param_wave_om].val[0]*xval));
	      sye[i] += pow(psr[0].wave_cos_err[j]*cos((j+1)*psr[0].param[param_wave_om].val[0]*xval),2);
	      sye[i] += pow(psr[0].wave_sine_err[j]*sin((j+1)*psr[0].param[param_wave_om].val[0]*xval),2);
	      taperVal = pow(1.0+pow((j+1)/(float)harmTaper,2),2);
		taperY1[i] += (psr[0].wave_cos[j]/taperVal*cos((j+1)*psr[0].param[param_wave_om].val[0]*xval)+psr[0].wave_sine[j]/taperVal*sin((j+1)*psr[0].param[param_wave_om].val[0]*xval));
	    }
	  sy2[i] = sy[i];
	  sye[i] = sqrt(sye[i]);
	}
    }
  if (psr[0].param[param_wave_om].paramSet[0]==1)
    {
      TKremovePoly_d(sx,sy2,npt,3);
      TKremovePoly_d(sx,taperY1,npt,3);
      for (i=0;i<npt;i++)
	{
	  if (psr[0].param[param_wave_om].paramSet[0]==1)
	    px[i] = (float)mjd2year(sx[i]+(double)psr[0].param[param_waveepoch].val[0]);
	  else
	    px[i] = (float)mjd2year(sx[i]);
	  
	  fx[i] = (float)sx[i];
	  fy[i] = (float)sy2[i];
	  fyTaper[i] = (float)taperY1[i];
	  fye1[i] = (float)(sy2[i] - sye[i]);
	  fye2[i] = (float)(sy2[i] + sye[i]);
	  //      printf("fye = %g %g\n",fye1[i],fye2[i]);
	  fprintf(fout_clkcurve,"%g %g %g %g %g\n",px[i],sx[i],sy[i],sy2[i],sx[i]+(double)psr[0].param[param_waveepoch].val[0]);
	  fprintf(fout_newclk,"%g %.15f\n",sx[i]+(double)psr[0].param[param_waveepoch].val[0],(double)(32.184L+taperY1[i]));
	}
    }
  else
    {
      for (i=0;i<npt;i++)
	{
	    px[i] = (float)mjd2year(sx[i]);
	}
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
	  cvTime[i] = (double)(t+psr[0].param[param_waveepoch].val[0]);
	  cvX[i] = (float)mjd2year(t+psr[0].param[param_waveepoch].val[0]);
	  cvY[i] = 0.0;
	  for (j=0;j<psr[0].nWhite;j++)
	    {
	      taperVal = pow(1.0+pow((j+1)/(float)harmTaper,2),2);
	      cvY[i]+=(psr[0].wave_cos[j]/taperVal*cos((j+1)*psr[0].param[param_wave_om].val[0]*t)+psr[0].wave_sine[j]/taperVal*sin((j+1)*psr[0].param[param_wave_om].val[0]*t));
	    }
	  for (j=0;j<nfit/2;j++)
	    {
	      taperVal = pow(1.0+pow((j+1)/(float)harmTaper,2),2);
	      matM[i][2*j] = cos((j+1)*psr[0].param[param_wave_om].val[0]*t)/taperVal;
	      matM[i][2*j+1] = sin((j+1)*psr[0].param[param_wave_om].val[0]*t)/taperVal;
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
	  fprintf(fout_clkpts,"%g %g %g %g\n",cvX[i],cvY[i],sqrt(m2[i][i]),cvTime[i]);
	}
    }
  else
    {
      double dx[psr[0].ifuncN];
      double dy[psr[0].ifuncN];
      double de[psr[0].ifuncN];
      double **cvm;
      double chisq;
      double out_p[3],out_e[3];
      cvm = (double **)malloc(sizeof(double *)*3);
      for (i=0;i<3;i++)
	cvm[i] = (double *)malloc(sizeof(double)*3);
   
      printf("In here - GEORGE\n");
      for (k=0;k<psr[0].ifuncN;k++)
	{
	  ex[k] = (float)mjd2year(psr[0].ifuncT[k]);
	  ey0[k] = (float)psr[0].ifuncV[k];
	  if (invertRes==1) ey0[k]*=-1;
	  dx[k] = (double)ex[k];
	  dy[k] = (double)ey0[k];
	  de[k] = psr[0].ifuncE[k];
	  printf("Using: %g %g %g\n",dx[k],dy[k],de[k]);
	}
      ne = psr[0].ifuncN;

      // Removing a weighted quadratic from the IFUNCS
      TKleastSquares_svd(dx,dy,de,ne,out_p,out_e,3,cvm,&chisq,TKfitPoly,1);
      for (i=0;i<ne;i++)
	{
	  printf("ifunc output: %g %g %g %g\n",ex[i],ey0[i],dy[i] - (out_p[0] + out_p[1]*dx[i] + out_p[2]*dx[i]*dx[i]),de[i]);
	  ey0[i] = dy[i] - (out_p[0] + out_p[1]*dx[i] + out_p[2]*dx[i]*dx[i]);
	}
      //TKremovePoly_f(ex,ey0,ne,3);
      for (i=0;i<3;i++)
	free(cvm[i]);
      free(cvm);
		      
      for (k=0;k<psr[0].ifuncN;k++)
	{
	  ey1[k] = (float)(ey0[k]-psr[0].ifuncE[k]);
	  ey2[k] = (float)(ey0[k]+psr[0].ifuncE[k]);
	  printf("Error has size %g value = %g\n",psr[0].ifuncE[k],ey0[k]);
	}
    }

  fclose(fout_clkpts);
  fclose(fout_clkcurve);
  fclose(fout_newclk);

  cpgbeg(0,grDev,1,1);
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
      cpgswin(px[0]-1,px[npt-1]+1,0,1);
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
	    if (showName==1)
	      cpgtext(px[0]-1-(px[npt-1]-px[0]+2)*0.10,dspanY[0]-0.8/(*npsr)*0.25,psr[p].name);
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
		//		cpgsah(2,45,0.3);
		//		cpgsch(0.5);cpgarro(dspanX[0],0.05,dspanX[1],0.05); cpgsch(1.4);
		//		cpgsch(0.5);cpgarro(dspanX[1],0.05,dspanX[0],0.05); cpgsch(1.4);
	      }
	  }
	cpgsci(1);
      }
      cpgsci(1);


      cpgsvp(0.1,0.95,0.13,0.50);
      cpgswin(px[0]-1,px[npt-1]+1,miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny));
      cpgsch(1);
      cpgbox("BNCTS",0,0,"BCNTS",0,0);
      cpglab("Year","Clock difference (\\gms)","");
    }
  else
    {
      cpgenv(px[0],px[npt-1],miny-0.1*(maxy-miny),maxy+0.1*(maxy-miny),0,1);
      cpglab("Year","Clock difference (sec)","");
    }

  cpgsch(1.4);
  //  cpgsci(2);
  //      cpgsls(4);cpgline(npt,px,fy);cpgsls(1);
  //  cpgsci(1); cpgsls(2);cpgline(npt,px,fyTaper);cpgsls(1);
  cpgsci(1);

  printf("GOT HERE\n");
  if (strcmp(overlay,"NULL")!=0)
    {
      FILE *fin;
      FILE *fout;
      double rx0[10000];
      double rx[10000],ry[10000],ry2[10000];
      char str[1024];

      fout = fopen("expected.dat","w");

      if (!(fin = fopen(overlay,"r")))
	{
	  printf("Unable to open file: >%s<\n",overlay);
	  exit(1);
	}
      while (!feof(fin))
	{
	  if (fgets(str,1024,fin)!=NULL)
	    {
	      if (str[0]!='#') { // Check for comment lines
		if (sscanf(str,"%lf %lf",&rx[nr],&ry[nr])==2)
		  {
		    //	      if (rx[nr] >= earliestTime && rx[nr] < latestTime)
		    //		{
		    rx0[nr] = (double)(rx[nr]-psr[0].param[param_pepoch].val[0]);
		    rx[nr] = mjd2year(rx[nr]);
		    ry[nr] *= overlayInvert;
		    ry[nr] /= 1.0e-6;
		    nr++;
		    //		}
		  }
	      }
	    }
	}
      fclose(fin);
      printf("Have read %d lines from the clock file\n",nr);
      // Remove quadratic      
      if (psr[0].ifuncN > 0)
	{
	  double sx[psr[0].ifuncN],sy[psr[0].ifuncN];
	  float fsx[psr[0].ifuncN],fsy[psr[0].ifuncN];
	  double param[3];
	  double chisq=0.0;

	  double dx[psr[0].ifuncN];
	  double dy[psr[0].ifuncN];
	  double de[psr[0].ifuncN];
	  double **cvm;
	  double out_p[3],out_e[3];

	  cvm = (double **)malloc(sizeof(double *)*3);
	  for (i=0;i<3;i++)
	    cvm[i] = (double *)malloc(sizeof(double)*3);

	  // Want to remove the same quadratic as in the pulsar clock function (i.e., same sampling)
	  for (i=0;i<psr[0].ifuncN;i++)
	    {
	      sx[i] = (double)(psr[0].ifuncT[i] - psr[0].param[param_pepoch].val[0]);
	      de[i] = psr[0].ifuncE[i];
	      sy[i] = 0;
	      for (j=0;j<nr-1;j++)
		{
		  if (sx[i] < rx0[j+1]+1 && sx[i] >= rx0[j]-1)
		    {
		      sy[i] = ry[j];
		      break;
		    }
		}
	      printf("overlay: %g %.20f\n",sx[i],sy[i]);
	    }

	  // NOTE: Not using error bars
	  // Should use the constraint weights
	  //
	  TKleastSquares_svd_noErr(sx,sy,psr[0].ifuncN, param, 3, TKfitPoly);  
	  //  TKleastSquares_svd(sx,sy,de,psr[0].ifuncN,param,out_e,3,cvm,&chisq,TKfitPoly,1);
		   //  TKleastSquares_svd(sx,sy,de,psr[0].ifuncN,param,out_e,3,cvm,&chisq,TKfitPoly,1);
	  //	  TKremovePoly_d(sx,sy,psr[0].ifuncN,3);
	  chisq=0.0;
	  for (i=0;i<3;i++)
	    free(cvm[i]);
	  free(cvm);

	  for (i=0;i<psr[0].ifuncN;i++)
	    {
	      printf("sxsy_val: %.15f %.15f %g %g %g %g %g\n",sx[i],sy[i],chisq,psr[0].ifuncV[i],sy[i],sy[i]-(param[0]+param[1]*sx[i]+param[2]*sx[i]*sx[i]),psr[0].ifuncE[i]);
	      
	      chisq+= pow((psr[0].ifuncV[i]-(sy[i]-(param[0]+param[1]*sx[i]+param[2]*sx[i]*sx[i])))/(psr[0].ifuncE[i]),2);
	    }
	  for (i=0;i<nr;i++)
	    {
	      //	      printf("sxsy = %.15f %.15f %.15f",rx[i],ry[i],(param[0]-param[1]*rx0[i]-param[2]*rx0[i]*rx0[i]));
	      ry[i] -= (param[0]+param[1]*rx0[i]+param[2]*rx0[i]*rx0[i]);
	      printf("plotVals: %.15f %.15f %.15f\n",rx[i],rx0[i],ry[i]);
	      //	      printf(" %.15f\n",ry[i]);
	    }
	  printf("Params %.15f %.15f %.15f\n",param[0],param[1],param[2]);
	  printf("Chisq of match between data and overlay = %g\n",chisq);
	  printf("Chisq/(npts) = %g\n",chisq/(double)psr[0].ifuncN);
	  //	  cpgpt(psr[0].ifuncN,fsx,fsy,19);
	}
      else
	TKremovePoly_d(rx,ry,nr,3);

      for (i=0;i<nr;i++)
	{
	  frx[i]=(float)rx[i];
	  fry[i]=(float)ry[i]; // Check not including a MINUS SIGN!!!
	  //	  printf("Expected: %g %g %g\n",frx[i],rx0[i],fry[i]);
	  fprintf(fout,"%g %g %g\n",rx[i],ry[i],rx0[i]+(double)psr[0].param[param_pepoch].val[0]);
	}
      fclose(fout);
      cpgsci(3);
      //      cpgsls(2); 
      //
      // This is the expectation line
      //
      cpgslw(4); cpgline(nr,frx,fry);  cpgslw(2);

      // Now overlay the error bar due to the correlation
      if (covarError == 1)
	{
	  float fre1[nr],fre2[nr];
	  double m,c,ival;
	  for (i=0;i<nr;i++)
	    {
	      // Find interpolated error bar at this time
	      for (k=0;k<psr[0].ifuncN-1;k++)
		{
		  if (rx0[i]+psr[0].param[param_pepoch].val[0] >= psr[0].ifuncT[k] &&
		      rx0[i]+psr[0].param[param_pepoch].val[0] < psr[0].ifuncT[k+1])
		    {
		      m = (psr[0].ifuncE[k]-psr[0].ifuncE[k+1])/(psr[0].ifuncT[k]-psr[0].ifuncT[k+1]);
		      c = psr[0].ifuncE[k]-m*psr[0].ifuncT[k];
		      
		      ival = m*(double)(rx0[i]+psr[0].param[param_pepoch].val[0])+c;
		      break;
		    }
		}
	      fre1[i] = fry[i] - ival;
	      fre2[i] = fry[i] + ival;
	    }
	  cpgsci(8); cpgsls(4);cpgline(nr,frx,fre1);
	  cpgline(nr,frx,fre2); cpgsls(1); cpgsci(1);
	}

      //      cpgsls(1);
      cpgsci(1);
      if (psr[0].param[param_wave_om].paramSet[1]==1)
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
  if (psr[0].param[param_ifunc].paramSet[0]==1)
    {
      float er1[psr[0].ifuncN],er2[psr[0].ifuncN];
      float er3[psr[0].ifuncN],er4[psr[0].ifuncN];
      for (i=0;i<psr[0].ifuncN;i++)
	{	  
	  if (covarError==1)
	    {
	      er1[i] = ey0[i] + whiteErr[i];
	      er2[i] = ey0[i] - whiteErr[i];
	      if (whiteErr[i] > psr[0].ifuncE[i]) {er3[i] = 0.0; er4[i] = 0.0;}
	      else
		{
		  er3[i] = sqrt(psr[0].ifuncE[i]*psr[0].ifuncE[i]-whiteErr[i]*whiteErr[i]);
		  er4[i] = -sqrt(psr[0].ifuncE[i]*psr[0].ifuncE[i]-whiteErr[i]*whiteErr[i]);
		}
	    }
	  else
	    {
	      er1[i] = ey0[i] + psr[0].ifuncE[i];
	      er2[i] = ey0[i] - psr[0].ifuncE[i];
	    }
	}
      printf("plotting %d points\n",npt);
      //      cpgsls(4);
      //      cpgline(npt,px,fye1);
      // cpgline(npt,px,fye2);
      cpgsls(1);


      // Plotting the error bar
            cpgerry(ne,ex,er1,er2,1);
            cpgpt(ne,ex,ey0,4);


      //      if (covarError==1)
      //	{
      //	  cpgsls(4); cpgline(ne,ex,er3);
      //	  cpgline(ne,ex,er4); cpgsls(1);
      //	}
    }
  else
    {
      // Covariance
      cpgsci(1); cpgpt(nptsClk,cvX,cvY,4); cpgerry(nptsClk,cvX,cvY1,cvY2,1); cpgsci(1);
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

char * plugVersionCheck = TEMPO2_h_VER;
