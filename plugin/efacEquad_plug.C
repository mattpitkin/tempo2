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
void calcEfacEquad(double *px,double *py,double *pe,int npts,double *efacRet,double *equadRet,int disp);

void calcEfacEquad2(double *px,double *py,double *pe,int npts,double *efacRet,double *equadRet,int disp,double correctEfac,double correctEquad,double minEquad,double maxEquad,double stepEquad,double minEfac,double maxEfac,double stepEfac,char *grDev);



#define FMAX(x,y) ((x<y)?y:x)
#define EPS1 0.001
#define EPS2 1.0e-8
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50
#define FREE_ARG char*
#define NR_END 1
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

int *ivector(long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void kstwo(double data1[], unsigned long n1, double data2[], unsigned long n2,
	   double *d, double *prob);
double probks(double alam);
double gaussFunc(double val);
void sort(unsigned long n, double arr[]);
void nrerror(char error_text[]);
double erff(double x);
void ksone(double data[], unsigned long n, double (*func)(double), double *d,
	   double *prob);
double gammp(double a, double x);
double gammln(double xx);
void gcf(double *gammcf, double a, double x, double *gln);
void gser(double *gamser, double a, double x, double *gln);

void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,j,k;
  double globalParameter;
  const char *CVS_verNum = "$Revision: 1.6 $";
  char flagID[128];
  int nFlag=0;
  char flagVal[100][1024];
  int found=0;
  double px[MAX_OBSN],py[MAX_OBSN],pe[MAX_OBSN];
  int plot=0;
  int npts=0;
  int lastUsed=0;
  char groupFile[1024];
  int group=0;
  FILE *fout;
  double efac,equad;
  int obsID[MAX_OBSN];
  int l;
  double minx,maxx;
  double dstep = 100;
  int    reqPoints = 10;
  int    usePoints;
  double xpos;
  double correctEfac=-1,correctEquad=-1;
  double maxEquad = 3;
  double minEquad = 0;
  double stepEquad = 0.1;
  double minEfac = 1;
  double maxEfac = 5;
  double stepEfac = 0.1;
  char grDev[128]="/xs";

  printf("Starting the plugin now\n");
  strcpy(flagID,"NULL");

  if (displayCVSversion == 1) CVSdisplayVersion((char *)"efacEquad.C",(char *)"plugin",CVS_verNum);

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: efacEquad\n");
  printf("Author:              G. Hobbs\n");
  printf("CVS Version:         $Revision: 1.6 $\n");
  printf(" --- type 'h' for help information\n");


  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[0],argv[i+1]); 
	  strcpy(timFile[0],argv[i+2]);
	}
      else if (strcmp(argv[i],"-reqPoints")==0)
	sscanf(argv[++i],"%d",&reqPoints);
      else if (strcmp(argv[i],"-dstep")==0)
	sscanf(argv[++i],"%lf",&dstep);
      else if (strcmp(argv[i],"-flag")==0)
	strcpy(flagID,argv[++i]);
      else if (strcmp(argv[i],"-plot")==0)
	plot=1;
      else if (strcmp(argv[i],"-correct")==0)
	{
	  sscanf(argv[++i],"%lf",&correctEfac);
	  sscanf(argv[++i],"%lf",&correctEquad);
	}
      else if (strcmp(argv[i],"-group")==0)
	{
	  strcpy(groupFile,argv[++i]);
	  group=1;
	}
      else if (strcasecmp(argv[i],"-minEfac")==0)
	sscanf(argv[++i],"%lf",&minEfac);
      else if (strcasecmp(argv[i],"-maxEfac")==0)
	sscanf(argv[++i],"%lf",&maxEfac);
      else if (strcasecmp(argv[i],"-stepEfac")==0)
	sscanf(argv[++i],"%lf",&stepEfac);
      else if (strcasecmp(argv[i],"-minEquad")==0)
	sscanf(argv[++i],"%lf",&minEquad);
      else if (strcasecmp(argv[i],"-maxEquad")==0)
	sscanf(argv[++i],"%lf",&maxEquad);
      else if (strcasecmp(argv[i],"-stepEquad")==0)
	sscanf(argv[++i],"%lf",&stepEquad);
      else if (strcasecmp(argv[i],"-grDev")==0)
	strcpy(grDev,argv[++i]);
    }
  if (strcmp(flagID,"NULL")==0)
    {
      printf("Must use the -flag option to define a particular flag\n");
      exit(1);
    }


  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);
  minx = maxx = 0.0;
  for (i=0;i<psr[0].nobs;i++)
    {
      if (psr[0].obsn[i].deleted==0)
	{
	  if (minx == 0 && maxx == 0)
	    {
	      minx = maxx = (double)psr[0].obsn[i].sat;
	    }
	  else 
	    {
	      if (minx > (double)psr[0].obsn[i].sat) minx = (double)psr[0].obsn[i].sat;
	      if (maxx < (double)psr[0].obsn[i].sat) maxx = (double)psr[0].obsn[i].sat;
	      if (psr[0].obsn[i].toaErr > maxEquad){
		printf("WARNING: maxEquad = %g is less than error bar size for observation %d of %g. Probably should increase EQUAD\n",maxEquad,i,psr[0].obsn[i].toaErr);
	      }
	    }
	}
    }

  // Turn off all fitting
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (j=0;j<psr[0].param[i].aSize;j++)
	{
	  if (psr[0].param[i].fitFlag[j]==1)
	    psr[0].param[i].fitFlag[j]=0;
	}
    }
  for (i=1;i<=psr[0].nJumps;i++)
    {
      if (psr[0].fitJump[i]==1)
	psr[0].fitJump[i]=0;
    }

  // 2. Turn on fitting for F0 and F1
  psr[0].param[param_f].fitFlag[0] = 1;
  psr[0].param[param_f].fitFlag[1] = 1;
  
  // 3. Add in constrained IFUNC values
  psr[0].param[param_ifunc].paramSet[0]=1;
  psr[0].param[param_ifunc].fitFlag[0]=1;
  psr[0].param[param_ifunc].val[0]=2;
  
  // Set up IFUNCS
  psr[0].ifuncN=0;
  xpos = minx-1;
  psr[0].ifuncT[psr[0].ifuncN] = xpos;
  psr[0].ifuncV[psr[0].ifuncN] = 0.0;
  psr[0].ifuncE[psr[0].ifuncN] = 0.0;
  psr[0].ifuncN++;
  lastUsed=0;
  // Assume time sorted points
  printf("Using dstep = %f\n",dstep);
  for (i=0;i<psr[0].nobs;i++)
    {
      if (i-lastUsed > reqPoints && (psr[0].obsn[i].sat - xpos) > dstep)
	{
	  xpos = minx-1;
	  lastUsed = i;
	  xpos = psr[0].obsn[i].sat + 1;
	  psr[0].ifuncT[psr[0].ifuncN] = xpos;
	  psr[0].ifuncV[psr[0].ifuncN] = 0.0;
	  psr[0].ifuncE[psr[0].ifuncN] = 0.0;
	  psr[0].ifuncN++;	  
	}
    }
  // Check last point
  if (xpos < psr[0].obsn[psr[0].nobs-1].sat)
    {
	  psr[0].ifuncT[psr[0].ifuncN] = psr[0].obsn[psr[0].nobs-1].sat + 1;
	  psr[0].ifuncV[psr[0].ifuncN] = 0.0;
	  psr[0].ifuncE[psr[0].ifuncN] = 0.0;
	  psr[0].ifuncN++;	  
    }

  psr[0].constraints[psr[0].nconstraints++]= constraint_ifunc_0;
  psr[0].constraints[psr[0].nconstraints++]= constraint_ifunc_1;
  psr[0].constraints[psr[0].nconstraints++]= constraint_ifunc_2;
  

  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      if (i==0) doFit(psr,*npsr,0);   /* Do the fitting     */
      else textOutput(psr,*npsr,globalParameter,0,0,1,(char *)"efacEquad_try.par");  /* Display the output */
      //      else textOutput(psr,*npsr,globalParameter,0,0,0,(char *)"");  /* Display the output */
    }
  printf("Step 1\n");
  

  // Determine how many different flags to process
  if (group==0)
    {
      for (i=0;i<psr[0].nobs;i++)
	{
	  for (k=0;k<psr[0].obsn[i].nFlags;k++)
	    {
	      if (strcmp(psr[0].obsn[i].flagID[k],flagID)==0)
		{
		  found=0;
		  for (j=0;j<nFlag;j++)
		    {
		      if (strcmp(psr[0].obsn[i].flagVal[k],flagVal[j])==0)
			{found=1; break;}
		    }
		  if (found==0)
		    {strcpy(flagVal[nFlag],psr[0].obsn[i].flagVal[k]); nFlag++;}
		}
	    }
	}
    }
  else
    {
      FILE *fin;
      char line[1024];

      nFlag=0;
      if (!(fin = fopen(groupFile,"r")))
	{
	  printf("Unable to open file %s\n",groupFile);
	  exit(1);
	}
      while (!feof(fin))
	{
	  if (fgets(line,1024,fin)!=NULL)
	    {
	      if (strlen(line)>1)
		{
		  strcpy(flagVal[nFlag],line);
		  flagVal[nFlag][strlen(flagVal[nFlag])-1]='\0';
		  printf("Got flag: >%s<\n",flagVal[nFlag]);
		  nFlag++;
		}
	    }
	}
      fclose(fin);
    }
  printf("Step 2\n");
  fout = fopen("efacEquad_output.dat","w");
  for (i=0;i<nFlag;i++)
    {
      printf("%d %s %s ",i+1,flagID,flagVal[i]);
      // Obtain the data set
      //      for (l=0;l<2;l++)
	{
	  npts=0;
	  for (j=0;j<psr[0].nobs;j++)
	    {
	      for (k=0;k<psr[0].obsn[j].nFlags;k++)
		{
		  //	      printf("Checking >%s< >%s< >%s<\n",flagVal[i],psr[0].obsn[j].flagID[k],psr[0].obsn[j].flagVal[k]);
		  //	      printf("Result %s\n",strstr(flagVal[i],psr[0].obsn[j].flagVal[k]));
		  if (psr[0].obsn[i].deleted==0)
		    {
		      //		      printf("Have %s %s %s %s\n",psr[0].obsn[j].flagID[k],flagID,psr[0].obsn[j].flagVal[k],flagVal[i]);
		      if (strcmp(psr[0].obsn[j].flagID[k],flagID)==0
			  && ((group==0 && strcmp(psr[0].obsn[j].flagVal[k],flagVal[i])==0)
			      || (group==1 && strcmp(flagVal[i],psr[0].obsn[j].flagVal[k])==0)))
			{
			  //			  		  printf("In here with %s %s\n",psr[0].obsn[j].flagVal[k],flagVal[i]);
			  // Check for outliers
			  /*		  if (l==0)
			    {
			      px[npts] = (double)(psr[0].obsn[j].sat-psr[0].param[param_pepoch].val[0]);
			      py[npts] = (double)psr[0].obsn[j].residual;
			      pe[npts] = (double)psr[0].obsn[j].toaErr*1.0e-6;
			      obsID[npts] = j;
			      npts++;
			    }
			    else if (l==1) */
			    {
			      double err;
			      px[npts] = (double)(psr[0].obsn[j].sat-psr[0].param[param_pepoch].val[0]);
			      py[npts] = (double)psr[0].obsn[j].residual;
			      pe[npts] = (double)psr[0].obsn[j].toaErr*1.0e-6;
			      //			      printf("Error = %g\n",pe[npts]);
			      obsID[npts] = j;
			      //			      err = sqrt(pe[i]*pe[i]+equad*1e-6*equad*1e-6)*efac;
			      //			      if (fabs(py[npts]/err) < 40)
			      npts++;
				//			      else
				//				printf("Removing outlier\n");
			    }
			}
		    }
		}
	    }
	  //	  calcEfacEquad(px,py,pe,npts,&efac,&equad,1);
	  //  printf("npts = %d ",npts);
  printf("Step 3\n");
	  calcEfacEquad2(px,py,pe,npts,&efac,&equad,plot,correctEfac,correctEquad,minEquad,maxEquad,stepEquad,minEfac,maxEfac,stepEfac,grDev);
  printf("Step 4\n");
	  //	  if (l==1){
	    printf("npts = %d ",npts);
	    printf("\n");
	    if (group==0)
	      {
		if (efac!=1.0)
		  fprintf(fout,"T2EFAC %s %s %g\n",flagID,flagVal[i],efac);
		if (equad!=0.0)
		  fprintf(fout,"T2EQUAD %s %s %g\n",flagID,flagVal[i],equad);
	      }
	    else if (group==1)
	      {
		char *pch;
		pch = strtok(flagVal[i]," ");
		
		while (pch != NULL)
		  {
		    fprintf(fout,"T2EFAC %s %s %g\n",flagID,pch,efac);
		    fprintf(fout,"T2EQUAD %s %s %g\n",flagID,pch,equad);
		    pch = strtok(NULL," ");
		  }
	      }
	    //	  }
	}
    }
  fclose(fout);
  return 0;
}

void calcEfacEquad(double *px,double *py,double *pe,int npts,double *efacRet,double *equadRet,int disp,char *grDev)
{
  int i;
  double nr[npts];
  double sx=0.0,sx2=0.0;
  double var;
  double equad = 0.0;
  double efac;
  double diff,bestDiff,bestEfac;
  double bestEquad,bestVar;
  double err;
  double efac1,efac2,equad1,equad2,var1,var2;
  double diff1,diff2;
  double maxEquad = 5;
  int t=0;
  FILE *fout;
  printf("Starting search\n");
  // Searching over a full grid of EFAC and EQUAD does not work
  // Must first hold EFAC fixed and then find the best EQUAD
  //    for (efac = 0.1;efac < 4;efac +=0.1)
  fout = fopen("equadData1.dat","w");
  efac = 1.0;
    {
      for (equad = 0.0;equad <= maxEquad;equad+=0.1)
	{
	  //	  efac = 1.0; equad=0;
	  // Calculate normalised residuals
	  sx = sx2 = 0.0;
	  for (i=0;i<npts;i++)
	    {
	      err = sqrt(pe[i]*pe[i]+equad*1e-6*equad*1e-6)*efac;
	      nr[i] = py[i]/err;
	      //	        printf("%g %g %g\n",py[i],err,nr[i]);
	      sx+=nr[i];
	      sx2+=(nr[i]*nr[i]);
	    }  
	  var = (sx2/(double)npts - pow(sx/(double)npts,2));
	  if (t==0)
	    {
	      bestEquad = equad;
	      bestEfac = efac;
	      bestDiff = diff = pow(var-1,2);
	      bestVar = var;
	    }
	  else
	    {
	      diff = pow(var-1.0,2);
	      if (diff < bestDiff)
		{
		  bestDiff = diff;
		  bestEquad = equad;
		  bestEfac = efac;
		  bestVar = var;
		}
	    }
	  t++;
	  //	  printf("%g %g %g %g %g %g %g %g %d\n",efac,equad,var,diff,bestDiff,bestEquad,sx,sx2,npts);
	  fprintf(fout,"%g %g\n",equad,var);
	}
    }
    equad = bestEquad;
    diff = bestDiff;

    fclose(fout);
    fout = fopen("efacData1.dat","w");
    for (efac = 0.1;efac < 4;efac+=0.1)
      {	
	//      for (equad = 0.0;equad < 3;equad+=0.1)
	{
	  //	  efac = 1.0; equad=0;
	  // Calculate normalised residuals
	  sx = sx2 = 0.0;
	  for (i=0;i<npts;i++)
	    {
	      err = sqrt(pe[i]*pe[i]+equad*1e-6*equad*1e-6)*efac;
	      nr[i] = py[i]/err;
	      //	        printf("%g %g %g\n",py[i],err,nr[i]);
	      sx+=nr[i];
	      sx2+=(nr[i]*nr[i]);
	    }  
	  var = (sx2/(double)npts - pow(sx/(double)npts,2));
	  diff = pow(var-1.0,2);
	  if (diff < bestDiff)
	    {
	      bestDiff = diff;
	      bestEquad = equad;
	      bestEfac = efac;
	      bestVar = var;
	    }
	  t++;
	  //	  printf("%g %g %g %g %g %g %g %g %d\n",efac,equad,var,diff,bestDiff,bestEquad,sx,sx2,npts);
	  fprintf(fout,"%g %g\n",efac,var);	  
	}
    }
    fclose(fout);
    efac1 = bestEfac;
    equad1 = bestEquad;
    var1 = bestVar;

    // Now try the other way around and see if EFAC is dominating
    fout = fopen("efacData2.dat","w");
    t=0;
    equad = 0.0;
    {
      for (efac = 0.1;efac < 4;efac+=0.1)
	{
	  //	  efac = 1.0; equad=0;
	  // Calculate normalised residuals
	  sx = sx2 = 0.0;
	  for (i=0;i<npts;i++)
	    {
	      err = sqrt(pe[i]*pe[i]+equad*1e-6*equad*1e-6)*efac;
	      nr[i] = py[i]/err;
	      //	        printf("%g %g %g\n",py[i],err,nr[i]);
	      sx+=nr[i];
	      sx2+=(nr[i]*nr[i]);
	    }  
	  var = (sx2/(double)npts - pow(sx/(double)npts,2));
	  if (t==0)
	    {
	      bestEquad = equad;
	      bestEfac = efac;
	      bestDiff = diff = pow(var-1,2);
	      bestVar = var;
	    }
	  else
	    {
	      diff = pow(var-1.0,2);
	      if (diff < bestDiff)
		{
		  bestDiff = diff;
		  bestEquad = equad;
		  bestEfac = efac;
		  bestVar = var;
		}
	    }
	  t++;
	  //	  printf("%g %g %g %g %g %g %g %g %d\n",efac,equad,var,diff,bestDiff,bestEquad,sx,sx2,npts);
	  fprintf(fout,"%g %g\n",efac,var);
	}
    }
    efac = bestEfac;
    diff = bestDiff;
    fclose(fout);
    fout = fopen("equadData2.dat","w");
    for (equad = 0.0;equad <= maxEquad;equad+=0.1)
      {	
	//      for (equad = 0.0;equad < 3;equad+=0.1)
	{
	  //	  efac = 1.0; equad=0;
	  // Calculate normalised residuals
	  sx = sx2 = 0.0;
	  for (i=0;i<npts;i++)
	    {
	      err = sqrt(pe[i]*pe[i]+equad*1e-6*equad*1e-6)*efac;
	      nr[i] = py[i]/err;
	      //	        printf("%g %g %g\n",py[i],err,nr[i]);
	      sx+=nr[i];
	      sx2+=(nr[i]*nr[i]);
	    }  
	  var = (sx2/(double)npts - pow(sx/(double)npts,2));
	  diff = pow(var-1.0,2);
	  if (diff < bestDiff)
	    {
	      bestDiff = diff;
	      bestEquad = equad;
	      bestEfac = efac;
	      bestVar = var;
	    }
	  t++;
	  //	  printf("%g %g %g %g %g %g %g %g %d\n",efac,equad,var,diff,bestDiff,bestEquad,sx,sx2,npts);
	  fprintf(fout,"%g %g\n",equad,var);	  
	}
    }
    fclose(fout);
    efac2 = bestEfac;
    equad2 = bestEquad;
    var2 = bestVar;
    diff1 = pow(var1-1.0,2);
    diff2 = pow(var2-1.0,2);
    if (diff2 < diff1) {
      bestVar = var2;
      bestEfac = efac2;
      bestEquad = equad2;
    }
    else
      {
	bestVar = var1;
	bestEfac = efac1;
	bestEquad = equad1;
      }
  if (disp==1) {
    printf("norm. var. = %g, equad = %g efac = %g ",bestVar,bestEquad,bestEfac);
    printf("norm. var1. = %g, equad1 = %g efac1 = %g ",var1,equad1,efac1);
    printf("norm. var2. = %g, equad2 = %g efac2 = %g ",var2,equad2,efac2);
  }
  //    exit(1);
  *efacRet =  bestEfac;
  *equadRet = bestEquad;
  printf("Complete search\n");
}

void calcEfacEquad2(double *px,double *py,double *pe,int npts,double *efacRet,double *equadRet,int disp,double correctEfac,double correctEquad,double minEquad,double maxEquad,double stepEquad,double minEfac,double maxEfac,double stepEfac,char *grDev)
{
  int i;
  double nr[npts];
  double sx=0.0,sx2=0.0;
  double var;
  double equad = 0.0;
  double efac = 0.0;
  double diff,bestDiff,bestEfac;
  double bestEquad,bestVar;
  double err;
  double efac1,efac2,equad1,equad2,var1,var2;
  double diff1,diff2;
  double maxV = 0;
  double minV = 0;
  float col;
  int t=0;
  float mx,my;
  char key;

  FILE *fout;
  FILE *fout2;
  FILE *fout3;
  FILE *fout4;

  float fx[8192*3],fy[8192*3],fv[8192*3];
  int n=0;
  double se=0;
  double sw;

  double w[npts];
  double meanErr=0.0;
  double meanVal=0.0;
  long seed = TKsetSeed();
  double d;
  double gauss[8192*4];
  double chisq;
  int nval;

  printf("Starting search\n");

  chisq = 0.0;
  for (i=0;i<npts;i++)
    chisq += pow(py[i]/pe[i],2);
  if (chisq/(double)npts < 1)
    {
      printf("Reduced chisq is %g which is less than 1\n",chisq/(double)npts);
      *efacRet = sqrt(chisq/(double)npts);
      *equadRet = 0.0;
      return;
    }
  printf("Starting %g %g\n",correctEfac,correctEquad);
    fout = fopen("gauss.dat","w");
    for (i=0;i<8192*4;i++)
      {
        gauss[i]=TKgaussDev(&seed);
        fprintf(fout,"%g\n",gauss[i]);
      }
    fclose(fout);
  // Searching over a full grid of EFAC and EQUAD does not work
  // Must first hold EFAC fixed and then find the best EQUAD
  //    for (efac = 0.1;efac < 4;efac +=0.1)
    //  printf("Starting2 %g %g\n",correctEfac,correctEquad);
  for (i=0;i<npts;i++)
    meanErr+=pe[i];
  meanErr/=(double)npts;
    printf("Starting3 %g %g\n",correctEfac,correctEquad);
  fout = fopen("equadEfac.dat","w");
  fout2 = fopen("norm.dat","w");
  fout3 = fopen("normCorrect.dat","w");
  fout4 = fopen("orig.dat","w");


  nval = (int)((maxEfac-minEfac)/stepEfac*(maxEquad-minEquad)/stepEquad+1);
  printf("nval = %d\n",nval);
  do {
    stepEquad *= 2.0;
    nval = (int)((maxEfac-minEfac)/stepEfac*(maxEquad-minEquad)/stepEquad+1);
    printf("nval = %d\n",nval);
  } while (nval > 8192*3);

  printf("Starting4 %g %g\n",correctEfac,correctEquad);
  printf("Number of grid points = %g\n",(maxEfac-minEfac)/stepEfac*(maxEquad-minEquad)/stepEquad);


  
  
  for (efac = minEfac;efac < maxEfac;efac+=stepEfac)
  //   efac = 2;
    {
      //      equad = 3.5;
      for (equad = minEquad;equad <= maxEquad;equad+=stepEquad)
	{
	  //	  efac = 1.0; equad=0;
	  // Calculate normalised residuals
	  sx = sx2 = se = 0.0;
	  meanVal = 0.0;
	  //	  printf("%g %g %d\n",efac,equad,npts);
	  //	  printf("place 1 %d\n",npts);
	  for (i=0;i<npts;i++)
	    {
	      //	      printf("Trying %g %g %g\n",pe[i],equad,efac);
	      err = sqrt(pe[i]*pe[i]+equad*1e-6*equad*1e-6)*efac;
	      //	      printf("Err = %g\n",err);
	      //err = sqrt(pow(pe[i]*efac,2)+equad*1e-6*equad*1e-6);
	      //	      err = sqrt(pe[i]*pe[i]+equad*1e-6*equad*1e-6)+pe[i]*efac;
	      //	      printf("py = %g %d %g %d\n",py[i],i,err,npts);
	      nr[i] = py[i]/err;
	      //	      printf("%g %g %g\n",py[i],err,nr[i]);
	      meanVal += nr[i];
	      w[i] = pow(pe[i]-meanErr,2);
	      fprintf(fout2,"efac = %g equad = %g %g\n",efac,equad,nr[i]);
	      //	      printf("Have: %g %g %g %g\n",efac,correctEfac,equad,correctEquad);
	      if (efac==1 && equad==0)
		  fprintf(fout4,"%g %g %g %g\n",px[i],py[i],err,nr[i]);
	      //	      printf("This pos\n");
	      if (fabs(efac-correctEfac) < 1e-5 && fabs(equad-correctEquad)<1e-5)
		  fprintf(fout3,"%g %g %g %g\n",px[i],py[i],err,nr[i]);
	      //	      printf("This pos2\n");
	      //  w[i] = 1.0;
	      sx+=nr[i];
	      sx2+=(nr[i]*nr[i]);
	      se+=1;
	      //	      printf("This pos3\n");
	    }  
	  //	  printf("Here with %g %g\n",efac,equad);
	  meanVal /= (double)npts;
	  for (i=0;i<npts;i++)
	    nr[i]-=meanVal;
	  //	  kstwo(nr-1,npts,gauss-1,8192*4,&d,&var);
	  ksone(nr-1,npts,gaussFunc,&d,&var);
	  if (fabs(efac-correctEfac) < 1e-5 && fabs(equad-correctEquad)<1e-5)
	    fprintf(fout3,"# %g\n",var);
	  if (efac==1 && equad==0)
	    fprintf(fout4,"# %g\n",var);
	  //	  printf("var = %g\n",var);
	  if (t==0)
	    {
	      bestEquad = equad;
	      bestEfac = efac;
	      //	      bestDiff = diff = pow(var-1,2);
	      bestDiff = diff = var;
	      bestVar = var;
	    }
	  else
	    {
	      //	      diff = pow(var-1.0,2);
	      diff = var;
	      //	      if (diff < bestDiff)
	      if (diff > bestDiff)
		{
		  bestDiff = diff;
		  bestEquad = equad;
		  bestEfac = efac;
		  bestVar = var;
		}
	    }
	  fx[n] = efac;
	  fy[n] = equad;
	  //	  fv[n] = log10(1.0/pow(var-1,2));
	  fv[n] = var; //log10(1.0/var); //pow(var-1,2));
	  if (t==0){
	    minV = maxV = fv[n];
	  }
	  if (fv[n] > maxV) maxV = fv[n];
	  
	  if (fv[n] < minV) minV = fv[n];
	  n++;
	  t++;
	  //	  printf("%g %g %g %g %g %g %g %g %d\n",efac,equad,var,diff,bestDiff,bestEquad,sx,sx2,npts);
	  fprintf(fout,"efac = %g equad = %g %g %g\n",efac,equad,var,1.0/pow(var-1.0,2));
	}
    }
  fclose(fout2);
  fclose(fout3);
  fclose(fout4);
  printf("Max = %g %g\n",maxV,minV);
  printf("Number of points = %d\n",npts);
  printf("T2EFAC = %g\n",bestEfac);
  printf("T2EQUAD = %g\n",bestEquad);
    equad = bestEquad;
    efac = bestEfac;
    diff = bestDiff;

    fclose(fout);
    if (disp==1) {
      printf("norm. var. = %g, equad = %g efac = %g ",bestVar,bestEquad,bestEfac);
      // Make the figure
      cpgbeg(0,grDev,1,1);
      cpgenv(minEfac-(maxEfac-minEfac)*0.1,maxEfac+(maxEfac-minEfac)*0.1,minEquad,maxEquad,0,1);
      cpglab("EFAC","EQUAD (\\gms)","");
      for (i=0;i<n;i++)
	{
	  //	  printf("Selecting %g\n",fv[i]/maxV);
	  col = (fv[i]-minV)/(maxV-minV);
	  //	  	  printf("%d %g col = %g\n",i,fv[i],col);
	  //col = fv[i];
	  //	  cpgscr(2,1-col,1-col,1-col);
	  cpgscr(2,col,col,col);
	  cpgsci(2);
	  cpgsch(3*col);
	  cpgpt(1,fx+i,fy+i,-15);
	}
      fx[0] = bestEfac;
      fy[0] = bestEquad;
      cpgsch(1.4); cpgsci(3);cpgpt(1,fx,fy,6); cpgsci(1);
      fx[0] = correctEfac;
      fy[0] = correctEquad;
      cpgsch(4); cpgsci(7);cpgpt(1,fx,fy,18); cpgsci(1); cpgsch(1);
      cpgbox("ABCTS",0,0,"ABCTS",0,0);
      cpgcurs(&mx,&my,&key);
      cpgend();
    }
  //    exit(1);
  *efacRet =  bestEfac;
  *equadRet = bestEquad;

  printf("Complete search\n");
}

char * plugVersionCheck = TEMPO2_h_VER;


void kstwo(double data1[], unsigned long n1, double data2[], unsigned long n2,
	double *d, double *prob)
{
	double probks(double alam);
	void sort(unsigned long n, double arr[]);
	unsigned long j1=1,j2=1;
	double d1,d2,dt,en1,en2,en,fn1=0.0,fn2=0.0;

	sort(n1,data1);
	sort(n2,data2);
	en1=n1;
	en2=n2;
	*d=0.0;
	while (j1 <= n1 && j2 <= n2) {
		if ((d1=data1[j1]) <= (d2=data2[j2])) fn1=j1++/en1;
		if (d2 <= d1) fn2=j2++/en2;
		if ((dt=fabs(fn2-fn1)) > *d) *d=dt;
	}
	en=sqrt(en1*en2/(en1+en2));
	*prob=probks((en+0.12+0.11/en)*(*d));
}
/* (C) Copr. 1986-92 Numerical Recipes Software #p21E6W)1.1&iE10(9p#. */

double probks(double alam)
{
	int j;
	double a2,fac=2.0,sum=0.0,term,termbf=0.0;

	a2 = -2.0*alam*alam;
	for (j=1;j<=100;j++) {
		term=fac*exp(a2*j*j);
		sum += term;
		if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
		fac = -fac;
		termbf=fabs(term);
	}
	return 1.0;
}


void sort(unsigned long n, double arr[])
{
	unsigned long i,ir=n,j,k,l=1;
	int jstack=0,*istack;
	double a,temp;

	istack=ivector(1,NSTACK);
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				for (i=j-1;i>=1;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
				}
				arr[i+1]=a;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j]);
			}
			arr[l]=arr[j];
			arr[j]=a;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}
/* (C) Copr. 1986-92 Numerical Recipes Software #p21E6W)1.1&iE10(9p#. */

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


void ksone(double data[], unsigned long n, double (*func)(double), double *d,
	double *prob)
{
	double probks(double alam);
	void sort(unsigned long n, double arr[]);
	unsigned long j;
	double dt,en,ff,fn,fo=0.0;

	sort(n,data);
	en=n;
	*d=0.0;
	for (j=1;j<=n;j++) {
		fn=j/en;
		ff=(*func)(data[j]);
		dt=FMAX(fabs(fo-ff),fabs(fn-ff));
		if (dt > *d) *d=dt;
		fo=fn;
	}
	en=sqrt(en);
	*prob=probks((en+0.12+0.11/en)*(*d));
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software #p21E6W)1.1&iE10(9p#. */

double erff(double x)
{
	double gammp(double a, double x);

	return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}
/* (C) Copr. 1986-92 Numerical Recipes Software #p21E6W)1.1&iE10(9p#. */

double gammp(double a, double x)
{
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	void nrerror(char error_text[]);
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammp");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software #p21E6W)1.1&iE10(9p#. */

#include <math.h>

void gcf(double *gammcf, double a, double x, double *gln)
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
/* (C) Copr. 1986-92 Numerical Recipes Software #p21E6W)1.1&iE10(9p#. */

#include <math.h>

void gser(double *gamser, double a, double x, double *gln)
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}
#undef ITMAX
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software #p21E6W)1.1&iE10(9p#. */

#include <math.h>

double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
/* (C) Copr. 1986-92 Numerical Recipes Software #p21E6W)1.1&iE10(9p#. */


double gaussFunc(double val)
{
  return 0.5*(1+erff(val/sqrt(2)));
}
