//  Copyright (C) 2006,2007,2008,2009,2010,2011 George Hobbs, Russell Edwards

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
#include "TKspectrum.h"

using namespace std;

bool write_debug_files=true;
double OMEGA0=0;
double GLOBAL_MEANSUB=0;
double GLOBAL_COSVAL=0;

void getSpectrum(pulsar *psr,double *px,double *py_r,double *py_i,int *nSpec,double toffset,double startOverlap,double endOverlap,double stepMJD,char* covarFuncFile);
void formCholeskyMatrixPlugin(double *c,double *resx,double *resy,double *rese,int np,double **uinv);
int calcSpectra_plugin(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY_R,double *specY_I,int nfit);
double psrangle(double centre_long,double centre_lat,double psr_long,double psr_lat);
void hdfunc(double x,double *p,int ma);
void hdfunc_offs(double x,double *p,int ma);
void hdfunc_meanSub(double x,double *p,int ma);
void hdfunc_removeCosine(double x,double *p,int ma);
void hdfunc_cosineSub(double x,double *p,int ma);
void cosineFunc(double x,double *p,int ma);

void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,p,n,p1,p2;
  double globalParameter;
  double startMJD,endMJD,stepMJD,t;
  char tstr[128]="";
  FILE *fout,*fin,*fout2;
  int  nSpec,nSpec1,nSpec2;
  double px1[4096],py_r1[4096],py_i1[4096];
  double px2[4096],py_r2[4096],py_i2[4096];
  double crossX[4096],crossY_r[4096],crossY_i[4096];
  double toverlap;
  
  double angle[4096];
  double a2zeta[4096],a2zeta_im[4096],a2zeta_e[4096];
  long double sum1,sum2,sum3,sum4,weight1;
  double cx,zeta;
  int    npair=0;
  char fname[4096];

  startMJD = 53000.0;
  endMJD   = 54806.0;
  stepMJD  = 100.0;

  char dummy[4096];

  double alpha1,fc1,amp1,pw1,alpha2,fc2,amp2,pw2;

  double guessGWamp = -1;
  double alpha_res = -13.0/3.0;
  double pg,modelPwr1,modelPwr2;
  double crossPowerErr[4096];
  double startOverlap,endOverlap;

  double red_alpha1,red_pwr1;
  double red_alpha2,red_pwr2;

  double toffset;
  double maxOverlap = 300; // Only consider pairs where the overlap is larger than this value
  char dcmFile[128];
  char covarFuncFile[128];
  char newname[128];

  int useRed=0;
  int nSpecNdof;
  int specStart;

  int giveHD=0;
  char hdFile[128];


  strcpy(dcmFile,"NULL");
  strcpy(covarFuncFile,"PSRJ");

  *npsr = 0;  

  printf("Graphical Interface: detectGWB\n");
  printf("Author:              G. Hobbs\n");
  printf("Version:             1.0\n");
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
      else if (strcmp(argv[i],"-gwamp")==0)
	sscanf(argv[++i],"%lf",&guessGWamp);
      else if (strcmp(argv[i],"-alpha_res")==0)
	sscanf(argv[++i],"%lf",&alpha_res);
      else if (strcmp(argv[i],"-t1")==0)
	sscanf(argv[++i],"%lf",&startMJD);
      else if (strcmp(argv[i],"-t2")==0)
	sscanf(argv[++i],"%lf",&endMJD);
      else if (strcmp(argv[i],"-dt")==0)
	sscanf(argv[++i],"%lf",&stepMJD);
      else if (strcmp(argv[i],"-dcf")==0)
	strcpy(covarFuncFile,argv[++i]);
      else if (strcasecmp(argv[i],"-usered")==0)
	useRed=1;
      else if (strcasecmp(argv[i],"-fast")==0)
	write_debug_files=false;
      else if (strcasecmp(argv[i],"-hd")==0)
	{
	  giveHD=1;
	  strcpy(hdFile,argv[++i]);
	}

    }

  if (giveHD==0)
    {
      toffset = 0.5*(startMJD+endMJD);

      if (guessGWamp == -1)
	{
	  printf("Must use the -gwamp command-line option to specify\n");
	  printf("a gravitational wave amplitude for use in the weighting\n");
	  exit(1);
	}
      
      //
      // Standard tempo2 processing to obtain post-fit timing residuals
      //
      readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
      readTimfile(psr,timFile,*npsr);         /* Load the arrival times    */
      preProcess(psr,*npsr,argc,argv);
      // Create some IFUNCs that will represent the low frequency signal
      // ensure that these IFUNCs are constrained not to include 
      // the quadratic etc.
      for (p=0;p<*npsr;p++)
	{
	  // Check if IFUNCs already exist
	  if (psr[p].ifuncN > 0)
	    {
	      printf("Pulsar: %s already has IFUNC parameters\n",psr[p].name);
	      exit(1);
	    }
	  n=0;
	  printf("Creating ifuncs from %g to %g in steps of %g\n",(double)startMJD,(double)endMJD,(double)stepMJD);
	  for (t = startMJD; t<=endMJD;t+=stepMJD)
	    {
	      if (t+stepMJD >= psr[p].obsn[0].sat && t-stepMJD <= psr[p].obsn[psr[p].nobs-1].sat)
		//	  if (t >= psr[p].obsn[0].sat && t <= psr[p].obsn[psr[p].nobs-1].sat)
		{
		  psr[p].ifuncT[n] = t;
		  psr[p].ifuncV[n] = 0.0;
		  psr[p].ifuncE[n] = 0.0;
		  n++;
		}
	    }
	  if (n < 2){
	    printf("\n\n\nERROR: For some reason we cannot determine the interpolated function. Please check the -t1, -t2, -dt options and also the start and finish arrival times in your .tim files.\n\n\n");
	    exit(1);
	  }
	  psr[p].ifuncN = n;
	  psr[p].param[param_ifunc].fitFlag[0] = 1;
	  psr[p].param[param_ifunc].paramSet[0] = 1;
	  psr[p].param[param_ifunc].val[0] = 2;
	  psr[p].constraints[psr[p].nconstraints++] = constraint_ifunc_0;
	  psr[p].constraints[psr[p].nconstraints++] = constraint_ifunc_1;
	  psr[p].constraints[psr[p].nconstraints++] = constraint_ifunc_2;
	  
	  /*      psr[p].constraints[psr[p].nconstraints++] = constraint_ifunc_year_sin;
		  psr[p].constraints[psr[p].nconstraints++] = constraint_ifunc_year_cos;
		  psr[p].constraints[psr[p].nconstraints++] = constraint_ifunc_year_xsin;
		  psr[p].constraints[psr[p].nconstraints++] = constraint_ifunc_year_xcos;
		  psr[p].constraints[psr[p].nconstraints++] = constraint_ifunc_year_sin2;
		  psr[p].constraints[psr[p].nconstraints++] = constraint_ifunc_year_cos2; */
	  sprintf(newname,"%s.new.par",psr[p].name);
	  if(write_debug_files)textOutput(psr+p,1,0,0,0,1,newname);
	}
      
      //
      // Do the standard pulsar fits and get record the IFUNCs
      //    
      printf("Currently npts = %d\n",psr[0].nobs);
      for (i=0;i<2;i++)                       /* Do two iterations for pre- and post-fit residuals*/
	{
	  formBatsAll(psr,*npsr);             /* Form the barycentric arrival times */
	  formResiduals(psr,*npsr,1);         /* Form the residuals                 */
	  //      if (i==0) doFitDCM(psr,dcmFile,covarFuncFile,*npsr,0);       /* Do the fitting     */
	  if (i==0) doFit(psr,*npsr,0);       /* Do the fitting     */
	  else if (write_debug_files)textOutput(psr,*npsr,globalParameter,0,0,0,tstr);  /* Display the output */
	}
      printf("Completed initial fit\n");
      for (p=0;p<*npsr;p++)
	{
	  if(write_debug_files){
	    sprintf(newname,"%s.afterfit.par",psr[p].name);
	    textOutput(psr+p,1,0,0,0,1,newname);
	  }
	  
	  sprintf(tstr,"%s.ifuncDGW",psr[p].name);
	  fout = fopen(tstr,"w");
	  for (i=0;i<psr[p].ifuncN;i++)
	    fprintf(fout,"%.2f %.10g %.10g\n",psr[p].ifuncT[i],psr[p].ifuncV[i],psr[p].ifuncE[i]);
	  
	  fclose(fout);
	}
      printf("Output afterfit.par and ifuncDGW files\n");
      // Now get a spectrum of each pulsar (should do this pairwise using the same
      // frequency channels -- for now assume identical sampling)
      /*  for (p=0;p<*npsr;p++)
	  {
	  OMEGA0 = (double)(2.0*M_PI/(psr[p].ifuncT[n-1]-psr[p].ifuncT[0]));
	  sprintf(tstr,"%s.specDGW",psr[p].name);
	  fout = fopen(tstr,"w");
	  getSpectrum(&psr[p],px1,py_r1,py_i1,&nSpec,toffset);
	  for (i=0;i<nSpec;i++)
	  if(write_debug_files)fprintf(fout,"%g %g %g\n",px1[i],py_r1[i],py_i1[i]);
	  fclose(fout);
	  }*/
      
      // Now do each pair of pulsars
      fout = fopen("hdcurve.dat","w");
      
      for (p1 = 0;p1 < *npsr ;p1++)
	{
	  for (p2 = p1+1;p2 < *npsr;p2++)
	    {
	      printf("Processing: %s %s\n",psr[p1].name,psr[p2].name);
	      // Find overlap region
	      if (psr[p1].ifuncT[0] > psr[p2].ifuncT[0])
		startOverlap = psr[p1].ifuncT[0];
	      else
		startOverlap = psr[p2].ifuncT[0];
	      
	      if (psr[p1].ifuncT[psr[p1].ifuncN-1] > psr[p2].ifuncT[psr[p2].ifuncN-1])
		endOverlap = psr[p2].ifuncT[psr[p2].ifuncN-1];
	      else
		endOverlap = psr[p1].ifuncT[psr[p1].ifuncN-1];
	      
	      
	      toverlap = endOverlap-startOverlap;
	      if (toverlap > maxOverlap)
		{
		  OMEGA0 = (double)(2.0*M_PI/toverlap);
		  
		  printf("Processing pair: %s--%s\n",psr[p1].name,psr[p2].name);
		  printf("Overlap time = %g\n",toverlap);
		  
		  getSpectrum(&psr[p1],px1,py_r1,py_i1,&nSpec1,toffset,startOverlap,endOverlap,stepMJD,covarFuncFile);
		  printf("Got initial spectrum\n");
		  getSpectrum(&psr[p2],px2,py_r2,py_i2,&nSpec2,toffset,startOverlap,endOverlap,stepMJD,covarFuncFile);
		  printf("Got second spectrum\n");
		  
		  if (nSpec1 < nSpec2) nSpec = nSpec1;
		  else nSpec = nSpec2;
		  
		  
		  sprintf(tstr,"%s-%s.crossSpecDGW",psr[p1].name,psr[p2].name);
		  if (write_debug_files) {
		    logmsg("Writing '%s'",tstr);
		    fout2 = fopen(tstr,"w");
		  }
		  
		  // Form the cross power spectrum
		  for (i=0;i<nSpec;i++)
		    {
		      crossX[i] = px1[i]; 
		      crossY_r[i] = (py_r1[i]*py_r2[i]+py_i1[i]*py_i2[i]); // /(toverlap/365.25);
		      crossY_i[i] = (py_i1[i]*py_r2[i]-py_i2[i]*py_r1[i]); // /(toverlap/365.25);
		    }

		  angle[npair]  = (double)psrangle(psr[p1].param[param_raj].val[0],
						   psr[p1].param[param_decj].val[0],
						   psr[p2].param[param_raj].val[0],
						   psr[p2].param[param_decj].val[0]);
		  for (i=0;i<nSpec;i++)
		    {
		      if (write_debug_files) fprintf(fout2,"%g %g %g %g %g %g %g %g %g %g %g %g\n",crossX[i],crossY_r[i],crossY_i[i],px1[i],py_r1[i]*py_r1[i]+py_i1[i]*py_i1[i],px2[i],py_r2[i]*py_r2[i]+py_i2[i]*py_i2[i], py_r1[i],py_i1[i],py_r2[i],py_i2[i],angle[npair]);
		    }
		  if (write_debug_files) fclose(fout2);

		  cx = (1.0-cos(angle[npair]*M_PI/180.0))/2.0;
		  zeta = 3.0/2.0*(cx*log(cx)-cx/6.0+1.0/3.0);
		  
		  sprintf(fname,"%s.model",psr[p1].name);
		  if (!(fin = fopen(fname,"r")))
		    {
		      printf("Unable to open file >%s<\n",fname);
		      exit(1);
		    }
		  fscanf(fin,"%s %s",dummy,dummy);
		  fscanf(fin,"%s %lf",dummy,&alpha1);
		  fscanf(fin,"%s %lf",dummy,&fc1);
		  fscanf(fin,"%s %lf",dummy,&amp1);
		  if (fscanf(fin,"%s %lf",dummy,&pw1)!=2)
		    {
		      printf("Unable to read whitenoise level >%s<\n",fname);
		      exit(1);
		    }
		  fclose(fin); 

		  sprintf(fname,"%s.model",psr[p2].name);
		  fin = fopen(fname,"r");
		  fscanf(fin,"%s %s",dummy,dummy);
		  fscanf(fin,"%s %lf",dummy,&alpha2);
		  fscanf(fin,"%s %lf",dummy,&fc2);
		  fscanf(fin,"%s %lf",dummy,&amp2);
		  fscanf(fin,"%s %lf",dummy,&pw2);
		  fclose(fin);
		  
		  if (useRed==1)
		    {
		      sprintf(fname,"%s.rednoise",psr[p1].name);
		      fin = fopen(fname,"r");
		      fscanf(fin,"%s %lf",dummy,&red_pwr1);
		      fscanf(fin,"%s %lf",dummy,&red_alpha1);
		      fclose(fin);
		      
		      sprintf(fname,"%s.rednoise",psr[p2].name);
		      fin = fopen(fname,"r");
		      fscanf(fin,"%s %lf",dummy,&red_pwr2);
		      fscanf(fin,"%s %lf",dummy,&red_alpha2);
		      fclose(fin);
		    }
		  
		  double pwe=0;
		  for (i=nSpec-10;i<nSpec;i++){
		    pwe+=py_r1[i]*py_r1[i]+py_i1[i]*py_i1[i];
		  }
		  pwe/=10.0;
		  logmsg("pw1=%lg pwE=%lg",pw1,pwe);
		  //		  pw1*=3.0;
		  pwe=0;
		  for (i=nSpec-10;i<nSpec;i++){
		    pwe+=py_r2[i]*py_r2[i]+py_i2[i]*py_i2[i];
		  }
		  pwe/=10.0;
		  logmsg("pw2=%lg pwe=%lg",pw2,pwe);
		  //		  pw2*=3.0;
		  //
		  nSpecNdof=nSpec;
		  specStart=0;
		  if(useRed==1)nSpecNdof=0;
		  for (i=0;i<nSpec;i++)
		    {
		      pg = guessGWamp*guessGWamp*pow((double)(crossX[i]*365.25),alpha_res)/12.0/M_PI/M_PI;
		      modelPwr1 = pg + pw1;
		      modelPwr2 = pg + pw2;
		      if (useRed==1)
			{
			  double pr1,pr2;
			  pr1 = red_pwr1*pow((double)(crossX[i]*365.25),red_alpha1);
			  pr2 = red_pwr2*pow((double)(crossX[i]*365.25),red_alpha2);
			  modelPwr1 += pr1;
			  modelPwr2 += pr2;
			  if (nSpecNdof==0){
			    if(pg > pr1+pw1 && pg > pr2+pw2){
			      nSpecNdof++;
			    } else {
			      specStart++;
			      }
			  } else {
			    if(pg > pr1 && pg > pr2){
			      nSpecNdof++;
			    } else{
			      break;
			    }
			  }
			}
		  //
		  // Should check this 0.5+zeta^2
		  //
		  //crossPowerErr[i] = sqrt((0.5+zeta*zeta)*modelPwr1*modelPwr2);
		      crossPowerErr[i] = sqrt((0.5*(1+zeta*zeta))*modelPwr1*modelPwr2);
		      //crossPowerErr[i] = sqrt(modelPwr1*modelPwr2);
		      //		  printf("In here with %s-%s %g %g %g %g %g %g\n",psr[p1].name,psr[p2].name,pg,pw1,pw2,guessGWamp,crossX[i],crossPowerErr[i]);
		    }
		  
		  nSpec = nSpecNdof+specStart;
		  
		  if (nSpecNdof > 0)
		    {
		      sum1=sum2=sum3=sum4=0.0L;
		      weight1 = -alpha_res;
		      sprintf(tstr,"%s-%s.weightDGW",psr[p1].name,psr[p2].name);
		      if(write_debug_files)fout2 = fopen(tstr,"w");
		      
		      for (i=specStart;i<nSpec;i++)
			{
			  sum1 += (crossY_r[i]*pow((double)(i+1),-1.0*weight1)/pow(crossPowerErr[i],2));
			  sum2 += (pow((double)(i+1),-2.0*weight1)/pow(crossPowerErr[i],2));
			  sum3 += (crossY_i[i]*pow((double)(i+1),-1.0*weight1)/pow(crossPowerErr[i],2));
			  sum4 += (1.0/pow(crossPowerErr[i],2)/pow((double)(i+1),2.0*weight1));
			  if(write_debug_files)fprintf(fout2,"%g %g %g %g %g\n",crossX[i],crossPowerErr[i],(pow((double)(i+1),-1.0*weight1)/pow(crossPowerErr[i],2)),(pow((double)(i+1),-2.0*weight1)/pow(crossPowerErr[i],2)),(1.0/pow(crossPowerErr[i],2)/pow((double)(i+1),2.0*weight1)));
			}
		      if (write_debug_files) fclose(fout2);
		      printf("Complete weighting\n");
		      a2zeta[npair]    = 12.0*M_PI*M_PI/pow(toverlap/365.25,weight1)*(sum1/sum2); // Equation 9 in Yardley et al.
		      printf("Step2\n");
		      a2zeta_im[npair] = 12.0*M_PI*M_PI/pow(toverlap/365.25,weight1)*(sum3/sum2); // Equation 14??
		      printf("Step3\n");
		      a2zeta_e[npair]  = 12.0*M_PI*M_PI/pow(toverlap/365.25,weight1)/sqrt(sum4); // Equation 12??
		      printf("Step4\n");
		      printf("v1 = %g\n",angle[npair]);
		      printf("v2 = %g\n",a2zeta[npair]);
		      
		      fprintf(fout,"%g %g %g %s %s %g %g %g %g %d %d %d %s %s\n",angle[npair],a2zeta[npair],a2zeta_e[npair],psr[p1].name,psr[p2].name,(double)sum1,(double)sum2,(double)weight1,(double)toverlap,nSpec,nSpec1,nSpec2,timFile[p1],timFile[p2]);
		      printf("Step5\n");
		      npair++;
		      printf("Trying next pulsar pair\n");
		    }
		}
	      printf("Completed processing: %s %s\n",psr[p1].name,psr[p2].name);

	    }
	}
      fclose(fout);
    }
  else {
    FILE *fin;
    char temp[1024];

    npair = 0;
    fin = fopen(hdFile,"r");
    while (!feof(fin))
      {
	if (fscanf(fin,"%lf %lf %lf %s %s %s %s %s %s %s",&angle[npair],&a2zeta[npair],&a2zeta_e[npair],temp,temp,temp,temp,temp,temp,temp)==10)
	  npair++;
      }
    fclose(fin);

  }

  // Fit for the amplitude of the HD curve and get its error
  // Currently not including errors in this fit
  {
    double p[2],e[2];
    double **cvm,**cvm1;
    double chisq;
    double storeA2zeta[4096];
    int nf=1;
    float fx[180],fy[180];
    double amp;
    
    double p1,e1,p2,e2,p3,e3,c1;
    double p4,p5,p6,e4,e5,e6;
    double p7,e7,p8,e8;
    FILE *fout;
    double meanCovar=0;



    printf("Calculating the significance\n");
    cvm = (double**)malloc(sizeof(double *)*2);
    for (i=0;i<2;i++)
      cvm[i] = (double *)malloc(sizeof(double)*2);

    cvm1 = (double**)malloc(sizeof(double *)*1);
    for (i=0;i<1;i++)
      cvm1[i] = (double *)malloc(sizeof(double)*1);
    printf("Allocated memory %d %d\n",npair,nf);
    TKleastSquares_svd(angle,a2zeta,a2zeta_e,npair,p,e,nf,cvm1,&chisq,hdfunc,1);
    printf("Done the fit\n");

    amp = p[0];
    printf("A2 = %g +/- %g\n",p[0],e[0]);
    if (p[0] > 0)
      printf("A = %g\n",sqrt(p[0]));
    printf("Significance = %g\n",p[0]/e[0]);
    printf("Chisq = %g, red chisq = %g\n",chisq,chisq/(double)(npair-1));
    p1 = p[0];
    e1 = e[0];
    c1 = chisq/(double)(npair-1);
    p2 = p[0];
    e2 = e[0]*sqrt(c1); // Increase error by chisq
    printf("2] A2 = %g +/- %g\n",p2,e2);
    if (p2 > 0)
      printf("2] A = %g\n",sqrt(p2));
    printf("2] Significance = %g\n",p2/e2);

    TKleastSquares_svd(angle,a2zeta,a2zeta_e,npair,p,e,nf,cvm1,&chisq,hdfunc,0);
    amp = p[0];
    printf("3] A2 = %g +/- %g\n",p[0],e[0]);
    if (p[0] > 0)
      printf("3] A = %g\n",sqrt(p[0]));
    printf("3] Significance = %g\n",p[0]/e[0]);
    printf("3] Chisq = %g, red chisq = %g\n",chisq,chisq/(double)(npair-1));
    p3 = p[0];
    e3 = e[0];

    // Now fit for an offset as well as the hd curve

    TKleastSquares_svd(angle,a2zeta,a2zeta_e,npair,p,e,2,cvm,&chisq,hdfunc_offs,1);
    printf("4] A2 = %g +/- %g\n",p[0],e[0]);
    printf("4] Significance = %g\n",p[0]/e[0]);
    printf("4] offset = %g +/- %g\n",p[1],e[1]);
    printf("4] chisq = %g\n",chisq);
    printf("4] chisq reduced = %g\n",chisq/(double)(npair-2));
    p4 = p[0];
    e4 = e[0];

    // Now calculate using the mean subtracted HD curve
    {
      double ctheta,cx;

      GLOBAL_MEANSUB=0;
      for (i=0;i<npair;i++)
	{
	  ctheta = cos(angle[i]*M_PI/180.0);
	  cx = (1.0-ctheta)/2.0;
	  GLOBAL_MEANSUB += (cx*log(cx)-cx/6.0+1.0/3.0)*3.0/2.0;
	}
      GLOBAL_MEANSUB /= (double)npair;
    }
    printf("Global mean subtracted value from HD curve is %g\n",GLOBAL_MEANSUB);


    TKleastSquares_svd(angle,a2zeta,a2zeta_e,npair,p,e,nf,cvm,&chisq,hdfunc_meanSub,1);
    printf("5] A2 = %g +/- %g\n",p[0],e[0]);
    printf("5] Significance = %g\n",p[0]/e[0]);
    printf("5] offset = %g +/- %g\n",p[1],e[1]);
    printf("5] chisq = %g\n",chisq);
    printf("5] chisq reduced = %g\n",chisq/(double)(npair-2));
    p5 = p[0];
    e5 = e[0];


    // Now remove the mean from the covariance estimates
    printf("Removing mean from the measured covariance estimates\n");
    for (i=0;i<npair;i++)
      meanCovar += a2zeta[i];
    meanCovar /= (double)npair;
    printf(" .. mean value = %g\n",meanCovar);
    for (i=0;i<npair;i++)
      {
	storeA2zeta[i] = a2zeta[i];
	a2zeta[i] -= meanCovar;
      }
    TKleastSquares_svd(angle,a2zeta,a2zeta_e,npair,p,e,nf,cvm,&chisq,hdfunc_meanSub,1);
    printf("6] A2 = %g +/- %g\n",p[0],e[0]);
    printf("6] Significance = %g\n",p[0]/e[0]);
    printf("6] offset = %g +/- %g\n",p[1],e[1]);
    printf("6] chisq = %g\n",chisq);
    printf("6] chisq reduced = %g\n",chisq/(double)(npair-2));
    p6 = p[0];
    e6 = e[0];


    // Now fit simultaneously for the A^2 and the amplitude of a cosine
    for (i=0;i<npair;i++)
      a2zeta[i] = storeA2zeta[i];
    TKleastSquares_svd(angle,a2zeta,a2zeta_e,npair,p,e,2,cvm,&chisq,hdfunc_removeCosine,1);
    printf("7] Fit for A^2 and cosine: A2 = %g +/- %g\n",p[0],e[0]);
    printf("7] cosine amplitude = %g +/- %g\n",p[1],e[1]);
    p7 = p[0];
    e7 = e[0];

    // Now calculate the cosine component of the actual HD curve
    // Now calculate using the mean subtracted HD curve
    {
      double ctheta,cx;
      double cpx[npair];
      double cpy[npair];
      
      GLOBAL_COSVAL=0;
      for (i=0;i<npair;i++)
	{
	  cpx[i] = angle[i];

	  ctheta = cos(angle[i]*M_PI/180.0);
	  cx = (1.0-ctheta)/2.0;
	  cpy[i] = (cx*log(cx)-cx/6.0+1.0/3.0)*3.0/2.0;
	}
      TKleastSquares_svd(cpx,cpy,a2zeta_e,npair,p,e,1,cvm,&chisq,cosineFunc,0);
      GLOBAL_COSVAL = p[0];
      printf("The fit of a cosine to the HD curve gives amplitude = %g +/- %g\n",p[0],e[0]);
    }


    TKleastSquares_svd(angle,a2zeta,a2zeta_e,npair,p,e,1,cvm,&chisq,hdfunc_cosineSub,1);
    printf("8] A2 = %g +/- %g\n",p[0],e[0]);
    printf("8] Significance = %g\n",p[0]/e[0]);
    printf("8] chisq = %g\n",chisq);
    printf("8] chisq reduced = %g\n",chisq/(double)(npair-2));
    p8 = p[0];
    e8 = e[0];



    

    fout = fopen("result.dat","w");
    fprintf(fout,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",p1,e1,p1/e1,c1,p2,e2,p2/e2,p3,e3,p3/e3,p4,e4,p4/e4,p5,e5,p5/e5,p6,e6,p6/e6,p7,e7,p7/e7,p8,e8,p8/e8,chisq);
    fclose(fout);

    free(cvm1[0]);
    free(cvm1);
    for (i=0;i<2;i++)
      free(cvm[i]);
    free(cvm);
  }

  return 0;
}

void getSpectrum(pulsar *psr,double *px,double *py_r,double *py_i,int *nSpec,double toffset,double startOverlap,double endOverlap,double stepMJD,char *covarFuncFile)
{


   /* We make a fake pulsar to allow us to re-use the 'standard' tempo2
	* spectral analysis code
	*/
   int nobs=0;
   int i;
   longdouble peopoch=toffset;
   pulsar* fakepsr = (pulsar*)malloc(sizeof(pulsar));

   printf("In get spectrum\n");
   fakepsr->obsn = (observation*) malloc(sizeof(observation)*psr->ifuncN);
   fakepsr->ToAextraCovar = NULL;
   fakepsr->robust=1;
   strcpy(fakepsr->name,psr->name);

   *nSpec = (int)((((endOverlap-startOverlap)/(double)(stepMJD))-1)/2.0);
   printf("nSpec = %d\n",*nSpec);
   fakepsr->param[param_pepoch].val=&peopoch;
   for(i=0; i < psr->ifuncN;i++){
      if (psr->ifuncT[i] >= startOverlap && psr->ifuncT[i] <= endOverlap){
		 fakepsr->obsn[nobs].sat=psr->ifuncT[i];
		 fakepsr->obsn[nobs].bat=psr->ifuncT[i];
		 fakepsr->obsn[nobs].bbat=psr->ifuncT[i];
		 fakepsr->obsn[nobs].residual=psr->ifuncV[i];
		 fakepsr->obsn[nobs].toaErr=psr->ifuncE[i]*1e6; // convert to us
		 nobs++;
	  }
   }
   fakepsr->nobs=nobs;
   printf("Set up spectrum %d\n",fakepsr->nobs);
   printf("nSpec = %d\n",*nSpec);
/* From here on in we copy the cholSpectra plugin
 */
   double **uinv;
   FILE *fin;
   char fname[128];
   int ndays=0;
   double resx[fakepsr->nobs],resy[fakepsr->nobs],rese[fakepsr->nobs];
   int ip[fakepsr->nobs];
   FILE *fout;

   //  printf("Calculating the spectrum\n");
   uinv = malloc_uinv(fakepsr->nobs);


   logmsg("Nspec=%d, nobs=%d",*nSpec,fakepsr->nobs);
   for (i=0;i<fakepsr->nobs;i++)
   {
	  resx[i] = (double)(fakepsr->obsn[i].sat-toffset);
	  resy[i] = (double)(fakepsr->obsn[i].residual);
	  rese[i] = (double)(fakepsr->obsn[i].toaErr*1.0e-6);
	  ip[i]=i;
	  printf("Handing Cholesky: %g %g %g %d\n",resx[i],resy[i],rese[i],ip[i]);
   }
   logmsg("Get Cholesky 'uinv' matrix from '%s'",covarFuncFile);
   getCholeskyMatrix(uinv,covarFuncFile,fakepsr,resx,resy,rese,fakepsr->nobs,0,ip);
   logdbg("Got uinv, now compute spectrum.");

   // Must calculate uinv for the pulsar
   printf("calculating spectra\n");
   printf("fake_psr.robust = %d\n",fakepsr->robust);
   calcSpectra_ri(uinv,resx,resy,fakepsr->nobs,px,py_r,py_i,*nSpec,psr);
   printf("complete calculating spectra\n");
   // Free uinv
   free_uinv(uinv);
   free(fakepsr->obsn);
   free(fakepsr);
}


double psrangle(double centre_long,double centre_lat,double psr_long,double psr_lat)
{
  double dlon,dlat,a,c;
  double deg2rad = M_PI/180.0;
  
  /* Apply the Haversine formula */
  dlon = (psr_long - centre_long);
  dlat = (psr_lat  - centre_lat);
  a = pow(sin(dlat/2.0),2) + cos(centre_lat) * 
    cos(psr_lat)*pow(sin(dlon/2.0),2);
  if (a==1)
    c = M_PI/deg2rad;
  else
    c = 2.0 * atan2(sqrt(a),sqrt(1.0-a))/deg2rad;  
  return c;
}



void hdfunc_meanSub(double x,double *p,int ma)
{
  double ctheta,cx;

  ctheta = cos(x*M_PI/180.0);
  cx = (1.0-ctheta)/2.0;
  p[0] = (cx*log(cx)-cx/6.0+1.0/3.0)*3.0/2.0-GLOBAL_MEANSUB;

}


void hdfunc_cosineSub(double x,double *p,int ma)
{
  double ctheta,cx;

  ctheta = cos(x*M_PI/180.0);
  cx = (1.0-ctheta)/2.0;
  p[0] = (cx*log(cx)-cx/6.0+1.0/3.0)*3.0/2.0-GLOBAL_COSVAL;

}


void hdfunc(double x,double *p,int ma)
{
  double ctheta,cx;

  ctheta = cos(x*M_PI/180.0);
  cx = (1.0-ctheta)/2.0;
  p[0] = (cx*log(cx)-cx/6.0+1.0/3.0)*3.0/2.0;

}

void hdfunc_offs(double x,double *p,int ma)
{
  double ctheta,cx;

  ctheta = cos(x*M_PI/180.0);
  cx = (1.0-ctheta)/2.0;
  p[0] = (cx*log(cx)-cx/6.0+1.0/3.0)*3.0/2.0;
  p[1] = 1.0;
}

void hdfunc_removeCosine(double x,double *p,int ma)
{
  double ctheta,cx;

  ctheta = cos(x*M_PI/180.0);
  cx = (1.0-ctheta)/2.0;
  p[0] = (cx*log(cx)-cx/6.0+1.0/3.0)*3.0/2.0;
  p[1] = ctheta;
}

void cosineFunc(double x,double *p,int ma)
{
  p[0] = cos(x*M_PI/180.0);
}








