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

void getSpectrum(pulsar *psr,double *px,double *py_r,double *py_i,int *nSpec,double toffset,double startOverlap,double endOverlap,double stepMJD,char *dofFile);
void formCholeskyMatrixPlugin(double *c,double *resx,double *resy,double *rese,int np,double **uinv);
int calcSpectra_plugin(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY_R,double *specY_I,int nfit);
double psrangle(double centre_long,double centre_lat,double psr_long,double psr_lat);
void hdfunc(double x,double *p,int ma);

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
  char dofFile[128]="";
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
  double maxOverlap = 500; // Only consider pairs where the overlap is larger than this value
  char dcmFile[128];
  char covarFuncFile[128];
  char newname[128];

  int useRed=0;
  int nSpecNdof;
  int specStart;

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
      else if (strcmp(argv[i],"-dof")==0)
	strcpy(dofFile,argv[++i]);
      else if (strcasecmp(argv[i],"-usered")==0)
	useRed=1;
      else if (strcasecmp(argv[i],"-fast")==0)
	write_debug_files=false;
    }

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
  printf("Starting 2 >%s<\n",dofFile);
  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr);         /* Load the arrival times    */
  printf("Here: %d\n",psr[0].nobs);
  preProcess(psr,*npsr,argc,argv);
  printf("Here 2: %d\n",psr[0].nobs);
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
	      getSpectrum(&psr[p1],px1,py_r1,py_i1,&nSpec1,toffset,startOverlap,endOverlap,stepMJD,dofFile);
	      getSpectrum(&psr[p2],px2,py_r2,py_i2,&nSpec2,toffset,startOverlap,endOverlap,stepMJD,dofFile);


	      if (nSpec1 < nSpec2) nSpec = nSpec1;
	      else nSpec = nSpec2;


	      sprintf(tstr,"%s-%s.crossSpecDGW",psr[p1].name,psr[p2].name);
	      if (write_debug_files) fout2 = fopen(tstr,"w");
	      
	      // Form the cross power spectrum
	      for (i=0;i<nSpec;i++)
		{
		  crossX[i] = px1[i]; 
		  crossY_r[i] = (py_r1[i]*py_r2[i]+py_i1[i]*py_i2[i]); // /(toverlap/365.25);
		  crossY_i[i] = (py_i1[i]*py_r2[i]-py_i2[i]*py_r1[i]); // /(toverlap/365.25);
		  if (write_debug_files) fprintf(fout2,"%g %g %g %g %g %g %g\n",crossX[i],crossY_r[i],crossY_i[i],px1[i],py_r1[i]*py_r1[i]+py_i1[i]*py_i1[i],px2[i],py_r2[i]*py_r2[i]+py_i2[i]*py_i2[i]);
		}
	      if (write_debug_files) fclose(fout2);
	      angle[npair]  = (double)psrangle(psr[p1].param[param_raj].val[0],
					       psr[p1].param[param_decj].val[0],
					       psr[p2].param[param_raj].val[0],
					       psr[p2].param[param_decj].val[0]);
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
		  //	      crossPowerErr[i] = sqrt((0.5+zeta*zeta)*modelPwr1*modelPwr2);
		  crossPowerErr[i] = sqrt(modelPwr1*modelPwr2);
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
		      if(write_debug_files)fprintf(fout2,"%g %g %g %g\n",crossPowerErr[i],(pow((double)(i+1),-1.0*weight1)/pow(crossPowerErr[i],2)),(pow((double)(i+1),-2.0*weight1)/pow(crossPowerErr[i],2)),(1.0/pow(crossPowerErr[i],2)/pow((double)(i+1),2.0*weight1)));
		    }
		  if (write_debug_files) fclose(fout2);
		  printf("Complete weighting\n");
		  a2zeta[npair]    = 12.0*M_PI*M_PI/pow(toverlap/365.25,weight1)*(sum1/sum2); // Equation 9 in Yardley et al.
		  printf("Step2\n");
		  a2zeta_im[npair] = 12.0*M_PI*M_PI/pow(toverlap/365.25,weight1)*(sum3/sum2);
		  printf("Step3\n");
		  a2zeta_e[npair]  = 12.0*M_PI*M_PI/pow(toverlap/365.25,weight1)/sqrt(sum4);
		  printf("Step4\n");
		  printf("v1 = %g\n",angle[npair]);
		  printf("v2 = %g\n",a2zeta[npair]);

		  fprintf(fout,"%g %g %g %s %s %g %g %g %g %d %d %d %s %s\n",angle[npair],a2zeta[npair],a2zeta_e[npair],psr[p1].name,psr[p2].name,(double)sum1,(double)sum2,(double)weight1,(double)toverlap,nSpec,nSpec1,nSpec2,timFile[p1],timFile[p2]);
		  printf("Step5\n");
		  npair++;
		  printf("Trying next pulsar pair\n");
		}
	    }
	}
    }
  fclose(fout);
  // Fit for the amplitude of the HD curve and get its error
  // Currently not including errors in this fit
  {
    double p[2],e[2];
    double **cvm;
    double chisq;
    int nf=1;
    float fx[180],fy[180];
    double amp;
    double p1,e1,p2,e2,p3,e3,c1;
    FILE *fout;


    cvm = (double**)malloc(sizeof(double *)*nf);
    for (i=0;i<nf;i++)
      cvm[i] = (double *)malloc(sizeof(double)*nf);

    TKleastSquares_svd(angle,a2zeta,a2zeta_e,npair,p,e,nf,cvm,&chisq,hdfunc,1);
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

    TKleastSquares_svd(angle,a2zeta,a2zeta_e,npair,p,e,nf,cvm,&chisq,hdfunc,0);
    amp = p[0];
    printf("3] A2 = %g +/- %g\n",p[0],e[0]);
    if (p[0] > 0)
      printf("3] A = %g\n",sqrt(p[0]));
    printf("3] Significance = %g\n",p[0]/e[0]);
    printf("3] Chisq = %g, red chisq = %g\n",chisq,chisq/(double)(npair-1));
    p3 = p[0];
    e3 = e[0];

    fout = fopen("result.dat","w");
    fprintf(fout,"%g %g %g %g %g %g %g %g %g %g\n",p1,e1,p1/e1,c1,p2,e2,p2/e2,p3,e3,p3/e3,chisq);
    fclose(fout);

    for (i=0;i<nf;i++)
      free(cvm[i]);
    free(cvm);
  }

  return 0;
}

void getSpectrum(pulsar *psr,double *px,double *py_r,double *py_i,int *nSpec,double toffset,double startOverlap,double endOverlap,double stepMJD,char *dofFile)
{
  double xv[psr->ifuncN],yv[psr->ifuncN],ev[psr->ifuncN];
  double **uinv;
  int    n = psr->ifuncN;
  int    i;
  FILE   *fin;
  char   fname[128];
  double escaleFactor,tt;
  double covarFunc[100000];
  int    ndays;
  int n2;

  // Have 8 originally
  //  *nSpec = 20;
    *nSpec = (int)((((endOverlap-startOverlap)/(double)(stepMJD))-1)/2.0);
  //  *nSpec = 10;
  //  *nSpec=20;

  //  *nSpec = 5;

  printf("Number of spectral channels = %d\n",*nSpec);
  uinv = (double **)malloc(sizeof(double *)*n);
  for (i=0;i<n;i++)
    uinv[i] = (double *)malloc(sizeof(double)*n);

  // Read in the covariance function
  sprintf(fname,"covarFunc.dat_%s",psr->name);
  fin = fopen(fname,"r");
  fscanf(fin,"%lf",&escaleFactor);
  ndays = 0;
  while (!feof(fin))
    {
      if (fscanf(fin,"%lf",&covarFunc[ndays])==1)
	ndays++;
    }
  fclose(fin);
  tt = covarFunc[ndays-1];
  printf("Size of covarFunc = %d, %g %g %g\n",ndays,endOverlap,startOverlap,(endOverlap-startOverlap)-ndays+1);

  for (i=ndays;i<(endOverlap-startOverlap)+1;i++) // Note that the IFUNC gridding can extend outside of the standard range
    {
      covarFunc[ndays++] = tt;
    }

  n2=0;
  for (i=0;i<n;i++)
    {
      if (psr->ifuncT[i] >= startOverlap && psr->ifuncT[i] <= endOverlap)
	{
	  xv[n2] = psr->ifuncT[i]-toffset;
	  yv[n2] = psr->ifuncV[i];
	  ev[n2] = psr->ifuncE[i];
	  n2++;
	}
    }
  printf("Forming Cholesky matrix with n = %d, n2 = %d, name = %s, start = %g, end = %g\n",n,n2,psr->name,startOverlap,endOverlap);
  formCholeskyMatrixPlugin(covarFunc,xv,yv,ev,n2,uinv);
  calcSpectra_plugin(uinv,xv,yv,n2,px,py_r,py_i,*nSpec);
  printf("dof file = %s\n",dofFile);
  if (strlen(dofFile) > 0) // Have an NDOF file
    {
      FILE *fin;
      char psrs[128];
      double ndofV;
      int found=0;
      if (!(fin = fopen(dofFile,"r")))
	{
	  printf("Unable to read file %s\n",dofFile);
	  exit(1);
	}
      while (!feof(fin))
	{
	  if (fscanf(fin,"%s %lf",psrs,&ndofV)==2)
	    {
	      if (strcmp(psrs,psr->name)==0)
		{
		  found=1;
		  printf("%s changing ndof from %d to %d\n",psrs,*nSpec,(int)(ndofV+0.5));
		  *nSpec = (int)(ndofV+0.5);
		}
	    }
	}
      fclose(fin);
      if (found==0)
	{
	  printf("Unable to find pulsar >%s< in file >%s<\n",psr->name,dofFile);
	  exit(1);
	}
	
    }
  for (i=0;i<n;i++)
    free(uinv[i]);
  free(uinv);
}

void formCholeskyMatrixPlugin(double *c,double *resx,double *resy,double *rese,int np,double **uinv)
{
  double **m,**u,sum;
  double *cholp;
  int i,j,k,ix,iy;
  double t0,cint,t;
  int t1,t2;
  int debug=0;

  //  printf("Getting the covariance matrix in doFit\n");
  m = (double **)malloc(sizeof(double *)*(np+1));
  u= (double **)malloc(sizeof(double *)*(np+1));
  cholp  = (double *)malloc(sizeof(double)*(np+1));  // Was ndays
  
  for (i=0;i<np+1;i++)
    {
      m[i] = (double *)malloc(sizeof(double)*(np+1));
      u[i] = (double *)malloc(sizeof(double)*(np+1));
    }
  //  printf("Allocated memory\n");
  for (ix=0;ix<np;ix++)
    {
      for (iy=0;iy<np;iy++)
	m[ix][iy] = fabs(resx[ix]-resx[iy]);
    }
  if (debug==1)
    {
      printf("First m = \n");
      for (i=0;i<5;i++)
	{ 
	  for (j=0;j<5;j++) printf("%10g ",m[i][j]); 
	  printf("\n");
	}

    }
  // Insert the covariance which depends only on the time difference.
  // Linearly interpolate between elements on the covariance function because
  // valid covariance matrix must have decreasing off diagonal elements.
  //  printf("Inserting into the covariance matrix\n");
  for (ix=0;ix<np;ix++)
    {
      for (iy=0;iy<np;iy++)
	{
	  t0 = m[ix][iy];
	  t1 = (int)floor(t0);
	  t2 = t1+1;
	  t  = t0-t1;
	  cint = c[t1]*(1-t)+c[t2]*t; // Linear interpolation
	  m[ix][iy] = cint;
	}
    }
  //  printf("Multiplying by errors\n");
  for (ix=0;ix<np;ix++)
    m[ix][ix]+=rese[ix]*rese[ix];
  if (debug==1)
    {
      printf("m = \n\n");
      for (i=0;i<5;i++)
	{ 
	  for (j=0;j<5;j++) printf("%10g ",m[i][j]); 
	  printf("\n");
	}
    }

  // Do the Cholesky
  //  printf("Cholesky decomposition\n");
  TKcholDecomposition(m,np,cholp);
  //  printf("Complete cholesky decomposition\n");
  // Now calculate uinv
  for (i=0;i<np;i++)
    {
      //      printf("i = %d ... %d\n",i,np);
      m[i][i] = 1.0/cholp[i];
      //      printf("s1\n");
      uinv[i][i] = m[i][i];
      //      printf("s2\n");
      for (j=0;j<i;j++)
      	uinv[i][j] = 0.0;
      //      printf("s3\n");
      for (j=i+1;j<np;j++)
	{
	  sum=0.0;
	  for (k=i;k<j;k++) sum-=m[j][k]*m[k][i];
	  m[j][i]=sum/cholp[j];
	  uinv[i][j] = m[j][i];
	}
      //      printf("s4\n");
    } 
  //  printf("Complete cholesky\n");
  if (debug==1)
    {
            printf("uinv = \n\n");
      for (i=0;i<5;i++)
	{ 
	  for (j=0;j<5;j++) printf("%10g ",uinv[i][j]); 
	  printf("\n");
	}
      printf("Completed inverting the matrix\n");
    }



  // Should free memory not required
  // (note: not freeing uinv)

  for (i=0;i<np+1;i++)
    {
      free(m[i]);
      free(u[i]);
    }
  free(m);
  free(u);
  free(cholp);
}

int calcSpectra_plugin(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY_R,double *specY_I,int nfit)
{
  int i,j,k;
  //  int nfit=nres/2-1;
  int nSpec;

  if (nfit < 0)
    nfit=nres/2-1;

  double sig[nres];
  double **newUinv;
  double **cvm;
  double chisq;
  pulsar *psr;
  int ip[nres];
  int nf = 3;
  double param[nf],error[nf];
  double v[nf];

  cvm = (double **)malloc(sizeof(double *)*nf);
  for (i=0;i<nf;i++)
    cvm[i] = (double *)malloc(sizeof(double)*nf);

  // Should fit independently to all frequencies
  for (i=0;i<nres;i++)
    {
      sig[i] = 1.0; // The errors are built into the uinv matrix
      ip[i] = 0;
    }
 
  for (k=0;k<nfit;k++)
    {
      //      printf("%5.2g\%\r",(double)k/(double)nfit*100.0);
      //      fflush(stdout);
      GLOBAL_OMEGA = OMEGA0*(k+1);
      //      TKleastSquares_svd_psr_dcm(resx,resy,sig,nres,param,error,3,cvm,&chisq,fitMeanSineFunc,0,psr,1.0e-40,ip,uinv);
      TKleastSquares_svd_psr_dcm(resx,resy,sig,nres,param,error,2,cvm,&chisq,fitCosSineFunc,0,psr,1.0e-40,ip,uinv);
      v[k] = (resx[nres-1]-resx[0])/365.25/2.0/pow(365.25*86400.0,2); 
      specX[k] = GLOBAL_OMEGA/2.0/M_PI;
      //            specY_R[k] = sqrt(v[k])*param[1];
      //            specY_I[k] = sqrt(v[k])*param[2];

      specY_R[k] = sqrt(v[k])*param[0];
      specY_I[k] = sqrt(v[k])*param[1];

      //specY_R[k] = (2.0/nres)*param[1]; //sqrt(v[k])*param[1];
      //      specY_I[k] = (2.0/nres)*param[2]; //sqrt(v[k])*param[2];
    }

  for (i=0;i<nf;i++)
    free(cvm[i]);
  free(cvm);
  //  printf("Complete spectra\n");
  return nfit;
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

void hdfunc(double x,double *p,int ma)
{
  double ctheta,cx;

  ctheta = cos(x*M_PI/180.0);
  cx = (1.0-ctheta)/2.0;
  p[0] = (cx*log(cx)-cx/6.0+1.0/3.0)*3.0/2.0;

}
