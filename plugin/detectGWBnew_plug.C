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
#include "TKfit.h"
#include "constraints.h"

using namespace std;

bool write_debug_files=true;
bool write_python_files=false;
double OMEGA0=0;

char notim=0;

double getSpectrum(pulsar *psr,double *px,double *py_r,double *py_i,int *nSpec,double toffset,double startOverlap,double endOverlap,double stepMJD,char* covarFuncFile,double t);
void formCholeskyMatrixPlugin(double *c,double *resx,double *resy,double *rese,int np,double **uinv);
int calcSpectra_plugin(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY_R,double *specY_I,int nfit);
double psrangle(double centre_long,double centre_lat,double psr_long,double psr_lat);
void hdfunc(double x,double *p,int ma);
int offsetToCM(pulsar* psr);



void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
   char parFile[MAX_PSR][MAX_FILELEN];
   char timFile[MAX_PSR][MAX_FILELEN];
   int i,p,n,p1,p2,j;
   double globalParameter;
   double startMJD,endMJD,stepMJD,t;
   char tstr[128]="";
   FILE *fout,*fin,*fout2,*pyfile;
   int  nSpec,nSpec1,nSpec2;
   double px1[4096],py_r1[4096],py_i1[4096];
   double px2[4096],py_r2[4096],py_i2[4096];
   double crossX[4096],crossY_r[4096],crossY_i[4096];
   double toverlap;
   double whitePSD[MAX_PSR];

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

   double dmStep=1e100;
   double toffset;
   double maxOverlap = 300; // Only consider pairs where the overlap is larger than this value
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
	  else if (strcmp(argv[i],"-notim")==0)
		 notim=1;
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
	  else if (strcasecmp(argv[i],"-dm")==0)
		 dmStep=stepMJD;
	  else if (strcasecmp(argv[i],"-fast")==0)
		 write_debug_files=false;
	  else if (strcasecmp(argv[i],"-py")==0)
		 write_python_files=true;

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
   readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
   readTimfile(psr,timFile,*npsr);         /* Load the arrival times    */
   printf("Here: %d\n",psr[0].nobs);
   preProcess(psr,*npsr,argc,argv);

   formBatsAll(psr,*npsr);             /* Form the barycentric arrival times */
   formResiduals(psr,*npsr,1);         /* Form the residuals                 */
   if (!notim){
	  printf("Here 2: %d\n",psr[0].nobs);
	  // Create some IFUNCs that will represent the low frequency signal
	  // ensure that these IFUNCs are constrained not to include 
	  // the quadratic etc.
	  for (p=0;p<*npsr;p++) {
		 n=0;

		 psr[p].param[param_dmmodel].fitFlag[0] = 1;
		 psr[p].param[param_dmmodel].paramSet[0] = 1;
		 psr[p].param[param_dmmodel].val[0] = 2;
		 psr[p].param[param_dmmodel].nLinkTo=0;
		 psr[p].param[param_dm].nLinkFrom=0;
		 psr[p].param[param_dmmodel].linkTo[(psr[p].param[param_dmmodel].nLinkTo)++] = param_dm;
		 psr[p].param[param_dm].linkFrom[(psr[p].param[param_dm].nLinkFrom)++]=param_dmmodel;

		 autosetDMCM(psr+p,dmStep,stepMJD,startMJD,endMJD,false);
		 psr[p].nconstraints = 0;
		 psr[p].fitMode=1;

		 psr[p].constraints[psr[p].nconstraints++] = constraint_dmmodel_mean;
		 psr[p].constraints[psr[p].nconstraints++] = constraint_dmmodel_cw_0;
		 //psr[p].constraints[psr[p].nconstraints++] = constraint_dmmodel_cw_1;
		 //psr[p].constraints[psr[p].nconstraints++] = constraint_dmmodel_cw_2;

		 sprintf(newname,"%s.new.par",psr[p].name);
		 if(write_debug_files)textOutput(psr+p,1,0,0,0,1,newname);
	  }


	  //
	  // Do the standard pulsar fits and get record the IFUNCs
	  //    
	  printf("Currently npts = %d\n",psr[0].nobs);
	  //doFitDCM(psr,"NULL",covarFuncFile,*npsr,0);       /* Do the fitting     */
	  doFit(psr,*npsr,0);       /* Do the fitting     */
	  //   if (write_debug_files){
	  //	  formBatsAll(psr,*npsr);             /* Form the barycentric arrival times */
	  //	  formResiduals(psr,*npsr,1);         /* Form the residuals                 */
	  //	  textOutput(psr,*npsr,globalParameter,0,0,0,tstr);  /* Display the output */
	  //  }

	  for (p=0;p<*npsr;p++)
	  {
		 if(write_debug_files){
			sprintf(newname,"%s.afterfit.par",psr[p].name);
			textOutput(psr+p,1,0,0,0,1,newname);
		 }
	  }
   } else {
	  printf("SKIP FITTING!!!! Using already existing Common Mode\n");
   }

   for (p=0;p<*npsr;p++){
	  double sss=0;
	  double eee;
	  sprintf(tstr,"%s.ifuncDGW",psr[p].name);
	  fout = fopen(tstr,"w");
	  for (i=0;i<psr[p].dmoffsCMnum;i++){
		 fprintf(fout,"%.2f %.10g %.10g\n",psr[p].dmoffsCM_mjd[i],psr[p].dmoffsCM[i],psr[p].dmoffsCM_error[i]);
		 eee=psr[p].dmoffsCM_error[i]/86400.0/365.25;
		 sss+=1.0/(eee*eee);
	  }
	  double tobs=(psr[p].dmoffsCM_mjd[psr[p].dmoffsCMnum-1]-psr[p].dmoffsCM_mjd[0])/365.25;
	  whitePSD[p]=2*tobs/sss;
	  logmsg("White %s = %g",psr[p].name,whitePSD[p]);

	  fclose(fout);
   }
   // Now get a spectrum of each pulsar (should do this pairwise using the same
   // frequency channels -- for now assume identical sampling)
   /*  for (p=0;p<*npsr;p++)
	   {
	   OMEGA0 = (double)(2.0*M_PI/(psr[p].dmoffsCM_mjd[n-1]-psr[p].dmoffsCM_mjd[0]));
	   sprintf(tstr,"%s.specDGW",psr[p].name);
	   fout = fopen(tstr,"w");
	   getSpectrum(&psr[p],px1,py_r1,py_i1,&nSpec,toffset);
	   for (i=0;i<nSpec;i++)
	   if(write_debug_files)fprintf(fout,"%g %g %g\n",px1[i],py_r1[i],py_i1[i]);
	   fclose(fout);
	   }*/

   // Now do each pair of pulsars
   fout = fopen("hdcurve.dat","w");
   if(write_python_files){
	  pyfile=fopen("GW.sum","w");
   }

   for (p1 = 0;p1 < *npsr ;p1++)
   {
	  for (p2 = p1+1;p2 < *npsr;p2++)
	  {
		 startOverlap=0;
		 endOverlap=0;
		 // Find overlap region
		 for (i=0; i<psr[p1].dmoffsCMnum && startOverlap==0; i++){
			for (j=0; j<psr[p2].dmoffsCMnum; j++){
			   if(fabs(psr[p1].dmoffsCM_mjd[i]-psr[p2].dmoffsCM_mjd[j])<0.5){
				  startOverlap=psr[p1].dmoffsCM_mjd[i];
				  break;
			   }
			}
		 }
		 for (i=psr[p1].dmoffsCMnum-1 ; i>=0 && endOverlap ==0; i--){
			for (j=psr[p2].dmoffsCMnum-1 ; j>=0 ; j--){
			   if(fabs(psr[p1].dmoffsCM_mjd[i]-psr[p2].dmoffsCM_mjd[j])<0.5){
				  endOverlap=psr[p1].dmoffsCM_mjd[i];
				  break;
			   }
			}
		 }

		 /* 
			if (psr[p1].dmoffsCM_mjd[0] > psr[p2].dmoffsCM_mjd[0])
			startOverlap = psr[p1].dmoffsCM_mjd[0];
			else
			startOverlap = psr[p2].dmoffsCM_mjd[0];

			if (psr[p1].dmoffsCM_mjd[psr[p1].dmoffsCMnum-1] > psr[p2].dmoffsCM_mjd[psr[p2].dmoffsCMnum-1])
			endOverlap = psr[p2].dmoffsCM_mjd[psr[p2].dmoffsCMnum-1];
			else
			endOverlap = psr[p1].dmoffsCM_mjd[psr[p1].dmoffsCMnum-1];
			*/

		 toverlap = endOverlap-startOverlap;
		 logmsg("%s %s %.1lf -> %.1lf overlap=%.1lf days (%.1lf yr)",psr[p1].name,psr[p2].name,startOverlap,endOverlap,toverlap,toverlap/365.25);
		 if (toverlap > maxOverlap)
		 {
			OMEGA0 = (double)(2.0*M_PI/toverlap);

			printf("Processing pair: %s--%s\n",psr[p1].name,psr[p2].name);
			nSpec1=0; // autoselect
			nSpec2=0; // autoselect
//			whitePSD[p1]= 
			   getSpectrum(&psr[p1],px1,py_r1,py_i1,&nSpec1,toffset,startOverlap,endOverlap,stepMJD,covarFuncFile,toverlap);
//			whitePSD[p2]=
			   getSpectrum(&psr[p2],px2,py_r2,py_i2,&nSpec2,toffset,startOverlap,endOverlap,stepMJD,covarFuncFile,toverlap);



			if (nSpec1 < nSpec2) nSpec = nSpec1;
			else nSpec = nSpec2;


			sprintf(tstr,"%s-%s.crossSpecDGW",psr[p1].name,psr[p2].name);
			if (write_debug_files) {
			   logmsg("Writing '%s'",tstr);
			   fout2 = fopen(tstr,"w");
			}
			angle[npair]  = (double)psrangle(psr[p1].param[param_raj].val[0],
				  psr[p1].param[param_decj].val[0],
				  psr[p2].param[param_raj].val[0],
				  psr[p2].param[param_decj].val[0]);

			if(write_python_files)fprintf(pyfile,"# % 14s % 14s % 8.3f % 12.4g % 12.4g\n",psr[p1].name,psr[p2].name,angle[npair],whitePSD[p1],whitePSD[p2]);

			// Form the cross power spectrum
			for (i=0;i<nSpec;i++)
			{
			   crossX[i] = px1[i]; 
			   crossY_r[i] = (py_r1[i]*py_r2[i]+py_i1[i]*py_i2[i]); // /(toverlap/365.25);
			   crossY_i[i] = (py_i1[i]*py_r2[i]-py_i2[i]*py_r1[i]); // /(toverlap/365.25);
			   if (write_debug_files) fprintf(fout2,"%g %g %g %g %g %g %g %g %g %g %g\n",crossX[i],crossY_r[i],crossY_i[i],px1[i],py_r1[i]*py_r1[i]+py_i1[i]*py_i1[i],px2[i],py_r2[i]*py_r2[i]+py_i2[i]*py_i2[i], py_r1[i],py_i1[i],py_r2[i],py_i2[i]);

			   if(write_python_files){
				  fprintf(pyfile,"% 12.8f % 12.4g % 12.4g % 12.4g % 12.4g\n",365.25*crossX[i],crossY_r[i],crossY_i[i],py_r1[i]*py_r1[i]+py_i1[i]*py_i1[i],py_r2[i]*py_r2[i]+py_i2[i]*py_i2[i]);
			   }
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
			if (!(fin = fopen(fname,"r")))
			{
			   printf("Unable to open file >%s<\n",fname);
			   exit(1);
			}

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

			logmsg("pwA=%lg pwB=%lg %lg",pw1,whitePSD[p1],whitePSD[p1]/pw1);
			pw1=whitePSD[p1];
			pw2=whitePSD[p2];

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

		 } else {
			printf("%f < %f",toverlap,maxOverlap);
		 }
	  }
   }
   fclose(fout);

   if(write_python_files){
	  fclose(pyfile);
   }
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
void fitPolyFunc(double x,double *v,int nfit,pulsar *psr,int ival,int ipsr)
{
   int i;
   v[0] = 1e-6; // Fit for mean
   v[1] = x*1e-10;
   v[2] = x*x*1e-20;
}

void fitMeanSineFunc(double x,double *v,int nfit,pulsar *psr,int ival,int ipsr)
{
   int i;
   v[0] = 1; // Fit for mean
   v[1] = cos(1e-5*x);
   v[2] = sin(1e-5*x);

}





double getSpectrum(pulsar *psr,double *px,double *py_r,double *py_i,int *nSpec,double toffset,double startOverlap,double endOverlap,double stepMJD,char *covarFuncFile,double T)
{
   /* We make a fake pulsar to allow us to re-use the 'standard' tempo2
	* spectral analysis code
	*/
   int i,j;
   FILE *f;

   longdouble peopoch=toffset;
   pulsar* fakepsr = (pulsar*)malloc(sizeof(pulsar));
   if (*nSpec < 1){
	  *nSpec = (int)((((endOverlap-startOverlap)/(double)(stepMJD))-1)/2.0);
   }

   const int ninObs = psr->dmoffsCMnum;
   int nfakeObs;
   double dt;
   if (notim==0){
	  nfakeObs = (int)(4*((endOverlap-startOverlap)/(double)(stepMJD)+0.5));
	  dt=(endOverlap-startOverlap)/(double)(nfakeObs-1); // this should be the same as stepMJD/4, but just in case.
   } else {
	  // if I simulated the CM then don't interpolate
	  nfakeObs = ninObs;
	  dt=stepMJD;
   }

   double *ifunc_in = (double*)malloc(sizeof(double)*ninObs);
   double *ifunc_out = (double*)malloc(sizeof(double)*nfakeObs);
   double **M = malloc_blas(nfakeObs,ninObs);
   double **Mt= malloc_blas(ninObs,nfakeObs);
   double **MC= malloc_blas(nfakeObs,ninObs);
   double **ifunc_CVM= malloc_uinv(ninObs);
   double ** CVM=malloc_uinv(nfakeObs);

   logmsg("dt=%lf stepMJD=%lf stepMJD/dt=%lf Nin=%d Nout=%d\n",dt,stepMJD,stepMJD/dt,ninObs,nfakeObs);

   fakepsr->param[param_pepoch].val=&peopoch;
   int k=0;

   for(i=0; i < nfakeObs ;i++){
	  ifunc_out[i]=0;
	  for(k=0;k < ninObs;k++){
		 M[i][k]=0;
		 Mt[k][i]=0;
	  }
   }
   // compute the matrix that translates the IFUNC to the finely sampled region.
   double t=startOverlap; // time
   k=0;
   for(i=0; i < nfakeObs ;i++){
	  while (k < ninObs){
		 if(t >=psr->dmoffsCM_mjd[k] && t < psr->dmoffsCM_mjd[k+1]){
			double v = (t-psr->dmoffsCM_mjd[k]) / (psr->dmoffsCM_mjd[k+1]-psr->dmoffsCM_mjd[k]); // fractional way through
			M[i][k] = 1-v;
			M[i][k+1] = v;
			Mt[k][i] = 1-v;
			Mt[k+1][i] = v;
			break;
		 } else{
			k+=1;
		 }
	  }
	  if (k == psr->dmoffsCMnum){
		 // we are "off the end"
		 M[i][k-1] = 1.0;
		 Mt[k-1][i] = 1.0;
	  }
	  t+=dt;
   }

   for (i=0;i<psr->dmoffsCMnum;i++){
	  ifunc_in[i]=psr->dmoffsCM[i];
   }

   fakepsr->obsn = (observation*) malloc(sizeof(observation)*nfakeObs);
   strcpy(fakepsr->name,psr->name);
   fakepsr->nobs=nfakeObs;
   fakepsr->nconstraints=0;

   // compute the interpolated function
   logmsg("M * ifunc");
   TKmultMatrixVec(M,ifunc_in,nfakeObs,ninObs,ifunc_out);
   logmsg("Get ifCVM");

   // Get the CVM for the input IFUNC.
   double r_chisq = psr->fitChisq / (float)(psr->fitNfree);
   int CMoffset=offsetToCM(psr);
   logdbg("CMoffset=%d\n",CMoffset);
   for (i=0;i<ninObs;i++){
	  for (j=0;j<=i;j++){
		 if(i==j){
			//ifunc_CVM[i][j]=psr->covar[i+CMoffset][j+CMoffset]*r_chisq;
			ifunc_CVM[i][i]=pow(psr->dmoffsCM_error[i],2);
		 } else {
			ifunc_CVM[i][j]=0;
			ifunc_CVM[j][i]=0;
		 }
	  }
   }

   // compute the covarince function of the new function.
   logmsg("MC=M * ifCVM");
   TKmultMatrix(M,ifunc_CVM,nfakeObs,ninObs,ninObs,MC);

   // finally compute the covariance function.
   logmsg("CVM = ifCVM * Mt");
   TKmultMatrix(MC,Mt,nfakeObs,ninObs,nfakeObs,CVM);

   // unfortunately we have to throw away the non-diagonal CVM.

   t=startOverlap;
   for(i=0; i < nfakeObs ;i++){
	  fakepsr->obsn[i].sat=t;
	  fakepsr->obsn[i].bat=t;
	  fakepsr->obsn[i].bbat=t;
	  fakepsr->obsn[i].residual=ifunc_out[i];
	  fakepsr->obsn[i].toaErr=sqrt(CVM[i][i]);
	  t+=dt;
   }

   logmsg("%lf %lf, %llf %llf, %s",startOverlap,endOverlap,fakepsr->obsn[0].sat,fakepsr->obsn[fakepsr->nobs-1].sat,fakepsr->name);
   /* From here on in we copy the cholSpectra plugin
   */
   double **uinv;
   FILE *fin;
   char fname[128];
   int ndays=0;
   double resx[fakepsr->nobs],resy[fakepsr->nobs],rese[fakepsr->nobs],sig[fakepsr->nobs];
   int ip[fakepsr->nobs];
   FILE *fout;

   //  printf("Calculating the spectrum\n");
   uinv = malloc_uinv(fakepsr->nobs);
   double sss=0;
   double eee;

   logmsg("Nspec=%d, nobs=%d",*nSpec,fakepsr->nobs);
   for (i=0;i<fakepsr->nobs;i++)
   {
	  resx[i] = (double)(fakepsr->obsn[i].sat-toffset);
	  resy[i] = (double)(fakepsr->obsn[i].residual);
	  rese[i] = fakepsr->obsn[i].toaErr;
	  sig[i]=1.0;
	  ip[i]=i;
	  eee=fakepsr->obsn[i].toaErr/1e6/86400.0/365.25;
	  sss+=1.0/(eee*eee);
   }
   fakepsr->ToAextraCovar=NULL;
   logmsg("Get Cholesky 'uinv' matrix from '%s'",covarFuncFile);
   getCholeskyMatrix(uinv,covarFuncFile,fakepsr,resx,resy,rese,fakepsr->nobs,0,ip);

   logmsg("Remove polynomial");
   double param[99];
   double error[99];
   double chisq;
   double** cvm = (double **)alloca(sizeof(double *)*99);
   for (i=0;i<99;i++){
	  cvm[i] = (double *)alloca(sizeof(double)*99);
	  param[i]=0;
   }

   TKleastSquares_svd_psr_dcm(resx,resy,sig,fakepsr->nobs,param,error,3,cvm,&chisq,fitPolyFunc,0,fakepsr,1.0e-40,ip,uinv);

   logmsg("Poly: %lg %lg %lg",param[0],param[1],param[2]);

   for (i=0;i<fakepsr->nobs;i++){
	  double x=resx[i];
	  resy[i]-=param[0]*1e-6+x*param[1]*1e-10+x*x*param[2]*1e-20;
   }
   logdbg("Got uinv, now compute spectrum.");

   // Must calculate uinv for the pulsar
   calcSpectra_ri_T(uinv,resx,resy,fakepsr->nobs,px,py_r,py_i,*nSpec,T,'N',fakepsr);

   double tspan = (resx[fakepsr->nobs-1]-resx[0])/365.25;
   // Free uinv
   free_blas(M);
   free_blas(Mt);
   free_blas(MC);
   free_blas(uinv);
   free_blas(ifunc_CVM);
   free_blas(CVM);

   free(fakepsr->obsn);
   free(fakepsr);
   return 2*tspan/sss;
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

int offsetToCM(pulsar* psr){
   int i,j,k;
   int n=0;
   n++;
   for (i=0;i<MAX_PARAMS;i++)
   {
	  for (k=0;k<psr->param[i].aSize;k++)
	  {
		 if (psr->param[i].fitFlag[k]==1) /* If we are fitting for this parameter */
		 {
			if (i!=param_start && i!=param_finish)
			{
			   if (i==param_wave_om)
			   {
				  if (psr->waveScale==2)
				  {
					 //                      for (j=0;j<psr->nWhite*2;j++)                           
					 for (j=0;j<psr->nWhite*4;j++)n++;
				  }
				  else
				  {
					 for (j=0;j<psr->nWhite*2;j++)n++;
				  }
			   }
			   else if (i==param_quad_om)
			   {
				  for (j=0;j<psr->nQuad*4;j++)n++;
			   }
			   else if (i==param_ifunc)
			   {
				  for (j=0;j<psr->ifuncN;j++)n++;
			   }
			   else if (i==param_gwsingle)
			   {
				  n+=4;
			   }
			   else if (i==param_dmmodel)
			   {
				  for (j=0;j<(int)psr->dmoffsDMnum;j++)n++;
				  return n;
			   }
			   else
				  n++;
			}
		 }
	  }
   }
}




