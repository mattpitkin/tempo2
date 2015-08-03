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
#include <inttypes.h>
#include <math.h>
#include "tempo2.h"
#include <cpgplot.h>
#include "T2toolkit.h"
#include "TKfit.h"

using namespace std;

double iterativeFit(pulsar* psr, double pb, double x,int *npsr);
void _itt(pulsar* psr,int *npsr);
void saveparams(pulsar* from, pulsar* to);

void help() /* Display help */
{
   /* This function should contain usage information about the plugin which should (in general) be accessed */
   /* by the user pressing 'h'                                                                              */
}
#define NIT 4

/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
   char parFile[MAX_PSR][MAX_FILELEN];
   char timFile[MAX_PSR][MAX_FILELEN];
   int i;
   *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

   double pb_start=7000;
   double pb_end=12000;
   double pb_step=1000;
   uint64_t n_m = 5;
   double m[]={2,5,10,15,20};
   double globalParameter;


   bool PLOTIT=false;
   bool OPT=false;

   printf("Graphical Interface: mjk\n");
   printf("Author:              Michael Keith\n");
   printf("Version:             v1.0\n");
   printf(" --- type 'h' for help information\n");

   /* Obtain all parameters from the command line */
   for (i=2;i<argc;i++)
   {
	  if (strcmp(argv[i],"-f")==0)
	  {
		 strcpy(parFile[0],argv[++i]);
		 strcpy(timFile[0],argv[++i]);
	  }

	  if (strcmp(argv[i],"-p")==0)PLOTIT=true;
	  if (strcmp(argv[i],"-o")==0)OPT=true;


	  if (strcmp(argv[i],"-pb")==0){
		 pb_start=atof(argv[++i]);
		 pb_end=atof(argv[++i]);
		 pb_step=atof(argv[++i]);
	  }
	  if (strcmp(argv[i],"-1")==0){n_m=1;
		 m[0]=10;
	  }
   }



   uint64_t n_pb = (pb_step+pb_end-pb_start)/pb_step;
   float xvals[n_pb];

   float map_rms[n_m][n_pb];
   float map_t0[n_m][n_pb];
   float map_ecc[n_m][n_pb];
   float map_f1[n_m][n_pb];
   float map_raj[n_m][n_pb];
   float map_decj[n_m][n_pb];

   double raj_zero = 2*M_PI*(20/24.0+32/60.0/24.0);
   double raj_rad,decj_rad;
   double rad2s = 24.0*3600.0/2.0/M_PI;
   double decj_zero = 2*M_PI*(41/360.0+27/60.0/360.0);
   double dec2as = 360.0*3600.0/2.0/M_PI;


   readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
   readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
   preProcess(psr,*npsr,argc,argv);

   psr->param[param_pb].fitFlag[0] = 0;
   psr->param[param_a1].fitFlag[0] = 0;
   uint64_t ntot=n_m*n_pb;
   double x;
   char newpar[1024];
   double best_rms=1e99;
   int best_ipb=-1;
   int best_im=-1;
   double rmsmin=9999999;
   logmsg("Running %d PB steps",n_pb);
   logmsg("Running %d Mf steps",n_m);
   for (int i_pb=0; i_pb<n_pb; i_pb++){
	  double pb=pb_start+pb_step*i_pb;
	  xvals[i_pb]=pb;
	  double pbHR=pb*24.0;
	  for (int i_m=0; i_m<n_m; i_m++){

		 readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
		 x=pow(m[i_m]*pbHR*pbHR*1.61812,1.0/3.0);
		 logmsg("PB: %lf, Mf: %lf, A1: %lf",pb,m[i_m],x);
		 psr->param[param_t0].fitFlag[0] = 1;
		 psr->param[param_om].fitFlag[0] = 1;
		 rmsmin=iterativeFit(psr,pb,x,npsr);
		 map_rms[i_m][i_pb] = psr->rmsPost*1e-3;
		 map_f1[i_m][i_pb] = psr->param[param_f].val[1]*1e13;
		 map_t0[i_m][i_pb] = psr->param[param_t0].val[0];
		 map_ecc[i_m][i_pb] = psr->param[param_ecc].val[0];
		 raj_rad =  psr->param[param_raj].val[0];
		 raj_rad -= raj_zero;
		 map_raj[i_m][i_pb] = raj_rad*rad2s;
		 decj_rad =  psr->param[param_decj].val[0];
		 decj_rad -= decj_zero;
		 map_decj[i_m][i_pb] = decj_rad*dec2as;

		 if(rmsmin < best_rms){
			best_rms=rmsmin;
			best_ipb=i_pb;
			best_im=i_m;
		 }

		 //sprintf(newpar,"best-%.0f%.0f.par",pb,m[i_m]);
		 //textOutput(psr,*npsr,globalParameter,0,0,1,newpar);/* Display the output */
	  }
   }

   logmsg("Best RMS: %lf",best_rms);
   printf("\n\n");
   printf("PB:\tMfunc\tT0:\tECC:\tF1:\tRAJ:\tDECJ:\n");
   double best_pb=xvals[best_ipb];
   double best_m=m[best_im];
   double best_t0= map_t0[best_im][best_ipb];
   double best_ecc = map_ecc[best_im][best_ipb];
   double best_f1 = map_f1[best_im][best_ipb];
   double best_raj = map_raj[best_im][best_ipb];
   double best_decj = map_decj[best_im][best_ipb];
   double best_f0=0;
   double best_glph=0;
   double best_glf0=0;
   double best_glf1=0;

   double inner_pb=xvals[best_ipb];
   double inner_m=m[best_im];
   double inner_t0= map_t0[best_im][best_ipb];
   double inner_ecc = map_ecc[best_im][best_ipb];
   double inner_f1 = map_f1[best_im][best_ipb];
   double inner_raj = map_raj[best_im][best_ipb];
   double inner_decj = map_decj[best_im][best_ipb];
   double inner_f0=0;
   double inner_glph=0;
   double inner_glf0=0;
   double inner_glf1=0;





   printf("%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf RYYY\n",xvals[best_ipb],m[best_im],map_t0[best_im][best_ipb],map_ecc[best_im][best_ipb],map_f1[best_im][best_ipb],map_raj[best_im][best_ipb],map_decj[best_im][best_ipb]);
   if(OPT){
	  double centre_m = m[best_im];
	  double delta_m=0.5;
	  double centre_pb=xvals[best_ipb];
	  double pb=centre_pb;
	  double delta_pb=pb_step/4.0;
	  double inner_rms=best_rms;
	  double prev_rms=best_rms;
	  double pb_R = centre_pb+pb_step;
	  double pb_L = centre_pb;
	  for(i=0; i<10; i++){
		 pb+=delta_pb;
		 double pbHR=pb*24.0;
		 readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
		 x=pow(centre_m*pbHR*pbHR*1.61812,1.0/3.0);
		 logmsg("PB: %lf, Mf: %lf, A1: %lf",pb,centre_m,x);
		 logmsg("L: %lf R: %lf dPb: %lf",pb_L,pb_R,delta_pb);
		 double rms=iterativeFit(psr,pb,x,npsr);
		 if(rms < inner_rms){
			inner_rms = rms;
			inner_pb=pb;
			inner_m=centre_m;
			inner_f1 = psr->param[param_f].val[1]*1e13;
			inner_t0 = psr->param[param_t0].val[0];
			inner_ecc = psr->param[param_ecc].val[0];
			raj_rad =  psr->param[param_raj].val[0];
			raj_rad -= raj_zero;
			inner_raj = raj_rad*rad2s;
			decj_rad =  psr->param[param_decj].val[0];
			decj_rad -= decj_zero;
			inner_decj = decj_rad*dec2as;
			inner_f0 = psr->param[param_f].val[0];
			inner_glph = psr->param[param_glph].val[0];
			inner_glf0 = psr->param[param_glf0].val[0];
			inner_glf1 = psr->param[param_glf1].val[0];
		 } 
		 if(rms > prev_rms){
			logmsg("No better");
			if(i==0){
			   logmsg("Reverse");
			   delta_pb = - delta_pb;
			   pb+=delta_pb;
			   pb_R=centre_pb;
			   pb_L=centre_pb-pb_step;
			} else {
			   if(delta_pb > 0){
				  pb_R=pb;
				  delta_pb = (pb_L-pb_R)/4.0;
			   } else {
				  pb_L=pb;
				  delta_pb = (pb_R-pb_L)/4.0;
			   }
			}
		 } else {
			logmsg("Better");
			if(delta_pb > 0){
			   // going Right.
			   pb_L = pb-delta_pb;
			} else {
			   // going Left.
			   pb_R = pb-delta_pb;
			}

		 }
		 prev_rms=rms;
		 if(inner_rms < best_rms){
			best_rms=inner_rms;
			best_pb=inner_pb;
			best_m=inner_m;
			best_f1=inner_f1;
			best_t0=inner_t0;
			best_ecc = inner_ecc;
			best_raj=inner_raj;
			best_decj=inner_decj;
			best_f0=inner_f0;
			best_glph=inner_glph;
			best_glf0=inner_glf0;
			best_glf1=inner_glf1;
		 }
	  }
   }

   logmsg("Best RMS: %lf",best_rms);
   printf("\n\n");
   printf("PB:\tMfunc\tT0:\tECC:\tF1:\tRAJ:\tDECJ:\n");
   printf("%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lg \t%lg \t%lg \t%lg RXXX\n",best_pb,best_m,best_t0,best_ecc,best_f1,best_raj,best_decj,best_f0,best_glph,best_glf0,best_glf1);

   FILE *fo;
   int i_pb,i_m;
   char fx[1024];
   char fn[1024];
   fo=fopen("out.asc","w");
   for(i_m=0; i_m<n_m; i_m++){
	  for (i_pb=0; i_pb < n_pb ; i_pb++){
		 fprintf(fo,"%g %g %g %g %g %g %g %g\n",xvals[i_pb],m[i_m],
			   map_rms[i_m][i_pb],map_t0[i_m][i_pb],
			map_ecc[i_m][i_pb], map_f1[i_m][i_pb], map_raj[i_m][i_pb], map_decj[i_m][i_pb]);
	  }
   }
   fclose(fo);




   if(PLOTIT){
	  textOutput(psr,1,0,0,0,1,"best.par");/* Display the output */
	  for (i=0;i<2;i++){
		 if(i==0)cpgopen("1/xs");
		 else cpgopen("out.ps/vcps");
		 cpgsch(0.5);
		 cpgsci(1);
		 float yy=0.99;
		 float yoo=0.15;
		 cpgsvp(0.1,0.95,yy-yoo,yy);
		 cpgswin(pb_start-pb_step*0.6,pb_end+pb_step*0.6,0.0,2.5);
		 cpgbox("ABCTS",1000,2,"ABCTSN",0.5,5);
		 cpglab("","RMS (s)","");
		 for(i_m=0; i_m<n_m; i_m++){
			cpgsci(i_m+1);
			cpgline(n_pb,xvals,map_rms[i_m]);
			cpgpt(n_pb,xvals,map_rms[i_m],0);
		 }

		 cpgsci(1);
		 yy-=yoo;
		 cpgsvp(0.1,0.95,yy-yoo,yy);
		 cpgswin(pb_start-pb_step*0.6,pb_end+pb_step*0.6,57950,58350);
		 cpgbox("ABCTS",1000,2,"ABCTS1N",100,5);
		 cpglab("","T0 (MJD)","");
		 for(i_m=0; i_m<n_m; i_m++){
			cpgsci(i_m+1);
			cpgline(n_pb,xvals,map_t0[i_m]);
			cpgpt(n_pb,xvals,map_t0[i_m],0);
		 }

		 cpgsci(1);
		 yy-=yoo;
		 cpgsvp(0.1,0.95,yy-yoo,yy);
		 cpgswin(pb_start-pb_step*0.6,pb_end+pb_step*0.6,0.71,0.99);
		 cpgbox("ABCTS",1000,2,"ABCTS1N",0.1,5);
		 cpglab("","Ecc","");
		 for(i_m=0; i_m<n_m; i_m++){
			cpgsci(i_m+1);
			cpgline(n_pb,xvals,map_ecc[i_m]);
			cpgpt(n_pb,xvals,map_ecc[i_m],0);
		 }

		 cpgsci(1);
		 yy-=yoo;
		 cpgsvp(0.1,0.95,yy-yoo,yy);
		 cpgswin(pb_start-pb_step*0.6,pb_end+pb_step*0.6,-9,-3.5);
		 cpgbox("ABCTS",1000,2,"ABCTS1N",2,2);
		 cpglab("","F1 (10^-13 /s/s)","");
		 for(i_m=0; i_m<n_m; i_m++){
			cpgsci(i_m+1);
			cpgline(n_pb,xvals,map_f1[i_m]);
			cpgpt(n_pb,xvals,map_f1[i_m],0);
		 }

		 cpgsci(1);
		 yy-=yoo;
		 cpgsvp(0.1,0.95,yy-yoo,yy);
		 cpgswin(pb_start-pb_step*0.6,pb_end+pb_step*0.6,13.05,13.3);
		 cpgbox("ABCTS",1000,2,"ABCTS1N",0.1,10);
		 cpglab("","RAJ (s)","");
		 for(i_m=0; i_m<n_m; i_m++){
			cpgsci(i_m+1);
			cpgline(n_pb,xvals,map_raj[i_m]);
			cpgpt(n_pb,xvals,map_raj[i_m],0);
		 }

		 cpgsci(1);
		 yy-=yoo;
		 cpgsvp(0.1,0.95,yy-yoo,yy);
		 cpgswin(pb_start-pb_step*0.6,pb_end+pb_step*0.6,24.35,24.65);
		 cpgbox("ABCTS1N",1000,2,"ABCTS1N",0.1,10);
		 cpglab("PB","DECJ (arcsec)","");
		 for(i_m=0; i_m<n_m; i_m++){
			cpgsci(i_m+1);
			cpgline(n_pb,xvals,map_decj[i_m]);
			cpgpt(n_pb,xvals,map_decj[i_m],0);
		 }
		 cpgclos();
	  }
   }

   return 0;
}


double iterativeFit(pulsar* psr, double pb, double x,int *npsr){
   double prev_rms=99999;
   double rmsmin=99999;
   int i;
   psr->param[param_pb].val[0]=pb;
   psr->param[param_a1].val[0]=x;
   psr->rmsPost=rmsmin;
   psr->rmsPre=rmsmin;
   saveparams(psr,psr+1);
   saveparams(psr,psr+4);
   psr->param[param_t0].fitFlag[0]=1;
   psr->param[param_om].fitFlag[0]=1;
   _itt(psr,npsr);

   saveparams(psr+4,psr);
   psr->param[param_t0].fitFlag[0]=0;
   psr->param[param_om].fitFlag[0]=1;
   _itt(psr,npsr);

   saveparams(psr,psr+4);
   psr->param[param_t0].fitFlag[0]=1;
   psr->param[param_om].fitFlag[0]=1;
   _itt(psr,npsr);

   saveparams(psr+4,psr);
   psr->param[param_t0].fitFlag[0]=1;
   psr->param[param_om].fitFlag[0]=0;
   _itt(psr,npsr);

   saveparams(psr,psr+4);
   psr->param[param_t0].fitFlag[0]=1;
   psr->param[param_om].fitFlag[0]=1;
   _itt(psr,npsr);

   saveparams(psr+4,psr);

   logmsg("PostRMS: %.2lf us",psr->rmsPost);

   //	  exit(1);
   return psr->rmsPost;
}

void _itt(pulsar* psr,int *npsr){
   int i;
   double prev_rms=999999;
   for (i=0;i<=NIT;i++)
   {
	  formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
	  formResiduals(psr,*npsr,1);    /* Form the residuals                 */
	  if (i<NIT) doFitAll(psr,*npsr,0);   /* Do the fitting     */
	  if (fabs(psr->param[param_ecc].val[0] - 0.5)>0.5){
		 // bad eccentricity
		 logmsg("Bad ECC");
		 saveparams(psr+1,psr);
		 break;
	  }
	  calcRMS(psr,0);
	  psr->rmsPre=prev_rms;
	  //logmsg("%d PreRMS: %.2lf us, PostRMS: %.2lf us",i,psr->rmsPre,psr->rmsPost);
	  if(i>0){
		 if(psr->rmsPost < (psr+4)->rmsPost)saveparams(psr,psr+4);
		 if(psr->rmsPost < 5000 && fabs(psr->rmsPost-psr->rmsPre)<0.1)break;
	  }
	  prev_rms=psr->rmsPost;
   }

}

void saveparams(pulsar* from, pulsar* to){
   int i=0;
   int j=0;
   to->rmsPost=from->rmsPost;
   to->rmsPre=from->rmsPre;
   for(i=0; i<MAX_PARAMS;i++){
	  for(j=0; j<from->param[i].aSize;j++){
		 to->param[i].val[j] = from->param[i].val[j];
		 to->param[i].err[j] = from->param[i].err[j];
		 to->param[i].fitFlag[j] = from->param[i].fitFlag[j];
		 to->param[i].paramSet[j] = from->param[i].paramSet[j];
	  }
   }
}
