/* GW sensitivity plugin */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "TKspectrum.h"
#include "T2toolkit.h"
#include "TKfit.h"
#include "GWsim.h"

using namespace std;

void doPlugin(pulsar *psr,int npsr,int doFitV,char parFile[MAX_PSR_VAL][MAX_FILELEN],char timFile[MAX_PSR_VAL][MAX_FILELEN]);
void getSensCurv(pulsar *psr,int npsr,double **resX,
		 double **resY,double **resE,
		 int *nObs,int doFitV);
int detectSource(pulsar *psr,int npsr,double **resX,
		 double **resY,double **resE);

void help() /* Display help */
{
  /* This function should contain usage information about the plugin which should (in general) be accessed */
  /* by the user pressing 'h'                                                                              */
}


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR_VAL][MAX_FILELEN];
  char timFile[MAX_PSR_VAL][MAX_FILELEN];
  int i;
  double globalParameter;
  int doFitV=0;

  *npsr = 0;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: GWsens\n");
  printf("Author:              author\n");
  printf("Version:             version number\n");
  printf(" --- type 'h' for help information\n");


  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[*npsr],argv[i+1]); 
	  strcpy(timFile[*npsr],argv[i+2]);
	  (*npsr)++;
	}
      else if (strcmp(argv[i],"-dofit")==0)
	doFitV=1;
    }

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr);         /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);    /* Form the residuals                 */
      if (i==0) doFit(psr,*npsr,0);   /* Do the fitting     */
      else textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
    }

  doPlugin(psr,*npsr,doFitV,parFile,timFile);

  return 0;
}

void doPlugin(pulsar *psr,int npsr,int doFitV,char parFile[MAX_PSR_VAL][MAX_FILELEN],char timFile[MAX_PSR_VAL][MAX_FILELEN])
{
  int i,j,k,p;

  int sourceDetected;
  int ndetected=0;
  int it,nit;
  gwSrc gw;
  gwSrc *gw_background;
  long seed = TKsetSeed();
  long double kp[MAX_PSR][3],dist[MAX_PSR];
  double gwfreq;
  double eplus_r,eplus_i,ecross_r,ecross_i;
  long double GWamp1 = 1.0e-12,GWamp2=0.0,GWamp;
  double amp;
  double **checkResY;
  double **resY;
  double **resX;
  double **resE;
  long double **sat0;
  double freq0,freq1;
  int addSS=1; // = 1 to add single source, = 0 to not add single source
  int addGWB=1; // = 1 to add GWB, = 0 otherwise
  int nfreq = 20;
  int ifreq;
  int finish;
  int totCount=0;
  long double gwres,toffset;
  char fname[100];
  FILE *f_result;
  FILE *fout;
  int ngw=1000;
  long double gwbA,alpha,flo,fhi;
  double **gwbRes;

  gw_background = (gwSrc *)malloc(sizeof(gwSrc)*ngw);
  flo = 1.0L/(30*365.25*86400.0L);
  fhi = 1.0L/(2.0*86400.0L);
  gwbA = 1.0e-19L;
  alpha = -2.0/3.0L;


  checkResY = (double **)malloc(MAX_PSR*sizeof(double *));
  gwbRes = (double **)malloc(MAX_PSR*sizeof(double *));
  resY = (double **)malloc(MAX_PSR*sizeof(double *));
  resX = (double **)malloc(MAX_PSR*sizeof(double *));
  resE = (double **)malloc(MAX_PSR*sizeof(double *));
  sat0 = (long double **)malloc(MAX_PSR*sizeof(long double *));
  for (i=0;i<MAX_PSR;i++)
    {
      checkResY[i] = (double *)malloc(MAX_OBSN*sizeof(double));
      gwbRes[i] = (double *)malloc(MAX_OBSN*sizeof(double));
      resY[i] = (double *)malloc(MAX_OBSN*sizeof(double));
      resX[i] = (double *)malloc(MAX_OBSN*sizeof(double));
      resE[i] = (double *)malloc(MAX_OBSN*sizeof(double));
      sat0[i] = (long double *)malloc(MAX_OBSN*sizeof(long double));
    }
  toffset = psr[0].param[param_pepoch].val[0];
  // Store residuals
  for (p=0;p<npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)
	{
	  gwbRes[p][i]=0.0;      

	  if (psr[p].obsn[i].deleted!=0)
	    {
	      printf("Must remove deleted points from the .tim file for psr %s\n",psr[p].name);
	      exit(1);
	    }
	  resY[p][i] = (double)psr[p].obsn[i].residual;
	  resX[p][i] = (double)(psr[p].obsn[i].sat - toffset);
	  resE[p][i] = (double)(psr[p].obsn[i].toaErr*1.0e-6);
	}
    }
  // Determine the idealised site arrival times
  for (j=0;j<5;j++)
    {
      for (p=0;p<npsr;p++)
	{
	  psr[p].nJumps = 0;
	  for(i=0;i<MAX_PARAMS;i++){
	    psr[p].param[i].nLinkTo = 0;
	    psr[p].param[i].nLinkFrom = 0;
	  }
	}
      readParfile(psr,parFile,timFile,npsr); /* Load the parameters       */
      formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,npsr,0);    /* Form the residuals                 */
      for (p=0;p<npsr;p++)
	{
	  for (i=0;i<psr[p].nobs;i++)
	    psr[p].obsn[i].sat -= (long double)psr[p].obsn[i].residual/86400.0L;
	}
    }

  for (p=0;p<npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)
	sat0[p][i] = psr[p].obsn[i].sat;
    }

  nit = 1000;

  f_result = fopen("sensitivity.dat","w");


  for (p=0;p<npsr;p++)
    {
      setupPulsar_GWsim(psr[p].param[param_raj].val[0],psr[p].param[param_decj].val[0],kp[p]);
      dist[p] = 0.91*3.08568e19;
    }

  // Do this procedure for all frequencies
  freq0 = 1.0/(30.0*86400.0*365.25); // 30 years
  freq1 = 1.0/(4*7*86400.0); // 4 weeks
  
  //  for (ifreq=0;ifreq<nfreq;ifreq++)
  for (ifreq=10;ifreq<nfreq;ifreq++)
    {
      //      GWamp1 = 1.0e-11;
      gwfreq = pow(10,log10(freq0)+(log10(freq1)-log10(freq0))*ifreq/(double)nfreq);
      GWamp1 = 0.1*gwfreq;
      GWamp2 = 1e-7*gwfreq;      
      totCount=0;

      printf("======================================================\n");
      printf("gwfreq = %g\n",gwfreq);
      // Loop around different amplitudes to find the 95% level
      do {
	finish=0;
	ndetected=0;
	GWamp = sqrt(GWamp1*GWamp2);
	//		GWamp = 7.0e-13;
	printf("GWamp = %Lg\n",GWamp);
	for (it=0;it<nit;it++)
	  {
	    if ((it+1)%100 == 0)
	      {
		printf("it: %d/%d\r",it+1,nit);
		fflush(stdout);
	      }
	    // For testing purposes we can create white data here
	    //      for (i=0;i<nObs[0];i++)
	    //	resY[0][i] = TKgaussDev(&seed); //*1.0e-6;
	    
	    // Step 2.1: Add GWB
	    if (addGWB==1)
	      {
	        double mean;
		GWbackground(gw_background,ngw,&seed,flo,fhi,gwbA,alpha,1);		
		for (p=0;p<npsr;p++)
		  {
		    mean = 0.0;
		    for (i=0;i<psr[p].nobs;i++)
		      {
			gwbRes[p][i] = 0.0;
			for (k=0;k<ngw;k++)
			  gwbRes[p][i]+=calculateResidualGW(kp[p],&gw_background[k],
							 (psr[p].obsn[i].sat-toffset)*86400.0L,dist[p]);
			mean+=gwbRes[p][i];
		      }
		    for (i=0;i<psr[p].nobs;i++)
		      gwbRes[p][i]-=mean/(double)(psr[p].nobs);
		  }
	      }
	    // Step 2.2: Add single source
	    /* Simulate a GW with a random polarisation and random position */
	    if (addSS==1)
	      {
		gw.phi_polar_g = 0.0; //M_PI/4L
		gw.theta_g     = acos((TKranDev(&seed)-0.5)*2); //puts grav wave on sky.
		gw.phi_g       = TKranDev(&seed)*2*M_PI;   
		gw.phase_g     = TKranDev(&seed)*2*M_PI; 
		//!!!!!!temporary marker for finding this location in file		
		gw.omega_g     = 2.0*M_PI*gwfreq; //((gwfreq - bin_width * (TKranDev(&seed) - 0.5))/86400.0); //somewhere in frequency bin
		
		gw.phi_bin = TKranDev(&seed)*2*M_PI;
		gw.theta_bin = gw.phase_g;
		gw.inc_bin = acos((TKranDev(&seed)-0.5)*2);
		
		//Calculate the polarization "unit vector"
		
		eplus_r =(-4*cos(gw.inc_bin)*cos(2*gw.theta_bin)*sin(2*gw.phi_bin) +
			  (3 + cos(2*gw.inc_bin))*cos(2*gw.phi_bin)*sin(2*gw.theta_bin))/2.;
		eplus_i =((3 + cos(2*gw.inc_bin))*cos(2*gw.phi_bin)*cos(2*gw.theta_bin) +
			  4*cos(gw.inc_bin)*sin(2*gw.phi_bin)*sin(2*gw.theta_bin))/2.;
		ecross_r =(4*cos(gw.inc_bin)*cos(2*gw.phi_bin)*cos(2*gw.theta_bin) +
			   (3 + cos(2*gw.inc_bin))*sin(2*gw.phi_bin)*sin(2*gw.theta_bin))/2.;
		ecross_i = ((3 + cos(2*gw.inc_bin))*cos(2*gw.theta_bin)*sin(2*gw.phi_bin) -
			    4*cos(gw.inc_bin)*cos(2*gw.phi_bin)*sin(2*gw.theta_bin))/2.;
		
		//Calculate the amplitude of the soon to be unit vector
		amp   = sqrt(pow(eplus_r,2) + pow(eplus_i,2) + pow(ecross_r,2) + pow(ecross_i,2));
		
		//Normalize
		eplus_r /= amp;
		eplus_i /= amp;
		ecross_r /= amp;
		ecross_i /= amp;
		
		//Set the complex GW amplitudes
		gw.aplus_g         = GWamp*eplus_r;
		gw.aplus_im_g    = GWamp*eplus_i;
		gw.across_g        = GWamp*ecross_r;
		gw.across_im_g  = GWamp*ecross_i;
		
		setupGW(&gw);  //same as GWsim
		for (p=0;p<npsr;p++)
		  {
		    for (i=0;i<psr[p].nobs;i++)
		      {
			gwres = calculateResidualGW(kp[p],&gw,(long double)resX[p][i]*86400.0L,dist[p]);
			//		if (doFitV==0)
			//		  checkResY[0][i]=resY[0][i]+calculateResidualGW(kp,&gw,(long double)resX[0][i]*86400.0L,dist);
			psr[p].obsn[i].sat = sat0[p][i]+((long double)(resY[p][i])+gwres+gwbRes[p][i])/86400.0L;
			//			psr[p].obsn[i].sat = sat0[p][i]+(gwbRes[p][i])/86400.0L; 
		    //		printf("gwres = %g %g\n",(double)resX[0][i],(double)gwres);
		      }
		  }
		for (p=0;p<npsr;p++)
		  {
		    psr[p].nJumps = 0;
		    for(i=0;i<MAX_PARAMS;i++){
		      psr[p].param[i].nLinkTo = 0;
		      psr[p].param[i].nLinkFrom = 0;
		    }
		  }
		readParfile(psr,parFile,timFile,npsr); /* Load the parameters       */
		formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
		formResiduals(psr,npsr,0);    /* Form the residuals                 */
		if (doFitV==1)
		  {
		    doFit(psr,npsr,0);   /* Do the fitting     */
		    formBatsAll(psr,npsr);         /* Form the barycentric arrival times */
		    formResiduals(psr,npsr,0);    /* Form the residuals                 */
		  }
		for (p=0;p<npsr;p++)
		  {
		    for (i=0;i<psr[p].nobs;i++)
		      {
			checkResY[p][i] = (double)psr[p].obsn[i].residual;
		      }
		  }
		
	      }
	    if (it==0)
	      {
		writeTim("sim0.tim",psr,"tempo2");
		
		fout = fopen("residuals.dat","w");
		for (i=0;i<psr[0].nobs;i++)
		  fprintf(fout,"%g %g\n",resX[0][i],checkResY[0][i]);
		
		fclose(fout);
		//	    exit(1);
	      }
	    // Step 2.3: Run 'black box 2' to determine whether we have detected a single source or not
	    sourceDetected = detectSource(psr,npsr,resX,checkResY,resE);
	    ndetected+=sourceDetected;
	    //      printf("Detected = %d/%d = %g percent\n",ndetected,it+1,(double)ndetected/(double)(it+1)*100.0);
	    
	    // Step 2.4: loop
	  }
	printf("Detected = %d/%d = %g percent\n",ndetected,it+1,(double)ndetected/(double)(it+1)*100.0);
	if ((double)ndetected/(double)(it+1)*100.0 > 94.0 && (double)ndetected/(double)(it+1)*100.0 < 96)
	  finish=1;
	else
	  {
	    if ((double)ndetected/(double)(it+1)*100.0 > 95.0)
	      GWamp1 = GWamp;
	    else
	      GWamp2 = GWamp;
	  }
	totCount++;   
      } while (finish==0 && totCount!=10);
      sprintf(fname,"sim_%.3g_%.3g.tim",(double)gwfreq,(double)GWamp);
      writeTim(fname,psr,"tempo2");
      if (finish==0)
	fprintf(f_result,"%g %g %g **\n",(double)gwfreq,(double)GWamp,(double)ndetected/(double)(it+1)*100.0);
      else
	fprintf(f_result,"%g %g %g\n",(double)gwfreq,(double)GWamp,(double)ndetected/(double)(it+1)*100.0);
      fflush(f_result);
    }
  // SHOULD FREE CHECKRESY
  fclose(f_result);
  free(gw_background);
}




int detectSource(pulsar *psr,int npsr,double **resX,
		 double **resY,
		 double **resE)
{
  double specX[MAX_OBSN],specY[MAX_OBSN],specE[MAX_OBSN];
  double splineX[MAX_OBSN],splineY[MAX_OBSN];
  double modelX[MAX_OBSN],modelY[MAX_OBSN]; // Model of the noise in the spectrum
  double threshold[MAX_OBSN];
  double falsep = 0.001;
  int nSpec,i,j,k,nSpline;
  int detection=0;
  static int time=1;
  FILE *fout;

  //  nSpline = 10;
  //  printf("Hardcoding nSpline = 30\n");
  if (time==1) printf("only considering 1 pulsar\n");

  // Get a power spectrum
  TKspectrum(resX[0],resY[0],resE[0],psr[0].nobs,0,0,0,0,0,2,1,1,specX,specY,&nSpec,0,0);
  if (time==1) fout = fopen("spectrum.dat","w");
  for (i=0;i<nSpec;i++)
    {
      specX[i] = log10(specX[i]);
	specY[i] = log10(specY[i]);
      if (time==1) fprintf(fout,"%g %g\n",specX[i],specY[i]);
    }
  if (time==1) fclose(fout);
  // Do spline fit
  /*  {
    double yd[MAX_OBSN][4],h;
    for (i=0;i<nSpline;i++)
      splineX[i] = specX[0]+(specX[nSpec-1]-specX[0])*i/(double)nSpline;
    TKcmonot(nSpec,specX,specY,yd); // Determines 'yd' for the spline fit
    TKspline_interpolate(nSpec,specX,specY,yd,splineX,splineY,nSpline); // Calculates the spline
    for (i=0;i<nSpline;i++)
      printf("spline %g %g\n",splineX[i],splineY[i]);
      } */
  if (time==1) printf("Not using spline: using cubic fit\n");
  {
    int ncoeff=3;
    if (time==1) printf("ncoeff is hardcoded to %d\n",ncoeff);
    double params[ncoeff];
    double val[ncoeff];
    //
    // TKfitPoly simply defines the type of polynomial that we wish to fit.
    // All fitting is done within TKleastSquares_svd_noErr
    //
    TKleastSquares_svd_noErr(specX,specY,nSpec,params,ncoeff,TKfitPoly);
    if (time==1) printf("Have fitted: %g %g %g %g\n",params[0],params[1],params[2],params[3]);
    if (time==1) fout = fopen("model.dat","w");
    for (i=0;i<nSpec;i++)
      {
	modelX[i] = specX[i];
	TKfitPoly(modelX[i],val,ncoeff);
	modelY[i] = 0.0;
	for (j=0;j<ncoeff;j++)
	  modelY[i] += val[j]*params[j];
	if (time==1) fprintf(fout,"%g %g\n",modelX[i],modelY[i]);
      }
    if (time==1) fclose(fout);
  }
  

  // Get the threshold values
  if (time==1) printf("Should be doing a simulation, but currently assuming chisq with 2dof\n");
  if (time==1) fopen("threshold.dat","w");
  for (i=0;i<nSpec;i++)
    {
      //            threshold[i] = log10(5.99*pow(10,modelY[i])*3.03);
      //            threshold[i] = log10(13.8*pow(10,modelY[i])*1.9);
      //      printf("Have %g\n",-2.0*log(falsep/(double)(nSpec+1)));
            threshold[i] = log10(-2.0*log(falsep/(double)(nSpec+1))*pow(10,modelY[i])); //*1.89);
      //    threshold[i] = log10(22.1*pow(10,modelY[i]));
      if (time==1) fprintf(fout,"%g %g\n",modelX[i],threshold[i]);
      // Check whether we have a detection
      if (specY[i] > threshold[i]) detection=1;
    }
  if (time==1) fclose(fout);
  time=2;
  return detection;
}
char * plugVersionCheck = TEMPO2_h_VER;
