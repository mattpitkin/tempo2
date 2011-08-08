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

/* Plugin to simulate an evolving GW source */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "GWsim.h"


using namespace std;

void ThetaEderivs(double t,double *vals,double *derivs);
void setup3C66B(double *e0,double *theta0,double *mjdOmega0,
		double *mjdObs0,double *mjdLast,int *nObs,
		double *psrdist,double *dist, double *omega0,
		double *mu,double *mc,double *phi);
double psrangle(double centre_long,double centre_lat,double psr_long,double psr_lat);
void setupTest(double *e0,double *theta0,double *mjdOmega0,
		double *mjdObs0,double *mjdLast,int *nObs,
		double *psrdist,double *dist, double *omega0,
		double *mu,double *mc,double *phi);
void RungeKuttaStep(double *y, double *dydx, int n, double *x, double htry, double eps,
	  double *yscal, double *hdid, double *hnext,
		    void (*derivs)(double, double *, double *));
void RungeKuttaCashKarp(double *y, double *dydx, int n, double x, double h, double *yout,
			double *yerr, void (*derivs)(double, double *, double *));

double constA0;
double const2;

       
#define SPEED_LIGHT 299792458.0 /* Speed of light (m/s)                       */
#define SOLAR_MASS  1.98892e30  /* Mass of Sun (kg)      */ 
#define BIG_G       6.673e-11   /* Gravitational constant           */
#define PCM         3.08568025e16        /* one parsec in meters    */
#define MAX_VAL     6000

/* Numerical recipes stuff */
#define MAXSTP 10000
#define TINY 1.0e-30
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define NR_END 1
#define FMAX(x,y) ((x<y)?y:x)
#define FREE_ARG char*

void ode(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	    double hmin, int *nok, int *nbad,
	    void (*derivs)(double, double [], double []),
	    void (*RungeKuttaStep)(double [], double [], int, double *, double, 
			 double, double [],
			 double *, double *, void (*)(double, double [], double [])));

int kmax,kount;
double *xp,**yp,dxsav;


void help() /* Display help */
{
  printf("evolveGW\n\n");
  printf("Required inputs:\n\n");
  printf("-f parfile timefile      Parameter and observation files\n");
  printf("-gwra                    Right ascension of GW source\n");
  printf("-gwdec                   Declination of GW source\n");
  printf("-e0                      Initial orbital eccentricity\n");
  printf("-period                  Initial orbital period (yr)\n");
  printf("-theta0                  Initial orbital phase (deg)\n");
  printf("-phi                     Phi (deg)\n");
  printf("-mc                      Chirp mass (M_solar)\n");
  printf("-gwdist                  GW source distance (pc)\n");
  printf("-epoch                   MJD of GW source parameters\n");
  printf("\n");
  printf("-psrdist                 Distance to pulsar (pc)\n");
  printf("\n\n\n");
  printf("Example for 3C66B source:\n");
  printf("tempo2 -gr GWevolve -f mypar.par mytim.tim -gwra 02:23:11.4 -gwdec 42:59:31 -e0 0.0001 -period 1.05 -theta0 60 -phi 45 -mc 1.3e10 -gwdist 80e6 -epoch 51981.0 -psrdist 1.0e3 \n");
}

long double calcAmp(gwSrc *gw);

/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  double globalParameter;
  double theta[MAX_VAL],e[MAX_VAL];
  double e0,theta0;
  double omega[MAX_VAL],omega0;
  double aVal,bVal,alpha;
  char gwRAs[1000],gwDECs[1000];
  double rPlusE[MAX_VAL],rCrossE[MAX_VAL];
  double rPlusP[MAX_VAL],rCrossP[MAX_VAL];
  double gwRA,gwDEC,psrRA,psrDEC;
  double rPlus,rCross;
  double period;
  long double res[MAX_VAL];
  double t[MAX_VAL],t0;
  double yvals[3];
  double h1,eps,hmin;
  double chi0;
  double mc;
  double dist;
  double mu;
  double phi;
  double resE[MAX_VAL],resP[MAX_VAL];
  double tdelay;
  int nok,nbad;
  int i;
  double toffset;
  double psrdist;
  double mjdOmega0;
  double mjdObs0;
  double mjdLast;
  int    nObs;
  FILE *fo1,*fo2,*fo3;


  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: singleSource\n");
  printf("Author:              G. Hobbs\n");
  printf("Version:             1.0\n");
  printf(" --- type 'h' for help information\n");

  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(parFile[0],argv[i+1]); 
	  strcpy(timFile[0],argv[i+2]); 
	}
      else if (strcmp(argv[i],"-gwra")==0)
	sscanf(argv[++i],"%s",gwRAs);
      else if (strcmp(argv[i],"-gwdec")==0)
	sscanf(argv[++i],"%s",gwDECs);
      else if (strcmp(argv[i],"-e0")==0)
	sscanf(argv[++i],"%lf",&e0);
      else if (strcmp(argv[i],"-theta0")==0)
	sscanf(argv[++i],"%lf",&theta0);
      else if (strcmp(argv[i],"-epoch")==0)
	sscanf(argv[++i],"%lf",&mjdOmega0);
      else if (strcmp(argv[i],"-period")==0)
	sscanf(argv[++i],"%lf",&period);
      else if (strcmp(argv[i],"-phi")==0)
	sscanf(argv[++i],"%lf",&phi);
      else if (strcmp(argv[i],"-gwdist")==0)
	sscanf(argv[++i],"%lf",&dist);
      else if (strcmp(argv[i],"-mc")==0)
	sscanf(argv[++i],"%lf",&mc);
      else if (strcmp(argv[i],"-psrdist")==0)
	sscanf(argv[++i],"%lf",&psrdist);
      else if (strcmp(argv[i],"-h")==0)
	{
	  help();
	  exit(1);
	}
    }  

  // Fix up units 
  theta0 = theta0*M_PI/180.0;
  psrdist *= PCM; 
  dist    *= PCM;
  omega0  = 2.0*M_PI/(period*365.25*86400.0);
  mc      *= SOLAR_MASS;
  phi     *= M_PI/180.0;
  {
    int h,m,d;
    double s;
    sscanf(gwRAs,"%d:%d:%lf",&h,&m,&s);
    if (gwRAs[0]=='-')
      gwRA = (h-m/60.0-s/60.0/60.0)*180.0/12.0;
    else
      gwRA = (h+m/60.0+s/60.0/60.0)*180.0/12.0;
    gwRA*=M_PI/180.0;

    sscanf(gwDECs,"%d:%d:%lf",&d,&m,&s);
    if (gwDECs[0]=='-')
      gwDEC = (d-m/60.0-s/60.0/60.0);
    else
      gwDEC = (d+m/60.0+s/60.0/60.0);
    gwDEC*=M_PI/180;
  }


  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr);         /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);
  psrRA = (double)psr[0].param[param_raj].val[0];
  psrDEC = (double)psr[0].param[param_decj].val[0];

  printf("psr position = %g %g\n",psrRA,psrDEC);
  printf("GW source position = %g %g\n",gwRA,gwDEC);
  mu = psrangle(psrRA,psrDEC,gwRA,gwDEC);
  printf("angle between pulsar and source = %g\n",mu);
  mu*=M_PI/180.0;
  //  mu = 81.5*M_PI/180.0;
  //  exit(1);
    

  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
    {
      formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,0);    /* Form the residuals                 */
      if (i==0) doFit(psr,*npsr,0);   /* Do the fitting     */
      else textOutput(psr,*npsr,globalParameter,0,0,0,"");  /* Display the output */
    }

  // ******************************************************
  // Setup source parameters
  //


      setup3C66B(&e0,&theta0,&mjdOmega0,&mjdObs0,&mjdLast,&nObs,&psrdist,&dist,
        	     &omega0,&mu,&mc,&phi);
  //    setupTest(&e0,&theta0,&mjdOmega0,&mjdObs0,&mjdLast,&nObs,&psrdist,&dist,
  //	     &omega0,&mu,&mc,&phi);

    nObs = psr[0].nobs;
    mjdObs0 = (double)psr[0].obsn[0].sat;
    mjdLast = (double)psr[0].obsn[psr[0].nobs-1].sat;
  // ******************************************************
  // Convert to geometrised units (units of time)
  // i.e. D_sec = D_m/c
  //      M_sec = M_kg G/c^3 
  // Note: keeping pulsar distance in SI units
  mc     = mc*BIG_G/pow(SPEED_LIGHT,3);
  dist   = dist/SPEED_LIGHT;

  tdelay = ((1.0-cos(mu))*psrdist/SPEED_LIGHT);
  printf("tdelay = %g\n",tdelay/365.25/86400.0);
  toffset= mjdObs0-mjdOmega0;

  printf("omega0 = %g\n",omega0*SPEED_LIGHT);
  printf("chirp mass = %g\n",mc);

  fo1 = fopen("out1.dat","w");
  fo2 = fopen("out2.dat","w");
  fo3 = fopen("out3.dat","w");


  /* Calculate A0 using initial conditions and equation 14 in Jenet et al. */
  constA0 = omega0*pow(e0,18.0/19.0)*pow(1.0-e0*e0,-3.0/2.0)*pow(1.0+121.0/304.0*e0*e0,1305.0/2299.0);
  printf("constA0 = %g\n",constA0);

  /* Calculate constant in the de/dt equation */
  chi0 = (1-e0*e0)*pow(e0,-12.0/19.0)*pow(1.0+121.0/304.0*e0*e0,-870.0/2299.0);
  const2 = -304.0/15.0*pow(mc,5.0/3.0)*pow(omega0,8.0/3.0)*pow(chi0,-4);
  printf("chi0 = %g\n",chi0);
  printf("const2 = %g\n",const2);

  /* Now solve for the coupled differential equations (11 and 12 in Jenet et al.)
   * to determine theta and e at time t
   */
  printf("Calculating Earth term\n");
  for (i=0;i<nObs;i++)
    {
      t[i] = (psr[0].obsn[i].sat - mjdOmega0)*86400.0;
      //      t[i] = toffset*86400.0+i*86400.0*((mjdLast-mjdObs0)/(double)nObs);
      t0 = 0.0;
      yvals[0] = theta0;
      yvals[1] = e0;
            
      eps=1.0e-12;    /* Set up properly */
      h1 = 1.0e-3;   /* Set up properly */
      hmin = 0.0;    /* Set up properly */
      ode(yvals,2,t0,t[i],eps,h1,hmin,&nok,&nbad,ThetaEderivs,RungeKuttaStep);
      //      printf("nok = %d, nbad = %d\n",nok,nbad);
      theta[i] = yvals[0];
      e[i]     = yvals[1];
      
      /* Now use value of theta and e to calculate omega */
      omega[i] = constA0*pow(e[i],-18.0/19.0)*pow(1.0-e[i]*e[i],3.0/2.0)*
	pow(1.0+121.0/304.0*e[i]*e[i],-1305.0/2299.0);
      
    }
  printf("Calculating residuals due to Earth term\n");
  /* Now calculate the timing residuals due to the Earth term */
  for (i=0;i<nObs;i++)
    {
      aVal = 2*e[i]*sin(theta[i])*cos(pow(theta[i],2))-
	0.5*sin(2.0*theta[i])*(1.0+e[i]*cos(theta[i]))*4.0;
      bVal = 2*cos(2*theta[i])+e[i]*cos(theta[i]);
      alpha = pow(mc,5.0/3.0)/dist/pow(omega[i],1.0/3.0)*
	sqrt(1.0-e[i]*e[i])/(1.0+e[i]*cos(theta[i]));
      rPlusE[i] = alpha*aVal;
      rCrossE[i] = alpha*bVal;
      resE[i] = 0.5*(1+cos(mu))*(rPlusE[i]*cos(2.0*phi)+rCrossE[i]*sin(2.0*phi));
      printf("ResE %g %g %g\n",(t[i]/86400.0/365.25)+16.375,resE[i],bVal);
      fprintf(fo1,"TimeE %g = %g %g %g %g\n",(t[i]/86400.0/365.25),theta[i],e[i],
	     omega[i],(double)resE[i]);
    }
  printf("Calculating pulsar term\n");
  /* Now do the same again to get the pulsar term */
  for (i=0;i<nObs;i++)
    {
      h1 = 1e-1; // Setup properly
      eps = 1e-5;
      t0 = 0.0;
      yvals[0] = theta0;
      yvals[1] = e0;
            
      ode(yvals,2,t0,(t[i]-tdelay),eps,h1,hmin,&nok,&nbad,ThetaEderivs,RungeKuttaStep);
      //      printf("nok = %d, nbad = %d\n",nok,nbad);
      theta[i] = yvals[0];
      e[i]     = yvals[1];
      
      /* Now use value of theta and e to calculate omega */
      omega[i] = constA0*pow(e[i],-18.0/19.0)*pow(1.0-e[i]*e[i],3.0/2.0)*
	pow(1.0+121.0/304.0*e[i]*e[i],-1305.0/2299.0);
      
    }
  printf("Calculating residuals due to pulsar term\n");
  /* Now calculate the timing residuals due to the Earth term */
  for (i=0;i<nObs;i++)
    {
      aVal = 2*e[i]*sin(theta[i])*cos(pow(theta[i],2))-
	0.5*sin(2.0*theta[i])*(1.0+e[i]*cos(theta[i]))*4.0;
      bVal = 2*cos(2*theta[i])+e[i]*cos(theta[i]);
      alpha = pow(mc,5.0/3.0)/dist/pow(omega[i],1.0/3.0)*
	sqrt(1.0-e[i]*e[i])/(1.0+e[i]*cos(theta[i]));
      rPlusP[i] = alpha*aVal;
      rCrossP[i] = alpha*bVal;
      resP[i] = 0.5*(1+cos(mu))*(rPlusP[i]*cos(2.0*phi)+rCrossP[i]*sin(2.0*phi));
      printf("ResP %g %g %g\n",(t[i]/86400.0/365.25)+16.375,resP[i],bVal);
      fprintf(fo2,"TimePSR %g = %g %g %g %g\n",t[i]/86400.0/365.25,theta[i],e[i],
	     omega[i],(double)resP[i]);
    }
  // Now calculate the final residuals 
  for (i=0;i<nObs;i++)
    {
      rPlus = rPlusE[i]-rPlusP[i];
      rCross = rCrossE[i]-rCrossP[i];
      res[i] = psr[0].obsn[i].residual+(long double)(resE[i]-resP[i]); //0.5*(1+cos(mu))*(rPlus*cos(2.0*phi)+rCross*sin(2.0*phi));
      psr[0].obsn[i].sat+=(long double)(resE[i]-resP[i])/86400.0L;
      fprintf(fo3,"residuals %g %g %g\n",(t[i]/86400.0/365.25)-toffset/365.25,(double)res[i],(double)psr[0].obsn[i].toaErr*1e-6);
    }
  fclose(fo1);
  fclose(fo2);
  fclose(fo3);
  writeTim("evolve.tim",psr,"tempo2");

  // Now write a pure GW data-set
  for (i=0;i<nObs;i++)
    {
      psr[0].obsn[i].sat-=(psr[0].obsn[i].residual/86400.0L);
    }  
  writeTim("evolveGW.tim",psr,"tempo2");

  return 0;
}

/* Define the derivatives of theta and e */
void ThetaEderivs(double t,double *vals,double *derivs)
{
  double theta = vals[0];
  double e = vals[1];

  if (e < 0)
    {
      printf("Big problem because eccentricty has gone negative when \n");
      printf("solving differential equations ... stopping ...\n");
      exit(1);
    }

  //  printf("In here with %g %g %g\n",t,vals[1],vals[2]);
  derivs[0] = constA0*pow(e,-18.0/19.0)*pow(1.0-e*e,3.0/2.0)*
    pow(1.0+121.0/304.0*e*e,-1305.0/2299.0)*
    pow(1.0+e*cos(theta),2)/pow(1-e*e,3.0/2.0);
  
  derivs[1] = const2*pow(e,-29.0/19.0)*pow(1.0-e*e,3.0/2.0)/
    pow(1.0+(121.0/304.0)*e*e,1181.0/2299.0);

}

void setup3C66B(double *e0,double *theta0,double *mjdOmega0,
		double *mjdObs0,double *mjdLast,int *nObs,
		double *psrdist,double *dist, double *omega0,
		double *mu,double *mc,double *phi)
{
  *e0 = 0.0001;            // Initial orbital eccentricity
  *theta0 = -80.0*M_PI/180.0; // WHAT DID RICK SET THIS TO?
  *mjdOmega0 = 51981.0;    // MJD of BHB parameters
  *mjdObs0   = 46436.6869; // MJD of first observation
  *mjdLast   = 48973.7153; // MJD of last observation
  *nObs      = 1000;       // Number of observations

  *psrdist = 1e3*PCM; // Pulsar distance
  *dist = 80e6*PCM; // 80 Mpc
  *omega0 = 2.0*M_PI/(1.05*365.25*86400.0);
  *mu = 81.5*M_PI/180.0;
  *mc   = 1.3e10*SOLAR_MASS;
  *phi = M_PI/4.0;
}

void setupTest(double *e0,double *theta0,double *mjdOmega0,
		double *mjdObs0,double *mjdLast,int *nObs,
		double *psrdist,double *dist, double *omega0,
		double *mu,double *mc,double *phi)
{
  *e0 = 0.00000001;            // Initial orbital eccentricity
  *theta0 = 90.0*M_PI/180.0;
  *mjdOmega0 = 46436.0;    // MJD of BHB parameters
  *mjdObs0   = 46436.6869; // MJD of first observation
  *mjdLast   = 53000.7153; // MJD of last observation
  *nObs      = 3000;       // Number of observations

  *psrdist = 1e3*PCM; // Pulsar distance
  *dist = 20e6*PCM; // 15 Mpc = Virgo cluster
  *omega0 = 2.0*M_PI/(10.0*365.25*86400.0);
  *mu = 85*M_PI/180.0;
  *mc   = 1.0e9*SOLAR_MASS;
  *phi = M_PI/4.0;
}
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

/* Based on odeint() routine in numerical recipes */
/* Modified by G. Hobbs to work within the tempo2 framework
 *  - uses double precision 
 *  - uses TEMPO2 TKroutines
 *  - arrays start from zero
 */
void ode(double *ystart, int nvar, double x1, double x2, double eps, double h1,
	 double hmin, int *nok, int *nbad,
	    void (*derivs)(double, double *, double *),
	    void (*RungeKuttaStep)(double *, double *, int, double *, double, 
			 double, double *,
			 double *, double *, void (*)(double, double *, double *)))
{
  int nstp,i;
  double xsav,x,hnext,hdid,h;
  double yscal[nvar],y[nvar],dydx[nvar];
  
  x=x1;
  h=TKsign_d(h1,x2-x1);
  *nok = (*nbad) = kount = 0;
  for (i=0;i<nvar;i++) y[i]=ystart[i];
  if (kmax > 0) xsav=x-dxsav*2.0;
  for (nstp=1;nstp<=10000;nstp++) {
    (*derivs)(x,y,dydx);
    for (i=0;i<nvar;i++)
      yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+1.0e-30;
    if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
      xp[(++kount)-1]=x;
      for (i=0;i<nvar;i++) yp[i][kount-1]=y[i];
      xsav=x;
    }
    if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
    (*RungeKuttaStep)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
    if (hdid == h) ++(*nok); else ++(*nbad);
    if ((x-x2)*(x2-x1) >= 0.0) {
      for (i=0;i<nvar;i++) ystart[i]=y[i];
      if (kmax) {
	xp[(++kount)-1]=x;
	for (i=0;i<nvar;i++) yp[i][kount-1]=y[i];
      }
      return;
    }
    if (fabs(hnext) <= hmin) {
      printf("Error in ode solver: step size too small\n");
      exit(1);
    }
    h=hnext;
  }
  printf("Error in ode solver: too many steps\n");
  exit(1);
}


/* Based on rkqs() routine in numerical recipes */
/* Modified by G. Hobbs to work within the tempo2 framework
 *  - uses double precision 
 *  - uses TEMPO2 TKroutines
 *  - arrays start from zero
 */

void RungeKuttaStep(double *y, double *dydx, int n, double *x, double htry, double eps,
	  double *yscal, double *hdid, double *hnext,
	  void (*derivs)(double, double *, double *))
{
	int i;
	double errmax,h,xnew,yerr[n],ytemp[n];

	h=htry;
	for (;;) {
	  RungeKuttaCashKarp(y,dydx,n,*x,h,ytemp,yerr,derivs);
	  errmax=0.0;
	  for (i=0;i<n;i++) errmax=TKretMax_d(errmax,fabs(yerr[i]/yscal[i]));
	  errmax /= eps;
	  if (errmax > 1.0) {
	    h=0.9*h*pow(errmax,-0.25);
	    if (h < 0.1*h) h *= 0.1;
	    xnew=(*x)+h;
	    if (xnew == *x) {
	      printf("Error in GWevolve_plug.C - routine: RungeKuttaStep\n");
	      exit(1);
	    }
	    continue;
	  } else {
	    if (errmax > 1.89e-4) *hnext=0.9*h*pow(errmax,-0.2);
	    else *hnext=5.0*h;
	    *x += (*hdid=h);
	    for (i=0;i<n;i++) y[i]=ytemp[i];
	    break;
	  }
	}
}

/* Based on rkck() routine in numerical recipes */
/* Modified by G. Hobbs to work within the tempo2 framework
 *  - uses double precision 
 *  - uses TEMPO2 TKroutines
 *  - arrays start from zero
 */

void RungeKuttaCashKarp(double *y, double *dydx, int n, double x, double h, double *yout,
	double *yerr, void (*derivs)(double, double *, double *))
{
	int i;
	static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
		b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
		dc5 = -277.0/14336.0;
	double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0,dc6=c6-0.25;
	double ak2[n],ak3[n],ak4[n],ak5[n],ak6[n],ytemp[n];

	for (i=0;i<n;i++) ytemp[i]=y[i]+b21*h*dydx[i];
	(*derivs)(x+a2*h,ytemp,ak2);
	for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	(*derivs)(x+a3*h,ytemp,ak3);
	for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	(*derivs)(x+a4*h,ytemp,ak4);
	for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	(*derivs)(x+a5*h,ytemp,ak5);
	for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	(*derivs)(x+a6*h,ytemp,ak6);
	for (i=0;i<n;i++) yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=0;i<n;i++) yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
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
char * plugVersionCheck = TEMPO2_h_VER;
