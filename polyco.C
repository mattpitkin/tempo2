#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
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
/* Routines to generate polycos               
 * Obtained using -polyco <mjd> flag in tempo2 
 *
 * polycos are used to predict a pulsar's parameters at a particular 
 * epoch.  A polyco file contains pulsar ephemerides over a short period of
 * time (typically hours) in the form of a simple polynomial expansion.
 * These polyco files are used for on-line folding of pulsar data while
 * observing.  Ephemerides are calculated on a day-by-day basis centred about
 * the transit time of the pulsar at the observatory. 
 *
 * In the original tempo program the polynomial ephemerides are written to
 * a file called `polyco.dat'.  Entries are listed sequentially within the 
 * file the following file format:
 *
 * On first line:
 * --------------
 * Pulsar Name                                   (pname)
 * Date                                          (date)
 * UTC (hhmmss.ss)                               (utprint)
 * TMID (MJD)                                    (nx fx)
 * DM                                            (dmpsr)
 * Doppler shift due to earth motion (10^-4)     (z)
 * Log_10 of fit rms residual in periods         (rms)
 * 
 * On second line:
 * ---------------
 * Reference Phase
 * Reference rotation frequency (F0)
 * Observatory number
 * Data span (minutes)
 * Number of coefficients 
 * Observing frequency (MHz)
 * Binary phase
 * On third line: (subsequent lines have three coefficients each up to NCOEFF
 * --------------
 * Coefficient 1 
 * Coefficient 2
 * Coefficient 3
 *
 * The pulse phase and frequency at time T are then caclulated as
 * DT = (T-TMID)*1440
 * PHASE = RPHASE + DT*60*F0 + COEFF(1) + DT*COEFF(2) + DT^2*COEFF(3) + ...
 * FREQ(Hz) = F0 + (1/60)*(COEFF(2) + 2*DT*COEFF(3) + 3*DT^2*COEFF(4) + ...)
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"
#include <time.h> 

void atimfake(pulsar *psr,int npsr,int tspan,int ncoeff,longdouble maxha,char *sitename,longdouble freq,
	      longdouble afmjd,longdouble *tmin,longdouble *tmidMJD,int *retTspan,int *nsets,longdouble *val);
void tzFit(pulsar *psr,int npsr,longdouble *tmin,double *doppler,double *rms,double *utc,
	   longdouble tmidMJD,int ncoeff, longdouble *coeff,char *binPhase,int nsets,longdouble afmjd,
	   char* sitename,int dspan,double obsFreq,char *date,longdouble *val,int trueDM,char* polyco_file);
void pcshft(longdouble a,longdouble b,longdouble *d,int n);
void chebpc(longdouble *c,longdouble *d, int n);

void polyco(pulsar *psr,int npsr,longdouble polyco_MJD1,longdouble polyco_MJD2,int nspan,int ncoeff,
	    longdouble maxha,char *sitename,longdouble freq,longdouble coeff[MAX_COEFF],int trueDM,char* polyco_file)
{
  longdouble tsid = 1.0L/1.002737909L;
  longdouble tmidMJD;
  double doppler,rms;
  int dspan;
  char binPhase[32];
  int    i;
  longdouble afmjd;
  longdouble val[800];
  char date[12];
  longdouble tmin[800];
  double globalParameter,utc;
  int nsets;
  const char *CVS_verNum = "$Revision: 1.8 $";

  /* buffers for three polyco outputs -- any input from the user is
   * prepended to the default file names */
  char fname1[256];
  char fname2[256];
  char fname3[256];
  strcpy(fname1,polyco_file);
  strcpy(fname2,polyco_file);
  strcpy(fname3,polyco_file);
  strcat(fname1,"polyco_new.dat");
  strcat(fname2,"newpolyco.dat");
  strcat(fname3,"polyco.tim");

  /* erase any existing files */
  fclose(fopen(fname1,"w"));
  fclose(fopen(fname2,"w"));

  /* special flag for barycenter polycos */
  bool bary = strcmp(sitename,"@")==0 || strcasecmp(sitename,"bat")==0;


  if (displayCVSversion == 1) CVSdisplayVersion("polyco.C","polyco()",CVS_verNum);

  binPhase[0]='\0';

  /* Set some defaults */
  psr[0].param[param_track].paramSet[0]=0;
  globalParameter=0; 

  // Zap the output files so we can append days to them later
  fclose(fopen("polyco_new.dat","w"));
  fclose(fopen("newpolyco.dat","w"));

  for (afmjd=polyco_MJD1; afmjd <= polyco_MJD2; afmjd+=tsid)
    {
      /* Create some arrival times spread throughout the day */
      atimfake(psr,npsr,nspan,ncoeff,maxha,sitename,freq,afmjd,tmin,&tmidMJD,&dspan,&nsets,val);
      for (i=0;i<psr[0].nobs;i++)
	{
	  psr->obsn[i].clockCorr=(!bary);
	  psr->obsn[i].delayCorr=(!bary);
	  psr->obsn[i].phaseOffset = 0.0;
	}
    /* Use reference phase TZR information unless at bary.*/
    if (strcmp(psr->tzrsite,"@")!=0) {
      psr->obsn[0].clockCorr=1;
      psr->obsn[0].delayCorr=1;
    }
    else {
      psr->obsn[0].clockCorr=0;
      psr->obsn[0].delayCorr=0;
    }
      writeTim(fname3,psr,"tempo2"); 
      formBatsAll(psr,npsr);
      formResiduals(psr,npsr,0);
      tzFit(psr,npsr,tmin,&doppler,&rms,&utc,tmidMJD,ncoeff,coeff,binPhase,nsets,afmjd,sitename,nspan,
	    freq,date,val,trueDM,polyco_file);
    }      
  if (psr->tempo1==0)
    printf("WARNING: Should probably use -tempo1 option\n");
}

/* Based on atimfake.f from original TEMPO */
void atimfake(pulsar *psr,int npsr,int tspan,int ncoeff,longdouble maxha,char *sitename,longdouble freq,
	      longdouble afmjd,longdouble *tmin,longdouble *tmidMJD,int *retTspan,int *nsets,longdouble *val)
{
  int i,k;
  longdouble rax,rajVal,mjd2;
  longdouble hlst,wait;
  longdouble x[33];
  longdouble solsid = 1.002737909; /* 1 solar day = 1.002737909 sidereal days */
                                   /* Note: slightly different in almanac p50 */
  int    ntoas=31;  
  int    nmjd;
  int    j;
  double alng;
  longdouble mjd1,val0=0.0L,bma,bpa;

  *retTspan = tspan;

  /* Get observatory coordinates - if barycenter do nothing */
  if (strcmp(sitename,"@")==0 || strcasecmp(sitename,"bat")==0) {
    alng = 0;
  }
  else {
    observatory *obs = getObservatory(sitename);
    /* See tzinit.f for if icoord = 0 : NOT IMPLEMENTED */
    alng = atan2(-obs->y, obs->x);
  }

  /* Conversion of UT into local sidereal time (WHERE DOES 0.154374 come from????) */
  /* Note: 0.154374 is the GST at MJD 0 = 3.716667309 h = 0.15486*24               */

  /* Local sidereal time */
  hlst = 24.0*fortran_mod((longdouble)(solsid*afmjd + 0.154374 - alng/(2*M_PI)),1.0);
  rajVal = psr->param[param_raj].val[0]/M_PI*12.0; /* RAJ in hours */

  /* Why +2 in the next expression ? */
  /* Convert to RAJ at CURRENT epoch */
  rax = (int)rajVal + ((int)((rajVal - ((int)rajVal))*60.0)+2)/60.0;

  /* local hour angle */
  wait = (rax - hlst)/solsid; /* solar to sidereal */

  /* Add a day if not in range */
  if (wait < -maxha) wait+= (24.0/solsid);

  mjd1 = afmjd + (wait-maxha)/24.0+tspan/(2.0*1440.0); /* Compute start time */
  nmjd = (int)mjd1;
  mjd1 = fortran_nint(48.0*(mjd1-(int)mjd1))/48.0; /* Round to nearest 1/2 hour, accurate precision not required here */
  *nsets = (int)((long double)(120.0L*maxha+tspan-1)/(long double)tspan);

  mjd2 = mjd1;
  /* MJD2 is rounded to 1.e-10 -- WHY DO THIS? */  
  /* IN TEMPO NOT BEING ROUNDED CORRECTLY (being cutoff instead) */
  {
    int ntmp1,ntmp2;
    
    ntmp1 = (int)(mjd2*1.0e8L);
    ntmp2 = (int)((mjd2*1.0e8L-ntmp1)*100.0L);
    mjd2 = ntmp2*1e-10L+ntmp1*1.0e-8L;
  }
  mjd1 = mjd2;
  /* This is the first part of the chebft numerical recipes routine 
   * for carrying out a Chebyshev fit 
   */
  {
    longdouble a,b;
    b = tspan/2.0+5.0;
    a = -b;
    bma = 0.5*(b-a);
    bpa = 0.5*(b+a);
  }

  psr->nobs = (ntoas*(*nsets))+1;
  /* First store the reference TOA */
  if (psr->param[param_tzrmjd].paramSet[0]!=1)
    {
      printf("ERROR: must set tzrmjd in the parameter file\n");
      exit(1);
    }
  psr->obsn[0].sat = psr->param[param_tzrmjd].val[0];
  psr->obsn[0].freq = psr->param[param_tzrfrq].val[0]; 
  psr->obsn[0].deleted = 0;
  psr->obsn[0].toaErr = 0;
  psr->obsn[0].efac=1.0;
  strcpy(psr->obsn[0].fname,"POLYCO_REF");
  strcpy(psr->obsn[0].telID,psr->tzrsite); 
  psr->obsn[0].nFlags = 0;

  for (k=1;k<=ntoas;k++)
    x[k-1]=cosl(M_PI*(k-0.5L)/ntoas)*bma+bpa;

  i = -1;
  for (j=0;j<*nsets;j++)
    {
      /* Part of the Chebyshev fit routine */
      mjd2 = mjd1 + (j)*tspan/1440.0L;

      for (k=0;k<ntoas;k++)
	{
	  i++;
	  if (i>800)
	    {
	      printf("ERROR: Tspan too small\n");
	      exit(1);
	    }
	  /* ************************************* */
	  val[i] = mjd2 + x[k]/1440.0L;
	  /* Do some more rounding */
	  {
	    int ntmp1,ntmp2;
	    
	    ntmp1 = (int)(val[i]*1.0e8L);
	    ntmp2 = (int)((val[i]*1.0e8L-ntmp1)*100.0L);
	    val[i] = ntmp2*1e-10L+ntmp1*1.0e-8L;
	  }	  	  
	  if (i==0) val0=val[i];
	  		  
	  psr->obsn[i+1].sat     = val[i]+nmjd; 
	  tmin[i+1]              = 1440.0*(val[i]-val0); 
	  psr->obsn[i+1].freq    = freq; 
	  psr->obsn[i+1].toaErr = 0.0;
	  psr->obsn[i+1].deleted = 0;
	  psr->obsn[i+1].efac=1.0;
	  strcpy(psr->obsn[i+1].fname,"POLYCO");
	  
	  /* HARDCODED TO PARKES ... */
	  strcpy(psr->obsn[i+1].telID,sitename);
	  psr->obsn[i+1].nFlags = 0;
	}


    }
  *tmidMJD = mjd2;
}


void tzFit(pulsar *psr,int npsr,longdouble *tmin,double *doppler,double *rms,double *utc,
	   longdouble tmidMJD,int ncoeff, longdouble *coeff,char *binPhase,int nsets,longdouble afmjd,
	   char* sitename,int tspan,double obsFreq,char *date,longdouble *val,int trueDM,char* polyco_file)
{
  FILE *fout,*fout2;
  int ntoas=32; 
  int nb;
  int j,i,k,iref;
  longdouble fac,sum,rphase,dm;
  longdouble c[31],d[31],ph[31],t[31],phase[800];
  longdouble b,a,rtime,sq,phtst,ct,ct2,phifac,phi0;
  struct tm *timePtr;
  time_t tm;
  char fname1[128];
  char fname2[128];

  strcpy(fname1,polyco_file);
  strcpy(fname2,polyco_file);
  strcat(fname1,"polyco_new.dat");
  strcat(fname2,"newpolyco.dat");
  fout = fopen(fname1,"a");
  fout2 = fopen(fname2,"a");

  /* Note: throughout we are ignoring the first 'TOA' */
  for (i=1;i<psr->nobs;i++)
    phase[i]=psr->obsn[i].phase;

  for (nb = 0;nb < nsets; nb++)
    {
      iref   = (nb+1)*(ntoas-1)-(int)((ntoas-1)/2.0);                            
      rtime  = tmin[iref];
      afmjd = (int)psr->obsn[iref].sat;
      tmidMJD = psr->obsn[iref].sat - (int)psr->obsn[iref].sat;

      tm = (time_t)((double)afmjd-40587.0)*SECDAY; /* Number of seconds since 01-Jan-70 */
      timePtr = gmtime(&tm);
      strftime(date,100,"%d-%b-%y",timePtr);
      if (date[0]=='0') date[0]=' ';
      rphase = psr->obsn[iref].phase;
      //      printf("rphase = %d %d %Lg\n",ntoas,iref,rphase);
      i=iref-(int)((ntoas-1)/2.0)-1; /* -1;*/ /* Maybe another -1 */
      for (j=1;j<ntoas;j++) 
	{      
	  i++;
	  t[j]  = tmin[i]-rtime; 
	  ph[j] = phase[i]-rphase-t[j]*psr->param[param_f].val[0]*60.0L;        
	}
      
      /* This bit of code is based on 'chebft' in numerical recipes for fitting
       * a Chebyshev polynomial
       */
      
      fac = 2.0/(ntoas-1);
      for (j=1;j<ntoas;j++)
	{
	  sum = 0.0;
	  for (k=1;k<ntoas;k++)
	    sum+=ph[k]*cos((M_PI*(j-1))*((k-0.5)/(longdouble)(ntoas-1.0)));
	  c[j-1] = fac*sum;
	}
      
      b = tspan/2.0+5.0;
      a = -b;
      
      chebpc(c,d,ncoeff);
      pcshft(a,b,d,ncoeff);
      for (i=0;i<ncoeff;i++)
	coeff[i] = d[i];
      sq=0.0;
            
      longdouble arg;
      for (j=1;j<ntoas;j++)
	{
	  phtst = d[0];
	  arg = t[j];
	  for (k=1;k<ncoeff;k++)
	    {
	      phtst+=d[k]*arg;
	      arg *= t[j];
	    }
	  sq+=(phtst-ph[j])*(phtst-ph[j]);
	}
      *rms = 1.0e6*sqrt(sq/(ntoas-1))/psr->param[param_f].val[0];
            
      /* Calculate UTC */
      {
	int nutsec,nuthrs,nutmin;
	longdouble uts;
	nutsec = (int)fortran_mod((longdouble)(86400.0*tmidMJD+0.005),86400.0);
	nuthrs = (int)(nutsec/3600.0);
	nutmin = (int)((nutsec-3600*nuthrs)/60.0);
	uts = nutsec-3600.0*nuthrs-60.0*nutmin;
	*utc = 10000.0*nuthrs+100.0*nutmin+uts;
      }
      
      /* Calculate binary phase */
      strcpy(binPhase," ");
      if (psr->param[param_pb].paramSet[0]==1) /* If binary pulsar */
	{
	  ct2 = psr->obsn[1].bat;   
	  ct = ct2 + (psr->obsn[iref].sat-psr->obsn[1].sat);
	  phifac=8.64e4/(psr->param[param_pb].val[0]*86400.0);
	  phi0=fortran_mod((ct-psr->param[param_t0].val[0])*phifac+90000.0,1.0); 
	  sprintf(binPhase,"%7.4Lf%9.4Lf",phi0,phifac);
	}
      *doppler = ((psr->obsn[iref].freqSSB-psr->obsn[iref].freq*1.0e6)/(psr->obsn[iref].freq*1.0e6));
      (*doppler)*=1.0e6/100.0; /* WHY THIS SCALING? */
      *rms     = log10(1.0e-6*(*rms)*psr->param[param_f].val[0]);
	  
      /* Output in original TEMPO format */
      if (psr->name[0]=='J') fprintf(fout,"%-10.10s ",(psr->name)+1);
      else fprintf(fout,"%-10.10s ",psr->name);
      fprintf(fout,"%9.9s",date);
      fprintf(fout,"%11.2f",*utc);
      fprintf(fout,"%20.11f",(double)(afmjd+tmidMJD));
      if (trueDM==0) dm = psr->param[param_dm].val[0];
      else dm = (psr->obsn[iref].tdis1+psr->obsn[iref].tdis2)*DM_CONST*1.0e-12*psr->obsn[iref].freqSSB*psr->obsn[iref].freqSSB;
	
      fprintf(fout,"%21.6f ",(double)dm);
      fprintf(fout,"%6.3f",*doppler);
      fprintf(fout,"%7.3f",*rms);
      fprintf(fout,"\n");
      fprintf(fout,"%20.6Lf",rphase);
      fprintf(fout,"%18.12f",(double)psr->param[param_f].val[0]);
      fprintf(fout,"%5s",sitename);
      fprintf(fout,"%5d",tspan);
      fprintf(fout,"%5d",ncoeff);
      fprintf(fout,"%10.3f",obsFreq);
      fprintf(fout,"%16s",binPhase);
      fprintf(fout,"\n");
      for (i=0;i<ncoeff;i++)
	{
	  fprintf(fout,"%25.17le",(double)coeff[i]);
	  if ((i+1)%3==0) fprintf(fout,"\n");
	}
      if ((i%3)!=0) fprintf(fout,"\n");
      
      fprintf(fout2,"TEMPO2: POLYCO TEMPO1 emulation\n");
      if (psr->name[0]=='J') fprintf(fout2,"%-10.10s\n",(psr->name)+1);
      else fprintf(fout2,"%-10.10s\n",psr->name);
      fprintf(fout2,"%-9.9s\n",date);
      fprintf(fout2,"%-11.2f\n",*utc);
      fprintf(fout2,"%-.20Lf\n",(afmjd+tmidMJD));
      fprintf(fout2,"%-25.10Lf\n",dm);
      fprintf(fout2,"%-10.7f\n",*doppler);
      fprintf(fout2,"%-7.3f\n",*rms);
      fprintf(fout2,"%-.15Lf\n",rphase);
      fprintf(fout2,"%-.20Lf\n",psr->param[param_f].val[0]);
      fprintf(fout2,"%-5s\n",sitename);
      fprintf(fout2,"%-5d\n",tspan);
      fprintf(fout2,"%-5d\n",ncoeff);
      fprintf(fout2,"%-10.3f\n",obsFreq);
      if (strlen(binPhase)<2) fprintf(fout2,"0.000000 0.00000\n");
      else fprintf(fout2,"%-16s\n",binPhase);
      
      for (i=0;i<ncoeff;i++)
	fprintf(fout2,"%-.30Le\n",coeff[i]);
    }
  fclose(fout);
  fclose(fout2);
}


/* Copied from chebpc.f */
/*
 * The following routines, chebpc and pcshft are taken from the numerical recipes
 * book.  They are used to get a polynomial approximation from Chebyshev coefficients
 *
 * It is not completely clear why Chebyshev polynomials are used for this:
 *   - They can avoid the use of extremely high numbers
 *   - Chebyshev polynomials are also useful if the series is calculated for large 'n',
 *     but then truncated to smaller 'n'. However, I don't think that this is done!
 */

void chebpc(longdouble *c,longdouble *d, int n)
{
  int j,k;
  longdouble dd[50],sv; /* 50 from NMAX in chebpc.f */
  for (j=0;j<n;j++)
    {
      d[j] = 0.0;
      dd[j] = 0.0;
    }
  d[0] = c[n-1];

  for (j=n-2;j>=1;j--)
    {
      for (k=n-1-j;k>=1;k--)
	{
	  sv = d[k];
	  d[k] = 2.0*d[k-1]-dd[k];
	  dd[k] = sv;
	}

      sv = d[0];
      d[0] = -dd[0]+c[j];
      dd[0] = sv;
    }
  for (j=n-1;j>=1;j--)
    d[j]=d[j-1]-dd[j];
  d[0] = -dd[0]+0.5*c[0];
} 


void pcshft(longdouble a,longdouble b,longdouble *d,int n)
{
  longdouble constant,fac;
  int j,k;
  constant = 2.0/(b-a);
  fac = constant;
  for (j=2;j<=n;j++)
    {
      d[j-1] = d[j-1]*fac;
      fac=fac*constant;
    }
  constant = 0.5*(a+b);
  for (j=1;j<n;j++)
    {
      for (k=n-1;k>=j;k--)
	d[k-1]=d[k-1]-constant*d[k];
    }
}
