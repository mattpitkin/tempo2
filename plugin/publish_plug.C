/* Plugin to produce output in publication quality format */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"

void dispParameter(int i,int k,pulsar *psr,FILE *fout,int err,double efac,int useCompare,char *compareFile,int *useDagger);
int nint_derived(double x);
int rnd8(double rval,double rerr,int ifac,char *cval,int *lv,char *cerr,
         int *le,char *msg);
void parseMinus(char *str);
void parseExp(char *str);
double fixRA(char *tstr,double err,char *valStr);
double fixDec(char *tstr,double err,char *valStr);

extern "C" int tempoOutput(int argc,char *argv[],pulsar *psr,int npsr) 
{  
  int i,k;
  FILE *fout;
  char name[500];
  double pval;
  double efac=1.0;
  char efacStr[500];
  int nohead=0;
  char compareFile[1024];
  int useCompare=0;
  int useDagger=0;
  char dist1[1024]="NULL";
  char dist2[1024]="NULL";
  double pdot;

  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-efac")==0)
        sscanf(argv[i+1],"%lf",&efac);
      else if (strcmp(argv[i],"-nohead")==0)
	nohead=1;
      else if (strcmp(argv[i],"-dist_tc")==0)
	strcpy(dist1,argv[++i]);
      else if (strcmp(argv[i],"-dist_cl")==0)
	strcpy(dist2,argv[++i]);
      else if (strcmp(argv[i],"-compare")==0)
	{
	  useCompare=1;
	  strcpy(compareFile,argv[++i]);
	}
    }
  char format[20]="";
  textOutput(psr,npsr,0,0,0,0,format);
  
  fout = fopen("table.tex","w");
  if (nohead==0){
    fprintf(fout,"\\documentclass{article}\n");
    fprintf(fout,"\\begin{document}\n");
  }
  if (useCompare==1)
    fprintf(fout,"\\begin{table*}\n");
  else
    fprintf(fout,"\\begin{table}\n");

  strcpy(name,psr[0].name); parseMinus(name);
  fprintf(fout,"\\caption{Parameters for PSR %s}\n",name);
  if (nohead==1){
    fprintf(fout,"\\resizebox{35pc}{!}{\\begin{minipage}{\\textwidth}\n");
  }
  if (useCompare==1)
    fprintf(fout,"\\begin{tabular}{llll}\n");
  else
    fprintf(fout,"\\begin{tabular}{ll}\n");
  fprintf(fout,"\\hline\\hline\n");
  if (useCompare==1)
    {
      fprintf(fout,"\\multicolumn{2}{c}{Fit and data-set} \\\\\n");
      fprintf(fout,"& This work & Previous & Change \\\\\n");
      fprintf(fout,"&           &          & ($\\sigma$) \\\\\n");
    }
  else
    fprintf(fout,"\\multicolumn{2}{c}{Fit and data-set} \\\\\n");
  
  fprintf(fout,"\\hline\n");
  strcpy(name,psr[0].name); parseMinus(name);
  // Check if a J already exists in the name
  if (psr[0].name[0]=='J')
    fprintf(fout,"Pulsar name\\dotfill & %s \\\\ \n",name);
  else
    fprintf(fout,"Pulsar name\\dotfill & J%s \\\\ \n",name);

  fprintf(fout,"MJD range\\dotfill & %7.1Lf---%7.1Lf \\\\ \n",
          psr[0].param[param_start].val[0],
	 psr[0].param[param_finish].val[0]);
  fprintf(fout,"Data span (yr)\\dotfill & %.2Lf \\\\ \n",
          (psr[0].param[param_finish].val[0]-psr[0].param[param_start].val[0])/365.25);
  fprintf(fout,"Number of TOAs\\dotfill & %d \\\\\n",psr[0].nFit);
  fprintf(fout,"Rms timing residual ($\\mu s$)\\dotfill & %.1Lf \\\\\n",
          psr[0].param[param_tres].val[0]);
  fprintf(fout,"Weighted fit\\dotfill & ");
  if (psr[0].fitMode==1) 
    {
      fprintf(fout," Y \\\\ \n");
      fprintf(fout,"Reduced $\\chi^2$ value \\dotfill & %.1f \\\\\n",
              psr[0].fitChisq/(double)psr[0].fitNfree);
    }
  else fprintf(fout," N \\\\ \n");

  /* *******************
     MEASURED Quantities
     ******************* */
  fprintf(fout,"\\hline\n");
  fprintf(fout,"\\multicolumn{2}{c}{Measured Quantities} \\\\ \n");
  fprintf(fout,"\\hline\n");
  
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[0].param[i].aSize;k++)
        {
          // Check if parameter is used and fitted for.  Also check if
          // it isn't START and FINISH because these cause a buffer
          // overflow as they don't have errors. Nevertheless, they
          // _can_ be fitted (somehow), which means they would pass
          // through this if-statement unless sorted out. The same
          // goes for any epoch, but since they cannot be fitted
          // anyway, this is no issue. (JPWV, 5.8.2010)
          if (psr[0].param[i].paramSet[k]==1 && 
              psr[0].param[i].fitFlag[k]==1 && 
              i!=param_start && i!=param_finish){
            dispParameter(i,k,psr,fout,1,efac,useCompare,compareFile,&useDagger);
          }
        }
    }
  
  /* **************
     SET Quantities
     ************** */
  // All non-fitted parameters (except prewhitening terms)
  fprintf(fout,"\\hline\n");
  fprintf(fout,"\\multicolumn{2}{c}{Set Quantities} \\\\ \n");
  fprintf(fout,"\\hline\n");
  
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[0].param[i].aSize;k++)
        {
          if (i!=param_track && i!=param_tres && 
              i!=param_tzrmjd && i!=param_tzrfrq
              && i!=param_start && i!=param_finish 
              && i!=param_waveepoch && i!=param_wave_om
              && i!=param_ephver && i!=param_dmmodel)
            {
              if (psr[0].param[i].paramSet[k]==1 && 
                  psr[0].param[i].fitFlag[k]==0)
                dispParameter(i,k,psr,fout,0,efac,useCompare,compareFile,&useDagger);
            }
        }
    }
   
  /* ******************
     PREWHITENING Terms 
     ****************** */
 if(psr[0].param[param_waveepoch].paramSet[0]==1){
   fprintf(fout,"\\hline\n");
   fprintf(fout,"\\multicolumn{2}{c}{Pre-whitening terms} \\\\ \n");
   fprintf(fout,"\\hline\n");

   // First waves reference epoch
   fprintf(fout,"Reference epoch for fitwaves\\dotfill & %Lg \\\\ \n",
           psr[0].param[param_waveepoch].val[0]);
   // Then fundamental wave frequency
   fprintf(fout,"Fundamental wave frequency, $\\omega_{\\rm pw}$ (yr$^{-1}$)\\dotfill");
   fprintf(fout," & %Lg \\\\ \n",psr[0].param[param_wave_om].val[0]*365.25);

   // Finally the amplitudes of the waves:
   int lv,le;
   char cval[500],cerr[500],msg[500],cval2[500],cerr2[500];
   for(i=0;i<psr[0].nWhite;i++){
     if(psr[0].wave_cos_err[i]!=0.0){
       rnd8((double)psr[0].wave_cos[i],psr[0].wave_cos_err[i]*efac,
            1,cval,&lv,cerr,&le,msg);
       rnd8((double)psr[0].wave_sine[i],psr[0].wave_sine_err[i]*efac,
            1,cval2,&lv,cerr2,&le,msg);
       fprintf(fout,"Wave %d: $A_{\\rm cos, %d}$; $A_{\\rm sin, %d}$\\dotfill & %s(%s); %s(%s) \\\\ \n",i+1,i+1,i+1,cval,cerr,cval2,cerr2);
     }else{
       // No errors to be printed out.
       fprintf(fout,"Wave %d: $A_{\\rm cos, %d}$; $A_{\\rm sin, %d}$\\dotfill & %lg; %lg \\\\ \n",i+1,i+1,i+1,psr[0].wave_cos[i],psr[0].wave_sine[i]);
     }
   }
 }

  /* ******************
     DERIVED Quantities
     ****************** */
  fprintf(fout,"\\hline\n");
  fprintf(fout,"\\multicolumn{2}{c}{Derived Quantities} \\\\\n");
  fprintf(fout,"\\hline\n"); 
  /* Characteristic age */
  pval = (double)(-psr[0].param[param_f].val[0]/2.0
                  /psr[0].param[param_f].val[1]/86400.0/365.25);
  fprintf(fout,"$\\log_{10}$(Characteristic age, yr) \\dotfill & %.2f \\\\\n",
          log10(pval));

  /* Surface magnetic field */
  pval = (double)(sqrt(-psr[0].param[param_f].val[1]/
                       pow(psr[0].param[param_f].val[0],3))*3.2e19);
  fprintf(fout,"$\\log_{10}$(Surface magnetic field strength, G) \\dotfill & %.2f \\\\\n",log10(pval));

  /* Edot */
  pdot = psr[0].param[param_f].val[1]/pow(psr[0].param[param_f].val[0],2);
  pval = (double)4.0/1.0e-15*M_PI*M_PI/pow(1.0/psr[0].param[param_f].val[0],3)*1e30*pdot;
  printf("Have %g %g\n",pdot,pval);
    //(sqrt(-psr[0].param[param_f].val[1]/
    //                       pow(psr[0].param[param_f].val[0],3))*3.2e19);
  fprintf(fout,"$\\log_{10}$(Edot, ergs/s) \\dotfill & %.2f \\\\\n",log10(fabs(pval)));


  /* Distance from parallax */
  if (psr[0].param[param_px].paramSet[0]==1){
    double pxdist,pxdistErr;
    int lv,le;
    char cval[500],cerr[500],msg[500];

    pxdist = (1.0/psr[0].param[param_px].val[0]*1000.0);
    pxdistErr = 1.0/powl(psr[0].param[param_px].val[0],2.0)*psr[0].param[param_px].err[0]*1000.0;
    rnd8((double)pxdist/1000.0,pxdistErr*efac/1000.0,1,cval,&lv,cerr,&le,msg);
    fprintf(fout,"Distance from parallax (kpc) \\dotfill & %s(%s) \\\\\n",cval,cerr);
  }

  /* Distances from command line */
  if (strcmp(dist1,"NULL")!=0)
    fprintf(fout,"TC93 distance (kpc) \\dotfill & %s\\\\\n",dist1);
  if (strcmp(dist2,"NULL")!=0)
    fprintf(fout,"CL02 distance (kpc) \\dotfill & %s\\\\\n",dist2);

  fprintf(fout,"\\hline\n");
  fprintf(fout,"\\multicolumn{2}{c}{Assumptions} \\\\\n");
  fprintf(fout,"\\hline\n"); 
  if (useCompare==0)
    {
      fprintf(fout,"Clock correction procedure\\dotfill & %s \\\\\n",psr[0].clock);
      fprintf(fout,"Solar system ephemeris model\\dotfill & %s \\\\\n",
	      psr[0].ephemeris);
      fprintf(fout,"Binary model\\dotfill & %s \\\\\n",psr[0].binaryModel);
    }
  else
    {
      char str[1024];
      FILE *fin;
      fin = fopen(compareFile,"r");
      while (!feof(fin))
	{
	  if (fscanf(fin,"%s",str)==1)
	    {
	      if (strcmp(str,"EPHEM")==0)
		{
		  fscanf(fin,"%s",str);
		  fprintf(fout,"Solar system ephemeris model\\dotfill & %s & %s \\\\\n",
			  psr[0].ephemeris,str);
		}
	      else if (strcmp(str,"CLK")==0)
		{
		  fscanf(fin,"%s",str);
		  fprintf(fout,"Clock correction procedure\\dotfill & %s & %s \\\\\n",psr[0].clock,str);
		}
	      else if (strcmp(str,"BINARY")==0)
		{
		  fscanf(fin,"%s",str);
		  fprintf(fout,"Binary model\\dotfill & %s & %s\\\\\n",psr[0].binaryModel,str);
		}
	    }
	}
      fclose(fin);
    }
  
  if (psr[0].calcShapiro==-1) 
    fprintf(fout,"Solar system Shapiro delay \\dotfill & N \\\\\n");
  if (psr[0].ipm!=1) 
    fprintf(fout,"Interplanetary medium delay \\dotfill & N \\\\\n");
  if (psr[0].units==TDB_UNITS) 
    fprintf(fout,"TDB units (tempo1 mode)\\dotfill & Y \\\\\n");
  if (psr[0].timeEphemeris==FB90_TIMEEPH) 
    fprintf(fout,"FB90 time ephemeris (tempo1 mode)\\dotfill & Y \\\\\n");
  if (psr[0].t2cMethod == T2C_TEMPO) 
    fprintf(fout,"T2C (tempo1 mode)\\dotfill & Y \\\\\n");
  if (psr[0].planetShapiro==0) 
    fprintf(fout,"Shapiro delay due to planets\\dotfill & N \\\\\n");
  if (psr[0].correctTroposphere==0) 
    fprintf(fout,"Tropospheric delay\\dotfill & N \\\\\n");
  if (psr[0].dilateFreq==0) 
    fprintf(fout,"Dilate frequency\\dotfill & N \\\\\n");
  if (psr[0].ne_sw!=NE_SW_DEFAULT) 
    fprintf(fout,
            "Electron density at 1 AU (cm$^{-3}$)\\dotfill & %.2f \\\\ \n",
            psr[0].ne_sw);
  if (psr[0].units==TDB_UNITS) /* TEMPO1 mode */
    fprintf(fout,"Model version number\\dotfill & %.2f \\\\ \n",2.0);
  else
    fprintf(fout,"Model version number\\dotfill & %.2f \\\\ \n",5.0);

  fprintf(fout,"\\hline\n");
  fprintf(fout,"\\end{tabular}\n");
  if (useDagger==1)
    fprintf(fout,"$\\dagger$ comparison made by updating the new parameters to the epoch of the previous parameters\n");
  if (nohead==1){
    fprintf(fout,"\\end{minipage}}\n");
  }

  if (efac==1)
    strcpy(efacStr,"");
  else if (efac==2)
    strcpy(efacStr,"twice");
  else if (efac==3)
    strcpy(efacStr,"three times");
  else
    sprintf(efacStr,"%f times",efac);

  printf("Suggest adding the following into your paper:\n\n");
  printf("Note: Figures in parentheses are %s the nominal 1$\\sigma$ \\textsc{tempo2} uncertainties in the least-significant digits quoted.\n",efacStr);

  if (useCompare==1)
    fprintf(fout,"\\end{table*}\n");
  else
    fprintf(fout,"\\end{table}\n");
  if (nohead==0)
    {      
      fprintf(fout,"\\end{document}\n");
    }
  fclose(fout);
}

void dispParameter(int i,int k,pulsar *psr,FILE *fout,int err,double efac,int useCompare,char *compareFile,int *useDagger)
{
  char label[1000];
  char valStr[1000];
  int lv,le;
  char cval[500],cerr[500],msg[500];

  if (i==param_raj)
    {
      strcpy(valStr,psr[0].rajStrPost);
      if (err==1)
	fixRA(psr[0].rajStrPost,psr[0].param[i].err[0]*efac,valStr);
    }
  else if (i==param_decj)
    {
      strcpy(valStr,psr[0].decjStrPost);
      if (err==1)
	fixDec(psr[0].decjStrPost,psr[0].param[i].err[0]*efac,valStr);
    }
  else
    {
      if (err==1)
        {
          rnd8((double)psr[0].param[i].val[k],psr[0].param[i].err[k]*efac,1,cval,&lv,cerr,&le,msg);
	  //	  printf("Trying %g %g %g %s %s\n",(double)psr[0].param[i].val[k],(double)psr[0].param[i].err[k],(double)efac,cval,cerr);
          sprintf(valStr,"%s(%s)",cval,cerr);
        }
      else
        sprintf(valStr,"%Lg",psr[0].param[i].val[k]);
   }
  parseMinus(valStr);
  parseExp(valStr);
  
  if (i==param_raj)
    strcpy(label,"Right ascension, $\\alpha$ (hh:mm:ss)");
  else if (i==param_decj)
    strcpy(label,"Declination, $\\delta$ (dd:mm:ss)");
  else if (i==param_f && k==0)
    strcpy(label,"Pulse frequency, $\\nu$ (s$^{-1}$)");
  else if (i==param_f && k==1)
    strcpy(label,"First derivative of pulse frequency, $\\dot{\\nu}$ (s$^{-2}$)");
  else if (i==param_f && k==2)
    strcpy(label,"Second derivative of pulse frequency, $\\ddot{\\nu}$ (s$^{-3}$)");
  else if (i==param_dm && k==0)
    strcpy(label,"Dispersion measure, DM (cm$^{-3}$pc)");
  else if (i==param_dm && k==1)
    strcpy(label,"First derivative of dispersion measure, $\\dot{DM}$ (cm$^{-3}$pc\\,yr$^{-1}$)");
  else if (i==param_pmra)
    strcpy(label,"Proper motion in right ascension, $\\mu_{\\alpha} \\cos \\delta$ (mas\\,yr$^{-1}$)");
  else if (i==param_pmdec)
    strcpy(label,"Proper motion in declination, $\\mu_{\\delta}$ (mas\\,yr$^{-1}$)");
  else if (i==param_px)
    strcpy(label,"Parallax, $\\pi$ (mas)");
  else if (i==param_pb)
    strcpy(label,"Orbital period, $P_b$ (d)");
  else if (i==param_t0)
    strcpy(label,"Epoch of periastron, $T_0$ (MJD)");
  else if (i==param_a1)
    strcpy(label,"Projected semi-major axis of orbit, $x$ (lt-s)");
  else if (i==param_om)
    strcpy(label,"Longitude of periastron, $\\omega_0$ (deg)");
  else if (i==param_omdot)
    strcpy(label,"Periastron advance, $\\dot{\\omega}$ (deg/yr)");
  else if (i==param_sini)
    strcpy(label,"Sine of inclination angle,$\\sin{i}$");
  else if (i==param_m2)
    strcpy(label,"Companion mass, $M_c$ ($M_\\odot$)");
  else if (i==param_ecc)
    strcpy(label,"Orbital eccentricity, $e$");
  else if (i==param_pepoch)
    strcpy(label,"Epoch of frequency determination (MJD)");
  else if (i==param_posepoch)
    strcpy(label,"Epoch of position determination (MJD)");
  else if (i==param_dmepoch)
    strcpy(label,"Epoch of dispersion measure determination (MJD)");
  else if (i==param_pbdot)
    strcpy(label,"First derivative of orbital period, $\\dot{P_b}$");
  else if (i==param_a1dot)
    strcpy(label,"First derivative of $x$, $\\dot{x}$ ($10^{-12}$)");
  else if (i==param_kom)
    strcpy(label,"Longitude of ascending node, $\\Omega$ (degrees)");
  else if (i== param_kin)
    strcpy(label,"Orbital inclination angle, $i$ (degrees)");
  else
    sprintf(label,"%s",psr[0].param[i].label[k]);

  if (useCompare==1)
    {
      FILE *fin;
      char catParam[128];
      char line[128];
      double catVal,catErr;
      char catStr[1024],catValStr[128];
      int nread;
      int match;
      double change;
      double raVal1,decVal1,raVal2,decVal2;
      double thisNewEpochVal;
      double diff,errDiff;
      double epoch1,epoch2,epoch2_p,epoch2_dm,epoch2_pos;
      int changeEpoch;
      char changeStr[128];

      epoch2_p = epoch2_dm = epoch2_pos = -1;
      match=0;
      changeEpoch=0;

      // Find the epochs
      fin = fopen(compareFile,"r");
      while (!feof(fin))
	{
	  fscanf(fin,"%s",catParam);
	  if (strcmp(catParam,"PEPOCH")==0)
	    fscanf(fin,"%lf",&epoch2_p);
	  else if (strcmp(catParam,"DMEPOCH")==0)
	    fscanf(fin,"%lf",&epoch2_dm);
	  else if (strcmp(catParam,"POSEPOCH")==0)
	    fscanf(fin,"%lf",&epoch2_pos);
	}
      fclose(fin);

      fin = fopen(compareFile,"r");
      while (!feof(fin))
	{
	  fgets(line,128,fin);
	  nread = sscanf(line,"%s %s %lf",catParam,catValStr,&catErr);
	  if (nread > 0)
	    {
	      if (strcmp(catParam,psr[0].param[i].shortlabel[k])==0)
		{
		  //		  printf("GOT A MATCH: %s\n",catParam);
		  match=1;
		  
		  if (nread==3)
		    {
		      // Check for RAJ
		      if (i==param_raj)
			{
			  // Error in arcsec
			  catErr *= M_PI/12.0/60.0/60.0;
			  raVal2 = fixRA(catValStr,catErr,catStr);			  
			  
			  //			  printf("raVal = %g %g %g\n",raVal2,catErr,(double)psr[0].param[param_raj].val[0]);
			}
		      else if (i==param_decj)
			{
			  catErr *= M_PI/180.0/60.0/60.0;
			  decVal2 = fixDec(catValStr,catErr,catStr);			  
			}
		      else
			{			  
			  sscanf(catValStr,"%lf",&catVal);
			  rnd8(catVal,catErr,1,cval,&lv,cerr,&le,msg);
			  sprintf(catStr,"%s(%s)",cval,cerr);
			}
		      // Calculate change
		      if (i==param_f && k==0) {
			epoch1 = (double)psr[0].param[param_pepoch].val[0];
			epoch2 = epoch2_p;
			thisNewEpochVal = (double)psr[0].param[i].val[k]+
			  (epoch2-epoch1)*SECDAY*psr[0].param[param_f].val[1];
			changeEpoch=1;
		      }
		      else if (i==param_raj)
			{
			  epoch1 = (double)psr[0].param[param_posepoch].val[0];
			  if (epoch2_pos==-1)
			    epoch2 = epoch2_p;
			  else
			    epoch2 = epoch2_pos;

			  thisNewEpochVal = (double)psr[0].param[param_raj].val[0]+
			    (epoch2-epoch1)*SECDAY*psr[0].param[param_pmra].val[0]/cos(psr[0].param[param_decj].val[0])*M_PI/180.0/1000.0/365.25/SECDAY/60.0/60.0;
			  changeEpoch=1;
			  catVal = raVal2;
			  printf("Testing %f %f %f\n",thisNewEpochVal,catVal,catErr);
			}
		      else if (i==param_decj)
			{
			  epoch1 = (double)psr[0].param[param_posepoch].val[0];
			  if (epoch2_pos==-1)
			    epoch2 = epoch2_p;
			  else
			    epoch2 = epoch2_pos;

			  thisNewEpochVal = (double)psr[0].param[param_decj].val[0]+
			    (epoch2-epoch1)*SECDAY*psr[0].param[param_pmdec].val[0]*M_PI/180.0/1000.0/365.25/SECDAY/60.0/60.0;
			  changeEpoch=1;
			  catVal = decVal2;
			  printf("Testing %.10f %.10f %.10f %f\n",(double)psr[0].param[param_decj].val[0],thisNewEpochVal,catVal,catErr);
			}
		      else
			thisNewEpochVal = (double)psr[0].param[i].val[k];
		    
		      diff = (thisNewEpochVal-catVal);
		      errDiff = sqrt(pow(psr[0].param[i].err[k],2) + pow(catErr,2));
		      change=diff/errDiff;
		      sprintf(changeStr,"%.1f",change);
		      parseMinus(changeStr);
		    }
		  else
		    sprintf(catStr,"%s",catValStr);
	  //	  printf("Trying %g %g %g %s %s\n",(double)psr[0].param[i].val[k],(double)psr[0].param[i].err[k],(double)efac,cval,cerr);


		  break;
		}
	    }
	}
      fclose(fin);
      parseMinus(catStr);
      parseExp(catStr);

      if (match==1)
	{
	  if (nread==3) 	    	   
	    {
	      if (changeEpoch==1)
		{
		  fprintf(fout,"%s\\dotfill & %s & %s & %s$^\\dagger$ \\\\ \n",label,valStr,catStr,changeStr);
		  *useDagger=1;
		}
	      else
		fprintf(fout,"%s\\dotfill & %s & %s & %s \\\\ \n",label,valStr,catStr,changeStr);
	    }
	  else
	    fprintf(fout,"%s\\dotfill & %s & %s & -- \\\\ \n",label,valStr,catStr);
	}
      else
	fprintf(fout,"%s\\dotfill & %s & -- & -- \\\\ \n",label,valStr);
    }
  else
    fprintf(fout,"%s\\dotfill & %s \\\\ \n",label,valStr);
}


double fixRA(char *tstr,double err,char *valStr)
{
  int dp,ierr,sym=0,k;
  double hr,min,sec;
  char disp[500];
  double ra;

  sscanf(tstr,"%lf:%lf:%lf",&hr,&min,&sec);
  ra = (hr+min/60.0+sec/3600.0)/12.0*M_PI;
  
  strcpy(disp,tstr);
  dp = nint_derived(log10(err*3600.0/M_PI*12.0));
  if (err*3600.0/M_PI*12.0/pow(10.0,(double)dp)<1.90)
    {
      if (dp < 0) dp --;
      else dp ++;
    }
  ierr = (int)(err*3600.0/M_PI*12.0/pow(10.0,(double)dp)+0.9999);
  if (dp > 0)
    ierr = (int)(err*3600.0/M_PI*12.0+0.9999);
  for (k=0;k<strlen(disp);k++)
    {
      if (disp[k] == ':') sym++;
      if (disp[k] == '.' && dp < 0) {disp[k-dp+1]='\0'; break;}
      if (disp[k] == '.' && dp >= 0) {disp[k]='\0'; break;}
    }
  sprintf(valStr,"%s(%d)",disp,ierr);

  return ra;
}

double fixDec(char *tstr,double err,char *valStr)
{
  int dp,ierr,sym=0,k;
  char disp[500];
  double deg,min,sec;
  double dec;

  //    strcpy(tstr,"10:11:12.3456789");
  //    err = 32/3600.0*M_PI/180.0;

  //  printf("fixDec = %s %g\n",tstr,err);
  sscanf(tstr,"%lf:%lf:%lf",&deg,&min,&sec);
  dec = (fabs(deg)+min/60.0+sec/3600.0)/180.0*M_PI;
  if (deg < 0) dec=-dec;

  strcpy(disp,tstr);
  dp = nint_derived(log10(err*3600.0/M_PI*180.0)-0.5);
  //  dp = (int)(log10(err*3600.0/M_PI*180.0));
  //  dp=-2;
  //  printf("dp = %d\n",dp);
    if (err*3600.0/M_PI*180.0/pow(10.0,(double)dp)<1.90)
      {
        if (dp < 0) dp --;
        else dp ++;
      }
  //  printf("dp here = %d\n",dp);
  ierr = (int)(err*3600.0/M_PI*180.0/pow(10.0,(double)dp)+0.9999);
  if (dp > 0)
    ierr = (int)(err*3600.0/M_PI*180.0+0.9999);
  //      printf("ierr = %d\n",ierr);
  //  printf("disp = %s %d\n",disp,dp);
  for (k=0;k<strlen(disp);k++)
    {
      if (disp[k] == ':') sym++;
      if (disp[k] == '.' && dp < 0) {disp[k-dp+1]='\0'; break;}
      if (disp[k] == '.' && dp >= 0) {disp[k]='\0'; break;}
    }
  sprintf(valStr,"%s(%d)",disp,ierr);
  //    printf("returning = %s %g\n",valStr,dec);
  return dec;
}

int rnd8(double rval,double rerr,int ifac,char *cval,int *lv,char *cerr,int *le,char *msg)
{
  double vv, ee, xv, xe;
  int ixv, ixe, iee, j, ivv, ise, irnd, ilim,ret,ret_lv,ret_le;
  char cexp[9], fmt[12],temp[20];

  ilim = 7;

  strcpy(cval," ");
  strcpy(cerr," ");
  strcpy(msg," ");
  ret = 0;
  ret_lv = 0;
  ret_le = 0;

  /* Get scale factors */
  vv = fabs(rval);
  ee = rerr*ifac;
  if (vv>0.0)
    xv=log10(vv);
  else
    xv=0.0;

  ixv=(int)fabs(xv);
  if (xv<0.0)ixv=-ixv-1;

  xe = log10(ee/2.0)+1.0e-10;
  ixe = (int)fabs(xe);
  if (xe<0.0)ixe=-ixe-1; 
  ise=ixe-ixv; /* Scale of error wrt value */

  ixv = -ixv; /* Reverse signs of scale factors */
  ixe = -ixe;
  ise = -ise;

  /* Check for encoding as integer */
  if (xv>=log10(2.0) && xv<(ilim+10) && ee>=2.0)
    {
      irnd=(int)xe;
      ivv=nint_derived(vv);
      if (irnd>1)
	{
	  irnd=nint_derived(pow(10.0,(double)irnd));
	  ivv = nint_derived(vv/irnd)*irnd;
	}
      if (ivv!=0)
	ret_lv=(int)(log10((double)ivv+1.0e-10)+2.0);

      if (rval<0.0) ivv=-ivv;
      sprintf(fmt,"%%%dd",ret_lv);
      sprintf(cval,fmt,ivv);
    }
  else /* Encode as real */
    {
      vv=rval*pow(10.0,(double)ixv); /* Scale for mantissa */
      ee=ee*pow(10.0,(double)ixe);   /* scale error */
      if (ixv<-ilim || ixv>3) /* Use exponent notation */
	{
	  if (ise<1) /* If err > 0.2, print val as fn.1 and scale error */
	    {
	      ee=ee*pow(10.0,1-ise);
	      ise=1; 
	    }
	  strcpy(cexp," ");
	  sprintf(cexp,"%d",-ixv); 
	  j=0; /* Strip off leading blanks */
	  do
	    {
	      j=strlen(cexp)-1;
	      if (j==0)
		{
		  strcpy(temp,cexp+1);
		  strcpy(cexp,temp);
		}
	    }
	  while (j==0);
	  if (vv<0) /* allow space for - sign */
	    sprintf(fmt,"%%%d.%df",ise+4,ise);
	  else
	    sprintf(fmt,"%%%d.%df",ise+3,ise);
	  sprintf(cval,fmt,vv);
	  if (cval[0]==' ') 
	    {
	      strcpy(temp,cval+1);
	      strcpy(cval,temp);
	    }
	  ret_lv = strlen(cval)-1;
	  strcat(cval,"E");
	  strcat(cval,cexp);
	}
      else
	{
	  if (ise<1)
	    {
	      if (ixv<1)
		{
		  ixv=1;
		  ise=ise-1;
		}
	      sprintf(fmt,"%%%d.%df",3+ixv,ixv);
	      ee=ee*pow(10.0,-ise);
	    }
	  else
	    {
	      if (ixv<0)ixv=0;
	      if (ixe<1)ixe=1;
	      sprintf(fmt,"%%%d.%df",3+ixv+ise,ixe);
	    }
	  sprintf(cval,fmt,rval);
	}
    }
  if (cval[0]==' ') {
    strcpy(temp,cval+1);
    strcpy(cval,temp);  /* For positive numbers */
  }
  ret_lv=strlen(cval)-1; 

  irnd = (int)log10(ee/2.0);
  if (irnd>1) /* Round error */
    {
      irnd=nint_derived(pow(10.0,irnd));
      iee=(int)(ee/irnd+0.999)*irnd;
    }
  else
    iee = (int)(ee+0.999);  /* Round error up */

  ee=iee;
  ret_le = (int)(log10(ee+0.999)+1.0);
  sprintf(fmt,"%%%dd",ret_le);
  sprintf(cerr,fmt,iee);


  *le = ret_le;
  *lv = ret_lv;


  return 0;
}
int nint_derived(double x){
  int i;
  if(x>0.){
    i=(int)(x+0.5);
  }
  else{
    i=(int)(x-0.5);
  }
  return(i);
}

void parseMinus(char *str)
{
  int i,j=0;
  char str2[500];

  strcpy(str2,str);
  for (i=0;i<strlen(str);i++)
    {
      if (str[i]=='-')
	{
	  str2[j++] = '$';
	  str2[j++] = '-';
	  str2[j++] = '$';
	}
      else if (str[i]=='e')
	str2[j++] = 'E';
      else
	str2[j++] = str[i];
    }
  str2[j] = '\0';
  strcpy(str,str2);
}

void parseExp(char *str)
{
  char str2[500],str3[500],str4[500],str5[500],str6[500];
  int i,j=0;
  char *ptr,*ptr2;

  ptr = strtok(str,"E");
  if (ptr != NULL)
    {
      strcpy(str2,ptr);
      ptr = strtok(NULL,"E");
      if (ptr != NULL)
	{
	  strcpy(str3,str2);
	  strcpy(str5,ptr);
	  for (i=0;i<strlen(str5);i++)
	    {
	      if (str5[i]=='(')
		break;
	    }
	  strcpy(str6,str5+i);
	  strcat(str3,str6);
	  ptr[i] = '\0';
	  strcat(str3,"$\\times 10^{");
	  /* Remove dollar signs from ptr (exponent) */
	  for (i=0;i<strlen(ptr);i++)
	    {
	      if (ptr[i]!='$')
		str4[j++] = ptr[i];
	    }
	  str4[j]='\0';
	  strcat(str3,str4);
	  strcat(str3,"}$");
	  strcpy(str,str3);
	}
    }
  

}
char * plugVersionCheck = TEMPO2_h_VER;
