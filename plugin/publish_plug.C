/* Plugin to produce output in publication quality format */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"

void dispParameter(int i,int k,pulsar *psr,FILE *fout,int err,double efac);
int nint_derived(double x);
int rnd8(double rval,double rerr,int ifac,char *cval,int *lv,char *cerr,
         int *le,char *msg);
void parseMinus(char *str);
void parseExp(char *str);

extern "C" int tempoOutput(int argc,char *argv[],pulsar *psr,int npsr) 
{  
  int i,k;
  FILE *fout;
  char name[500];
  double pval;
  double efac=1.0;
  char efacStr[500];
  
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-efac")==0)
        sscanf(argv[i+1],"%lf",&efac);
    }
  char format[20]="";
  textOutput(psr,npsr,0,0,0,0,format);
  
  fout = fopen("table.tex","w");
  fprintf(fout,"\\documentclass{article}\n");
  fprintf(fout,"\\begin{document}\n");

  fprintf(fout,"\\begin{table}\n");
  fprintf(fout,"\\begin{tabular}{ll}\n");
  fprintf(fout,"\\hline\\hline\n");
  fprintf(fout,"\\multicolumn{2}{c}{Fit and data-set} \\\\\n");
  fprintf(fout,"\\hline\n");
  strcpy(name,psr[0].name); parseMinus(name);
  fprintf(fout,"Pulsar name\\dotfill & J%s \\\\ \n",name);
  fprintf(fout,"MJD range\\dotfill & %7.1Lf---%7.1Lf \\\\ \n",
          psr[0].param[param_start].val[0],
	 psr[0].param[param_finish].val[0]);
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
            dispParameter(i,k,psr,fout,1,efac);
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
              && i!=param_ephver)
            {
              if (psr[0].param[i].paramSet[k]==1 && 
                  psr[0].param[i].fitFlag[k]==0)
                dispParameter(i,k,psr,fout,0,efac);
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

  fprintf(fout,"\\hline\n");
  fprintf(fout,"\\multicolumn{2}{c}{Assumptions} \\\\\n");
  fprintf(fout,"\\hline\n"); 
  fprintf(fout,"Clock correction procedure\\dotfill & %s \\\\\n",psr[0].clock);
  fprintf(fout,"Solar system ephemeris model\\dotfill & %s \\\\\n",
          psr[0].ephemeris);
  fprintf(fout,"Binary model\\dotfill & %s \\\\\n",psr[0].binaryModel);
  
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
  if (efac==1)
    strcpy(efacStr,"");
  else if (efac==2)
    strcpy(efacStr,"twice");
  else if (efac==3)
    strcpy(efacStr,"three times");
  else
    sprintf(efacStr,"%f times",efac);

  fprintf(fout,"Note: Figures in parentheses are %s the nominal 1$\\sigma$ \\textsc{tempo2} uncertainties in the least-significant digits quoted.\n",efacStr);
  fprintf(fout,"\\end{table}\n");

  fprintf(fout,"\\end{document}\n");
  fclose(fout);
}


void dispParameter(int i,int k,pulsar *psr,FILE *fout,int err,double efac)
{
  char label[1000];
  char valStr[1000];
  int lv,le;
  char cval[500],cerr[500],msg[500];

  if (i==param_raj)
    {
      strcpy(valStr,psr[0].rajStrPost);
      if (err==1)
        {
          int dp,ierr,sym,k;
          double err;
          char disp[500];
          
          strcpy(disp,psr[0].rajStrPost);
          err = psr[0].param[i].err[0]*(efac);
          dp = nint_derived(log10(err*3600.0/M_PI*12.0));
          if (err*3600.0/M_PI*12.0/pow(10.0,(double)dp)<1.90)
            {
              if (dp < 0) dp --;
              else dp ++;
            }
          ierr = (int)(err*3600.0/M_PI*12.0/pow(10.0,(double)dp)+0.9999);
          for (k=0;k<strlen(disp);k++)
            {
              if (disp[k] == ':') sym++;
              if (disp[k] == '.' && dp < 0) {disp[k-dp+1]='\0'; break;}
            }
          sprintf(valStr,"%s(%d)",disp,ierr);
        }
    }
  else if (i==param_decj)
    {
      strcpy(valStr,psr[0].decjStrPost);
      if (err==1)
        {
          int dp,ierr,sym,k;
          double err;
          char disp[500];
          
          strcpy(disp,psr[0].decjStrPost);
          err = psr[0].param[i].err[0]*(efac);
          dp = nint_derived(log10(err*3600.0/M_PI*180.0));
          if (err*3600.0/M_PI*180.0/pow(10.0,(double)dp)<1.90)
            {
              if (dp < 0) dp --;
              else dp ++;
            }
          ierr = (int)(err*3600.0/M_PI*180.0/pow(10.0,(double)dp)+0.9999);
          for (k=0;k<strlen(disp);k++)
            {
              if (disp[k] == ':') sym++;
              if (disp[k] == '.' && dp < 0) {disp[k-dp+1]='\0'; break;}
            }
          sprintf(valStr,"%s(%d)",disp,ierr);
        }
    }
  else
    {
      if (err==1)
        {
          rnd8((double)psr[0].param[i].val[k],psr[0].param[i].err[k]*efac,1,cval,&lv,cerr,&le,msg);
          sprintf(valStr,"%s(%s)",cval,cerr);
        }
      else
        sprintf(valStr,"%Lg",psr[0].param[i].val[k]);
    }
  parseMinus(valStr);
  parseExp(valStr);
  
  if (i==param_raj)
    strcpy(label,"Right ascension, $\\alpha$");
  else if (i==param_decj)
    strcpy(label,"Declination, $\\delta$");
  else if (i==param_f && k==0)
    strcpy(label,"Pulse frequency, $\\nu$ (s$^{-1}$)");
  else if (i==param_f && k==1)
    strcpy(label,"First derivative of pulse frequency, $\\dot{\\nu}$ (s$^{-2}$)");
  else if (i==param_f && k==2)
    strcpy(label,"Second derivative of pulse frequency, $\\ddot{\\nu}$ (s$^{-3}$)");
  else if (i==param_dm && k==0)
    strcpy(label,"Dispersion measure, $DM$ (cm$^{-3}$pc)");
  else if (i==param_dm && k==1)
    strcpy(label,"First derivative of dispersion measure, $\\dot{DM}$ (cm$^{-3}$pc\\,yr$^{-1}$)");
  else if (i==param_pmra)
    strcpy(label,"Proper motion in right ascension, $\\mu_{\\alpha}$ (mas\\,yr$^{-1}$)");
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

  fprintf(fout,"%s\\dotfill & %s \\\\ \n",label,valStr);
  
}


int rnd8(double rval,double rerr,int ifac,char *cval,int *lv,char *cerr,int *le,char *msg)
{
  double vv, ee, xv, xe;
  int ixv, ixe, iee, j, ivv, ise, irnd, ilim,ret,ret_lv,ret_le;
  char cexp[9], fmt[12];

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
	      if (j==0)strcpy(cexp,cexp+1);
	    }
	  while (j==0);
	  if (vv<0) /* allow space for - sign */
	    sprintf(fmt,"%%%d.%df",ise+4,ise);
	  else
	    sprintf(fmt,"%%%d.%df",ise+3,ise);
	  sprintf(cval,fmt,vv);
	  if (cval[0]==' ') strcpy(cval,cval+1);
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
  if (cval[0]==' ') strcpy(cval,cval+1);  /* For positive numbers */
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
