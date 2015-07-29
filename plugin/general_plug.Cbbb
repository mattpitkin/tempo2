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

/* This plugin allows the user to define the output format 
 * The output format can be defined on the command line or in a file 
 *
 * On the command line: use -s "....." option
 * In a file: use -file filename
 *
 * To print the value of a particular parameter use e.g. {F0_v}, {DM_v}, {RAJ_v} etc.
 * To print the error of a parameter use e.g. {F0_e}, {DM_e}
 * To print the label for a parameter use {F0_l}, {DM_l}
 * To display all parameters that are set use {ALL_l} {ALL_v} {ALL_e}
 * To display a parameter in "publication" mode us {F0_p} or {ALL_p}
 *
 * Other commands
 *
 * {ERRMULT x} multiple all errors by x
 * {RAD}       display RAJ and DECJ in radians
 * {NORAD}     display RAJ and DECJ in hms/dms respectively
 * {FORMAT c}  sets the format for displaying values and errors (printf-type format)
 * {TAB x}     move the cursor to position x
 * {NULL c}    set the null string to display if a parameter is not set
 * {UNITS}     choose default units for output
 *
 * All other characters are printed 
 * \n newline
 * \{ print a { character
 * \} print a } character
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"

int nint_derived(double x);
int rnd8(double rval,double rerr,int ifac,char *cval,int *lv,char *cerr,int *le,char *msg);
void parseLine(pulsar *psr,char *line,double *errMult,char *null,char *format,int *rad);

extern "C" int tempoOutput(int argc,char *argv[],pulsar *psr,int npsr) 
{  
  FILE *fin;
  char filename[100];
  char line[500];
  char var[500],type[100],cline[1000];
  char null[500]="*";
  char format[500]="%.14Lg"; 
  char cval[500],cerr[500],msg[500];
  int  i,j,varN,file=-1,rad=0;
  double errMult=1;

  for (i=0;i<argc;i++)
    {
      if (strcasecmp(argv[i],"-file")==0)
	{
	  file=1;
	  strcpy(filename,argv[i+1]);
	}
      else if (strcasecmp(argv[i],"-s")==0)
	{
	  file=0;
	  strcpy(cline,argv[i+1]);
	}
    }
  if (file==-1)
    {
      printf("must use -s or -file option\n");
      exit(1);
    }
  if (file==1) /* Read from file */
    {
      if (!(fin = fopen(filename,"r")))
	{
	  printf("Unable to open file %s\n",filename);
	  exit(1);
	}
    
      while (!feof(fin))
	{
	  if (fgets(line,500,fin)!=NULL)
	    {
	      /* Now parse line */
	      parseLine(psr,line,&errMult,null,format,&rad);
	    }
	}
    }
  else if (file==0) /* Read from string */
    parseLine(psr,cline,&errMult,null,format,&rad);
}

void parseLine(pulsar *psr,char *line,double *errMult,char *null,char *format,int *rad)
{
  int i,j,k,l;
  char var[500],type[100];
  int varN,varA,lv,le,parameter=0,parameterA=0,end=1;
  char cval[500],cerr[500],msg[500];
  char disp[1000];
  int pos=0;
  int units=0;
  double unitVal=1;

  do {
    pos=0;
    for (i=0;i<strlen(line);i++)
      {	  
	if (line[i]=='\\' && (line[i+1]=='{' || line[i+1]=='}'))
	  { /* Do nothing */ }
	else if (line[i]=='\\' && line[i+1]=='n')
	  {
	    i++;
	    printf("\n");
	  }
	else if (line[i]=='{' && line[i-1]!='\\') /* Have command */
	  {
	    for (j=i+1;j<strlen(line);j++)
	      {
		if (line[j]=='}')
		  break;
	      }
	    strcpy(var,line+i+1);
	    var[j-i-1]='\0';
	    if (strstr(var,"FORMAT")!=NULL || strstr(var,"format")!=NULL)
	      {
		strcpy(type,var+6);
		strcpy(format,type);
	      }
	    else if (strstr(var,"TAB")!=NULL)
	      {
		int reqSpace,k;
		sscanf(var+4,"%d",&reqSpace);
		if (reqSpace>pos)
		  {
		    for (k=pos;k<reqSpace;k++)
		      printf(" ");
		    pos=reqSpace;
		  }
	      }
	    else if (strstr(var,"UNITS")!=NULL)
	      units=1;
	    else if (strstr(var,"NORAD")!=NULL)
	      *rad=0;
	    else if (strstr(var,"RAD")!=NULL)
	      *rad=1;
	    else if (strstr(var,"ERRMULT")!=NULL || strstr(var,"errmult")!=NULL)
	      sscanf(var+7,"%lf",errMult);
	    else if (strstr(var,"NULL")!=NULL || strstr(var,"null")!=NULL)
	      sscanf(var+5,"%s",null);
	    else if (strcasecmp(var,"name")==0)
	      {
		sprintf(disp,"%s",psr[0].name); 
		printf("%s",disp);
		pos+=strlen(disp);		
	      }
	    else if (strcasecmp(var,"nobs")==0)
	      {
		sprintf(disp,"%d",psr[0].nobs); 
		printf("%s",disp);
		pos+=strlen(disp);		
	      }
	    else if (strcasecmp(var,"prerms")==0)
	      {
		sprintf(disp,"%.3g",psr[0].rmsPre); 
		printf("%s",disp);
		pos+=strlen(disp);		
	      }
	    else if (strcasecmp(var,"postrms")==0)
	      {
		sprintf(disp,"%.3g",psr[0].rmsPost); 
		printf("%s",disp);
		pos+=strlen(disp);		
	      }
	    else if (strcasecmp(var,"tspan")==0)
	      {
		double start=55000,end=-1;
		for (l=0;l<psr[0].nobs;l++)
		  {
		    if (start > psr[0].obsn[l].sat)
		      start = (double)psr[0].obsn[l].sat;
		    if (end < psr[0].obsn[l].sat)
		      end = (double)psr[0].obsn[l].sat;
		  }
		sprintf(disp,"%.3f",end-start); 
		printf("%s",disp);
		pos+=strlen(disp);				
	      }

	    var[j-i-3]='\0';
	    varN=-1;
	    
	    strcpy(type,line+i+1+j-i-2);
	    i=j;
	    type[1]='\0';
	    for (j=0;j<MAX_PARAMS;j++)
	      {
		for (k=0;k<psr[0].param[j].aSize;k++)
		  {
		    if (strcasecmp(var,psr[0].param[j].shortlabel[k])==0)
		      {
			varN=j;
			varA=k;
			j=MAX_PARAMS;
			break;
		      }
		  }
	      }
	    if (strcasecmp(var,"ALL")==0)
	      {
		if (psr[0].param[parameter].paramSet[parameterA]==1)
		  {
		    varN=parameter;
		    varA=parameterA;
		  }
		end=0;
	      }
	    if (varN>=0)
	      {
		unitVal=1;
		if (units==1 && varN==param_f && varA==1)      unitVal=1e15;
		else if (units==1 && varN==param_f && varA==2) unitVal=1e24;

		if (psr[0].param[varN].paramSet[varA]==0)
		  {
		    printf("%s",null);
		    pos+=strlen(null);
		  }
		else
		{
		  if (strcasecmp(type,"v")==0) /* Display value */
		    {
		      if (varN==param_raj && *rad==0)
			{printf("%s",psr[0].rajStrPost); pos+=strlen(psr[0].rajStrPost);}
		      else if (varN==param_decj && *rad==0)
			{printf("%s",psr[0].decjStrPost); pos+=strlen(psr[0].decjStrPost);}
		      else
			{
			  sprintf(disp,format,psr[0].param[varN].val[varA]*unitVal); 
			  printf("%s",disp);
			  pos+=strlen(disp);
			}
		    }
		    if (strcasecmp(type,"l")==0) /* Display label */
		      {
			sprintf(disp,"%-10.10s",psr[0].param[varN].shortlabel[varA]); 
			printf("%s",disp);
			pos+=strlen(disp);
		      }
		    if (strcasecmp(type,"f")==0) /* Display full label */
		      {
			sprintf(disp,"%s",psr[0].param[varN].label[varA]); 
			printf("%s",disp);
			pos+=strlen(disp);
		      }
		  else if (strcasecmp(type,"e")==0) /* Display error */
		    {
		      if (varN==param_raj && *rad==0)
			{
			  sprintf(disp,format,psr[0].param[varN].err[0]*(*errMult)*12.0*60.0*60.0/M_PI);
			  printf("%s",disp);
			  pos+=strlen(disp);
			}
		      else if (varN==param_decj && *rad==0)
			{
			  sprintf(disp,format,psr[0].param[varN].err[0]*(*errMult)*180.0*60.0*60.0/M_PI);
			  printf("%s",disp);
			  pos+=strlen(disp);
			}
		      else
			{
			  sprintf(disp,format,psr[0].param[varN].err[varA]*(*errMult)*unitVal);
			  printf("%s",disp);
			  pos+=strlen(disp);
			}
		    }
		  else if (strcasecmp(type,"p")==0) /* Publication format */
		    {
		      if (psr[0].param[varN].err[varA]>0)
			{
			  rnd8((double)psr[0].param[varN].val[varA]*unitVal,psr[0].param[varN].err[varA]*(*errMult)*unitVal,1,cval,&lv,cerr,&le,msg);
			  /* More complicated if raj, decj in string format */
			  if (varN==param_raj && *rad==0)
			    {
			      int dp,ierr,sym,k;
			      double err;
			      strcpy(disp,psr[0].rajStrPost);
			      err = psr[0].param[varN].err[0]*(*errMult);
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
			      sprintf(disp,"%s(%d)",disp,ierr);
			      printf("%s",disp);
			      pos+=strlen(disp);
			    }
			  else if (varN==param_decj && *rad==0)
			    {
			      int dp,ierr,sym,k;
			      double err;
			      strcpy(disp,psr[0].decjStrPost);
			      err = psr[0].param[varN].err[0]*(*errMult);
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
			      sprintf(disp,"%s(%d)",disp,ierr);
			      printf("%s",disp);
			      pos+=strlen(disp);
			    }
			  else
			    {
			      printf("%s(%s)",cval,cerr);
			      pos+=strlen(cval)+strlen(cerr)+2;
			    }
			}
		      else
			{
			  if (varN==param_raj && *rad==0)
			    {printf("%s",psr[0].rajStrPost); pos+=strlen(psr[0].rajStrPost);}
			  else if (varN==param_decj && *rad==0)
			    {printf("%s",psr[0].decjStrPost); pos+=strlen(psr[0].decjStrPost);}
			  else
			    {
			      sprintf(disp,format,psr[0].param[varN].val[varA]*unitVal);				
			      printf("%s",disp);
			      pos+=strlen(disp);
			    }
			}
		    }
		}
	      }
	  }
	else
	  {
	    printf("%c",line[i]);	      
	    pos++;
	  }
      }
    if (end==0) 
      {
	do {
	  parameterA++;
	  if (parameterA == psr[0].param[parameter].aSize)
	    {
	      parameterA=0;
	      parameter++;
	    }
	} while (parameter<MAX_PARAMS && psr[0].param[parameter].paramSet[parameterA]!=1);
	if (parameter==MAX_PARAMS) end=1;

	/*	for (i=parameter+1;i<MAX_PARAMS;i++)
		{
		if (psr[0].param[i].paramSet[0]==1)
		{
		parameter=i;
		break;
		}
		}
		if (i==MAX_PARAMS) end=1; */
      }
  } while (end==0 && parameter<MAX_PARAMS);
  printf("chisq = %f\n",psr[0].fitChisq);
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
      if (ixv<-ilim || ixv>ilim) /* Use exponent notation */
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

char * plugVersionCheck = TEMPO2_h_VER;
