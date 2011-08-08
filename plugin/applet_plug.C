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

/* This plugin allows the user to define the output format for listing observational based
 * values
 *
 * On the command line: use -s "....." option
 * In a file: use -file filename
 *
 * options {resPre} - prefit residuals
 *         {resPost} - postfit residuals
 * 
 * 
 * 
 *
 * Other commands
 *
 * {ERRMULT x} multiple all errors by x
 * {FORMAT c}  sets the format for displaying values and errors (printf-type format)
 * {TAB x}     move the cursor to position x
 * {NULL c}    set the null string to display if a parameter is not set
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
void parseLine(pulsar *psr,char *line,double *errMult,char *null,char *format,char *dformat,int *rad, FILE *fout);
double fortranMod(double a,double p);
void help() /* Display help */
{
	  /* This function should contain usage information about the plugin which should (in general) be accessed */
	  /* by the user pressing 'h'                                                                              */
}


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr)
{
FILE *fin, *fout;
  char filename[100], filenameout[100];
  char line[500];
  char var[500],type[100],cline[1000];
  char null[500]="*";
  char format[500]="%.20Lg";
  char dformat[500]="%.20lg";
  char cval[500],cerr[500],msg[500];
  char negative,header;
  int  i,j,k,varN,file=-1,fileout=0,rad=0;
  double errMult=1;
  double mjd;


  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  double globalParameter;

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: name\n");
  printf("Author:              author\n");
  printf("Version:             version number\n");
  printf(" --- type 'h' for help information\n");


  /* Obtain the .par and the .tim file from the command line */
  if (argc==4) /* Only provided .tim name */
  {
	  strcpy(timFile[0],argv[3]);
	  strcpy(parFile[0],argv[3]);
	  parFile[0][strlen(parFile[0])-3] = '\0';
	  strcat(parFile[0],"par");
  }

  for (i=0;i<argc;i++)
  {
	  if (strcmp(argv[i],"-f")==0)
	  {
		  strcpy(parFile[0],argv[i+1]);
		  strcpy(timFile[0],argv[i+2]);
	  }
  }

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);

  for (i=0;i<argc;i++)
  {
	  if ((strcasecmp(argv[i],"-jp")==0) || (strcasecmp(argv[i],"-jn")==0)){
		  negative = (strcasecmp(argv[i],"-jn")==0);
		  psr[0].phaseJump[psr[0].nPhaseJump] = atof(argv[++i]);
		  if (negative) psr[0].phaseJumpDir[psr[0].nPhaseJump] = -1;
		  else psr[0].phaseJumpDir[psr[0].nPhaseJump] = 1;
		  /*		  for (j=0;j<psr[0].nobs;j++)
				  {
				  k=psr[0].nPhaseJump;
				  if (psr[0].obsn[j].sat > psr[0].phaseJump[k])
				  {
				  psr[0].obsn[j].residual += (double)psr[0].phaseJumpDir[k]/psr[0].param[param_f].val[0];
				  psr[0].obsn[j].prefitResidual += (double)psr[0].phaseJumpDir[k]/psr[0].param[param_f].val[0];
				  }
				  }*/
		  psr[0].nPhaseJump++;
	  }

	  if (strcasecmp(argv[i],"-start")==0){
		  psr[0].param[param_start].val[0]=atof(argv[++i]);
		  psr[0].param[param_start].paramSet[0]=1;
		  psr[0].param[param_start].fitFlag[0]=1;
	  }

	  if (strcasecmp(argv[i],"-finish")==0){
		  psr[0].param[param_finish].val[0]=atof(argv[++i]);
		  psr[0].param[param_finish].paramSet[0]=1;
		  psr[0].param[param_finish].fitFlag[0]=1;
	  }
  }
  printf("Inserted %d phase jumps\n",psr[0].nPhaseJump);

  header=0;

  for (i=0;i<2;i++)                   /* Do two iterations for pre- and post-fit residuals*/
  {
	  formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
	  formResiduals(psr,*npsr,1);    /* Form the residuals                 */
	  if (i==0) doFit(psr,*npsr,0);   /* Do the fitting     */
	  else textOutput(psr,*npsr,0,0,0,1,"out.par");  /* Display the output */
  }

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
      if (strcasecmp(argv[i],"-outfile")==0)
	{
	  fileout=1;
	  strcpy(filenameout,argv[i+1]);
	}
      if (strcasecmp(argv[i],"-head")==0)
      {
	      header=1;
      }

   }
  if (file==-1)
    {
      printf("must use -s or -file option\n");
      exit(1);
    }
  if (fileout==1) /* Write to file */
    {
      if (!(fout = fopen(filenameout,"w")))
	{
	  printf("Unable to open file %s\n",filenameout);
	  exit(1);
	}
    }
  else if (fileout==0) /* Write to standard out */
    {
      fout = stdout;
    }  

  
  if(header){
	  fprintf(fout,"HDR\tP0\t%llf\n",1.0/psr[0].param[param_f].val[0]);
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
	      parseLine(psr,line,&errMult,null,format,dformat,&rad,fout);
	    }
	}
    }
  else if (file==0) /* Read from string */
    parseLine(psr,cline,&errMult,null,format,dformat,&rad,fout);
}

void parseLine(pulsar *psr,char *line,double *errMult,char *null,char *format,char *dformat,int *rad,FILE *fout)
{
  int i,j;
  char var[500],type[100];
  int varN,lv,le,parameter=0,end=1;
  char cval[500],cerr[500],msg[500];
  char disp[1000];
  int pos=0;
  int sub1=0,first=0;

  varN=-1;
  do {
    pos=0;
    end=0;
    if (varN==psr[0].nobs-2) end=1;
    varN++;

    if (psr[0].obsn[varN].deleted==0)
      {
	if (first==0) first = varN; /* Record first residual */

	for (i=0;i<strlen(line);i++)
	  {	  
	    if (line[i]=='\\' && (line[i+1]=='{' || line[i+1]=='}'))
	      { /* Do nothing */ }
	    else if (line[i]=='\\' && line[i+1]=='n')
	      {
		i++;
		fprintf(fout,"\n");
	      }
	    else if (line[i]=='\\' && line[i+1]=='t')
	      {
		i++;
		sprintf(disp,"\t");
		fprintf(fout,"%s",disp);
		pos+=strlen(disp);
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
		else if (strstr(var,"SUB1")!=NULL || strstr(var,"sub1")!=NULL) /* Subtract first residual from remainder */
		    sub1 = 1;
		else if (strstr(var,"TAB")!=NULL)
		  {
		    int reqSpace,k;
		    sscanf(var+4,"%d",&reqSpace);
		    if (reqSpace>pos)
		      {
			for (k=pos;k<reqSpace;k++)
			  fprintf(fout," ");
			pos=reqSpace;
		      }
		  }
		else if (strstr(var,"NORAD")!=NULL)
		  *rad=0;
		else if (strstr(var,"RAD")!=NULL)
		  *rad=1;
		else if (strstr(var,"ERRMULT")!=NULL || strstr(var,"errmult")!=NULL)
		  sscanf(var+7,"%lf",errMult);
		else if (strstr(var,"NULL")!=NULL || strstr(var,"null")!=NULL)
		  sscanf(var+5,"%s",null);
		
		/*	    var[j-i-3]='\0'; */
		
		
		/*	    strcpy(type,line+i+1+j-i-2);*/
		i=j;
		type[1]='\0';
		
		if (varN>=0)
		  {
		    if (strcasecmp(var,"bat")==0) /* barycentric arrival time */
		      {
			sprintf(disp,format,psr[0].obsn[varN].bat);
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		    if (strcasecmp(var,"bbat")==0) /* binary barycentric arrival time */
		      {
			sprintf(disp,format,psr[0].obsn[varN].bbat);
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		    if (strcasecmp(var,"clock4")==0) /* 4th clock correction */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].correctionsTT[4].correction); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		    if (strcasecmp(var,"ut1")==0) /* UT1 correction */
		      {
			sprintf(disp,dformat,(double)psr[0].obsn[varN].correctionUT1);
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		    if (strcasecmp(var,"clock3")==0) /* 3th clock correction */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].correctionsTT[3].correction);  
			fprintf(fout,"%s",disp);
			pos+=strlen(disp); 
		      }
		    if (strcasecmp(var,"clock2")==0) /* 2th clock correction */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].correctionsTT[2].correction); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		    if (strcasecmp(var,"clock1")==0) /* 1th clock correction */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].correctionsTT[1].correction); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		    if (strcasecmp(var,"del")==0) 
		      {
			sprintf(disp,dformat,(double)psr[0].obsn[varN].deleted); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		    if (strcasecmp(var,"file")==0)
		      {
			sprintf(disp,"%s",psr[0].obsn[varN].fname); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		    if (strcasecmp(var,"clock0")==0) /* 0th clock correction */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].correctionsTT[0].correction); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		    if (strcasecmp(var,"sat")==0) /* site arrival time */
		      {
			sprintf(disp,format,psr[0].obsn[varN].sat); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		    else if (strcasecmp(var,"shapiro")==0) /* Solar shapiro delay */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].shapiroDelaySun); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"shapiroJ")==0) /* Solar shapiro delay */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].shapiroDelayJupiter); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"shapiroS")==0) /* Solar shapiro delay */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].shapiroDelaySaturn); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"shapiroV")==0) /* Solar shapiro delay */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].shapiroDelayVenus); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"shapiroU")==0) /* Solar shapiro delay */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].shapiroDelayUranus); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"shapiroN")==0) /* Solar shapiro delay */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].shapiroDelayNeptune); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"roemer")==0) /* Roemer delay */
		      {
			sprintf(disp,format,psr[0].obsn[varN].roemer); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"flags")==0) 
		      {
			int l;
			for (l=0;l<psr[0].obsn[varN].nFlags;l++)
			  {
			    sprintf(disp,"%s %s ",psr[0].obsn[varN].flagID[l],psr[0].obsn[varN].flagVal[l]); 
			    fprintf(fout,"%s",disp);
			    pos+=strlen(disp);			
			  }
		      }
		    else if (strcasecmp(var,"tropo")==0)
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].troposphericDelay); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"tt")==0) /* Correction to TT */
		      {
			sprintf(disp,dformat,getCorrectionTT(psr[0].obsn+varN)/SECDAY); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"tt2tb")==0) /* Correction from TT to tb */
		      {
			sprintf(disp,format,psr[0].obsn[varN].correctionTT_TB); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"earth_ssb")==0) /* Vector from Earth to SSB */
		      {
			sprintf(disp,dformat,sqrt(pow(psr[0].obsn[varN].earth_ssb[0],2)+
						 pow(psr[0].obsn[varN].earth_ssb[0],2)+
						 pow(psr[0].obsn[varN].earth_ssb[0],2))); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"earth_ssb1")==0) /* x from Earth to SSB */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].earth_ssb[0]); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"earth_ssb2")==0) /* y from Earth to SSB */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].earth_ssb[1]); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"earth_ssb3")==0) /* z from Earth to SSB */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].earth_ssb[2]); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"sun_earth1")==0) 
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].sun_earth[0]); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"sun_earth2")==0)
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].sun_earth[1]); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"sun_earth3")==0) 
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].sun_earth[2]); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    /* partial derivative w.r.t any parameter (post-fit) */
		    else if (strncasecmp(var,"d_",2)==0)
		    {
		      int iparam, iindex;
		      for (iparam=0; iparam < MAX_PARAMS; iparam++)
		      {
			for (iindex=0; iindex < psr->param[iparam].aSize; iindex++)
			  if (!strcasecmp(psr->param[iparam].shortlabel[iindex], var+2))
			    break;
			if (iindex < psr->param[iparam].aSize)
			  break;
		      }
		      
		      if (iparam==MAX_PARAMS)
		      {
			fprintf(stderr, "Parameter %s not found!\n",
				var+2);
			exit(1);
		      }

		      sprintf(disp,dformat, 
			      getParamDeriv(psr,varN,psr[0].obsn[varN].bbat-psr[0].param[param_pepoch].val[0], iparam, iindex)); 
		      fprintf(fout,"%s",disp);
		      pos+=strlen(disp);			
		     
		    }
		    else if (strcasecmp(var,"ism")==0) /* ISM delay */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].tdis1); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"elev")==0) /* Sourve elevation neglecting proper motion */
		      {
			observatory *obs;
			double source_elevation;

			obs = getObservatory(psr[0].obsn[varN].telID);
			// get source elevation neglecting proper motion
			source_elevation = asin(dotproduct(psr[0].obsn[varN].zenith,
							   psr[0].posPulsar)
						/ obs->height_grs80);
			sprintf(disp,dformat,source_elevation*180.0/M_PI); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"posPulsar")==0)
		      {
			sprintf(disp,"%f %f %f",psr[0].posPulsar[0],psr[0].posPulsar[1],psr[0].posPulsar[2]); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"zenith")==0)
		      {
			sprintf(disp,"%f %f %f",psr[0].obsn[varN].zenith[0],psr[0].obsn[varN].zenith[1],
				psr[0].obsn[varN].zenith[2]); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"npulse")==0) /* Pulse number */
		      {
			longdouble phase;
			sprintf(disp,"%.0Lf",psr[0].obsn[varN].phase); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"clock")==0)
		      {
			sprintf(disp,format,getCorrectionTT(psr[0].obsn+varN)+psr[0].obsn[varN].correctionTT_TB-psr[0].obsn[varN].correctionsTT[0].correction-psr[0].obsn[varN].correctionsTT[1].correction-psr[0].obsn[varN].correctionsTT[2].correction); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"clk_corr1")==0)
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].correctionsTT[0].correction); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"clk_corr2")==0)
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].correctionsTT[1].correction); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"clk2")==0)
		      {
			sprintf(disp,dformat,(psr[0].obsn[varN].correctionsTT[0].correction+psr[0].obsn[varN].correctionsTT[1].correction+psr[0].obsn[varN].correctionsTT[2].correction)/SECDAY); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }		    
		    else if (strcasecmp(var,"clkchain")==0)
		      {
			int l;
			strcpy(disp,"");
			for (l=0;l<psr[0].obsn[varN].nclock_correction;l++)
			  {
			    if (l!=0)
			      strcat(disp,"->");
			    strcat(disp,psr[0].obsn[varN].correctionsTT[l].corrects_to);
			  }			    
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }		    
		    else if (strcasecmp(var,"ipm")==0) /* Interplanetary medium delay */
		      {
			sprintf(disp,dformat,psr[0].obsn[varN].tdis2); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    else if (strcasecmp(var,"fmjdu")==0) 
		      {
			sprintf(disp,format,psr[0].obsn[varN].sat+(psr[0].obsn[varN].correctionsTT[0].correction+psr[0].obsn[varN].correctionsTT[1].correction+psr[0].obsn[varN].correctionsTT[2].correction)/SECDAY); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);			
		      }
		    if (strcasecmp(var,"x")==0) /* BAT - epoch */
		      {
			sprintf(disp,format,psr[0].obsn[varN].bat - psr[0].param[param_pepoch].val[0]); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		    if (strcasecmp(var,"freq")==0) /* observing frequency */
		      {
			sprintf(disp,format,(longdouble)psr[0].obsn[varN].freq); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		    if (strcasecmp(var,"pre")==0) /* prefit residual */
		      {
			sprintf(disp,format,(longdouble)(psr[0].obsn[varN].prefitResidual-sub1*psr[0].obsn[first].prefitResidual)); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		    if (strcasecmp(var,"pre_phase")==0) /* prefit residual in phase */
		      {
			sprintf(disp,format,(longdouble)(psr[0].obsn[varN].prefitResidual-sub1*psr[0].obsn[first].prefitResidual)*psr[0].param[param_f].val[0]); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		    if (strcasecmp(var,"post_phase")==0) /* postfit residual in phase */
		      {
			sprintf(disp,format,(longdouble)(psr[0].obsn[varN].residual-sub1*psr[0].obsn[first].residual)*psr[0].param[param_f].val[0]); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		    if (strcasecmp(var,"post")==0) /* postfit residual */
		      {
			sprintf(disp,format,(longdouble)(psr[0].obsn[varN].residual-sub1*psr[0].obsn[first].residual)); 
			fprintf(fout,"%s",disp);
			pos+=strlen(disp);
		      }
		if (strcasecmp(var,"err")==0) /* toaErr */
		  {
		    sprintf(disp,format,(longdouble)psr[0].obsn[varN].toaErr); 
		    fprintf(fout,"%s",disp);
		    pos+=strlen(disp);
		  }
		if (strcasecmp(var,"binphase")==0) /* binary phase */
		  {
		    double pbdot=0.0;
		    double tpb = (psr[0].obsn[varN].bat-psr[0].param[param_t0].val[0])*86400.0/(psr[0].param[param_pb].val[0]*SECDAY);
		    double phase;
		    phase = fortranMod(tpb+1000000.0,1.0);
		    if (phase < 0.0) phase+=1.0; 


		    sprintf(disp,format,(longdouble)phase); 
		    fprintf(fout,"%s",disp);
		    pos+=strlen(disp);
		  }
		
	      }
	      }
	else
	  {
	    fprintf(fout,"%c",line[i]);	      
	    pos++;
	  }
	  }
      }
  } while (end==0);
}

int rnd8(double rval,double rerr,int ifac,char *cval,int *lv,char *cerr,int *le,char *msg)
{
  double vv, ee, xv, xe;
  int ixv, ixe, iee, j, ivv, ise, irnd, ilim,ret,ret_lv,ret_le;
  char cexp[9], fmt[12];

  ilim = 4;

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
	  irnd=nint_derived(pow(10.0,irnd));
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
      vv=rval*pow(10.0,ixv); /* Scale for mantissa */
      ee=ee*pow(10.0,ixe);   /* scale error */
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


double fortranMod(double a,double p)
{
  double ret;

  ret = a - (int)(a/p)*p;
  return ret;
}
char * plugVersionCheck = TEMPO2_h_VER;
