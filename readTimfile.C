#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"

/* ******************************************** */
/* read_timfile                                 */
/* Author:  G. Hobbs (02 May 2003)              */
/* Purpose: Reads a .tim file and fills         */
/*          observation structure               */
/* Inputs:  timFile - filename, obsn - structure*/
/*          of observations, nObs number of     */
/*          lines in .tim file                  */
/* Outputs: Fills obsn and nObs                 */
/*                                              */
/* Notes:   Default is to use the new .tim      */
/*          format file, but can read the old   */
/*          Parkes and Jodrell tempo format.    */ 
/*                                              */
/* Changes:                                     */
/* ******************************************** */
void readTim(char *timname,pulsar *psr,int *jumpVal);
void removeCR2(char *str);

void readTimfile(pulsar *psr,char timFile[][MAX_FILELEN],int npsr)
{
  int p,i;
  int jumpVal=0;
  FILE *fin;
  if (debugFlag==1) printf("In reading tim file\n");

  for (p=0;p<npsr;p++)
    {

      psr[p].nobs=0;
      /*      jumpVal=psr[p].nJumps;
      for (i=0;i<jumpVal;i++)
      psr[p].jumpVal[i]=0.0; */
      jumpVal=0;

      if (psr[0].jboFormat!=0)
	readJBO_bat(timFile[p],&psr[p],p);
      else
	readTim(timFile[p],&psr[p],&jumpVal);
      if (debugFlag==1) printf("Checking for deleted points >%s<\n",psr[0].deleteFileName);

      /* Check for deleted points in separate file */
      if (strcmp(psr[0].deleteFileName,"NONE")!=0)
	{
	  if (debugFlag==1) printf("In checking for deleted points\n");
	  if (!(fin = fopen(psr[0].deleteFileName,"r")))
	    printf("Warning: unable to open %s\n",psr[0].deleteFileName);
	  else
	    {
	      longdouble input;
	      while (!feof(fin))
		{
		  char inputstr[1024];
		  if (fscanf(fin,"%s",inputstr)==1)
		    {
		      input = parse_longdouble(inputstr);
		      for (i=0;i<psr[p].nobs;i++)
			{
			  if (fabs(psr[p].obsn[i].sat-input)<0.001)
			    psr[p].obsn[i].deleted=1;
			}
		    }
		}
	    }
	}
    }
  if (debugFlag==1) printf("Leaving readTimfile\n");  
}

void readTim(char *timname,pulsar *psr,int *jumpVal)
{
  FILE *fin;
  char profileDir[MAX_STRLEN]="";
  char tt[MAX_STRLEN];
  int nread=0,nread2,nObs=0,i,k;
  char firstWord[1000],line[1000]="",dummy[1000];
  char param1[100];//,param2[100],param3[100],param4[100],param5[100];
  //char param6[100],param7[100],param8[100];
  //  int val1;
  int format,endit=0;
  int valid;
  int skip=0;
  static double global_efac=1.0;
  static double efac=1.0;
  static double eset=-1.0;
  static double emax=-1.0;  /* Maximum error used */
  static double efloor=-1.0;
  static double emin=-1.0;
  static double equad = 0.0;
  static double fmax=-1.0;
  static double fmin=-1.0;
  static double sigma=0.0;  /* Forced setting of TOA error */
  static double time=0.0;
  int    add=0;
  static int infoNum=-1;
  char oldLine[1000];
  char global_efacFlag[100][1000];
  char global_efacFlagID[100][1000];
  double global_efacFlagVal[100];
  int nglobal_efacFlag=0;
  char efacFlag[100][1000];
  char efacFlagID[100][1000];
  double efacFlagVal[100];
  int  nefacFlag=0;

  char equadFlag[100][1000];
  char equadFlagID[100][1000];
  double equadFlagVal[100];
  int  nequadFlag=0;
  //  printf("Reading file %s %d %d\n",timname,nObs,psr[0].nobs);
  /* Attempt to open .tim file, exit if not possible */
  if (!(fin = fopen(timname,"r")))
    {
      printf("ERROR [FILE4]: Unable to open timfile >%s<\n",timname);
      exit(1);
    }
  /* Have successfully opened the file */
  
  /* Attempt to determine the file format */
  if (psr->fixedFormat > 0)
    {
      int k;
      format=1;
      for (k=0;k<psr->fixedFormat;k++)
	fgets(oldLine,1000,fin);
    }
  else
    {
      format=1;
      
      while (!feof(fin))
	{
	  fscanf(fin,"%s",firstWord);
	  if (strcmp(firstWord,"FORMAT")==0) /* Tempo2 format */
	    format = 0;
	  else if (strcmp(firstWord,"HEAD")==0) /* .tpo format */
	    format = 2;
	}
      fclose(fin);
      fin = fopen(timname,"r");  /* Could use rewind instead */
    }

  if (format==2)
    {
      do {
	fgets(line,1000,fin); removeCR2(line);
      } while (strcasecmp(line,"TOAS")!=0 && feof(fin)==0);
      format=1;

    }

  while (!feof(fin) && endit==0)
    {      
      valid=0;
      if (format==0) /* New Tempo2 (pseudo-)free format file format */
	{
	  if (nread<5 && strlen(line)>1)
	    {
	      line[strlen(line)-1]='\0';
	      /* printf("Warning: >%s< Not a complete line in .tim file. Ignoring line\n",line); */
	    }
	  nObs = psr->nobs;
	  if (fgets(line,1000,fin)==NULL)             /* Read line from .tim file */
	    valid=-2;
	  //	  printf("Have read %s\n",line);

	  /* Remove \n character at end */
	  if (line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0';
	  if (line[0]!='C' && line[0]!='#') /* Ignore comment lines     */
	    {
	      /* Read default columns */
	      if (strlen(line)>0)
		{
		  char sat_str[1024];
		  psr->obsn[nObs].deleted=0;

		  if (line[0]=='I')
		    {
		      nread = sscanf(line,"%s %s %lf %s %lf %s",dummy,psr->obsn[nObs].fname,
				     &(psr->obsn[nObs].freq),sat_str,
				     &(psr->obsn[nObs].toaErr),psr->obsn[nObs].telID);	      
		      psr->obsn[nObs].deleted=1;
		    }		  
		  else
		    {
		      nread = sscanf(line,"%s %lf %s %lf %s",psr->obsn[nObs].fname,
				 &(psr->obsn[nObs].freq),sat_str,
				     &(psr->obsn[nObs].toaErr),psr->obsn[nObs].telID);
		      if (strlen(profileDir)>0)
			{
			  sprintf(tt,"%s/%s",profileDir,psr->obsn[nObs].fname);
			  strcpy(psr->obsn[nObs].fname,tt);
			}
		    }
		  psr->obsn[nObs].sat = parse_longdouble(sat_str);
		  psr->obsn[nObs].phaseOffset = 0.0;
		  /* Read the rest of the line */		  
		  psr->obsn[nObs].nFlags = 0;
    
		  /*		  strcpy(psr->obsn[nObs].flagID[0],"FLAGID");
				  strcpy(psr->obsn[nObs].flagVal[0],"FLAGVAL"); */
		  for (i=0;i<(int)strlen(line)-1;i++)
		    {
		      if (line[i]=='-' && (line[i+1] < 48 || line[i+1] > 57))
			{
			  strcpy(oldLine,line);
			  if (strchr(line+i,' ')!=NULL)
			    {
			      strcpy(strchr(line+i,' '),"");
			      strcpy(psr->obsn[nObs].flagID[psr->obsn[nObs].nFlags],line+i);
			      i+=strlen(line+i)+1;
			      strcpy(line,oldLine);
			      if (strchr(line+i,' ')!=NULL)
				{
				  int j;
				  for (j=i;j<strlen(line);j++)
				    {
				      if (line[j]!=' ')
					{
					  i+=(j-i);
					  break;
					}
				    }
				  strcpy(strchr(line+i,' '),"");
				}
			      strcpy(psr->obsn[nObs].flagVal[psr->obsn[nObs].nFlags],line+i);
			      i+=strlen(line+i);
			      strcpy(line,oldLine); 
			      psr->obsn[nObs].nFlags++;
			    }
			}
		    }
		}
	      else
		nread=0;
	      
	      /* Now check whether this is a valid observation */
	      if (nread>=5 && valid==0) valid=1;
	    }
	  else
	    nread=5;	      
	}
      else if (format==1) /* Parkes, princeton or ITOA file format */
	{
	  if (fgets(line,1000,fin)==NULL)             /* Read line from .tim file */
	    valid=-2;
	  else
	    {	    
	      /* Remove \n character at end */
	      if (line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0';
	      nObs = psr->nobs;
	      psr->obsn[nObs].nFlags = 0;
	      nread = sscanf(line,"%s",param1);
	      add=0;

	      if (strcmp(param1,"C")==0 || strcmp(param1,"c")==0) /* Comment line */
		{
		  valid=0;
		  /*		  add=1; */ /* NOW IGNORE COMMENTED LINES COMPLETELY */
		  /*		  psr->obsn[nObs].deleted=1; */
		}
	      else
		{
		  psr->obsn[nObs].deleted=0;
		  
		  if (strcasecmp(param1,"INCLUDE")!=0)
		    {
		      if (line[0+add]==' ')      /* Parkes format */
			{
			  valid=1;
			  if (strlen(line)+add < 79) valid=-2;
			  else
			    {
			      strcpy(psr->obsn[nObs].fname,line+1+add); psr->obsn[nObs].fname[25]='\0';
			      strcpy(param1,line+25+add); param1[9]='\0'; if (strlen(param1)<2) valid=-2; 
			      
			      if (sscanf(param1,"%lf",&(psr->obsn[nObs].freq))!=1) valid=-2;
			      strcpy(param1,line+34+add); param1[21]='\0';
			      psr->obsn[nObs].sat = parse_longdouble(param1);
			      strcpy(param1,line+55+add); param1[8]='\0';
			      if (sscanf(param1,"%lf",&(psr->obsn[nObs].phaseOffset))!=1) valid=-2;
			      strcpy(param1,line+63+add); param1[8]='\0';
			      if (sscanf(param1,"%lf",&(psr->obsn[nObs].toaErr))!=1) valid=-2;
			      sscanf(line+79+add,"%s",psr->obsn[nObs].telID); psr->obsn[nObs].telID[1]='\0';
			    }
			}
		      else if (line[1+add]==' ') /* Princeton format */
			{
			  double dmoffset;
			  valid=1;
			  sscanf(line+add,"%s",psr->obsn[nObs].telID); psr->obsn[nObs].telID[1]='\0';	
			  strcpy(psr->obsn[nObs].fname,"NOT SET");
			  strcpy(param1,line+15+add); param1[9]='\0';
			  if (sscanf(param1,"%lf",&(psr->obsn[nObs].freq))!=1) valid=-2;
			  strcpy(param1,line+25+add); param1[20]='\0';
			  psr->obsn[nObs].sat = parse_longdouble(param1);
			  strcpy(param1,line+45+add); param1[9]='\0';
			  if (sscanf(param1,"%lf",&(psr->obsn[nObs].toaErr))!=1) valid=-2;
			  strcpy(param1,line+68+add); param1[10]='\0'; /* SHOULD BE DM OFFSET */
			  if (sscanf(param1,"%lf",&dmoffset)==1)
			    {
			      strcpy(psr->obsn[nObs].flagID[psr->obsn[nObs].nFlags],"-dmo");
			      sprintf(psr->obsn[nObs].flagVal[psr->obsn[nObs].nFlags],"%.10g",(double)(dmoffset)); 
			      psr->obsn[nObs].nFlags++;
			    }
			}
		      else if (strlen(line)>14 && line[14+add]=='.')  /* ITOA format */
			{
			  valid=1;
			  strcpy(psr->obsn[nObs].fname,"NOT SET");		      
			  strcpy(param1,line+9+add); param1[19]='\0';
			  psr->obsn[nObs].sat = parse_longdouble(param1);
			  strcpy(param1,line+28+add); param1[6]='\0';
			  if (sscanf(param1,"%lf",&(psr->obsn[nObs].toaErr))!=1) valid=-2;
			  strcpy(param1,line+34+add); param1[11]='\0';
			  if (sscanf(param1,"%lf",&(psr->obsn[nObs].freq))!=1) valid=-2;
			  strcpy(param1,line+45+add); param1[10]='\0'; /* SHOULD BE DM OFFSET */
			  if (sscanf(param1,"%lf",&(psr->obsn[nObs].phaseOffset))!=1) valid=-2;
			  sscanf(line+57+add,"%s",psr->obsn[nObs].telID); psr->obsn[nObs].telID[2]='\0';
			}
		    }
		}
	    }
	} // end of "else if (format == 1)"
      if (valid==1)  /* Check for validity of this TOA */
	{
	  /* The above comment is misleading. "valid==1" implies that the TOA _is_ valid. 
	     However, below, the observation undergoes a series of actions that 
	     - Should be performed on all observations 
	     - Are totally independent of the file format.
	     
	     Just so you know. 
	  */
	  if (psr->obsn[nObs].telID[0]=='@'            /* Have barycentric arrival time */
	      || strcasecmp(psr->obsn[nObs].telID,"bat")==0) 
	    {
	      psr->obsn[nObs].clockCorr=0;  /* therefore don't do clock corrections */
	      psr->obsn[nObs].delayCorr=0;
	    }
	  else if (strcmp(psr->obsn[nObs].telID,"STL")==0)
	    {
	      psr->obsn[nObs].clockCorr=0;  /* don't do clock corrections */
	      psr->obsn[nObs].delayCorr=1;
	    }
	  else
	    {
	      psr->obsn[nObs].clockCorr=1;
	      psr->obsn[nObs].delayCorr=1;
	    }	  
	
	  /* Check for conditionals */


	  /* Note ordering: EQUAD BEFORE EFAC */
	  /* Joris is going to re-write the equad business now. So stay tuned. 
	     In case everything goes wrong, you can put the following two lines back again.
	     if (equad!=0.0)
	     psr->obsn[nObs].toaErr = sqrt(pow(psr->obsn[nObs].toaErr,2)+equad*equad); */
	  
	  psr->obsn[nObs].toaErr = sqrt(pow(psr->obsn[nObs].toaErr,2)+equad*equad); 
	  psr->obsn[nObs].equad = equad;

	  for (i=0;i<nequadFlag;i++)
	    {
	      for (k=0;k<psr->obsn[nObs].nFlags;k++)
		{
		  if (strcmp(psr->obsn[nObs].flagID[k],equadFlagID[i])==0)
		    {
		      if (strcmp(psr->obsn[nObs].flagVal[k],equadFlag[i])==0){
			psr->obsn[nObs].toaErr = 
			  sqrt(pow(psr->obsn[nObs].toaErr,2)+pow(equadFlagVal[i],2));
			psr->obsn[nObs].equad += equadFlagVal[i];
		      }
		    }
		}
	    }
	  //printf("%lg  ",psr->obsn[nObs].equad);

	  // Furthermore, I suspect that the "EQUAD on certain flags" should be put before the
	  // efac adding, but well, who am I anyway?
	  if (sigma!=0.0)
	    psr->obsn[nObs].toaErr = sigma;
	  psr->obsn[nObs].toaErr *= (efac*global_efac);
	  psr->obsn[nObs].efac = (efac*global_efac);
	  /*	  for (i=0;i<nequadFlag;i++)
	    {
	      for (k=0;k<psr->obsn[nObs].nFlags;k++)
		{
		  if (strcmp(psr->obsn[nObs].flagID[k],equadFlagID[i])==0)
		    {
		      if (strcmp(psr->obsn[nObs].flagVal[k],equadFlag[i])==0)
			psr->obsn[nObs].toaErr = sqrt(pow(psr->obsn[nObs].toaErr,2)+pow(equadFlagVal[i],2));
		    }
		}
		}*/
	  for (i=0;i<nefacFlag;i++)
	    {
	      for (k=0;k<psr->obsn[nObs].nFlags;k++)
		{
		  if (strcmp(psr->obsn[nObs].flagID[k],efacFlagID[i])==0)
		    {
		      if (strcmp(psr->obsn[nObs].flagVal[k],efacFlag[i])==0){
			psr->obsn[nObs].toaErr *= efacFlagVal[i];
			psr->obsn[nObs].efac *= efacFlagVal[i];
		      }
		    }
		}
	    }
	  if (eset>-1)
	    psr->obsn[nObs].toaErr = eset;
  
	  psr->obsn[nObs].jump = *jumpVal;
	
	  if (time!=0.0) psr->obsn[nObs].sat += time/60.0/60.0/24.0;
	  if (efloor!=-1 && psr->obsn[nObs].toaErr < efloor) psr->obsn[nObs].toaErr = efloor;
	  if (emax!=-1 && psr->obsn[nObs].toaErr > emax) psr->obsn[nObs].deleted = 1;
	  if (emin!=-1 && psr->obsn[nObs].toaErr < emin) psr->obsn[nObs].deleted = 1;
	  if (fmax!=-1 && psr->obsn[nObs].freq > fmax) valid=0;
	  if (fmin!=-1 && psr->obsn[nObs].freq < fmin) valid=0;
	  if (infoNum!=-1)
	    {
	       strcpy(psr->obsn[nObs].flagID[psr->obsn[nObs].nFlags],"-i");
	       sprintf(psr->obsn[nObs].flagVal[psr->obsn[nObs].nFlags],"%d",infoNum);
	       psr->obsn[nObs].nFlags++;
	    }
	  if (skip==1) valid=0;
	  
	  if (valid==1)(psr->nobs)++;
	  if (psr->nobs > MAX_OBSN-2)
	    {
	      fprintf(stderr, "Too many TOAs! Change tempo.h and recompile!\n");
	      exit(1);
	    }
	}
      else if (valid!=-2) // This means: if the line does contain information, but not an observation
	{
	  nread2 = sscanf(line,"%s",param1);
	  if (skip==0)
	    {
	      if (strcasecmp(param1,"END")==0)
		endit=1;	
	      /* Global error multiplying factor */
	      else if (strcasecmp(param1,"PROFILE_DIR")==0)
		sscanf(line,"%s %s",param1,profileDir);
	      else if (strcasecmp(param1,"GLOBAL_EFAC")==0)
		sscanf(line,"%s %lf",param1,&global_efac);
	      else if (strcasecmp(param1,"EFAC")==0)  /* Error multiplying factor */    
		sscanf(line,"%s %lf",param1,&efac);
	      else if (strcasecmp(param1,"EFLOOR")==0)  /* Minimum error                        */    
		sscanf(line,"%s %lf",param1,&efloor);
	      else if (strcasecmp(param1,"T2EFAC")==0) /* EFAC for given flag                    */
		{
		  sscanf(line,"%s %s %s %lf",param1,efacFlagID[nefacFlag],efacFlag[nefacFlag],
			 &efacFlagVal[nefacFlag]);
		  nefacFlag++;
		}
	      else if (strcasecmp(param1,"T2EQUAD")==0) /* EQUAD for given flag                  */
		{
		  sscanf(line,"%s %s %s %lf",param1,equadFlagID[nequadFlag],equadFlag[nequadFlag],
			 &equadFlagVal[nequadFlag]);
		  nequadFlag++;
		}
	      else if (strcasecmp(param1,"EMAX")==0)  /* Maximum error                           */
		sscanf(line,"%s %lf",param1,&emax);
	      else if (strcasecmp(param1,"EMIN")==0)  /* Minimum error                           */
		sscanf(line,"%s %lf",param1,&emin);
	      else if (strcasecmp(param1,"ESET")==0)  /* Set all errors to given amount          */
		sscanf(line,"%s %lf",param1,&eset);
	      else if (strcasecmp(param1,"FMAX")==0)  /* Maximum observing frequency             */
		sscanf(line,"%s %lf",param1,&fmax);
	      else if (strcasecmp(param1,"FMIN")==0)  /* Minimum observing frequency             */
		sscanf(line,"%s %lf",param1,&fmin);
	      else if (strcasecmp(param1,"INFO")==0)  /* Highlighting flag                       */
		{
		  sscanf(line,"%s %d",param1,&infoNum);
		}
	      else if (strcasecmp(param1,"PHASE")==0) /* Add phase jump                          */
		{
		  psr->phaseJump[psr->nPhaseJump] = psr->obsn[nObs-1].sat;		  
		  sscanf(line,"%s %d",param1,&psr->phaseJumpDir[psr->nPhaseJump]);
		  psr->nPhaseJump++;		  
		}
	      else if (strcasecmp(param1,"EQUAD")==0) /* Error to add in quadrature               */
		sscanf(line,"%s %lf",param1,&equad);
	      else if (strcasecmp(param1,"SIGMA")==0) /* Set all errors to constant value         */
		sscanf(line,"%s %lf",param1,&sigma);
	      else if (strcasecmp(param1,"TIME")==0)  /* Add a constant time to all arrival times */
		{
		  double dtime;
		  sscanf(line,"%s %lf",param1,&dtime);
		  time+=dtime;
		}
	      else if (strcasecmp(param1,"MODE")==0) /* Fit with errors */
		{
		  sscanf(line,"%s %d",param1,&(psr->fitMode));
		  displayMsg(1,"TIM1","Please place MODE flags in the parameter file","",psr->noWarnings);
		}
	      else if (strcasecmp(param1,"INCLUDE")==0) /* Include another .tim file */
		{
		  char newtim[MAX_FILELEN];
		  if (sscanf(line,"%s %s",param1,newtim)==2)
		    readTim(newtim,psr,jumpVal);
		  else
		    printf("Unable to parse INCLUDE line >%s<\n",line);
		  strcpy(param1,"");
		}
	      else if (strcasecmp(param1,"SKIP")==0) /* Skip data */
		skip=1;
	      else if (strcasecmp(param1,"JUMP")==0) /* JUMP */
		(*jumpVal)++;
	      else if (strcasecmp(param1,"NOSKIP")==0) /* Stop skipping data */
		skip=0;
	    }
	  else
	    {
	      if (strcasecmp(param1,"NOSKIP")==0) /* Stop skipping data */
		skip=0;
	    }
	} // end of "else if (valid != -2)"
    } // end of "while (!feof(fin) && endit == 0)"
  fclose(fin);
}


/* Write an arrival time file in the tempo2 format */
void writeTim(char *timname,pulsar *psr,char *fileFormat)
{
  FILE *fout;
  int i,j;
  char name[1000];
  double current_efac;
  double current_equad;
  double interim_error;

  fout = fopen(timname,"w");
  if (strcmp(fileFormat,"tempo2")==0) fprintf(fout,"FORMAT 1\n");
  if (psr->fitMode==1) fprintf(fout,"MODE 1\n");    
  current_efac = psr[0].obsn[0].efac;
  current_equad = psr[0].obsn[0].equad;
  if(current_equad!=0.0)
    fprintf(fout,"EQUAD %lf\n",current_equad);
  if(current_efac!=1.0)
    fprintf(fout,"EFAC %lf\n",current_efac);

  for (i=0;i<psr->nobs;i++)
    {
      if((double)psr[0].obsn[i].equad !=  current_equad){
	fprintf(fout,"EQUAD %lf\n",psr[0].obsn[i].equad);
	current_equad = psr[0].obsn[i].equad;
      }
      if(psr[0].obsn[i].efac!=current_efac){
	fprintf(fout,"EFAC %lf\n",psr[0].obsn[i].efac);
	current_efac = psr[0].obsn[i].efac;
      }
      if (psr[0].obsn[i].deleted==1)
	fprintf(fout,"C");

      if (strcmp(fileFormat,"tempo2")==0) 
	{
	  sscanf(psr->obsn[i].fname,"%s",name);
	  interim_error = psr->obsn[i].toaErr/current_efac;
	  interim_error = sqrt(pow(interim_error,2.0)-(pow(current_equad,2.0)));
	  fprintf(fout," %s %.8f %.17Lf %.5f %s ", name,psr->obsn[i].freq,
		  // psr->obsn[i].sat, (psr->obsn[i].toaErr/current_efac),psr->obsn[i].telID);
		  psr->obsn[i].sat, interim_error,psr->obsn[i].telID);
	  if(interim_error<0.0L){
	    printf("ERROR - TOAerror < 0!!\n");
	    exit(1);
	  }

	  /* Now add flags */
	  for (j=0;j<psr->obsn[i].nFlags;j++)
	    {
	      if (strcasecmp(psr->obsn[i].flagID[j],"-selectAll")!=0)
		fprintf(fout,"%s %s ",psr->obsn[i].flagID[j],psr->obsn[i].flagVal[j]);
	    }
	  fprintf(fout,"\n");
	}
      else
	{
	  if (strlen(psr[0].obsn[i].fname) > 22)
	    psr[0].obsn[i].fname[22]='\0';
	    
	  fprintf(fout," %-25.25s%8.3f  %.13Lf    0.00 %7.2f        %s",
		  psr[0].obsn[i].fname,(double)psr[0].obsn[i].freq,
		  psr[0].obsn[i].sat,(double)psr[0].obsn[i].toaErr,psr[0].obsn[i].telID);
	  fprintf(fout,"\n");
	}

      //      for (j=0;j<psr->nPhaseJump;j++)
      //	{
      //	  if (i<psr->nobs && psr->obsn[i].sat < psr->phaseJump[j] && 
      //	      psr->obsn[i+1].sat >= psr->phaseJump[j])
      //	    fprintf(fout,"PHASE %d\n",psr->phaseJumpDir[j]);
      //	}
      
    }
  fclose(fout);
}

/* Removes newline at end of string */
void removeCR2(char *str)
{
  if (str[strlen(str)-1]=='\n') str[strlen(str)-1] = '\0';
}
