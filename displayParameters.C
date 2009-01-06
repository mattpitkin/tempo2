#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"

/* ******************************************** */
/* displayParameters                            */
/* Author:  G. Hobbs (13 June 2003)             */
/* Purpose: Displays the values of all params   */
/*                                              */
/* Notes:                                       */
/* Changes:                                     */
/* ******************************************** */
void displayParameters(int pos,char timFile[][MAX_FILELEN],char parFile[][MAX_FILELEN],pulsar *psr,int npsr)
{
  int i, ic, p;
  
  printf("\n\n");
  if (pos==1) printf("Parameter values after getInputs from command line\n");
  if (pos==2) printf("Parameter values after reading of .par file\n");
  if (pos==3) printf("Parameter values after reading of .tim file\n");
  if (pos==4) printf("After obtaining correction between observatory UTC and UTC(NIST)\n");
  if (pos==5) printf("After calculating correction between UTC and TAI\n");
  if (pos==6) printf("After calculating correction between TAI and TB\n");
  if (pos==7) printf("After calculating vector pointing at pulsar\n");

  printf("\n\n");
  printf("lobal definitions:\n\n");
  printf("MAX_FREQ_DERIVATIVES = %d\n",MAX_FREQ_DERIVATIVES);
  printf("MAX_DM_DERIVATIVES   = %d\n",MAX_DM_DERIVATIVES);
  printf("MAX_FILELEN          = %d\n",MAX_FILELEN);
  printf("\n\n");
  for (p=0;p<npsr;p++)
    {
      printf("Command line inputs:\n\n");
      printf("Pulse arrival time file[0] = %s\n",timFile[0]);
      printf("Pulsar parameter file[0]   = %s\n",parFile[0]);
      
      printf("\n\n");
      printf("Pulsar Parameters:\n\n");
      printf("Pulsar name          = %s\n",psr[p].name);
      if (pos>2) /* Have .tim file data */
	{
	  printf("\n\nObservations:\n\n");
	  printf("Number of observations = %d\n",psr[p].nobs);
	  printf("--------------------------------------------------------------------------------\n");
	  printf("   #  | SAT                   | Clock corrections(s)     |\n");
	  printf("--------------------------------------------------------------------------------\n");
	  for (i=0;i<psr[p].nobs;i++)
	  {
	    printf("|%4d | %s |", i+1, print_longdouble(psr[p].obsn[i].sat).c_str());
	    for (ic=0; ic < psr[p].obsn[i].nclock_correction; ic++)
	      printf(" %-14.7g(->%s)", 
		     (double)psr[p].obsn[i].correctionsTT[ic].correction,
		     psr[p].obsn[i].correctionsTT[ic].corrects_to);
	    printf("\n");
	  }
 	  printf("--------------------------------------------------------------------------------\n");
	  printf("   #   TAI->TB   | SAT->UT1      |\n");
	  printf("--------------------------------------------------------------------------------\n");
	  for (i=0;i<psr[p].nobs;i++)
	    printf("|%4d | %-10.9g | %-16.9e |\n",
		   i+1,
		   (double)psr[p].obsn[i].correctionTT_TB,
		   (double)psr[p].obsn[i].correctionUT1);
	  printf("--------------------------------------------------------------------------------\n");
	}
      if (pos>6)  /* Ephemeris information */
	{
	  printf("\n\nEphemeris information:\n\n");
	  printf("3-vector pointing at pulsar       = (%-13.13f,%-13.13f,%13.13f)\n",
		 psr[p].posPulsar[0],psr[p].posPulsar[1],psr[p].posPulsar[2]);
	  printf("3-vector giving pulsar's velocity = (%-13.8e,%13.8e,%13.8e)\n",
		 psr[p].velPulsar[0],psr[p].velPulsar[1],psr[p].velPulsar[2]);
	}
      if (pos>7)
	{
	  printf("\n\n");
	  printf("--------------------------------------------------------------------------------------------------------------\n");
	  printf("   JD            RCS(1)   RCS(2)     RCS(3)      RCB(1)     RCB(2)     RCB(3)    RBE(1)    RBE(2)    RBE(3)\n");
	  printf("--------------------------------------------------------------------------------------------------------------\n");
	  for (i=0;i<psr[p].nobs;i++)
	    {
	      printf("%-.5f %10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f\n",
		     (double)(psr[p].obsn[i].sat + getCorrectionTT(psr[p].obsn+i)/SECDAY +
		     psr[p].obsn[i].correctionTT_TB/SECDAY 
		     + 2400000.5),
		     (double)psr[p].obsn[i].sun_ssb[0],
		     (double)psr[p].obsn[i].sun_ssb[1],
		     (double)psr[p].obsn[i].sun_ssb[2],
		     (double)psr[p].obsn[i].earthMoonBary_ssb[0],
		     (double)psr[p].obsn[i].earthMoonBary_ssb[1],
		     (double)psr[p].obsn[i].earthMoonBary_ssb[2],
		     (double)psr[p].obsn[i].earthMoonBary_earth[0],
		     (double)psr[p].obsn[i].earthMoonBary_earth[1],
		     (double)psr[p].obsn[i].earthMoonBary_earth[2]);
	    }
	  printf("--------------------------------------------------------------------------------------------------------------\n");
	  
	  printf("\n\n");
	  printf("--------------------------------------------------------------------------------------------------------------\n");
	  printf("   JD            RCS(4)   RCS(5)     RCS(6)      RCB(4)     RCB(5)     RCB(6)    RBE(4)    RBE(5)    RBE(6)\n");
	  printf("--------------------------------------------------------------------------------------------------------------\n");
	  for (i=0;i<psr[p].nobs;i++)
	    {
	      printf("%-.5lf %10.3e %-10.3e %-10.3e %-10.3e %-10.3e %-10.3e %-10.3e %-10.3e %-10.3e\n",
		     (double)(psr[p].obsn[i].sat + 
		     getCorrectionTT(psr[p].obsn+i)/SECDAY + 
		      psr[p].obsn[i].correctionTT_TB/SECDAY 
		     + 2400000.5),
		     (double)psr[p].obsn[i].sun_ssb[3],
		     (double)psr[p].obsn[i].sun_ssb[4],
		     (double)psr[p].obsn[i].sun_ssb[5],
		     (double)psr[p].obsn[i].earthMoonBary_ssb[3],
		     (double)psr[p].obsn[i].earthMoonBary_ssb[4],
		     (double)psr[p].obsn[i].earthMoonBary_ssb[5],
		     (double)psr[p].obsn[i].earthMoonBary_earth[3],
		     (double)psr[p].obsn[i].earthMoonBary_earth[4],
		     (double)psr[p].obsn[i].earthMoonBary_earth[5]);
	    }
	  printf("--------------------------------------------------------------------------------------------------------------\n");
	  printf("MJD         OBS->EARTH1   OBS->EARTH2  OBS->EARTH3  SITEVEL1    SITEVEL2     SITEVEL3\n");
	  printf("--------------------------------------------------------------------------------------------------------------\n");

	  for (i=0;i<psr[p].nobs;i++)
	    {
	      printf("%-.5lf %12.5e %-12.5e %-12.5e %-12.5e %-12.5e %-12.5e\n",
		     (double)psr[p].obsn[i].sat + 
		     getCorrectionTT(psr[p].obsn+i)/SECDAY,
		     psr[p].obsn[i].observatory_earth[0],psr[p].obsn[i].observatory_earth[1],
		     psr[p].obsn[i].observatory_earth[2],
		     psr[p].obsn[i].siteVel[0],psr[p].obsn[i].siteVel[1],psr[p].obsn[i].siteVel[2]);
	    }
	  printf("--------------------------------------------------------------------------------------------------------------\n");

	}
      if (pos>9) /* Extra Delays */
	{
	  printf("\n\nExtra Delays:\n\n");
	  
	  printf("--------------------------------------------------------------------------------------------------------------\n"); 
	  printf("MJD          ShapiroSun          TDIS               Roemer           BAT\n");
	  printf("--------------------------------------------------------------------------------------------------------------\n");
	  for (i=0;i<psr[p].nobs;i++)
	    {
	      printf("%s %.12e %.12e %.10Lf %.10f\n",print_longdouble(psr[p].obsn[i].sat).c_str(),psr[p].obsn[i].shapiroDelaySun,
		     psr[p].obsn[i].tdis1+psr[p].obsn[i].tdis2,psr[p].obsn[i].roemer,(double)psr[p].obsn[i].bat);
	    }
	  printf("--------------------------------------------------------------------------------------------------------------\n");
    }
    }

}

