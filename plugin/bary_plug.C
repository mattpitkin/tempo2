//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

// Plugin to calculate barycentric arrival times given a site arrival
// time and a telescope site code or coordinate
//

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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "dynarr.h"

using namespace std;

#define GRS80_A 6378137.0           /* semi-major axis (m) */
#define GRS80_F 1.0/298.257222101   /* flattening */

void ITRF_to_GRS80(observatory *obs);

void help() /* Display help */
{
}

static DynamicArray observatories;

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  int i,j;
  double globalParameter;
  int defineMJD=0;
  long double mjd,bat;
  int defineTelCoord=0;
  double telx,tely,telz;
  int definePsrCoord=0;
  double raj, decj;
  double roemer,rcos1,shapiro;
  double rca[3];
  observatory *obs;
  observatory newObs;
  double clockCor;

  const char *CVS_verNum = "$Revision: 1.1 $";


  if (displayCVSversion == 1) CVSdisplayVersion((char *)"bary.C",(char *)"plugin",CVS_verNum);

  *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: bary\n");
  printf("Author:              G. Hobbs, D. Kaplan\n");
  printf("CVS Version:         $Revision: 1.1 $\n");
  printf(" --- type 'h' for help information\n");


  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-mjd")==0)
	{
	  defineMJD = 1;
	  sscanf(argv[++i],"%Lf,",&mjd);
	}
      else if (strcmp(argv[i],"-coord")==0)
	{
	  defineTelCoord = 1;
	  sscanf(argv[++i],"%lf",&telx);
	  sscanf(argv[++i],"%lf",&tely);
	  sscanf(argv[++i],"%lf",&telz);
	}
      else if (strcmp(argv[i],"-psr")==0)
	{
	  definePsrCoord = 1;
	  sscanf(argv[++i],"%lf",&raj);
	  sscanf(argv[++i],"%lf",&decj);
	}
    }

  if (defineMJD == 0)
    {
      printf("To use this plugin you have to define the site arrival time on the command line using the -mjd option\n");
      exit(1);
    }

  if (defineTelCoord==0)
    {
      printf("To use this plugin you have to define the telescope geocentric coordinates on the command line using the -psr raj decj option (in radians)\n");
      exit(1);
    }

  if (definePsrCoord==0)
    {
      printf("To use this plugin you have to define the telescope geocentric coordinates on the command line using the -coord x y z option\n");
      exit(1);
    }

  // Set up pulsar structure
  psr[0].param[param_raj].val[0] = (long double)raj;
  psr[0].param[param_raj].paramSet[0] = 1;
  psr[0].param[param_decj].val[0] = (long double)decj;
  psr[0].param[param_decj].paramSet[0] = 1;
  psr[0].obsn[0].sat = mjd;
  psr[0].obsn[0].clockCorr=1;
  psr[0].nobs = 1;
  psr[0].units = SI_UNITS;
  psr[0].timeEphemeris = FB90_TIMEEPH;
  psr[0].obsn[0].delayCorr = 1;
  psr[0].t2cMethod = T2C_IAU2000B;
  strcpy(psr[0].JPL_EPHEMERIS,getenv(TEMPO2_ENVIRON));
  strcpy(psr[0].ephemeris,"DE405");
  strcat(psr[0].JPL_EPHEMERIS,"/ephemeris/DE405.1950.2050");
  strcpy(psr[0].obsn[0].telID,"Unset");
  printf("Warning: defaulting to DE405 ephemeris\n");

  // Clock corrections
  //  toa2utc(psr,1);
  clockCor = getCorrection(&(psr[0].obsn[0]),(char *)"UTC",(char *)"TT(TAI)",1);
  psr[0].obsn[0].nclock_correction = 1;
  psr[0].obsn[0].correctionsTT[0].correction = clockCor;

  printf("Clock correction = %g\n",clockCor);
  // Form a vector pointing at the pulsar
  vectorPulsar(psr,1);

  // Read the Solar System ephemeris
  readEphemeris(psr,1,0);/* 2. Read the ephemeris */

  // Set up the observatory to centre-of-earth vector
  DynamicArray_init(&observatories, sizeof(observatory));
  newObs.x = telx;
  newObs.y = tely;
  newObs.z = telz;
  ITRF_to_GRS80(&newObs);

  //  psr[0].obsn[0].observatory_earth[0] = telx;
  //  psr[0].obsn[0].observatory_earth[1] = tely;
  //  psr[0].obsn[0].observatory_earth[2] = telz;

  double trs[3], zenith_trs[3];
  trs[0]=newObs.x;
  trs[1]=newObs.y;
  trs[2]=newObs.z;
  zenith_trs[0]= newObs.height_grs80 
    * cos(newObs.longitude_grs80) * cos(newObs.latitude_grs80);
  zenith_trs[1]= newObs.height_grs80 
    * sin(newObs.longitude_grs80) * cos(newObs.latitude_grs80);
  zenith_trs[2] = newObs.height_grs80*sin(newObs.latitude_grs80);
  long double utc = psr[0].obsn[0].sat;
  //  if (psr[p].obsn[i].clockCorr!=0 && psr[p].obsn[i].clockCorr!=2)
  //    utc += getCorrection(psr[p].obsn+i, 
  //			 psr[p].clockFromOverride,
  //			 "UTC", psr[p].noWarnings)/SECDAY;
  utc=psr[0].obsn[0].sat;
  get_obsCoord_IAU2000B(trs, zenith_trs,
			psr[0].obsn[0].sat,
			//			+getCorrectionTT(psr[0].obsn)/SECDAY,
			utc,
			psr[0].obsn[0].observatory_earth,
			psr[0].obsn[0].zenith,
			psr[0].obsn[0].siteVel);
  printf("obs eqrth = %g %g %g\n",psr[0].obsn[0].observatory_earth[0],psr[0].obsn[0].observatory_earth[1],
	 psr[0].obsn[0].observatory_earth[2]);

  



  // Roemer delay (ignoring parallax etc.)
  for (j=0;j<3;j++)
    rca[j] = psr[0].obsn[0].earth_ssb[j] + psr[0].obsn[0].observatory_earth[j]; 

  rcos1 = dotproduct(psr[0].posPulsar,rca);
  roemer = rcos1;

  // Shapiro delay
  shapiro_delay(psr,1,0,0,0,0); // Ignoring pulsar velocity term
  shapiro = psr[0].obsn[0].shapiroDelaySun;

  // Calculate TT to TB correction
  tt2tb(psr,1);


  bat = mjd + (long double)(roemer-shapiro+clockCor+psr[0].obsn[0].correctionTT_TB)/86400.0L;

  printf("\n\n");
  printf("site arrival time        = %.15Lf\n",mjd);
  printf("site to TT correction (s)= %g\n",clockCor);
  printf("TT to TB correction      = %g\n",(double)psr[0].obsn[0].correctionTT_TB);
  printf("Roemer delay (s)         = %g\n",roemer);
  printf("Shapiro delay (s)        = %g\n",shapiro);
  printf("barycentric arrival time = %.15Lf\n",bat);
  return 0;
}

char * plugVersionCheck = (char *)TEMPO2_h_VER;

// Geocentric to geodetic.
// Uses Vermeille (2004)'s method:
//http://www.springerlink.com/app/home/contribution.asp?wasp=08ea5d2c4c62464789a7961196d84ab5&referrer=parent&backto=issue,11,18;journal,9,85;linkingpublicationresults,1:100435,1
void ITRF_to_GRS80(observatory *obs)
{
  double p = (obs->x*obs->x + obs->y*obs->y)/ (GRS80_A*GRS80_A);
  double esq = GRS80_F*(2.0-GRS80_F);
  double q = (1.0-esq)/(GRS80_A*GRS80_A)*obs->z*obs->z;
  double r = (p+q-esq*esq)/6.0;
  double s = esq*esq*p*q/(4*r*r*r);
  double t = pow(1.0+s+sqrt(s*(2.0+s)), 1.0/3.0);
  double u = r*(1.0+t+1.0/t);
  double v = sqrt(u*u+esq*esq*q);
  double w = esq*(u+v-q)/(2.0*v);
  double k = sqrt(u+v+w*w)-w;
  double D = k*sqrt(obs->x*obs->x+obs->y*obs->y)/(k+esq);
  
  obs->height_grs80 = (k+esq-1.0)/k * sqrt(D*D+obs->z*obs->z);
  obs->latitude_grs80 = 2.0*atan2(obs->z, D+sqrt(D*D+obs->z*obs->z));
  if (obs->y >= 0.0)
    obs->longitude_grs80 =
      0.5*M_PI - 2.0*atan2(obs->x, sqrt(obs->x*obs->x+obs->y*obs->y)+obs->y);
  else
    obs->longitude_grs80 = 
      -0.5*M_PI + 2.0*atan2(obs->x, sqrt(obs->x*obs->x+obs->y*obs->y)-obs->y);
}
