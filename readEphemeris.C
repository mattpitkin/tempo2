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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tempo2.h"
#include "jpleph.h"
#include "ifteph.h"

#define MAX_SHOTS 10000

double random2(long *idum);
double gasdev(long *idum);

/* ******************************************** */
/* readEphemeris                                */
/* Author:  G. Hobbs (05 May 2003)              */
/* Purpose: Reads and interpolates the JPL      */
/*          planetary ephemeris                 */
/* Inputs:                                      */
/* Outputs:                                     */
/*                                              */
/* Notes: This routine does not give exactly    */
/*        the same values as the old tempo      */
/*        THIS MUST BE CHECKED                  */
/* Changes:                                     */
/* ******************************************** */
void readEphemeris(pulsar *psr,int npsr,int addEphemNoise)
{
  double jd_teph[2]; /* Julian date = MJD - 2400000.5, in Eph frame (~=TDB) */
  longdouble jd;
  void *ephem;
  char nams[400][6];
  double emrat;
  double vals[400];
  double one_au;
  int err_code;
  int i,p;
  const char *CVS_verNum = "$Revision: 1.12 $";

  if (displayCVSversion == 1) CVSdisplayVersion("readEphemeris.C","readEphemeris()",CVS_verNum);

  if (addEphemNoise!=0)
    printf("WARNING: Adding noise into ephemeris %d\n",addEphemNoise); 

  //  printf("Reading the ephemeris\n");
  for (p=0;p<npsr;p++)
    {
      /* Must convert MJD to JD */
      /*      printf("Ephemeris reading = %s\n",psr[p].JPL_EPHEMERIS); */
      //      printf("Ephemeris being used = %s\n",psr[p].JPL_EPHEMERIS);
      ephem = jpl_init_ephemeris(psr[p].JPL_EPHEMERIS, nams, vals);
      
      if( !ephem)
	{
	  printf( "ERROR [EPHEM1]: Ephemeris file '%s' not loaded\n", psr[p].JPL_EPHEMERIS);
	  exit(1);
	}
      
      one_au = jpl_get_double(ephem, JPL_EPHEM_AU_IN_KM) * 1000.0;
      emrat  = jpl_get_double(ephem, JPL_EPHEM_EARTH_MOON_RATIO);
      if (debugFlag) printf("one_au = %g, emrat = %g\n",one_au,emrat);
      /* Get "Ephemeris AU" in SI m instead of Ephemeris m */
      if (psr[p].units == SI_UNITS)
	      one_au = (double)one_au*IFTE_K; 
      /* redwards shapiro delay checking hack */
#if 0
      {
	double radius[] = {3.440e6, 6.052e6, 6.378e6,
			  3.397e6, 7.149e7, 6.027e7, 2.556e7,
			  2.477e7, 1.150e6, 1.738e6, 6.95e8};
	double mass[] = {3.30e23, 4.87e24, 5.97e24,
			6.42e23, 1.90e27, 5.68e26, 8.68e25,
			1.02e26, 1.27e22, 7.35e22, 1.99e30};
	double shap_min[11], shap_max[11], shap2_min[11], shap2_max[11];
	double jd, jd0=2433265.5, jd1=2469807.5;
	int i;
	for (i=0; i < 11; i++)
	{
	  shap_min[i] = shap2_min[i] = 1e6;
	  shap_max[i] = shap2_max[i] = -1e6;
	}
	for (jd=jd0;jd < jd1; jd+=10.0)
	{
	  printf("%.12lf", jd);
	  for (i=0; i < 11; i++)
	  {
	    if (i!=2)
	    {
	    double r_schwarz = 6.67e-11*mass[i]/9e16;
	    double r_p[6];
	    r_p[0]=r_p[1]=r_p[2] = 0.0;
	    jpl_pleph(ephem, (double)jd, i+1, 3, r_p, 0);
	    double r = sqrt(r_p[0]*r_p[0]+
			    r_p[1]*r_p[1]+
			    r_p[2]*r_p[2]) * one_au * 1000.0;
	    double theta = radius[i]/ r;
	    double shaphi = -2.0*r_schwarz / 3e8 * log(r*(1-cos(theta)));
	    double shaplo = -2.0*r_schwarz / 3e8 * log(r*2.0);
	    double shap2hi = r_schwarz * r_schwarz / 3e8
	      * 4.0 / (r * theta*theta);
	    double shap2lo = 0.0;
// 	    printf (" %lg %lg %lg", theta, shaphi, shap2hi);
	    if (shaplo < shap_min[i]) 
	      shap_min[i] = shaplo;
	    if (shaphi > shap_max[i])
	      shap_max[i] = shaphi;
	    if (shap2lo < shap2_min[i])
	      shap2_min[i] = shap2lo;
	    if (shap2hi > shap2_max[i])
	      shap2_max[i] = shap2hi;
	  }
	  }
// 	  printf(" SHAP\n");
	}
	for (i=0; i < 11; i++)
	{
	  fprintf(stderr, "%d %lg %lg\n", i, shap_max[i]-shap_min[i],
		  shap2_max[i]-shap2_min[i]);
	}
// 	exit(1);
      }
#endif
      for (i=0;i<psr[p].nobs;i++)
	{	  
	  // If the arrival time is from a satellite whose coordinates are known in the
	  // barycentric reference frame then set these values directly without reading
	  // the ephemeris
	  /* Note, interpolation takes place within the JPL reader.  */
	  /* JPL ephemeris is based on ephemeris time */
	  jd = psr[p].obsn[i].sat + getCorrectionTT(psr[p].obsn+i)/SECDAY + 
	    psr[p].obsn[i].correctionTT_Teph/SECDAY+2400000.5; 
	  jd_teph[0] = (double)((int)jd); /* 2452620.0; */
	  jd_teph[1] = (double)(jd - (int)jd); /* 0.08342753346369;  */
	  /*	  if (psr[0].obsn[i].deleted==0) printf("Giving ephemeris reader: %.14Lf %.14f %.14f\n",jd,jd_teph[0],jd_teph[1]); */
	  
	  /* Convert to TDB if necessary */ 
	  // 	  if (psr[p].units == SI_UNITS)
	  // 	    jd_teph = (jd_teph - IFTE_JD0)/IFTE_K + IFTE_JD0;
	  // 	  /* Convert TDB to T_eph */
	  // 	  jd_teph += IFTE_TEPH0/86400.0;
	  
	  /* Calculate position of Sun from solar-system BC (SSB) */
	  err_code = jpl_pleph(ephem, jd_teph, 11, 12, psr[p].obsn[i].sun_ssb, 1);
	  /* Convert to sec and lt-s/s from AU and AU/day*/
	  psr[p].obsn[i].sun_ssb[0]*=one_au/SPEED_LIGHT; psr[p].obsn[i].sun_ssb[1]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].sun_ssb[2]*=one_au/SPEED_LIGHT;
	  psr[p].obsn[i].sun_ssb[3]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].sun_ssb[4]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].sun_ssb[5]*=one_au/86400.0/SPEED_LIGHT; 
	  /* Convert to ecliptic coordinates if necessary */
	  if (psr[p].eclCoord==1) equ2ecl(psr[p].obsn[i].sun_ssb);
	  
	  /* Calculate position of Sun from Earth */
	  err_code = jpl_pleph(ephem, jd_teph, 11, 3, psr[p].obsn[i].sun_earth, 1);
	  /* Convert to sec and lt-s/s from AU and AU/day*/
	  psr[p].obsn[i].sun_earth[0]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].sun_earth[1]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].sun_earth[2]*=one_au/SPEED_LIGHT;
	  psr[p].obsn[i].sun_earth[3]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].sun_earth[4]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].sun_earth[5]*=one_au/86400.0/SPEED_LIGHT; 
	  if (psr[p].eclCoord==1) equ2ecl(psr[p].obsn[i].sun_earth);
	  
	      
	  /* Calculate position of all planets from solar-system BC (SSB),
	     if necessary */
	  int iplanet;
	  for (iplanet=0; iplanet < 9; iplanet++)
	    {
	      if (psr[p].param[param_dmassplanet].paramSet[iplanet])
		{
		  err_code = jpl_pleph(ephem, jd_teph, iplanet+1, 12, 
				       psr[p].obsn[i].planet_ssb[iplanet], 1);
		  /* Convert to sec and lt-s/s from AU and AU/day*/
		  psr[p].obsn[i].planet_ssb[iplanet][0]*=one_au/SPEED_LIGHT; 
		  psr[p].obsn[i].planet_ssb[iplanet][1]*=one_au/SPEED_LIGHT; 
		  psr[p].obsn[i].planet_ssb[iplanet][2]*=one_au/SPEED_LIGHT;
		  psr[p].obsn[i].planet_ssb[iplanet][3]*=one_au/86400.0/SPEED_LIGHT; 
		  psr[p].obsn[i].planet_ssb[iplanet][4]*=one_au/86400.0/SPEED_LIGHT; 
		  psr[p].obsn[i].planet_ssb[iplanet][5]*=one_au/86400.0/SPEED_LIGHT; 
		  if (psr[p].eclCoord==1) 
		    equ2ecl(psr[p].obsn[i].planet_ssb[iplanet]);
		}
	    }
	  /* Actually should determine this position at the time the pulse passes the pulsar */
	  err_code = jpl_pleph(ephem, jd_teph, 5, 3, psr[p].obsn[i].jupiter_earth, 1);
	  /* Convert to sec and lt-s/s from AU and AU/day*/
	  psr[p].obsn[i].jupiter_earth[0]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].jupiter_earth[1]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].jupiter_earth[2]*=one_au/SPEED_LIGHT;
	  psr[p].obsn[i].jupiter_earth[3]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].jupiter_earth[4]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].jupiter_earth[5]*=one_au/86400.0/SPEED_LIGHT; 
	  /* Convert to ecliptic coordinates if necessary */
	  if (psr[p].eclCoord==1) equ2ecl(psr[p].obsn[i].jupiter_earth);

	  /* Actually should determine this position at the time the pulse passes the pulsar */
	  err_code = jpl_pleph(ephem, jd_teph, 6, 3, psr[p].obsn[i].saturn_earth, 1);
	  /* Convert to sec and lt-s/s from AU and AU/day*/
	  psr[p].obsn[i].saturn_earth[0]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].saturn_earth[1]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].saturn_earth[2]*=one_au/SPEED_LIGHT;
	  psr[p].obsn[i].saturn_earth[3]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].saturn_earth[4]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].saturn_earth[5]*=one_au/86400.0/SPEED_LIGHT; 
	  if (psr[p].eclCoord==1) equ2ecl(psr[p].obsn[i].saturn_earth);

	  err_code = jpl_pleph(ephem, jd_teph, 2, 3, psr[p].obsn[i].venus_earth, 1);
	  /* Convert to sec and lt-s/s from AU and AU/day*/
	  psr[p].obsn[i].venus_earth[0]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].venus_earth[1]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].venus_earth[2]*=one_au/SPEED_LIGHT;
	  psr[p].obsn[i].venus_earth[3]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].venus_earth[4]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].venus_earth[5]*=one_au/86400.0/SPEED_LIGHT; 
	  if (psr[p].eclCoord==1) equ2ecl(psr[p].obsn[i].venus_earth);

	  /* Actually should determine this position at the time the pulse passes the pulsar */
	  err_code = jpl_pleph(ephem, jd_teph, 7, 3, psr[p].obsn[i].uranus_earth, 1);
	  /* Convert to sec and lt-s/s from AU and AU/day*/
	  psr[p].obsn[i].uranus_earth[0]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].uranus_earth[1]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].uranus_earth[2]*=one_au/SPEED_LIGHT;
	  psr[p].obsn[i].uranus_earth[3]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].uranus_earth[4]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].uranus_earth[5]*=one_au/86400.0/SPEED_LIGHT; 
	  /* Convert to ecliptic coordinates if necessary */
	  if (psr[p].eclCoord==1) equ2ecl(psr[p].obsn[i].uranus_earth);

	  /* Actually should determine this position at the time the pulse passes the pulsar */
	  err_code = jpl_pleph(ephem, jd_teph, 8, 3, psr[p].obsn[i].neptune_earth, 1);
	  /* Convert to sec and lt-s/s from AU and AU/day*/
	  psr[p].obsn[i].neptune_earth[0]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].neptune_earth[1]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].neptune_earth[2]*=one_au/SPEED_LIGHT;
	  psr[p].obsn[i].neptune_earth[3]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].neptune_earth[4]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].neptune_earth[5]*=one_au/86400.0/SPEED_LIGHT; 
	  /* Convert to ecliptic coordinates if necessary */
	  if (psr[p].eclCoord==1) equ2ecl(psr[p].obsn[i].neptune_earth);

	  /* Calculate position of centre of Earth from the SSB */
	  err_code = jpl_pleph(ephem, jd_teph, 3, 12, psr[p].obsn[i].earth_ssb, 1);
	  /* Convert to sec and lt-s/s from AU and AU/day*/
	  psr[p].obsn[i].earth_ssb[0]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].earth_ssb[1]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].earth_ssb[2]*=one_au/SPEED_LIGHT;
	  psr[p].obsn[i].earth_ssb[3]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].earth_ssb[4]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].earth_ssb[5]*=one_au/86400.0/SPEED_LIGHT;

	  /* Redwards adjust SSB position for reflex motion due to error
	     in planetary masses */
	  for (iplanet=0; iplanet < 9; iplanet++)
	  {
	    if (psr[p].param[param_dmassplanet].paramSet[iplanet])
	    {
	      //	      	      printf("NEW PLANET DMASS: %g %g %g\n",(double)psr[p].param[param_dmassplanet].val[4],
	      //	      		     (double)psr[p].obsn[i].earth_ssb[0],(double)(psr[p].param[param_dmassplanet].val[4] *
	      //	      								  psr[p].obsn[i].planet_ssb[iplanet][0]));
		      for (int icomp=0; icomp < 6; icomp++)
			psr[p].obsn[i].earth_ssb[icomp] -= 
			  psr[p].param[param_dmassplanet].val[iplanet] *
			  psr[p].obsn[i].planet_ssb[iplanet][icomp];
	    }
	  }
;
	  /* Convert to ecliptic coordinates if necessary */
	  if (psr[p].eclCoord==1) 
	    {
	      equ2ecl(psr[p].obsn[i].earth_ssb);	      
	      equ2ecl(psr[p].obsn[i].earth_ssb+3);
	    }
	  
	  /* Calculate position of Earth-Moon BC from the SSB */
	  err_code = jpl_pleph(ephem, jd_teph, 13, 12, psr[p].obsn[i].earthMoonBary_ssb, 1);
	  /* Convert to sec and lt-s/s from AU and AU/day*/
	  psr[p].obsn[i].earthMoonBary_ssb[0]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].earthMoonBary_ssb[1]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].earthMoonBary_ssb[2]*=one_au/SPEED_LIGHT;
	  psr[p].obsn[i].earthMoonBary_ssb[3]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].earthMoonBary_ssb[4]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].earthMoonBary_ssb[5]*=one_au/86400.0/SPEED_LIGHT;
	  /* Convert to ecliptic coordinates if necessary */
	  if (psr[p].eclCoord==1) equ2ecl(psr[p].obsn[i].earthMoonBary_ssb);
	  

	  /* Calculate position of Moon from the Earth */

	  err_code = jpl_pleph(ephem, jd_teph, 13, 3, psr[p].obsn[i].earthMoonBary_earth, 1);
	  /* Convert to sec and lt-s/s from AU and AU/day*/
	  psr[p].obsn[i].earthMoonBary_earth[0]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].earthMoonBary_earth[1]*=one_au/SPEED_LIGHT; 
	  psr[p].obsn[i].earthMoonBary_earth[2]*=one_au/SPEED_LIGHT;
	  psr[p].obsn[i].earthMoonBary_earth[3]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].earthMoonBary_earth[4]*=one_au/86400.0/SPEED_LIGHT; 
	  psr[p].obsn[i].earthMoonBary_earth[5]*=one_au/86400.0/SPEED_LIGHT;

	  /* 	  err_code = jpl_pleph(ephem, jd_teph, 10, 3, psr[p].obsn[i].earthMoonBary_earth, 1); */
	  /*	  psr[p].obsn[i].earthMoonBary_earth[0]*=-one_au/(1.0+emrat)/SPEED_LIGHT; 
	  psr[p].obsn[i].earthMoonBary_earth[1]*=-one_au/(1.0+emrat)/SPEED_LIGHT; 
	  psr[p].obsn[i].earthMoonBary_earth[2]*=-one_au/(1.0+emrat)/SPEED_LIGHT;
	  psr[p].obsn[i].earthMoonBary_earth[3]*=-one_au/86400.0/(1.0+emrat)/SPEED_LIGHT; 
	  psr[p].obsn[i].earthMoonBary_earth[4]*=-one_au/86400.0/(1.0+emrat)/SPEED_LIGHT; 
	  psr[p].obsn[i].earthMoonBary_earth[5]*=-one_au/86400.0/(1.0+emrat)/SPEED_LIGHT; */
	  /* Convert to ecliptic coordinates if necessary */
	  if (psr[p].eclCoord==1) equ2ecl(psr[p].obsn[i].earthMoonBary_earth);

	  /* Calculate obsn[i].nutations */
	  err_code = jpl_pleph(ephem, jd_teph, 14, 3, psr[p].obsn[i].nutations, 1);
	  /* Note interpolation done in the JPL reading routines */
	  //	  printf("Checking for interpolation %d\n",psr[p].nTelDX);
	  if (psr[p].nTelDX > 0)
	    {
	      int k;
	      double m,c,ival;
	      ival=0.0;
	      for (k=0;k<psr[p].nTelDX-1;k++)
		{
		  if ((double)psr[p].obsn[i].sat >= psr[p].telDX_t[k] &&
		      (double)psr[p].obsn[i].sat < psr[p].telDX_t[k+1])
		    {
		      m = (psr[p].telDX_v[k]-psr[p].telDX_v[k+1])/(psr[p].telDX_t[k]-psr[p].telDX_t[k+1]);
		      c = psr[p].telDX_v[k]-m*psr[p].telDX_t[k];
		      
		      ival = m*(double)psr[p].obsn[i].sat+c;
		      break;
		    }
		}
	      psr[p].obsn[i].earth_ssb[0] += ival; // Change x-position
	    }
	  if (psr[p].nTelDY > 0)
	    {
	      int k;
	      double m,c,ival;
	      ival=0.0;
	      //	      printf("Using interpolation function for telDY\n");
	      for (k=0;k<psr[p].nTelDY-1;k++)
		{
		  if ((double)psr[p].obsn[i].sat >= psr[p].telDY_t[k] &&
		      (double)psr[p].obsn[i].sat < psr[p].telDY_t[k+1])
		    {
		      m = (psr[p].telDY_v[k]-psr[p].telDY_v[k+1])/(psr[p].telDY_t[k]-psr[p].telDY_t[k+1]);
		      c = psr[p].telDY_v[k]-m*psr[p].telDY_t[k];
		      
		      ival = m*(double)psr[p].obsn[i].sat+c;
		      break;
		    }
		}
	      psr[p].obsn[i].earth_ssb[1] += ival; // Change y-position
	    }
	  if (psr[p].nTelDZ > 0)
	    {
	      int k;
	      double m,c,ival;
	      ival=0.0;
	      //	      printf("Using interpolation function for telDZ\n");
	      for (k=0;k<psr[p].nTelDZ-1;k++)
		{
		  if ((double)psr[p].obsn[i].sat >= psr[p].telDZ_t[k] &&
		      (double)psr[p].obsn[i].sat < psr[p].telDZ_t[k+1])
		    {
		      m = (psr[p].telDZ_v[k]-psr[p].telDZ_v[k+1])/(psr[p].telDZ_t[k]-psr[p].telDZ_t[k+1]);
		      c = psr[p].telDZ_v[k]-m*psr[p].telDZ_t[k];
		      
		      ival = m*(double)psr[p].obsn[i].sat+c;
		      break;
		    }
		}
	      psr[p].obsn[i].earth_ssb[2] += ival; // Change z-position
	    }
	  if (strcmp(psr[p].obsn[i].telID,"STL_BAT")==0)
	    {
	      int k;
	      for (k=0;k<psr[p].obsn[i].nFlags;k++)
		{
		  if (strcmp(psr[p].obsn[i].flagID[k],"-telx")==0){
		    sscanf(psr[p].obsn[i].flagVal[k],"%lf",&psr[p].obsn[i].earth_ssb[0]);
		  }
		  if (strcmp(psr[p].obsn[i].flagID[k],"-tely")==0){
		    sscanf(psr[p].obsn[i].flagVal[k],"%lf",&psr[p].obsn[i].earth_ssb[1]);
		  }
		  if (strcmp(psr[p].obsn[i].flagID[k],"-telz")==0){
		    sscanf(psr[p].obsn[i].flagVal[k],"%lf",&psr[p].obsn[i].earth_ssb[2]);
		  }
		}
	      for (k=3;k<6;k++) psr[p].obsn[i].earth_ssb[k]=0; // Should set if we know the velocity
	    }
	}
      jpl_close_ephemeris(ephem); 
    }
}

/* Based on ran1.f in original fortran */

double random2(long *idum)
{
  int j;
  static longdouble r[100],result;
  long m1=259100;
  long ia1=7141;
  long ic1=54773;
  long m2=134456;
  long ia2=8121;
  long ic2=28411;
  long m3=243000;
  long ia3=4561;
  long ic3=51349;
  longdouble rm1,rm2;

  static int iff=0;
  static int ix1=0;
  static int ix2=0;
  static int ix3=0;

  rm1=1./m1;
  rm2=1./m2;
    
  if(*idum < 0 || iff == 0) 
    {
      iff=1;
      ix1=(int)fortran_mod(ic1-(*idum),m1);
      ix1=(int)fortran_mod(ia1*ix1+ic1,m1);
      ix2=(int)fortran_mod(ix1,m2);
      ix1=(int)fortran_mod(ia1*ix1+ic1,m1);
      ix3=(int)fortran_mod(ix1,m3);

      for (j=0;j<97;j++)
	{
	  ix1=(int)fortran_mod(ia1*ix1+ic1,m1);
	  ix2=(int)fortran_mod(ia2*ix2+ic2,m2);
	  r[j]=(ix1+ix2*rm2)*rm1;
	}
      *idum=1;
    }
  
  ix1=(int)fortran_mod(ia1*ix1+ic1,m1);
  ix2=(int)fortran_mod(ia2*ix2+ic2,m2);
  ix3=(int)fortran_mod(ia3*ix3+ic3,m3);
  j=(97*ix3)/m3;
  if(j>96 || j<0) {printf("Problem in bootstrap.C\n"); exit(1);}

  result = r[j];
  r[j]=(ix1+ix2*rm2)*rm1;
  return result;
}

double gasdev(long *idum)
{
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*random2(idum)-1.0;
      v2=2.0*random2(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}
