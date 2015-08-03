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

// redwards code for computing tropospheric delays

#include <math.h>
#include <strings.h>
#include "tempo2.h"
#include <stdio.h>
#include <stdlib.h>
#include <glob.h>
#include <string.h>
#include "tabulatedfunction.h"

// Note, these have been tested & differ at the 10^-3.5 (absolute) level
// from the results of the CALC/SOLVE routine. A major part of that is
// due to the incorrect length of year in the NMF formula, combined with
// a different choice of epoch, but smaller discrepancies remain and
// ought to be chased up if these routines are used for VLBI.
// for tempo2 there's no problem, we only need the mapping function to 0.1
double
NMF_hydrostatic(double utc_mjd,
		double site_latitude, double site_height,
		double source_elevation)
{
  // averages as a function of latitude
  double avgs_a[] = {1.2769934e-3, 
		    1.2683230e-3, 
		    1.2465397e-3, 
		    1.2196049e-3, 
		    1.2045996e-3};
  double avgs_b[] = {2.9153695e-3, 
		    2.9152299e-3, 
		    2.9288445e-3, 
		    2.9022565e-3, 
		    2.9024912e-3};
  double avgs_c[] = {62.610505e-3, 
		    62.837393e-3,
		    63.721774e-3,
		    63.824265e-3, 
		    64.258455e-3};
  double amps_a[] = {0.0,   
		    1.2709626e-5, 
		    2.6523662e-5, 
		    3.4000452e-5, 
		    4.1202191e-5};
  double amps_b[] = {0.0,
		    2.1414979e-5, 
		    3.0160779e-5, 
		    7.2562722e-5, 
		    11.723375e-5};
  double amps_c[] = {0.0,
		    0.0128400e-5, 
		    4.3497037e-5, 
		    84.795348e-5, 
		     170.37206e-5};
  // Height correction coefficients
  double a_h = 2.53e-5, b_h = 5.49e-3, c_h = 1.14e-3;
  //
  int ilat1, ilat2; // nearest 2 latitude marks
  double absolute_latitude;
  double frac=0.0;      // fraction of distance from lat 1 to lat 2
  double a, b, c;
  double cos_phase;
  double basic, height_correction;
  double sin_elevation=sin(source_elevation);

  // work out which latitude(s) to use
  absolute_latitude = fabs(site_latitude)*180.0/M_PI;
  ilat1 = (int)floor(absolute_latitude / 15.0)-1;
  if (ilat1 < 0)
    ilat1 = ilat2 = 0;
  else if (ilat1>=4)
    ilat1 = ilat2 = 4;
  else
  {
    frac = absolute_latitude / 15.0 - 1.0 - ilat1;
    ilat2=ilat1+1;
  }

  // work out the phase wrt DOY 28 of 2005 (or any other year), mjd 53398
  cos_phase = cos( (utc_mjd - 53398.0)
		   * 2.0*M_PI / 365.25)
    * (site_latitude < 0.0 ? -1.0 : 1.0); // Invert seasons if needed

  // get interpolated values for a, b, c
  if (ilat1==ilat2)
  {
    a = avgs_a[ilat1] - amps_a[ilat1] * cos_phase;
    b = avgs_b[ilat1] - amps_b[ilat1] * cos_phase;
    c = avgs_c[ilat1] - amps_c[ilat1] * cos_phase;
  }
  else
  {
    a = frac*avgs_a[ilat2] + (1.0-frac)*avgs_a[ilat1]
      - (frac*amps_a[ilat2] + (1.0-frac)*amps_a[ilat1]) * cos_phase;
    b = frac*avgs_b[ilat2] + (1.0-frac)*avgs_b[ilat1]
      - (frac*amps_b[ilat2] + (1.0-frac)*amps_b[ilat1]) * cos_phase;
    c = frac*avgs_c[ilat2] + (1.0-frac)*avgs_c[ilat1]
      - (frac*amps_c[ilat2] + (1.0-frac)*amps_c[ilat1]) * cos_phase;
  }

  // get the basic mapping function value (no height correction)
  basic = (1.0+a/(1.0+b/(1.0+c))) /
    (sin_elevation+a/(sin_elevation+b/(sin_elevation+c)));
  // get height correction
  // how well do we need to know this ???
  // max value is about 0.02 * height(km)
  // we need m.f. to about 0.1 .. so only need height to 5 km.
  // ellipsoid vs sphere = 20 km or so 

  // source elevation ... sin(elevation) = cos(zenith) 
  // where zenith is the angle between the obs-geocentre and source-geocentre
  // vectors ... roughly! actually this is not quite true as the zenith
  // is not aligned with the obs-geocentre vector due to ellipticity of earth.
  // Max deviation is 11.5' at around 45 deg.
  // At 5 deg elevation this amounts to about 0.5 in the mapping function,
  // which is too much! So, we need the geodetic latitude in order to
  // compute the elevation ... then we need to transform the zenith
  // vector into the celestial frame where we have the source direction
  // then finally we can do our dot product
// geoid vs ellipsoid = 100 m
  // so yes we do need the height above the reference ellipsoid
  height_correction = site_height * 1.0e-3 * 
    (1.0/sin_elevation -
     (1.0+a_h/(1.0+b_h/(1.0+c_h))) /
     (sin_elevation+a_h/(sin_elevation+b_h/(sin_elevation+c_h))));

  return basic + height_correction;
}

double
NMF_wet(double site_latitude, double source_elevation)
{
  double as[] = {5.8021897e-4, 
	       5.6794847e-4, 
	       5.8118019e-4, 
	       5.9727542e-4, 
	       6.1641693e-4};
  double bs[] = {1.4275268e-3, 
	       1.5138625e-3, 
	       1.4572752e-3, 
	       1.5007428e-3, 
	       1.7599082e-3};
  double cs[] = {4.3472961e-2, 
	       4.6729510e-2, 
	       4.3908931e-2, 
	       4.4626982e-2, 
	       5.4736038e-2};

  int ilat1, ilat2; // nearest 2 latitude marks
  double absolute_latitude;
  double frac=0.0;      // fraction of distance from lat 1 to lat 2
  double sin_elevation=sin(source_elevation);
  double a,b,c;
  // work out which latitude(s) to use
  absolute_latitude = fabs(site_latitude)*180.0/M_PI;
  ilat1 = (int)floor(absolute_latitude / 15.0)-1;
  if (ilat1 < 0)
    ilat1 = ilat2 = 0;
  else if (ilat1>=4)
    ilat1 = ilat2 = 4;
  else
  {
    frac = absolute_latitude / 15.0 - 1.0 - ilat1;
    ilat2=ilat1+1;
  }
  // get interpolated values for a, b, c
  if (ilat1==ilat2)
  {
    a = as[ilat1];
    b = bs[ilat1];
    c = cs[ilat1];
  }
  else
  {
    a = frac*as[ilat2] + (1.0-frac)*as[ilat1];
    b = frac*bs[ilat2] + (1.0-frac)*bs[ilat1];
    c = frac*bs[ilat2] + (1.0-frac)*cs[ilat1];
  }

  return (1.0+a/(1.0+b/(1.0+c))) /
    (sin_elevation+a/(sin_elevation+b/(sin_elevation+c)));
}

// Tables of Surface Atmpospheric Pressures and Zenith Wet Delays

typedef struct
{
  TabulatedFunction table;
  char siteName[256];
} MeteorologyFunction;

static int meteorology_tables_initialized = 0;
DynamicArray zenithWetDelayTables;
DynamicArray surfaceAtmosphericPressureTables;

void
MeteorologyFunction_load(MeteorologyFunction *func, char *fileName)
{
  TabulatedFunction_load(&func->table, fileName);
  if (sscanf(func->table.header_line+1, // skip # 
	     "%s", &func->siteName)!=1)
  {
    fprintf(stderr, 
	    "Error parsing meterology file %s: first line must be of form # site_name\n",
	    fileName);
    exit(1);
  }
}

double
MeteorologyFunction_getValue(MeteorologyFunction *func,
				      double mjd)
{
  return TabulatedFunction_getValue(&func->table, mjd);
}

double
MeteorologyFunction_getStartMJD(MeteorologyFunction *func)
{
  return TabulatedFunction_getStartX(&func->table);
}
double
MeteorologyFunction_getEndMJD(MeteorologyFunction *func)
{
  return TabulatedFunction_getEndX(&func->table);
}

void
initialize_meteorology_table(int dispWarnings,
			      char *path, char *extension,
			      DynamicArray *tables,
			      char *description)
{
  glob_t g;
  char pattern[1024];
  char **pfname;
  int globRet;

  MeteorologyFunction func;
  DynamicArray_init(tables, sizeof(MeteorologyFunction));

  /* load all  files in specified spot */
  sprintf(pattern, "%s/%s/*.%s", getenv(TEMPO2_ENVIRON), path, extension);
  globRet = glob(pattern, 0, NULL, &g);
  if (globRet == GLOB_NOSPACE)
    { printf("Out of memory in tropo.C\n"); exit(1);}
#ifdef GLOB_ABORTED
  else if (globRet == GLOB_ABORTED)
    { printf("Read error in tropo.C\n"); exit(1); }
#endif
#ifdef GLOB_NOMATCH
  else if (globRet == GLOB_NOMATCH)
    { if (dispWarnings==0)
      printf("No .%s files in $TEMPO2/%s\n", extension, path); return; }
#endif

  for (pfname = g.gl_pathv; *pfname != NULL; pfname++)
  {
    MeteorologyFunction_load(&func, *pfname);
    DynamicArray_push_back(tables, &func);
    if (dispWarnings==0)
      printf("Loaded %s for site %s from %s\n", 
	     description, func.siteName, func.table.fileName);
  }
}
  
void
initialize_meteorology_tables(int dispWarnings)
{
  initialize_meteorology_table(dispWarnings,
			       "atmosphere", "zwd", 
			       &zenithWetDelayTables,
			       "zenith wet delays");
  initialize_meteorology_table(dispWarnings,
			       "atmosphere", "sap", 
			       &surfaceAtmosphericPressureTables,
			       "surface atmospheric pressures");
  meteorology_tables_initialized = 1;
}

double getMeteorologicalValue(DynamicArray *tables,
			      char *siteName, double mjd, int warnings)
{
  if (meteorology_tables_initialized != 1)
    initialize_meteorology_tables(warnings);

  // search for that site name in that MJD range
  size_t itab;
  MeteorologyFunction *func;
  for (itab=0; itab < tables->nelem; itab++)
  {
    func = ((MeteorologyFunction *)tables->data) + itab;
    if (!strcasecmp(func->siteName, siteName)
	&& MeteorologyFunction_getStartMJD(func) <= mjd
	&& MeteorologyFunction_getEndMJD(func) >= mjd)
      break;
  }
  if (itab < tables->nelem)
    return MeteorologyFunction_getValue(func, mjd);
  else
    return 0.0;
}

double getZenithWetDelay(char *siteName, double mjd, int warnings)
{
  double zwd = getMeteorologicalValue(&zenithWetDelayTables, siteName, mjd,
				      warnings);
  if (zwd==0.0)
  {
    char msg[1000],msg2[1000];
    sprintf(msg,"Assume zero zenith wet delay (no data) for site %s  at MJD",
	    siteName);
    sprintf(msg2,"%.1f",mjd);
    displayMsg(1,"TROP1",msg,msg2,warnings);
  }

  return zwd;
}

double getSurfaceAtmosphericPressure(char *siteName, double mjd, int warnings)
{
  double sap = getMeteorologicalValue(&surfaceAtmosphericPressureTables, 
				      siteName, mjd,
				      warnings);
  if (sap==0.0)
  {
    char msg[1000],msg2[1000];
    sprintf(msg,"Assume standard atmospheric pressure (no data) for site %s at MJD",
	    siteName);
    sprintf(msg2,"%.1f",mjd);
    displayMsg(1,"TROP2",msg,msg2,warnings);
    sap = 101.325;
  }

  return sap;
}


void
compute_tropospheric_delays(pulsar *psr,int npsr)
{
  double zenith_delay_hydrostatic, zenith_delay_wet;
  double mapping_hydrostatic, mapping_wet;
  int i, p;
  observatory *obs;
  double source_elevation;
  double pressure;
  const char *CVS_verNum = "$Revision: 1.9 $";

  if (displayCVSversion == 1) CVSdisplayVersion("tropo.C","computer_tropospheric_delays()",CVS_verNum);

  // test code
#if 0
  {
    int ilat, iel;
    double lat, el, h, w;
    for (ilat=0;ilat<=11;ilat++)
    {
      lat = (ilat * 8.0)* 3.14159265358979/180.0;
      for (iel=0; iel<=100; iel++)
      {
	el = (5.0+iel*0.84) * 3.14159265358979/180.0;
	h = NMF_hydrostatic(50000.0, lat, 10000.0, el);
	w = NMF_wet(lat, el);
	printf("%.15lg %.15lg %.15lg %.15lg NMF\n",
	       el, h, w, lat);
      }
    }
    exit(1);
  }
#endif

  for (p=0;p<npsr;p++)
  {
    bool warned = false;
    for (i=0;i<psr[p].nobs;i++)
    { 
      psr[p].obsn[i].troposphericDelay = 0.0;
      if (psr[p].obsn[i].delayCorr!=0 && psr[p].correctTroposphere!=0 && psr[p].obsn[i].zenith[2]==0.0 && !warned)
      {
	logdbg( "WARNING: Tropospheric delay correction not possible with T2CMETHOD TEMPO %d %lf", i, psr[p].obsn[i].zenith[2]);
	warned = true;
      }
      else if (strcmp(psr[p].obsn[i].telID,"STL_FBAT")==0)
	{
	  logdbg("Not correcting for tropospheric delay");
	  warned=true;
	}
      else if (psr[p].obsn[i].delayCorr!=0 && psr[p].correctTroposphere!=0)
      {
	obs = getObservatory(psr[p].obsn[i].telID);
	/* Check whether the observation represents the COE */
	if (strcasecmp(obs->name,"COE")!=0)
	  {
	    // get source elevation neglecting proper motion
	    source_elevation = asin(dotproduct(psr[p].obsn[i].zenith,
					       psr[p].posPulsar)
				    / obs->height_grs80);
	    // get surface atmospheric pressure
	    pressure = 
	      getSurfaceAtmosphericPressure(obs->code, psr[p].obsn[i].sat, 
					    psr[p].noWarnings);
	    
	    // ------------------- Hydrostatic delay
	    // Zenith delay from Davies et al (Radio Sci 20 1593)
	    zenith_delay_hydrostatic = 0.02268 * pressure 
	      / (SPEED_LIGHT * (1.0-0.00266*cos(obs->latitude_grs80)
				-2.8e-7*obs->height_grs80));
	    // mapping function
	    mapping_hydrostatic = 
	      NMF_hydrostatic(psr[p].obsn[i].sat
			      + getCorrection(psr[p].obsn+i, 
					      psr[p].clockFromOverride,
					      "UTC", psr[p].noWarnings)/SECDAY,
			      obs->latitude_grs80, obs->height_grs80,
			      source_elevation);
	    // ------------------ Wet delay
	    zenith_delay_wet = 
	      getZenithWetDelay(obs->code, psr[p].obsn[i].sat, 
				psr[p].noWarnings);
	    mapping_wet = NMF_wet(obs->latitude_grs80, source_elevation);
	    
	    psr[p].obsn[i].troposphericDelay = 
	      zenith_delay_hydrostatic * mapping_hydrostatic
	      + zenith_delay_wet * mapping_wet;
	  }
	else
	  psr[p].obsn[i].troposphericDelay = 0.0;
      }
    }
  }
}

  
