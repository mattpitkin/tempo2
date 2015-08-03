
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

/*-*-C-*-*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <glob.h>
#include <math.h>

#include "tempo2.h"
#include "dynarr.h"

/* static database of observatories */

static DynamicArray observatories;
typedef struct
{
  char code[100];
  DynamicArray aliases;
} ObservatoryAliasList;

static DynamicArray observatoryAliasLists;

static int observatories_initialised = 0;



//void convertGeocentric(observatory newObs);
void GRS80_to_ITRF(observatory *obs);
void ITRF_to_GRS80(observatory *obs);

double fang(int i,double f);

void
readObservatoryFile(char *fname)
{
  FILE *f = fopen(fname, "r");
  char line[1024];
  observatory newObs;
  int nread, iline=0;
  const char *CVS_verNum = "$Revision: 1.9 $";

  if (displayCVSversion == 1) CVSdisplayVersion("observatory.C","observatory()",CVS_verNum);

  if (!f)
  {
    fprintf(stderr, "Unable to open observatory data file %s : %s\n",
	    fname, strerror(errno));
    exit(1);
  }

  while (fgets(line, 1024, f))
  {
    iline++;
    if (line[0]!='#')
    {
      nread = sscanf(line, "%lf %lf %lf %s %s %s",
		     &newObs.x, &newObs.y, &newObs.z, 
		     newObs.name, newObs.code, newObs.clock_name);
      if (nread > 0)
      {
	 /* Most likely given in geodetic coordinates */
	if (fabs(newObs.z) < 10000.0 && fabs(newObs.z) > 1.0)
	  {
	    newObs.latitude_grs80 = fang(1, newObs.x);
	    newObs.longitude_grs80 = -fang(1, newObs.y); // W -ve
	    newObs.height_grs80 = newObs.z;
	    GRS80_to_ITRF(&newObs);
	    printf("\n\n");
	    printf("ERROR: The observatory coordinates for line %d in %s\n",
		   iline, fname);
	    printf(" appear to be in TEMPO geodetic format. Tempo2 requires cartesian geocentric\n"); 
	    printf(" coordinates w.r.t the ITRF.\n");
	    printf(" ASSUMING your coordinates are geodetic in the ITRF realisation of the GRS80\n");
	    printf(" reference ellipsoid (unlikely to be precisely true), then the geocentric\n");
	    printf(" coordinates should be specified as follows:\n");

	    printf(" %.2lf %.2lf %.2lf\n", newObs.x, newObs.y, newObs.z);
	    exit(1);
	  } 
	if (nread != 6 && nread != 5)
	{
	  fprintf(stderr, "Error parsing line %d of file %s:\n%s\n",
		  iline, fname, line);
	  exit(1);
	}
	if (nread == 5)
	  sprintf(newObs.clock_name, "UTC(%s)", newObs.code);
	ITRF_to_GRS80(&newObs);
	DynamicArray_push_back(&observatories, &newObs);
#if 0
	{
	  // test geodetic<->geocentric
	  do
	  {
	    printf("GC %.20lg %.20lg %.20lg\n", newObs.x, newObs.y, newObs.z);
	    GRS80_to_ITRF(&newObs);
	    printf("GD %.20lg %.20lg %.20lg\n",
		   newObs.latitude_grs80,
		   newObs.longitude_grs80,
		   newObs.height_grs80);
	    ITRF_to_GRS80(&newObs);
	  } while (1);
	}
#endif
      }
    }
  }

    logdbg( "readObservatoryFile: %d entries\n", observatories.nelem);

  fclose(f);  /* Added by GH */
}

void
readAliases(char *fname)
{ 
  FILE *f;
  char line[1024], alias[128];
  ObservatoryAliasList list;
  int ichar, nread;


  f = fopen(fname, "r");
  if (!f)
  {
    fprintf(stderr, "Unable to open observatory aliases file %s : %s\n",
	    fname, strerror(errno));
    exit(1);
  }
 
  DynamicArray_init(&observatoryAliasLists, sizeof(ObservatoryAliasList));

  while (fgets(line, 1024, f))
  {
    if (line[0]!='#')
    {
      if (sscanf(line, "%s %n", list.code, &ichar)==1)
      {
	DynamicArray_init(&list.aliases, 128);      
	while (sscanf(line+ichar, "%s %n", alias, &nread)==1)
	{
	  DynamicArray_push_back(&list.aliases, alias);
	  ichar += nread;
	}
	DynamicArray_push_back(&observatoryAliasLists, &list);
      }
    }
  }

  fclose(f);

}



void
initObservatories()
{
  char fname[1024];
  glob_t g;
  char pattern[1024];
  char **pfname;
  int ret;

  logdbg("In initObservatories");
  DynamicArray_init(&observatories, sizeof(observatory));
	
  /* load observatories.dat first */
  sprintf(fname, "%s/observatory/observatories.dat", getenv(TEMPO2_ENVIRON));

  logdbg("Reading observatories file >%s<",fname);
  readObservatoryFile(fname);
  logdbg("Reading other .dat files");
  /* now load other .dat files in that directory */
  sprintf(pattern, "%s/observatory/*.dat", getenv(TEMPO2_ENVIRON));
  logdbg("Reading individual files >%s<",pattern);

  ret = glob(pattern, 0, NULL, &g);
  logdbg("Return from glob = %d",ret);
  logdbg("Glob success");
  for (pfname = g.gl_pathv; *pfname != NULL; pfname++)
  {
    if (strcmp(*pfname, fname))
    {
      logdbg("Reading >%s<",*pfname);
      readObservatoryFile(*pfname);
    }
  }
  globfree(&g);
  logdbg("Reading aliases file");
  /* and load aliases list */
  sprintf(fname, "%s/observatory/aliases", getenv(TEMPO2_ENVIRON));
  readAliases(fname);

  observatories_initialised = 1;
  logdbg("Leaving initObservatories");
		   
}

void lookup_observatory_alias(char *incode, char *outcode)
{
  int ilist, ialias;
  char *alias;
  ObservatoryAliasList *list;

  //  logdbg("In lookup_observatory_alias with >%s< >%d<",incode,observatories_initialised);
  if (observatories_initialised != 1)
    initObservatories();
  //  logdbg("Initialised observatories ",incode);
  for (ilist=0; ilist < observatoryAliasLists.nelem; ilist++)
  {
    list = ((ObservatoryAliasList*)observatoryAliasLists.data)+ilist;
    for (ialias=0; ialias < list->aliases.nelem; ialias++)
    {
      alias = (char *)list->aliases.data + ialias*128;
      if (!strcasecmp(alias, incode))
	break;
    }
    if (ialias < list->aliases.nelem)
      break;
  }
  if (ilist < observatoryAliasLists.nelem && ialias < list->aliases.nelem)
  {
    //    if (debugFlag)
    //      fprintf (stderr, "Copying alias = '%s' into outcode\n", list->code);
    strcpy(outcode, list->code);
  }
  else if (outcode!=incode) {
    //      logdbg("Copying incode = '%s' over outcode\n", incode);
    strcpy(outcode, incode);
  }

  //  logdbg("Leaving lookup_observatory_alias");
}

observatory *
getObservatory(char *code)
{
  int io;
  observatory *obs;
  char alias[1024];

  //  logdbg("In getObservatory code='%s'", code);
  /* Replace site alias if necessary */
  lookup_observatory_alias(code, alias);
  //  logdbg("Checked observatory alias='%s'", alias);
  if (observatories_initialised != 1)
    initObservatories();
  //  logdbg("Initialised observatory");

  for (io=0; io < (int)observatories.nelem; io++)
  {
    obs = ((observatory*)observatories.data)+io;

    //    if (debugFlag)
    //      fprintf (stderr, "getObservatory: compare alias"
    //	       " with code='%s' or name='%s')\n", obs->code, obs->name);

    if (!strcasecmp(obs->code, alias) || !strcasecmp(obs->name, alias))
      {
	//	logdbg("leaving getObservatory");
  
	return obs;
      }
  }
  fprintf(stderr, "Observatory code '%s' not found!\n", code);
  exit(1);
}


// redwards commented out; see GRS80_to_itrf below
#if 0 
/* From setup.f in original tempo:
 *  old approach is IAU 1964 power series.  Tests show this is
 *  good to 5.10^-10 for 1964 figure.  But series misapplied
 *  below (erad series should be in alat, not hlt)
 *  so actual error is 1.d-5.  Furthermore, want flexibility
 *  of going to any figure to use full equation approach.
 *  see AA K11.
 *  IAU 1976 flattening f, equatorial radius a
 */
void convertGeocentric(observatory newObs)
{
  double alat,elev,along;
  double alng,aa_f,aa_a,aa_c,aa_arcf,aa_arsf,hlt,erad;
  double erad2,elev2,along2,alat2,hrd,ault;

  ault = 499.004786;

  alat  = newObs.x;
  along = newObs.y;
  elev  = newObs.z;
 
  printf("An approximate conversion to geocentric coordinates is: \n\n");
  alat = fang(1,alat);
  alng = fang(1,along);
  
  aa_f = 1.0/298.257;
  aa_a = 6378140.0;
  aa_c = 1.0/sqrt(1.0+(-2.0+aa_f)*aa_f*sin(alat)*sin(alat));
  aa_arcf = (aa_a*aa_c+elev)*cos(alat);
  aa_arsf = (aa_a*(1.0-aa_f)*(1.0-aa_f)*aa_c+elev)*sin(alat);
  hlt = atan2(aa_arsf,aa_arcf);
  erad = sqrt(aa_arcf*aa_arcf+aa_arsf*aa_arsf);

  hrd=erad/(2.99792458e8*ault);

  /* Convert to geocentric coords */
  
  erad2  = (2.99792458e8*ault)*hrd;
  elev2  = erad2*sin(hlt);
  along2 = sqrt(((pow(erad2,2)-pow(elev2,2))*pow(tan(alng),2))/(1.0+pow(tan(alng),2)));
  alat2  = sqrt(pow(erad2,2)-pow(elev2,2)-pow(along2,2));
  
  printf(" %11.3f   %11.3f      %11.3f      %-17.17s   %s\n",along2,alat2,elev2,newObs.name, newObs.code);
}
#endif
/* From ang.f: Converts number to radians or to fraction of 2pi */
double fang(int i,double f)
{
  int ia,ib;
  double g,ang;

  ia = (int)(f/1.0e4);
  ib = (int)((f-ia*1.0e4)/1.0e2);
  g = f-ia*1.0e4-ib*1.0e2;
  ang = (ia+(ib+g/6.0e1)/6.0e1)/36.0e1;
  if (i>2) ang=ang*15.0;
  if (i==1 || i==3) ang=ang*M_PI*2.0;

  return ang;
}

// redwards functions for handling geodetic coordinates, which are
// needed for atomspheric corrections, and also may be inadvertantly
// provided on input
#define GRS80_A 6378137.0           /* semi-major axis (m) */
#define GRS80_F 1.0/298.257222101   /* flattening */


// Geocentric to geodetic.
// Uses Vermeille (2004)'s method:
//http://www.springerlink.com/app/home/contribution.asp?wasp=08ea5d2c4c62464789a7961196d84ab5&referrer=parent&backto=issue,11,18;journal,9,85;linkingpublicationresults,1:100435,1
void
ITRF_to_GRS80(observatory *obs)
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

// Geodetic to geocentric: standard formula
void
GRS80_to_ITRF(observatory *obs)
{
  double esq = GRS80_F * (2.0 - GRS80_F);
  double N = GRS80_A / sqrt(1.0-esq*sin(obs->latitude_grs80)*sin(obs->latitude_grs80));
  obs->x = (N+obs->height_grs80)*cos(obs->latitude_grs80)*cos(obs->longitude_grs80);
  obs->y = (N+obs->height_grs80)*cos(obs->latitude_grs80)*sin(obs->longitude_grs80);
  obs->z = (N*(1.0-esq)+obs->height_grs80)*sin(obs->latitude_grs80);
}
 
