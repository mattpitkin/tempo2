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

/* Solar wind model developed by You Xiaopeng, Bill Coles and George Hobbs */
/* This model is written up in You et al. (2007)                           */

// Software for dealing with the solar wind for pulsar observations
// G. Hobbs, based on code by W. Coles (and following input from Xiaopeng You)
//
// based on mcl2.f by W. Coles
//

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "tempo2.h"

#define MAX_CURRENT 500

void mcl2(double eclon,double eclat,int iyr,int iday,double secs,double vel,
	  double helat[36],double crlon[36],double rots[36],double *helate,
	  double *crlne,double *rote,double *elong,double *beta,double *dlon,double *delng,int *nLineOfSight);
double amod(double a,double p);
double elsun2(int iyr,int iday,double secs,double *gst,double *sra,double *sdec);
void   outputResults(double *crlon,double *helat);
void   calcRotN(double crlne,double rote,int *irot1,int *irot2,double *bcrlon);
int    readCurrentSheet(char *fname,float *currentLon,float *currentLat);
double findAngle(double lon1,double lat1,double lon2,double lat2);
void mjd2date(int mjd,int *iyr,int *yy, int *mm, int *dd, int *iday); /*MJD to date, to iyr=year-1900, to iday (number of the day in a year) */
void convertEcliptic(double raj,double decj,double *elong,double *elat);

using namespace std;

 
double solarWindModel(pulsar psr,int iobs)
{
  const char *CVS_verNum = "$Revision: 1.11 $";
  int i,j;
  double deg2rad = M_PI/180.0;
  long double mjd;

  // Inputs
  
  double eclon;  // Ecliptic longitude in rad
  double eclat;  // Ecliptic latitude in rad
  int iyr,iday,yy,mm,dd;
  double secs;
  double vel;

  double angleFastSlow = 20; // Use 20 degrees as the separator between the fast and slow winds

  // Outputs
  double helat[36];
  double crlon[36];
  double rots[36];
  double helate;
  double crlne;
  double rote;
  double elong;
  double beta;
  double dlon;
  double delng;

  int irot1,irot2;
  double bcrlon;
  int nLineOfSight;

  // SWO files
  float currentLon[MAX_CURRENT];
  float currentLat[MAX_CURRENT];
  int   nCurrent;
  int fileNum;
  double closestAngle,angle;
  int iAngle;


  double fast_ne = 2.5; // From 3. This change was suggested by Bill Coles by phone to G. Hobbs
  double slow_ne = 10;
  double integral=0;
  double DM_sun;
  char fname[1024];

  int nFast,nSlow;

  if (displayCVSversion == 1) CVSdisplayVersion("sw_delay.C","solarWindModel()",CVS_verNum);

  if(psr.eclCoord==1)
    {
      eclat=psr.param[param_decj].val[0];
      eclon=psr.param[param_raj].val[0];
    }
  else if(psr.eclCoord==0)
    convertEcliptic(psr.param[param_raj].val[0],psr.param[param_decj].val[0],&eclon,&eclat);

  //  printf("%d, %30.20Lf\n",iobs, psr.obsn[iobs].sat);
  mjd=psr.obsn[iobs].sat;

  //  convertMJD(mjd,&iyr,&iday,&secs);
  mjd2date((int)mjd,&iyr,&yy, &mm, &dd, &iday); /*MJD to date, to iyr=year-1900, to iday (number of the day in a year) */
  secs=(mjd-(int)mjd)*86400;
  vel = 400;

  mcl2(eclon,eclat,iyr,iday,secs,vel,
       helat,crlon,rots,&helate,&crlne,&rote,&elong,&beta,&dlon,&delng,&nLineOfSight);
  //  for (i=0;i<36;i++)
  //    printf("%d %g %g %g %g %g %g %g %g %g %g\n",i,crlon[i],rots[i],helat[i],helate,crlne,rote,elong,beta,dlon,delng);
  
  // Now into Xiaopeng's code ...
  // calculate elsun
  // calculate rotN
  calcRotN(crlne,rote,&irot1,&irot2,&bcrlon);
  printf("rotN output = %d %d %g\n",irot1,irot2,bcrlon);
  //  outputResults(crlon,helat);

  // Go through the line-of-sight track
  integral = 0;
  nFast = 0;
  nSlow = 0;
  for (i=0;i<nLineOfSight;i++)
    {
      //      printf("Selecting: %g %g %g\n",crlon[i],helat[i],rots[i]);
      fileNum = (int)rots[i];
      sprintf(fname,"%s/solarWindModel/CS%d.txt",getenv(TEMPO2_ENVIRON),fileNum);
      nCurrent = readCurrentSheet(fname,currentLon,currentLat);      
      // Now check the closest angle to the current sheet for the current position
      closestAngle = findAngle(crlon[i],helat[i],currentLon[0]*deg2rad,currentLat[0]*deg2rad)/deg2rad;
      iAngle=0;
      for (j=1;j<nCurrent;j++)
	{
	  angle = findAngle(crlon[i],helat[i],currentLon[j]*deg2rad,currentLat[j]*deg2rad)/deg2rad;
	  if (angle < closestAngle)
	    {
	      closestAngle = angle;
	      iAngle = j;
	    } 
	}
      //      printf("Closest = %g %d\n",closestAngle,iAngle);
      if (closestAngle < angleFastSlow)
	{
	  integral += slow_ne;
	  nSlow++;
	}
      else
	{
	  integral += fast_ne;
	  nFast++;
	}
    }
  DM_sun = integral*5*(M_PI/180)*4.85e-6/sin(elong);
  printf("Integral = %.5Lg %g %g %d %d %g %g\n",mjd,integral,DM_sun,nSlow,nFast,elong,elong/deg2rad);

  return(DM_sun);
}



// Angle in radians
double findAngle(double lon1,double lat1,double lon2,double lat2)
{
  double dlon,dlat,a,c;
  double deg2rad = M_PI/180.0;

  /* Apply the Haversine formula */
  dlon = (lon2 - lon1);
  dlat = (lat2  - lat1);
  a = pow(sin(dlat/2.0),2) + cos(lat1) *
    cos(lat2)*pow(sin(dlon/2.0),2);
  if (a==1)
    c = M_PI;
  else
    c = 2.0 * atan2(sqrt(a),sqrt(1.0-a));
  return c;
}

int readCurrentSheet(char *fname,float *currentLon,float *currentLat)
{
  FILE *fin;
  int n=0;
  if (!(fin = fopen(fname,"r")))
    {
      printf("Unable to open file %s\n",fname);
      exit(1);
    }
  while (!feof(fin))
    {
      if (fscanf(fin,"%f %f",&currentLon[n],&currentLat[n])==2)
	n++;
    }
  fclose(fin);
  return n;
}

void outputResults(double *crlon,double *helat)
{
  float psrx[36],psry[36];
  int i;

  for (i=0;i<36;i++)
    {
      psrx[i] = crlon[i]*180.0/M_PI;
      if (psrx[i] <0) psrx[i] += 360.0;
      if (psrx[i] >360.0) psrx[i]-=360.0;

      psry[i] = helat[i]*180.0/M_PI;
      printf("output: %g %g\n",psrx[i],psry[i]);
    }
}

//	eclon= ecliptic longitude in rad.                       (input)
//	eclat= ecliptic latitude in rad.                        (input)
//	iyr= last two digits (year-1900) integer*4		(input)
//	iday= the day's number at that year   UT integer*4      (input)
//	secs= the seconds of that day UT      real*4            (input)
//	vel= solar wind velocity in km/sec                      (input)
//
//	helat= heliographic latitude of stream footpoint in rad.(output)
//		real*4(36) array
//	crlon= carrington longitude of stream footpoint in rad. (output)
//		real*4(36) array
//	rots= rotation number of stream footpoint               (output)
//		real*4(36) array
//	helate= heliographic latitude of sub-earth point in rad.(output)
//	crlne= carrington longitude of sub-earth point in rad.  (output)
//	rote= rotation number of sub-earth point *100           (output)
//	elong= elongation of source in rad.                     (output)
//	beta= heliographic latitude of the source
//                       as viewed from earth in rad.		(output)
//	dlon= angle at the sun between earth and point of
//                       closest approach projected onto the
//                       heliographic equator in rad.		(output)
//	delng= the elongation projected on the
//                          equatorial plane in rad.             (output)
//
//     !! BEWARE all angles are in RADIANS on output  !!  

void mcl2(double eclon,double eclat,int iyr,int iday,double secs,double vel,
	  double helat[36],double crlon[36],double rots[36],double *helate,
	  double *crlne,double *rote,double *elong,double *beta,double *dlon,double *delng,int *nLineOfSight)
{
  int i,j;
  int last=0;

  double rad2deg = 180.0/M_PI;
  double zlonan;
  double zkl;
  double ya;
  double gs,ra,de;
  double delta_PA,cl;
  double dp;

  int lp;
  double cl0;
  double inc;
  double theta;
  int finished=0;
  double clan;
  double slan;
  double sl,snx,sny;
  double b1b2,cep,blp,ble;
  double zlagtm;
  int ileap;
  double zhr,zle1,zle2,zlpp;
  double delcrle,rote1,rote2,delcrl;

  double x[3][3];
  double a[3],b[3],c[3],e[2],h[3];
  double pe1[36],pe2[36],pe3[36];

  zlonan = 1.28573 + 7.89327 * (iyr + iday/365.+ 50. )/ 32400.;
  zkl    = cos(eclat);
  // ya = ecliptic longitude of the sun
  ya     = elsun2(iyr,iday,secs,&gs,&ra,&de); // Check inputs and outputs here
  delta_PA = atan(-0.12722*cos(ya-zlonan)); // Check atan or atan2
  // cl = distance from Earth to the point P of closest approach
  cl = cos(ya-eclon)*zkl;
  dp = sqrt(1.0-cl*cl);
  // elong is the true elongation in rad
  *elong = atan2(dp,cl);
  
  // Determine latitude and longitude values for
  // points along the line of sight from the -88 degree to the earth
  // (at -(90-elong) deg) in 5 deg steps

  lp = 0;
  cl0 = cl;
  inc = 5.0/rad2deg;
  theta = M_PI/2.0 - inc;
  cl = cl0 + dp*tan(theta);
  do {
    clan = cos(zlonan);
    slan = sin(zlonan);
    
    //    x is a rotation matrix ( converts ecliptic at
    //   heliographic rest coordinates).
    //		[  cos(lonan)            sin(lonan)             0       ]
    //          [ -cos(b0)sin(lonan)     cos(b0)cos(lonan)    sin(b0)   ]
    //		[  sin(b0)sin(lonan)    -sin(b0)cos(lonan)    cos(b0)   ]


    x[0][0]=clan;
    x[0][1]=slan;
    x[0][2] = 0.0;
    x[1][0]=-0.9920049497*slan;
    x[1][1]=0.9920049497*clan;
    x[1][2] = 0.1262; // sin(b0)
    x[2][0]=0.1261989691*slan;
    x[2][1]=-0.1261989691*clan;
    x[2][2] = 0.9920; // cos(b0) b0 is the tilt of the sun's spin axis

    sl = sqrt((1.0-cl0*cl0)+(pow(cl0-cl,2)));
    //    sl= distance from the sun to the point P
    //	special calculation if the elongation is more than 73 deg
    snx = -cos(ya);
    sny = -sin(ya);
    // snx and sny are the ecliptic coordinates of the earth as
    // viewed from the Sun
    //		a is a unit vector from the earth to the source
    //		a(1) = cos(eclat)*cos(eclon)
    //          a(2) = cos(eclat)*sin(eclon)
    //		a(3) = sin(eclat)
    //	a are cartesian ecliptic coordinates
    a[0] = cos(eclon)*zkl;
    a[1] = sin(eclon)*zkl;
    a[2] = sin(eclat);


    // be are the ecliptic coordinates of the points of
    // closest approach, with the sun as the origin
    b[0] = a[0]*cl+snx;
    b[1] = a[1]*cl+sny;
    b[2] = a[2]*cl;

    if (lp > 35) {
      printf("ERROR: lp cannot be %d\n",lp);
      exit(1);
    }
    pe1[lp] = b[0];
    pe2[lp] = b[1];
    pe3[lp] = b[2];

    // c is a unit vector pointing from the earth to the
    // source in heliographic coordinates
    for (i=0;i<3;i++)
      {
	c[i]=0.0;
	for (j=0;j<3;j++)
	  c[i]+=x[i][j]*a[j];
      }

    // beta is the latitude of the source above the equatorial as
    // viewed from earth in heliographic coordinates
    *beta = atan2(c[2],sqrt(c[0]*c[0]+c[1]*c[1]));

    // a are now the heliographic coordinates of the point of closest
    // approach
    for (i=0;i<3;i++)
      {
	a[i] = 0.0;
	for (j=0;j<3;j++)
	  a[i] = a[i]+x[i][j]*b[j];
      }

    // helat = heliographic latitude of point P
    helat[lp] = atan2(a[2],sqrt(a[0]*a[0]+a[1]*a[1]));

    // b are now the heliographic coordinates of the earth
    b[0] = x[0][0]*snx+x[0][1]*sny;
    b[1] = x[1][0]*snx+x[1][1]*sny;
    b[2] = x[2][0]*snx+x[2][1]*sny;

    b1b2 = sqrt(b[0]*b[0]+b[1]*b[1]);
    e[0] = snx;
    e[1] = sny;
    
    h[0] = b[0];
    h[1] = b[1];
    h[2] = b[2];

    // heliographic latitude of the earth at the time of observation
    *helate = atan2(b[2],b1b2);
    cep = ((-b[0]*c[0]-b[1]*c[1])/b1b2)/sqrt(c[0]*c[0]+c[1]*c[1]);
    if (cep > 1) cep=1;

    // the equatorial projection of the elongation
    *delng = atan2(sqrt(1.-cep*cep),cep);
    blp = atan2(a[1],a[0]);
    if (blp < 0) blp+=2.0*M_PI;
    ble = atan2(b[1],b[0]);

    // ble is the longitude of the earth
    if (ble < 0) ble+=2.0*M_PI;

    // propagation time frmo the sun to the scattering point
    zlagtm = 1.496e8*sl/3600.0/vel;
    ileap = (iyr-69)/4; // CHECK THIS
    // number of hours that elapsed since the reference date: dec. 7th 1969, ut 14.77 hr.
    zhr = 24.*((iyr-69)*365.+ileap+iday-341.0) + secs/3600. -14.77;

    zle1 = 228.42-zhr/1.81835;
    //		this is an approximate algorithm
    //		for zle, independent of a calculation of the earth's
    //		position.  It assumes constant angular speed of earth in
    //		orbit; (1.818735 is the hours per degree of the
    //		sub earth point, ie omegasun - omegaearthorbit)
    //	The accurate code is below - difference of up to
    //	2 deg due to eccentricity of earth's orbit

    zle2 = 228.42 - zhr/1.692002 + ble*rad2deg;
    // 		carrington long of earth was 228.42 on reference date/time
    *crlne=360.+amod(zle2,360.); // WHAT IS AMOD???
    if (*crlne >= 360.0) *crlne = amod(*crlne,360.0);
    
    // crlne - carrington longitude of sub-earth point at the time of the
    // observation in radians
    delcrle = zlagtm/1.692002;
    rote1 = 1556.-zle1/360.0;
    rote2 = 1556.-zle2/360.0;
    *rote = (int)(rote1) + (rote2 - (int)rote2);

    if (rote1 - *rote > (4.0/360.0)) *rote = *rote + 1;
    if (rote1 - *rote < (-4.0/360.0)) *rote = *rote - 1;

    // heliographic longitude difference between the earth and P in radians
    *dlon = blp-ble;
    //		1.692002 is the hours per degree in fixed heliographic frame.
    //		crlon is longitude of stream footpoint in the rotating system
    delcrl = zlagtm/1.692002 + (*dlon)*rad2deg;
    zlpp = zle2+delcrl;

    crlon[lp] = 360.0 + amod(zlpp,360.0);
    if (crlon[lp] >= 360.0) crlon[lp] = amod(crlon[lp],360.0);

    if (crlon[lp] >= *crlne)
      rots[lp] = (int)(*rote) + (1.0-crlon[lp]/360.0);
    else
      {
	rots[lp] = (int)(*rote) - (crlon[lp]/360.0);
	if (crlon[lp]/360.0 < 0.01) rots[lp] = (int)(*rote)-0.01;
      }
    crlon[lp] = crlon[lp]/rad2deg;
    *crlne = *crlne/rad2deg;

    theta = theta - inc;
    cl = cl0 + dp*tan(theta);
    
    if (last) finished=1;
    if (cl < 0.0)
      {
	cl = 0.0;
	last = 1;
      }
    lp++; // Moved from the Fortran to ensure starting at 0
  } while (finished==0);
  *nLineOfSight = lp;

  for (i=lp;i<36;i++)
    {
      helat[i] = helat[lp-1];
      crlon[i] = crlon[lp-1];
      rots[i]  = rots[lp-1];
    }
}

double amod(double a,double p)
{
  double out;
  if (p==0) return 0;
  out = a-(int)(a/p)*p;
  return out;
}

// Useful subroutines

double elsun2(int iyr,int iday,double secs,double *gst,double *sra,double *sdec)
{
// ----------------------------------------------------------------------
// 	elsun : function to calculate sidereal time and position of the
// 	sun good for year 1901 through 2099 (only last two digits of the
// 	year are used, 1900 is suppressed). accuracy 0.006 degree.
// 	Input is iyr, iday(integers), and secs, defining universal time.
// 	Output is greenwich mean side real time (gst) in degrees;
// 	longitude along ecliptic (elsun2), and apparent right ascension
// 	and declination (sra,sdec) of the sun, all in radians.
// -----------------------------------------------------------------

  double dj,fday;
  double rad = 57.29578;
  int idd;
  double t,vl,g;
  double elsun;
  double obliq,slp,sind,cosd,cot,elsun2;
  //
  if (iyr < 1  || iyr > 199) {
    printf("ERROR in elsun2\n");
    exit(1);
  }

  fday= (float)(secs)/86400.0;
  idd= 365*iyr + (iyr-1)/4 + iday;
  dj= idd + fday - 0.5;
  t= dj/36525.0;
  vl= amod(279.696678+0.9856473354*dj,360.0);
  *gst= amod(279.690983+0.9856473354*dj+360.*fday+180.,360.);
  g= amod(358.475845+0.985600267*dj,360.0)/rad;
  elsun= vl+(1.91946-0.004789*t)*sin(g)+0.020094*sin(2.*g);
  obliq= (23.45229- 0.0130125*t)/rad;
  slp= (elsun-0.005686)/rad;
  sind= sin(obliq)*sin(slp);
  cosd= sqrt(1.-sind*sind);
  *sdec= rad*atan(sind/cosd);
  cot= cos(obliq)/sin(obliq);
  *sra= 180.-rad*atan2(sind/cosd*cot,-cos(slp)/cosd);
  elsun2= elsun/rad;
  *sdec= *sdec/rad;
  *sra= *sra/rad;
  return elsun2;
}

// From xiaopeng's code - should check
void calcRotN(double crlne,double rote,int *irot1,int *irot2,double *bcrlon)
{
  *bcrlon = (double)((int)((0.5+crlne*180.0/M_PI))+180.0);

  if (crlne <= M_PI)
    {
      *irot1 = (int)rote;
      *irot2 = (int)rote+1;
    }
  else
    {
      *irot1 = (int)(rote)-1;      
      *irot2 = (int)rote;
      if (*bcrlon > 360.0)
	*bcrlon = *bcrlon - 360.0;
    }
}

void mjd2date(int mjd,int *iyr,int *yy, int *mm, int *dd, int *iday) /*MJD to date, to iyr=year-1900, to iday (number of the day in a year) */
{
  int y1,m1,k;
  y1=(int)(((double)mjd-15078.2)/365.25);
  m1=(int)(((double)mjd-14956.1-(int)((double)y1*365.25))/30.6001);
  *dd=mjd-14956-(int)((double)y1*365.25)-(int)((double)m1*30.6001);
  k=0;
  if( m1==14 || m1==15 ) k=1;
  *mm=m1-1-k*12;
  int mk=m1-1-k*12;
  *iyr=y1+k;
  *yy=*iyr+1900;
  if((*yy%4==0 && *yy%100!=0) || (*yy%400==0))
    {
      if(*mm==1)*iday=    (*dd);
      else if(*mm==2)*iday=31 +(*dd);
      else if(*mm==3)*iday=60 +(*dd);
      else if(*mm==4)*iday=91 +(*dd);
      else if(*mm==5)*iday=121+(*dd);
      else if(*mm==6)*iday=152+(*dd);
      else if(*mm==7)*iday=182+(*dd);
      else if(*mm==8)*iday=213+(*dd);
      else if(*mm==9)*iday=244+(*dd);
      else if(*mm==10)*iday=274+(*dd);
      else if(*mm==11)*iday=305+(*dd);
      else if(*mm==12)*iday=335+(*dd);
      else printf("date error\n");
    }
  else
    {
      if(*mm==1)*iday=    (*dd);
      else if(*mm==2)*iday=31 +(*dd);
      else if(*mm==3)*iday=59 +(*dd);
      else if(*mm==4)*iday=90 +(*dd);
      else if(*mm==5)*iday=120+(*dd);
      else if(*mm==6)*iday=151+(*dd);
      else if(*mm==7)*iday=181+(*dd);
      else if(*mm==8)*iday=212+(*dd);
      else if(*mm==9)*iday=243+(*dd);
      else if(*mm==10)*iday=273+(*dd);
      else if(*mm==11)*iday=304+(*dd);
      else if(*mm==12)*iday=334+(*dd);
      else printf("date error\n");
    }
}


void convertEcliptic(double raj,double decj,double *elong,double *elat)
{
  double sinb,beta,x,y,lambdap,lambda;
  double deg2rad = M_PI/180.0;
  double epsilon = 23.439292*deg2rad;
  /*  double epsilon = 23.441884*deg2rad;*/

  sinb = sin(decj)*cos(epsilon)-cos(decj)*sin(epsilon)*sin(raj);
  beta = asin(sinb);
  y = sin(raj)*cos(epsilon)+tan(decj)*sin(epsilon);
  x = cos(raj);
  
  lambdap = atan2(y,x);
  if (lambdap<0) lambda=lambdap+2*M_PI;
  else lambda = lambdap;

  *elong = lambda;
  *elat  = beta;
}
