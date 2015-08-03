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
#include <errno.h>
#include "tempo2.h"

double get_precessionMatrix(double prn[3][3],double mjd,double delp,double dele);
void remove_white(char *str);
double lmst(double mjd,double olong,double *tsid,double *tsid_der);
double ang(int i,double f);

extern "C" long iau_cp_(double *, double *);
extern "C" long iau_pom00_(double *, double *, double *, double *);
extern "C" long iau_rxp_(double *, double *, double *);
extern "C" long iau_sxp_(double *, double *, double *);
extern "C" long iau_pxp_(double *, double *, double *);
extern "C" long iau_c2t00b_(double *, double *, double *, double *, 
			    double *, double *, double *);
extern "C" long iau_trxpv_(double *, double *, double *);
extern "C" long iau_trxp_(double *, double *, double *);

void 
get_obsCoord_IAU2000B(double observatory_trs[3],
		      double zenith_trs[3],
		      longdouble tt_mjd, longdouble utc_mjd,
		      double observatory_crs[3],
		      double zenith_crs[3],
		      double observatory_velocity_crs[3],char *eopcFile)
{
  double trs[2][3], crs[2][3], north[3], pole_itrs[3], omega_itrs[3];
  double t2c[3][3], polarmotion[3][3];
  double dut1, dut1dot, eradot, xp, yp; 
  longdouble tt_jd = tt_mjd + (longdouble)2400000.5, ut1_jd;
  double tt_jd1 = (int)tt_jd, tt_jd2 = tt_jd-tt_jd1;
  double ut1_jd1, ut1_jd2;
  double sprime=0.0;
  double one_on_c = 1.0/SPEED_LIGHT;
  const char *CVS_verNum = "$Revision: 1.12 $";

  if (displayCVSversion == 1) CVSdisplayVersion("get_obsCoord.C","get_obsCoord_IAU2000B()",CVS_verNum);

  // Get Earth orientation parameters
  get_EOP((double)utc_mjd, &xp, &yp, &dut1, &dut1dot, 2,eopcFile);
  ut1_jd = utc_mjd + dut1/86400.0 + (longdouble)2400000.5;
  ut1_jd1 = (int)ut1_jd;
  ut1_jd2 = ut1_jd-ut1_jd1;

  // Stick the site position vector in
  iau_cp_(observatory_trs, trs[0]);

  // Work out site velocity in CRS... 
  // first we need to know the angular velocity vector
  iau_pom00_(&xp, &yp, &sprime, polarmotion[0]); // polar motion matrix
  north[0] = north[1] = 0.0; north[2] = 1.0; // Vector to +ve pole
  iau_rxp_(polarmotion[0], north, pole_itrs); // Spin pole in ITRS 
  eradot = 2.0*M_PI*1.00273781191135448*(1.0+dut1dot)/86400.0;
  iau_sxp_(&eradot, pole_itrs, omega_itrs); // Angular velocity in ITRS (rad/s)
  iau_pxp_(omega_itrs, trs[0], trs[1]); // Tangential velocity (m/s)

  // Get the Celestial->Terrestrial matrix
  iau_c2t00b_(&tt_jd1, &tt_jd2, &ut1_jd1, &ut1_jd2, &xp, &yp, t2c[0]);

  // Multiply the itrs position/velocity vector by its transpose (=inverse)
  // to transform trs->crs
  iau_trxpv_(t2c[0], trs[0], crs[0]);
  iau_trxp_(t2c[0], zenith_trs, zenith_crs);

  // Convert to light seconds / dimensionless units and copy to output
  iau_sxp_(&one_on_c, crs[0], observatory_crs);
  iau_sxp_(&one_on_c, crs[1], observatory_velocity_crs);
}



/* NOT COMPLETE REQUIRES NUTATIONS TO BE KNOWN: Must read the ephemeris first */

/* ******************************************** */
/* get_obsCoord                                 */
/* Author:  G. Hobbs (06 May 2003)              */
/* Purpose: Calculates the position of the      */
/*          observatory site relative to the    */
/*          Earth's centre                      */
/*                                              */
/* Inputs:                                      */
/* Outputs:                                     */
/*                                              */
/* Notes: based on obsite.f                     */
/*        NONE of this has been checked         */
/*                                              */
/*   Requires nutations to be set               */
/*   Results are not quite the same as TEMPO    */
/*   THIS MUST BE CHECKED                       */
/* Changes:                                     */
/* ******************************************** */

void get_obsCoord(pulsar *psr,int npsr)
{
  double ph,eeq[3];
  int i,j,k; 
  double prn[3][3];
  double pc,oblq,toblq;
  double siteCoord[3];
  double erad,hlt,alng,hrd,tsid,sdd,speed,sitera;
  int p;
  observatory *obs;
  const char *CVS_verNum = "$Revision: 1.12 $";

  if (displayCVSversion == 1) CVSdisplayVersion("get_obsCoord.C","get_obsCoord()",CVS_verNum);

  for (p=0;p<npsr;p++)
    {
       for (i=0;i<psr[p].nobs;i++)
	{ 
	  if (psr[p].obsn[i].delayCorr!=0)
	    {	     
	      if (strcmp(psr[p].obsn[i].telID,"STL")==0 ||
		  strcmp(psr[p].obsn[i].telID,"STL_FBAT")==0) // Satellite
		{
		  psr[p].obsn[i].siteVel[0] = 0.0;
		  psr[p].obsn[i].siteVel[1] = 0.0;
		  psr[p].obsn[i].siteVel[2] = 0.0;		  

		  // Set from the TELX, TELY, TELZ parameters
		  if (psr[p].param[param_telx].paramSet[0] == 1)
		    {
		      longdouble arg,deltaT,t0,pos;

		      if (psr[p].param[param_telEpoch].paramSet[0]==1)
			t0 = psr[p].param[param_telEpoch].val[0];
		      else
			t0 = 0.0L;

		      deltaT = (psr[p].obsn[i].sat - t0)*SECDAY;		
		      arg = deltaT;
		      pos = psr[p].param[param_telx].val[0];
		      if (psr[p].param[param_telx].paramSet[1] == 1) {
			pos += psr[p].param[param_telx].val[1]*arg;
			if (psr[p].param[param_telEpoch].paramSet[0]==0)
			  {
			    printf("ERROR: Using telescope velocity without setting telEpoch\n");
			    exit(1);
			  }
		      }
		      arg *= deltaT; if (psr[p].param[param_telx].paramSet[2] == 1) pos += (0.5*psr[p].param[param_telx].val[2]*arg);
		      arg *= deltaT; if (psr[p].param[param_telx].paramSet[3] == 1) pos += (1.0L/6.0L*psr[p].param[param_telx].val[3]*arg);
		      psr[p].obsn[i].observatory_earth[0] = (double)pos;
		      //logdbg("Setting x to %g",(double)psr[p].obsn[i].observatory_earth[0]);
		    }
		  if (psr[p].param[param_tely].paramSet[0] == 1)
		    {
		      longdouble arg,deltaT,t0,pos;

		      if (psr[p].param[param_telEpoch].paramSet[0]==1)
			t0 = psr[p].param[param_telEpoch].val[0];
		      else
			t0 = 0.0L;

		      deltaT = (psr[p].obsn[i].sat - t0)*SECDAY;		
		      arg = deltaT;
		      pos = psr[p].param[param_tely].val[0];
		      if (psr[p].param[param_tely].paramSet[1] == 1) pos += psr[p].param[param_tely].val[1]*arg;
		      arg *= deltaT; if (psr[p].param[param_tely].paramSet[2] == 1) pos += (0.5*psr[p].param[param_tely].val[2]*arg);
		      arg *= deltaT; if (psr[p].param[param_tely].paramSet[3] == 1) pos += (1.0L/6.0L*psr[p].param[param_tely].val[3]*arg);
		      psr[p].obsn[i].observatory_earth[1] = (double)pos;
		      //		      printf("Setting y to %g\n",(double)psr[p].obsn[i].observatory_earth[1]);
		    }
		  if (psr[p].param[param_telz].paramSet[0] == 1)
		    {
		      longdouble arg,deltaT,t0,pos;

		      if (psr[p].param[param_telEpoch].paramSet[0]==1)
			t0 = psr[p].param[param_telEpoch].val[0];
		      else
			t0 = 0.0L;

		      deltaT = (psr[p].obsn[i].sat - t0)*SECDAY;		
		      arg = deltaT;
		      pos = psr[p].param[param_telz].val[0];
		      if (psr[p].param[param_telz].paramSet[1] == 1) pos += psr[p].param[param_telz].val[1]*arg;
		      arg *= deltaT; if (psr[p].param[param_telz].paramSet[2] == 1) pos += (0.5*psr[p].param[param_telz].val[2]*arg);
		      arg *= deltaT; if (psr[p].param[param_telz].paramSet[3] == 1) pos += (1.0L/6.0L*psr[p].param[param_telz].val[3]*arg);
		      psr[p].obsn[i].observatory_earth[2] = (double)pos;
		      //		      printf("Setting z to %g\n",(double)psr[p].obsn[i].observatory_earth[2]);
		    }



		  //		  printf("Setting observatory coordinates for TOA %d\n",i);
		  // Now check flags to obtain the telescope coordinates for this time
		  for (k=0;k<psr[p].obsn[i].nFlags;k++)
		    {
		      //
		      // NOTE: For a STL_FBAT setting OBSERVATORY->BAT not OBSERVATORY->EARTH
		      // The terminology is misleading
		      //
		      if (strcmp(psr[p].obsn[i].flagID[k],"-telx")==0){
			sscanf(psr[p].obsn[i].flagVal[k],"%lf",&psr[p].obsn[i].observatory_earth[0]);
			if (strcmp(psr[p].obsn[i].telID,"STL")==0) psr[p].obsn[i].observatory_earth[0]/=SPEED_LIGHT;
		      }
		      if (strcmp(psr[p].obsn[i].flagID[k],"-tely")==0){
			sscanf(psr[p].obsn[i].flagVal[k],"%lf",&psr[p].obsn[i].observatory_earth[1]);
			if (strcmp(psr[p].obsn[i].telID,"STL")==0) psr[p].obsn[i].observatory_earth[1]/=SPEED_LIGHT;
		      }
		      if (strcmp(psr[p].obsn[i].flagID[k],"-telz")==0){
			sscanf(psr[p].obsn[i].flagVal[k],"%lf",&psr[p].obsn[i].observatory_earth[2]);		      
			if (strcmp(psr[p].obsn[i].telID,"STL")==0) psr[p].obsn[i].observatory_earth[2]/=SPEED_LIGHT;
		      }
			
		    }
		  // Check for offsets in satellite position
		  if (psr[p].nTelDX > 0)
		    {
		      int k;
		      double m,c,ival;
		      ival=0.0;
		      //	      printf("Using interpolation function for telDX\n");
		      for (k=0;k<psr[p].nTelDX-1;k++)
			{
			  if ((double)psr[p].obsn[i].sat >= psr[p].telDX_t[k] &&
			      (double)psr[p].obsn[i].sat < psr[p].telDX_t[k+1])
			    {
			      m = (psr[p].telDX_v[k]-psr[p].telDX_v[k+1])/(psr[p].telDX_t[k]-psr[p].telDX_t[k+1]);
			      c = psr[p].telDX_v[k]-m*psr[p].telDX_t[k];
			      
			      if (psr[p].param[param_tel_dx].val[0] == 0 ||
				  psr[p].param[param_tel_dx].val[0] == 1)
				ival = m*(double)psr[p].obsn[i].sat+c;
			      else if (psr[p].param[param_tel_dx].val[0] == 2)
				ival = psr[p].telDX_v[k];
			      break;
			    }
			}
		      //		      printf("xpos add (stl) = %g %d\n",ival,i);
		      psr[p].obsn[i].observatory_earth[0] += ival; // Change x-position
		    }		  
		  // Check for offsets in satellite position
		  if (psr[p].nTelDY > 0)
		    {
		      int k;
		      double m,c,ival;
		      ival=0.0;
		      //	      printf("Using interpolation function for telDX\n");
		      for (k=0;k<psr[p].nTelDY-1;k++)
			{
			  if ((double)psr[p].obsn[i].sat >= psr[p].telDY_t[k] &&
			      (double)psr[p].obsn[i].sat < psr[p].telDY_t[k+1])
			    {
			      m = (psr[p].telDY_v[k]-psr[p].telDY_v[k+1])/(psr[p].telDY_t[k]-psr[p].telDY_t[k+1]);
			      c = psr[p].telDY_v[k]-m*psr[p].telDY_t[k];
			      
			      if (psr[p].param[param_tel_dy].val[0] == 0 ||
				  psr[p].param[param_tel_dy].val[0] == 1)
				ival = m*(double)psr[p].obsn[i].sat+c;
			      else if (psr[p].param[param_tel_dy].val[0] == 2)
				ival = psr[p].telDY_v[k];

			      break;
			    }
			}
		      //		      printf("ypos add (stl) = %g %d\n",ival,i);
		      psr[p].obsn[i].observatory_earth[1] += ival; // Change y-position
		    }		  
		  // Check for offsets in satellite position
		  if (psr[p].nTelDZ > 0)
		    {
		      int k;
		      double m,c,ival;
		      ival=0.0;
		      //	      printf("Using interpolation function for telDX\n");
		      for (k=0;k<psr[p].nTelDZ-1;k++)
			{
			  if ((double)psr[p].obsn[i].sat >= psr[p].telDZ_t[k] &&
			      (double)psr[p].obsn[i].sat < psr[p].telDZ_t[k+1])
			    {
			      m = (psr[p].telDZ_v[k]-psr[p].telDZ_v[k+1])/(psr[p].telDZ_t[k]-psr[p].telDZ_t[k+1]);
			      c = psr[p].telDZ_v[k]-m*psr[p].telDZ_t[k];
			      
			      if (psr[p].param[param_tel_dz].val[0] == 0 ||
				  psr[p].param[param_tel_dz].val[0] == 1)
				ival = m*(double)psr[p].obsn[i].sat+c;
			      else if (psr[p].param[param_tel_dz].val[0] == 2)
				ival = psr[p].telDZ_v[k];
			      break;
			    }
			}
		      //		      printf("zpos add (stl) = %g %d\n",ival,i);
		      psr[p].obsn[i].observatory_earth[2] += ival; // Change z-position
		    }		  
		}
	      else if (strcmp(psr[p].obsn[i].telID,"STL_BAT")==0) // Satellite in barycentric coordinates
		{
		  psr[p].obsn[i].siteVel[0] = 0.0;
		  psr[p].obsn[i].siteVel[1] = 0.0;
		  psr[p].obsn[i].siteVel[2] = 0.0;
		  psr[p].obsn[i].observatory_earth[0] = 0.0;
		  psr[p].obsn[i].observatory_earth[1] = 0.0;
		  psr[p].obsn[i].observatory_earth[2] = 0.0;
		  psr[p].correctTroposphere = 0;
		}
	      else
		{
		  obs = getObservatory(psr[p].obsn[i].telID);
		  // Check for override:
		  if (strcmp(psr[p].obsn[i].telID,"IMAG")==0)
		    {
		      //
		      // Must check what is happening to other parameters - such as velocities
		      //
		      for (k=0;k<psr[p].obsn[i].nFlags;k++)
			{
			  if (strcmp(psr[p].obsn[i].flagID[k],"-telx")==0){
			    sscanf(psr[p].obsn[i].flagVal[k],"%lf",&obs->x);
			  }
			  if (strcmp(psr[p].obsn[i].flagID[k],"-tely")==0){
			    sscanf(psr[p].obsn[i].flagVal[k],"%lf",&obs->y);
			  }
			  if (strcmp(psr[p].obsn[i].flagID[k],"-telz")==0){
			    sscanf(psr[p].obsn[i].flagVal[k],"%lf",&obs->z);		      
			  }
			}		
		    }
		  // New way
		  if (psr[p].t2cMethod == T2C_IAU2000B)
		    {
		      double trs[3], zenith_trs[3];
		      trs[0]=obs->x;
		      trs[1]=obs->y;
		      trs[2]=obs->z;
		      zenith_trs[0]= obs->height_grs80 
			* cos(obs->longitude_grs80) * cos(obs->latitude_grs80);
		      zenith_trs[1]= obs->height_grs80 
			* sin(obs->longitude_grs80) * cos(obs->latitude_grs80);
		      zenith_trs[2] = obs->height_grs80*sin(obs->latitude_grs80);
		      long double utc = psr[p].obsn[i].sat;
		      if (psr[p].obsn[i].clockCorr!=0 && psr[p].obsn[i].clockCorr!=2)
			utc += getCorrection(psr[p].obsn+i, 
					     psr[p].clockFromOverride,
					     "UTC", psr[p].noWarnings)/SECDAY;
		      get_obsCoord_IAU2000B(trs, zenith_trs,
					    psr[p].obsn[i].sat
					    +getCorrectionTT(psr[p].obsn+i)/SECDAY,
					    utc,
					    psr[p].obsn[i].observatory_earth,
					    psr[p].obsn[i].zenith,
					    psr[p].obsn[i].siteVel,psr[p].eopc04_file);
          if( psr[p].eclCoord == 1 ){
            equ2ecl( psr[p].obsn[i].observatory_earth );
            equ2ecl( psr[p].obsn[i].siteVel );
          }
		    }
		  else {
		    psr[p].obsn[i].zenith[0]=psr[p].obsn[i].zenith[1]=psr[p].obsn[i].zenith[2]=0.0; // Only calc'd by IAU code
		    //	      if (	psr[p].obsn[i].zenith[2]==0.0)
		    erad = sqrt(obs->x*obs->x+obs->y*obs->y+obs->z*obs->z);//height(m)
		    hlt  = asin(obs->z/erad); // latitude
		    alng = atan2(-obs->y,obs->x); // longitude
		    hrd  = erad/(2.99792458e8*499.004786); // height (AU)
		    siteCoord[0] = hrd * cos(hlt) * 499.004786; // dist from axis (lt-sec)
		    siteCoord[1] = siteCoord[0]*tan(hlt); // z (lt-sec)
		    siteCoord[2] = alng; // longitude
		    
		    /* PC,PS equatorial and meridional components of nutations of longitude */      
		    toblq = (psr[p].obsn[i].sat+2400000.5-2451545.0)/36525.0;
		    oblq = (((1.813e-3*toblq-5.9e-4)*toblq-4.6815e1)*toblq +84381.448)/3600.0;
		    
		    pc = cos(oblq*M_PI/180.0+psr[p].obsn[i].nutations[1])*psr[p].obsn[i].nutations[0];
		    
		    /* TSID = sidereal time (lmst, timcalc -> obsite) */
		    
		    lmst(psr[p].obsn[i].sat+psr[p].obsn[i].correctionUT1/SECDAY
			 ,0.0,&tsid,&sdd);
		    tsid*=2.0*M_PI;
		    
		    /* Compute the local, true sidereal time */
		    ph = tsid+pc-siteCoord[2];  
		    /* Get X-Y-Z coordinates taking out earth rotation */
		    eeq[0] = siteCoord[0]*cos(ph); 
		    eeq[1] = siteCoord[0]*sin(ph);
		    eeq[2] = siteCoord[1];
		    
		    /* Now obtain PRN -- the precession matrix */
		    get_precessionMatrix(prn,(double)psr[p].obsn[i].sat
					 +psr[p].obsn[i].correctionUT1/SECDAY
					 ,psr[p].obsn[i].nutations[0],
					 psr[p].obsn[i].nutations[1]);
		    /* Calculate the position after precession/nutation*/
		    for (j=0;j<3;j++)
		      psr[p].obsn[i].observatory_earth[j]=prn[j][0]*eeq[0]+prn[j][1]*eeq[1]+prn[j][2]*eeq[2]; 
		    /* Rotate vector if we are working in ecliptic coordinates */
		    if (psr[p].eclCoord==1) equ2ecl(psr[p].obsn[i].observatory_earth);
		    
		    /* Calculate observatory velocity w.r.t. geocentre (1950.0)!!!! <<<<<<< */
		    speed = 2.0*M_PI*siteCoord[0]/(86400.0/1.00273);
		    sitera = 0.0;
		    if (speed>1.0e-10) sitera = atan2(psr[p].obsn[i].observatory_earth[1],
						      psr[p].obsn[i].observatory_earth[0]);
		    psr[p].obsn[i].siteVel[0] = -sin(sitera)*speed;
		    psr[p].obsn[i].siteVel[1] =  cos(sitera)*speed;
		    psr[p].obsn[i].siteVel[2] =  0.0;
		    if (psr[p].eclCoord==1) equ2ecl(psr[p].obsn[i].siteVel);
		    
		    /* Technically if using TDB these coordinates should be 
		       transformed to that frame. In practise it doesn't matter,
		       we're only talking about 0.3 ns, or 2.5e-14 in v/c */
		    
		    /* hack to transform for faked values of "K" */
		    //  	      vectorscale(psr[p].obsn[i].observatory_earth, 1.15505);
		    //  	      vectorscale(psr[p].obsn[i].siteVel, 1.15505);
		  }
		}
	    }
	  else if (i==0)
	    printf("Delay correction turned off for psr %d\n", p);
	}
    }
}


/* Remove white space and new line characters from a line. Leaving one white space is multiple spaces are given */
/* All white spaces are removed at the start of a line: */
/*   fred    blogs joe   kate     -> */
/* fred blogs joe kate               */
void remove_white(char *str)
{
  int i;
  int time=1;

  for (i=0;i<(int)strlen(str);i++)
    {
      if ((str[i]==' ' || str[i]=='\t') && time==1)
	{
	  strcpy(str+i,str+i+1);
	  i--;
	}
      else if (time==2 && (str[i]==' ' || str[i]=='\t') && (str[i-1]==' ' || str[i-1]=='\t')) 
	{
	  strcpy(str+i,str+i+1);
	  i--;
	}	       
      else if (str[i]=='\n')
	{
	  strcpy(str+i,str+i+1);
	  i--;
	}
      else if (str[i]!=' ' && str[i]!='\t')
	  time=2;
    }
}

/* ******************************************** */
/* lmst                                         */
/* Author:  G. Hobbs (06 May 2003)              */
/* Purpose: Calculates the local mean sidereal  */
/*          time (tsid) and its derivative      */
/*                                              */
/* Inputs:                                      */
/* Outputs: tsid, tsid_der                      */
/*                                              */
/* Notes: based on lmst.f                       */
/*        NONE of this has been checked         */
/* Changes:                                     */
/* ******************************************** */

double lmst(double mjd,double olong,double *tsid,double *tsid_der)
{
  double xlst,sdd;
  double gmst0;
  double a = 24110.54841;
  double b = 8640184.812866;
  double c = 0.093104;
  double d = -6.2e-6;
  double bprime,cprime,dprime;
  double tu0,fmjdu1,dtu,tu,seconds_per_jc,gst;
  int nmjdu1;

  nmjdu1 = (int)mjd;
  fmjdu1 = mjd - nmjdu1;

  tu0 = ((double)(nmjdu1-51545)+0.5)/3.6525e4;
  dtu  =fmjdu1/3.6525e4;
  tu = tu0+dtu;
  gmst0 = (a + tu0*(b+tu0*(c+tu0*d)))/86400.0;
  seconds_per_jc = 86400.0*36525.0;

  bprime = 1.0 + b/seconds_per_jc;
  cprime = 2.0 * c/seconds_per_jc;
  dprime = 3.0 * d/seconds_per_jc;

  sdd = bprime+tu*(cprime+tu*dprime);

  gst = gmst0 + dtu*(seconds_per_jc + b + c*(tu+tu0) + d*(tu*tu+tu*tu0+tu0*tu0))/86400;
  xlst = gst - olong/360.0;
  xlst = fortran_mod(xlst,1.0);

  if (xlst<0.0)xlst=xlst+1.0;

  *tsid = xlst;
  *tsid_der = sdd;
  return 0.0;
}

/* Copied from the Fortran, no checks have been made */

double get_precessionMatrix(double prn[3][3],double mjd,double delp,double dele)
{
  int i,j;
  /* For precession */
  double t,zeta,dzeta,z,dz,theta,dtheta;
  double par_zeta[3] = {2306.2181, 0.30188, 0.017998};
  double par_z[3] = {2306.2181, 1.09468, 0.018203};
  double par_theta[3] = {2004.3109, -0.42665, -0.041833};
  double seconds_per_rad = 3600.0*180.0/M_PI;
  
  double czeta,szeta,dczeta,dszeta,cz,sz,dcz,dsz,ctheta,stheta,dctheta,dstheta;

  double nut[3][3],prc[3][3],dprecess[3][3];

  /* For nutation */
  double dt;
  double eps=OBLQ*M_PI/180.0;
  double ceps;
  double seps;

  /* Precession -- copied from precession.f */
  /*  printf("PRECESS %f\n",mjd); */
  t = (mjd-51544.5)/36525.0;

  zeta=t*(par_zeta[0]+t*(par_zeta[1]+t*par_zeta[2]))/seconds_per_rad;
  dzeta=(par_zeta[0]+t*(2.0*par_zeta[1]+t*3.0*par_zeta[2]))/seconds_per_rad/36525.0;
  z=t*(par_z[0]+t*(par_z[1]+t*par_z[2]))/seconds_per_rad;
  dz=(par_z[0]+t*(2.0*par_z[1]+t*3.0*par_z[2]))/seconds_per_rad/36525.0;
  theta=t*(par_theta[0]+t*(par_theta[1]+t*par_theta[2]))/seconds_per_rad;
  dtheta=(par_theta[0]+t*(2.0*par_theta[1]+t*3.0*par_theta[2]))/seconds_per_rad/36525.0;

  czeta = cos(zeta);
  szeta = sin(zeta);
  dczeta = -szeta*dzeta;
  dszeta = czeta*dzeta;
  cz = cos(z);
  sz = sin(z);
  dcz = -sz*dz;
  dsz = cz*dz;
  ctheta = cos(theta);
  stheta = sin(theta);
  dctheta = -stheta*dtheta;
  dstheta = ctheta*dtheta;

  prc[0][0] = czeta*ctheta*cz - szeta*sz;
  prc[1][0] = czeta*ctheta*sz + szeta*cz;
  prc[2][0] = czeta*stheta;
  prc[0][1] = -szeta*ctheta*cz - czeta*sz;
  prc[1][1] = -szeta*ctheta*sz + czeta*cz;
  prc[2][1] = -szeta*stheta;
  prc[0][2] = -stheta*cz;
  prc[1][2] = -stheta*sz;
  prc[2][2] = ctheta;

  dprecess[0][0] = dczeta*ctheta*cz + czeta*dctheta*cz + czeta*ctheta*dcz - dszeta*sz - szeta*dsz;
  dprecess[1][0] = dczeta*ctheta*sz + czeta*dctheta*sz + czeta*ctheta*dsz + dszeta*cz + szeta*dcz;
  dprecess[2][0] = dczeta*stheta + czeta*dstheta;
  dprecess[0][1] = - dszeta*ctheta*cz - szeta*dctheta*cz - szeta*ctheta*dcz - dczeta*sz - czeta*dsz;
  dprecess[1][1] = - dszeta*ctheta*sz - szeta*dctheta*sz - szeta*ctheta*dsz + dczeta*cz + czeta*dcz;
  dprecess[2][1] = - dszeta*stheta - szeta*dstheta;
  dprecess[0][2] = - dstheta*cz - stheta*dcz;
  dprecess[1][2] = - dstheta*sz - stheta*dsz;
  dprecess[2][2] = dctheta;

  /* END OF PRECESSION.F */

  /* Start of NUTATION.F */
  ceps=cos(eps);
  seps=sin(eps);

  nut[0][0]=1.0;
  nut[0][1]=-delp*ceps;
  nut[0][2]=-delp*seps;
  nut[1][0]=-nut[0][1];
  nut[1][1]=1.0;
  nut[1][2]=-dele;
  nut[2][0]=-nut[0][2];
  nut[2][1]=-nut[1][2];
  nut[2][2]=1.0;

  dt = delp*cos(eps+dele);
 
  /* From PRCNUT.f */

  for (i=0;i<3;i++)
    {
      for (j=0;j<3;j++)
	  prn[j][i] = nut[i][0]*prc[0][j] + nut[i][1]*prc[1][j] + nut[i][2]*prc[2][j];
    }
  return 0.0;
}




/* COPIED FROM TEMPO, ANG.F */
double ang(int i,double f)
{
  int ia,ib;
  double g,ang;

  ia = (int)floor(f/1.0e4);
  ib = (int)floor((f-ia*1e4)/1.0e2);
  g = f-ia*1.0e4-ib*1.0e2;
  ang = (ia+(ib+g/6.0e1)/6.0e1)/36.0e1;
  if (i>2) ang=ang*15.0;
  if (i==1 || i==3) ang=ang*M_PI*2;
  return ang;

}
