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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"
#include "ifunc.h"
#include "shapelet.h"
#include "t2fit_solarwind.h"


void dm_delays(pulsar *psr,int npsr,int p,int i,double delt,double dt_SSB)
{
    double ctheta,freqf,pospos,r,voverc,dmval,dmval2;
    double rsa[3],pos[3],vobs[3],yrs;
    double dmDot;
    longdouble dt;
    int j,k;
    const char *CVS_verNum = "$Id$";
    if (displayCVSversion == 1) CVSdisplayVersion("dm_delays.C","dm_delays()",CVS_verNum);
    logdbg("dm_delays with pulsar %d; number of obs = %d",p,psr[p].nobs);

    if (psr[p].obsn[i].delayCorr==0) /* No correction */
    {
        psr[p].obsn[i].tdis1 = 0.0;
        psr[p].obsn[i].tdis2 = 0.0;
    }
    else
    {
        for (k=0;k<3;k++)
            pos[k] = psr[p].posPulsar[k]+delt*psr[p].velPulsar[k];      
        pospos = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
        for (k=0;k<3;k++)
            pos[k] /= pospos;
        //      logdbg("with pos = %f %f %f",pos[0],pos[1],pos[2]);      

        /* Calculate position of observatory from Sun */
        for (j=0;j<3;j++)
        {
            /*	  rsa[j] = -psr[p].obsn[i].sun_ssb[j] + psr[p].obsn[i].earthMoonBary_ssb[j] -
                  psr[p].obsn[i].earthMoonBary_earth[j] + psr[p].obsn[i].observatory_earth[j]; */
            if (strcmp(psr[p].obsn[i].telID,"STL_FBAT")==0)
                rsa[j] = -psr[p].obsn[i].sun_ssb[j] + psr[p].obsn[i].observatory_earth[j];
            else
                rsa[j] = -psr[p].obsn[i].sun_ssb[j] + psr[p].obsn[i].earth_ssb[j] + psr[p].obsn[i].observatory_earth[j];
        }
        //      logdbg("with rsa = %f %f %f",rsa[0],rsa[1],rsa[2]);            
        //      printf("with rsa = %s %f %f %f\n",psr[p].obsn[i].fname,rsa[0],rsa[1],rsa[2]);            
        /* What about Sun from SSB? */
        for (j=0;j<3;j++)
        {
            if (strcmp(psr[p].obsn[i].telID,"STL_FBAT")==0)
            {
                if (psr[p].param[param_telx].paramSet[1] == 1 &&
                        psr[p].param[param_tely].paramSet[1] == 1 &&
                        psr[p].param[param_telz].paramSet[1] == 1)
                {
                    if (j==0)      vobs[j] = (double)psr[p].param[param_telx].val[1];
                    else if (j==1) vobs[j] = (double)psr[p].param[param_tely].val[1];
                    else if (j==2) vobs[j] = (double)psr[p].param[param_telz].val[1];
                    //		  printf("Setting vobs = %g\n",vobs[j]);
                }
                else
                {
                    logerr("setting vobs = 0, this can give a large error!");
                    vobs[j] = 0.0;
                }
            }
            else
                vobs[j] = psr[p].obsn[i].earth_ssb[j+3] + psr[p].obsn[i].siteVel[j];
        }
        /*	vobs[j] = psr[p].obsn[i].earthMoonBary_ssb[j+3] - psr[p].obsn[i].earthMoonBary_earth[j+3] + 
            psr[p].obsn[i].siteVel[j];*/
        //      logdbg("with vobs = %f %f %f",vobs[0],vobs[1],vobs[2]);      
        // printf("with vobs = %s %g %g %g\n",psr[p].obsn[i].fname,vobs[0],vobs[1],vobs[2]);      
        r = sqrt(dotproduct(rsa,rsa));
        ctheta = dotproduct(pos,rsa)/r;
        voverc = dotproduct(pos,vobs);
        freqf = psr[p].obsn[i].freq*1.0e6*(1.0-voverc);   /* Converting frequency in to Hz */
        //      printf("freqf = %g %g %g %g\n",freqf,psr[p].obsn[i].freq,voverc,psr[p].obsn[i].einsteinRate);
        /* Transform freq due to Einstein delay ! */
        if (psr[p].dilateFreq && freqf > 0 && psr[p].obsn[i].einsteinRate != 0.0)
            freqf /= psr[p].obsn[i].einsteinRate;

        //      logdbg("Transforming frequency due to Einstein delay");      

        psr[p].obsn[i].freqSSB = freqf; /* Record observing frequency in barycentric frame */
        //      logdbg("set freqSSB");      
        dt = psr[p].obsn[i].sat - psr[p].param[param_dmepoch].val[0];
        yrs = dt/365.25;
        dmDot=0.0;

        longdouble arg = 1.0;
        // redwards changed to avoid using slow calls to pow
        // Note this is not a Taylor Series - the Edwards paper says that it is!
        // 2020-06 MJK: Added support for having a Taylor series!
        double series_fac=1.0;
        for (k=1;k<9;k++)
        {
            arg *= yrs;
            if (psr->dm_series_type == series_taylor_pn) {
                series_fac *= k;
            }
            if (psr[p].param[param_dm].paramSet[k]==1)
                dmDot+=(double)(psr[p].param[param_dm].val[k]*arg/series_fac); 
        }
        //      logdbg("calculated dmDot %Lg",psr[p].param[param_dm].val[0]);
        dmval = psr[p].param[param_dm].val[0]; //+dmDot;

        long double cmval=0;
        //long double cmDot;
        arg=1.0;
        for (k=1;k<9;k++)
        {
            arg *= yrs;
            if (psr[p].param[param_cm].paramSet[k]==1)
                cmval+=(double)(psr[p].param[param_cm].val[k]*arg); 
        }


        //      logdbg("calculating dmval");      
        // NOT DONE ANYMORE:      dmval += psr[p].obsn[i].phaseOffset;  /* In completely the wrong place - phaseoffset is actually DM offset */



        /* Are we using DM values using DMMODEL parameter? */
        if (psr[p].param[param_dmmodel].paramSet[0] == 1)
        {
            static int t=0xFFFF;
            dmval=0;
            //double meanDM=(double)psr[p].param[param_dm].val[0];
            double meanDM=getParameterValue(&psr[p],param_dmmodel,0);

            if (psr[p].param[param_dm].paramSet[0] == 1 && psr->param[param_dmmodel].linkTo[0] != param_dm)
            {
                if(t & 0x000F){
                    logerr("Using DMMODEL value for dispersion measure instead of the DM parameter");
                    logerr("Please set 'DMMODEL DM' in the par file rather than specifying a value");
                    t= t & 0xFFF0;
                }
                //	      psr[p].param[param_dmmodel].val[0]=psr[p].param[param_dm].val[0];
            }

            // set the mean DMOFF to zero.
            if(t&0x00F0){
                double mean_dmoff=0;
                for (k=0;k<psr[p].dmoffsDMnum;k++){
                    mean_dmoff+=psr[p].dmoffsDM[k];
                }
                mean_dmoff /= (double)psr[p].dmoffsDMnum;
                if(mean_dmoff > 1e-12){
                    printf("WARNING: Mean DMOFF = %lg\n",mean_dmoff);
                    /*			  for(k=0;i < psr[p].nconstraints; k++){
                                  if (psr[p].constraints[k]==constraint_dmmodel_mean){
                                  printf("         Removing mean to keep DMMODEL constraint correct.\n");
                                  for (k=0;k<psr[p].dmoffsNum;k++){
                                  psr[p].dmoffsDM[k] -= mean_dmoff;
                                  }
                                  break;
                                  }
                                  }*/
                }
                t = t & 0xFF0F;
            }

            //dmoffs is basically an ifunc
            dmval = meanDM + ifunc(psr[p].dmoffsDM_mjd,psr[p].dmoffsDM,(double)psr[p].obsn[i].sat,psr[p].dmoffsDMnum);

        }
        else
        {
            /* Set using DMX ranges */
            for (k=0; k<psr[p].ndmx; k++) 
            {
                if ((psr[p].obsn[i].sat > psr[p].param[param_dmxr1].val[k])
                        && (psr[p].obsn[i].sat < psr[p].param[param_dmxr2].val[k]))
                    dmval += psr[p].param[param_dmx].val[k];
            }

            /* Check for flag to set dm to exact value */
            for (k=0;k<psr[p].obsn[i].nFlags;k++)
            {
                if (strcmp(psr[p].obsn[i].flagID[k],"-dm")==0) /* Have DM correction */
                {
                    sscanf(psr[p].obsn[i].flagVal[k],"%lf",&dmval);
                    // The dmOffset value is never set and therefore the following line only adds confusion.
                    //                                                   JPWV, 08.05.2014
                    //dmval+=psr[p].dmOffset;
                }
                if (strcmp(psr[p].obsn[i].flagID[k],"-dmo")==0) /* Have DM correction */
                {
                    sscanf(psr[p].obsn[i].flagVal[k],"%lf",&dmval2);
                    // The dmOffset value is never set and therefore the following line only adds confusion.
                    //                                                   JPWV, 08.05.2014
                    //dmval+=(dmval2+psr[p].dmOffset);
                    dmval += dmval2; // JPWV, 08.05.2014
                }
            }
            //	  logdbg("Looked for flags");      
        }



        // IS there a TN DM Shapelet parameter to deal with?
        for (int iTNShape=0; iTNShape < psr[p].nTNShapeletEvents; ++iTNShape) {
            // We only want to do DM for ones with a spectral index is equal to 2
            if (psr[p].TNShapeletEvFScale[iTNShape] == 2) {
                double t = (double)psr[p].obsn[i].bat;
                if (t==0) t = (double)psr[p].obsn[i].sat; // I guess we approximate it?
                double shapedm = evaluateShapelet(psr->TNShapeletEvN[iTNShape],
                        psr->TNShapeletEvPos[iTNShape],
                        psr->TNShapeletEvWidth[iTNShape],
                        psr->TNShapeletEvCoef[iTNShape],
                        t);
                dmval += shapedm;
            }
        }

        // Add in derivative term
        dmval += dmDot;

        // Add in annual terms
        //

        if (psr[p].param[param_dm_sin1yr].paramSet[0]==1){
            dmval += psr[p].param[param_dm_sin1yr].val[0]*sin(2*M_PI/(365.25)*dt);
        }
        if (psr[p].param[param_dm_cos1yr].paramSet[0]==1){
            dmval += psr[p].param[param_dm_cos1yr].val[0]*cos(2*M_PI/(365.25)*dt);
        }


        if (freqf<=1) /* Have infinitive frequency */
            psr[p].obsn[i].tdis1 = 0.0;
        else{
            psr[p].obsn[i].tdis1 = dmval/DM_CONST/1.0e-12/freqf/freqf; 
            /* The lines below were used to test scattering variation instead of DM variation.
            */
            //	psr[p].obsn[i].tdis1 = psr[p].param[param_dm].val[0]/DM_CONST/1.0e-12/freqf/freqf; 
            //	psr[p].obsn[i].tdis1 += dmval/pow(freqf/1e6,4.4); 
        }
        // add in chromatic noise derivatives

        double gam;

        if (psr[p].TNChromIdx != 0)
        {
            if (freqf>=1)
            {


                gam=-psr[p].TNChromIdx;
                //fprintf(stderr, "TNChrom %.3e\n", psr[p].TNChromIdx);
                //exit(0);

                //psr[p].obsn[i].tdis1 += cmval/DM_CONST/1e-12/freqf/freqf;//*powl(freqf/1e6,gam);
                psr[p].obsn[i].tdis1 += cmval/DM_CONST*powl(freqf/1e6,gam);

            }
        }

        //      logdbg("calculate tdis1");      
        /* Add frequency dependent delay term */
        if (psr[p].param[param_fddc].paramSet[0]==1 && freqf>1)
            psr[p].obsn[i].tdis1 += (double)(psr[p].param[param_fddc].val[0]/pow(freqf*1.0e-6,(double)psr[p].param[param_fddi].val[0]));
        for (k=0; k<psr[p].param[param_fd].aSize; k++) {
            if (psr[p].param[param_fd].paramSet[k]==1 && freqf>1)
                psr[p].obsn[i].tdis1 += (double)(psr[p].param[param_fd].val[k] * 
                        pow(log(freqf/1e9),k+1));
        }
        if (freqf<=1.0 || psr[p].ipm==0){
            psr[p].obsn[i].tdis2 = 0.0;
            psr[p].obsn[i].spherical_solar_wind = 0.0; // there is no effect of the solar wind
        }
        else if (psr[p].swm==0) // Simple tempo2 solar wind model
        {
            /* TEMPO1: Note: the following line simplifies to 2e14 theta/r/sin(theta)/2/freqf^2 */
            /* psr[p].obsn[i].tdis2 = 2.0e14*acos(ctheta)/r/sqrt(1.0-ctheta*ctheta)/2.0/freqf/freqf; 
               printf("tdis %f %f %f %f %f %g %g\n",1.0e6,AU_DIST,SPEED_LIGHT,DM_CONST_SI,psr[p].ne_sw,1.0e6*AU_DIST*AU_DIST/SPEED_LIGHT/DM_CONST_SI*psr[p].ne_sw*
               acos(ctheta)/r/sqrt(1.0-ctheta*ctheta)/freqf/freqf,psr[p].obsn[i].tdis2); */

            // this is the conversion from ne_sw to solar wind contribution to the residual.
            psr[p].obsn[i].spherical_solar_wind = 1.0e6*AU_DIST*AU_DIST/SPEED_LIGHT/DM_CONST_SI*acos(ctheta)/r/sqrt(1.0-ctheta*ctheta)/freqf/freqf;

            psr[p].obsn[i].tdis2 = psr[p].ne_sw * psr[p].obsn[i].spherical_solar_wind; 

            if (psr[p].param[param_ne_sw_sin].paramSet[0]) {
                double ne_sw_sin = -cos(2*M_PI*(psr[p].obsn[i].sat - SOLAR_CYCLE_EPOCH)/SOLAR_CYCLE_PERIOD) * psr[p].param[param_ne_sw_sin].val[0];
                psr[p].obsn[i].tdis2 += ne_sw_sin * psr[p].obsn[i].spherical_solar_wind; 
            }

            if (psr[p].ne_sw_ifuncN > 0) {
                // we have an 'ifunc-like' solar wind model
                //
                double ne_sw_ifunc = ifunc(psr[p].ne_sw_ifuncT,psr[p].ne_sw_ifuncV,static_cast<double>(psr[p].obsn[i].sat),psr[p].ne_sw_ifuncN);
                psr[p].obsn[i].tdis2 += ne_sw_ifunc * psr[p].obsn[i].spherical_solar_wind;

            }
        }
        else if (psr[p].swm==1) // More complex solar wind model introduced by Xiaopeng You and William Coles
        {
            //	  if(acos(ctheta)*180/3.14159265>120)
            //	  printf("SW: calling observation %d\n",i);
            psr[p].obsn[i].tdis2 = (double)solarWindModel(psr+p,i)/DM_CONST/1.0e-12/freqf/freqf;
            //	    printf("Returning: %g\n",(double)psr[p].obsn[i].tdis2);
            //	  else
            //	    psr[p].obsn[i].tdis2 = 1.0e6*AU_DIST*AU_DIST/SPEED_LIGHT/DM_CONST_SI*psr[p].ne_sw*acos(ctheta)/r/sqrt(1.0-ctheta*ctheta)/freqf/freqf; 	    
        }

        //	  logdbg("[%d/%d] with tdis2 = %g, freqf = %g %g %g %d %g",i,
        //	 psr[p].nobs,(double)psr[p].obsn[i].tdis2,(double)freqf,psr[p].obsn[i].freq,
        //	 (1.0-voverc),psr[p].dilateFreq,psr[p].obsn[i].einsteinRate);      

        /* psr[p].obsn[i].tdis = (psr[p].param[param_dm].val/2.41e-16 +
           2.0e14*acos(ctheta)/r/sqrt(1.0-ctheta*ctheta)/2.0)/freqf/freqf; */
        /* printf("Dispersion delay = %g\n",(double)(psr[p].obsn[i].tdis1+psr[p].obsn[i].tdis2));  */

    }
    //        logdbg("Exiting dm_delays with pulsar %d; number of obs = %d",p,psr[p].nobs);

}

