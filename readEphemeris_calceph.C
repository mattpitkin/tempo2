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
#include "ifteph.h"

void convertUnits(double *val,int units,int eclCoord);

#ifdef HAVE_CALCEPH
#include <calceph.h>

/*
 * readEphemeris_calceph
 *
 * ephemeris reading routine making use of the calceph library
 *
 * This routine sets:
 *  psr[p].obsn[i].sun_ssb[0->5]
 *  psr[p].obsn[i].sun_earth[0->5]
 *  psr[p].obsn[i].planet_ssb[0->9][0->5]
 *  psr[p].obsn[i].jupiter_earth[0->5]
 *  psr[p].obsn[i].saturn_earth[0->5]
 *  psr[p].obsn[i].venus_earth[0->5]
 *  psr[p].obsn[i].uranus_earth[0->5]
 *  psr[p].obsn[i].neptune_earth[0->5]
 *  psr[p].obsn[i].earth_ssb[0->5]
 * 
 * Also
 *  update SSB position for reflex motion due to error in planetary masses
 *  convert to ecliptic coordinates if requested
 *  set earth_ssb with TEL_DX, TEL_DY, TEL_DZ parameters are required
 */

void readEphemeris_calceph(pulsar *psr,int npsr)
{
    t_calcephbin *eph;
    int i,p;
    longdouble jd;
    double jd0,jd1;


    for (p=0;p<npsr;p++)
    {
        eph = calceph_open(psr[p].JPL_EPHEMERIS);
        if (eph) {
            printf("Successfully opened ephemeris >%s<\n",psr[p].JPL_EPHEMERIS);
        } else {
            printf("Error: unable to open ephemeris >%s< for pulsar >%s<\n",psr[p].JPL_EPHEMERIS, psr[p].name);
            exit(1);
        }
        // Now read the ephemeris for each observation
        for (i=0;i<psr[p].nobs;i++)
        {
            jd = psr[p].obsn[i].sat + getCorrectionTT(psr[p].obsn+i)/SECDAY + 
                psr[p].obsn[i].correctionTT_Teph/SECDAY+2400000.5; 
            jd0 = (double)((int)jd);
            jd1 = (double)(jd-(int)jd);

            // Calculate the Earth to SSB vector
            calceph_compute_unit(eph,jd0,jd1,3,12,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,psr[p].obsn[i].earth_ssb);
            convertUnits(psr[p].obsn[i].earth_ssb,psr[p].units, psr[p].eclCoord);

            // Calculate the Sun to SSB vector
            calceph_compute_unit(eph,jd0,jd1,11,12,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,psr[p].obsn[i].sun_ssb);
            convertUnits(psr[p].obsn[i].sun_ssb,psr[p].units,psr[p].eclCoord);




            int iplanet;

            for (iplanet=0; iplanet < 9; iplanet++)
            {
                if ((psr[p].param[param_dphaseplanet].paramSet[iplanet] == 1) || (psr[p].param[param_dmassplanet].paramSet[iplanet] == 1))
                {
                    //err_code = jpl_pleph(ephem, jd, iplanet+1, 12, 
                    //		       psr[p].obsn[i].planet_ssb[iplanet], 1);
                    // Convert to sec and lt-s/s from AU and AU/day

                    calceph_compute_unit(eph,jd0,jd1,iplanet+1,12,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,psr[p].obsn[i].planet_ssb[iplanet]);

                    convertUnits(psr[p].obsn[i].planet_ssb[iplanet],psr[p].units,psr[p].eclCoord);


                }
            }


            for (iplanet=0; iplanet < 9; iplanet++)
            {
                if ((psr[p].param[param_dphaseplanet].paramSet[iplanet] == 1)  || (psr[p].param[param_dmassplanet].paramSet[iplanet] == 1))
                {
                    // err_code = jpl_pleph(ephem, jd+1e-6, iplanet+1, 12, 
                    //		       psr[p].obsn[i].planet_ssb_tmr[iplanet], 1);

                    calceph_compute_unit(eph,jd0+1,jd1,iplanet+1,12,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,psr[p].obsn[i].planet_ssb_tmr[iplanet]);


                    convertUnits(psr[p].obsn[i].planet_ssb_tmr[iplanet],psr[p].units,psr[p].eclCoord);



                }
            }



            //  fprintf(planetfile, "%.5Le ", psr[p].obsn[i].sat);


            for(iplanet=0;iplanet<9;iplanet++)
            {
                psr[p].obsn[i].planet_ssb_derv[iplanet][0]=   psr[p].obsn[i].planet_ssb_tmr[iplanet][0]- psr[p].obsn[i].planet_ssb[iplanet][0];
                psr[p].obsn[i].planet_ssb_derv[iplanet][1]=  psr[p].obsn[i].planet_ssb_tmr[iplanet][1]- psr[p].obsn[i].planet_ssb[iplanet][1];
                psr[p].obsn[i].planet_ssb_derv[iplanet][2]=   psr[p].obsn[i].planet_ssb_tmr[iplanet][2]- psr[p].obsn[i].planet_ssb[iplanet][2];
                psr[p].obsn[i].planet_ssb_derv[iplanet][3]=   psr[p].obsn[i].planet_ssb_tmr[iplanet][3]- psr[p].obsn[i].planet_ssb[iplanet][3];
                psr[p].obsn[i].planet_ssb_derv[iplanet][4]=   psr[p].obsn[i].planet_ssb_tmr[iplanet][4]- psr[p].obsn[i].planet_ssb[iplanet][4];
                psr[p].obsn[i].planet_ssb_derv[iplanet][5]=   psr[p].obsn[i].planet_ssb_tmr[iplanet][5]- psr[p].obsn[i].planet_ssb[iplanet][5]; 
                //fprintf(planetfile, "%.5e %.5e %.5e ",psr[p].obsn[i].planet_ssb_tmr[iplanet][1], psr[p].obsn[i].planet_ssb[iplanet][1], psr[p].obsn[i].planet_ssb_derv[iplanet][1] );
            }
            //fprintf(planetfile, "\n");








            // calculate Jupiter to Earth vector
            calceph_compute_unit(eph,jd0,jd1,5,3,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,psr[p].obsn[i].jupiter_earth);
            convertUnits(psr[p].obsn[i].jupiter_earth,psr[p].units,psr[p].eclCoord);

            // calculate Saturn to Earth vector
            calceph_compute_unit(eph,jd0,jd1,6,3,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,psr[p].obsn[i].saturn_earth);
            convertUnits(psr[p].obsn[i].saturn_earth,psr[p].units,psr[p].eclCoord);

            // calculate Uranus to Earth vector
            calceph_compute_unit(eph,jd0,jd1,7,3,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,psr[p].obsn[i].uranus_earth);
            convertUnits(psr[p].obsn[i].uranus_earth,psr[p].units,psr[p].eclCoord);

            // Neptune-Earth
            calceph_compute_unit(eph,jd0,jd1,8,3,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,psr[p].obsn[i].neptune_earth);
            convertUnits(psr[p].obsn[i].neptune_earth,psr[p].units,psr[p].eclCoord);

            // Venus-earth
            calceph_compute_unit(eph,jd0,jd1,2,3,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,psr[p].obsn[i].venus_earth);
            convertUnits(psr[p].obsn[i].venus_earth,psr[p].units,psr[p].eclCoord);

            // Earth-moon bary to SSB
            calceph_compute_unit(eph,jd0,jd1,13,12,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,psr[p].obsn[i].earthMoonBary_ssb);
            convertUnits(psr[p].obsn[i].earthMoonBary_ssb,psr[p].units,psr[p].eclCoord);

            //Earth-moon bary to earth
            calceph_compute_unit(eph,jd0,jd1,13,3,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,psr[p].obsn[i].earthMoonBary_earth);
            convertUnits(psr[p].obsn[i].earthMoonBary_earth,psr[p].units,psr[p].eclCoord);


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

        }	
        calceph_close(eph);
    }
}


void tt2tb_calceph(pulsar *psr, int npsr)
{
    t_calcephbin *eph;
    int i,p;
    long double jd;
    double jd0,jd1;
    double ttcorr[6];

    for (p=0;p<npsr;p++)
    {
        eph = calceph_open(psr[p].JPL_EPHEMERIS);
        if (eph) {
            printf("Successfully opened ephemeris >%s<\n",psr[p].JPL_EPHEMERIS);
        } else {
            printf("Error: unable to open ephemeris >%s< for pulsar >%s<\n",psr[p].JPL_EPHEMERIS,psr[p].name);
            exit(1);
        }
        // Now read the ephemeris for each observation
        for (i=0;i<psr[p].nobs;i++)
        {
            jd = psr[p].obsn[i].sat + getCorrectionTT(psr[p].obsn+i)/SECDAY + 
                psr[p].obsn[i].correctionTT_Teph/SECDAY+2400000.5; 
            jd0 = (double)((int)jd);
            jd1 = (double)(jd-(int)jd);

            // Calculate the Earth to SSB vector
            //calceph_compute_unit(eph,jd0,jd1,3,12,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,psr[p].obsn[i].earth_ssb);
            calceph_compute(eph,jd0,jd1,16,0,ttcorr);



            psr[p].obsn[i].correctionTT_calcEph=-ttcorr[0]; 

            //fprintf(stderr, "%.8e\n", ttcorr[0]);

            //  convertUnits(psr[p].obsn[i].earth_ssb,psr[p].units);

            // Calculate the Sun to SSB vector
            //calceph_compute_unit(eph,jd0,jd1,11,12,CALCEPH_UNIT_KM|CALCEPH_UNIT_SEC,psr[p].obsn[i].sun_ssb);
            // convertUnits(psr[p].obsn[i].sun_ssb,psr[p].units);
        }	
        calceph_close(eph);


    }
    //exit(0);
    return;
}
#else
void readEphemeris_calceph(pulsar *psr,int npsr)
{
    printf("ERROR: unable to use calceph library routines as library not installed\n");
    exit(1);
}
void tt2tb_calceph(pulsar *psr, int npsr)
{
    printf("ERROR: unable to use calceph library routines as library not installed\n");
    exit(1);
}
#endif

void convertUnits(double *val,int units,int eclCoord)
{
    double scale=1;
    if (units == SI_UNITS)
        scale = IFTE_K;
    val[0] = val[0]*1000.0*scale/SPEED_LIGHT;
    val[1] = val[1]*1000.0*scale/SPEED_LIGHT;
    val[2] = val[2]*1000.0*scale/SPEED_LIGHT;
    val[3] = val[3]*1000.0*scale/SPEED_LIGHT;
    val[4] = val[4]*1000.0*scale/SPEED_LIGHT;
    val[5] = val[5]*1000.0*scale/SPEED_LIGHT;
    if (eclCoord==1) equ2ecl(val);
    if (eclCoord==1) equ2ecl(val+3);
}
