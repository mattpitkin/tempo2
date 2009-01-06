#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"

/* ******************************************** */
/* vectorPulsar                                 */
/* Author:  G. Hobbs (06 May 2003)              */
/* Purpose: Forms a 3-vector pointing at the    */
/*          pulsar                              */
/*                                              */
/* Inputs:                                      */
/* Outputs: Fills posPulsar and velPulsar       */
/*                                              */
/* Notes: based on pos_pmconv.f                 */
/*        NONE of this has been checked         */
/*        pospm_conv.f replaced space_motion.f  */
/*                                              */
/*        Do we have to change RAJ and DECJ     */
/*        to particular epoch? (NOT DONE)       */
/*        What about taking proper motion into  */
/*        account when calculating RAJ and DECJ */
/*        (NOT DONE)                            */
/*                                              */
/* Changes:                                     */
/* ******************************************** */
void vectorPulsar(pulsar *psr,int npsr)
{
  double ca,sa,cd,sd,convert,dec;
  double alpha,delta;
  int p;

  for (p=0;p<npsr;p++)
    {
      alpha = psr[p].param[param_raj].val[0];
      delta = psr[p].param[param_decj].val[0];
      ca = cos(alpha);
      sa = sin(alpha);
      cd = cos(delta);
      sd = sin(delta);
      
      psr[p].posPulsar[0] = ca*cd;
      psr[p].posPulsar[1] = sa*cd;
      psr[p].posPulsar[2] = sd;
      
      /* Conversion from mas/yr to radians of arc per century */
      convert = 1.0/1000.0/60.0/60.0*M_PI/180.0*100.0;
      dec = (double)psr[p].param[param_decj].val[0];
      psr[p].velPulsar[0] = convert*(-psr[p].param[param_pmra].val[0]/cos(dec)*sa*cd - 
				     psr[p].param[param_pmdec].val[0]*ca*sd);
      psr[p].velPulsar[1] = convert*(psr[p].param[param_pmra].val[0]/cos(dec)*ca*cd 
				     - psr[p].param[param_pmdec].val[0]*sa*sd);
      psr[p].velPulsar[2] = convert*(psr[p].param[param_pmdec].val[0]*cd);
    }
}


