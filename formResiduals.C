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
#include <string.h>
#include <math.h>
#include "tempo2.h"

/* Form the timing residuals from the timing model and the barycentric arrival times */
void residualTracking(pulsar *psr);

void formResiduals(pulsar *psr,int npsr,int removeMean)
{
   longdouble residual;  /* Residual in phase */
   longdouble ddnprd;
   longdouble dphase;
   longdouble nphase,phase5,phase2,phase3,phase4,lastResidual=0,priorResidual=0,ppRes=0;
   longdouble lastBat=0.0,priorBat=0.0,ppBat=0.0;
   longdouble phaseJ,phaseW;
   longdouble ftpd,fct,ff0,phaseint;
   longdouble torb,deltaT,dt00=0.0,phas1=0.0;
   longdouble mean,ct00=0.0;
   int nmean;
   int ntpd,nf0;
   int i,p,k;
   int time=0;
   int ntrk=0;
   int gotit=0;

   for (p=0;p<npsr;p++)
     {
       mean = 0.0;
       nmean = 0;

       for (i=0;i<psr[p].nobs;i++)
	 {
	   ddnprd = 0.0; /* ALWAYS ?? */
	   dphase = 0.0; /*    "      */
	   
	   torb = 0.0; 

	   /* Binary parameters */
	   if (psr[p].param[param_pb].paramSet[0]==1) 
	     {
	       if (strcmp(psr[p].binaryModel,"BT")==0)         torb = BTmodel(psr,p,i,-1);
	       else if (strcmp(psr[p].binaryModel,"BTJ")==0)   torb = BTJmodel(psr,p,i,-1,0);
	       else if (strcmp(psr[p].binaryModel,"ELL1")==0)  torb = ELL1model(psr,p,i,-1);
	       else if (strcmp(psr[p].binaryModel,"DD")==0)    torb = DDmodel(psr,p,i,-1);
	       else if (strcmp(psr[p].binaryModel,"DDK")==0)   torb = DDKmodel(psr,p,i,-1);
	       else if (strcmp(psr[p].binaryModel,"DDS")==0)   torb = DDSmodel(psr,p,i,-1);
	       else if (strcmp(psr[p].binaryModel,"MSS")==0)   torb = MSSmodel(psr,p,i,-1);
	       else if (strcmp(psr[p].binaryModel,"DDGR")==0)  torb = DDGRmodel(psr,p,i,-1);
	       else if (strcmp(psr[p].binaryModel,"T2")==0)    torb = T2model(psr,p,i,-1,0);
	       else {printf("Warning: Unknown binary model '%s'\n",psr[p].binaryModel); exit(1);}
	     }
	   psr[p].obsn[i].torb = torb; // save for plotting etc
	   deltaT = (psr[p].obsn[i].bbat-psr[p].param[param_pepoch].val[0])*86400.0+torb;
	   fct = psr[p].obsn[i].bbat-(int)psr[p].obsn[i].bbat - (psr[p].param[param_pepoch].val[0] - 
							       (int)psr[p].param[param_pepoch].val[0]);	   
	   ftpd = fct+torb/86400.0;	   
	   psr[p].obsn[i].pet = psr[p].obsn[i].bbat - torb/SECDAY;

	   /* What is the extra -ntpd*(int)(psr.f0) term on the end??? */
	   /* The extra (not included integral parts (int)f0.(int)bat) are not necessary! */
	   /* ntpd = (int)psr[p].obsn[i].bat-(int)psr[p].param[param_pepoch].val[0]; */
	   /*	   phase2 = deltaT*(psr[p].param[param_f0].val)-(ntpd*(int)(psr[p].param[param_f0].val)*86400.0); */

	   nf0  = (int)psr[p].param[param_f].val[0];
	   ntpd = ((int)psr[p].obsn[i].bbat-(int)psr[p].param[param_pepoch].val[0]);
	   ff0  = psr[p].param[param_f].val[0]-(int)psr[p].param[param_f].val[0];
	   phase2 = (nf0*ftpd+ntpd*ff0+ftpd*ff0); 
	   phase2 *= 86400.0;
	   /* redwards, all these calls to pow are slow & imprecise. changed */
	   longdouble arg = deltaT*deltaT;
	   phase3 = 0.5*psr[p].param[param_f].val[1]*arg; 
	   arg *= deltaT; if (psr[p].param[param_f].paramSet[2]==1) phase3 += (psr[p].param[param_f].val[2]/6.0L)*arg;
	   arg *= deltaT; if (psr[p].param[param_f].paramSet[3]==1) phase3 += (psr[p].param[param_f].val[3]/24.0L)*arg;
	   arg *= deltaT; if (psr[p].param[param_f].paramSet[4]==1) phase3 += (psr[p].param[param_f].val[4]/120.0L)*arg;
	   arg *= deltaT; if (psr[p].param[param_f].paramSet[5]==1) phase3 += (psr[p].param[param_f].val[5]/720.0L)*arg;
	   arg *= deltaT; if (psr[p].param[param_f].paramSet[6]==1) phase3 += (psr[p].param[param_f].val[6]/5040.0L)*arg;
	   arg *= deltaT; if (psr[p].param[param_f].paramSet[7]==1) phase3 += (psr[p].param[param_f].val[7]/40320.0L)*arg;
	   arg *= deltaT; if (psr[p].param[param_f].paramSet[8]==1) phase3 += (psr[p].param[param_f].val[8]/362880.0L)*arg;
	   arg *= deltaT; if (psr[p].param[param_f].paramSet[9]==1) phase3 += (psr[p].param[param_f].val[9]/3628800.0L)*arg; 
	   arg *= deltaT; if (psr[p].param[param_f].paramSet[10]==1) phase3 += (psr[p].param[param_f].val[10]/39916800.0L)*arg; 
	   arg *= deltaT; if (psr[p].param[param_f].paramSet[11]==1) phase3 += (psr[p].param[param_f].val[11]/479001600.0L)*arg; 
	   arg *= deltaT; if (psr[p].param[param_f].paramSet[12]==1) phase3 += (psr[p].param[param_f].val[12]/6227020800.0L)*arg; 

	   /* Must check glitch parameters */
	   phase4 = 0.0; 
	   for (k=0;k<psr[p].param[param_glep].aSize;k++)
	     {
	       if (psr[p].param[param_glep].paramSet[k]==1)
		 {
		   longdouble tp,dt1,expf,tgl;
		   tp = (ntpd+ftpd)*86400.0;	
		   tgl = (psr[p].param[param_glep].val[k] - psr[p].param[param_pepoch].val[0])*86400.0;
		   /* Only has effect after the glitch epoch */
		   if (tp >= tgl) 
		     {
		       dt1 = tp-tgl;
		       if (psr[p].param[param_gltd].val[k]!=0.0)
			 expf = exp(-dt1/86400.0/psr[p].param[param_gltd].val[k]);
		       else
			 expf = 1.0;
		       
		       // What happens if GLF2 (or GLF1) is not set?
		       phase4+=psr[p].param[param_glph].val[k]+
			 psr[p].param[param_glf0].val[k]*dt1 + 
			 0.5*psr[p].param[param_glf1].val[k]*dt1*dt1 +
			 1.0/6.0*psr[p].param[param_glf2].val[k]*dt1*dt1*dt1 +
			 psr[p].param[param_glf0d].val[k]*
			 psr[p].param[param_gltd].val[k]*86400.0*(1.0-expf);
		     }
		 }
	     }
	   
	   /* Add in extra phase due to jumps */
	   phaseJ = 0.0;
	   for (k=1;k<=psr[p].nJumps;k++)	    
	     {
	       if (psr[p].obsn[i].jump==k)
		 phaseJ+=psr[p].jumpVal[k]*psr[p].param[param_f].val[0];
	     }

	   /* Add in extra phase due to whitening procedures */
	   phaseW = 0.0;
	   if (psr[p].param[param_wave_om].paramSet[0] == 1)
	     {
	       double om;  /* Fundamental frequency */
	       double dt;  /* Change in time from pepoch */

	       dt = psr[p].obsn[i].bbat - psr[p].param[param_pepoch].val[0];
	       om = psr[p].param[param_wave_om].val[0];
	       for (k=0;k<psr[p].nWhite;k++)
		 {
		   phaseW += psr[p].wave_sine[k]*sin(om*(k+1)*dt)*psr[p].param[param_f].val[0] +
		     psr[p].wave_cos[k]*cos(om*(k+1)*dt)*psr[p].param[param_f].val[0];
		 }
	     }

	   phase5 = phase2+phase3+phase4+phaseJ+phaseW;

	   if (psr[p].obsn[i].nFlags>0) /* Look for extra factor to add to residuals */
	     {
	       int k;
	       longdouble extra;

	       for (k=0;k<psr[p].obsn[i].nFlags;k++)
		 {
		   if (strcmp(psr[p].obsn[i].flagID[k],"-radd")==0)
		     {
		       sscanf(psr[p].obsn[i].flagVal[k],"%Lf",&extra);
		       /* psr[p].obsn[i].residual+=extra; */
		       phase5+=(extra*psr[p].param[param_f].val[0]);
		     }
		 }
	     }
	   
	   /* Set the first residual to be zero */
	   /* in this case fortran_mod) returns [phase5 - (int)(phase5)] 
            * ie. returns the fractional part of phase5
            */
	   if (psr[p].param[param_iperharm].paramSet[0]==1)
	     {  /* This code has been added for observations taken at the wrong rotational period */
	       phase5 = fortran_mod(phase5,1.0/psr[p].param[param_iperharm].val[0]);
	       if (phase5 >= 1.0/2.0/psr[p].param[param_iperharm].val[0])
		 {
		   while (phase5 >= 1.0/2.0/psr[p].param[param_iperharm].val[0])
		     phase5 -= 1.0/psr[p].param[param_iperharm].val[0];
		 }
	       else if (phase5 < -1.0/2.0/psr[p].param[param_iperharm].val[0])
		 {
		   while (phase5 < -1.0/2.0/psr[p].param[param_iperharm].val[0])
		     phase5 += 1.0/psr[p].param[param_iperharm].val[0];
		 }
	       if (i==0) phas1  = phase5;
	       phase5 = phase5-phas1; 	       
	     }
	   else
	     {
	       if (i==0) phas1  = fortran_mod((phase5),1.0L); 
	       phase5 = phase5-phas1; 
	     }
	   
	   /* phase5 - nphase is also the fractional part of phase5
            * -- this is different to using fortran_mod() as above as the nearest integer
            *    to a negative or a positive number is calculated differently. 
	    */
	   nphase = (longdouble)fortran_nint(phase5);

	   residual = phase5 - nphase + dphase+ddnprd;
	   psr[p].obsn[i].nphase = nphase;
	   /* residual = residual in phase */
	   if (psr[p].obsn[i].deleted!=1)
	     {
	       gotit = 1;
	     }

	   if (psr[p].param[param_track].paramSet[0]==1 && psr[p].param[param_track].val[0] != 0 && gotit==1 && time==1)
	     {
	       residual+=ntrk;
	       // Note that this requires that the points be in time order
	       if (psr[p].obsn[i].bbat - ct00 < 0.0L && fabs(psr[p].obsn[i].bbat-ct00) > 1)
		 {
		   printf("ERROR: Points must be in time order for tracking to work\n");
		   printf("Pulsar = %s\n",psr[p].name);
		   printf("Observation %d (%s) has BBAT = %.5g\n",i,psr[p].obsn[i].fname,(double)psr[p].obsn[i].bbat);
		   printf("Observation %d (%s) has BBAT = %.5g\n",i-1,psr[p].obsn[i-1].fname,(double)psr[p].obsn[i-1].bbat);
		   printf("Difference = %Lg %Lg %Lg\n",psr[p].obsn[i].bbat, ct00,psr[p].obsn[i].bbat - ct00);
		   exit(1);
		 }
	       if (psr[p].obsn[i].bbat - ct00 < fabs(psr[p].param[param_track].val[0]))
		 {
		   if (fabs(residual+1.0-dt00) < fabs(residual-dt00))
		     {
		       residual+=1.0;
		       ntrk+=1;
		     } 
		   else if (fabs(residual-1.0-dt00) < fabs(residual-dt00))
		     {
		       residual-=1.0;
		       ntrk-=1;
		     } 
		 }
	     }
	   if (psr[p].param[param_track].val < 0) // Do extra tracking
	     {
	     }

	   dt00  = residual;
	   ct00 = psr[p].obsn[i].bbat;
	   
	   if (gotit==1 && time==0)
	     time=1;
	   psr[p].obsn[i].residual = residual/psr[p].param[param_f].val[0];

	   ppRes = priorResidual;
	   priorResidual=lastResidual;
	   ppBat = priorBat;
	   priorBat = lastBat;
	   lastResidual=residual;
	   lastBat = psr[p].obsn[i].bat;

 	   /* Check for phase offsets */
	   for (k=0;k<psr[p].nPhaseJump;k++)
	     {
	       //	       if (psr[p].obsn[i].sat > psr[p].phaseJump[k])
	       //	       printf("Comparing with %g\n",(double)psr[p].obsn[psr[p].phaseJumpID[k]].sat);
	       if (psr[p].obsn[i].sat > psr[p].obsn[psr[p].phaseJumpID[k]].sat)
		 psr[p].obsn[i].residual += (double)psr[p].phaseJumpDir[k]/psr[p].param[param_f].val[0];
	     }

	   phaseint = nf0*ntpd*86400.0;
	   /*	   psr[p].obsn[i].phase = nint(phase5)+phaseint; */
	   /* dt in resid.f */
	   /*	   psr[p].obsn[i].phase = phaseint+nphase+phase5-nphase+dphase+ddnprd; */	   
	   psr[p].obsn[i].phase = phaseint+phase5;  
	   if (psr[p].obsn[i].deleted!=1)
	     {
	       mean+=psr[p].obsn[i].residual;
	       nmean++;
	     }
	 }
       if (removeMean==1)
	 {
	   mean/=(long double)nmean;
	   for (i=0;i<psr[p].nobs;i++)
	     psr[p].obsn[i].residual-=mean;
	 }

     }
}

