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

/* Template for a tempo2 plugin */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include "tempo2.h"
#include "TKmatrix.h"
#include "T2toolkit.h"
#include "t2fit.h"
#include "enum_str.h"
#include "cholesky.h"
#include "choleskyRoutines.h"


void help() /* Display help */
{
   /* This function should contain usage information about the plugin which should (in general) be accessed */
   /* by the user pressing 'h'                                                                              */
}

/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
   char parFile[MAX_PSR][MAX_FILELEN];
   char timFile[MAX_PSR][MAX_FILELEN];
   bool vary_y = true;
   bool vary_M = true;
   int i;
   int nloop=16;
   *npsr = 1;  /* For a graphical interface that only shows results for one pulsar */

   printf("Graphical Interface: fitStability\n");
   printf("Author:              Michael Keith\n");
   printf("Version:             v1.0\n");
   printf(" --- type 'h' for help information\n");

   /* Obtain all parameters from the command line */
   for (i=1;i<argc;i++)
   {
	  if (strcmp(argv[i],"-l")==0){
          nloop= atoi(argv[++i]);
      } else if (strcmp(argv[i],"-y")==0){
          vary_M=false;
      } else if (strcmp(argv[i],"-M")==0){
          vary_y=false;
      } else if (strcmp(argv[i],"-f")==0) {
		 strcpy(parFile[0],argv[++i]);
		 strcpy(timFile[0],argv[++i]);
	  }
   }



   long seed = TKsetSeed();

   readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
   readTimfile(psr,timFile,*npsr); /* Load the arrival times    */
   preProcess(psr,*npsr,argc,argv);


   formBatsAll(psr,*npsr);                /* Form Barycentric arrival times */
   formResiduals(psr,*npsr,1);       /* Form residuals */
   t2Fit(psr,*npsr,covarFuncFile);

   copyPSR(psr,0,1); // make a rough copy of the pre-fit state. This doesn't do a good job as it is a bit out of date.

   if (vary_y) { 
       // We add white noise to the 'pre-fit' to make it fair when we 
       // perturb the data later.
       for (int iobs =0 ; iobs < psr[0].nobs; ++ iobs) {
           longdouble eta = TKgaussDev(&seed)*psr[1].obsn[iobs].toaErr;
           psr[0].obsn[iobs].sat = psr[1].obsn[iobs].sat + eta/86400.0e6;
           psr[0].obsn[iobs].toaErr = psr[1].obsn[iobs].toaErr*sqrt(2);
       }
   }

   formBatsAll(psr,*npsr);                /* Form Barycentric arrival times */
   formResiduals(psr,*npsr,1);       /* Form residuals */
   t2Fit(psr,*npsr,covarFuncFile);



   double globalParameter=0;
   textOutput(psr,*npsr,globalParameter,0,0,0,"");


   FitInfo orig_fitinfo;
   memcpy(&orig_fitinfo,&(psr[0].fitinfo),sizeof(FitInfo));

   const int nparam = psr[0].fitinfo.nParams;
   double **parameter_cvm = malloc_uinv(nparam); // allocate parameter_cvm
   double **parameter_uinv = malloc_uinv(nparam); // allocate parameter_uinv

   for (int iparam =0; iparam < nparam ; ++iparam) {
       memcpy(parameter_cvm[iparam],psr[0].covar[iparam],nparam*sizeof(double));
   }
   cholesky_formLinv(parameter_uinv,parameter_cvm,nparam);

   double *ssq_deviations = new double[nparam];
   for (int iparam =0; iparam < nparam ; ++iparam) {
       ssq_deviations[iparam]=0;
   }

   double total_ssq=0;
   double total_ssqW=0;
   for (int iloop =0; iloop < nloop; ++iloop) {
       copyPSR(psr,1,0); // copy back over the original...

       if (vary_y) { 
           // perturb the input data...
           for (int iobs =0 ; iobs < psr[0].nobs; ++ iobs) {
               longdouble eta = TKgaussDev(&seed)*psr[1].obsn[iobs].toaErr; // make sure to use original error.
               psr[0].obsn[iobs].sat = psr[1].obsn[iobs].sat + eta/86400.0e6;
               psr[0].obsn[iobs].toaErr = psr[1].obsn[iobs].toaErr*sqrt(2);
           }
       }

       // perturb the initial parameters

       double *white_noise=new double[nparam];
       double *param_noise=new double[nparam];

       for (int iparam =0; iparam < nparam ; ++iparam) {
           white_noise[iparam] = TKgaussDev(&seed);
           param_noise[iparam] = 0;
       }
       if (vary_M) {
           // perturb the parameters...
           for (int iparam =0; iparam < nparam ; ++iparam) {
               double numerator=white_noise[iparam];
               for (int jparam =0; jparam < iparam ; ++jparam) {
                   numerator -= parameter_uinv[iparam][jparam] * param_noise[jparam];
                   //               logmsg("%d %d %lg %lg\n",iparam,jparam,parameter_uinv[jparam][iparam],parameter_uinv[iparam][jparam]);
               }
               param_noise[iparam] = numerator / parameter_uinv[iparam][iparam];
           }

           for (int iparam =0; iparam < nparam ; ++iparam) {
               param_label p = orig_fitinfo.paramIndex[iparam];
               int k = orig_fitinfo.paramCounters[iparam];
               switch(p) {
                   case param_red_sin:
                   case param_red_cos:
                   case param_red_dm_sin:
                   case param_red_dm_cos:
                   case param_band_red_sin:
                   case param_band_red_cos:
                   case param_group_red_sin:
                   case param_group_red_cos:
                   case param_jitter:
//                   case param_ZERO:
                       param_noise[iparam]=0; // these parameters don't get stored...
                       break;
                   default:
                       orig_fitinfo.updateFunctions[iparam](psr,0,p,k,param_noise[iparam],0);
                       break;
               }
           }
       }


       // redo the fit
       formBatsAll(psr,*npsr);                /* Form Barycentric arrival times */
       formResiduals(psr,*npsr,1);       /* Form residuals */
       t2Fit(psr,*npsr,covarFuncFile);

       double sumsq=0;
       double* deviations = new double[nparam];
       double* white_deviations = new double[nparam];

       for (int iparam =0; iparam < nparam ; ++iparam) {
           double abs_deviation = psr[0].fitinfo.output.parameterEstimates[iparam] - orig_fitinfo.output.parameterEstimates[iparam] + param_noise[iparam];
           deviations[iparam] = abs_deviation;
           double sig_deviation = abs_deviation / sqrt(parameter_cvm[iparam][iparam]);

           sumsq += sig_deviation*sig_deviation;
           ssq_deviations[iparam] += sig_deviation*sig_deviation;
       }
       total_ssq += sumsq;




       TKmultMatrixVec(parameter_uinv,deviations,nparam,nparam,white_deviations);

       double sumsqW=0;
       for (int iparam =0; iparam < nparam ; ++iparam) {
           sumsqW += white_deviations[iparam]*white_deviations[iparam];
       }
       total_ssqW += sumsqW;
       logmsg("Iteration %02d: RMS white deviation %lg, RMS devitation: %lg chisq: %lg",iloop,sqrt(sumsqW/double(nparam)),sqrt(sumsq/double(nparam)),psr[0].fitChisq/(double)psr[0].fitNfree);

       delete[] deviations;
       delete[] white_deviations;

   }






   logmsg("RMS Deviations per parameter");
   for (int iparam =0; iparam < nparam ; ++iparam) {
       param_label p1 = psr[0].fitinfo.paramIndex[iparam];
       int c1= psr[0].fitinfo.paramCounters[iparam];
       char n1[80];
       if (c1 < psr[0].param[p1].aSize && strlen(psr[0].param[p1].shortlabel[c1])>0) {
           strncpy(n1,psr[0].param[p1].shortlabel[c1],80);
       } else {
           snprintf(n1,80,"%s(%d)",label_str[p1],c1);
       }
       logmsg("% 2d %20s %.4f",iparam,n1,sqrt(ssq_deviations[iparam]/(double)nloop));
   }
   logmsg("");
   logmsg("Total RMS Deviation: %lg",sqrt(total_ssq/(double)(nloop*nparam)));
   logmsg("Total RMS White Deviation: %lg",sqrt(total_ssqW/(double)(nloop*nparam)));

   // free memory
   delete[] ssq_deviations;
   free_uinv(parameter_cvm);
   free_uinv(parameter_uinv);
}
