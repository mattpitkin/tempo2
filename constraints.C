#include "constraints.h"
#include "TKfit.h"
#include "T2toolkit.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#define CONSTRAINT_WEIGHTS

std::string get_constraint_name(enum constraint c){
#ifdef CONSTRAINT_WEIGHTS
		printf("[WT]");
#endif
	switch(c){
		case constraint_dmmodel_mean:
			return "DMMODEL mean(DM) = 0";
		case constraint_dmmodel_dm1:
			return "DMMODEL mean(DM1) = 0";
		case constraint_dmmodel_cw_0:
			return "DMMODEL mean(C) = 0";
		case constraint_dmmodel_cw_1:
			return "DMMODEL linear(C) = 0";
		case constraint_dmmodel_cw_2:
			return "DMMODEL quadratic(C) = 0";
		case constraint_dmmodel_cw_3:
			return "DMMODEL cubic(C) = 0";
		case constraint_ifunc_0:
			return "IFUNC mean(C) = 0";
		case constraint_ifunc_1:
			return "IFUNC linear(C) = 0";
		case constraint_ifunc_2:
			return "IFUNC quadratic(C) = 0";
		case constraint_quad_ifunc_p_0:
			return "QIFUNC_p mean(C) = 0";
		case constraint_quad_ifunc_p_1:
			return "QIFUNC_p linear(C) = 0";
		case constraint_quad_ifunc_p_2:
			return "QIFUNC_p quadratic(C) = 0";
		case constraint_quad_ifunc_c_0:
			return "QIFUNC_c mean(C) = 0";
		case constraint_quad_ifunc_c_1:
			return "QIFUNC_c linear(C) = 0";
		case constraint_quad_ifunc_c_2:
			return "QIFUNC_c quadratic(C) = 0";
		case constraint_tel_dx_0:
			return "TEL_DX mean(C) = 0";
		case constraint_tel_dx_1:
			return "TEL_DX linear(C) = 0";
		case constraint_tel_dx_2:
			return "TEL_DX quadratic(C) = 0";
		case constraint_tel_dy_0:
			return "TEL_DY mean(C) = 0";
		case constraint_tel_dy_1:
			return "TEL_DY linear(C) = 0";
		case constraint_tel_dy_2:
			return "TEL_DY quadratic(C) = 0";
		case constraint_tel_dz_0:
			return "TEL_DZ mean(C) = 0";
		case constraint_tel_dz_1:
			return "TEL_DZ linear(C) = 0";
		case constraint_tel_dz_2:
			return "TEL_DZ quadratic(C) = 0";
		case constraint_dmmodel_cw_year_sin:
			return "DMMODEL_YEAR C.sin(t) = 0";
		case constraint_dmmodel_cw_year_cos:
			return "DMMODEL_YEAR C.cos(t) = 0";
		case constraint_dmmodel_cw_year_xsin:
			return "DMMODEL_YEAR C.t.sin(t) = 0";
		case constraint_dmmodel_cw_year_xcos:
			return "DMMODEL_YEAR C.t.cos(t) = 0";
		case constraint_dmmodel_cw_px:
			return "DMMODEL_YEAR C:PX = 0";
		case constraint_dmmodel_cw_year_sin2:
			return "DMMODEL_YEAR C.sin(2t) = 0";
		case constraint_dmmodel_cw_year_cos2:
			return "DMMODEL_YEAR C.cos(2t) = 0";
		case constraint_ifunc_year_sin:
			return "IFUNC_YEAR C.sin(t) = 0";
		case constraint_ifunc_year_cos:
			return "IFUNC_YEAR C.cos(t) = 0";
		case constraint_ifunc_year_xsin:
			return "IFUNC_YEAR C.t.sin(t) = 0";
		case constraint_ifunc_year_xcos:
			return "IFUNC_YEAR C.t.cos(t) = 0";
		case constraint_ifunc_year_sin2:
			return "IFUNC_YEAR C.sin(2t) = 0";
		case constraint_ifunc_year_cos2:
			return "IFUNC_YEAR C.cos(2t) = 0";
		case constraint_qifunc_p_year_sin:
			return "QIFUNC_p_YEAR C.sin(t) = 0";
		case constraint_qifunc_p_year_cos:
			return "QIFUNC_p_YEAR C.cos(t) = 0";
		case constraint_qifunc_p_year_xsin:
			return "QIFUNC_p_YEAR C.t.sin(t) = 0";
		case constraint_qifunc_p_year_xcos:
			return "QIFUNC_p_YEAR C.t.cos(t) = 0";
		case constraint_qifunc_p_year_sin2:
			return "QIFUNC_p_YEAR C.sin(2t) = 0";
		case constraint_qifunc_p_year_cos2:
			return "QIFUNC_p_YEAR C.cos(2t) = 0";
		case constraint_qifunc_c_year_sin:
			return "QIFUNC_c_YEAR C.sin(t) = 0";
		case constraint_qifunc_c_year_cos:
			return "QIFUNC_c_YEAR C.cos(t) = 0";
		case constraint_qifunc_c_year_xsin:
			return "QIFUNC_c_YEAR C.t.sin(t) = 0";
		case constraint_qifunc_c_year_xcos:
			return "QIFUNC_c_YEAR C.t.cos(t) = 0";
		case constraint_qifunc_c_year_sin2:
			return "QIFUNC_c_YEAR C.sin(2t) = 0";
		case constraint_qifunc_c_year_cos2:
			return "QIFUNC_c_YEAR C.cos(2t) = 0";

		default:
			return "UNKNOWN!";
	}
}



void matrixDMConstraintWeights(pulsar *psr){
		int i,j,k;
		int nobs=0;
		int nfit=psr->dmoffsDMnum+psr->dmoffsCMnum;
		int nDM=psr->dmoffsDMnum;

		double x;


		if (psr->param[param_dmmodel].fitFlag[0]==1)
		  {
		    logdbg("Getting DM constraints for %s",psr->name);
		    // find out how many obs we have.
		    for(i=0; i < psr->nobs; i++){
		      if (psr->obsn[i].deleted==0)
			{
			  char okay=1;
			  /* Check for START and FINISH flags */
			  if (psr->param[param_start].paramSet[0]==1 && psr->param[param_start].fitFlag[0]==1 &&
			      (psr->param[param_start].val[0] > psr->obsn[i].sat))
			    okay=0;
			  if (psr->param[param_finish].paramSet[0]==1 && psr->param[param_finish].fitFlag[0]==1 &&
			      psr->param[param_finish].val[0] < psr->obsn[i].sat)
			    okay=0;
			  if (okay==1)
			    {
			      nobs++;
			    }
			}
		    }
		// originally was sizeof(double*)*nobs.
		double ** designMatrix=(double**)malloc_uinv(nobs);
		double *e=(double*)malloc(sizeof(double)*nobs);

		nobs=0;
		//		printf("c0: %d %s\n",psr->nobs,psr->name);
		for(i=0; i < psr->nobs; i++){
		  // Check for "ok" and start/finish coppied from dofit.
		  // Needs to be the same so that the weighting is the same as the fit,
		  // but in practice probably doesn't matter.
		  //		  printf("c1 checking: %d %d %s\n",psr->obsn[i].deleted,psr->nobs,psr->name);
		  if (psr->obsn[i].deleted==0)
		    {
		      char okay=1;
		      
		      /* Check for START and FINISH flags */
		      if (psr->param[param_start].paramSet[0]==1 && psr->param[param_start].fitFlag[0]==1 &&
			  (psr->param[param_start].val[0] > psr->obsn[i].sat))
			okay=0;
		      if (psr->param[param_finish].paramSet[0]==1 && psr->param[param_finish].fitFlag[0]==1 &&
			  psr->param[param_finish].val[0] < psr->obsn[i].sat)
			okay=0;

		      //		  printf("c2 checking: %d\n",okay);
		      if (okay==1)
			{
			  x   = (double)(psr->obsn[i].bbat-psr->param[param_pepoch].val[0]);
			  double sig=psr->obsn[i].toaErr*1e-6;
			  double dmf = 1.0/(DM_CONST*powl(psr->obsn[i].freqSSB/1.0e6,2));
			  for (k=0; k < nfit;k++){
			    
			    if(k<nDM)
			      designMatrix[nobs][k]=dmf*getParamDeriv(psr,i,x,param_dmmodel,k)/sig;
			    else
			      designMatrix[nobs][k]=getParamDeriv(psr,i,x,param_dmmodel,k)/sig;
			  }
			  nobs++;
			  
			}
		    }
		} 
		printf("Calling TKleastSquares %d %d\n",nobs,nfit); fflush(stdout);
		TKleastSquares(NULL,NULL,designMatrix,designMatrix,nobs,nfit,1e-20,0,NULL,e,NULL);
		printf("Returning from calling TKleastSquares\n"); fflush(stdout);
		double sum_wDM=0;
		double sum_wCM=0;
		for (i=0;i<nfit;i++)
		  {
		    double sum=0.0;
		    if(i < nDM){
		      psr->dmoffsDM_weight[i]=1.0/e[i]/e[i];
		      printf("constraints.C here: %d %g\n",i,e[i]);
		      sum_wDM+=psr->dmoffsDM_weight[i];
		    } else {
		      psr->dmoffsCM_weight[i-nDM]=1.0/e[i]/e[i];
		      sum_wCM+=psr->dmoffsCM_weight[i-nDM];
		    }
		    
		    
		    
		  }
		//normalise the weights
		for (i=0;i<psr->dmoffsDMnum;i++)
		  psr->dmoffsDM_weight[i]/=sum_wDM;
		for (i=0;i<psr->dmoffsCMnum;i++)
		  psr->dmoffsCM_weight[i]/=sum_wCM;
		
		// free everything .
		free_uinv(designMatrix);
		free(e);
		  }
}


/*
 * Derive the weighting functions for the constraints, based upon the ToA errors.
 *
 */
void computeConstraintWeights(pulsar *psr){
  // printf("GH: in computeConstraintWeights with %s\n",psr->name);
   for (int k=0; k < psr->ifuncN; k++){
	  psr->ifunc_weights[k]=1.0/(double)psr->ifuncN;
   }   
   for (int k=0; k < psr->dmoffsDMnum; k++)
	  psr->dmoffsDM_weight[k]=1.0/(double)psr->dmoffsDMnum;

   for (int k=0; k < psr->dmoffsCMnum; k++)
	  psr->dmoffsCM_weight[k]=1.0/(double)psr->dmoffsCMnum;

#ifdef CONSTRAINT_WEIGHTS
   if(psr->dmoffsDMnum>0 || psr->dmoffsCMnum>0) {
	  matrixDMConstraintWeights(psr);
   }
   /*
	* Derive weights for ifuncs
	*/
   if(psr->ifuncN>0) {
     //  printf("GH: in computeConstraintWeights part 2 with %s\n",psr->name);
	  for (int i=0; i < psr->nobs; i++){
		 if (psr->obsn[i].deleted==0){
			// compute the weight that this ToA applies to each IFUNC
			for (int k=0;k<psr->ifuncN-1;k++)
			{
			   if ((double)psr->obsn[i].sat >= psr->ifuncT[k] &&
					 (double)psr->obsn[i].sat < psr->ifuncT[k+1])
			   {
				  double w=1.0/pow(psr->obsn[i].toaErr,2);
				  double dt=(psr->ifuncT[k+1]-psr->ifuncT[k]);
				  double t=(double)psr->obsn[i].sat-psr->ifuncT[k];
				  psr->ifunc_weights[k+1]+=w*t/dt;
				  psr->ifunc_weights[k]+=w*(1.0-t/dt);
				  break;
			   }
			}

		 }
	  }
	  // Normalise the weights
	  double sum=0;
	  for (int k=0; k < psr->ifuncN; k++){
		 sum+=psr->ifunc_weights[k];
	  }
	  for (int k=0; k < psr->ifuncN; k++){
		 psr->ifunc_weights[k]/=sum;
		 //  printf("GH: in computeConstraintWeights with %s setting %g\n",psr->name,psr->ifunc_weights[k]);
	  }
   }
#endif
   return;
}

double consFunc_dmmodel_mean(pulsar *psr,int i,int k){
   /*
	* Only operate on param=dmmodel and when fit parameter is 
	* one of the frequency dependant parts (i.e. first dmoffsNum)
	*/

   if(i==param_dmmodel && k < psr->dmoffsDMnum){
	  return psr->dmoffsDM_weight[k];
   } else return 0;
}

double consFunc_dmmodel_dm1(pulsar *psr,int i,int k){
   /*
	* Only operate on param=dmmodel and when fit parameter is 
	* one of the frequency dependant parts (i.e. first dmoffsNum)
	*/
  long double epoch = psr->param[param_dmepoch].val[0];
   int nDM=psr->dmoffsDMnum;
   printf("WE ARE IN HERE\n");
   if(i==param_dmmodel && k < nDM){
     //     k-=nDM;
     return psr->dmoffsDM_weight[k]*(psr->dmoffsDM_mjd[k]-epoch);
   } else return 0;
}

double consFunc_dmmodel_cw(pulsar *psr,int i,int k,int order){
   /*
	* Only operate on param=dmmodel and when fit parameter is 
	* one of the frequency independant parts (i.e. last dmoffsNum).
	*/
   int nDM=psr->dmoffsDMnum;
   if(i==param_dmmodel && k >= nDM){
	  k-=nDM;
	  long double epoch = psr->param[param_pepoch].val[0];
	  long double w=psr->dmoffsCM_weight[k];
	  return w*pow(psr->dmoffsCM_mjd[k]-epoch,order);
   } else return 0;

}

double consFunc_dmmodel_cw_year(pulsar *psr,int i,int k,int order){
   /*
	* Only operate on param=dmmodel and when fit parameter is 
	* one of the frequency independant parts (i.e. last dmoffsNum).
	*/
   int nDM=psr->dmoffsDMnum;
   if(i==param_dmmodel && k >= nDM){
	  k-=nDM;
	  double rc;
	  double s[3];
	  long double epoch = psr->param[param_pepoch].val[0];
	  long double t = psr->dmoffsCM_mjd[k]-epoch;
	  long double x = 2.0*M_PI*t/365.25;
	  long double w=psr->dmoffsCM_weight[k];
	  switch (order){
		 case 0:
			return w*sin(x);
		 case 1:
			return w*cos(x);
		 case 2:
			return w*x*sin(x);
		 case 3:
			return w*x*cos(x);
		 case 4:
			return w*sin(2*x);
		 case 5:
			return w*cos(2*x);
		 case 6:
			t = psr->dmoffsCM_mjd[k]-52995.0;
			x=2*M_PI*t/365.25;
			s[0]=-500*sin(x);
			s[1]=450*cos(x);
			s[2]=200*cos(x);
			rc=dotproduct(psr[0].posPulsar,s) / sqrt(dotproduct(psr[0].posPulsar,psr[0].posPulsar)*dotproduct(s,s));
			return w*rc;
		 default:
			return 0;
	  }
   } else return 0;

}


double consFunc_tel_dx(pulsar *psr,int i,int k,int order){
   if(i==param_tel_dx){
	  long double epoch = psr->param[param_pepoch].val[0];
	  return pow(psr->telDX_t[k]-epoch,order);
   }
   else return 0;
}

double consFunc_tel_dy(pulsar *psr,int i,int k,int order){
   if(i==param_tel_dy){
	  long double epoch = psr->param[param_pepoch].val[0];
	  return pow(psr->telDY_t[k]-epoch,order);
   }
   else return 0;
}

double consFunc_tel_dz(pulsar *psr,int i,int k,int order){
   if(i==param_tel_dz){
	  //    printf("In contraint with k = %d, order = %d\n",k,order); 
	  long double epoch = psr->param[param_pepoch].val[0];
	  return pow(psr->telDZ_t[k]-epoch,order);
   }
   else return 0;
}

double consFunc_ifunc(pulsar *psr,int i,int k,int order){
   /*
	* Only operate on param=ifunc and when fit parameter is 
	* one of the frequency independant parts (i.e. last ifuncN).
	*/
   if(i==param_ifunc){
	  long double epoch = psr->param[param_pepoch].val[0];
	  //	  printf("IFUNC: CONSTRAINT: %s %g %g %d\n",psr->name,psr->ifunc_weights[k],(double)epoch,order);
	  return psr->ifunc_weights[k]*pow(psr->ifuncT[k]-epoch,order);

   }
   else return 0;
}

double consFunc_qifunc_p_year(pulsar *psr,int i,int k,int order){
   /*
	* Only operate on param=dmmodel and when fit parameter is 
	* one of the frequency independant parts (i.e. last dmoffsNum).
	*/
   if(i==param_quad_ifunc_p && k < psr->quad_ifuncN_p){
	  long double epoch = psr->param[param_pepoch].val[0];
	  long double t = psr->quad_ifuncT_p[k%psr->quad_ifuncN_p]-epoch;
	  long double x = 2.0*M_PI*t/365.25;
	  switch (order){
		 case 0:
			return sin(x);
		 case 1:
			return cos(x);
		 case 2:
			return x*sin(x);
		 case 3:
			return x*cos(x);
		 case 4:
			return sin(2*x);
		 case 5:
			return cos(2*x);

		 default:
			return 0;
	  }
   } else return 0;
}

double consFunc_qifunc_c_year(pulsar *psr,int i,int k,int order){
   /*
	* Only operate on param=dmmodel and when fit parameter is 
	* one of the frequency independant parts (i.e. last dmoffsNum).
	*/
   if(i==param_quad_ifunc_c && k < psr->quad_ifuncN_c){
	  long double epoch = psr->param[param_pepoch].val[0];
	  long double t = psr->quad_ifuncT_c[k%psr->quad_ifuncN_c]-epoch;
	  long double x = 2.0*M_PI*t/365.25;
	  switch (order){
		 case 0:
			return sin(x);
		 case 1:
			return cos(x);
		 case 2:
			return x*sin(x);
		 case 3:
			return x*cos(x);
		 case 4:
			return sin(2*x);
		 case 5:
			return cos(2*x);
		 default:
			return 0;
	  }
   } else return 0;
}

double consFunc_ifunc_year(pulsar *psr,int i,int k,int order){
   /*
	* Only operate on param=dmmodel and when fit parameter is 
	* one of the frequency independant parts (i.e. last dmoffsNum).
	*/
   if(i==param_ifunc && k < psr->ifuncN){
	  long double epoch = psr->param[param_pepoch].val[0];
	  long double t = psr->ifuncT[k%psr->ifuncN]-epoch;
	  long double x = 2.0*M_PI*t/365.25;
	  switch (order){
		 case 0:
			return psr->ifunc_weights[k]*sin(x);
		 case 1:
			return psr->ifunc_weights[k]*cos(x);
		 case 2:
			return psr->ifunc_weights[k]*x*sin(x);
		 case 3:
			return psr->ifunc_weights[k]*x*cos(x);
		 case 4:
			return psr->ifunc_weights[k]*sin(2*x);
		 case 5:
			return psr->ifunc_weights[k]*cos(2*x);

		 default:
			return 0;
	  }
   } else return 0;
}


double consFunc_quad_ifunc_p(pulsar *psr,int i,int k,int order){
   /*
	* Only operate on param=ifunc and when fit parameter is 
	* one of the frequency independant parts (i.e. last ifuncN).
	*/
   if(i==param_quad_ifunc_p){
	  long double epoch = psr->param[param_pepoch].val[0];
	  return pow(psr->quad_ifuncT_p[k]-epoch,order);
   }
   else return 0;
}
double consFunc_quad_ifunc_c(pulsar *psr,int i,int k,int order){
   /*
	* Only operate on param=ifunc and when fit parameter is 
	* one of the frequency independant parts (i.e. last ifuncN).
	*/
   if(i==param_quad_ifunc_c){
	  long double epoch = psr->param[param_pepoch].val[0];
	  return pow(psr->quad_ifuncT_c[k]-epoch,order);
   }
   else return 0;
}


void autoConstraints(pulsar* psr_array, int ipsr,int npsr){
   pulsar* psr = psr_array+ipsr;
   psr->nconstraints=0;

   char dmmodel = psr->param[param_dmmodel].fitFlag[0] && psr->dmoffsDMnum > 1;
   char cmmodel = psr->param[param_dmmodel].fitFlag[0] && psr->dmoffsCMnum > 1;
   // If these are global fits, only add the constraints to the first pulsar.
   char ifunc = psr->param[param_ifunc].fitFlag[0]==1 || (ipsr==0 && psr->param[param_ifunc].fitFlag[0]==2);
   char qifunc_p = psr->param[param_quad_ifunc_p].fitFlag[0] ==1 || (ipsr==0 && psr->param[param_quad_ifunc_p].fitFlag[0] ==2);
   char qifunc_c = psr->param[param_quad_ifunc_c].fitFlag[0] ==1 || (ipsr==0 && psr->param[param_quad_ifunc_c].fitFlag[0] ==2);
   char tel_dx = psr->param[param_tel_dx].fitFlag[0] ==1 || (ipsr==0 && psr->param[param_tel_dx].fitFlag[0] ==2);
   char tel_dy = psr->param[param_tel_dy].fitFlag[0] ==1 || (ipsr==0 && psr->param[param_tel_dy].fitFlag[0] ==2);
   char tel_dz = psr->param[param_tel_dz].fitFlag[0] ==1 || (ipsr==0 && psr->param[param_tel_dz].fitFlag[0] ==2);
   

   // offset
   if (cmmodel) psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_0;
   if (dmmodel) psr->constraints[psr->nconstraints++] = constraint_dmmodel_mean;
   if (ifunc) psr->constraints[psr->nconstraints++] = constraint_ifunc_0;
   if (qifunc_p) psr->constraints[psr->nconstraints++] = constraint_quad_ifunc_p_0;
   if (qifunc_c) psr->constraints[psr->nconstraints++] = constraint_quad_ifunc_c_0;
   if (tel_dx) psr->constraints[psr->nconstraints++] = constraint_tel_dx_0;
   if (tel_dy) psr->constraints[psr->nconstraints++] = constraint_tel_dy_0;
   if (tel_dz) psr->constraints[psr->nconstraints++] = constraint_tel_dz_0;

   // F0
   if(psr->param[param_f].fitFlag[0]){
	  if(dmmodel)psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_1;
	  if (ifunc) psr->constraints[psr->nconstraints++] = constraint_ifunc_1;
	  if (qifunc_p) psr->constraints[psr->nconstraints++] = constraint_quad_ifunc_p_1;
	  if (qifunc_c) psr->constraints[psr->nconstraints++] = constraint_quad_ifunc_c_1;
	  if (tel_dx) psr->constraints[psr->nconstraints++] = constraint_tel_dx_1;
	  if (tel_dy) psr->constraints[psr->nconstraints++] = constraint_tel_dy_1;
	  if (tel_dz) psr->constraints[psr->nconstraints++] = constraint_tel_dz_1;
   }
   //F1
   if(psr->param[param_f].fitFlag[1]){
	  if(dmmodel)psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_2;
	  if (ifunc) psr->constraints[psr->nconstraints++] = constraint_ifunc_2;
	  if (qifunc_p) psr->constraints[psr->nconstraints++] = constraint_quad_ifunc_p_2;
	  if (qifunc_c) psr->constraints[psr->nconstraints++] = constraint_quad_ifunc_c_2;
	  if (tel_dx) psr->constraints[psr->nconstraints++] = constraint_tel_dx_2;
	  if (tel_dy) psr->constraints[psr->nconstraints++] = constraint_tel_dy_2;
	  if (tel_dz) psr->constraints[psr->nconstraints++] = constraint_tel_dz_2;

   }
   //F2
   if(psr->param[param_f].fitFlag[2]){
	  if(dmmodel)psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_3;
   }

   // RA/DEC
   if (psr->param[param_raj].fitFlag[0] || psr->param[param_decj].fitFlag[0]){
	  if(dmmodel){
		 psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_year_sin;
		 psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_year_cos;
	  }
	  if (ifunc){
		 psr->constraints[psr->nconstraints++] = constraint_ifunc_year_sin;
		 psr->constraints[psr->nconstraints++] = constraint_ifunc_year_cos;
	  }
	  if (qifunc_p){
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_p_year_sin;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_p_year_cos;
	  }
	  if (qifunc_c){
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_c_year_sin;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_c_year_cos;
	  }
   }
   if (psr->param[param_pmra].fitFlag[0] || psr->param[param_pmdec].fitFlag[0]){
	  if(dmmodel){
		 psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_year_xsin;
		 psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_year_xcos;
	  }
	  if (ifunc){
		 psr->constraints[psr->nconstraints++] = constraint_ifunc_year_xsin;
		 psr->constraints[psr->nconstraints++] = constraint_ifunc_year_xcos;
	  }
	  if (qifunc_p){
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_p_year_xsin;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_p_year_xcos;
	  }
	  if (qifunc_c){
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_c_year_xsin;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_c_year_xcos;
	  }
   }
   if (psr->param[param_px].fitFlag[0]){
	  if(dmmodel){
		 psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_year_sin2;
		 psr->constraints[psr->nconstraints++] = constraint_dmmodel_cw_year_cos2;
	  }
	  if (ifunc){
		 psr->constraints[psr->nconstraints++] = constraint_ifunc_year_sin2;
		 psr->constraints[psr->nconstraints++] = constraint_ifunc_year_cos2;
	  }
	  //
	  // Dustin Madison argues that constraints are not needed on PX
	  /*	  if (qifunc_p){
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_p_year_sin2;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_p_year_cos2;
	  }
	  if (qifunc_c){
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_c_year_sin2;
		 psr->constraints[psr->nconstraints++] = constraint_qifunc_c_year_cos2;
		 } */
   }


}



void autosetDMCM(pulsar* psr, double dmstep,double cmstep, double start, double end, bool fixCMgrid){
   int i,j;
   double mjd;
   bool ok;
   double psrstart,psrend;

   double D[MAX_IFUNC];
   double C[MAX_IFUNC];

   psr->dmoffsDMnum=0;
   psrstart=(double)(psr->obsn[0].sat)-cmstep/2.0;
   psrend=(double)(psr->obsn[psr->nobs-1].sat)+cmstep/2.0;
   i=0;
   mjd=start;
   while (mjd <= end){
	  if (mjd > psrstart && mjd < psrend){
		 psr->dmoffsCM_mjd[i]=mjd;
		 psr->dmoffsCM[i]=0;
		 i++;
	  }
	  mjd+=cmstep;
   }
   psr->dmoffsCMnum=i;
   if (psr->dmoffsCMnum > MAX_IFUNC){
	  logerr("too many CM steps!");
	  exit(1);
   }


   ok=fixCMgrid; // only do the CM if we are not fixing the grid
   while (!ok){
	  matrixDMConstraintWeights(psr);
	  ok=true;
	  double threshold = 0.05/(double)psr->dmoffsCMnum;
	  j=0;
	  mjd=start;
	  bool bb=false;
	  for (j=0; j < psr->dmoffsCMnum; j++){
		 C[j]=psr->dmoffsCM_weight[j];
	  }
	  double min=TKfindMin_d(C,psr->dmoffsCMnum);
	  if(min < threshold){
		 i=0;
		 logmsg("min=%lg / %lg",min,threshold);
		 for (j=0; j < psr->dmoffsCMnum; j++){
			psr->dmoffsCM_mjd[i]=psr->dmoffsCM_mjd[j];
			psr->dmoffsCM[i]=0;
			if(psr->dmoffsCM_weight[j]>min){
			   i++;
			   bb=false;
			} else{
			   if (bb) i++; // don't delete multiple points in a row
			   else logmsg("CM Skip %lf %lg",psr->dmoffsCM_mjd[j],psr->dmoffsCM_weight[i]);
			   ok=false;
			   bb=true;
			}
			mjd+=cmstep;
		 }
		 psr->dmoffsCMnum=i;
	  }
   }

   logmsg("psr: %s ncm=%d",psr->name,psr->dmoffsCMnum);


   psrstart=(double)(psr->obsn[0].sat)-dmstep/2.0;
   psrend=(double)(psr->obsn[psr->nobs-1].sat)+dmstep/2.0;

   i=0;
   mjd=start;
   while (mjd <= end){
	  if (mjd > psrstart && mjd < psrend){
		 psr->dmoffsDM_mjd[i]=mjd;
		 psr->dmoffsDM[i]=0;
		 i++;
	  }
	  mjd+=dmstep;
   }
   psr->dmoffsDMnum=i;
   if (psr->dmoffsDMnum > MAX_IFUNC){
	  logerr("too many DM steps!");
	  exit(1);
   }


   ok=false;
   while (!ok){
	  matrixDMConstraintWeights(psr);
	  ok=true;
	  double threshold = 0.01/(double)psr->dmoffsDMnum;
	  j=0;
	  mjd=start;
	  bool bb=false;
	  for (j=0; j < psr->dmoffsDMnum; j++){
		 D[j]=psr->dmoffsDM_weight[j];
	  }
	  double min=TKfindMin_d(D,psr->dmoffsDMnum);
	  if(min < threshold){
		 i=0;
		 logmsg("min=%lg / %lg",min,threshold);
		 for (j=0; j < psr->dmoffsDMnum; j++){
			psr->dmoffsDM_mjd[i]=psr->dmoffsDM_mjd[j];
			psr->dmoffsDM[i]=0;
			if(psr->dmoffsDM_weight[j]>min){
			   i++;
			   bb=false;
			} else{
			   if (bb) i++; // don't delete multiple points in a row
			   else logmsg("DM Skip %lf %lg",psr->dmoffsDM_mjd[j],psr->dmoffsDM_weight[i]);
			   ok=false;
			   bb=true;
			}
			mjd+=cmstep;
		 }
		 psr->dmoffsDMnum=i;
	  }
   }

   logmsg("psr: %s ndm=%d",psr->name,psr->dmoffsDMnum);

   ok=fixCMgrid; // only do the CM if we are not fixing the grid
   while (!ok){
	  matrixDMConstraintWeights(psr);
	  ok=true;
	  double threshold = 0.05/(double)psr->dmoffsCMnum;
	  j=0;
	  mjd=start;
	  bool bb=false;
	  for (j=0; j < psr->dmoffsCMnum; j++){
		 C[j]=psr->dmoffsCM_weight[j];
	  }
	  double min=TKfindMin_d(C,psr->dmoffsCMnum);
	  if(min < threshold){
		 i=0;
		 logmsg("min=%lg / %lg",min,threshold);
		 for (j=0; j < psr->dmoffsCMnum; j++){
			psr->dmoffsCM_mjd[i]=psr->dmoffsCM_mjd[j];
			psr->dmoffsCM[i]=0;
			if(psr->dmoffsCM_weight[j]>min){
			   i++;
			   bb=false;
			} else{
			   if (bb) i++; // don't delete multiple points in a row
			   else logmsg("CM Skip %lf %lg",psr->dmoffsCM_mjd[j],psr->dmoffsCM_weight[i]);
			   ok=false;
			   bb=true;
			}
			mjd+=cmstep;
		 }
		 psr->dmoffsCMnum=i;
	  }
   }

   logmsg("psr: %s ncm=%d",psr->name,psr->dmoffsCMnum);



}

