#include "constraints.h"
#include "TKfit.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define CONSTRAINT_WEIGHTS

std::string get_constraint_name(enum constraint c){
#ifdef CONSTRAINT_WEIGHTS
		printf("[WT]");
#endif
	switch(c){
		case constraint_dmmodel_mean:
			return "DMMODEL mean(DM) = 0";
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

		default:
			return "UNKNOWN!";
	}
}



void matrixDMConstraintWeights(pulsar *psr,double** uinv){
		int i,j,k;
		int nobs=0;
		int nfit=psr->dmoffsDMnum+psr->dmoffsCMnum;
		int nDM=psr->dmoffsDMnum;

		double x;


		logdbg("Getting DM constraints for %s",psr->name);
		if(uinv!=NULL)logdbg("Using UINV for weighting");


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

		double* _designMatrix=(double*)malloc(sizeof(double)*nobs*nfit);
		double* _v=(double*)malloc(sizeof(double)*nfit*nfit);
		double* _u=(double*)malloc(sizeof(double)*nobs*nfit);
		double* _uout=(double*)malloc(sizeof(double)*nobs*nfit);

		double ** designMatrix=(double**)malloc(sizeof(double*)*nobs);
		double ** u=(double**)malloc(sizeof(double*)*nobs);
		double ** uout=(double**)malloc(sizeof(double*)*nobs);
		double ** v=(double**)malloc(sizeof(double*)*nfit);
		double  * w=(double*)malloc(sizeof(double)*nfit);
		for (k=0; k < nfit;k++){
				v[k]=_v+nfit*k;
		}
		for (k=0; k < nobs;k++){
		   designMatrix[k]=_designMatrix+nfit*k;
		   u[k]=_u+nfit*k;
		   uout[k]=_uout+nfit*k;
		}

		
		nobs=0;
		for(i=0; i < psr->nobs; i++){
				// Check for "ok" and start/finish coppied from dofit.
				// Needs to be the same so that the weighting is the same as the fit,
				// but in practice probably doesn't matter.
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
								x   = (double)(psr->obsn[i].bbat-psr->param[param_pepoch].val[0]);
								double sig;
								if(uinv==NULL) sig=psr->obsn[i].toaErr*1e-6;
								else sig=1;
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

		if(uinv!=NULL){
		   // if we are using a covariance matrix.
		   // do pre-whitening of data and model
			  TKmultMatrix(uinv,designMatrix,nobs,nfit,uout);
			  for (i=0;i<nobs;i++)
			  {
				 for (j=0;j<nfit;j++)
				 {
					designMatrix[i][j] = uout[i][j];
				 }
			  }
		}


		TKsingularValueDecomposition_lsq(designMatrix,nobs,nfit,v,w,u);
		double *wt = (double*)malloc(sizeof(double)*nfit);
		for (i=0;i<nfit;i++)
		{
		   if (w[i]!=0) wt[i] = 1.0/w[i]/w[i];
		   else wt[i] = 0.0;
		}
		double sum_wDM=0;
		double sum_wCM=0;
		for (i=0;i<nfit;i++)
		{
		   double sum=0.0;
		   for (k=0;k<nfit;k++)
			  sum+=v[i][k]*v[i][k]*wt[k];
		   if(i < nDM){
			  psr->dmoffsDM_weight[i]=1.0/(sum);
			  sum_wDM+=psr->dmoffsDM_weight[i];
		   } else {
			  psr->dmoffsCM_weight[i-nDM]=1.0/(sum);
			  sum_wCM+=psr->dmoffsCM_weight[i-nDM];
		   }



		}

		//normalise the weights
		for (i=0;i<psr->dmoffsDMnum;i++)
		   psr->dmoffsDM_weight[i]/=sum_wDM;
		for (i=0;i<psr->dmoffsCMnum;i++)
		   psr->dmoffsCM_weight[i]/=sum_wCM;



		// free everything .
		free(u);
		free(uout);
		free(designMatrix);
		free(v);
		free(w);
		free(wt);
		free(_u);
		free(_uout);
		free(_designMatrix);
		free(_v);

}


/*
 * Derive the weighting functions for the constraints, based upon the ToA errors.
 *
 */
void computeConstraintWeights(pulsar *psr, double** uinv){
   for (int k=0; k < psr->ifuncN; k++){
	  psr->ifunc_weights[k]=1.0/(double)psr->ifuncN;
   }
   for (int k=0; k < psr->dmoffsDMnum; k++)
	  psr->dmoffsDM_weight[k]=1.0/(double)psr->dmoffsDMnum;

   for (int k=0; k < psr->dmoffsCMnum; k++)
	  psr->dmoffsCM_weight[k]=1.0/(double)psr->dmoffsCMnum;

#ifdef CONSTRAINT_WEIGHTS
   if(psr->dmoffsDMnum>0 || psr->dmoffsCMnum>0) {
	  matrixDMConstraintWeights(psr,uinv);
   }
   /*
	* Derive weights for ifuncs
	*/
   if(psr->ifuncN>0) {
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
	  return psr->ifunc_weights[k]*pow(psr->ifuncT[k]-epoch,order);

   }
   else return 0;
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






void autosetDMCM(pulsar* psr, double dmstep,double cmstep, double start, double end, bool fixCMgrid){
   int i,j;
   double mjd;
   bool ok;

   
   double psrstart=(double)(psr->obsn[0].sat)-dmstep/2.0;
   double psrend=(double)(psr->obsn[psr->nobs-1].sat)+dmstep/2.0;

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


   ok=false;
   while (!ok){
	  matrixDMConstraintWeights(psr,NULL);
	  ok=true;
	  double threshold = 0.1/(double)psr->dmoffsDMnum;
	  i=0;
	  j=0;
	  mjd=start;
	  for (j=0; j < psr->dmoffsDMnum; j++){
		 psr->dmoffsDM_mjd[i]=psr->dmoffsDM_mjd[j];
		 psr->dmoffsDM[i]=0;
		 if(!ok || psr->dmoffsDM_weight[j]>threshold)i++;
		 else{
			logmsg("Skip %lf %lg",psr->dmoffsDM_mjd[j],psr->dmoffsDM_weight[i]);
			ok=false;
		 }
		 mjd+=dmstep;
	  }
	  psr->dmoffsDMnum=i;
   }

   logmsg("psr: %s ndm=%d",psr->name,psr->dmoffsDMnum);
   ok=fixCMgrid; // only do the CM if we are not fixing the grid
   while (!ok){
	  matrixDMConstraintWeights(psr,NULL);
	  ok=true;
	  double threshold = 0.05/(double)psr->dmoffsCMnum;
	  i=0;
	  j=0;
	  mjd=start;
	  bool bb=false;
	  for (j=0; j < psr->dmoffsCMnum; j++){
		 psr->dmoffsCM_mjd[i]=psr->dmoffsCM_mjd[j];
		 psr->dmoffsCM[i]=0;
		 if(psr->dmoffsCM_weight[j]>threshold){
			i++;
			bb=false;
		 } else{
			if (bb) i++; // don't delete multiple points in a row
			else logmsg("Skip %lf %lg",psr->dmoffsCM_mjd[j],psr->dmoffsCM_weight[i]);
			ok=false;
			bb=true;
		 }
		 mjd+=cmstep;
	  }
	  psr->dmoffsCMnum=i;
   }
   logmsg("psr: %s ncm=%d",psr->name,psr->dmoffsCMnum);
}

