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


//double getParamDeriv(pulsar *psr,int ipos,double x,int i,int k);

void matrixDMConstraintWeights(pulsar *psr){
		int i,k;
		int nobs=0;
		int nfit=psr->dmoffsNum*2;
		double x;
		double ** fitMatrix=(double**)malloc(sizeof(double*)*psr->nobs);
		double ** u=(double**)malloc(sizeof(double*)*psr->nobs);
		double ** v=(double**)malloc(sizeof(double*)*nfit);
		double  * w=(double*)malloc(sizeof(double)*nfit);
		for (k=0; k < nfit;k++){
				v[k]=(double*)malloc(sizeof(double)*nfit);
		}
		for(i=0; i < psr->nobs; i++){
				// Check for "ok" and start/finish coppied from dofit.
				// Needs to be the same so that the weighting is the same as the fit,
				// but in practice probably doesn't matter.
				if (psr->obsn[i].deleted==0)
				{
						fitMatrix[nobs]=(double*)malloc(sizeof(double)*nfit);
						u[nobs]=(double*)malloc(sizeof(double)*nfit);
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
								double sig=psr->obsn[i].toaErr*1e-6;
								double dmf = 1.0/(DM_CONST*powl(psr->obsn[i].freqSSB/1.0e6,2));
								for (k=0; k < nfit;k++){
										fitMatrix[nobs][k]=getParamDeriv(psr,i,x,param_dmmodel,k)/sig;
										if(k<psr->dmoffsNum)fitMatrix[nobs][k]*=dmf;
								}
								nobs++;

						}
				}
		}
		TKsingularValueDecomposition_lsq(fitMatrix,nobs,nfit,v,w,u);
		double *wt = (double*)malloc(sizeof(double)*nfit);
		double *e = (double*)malloc(sizeof(double)*nfit);
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
				e[i] = sqrt(sum);
				psr->dmoffsCWeights[i]=1.0/(sum);
				if((int)(i/psr->dmoffsNum)==0)
						sum_wDM+=psr->dmoffsCWeights[i];
				else
						sum_wCM+=psr->dmoffsCWeights[i];

		}
		//normalise the weights
		sum_wDM/=(double)psr->dmoffsNum;
		sum_wCM/=(double)psr->dmoffsNum;
		for (i=0;i<nfit;i++){
				if((int)(i/psr->dmoffsNum)==0)
						psr->dmoffsCWeights[i]/=sum_wDM;
				else
						psr->dmoffsCWeights[i]/=sum_wCM;
		}


		// free everything .
		for(i=0; i < nobs; i++){
				free(fitMatrix[i]);
				free(u[i]);
		}
		for (k=0; k < nfit;k++)free(v[k]);
		free(u);
		free(fitMatrix);
		free(v);
		free(w);
		free(wt);
		free(e);

}


/*
 * Derive the weighting functions for the constraints, based upon the ToA errors.
 *
 */
void computeConstraintWeights(pulsar *psr, int npsr){
		for (int p=0; p<npsr; p++){
				for (int k=0; k < psr[p].ifuncN; k++){
						psr[p].ifunc_weights[k]=1.0;
				}
				for (int k=0; k < psr[p].dmoffsNum*2; k++){
						psr[p].dmoffsCWeights[k]=1.0;
				}

#ifdef CONSTRAINT_WEIGHTS
				if(psr[p].dmoffsNum>0) {
						matrixDMConstraintWeights(psr+p);
				}
				/*
				 * Derive weights for ifuncs
				 */
				if(psr[p].ifuncN>0) {
						for (int i=0; i < psr[p].nobs; i++){
								if (psr[p].obsn[i].deleted==0){
										// compute the weight that this ToA applies to each IFUNC
										for (int k=0;k<psr[p].ifuncN-1;k++)
										{
												if ((double)psr[p].obsn[i].sat >= psr[p].ifuncT[k] &&
																(double)psr[p].obsn[i].sat < psr[p].ifuncT[k+1])
												{
														double w=1.0/pow(psr[p].obsn[i].toaErr,2);
														double dt=(psr[p].ifuncT[k+1]-psr[p].ifuncT[k]);
														double t=(double)psr[p].obsn[i].sat-psr[p].ifuncT[k];
														psr[p].ifunc_weights[k+1]+=w*t/dt;
														psr[p].ifunc_weights[k]+=w*(1.0-t/dt);
														break;
												}
										}

								}
						}
						// Normalise the weights
						double sum=0;
						for (int k=0; k < psr[p].ifuncN; k++){
								sum+=psr[p].ifunc_weights[k];
						}
						sum/=(float)psr[p].ifuncN;
						for (int k=0; k < psr[p].ifuncN; k++){
								psr[p].ifunc_weights[k]/=sum;
						}
				}
#endif
		}
		return;
}

double consFunc_dmmodel_mean(pulsar *psr,int i,int k){
		/*
		 * Only operate on param=dmmodel and when fit parameter is 
		 * one of the frequency dependant parts (i.e. first dmoffsNum)
		 */

		if(i==param_dmmodel && k < psr->dmoffsNum){
				return psr->dmoffsCWeights[k];
		} else return 0;
}
double consFunc_dmmodel_cw(pulsar *psr,int i,int k,int order){
		/*
		 * Only operate on param=dmmodel and when fit parameter is 
		 * one of the frequency independant parts (i.e. last dmoffsNum).
		 */
		if(i==param_dmmodel && k >= psr->dmoffsNum){
				long double epoch = psr->param[param_pepoch].val[0];
				long double w=psr->dmoffsCWeights[k];
				return w*pow(psr->dmoffsMJD[k%psr->dmoffsNum]-epoch,order);
		} else return 0;

}

double consFunc_dmmodel_cw_year(pulsar *psr,int i,int k,int order){
		/*
		 * Only operate on param=dmmodel and when fit parameter is 
		 * one of the frequency independant parts (i.e. last dmoffsNum).
		 */
		if(i==param_dmmodel && k >= psr->dmoffsNum){
				long double epoch = psr->param[param_pepoch].val[0];
				long double t = psr->dmoffsMJD[k%psr->dmoffsNum]-epoch;
				long double x = 2.0*M_PI*t/365.25;
				long double w=psr->dmoffsCWeights[k];
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
