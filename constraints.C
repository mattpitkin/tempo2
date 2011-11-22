#include "constraints.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

std::string get_constraint_name(enum constraint c){
	switch(c){
		case constraint_dmmodel_mean:
			return "DMMODEL mean(DM) = 0";
		case constraint_dmmodel_cw_0:
			return "DMMODEL mean(C) = 0";
		case constraint_dmmodel_cw_1:
			return "DMMODEL linear(C) = 0";
		case constraint_dmmodel_cw_2:
			return "DMMODEL quadratic(C) = 0";
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


		default:
			return "UNKNOWN!";
	}
}

double consFunc_dmmodel_mean(pulsar *psr,int i,int k){
	/*
	 * Only operate on param=dmmodel and when fit parameter is 
	 * one of the frequency dependant parts (i.e. first dmoffsNum)
	 */

	if(i==param_dmmodel && k < psr->dmoffsNum){
		return 1.0;
	} else return 0;
}
double consFunc_dmmodel_cw(pulsar *psr,int i,int k,int order){
	/*
	 * Only operate on param=dmmodel and when fit parameter is 
	 * one of the frequency independant parts (i.e. last dmoffsNum).
	 */
	if(i==param_dmmodel && k >= psr->dmoffsNum){
		long double epoch = psr->param[param_pepoch].val[0];
		return pow(psr->dmoffsMJD[k%psr->dmoffsNum]-epoch,order);
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
  printf("In here with %d %d %d %s\n",i,k,order,psr->name);
  if(i==param_ifunc){
    long double epoch = psr->param[param_pepoch].val[0];
    return pow(psr->ifuncT[k]-epoch,order);
  }
  else return 0;
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
