#include "constraints.h"
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


double consFunc_ifunc(pulsar *psr,int i,int k,int order){
	/*
	 * Only operate on param=ifunc and when fit parameter is 
	 * one of the frequency independant parts (i.e. last ifuncN).
	 */
	if(i==param_ifunc){
		long double epoch = psr->param[param_pepoch].val[0];
		return pow(psr->ifuncT[k]-epoch,order);
	} else return 0;

}



