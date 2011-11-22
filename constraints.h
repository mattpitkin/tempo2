#include <string.h>
#include "tempo2.h"

std::string get_constraint_name(enum constraint c);

/*
 * Constraints must be specified in the constraint enum in tempo2.h
 *
 */

double consFunc_dmmodel_mean(pulsar *psr,int i,int k);
double consFunc_dmmodel_cw(pulsar *psr,int i,int k, int order);
double consFunc_dmmodel_cw_year(pulsar *psr,int i,int k);


double consFunc_ifunc(pulsar *psr,int i,int k, int order);
double consFunc_tel_dx(pulsar *psr,int i,int k, int order);
double consFunc_tel_dy(pulsar *psr,int i,int k, int order);
double consFunc_tel_dz(pulsar *psr,int i,int k, int order);

double consFunc_quad_ifunc_p(pulsar *psr,int i,int k, int order);
double consFunc_quad_ifunc_c(pulsar *psr,int i,int k, int order);

