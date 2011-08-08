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



