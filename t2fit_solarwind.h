#include "tempo2.h"

// These definitely need refinement! I just pulled them from thin air. MJK
#define SOLAR_CYCLE_EPOCH 58818.0
#define SOLAR_CYCLE_PERIOD 4018.0


double t2FitFunc_ne_sw(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
double t2FitFunc_ne_sw_sin(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
double t2FitFunc_ne_sw_ifunc(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);



void t2UpdateFunc_ne_sw(pulsar *psr, int ipsr ,param_label label,int k, double val, double error);
void t2UpdateFunc_ne_sw_ifunc(pulsar *psr, int ipsr ,param_label label,int k, double val, double error);


double constraint_ne_sw_ifunc_function(pulsar *psr,int ipsr, int iconstraint,int iparam,int constraintk,int k,void* extra);


