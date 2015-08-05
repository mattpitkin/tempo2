#include "tempo2.h"
double t2FitFunc_dmmodelDM(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_dmmodelDM(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);
double t2FitFunc_dmmodelCM(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_dmmodelCM(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);
