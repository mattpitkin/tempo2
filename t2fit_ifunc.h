#include "tempo2.h"
double ifunc(const double* mjd, const double t,const int N, const int k);
double sinfunc(const double *T, const double t, const int k);


double t2FitFunc_sifunc(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
double t2FitFunc_ifunc(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_ifunc(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);
