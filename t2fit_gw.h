#include <tempo2.h>

double t2FitFunc_gwm_amp(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
double t2FitFunc_gwb_amp(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
double t2FitFunc_gwcs_amp(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
double t2FitFunc_quad_om(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
double t2FitFunc_gwsingle(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);


void t2UpdateFunc_gwsingle(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);
void t2UpdateFunc_quad_om(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) ;
