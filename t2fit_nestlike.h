#include <tempo2.h>

double t2FitFunc_nestlike_red(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_nestlike_red(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_nestlike_red_dm(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_nestlike_red_dm(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);


double t2FitFunc_nestlike_red_chrom(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_nestlike_red_chrom(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);


double t2FitFunc_nestlike_jitter(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);



void t2UpdateFunc_nestlike_jitter(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);


double t2FitFunc_nestlike_band(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);

void t2UpdateFunc_nestlike_band(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);


double t2FitFunc_nestlike_group(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_nestlike_group(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);


double t2FitFunc_nestlike_shape_red(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
double t2FitFunc_nestlike_shape_dm(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);


void t2UpdateFunc_nestlike_shape(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);
