#include <tempo2.h>
#include "t2fit_position.h"
#include "t2fit_fitwaves.h"
#include "t2fit_glitch.h"
#include "t2fit_ifunc.h"
#include "t2fit_dmmodel.h"
#include "t2fit_dmother.h"

void t2UpdateFunc_simpleAdd(pulsar *psr, int ipsr ,param_label label,int k, double val, double error);
void t2UpdateFunc_simpleMinus(pulsar *psr, int ipsr ,param_label label,int k, double val, double error);

double t2FitFunc_zero(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_zero(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_stdFreq(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_stdFreq(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_binaryModels(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_binaryModels(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_planet(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_planet(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_stdDm(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
// dm is "simple"


double t2FitFunc_stdGravWav(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_stdGravWav(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);


double t2FitFunc_telPos(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_telPos(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);


double t2FitFunc_ifunc(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_ifunc(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);


double t2FitFunc_jump(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_jump(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_notImplemented(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_notImplemented(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

