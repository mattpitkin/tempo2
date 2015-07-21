#include <tempo2.h>
#include "t2fit_positionFitFuncs.h"

double t2FitFunc_zero(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_zero(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_stdFreq(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_stdFreq(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_binaryModels(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_binaryModels(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_planet(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_planet(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_stdDm(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_stdDm(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_stdGlitch(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_stdGlitch(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_stdGravWav(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_stdGravWav(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);


double t2FitFunc_telPos(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_telPos(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_fitwaves(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_fitwaves(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_ifunc(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_ifunc(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_miscDm(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_miscDm(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);


