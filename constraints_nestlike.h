#include <tempo2.h>

double constraints_nestlike_red(pulsar *psr,int ipsr, int iconstraint,int iparam,int constraintk,int k,void* special);
double constraints_nestlike_red_dm(pulsar *psr,int ipsr, int iconstraint,int iparam,int constraintk,int k,void* special);

double constraints_nestlike_red_chrom(pulsar *psr,int ipsr, int iconstraint,int iparam,int constraintk,int k,void* special);

double constraints_nestlike_jitter(pulsar *psr,int ipsr, int iconstraint,int iparam,int constraintk,int k,void* special);

double constraints_nestlike_band(pulsar *psr,int ipsr, int iconstraint,int iparam,int constraintk,int k,void* special);

double constraints_nestlike_group(pulsar *psr,int ipsr, int iconstraint,int iparam,int constraintk,int k,void* special);
