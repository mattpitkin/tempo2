#include <tempo2pred.h>
#include "t2pred_loops.h"

void T2Predictor_GetPhase_array_ld(const T2Predictor *t2p, long double *mjd, const int nmjd, long double freq, double* out) {
    int i;
    for (i=0; i < nmjd; ++i){
        out[i] = T2Predictor_GetPhase(t2p,mjd[i],freq);
    }
}
void T2Predictor_GetFrequency_array_ld(const T2Predictor *t2p, long double *mjd, const int nmjd, long double freq, double* out) {
    int i;
    for (i=0; i < nmjd; ++i){
        out[i] = T2Predictor_GetFrequency(t2p,mjd[i],freq);
    }
}



void T2Predictor_GetPhase_array(const T2Predictor *t2p, double refmjd, double *mjd, const int nmjd, double freq, double* out) {
    int i;
    for (i=0; i < nmjd; ++i){
        out[i] = T2Predictor_GetPhase(t2p,mjd[i],freq);
    }
}
void T2Predictor_GetFrequency_array(const T2Predictor *t2p, double refmjd, double *mjd, const int nmjd, double freq, double* out) {
    int i;
    for (i=0; i < nmjd; ++i){
        out[i] = T2Predictor_GetFrequency(t2p,mjd[i],freq);
    }
}


