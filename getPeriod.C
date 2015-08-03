#include <cstdio>
#include <cstdlib>
#include <tempo2pred.h>

int main (int argc, char** argv){
    T2Predictor pred;

    char* pred_fname = argv[1];
    double mjd=atof(argv[2]);

    T2Predictor_Read(&pred, pred_fname);

    const double psr_freq = (double)T2Predictor_GetFrequency(&pred,mjd,1400);

    printf("%lf %24.20lf %24.20lf\n",mjd,psr_freq,1.0/psr_freq);

    return 0;
}
