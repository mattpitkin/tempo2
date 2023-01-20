#include <math.h>
#include "shapelet.h"
#include "TKlog.h"

/**
 * iter_factorial copied from TempoNestUtilities.cpp
 */
double iter_factorial(unsigned int n)
{
    double ret = 1;
    for(unsigned int i = 1; i <= n; ++i)
        ret *= i;
    return ret;
}

/**
 * Copied from TempoNestLikeFuncs.cpp
 *
 * With some minor cleanup
 */
void TNothpl(int n,double x,double *pl){
    double a=2.0;
    double b=0.0;
    double c=1.0;
    double y0=1.0;
    double y1=2.0*x;
    pl[0]=1.0;
    pl[1]=2.0*x;
    for(int k=2;k<n;k++){
        double c=2.0*(k-1.0);
        y0=y0/sqrt(double(k*1.0));
        y1=y1/sqrt(double(k*1.0));
        double yn=(a*x+b)*y1-c*y0;
        yn=yn;
        pl[k]=yn;
        y0=y1;
        y1=yn;
    }
}

/*
 * Avoid a load of arrays and such - just do it iteratively.
 * Should return same values as above in sequence.
 */
double iterative_TNothpl(int n, double x, double &y0, double &y1) {
    if (n==0) {
        y0 = 1.0;
        return y0;
    }
    else if (n==1) {
        y1 = 2.0*x;
        return y1;
    } else {
        double a=2.0;
        double b=0.0;
        double sqrt_k = 1.0/sqrt((double)n);
        double c=2.0*(n-1.0);
        y0 *= sqrt_k;
        y1 *= sqrt_k;
        double yn = (a*x+b)*y1-c*y0;
        y0 = y1;
        y1 = yn;
        return yn;
    }
}


double evaluateShapelet(int ncoeff, double pos, double width, double* coef, double bat) {

    double result = 0;
    double H = (bat-pos)/(sqrt(2.0)*width);
    double y0,y1;


    for(int c=0; c < ncoeff; ++c){
        double shapeNorm=1.0/sqrt(sqrt(2.0*M_PI)*pow(2.0,c)*iter_factorial(c));
        result += shapeNorm * coef[c] * iterative_TNothpl(c, H, y0,y1);
    }
    return result * exp(-0.5*pow((bat-pos)/width, 2));
}
