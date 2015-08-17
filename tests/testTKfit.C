#include <gtest/gtest.h>
#include <src/gtest_main.cc>
#include <cfloat>
#include <cmath>
#include "TKfit.h"
#include "T2accel.h"
#define FIT_EPSILON 1.01*DBL_EPSILON


TEST(testTKpoly, trivial){
    const int m=2;
    double x[]={-1.0,0.0,1.0};
    double y[]={-1.0,0.0,1.0};
    double p[m];
    useT2accel=1;
    TKfindPoly_d(x,y,3,m,p);
    ASSERT_LT(fabs(p[0]),FIT_EPSILON);
    ASSERT_LT(fabs(p[1]-1.0),FIT_EPSILON);
}

TEST(testTKpoly, quadratic){
    const int m=3;
    const int n=1024;
    const double A = 1.0;
    const double B = 100.0;
    const double C = 1123.123;
    double x[n];
    double y[n];
    double Y[n];
    double p[m];
    double q[m];
    q[0]=q[1]=q[2]=0;
    long double d = 1.0L/static_cast<long double>(n);
    for (int i = 0; i < n; i++){
        long double ix = static_cast<long double>(i-n/2) * d;
        x[i] = static_cast<double>(ix);
        y[i] = A + ix*B + ix*ix*C;
        Y[i] = A + ix*B + ix*ix*C;
    }
    useT2accel=1;
    int itr=1;
    while(itr < 10){
        TKfindPoly_d(x,y,n,m,p);

        q[0] += p[0];
        q[1] += p[1];
        q[2] += p[2];
        if (
                (fabs(q[0] - A) < FIT_EPSILON) && 
                (fabs(q[1] - B) < FIT_EPSILON) && 
                (fabs(q[2] - C) < FIT_EPSILON)
           ) break;

        for (int i = 0; i < n; i++){
            long double ix = static_cast<long double>(i-n/2) * d;
            y[i] = Y[i] - (q[0] + ix*q[1] + ix*ix*q[2]);
        }
        itr++;

    }
    ASSERT_LT(itr,3);
}

