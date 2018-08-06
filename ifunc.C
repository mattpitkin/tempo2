#include <assert.h>
#include "tempo2.h"
#include "ifunc.h"

double ifunc(const double* mjd, const double t,const int N, const int k){
    double yoffs[MAX_IFUNC];
    logdbg("ifunc %lg => %lg : %lg %d %d",mjd[0],mjd[N-1],t,N,k);
    assert(k<N);

    // set up a mask that only has 1 for our chosen fit parameter.
    for (int ioff =0;ioff<N;ioff++){
        if (ioff==k){
            yoffs[ioff]=1;
        } else {
            yoffs[ioff]=0;
        }
    }

    // call the function that sums over all.
    return ifunc(mjd,yoffs,t,N);
}

double sinfunc(const double *T, const double bat, const int k){
    double dt = bat - T[k];
    double tt = M_PI/(T[1] - T[0])*dt;
    if (tt ==0) {
        return 1.0;
    } else {
        return sin(tt)/tt;
    }
}

double ifunc(const double* mjd, const double* yoffs, const double t,const int N){
    if (t < mjd[0]){
        // we are before the first jump
        // so our gradient is just the zeroth offset.
        return yoffs[0];
    } else if(t > mjd[N-1]){
        return yoffs[N-1];
    } else{
        // find the pair we are between...
        for (int ioff =0;ioff<N;ioff++){
            if(t >= mjd[ioff] && t < mjd[ioff+1]){
                double x1 = mjd[ioff];
                double x2 = mjd[ioff+1];
                double x = (t-x1)/(x2-x1);
                double y1=yoffs[ioff];
                double y2=yoffs[ioff+1];
                return (y2-y1)*x + y1;
            }
        }
    }
    logerr("Shouldn't get here");
    assert(0);
    return 0;
}

