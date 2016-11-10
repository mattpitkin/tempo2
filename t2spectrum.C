#include "tempo2.h"
#include "TKspectrum.h"
#include "TKfit.h"
#include <math.h>

double GLOBAL_OMEGA;

// Spectral analysis using covariance matrix
// note: uinv array must start from 0, not 1
// NEW FEATURE:
// set nfit < 0 to automatically set it to nres/2-1
int calcSpectra_ri(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY_R,double *specY_I,int nfit,pulsar* psr) {

    //  pulsar *psr=NULL;
//    return calcSpectra_ri_T(uinv,resx,resy,nres,specX,specY_R,specY_I,nfit,(resx[nres-1]-resx[0]),'N',psr);
    return calcSpectraErr_complex(uinv,resx,resy,nres,specX,specY_R,specY_I,NULL,nfit);
}




int calcSpectra_ri_T(double **uinv,double *resx,double *resy,int nres,double *specX,double *specY_R,double *specY_I,int nfit,double T,char fitfuncMode, pulsar* psr) {
    int i,k;
// UNUSED VARIABLE //     int nSpec;

    if (nfit < 0)
        nfit=nres/2-1;

    double v[nfit];
    double chisq;
    int ip[nres];
    double param[3];

    // to allow for computing the spectrum of IFUNC/CM etc we need different fit functions.
    void (*FIT_FUNC)(double, double [], int,pulsar *,int,int); // the fit function we will use
    FIT_FUNC=fitMeanSineFunc; // this is the standard periodogram/
    if(fitfuncMode=='I'){
        if (psr==NULL){
            logerr("ERROR: cannot use IFUNC spectral estimation without a real pulsar\nCheck call to calcSpectra_ri_T()");
            return -1;
        }
        FIT_FUNC=fitMeanSineFunc_IFUNC; // this accounts for smoothing of the CM/IFUNC
    }
    double binfactor = (double)nres/(double)(nres-1);
    if(fitfuncMode=='T'){
        binfactor=1.0;
    }

    // Should fit independently to all frequencies
    for (i=0;i<nres;i++)
    {
        ip[i] = i;
    }
    logmsg("Computing %d spectral channels",nfit);
    for (k=0;k<nfit;k++)
    {
        GLOBAL_OMEGA = 2.0*M_PI/(T*binfactor)*(k+1);
        printf("In here k = %d\n",k);
        TKleastSquares_single_pulsar(resx,resy,nres,param,NULL,3,NULL,&chisq,FIT_FUNC,psr,1.0e-40,ip,1,uinv);
        printf("Out here k = %d\n",k);
        v[k] = (resx[nres-1]-resx[0])/365.25/2.0/pow(365.25*86400.0,2); 
        specX[k] = GLOBAL_OMEGA/2.0/M_PI;
        specY_R[k] = sqrt(v[k])*param[1];
        specY_I[k] = sqrt(v[k])*param[2];
    }
    return nfit;
}

// Fit for mean and sine and cosine terms at a specified frequency (G_OMEGA)
// The psr and ival parameters are ignored
void fitMeanSineFunc(double x,double *v,int nfit,pulsar *psr,int ival,int ipsr)
{
// UNUSED VARIABLE //     int i;
    v[0] = 1; // Fit for mean
    v[1] = cos(GLOBAL_OMEGA*x);
    v[2] = sin(GLOBAL_OMEGA*x);

}

/*
 * This only works for residuals faked from an IFUNC.
 * Note that ival must be in order.
 */
void fitMeanSineFunc_IFUNC(double x,double *v,int nfit,pulsar *psr,int ival,int ipsr)
{
    double m,c; // for the straight line.
    double x0,x1;
    double i0,i1;
    double w = GLOBAL_OMEGA;
    double A=0;

    //   logmsg("%lf %lf %d %d",x,w,ival,psr->nobs);
    v[0] = 0; // Fit for mean
    v[1]=0; // will be cos
    v[2]=0; // will be sin

    x0 = x;
    if (ival< psr->nobs-1){
        x1 = x0 + (psr->obsn[ival+1].sat-psr->obsn[ival].sat); // this is the next x-val, under the assumptions.
        m=-1/pow(x0-x1,2);
        c=x1/pow(x0-x1,2);

        //logmsg("R) x0=%lf x1=%lf m=%lg c=%lg",x0,x1,m,c);
        i0=m*x0*x0/2.0 + c*x0;
        i1=m*x1*x1/2.0 + c*x1;
        v[0]+=i1-i0;

        // integrate the right side of the triangle.
        i0=(m*cos(w*x0) - w*(c+m*x0)*sin(w*x0))/w/w; // integral of (m*x+c).cos(w*x)
        i1=(m*cos(w*x1) - w*(c+m*x1)*sin(w*x1))/w/w; // integral of (m*x+c).sin(w*x)
        v[1]+=i1-i0;

        i0=(m*sin(w*x0) - w*(c+m*x0)*cos(w*x0))/w/w;
        i1=(m*sin(w*x1) - w*(c+m*x1)*cos(w*x1))/w/w;
        v[2]+=i1-i0;
        A+=0.5; // area
    }
    if (ival > 0){
        x1 = x0 + (psr->obsn[ival-1].sat-psr->obsn[ival].sat); // this is the next x-val, under the assumptions.
        m=1/pow(x0-x1,2);
        c=-x1/pow(x0-x1,2);

        //	  logmsg("L) x0=%lf x1=%lf m=%lg c=%lg",x0,x1,m,c);
        i0=m*x0*x0/2.0 + c*x0;
        i1=m*x1*x1/2.0 + c*x1;
        v[0]+=i0-i1;

        // integrate the right side of the triangle.
        i0=(m*cos(w*x0) - w*(c+m*x0)*sin(w*x0))/w/w;
        i1=(m*cos(w*x1) - w*(c+m*x1)*sin(w*x1))/w/w;
        v[1]+=i0-i1;

        i0=(m*sin(w*x0) - w*(c+m*x0)*cos(w*x0))/w/w;
        i1=(m*sin(w*x1) - w*(c+m*x1)*cos(w*x1))/w/w;
        v[2]+=i0-i1;
        A+=0.5; // area
    }
    v[0]/=A; // normalise to 1 in case we only had half the triangle.
    v[1]/=A; // normalise to 1 in case we only had half the triangle.
    v[2]/=A;

    //  double old[3];
    //  fitMeanSineFunc(x, old,nfit,psr,ival);
    //   logmsg("%lg %lg --  %lg %lg -- %lg %lg\n",v[0],old[0],v[1],old[1],v[2],old[2]);
}


// Fit for mean and sine and cosine terms at a specified frequency (G_OMEGA)
// The psr and ival parameters are ignored
void fitCosSineFunc(double x,double *v,int nfit,pulsar *psr,int ival,int ipsr)
{
// UNUSED VARIABLE //     int i;
    v[0] = cos(GLOBAL_OMEGA*x);
    v[1] = sin(GLOBAL_OMEGA*x);

}
