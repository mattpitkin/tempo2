#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

/*
 *    This file is part of TEMPO2. 
 * 
 *    TEMPO2 is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    TEMPO2 is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
 */

/*
 *    If you use TEMPO2 then please acknowledge it by citing 
 *    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
 *    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
 *    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
 *    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
 *    timing model.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <string.h>
#include <dlfcn.h>
#include "tempo2.h"
#include "TKfit.h"
#include "t2fit.h"
#include "TKsvd.h"
#include "T2accel.h"
void globalFITfuncs(double x,double afunc[],int ma,pulsar *psr,int ipos,int ipsr);

void updateGlobalParameters(pulsar* psr,int npsr, double* val,double* error);

int getNparams(pulsar *psr,int offset);
int getNglobal(pulsar *psr,int npsr);
double getConstraintDeriv(pulsar *psr,int ipos,int i,int k);

void FITfuncs(double x,double afunc[],int ma,pulsar *psr,int ipos,int ipsr);
void updateParameters(pulsar *psr,int p,double *val,double *error);


///////////////////Functions for TempoNest maximum likelihood Fitting////////////////

#ifdef HAVE_LAPACK
#ifdef HAVE_BLAS
void getTempoNestMaxLike(pulsar *pulse, int npsr);

void dgesvd_ftoc(double *in, double **out, int rows, int cols);
double* dgesvd_ctof(double **in, int rows, int cols);
void dgesvd(double **A, int m, int n, double *S, double **U, double **VT);
void dgemv(double **A, double *vecin,double *vecout,int rowa, int cola, char AT);
double *dgemv_ctof(double **in, int rows, int cols);
void dgemv_ftoc(double *in, double **out, int rows, int cols);
void dgemm(double **A, double **B,double **C,int rowa, int cola, int rowb, int colb, char AT, char BT);
double *dgemm_ctof(double **in, int rows, int cols);
void dgemm_ftoc(double *in, double **out, int rows, int cols);
void dpotri(double **A, int msize);
double *dpotri_ctof(double **in, int rows, int cols);
void dpotri_ftoc(double *in, double **out, int rows, int cols);
void dpotrf(double **A, int msize, double &det);
double *dpotrf_ctof(double **in, int rows, int cols);
void dpotrf_ftoc(double *in, double **out, int rows, int cols);

extern "C" void dpotrf_(char *UPLO, int *msize, double *a, int *lda, int *info);
extern "C" void dpotri_(char *UPLO, int *msize, double *a, int *lda, int *info);
extern "C" void dgesvd_(char *jobu, char *jobvt, int *m, int *n,
        double *a, int *lda, double *s, double *u,
        int *ldu, double *vt, int *ldvt, double *work,
        int *lwork, int *info);
extern "C" void dgemv_(char *jobu, int *m, int *n,
        double *alpha, double *a, int *lda,
        double *x, int *incx, double *beta, double *y, int *incy);
extern "C" void dgemm_(char *jobu, char *jobvt, int *m, int *n,
        int *k, double *alpha, double *a, int *lda,
        double *b, int *ldb, double *beta, double *c,
        int *ldc);
#endif
#endif
/////////////////End Functions for TempoNest maximum likelihood Fitting/////////////////////


extern "C" char* pluginFitFunc(pulsar *psr,int npsr,const char *covarFuncFile) {

    printf("\n\nCALLING getTempoNextMaxLike!\n\n");
#ifdef HAVE_LAPACK
#ifdef HAVE_BLAS
    getTempoNestMaxLike(psr, npsr);
    return 0;
#else
    printf("LAPACK required for use of temponest noise parameters\n");
    return 0;
#endif
#endif

}


int getNglobal(pulsar *psr,int npsr){
    int nGlobal=0;
    int i,k;
    // Add global parameters
    for (i=0;i<MAX_PARAMS;i++)
    {
        for (k=0;k<psr[0].param[i].aSize;k++)
        {
            if (psr[0].param[i].fitFlag[k]==2) {
                {
                    if (i!=param_wave_om && i!= param_ifunc && i!=param_quad_om &&
                            i!=param_tel_dx && i!= param_tel_dy && i!=param_tel_dz &&
                            i!=param_quad_ifunc_p && i!=param_quad_ifunc_c && i!=param_gwsingle && i!=param_wave_dm)
                    {
                        psr->fitParamI[nGlobal]  = i;
                        psr->fitParamK[nGlobal]  = k;

                        nGlobal++;
                    }
                }
            }
        }
    }
    /* Add extra parameters for sinusoidal whitening */
    if (psr[0].param[param_wave_om].fitFlag[0]==2)
    {
        for (i=0;i<psr->nWhite*2;i++)
        {psr->fitParamI[nGlobal+i]  = param_wave_om; psr->fitParamK[nGlobal+i]  = i;}
        nGlobal+=psr[0].nWhite*2;
    }
    if (psr[0].param[param_wave_dm].fitFlag[0]==2)
    {
        for (i=0;i<psr->nWhite_dm*2;i++)
        {psr->fitParamI[nGlobal+i]  = param_wave_dm; psr->fitParamK[nGlobal+i]  = i;}
        nGlobal+=psr[0].nWhite_dm*2;
    }





    if (psr[0].param[param_ifunc].fitFlag[0]==2){
        if (psr[0].param[param_ifunc].val[0] == 0)
        {
            for (i=0;i<psr->ifuncN-1;i++)
            {psr->fitParamI[nGlobal+i]  = param_ifunc; psr->fitParamK[nGlobal+i]  = i;}
            nGlobal+=psr[0].ifuncN-1;
        }
        else
        {
            for (i=0;i<psr->ifuncN;i++)
            {psr->fitParamI[nGlobal+i]  = param_ifunc; psr->fitParamK[nGlobal+i]  = i;}
            nGlobal+=psr[0].ifuncN;
        }
    }
    if (psr[0].param[param_quad_om].fitFlag[0]==2){
        for (i=0;i<psr->nQuad*4;i++)
        {psr->fitParamI[nGlobal+i]  = param_quad_om; psr->fitParamK[nGlobal+i]  = i;}
        nGlobal+=psr[0].nQuad*4;
    }

    if (psr[0].param[param_tel_dx].fitFlag[0]==2)
    {
        if (psr->param[param_tel_dx].val[0] < 2){
            for (i=0;i<psr->nTelDX;i++)
            {psr->fitParamI[nGlobal+i]  = param_tel_dx; psr->fitParamK[nGlobal+i]  = i;}
            nGlobal+=(psr[0].nTelDX);
        } else if (psr->param[param_tel_dx].val[0] == 2){
            for (i=0;i<psr->nTelDX-1;i++)
            {psr->fitParamI[nGlobal+i]  = param_tel_dx; psr->fitParamK[nGlobal+i]  = i;}
            nGlobal+=(psr[0].nTelDX-1);
        }
    }
    if (psr[0].param[param_tel_dy].fitFlag[0]==2)
    {     
        if (psr->param[param_tel_dy].val[0] < 2){
            for (i=0;i<psr->nTelDY;i++)
            {psr->fitParamI[nGlobal+i]  = param_tel_dy; psr->fitParamK[nGlobal+i]  = i;}
            nGlobal+=(psr[0].nTelDY);
        } else if (psr->param[param_tel_dy].val[0] < 2) {
            for (i=0;i<psr->nTelDY-1;i++)
            {psr->fitParamI[nGlobal+i]  = param_tel_dy; psr->fitParamK[nGlobal+i]  = i;}
            nGlobal+=(psr[0].nTelDY-1);
        }
    }
    if (psr[0].param[param_tel_dz].fitFlag[0]==2)
    {
        if (psr->param[param_tel_dz].val[0] < 2){
            for (i=0;i<psr->nTelDZ;i++)
            {psr->fitParamI[nGlobal+i]  = param_tel_dz; psr->fitParamK[nGlobal+i]  = i;}

            nGlobal+=(psr[0].nTelDZ);
        } else if (psr->param[param_tel_dz].val[0] == 2){
            for (i=0;i<psr->nTelDZ-1;i++)
            {psr->fitParamI[nGlobal+i]  = param_tel_dz; psr->fitParamK[nGlobal+i]  = i;}

            nGlobal+=(psr[0].nTelDZ-1);
        }
    }
    if (psr->param[param_quad_ifunc_p].fitFlag[0]==2)
    {
        for (i=0;i<psr->quad_ifuncN_p;i++)
        {psr->fitParamI[nGlobal+i]  = param_quad_ifunc_p; psr->fitParamK[nGlobal+i]  = i;}

        nGlobal+=(psr->quad_ifuncN_p);
    }
    if (psr->param[param_quad_ifunc_c].fitFlag[0]==2)
    {
        for (i=0;i<psr->quad_ifuncN_c;i++)
        {psr->fitParamI[nGlobal+i]  = param_quad_ifunc_c; psr->fitParamK[nGlobal+i]  = i;}
        nGlobal+=(psr->quad_ifuncN_c);
    }
    if (psr[0].param[param_gwsingle].fitFlag[0]==2)
    {
        for (i=0;i<4;i++)
        {psr->fitParamI[nGlobal+i]  = param_gwsingle; psr->fitParamK[nGlobal+i]  = i;}

        nGlobal+=(4);
    }
    return nGlobal;

}

int getNparams(pulsar *psr,int offset)
{
    int npol;
    int i,k;

    npol = 1;
    for (i=0;i<MAX_PARAMS;i++)
    {
        for (k=0;k<psr->param[i].aSize;k++)
        {
            if (psr->param[i].paramSet[k]==1 && psr->param[i].fitFlag[k]==1) {
                if (i!=param_start && i!=param_finish && i!=param_dmmodel && i!=param_gwsingle)
                {
                    psr->fitParamI[npol+offset]  = i;
                    psr->fitParamK[npol+offset]  = k;
                    npol++;
                }
            }
        }
    }
    /* Add extra parameters for jumps */
    for (i=1;i<=psr->nJumps;i++)
    {
        if (psr->fitJump[i]==1)
        {
            psr->fitParamI[npol+offset]  = -1;
            psr->fitParamK[npol+offset]  = 0;
            npol++;
        }
    }
    /* Add extra parameters for sinusoidal whitening */
    if (psr->param[param_wave_om].fitFlag[0]==1)
    {
        printf("waveScale at this point = %d\n",static_cast<int>(psr->waveScale));
        if (psr->waveScale==1)
        {
            for (i=0;i<psr->nWhite*2-1;i++)
            {psr->fitParamI[npol+i+offset]  = param_wave_om; psr->fitParamK[npol+i+offset]  = i;}
            npol+=psr->nWhite*2-1;

        }
        else if (psr->waveScale==2)
        {
            for (i=0;i<psr->nWhite*4-1;i++)
            {psr->fitParamI[npol+i+offset]  = param_wave_om; psr->fitParamK[npol+i+offset]  = i;}

            npol+=psr->nWhite*4-1;
        }
        else
        {
            for (i=0;i<psr->nWhite*2-1;i++)
            {psr->fitParamI[npol+i+offset]  = param_wave_om; psr->fitParamK[npol+i+offset]  = i;}
            npol+=psr->nWhite*2-1;      
        }
    }
    if (psr->param[param_wave_dm].fitFlag[0]==1)
    {

        for (i=0;i<psr->nWhite_dm*2-1;i++)
        {psr->fitParamI[npol+i+offset]  = param_wave_dm; psr->fitParamK[npol+i+offset]  = i;}
        npol+=psr->nWhite_dm*2-1;

    }


    if (psr->param[param_quad_om].fitFlag[0]==1)
    {
        for (i=0;i<psr->nQuad*4-1;i++)
        {psr->fitParamI[npol+i+offset]  = param_quad_om; psr->fitParamK[npol+i+offset]  = i;}

        npol+=(psr->nQuad*4)-1;
    }
    if (psr->param[param_ifunc].fitFlag[0]==1)
    {
        if (psr->param[param_ifunc].val[0] == 0)
        {
            for (i=0;i<psr->ifuncN-1;i++)
            {psr->fitParamI[npol+i+offset]  = param_ifunc; psr->fitParamK[npol+i+offset]  = i;}

            npol+=(psr->ifuncN-1);
        }
        else
        {
            for (i=0;i<psr->ifuncN-1;i++)
            {psr->fitParamI[npol+i+offset]  = param_ifunc; psr->fitParamK[npol+i+offset]  = i;}

            npol+=(psr->ifuncN-1);
        }
    }
    if (psr->param[param_clk_offs].fitFlag[0]==1)
    {
        for (i=0;i<psr->clkOffsN-1;i++)
        {psr->fitParamI[npol+i+offset]  = param_clk_offs; psr->fitParamK[npol+i+offset]  = i;}

        npol+=(psr->clkOffsN-1);
    }
    if (psr->param[param_tel_dx].fitFlag[0]==1 && psr->param[param_tel_dx].val[0] < 2)
    {
        for (i=0;i<psr->nTelDX-1;i++)
        {psr->fitParamI[npol+i+offset]  = param_tel_dx; psr->fitParamK[npol+i+offset]  = i;}

        npol+=(psr->nTelDX-1);
    }
    else if (psr->param[param_tel_dx].fitFlag[0]==1 && psr->param[param_tel_dx].val[0] == 2)
    {
        for (i=0;i<psr->nTelDX-2;i++)
        {psr->fitParamI[npol+i+offset]  = param_tel_dx; psr->fitParamK[npol+i+offset]  = i;}
        npol+=(psr->nTelDX-2);
    }
    if (psr->param[param_tel_dy].fitFlag[0]==1 && psr->param[param_tel_dy].val[0] < 2)
    {
        for (i=0;i<psr->nTelDY-1;i++)
        {psr->fitParamI[npol+i+offset]  = param_tel_dy; psr->fitParamK[npol+i+offset]  = i;}

        npol+=(psr->nTelDY-1);
    }
    else if (psr->param[param_tel_dy].fitFlag[0]==1 && psr->param[param_tel_dy].val[0] == 2)
    {
        for (i=0;i<psr->nTelDY-2;i++)
        {psr->fitParamI[npol+i+offset]  = param_tel_dy; psr->fitParamK[npol+i+offset]  = i;}


        npol+=(psr->nTelDY-2);
    }
    if (psr->param[param_tel_dz].fitFlag[0]==1 && psr->param[param_tel_dz].val[0] < 2)
    {
        for (i=0;i<psr->nTelDZ-1;i++)
        {psr->fitParamI[npol+i+offset]  = param_tel_dz; psr->fitParamK[npol+i+offset]  = i;}

        npol+=(psr->nTelDZ-1);
    }
    else if (psr->param[param_tel_dz].fitFlag[0]==1 && psr->param[param_tel_dz].val[0] == 2)
    {
        for (i=0;i<psr->nTelDZ-2;i++)
        {psr->fitParamI[npol+i+offset]  = param_tel_dz; psr->fitParamK[npol+i+offset]  = i;}

        npol+=(psr->nTelDZ-2);
    }
    if (psr->param[param_quad_ifunc_p].fitFlag[0]==1)
    {
        for (i=0;i<psr->quad_ifuncN_p-1;i++)
        {psr->fitParamI[npol+i+offset]  = param_quad_ifunc_p; psr->fitParamK[npol+i+offset]  = i;}

        npol+=(psr->quad_ifuncN_p-1);
    }
    if (psr->param[param_quad_ifunc_c].fitFlag[0]==1)
    {
        for (i=0;i<psr->quad_ifuncN_c-1;i++)
        {psr->fitParamI[npol+i+offset]  = param_quad_ifunc_c; psr->fitParamK[npol+i+offset]  = i;}

        npol+=(psr->quad_ifuncN_c-1);
    }
    /* Add extra parameters for DMMODEL fitting */
    if (psr->param[param_dmmodel].fitFlag[0]==1){
        for (i=0;i<psr->dmoffsDMnum;i++)
        {psr->fitParamI[npol+i+offset]  = -2; psr->fitParamK[npol+i+offset]  = i;}

        npol+=psr->dmoffsDMnum;
        for (i=0;i<psr->dmoffsCMnum;i++)
        {psr->fitParamI[npol+i+offset]  = -3; psr->fitParamK[npol+i+offset]  = i;}

        npol+=psr->dmoffsCMnum;
    }
    /* Add extra parameters for GW single source fitting */
    if (psr->param[param_gwsingle].fitFlag[0]==1)
    {
        for (i=0;i<4;i++)
        {psr->fitParamI[npol+i+offset]  = param_gwsingle; psr->fitParamK[npol+i+offset]  = i;}
        npol+=4; 
    }
    return npol;
}


void globalFITfuncs(double x,double afunc[],int ma,pulsar *psr,int ipos,int p){
    int i,j,k;
    int nglobal = psr[p].nGlobal;
    for (i=0;i<ma;i++) afunc[i]=0.0;

    // add global parameters.
    // Global fit
    int c=0;
    int kk;
    for (i=0;i<MAX_PARAMS;i++)
    {
        for (k=0;k<psr[p].param[i].aSize;k++)
        {
            if (psr[p].param[i].fitFlag[k] == 2)
            {
                //          afunc[c] = dotproduct(psr[p].posPulsar,psr[p].obsn[ipos].planet_ssb[4]);
                if (i==param_wave_om)
                {
                    for (kk=0;kk<2*psr[0].nWhite;kk++)
                    {
                        afunc[c]= getParamDeriv(&psr[p],ipos,x+(double)psr[p].param[param_pepoch].val[0] - (double)psr[0].param[param_waveepoch].val[0],i,kk);
                        c++;
                    }
                }

                else if (i==param_wave_dm)
                {
                    for (kk=0;kk<2*psr[0].nWhite_dm;kk++)
                    {
                        afunc[c]= getParamDeriv(&psr[p],ipos,x+(double)psr[p].param[param_pepoch].val[0] - (double)psr[0].param[param_waveepoch_dm].val[0],i,kk);
                        c++;
                    }
                }

                else if(i==param_ifunc){
                    if (psr[p].param[param_ifunc].val[0] == 0)
                    {
                        for (j=0;j<psr[p].ifuncN-1;j++)
                        {
                            afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
                            //                printf("ifc=%d %d %g\n",counter,c,afunc[c]);
                            c++;
                        }
                    }
                    else
                    {
                        for (j=0;j<psr[p].ifuncN;j++)
                        {
                            afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
                            //                printf("ifc=%d %d %g\n",counter,c,afunc[c]);
                            c++;
                        }
                    }

                }
                else if(i==param_quad_ifunc_p){
                    for (j=0;j<psr[p].quad_ifuncN_p;j++)
                    {
                        afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
                        c++;
                    }

                }
                else if(i==param_quad_ifunc_c){
                    for (j=0;j<psr[p].quad_ifuncN_c;j++)
                    {
                        afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
                        c++;
                    }

                }
                else if(i==param_quad_om){
                    for (j=0;j<psr[p].nQuad*4;j++)
                    {
                        afunc[c] = getParamDeriv(&psr[p],ipos,x+(double)psr[p].param[param_pepoch].val[0],i,j);
                        //                printf("ifc=%d %d %g\n",counter,c,afunc[c]);
                        c++;
                    }
                }
                else if(i==param_tel_dx){
                    if (psr[p].param[i].val[0]<2){
                        for (j=0;j<psr[p].nTelDX;j++)
                        {
                            afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
                            c++;
                        }
                    } else {
                        for (j=0;j<psr[p].nTelDX-1;j++)
                        {
                            afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
                            c++;
                        }
                    }
                }
                else if(i==param_tel_dy){
                    if (psr[p].param[i].val[0]<2){
                        for (j=0;j<psr[p].nTelDY;j++)
                        {
                            afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
                            c++;
                        }
                    } else {
                        for (j=0;j<psr[p].nTelDY-1;j++)
                        {
                            afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
                            c++;
                        }
                    }
                }
                else if(i==param_tel_dz){
                    if (psr[p].param[i].val[0]<2){
                        for (j=0;j<psr[p].nTelDZ;j++)
                        {
                            afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
                            c++;
                        }
                    }else{
                        for (j=0;j<psr[p].nTelDZ-1;j++)
                        {
                            afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
                            c++;
                        }
                    }
                }
                else if (i==param_gwsingle)
                {
                    afunc[c] = getParamDeriv(&psr[p],ipos,x,i,0);
                    c++;
                    afunc[c] = getParamDeriv(&psr[p],ipos,x,i,1);
                    c++;
                    afunc[c] = getParamDeriv(&psr[p],ipos,x,i,2);
                    c++;
                    afunc[c] = getParamDeriv(&psr[p],ipos,x,i,3);
                    c++;
                }
                else
                {
                    afunc[c] = getParamDeriv(&psr[p],ipos,x,i,k);
                    c++;
                }
            }
        }

    }

    FITfuncs(x,afunc+nglobal,ma-nglobal,psr,ipos,p); // the non-global parameters.
}


void FITfuncs(double x,double afunc[],int ma,pulsar *psr,int ipos,int ipsr)
{

    int i,n=0,k,j,l,found;
    psr+=ipsr; // this avoids having to change the code now we have unified with global fitting.

    /*
     * NEW NOTICE: As of July 2015 the constraints have moved out of this function into a 
     * constraints matrix! You can now continue to assume that ipos is always less than psr->nobs
     *
     * *** Old notice below can be disregarded ***
     * OLD NOTICE: To allow for constraints, ipos may be larger than psr->nobs!
     * 
     * Therefore, if you modify this function, make sure to check (ipos < psr->nobs) before
     * using ipos, and call getParamDeriv to get the derivative since this function checks
     * for constraints and calls the appropriate constraint function/
     * 
     * M. Keith August 2011
     * *** end of old notice ***
     *
     */

    if(ipos < psr->nobs)afunc[n++] = 1;  /* Always fit for an arbitrary offset (unless this obs is a constraint!)*/
    else afunc[n++] = 0;
    /* See what we are fitting for */
    for (i=0;i<MAX_PARAMS;i++)
    {
        for (k=0;k<psr->param[i].aSize;k++)
        {
            if (psr->param[i].paramSet[k]==1 && psr->param[i].fitFlag[k]==1) /* If we are fitting for this parameter */
            {
                if (i!=param_start && i!=param_finish)
                {

                    logdbg("Fitting for %d (%s)",i,psr->param[i].label[k]);
                    if (i==param_wave_om)
                    {
                        if (psr->waveScale==2)
                        {
                            // 			  for (j=0;j<psr->nWhite*2;j++) 			  
                            for (j=0;j<psr->nWhite*4;j++)
                                afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
                        }
                        else
                        {
                            for (j=0;j<psr->nWhite*2;j++)
                                afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
                        }
                    }
                    else if (i==param_wave_dm)
                    {


                        fprintf(stderr, "here\n");
                        exit(0);
                        for (j=0;j<psr->nWhite_dm*2;j++)
                            afunc[n++] = getParamDeriv(psr,ipos,x,i,j);

                    }
                    else if (i==param_quad_om)
                    {
                        for (j=0;j<psr->nQuad*4;j++)
                            afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
                    }
                    else if (i==param_ifunc)
                    {
                        if (psr->param[param_ifunc].val[0] == 0)
                        {
                            for (j=0;j<psr->ifuncN-1;j++)
                                afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
                        }
                        else
                        {
                            for (j=0;j<psr->ifuncN;j++)
                                afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
                        }
                    }
                    else if (i==param_ifunc)
                    {
                        for (j=0;j<psr->clkOffsN;j++)
                            afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
                    }
                    else if (i==param_tel_dx)
                    {
                        if (psr->param[param_tel_dx].val[0]<2)
                        {
                            for (j=0;j<psr->nTelDX;j++)
                                afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
                        }
                        else if (psr->param[param_tel_dx].val[0]==2)
                        {
                            for (j=0;j<psr->nTelDX-1;j++)
                                afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
                        }
                    }
                    else if (i==param_tel_dy)
                    {
                        if (psr->param[param_tel_dy].val[0]<2)
                        {
                            for (j=0;j<psr->nTelDY;j++)
                                afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
                        }
                        else if (psr->param[param_tel_dy].val[0] == 2)
                        {
                            for (j=0;j<psr->nTelDY-1;j++)
                                afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
                        }
                    }
                    else if (i==param_tel_dz)
                    {
                        if (psr->param[param_tel_dz].val[0]<2)
                        {
                            for (j=0;j<psr->nTelDZ;j++)
                                afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
                        }
                        else if (psr->param[param_tel_dz].val[0] == 2)
                        {
                            for (j=0;j<psr->nTelDZ-1;j++)
                                afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
                        }
                    }
                    else if (i==param_quad_ifunc_p)
                    {
                        for (j=0;j<psr->quad_ifuncN_p;j++)
                            afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
                    }
                    else if (i==param_quad_ifunc_c)
                    {
                        for (j=0;j<psr->quad_ifuncN_c;j++)
                            afunc[n++] = getParamDeriv(psr,ipos,x,i,j);
                    }
                    else if (i==param_gwsingle)
                    {
                        afunc[n++] = getParamDeriv(psr,ipos,x,i,0);
                        afunc[n++] = getParamDeriv(psr,ipos,x,i,1);
                        afunc[n++] = getParamDeriv(psr,ipos,x,i,2);
                        afunc[n++] = getParamDeriv(psr,ipos,x,i,3);
                    }
                    else if (i==param_dmmodel)
                    {		      

                        double dmf = 1;
                        if (ipos < psr->nobs)dmf = 1.0/(DM_CONST*powl(psr->obsn[ipos].freqSSB/1.0e6,2));

                        for (j=0;j<(int)psr->dmoffsDMnum;j++)
                            afunc[n++] = dmf*getParamDeriv(psr,ipos,x,i,j);
                        for (j=0;j<(int)psr->dmoffsCMnum;j++)
                            afunc[n++] =     getParamDeriv(psr,ipos,x,i,j+psr->dmoffsDMnum);
                    }
                    else
                    {
                        afunc[n++] = getParamDeriv(psr,ipos,x,i,k);
                        //		      printf("getParamDeriv: n = %d, result = %g, %d %d %d %g %s\n",n-1,afunc[n-1],ipos,i,k,x,psr->param[i].shortlabel[k]);
                    }
                }
            }
        }
    }
    /* JUMPS */

    for (i=1;i<=psr->nJumps;i++)
    {
        if (psr->fitJump[i]==1)
        {
            found = 0;
            if(ipos < psr->nobs){
                for (l=0;l<psr->obsn[ipos].obsNjump;l++)
                {
                    if (psr->obsn[ipos].jump[l]==i)
                    {
                        found = 1;
                        break;
                    }
                    //	      else
                    //		found = 0.0;
                }
            }
            afunc[n++] = found;

        }
    } 
    if (n!=ma) { 
        printf("Problem in fitting routine n = %d, ma = %d\n",n,ma);
    }
}


void updateGlobalParameters(pulsar* psr,int npsr, double* val,double* error){
    int offset=0;
    int i,j,k,p;
    // update global parameters
    for (i=0;i<MAX_PARAMS;i++)
    {
        for (k=0;k<psr[0].param[i].aSize;k++)
        {
            if (psr[0].param[i].fitFlag[k] == 2)
            {
                //			printf("Have global parameter %d %d %d %g %g\n",i,param_wave_om,param_ifunc,val[offset],error[offset]);
                if (i==param_wave_om)
                {
                    int kk;
                    for (kk=0;kk<psr[0].nWhite;kk++)
                    {
                        for (p=0;p<npsr;p++)
                        {
                            psr[p].wave_cos[kk]  -= val[offset];
                            psr[p].wave_cos_err[kk] = error[offset];
                            if (p==0) printf("Have wave %d %d %g %g\n",offset,kk,val[offset],error[offset]);
                        }
                        offset++;
                        for (p=0;p<npsr;p++)
                        {
                            psr[p].wave_sine[kk] -= val[offset];
                            psr[p].wave_sine_err[kk] = error[offset];
                        }
                        offset++;
                    }
                    offset--;
                }
                else if (i==param_wave_dm) /* Whitening procedure using sinusoids */
                {
                    int k;
                    //exit(0);


                    for (k=0;k<psr[p].nWhite_dm;k++)
                    {
                        //fprintf(stderr, "%.3e\n", val[j]);
                        psr[p].wave_cos_dm[k]  -= val[j]; 
                        psr[p].wave_cos_dm_err[k] = error[j]; j++;
                        psr[p].wave_sine_dm[k] -= val[j]; 
                        psr[p].wave_sine_dm_err[k] = error[j]; j++;	      
                    }
                    j--;  

                    //exit(0);
                }




                else if (i==param_gwm_amp)
                {
                    //			  printf("In here with offset = %d\n",offset);
                    for (p=0;p<npsr;p++)
                    {			      
                        if (psr[p].param[param_gwm_amp].paramSet[1]==1)
                        {
                            psr[p].param[i].val[k] -= val[offset];
                            psr[p].param[i].err[k] = error[offset];
                        }
                        else
                        {
                            psr[p].param[i].val[0] -= val[offset];
                            psr[p].param[i].err[0] = error[offset];
                        }
                    }
                    //			  printf("Setting %d %g %g\n",offset,val[offset],error[offset]);
                    //			  offset++;

                }

                else if (i==param_gwb_amp)
                {
                    //			  printf("In here with offset = %d\n",offset);
                    for (p=0;p<npsr;p++)
                    {			      
                        if (psr[p].param[param_gwb_amp].paramSet[1]==1)
                        {
                            psr[p].param[i].val[k] -= val[offset];
                            psr[p].param[i].err[k] = error[offset];
                        }
                        else
                        {
                            psr[p].param[i].val[0] -= val[offset];
                            psr[p].param[i].err[0] = error[offset];
                        }
                    }
                    //			  printf("Setting %d %g %g\n",offset,val[offset],error[offset]);
                    //			  offset++;

                }


                else if(i==param_ifunc) {
                    printf("Updating %d point\n",psr[0].ifuncN);
                    if (psr[0].param[param_ifunc].val[0]==0)
                    {
                        for (j=0;j<psr[0].ifuncN-1;j++)
                        {
                            printf("Updating %d %g\n",offset,val[offset]);
                            for (p=0;p<npsr;p++)
                            {
                                psr[p].ifuncV[j]-=val[offset];
                                psr[p].ifuncE[j]=error[offset];
                            }
                            offset++;
                        }
                    }
                    else
                    {
                        for (j=0;j<psr[0].ifuncN;j++)
                        {
                            printf("Updating %d %g\n",offset,val[offset]);
                            for (p=0;p<npsr;p++)
                            {
                                psr[p].ifuncV[j]-=val[offset];
                                psr[p].ifuncE[j]=error[offset];
                            }
                            offset++;
                        }
                    }
                    offset--;
                }
                else if (i==param_quad_ifunc_p)
                {
                    for (j=0;j<psr[0].quad_ifuncN_p;j++)
                    {
                        printf("Updating %g\n",val[offset]);
                        for (p=0;p<npsr;p++)
                        {
                            psr[p].quad_ifuncV_p[j]-=val[offset];
                            psr[p].quad_ifuncE_p[j]=error[offset];
                        }
                        offset++;
                    }
                    offset--;
                }
                else if (i==param_quad_ifunc_c)
                {
                    for (j=0;j<psr[0].quad_ifuncN_c;j++)
                    {
                        printf("Updating %g\n",val[offset]);
                        for (p=0;p<npsr;p++)
                        {
                            psr[p].quad_ifuncV_c[j]-=val[offset];
                            psr[p].quad_ifuncE_c[j]=error[offset];
                        }
                        offset++;
                    }
                    offset--;
                }
                else if(i==param_quad_om) {
                    printf("Updating %d point\n",psr[0].ifuncN);
                    for (j=0;j<psr[0].nQuad;j++)
                    {
                        printf("Updating %d %g\n",offset,val[offset]);
                        for (p=0;p<npsr;p++)
                        {
                            psr[p].quad_aplus_r[j]    -= val[offset];
                            psr[p].quad_aplus_i[j]    -= val[offset+1];
                            psr[p].quad_across_r[j]   -= val[offset+2];
                            psr[p].quad_across_i[j]   -= val[offset+3];
                            psr[p].quad_aplus_r_e[j]   = error[offset];
                            psr[p].quad_aplus_i_e[j]   = error[offset+1];
                            psr[p].quad_across_r_e[j]  = error[offset+2];
                            psr[p].quad_across_i_e[j]  = error[offset+3];
                        }
                        offset+=4;
                    }
                    offset--;
                }
                else if (i==param_tel_dx)
                {
                    if (psr[0].param[i].val[0]<2){
                        for (j=0;j<psr[0].nTelDX;j++)
                        {
                            printf("Setting a: %d %g\n",j,val[offset]);
                            for (p=0;p<npsr;p++)
                            {
                                psr[p].telDX_v[j]-=val[offset];
                                psr[p].telDX_e[j]=error[offset];
                            }
                            offset++;
                        }
                        offset--;
                    }  else {
                        for (j=0;j<psr[0].nTelDX-1;j++)
                        {
                            printf("Setting b: %d %g\n",j,val[offset]);
                            for (p=0;p<npsr;p++)
                            {
                                psr[p].telDX_v[j]-=val[offset];
                                psr[p].telDX_e[j]=error[offset];
                            }
                            offset++;
                        }
                        offset--;

                    }
                }
                else if (i==param_tel_dy)
                {
                    if (psr[0].param[i].val[0]<2){
                        for (j=0;j<psr[0].nTelDY;j++)
                        {
                            printf("Setting c: %d %g\n",j,val[offset]);
                            for (p=0;p<npsr;p++)
                            {
                                psr[p].telDY_v[j]-=val[offset];
                                psr[p].telDY_e[j]=error[offset];
                            }
                            offset++;
                        }
                        offset--;
                    } else {
                        for (j=0;j<psr[0].nTelDY-1;j++)
                        {
                            printf("Setting d: %d %g\n",j,val[offset]);
                            for (p=0;p<npsr;p++)
                            {
                                psr[p].telDY_v[j]-=val[offset];
                                psr[p].telDY_e[j]=error[offset];
                            }
                            offset++;
                        }
                        offset--;
                    }
                }
                else if (i==param_tel_dz)
                {
                    if (psr[0].param[i].val[0]<2){
                        for (j=0;j<psr[0].nTelDZ;j++)
                        {
                            printf("Setting e: %d %g\n",j,val[offset]);
                            for (p=0;p<npsr;p++)
                            {
                                psr[p].telDZ_v[j]-=val[offset];
                                psr[p].telDZ_e[j]=error[offset];
                            }
                            offset++;
                        }
                        offset--;
                    } else {
                        for (j=0;j<psr[0].nTelDZ-1;j++)
                        {
                            for (p=0;p<npsr;p++)
                            {
                                psr[p].telDZ_v[j]-=val[offset];
                                psr[p].telDZ_e[j]=error[offset];
                            }
                            offset++;
                        }
                        offset--;

                    }
                }
                else if (i==param_gwsingle)
                {
                    for (p=0;p<npsr;p++)
                    {
                        psr[p].gwsrc_aplus_r -= val[offset];
                        psr[p].gwsrc_across_r -= val[offset+1];
                        psr[p].gwsrc_aplus_r_e = error[offset];
                        psr[p].gwsrc_across_r_e = error[offset+1];
                        psr[p].gwsrc_aplus_i -= val[offset+2];
                        psr[p].gwsrc_across_i -= val[offset+3];
                        psr[p].gwsrc_aplus_i_e = error[offset+2];
                        psr[p].gwsrc_across_i_e = error[offset+3];
                    }
                    offset+=3;
                }
                else
                {
                    for (p=0;p<npsr;p++)
                    {
                        if (i==param_telx || i==param_tely || i==param_telz)
                            psr[p].param[i].val[k] -= val[offset];
                        else
                            psr[p].param[i].val[k] += val[offset];
                        psr[p].param[i].err[k] = error[offset];
                    }
                }
                offset++;
            }
        }
    }


}
void updateParameters(pulsar *psr,int p,double *val,double *error)
{
    int i,j,k;
    logdbg("Updating parameters");
    psr[p].offset = val[0];
    psr[p].offset_e = error[0];
    j=1;
    for (i=0;i<MAX_PARAMS;i++)
    {
        for (k=0;k<psr[p].param[i].aSize;k++)
        {
            if (psr[p].param[i].paramSet[k]==1 && psr[p].param[i].fitFlag[k]==1 && (i!=param_start && i!=param_finish))
            {
                if (i==param_f) 
                {
                    if (k==0)
                    {
                        psr[p].param[param_f].val[k] *= (1.0-val[j]/psr[p].param[param_f].val[0]);	
                        psr[p].param[param_f].err[k]  = error[j];
                    }
                    else 
                    {
                        longdouble scale;
                        scale=longdouble(1.0);
                        if (k==2)      scale=1.0e9L;
                        else if (k>2 && k<10)  scale=1.0e18L;
                        else if (k>9) scale=1.0e23L;

                        psr[p].param[param_f].val[k] = psr[p].param[param_f].val[k] - 
                            (psr[p].param[param_f].val[0]*(val[j]/pow(24.0*3600.0,k+1))/scale);
                        psr[p].param[param_f].err[k] = error[j]/(pow(24.0*3600.0,k+1))/scale*
                            psr[p].param[param_f].val[0];
                    }
                }
                else if (i==param_dm || i==param_px || i==param_fddc || i==param_fddi || i==param_dmassplanet || i==param_dmx || i==param_fd || i==param_dm_sin1yr || i==param_dm_cos1yr)
                {
                    psr[p].param[i].val[k] += val[j];
                    psr[p].param[i].err[k]  = error[j];
                    // The following lines break the -dmo tim-file option and have therefore been disabled. 
                    // As far as I know, they don't really have an effect anyway.
                    //                                               JPWV, 08.05.2014
                    // This is slow - should be a better approach
                    // if (i==param_dm){
                    //   psr[p].dmOffset+=val[j];
                    // }
                }
                else if (i==param_dshk)
                {
                    psr[p].param[i].val[k] += val[j];
                    psr[p].param[i].err[k]  = error[j];
                }
                else if (i==param_pmrv)
                {
                    psr[p].param[i].val[k] += 10.0*val[j]*360.0*60.0*60.0/(2.0*M_PI);
                    psr[p].param[i].err[k]  = 10.0*error[j]*360.0*60.0*60.0/(2.0*M_PI);
                }
                else if (i==param_glph) /* Glitch phase */
                {
                    psr[p].param[i].val[k] -= val[j];     
                    psr[p].param[i].err[k]  = error[j];   
                }
                else if (i==param_glf0d) /* Glitch */
                {
                    psr[p].param[i].val[k] -= val[j]; 
                    psr[p].param[i].err[k]  = error[j]; 
                }
                else if (i==param_gltd) /* Glitch time delay */
                {
                    psr[p].param[i].val[k] -= val[j]; 
                    psr[p].param[i].err[k]  = error[j]; 
                }
                else if (i==param_glf0) /* Glitch permanent pulse frequency increment */
                {
                    psr[p].param[i].val[k] -= val[j]; 
                    psr[p].param[i].err[k]  = error[j];                        
                }
                else if (i==param_glf1) /* Glitch permanent pulse frequency deriv. increment */
                {
                    psr[p].param[i].val[k] -= val[j]; //*psr[p].param[param_f].val[0]; 
                    psr[p].param[i].err[k]  = error[j];                        
                }
                else if (i==param_glf2) /* Glitch permanent pulse frequency second deriv. increment */
                {
                    psr[p].param[i].val[k] -= val[j]*1.0e-27; //*psr[p].param[param_f].val[0]; 
                    psr[p].param[i].err[k]  = error[j]*1.0e-27;                        
                }
                else if (i==param_telx || i==param_tely || i==param_telz)
                {
                    psr[p].param[i].val[k] -= val[j];
                    psr[p].param[i].err[k] = error[j];
                }	      
                else if (i==param_raj)
                {
                    char retstr[100];
                    psr[p].param[param_raj].val[k] += val[j];
                    psr[p].param[param_raj].err[k] = error[j];

                    /* Must obtain this in hms form */
                    turn_hms(psr[p].param[param_raj].val[k]/(2.0*M_PI), retstr);
                    strcpy(psr[p].rajStrPost,retstr);
                }
                else if (i==param_decj)
                {
                    char retstr[100];
                    psr[p].param[param_decj].val[k] += val[j];
                    psr[p].param[param_decj].err[k] = error[j];
                    /* Must obtain this in dms form */
                    turn_dms(psr[p].param[param_decj].val[k]/(2.0*M_PI), retstr);
                    strcpy(psr[p].decjStrPost,retstr);
                }
                else if (i==param_pmra) /* Return in radian/sec */
                {
                    psr[p].param[param_pmra].val[k] += val[j]*180.0/M_PI*60.0*60.0*
                        1000.0*SECDAY*365.25/24.0/3600.0*cos(psr[p].param[param_decj].val[0]);
                    psr[p].param[param_pmra].err[k] = error[j]*180.0/M_PI*60.0*60.0*
                        1000.0*SECDAY*365.25/24.0/3600.0*cos(psr[p].param[param_decj].val[0]);
                }
                else if (i==param_pmdec) /* Return in radian/sec */
                {
                    psr[p].param[param_pmdec].val[k] += val[j]*180.0/M_PI*60.0*60.0*1000.0*
                        SECDAY*365.25/24.0/3600.0;
                    psr[p].param[param_pmdec].err[k] = error[j]*180.0/M_PI*60.0*60.0*1000.0*
                        SECDAY*365.25/24.0/3600.0;
                }	       
                else if (i==param_wave_om) /* Whitening procedure using sinusoids */
                {

                    /*
                       int k;
                       for (k=0;k<psr[p].nWhite;k++)
                       {
                       psr[p].wave_cos[k]  -= val[j]; 
                       psr[p].wave_cos_err[k] = error[j]; j++;
                       psr[p].wave_sine[k] -= val[j]; 
                       psr[p].wave_sine_err[k] = error[j]; j++;	      
                       }
                       if (psr->waveScale==2) // Ignore the non-frequency derivative terms
                       {
                       for (k=0;k<psr[p].nWhite;k++)
                       {
                       printf("Ignoring cos %g %g\n",val[j],error[j]); j++;
                       printf("Ignoring sin %g %g\n",val[j],error[j]); j++;
                       }
                    //		      j+=psr->nWhite*2; 
                    }
                    j--;

*/



                    if (psr->waveScale != 2)
                    {

                        for (k=0;k<psr[p].nWhite;k++)
                        {
                            fprintf(stderr, "%d %.3e \n", j, val[j]);
                            psr[p].wave_cos[k]  -= val[j]; 
                            psr[p].wave_cos_err[k] = error[j]; j++;
                            psr[p].wave_sine[k] -= val[j]; 
                            psr[p].wave_sine_err[k] = error[j]; j++;	      
                        }
                    }
                    else if(psr->waveScale==2) // Ignore the non-frequency derivative terms
                    {


                        for (k=0;k<psr[p].nWhite;k++)
                        {
                            //printf("Ignoring cos %g %g\n",val[j],error[j]); j++;
                            //printf("Ignoring sin %g %g\n",val[j],error[j]); j++;
                            psr[p].wave_cos_dm[k]-= val[j];
                            psr[p].wave_cos_dm_err[k]= error[j]; j++;
                            psr[p].wave_sine_dm[k] -= val[j];
                            psr[p].wave_sine_dm_err[k]= error[j]; j++;
                        }
                        // now do timing noise terms
                        for (k=0;k<psr[p].nWhite;k++)
                        {
                            psr[p].wave_cos[k]  -= val[j]; 
                            psr[p].wave_cos_err[k] = error[j]; j++;
                            psr[p].wave_sine[k] -= val[j]; 
                            psr[p].wave_sine_err[k] = error[j]; j++;	      
                        }

                        //		      j+=psr->nWhite*2; 
                    }
                    j--;

                }	


                else if (i==param_wave_dm) /* Whitening procedure using sinusoids */
                {
                    int k;
                    //exit(0);


                    for (k=0;k<psr[p].nWhite_dm;k++)
                    {
                        //fprintf(stderr, "%.3e\n", val[j]);
                        psr[p].wave_cos_dm[k]  -= val[j]; 
                        psr[p].wave_cos_dm_err[k] = error[j]; j++;
                        psr[p].wave_sine_dm[k] -= val[j]; 
                        psr[p].wave_sine_dm_err[k] = error[j]; j++;	      
                    }
                    j--;  

                    //exit(0);
                }


                else if (i==param_ifunc) 
                {
                    int k;
                    for (k=0;k<psr[p].ifuncN;k++)
                    {
                        psr[p].ifuncV[k] -= val[j];
                        psr[p].ifuncE[k] = error[j];
                        j++;
                    }
                }		  
                else if (i==param_clk_offs) 
                {
                    int k;
                    for (k=0;k<psr[p].clkOffsN;k++)
                    {
                        psr[p].clk_offsV[k] -= val[j];
                        psr[p].clk_offsE[k] = error[j];
                        j++;
                    }
                }		  
                else if (i==param_tel_dx) 
                {
                    int k;
                    if (psr[p].param[param_tel_dx].val[0] == -1)
                    {
                        printf("Updating with %g %g\n",val[j],error[j]);
                        psr[p].telDX_v[0] -= val[j]; //*SPEED_LIGHT/1000.0;
                        psr[p].telDX_e[0] = error[j]; //*SPEED_LIGHT/1000.0;
                    }
                    else if (psr[p].param[param_tel_dx].val[0] < 2)
                    {
                        for (k=0;k<psr[p].nTelDX;k++)
                        {
                            psr[p].telDX_v[k] -= val[j];
                            psr[p].telDX_e[k] = error[j];
                            j++;
                        }
                    }
                    else
                    {
                        for (k=0;k<psr[p].nTelDX-1;k++)
                        {
                            psr[p].telDX_v[k] -= val[j];
                            psr[p].telDX_e[k] = error[j];
                            j++;
                        }
                    }
                }
                else if (i==param_tel_dy) 
                {
                    int k;
                    if (psr[p].param[param_tel_dy].val[0] == -1)
                    {
                        psr[p].telDY_v[0] -= val[j]; //*SPEED_LIGHT/1000.0;
                        psr[p].telDY_e[0] = error[j]; //*SPEED_LIGHT/1000.0;
                    }
                    else if (psr[p].param[param_tel_dy].val[0] < 2)
                    {
                        for (k=0;k<psr[p].nTelDY;k++)
                        {
                            psr[p].telDY_v[k] -= val[j];
                            psr[p].telDY_e[k] = error[j];
                            j++;
                        }
                    }
                    else
                    {
                        for (k=0;k<psr[p].nTelDY-1;k++)
                        {
                            psr[p].telDY_v[k] -= val[j];
                            psr[p].telDY_e[k] = error[j];
                            j++;
                        }
                    }


                }		  
                else if (i==param_tel_dz) 
                {
                    int k;
                    if (psr[p].param[param_tel_dz].val[0] == -1)
                    {
                        psr[p].telDZ_v[0] -= val[j]; //*SPEED_LIGHT/1000.0;
                        psr[p].telDZ_e[0] = error[j]; //*SPEED_LIGHT/1000.0;
                    }
                    else if (psr[p].param[param_tel_dz].val[0] < 2)
                    {
                        for (k=0;k<psr[p].nTelDZ;k++)
                        {
                            psr[p].telDZ_v[k] -= val[j];
                            psr[p].telDZ_e[k] = error[j];
                            j++;
                        }
                    }
                    else
                    {
                        for (k=0;k<psr[p].nTelDZ-1;k++)
                        {
                            psr[p].telDZ_v[k] -= val[j];
                            psr[p].telDZ_e[k] = error[j];
                            j++;
                        }
                    }
                }		  
                else if (i==param_quad_ifunc_p) 
                {
                    int k;
                    for (k=0;k<psr[p].quad_ifuncN_p;k++)
                    {
                        psr[p].quad_ifuncV_p[k] -= val[j];
                        psr[p].quad_ifuncE_p[k] = error[j];
                        j++;
                    }
                }		  
                else if (i==param_quad_ifunc_c) 
                {
                    int k;
                    for (k=0;k<psr[p].quad_ifuncN_c;k++)
                    {
                        psr[p].quad_ifuncV_c[k] -= val[j];
                        psr[p].quad_ifuncE_c[k] = error[j];
                        j++;
                    }
                }		  
                else if (i==param_gwsingle)
                {
                    printf("%d GW SINGLE -- in here %g %g %g %g\n",j,val[j],error[j],val[j+1],error[j+1]);
                    psr[p].gwsrc_aplus_r -= val[j];
                    j++;
                    psr[p].gwsrc_across_r -= val[j];
                    j++;
                    psr[p].gwsrc_aplus_i -= val[j];
                    j++;
                    psr[p].gwsrc_across_i -= val[j];
                    printf("Now have: %g %g\n",psr[p].gwsrc_aplus_r,psr[p].gwsrc_across_r);
                }
                else if (i==param_gwm_amp)
                {
                    printf("Should I be in here???\n");
                    if (psr[p].param[param_gwm_amp].paramSet[1]==1)
                    {
                        psr[p].param[i].val[k] -= val[j];
                        psr[p].param[i].err[k] = error[j];
                        j++;
                    }
                    else
                    {
                        psr[p].param[i].val[0] -= val[j];
                        psr[p].param[i].err[0] = error[j];
                        j++;
                    }
                }

                else if (i==param_gwb_amp)
                {
                    printf("Should I be in here???\n");
                    if (psr[p].param[param_gwb_amp].paramSet[1]==1)
                    {
                        psr[p].param[i].val[k] -= val[j];
                        psr[p].param[i].err[k] = error[j];
                        j++;
                    }
                    else
                    {
                        psr[p].param[i].val[0] -= val[j];
                        psr[p].param[i].err[0] = error[j];
                        j++;
                    }
                }
                else if (i==param_brake)
                {
                    psr[p].param[i].val[0] -= val[j];
                    psr[p].param[i].err[0] = error[j];

                }


                else if (i==param_dmmodel)
                {
                    for (int k=0;k<psr[p].dmoffsDMnum;k++)
                    {
                        psr[p].dmoffsDM[k] += val[j];
                        psr[p].dmoffsDM_error[k] = error[j];
                        j++;
                    }
                    for (int k=0;k<psr[p].dmoffsCMnum;k++){

                        psr[p].dmoffsCM[k] = val[j];
                        psr[p].dmoffsCM_error[k] =  error[j];
                        j++;
                    }
                    j--;
                }
                else if (i==param_start)
                {
                }
                else if (i==param_finish)
                {
                }
                else if (strcmp(psr[p].binaryModel,"BT")==0)
                    updateBT(&psr[p],val[j],error[j],i);
                else if (strcmp(psr[p].binaryModel,"BTJ")==0)
                    updateBTJ(&psr[p],val[j],error[j],i,k);
                else if (strcmp(psr[p].binaryModel,"BTX")==0)
                    updateBTX(&psr[p],val[j],error[j],i,k);
                else if (strcmp(psr[p].binaryModel,"ELL1")==0)
                    updateELL1(&psr[p],val[j],error[j],i,k);
                else if (strcmp(psr[p].binaryModel,"DD")==0)
                    updateDD(&psr[p],val[j],error[j],i);
                else if (strcmp(psr[p].binaryModel,"DDK")==0)
                    updateDDK(&psr[p],val[j],error[j],i);
                else if (strcmp(psr[p].binaryModel,"DDS")==0)
                    updateDDS(&psr[p],val[j],error[j],i);
                else if (strcmp(psr[p].binaryModel,"DDGR")==0)
                    updateDDGR(&psr[p],val[j],error[j],i);
                else if (strcmp(psr[p].binaryModel,"MSS")==0)
                    updateMSS(&psr[p],val[j],error[j],i);
                else if (strcmp(psr[p].binaryModel,"T2")==0)
                    updateT2(&psr[p],val[j],error[j],i,k);
                else if (strcmp(psr[p].binaryModel,"T2-PTA")==0)
                    updateT2_PTA(&psr[p],val[j],error[j],i,k);
                else if (strcmp(psr[p].binaryModel,"DDH")==0)
                    updateDDH(&psr[p],val[j],error[j],i);
                else if (strcmp(psr[p].binaryModel,"ELL1H")==0)
                    updateELL1H(&psr[p],val[j],error[j],i);

                j++; /* Increment position in fit list */

            }
        }
    }
    if (strcmp(psr[p].binaryModel,"DDGR")==0) 
        DDGRmodel(psr,0,0,-2);  /* Update GR parameters */	  

    logdbg("Updating jumps; nJumps = %d",psr[p].nJumps);
    /* Now check jumps */
    for (i=1;i<=psr[p].nJumps;i++)
    {
        logdbg("%d fitJump = %d",i,psr[p].fitJump[i]);
        if (psr[p].fitJump[i]==1)
        {
            logdbg("%d Jump changes",i);
            logdbg("value = %g",(double)val[j]);
            logdbg("error = %g",(double)error[j]);
            psr[p].jumpVal[i] += -val[j];
            psr[p].jumpValErr[i] = error[j];
            j++;
        }
        /*	      printf("Have jumps %g %g\n",(double)val[j],error[j][j]); */
    }
    logdbg("Complete updating parameters");
}


#ifdef HAVE_LAPACK
#ifdef HAVE_BLAS



void othpl(int n,double x,double *pl){


    double a=2.0;
    double b=0.0;
    double y0=1.0;
    double y1=2.0*x;
    pl[0]=1.0;
    pl[1]=2.0*x;



    for(int k=2;k<n;k++){

        double c=2.0*(k-1.0);
        double yn=(a*x+b)*y1-c*y0;
        pl[k]=yn;
        y0=y1;
        y1=yn;

    }



}




void getTempoNestMaxLike(pulsar *pulse, int npsr){

    int subDM=pulse->TNsubtractDM;
    int subRed=pulse->TNsubtractRed;

    pulse->TNsubtractDM=0;
    pulse->TNsubtractRed=0;


    formBatsAll(pulse,npsr);       /* Form Barycentric arrival times */
    formResiduals(pulse,npsr,1);       /* Form residuals */


    pulse->TNsubtractDM=subDM;
    pulse->TNsubtractRed=subRed;

    /////////////////////////////////////////////////////////////////////////////////////////////  
    /////////////////////////Form the Design Matrix////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////  

    int numtofit = 1;
    int numTime=0;
    int numJumps=0;
    for(int i=0; i<MAX_PARAMS; i++) {
        for(int k=0; k<pulse->param[i].aSize; k++) {
            if(pulse->param[i].fitFlag[k]==1) {
                if(i!=param_start && i!=param_finish) {
                    numtofit++;
                    numTime++;
                }
            } // if fitFlag
        } // for k
    } // for i

    /* Add extra parameters for jumps */
    for(int i=1; i<=pulse->nJumps; i++) {
        if(pulse->fitJump[i]==1){
            numtofit++;
            numJumps++;
        }
    }


    double **TNDM=new double*[pulse->nobs];
    for(int i=0;i<pulse->nobs;i++){
        TNDM[i]=new double[numtofit];
    }


    double pdParamDeriv[MAX_PARAMS];
    double *TNDMScale=new double[numtofit];
    for(int j=0; j<numtofit; j++) {
        TNDMScale[j]=0;
    }

    for(int i=0; i < pulse->nobs; i++) {
        FITfuncs(pulse[0].obsn[i].bat - pulse[0].param[param_pepoch].val[0], pdParamDeriv, numtofit, pulse, i,0);
        for(int j=0; j<numtofit; j++) {
            TNDM[i][j]=pdParamDeriv[j];
            //			printf("PD: %i %i %g \n", i, j, pdParamDeriv[j]);
            TNDMScale[j]+=TNDM[i][j]*TNDM[i][j];
        } 
    } 


    for(int j=0; j<numtofit; j++) {
        TNDMScale[j]=sqrt(TNDMScale[j]);
    }


    for(int i=0; i < pulse->nobs; i++) {
        for(int j=0; j<numtofit; j++) {
            TNDM[i][j] = TNDM[i][j]/TNDMScale[j];
        }
    }


    int useOrthogonal=pulse->useTNOrth;	


    double *S;
    double **U;
    double **VT;
    double **V;

    if(useOrthogonal==1){


        S = new double[numtofit];
        U = new double*[pulse->nobs];
        for(int k=0; k < pulse->nobs; k++){
            U[k] = new double[pulse->nobs];
        }
        VT = new double*[numtofit]; 
        for (int k=0; k<numtofit; k++) VT[k] = new double[numtofit];

        dgesvd(TNDM,pulse->nobs, numtofit, S, U, VT);


        V=new double*[numtofit];


        for(int i=0;i<numtofit;i++){
            V[i]=new double[numtofit];
        }

        for(int j=0;j < numtofit;j++){
            for(int k=0;k < numtofit;k++){
                V[j][k]=VT[k][j];
            }
        }

        for(int j=0;j<pulse->nobs;j++){
            for(int k=0;k < numtofit;k++){
                TNDM[j][k]=U[j][k];	
            }	
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////Noise Hyperparameters//////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////


    double start,end;
    int go=0;
    for (int i=0;i<pulse->nobs;i++)
    {
        if (pulse->obsn[i].deleted==0)
        {
            if (go==0)
            {
                go = 1;
                start = (double)pulse->obsn[i].bat;
                end  = start;
            }
            else
            {
                if (start > (double)pulse->obsn[i].bat)
                    start = (double)pulse->obsn[i].bat;
                if (end < (double)pulse->obsn[i].bat)
                    end = (double)pulse->obsn[i].bat;
            }
        }
    }

    double maxtspan=1*(end-start);
    double averageTSamp=2*maxtspan/pulse->nobs;

    double **DMEventInfo;


    double RedAmp=0;
    double RedIndex=0;
    double RedFLow=0;
    double RedCorner=0;
    double DMAmp=0;
    double DMIndex=0;
    int FitRedCoeff=0;
    int FitDMCoeff=0;
    int DMEventCoeff=0;
    int totCoeff=0;
    int DMEventQuadTerms=0;

    if(pulse->TNRedAmp != 0 && pulse->TNRedGam != 0){
        RedAmp=pulse->TNRedAmp;
        RedIndex=pulse->TNRedGam;
        FitRedCoeff=2*pulse->TNRedC;
        RedFLow=pow(10.0, pulse->TNRedFLow);
        RedCorner=pulse->TNRedCorner/maxtspan;

        for(int i = 0; i < 200; i++){
            pulse->TNRedCoeffs[i] = 0;
        }


        printf("\nIncluding Red noise with %i Frequencies, %g Log_10 Amplitude, %g Spectral Index\n", FitRedCoeff/2,RedAmp,RedIndex);
    }
    if(pulse->TNDMAmp != 0 && pulse->TNDMGam != 0){
        DMAmp=pulse->TNDMAmp;
        DMIndex=pulse->TNDMGam;
        FitDMCoeff=2*pulse->TNDMC;


        for(int i = 0; i < 200; i++){
            pulse->TNDMCoeffs[i] = 0;
        }



        printf("\nIncluding DM Variations with %i Frequencies, %g Log_10 Amplitude, %g Spectral Index\n", FitDMCoeff/2,DMAmp,DMIndex);
    }

    if(pulse->nDMEvents > 0){
        DMEventInfo=new double*[pulse->nDMEvents];
        for(int i=0; i < pulse->nDMEvents; i++){


            printf("\nIncluding DM Event %i : %g Start, %g  Length, %g Log_10 Amp, %g Spectral Index\n", i,pulse->TNDMEvStart[i], pulse->TNDMEvLength[i], pulse->TNDMEvAmp[i], pulse->TNDMEvGam[i]);

            DMEventInfo[i]=new double[4];
            DMEventInfo[i][0]=pulse->TNDMEvStart[i]; //Start time
            DMEventInfo[i][1]=pulse->TNDMEvLength[i]; //Stop Time
            DMEventInfo[i][2]=pow(10.0, pulse->TNDMEvAmp[i]); //Amplitude
            DMEventInfo[i][3]=pulse->TNDMEvGam[i]; //SpectralIndex
            DMEventCoeff+=2*int(DMEventInfo[i][1]/averageTSamp);

            if(pulse->TNDMEvOff[i]==1){printf(" Including DMEvent offset \n"); DMEventQuadTerms++;}
            if(pulse->TNDMEvLin[i]==1){printf(" Including DMEvent Linear Term \n");DMEventQuadTerms++;}
            if(pulse->TNDMEvQuad[i]==1){printf(" Including DMEvent Quadratic Term \n");DMEventQuadTerms++;}

        }
    }


    // count ECORR values
    if(pulse->nTNECORR > 0){
        for(int i=0; i<pulse->nTNECORR; i++){
            printf("\nIncluding ECORR value for backend %s: %g mus", \
                    pulse->TNECORRFlagVal[i], pulse->TNECORRVal[i]);

        }
    }


    // find number of epochs (default dt= 10 s)
    int *Processed = new int[pulse->nobs];

    // initialize processed flags
    for (int i=0;i<pulse->nobs;i++){
        Processed[i] = 1;
    }

    // make sure we only process the epochs with the chosen flags
    for (int i=0;i<pulse->nobs;i++){
        for (int j=0;j<pulse->obsn[i].nFlags;j++){
            for (int k=0;k<pulse->nTNECORR;k++){
                if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k])==0){
                    if (strcmp(pulse->obsn[i].flagVal[j],pulse->TNECORRFlagVal[k])==0){
                        Processed[i] = 0;
                    }
                }
            }
        }
    }

    double dt = 10.0 / SECDAY;
    double satmin;
    double satmax;
    int nepoch = 0;
    int in = 0;
    int allProcessed = 0;
    while (!allProcessed){
        for (int i=0;i<pulse->nobs;i++){
            if (Processed[i]==0){
                satmin = (double)pulse->obsn[i].bat - dt;
                satmax = (double)pulse->obsn[i].bat + dt;
                break;
            }
        }
        for (int i=0;i<pulse->nobs;i++){
            for (int j=0;j<pulse->obsn[i].nFlags;j++){
                for (int k=0;k<pulse->nTNECORR;k++){
                    if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k])==0){
                        if (strcmp(pulse->obsn[i].flagVal[j],pulse->TNECORRFlagVal[k])==0){
                            if ((double)pulse->obsn[i].bat > satmin && \
                                    (double)pulse->obsn[i].bat < satmax){
                                Processed[i] = 1;
                                in++;
                            }
                        }
                    }
                }
            }
        }
        if (in != 0){
            nepoch++;
            in = 0;
        }
        allProcessed = 1;
        for (int i=0;i<pulse->nobs;i++){
            if (Processed[i]==0){
                allProcessed = 0;
                break;
            }
        }
    }


    if (nepoch > 0){
        printf("\n\nUsing %d epochs for PSR %s\n\n", nepoch, pulse->name);
    }
    // Total coefficients in F-matrix (include Jitter matrix here if present)
    totCoeff += FitRedCoeff;
    totCoeff += FitDMCoeff;
    totCoeff += DMEventCoeff;
    totCoeff += nepoch;

    if(pulse->TNBandDMAmp != 0 && pulse->TNBandDMGam != 0){
        printf("Including Band DM Noise: Amp %g   Index %g \n", pulse->TNBandDMAmp, pulse->TNBandDMGam);
        totCoeff += 6*pulse->TNBandDMC;
    }

    for(int i =0; i < pulse->nTNBandNoise; i++){
        printf("Including Band Noise between %g and %g MHz: Amp %g   Index %g \n", pulse->TNBandNoiseLF[i], pulse->TNBandNoiseHF[i], pulse->TNBandNoiseAmp[i], pulse->TNBandNoiseGam[i]);
        totCoeff += 2*pulse->TNBandNoiseC[i];
    }

    for(int i =0; i < pulse->nTNGroupNoise; i++){
        printf("Including Group Noise : Amp %g   Index %g \n", pulse->TNGroupNoiseAmp[i], pulse->TNGroupNoiseGam[i]);
        totCoeff += 2*pulse->TNGroupNoiseC[i];
    }

    printf("\nTotal number of coefficients: %d\n",totCoeff);

    double **FMatrix=new double*[pulse->nobs];
    for(int i=0;i<pulse->nobs;i++){
        FMatrix[i]=new double[totCoeff];
        for(int j=0;j<totCoeff;j++){
            FMatrix[i][j] = 0;
        }
    }


    double *freqs = new double[totCoeff];

    double *DMVec=new double[pulse->nobs];
    double DMKappa = 2.410*pow(10.0,-16);
    int startpos=0;
    double freqdet=0;

    double *powercoeff=new double[totCoeff];
    for(int o=0;o<totCoeff; o++){
        powercoeff[o]=0;
    }

    double f1yr = 1.0/3.16e7;


    //////////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////Red Noise///////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    RedAmp=pow(10.0, RedAmp);

    for (int i=0; i<FitRedCoeff/2; i++){

        freqs[startpos+i]=RedFLow*((double)(i+1.0))/(maxtspan);
        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

        double rho = pow((1+(pow((1.0/365.25)/RedCorner,RedIndex/2))),2)*(RedAmp*RedAmp/12.0/(M_PI*M_PI))/pow((1+(pow(freqs[i]/RedCorner,RedIndex/2))),2)/(maxtspan*24*60*60)*pow(f1yr,-3.0);
        powercoeff[i]+= rho;
        powercoeff[i+FitRedCoeff/2]+= rho;
        //printf("T2 RedC: %i %g %g %g %g %g \n", i, rho, RedAmp, RedIndex, freqs[i], maxtspan);
    }


    for(int i=0;i<FitRedCoeff/2;i++){
        for(int k=0;k<pulse->nobs;k++){
            double time=(double)pulse->obsn[k].bat;
            FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);

        }
    }

    for(int i=0;i<FitRedCoeff/2;i++){
        for(int k=0;k<pulse->nobs;k++){
            double time=(double)pulse->obsn[k].bat;
            FMatrix[k][i+FitRedCoeff/2]=sin(2*M_PI*freqs[i]*time);
        }
    }


    //////////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////DM Variations//////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////


    startpos=FitRedCoeff;
    DMAmp=pow(10.0, DMAmp);

    for (int i=0; i<FitDMCoeff/2; i++){

        freqs[startpos+i]=((double)(i+1.0))/maxtspan;
        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];

        double rho = (DMAmp*DMAmp)*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-DMIndex))/(maxtspan*24*60*60);	
        powercoeff[startpos+i]+=rho;
        powercoeff[startpos+i+FitDMCoeff/2]+=rho;
    }



    for(int o=0;o<pulse->nobs; o++){
        DMVec[o]=1.0/(DMKappa*pow((double)pulse->obsn[o].freqSSB,2));
    }

    for(int i=0;i<FitDMCoeff/2;i++){
        for(int k=0;k<pulse->nobs;k++){
            double time=(double)pulse->obsn[k].bat;
            FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
        }
    }

    for(int i=0;i<FitDMCoeff/2;i++){
        for(int k=0;k<pulse->nobs;k++){
            double time=(double)pulse->obsn[k].bat;
            FMatrix[k][startpos+i+FitDMCoeff/2]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
        }
    }



    //////////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////DM Events//////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////	

    startpos+=FitDMCoeff;

    if(pulse->nDMEvents > 0){

        for(int i =0; i < pulse->nDMEvents; i++){

            double DMamp=DMEventInfo[i][2];
            double DMindex=DMEventInfo[i][3];

            double Tspan = DMEventInfo[i][1];
            double f1yr = 1.0/3.16e7;
            int DMEventnumCoeff=int(Tspan/averageTSamp);

            for (int c=0; c<DMEventnumCoeff; c++){
                freqs[startpos+c]=((double)(c+1))/Tspan;
                freqs[startpos+c+DMEventnumCoeff]=freqs[startpos+c];

                double rho = (DMamp*DMamp)*pow(f1yr,(-3)) * pow(freqs[startpos+c]*365.25,(-DMindex))/(maxtspan*24*60*60);
                powercoeff[startpos+c]+=rho;
                powercoeff[startpos+c+DMEventnumCoeff]+=rho;
                freqdet=freqdet+2*log(powercoeff[startpos+c]);
            }


            for(int c=0;c<DMEventnumCoeff;c++){
                for(int k=0;k<pulse->nobs;k++){
                    double time=(double)pulse->obsn[k].bat;
                    if(time < DMEventInfo[i][0]+Tspan && time > DMEventInfo[i][0]){
                        FMatrix[k][startpos+c]=cos(2*M_PI*freqs[startpos+c]*time)*DMVec[k];
                        FMatrix[k][startpos+c+DMEventnumCoeff]=sin(2*M_PI*freqs[startpos+c]*time)*DMVec[k];
                    }
                    else{
                        FMatrix[k][startpos+c]=0;
                        FMatrix[k][startpos+c+DMEventnumCoeff]=0;
                    }
                }
            }

            startpos+=2*DMEventnumCoeff;

        }
    }

    // shapelet events
    int ShapeEventTerms=0;
    if(pulse->nTNShapeletEvents > 0){
        for(int i =0; i < pulse->nTNShapeletEvents; i++){
            ShapeEventTerms += pulse->TNShapeletEvN[i];
        }	

    }
    //printf("Total event coeffs: %i \n", ShapeEventTerms);


    startpos += DMEventCoeff;

    //////////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////DM Band Noise//////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    if(pulse->TNBandDMAmp != 0 && pulse->TNBandDMGam != 0){

        double BandDMAmp=pow(10.0, pulse->TNBandDMAmp);



        double startfreq=0;
        double stopfreq=0;

        /////////////////////////////////50CM/////////////////////////////////////////////////////////////////

        startfreq = 0;
        stopfreq=1000;

        for (int i=0; i<pulse->TNBandDMC; i++){

            freqs[startpos+i]=((double)(i+1.0))/maxtspan;
            freqs[startpos+i+pulse->TNBandDMC]=freqs[startpos+i];

            double rho = (BandDMAmp*BandDMAmp)*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-pulse->TNBandDMGam))/(maxtspan*24*60*60);	
            powercoeff[startpos+i]+=rho;
            powercoeff[startpos+i+pulse->TNBandDMC]+=rho;
        }




        for(int i=0;i<pulse->TNBandDMC;i++){
            for(int k=0;k<pulse->nobs;k++){
                if(pulse->obsn[k].freq > startfreq && pulse->obsn[k].freq < stopfreq){
                    double time=(double)pulse->obsn[k].bat;
                    FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                }
                else{	
                    FMatrix[k][startpos+i]=0;
                }
            }
        }

        for(int i=0;i<pulse->TNBandDMC;i++){
            for(int k=0;k<pulse->nobs;k++){
                if(pulse->obsn[k].freq > startfreq && pulse->obsn[k].freq < stopfreq){		
                    double time=(double)pulse->obsn[k].bat;
                    FMatrix[k][startpos+i+pulse->TNBandDMC]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                }
                else{
                    FMatrix[k][startpos+i+pulse->TNBandDMC]=0;
                }
            }
        }


        startpos=startpos+2*pulse->TNBandDMC;

        /////////////////////////////////20CM/////////////////////////////////////////////////////////////////

        startfreq = 1000;
        stopfreq=1800;

        for (int i=0; i<pulse->TNBandDMC; i++){

            freqs[startpos+i]=((double)(i+1.0))/maxtspan;
            freqs[startpos+i+pulse->TNBandDMC]=freqs[startpos+i];

            double rho = (BandDMAmp*BandDMAmp)*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-pulse->TNBandDMGam))/(maxtspan*24*60*60);	
            powercoeff[startpos+i]+=rho;
            powercoeff[startpos+i+pulse->TNBandDMC]+=rho;
        }




        for(int i=0;i<pulse->TNBandDMC;i++){
            for(int k=0;k<pulse->nobs;k++){
                if(pulse->obsn[k].freq > startfreq && pulse->obsn[k].freq < stopfreq){
                    double time=(double)pulse->obsn[k].bat;
                    FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                }
                else{	
                    FMatrix[k][startpos+i]=0;
                }
            }
        }

        for(int i=0;i<pulse->TNBandDMC;i++){
            for(int k=0;k<pulse->nobs;k++){
                if(pulse->obsn[k].freq > startfreq && pulse->obsn[k].freq < stopfreq){		
                    double time=(double)pulse->obsn[k].bat;
                    FMatrix[k][startpos+i+pulse->TNBandDMC]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                }
                else{
                    FMatrix[k][startpos+i+pulse->TNBandDMC]=0;
                }
            }
        }


        startpos=startpos+2*pulse->TNBandDMC;

        /////////////////////////////////10CM/////////////////////////////////////////////////////////////////


        startfreq = 1800;
        stopfreq=10000;

        for (int i=0; i<pulse->TNBandDMC; i++){

            freqs[startpos+i]=((double)(i+1.0))/maxtspan;
            freqs[startpos+i+pulse->TNBandDMC]=freqs[startpos+i];

            double rho = (BandDMAmp*BandDMAmp)*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-pulse->TNBandDMGam))/(maxtspan*24*60*60);	
            powercoeff[startpos+i]+=rho;
            powercoeff[startpos+i+pulse->TNBandDMC]+=rho;
        }




        for(int i=0;i<pulse->TNBandDMC;i++){
            for(int k=0;k<pulse->nobs;k++){
                if(pulse->obsn[k].freq > startfreq && pulse->obsn[k].freq < stopfreq){
                    double time=(double)pulse->obsn[k].bat;
                    FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                }
                else{	
                    FMatrix[k][startpos+i]=0;
                }
            }
        }

        for(int i=0;i<pulse->TNBandDMC;i++){
            for(int k=0;k<pulse->nobs;k++){
                if(pulse->obsn[k].freq > startfreq && pulse->obsn[k].freq < stopfreq){		
                    double time=(double)pulse->obsn[k].bat;
                    FMatrix[k][startpos+i+pulse->TNBandDMC]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                }
                else{
                    FMatrix[k][startpos+i+pulse->TNBandDMC]=0;
                }
            }
        }


        startpos=startpos+2*pulse->TNBandDMC;
    }


    //////////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////Band Noise//////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    int totalBandNoiseCoeff=0;
    for(int g =0; g < pulse->nTNBandNoise; g++){

        double BandLF = pulse->TNBandNoiseLF[g];
        double BandHF = pulse->TNBandNoiseHF[g];
        double BandAmp=pow(10.0, pulse->TNBandNoiseAmp[g]);
        double BandSpec=pulse->TNBandNoiseGam[g];
        int BandC=pulse->TNBandNoiseC[g];

        totalBandNoiseCoeff+=2*BandC;



        for (int i=0; i<BandC; i++){

            freqs[startpos+i]=((double)(i+1.0))/maxtspan;
            freqs[startpos+i+BandC]=freqs[startpos+i];

            double rho = (BandAmp*BandAmp/12.0/(M_PI*M_PI))*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-BandSpec))/(maxtspan*24*60*60);	
            powercoeff[startpos+i]+=rho;
            powercoeff[startpos+i+BandC]+=rho;
        }




        for(int i=0;i<BandC;i++){
            for(int k=0; k<pulse->nobs; k++){
                if(pulse->obsn[k].freq > BandLF && pulse->obsn[k].freq < BandHF){
                    double time=(double)pulse->obsn[k].bat;
                    FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time);
                    FMatrix[k][startpos+i+BandC]=sin(2*M_PI*freqs[startpos+i]*time);

                }
                else{
                    FMatrix[k][startpos+i]=0;
                    FMatrix[k][startpos+i+BandC]=0;
                }

            }
        }




        startpos=startpos+2*BandC;
    }


    //////////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////Group Noise//////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    int totalGroupCoeff=0;
    for(int g =0; g < pulse->nTNGroupNoise; g++){

        double GroupAmp=pow(10.0, pulse->TNGroupNoiseAmp[g]);
        double GroupSpec=pulse->TNGroupNoiseGam[g];
        int GroupC=pulse->TNGroupNoiseC[g];

        totalGroupCoeff+=2*GroupC;



        for (int i=0; i<GroupC; i++){

            freqs[startpos+i]=((double)(i+1.0))/maxtspan;
            freqs[startpos+i+GroupC]=freqs[startpos+i];

            double rho = (GroupAmp*GroupAmp/12.0/(M_PI*M_PI))*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-GroupSpec))/(maxtspan*24*60*60);	
            powercoeff[startpos+i]+=rho;
            powercoeff[startpos+i+GroupC]+=rho;
        }




        for(int i=0;i<GroupC;i++){
            for(int k=0; k<pulse->nobs; k++){
                int set=0;
                for (int j=0; j < pulse->obsn[k].nFlags; j++){

                    //Check Group Noise Flag

                    if (strcmp(pulse->obsn[k].flagID[j],pulse->TNGroupNoiseFlagID[g])==0){
                        if (strcmp(pulse->obsn[k].flagVal[j],pulse->TNGroupNoiseFlagVal[g])==0){
                            double time=(double)pulse->obsn[k].bat;
                            FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time);
                            FMatrix[k][startpos+i+GroupC]=sin(2*M_PI*freqs[startpos+i]*time);
                            set = 1;

                        }
                        else{
                            if(set == 0){	
                                FMatrix[k][startpos+i]=0;
                                FMatrix[k][startpos+i+GroupC]=0;
                            }
                        }
                    }
                    else{
                        if(set == 0){
                            FMatrix[k][startpos+i]=0;
                            FMatrix[k][startpos+i+GroupC]=0;
                        }
                    }	

                }
            }
        }




        startpos=startpos+2*GroupC;
    }



    //////////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////U Jitter Matrix////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////


    //// construct jitter matrix //


    //printf("\nStarting Coefficients at %d for Jitter:\n",startpos);

    // initialize processed flags
    for (int i=0;i<pulse->nobs;i++){
        Processed[i] = 1;
    }

    // make sure we only process the epochs with the chosen flags
    for (int i=0;i<pulse->nobs;i++){
        for (int j=0;j<pulse->obsn[i].nFlags;j++){
            for (int k=0;k<pulse->nTNECORR;k++){
                if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k])==0){
                    if (strcmp(pulse->obsn[i].flagVal[j],pulse->TNECORRFlagVal[k])==0){
                        Processed[i] = 0;
                    }
                }
            }
        }
    }


    nepoch = 0;
    in = 0;
    allProcessed = 0;
    while (!allProcessed){
        for (int i=0;i<pulse->nobs;i++){
            if (Processed[i]==0){
                satmin = (double)pulse->obsn[i].bat - dt;
                satmax = (double)pulse->obsn[i].bat + dt;
                break;
            }
        }
        for (int i=0;i<pulse->nobs;i++){
            for (int j=0;j<pulse->obsn[i].nFlags;j++){
                for (int k=0;k<pulse->nTNECORR;k++){
                    if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k])==0){
                        if (strcmp(pulse->obsn[i].flagVal[j],pulse->TNECORRFlagVal[k])==0){
                            if ((double)pulse->obsn[i].bat > satmin && \
                                    (double)pulse->obsn[i].bat < satmax){
                                Processed[i] = 1;
                                FMatrix[i][startpos+nepoch] = 1.0;
                                powercoeff[startpos+nepoch] = pow(pulse->TNECORRVal[k], \
                                        2.0)*1e-12;
                                in++;
                            }
                            else{
                                FMatrix[i][startpos+nepoch] = 0.0;
                            }
                        }
                    }
                }
            }
        }
        if (in != 0){
            nepoch++;
            in = 0;
        }
        allProcessed = 1;
        for (int i=0;i<pulse->nobs;i++){
            if (Processed[i]==0){
                allProcessed = 0;
                break;
            }
        }
    }


    /*for (int i;i<totCoeff;i++){
      printf("powercoeff %d: %g\n",i,powercoeff[i]);
      }*/

    // get Phi matrix
    double **PPFM=new double*[totCoeff];
    for(int i=0;i<totCoeff;i++){
        PPFM[i]=new double[totCoeff];
        for(int j=0;j<totCoeff;j++){
            PPFM[i][j]=0;
        }
    }


    for(int c1=0; c1<totCoeff; c1++){
        //		printf("PPFM %i %g \n", c1, powercoeff[c1]);
        PPFM[c1][c1]=1.0/powercoeff[c1];
    }

    //////////////////////////////////////////////////////////////////////////////////////////  
    //////////////////////Form Total Matrix///////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    int totalsize=numtofit+totCoeff+DMEventQuadTerms+ShapeEventTerms;
    int totalnoisesize=numtofit+FitRedCoeff+FitDMCoeff+totalBandNoiseCoeff+6*pulse->TNBandDMC+totalGroupCoeff;

    double **TotalMatrix=new double*[pulse->nobs];
    for(int i =0;i<pulse->nobs;i++){
        TotalMatrix[i]=new double[totalsize];
        for(int j =0;j<totalsize; j++){
            TotalMatrix[i][j]=0;
        }
    }


    for(int i =0;i<pulse->nobs;i++){
        int startpoint=0;
        double time=((double)pulse->obsn[i].bat);

        for(int j =0;j<numtofit; j++){
            TotalMatrix[i][j]=TNDM[i][j];
        }

        startpoint+=numtofit;		
        for(int j =0;j<totCoeff; j++){
            TotalMatrix[i][j+startpoint]=FMatrix[i][j];
        }

        startpoint += totCoeff;
        int DMEvterms=0;
        for(int e =0; e < pulse->nDMEvents; e++){
            if(time < DMEventInfo[e][0]+DMEventInfo[e][1] && time > DMEventInfo[e][0]){
                if(pulse->TNDMEvOff[e]==1){	
                    TotalMatrix[i][startpoint+DMEvterms]=DMVec[i];
                    DMEvterms++;
                }
                if(pulse->TNDMEvLin[e]==1){	
                    TotalMatrix[i][startpoint+DMEvterms]=(time-DMEventInfo[e][0])*DMVec[i];
                    DMEvterms++;
                }
                if(pulse->TNDMEvQuad[e]==1){	
                    TotalMatrix[i][startpoint+DMEvterms]=pow((time-DMEventInfo[e][0]),2)*DMVec[i];
                    DMEvterms++;
                }
            }
        }

        startpoint += DMEventQuadTerms;
        int ShapeEvterms=0;
        for(int e =0; e < pulse->nTNShapeletEvents; e++){

            double *shapeVec =  new double[pulse->TNShapeletEvN[e]];
            double HVal=(time-pulse->TNShapeletEvPos[e])/(sqrt(2.0)*pulse->TNShapeletEvWidth[e]);
            othpl(pulse->TNShapeletEvN[e],HVal,shapeVec);

            for(int s=0; s < pulse->TNShapeletEvN[e]; s++){

                double NormTerm=1.0/sqrt(sqrt(2.0*M_PI)*pow(2.0,s));
                TotalMatrix[i][startpoint+ShapeEvterms] = NormTerm*shapeVec[s]*exp(-0.5*pow((time-pulse->TNShapeletEvPos[e])/pulse->TNShapeletEvWidth[e], 2))*pow(DMVec[i], pulse->TNShapeletEvFScale[e]/2.0);

                ShapeEvterms++;

            }
            delete[] shapeVec;
        }

    }

    //////////////////////////////////////////////////////////////////////////////////////////  
    //////////////////////Get Residuals Vector////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    double *Resvec=new double[pulse->nobs];
    for(int o=0;o<pulse->nobs; o++){
        Resvec[o]=(double)pulse->obsn[o].residual;
        pulse->obsn[o].prefitResidual = pulse->obsn[o].residual;
    }




    //////////////////////////////////////////////////////////////////////////////////////////  
    //////////////////////Get White Noise Vector//////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    double *Noise=new double[pulse->nobs];
    for(int o=0;o<pulse->nobs; o++){
        Noise[o]=pow(((pulse->obsn[o].toaErr)*pow(10.0,-6)),2);
    }


    //////////////////////////////////////////////////////////////////////////////////////////  
    //////////////////////Do Algebra//////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    double **NG = new double*[pulse->nobs]; for (int k=0; k<pulse->nobs; k++) NG[k] = new double[totalsize];
    double **GNG = new double*[totalsize]; for (int k=0; k<totalsize; k++) GNG[k] = new double[totalsize];



    for(int i=0;i<pulse->nobs;i++){
        for(int j=0;j<totalsize; j++){

            NG[i][j]=TotalMatrix[i][j]/Noise[i];

        }
    }


    dgemm(TotalMatrix, NG,GNG,pulse->nobs, totalsize,pulse->nobs, totalsize, 'T','N');


    for(int j =0;j<totCoeff; j++){
        GNG[numtofit+j][numtofit+j]+=PPFM[j][j];
    }


    double tdet=0;
    dpotrf(GNG, totalsize, tdet);
    dpotri(GNG,totalsize);
    //	printf("logdet GNG = %g\n",tdet);

    double *dG=new double[totalsize];
    dgemv(NG,Resvec,dG,pulse->nobs,totalsize,'T');

    double *maxcoeff=new double[totalsize];
    dgemv(GNG,dG,maxcoeff,totalsize,totalsize,'N');

    longdouble *Errorvec=new longdouble[totalsize];

    for(int i =0; i < totalsize; i++){
        Errorvec[i]=pow(GNG[i][i], 0.5);
    }


    double *Scoeff;
    long double *Serr;

    if(useOrthogonal == 1){

        Scoeff=new double[numtofit];
        Serr=new long double[numtofit];
        for(int i =0; i < numtofit; i++){
            if(S[i] >= 0){
                Scoeff[i]=maxcoeff[i]/S[i];
                Serr[i]=Errorvec[i]/S[i];
            }
            else{
                Scoeff[i]=0;
                Serr[i]=0;
            }
        }
    }

    double *TempoCoeff = new double[numtofit];
    double *TempoErr =  new double[numtofit];


    if(useOrthogonal==1){
        dgemv(V,Scoeff,TempoCoeff,numtofit,numtofit, 'N');
    }

    for(int i=0;i<numtofit; i++){
        if(useOrthogonal==1){
            long double errsum=0;
            for(int j=0;j<numtofit; j++){
                errsum += pow(V[i][j]*Serr[j],2);
            }
            TempoCoeff[i]=TempoCoeff[i]/TNDMScale[i];
            TempoErr[i]=pow(errsum,0.5)/TNDMScale[i];
        }
        if(useOrthogonal==0){
            TempoCoeff[i]=maxcoeff[i]/TNDMScale[i];
            TempoErr[i]=((double)Errorvec[i])/TNDMScale[i];
        }
    }
    updateParameters(pulse,0,TempoCoeff,TempoErr);


    //////////////////Get Red noise and DM Coeffs and errors/////////////////////

    int redcounter = 0;
    int dmcounter = 0;
    for(int j=0;j<totalsize; j++){
        if(j>=numtofit && j < numtofit+FitRedCoeff){
            //pulse->TNRedCoeffs[redcounter] = maxcoeff[j];
            //pulse->TNRedCoeffs[redcounter+100] = Errorvec[j];
            //printf("Red Coeff: %i %g +/- %g \n", redcounter, pulse->TNRedCoeffs[redcounter], pulse->TNRedCoeffs[redcounter+100]);
            redcounter++;
        }

        if(j>=FitRedCoeff+numtofit && j < numtofit+FitRedCoeff+FitDMCoeff){
            //pulse->TNDMCoeffs[dmcounter] = maxcoeff[j];
            //pulse->TNDMCoeffs[dmcounter+100] = Errorvec[j];
            //printf("DM Coeff: %i %g +/- %g \n", dmcounter, pulse->TNDMCoeffs[dmcounter], pulse->TNDMCoeffs[dmcounter+100]);
            dmcounter++;
        }
    }

    double chisq=0;
    printf("subtractRed is %i \n", pulse->TNsubtractRed);
    printf("subtractDM is %i \n", pulse->TNsubtractDM);

    for(int i=0;i<pulse->nobs;i++){
        double dsum=0;
        double redsum=0;
        double rederr=0;
        double dmsum=0;
        double dmerr=0;
        double shapesum=0;
        double shapeerr=0;
        for(int j=0;j<totalsize; j++){
            //			if(i==0)printf("Max coeff: %i %g \n", j,maxcoeff[j]);		
            dsum=dsum+TotalMatrix[i][j]*maxcoeff[j];

            if(j>=numtofit && j < numtofit+FitRedCoeff){
                //if(i==20){
                //printf("TM: %i %g %g \n", j-numtofit, TotalMatrix[i][j], maxcoeff[j]);
                //}
                redsum+=TotalMatrix[i][j]*maxcoeff[j];
                rederr+=pow(TotalMatrix[i][j]*Errorvec[j],2);
            }
            if(j>=FitRedCoeff+numtofit && j < numtofit+FitRedCoeff+FitDMCoeff){
                //if(i==20){
                //printf("TM: %i %g %g \n", j-numtofit, TotalMatrix[i][j], maxcoeff[j]);
                //}

                dmsum+=TotalMatrix[i][j]*maxcoeff[j];
                dmerr+=pow(TotalMatrix[i][j]*Errorvec[j],2);
            }

            if(j>=numtofit+FitRedCoeff+FitDMCoeff && j < totalnoisesize){
                redsum+=TotalMatrix[i][j]*maxcoeff[j];
                rederr+=pow(TotalMatrix[i][j]*Errorvec[j],2);
            }		
            if(j>=numtofit+totCoeff && j < numtofit+totCoeff+ShapeEventTerms){
                //					printf("Shapeevent terms %i %i %g %g \n", i, j, TotalMatrix[i][j]*maxcoeff[j], pow(TotalMatrix[i][j]*Errorvec[j],2));
                dmsum+=TotalMatrix[i][j]*maxcoeff[j];
                dmerr+=pow(TotalMatrix[i][j]*Errorvec[j],2);
                shapesum+=TotalMatrix[i][j]*maxcoeff[j];
                shapeerr+=pow(TotalMatrix[i][j]*Errorvec[j],2);
            }		


        }

        //		if(fabs(shapesum) > pow(10.0, -10)){printf("Shapeevent terms %i %.10g %g %g \n", i, (double)pulse->obsn[i].bat, shapesum/DMVec[i], sqrt(shapeerr)/DMVec[i]);}


        longdouble yrs = (pulse->obsn[i].bat - pulse->param[param_dmepoch].val[0])/365.25;
        longdouble arg = 1.0;
        double dmDot=0;
        double dmDotErr=0;
        for (int d=0;d<9;d++){
            if(d>0){
                arg *= yrs;
            }
            if (pulse->param[param_dm].paramSet[d]==1){
                if(d>0){
                    dmDot+=(double)(pulse->param[param_dm].val[d]*arg);
                }
                dmDotErr+=pow((double)(pulse->param[param_dm].err[d]*arg),2);
            }
        }

        double pDotErr=0;
        pDotErr+=pow(TempoErr[0],2);
        if(pulse->param[param_f].paramSet[0]==1){
            arg=((pulse->obsn[i].bat - pulse->param[param_pepoch].val[0])/ \
                    pulse[0].param[param_f].val[0])*86400.0;
            pDotErr+=pow((double)(pulse->param[param_f].err[0]*arg),2);

        }

        if(pulse->param[param_f].paramSet[1]==1){
            arg=0.5*pow((double)(pulse->obsn[i].bat - pulse->param[param_pepoch].val[0]), 2);
            longdouble argerr = (pulse->param[param_f].err[1]/ \
                    pulse[0].param[param_f].val[0])*86400.0*86400;
            pDotErr+=pow((double)(argerr*arg),2);
        }


        chisq+=(Resvec[i]-dsum)*(Resvec[i]-dsum)/(Noise[i]);
        pulse->obsn[i].TNRedSignal=redsum;
        pulse->obsn[i].TNRedErr=pow(rederr+pDotErr,0.5);

        pulse->obsn[i].TNDMSignal=dmsum;
        pulse->obsn[i].TNDMErr=pow(dmerr/pow(DMVec[i],2) + dmDotErr,0.5);
    }

    //////////////////////////////////////////////////////////////////////////////////////////  
    //////////////////////Calculate GLS Chi-squared///////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    for (int i=0; i<pulse->nobs;i++){
        for (int j=0; j<numtofit;j++){
            Resvec[i] -= TotalMatrix[i][j]*maxcoeff[j];
        }
    }

    int nother = totalsize-numtofit;
    double **TMat = new double*[pulse->nobs];
    for (int k=0;k<pulse->nobs;k++){
        TMat[k] = new double[nother];
    }

    for(int i=0;i<pulse->nobs;i++){
        for(int j=0;j<nother; j++){
            TMat[i][j] = TotalMatrix[i][j+numtofit];
        }
    }

    double **NT = new double*[pulse->nobs];
    for (int k=0;k<pulse->nobs;k++){
        NT[k] = new double[nother];
    }

    double **TNT = new double*[nother]; 
    for (int k=0; k<nother; k++){
        TNT[k] = new double[nother];
    }

    for(int i=0;i<pulse->nobs;i++){
        for(int j=0;j<nother; j++){
            NT[i][j] = TMat[i][j]/Noise[i];
        }
    }

    dgemm(TMat, NT, TNT, pulse->nobs, nother, pulse->nobs, nother, 'T', 'N');

    for(int j =0;j<totCoeff; j++){
        TNT[j][j] += PPFM[j][j];
    }

    tdet = 0.0;
    dpotrf(TNT, nother, tdet);
    dpotri(TNT, nother);

    double *dT = new double[nother];
    dgemv(NT, Resvec, dT, pulse->nobs, nother, 'T');

    double *Sigmad = new double[nother];
    dgemv(TNT, dT, Sigmad, nother, nother, 'N');

    chisq = 0.0;
    double timesq=0;
    for (int i=0; i<pulse->nobs; i++){
        timesq+=Resvec[i]*Resvec[i]/Noise[i];
        chisq += Resvec[i]*Resvec[i]/Noise[i];
    }


    double freqsq=0;
    for (int k=0; k<nother; k++){
        freqsq+=dT[k]*Sigmad[k];
        chisq -= dT[k]*Sigmad[k];
    }

    pulse->fitChisq = chisq; 
    pulse->fitNfree = pulse->nobs-numtofit;
    pulse->nFit=pulse->nobs;
    pulse->nParam = numtofit;


    if(useOrthogonal == 1){
        for (int j = 0; j < pulse->nobs; j++){
            delete[]U[j];
        }
        delete[]U;

        delete[]S;
        delete[]Scoeff;	
        delete[]Serr;	

        for (int j = 0; j < numtofit; j++){
            delete[]VT[j];
            delete[]V[j];
        }

        delete[]V; 
        delete[]VT;
    }
    delete[] TempoErr;
    delete[] TempoCoeff;
    delete[] DMVec;
    delete[] dT;
    delete[] Sigmad;
    for(int i=0;i<pulse->nobs;i++){
        delete[] TNDM[i];
    }
    delete[] TNDM;
    for(int i=0; i<pulse->nobs; i++){
        delete[] TMat[i];
        delete[] NT[i];
    }
    delete[] TMat;
    delete[] NT;
    for(int i=0; i<nother; i++){
        delete[] TNT[i];
    }
    delete[] TNT;

    delete[] dG;
    delete[] maxcoeff;
    delete[] Errorvec;
    delete[] Resvec;
    delete[] Noise;

    for (int k=0; k<pulse->nobs; k++){
        delete[] NG[k];
    }
    delete[] NG;
    for (int k=0; k<totalsize; k++){
        delete[] GNG[k];
    }
    delete[] GNG;

    for(int i =0;i<pulse->nobs;i++){
        delete[] TotalMatrix[i];
    }
    delete[] TotalMatrix;

    for(int i=0;i<totCoeff;i++){
        delete[] PPFM[i];
    }
    delete[] PPFM;

    delete[] powercoeff;
    delete[] freqs;
    for(int i=0;i<pulse->nobs;i++){
        delete[] FMatrix[i];
    }
    delete[] FMatrix;
    delete[]Processed;

}


void dgesvd(double **A, int m, int n, double *S, double **U, double **VT)
{
    char jobu, jobvt;
    int lda, ldu, ldvt, lwork, info;
    double *a, *u, *vt, *work;


    jobu = 'A'; /* Specifies options for computing U.
A: all M columns of U are returned in array U;
S: the first min(m,n) columns of U (the left
singular vectors) are returned in the array U;
O: the first min(m,n) columns of U (the left
singular vectors) are overwritten on the array A;
N: no columns of U (no left singular vectors) are
computed. */

    jobvt = 'A'; /* Specifies options for computing VT.
A: all N rows of V**T are returned in the array
VT;
S: the first min(m,n) rows of V**T (the right
singular vectors) are returned in the array VT;
O: the first min(m,n) rows of V**T (the right
singular vectors) are overwritten on the array A;
N: no rows of V**T (no right singular vectors) are
computed. */

    lda = m; // The leading dimension of the matrix a.
    a = dgesvd_ctof(A, lda, n); /* Convert the matrix A from double pointer
                                   C form to single pointer Fortran form. */

    ldu = m;



    ldu = m; // Left singular vector matrix
    u = new double[ldu*ldu];

    ldvt = n; // Right singular vector matrix
    vt = new double[ldvt*n];

    int LMAX=100000;

    work = new double[LMAX];
    lwork = -1; // Set up the work array, larger than needed.

    // 	printf("parm 11 %i %i\n",ldu,ldvt);
    dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, S, u,&ldu, vt, &ldvt, work, &lwork, &info);

    lwork = std::min(LMAX,int(work[0]));

    dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, S, u,&ldu, vt, &ldvt, work, &lwork, &info);
    // 	printf("parm 11 out %i %i\n",ldu,ldvt);
    dgesvd_ftoc(u, U, ldu, ldu);
    dgesvd_ftoc(vt, VT, ldvt, n);

    delete a;
    delete u;
    delete vt;
    delete work;
}


double* dgesvd_ctof(double **in, int rows, int cols)
{
    double *out;
    int i, j;

    out = new double[rows*cols];
    for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
    return(out);
}


void dgesvd_ftoc(double *in, double **out, int rows, int cols)
{
    int i, j;

    for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}


void dgemv(double **A, double *vecin,double *vecout,int rowa, int cola, char AT)
{

    double *a;

    double alpha=1;
    double beta=0;
    int incX=1;
    int incY=1;

    a = dgemv_ctof(A, rowa, cola); 

    dgemv_(&AT, &rowa, &cola, &alpha, a, &rowa, vecin, &incX, &beta, vecout, &incY);

    delete a;
}


double* dgemv_ctof(double **in, int rows, int cols)
{
    double *out;
    int i, j;

    out = new double[rows*cols];
    for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
    return(out);
}


void dgemv_ftoc(double *in, double **out, int rows, int cols)
{
    int i, j;

    for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}

void dgemm(double **A, double **B,double **C,int rowa, int cola, int rowb, int colb, char AT, char BT)
{

    int M,N,K;
    double *a, *b, *c;

    double alpha=1;
    double beta=0;


    if(AT == 'N'){
        M=rowa;
        K=cola;
    }
    else if(AT == 'T'){
        M=cola;
        K=rowa;
    }

    if(BT == 'N'){
        N=colb;
    }
    else if(BT == 'T'){
        N=rowa;
    }	
    /*
       (TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)*/


    a = dgemm_ctof(A, rowa, cola); 
    b = dgemm_ctof(B, rowb, colb);
    c = new double[M*N];

    dgemm_(&AT, &BT, &M, &N, &K, &alpha, a, &rowa,b, &rowb, &beta, c, &M);
    dgemm_ftoc(c, C, M, N);


    delete a;
    delete b;
    delete c;

}


double* dgemm_ctof(double **in, int rows, int cols)
{
    double *out;
    int i, j;

    out = new double[rows*cols];
    for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
    return(out);
}


void dgemm_ftoc(double *in, double **out, int rows, int cols)
{
    int i, j;

    for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}

void dpotri(double **A, int msize)
{

    int info;
    double *a;
    char UPLO='L';

    a = dpotri_ctof(A, msize, msize); 

    dpotri_(&UPLO, &msize, a, &msize, &info);
    dpotri_ftoc(a, A, msize, msize);

    for(int i=0;i<msize;i++){
        for(int j=0;j<i;j++){
            A[j][i]=A[i][j];
        }
    }


    delete a;

}


double* dpotri_ctof(double **in, int rows, int cols)
{
    double *out;
    int i, j;

    out = new double[rows*cols];
    for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
    return(out);
}


void dpotri_ftoc(double *in, double **out, int rows, int cols)
{
    int i, j;

    for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}

void dpotrf(double **A, int msize, double &det)
{

    int info;
    double *a;
    char UPLO='L';

    a = dpotrf_ctof(A, msize, msize); 

    dpotrf_(&UPLO, &msize, a, &msize, &info);
    dpotrf_ftoc(a, A, msize, msize);

    det=0;
    for(int i=0;i<msize;i++){
        det+=log(A[i][i]);
    }

    det=det*2;

    //printf("info: %i \n", info);
    delete a;

}


double* dpotrf_ctof(double **in, int rows, int cols)
{
    double *out;
    int i, j;

    out = new double[rows*cols];
    for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
    return(out);
}


void dpotrf_ftoc(double *in, double **out, int rows, int cols)
{
    int i, j;

    for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}

#endif
#endif
