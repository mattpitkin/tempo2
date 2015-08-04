#include <gtest/gtest.h>
#include <src/gtest_main.cc>
#include "tempo2.h"

TEST(testEndToEnd,checkIdeal){
    pulsar _psr;
    pulsar *psr = &_psr;
    MAX_PSR=1;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
    char** argv=NULL;
    int npsr=1;
    initialise(psr,0); /* Initialise all */

    strcpy(parFile[0],"test3.par");
    strcpy(timFile[0],"test3.tim");

    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);
    preProcessSimple(psr);

    psr->noWarnings=2;
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        ASSERT_LT(static_cast<double>(fabsl(psr->obsn[iobs].residual)),1e-10l) << "test3.par test3.tim do not give idealised ToAs";
    }

    doFitAll(psr,npsr,"NULL");

    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        ASSERT_LT(static_cast<double>(fabsl(psr->obsn[iobs].residual)),1e-10l) << "Fitting has caused error";
    }

}



TEST(testEndToEnd,checkFit){
    pulsar _psr;
    pulsar *psr = &_psr;
    MAX_PSR=1;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
    char** argv=NULL;
    int npsr=1;
    initialise(psr,0); /* Initialise all */

    strcpy(parFile[0],"test3b.par");
    strcpy(timFile[0],"test3.tim");

    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);
    preProcessSimple(psr);

    psr->noWarnings=2;
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    doFitAll(psr,npsr,"NULL");
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    doFitAll(psr,npsr,"NULL");

    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        ASSERT_LT(static_cast<double>(fabsl(psr->obsn[iobs].residual)),1e-10l) << "Fitting has caused error";
    }

}

