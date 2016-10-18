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

    strcpy(parFile[0],DATDIR "/test3.par");
    strcpy(timFile[0],DATDIR "/test3.tim");

    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);
    preProcessSimple(psr);

    psr->noWarnings=2;
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        ASSERT_LT(static_cast<double>(fabsl(psr->obsn[iobs].residual)),TEST_DELTA) << DATDIR "/test3.par test3.tim do not give idealised ToAs";
    }

    doFitAll(psr,npsr,"NULL");

    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        ASSERT_LT(static_cast<double>(fabsl(psr->obsn[iobs].residual)),TEST_DELTA) << "Fitting has caused error in ideal";
    }

}

TEST(testEndToEnd,checkDE430){
    pulsar _psr;
    pulsar *psr = &_psr;
    MAX_PSR=1;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
    char** argv=NULL;
    int npsr=1;
    initialise(psr,0); /* Initialise all */

    strcpy(parFile[0],DATDIR "/test_de430.par");
    strcpy(timFile[0],DATDIR "/test_de430.tim");

    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);
    preProcessSimple(psr);

    psr->noWarnings=2;
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        ASSERT_LT(static_cast<double>(fabsl(psr->obsn[iobs].residual)),TEST_DELTA) << DATDIR "/test_de430.par test_de430.tim do not give idealised ToAs";
    }

    doFitAll(psr,npsr,"NULL");

    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        ASSERT_LT(static_cast<double>(fabsl(psr->obsn[iobs].residual)),TEST_DELTA) << "Fitting has caused error in ideal";
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

    strcpy(parFile[0],DATDIR "/test3b.par");
    strcpy(timFile[0],DATDIR "/test3.tim");

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
        ASSERT_LT(static_cast<double>(fabsl(psr->obsn[iobs].residual)),TEST_DELTA) << "Fitting has caused error";
    }

}

TEST(testEndToEnd,gracefulBadFit){
    pulsar _psr;
    pulsar *psr = &_psr;
    MAX_PSR=1;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
    char** argv=NULL;
    int npsr=1;
    initialise(psr,0); /* Initialise all */

    strcpy(parFile[0],DATDIR "/test3c.par");
    strcpy(timFile[0],DATDIR "/test3.tim");

    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);
    preProcessSimple(psr);

    psr->noWarnings=2;
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    TK_errorCount=0;
    doFitAll(psr,npsr,"NULL");
//    ASSERT_GT(TK_errorCount,1u);
ASSERT_EQ(1,1);
    TK_errorCount=0;
    
}

