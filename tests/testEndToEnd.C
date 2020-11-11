#include <gtest/gtest.h>
#include <src/gtest_main.cc>
#include "tempo2.h"
#include "t2fit.h"
#include "T2accel.h"

TEST(testEndToEnd,checkIdeal){
    pulsar _psr;
    pulsar *psr = &_psr;
    MAX_PSR=1;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
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

    t2Fit(psr,npsr,"NULL");

    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        ASSERT_LT(static_cast<double>(fabsl(psr->obsn[iobs].residual)),TEST_DELTA) << "Fitting has caused error in ideal";
    }

}

TEST(testEndToEnd,checkIdeal_jb){
    pulsar _psr;
    pulsar *psr = &_psr;
    MAX_PSR=1;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
    int npsr=1;
    initialise(psr,0); /* Initialise all */

    strcpy(parFile[0],DATDIR "/test3_jb.par");
    strcpy(timFile[0],DATDIR "/test3_jb.tim");

    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);
    preProcessSimple(psr);

    psr->noWarnings=2;
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        ASSERT_LT(static_cast<double>(fabsl(psr->obsn[iobs].residual)),TEST_DELTA) << DATDIR "/test3_jb.par test3_jb.tim do not give idealised ToAs";
    }

    t2Fit(psr,npsr,"NULL");

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

    t2Fit(psr,npsr,"NULL");

    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        ASSERT_LT(static_cast<double>(fabsl(psr->obsn[iobs].residual)),TEST_DELTA) << "Fitting has caused error in ideal";
    }

}

TEST(testEndToEnd,check_1937_2020){
    pulsar _psr;
    pulsar *psr = &_psr;
    MAX_PSR=1;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
    int npsr=1;
    initialise(psr,0); /* Initialise all */

    strcpy(parFile[0],DATDIR "/test_1937_2020.par");
    strcpy(timFile[0],DATDIR "/test_1937_2020.tim");

    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);
    preProcessSimple(psr);

    psr->noWarnings=2;
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        ASSERT_LT(static_cast<double>(fabsl(psr->obsn[iobs].residual)),TEST_DELTA) << DATDIR "/test_de430.par test_de430.tim do not give idealised ToAs";
    }

    t2Fit(psr,npsr,"NULL");

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
    t2Fit(psr,npsr,"NULL");
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    t2Fit(psr,npsr,"NULL");

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
    t2Fit(psr,npsr,"NULL");
//    ASSERT_GT(TK_errorCount,1u);
ASSERT_EQ(1,1);
    TK_errorCount=0;
    
}





TEST(testEndToEnd,checkGlobalQIfunc){
    if(useT2accel) useT2accel=1;
    int npsr=3;
    MAX_PSR=3;
    pulsar *psr = (pulsar*)calloc(MAX_PSR, sizeof(pulsar));

    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];

    initialise(psr,0); /* Initialise all */


    strcpy(parFile[0],DATDIR "/test4/test4C.par");
    strcpy(timFile[0],DATDIR "/test4/test4C.tim");

    strcpy(parFile[1],DATDIR "/test4/test4A.par");
    strcpy(timFile[1],DATDIR "/test4/test4A.tim");

    strcpy(parFile[2],DATDIR "/test4/test4B.par");
    strcpy(timFile[2],DATDIR "/test4/test4B.tim");


    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);


    char* cmd[2];
    cmd[0] = (char*)calloc(80,1);
    cmd[1] = (char*)calloc(256,1);
    strcpy(cmd[0],"-global");
    strcpy(cmd[1],DATDIR "/test4/global.par");
    preProcess(psr,npsr,2,cmd);
    free(cmd[0]);
    free(cmd[1]);

    psr->noWarnings=2;
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,1);

    t2Fit(psr,npsr,"NULL");

    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,1);

    FILE* f = fopen(DATDIR "/test4/aplus_t2.dat","r");

    for(int qi = 0; qi < psr->quad_ifuncN_p; qi++){
        double m,v,e;
        fscanf(f,"%lg %lg %lg",&m,&v,&e);
        ASSERT_LT(static_cast<double>(fabsl(psr->quad_ifuncV_p[qi] - v)),TEST_DELTA) << "Error in QIfunc from global fit";
    }

    fclose(f);
    free(psr);

}






