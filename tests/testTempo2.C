#include <gtest/gtest.h>
#include <src/gtest_main.cc>
#include "gtest_extensions.h"
#include "tempo2.h"
#ifndef DATDIR
#define DATDIR .
#endif
#define NANOSEC 1e-9

#ifdef LONGDOUBLE_IS_FLOAT128
#include <quadmath.h>
#endif
TEST(testTempo2h, maxpsrset){
    ASSERT_GT(MAX_PSR,0);
}

TEST(testLongDouble, precision){
#ifdef LONGDOUBLE_IS_FLOAT128
    NOTE("Tempo2 longdoubles are 128-bit IEEE 754 binary128 implemented as __float128 in GCC\n");
    NOTE("longdouble has %d bits in mantisa, and exponent between %d and %d\n",FLT128_MANT_DIG,FLT128_MIN_EXP,FLT128_MAX_EXP);
#endif

#ifdef LONGDOUBLE_IS_LD
    NOTE("Tempo2 longdoubles are GCC long double, probably 80-bit x86 Extended Precision\n");
    NOTE("longdouble has %d bits in mantisa, and exponent between %d and %d\n",LDBL_MANT_DIG,LDBL_MIN_EXP,LDBL_MAX_EXP);
#endif

    EXPECT_EQ(sizeof(longdouble),16);
#ifdef LONGDOUBLE_IS_FLOAT128
    EXPECT_GT(FLT128_DIG,32);
#else
    EXPECT_GT(LDBL_DIG,17);
#endif
}

TEST(testLongDouble, printAndParse){
    longdouble ld = longdouble(123.0);
    char sb[1024];
    ld_sprintf(sb,"%.1Lf",ld);
    ASSERT_STREQ(sb,"123.0");
    ld = longdouble(50000.12345678912345);
    ld_sprintf(sb,"%.14Lf",ld);

    ASSERT_STREQ("50000.12345678912345",sb);

    longdouble ld2 = parse_longdouble(sb);
    ASSERT_EQ(ld,ld2);

    unsigned uu = 2147483649;
    long long int ll = 2147483649L;
    long long int ll2 = 17179869185L;

    ld_sprintf(sb,"%% %d %.1f %u %.1lf %lld %lld %s %.1Lf %% %c",-1,0.1,uu,0.1,ll,ll2,"t",ld,'x');
    ASSERT_STREQ("% -1 0.1 2147483649 0.1 2147483649 17179869185 t 50000.1 % x",sb);
}

TEST(testFormResiduals, basicBATs){
    pulsar _psr;
    pulsar *psr = &_psr;
    MAX_PSR=1;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
    int npsr=1;
    initialise(psr,0); /* Initialise all */

    strcpy(parFile[0],DATDIR "/test1.par");
    strcpy(timFile[0],DATDIR "/test1.tim");

    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);
    psr->noWarnings=2;
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    ASSERT_LT(static_cast<double>(fabsl(psr[0].obsn[1].residual)),NANOSEC);
    ASSERT_LT(static_cast<double>(fabsl(psr[0].obsn[2].residual-longdouble(0.4))),NANOSEC);
}




TEST(testFormResiduals, subtractBATs){
    pulsar _psr;
    pulsar *psr = &_psr;
    MAX_PSR=1;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
    int npsr=1;
    initialise(psr,0); /* Initialise all */

    strcpy(parFile[0],DATDIR "/test1.par");
    strcpy(timFile[0],DATDIR "/test1.tim");

    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);

    for(int iobs = 0; iobs < psr->nobs; iobs++){
        psr[0].obsn[iobs].bat = psr[0].obsn[iobs].sat;
        psr[0].obsn[iobs].bbat = psr[0].obsn[iobs].sat;
        psr[0].obsn[iobs].delayCorr=0;
    }
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        psr[0].obsn[iobs].sat -= psr[0].obsn[iobs].residual/SECDAYl;
        psr[0].obsn[iobs].bat = psr[0].obsn[iobs].sat;
        psr[0].obsn[iobs].bbat = psr[0].obsn[iobs].sat;
    }
    formResiduals(psr,npsr,0);

    for(int iobs = 0; iobs < psr->nobs; iobs++){
        EXPECT_LT(static_cast<double>(fabsl(psr[0].obsn[iobs].residual)),NANOSEC) << "Precision lost in formResiduals (s)";
    }

    for(int iobs = 1; iobs < psr->nobs; iobs++){
        psr[0].obsn[iobs].sat += longdouble(4e-9)/SECDAYl;
        psr[0].obsn[iobs].bat = psr[0].obsn[iobs].sat;
        psr[0].obsn[iobs].bbat = psr[0].obsn[iobs].sat;
    }
    formResiduals(psr,npsr,0);

    for(int iobs = 1; iobs < psr->nobs; iobs++){
        EXPECT_NEAR(longdouble(4.0e-9),fabsl(psr[0].obsn[iobs].residual),longdouble(NANOSEC)) << "Precision lost in formResiduals (s)";
    }
}


TEST(testFormResiduals, subtractSATs){
    pulsar _psr;
    pulsar *psr = &_psr;
    MAX_PSR=1;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
    int npsr=1;
    initialise(psr,0); /* Initialise all */

    strcpy(parFile[0],DATDIR "/test1.par");
    strcpy(timFile[0],DATDIR "/test2.tim");

    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);
    psr->noWarnings=2;
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        psr[0].obsn[iobs].sat -= psr[0].obsn[iobs].residual/SECDAYl;
    }
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        psr[0].obsn[iobs].sat -= psr[0].obsn[iobs].residual/SECDAYl;
    }
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        psr[0].obsn[iobs].sat -= psr[0].obsn[iobs].residual/SECDAYl;
    }
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);

    for(int iobs = 0; iobs < psr->nobs; iobs++){
        EXPECT_LT(static_cast<double>(fabsl(psr[0].obsn[iobs].residual)),NANOSEC) << "Precision lost in formBats or formResiduals";
    }
}

TEST(testFormBats, offsetSATs){
    pulsar _psr;
    pulsar *psr = &_psr;
    MAX_PSR=1;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
    int npsr=1;
    initialise(psr,0); /* Initialise all */

    strcpy(parFile[0],DATDIR "/test1.par");
    strcpy(timFile[0],DATDIR "/test2.tim");

    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);
    psr->noWarnings=2;
    formBatsAll(psr,npsr);

    for(int iobs = 0; iobs < psr->nobs; iobs++){
        psr->obsn[iobs].prefitResidual = psr->obsn[iobs].bat;
        psr->obsn[iobs].sat += longdouble(4e-9)/SECDAYl;
    }
    formBatsAll(psr,npsr);

    for(int iobs = 0; iobs < psr->nobs; iobs++){
        EXPECT_LT(static_cast<double>(fabsl(psr->obsn[iobs].bat - psr->obsn[iobs].prefitResidual - longdouble(4e-9)/SECDAYl)),NANOSEC) << "Precision lost in formBats";
    }

}
