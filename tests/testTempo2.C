#include <gtest/gtest.h>
#include <src/gtest_main.cc>
#include "tempo2.h"

#ifdef LONGDOUBLE_IS_FLOAT128
#include <quadmath.h>
#endif
TEST(testTempo2, maxpsrset){
    ASSERT_GT(MAX_PSR,0);
}

TEST(testTempo2, sizeoflongdouble){
    EXPECT_EQ(sizeof(longdouble),16);
#ifdef LONGDOUBLE_IS_FLOAT128
    EXPECT_GT(FLT128_DIG,32);
#else
    EXPECT_GT(LDBL_DIG,17);
#endif
}

TEST(testTempo2, printlongdouble){
    longdouble ld = longdouble(123.0);
    char sb[1024];
    ld_sprintf(sb,"%.1Lf",ld);
    ASSERT_STREQ(sb,"123.0");
    ld = longdouble(0.123456789012345678901234567895);
    ld_sprintf(sb,"%.30Lf",ld);
    ASSERT_STREQ(sb,"0.123456789012345678901234567895");
    longdouble ld2 = parse_longdouble(sb);
    ASSERT_EQ(ld2,ld);
}
