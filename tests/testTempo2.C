#include <gtest/gtest.h>
#include <src/gtest_main.cc>
#include <quadmath.h>
#include "tempo2.h"

TEST(testTempo2, maxpsrset){
    ASSERT_GT(MAX_PSR,0);
}

TEST(testTempo2, sizeoflongdouble){
    EXPECT_EQ(sizeof(longdouble),16);
    EXPECT_GT(FLT128_DIG,32);
    EXPECT_GT(LDBL_DIG,32);
    ASSERT_GT(sizeof(longdouble),12);
}

