#include <gtest/gtest.h>
#include <src/gtest_main.cc>
#include "tempo2.h"

TEST(testTempo2, maxpsrset){
    ASSERT_GT(MAX_PSR,0);
}

TEST(testTempo2, sizeoflongdouble){
    EXPECT_EQ(sizeof(longdouble),16);
    ASSERT_GT(sizeof(longdouble),12);
}

