AM_DEFAULT_SOURCE_EXT=.C

test_programs=testTKfit testTempo2 testEndToEnd

noinst_LTLIBRARIES=libgtest.la
check_PROGRAMS=$(test_programs)

TESTS=$(test_programs)



CPPFLAGS=-Igtest-1.7.0/include -Igtest-1.7.0

libgtest_la_SOURCES=gtest-1.7.0/src/gtest-all.cc

LDADD=libgtest.la ../libtempo2.la



