AM_DEFAULT_SOURCE_EXT=.C

test_programs=testTKfit testTempo2 testEndToEnd
LOG_DRIVER=$(srcdir)/custom-test-driver

noinst_LTLIBRARIES=libgtest.la
check_PROGRAMS=$(test_programs)

TESTS=$(test_programs)

gtest_src=gtest-1.7.0/src/gtest-death-test.cc gtest-1.7.0/src/gtest-filepath.cc gtest-1.7.0/src/gtest-internal-inl.h gtest-1.7.0/src/gtest-port.cc gtest-1.7.0/src/gtest-printers.cc gtest-1.7.0/src/gtest-test-part.cc gtest-1.7.0/src/gtest-typed-test.cc gtest-1.7.0/src/gtest.cc gtest-1.7.0/src/gtest_main.cc

gtest_hdr=gtest-1.7.0/include/gtest/gtest-death-test.h gtest-1.7.0/include/gtest/gtest-message.h gtest-1.7.0/include/gtest/gtest-param-test.h gtest-1.7.0/include/gtest/gtest-printers.h gtest-1.7.0/include/gtest/gtest-spi.h gtest-1.7.0/include/gtest/gtest-test-part.h gtest-1.7.0/include/gtest/gtest-typed-test.h gtest-1.7.0/include/gtest/gtest.h gtest-1.7.0/include/gtest/gtest_pred_impl.h gtest-1.7.0/include/gtest/gtest_prod.h gtest-1.7.0/include/gtest/internal/gtest-death-test-internal.h gtest-1.7.0/include/gtest/internal/gtest-filepath.h gtest-1.7.0/include/gtest/internal/gtest-internal.h gtest-1.7.0/include/gtest/internal/gtest-linked_ptr.h gtest-1.7.0/include/gtest/internal/gtest-param-util-generated.h gtest-1.7.0/include/gtest/internal/gtest-param-util.h gtest-1.7.0/include/gtest/internal/gtest-port.h gtest-1.7.0/include/gtest/internal/gtest-string.h gtest-1.7.0/include/gtest/internal/gtest-tuple.h gtest-1.7.0/include/gtest/internal/gtest-type-util.h

EXTRA_DIST=test_data custom-test-driver $(gtest_src) $(gtest_hdr)
AM_CPPFLAGS = -DTEST_DELTA=1e-9l -DDATDIR=\"$(srcdir)/test_data\" @PTHREAD_CFLAGS@

CPPFLAGS=-I$(srcdir)/gtest-1.7.0/include -I$(srcdir)/gtest-1.7.0 -I$(top_srcdir)



libgtest_la_SOURCES=gtest-1.7.0/src/gtest-all.cc
testTKfit_LDADD=libgtest.la ../libt2toolkit.la @PTHREAD_LIBS@ @PTHREAD_CFLAGS@ $(LIBADD_DLOPEN)
testTempo2_LDADD=libgtest.la ../libtempo2.la @PTHREAD_LIBS@ @PTHREAD_CFLAGS@ $(LIBADD_DLOPEN)
testEndToEnd_LDADD=libgtest.la ../libtempo2.la @PTHREAD_LIBS@ @PTHREAD_CFLAGS@ $(LIBADD_DLOPEN)

