
if BUILTIN_MPACK
noinst_LTLIBRARIES=libmpack_lite.la

mblaslite_objects=
mlapacklite_objects=mlapack-lite/reference/Rlasdt.double.lo
#mlapack-lite/reference/iMlaenv.double.lo
include mpack_lite.double.inc

if HAVE_FLOAT128
include mpack_lite.__float128.inc
endif


libmpack_lite_la_DEPENDENCIES=$(mblaslite_objects) $(mlapacklite_objects)
libmpack_lite_la_LIBADD=$(mblaslite_objects) $(mlapacklite_objects)
libmpack_lite_la_SOURCES=empty.cpp

LTFILES=$(mblaslite_objects) $(mlapacklite_objects)
CLEANFILES=$(mlapacklite_objects) $(mblaslite_objects) $(LTFILES:.lo=.o)

MPACK_BUILTIN_CPPFLAGS=-Wno-fatal-errors -I$(top_srcdir)/mpack_lite/include -O3
MPACK_LT=$(LIBTOOL) $(AM_V_lt) --tag=CXX $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(CXX) $(MPACK_BUILTIN_CPPFLAGS)

endif

%.double.lo: %.cpp
	$(MPACK_LT) -D___MPACK_BUILD_WITH_DOUBLE___ -finline-functions -c $^ -o $@


%.__float128.lo: %.cpp
	$(MPACK_LT) -D___MPACK_BUILD_WITH___FLOAT128___ -finline-functions -c $^ -o $@


