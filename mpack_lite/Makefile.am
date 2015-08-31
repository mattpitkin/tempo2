
if BUILTIN_MPACK
noinst_LTLIBRARIES=libmpack_lite.la

mblaslite_objects=
mlapacklite_objects=
include mpack_lite.double.inc

if HAVE_FLOAT128
include mpack_lite.__float128.inc
endif


libmpack_lite_la_DEPENDENCIES=$(mblaslite_objects) $(mlapacklite_objects)
libmpack_lite_la_LIBADD=$(mblaslite_objects) $(mlapacklite_objects)
libmpack_lite_la_SOURCES=empty.cpp

endif

%.double.lo: %.cpp
	$(LTCXXCOMPILE) -I$(top_srcdir)/mpack_lite/include -D___MPACK_BUILD_WITH_DOUBLE___ -finline-functions -c $^ -o $@


%.__float128.lo: %.cpp
	$(LTCXXCOMPILE) -I$(top_srcdir)/mpack_lite/include -D___MPACK_BUILD_WITH___FLOAT128___ -finline-functions -c $^ -o $@

