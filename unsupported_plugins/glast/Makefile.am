plugindir=@TEMPO2_PLUG_INSTALL@

plugin_LTLIBRARIES=glast_@T2ARCH@_plug.la


###
#
## This sets the includes and options that are needed to build nice modules with libtool
#
####
LDFLAGS+= -shared -module -avoid-version -shrext .t2
AM_CPPFLAGS=-I@top_srcdir@/


pluginINSTALLS=glast-install



glast_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h glast_plug.C 
glast_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
glast-install: glast_@T2ARCH@_plug.la
	$(INSTALL) .libs/glast_@T2ARCH@_plug.t2 $(plugindir)


install: $(pluginINSTALLS)
