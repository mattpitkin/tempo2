plugindir=@TEMPO2_PLUG_INSTALL@

plugin_LTLIBRARIES=ppta_dr1_@T2ARCH@_splug.la


###
#
## This sets the includes and options that are needed to build nice modules with libtool
#
####
LDFLAGS+= -shared -module -avoid-version -shrext .t2
AM_CPPFLAGS=-I@top_srcdir@/


pluginINSTALLS=ppta_dr1-install



ppta_dr1_@T2ARCH@_splug_la_SOURCES=@top_srcdir@/tempo2.h ppta_dr1_splug.C 
ppta_dr1_@T2ARCH@_splug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
ppta_dr1-install: ppta_dr1_@T2ARCH@_splug.la
	$(INSTALL) .libs/ppta_dr1_@T2ARCH@_splug.t2 $(plugindir)


install: $(pluginINSTALLS)
