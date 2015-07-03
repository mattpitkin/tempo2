plugindir=@TEMPO2_PLUG_INSTALL@

plugin_LTLIBRARIES=addDmVar_@T2ARCH@_plug.la addGWB_@T2ARCH@_plug.la addGaussian_@T2ARCH@_plug.la add2state_@T2ARCH@_plug.la formIdeal_@T2ARCH@_plug.la createRealisation_@T2ARCH@_plug.la addRedNoise_@T2ARCH@_plug.la addCGW_@T2ARCH@_plug.la addInstrumentalIssues_@T2ARCH@_plug.la addOutliers_@T2ARCH@_plug.la

bin_PROGRAMS=read_toasim
lib_LTLIBRARIES=libtoasim.la


###
#
## This sets the includes and options that are needed to build nice modules with libtool
#
####
AM_CPPFLAGS=-I@top_srcdir@/

CFLAGS+=@FFTW_CFLAGS@

pluginINSTALLS=addDmVar-install addGWB-install addGaussian-install add2state-install addInstrumentalIssues-install formIdeal-install createRealisation-install addRedNoise-install addCGW-install addOutliers-install


TOASIM_src=toasim.h toasim.c toasim_compute.c toasim_io.c makeRedNoise.C makeRedNoise.h
libtoasim_la_SOURCES=$(TOASIM_src)

read_toasim_SOURCES=read_toasim.c
read_toasim_LDADD=libtoasim.la -lfftw3f ../../libtempo2.la

plugin_LDFLAGS= -shared -module -avoid-version -shrext .t2 -lfftw3f
addDmVar_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h addDmVar_plug.C $(TOASIM_src)
addDmVar_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
addDmVar_@T2ARCH@_plug_la_LDFLAGS=$(plugin_LDFLAGS)

addDmVar-install: addDmVar_@T2ARCH@_plug.la
	$(INSTALL) .libs/addDmVar_@T2ARCH@_plug.t2 $(plugindir)

addRedNoise_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h addRedNoise_plug.C $(TOASIM_src)
addRedNoise_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
addRedNoise_@T2ARCH@_plug_la_LDFLAGS=$(plugin_LDFLAGS)
addRedNoise-install: addRedNoise_@T2ARCH@_plug.la
	$(INSTALL) .libs/addRedNoise_@T2ARCH@_plug.t2 $(plugindir)


addGWB_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h addGWB_plug.C $(TOASIM_src)
addGWB_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
addGWB_@T2ARCH@_plug_la_LDFLAGS=$(plugin_LDFLAGS)
addGWB-install: addGWB_@T2ARCH@_plug.la
	$(INSTALL) .libs/addGWB_@T2ARCH@_plug.t2 $(plugindir)

addGaussian_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h addGaussian_plug.C $(TOASIM_src)
addGaussian_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
addGaussian_@T2ARCH@_plug_la_LDFLAGS=$(plugin_LDFLAGS)
addGaussian-install: addGaussian_@T2ARCH@_plug.la
	$(INSTALL) .libs/addGaussian_@T2ARCH@_plug.t2 $(plugindir)

addOutliers_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h addOutliers_plug.C $(TOASIM_src)
addOutliers_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
addOutliers_@T2ARCH@_plug_la_LDFLAGS=$(plugin_LDFLAGS)
addOutliers-install: addOutliers_@T2ARCH@_plug.la
	$(INSTALL) .libs/addOutliers_@T2ARCH@_plug.t2 $(plugindir)



add2state_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h add2state_plug.C $(TOASIM_src)
add2state_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
add2state_@T2ARCH@_plug_la_LDFLAGS=$(plugin_LDFLAGS)
add2state-install: add2state_@T2ARCH@_plug.la
	$(INSTALL) .libs/add2state_@T2ARCH@_plug.t2 $(plugindir)

addInstrumentalIssues_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h addInstrumentalIssues_plug.C $(TOASIM_src)
addInstrumentalIssues_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
addInstrumentalIssues_@T2ARCH@_plug_la_LDFLAGS=$(plugin_LDFLAGS)
addInstrumentalIssues-install: addInstrumentalIssues_@T2ARCH@_plug.la
	$(INSTALL) .libs/addInstrumentalIssues_@T2ARCH@_plug.t2 $(plugindir)

addCGW_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h addCGW_plug.C $(TOASIM_src)
addCGW_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
addCGW_@T2ARCH@_plug_la_LDFLAGS=$(plugin_LDFLAGS)
addCGW-install: addCGW_@T2ARCH@_plug.la
	$(INSTALL) .libs/addCGW_@T2ARCH@_plug.t2 $(plugindir)


formIdeal_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h formIdeal_plug.C $(TOASIM_src)
formIdeal_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
formIdeal_@T2ARCH@_plug_la_LDFLAGS=$(plugin_LDFLAGS)
formIdeal-install: formIdeal_@T2ARCH@_plug.la
	$(INSTALL) .libs/formIdeal_@T2ARCH@_plug.t2 $(plugindir)

createRealisation_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h createRealisation_plug.C $(TOASIM_src)
createRealisation_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
createRealisation_@T2ARCH@_plug_la_LDFLAGS=$(plugin_LDFLAGS)
createRealisation-install: createRealisation_@T2ARCH@_plug.la
	$(INSTALL) .libs/createRealisation_@T2ARCH@_plug.t2 $(plugindir)



install: $(pluginINSTALLS) .libs/libtoasim.la read_toasim
	$(INSTALL) .libs/libtoasim.la $(libdir)
	$(INSTALL) .libs/read_toasim $(bindir)
