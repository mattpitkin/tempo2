plugindir=@TEMPO2_PLUG_INSTALL@

plugin_LTLIBRARIES=addDmVar_@T2ARCH@_plug.la addGWB_@T2ARCH@_plug.la addGaussian_@T2ARCH@_plug.la formIdeal_@T2ARCH@_plug.la createRealisation_@T2ARCH@_plug.la addRedNoise_@T2ARCH@_plug.la


###
#
## This sets the includes and options that are needed to build nice modules with libtool
#
####
LDFLAGS+= -shared -module -avoid-version -shrext .t2 -lfftw3f
INCLUDES=-I@top_srcdir@/

CFLAGS+=@FFTW_CFLAGS@

pluginINSTALLS=addDmVar-install addGWB-install addGaussian-install formIdeal-install createRealisation-install addRedNoise-install

TOASIM_src=toasim.h toasim.c toasim_compute.c toasim_io.c makeRedNoise.C makeRedNoise.h


addDmVar_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h addDmVar_plug.C $(TOASIM_src)
addDmVar_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
addDmVar-install: addDmVar_@T2ARCH@_plug.la
	$(INSTALL) .libs/addDmVar_@T2ARCH@_plug.t2 $(plugindir)

addRedNoise_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h addRedNoise_plug.C $(TOASIM_src)
addRedNoise_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
addRedNoise-install: addRedNoise_@T2ARCH@_plug.la
	$(INSTALL) .libs/addRedNoise_@T2ARCH@_plug.t2 $(plugindir)


addGWB_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h addGWB_plug.C $(TOASIM_src)
addGWB_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
addGWB-install: addGWB_@T2ARCH@_plug.la
	$(INSTALL) .libs/addGWB_@T2ARCH@_plug.t2 $(plugindir)

addGaussian_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h addGaussian_plug.C $(TOASIM_src)
addGaussian_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
addGaussian-install: addGaussian_@T2ARCH@_plug.la
	$(INSTALL) .libs/addGaussian_@T2ARCH@_plug.t2 $(plugindir)


formIdeal_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h formIdeal_plug.C $(TOASIM_src)
formIdeal_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
formIdeal-install: formIdeal_@T2ARCH@_plug.la
	$(INSTALL) .libs/formIdeal_@T2ARCH@_plug.t2 $(plugindir)

createRealisation_@T2ARCH@_plug_la_SOURCES=@top_srcdir@/tempo2.h createRealisation_plug.C $(TOASIM_src)
createRealisation_@T2ARCH@_plug_la_LIBADD=../../libtempo2.la ../../sofa/libsofa.la
createRealisation-install: createRealisation_@T2ARCH@_plug.la
	$(INSTALL) .libs/createRealisation_@T2ARCH@_plug.t2 $(plugindir)



install: $(pluginINSTALLS)
