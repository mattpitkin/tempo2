plugindir=@TEMPO2_PLUG_INSTALL@
###
#
# This sets the includes and options that are needed to build nice modules with libtool
#
###
LDFLAGS+= -shared -module -avoid-version -shrext .t2
AM_CPPFLAGS=-I@top_srcdir@/
# AM_CPPFLAGS += @PGPLOT_CFLAGS@ @FFTW_CFLAGS@ @CFITSIO_CFLAGS@

##
# These values are set by the includes below, so don't change them from being blank!
##
PLUGINSTALLS=
PLUGINS=
plugin_LTLIBRARIES=

#
# Plugins that need no extra libs!
#
#
include build_settings/vanilla.am

#
# Now for plugins that need PGPLOT
#
if HAVE_PGPLOT

include build_settings/pgplot.am

if HAVE_FFTW3
include build_settings/pgplot_fftw.am
endif

if HAVE_CFITSIO
include build_settings/pgplot_cfitsio.am
endif

#END OF PGPLOT
endif


if HAVE_CFITSIO
include build_settings/cfitsio.am
endif

#
# GSL plugins
#
if HAVE_GSL
include build_settings/gsl.am
endif


.PHONY: $(PLUGINS) $(PLUGINSTALLS) plugins plugins-install

plugins: $(PLUGINS)
plugins-install: plugdir $(PLUGINSTALLS)


plugdir:
	mkdir -p $(plugindir)
