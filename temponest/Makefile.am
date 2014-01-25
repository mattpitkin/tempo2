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

# Only build TempoNest if MultiNest is installed (other samplers in the future?)
if HAVE_MULTINEST

include temponest.am

endif


.PHONY: $(PLUGINS) $(PLUGINSTALLS) temponest temponest-install

temponest: $(PLUGINS)
temponest-install: plugdir $(PLUGINSTALLS)


plugdir:
	mkdir -p $(plugindir)
