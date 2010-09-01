plugindir=@TEMPO2_DIR@/plugins
plugin_LTLIBRARIES=
###
#
# This sets the includes and options that are needed to build nice modules with libtool
#
###
LDFLAGS+= -shared -module -avoid-version -shrext .t2
LDADD=../libtempo2.la
INCLUDES=-I@top_srcdir@/

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
endif

AM_CPPFLAGS = @PGPLOT_CFLAGS@

#
# GSL plugins
#
if HAVE_GSL
include build_settings/gsl.am
endif


