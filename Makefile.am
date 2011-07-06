
ACLOCAL_AMFLAGS = -I autoconf

SUBDIRS=sofa . #plugin

bin_PROGRAMS=tempo2
lib_LTLIBRARIES=libtempo2.la libtempo2pred.la

#It'd be better to sort out which headers are neaded by which target, but I'm to busy to do that right now.
headers=GWsim.h T2toolkit.h TKfit.h TKspectrum.h dynarr.h ifteph.h jpl_int.h jpleph.h read_fortran.h read_fortran2.h tabulatedfunction.h tempo2.h tempo2Util.h tempo2pred.h tempo2pred_int.h

libtempo2_la_SOURCES = $(headers) preProcess.C calculate_bclt.C T2toolkit.C TKfit.C doFit.C readTimfile.C formBats.C initialise.C shapiro_delay.C textOutput.C displayParameters.C formResiduals.C jpleph.c toa2utc.C tai2tt.C dm_delays.C getInputs.C readEphemeris.C tt2tdb.C utc2tai.C tai2ut1.C  get_obsCoord.C readParfile.C units.C vectorPulsar.C tempo2Util.C polyco.C dynarr.C clkcorr.C observatory.C storePrecision.C readJBO_bat.C MSSmodel.C BTmodel.C BTJmodel.C ELL1model.C DDGRmodel.C DDmodel.C DDKmodel.C DDSmodel.C T2model.C ifteph.C eop.C tropo.C bootstrap.C cheby2d_int.C cheby2d.c tabulatedfunction.C t1polyco.c tempo2pred.c secularMotion.C sw_delay.C global.C preProcessSimple.C
libtempo2_la_LIBADD = sofa/libsofa.la
libtempo2_la_CFLAGS = -fPIC
#libtempo2_la_LDFLAGS = -static

libtempo2_la_includedir=$(includedir)
libtempo2_la_include_HEADERS = tempo2.h

libtempo2pred_la_SOURCES = $(headers) tempo2pred.c cheby2d.c t1polyco.c
libtempo2pred_la_CFLAGS = -fPIC
#libtempo2pred_la_LDFLAGS = -static

libtempo2pred_la_includedir=$(includedir)
libtempo2pred_la_include_HEADERS = tempo2pred.h tempo2pred_int.h



tempo2_SOURCES = $(headers) tempo2.C 
tempo2_LDADD = libtempo2.la sofa/libsofa.la @FLIBS@ $(LIBADD_DLOPEN)


# Specify to include the T2runtime in distributions
EXTRA_DIST = T2runtime unsupported_plugins

DIST_SUBDIRS=plugin

	

# This clears out any CVS files that might get in our distribution
dist-hook:
	rm -rf `find $(distdir)/T2runtime -name CVS`
	rm -rf $(distdir)/T2runtime/plugins/*`
	rm -rf $(distdir)/unsupported_plugins/built_plugins/*
	rm -rf `find $(distdir)/unsupported_plugins/ -name CVS`
	rm -rf `find $(distdir)/unsupported_plugins/ -name '*.o'`



plugins-install: all
	$(MAKE) -C plugin plugins-install

plugins: all
	$(MAKE) -C plugin


unsupported:
	$(MAKE) -C unsupported_plugins
	
unsupported-install: unsupported
	$(MAKE) -C unsupported_plugins install

.PHONY: plugins-install plugins unsupported_plugins

install-data-hook:
	@echo ""
	@echo ""
	@echo "*** NOTICE ***"
	@echo "Plugins are no longer installed by default"
	@echo "Use: 'make plugins-install' to install them"
	@echo ""
	@echo ""
