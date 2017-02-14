
ACLOCAL_AMFLAGS = -I autoconf

SUBDIRS=sofa . tests
DIST_SUBDIRS=sofa plugin unsupported_plugins tests

bin_PROGRAMS=tempo2 getPeriod
lib_LTLIBRARIES=libtempo2.la libtempo2pred.la libt2toolkit.la

BUILT_SOURCES=TKcholesky.h TKlongdouble.h enum_str.h enum_str.C


# VERSION 2016.10.0
# c=2
# r=0
# a=0
libtempo2_la_LDFLAGS=-version-info 2:0:0
libt2toolkit_la_LDFLAGS=-version-info 1:1:0
libtempo2pred_la_LDFLAGS=-version-info 1:1:0



# Source files used by plugins that are not part of libtempo2
xtra_dist_sources=choleskyRoutines.C choleskyRoutines.h

FITFUNCS=t2fit_stdFitFuncs.C t2fit_position.C t2fit_fitwaves.C t2fit_glitch.C t2fit_ifunc.C t2fit_dmmodel.C t2fit_dmother.C t2fit_nestlike.C constraints_nestlike.C

#It'd be better to sort out which headers are neaded by which target, but I'm to busy to do that right now.
headers=GWsim.h T2toolkit.h TKfit.h TKrobust.h TKspectrum.h dynarr.h ifteph.h jpl_int.h jpleph.h read_fortran.h read_fortran2.h tabulatedfunction.h tempo2.h tempo2Util.h tempo2pred.h tempo2pred_int.h constraints.h choleskyRoutines.h T2accel.h TKmatrix.h TKsvd.h cholesky.h TKlog.h t2fit.h $(FITFUNCS:.C=.h) $(BUILT_SOURCES) TKlongdouble.float128.h TKlongdouble.ld.h

libtempo2_la_SOURCES = $(headers) preProcess.C calculate_bclt.C T2toolkit.C TKfit.C TKrobust.C doFit.C readTimfile.C formBats.C initialise.C shapiro_delay.C textOutput.C displayParameters.C formResiduals.C jpleph.c toa2utc.C tai2tt.C dm_delays.C getInputs.C readEphemeris.C tt2tdb.C utc2tai.C tai2ut1.C  get_obsCoord.C readParfile.C units.C vectorPulsar.C tempo2Util.C polyco.C dynarr.C clkcorr.C observatory.C storePrecision.C readJBO_bat.C MSSmodel.C BTmodel.C BTJmodel.C BTXmodel.C ELL1model.C DDGRmodel.C DDmodel.C DDKmodel.C DDSmodel.C DDHmodel.C ELL1Hmodel.C T2model.C ifteph.C eop.C tropo.C bootstrap.C cheby2d_int.C cheby2d.c tabulatedfunction.C t1polyco.c tempo2pred.c secularMotion.C sw_delay.C global.C preProcessSimple.C constraints.C cholesky.C choleskyRoutines.C TKspectrum.C T2accel.C TKsvd.C TKmatrix.C sortToAs.C T2-PTAmodel.C GWsim.C readEphemeris_calceph.C TKlog.C choleskyAutomatic.C t2fit.C TKlongdouble.C t2spectrum.C updateEpoch.C $(FITFUNCS)

libt2toolkit_la_SOURCES=T2toolkit.C T2toolkit.h TKfit.C TKrobust.C TKfit.h TKmatrix.C TKmatrix.h TKsvd.C TKsvd.h  TKlog.C TKlog.h T2accel.C T2accel.h cholesky.C cholesky.h choleskyRoutines.C t2toolkit_global.C TKlongdouble.C TKspectrum.C
libt2toolkit_la_CFLAGS=-fPIC
libt2toolkit_la_LIBADD=

libtempo2_la_LIBADD = sofa/libsofa.la
libtempo2_la_CFLAGS = -fPIC 
#libtempo2_la_LDFLAGS = -static

if HAVE_LAPACK
libtempo2_la_LIBADD+=@LIBS_LAPACK@
libt2toolkit_la_LIBADD+=@LIBS_LAPACK@
endif

if HAVE_BLAS
libtempo2_la_LIBADD+=@BLAS_LIBS@
libt2toolkit_la_LIBADD+=@BLAS_LIBS@
endif

if HAVE_CALCEPH
libtempo2_la_LIBADD+=@CALCEPH_LIBS@
endif

libtempo2_la_includedir=$(includedir)
libtempo2_la_include_HEADERS = tempo2.h GWsim.h t2fit.h

libt2toolkit_la_includedir=$(includedir)
libt2toolkit_la_include_HEADERS = TKfit.h TKrobust.h TKmatrix.h T2toolkit.h TKsvd.h TKcholesky.h TKlog.h TKlongdouble.h

TKcholesky.h: cholesky.h
	cp $(srcdir)/cholesky.h TKcholesky.h

libtempo2pred_la_SOURCES = $(headers) tempo2pred.c cheby2d.c t1polyco.c
libtempo2pred_la_CFLAGS = -fPIC
#libtempo2pred_la_LDFLAGS = -static
libtempo2pred_la_LIBADD =

libtempo2pred_la_includedir=$(includedir)
libtempo2pred_la_include_HEADERS = tempo2pred.h tempo2pred_int.h


tempo2_SOURCES = $(headers) tempo2.C 
tempo2_LDADD = libtempo2.la sofa/libsofa.la @FLIBS@ $(LIBADD_DLOPEN)
tempo2_CPPFLAGS=

if HAVE_QDINSTALL
tempo2_LDADD+=@QDINSTALL_LIBS@
tempo2_CPPFLAGS+=@QDINSTALL_CFLAGS@
endif



getPeriod_SOURCES= getPeriod.C
getPeriod_LDADD = libtempo2pred.la


# Specify to include the T2runtime in distributions
EXTRA_DIST = T2runtime $(xtra_dist_sources) documentation extract_enums.sh


# This clears out any CVS files that might get in our distribution
dist-hook:
	rm -rf `find $(distdir)/T2runtime -name CVS`
	rm -rf $(distdir)/T2runtime/plugins/*
	rm -rf $(distdir)/unsupported_plugins/built_plugins/*
	rm -rf `find $(distdir)/unsupported_plugins/ -name CVS`
	rm -rf `find $(distdir)/unsupported_plugins/ -name '*.o'`
	rm -rf $(distdir)/documentation/doxygen



plugins: all
	$(MAKE) -C plugin

plugins-clean:
	$(MAKE) -C plugin clean

plugins-install: plugins
	mkdir -p $(TEMPO2)/plugins
	$(MAKE) -C plugin plugins-install


unsupported: all
	$(MAKE) -C unsupported_plugins

unsupported-install: unsupported
	mkdir -p $(TEMPO2)/plugins
	$(MAKE) -C unsupported_plugins install

unsupported-clean:
	$(MAKE) -C unsupported_plugins clean

temponest: all
	$(MAKE) -C temponest

temponest-install: temponest
	$(MAKE) -C temponest install

temponest-clean:
	$(MAKE) -C temponest clean

core: all
core-install: install
core-clean: clean

complete: core unsupported plugins
	@echo "Built core, plugins, unsupported"
complete-install: core-install unsupported-install plugins-install
	@echo "Installed tempo2, plugins, unsupported"
complete-clean: plugins-clean unsupported-clean core-clean
	@echo "Clean core, plugins, unsupported"

#check: all
#	$(MAKE) -C tests check

if HAVE_FLOAT128
TKlongdouble.h: TKlongdouble.float128.h
	cp $(srcdir)/TKlongdouble.float128.h $@
else
TKlongdouble.h: TKlongdouble.ld.h
	cp $(srcdir)/TKlongdouble.ld.h $@
endif

if HAVE_DOXYGEN
doc: $(libtempo2_la_SOURCES)
	doxygen
	make -C documentation/doxygen/latex
	cp documentation/doxygen/latex/refman.pdf documentation/tempo2_documentation.pdf
	
else
doc:
	echo "Sorry, you need doxygen to build the documentation"
endif

enum_str.h: enum_str.C

enum_str.C: extract_enums.sh tempo2.h
	$(srcdir)/extract_enums.sh < $(srcdir)/tempo2.h


CLEANFILES=$(BUILT_SOURCES) documentation/doxygen

.PHONY: plugins-install plugins plugins-clean \
   unsupported unsupported-install unsupported-clean \
   complete complete-install complete-clean \
   temponest temponest-install \
   core core-install core-clean test doc


install-data-hook:
	@echo ""
	@echo ""
	@echo "*** NOTICE ***"
	@echo "Plugins are no longer installed by default"
	@echo "Use: 'make plugins-install' to install them"
	@echo "Use: 'make temponest-install' to install temponest"
	@echo ""
	@echo ""
