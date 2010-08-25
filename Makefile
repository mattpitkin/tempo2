include Makefile.config
include Makefile.common

all: slim plugins

slim: tempo2 libtempo2pred.a

plugins: $(TEMPO2)/plugins/.exists
	cd plugin && $(MAKE)

ifeq ($(origin TEMPO2), undefined)
install sliminstall:
	@echo "ERROR: TEMPO2 environment variable not set!"
	exit 1
else
install: all sliminstall
	cd plugin && $(MAKE) install

sliminstall: dirs
	cp -p libtempo2pred.a libtempo2.a sofa/libsofa.a $(TEMPO2)/lib
	cp -p tempo2.h tempo2pred.h tempo2pred_int.h $(TEMPO2)/include
endif

dirs: $(TEMPO2)/lib/.exists $(TEMPO2)/include/.exists

# $* expands into the implicit rule stem
# $@ expands into the target name
%/.exists:
	mkdir -p $*
	touch $@

clean:
	rm -f $(LIBRARY) tempo2 tempo2.o $(LIB_OBJS)
	cd sofa && $(MAKE) clean
	cd plugin && $(MAKE) clean

sofa/libsofa.a:
	cd sofa && $(MAKE) libsofa.a


LIB_SRCS = preProcess.C preProcessSimple.C calculate_bclt.C T2toolkit.C TKfit.C doFit.C readTimfile.C formBats.C initialise.C shapiro_delay.C textOutput.C displayParameters.C formResiduals.C jpleph.c toa2utc.C tai2tt.C dm_delays.C getInputs.C readEphemeris.C tt2tdb.C utc2tai.C tai2ut1.C  get_obsCoord.C readParfile.C units.C vectorPulsar.C tempo2Util.C polyco.C dynarr.C clkcorr.C observatory.C storePrecision.C readJBO_bat.C MSSmodel.C BTmodel.C BTJmodel.C ELL1model.C DDGRmodel.C DDmodel.C DDKmodel.C DDSmodel.C T2model.C ifteph.C eop.C tropo.C bootstrap.C cheby2d_int.C cheby2d.c tabulatedfunction.C t1polyco.c tempo2pred.c secularMotion.C sw_delay.C global.C

LIB_OBJS = $(patsubst %.c,%.o,${patsubst %.f,%.o,${patsubst %.C,%.o,$(LIB_SRCS)}})

LIBRARY = libtempo2.a

$(LIBRARY): $(LIB_OBJS)
	$(AR) $(ARFLAGS) $(LIBRARY) $(LIB_OBJS)

PRED_LIBRARY = libtempo2pred.a
PRED_LIB_OBJS = tempo2pred.o cheby2d.o t1polyco.o

$(PRED_LIBRARY): $(PRED_LIB_OBJS)
	$(AR) $(ARFLAGS) $(PRED_LIBRARY) $(PRED_LIB_OBJS)


tempo2.o: tempo2.C

tempo2: tempo2.o $(LIBRARY) sofa/libsofa.a
	$(CXX) $(CXXFLAGS) -o tempo2 tempo2.o -L. -ltempo2 -Lsofa -lsofa $(LDFLAGS) $(LF) $(LDSO) -g


