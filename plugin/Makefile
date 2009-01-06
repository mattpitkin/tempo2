include ../Makefile.config
include ../Makefile.common

#Note, if you get runtime linkage problems when you make a new plugin,
# be sure to add calls to the offending functions in 
#tempo2.C:thwart_annoying_dynamic_library_stuff()

# Plugins involving no extra libraries
PLUGIN_SRCS := general2_plug.C general_plug.C matrix_plug.C transform_plug.C fake_plug.C

# Plugins involving PGPLOT
ifeq ($(origin LCPGPLOT), undefined)
WARNING += "Warning, LCPGPLOT not defined in Makefile.config; not building PGPLOT plugins"
else
PLUGIN_SRCS += basic_plug.C delays_plug.C plk_plug.C splk_plug.C  
CXXFLAGS += $(ICPGPLOT)
LDFLAGS += $(LCPGPLOT)
endif

# Plugins that are currently not compiling or working correctly
BROKEN_SRCS = compareTempo_plug.C  spectral_plug.C tnoise_plug.C

PLUGIN_LIBS := $(patsubst %_plug.C,%_$(ARCH)_plug.so,$(PLUGIN_SRCS))

CXXFLAGS += -I..

%_$(ARCH)_plug.so: %_plug.o
	$(DSO) $@ $< $(LDFLAGS) 

all: $(PLUGIN_LIBS)
	@echo $(WARNING)

clean:
	rm -f $(PLUGIN_LIBS)

ifeq ($(origin TEMPO2), undefined)
install:
	@echo "ERROR: You must set env. variable TEMPO2 before trying to install plugins!"
	exit 1
else
install:
	cp $(PLUGIN_LIBS) $(TEMPO2)/plugins
endif
