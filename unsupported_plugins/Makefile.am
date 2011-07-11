SUBD=toasim glast

.PHONY: $(SUBD) subdirs

subdirs: $(SUBD)

$(SUBD):
	mkdir -p $@/../built_plugins
	@echo "BUILD $@"
	@if test -e $@/Makefile ; then $(MAKE) -C $@ ; else if test -e $@/build.sh ; then cd $@ ; ./build.sh ; else echo "Don't know how to build $@" ; fi ; fi


install:
	$(INSTALL) built_plugins/*.t2 @TEMPO2_PLUG_INSTALL@


EXTRA_DIST=glast toasim

