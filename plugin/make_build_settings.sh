#!/bin/sh
list=$1
libextra=$2
incextra=$3

/bin/echo -n "plugin_LTLIBRARIES+="
for lib in `grep -v -e '^#' $list` ; do
	/bin/echo -n " "`/bin/echo $lib | sed -e 's:\(.*\)_plug\..*:\1_@T2ARCH@_plug.la:'`
done

/bin/echo ""

for lib in `grep -v -e '^#' $list` ; do
	/bin/echo $lib | sed -e 's:\(.*\)_plug\..*:\1_@T2ARCH@_plug_la_SOURCES=../tempo2.h '$lib':'
	if [ -n "$libextra" ] ; then
		/bin/echo $lib | sed -e 's:\(.*\)_plug\..*:\1_@T2ARCH@_plug_la_LIBADD='$libextra':'
	fi
	if [ -n "$incextra" ] ; then
		/bin/echo $lib | sed -e 's:\(.*\)_plug\..*:\1_@T2ARCH@_plug_la_INCLUDES='$incextra':'
	fi
done
