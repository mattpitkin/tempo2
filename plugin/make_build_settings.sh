#!/bin/bash
list=$1
libextra=$2
incextra=$3


/bin/echo -n "plugin_LTLIBRARIES+="
for lib in `grep -v -e '^#' $list` ; do
	if [ -e $lib ] ; then
		/bin/echo -n " "`/bin/echo $lib | sed -e 's:\(.*\)_\(s\{,1\}plug\)\..*:\1_@T2ARCH@_\2.la:'`
	else
		echo "WARNING: File plugin/$lib does not exist... cannot build plugin" >&2
	fi
done

/bin/echo ""

for lib in `grep -v -e '^#' $list` ; do
	if [ -e $lib ] ; then
		/bin/echo $lib | sed -e 's:\(.*\)_\(s\{,1\}plug\)\..*:\1_@T2ARCH@_\2_la_SOURCES=../tempo2.h '$lib':'
		if [ -n "$libextra" ] ; then
			/bin/echo $lib | sed -e 's:\(.*\)_\(s\{,1\}plug\)\..*:\1_@T2ARCH@_\2_la_LIBADD='$libextra':'
		fi
		if [ -n "$incextra" ] ; then
			/bin/echo $lib | sed -e 's:\(.*\)_\(s\{,1\}plug\)\..*:\1_@T2ARCH@_\2_la_INCLUDES='$incextra':'
		fi
	fi
done
