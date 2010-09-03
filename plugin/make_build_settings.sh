#!/bin/bash
list=$1
libextra=$2
incextra=$3


/bin/echo -n "plugin_LTLIBRARIES+="
for lib in `grep -v -e '^#' $list` ; do
	if [ -e $lib ] ; then
		/bin/echo -n " "`/bin/echo $lib | sed -e 's:\(.*\)_\(s\{0,1\}plug\)\..*:\1_@T2ARCH@_\2.la:'`
	else
		echo "WARNING: File plugin/$lib does not exist... this plugin will not be built" >&2
	fi
done

/bin/echo ""

for lib in `grep -v -e '^#' $list` ; do
	if [ -e $lib ] ; then
		/bin/echo $lib | sed -e 's:\(.*\)_\(s\{0,1\}plug\)\..*:\1_@T2ARCH@_\2_la_SOURCES=../tempo2.h '"$lib"':'
		/bin/echo $lib | sed -e 's:\(.*\)_\(s\{0,1\}plug\)\..*:\1_@T2ARCH@_\2_la_LIBADD=../libtempo2.la ../sofa/libsofa.la '"$libextra"':'
		if [ -n "$incextra" ] ; then
			/bin/echo $lib | sed -e 's:\(.*\)_\(s\{0,1\}plug\)\..*:\1_@T2ARCH@_\2_la_INCLUDES='"$incextra"':'
		fi
	fi
done

