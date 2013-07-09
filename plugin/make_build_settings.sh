#!/bin/bash
list=$1
libextra=$2
incextra=$3

installtargets=""
buildtargets=""

/bin/echo -n "plugin_LTLIBRARIES+="
for lib in `grep -v -e '^#' $list` ; do
	if [ -e $lib ] ; then
		plugname=`/bin/echo $lib | sed -e 's:\(.*\)_\(s\{0,1\}plug\)\..*:\1:'`
		installtargets="$installtargets $plugname-install"
		buildtargets="$buildtargets $plugname"
		/bin/echo -n " "`/bin/echo $lib | sed -e 's:\(.*\)_\(s\{0,1\}plug\)\..*:\1_@T2ARCH@_\2.la:'`
	else
		echo "WARNING: File plugin/$lib does not exist... this plugin will not be built" >&2
	fi
done

/bin/echo ""

/bin/echo "PLUGINSTALLS+=$installtargets"
/bin/echo "PLUGINS+=$buildtargets"

for lib in `grep -v -e '^#' $list` ; do
	if [ -e $lib ] ; then
		/bin/echo $lib | sed -e 's:\(.*\)_\(s\{0,1\}plug\)\..*:\1_@T2ARCH@_\2_la_SOURCES=../tempo2.h '"$lib"':'
		/bin/echo $lib | sed -e 's:\(.*\)_\(s\{0,1\}plug\)\..*:\1_@T2ARCH@_\2_la_LIBADD=../libtempo2.la ../sofa/libsofa.la '"$libextra"':'
		if [ -n "$incextra" ] ; then
			/bin/echo $lib | sed -e 's:\(.*\)_\(s\{0,1\}plug\)\..*:\1_@T2ARCH@_\2_la_CPPFLAGS='"$incextra"':'
		fi
		plugname=`/bin/echo $lib | sed -e 's:\(.*\)_\(s\{0,1\}plug\)\..*:\1:'`
		plugtype=`/bin/echo $lib | sed -e 's:\(.*\)_\(s\{0,1\}plug\)\..*:\2:'`
		/bin/echo "$plugname-install: ${plugname}_@T2ARCH@_${plugtype}.la plugdir"
		/bin/echo "	\$(INSTALL) .libs/${plugname}_@T2ARCH@_${plugtype}.t2 \$(plugindir)"
		/bin/echo "$plugname: ${plugname}_@T2ARCH@_${plugtype}.la"
	fi
done

