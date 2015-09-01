#!/bin/sh
objects=""
flt=$1
dir=`dirname $0`
cd $dir
if test -z "$flt" ; then
    exit 1
fi
while read alg ; do
for pack in mblas mlapack ; do
    oo=""
    oo=`find ${pack}-lite/optimized/$flt/ -name "$alg*.cpp" | sed -e "s:.cpp:.${flt}.lo:g" | grep -v Rgemm_omp`
    if test -z "$oo" ; then
        oo=`find ${pack}-lite/reference/ -name "$alg*.cpp" | sed -e "s:.cpp:.${flt}.lo:g"`
    fi
    if test -n "$oo" ; then
        oo=`echo $oo`
        objects="$objects $oo"
        break
    fi
done
done < mpack-routines

echo "# ** $flt **"

echo "mblaslite_objects+=$objects"
