#!/bin/sh

flt=$1
dir=`dirname $0`
cd $dir
if test -z "$flt" ; then
    exit 1
fi
while read alg ; do
    oo=`find mblas-lite/optimized/$flt/ -name "$alg*.cpp" | sed -e "s:.cpp:.${flt}.lo:g"`
    if test -z "$opt" ; then
        oo=`find mblas-lite/reference/ -name "$alg*.cpp" | sed -e "s:.cpp:.${flt}.lo:g"`
    fi
    objects="$objects $oo"
    objects=`echo $objects`

done < mblas-routines

echo "# ** $flt **"

echo "mblaslite_objects+=$objects"
