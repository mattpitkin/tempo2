#!/bin/bash -e

tempo2_src=../..
CXXFLAGS='-O3 -march=native'

for i in toasim*.c ; do
	echo $CC -c $i
	$CC -c $i
done
echo $CXX -fPIC -c $tempo2_src/T2toolkit.C $tempo2_src/TKfit.C $CXXFLAGS
$CXX -fPIC -c $tempo2_src/T2toolkit.C $tempo2_src/TKfit.C $CXXFLAGS
for plug in *_plug.C; do
	name=`basename $plug _plug.C`
	echo $CXX -fPIC -c ${name}_plug.C  -I$tempo2_src $CXXFLAGS
	$CXX -fPIC -c ${name}_plug.C  -I$tempo2_src -O3 -march=native
	$CXX -shared -o ${name}_${LOGIN_ARCH}_plug.t2 ${name}_plug.o T2toolkit.o TKfit.o toasim*.o -ltempo2 -lsofa $LDFLAGS 
	echo $CXX -shared -o $PSRHOME/people/mkeith/t2plugin/${name}_${LOGIN_ARCH}_plug.t2 ${name}_plug.o T2toolkit.o TKfit.o toasim*.o -ltempo2 -lsofa $LDFLAGS; rm ${name}_plug.o
done
echo mv *${LOGIN_ARCH}_plug.t2 ../built_plugins/
mv *${LOGIN_ARCH}_plug.t2 ../built_plugins/
