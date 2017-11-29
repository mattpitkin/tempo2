#!/bin/bash

#ARCH=linux-gnu
ARCH=darwin15

CXX=g++-7

cmd="$CXX -o $TEMPO2/plugins/mjkbayes_${ARCH}_plug.t2 mjkbayes_plug.C -ltempo2 -I../.. $LDFLAGS $CFLAGS -shared -fPIC -O2 /Users/mkeith/Code/tempo2/TempoNest/MultiNest/libnest3.a -lgfortran -llapack -lmpi -lmpi_mpifh"

echo $cmd
$cmd
