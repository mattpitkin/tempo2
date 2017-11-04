#!/bin/bash

ARCH=linux-gnu

CXX=g++

cmd="$CXX -o $TEMPO2/plugins/mjkbayes_${ARCH}_plug.t2 mjkbayes_plug.C -ltempo2 -I../.. $LDFLAGS $CFLAGS -shared -fPIC -O2 -lmultinest_mpi"

echo $cmd
$cmd
