#!/bin/bash

if [[ "$HOSTNAME" == "lindisfarne.jb.man.ac.uk" ]] ; then
ARCH=darwin15

CXX=g++-7 

cmd="$CXX  -o $TEMPO2/plugins/mjkbayes_${ARCH}_plug.t2 mjkbayes_plug.C -ltempo2 -I../.. $LDFLAGS $CFLAGS -shared -fPIC -O2 /Users/mkeith/Code/tempo2/TempoNest/MultiNest/libnest3.a -lgfortran -llapack -lmpi -lmpi_mpifh  -Wfatal-errors"

fi

if [[ "$HOSTNAME" == "hulk.jb.man.ac.uk" ]] ; then

ARCH=$LOGIN_ARCH
cmd="$CXX --std=gnu++11 -o $PSR_EXPERIMENTAL/share/tempo2-dev/mjkbayes_${ARCH}_plug.t2 mjkbayes_plug.C -ltempo2 -I../.. $LDFLAGS $CFLAGS -shared -fPIC -O2 -lmultinest_mpi -lgfortran -lopenblas -lmpi -lmpi_mpifh"

fi


if [[ "$HOSTNAME" == "ubuntu" ]] ; then
    ARCH="linux-gnu"
    CXX=g++
    cmd="$CXX --std=gnu++11 -o $TEMPO2/plugins/mjkbayes_${ARCH}_plug.t2 mjkbayes_plug.C -ltempo2 -I../.. -I/usr/include/mpi $LDFLAGS $CFLAGS -shared -fPIC -O2 -lmultinest_mpi -lgfortran -lmpi -lmpi_mpifh -lmpi_cxx"

fi

echo $cmd
$cmd


