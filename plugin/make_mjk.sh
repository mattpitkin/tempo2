#!/bin/bash


cmd="$CXX -o $TEMPO2/plugins/mjk_darwin15_plug.t2 mjk2_plug.C $CFLAGS $LDFLAGS -ltempo2 -I.. -shared -fPIC -O2" 

echo $cmd
$cmd
