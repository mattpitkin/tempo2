#!/bin/bash


cmd="$CXX -o $TEMPO2/plugins/mjk_darwin13_plug.t2 mjk_plug.C -ltempo2 -lpgplot -lcpgplot -I.. -shared -fPIC -O2" 

echo $cmd
$cmd
