#!/bin/bash

plug=$1

cmd="$CXX -o $TEMPO2/plugins/${plug}_darwin13_plug.t2 ${plug}_plug.C -ltempo2 -lpgplot -lcpgplot -I.. -shared -fPIC -O2" 

echo $cmd
$cmd
