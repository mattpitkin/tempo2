#!/bin/bash
tempo2_src=../../

echo g++ -fPIC -c glast_plug.C -I${PGPLOT_DIR}  -I$tempo2_src
g++ -fPIC -c glast_plug.C -I${PGPLOT_DIR} -I$tempo2_src
echo g++ -fPIC -shared -o ../built_plugins/glast_${LOGIN_ARCH}_plug.t2 glast_plug.o $LDFLAGS -lcpgplot -lpgplot -L/usr/X11R6/lib -lX11 -lpng -lz -lgfortran
g++ -fPIC -shared -o ../built_plugins/glast_${LOGIN_ARCH}_plug.t2 glast_plug.o $LDFLAGS -lcpgplot -lpgplot -L/usr/X11R6/lib -lX11 -lpng -lz -lgfortran
rm glast_plug.o
