#!/bin/bash -e
tempo2_src=../../

echo "g++ -fPIC -c ppta_dr1_splug.C -I${PGPLOT_DIR} -I$tempo2_src"
g++ -fPIC -c ppta_dr1_splug.C -I${PGPLOT_DIR} -I$tempo2_src
echo "g++ -fPIC -shared -o ../built_plugins/ppta_dr1_${LOGIN_ARCH}_splug.t2 ppta_dr1_splug.o $LDFLAGS"
g++ -fPIC -shared -o ../built_plugins/ppta_dr1_${LOGIN_ARCH}_splug.t2 ppta_dr1_splug.o $LDFLAGS
rm ppta_dr1_splug.o
