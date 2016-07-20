[![Build Status](https://drone.io/bitbucket.org/psrsoft/tempo2/status.png)](https://drone.io/bitbucket.org/psrsoft/tempo2/latest)

Git INSTALLATION README
=======================

0. Contents
	1. What this package is    
	2. Quick Guide    
	3. Requirements    
	4. Detailed instalation guide    
	5. Plugins    
	6. Changes from old makefile    
	7. Installation troubleshooting    
	8. Bugs and feature requests


1. What this package is
------------------------
You (or someone else) have checked out tempo2 from the Git
(https://bitbucket.org/psrsoft/tempo2)

This is the best way to get the latest/cutting edge version, and develop your
own additions to the tempo2 code or via plugins.

For more information on tempo2 see:
http://www.atnf.csiro.au/research/pulsar/tempo2/

This requires the gnu autotools. If you don't have or don't want to install 
autotools, we recommend you install the latest distributed release from
http://www.atnf.csiro.au/research/pulsar/tempo2/
or use PSRSOFT to install tempo2:
http://www.pulsarastronomy.net/wiki/Software/PSRSoft

2. Quick Guide
---------------
Bootstrap the build system:

    ./bootstrap

setup the tempo2 runtime dir

    cp -r T2runtime /usr/share/tempo2/
    export TEMPO2=/usr/share/tempo2/

Configure:

    ./configure [[--prefix=/your/install/path]]

use --prefix to set the path you want to install the binaries and libraries

Make and install...

    make && make install

You will probably want to build the default plugins (plk, etc). Do this with:

    make plugins && make plugins-install


And you're done.

3. Requirements
---------------
Tempo2 requires the following:

 - A fortran 77 compiler (tested with gfortran).
 - A C compiler (tested with gcc).

Plugins may have other requirements, notably PGPLOT.

5. Plugins
----------
The bootstrap command will create suitible makefiles for the default set of
plugins. This is controled by the contents of the files in
./plugin/plugin_lists/

 - vanilla.plugins lists plugins to install which have no dependancies.
 - pgplot.plugins lists plugins to install that are dependant on PGPLOT.
 - gsl.plugins lists plugins to install that are dependant on the GSL.

5.1 Building your own plugin
----------------------------
The easiest way to compile your own plugins is:

    g++ {$CFLAGS} {$LDFLAGS} -fPIC -shared -o {$TEMPO2}/plugins/{$PLG_NAME}_{$LOGIN_ARCH}_plug.t2 {$SRCLIST}

where:

 - `{$PLG_NAME}` is the name of your plugin
 - `{$SRCLIST}` is your plugin's source code.
 - `{$LOGIN_ARCH}` is the result of `` `uname` `` (usualy Linux).
 - `{$CFLAGS}` are the compiler flags your plugin needs... remeber to add a -I option to point to the location of tempo2.h
 - `{$LDFLAGS}` are any linking options you need, e.g. pgplot, etc.
 - `{$TEMPO2}` is the tempo2 runtime dir

For example, to compile a basic plugin called '`foo`' on linux, you might do

    g++ -I/usr/src/tempo2 -fPIC -shared -o $TEMPO2/plugins/foo_{$LOGIN_ARCH}_plug.t2 foo_plug.C


5.2 Adding a new plugin to the default build list
-------------------------------------------------
If your plugin has dependances that are already covered by the lists above,
just add the name to the appropriate list, and name your plugin source file as:

name_plug.C

6. Changes from the old Make system.
------------------------------------
At the start of 2010, tempo2 moved over to an autotools based make system,
replacing the old hand written makefiles. This may confuse some people!

Important notes:

 - Tempo2 plugins now have a .t2 extention, rather than the old .so
    This is to ensure reduce confusion on MacOSx and to allow the old
    make system and the new make system to co-exist for a while.
 - Any 3rd party plugins will still work as before. Indeed, to update
    a plugin, just change the .so extention to a .t2 extention.
    e.g. mv general_Linux_plug.so general_Linux_plug.t2

7. Installation Troubleshooting
-------------------------------

7.1 Can't find PGPLOT
---------------------
Download pgplot from:
http://www.astro.caltech.edu/~tjp/pgplot/

Or use PSRSOFT to manage the installation.
http://www.pulsarastronomy.net/wiki/Software/PSRSoft

If you have pgplot installed, but it is not detected by the configure script, check:

 - You have got at least libpgplot.a and libcpgplot.a in your LDFLAGS
 - Check you have `$PGPLOT_DIR` pointing to the folder with grfont.dat and rgb.txt
 - Check that you have `$F77` set to the same compiler that compiled PGPLOT
    (e.g. setenv F77 gfortran, if you used gfortran for PGPLOT)


7.2 Incompatible C and Fortran compilers
----------------------------------------
Check that you are using the same build of gcc and gfortran (or whatever compiler you are using).

Note that on MacOSX there is often an issue where the default compiler is incompatible with gfortran.
The gfortran compatible version is often called gcc-4 and gxx-4 or similar. Use this with:

    export CC=gcc-4
    export CXX=g++-4

and reconfigure.


8. Bugs and feature requests
-----------------------------
Please submit bug reports here: https://bitbucket.org/psrsoft/tempo2/issues/new

Note that it is very helpful if you can upload a small example demonstrating the bug!