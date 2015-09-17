Directory structure                    {#dirs}
===================

The tempo2 directory structure:

~~~~{.f}
.
+-- autoconf.boot
+-- documentation
+-- mpack_lite
+-- plugin
+-- sofa
+-- T2runtime
+-+ tests
  +-- gtest-1.7.0
  +-- test_data
+-- unsupported_plugins
~~~~

### autoconf.boot
This directory contains the .m4 files used by autoconf to build the configure script. It is copied to autoconf/ by the bootstrap script.

### documentation
Includes this documentation

### mpack_lite
Source code for multi-precision lapack/blas. This is a subset of the mplapack package from http://mplapack.sourceforge.net/

### plugin
Source code for plugins

### sofa
Source code for the 3rd party fortran SOFA library.

### T2runtime
This directory contains the runtime files for tempo2, i.e. the contents of this directory should be reached at $TEMPO2
This includes the clock correction files, observatory parameters and earth ephemerdies, etc.

### tests
Source code for the unit tests, and the gtest library. Also contains a number of data files in the test_data subdirectory used by the tests.

### unsupported_plugins
Source code for other plugins that are for whatever reason not part of the main distribution.

