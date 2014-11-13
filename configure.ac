#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([Tempo2],[2014.11.1],[george.hobbs@csiro.au])

AC_ARG_VAR([TEMPO2], [Default TEMPO2 installation directory])
AC_PREFIX_DEFAULT([${TEMPO2:-"/usr/local"}])

if test -n "${CXXFLAGS}"; then
	user_set_cxxflags=yes
fi

AC_CONFIG_MACRO_DIR([autoconf])
AC_CONFIG_AUX_DIR([config])

# Enable large-file support
AC_SYS_LARGEFILE

# Determine the machine endian
AC_C_BIGENDIAN

# Checks for compilers
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_CC
AC_PROG_CPP
AM_PROG_CC_C_O
AC_PROG_F77
AC_PROG_F77_C_O
AC_F77_LIBRARY_LDFLAGS
AC_PROG_FC

if test -n "${CXXFLAGS}"; then
	user_set_cxxflags=yes
fi
AC_PROG_CXX
if test X$user_set_cxxflags != Xyes; then
	CXXFLAGS="-O2"
fi


AC_PROG_INSTALL
AC_PROG_LIBTOOL



if test $can_build_shared == no ; then
	AC_WARN("Libtool cannot build shared libraries - ")
fi

# Check whether the 'cos' function can be found in library 'm'
# (standard C math library).
# Math library is required for the tempo2 predictor library
AC_CHECK_LIB(m, cos)

# look for third-party packages in various places (or not)
SWIN_LOCAL

AC_MSG_CHECKING([for \$TEMPO2])
if test -n "$TEMPO2" ; then
	TEMPO2_DIR=$TEMPO2
	AC_SUBST(TEMPO2_DIR)
	have_tempo2_dir=yes
else
	have_tempo2_dir=no
fi

AC_MSG_RESULT($have_tempo2_dir)

if test X$have_tempo2_dir != Xyes; then
	AC_WARN("You must set the $TEMPO2 variable before configuring.")
fi

AC_ARG_WITH([tempo2-plug-dir],
              AC_HELP_STRING([--with-tempo2-plug-dir=DIR],
                             [Specify DIR to install plugins]))


AC_MSG_CHECKING([for plugin install path])
if test -z "$with_tempo2_plug_dir" ; then
	TEMPO2_PLUG_INSTALL=$TEMPO2/plugins/
else
	TEMPO2_PLUG_INSTALL=$with_tempo2_plug_dir
	AC_DEFINE_UNQUOTED([TEMPO2_CONFIGURE_PLUG], "$with_tempo2_plug_dir",["Tempo2 plugin path"])
fi

AC_SUBST(TEMPO2_PLUG_INSTALL)
AC_MSG_RESULT($TEMPO2_PLUG_INSTALL)

MJK_FORTRAN_TEST

if test X$fortran_c_links != Xyes ; then
	AC_MSG_ERROR(Your C and F77 compilers do not cross link (\$CC=$CC \$F77=$F77))
fi

LT_LIB_DLLOAD

SWIN_LIB_PGPLOT
SWIN_LIB_FFTW
SWIN_LIB_CFITSIO


# No longer need openmp.
#AX_OPENMP([have_openmp=yes],[have_openmp=no])
#if test X$have_openmp == Xyes; then
#	CFLAGS="$CFLAGS $OPENMP_CFLAGS -DHAVE_OPENMP"
#	CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS -DHAVE_OPENMP"
#fi



#GSL basic stuff

AX_PATH_GSL(1.11.0, have_gsl=yes,have_gsl=no)
AM_CONDITIONAL(HAVE_GSL, [test x"$have_gsl" = xyes])

AC_F77_WRAPPERS

ACX_BLAS(have_blas=yes,have_blas=no)
AC_MSG_CHECKING([for BLAS])
AC_MSG_RESULT($have_blas)
AM_CONDITIONAL(HAVE_BLAS, [test x$have_blas = xyes])

# Test for calceph 
ACX_CALCEPH(have_calceph=yes, have_calceph=no)
AC_MSG_CHECKING([for CALCEPH])
AC_MSG_RESULT($have_calceph)
AM_CONDITIONAL(HAVE_CALCEPH, [test x"$have_calceph" = xyes])

# test for LAPACK acceleration
ACX_LAPACK(have_lapack=yes, have_lapack=no)
AC_MSG_CHECKING([for LAPACK])
AC_MSG_RESULT($have_lapack)
AM_CONDITIONAL(HAVE_LAPACK, [test x"$have_lapack" = xyes])

# TempoNest addition:
# Test for MultiNest and CUDA/CULA installation
SWIN_LIB_MULTINEST
SWIN_LIB_CULA


AC_MSG_CHECKING([for LOGIN_ARCH variable])
if test -n "$LOGIN_ARCH" ; then
	t2arch=$LOGIN_ARCH
	else if test -n "$OSTYPE" ; then
		AC_MSG_RESULT(no)
		AC_MSG_CHECKING([for OSTYPE variable])
		t2arch=$OSTYPE
		else if test -n $OSTYPE ; then
		AC_MSG_RESULT(no)
		AC_MSG_CHECKING([for uname])
			t2arch=`uname`
		fi
	fi
fi


AC_MSG_RESULT($t2arch)

AC_DEFINE_UNQUOTED([TEMPO2_ARCH], "$t2arch",["Tempo2 ARCH variable"])
T2ARCH=$t2arch

AC_SUBST(T2ARCH)


AM_INIT_AUTOMAKE([foreign subdir-objects])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile plugin/Makefile sofa/Makefile
		 unsupported_plugins/Makefile
		 unsupported_plugins/toasim/Makefile
		 unsupported_plugins/ppta/Makefile
		 unsupported_plugins/glast/Makefile
		 temponest/Makefile])

AC_OUTPUT

AC_MSG_NOTICE(Fortran compiler is $F77 $FFLAGS)
AC_MSG_NOTICE(C compiler is $CC $CFLAGS)
AC_MSG_NOTICE(C++ compiler is $CXX $CXXFLAGS)
AC_MSG_NOTICE(Global linking flags are $LDFLAGS)


# User warnings
if test X$have_pgplot != Xyes; then
	AC_MSG_NOTICE(Plugins: PGPLOT was not detected - plugins requiring pgplot will not be compiled)
else
	AC_MSG_NOTICE(PGPLOT compile flags are "$PGPLOT_CFLAGS")
	AC_MSG_NOTICE(PGPLOT linking flags are "$PGPLOT_LIBS")
fi

if test X$have_gsl != Xyes; then
	AC_MSG_NOTICE(Plugins: GSL was not found - plugins requiring GSL will not be compiled)
fi

if test X$have_fftw3 != Xyes; then
	AC_MSG_NOTICE(Plugins: FFTW3 was not found - plugins requiring fftw will not be compiled)
fi

if test X$have_cfitsio != Xyes; then
	AC_MSG_NOTICE(Plugins: cfitsio was not found - plugins requiring cfitsio will not be compiled)
fi

if test X$have_lapack = Xyes; then
	AC_MSG_NOTICE(Using LAPACK library for Cholesky decomposition)
fi

if test X$have_blas = Xyes; then
	AC_MSG_NOTICE(Using BLAS library for matrix multiplication)
fi

if test X$have_calceph = Xyes; then
	AC_MSG_NOTICE(Able to use calceph library for ephemeris reading)
else
	AC_MSG_NOTICE(Calceph library was not found)
fi

# TempoNest addition:
if test X$have_multinest != Xyes; then
	AC_MSG_NOTICE(MultiNest was not found - TempoNest can not be built)

    else if test X$have_cula = Xyes; then
        AC_MSG_NOTICE(TempoNest will be built with CUDA/CULA GPU support)
    else
        AC_MSG_NOTICE(TempoNest will be built without GPU support)
    fi
fi
if test $can_build_shared == no ; then
	echo ""
	AC_WARN("Libtool cannot build shared libraries - no plugins can be built!")
	echo ""
else
	echo ""
	echo "NOTICE: plugins can be built with 'make plugins && make plugins-install'"
    if test X$have_multinest = Xyes; then
        echo "        temponest can be built with 'make temponest && make temponest-install'"

	   echo ""
	   echo "Debugging TempoNest configure, library summary:"
	   echo " * MULTINEST:    $MULTINEST_LIBS   $MULTINEST_CFLAGS"
	   echo " * CULA:         $CULA_LIBS        $CULA_CFLAGS"
	   echo " * GSL:          $GSL_LIBS         $GSL_CFLAGS"
	   echo " * BLAS:         $BLAS_LIBS        $BLAS_CFLAGS"
	   echo " * CALCEPH:      $CALCEPH_LIBS     $CALCEPH_CFLAGS"
	   echo " * LAPACK:       $LIBS_LAPACK"
    fi
    echo ""
fi

echo ""
echo "****"
echo "Ok! Good to start compiling... try make && make install'"
echo "****"



