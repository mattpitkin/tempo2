#
# SWIN_LIB_QDINSTALL([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
#
# This m4 macro checks availability of the high precision qdinstall libraries.
#
# QDINSTALL_CFLAGS - autoconfig variable with flags required for compiling
# QDINSTALL_LIBS   - autoconfig variable with flags required for linking
# HAVE_QDINSTALL   - automake conditional
# HAVE_QDINSTALL   - pre-processor macro in config.h
#
# Some environment variables are required for the QDINSTALL libraries to
# function, if they are not installed in standard locations
#
# $QDINSTALL need to be pointing to the installation
# directory
#
# This macro tries to link a test program in the following way
#
#    -L$(QDINSTALL)/lib  -lmblas_qd -lqdinstall_qd -lmblas_dd -lqdinstall_dd  -I$(QDINSTALL)/include
#
#
# ----------------------------------------------------------
AC_DEFUN([SWIN_LIB_QDINSTALL],
[
  AC_PROVIDE([SWIN_LIB_QDINSTALL])
  AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
  AC_REQUIRE([ACX_BLAS])
  AC_REQUIRE([ACX_LAPACK])

  AC_MSG_CHECKING([for QDINSTALL installation])

  QDINSTALL_CFLAGS=""
  QDINSTALL_LIBS=""

  QDINSTALL_LIB="-lqd"

  if test x"$QDINSTALL" != x; then
    QDINSTALL_LIBS="-L$QDINSTALL/lib"
  fi

  if test x"$QDINSTALL" != x; then
    QDINSTALL_CFLAGS="-I$QDINSTALL/include"
  fi

  QDINSTALL_LIBS="$QDINSTALL_LIBS $QDINSTALL_LIB"

  AC_LANG_PUSH(C++)
  ac_save_CFLAGS="$CFLAGS"
  ac_save_CXXFLAGS="$CXXFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $QDINSTALL_LIBS"
  CFLAGS="$ac_save_CFLAGS $QDINSTALL_CFLAGS"
  CXXFLAGS="$ac_save_CXXFLAGS $QDINSTALL_CFLAGS"

  # test compilation of simple program
  AC_TRY_LINK([#include "qd/qd_real.h"],[qd_real a; fpu_fix_start(0)],
              have_qdinstall=yes, have_qdinstall=no)


  AC_MSG_RESULT($have_qdinstall)

  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"

  if test x"$have_qdinstall" = xyes; then
    AC_DEFINE([HAVE_QDINSTALL], [1], [Define to 1 if you have the QDINSTALL library])
    [$1]
  else
    AC_MSG_NOTICE([Will compile without QD code. Might cause problems if you use tempo2nest])
    QDINSTALL_CFLAGS=""
    QDINSTALL_LIBS=""
    [$2]
  fi

  AC_SUBST(QDINSTALL_CFLAGS)
  AC_SUBST(QDINSTALL_LIBS)
  AM_CONDITIONAL(HAVE_QDINSTALL, [test x"$have_qdinstall" = xyes])

  AC_LANG_POP(C++)
])

