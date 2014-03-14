dnl @synopsis ACX_CALCEPH([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl

AC_DEFUN([ACX_CALCEPH], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
acx_calceph_ok=no

AC_ARG_WITH(calceph,
	[AC_HELP_STRING([--with-calceph=<lib>], [use CALCEPH library <lib>])])
case $with_calceph in
	yes | "") ;;
	no) acx_calceph_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) CALCEPH_LIBS="$with_calceph" ;;
	*) CALCEPH_LIBS="-l$with_calceph" ;;
esac


AC_CHECK_FUNC([calceph_close])

acx_calceph_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check CALCEPH_LIBS environment variable
if test $acx_calceph_ok = no; then
if test "x$CALCEPH_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$CALCEPH_LIBS $LIBS"
	AC_MSG_CHECKING([for $calceph_close in $CALCEPH_LIBS])
#	AC_TRY_LINK_FUNC($calceph_close, [acx_calceph_ok=yes], [CALCEPH_LIBS=""])
#	AC_MSG_RESULT($acx_calceph_ok)
acx_calceph_ok=yes 
	LIBS="$save_LIBS"
        CALCEPH_LIBS="-lcalceph"
fi
fi

# Generic CALCEPH library?
if test $acx_calceph_ok = no; then
	AC_CHECK_LIB(calceph, $calceph_close, [acx_calceph_ok=yes; CALCEPH_LIBS="-lcalceph"])
fi

AC_SUBST(CALCEPH_LIBS)

LIBS="$acx_calceph_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_calceph_ok" = xyes; then
   AC_DEFINE(HAVE_CALCEPH,1,[Define if you have the calceph library.])
   $1
else
        acx_calceph_ok=no
        $2
fi
])dnl ACX_CALCEPH
