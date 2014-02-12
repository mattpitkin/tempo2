#
# SWIN_LIB_MULTINEST([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
#
# This m4 macro checks availability of the package MultiNest
# by F. Feroz, M.P. hobson, and M. Bridges
#
# MULTINEST_CFLAGS - autoconfig variable with flags required for compiling
# MULTINEST_LIBS   - autoconfig variable with flags required for linking
# HAVE_MULTINEST   - automake conditional
# HAVE_MULTINEST   - pre-processor macro in config.h
#
# This macro tries to link a test program, by using something like
#
#    -L$MULTINEST_DIR -lnest3 -llapack -blas
#
# Notice that the environment variable MULTINEST_DIR is used. In the case the
# library 'libnest' is not installed in a default location, let MULTINEST_DIR
# point to the location of libnest
#
#  MULTINEST_LIBS="$MULTINEST_LIBS -lnest3 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread"
#  MULTINEST_LIBS="$MULTINEST_LIBS -lnest3 -llapack"
# ----------------------------------------------------------
AC_DEFUN([SWIN_LIB_MULTINEST],
[
  AC_PROVIDE([SWIN_LIB_MULTINEST])
  AC_REQUIRE([ACX_BLAS])
  AC_REQUIRE([ACX_LAPACK])
  AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])

  MULTINEST_CFLAGS=""
  MULTINEST_LIBS=""

  if test x"$MULTINEST_DIR" != x; then
    MULTINEST_CFLAGS="-I$MULTINEST_DIR"
    MULTINEST_LIBS="-L$MULTINEST_DIR"
  fi

  MULTINEST_LIBS="$MULTINEST_LIBS $BLAS_LIBS $LIBS_LAPACK -lnest3"

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $MULTINEST_LIBS"
  CFLAGS="$ac_save_CFLAGS $MULTINEST_CFLAGS"

  AC_CHECK_LIB(nest3, __nested_MOD_nestrun,
            [

                compile_multinest=yes
            ],
            [
                AC_CHECK_LIB(nest3, nested_mp_nestrun_,
                    [
                        compile_multinest=yes
                    ])
            ])

  if test x"$compile_multinest" = xyes; then
      AC_MSG_NOTICE(Running MultiNest to confirm runtime linking...)
      AC_LANG_PUSH(C++)
      AC_TRY_RUN(
          [
#ifdef __INTEL_COMPILER 			/* if the MultiNest library was compiled with ifort */
        #define NESTRUN nested_mp_nestrun_
#elif defined __GNUC__ 				/* if the MultiNest library was compiled with gfortran */
           #define NESTRUN __nested_MOD_nestrun
#else
           #error Do not know how to link to Fortran libraries, check symbol table for your platform (nm libnest3.a | grep nestrun) & edit example_eggbox_C++/eggbox.cc
#endif

    extern "C" {
        void NESTRUN(int &IS, int &mmodal, int &ceff, int &nlive, double &tol, double &efr, int &ndims,
            int &nPar, int &nClsPar, int &maxModes, int &updInt, double &Ztol, char *root, int &seed,
            int *pWrap, int &fb, int &resume, int &outfile, int &initMPI, double &logZero, int &maxiter,
            void (*Loglike)(double *Cube, int &n_dim, int &n_par, double &lnew, void *),
            void (*dumper)(int &, int &, int &, double **, double **, double **, double &, double &, double &, void *),
            void *context, int &root_len);
    }


#include <stdio.h>
#include <string.h>

    void loglik(double *pd, int &ndims, int &ndim, double &lnew, void *context) {
        lnew = -5 * (pd[0]*pd[0] + pd[1]*pd[1]);
        return;
    } /* loglik */

    void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &logZerr, void *context) {
        return;
    } /* dumper */

    int main(int argc, const char **argv) {
        int IS = 1, mmodal=1, ceff = 0, nlive = 10, ndims=2, nPar=2, nClsPar=2, maxModes=2, updInt=50,
            seed=-1, i, fb=1, outfile=0, initMPI=1, resume=0, maxiter=20, pWrap[2];
        double tol=0.5, efr=0.8, Ztol=-1E90, logZero = -1.0e99;
        void *context=NULL;
        char root[]="./test";
        int root_len = strlen(root);

        for(i=0; i<2; i++) { pWrap[i]=0; }

        NESTRUN(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, &loglik, &dumper, context, root_len);

        return 0;
    } /* main */
          ],
          [
              AC_MSG_NOTICE(Running MultiNest successful)
              have_multinest=yes
          ],
          [
              AC_MSG_NOTICE(Running MultiNest failed)
          ])
      AC_LANG_POP(C++)
  fi

  AC_MSG_CHECKING([for MultiNest installation])
  AC_MSG_RESULT($have_multinest)

  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"

  if test x"$have_multinest" = xyes; then
    AC_DEFINE([HAVE_MULTINEST], [1], [Define to 1 if you have the MULTINEST library])
    [$1]
  else
    AC_MSG_WARN([MultiNest code will not be compiled. This only affects MultiNest plugins.])

    if test x"$compile_multinest" = xyes; then
      AC_MSG_WARN([***************************************************************])
      AC_MSG_WARN([MultiNest code can be compiled, but it cannot be run. This most])
      AC_MSG_WARN([likely means that the MultiNest library libnest3.so is not in the])
      AC_MSG_WARN([path. Please set (DY)LD_LIBRARY_PATH correctly.])
      AC_MSG_WARN([***************************************************************])
    else
      if test x"$MULTINEST_DIR" = x; then
        AC_MSG_NOTICE([MultiNest plugins will not be compiled. If you need this, set MULTINEST_DIR.])
# Reduced warning message size
#        AC_MSG_WARN([***************************************************************])
#        AC_MSG_WARN([Please set the MULTINEST_DIR environment variable, or set the LIBRARY_PATH variable])
#        AC_MSG_WARN([       Note: LIBRARY_PATH is used at linking stage, (DY)LD_LIBRARY_PATH is used at runtime])
#        AC_MSG_WARN([***************************************************************])
      fi
    fi

    MULTINEST_CFLAGS=""
    MULTINEST_LIBS=""
    [$2]
  fi

  AC_SUBST(MULTINEST_CFLAGS)
  AC_SUBST(MULTINEST_LIBS)
  AM_CONDITIONAL(HAVE_MULTINEST, [test x"$have_multinest" = xyes])
])

