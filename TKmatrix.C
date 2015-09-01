#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <map>
#include "TKlongdouble.h"
#include "TKmatrix.h"
#include "TKlog.h"
#include "T2accel.h"

#include "TKmatrix.cpp"

/// C functions

#define TK_POSTFIX _d
#define TKmatrix_D double
#include "TKmatrix.detail.cpp"
#undef TKmatrix_D
#undef TK_POSTFIX

#define TK_POSTFIX _f
#define TKmatrix_D float
#include "TKmatrix.detail.cpp"
#undef TKmatrix_D
#undef TK_POSTFIX

#define TK_POSTFIX _L
#define TKmatrix_D longdouble
#include "TKmatrix.detail.cpp"
#undef TKmatrix_D
#undef TK_POSTFIX


