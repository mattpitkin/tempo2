#ifndef TKlongdouble_h
#define TKlongdouble_h
#define USE_BUILTIN_LONGDOUBLE

#include <math.h>

#ifdef sun
#include <sunmath.h> 
//* there is no such file ! J. Wang */
// Note: you sometimes need a compile line e.g. -I/opt/SUNWspro/WS6U2/include/cc
#endif

#define longdouble(a) a##L
#define LD_PI M_PI

#ifdef USE_BUILTIN_LONGDOUBLE
#define LONGDOUBLE_IS_IEEE754
typedef long double longdouble;
#define LONGDOUBLE_ONE 1.0L
/* OSes/architectures lacking built-in double; use "qd" library */
#else
#define LONGDOUBLE_IS_DDREAL
#include "dd.h"
typedef dd_real longdouble;
#define LONGDOUBLE_ONE "1.0"
dd_real pow(const dd_real &a, const dd_real &b);
//operator float(const dd_real &a);
#endif

/* function to get longdouble as string (%g style). This returns
   std::string, from which you can get a normal C char * like this:
   print_longdouble(x).c_str(). Unfortunately we can't just return char *
   directly as it would end up pointing to de-allocated memory, and we
   can't do it statically in case you want to say call this twice for
   one printf */
#ifdef __cplusplus
#include <string>
std::string print_longdouble(const longdouble &ld);
#endif

#ifdef __cplusplus
extern "C" {
#endif

    /* function to parse a string as longdouble */
    longdouble parse_longdouble(const char *str);


#ifdef __cplusplus
}
#endif

#define ld_printf printf
#define ld_fprintf fprintf
#define ld_sprintf sprintf

#endif
