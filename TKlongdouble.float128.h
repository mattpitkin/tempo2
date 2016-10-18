#ifndef TKlongdouble_h
#define TKlongdouble_h
#define USE_BUILTIN_LONGDOUBLE
#define LONGDOUBLE_IS_FLOAT128

#include <math.h>
#include <quadmath.h>

typedef __float128 longdouble;
#define LONGDOUBLE_ONE 1.0Q

#define longdouble(a) a##Q
#define FMT_LD "Q"

#define LD_PI M_PIq
#define cosl cosq
#define sinl sinq
#define floorl floorq
#define fabsl fabsq
#define powl powq

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

/* function to parse a string as longdouble */

#ifdef __cplusplus

inline __float128 cos(__float128 v) {
    return cosq(v);
}

inline __float128 sin(__float128 v) {
    return sinq(v);
}

inline __float128 fabs(__float128 v) {
    return fabsq(v);
}
inline __float128 sqrt(__float128 v) {
    return sqrtq(v);
}
inline __float128 asin(__float128 v) {
    return asinq(v);
}

inline __float128 atan2(__float128 v, __float128 v2) {
    return atan2q(v,v2);
}

inline __float128 pow(__float128 v, double p) {
    return powq(v,p);
}
inline __float128 exp(__float128 v) {
    return expq(v);
}
inline __float128 log10(__float128 v) {
    return log10q(v);
}
inline __float128 log(__float128 v) {
    return logq(v);
}
inline __float128 tan(__float128 v) {
    return tanq(v);
}


extern "C" {
#endif // c++

    longdouble parse_longdouble(const char *str);

    int ld_printf(const char *__format, ...);
    int ld_fprintf(FILE* __stream, const char *__format, ...);
    int ld_sprintf(char* __str, const char *__format, ...);

#ifdef __cplusplus
}
#endif

#endif
