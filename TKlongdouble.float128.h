#define USE_BUILTIN_LONGDOUBLE
#include <quadmath.h>

typedef __float128 longdouble;
#define LONGDOUBLE_ONE 1.0Q

#define longdouble(a)


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
longdouble parse_longdouble(const char *str);


