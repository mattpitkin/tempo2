#ifndef _TKlog_h
#define _TKlog_h
#define TK_MAX_ERRORS 16
#define TK_MAX_ERROR_LEN 128

#ifdef __cplusplus
#include <cstdio>
#include <ctime>
extern "C" {
#else
#include <stdio.h>
#include <time.h>
#endif
    extern int debugFlag;   /* Global = 1 if debug mode is running */
    extern int writeResiduals;   /* Global. Bit1=prefit, bit2=designmatrix, bit3=postfit. Indicate we are writing out post-fit residuals */
    extern int tcheck;   /* Global = 1 if time check message should be printed is running */
    extern clock_t timer_clk;
    extern unsigned TK_errorCount;
    extern unsigned TK_warnCount;
    extern char TK_errorlog[TK_MAX_ERRORS][TK_MAX_ERROR_LEN];
    extern char TK_warnlog[TK_MAX_ERRORS][TK_MAX_ERROR_LEN];
    int logerr_check();
    void _TKchklog(FILE*, const char*, ...);

#ifdef __cplusplus
}
#endif

/* define some functions for log message 
 * M.Keith 2012 - let me know if this fails to compile anywhere.
 * mkeith@pulsarastronomy.net
 **/
#ifndef LOG_OUTFILE
#define LOG_OUTFILE stdout
#endif
#define RESETCOLOR "\033[0m"
#define WARNCOLOR  RESETCOLOR "\033[0;35m"
#define BOLDCOLOR  RESETCOLOR "\033[1m"
#define ERRORCOLOR RESETCOLOR "\033[1;31m"
#define WHERESTR  "[%s:%d] "
#define WHEREARG  __FILE__, __LINE__
#define ENDL "\n"
#define WHEREERR ERRORCOLOR "***ERROR***\n [%s:%d] " RESETCOLOR
#define WHEREWARN BOLDCOLOR "[%s:%d] " WARNCOLOR "Warning: " RESETCOLOR
#define ENDERR "\n***!!!!!***"
#define WHERETCHK "[%s:%d] T=%.2f s: "
#define _LOG(_fmt,...) _TKchklog(LOG_OUTFILE,_fmt,##__VA_ARGS__)
//fprintf(LOG_OUTFILE,_TKchklog(_fmt),##__VA_ARGS__)
#define logmsg(_fmt, ...) _LOG(WHERESTR _fmt ENDL, WHEREARG,##__VA_ARGS__)
#define logdbg(_fmt, ...)  if(debugFlag)logmsg(_fmt,##__VA_ARGS__)
#define logerr(_fmt, ...) do{TK_STORE_ERROR(_fmt,##__VA_ARGS__); _LOG(WHEREERR _fmt ENDERR ENDL,  WHEREARG,##__VA_ARGS__);}while(0)
#define logwarn(_fmt, ...) do{TK_STORE_WARNING(_fmt,##__VA_ARGS__); _LOG(WHEREWARN _fmt ENDL,  WHEREARG,##__VA_ARGS__);}while(0)
#define logtchk(_fmt, ...) if(tcheck)_LOG(WHERETCHK _fmt ENDL, WHEREARG,(clock()-timer_clk)/(float)CLOCKS_PER_SEC,##__VA_ARGS__)
#define TK_STORE_ERROR(_fmt,...) if(TK_errorCount < TK_MAX_ERRORS)snprintf(TK_errorlog[TK_errorCount],TK_MAX_ERROR_LEN, _fmt,##__VA_ARGS__); ++TK_errorCount
#define TK_STORE_WARNING(_fmt,...) if(TK_warnCount < TK_MAX_ERRORS)snprintf(TK_warnlog[TK_warnCount],TK_MAX_ERROR_LEN, _fmt,##__VA_ARGS__); ++TK_warnCount

#ifdef __GNUC__
#define DEPRECATED __attribute__ ((deprecated))
#else
#define DEPRECATED
#endif

#endif //_TKlog_h
