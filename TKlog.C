#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "TKlog.h"
int debugFlag = 0;
unsigned TK_errorCount = 0;
unsigned TK_warnCount = 0;
int writeResiduals=0;
int tcheck = 0;
clock_t timer_clk = 0;
char TK_errorlog[TK_MAX_ERRORS][TK_MAX_ERROR_LEN];
char TK_warnlog[TK_MAX_ERRORS][TK_MAX_ERROR_LEN];

int logerr_check(){
    const unsigned count=TK_errorCount;
    TK_errorCount=0;

    if(TK_warnCount){
        _LOG(WARNCOLOR "Notice:" RESETCOLOR " There were %u warnings. Sumaries are shown below, check logs for full details.\n",TK_warnCount);
        for (unsigned i=0; i < TK_warnCount; ++i){
            _LOG(BOLDCOLOR "Warning #%u:" RESETCOLOR " %s\n",i+1,TK_warnlog[i]);
            if(i > TK_MAX_ERROR_LEN){
                _LOG(WARNCOLOR "   ... too many warnings\n" RESETCOLOR);
                break;
            }
        }
        TK_warnCount=0;
    }

    if(count){
        _LOG("\n\n" ERRORCOLOR "Notice:" RESETCOLOR " There were %u errors! Sumaries are shown below, check logs for full details.\n",count);
        for (unsigned i=0; i < count; ++i){
            _LOG(ERRORCOLOR "Error #%u" RESETCOLOR " = " BOLDCOLOR "%s" RESETCOLOR ENDL,i+1,TK_errorlog[i]);
            if(i > TK_MAX_ERROR_LEN){
                _LOG(ERRORCOLOR "  ... too many errors\n" RESETCOLOR);
                break;
            }

        }
    }
    return count;
}
