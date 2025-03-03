/**
 * This file has the C++ code for the tempo2 predictors. This is only used 
 * interally by tempo2 and not in the libtempo2pred library.
 */

#include "tempo2.h"
#include "tempo2pred_int.h"



void ChebyModel_Construct(ChebyModel *cm, const pulsar *psr);
void ChebyModel_Test(ChebyModel *cm, const pulsar *psr, int nmjd, int nfreq,
        long double *residualRMS, long double *residualMAV);
void ChebyModelSet_Construct(ChebyModelSet *cms, const pulsar *psr,
        const char *sitename,
        long double mjd_start, long double mjd_end,
        long double segment_length, long double overlap,
        long double freq_start,long double freq_end,
        int nmjdcoeff, int nfreqcoeff);
void ChebyModelSet_Test(ChebyModelSet *cms, const pulsar *psr,
        int nmjd, int nfreq,
        long double *residualRMS, long double *residualMAV);