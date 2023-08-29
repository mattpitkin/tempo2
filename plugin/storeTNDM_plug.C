//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

/*
 *    This file is part of TEMPO2. 
 * 
 *    TEMPO2 is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    TEMPO2 is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
 */

/*
 *    If you use TEMPO2 then please acknowledge it by citing 
 *    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
 *    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
 *    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
 *    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
 *    timing model.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"

extern "C" int tempoOutput(int argc,char *argv[],pulsar *psr,int npsr) 
{
    FILE* out= fopen("tndm_offset.txt","w");

    for (int iobs = 0 ; iobs < psr->nobs; ++iobs) {
        double freq = (double)psr[0].obsn[iobs].freqSSB;
        double dm_offset = psr[0].obsn[iobs].TNDMSignal*(DM_CONST*pow(freq/1e6,2));

        int flagid=psr[0].obsn[iobs].nFlags;
        for(int k=0; k < psr[0].obsn[iobs].nFlags ; k++){
            if(strcmp(psr[0].obsn[iobs].flagID[k],"-dmo")==0){
                flagid=k;
                break;
            }
        }

        strcpy(psr[0].obsn[iobs].flagID[flagid],"-dmo");
        sprintf(psr[0].obsn[iobs].flagVal[flagid],"%lg",dm_offset);

        if (flagid==psr[0].obsn[iobs].nFlags) {
            psr[0].obsn[iobs].nFlags++;
        }
        if ( psr[0].obsn[iobs].deleted ){
            continue;
        }
        ld_fprintf(out,"%.15Lf %lg %lg\n",psr[0].obsn[iobs].sat, dm_offset, psr[0].obsn[iobs].TNDMSignal);
    }
    void writeTim(const char *timname,pulsar *psr,const char *fileFormat);

    writeTim("tndm.tim",psr,"tempo2");

    fclose(out);
    return 0;
}
