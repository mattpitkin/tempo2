//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russel Edwards

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
#include <tempo2.h>
#include <string.h>
#include <math.h>
#include "TKfit.h"

void updateDMvals(pulsar *psr,int p);


extern "C" int pluginFitFunc(pulsar *psr,int npsr,int writeModel) 
{
	int p,i,k;
	int count=0;
	int flags[MAX_PSR];
	longdouble* preobs[MAX_PSR];
	printf(" <DMMODEL> Reset DM component of error to 0\n");

	for (p=0; p < npsr; p++){
		if (psr[p].param[param_dmmodel].fitFlag[0]==1){
			flags[p]=1;
			count++;
			for(i=0; i < psr[p].nobs; i++){
				psr[p].obsn[i].toaErr=psr[p].obsn[i].origErr;
				psr[p].obsn[i].toaDMErr=0;
			}
		} else flags[p]=0;
	}
	strcpy(psr[0].fitFunc,"default");

	printf(" <DMMODEL> First fit...\n");
	if (strcmp(dcmFile,"NULL")==0 && strcmp(covarFuncFile,"NULL")==0)
		doFit(psr,npsr,writeModel);
	else
		doFitDCM(psr,dcmFile,covarFuncFile,npsr,0);

	if(count){

		printf(" <DMMODEL> Disable DMMODEL\n");
		// we also save the pre-fit residuals so it appears as if only one fit has happened
		for (p=0; p < npsr; p++){
			preobs[p] = (longdouble*)malloc(sizeof(longdouble)*psr[p].nobs);
			for(i=0; i < psr[p].nobs; i++){
				preobs[p][i]=psr[p].obsn[i].prefitResidual;
			}
			psr[p].param[param_dmmodel].fitFlag[0]=0;
		}
			

		formBatsAll(psr,npsr);
		formResiduals(psr,npsr,0);

		printf(" <DMMODEL> Second Fit...\n");
		if (strcmp(dcmFile,"NULL")==0 && strcmp(covarFuncFile,"NULL")==0)
			doFit(psr,npsr,writeModel);
		else
			doFitDCM(psr,dcmFile,covarFuncFile,npsr,0);


			printf(" <DMMODEL> Enable DMMODEL\n");
		for (p=0; p < npsr; p++){
			psr[p].param[param_dmmodel].fitFlag[0]=flags[p];
			// restore the pre-fit residuals
			for(i=0; i < psr[p].nobs; i++){
				psr[p].obsn[i].prefitResidual=preobs[p][i];
			}
		}
	}
	strcpy(psr[0].fitFunc,"dmmodel");
	return 0;

}
