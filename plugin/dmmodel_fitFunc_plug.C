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
	char flags[MAX_PSR];
	char Jflags[MAX_PSR][MAX_JUMPS];
	observation* preobs[MAX_PSR];
	printf("TEST\n");

	for (p=0; p < npsr; p++){
		if (psr[p].param[param_dmmodel].fitFlag[0]==1){
			flags[p]=1;
			count++;
			for(i=0; i < psr[p].nobs; i++){
				psr[p].obsn[i].toaErr=psr[p].obsn[i].origErr;
			}
			
		} else flags[p]=0;
	}
	strcpy(psr[0].fitFunc,"default");

	printf("ONE\n");
	if (strcmp(dcmFile,"NULL")==0 && strcmp(covarFuncFile,"NULL")==0)
		doFit(psr,npsr,writeModel);
	else
                doFitDCM(psr,dcmFile,covarFuncFile,npsr,0);



	if(count){

		formBatsAll(psr,npsr);
		formResiduals(psr,npsr,0);



		for (p=0; p < npsr; p++){
			preobs[p] = (observation*)malloc(sizeof(observation)*psr[p].nobs);
			memcpy(preobs[p],psr[p].obsn,sizeof(observation)*psr[p].nobs);
			double sum=0;
			for(i=0; i < psr[p].dmoffsNum; i++){
				sum+=psr[p].dmoffsDM[i];
			}
			sum/= (double)psr[p].dmoffsNum;
			for(i=0; i < psr[p].dmoffsNum; i++){
				psr[p].dmoffsDM[i]-=sum;
			}
			psr[p].param[param_dmmodel].val[0]+=(long double)sum;
		}
		for (p=0; p < npsr; p++){
			for(i=0; i < psr[p].nJumps; i++){
				Jflags[p][i]=psr[p].fitJump[i];
				psr[p].fitJump[i]=0;
			}
		}
		printf("TWO\n");
		if (strcmp(dcmFile,"NULL")==0 && strcmp(covarFuncFile,"NULL")==0)
			doFit(psr,npsr,writeModel);
		else
			doFitDCM(psr,dcmFile,covarFuncFile,npsr,0);


		formBatsAll(psr,npsr);
		formResiduals(psr,npsr,0);

		//	textOutput(psr,npsr,0.0,0,0,0,"");
		for (p=0; p < npsr; p++){
			if (psr[p].param[param_dmmodel].fitFlag[0]==1){
				updateDMvals(psr,p);
				double os=0;
				double ds=0;
				double ns=0;
				for(i=0; i < psr[p].nobs; i++){
					ds+=psr[p].obsn[i].toaDMErr;
					os+=psr[p].obsn[i].origErr;
					ns+=psr[p].obsn[i].toaErr;
				}
				printf("Mean orig  TOA error = %f\n",os/(double)psr[p].nobs);
				printf("Mean DM    TOA error = %f\n",ds/(double)psr[p].nobs);
				printf("Mean total TOA error = %f\n",ns/(double)psr[p].nobs);
			}
		}
		for (p=0; p < npsr; p++){
			psr[p].param[param_dmmodel].fitFlag[0]=0;
			double sum=0;
			for(i=0; i < psr[p].dmoffsNum; i++){
				sum+=psr[p].dmoffsDM[i];
			}
			sum/= (double)psr[p].dmoffsNum;
			for(i=0; i < psr[p].dmoffsNum; i++){
				psr[p].dmoffsDM[i]-=sum;
			}
			psr[p].param[param_dmmodel].val[0]+=(long double)sum;
		}




		printf("THREE\n");
		if (strcmp(dcmFile,"NULL")==0 && strcmp(covarFuncFile,"NULL")==0)
			doFit(psr,npsr,writeModel);
		else
			doFitDCM(psr,dcmFile,covarFuncFile,npsr,0);



//		textOutput(psr,npsr,0.0,0,0,0,"");
		for (p=0; p < npsr; p++){
			psr[p].param[param_dmmodel].fitFlag[0]=flags[p];
			memcpy(psr[p].obsn,preobs[p],sizeof(observation)*psr[p].nobs);
			updateDMvals(psr,p);
			free(preobs[p]);
			for(i=0; i < psr[p].nJumps; i++){
				psr[p].fitJump[i]=Jflags[p][i];
			}

		}
	}
	strcpy(psr[0].fitFunc,"dmmodel");
	return 0;

}
