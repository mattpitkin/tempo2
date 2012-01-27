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
void getFitLabels(pulsar* psr, int p,char** ret);


extern "C" int pluginFitFunc(pulsar *psr,int npsr,int writeModel) 
{
	int p,i,k,j;
	int count=0;
	int flags[MAX_PSR];
	longdouble* preobs[MAX_PSR];
	FILE *dmf;
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


	for (p=0; p < npsr; p++){
		char** labels = (char**)malloc(sizeof(char*)*MAX_PARAMS);
		for (i=0;i<MAX_PARAMS;i++){
			labels[i]= (char*)malloc(80);
			strcpy(labels[i],"UNK_PARAM");
		}
		getFitLabels(psr,p,labels);
		double** cvm=psr[p].covar;
		int npol = psr[p].nFit - psr[p].fitNfree;
		bool warn=false;
		for (i=0;i<npol;i++){
			for (j=0;j<i;j++){
				double cv=fabs(cvm[i][j]/sqrt((cvm[j][j])*(cvm[i][i])));
				if(cv > 0.5){
					if(!warn){
						printf(" <DMMODEL> Warning: highly covariant parameters in fit!\n");
						warn=true;
					}
					printf("  % 15s % 15s %lg\n",labels[i],labels[j],cv);
				}
			}
		}
		for (i=0;i<MAX_PARAMS;i++)free(labels[i]);
		free(labels);
	}


	if(count){


		{
			p=0;
			dmf=fopen("dm.model","w");
			for (i=0;i<(int)psr[p].dmoffsNum;i++){
				fprintf(dmf,"%lf %.15lg %.15lg %.15lg %.15lg\n",(double)psr[p].dmoffsMJD[i],(double)psr[p].dmoffsDM[i],(double)psr[p].dmoffsDMe[i],(double)psr[p].dmoffsOffset[i],(double)psr[p].dmoffsError[i]);
			}
			fclose(dmf);
			dmf=fopen("dm.toas","w");
			double dmval=0;
			double meanDM=0;
			double m,c;
			for (i=0;i<(int)psr[p].nobs;i++){
				if ((double)psr[p].obsn[i].sat < psr[p].dmoffsMJD[0])
					dmval = meanDM + psr[p].dmoffsDM[0];
				else if ((double)psr[p].obsn[i].sat > psr[p].dmoffsMJD[psr[p].dmoffsNum-1])
					dmval = meanDM+psr[p].dmoffsDM[psr[p].dmoffsNum-1];
				else
				{
					for (k=0;k<psr[p].dmoffsNum-1;k++)
					{
						if ((double)psr[p].obsn[i].sat >= psr[p].dmoffsMJD[k] &&
								(double)psr[p].obsn[i].sat < psr[p].dmoffsMJD[k+1])
						{
							// Do linear interpolation
							// Note: this is also used at various points in the 
							// code (e.g., textOutput.C - if any changes are 
							// made here then these changes should be made throughout!
							m = (psr[p].dmoffsDM[k]-psr[p].dmoffsDM[k+1])/(psr[p].dmoffsMJD[k]-psr[p].dmoffsMJD[k+1]);
							c = psr[p].dmoffsDM[k]-m*psr[p].dmoffsMJD[k];
							dmval = m*(double)psr[p].obsn[i].sat+c + meanDM;
							break;
						}
					}
				}
				double freqf= psr[p].obsn[i].freq*1.0e6;
				double delay= dmval/DM_CONST/1.0e-12/freqf/freqf;
				fprintf(dmf,"%lf %.10lg %.10lg %.10lg %.10lg\n",(double)psr[p].obsn[i].sat,dmval,delay,(double)psr[p].obsn[i].toaErr*1e-6,freqf);
			}

		}
		fclose(dmf);


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



void getFitLabels(pulsar* psr, int p,char** ret){
	int i,j,k;
	int n=0;
	sprintf(ret[n++],"ZERO_OFF");
	for (i=0;i<MAX_PARAMS;i++)
	{
		for (k=0;k<psr[p].param[i].aSize;k++)
		{
			if (psr[p].param[i].fitFlag[k]==1) /* If we are fitting for this parameter */
			{
				if (i!=param_start && i!=param_finish)
				{
					if (i==param_wave_om)
					{
						if (psr[p].waveScale==2)
						{
							//                      for (j=0;j<psr[p].nWhite*2;j++)                           
							for (j=0;j<psr[p].nWhite*4;j++)
								sprintf(ret[n++],"%s_%02d",psr[p].param[i].shortlabel[k],j);
						}
						else
						{
							for (j=0;j<psr[p].nWhite*2;j++)
								sprintf(ret[n++],"%s_%02d",psr[p].param[i].shortlabel[k],j);
						}
					}
					else if (i==param_quad_om)
					{
						for (j=0;j<psr[p].nQuad*4;j++)
							sprintf(ret[n++],"%s_%02d",psr[p].param[i].shortlabel[k],j);
					}
					else if (i==param_ifunc)
					{
						for (j=0;j<psr[p].ifuncN;j++)
							sprintf(ret[n++],"%s_%02d",psr[p].param[i].shortlabel[k],j);
					}
					else if (i==param_gwsingle)
					{
						sprintf(ret[n++],"%s_%02d",psr[p].param[i].shortlabel[k],0);
						sprintf(ret[n++],"%s_%02d",psr[p].param[i].shortlabel[k],1);
						sprintf(ret[n++],"%s_%02d",psr[p].param[i].shortlabel[k],2);
						sprintf(ret[n++],"%s_%02d",psr[p].param[i].shortlabel[k],3);
					}
					else if (i==param_dmmodel)
					{
						for (j=0;j<(int)psr[p].dmoffsNum;j++)
						{
							sprintf(ret[n++],"%s_D_%02d",psr[p].param[i].shortlabel[k],j);
							sprintf(ret[n++],"%s_C_%02d",psr[p].param[i].shortlabel[k],j);
						}
					}
					else
						sprintf(ret[n++],"%s",psr[p].param[i].shortlabel[k]);
				}
			}
		}
	}
	for (i=1;i<=psr[p].nJumps;i++)
	{
		if (psr[p].fitJump[i]==1)
		{
			sprintf(ret[n++],"%s",psr[p].jumpStr[i]);
		}
	}

}
