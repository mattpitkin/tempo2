//  Copyright (C) 2013 Lindley Lentati

/*
*    This file is part of TempoNest 
* 
*    TempoNest is free software: you can redistribute it and/or modify 
*    it under the terms of the GNU General Public License as published by 
*    the Free Software Foundation, either version 3 of the License, or 
*    (at your option) any later version. 
*    TempoNest  is distributed in the hope that it will be useful, 
*    but WITHOUT ANY WARRANTY; without even the implied warranty of 
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
*    GNU General Public License for more details. 
*    You should have received a copy of the GNU General Public License 
*    along with TempoNest.  If not, see <http://www.gnu.org/licenses/>. 
*/

/*
*    If you use TempoNest and as a byproduct both Tempo2 and MultiNest
*    then please acknowledge it by citing Lentati L., Alexander P., Hobson M. P. (2013) for TempoNest,
*    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
*    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
*    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
*    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
*    timing model and MultiNest Papers here.
*/


#include <stdio.h>
#include <vector>
#include <gsl/gsl_sf_gamma.h>
#include "dgemm.h"
#include "dgemv.h"
#include "dpotri.h"
#include "dpotrf.h"
#include "dpotrs.h"
#include "tempo2.h"
#include "TempoNest.h"
#include "dgesvd.h"


void WhiteLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double *EQUAD;
	int pcount=0;
	
	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];

	double NLphase=0;
	for(int p=0;p<ndim;p++){

		Cube[p]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[p]+((MNStruct *)context)->Dpriors[p][0];
	}

	if(((MNStruct *)context)->doLinear==0){
	
		for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
			LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
		}

		NLphase=(double)LDparams[0];
		pcount++;
		for(int p=1;p<((MNStruct *)context)->numFitTiming;p++){
			((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
			pcount++;
		}
		for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
			((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
			pcount++;
		}
		
//		fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
//		formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
//		
//		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
//			Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+phase;
//		}
	
	}
	else if(((MNStruct *)context)->doLinear==1){
	
		for(int p=0;p < numfit; p++){
			Fitparams[p]=Cube[p];
			//printf("FitP: %i %g \n",p,Cube[p]);
			pcount++;
		}
		
		double *Fitvec=new double[((MNStruct *)context)->pulse->nobs];

		dgemv(((MNStruct *)context)->DMatrix,Fitparams,Fitvec,((MNStruct *)context)->pulse->nobs,numfit,'N');
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
			//printf("FitVec: %i %g \n",o,Fitvec[o]);
		}
		
		delete[] Fitvec;
	}

	
	double **Steps;
	if(((MNStruct *)context)->incStep > 0){
		
		Steps=new double*[((MNStruct *)context)->incStep];
		
		for(int i = 0; i < ((MNStruct *)context)->incStep; i++){
			Steps[i]=new double[2];
			Steps[i][0] = Cube[pcount];
			pcount++;
			Steps[i][1] = Cube[pcount];
			pcount++;
		}
	}
		

	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EFAC[o]=1;
		}
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EFAC[o]=Cube[pcount];
		}
		pcount++;
		
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int p=0;p< ((MNStruct *)context)->systemcount; p++){
			EFAC[p]=Cube[pcount];
			pcount++;
		}
	}				

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=0;
		}
	}
	else if(((MNStruct *)context)->numFitEQUAD == 1){
		EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=pow(10.0,2*Cube[pcount]);
		}
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEQUAD > 1){
                EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
                        EQUAD[o]=pow(10.0,2*Cube[pcount]);
			pcount++;
                }
        }

        if(((MNStruct *)context)->doLinear==0){

                fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
                formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                          Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+NLphase;
                }
        }


	if(((MNStruct *)context)->incStep > 0){
		for(int i = 0; i < ((MNStruct *)context)->incStep; i++){
			double StepAmp = Steps[i][0];
			double StepTime = Steps[i][1];

			for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
				double time=(double)((MNStruct *)context)->pulse->obsn[o1].bat ;
				if( time > StepTime){
					Resvec[o1] += StepAmp;
				}
			}
		}
	}





	double Chisq=0;
	double noiseval=0;
	double detN=0;
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){

                if(((MNStruct *)context)->numFitEQUAD < 2 && ((MNStruct *)context)->numFitEFAC < 2){
			if(((MNStruct *)context)->whitemodel == 0){
				noiseval=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[0],2) + EQUAD[0];
			}
			else if(((MNStruct *)context)->whitemodel == 1){
				noiseval= EFAC[0]*EFAC[0]*(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[0]);
			}
                }
                else if(((MNStruct *)context)->numFitEQUAD > 1 || ((MNStruct *)context)->numFitEFAC > 1){
			if(((MNStruct *)context)->whitemodel == 0){
                       	 	noiseval=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD[((MNStruct *)context)->sysFlags[o]];
			}
			else if(((MNStruct *)context)->whitemodel == 1){
                       	 	noiseval=EFAC[((MNStruct *)context)->sysFlags[o]]*EFAC[((MNStruct *)context)->sysFlags[o]]*(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o]]);
			}
                }
	
	

		Chisq += pow(Resvec[o],2)/noiseval;
		detN += log(noiseval);
	}
	
	//printf("White: %g %g \n", detN, Chisq);

	if(isnan(detN) || isinf(detN) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
	}
	else{
		lnew = -0.5*(((MNStruct *)context)->pulse->nobs*log(2*M_PI) + detN + Chisq);	
	}

	delete[] EFAC;
	delete[] Resvec;
	if(((MNStruct *)context)->incStep > 0){
		for(int i = 0; i < ((MNStruct *)context)->incStep; i++){
			delete[] Steps[i];
		}
		delete[] Steps;
	}


}


void WhiteMarginLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double **EFAC;
	double *EQUAD;
	int pcount=0;
	
	int totdims = ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	totdims += ((MNStruct *)context)->numFitEFAC*((MNStruct *)context)->EPolTerms + ((MNStruct *)context)->numFitEQUAD;
	totdims +=2*((MNStruct *)context)->incStep;
	if(((MNStruct *)context)->incRED==2)totdims+=((MNStruct *)context)->numFitRedCoeff;
	if(((MNStruct *)context)->incDM==2)totdims+=((MNStruct *)context)->numFitDMCoeff;
	if(((MNStruct *)context)->varyRedCoeff==1)totdims+=2;
	if(((MNStruct *)context)->varyDMCoeff==1)totdims+=2;
   	if(((MNStruct *)context)->incRED==3)totdims+=2;
    if(((MNStruct *)context)->incDM==3)totdims+=2;
	if(((MNStruct *)context)->yearlyDM == 1)totdims+=2;
 	if(((MNStruct *)context)->incsinusoid == 1)totdims+=3;   
   	if(((MNStruct *)context)->incGWB == 1)totdims+=1; 
	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];
	int fitcount=0;
	for(int p=0;p<totdims;p++){

		if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){
		Cube[pcount]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[pcount]+((MNStruct *)context)->Dpriors[p][0];

		pcount++;
		}
	}
	pcount=0;
	if(((MNStruct *)context)->doLinear==0){
	
		for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
			if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){

				LDparams[p]=Cube[fitcount]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
				fitcount++;
			}
			else if(((MNStruct *)context)->Dpriors[p][1] == ((MNStruct *)context)->Dpriors[p][0]){
				LDparams[p]=((MNStruct *)context)->Dpriors[p][0]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
			}


		}
		pcount=0;
		double phase=(double)LDparams[0];
		pcount++;
		for(int p=1;p<((MNStruct *)context)->numFitTiming;p++){
			((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
			pcount++;
		}
		for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
			((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
			pcount++;
		}
		
		fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
		formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
		
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+phase;

		}
	
	}
	else if(((MNStruct *)context)->doLinear==1){
		fitcount=0;

		
		for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
			if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){
				Fitparams[p]=Cube[fitcount];
				fitcount++;
			}
			else if(((MNStruct *)context)->Dpriors[p][1] == ((MNStruct *)context)->Dpriors[p][0]){
				Fitparams[p]=0;
			}
		}
		double *Fitvec=new double[((MNStruct *)context)->pulse->nobs];

		dgemv(((MNStruct *)context)->DMatrix,Fitparams,Fitvec,((MNStruct *)context)->pulse->nobs,numfit,'N');
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
		}
		
		delete[] Fitvec;
	}
	pcount=fitcount;

	if(((MNStruct *)context)->incStep > 0){
		for(int i = 0; i < ((MNStruct *)context)->incStep; i++){
			double StepAmp = Cube[pcount];
			pcount++;
			double StepTime = Cube[pcount];
			pcount++;
			for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
				if(((MNStruct *)context)->pulse->obsn[o1].bat > StepTime){
					Resvec[o1] += StepAmp;
				}
			}
		}
	}	
	
	
/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get White Noise vector///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////  


	double equadpriorterm=0;
	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double*[((MNStruct *)context)->EPolTerms];
		for(int n=1; n <=((MNStruct *)context)->EPolTerms; n++){
			EFAC[n-1]=new double[((MNStruct *)context)->systemcount];
			if(n==1){
				for(int o=0;o<((MNStruct *)context)->systemcount; o++){
					EFAC[n-1][o]=1;
				}
			}
			else{
                                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
                                       EFAC[n-1][o]=0;
                            }
			}
		}
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double*[((MNStruct *)context)->EPolTerms];
		for(int n=1; n <=((MNStruct *)context)->EPolTerms; n++){
			
			EFAC[n-1]=new double[((MNStruct *)context)->systemcount];
			if(n==1){
				for(int o=0;o<((MNStruct *)context)->systemcount; o++){
					
					EFAC[n-1][o]=Cube[pcount];
				}
				pcount++;
			}
			else{
                    for(int o=0;o<((MNStruct *)context)->systemcount; o++){

                            EFAC[n-1][o]=pow(10.0,Cube[pcount]);
                    }
                    pcount++;
            }
		}
		
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double*[((MNStruct *)context)->EPolTerms];
		for(int n=1; n <=((MNStruct *)context)->EPolTerms; n++){
			EFAC[n-1]=new double[((MNStruct *)context)->systemcount];
			if(n==1){
				for(int p=0;p< ((MNStruct *)context)->systemcount; p++){
					EFAC[n-1][p]=Cube[pcount];
					pcount++;
				}
			}
			else{
                                for(int p=0;p< ((MNStruct *)context)->systemcount; p++){
                                        EFAC[n-1][p]=pow(10.0,Cube[pcount]);
                                        pcount++;
                                }
                        }
		}
	}	

		

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=0;
		}
	}
	else if(((MNStruct *)context)->numFitEQUAD == 1){
		EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=pow(10.0,2*Cube[pcount]);
			equadpriorterm+=log(pow(10.0,Cube[pcount]));
		}
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEQUAD > 1){
        EQUAD=new double[((MNStruct *)context)->systemcount];
        for(int o=0;o<((MNStruct *)context)->systemcount; o++){
            EQUAD[o]=pow(10.0,2*Cube[pcount]);
	    equadpriorterm+=log(pow(10.0,Cube[pcount]));

			pcount++;
        }
    }




	double *Noise;	
	double *BATvec;
	Noise=new double[((MNStruct *)context)->pulse->nobs];
	BATvec=new double[((MNStruct *)context)->pulse->nobs];
	
	
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
	}
		
		
	if(((MNStruct *)context)->whitemodel == 0){
	
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			double EFACterm=0;
			for(int n=1; n <=((MNStruct *)context)->EPolTerms; n++){
				EFACterm=EFACterm + pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))/pow(pow(10.0,-7),n-1),n)*EFAC[n-1][((MNStruct *)context)->sysFlags[o]];
			}	

			Noise[o]= (pow(EFACterm,2) + EQUAD[((MNStruct *)context)->sysFlags[o]]);

		}
		
	}
	else if(((MNStruct *)context)->whitemodel == 1){
	
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){

			Noise[o]=(EFAC[0][((MNStruct *)context)->sysFlags[o]]*EFAC[0][((MNStruct *)context)->sysFlags[o]]*(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o]]));
		}
		
	}

	
	
/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get Time domain likelihood//////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////  


	double tdet=0;
	double timelike=0;

	for(int o=0; o<((MNStruct *)context)->pulse->nobs; o++){
		//printf("tlike %i %g %g %g \n", o, Resvec[o], Noise[o], timelike);
		timelike+=Resvec[o]*Resvec[o]/Noise[o];
		tdet += log(Noise[o]);
	}
		
		
/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Form the Design Matrix////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////  

	int TimetoMargin=0;
	for(int i =0; i < ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps; i++){
		if(((MNStruct *)context)->LDpriors[i][2]==1)TimetoMargin++;
	}


	double **TNDM=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		TNDM[i]=new double[TimetoMargin];
	}

	getCustomDMatrixLike(context, TNDM);
	

	double* S = new double[TimetoMargin];
	double** U = new double*[((MNStruct *)context)->pulse->nobs];
	for(int k=0; k < ((MNStruct *)context)->pulse->nobs; k++){
		U[k] = new double[((MNStruct *)context)->pulse->nobs];
	}
	double** VT = new double*[TimetoMargin]; 
	for (int k=0; k<TimetoMargin; k++) VT[k] = new double[TimetoMargin];

	dgesvd(TNDM,((MNStruct *)context)->pulse->nobs, TimetoMargin, S, U, VT);

	delete[]S;	

	for (int j = 0; j < TimetoMargin; j++){
		delete[]VT[j];
	}
	
	delete[]VT;
	
		
	

	for(int j=0;j<((MNStruct *)context)->pulse->nobs;j++){
		for(int k=0;k < TimetoMargin;k++){
				TNDM[j][k]=U[j][k];
		}
	}

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[]U[j];
	}
	delete[]U;


//////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////////////Do Noise Stuff////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////  


	if(((MNStruct *)context)->incsinusoid == 1){
	
		double start,end;
		int go=0;
		for (int i=0;i<((MNStruct *)context)->pulse->nobs;i++)
		  {
			if (((MNStruct *)context)->pulse->obsn[i].deleted==0)
			  {
			if (go==0)
			  {
				go = 1;
				start = (double)((MNStruct *)context)->pulse->obsn[i].bat;
				end  = start;
			  }
			else
			  {
				if (start > (double)((MNStruct *)context)->pulse->obsn[i].bat)
				  start = (double)((MNStruct *)context)->pulse->obsn[i].bat;
				if (end < (double)((MNStruct *)context)->pulse->obsn[i].bat)
				  end = (double)((MNStruct *)context)->pulse->obsn[i].bat;
			  }
			  }
		  }

		double maxtspan=1*(end-start);
		
		double sineamp=pow(10.0,Cube[pcount]);
		pcount++;
		double sinephase=Cube[pcount];
		pcount++;
		double sinefreq=pow(10.0,Cube[pcount])/maxtspan;
		pcount++;		
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]-= sineamp*sin(2*M_PI*sinefreq*(double)((MNStruct *)context)->pulse->obsn[o].bat + sinephase);
		}
	}



	if(((MNStruct *)context)->yearlyDM == 1){
		double DMKappa = 2.410*pow(10.0,-16);
		double yearlyamp=pow(10.0,Cube[pcount]);
		pcount++;
		double yearlyphase=Cube[pcount];
		pcount++;
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]-= yearlyamp*sin((2*M_PI/365.25)*(double)((MNStruct *)context)->pulse->obsn[o].bat + yearlyphase)/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));;
		}
	}


//////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////Form Total Matrices////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

	int totalsize=TimetoMargin;
	double **TotalMatrix=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i =0;i<((MNStruct *)context)->pulse->nobs;i++){
		TotalMatrix[i]=new double[totalsize];
		for(int j =0;j<totalsize; j++){
			TotalMatrix[i][j]=0;
		}
	}
	

	for(int i =0;i<((MNStruct *)context)->pulse->nobs;i++){
		for(int j =0;j<TimetoMargin; j++){
			TotalMatrix[i][j]=TNDM[i][j];
		}
	}

//////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////Do Algebra/////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

	double *NTd = new double[totalsize];
	double **NT=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		NT[i]=new double[totalsize];
	}

	double **TNT=new double*[totalsize];
	for(int i=0;i<totalsize;i++){
		TNT[i]=new double[totalsize];
	}
	
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		for(int j=0;j<totalsize;j++){
			//if(j<32)printf("NT: %i %i %g %g \n", j, i, TotalMatrix[i][j],1.0/Noise[i]);
			NT[i][j]=TotalMatrix[i][j]/Noise[i];
		}
	}


	dgemm(TotalMatrix, NT , TNT, ((MNStruct *)context)->pulse->nobs, totalsize, ((MNStruct *)context)->pulse->nobs, totalsize, 'T', 'N');

	dgemv(NT,Resvec,NTd,((MNStruct *)context)->pulse->nobs,totalsize,'T');


    double jointdet=0;
    double freqlike=0;
    double *WorkCoeff = new double[totalsize];
    for(int o1=0;o1<totalsize; o1++){
    		
            WorkCoeff[o1]=NTd[o1];
    }

	int globalinfo=0;
	int info=0;
	
    dpotrfInfo(TNT, totalsize, jointdet, info);
    if(info != 0)globalinfo=1;

    info=0;
    dpotrsInfo(TNT, WorkCoeff, totalsize, info);
    if(info != 0)globalinfo=1;
    
    for(int j=0;j<totalsize;j++){
            freqlike += NTd[j]*WorkCoeff[j];
            //printf("CPU copy %i %g %g \n", j,NTd[j], WorkCoeff[j]);
    }
	
	lnew=-0.5*(tdet+jointdet+timelike-freqlike)+ equadpriorterm;

	if(isnan(lnew) || isinf(lnew) || globalinfo != 0){

		lnew=-pow(10.0,20);
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);
// 		printf("Like: %g %g %g \n",lnew,Chisq,covdet);
		
	}

	
	delete[] WorkCoeff;
	delete[] EFAC;
	delete[] EQUAD;
	delete[] NTd;
	delete[] Noise;
	delete[] Resvec;
	
	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[]TotalMatrix[j];
		delete[]TNDM[j];
		delete[]NT[j];
	}
	delete[]TotalMatrix;
	delete[]TNDM;
	delete[]NT;

	for (int j = 0; j < totalsize; j++){
		delete[]TNT[j];
	}
	delete[]TNT;


	
//	printf("CPUChisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);

// 	if(isinf(lnew) || isinf(jointdet) || isinf(tdet) || isinf(freqdet) || isinf(timelike) || isinf(freqlike)){
 	//printf("Chisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
// 	}
	


}


void LRedLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;
	double phase=0;
	double *EFAC;
	double *EQUAD;
	int pcount=0;
	int bad=0;
	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];

	for(int p=0;p<ndim;p++){

		Cube[p]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[p]+((MNStruct *)context)->Dpriors[p][0];
	}

	if(((MNStruct *)context)->doLinear==0){
	
		for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
			LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
		}

		phase=(double)LDparams[0];
		pcount++;
		for(int p=1;p<((MNStruct *)context)->numFitTiming;p++){
			((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
			pcount++;
		}
		for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
			((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
			pcount++;
		}
		
		if(((MNStruct *)context)->pulse->param[param_dmmodel].fitFlag[0] == 1){
			int DMnum=((MNStruct *)context)->pulse[0].dmoffsDMnum;
			for(int i =0; i < DMnum; i++){
				((MNStruct *)context)->pulse[0].dmoffsDM[i]=Cube[ndim-DMnum+i];
			}
		}
		
		
		fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
		formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
		
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+phase;
			
		}
// 		printf("Phase: %g \n", phase);
	
	}
	else if(((MNStruct *)context)->doLinear==1){
	
		for(int p=0;p < numfit; p++){
			Fitparams[p]=Cube[p];
			pcount++;
		}
		
		double *Fitvec=new double[((MNStruct *)context)->pulse->nobs];

		dgemv(((MNStruct *)context)->DMatrix,Fitparams,Fitvec,((MNStruct *)context)->pulse->nobs,numfit,'N');
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
		}
		
		delete[] Fitvec;
	}

	if(((MNStruct *)context)->incStep > 0){
		for(int i = 0; i < ((MNStruct *)context)->incStep; i++){
			double StepAmp = Cube[pcount];
			pcount++;
			double StepTime = Cube[pcount];
			pcount++;
			for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
				double time = (double)((MNStruct *)context)->pulse->obsn[o1].bat;
				if(time > StepTime){
					Resvec[o1] += StepAmp;
				}
			}
		}
	}
		

	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EFAC[o]=1;
		}
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EFAC[o]=Cube[pcount];
		}
		pcount++;
		
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int p=0;p< ((MNStruct *)context)->systemcount; p++){
			EFAC[p]=Cube[pcount];
			pcount++;
		}
	}				

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=0;
		}
	}
	else if(((MNStruct *)context)->numFitEQUAD == 1){
		EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=pow(10.0,2*Cube[pcount]);
		}
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEQUAD > 1){
                EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
                        EQUAD[o]=pow(10.0,2*Cube[pcount]);
			pcount++;
                }
        }

	int FitRedCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
	int FitDMCoeff=2*(((MNStruct *)context)->numFitDMCoeff);
	if(((MNStruct *)context)->incFloatDM != 0)FitDMCoeff+=2*((MNStruct *)context)->incFloatDM;
	if(((MNStruct *)context)->incFloatRed != 0)FitRedCoeff+=2*((MNStruct *)context)->incFloatRed;
    	int totCoeff=0;
    	if(((MNStruct *)context)->incRED != 0)totCoeff+=FitRedCoeff;
    	if(((MNStruct *)context)->incDM != 0)totCoeff+=FitDMCoeff;

	double *WorkNoise=new double[((MNStruct *)context)->pulse->nobs];
	double *powercoeff=new double[totCoeff];

	double tdet=0;
	double timelike=0;
	double timelike2=0;

	if(((MNStruct *)context)->whitemodel == 0){
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			
			WorkNoise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD[((MNStruct *)context)->sysFlags[o]];
			
			tdet=tdet+log(WorkNoise[o]);
			WorkNoise[o]=1.0/WorkNoise[o];
			timelike=timelike+pow(Resvec[o],2)*WorkNoise[o];
			timelike2=timelike2+pow((double)((MNStruct *)context)->pulse->obsn[o].residual,2)*WorkNoise[o];

		}
	}
	else if(((MNStruct *)context)->whitemodel == 1){
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			
			WorkNoise[o]=EFAC[((MNStruct *)context)->sysFlags[o]]*EFAC[((MNStruct *)context)->sysFlags[o]]*(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o]]);
			
			tdet=tdet+log(WorkNoise[o]);
			WorkNoise[o]=1.0/WorkNoise[o];
			timelike=timelike+pow(Resvec[o],2)*WorkNoise[o];
			timelike2=timelike2+pow((double)((MNStruct *)context)->pulse->obsn[o].residual,2)*WorkNoise[o];

		}
	}

	
	double *NFd = new double[totCoeff];
	double **FMatrix=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		FMatrix[i]=new double[totCoeff];
	}

	double **NF=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		NF[i]=new double[totCoeff];
	}

	double **FNF=new double*[totCoeff];
	for(int i=0;i<totCoeff;i++){
		FNF[i]=new double[totCoeff];
	}





	double start,end;
	int go=0;
	for (int i=0;i<((MNStruct *)context)->pulse->nobs;i++)
	  {
	    if (((MNStruct *)context)->pulse->obsn[i].deleted==0)
	      {
		if (go==0)
		  {
		    go = 1;
		    start = (double)((MNStruct *)context)->pulse->obsn[i].bat;
		    end  = start;
		  }
		else
		  {
		    if (start > (double)((MNStruct *)context)->pulse->obsn[i].bat)
		      start = (double)((MNStruct *)context)->pulse->obsn[i].bat;
		    if (end < (double)((MNStruct *)context)->pulse->obsn[i].bat)
		      end = (double)((MNStruct *)context)->pulse->obsn[i].bat;
		  }
	      }
	  }
// 	printf("Total time span = %.6f days = %.6f years\n",end-start,(end-start)/365.25);
	double maxtspan=end-start;

        double *freqs = new double[totCoeff];

        double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
        double DMKappa = 2.410*pow(10.0,-16);
        int startpos=0;
        double freqdet=0;
        if(((MNStruct *)context)->incRED==2){
                for (int i=0; i<FitRedCoeff/2; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
                        freqs[i+FitRedCoeff/2]=freqs[i];

                        powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[i+FitRedCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);
                        pcount++;
                }


                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);
                        }
                }

                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i+FitRedCoeff/2]=sin(2*M_PI*freqs[i]*time);
                        }
                }



                startpos=FitRedCoeff;

        }
   else if(((MNStruct *)context)->incRED==3){

                double redamp=Cube[pcount];
                pcount++;
                double redindex=Cube[pcount];
                pcount++;
// 		printf("red: %g %g \n", redamp, redindex);
                 redamp=pow(10.0, redamp);

                freqdet=0;
                 for (int i=0; i<FitRedCoeff/2; i++){

                        freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

                        powercoeff[i]=redamp*redamp*pow((freqs[i]*365.25),-1.0*redindex)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[i+FitRedCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);


                 }

                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);
                        }
                }

                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i+FitRedCoeff/2]=sin(2*M_PI*freqs[i]*time);
                        }
                }


                startpos=FitRedCoeff;

        }


       if(((MNStruct *)context)->incDM==2){

                for (int i=0; i<FitDMCoeff/2; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[startpos+i]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                        pcount++;
                }

                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }

                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                        }
                }

                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i+FitDMCoeff/2]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                        }
                }



        }
        else if(((MNStruct *)context)->incDM==3){
                double DMamp=Cube[pcount];
                pcount++;
                double DMindex=Cube[pcount];
                pcount++;

                DMamp=pow(10.0, DMamp);

                 for (int i=0; i<FitDMCoeff/2; i++){
                        freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[startpos+i]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];

                        powercoeff[startpos+i]=DMamp*DMamp*pow((freqs[startpos+i]*365.25),-1.0*DMindex)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);


                 }
                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }


                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                        }
                }

                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i+FitDMCoeff/2]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                        }
                }



        }





// 	makeFourier(((MNStruct *)context)->pulse, FitCoeff, FMatrix);

	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		for(int j=0;j<totCoeff;j++){
// 			printf("%i %i %g %g \n",i,j,WorkNoise[i],FMatrix[i][j]);
			NF[i][j]=WorkNoise[i]*FMatrix[i][j];
		}
	}
	dgemv(NF,Resvec,NFd,((MNStruct *)context)->pulse->nobs,totCoeff,'T');
	dgemm(FMatrix, NF , FNF, ((MNStruct *)context)->pulse->nobs, totCoeff, ((MNStruct *)context)->pulse->nobs, totCoeff, 'T', 'N');


	double **PPFM=new double*[totCoeff];
	for(int i=0;i<totCoeff;i++){
		PPFM[i]=new double[totCoeff];
		for(int j=0;j<totCoeff;j++){
			PPFM[i][j]=0;
		}
	}


	for(int c1=0; c1<totCoeff; c1++){
		PPFM[c1][c1]=1.0/powercoeff[c1];
	}



	for(int j=0;j<totCoeff;j++){
		for(int k=0;k<totCoeff;k++){
			PPFM[j][k]=PPFM[j][k]+FNF[j][k];
		}
	}

        double jointdet=0;
        double freqlike=0;
       double *WorkCoeff = new double[totCoeff];
       for(int o1=0;o1<totCoeff; o1++){
                WorkCoeff[o1]=NFd[o1];
        }


	    int info=0;
        dpotrfInfo(PPFM, totCoeff, jointdet,info);
        dpotrs(PPFM, WorkCoeff, totCoeff);
        for(int j=0;j<totCoeff;j++){
                freqlike += NFd[j]*WorkCoeff[j];
        }
	
	lnew=-0.5*(jointdet+tdet+freqdet+timelike-freqlike);

	if(isnan(lnew) || isinf(lnew)){

		lnew=-pow(10.0,200);
		
	}
	
	//printf("CPULike: %g %g %g %g %g %g \n", lnew, jointdet, tdet, freqdet, timelike, freqlike);
	delete[] DMVec;
	delete[] WorkCoeff;
	delete[] EFAC;
	delete[] EQUAD;
	delete[] WorkNoise;
	delete[] powercoeff;
	delete[] Resvec;
	delete[] NFd;
	delete[] freqs;

	for (int j = 0; j < totCoeff; j++){
		delete[] PPFM[j];
	}
	delete[] PPFM;

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[] NF[j];
	}
	delete[] NF;

	for (int j = 0; j < totCoeff; j++){
		delete[] FNF[j];
	}
	delete[] FNF;

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[] FMatrix[j];
	}
	delete[] FMatrix;

}



void NewLRedMarginLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double **EFAC;
	double *EQUAD;
	int pcount=0;
	
	int totdims = ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	totdims += ((MNStruct *)context)->numFitEFAC*((MNStruct *)context)->EPolTerms + ((MNStruct *)context)->numFitEQUAD;
	totdims +=2*((MNStruct *)context)->incStep;
	if(((MNStruct *)context)->incRED==2)totdims+=((MNStruct *)context)->numFitRedCoeff;
	if(((MNStruct *)context)->incDM==2)totdims+=((MNStruct *)context)->numFitDMCoeff;
	if(((MNStruct *)context)->varyRedCoeff==1)totdims+=2;
	if(((MNStruct *)context)->varyDMCoeff==1)totdims+=2;
   	if(((MNStruct *)context)->incRED==3)totdims+=2;
    if(((MNStruct *)context)->incDM==3)totdims+=2;
	if(((MNStruct *)context)->yearlyDM == 1)totdims+=2;
 	if(((MNStruct *)context)->incsinusoid == 1)totdims+=3;   
   	if(((MNStruct *)context)->incGWB == 1)totdims+=1; 
	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];
	int fitcount=0;
	for(int p=0;p<totdims;p++){

		if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){
		Cube[pcount]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[pcount]+((MNStruct *)context)->Dpriors[p][0];

		pcount++;
		}
	}
	pcount=0;
	if(((MNStruct *)context)->doLinear==0){
	
		for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
			if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){

				LDparams[p]=Cube[fitcount]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
				fitcount++;
			}
			else if(((MNStruct *)context)->Dpriors[p][1] == ((MNStruct *)context)->Dpriors[p][0]){
				LDparams[p]=((MNStruct *)context)->Dpriors[p][0]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
			}


		}
		pcount=0;
		double phase=(double)LDparams[0];
		pcount++;
		for(int p=1;p<((MNStruct *)context)->numFitTiming;p++){
			((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
			pcount++;
		}
		for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
			((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
			pcount++;
		}
		
		fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
		formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
		
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+phase;

		}
	
	}
	else if(((MNStruct *)context)->doLinear==1){
		fitcount=0;

		
		for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
			if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){
				Fitparams[p]=Cube[fitcount];
				fitcount++;
			}
			else if(((MNStruct *)context)->Dpriors[p][1] == ((MNStruct *)context)->Dpriors[p][0]){
				Fitparams[p]=0;
			}
		}
		double *Fitvec=new double[((MNStruct *)context)->pulse->nobs];

		dgemv(((MNStruct *)context)->DMatrix,Fitparams,Fitvec,((MNStruct *)context)->pulse->nobs,numfit,'N');
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
		}
		
		delete[] Fitvec;
	}
	pcount=fitcount;

	if(((MNStruct *)context)->incStep > 0){
		for(int i = 0; i < ((MNStruct *)context)->incStep; i++){
			double StepAmp = Cube[pcount];
			pcount++;
			double StepTime = Cube[pcount];
			pcount++;
			for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
				if(((MNStruct *)context)->pulse->obsn[o1].bat > StepTime){
					Resvec[o1] += StepAmp;
				}
			}
		}
	}	
	
	
/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get White Noise vector///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////  


	double equadpriorterm=0;
	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double*[((MNStruct *)context)->EPolTerms];
		for(int n=1; n <=((MNStruct *)context)->EPolTerms; n++){
			EFAC[n-1]=new double[((MNStruct *)context)->systemcount];
			if(n==1){
				for(int o=0;o<((MNStruct *)context)->systemcount; o++){
					EFAC[n-1][o]=1;
				}
			}
			else{
                                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
                                       EFAC[n-1][o]=0;
                            }
			}
		}
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double*[((MNStruct *)context)->EPolTerms];
		for(int n=1; n <=((MNStruct *)context)->EPolTerms; n++){
			
			EFAC[n-1]=new double[((MNStruct *)context)->systemcount];
			if(n==1){
				for(int o=0;o<((MNStruct *)context)->systemcount; o++){
					
					EFAC[n-1][o]=Cube[pcount];
				}
				pcount++;
			}
			else{
                    for(int o=0;o<((MNStruct *)context)->systemcount; o++){

                            EFAC[n-1][o]=pow(10.0,Cube[pcount]);
                    }
                    pcount++;
            }
		}
		
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double*[((MNStruct *)context)->EPolTerms];
		for(int n=1; n <=((MNStruct *)context)->EPolTerms; n++){
			EFAC[n-1]=new double[((MNStruct *)context)->systemcount];
			if(n==1){
				for(int p=0;p< ((MNStruct *)context)->systemcount; p++){
					EFAC[n-1][p]=Cube[pcount];
					pcount++;
				}
			}
			else{
                                for(int p=0;p< ((MNStruct *)context)->systemcount; p++){
                                        EFAC[n-1][p]=pow(10.0,Cube[pcount]);
                                        pcount++;
                                }
                        }
		}
	}	

		

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=0;
		}
	}
	else if(((MNStruct *)context)->numFitEQUAD == 1){
		EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=pow(10.0,2*Cube[pcount]);
			equadpriorterm+=log(pow(10.0,Cube[pcount]));
		}
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEQUAD > 1){
        EQUAD=new double[((MNStruct *)context)->systemcount];
        for(int o=0;o<((MNStruct *)context)->systemcount; o++){
            EQUAD[o]=pow(10.0,2*Cube[pcount]);
	    equadpriorterm+=log(pow(10.0,Cube[pcount]));

			pcount++;
        }
    }




	double *Noise;	
	double *BATvec;
	Noise=new double[((MNStruct *)context)->pulse->nobs];
	BATvec=new double[((MNStruct *)context)->pulse->nobs];
	
	
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
	}
		
		
	if(((MNStruct *)context)->whitemodel == 0){
	
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			double EFACterm=0;
			for(int n=1; n <=((MNStruct *)context)->EPolTerms; n++){
				EFACterm=EFACterm + pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))/pow(pow(10.0,-7),n-1),n)*EFAC[n-1][((MNStruct *)context)->sysFlags[o]];
			}	

			Noise[o]= (pow(EFACterm,2) + EQUAD[((MNStruct *)context)->sysFlags[o]]);

		}
		
	}
	else if(((MNStruct *)context)->whitemodel == 1){
	
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){

			Noise[o]=(EFAC[0][((MNStruct *)context)->sysFlags[o]]*EFAC[0][((MNStruct *)context)->sysFlags[o]]*(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o]]));
		}
		
	}

	
	
/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get Time domain likelihood//////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////  


	double tdet=0;
	double timelike=0;

	for(int o=0; o<((MNStruct *)context)->pulse->nobs; o++){
		//printf("tlike %i %g %g %g \n", o, Resvec[o], Noise[o], timelike);
		timelike+=Resvec[o]*Resvec[o]/Noise[o];
		tdet += log(Noise[o]);
	}
		
		
/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Form the Design Matrix////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////  

	int TimetoMargin=0;
	for(int i =0; i < ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps; i++){
		if(((MNStruct *)context)->LDpriors[i][2]==1)TimetoMargin++;
	}


	double **TNDM=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		TNDM[i]=new double[TimetoMargin];
	}

	getCustomDMatrixLike(context, TNDM);
	

	double* S = new double[TimetoMargin];
	double** U = new double*[((MNStruct *)context)->pulse->nobs];
	for(int k=0; k < ((MNStruct *)context)->pulse->nobs; k++){
		U[k] = new double[((MNStruct *)context)->pulse->nobs];
	}
	double** VT = new double*[TimetoMargin]; 
	for (int k=0; k<TimetoMargin; k++) VT[k] = new double[TimetoMargin];

	dgesvd(TNDM,((MNStruct *)context)->pulse->nobs, TimetoMargin, S, U, VT);

	delete[]S;	

	for (int j = 0; j < TimetoMargin; j++){
		delete[]VT[j];
	}
	
	delete[]VT;
	
		
	

	for(int j=0;j<((MNStruct *)context)->pulse->nobs;j++){
		for(int k=0;k < TimetoMargin;k++){
				TNDM[j][k]=U[j][k];
		}
	}

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[]U[j];
	}
	delete[]U;


//////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////////////Form the Power Spectrum//////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////  

	int FitRedCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
	int FitDMCoeff=2*(((MNStruct *)context)->numFitDMCoeff);

	if(((MNStruct *)context)->incFloatDM != 0)FitDMCoeff+=2*((MNStruct *)context)->incFloatDM;
	if(((MNStruct *)context)->incFloatRed != 0)FitRedCoeff+=2*((MNStruct *)context)->incFloatRed;


    int totCoeff=0;
    if(((MNStruct *)context)->incRED != 0)totCoeff+=FitRedCoeff;
    if(((MNStruct *)context)->incDM != 0)totCoeff+=FitDMCoeff;

	double *powercoeff=new double[totCoeff];
	for(int o=0;o<totCoeff; o++){
		powercoeff[o]=0;
	}
	
	double **FMatrix=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		FMatrix[i]=new double[totCoeff];
	}


	double priorterm=0;
	bool uniformprior=0;
	double start,end;
	int go=0;
	for (int i=0;i<((MNStruct *)context)->pulse->nobs;i++)
	  {
	    if (((MNStruct *)context)->pulse->obsn[i].deleted==0)
	      {
		if (go==0)
		  {
		    go = 1;
		    start = (double)((MNStruct *)context)->pulse->obsn[i].bat;
		    end  = start;
		  }
		else
		  {
		    if (start > (double)((MNStruct *)context)->pulse->obsn[i].bat)
		      start = (double)((MNStruct *)context)->pulse->obsn[i].bat;
		    if (end < (double)((MNStruct *)context)->pulse->obsn[i].bat)
		      end = (double)((MNStruct *)context)->pulse->obsn[i].bat;
		  }
	      }
	  }

	double maxtspan=1*(end-start);


   double *freqs = new double[totCoeff];

    double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
    double DMKappa = 2.410*pow(10.0,-16);
    int startpos=0;
    double freqdet=0;
    double GWBAmpPrior=0;
    if(((MNStruct *)context)->incRED==2){

    
		for (int i=0; i<FitRedCoeff/2; i++){
			int pnum=pcount;
			double pc=Cube[pcount];
			freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
			freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

			powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);
			powercoeff[i+FitRedCoeff/2]=powercoeff[i];
			freqdet=freqdet+2*log(powercoeff[i]);
			pcount++;
		}
		
		
        for(int i=0;i<FitRedCoeff/2;i++){
                for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                        double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                        FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);

                }
        }

        for(int i=0;i<FitRedCoeff/2;i++){
                for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                        double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                        FMatrix[k][i+FitRedCoeff/2]=sin(2*M_PI*freqs[i]*time);
                }
        }
	            
	    startpos=FitRedCoeff;

    }
   else if(((MNStruct *)context)->incRED==3){

		freqdet=0;
		
		for(int pl = 0; pl < ((MNStruct *)context)->numFitRedPL; pl ++){
			double redamp=Cube[pcount];
			pcount++;
			double redindex=Cube[pcount];
			pcount++;
	
   			double Tspan = maxtspan;
			double f1yr = 1.0/3.16e7;
    			
			
			redamp=pow(10.0, redamp);
			if(uniformprior==1)priorterm+=log(redamp);


			double Agw=redamp;
	
			for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed ; i++){
	
				freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
				freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
				
 				double rho = (Agw*Agw/12.0/(M_PI*M_PI))*pow(f1yr,(-3)) * pow(freqs[i]*365.25,(-redindex))/(maxtspan*24*60*60);

 				powercoeff[i]+= rho;
 				powercoeff[i+FitRedCoeff/2]+= rho;

			}
		}
		
		int coefftovary=0;
		double amptovary=0.0;
		if(((MNStruct *)context)->varyRedCoeff==1){
			coefftovary=int(pow(10.0,Cube[pcount]))-1;
			pcount++;
			amptovary=pow(10.0,Cube[pcount])/(maxtspan*24*60*60);
			pcount++;

			powercoeff[coefftovary]=amptovary;
			powercoeff[coefftovary+FitRedCoeff/2]=amptovary;	
		}		
		
		double GWBAmp=0;
		if(((MNStruct *)context)->incGWB==1){
			GWBAmp=pow(10.0,Cube[pcount]);
			pcount++;
			GWBAmpPrior=log(GWBAmp);
			double Tspan = maxtspan;
			double f1yr = 1.0/3.16e7;
			for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed ; i++){
				double rho = (GWBAmp*GWBAmp/12.0/(M_PI*M_PI))*pow(f1yr,(-3)) * pow(freqs[i]*365.25,(-4.333))/(maxtspan*24*60*60);	
				powercoeff[i]+= rho;
				powercoeff[i+FitRedCoeff/2]+= rho;
			}
		}
		for (int i=0; i<FitRedCoeff/2; i++){
			freqdet=freqdet+2*log(powercoeff[i]);
		}
		
        for(int i=0;i<FitRedCoeff/2;i++){
                for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                        double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                        FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);

                }
        }

        for(int i=0;i<FitRedCoeff/2;i++){
                for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                        double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                        FMatrix[k][i+FitRedCoeff/2]=sin(2*M_PI*freqs[i]*time);
                }
        }

        startpos=FitRedCoeff;

    }


	if(((MNStruct *)context)->incsinusoid == 1){
		double sineamp=pow(10.0,Cube[pcount]);
		pcount++;
		double sinephase=Cube[pcount];
		pcount++;
		double sinefreq=pow(10.0,Cube[pcount])/maxtspan;
		pcount++;		
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]-= sineamp*sin(2*M_PI*sinefreq*(double)((MNStruct *)context)->pulse->obsn[o].bat + sinephase);
		}
	}





       if(((MNStruct *)context)->incDM==2){

			for (int i=0; i<FitDMCoeff/2; i++){
				int pnum=pcount;
				double pc=Cube[pcount];
				freqs[startpos+i]=((MNStruct *)context)->sampleFreq[startpos/2 - ((MNStruct *)context)->incFloatRed+i]/maxtspan;
				freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
	
				powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
				powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
				freqdet=freqdet+2*log(powercoeff[startpos+i]);
				pcount++;
			}
           	 


                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }
                
            for(int i=0;i<FitDMCoeff/2;i++){
                for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                        double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                        FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                }
        }

        for(int i=0;i<FitDMCoeff/2;i++){
                for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                        double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                        FMatrix[k][startpos+i+FitDMCoeff/2]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                }
        }

        }
        else if(((MNStruct *)context)->incDM==3){

		for(int pl = 0; pl < ((MNStruct *)context)->numFitDMPL; pl ++){
			double DMamp=Cube[pcount];
			pcount++;
			double DMindex=Cube[pcount];
			pcount++;
			
   			double Tspan = maxtspan;
			double f1yr = 1.0/3.16e7;
    			

			DMamp=pow(10.0, DMamp);
			if(uniformprior==1)priorterm+=log(DMamp);
			for (int i=0; i<FitDMCoeff/2; i++){
	
				freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[startpos/2 - ((MNStruct *)context)->incFloatRed +i]/maxtspan;
				freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
				
 				double rho = (DMamp*DMamp)*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-DMindex))/(maxtspan*24*60*60);	
 				powercoeff[startpos+i]+=rho;
 				powercoeff[startpos+i+FitDMCoeff/2]+=rho;
			}
		}
		
		
		int coefftovary=0;
		double amptovary=0.0;
		if(((MNStruct *)context)->varyDMCoeff==1){
			coefftovary=int(pow(10.0,Cube[pcount]))-1;
			pcount++;
			amptovary=pow(10.0,Cube[pcount])/(maxtspan*24*60*60);
			pcount++;

			powercoeff[startpos+coefftovary]=amptovary;
			powercoeff[startpos+coefftovary+FitDMCoeff/2]=amptovary;	
		}	
			
		
		for (int i=0; i<FitDMCoeff/2; i++){
			freqdet=freqdet+2*log(powercoeff[startpos+i]);
		}


        for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
        }

    }
    
        for(int i=0;i<FitDMCoeff/2;i++){
                for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                        double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                        FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                }
        }

        for(int i=0;i<FitDMCoeff/2;i++){
                for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                        double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                        FMatrix[k][startpos+i+FitDMCoeff/2]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                }
        }
    
	if(((MNStruct *)context)->yearlyDM == 1){
		double yearlyamp=pow(10.0,Cube[pcount]);
		pcount++;
		double yearlyphase=Cube[pcount];
		pcount++;
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]-= yearlyamp*sin((2*M_PI/365.25)*(double)((MNStruct *)context)->pulse->obsn[o].bat + yearlyphase)*DMVec[o];
		}
	}


//////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////Form Total Matrices////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

	int totalsize=TimetoMargin+totCoeff;
	double **TotalMatrix=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i =0;i<((MNStruct *)context)->pulse->nobs;i++){
		TotalMatrix[i]=new double[totalsize];
		for(int j =0;j<totalsize; j++){
			TotalMatrix[i][j]=0;
		}
	}
	

	for(int i =0;i<((MNStruct *)context)->pulse->nobs;i++){
		for(int j =0;j<TimetoMargin; j++){
			TotalMatrix[i][j]=TNDM[i][j];
		}
		
		for(int j =0;j<totCoeff; j++){
			TotalMatrix[i][j+TimetoMargin]=FMatrix[i][j];
		}
	}

//////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////Do Algebra/////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

	double *NTd = new double[totalsize];
	double **NT=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		NT[i]=new double[totalsize];
	}

	double **TNT=new double*[totalsize];
	for(int i=0;i<totalsize;i++){
		TNT[i]=new double[totalsize];
	}
	
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		for(int j=0;j<totalsize;j++){
			//if(j<32)printf("NT: %i %i %g %g \n", j, i, TotalMatrix[i][j],1.0/Noise[i]);
			NT[i][j]=TotalMatrix[i][j]/Noise[i];
		}
	}


	dgemm(TotalMatrix, NT , TNT, ((MNStruct *)context)->pulse->nobs, totalsize, ((MNStruct *)context)->pulse->nobs, totalsize, 'T', 'N');

	dgemv(NT,Resvec,NTd,((MNStruct *)context)->pulse->nobs,totalsize,'T');

	for(int j=0;j<totCoeff;j++){
			TNT[TimetoMargin+j][TimetoMargin+j] += 1.0/powercoeff[j];
	}

    double jointdet=0;
    double freqlike=0;
    double *WorkCoeff = new double[totalsize];
    for(int o1=0;o1<totalsize; o1++){
    		
            WorkCoeff[o1]=NTd[o1];
    }

	int globalinfo=0;
	int info=0;
	
    dpotrfInfo(TNT, totalsize, jointdet, info);
    if(info != 0)globalinfo=1;

    info=0;
    dpotrsInfo(TNT, WorkCoeff, totalsize, info);
    if(info != 0)globalinfo=1;
    
    for(int j=0;j<totalsize;j++){
            freqlike += NTd[j]*WorkCoeff[j];
            //printf("CPU copy %i %g %g \n", j,NTd[j], WorkCoeff[j]);
    }
	
	lnew=-0.5*(tdet+jointdet+freqdet+timelike-freqlike)+ equadpriorterm + GWBAmpPrior;

	if(isnan(lnew) || isinf(lnew) || globalinfo != 0){

		lnew=-pow(10.0,20);
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);
// 		printf("Like: %g %g %g \n",lnew,Chisq,covdet);
		
	}

	delete[] DMVec;
	delete[] WorkCoeff;
	delete[] EFAC;
	delete[] EQUAD;
	delete[] powercoeff;
	delete[] NTd;
	delete[] freqs;
	delete[] Noise;
	delete[] Resvec;
	
	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[]TotalMatrix[j];
		delete[]TNDM[j];
		delete[]FMatrix[j];
		delete[]NT[j];
	}
	delete[]TotalMatrix;
	delete[]TNDM;
	delete[]FMatrix;
	delete[]NT;

	for (int j = 0; j < totalsize; j++){
		delete[]TNT[j];
	}
	delete[]TNT;


	
//	printf("CPUChisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);

// 	if(isinf(lnew) || isinf(jointdet) || isinf(tdet) || isinf(freqdet) || isinf(timelike) || isinf(freqlike)){
 	//printf("Chisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
// 	}

}



