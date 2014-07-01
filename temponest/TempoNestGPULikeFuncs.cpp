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

double *TNFactorialList = new double[21];



extern "C" void WhiteMarginGPUWrapper_(void *context, double *TNDMVec, double *resvec, double *BATvec, double *Noise, int N, int D, int T,int NTime, int NJumps, double *likevals);



extern "C" void LRedGPUWrapper_(double *Freqs, double *resvec, double *BATvec, double *DMVec, double *Noise, double **FNF, double *NFd, int N, int RF, int DMF, int F, int incRED, int incDM);
extern "C" void NewLRedMarginGPUWrapper_(void *context, double *TNDMVec, double *Freqs, double *powercoeff, double *resvec, double *BATvec, double *DMVec, double *Noise, int N, int RF,int DMF, int D, int F, int T, int incRED, int incDM, int numFitTiming, int numFitJumps, double *likevals);




void store_factorial(){

	for(int k=0; k <=20; k++){
		TNFactorialList[k]=iter_factorial(k);
	}
}

void WhiteMarginGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double **EFAC;
	double *EQUAD;
	int pcount=0;
	
	int totdims = ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	totdims += ((MNStruct *)context)->numFitEFAC*((MNStruct *)context)->EPolTerms + ((MNStruct *)context)->numFitEQUAD;
	totdims +=2*((MNStruct *)context)->incStep;

	if(((MNStruct *)context)->yearlyDM == 1)totdims+=2;
 	if(((MNStruct *)context)->incsinusoid == 1)totdims+=3;   

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

			Noise[o]= 1.0/(pow(EFACterm,2) + EQUAD[((MNStruct *)context)->sysFlags[o]]);

		}
		
	}
	else if(((MNStruct *)context)->whitemodel == 1){
	
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){

			Noise[o]=1.0/(EFAC[0][((MNStruct *)context)->sysFlags[o]]*EFAC[0][((MNStruct *)context)->sysFlags[o]]*(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o]]));
		}
		
	}
	


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
			Resvec[o]-= yearlyamp*sin((2*M_PI/365.25)*(double)((MNStruct *)context)->pulse->obsn[o].bat + yearlyphase)/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
		}
	}


	
/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get Time domain likelihood//////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////  


	double tdet=0;
	double timelike=0;

	for(int o=0; o<((MNStruct *)context)->pulse->nobs; o++){
		timelike+=Resvec[o]*Resvec[o]*Noise[o];
		tdet -= log(Noise[o]);
	}

    
//////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////get TNDMVec////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////   

    
	int TimetoMargin=0;
	for(int i =0; i < ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps; i++){
		if(((MNStruct *)context)->LDpriors[i][2]==1)TimetoMargin++;
	}
	
	double **TNDM=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		TNDM[i]=new double[TimetoMargin];
	}
	double *TNDMVec=new double[((MNStruct *)context)->pulse->nobs*TimetoMargin];
	
	if(TimetoMargin != ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps){
	
		getCustomDMatrixLike(context, TNDM);
	
		for(int g=0;g<TimetoMargin; g++){
			for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){

				TNDMVec[g*((MNStruct *)context)->pulse->nobs + o]=TNDM[o][g];
			}
		}
	}
	

    
    

//////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////Call GPU Code/////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////    
    
    
	double *likevals=new double[2];
	int totalsize=TimetoMargin;
	WhiteMarginGPUWrapper_(context, TNDMVec, Resvec, BATvec, Noise,  ((MNStruct *)context)->pulse->nobs, TimetoMargin, totalsize, ((MNStruct *)context)->numFitTiming, ((MNStruct *)context)->numFitJumps, likevals);
    
    
//////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////calculate likelihood///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////    

	double jointdet=likevals[0];
	double freqlike=likevals[1];
    
	
	lnew=-0.5*(tdet+jointdet+timelike-freqlike) + equadpriorterm;

	if(isnan(lnew) || isinf(lnew)){

		lnew=-pow(10.0,20);
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);
// 		printf("Like: %g %g %g \n",lnew,Chisq,covdet);
		
	}


	delete[] EFAC;
	delete[] EQUAD;
	delete[] Noise;
	delete[] Resvec;
	delete[] BATvec;
	delete[] likevals;
	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[]TNDM[j];
	}
	delete[]TNDM;
	delete[]TNDMVec;

	
	//printf("GPUChisq: %.8g %.8g %.8g %.8g %.8g %.8g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);

// 	if(isinf(lnew) || isinf(jointdet) || isinf(tdet) || isinf(freqdet) || isinf(timelike) || isinf(freqlike)){
 	//printf("Chisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
// 	}


	

}



void LRedGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
	//printf("hereNM");
	clock_t startClock,endClock;

	double *EFAC;
	double *EQUAD;
	int pcount=0;

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

        int totCoeff=0;
        if(((MNStruct *)context)->incRED != 0)totCoeff+=FitRedCoeff;
        if(((MNStruct *)context)->incDM != 0)totCoeff+=FitDMCoeff;

        double *powercoeff=new double[totCoeff];
        for(int o=0;o<totCoeff; o++){
                powercoeff[o]=0;
        }

	double *WorkNoise=new double[((MNStruct *)context)->pulse->nobs];

	double tdet=0;
	double timelike=0;



	double *BATvec=new double[((MNStruct *)context)->pulse->nobs];

	if(((MNStruct *)context)->whitemodel == 0){
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			//	printf("Noise %i %g %g %g\n",m1,Noise[m1],EFAC[flagList[m1]],EQUAD);
				WorkNoise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD[((MNStruct *)context)->sysFlags[o]];
				
				tdet=tdet+log(WorkNoise[o]);
				WorkNoise[o]=1.0/WorkNoise[o];
				timelike=timelike+pow(Resvec[o],2)*WorkNoise[o];
	
				BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
	
	
		}
	}
	else if(((MNStruct *)context)->whitemodel == 1){
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			//	printf("Noise %i %g %g %g\n",m1,Noise[m1],EFAC[flagList[m1]],EQUAD);
				WorkNoise[o]=EFAC[((MNStruct *)context)->sysFlags[o]]*EFAC[((MNStruct *)context)->sysFlags[o]]*(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o]]);
				
				tdet=tdet+log(WorkNoise[o]);
				WorkNoise[o]=1.0/WorkNoise[o];
				timelike=timelike+pow(Resvec[o],2)*WorkNoise[o];
	
				BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
	
	
		}
	}

	double *NFd = new double[totCoeff];
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

	double maxtspan=end-start;


	double *freqs = new double[totCoeff];

        double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
        double DMKappa = 2.410*pow(10.0,-16);
        int startpos=0;
        double freqdet=0;
	

        if(((MNStruct *)context)->incRED==2){
        
        	if(((MNStruct *)context)->incFloatRed == 0){
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
            }
            else if(((MNStruct *)context)->incFloatRed >0){

                for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed ; i++){
                
                		int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

                        powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[i+FitRedCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);
                        pcount++;
                }
                
                 for (int i=FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed; i<FitRedCoeff/2; i++){
                		//printf("Freq: %g \n", Cube[pcount]);
                		
                        freqs[startpos+i]=Cube[pcount]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
                         pcount++;
                         
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        pcount++;
                        

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitRedCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                       
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

                
   			   redamp=pow(10.0, redamp);

				for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed ; i++){

		                freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
		                freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
		                
		                double PLcomp=redamp*redamp*pow((freqs[i]*365.25),-1.0*redindex)/(maxtspan*24*60*60);

		                powercoeff[i]+= PLcomp;
		                powercoeff[i+FitRedCoeff/2]+= PLcomp;
				}
			}
				
				
			for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed ; i++){
				freqdet=freqdet+2*log(powercoeff[i]);
		//		printf("%i %g %g \n",i,powercoeff[i], freqdet);
			}

                 for (int i=FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed; i<FitRedCoeff/2; i++){
                                
                    //    Cube[pcount]=floor(Cube[pcount]);
                        freqs[startpos+i]=Cube[pcount]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
                        pcount++;

                        int pnum=pcount;
                        double pc=Cube[pcount];
                        pcount++;


                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitRedCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);

                }

                startpos=FitRedCoeff;

        }
// 		printf("DM\n");
        double nlist[((MNStruct *)context)->incFloatDM][2];
	if(((MNStruct *)context)->incDM==2){
        
        	if(((MNStruct *)context)->incFloatDM == 0){

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
            }
            else if(((MNStruct *)context)->incFloatDM >0){

                for (int i=0; i<FitDMCoeff/2 - ((MNStruct *)context)->incFloatDM ; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=((MNStruct *)context)->sampleFreq[startpos/2 - ((MNStruct *)context)->incFloatRed +i]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                        pcount++;
                }
                
                for (int i=FitDMCoeff/2 - ((MNStruct *)context)->incFloatDM; i<FitDMCoeff/2; i++){

                        freqs[startpos+i]=Cube[pcount]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];

                         pcount++;
                         
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        pcount++;
                        
						
                        
                        

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                       
                }




                
            }            

                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }


        }
        else if(((MNStruct *)context)->incDM==3){
        
        
				for(int pl = 0; pl < ((MNStruct *)context)->numFitDMPL; pl ++){
			
		            double DMamp=Cube[pcount];
		            pcount++;
		            double DMindex=Cube[pcount];
		            pcount++;

		            DMamp=pow(10.0, DMamp);

					for (int i=0; i<FitDMCoeff/2 - ((MNStruct *)context)->incFloatDM ; i++){

				            freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[startpos/2 - ((MNStruct *)context)->incFloatRed +i]/maxtspan;
		                    freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
		                    
		                     double PLcomp=DMamp*DMamp*pow((freqs[startpos+i]*365.25),-1.0*DMindex)/(maxtspan*24*60*60);
				            
				            powercoeff[startpos+i]+=PLcomp;
		                    powercoeff[startpos+i+FitDMCoeff/2]+=PLcomp;
					}
				}
			
				for (int i=0; i<FitDMCoeff/2 - ((MNStruct *)context)->incFloatDM ; i++){
					freqdet=freqdet+2*log(powercoeff[startpos+i]);
//					printf("%i %g %g \n", i, powercoeff[startpos+i], freqdet);
				}
			



                 for (int i= FitDMCoeff/2 - ((MNStruct *)context)->incFloatDM ; i<FitDMCoeff/2; i++){

						//Cube[pcount]=floor(Cube[pcount]);
                        freqs[startpos+i]=Cube[pcount]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
                        pcount++;

                        int pnum=pcount;
                        double pc=Cube[pcount];
                        pcount++;


                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);

                }

                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }

        }
	LRedGPUWrapper_(freqs, Resvec, BATvec, DMVec, WorkNoise, FNF, NFd, ((MNStruct *)context)->pulse->nobs, FitRedCoeff, FitDMCoeff, totCoeff,((MNStruct *)context)->incRED, ((MNStruct *)context)->incDM);



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
        dpotrfInfo(PPFM, totCoeff, jointdet, info);
        dpotrs(PPFM, WorkCoeff, totCoeff);
        for(int j=0;j<totCoeff;j++){
                freqlike += NFd[j]*WorkCoeff[j];
        }
        lnew=-0.5*(((double)((MNStruct *)context)->pulse->nobs)*log(2.0*M_PI) + tdet+jointdet+freqdet+timelike-freqlike);


	if(isnan(lnew) || isinf(lnew)){

		lnew=-pow(10.0,200);
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);

		
	}
 	//printf("Like: %g %g %g %g %g %g\n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
	 //printf("CPULIKE: %g %g %g %g %g %g \n", lnew, jointdet,tdet,freqdet,timelike,freqlike);

	delete[] EFAC;
	delete[] EQUAD;
	delete[] WorkNoise;
	delete[] powercoeff;
	delete[] Resvec;
	delete[] BATvec;
	delete[] NFd;
	delete[] freqs;
	delete[] DMVec;
	delete[] WorkCoeff;

	for (int j = 0; j < totCoeff; j++){
		delete[] PPFM[j];
	}
	delete[] PPFM;



	for (int j = 0; j < totCoeff; j++){
		delete[] FNF[j];
	}
	delete[] FNF;


}




void NewLRedMarginGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
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

			Noise[o]= 1.0/(pow(EFACterm,2) + EQUAD[((MNStruct *)context)->sysFlags[o]]);

		}
		
	}
	else if(((MNStruct *)context)->whitemodel == 1){
	
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){

			Noise[o]=1.0/(EFAC[0][((MNStruct *)context)->sysFlags[o]]*EFAC[0][((MNStruct *)context)->sysFlags[o]]*(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o]]));
		}
		
	}
	


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
    
	if(((MNStruct *)context)->yearlyDM == 1){
		double yearlyamp=pow(10.0,Cube[pcount]);
		pcount++;
		double yearlyphase=Cube[pcount];
		pcount++;
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]-= yearlyamp*sin((2*M_PI/365.25)*(double)((MNStruct *)context)->pulse->obsn[o].bat + yearlyphase)*DMVec[o];
		}
	}


	
/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get Time domain likelihood//////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////  


	double tdet=0;
	double timelike=0;

	for(int o=0; o<((MNStruct *)context)->pulse->nobs; o++){
		timelike+=Resvec[o]*Resvec[o]*Noise[o];
		tdet -= log(Noise[o]);
	}

    
//////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////get TNDMVec////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////   

    
	int TimetoMargin=0;
	for(int i =0; i < ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps; i++){
		if(((MNStruct *)context)->LDpriors[i][2]==1)TimetoMargin++;
	}
	
	double **TNDM=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		TNDM[i]=new double[TimetoMargin];
	}
	double *TNDMVec=new double[((MNStruct *)context)->pulse->nobs*TimetoMargin];
	
	if(TimetoMargin != ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps){
	
		getCustomDMatrixLike(context, TNDM);
	
		for(int g=0;g<TimetoMargin; g++){
			for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){

				TNDMVec[g*((MNStruct *)context)->pulse->nobs + o]=TNDM[o][g];
			}
		}
	}
	

    
    

//////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////Call GPU Code/////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////    
    
    
	double *likevals=new double[2];
	int totalsize=TimetoMargin+totCoeff;
	NewLRedMarginGPUWrapper_(context, TNDMVec, freqs, powercoeff, Resvec, BATvec, DMVec, Noise,  ((MNStruct *)context)->pulse->nobs, FitRedCoeff,FitDMCoeff, TimetoMargin, totCoeff, totalsize, ((MNStruct *)context)->incRED,((MNStruct *)context)->incDM, ((MNStruct *)context)->numFitTiming, ((MNStruct *)context)->numFitJumps, likevals);
    
    
//////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////calculate likelihood///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////    

	double jointdet=likevals[0];
	double freqlike=likevals[1];
    
	
	lnew=-0.5*(tdet+jointdet+freqdet+timelike-freqlike) + equadpriorterm + GWBAmpPrior;

	if(isnan(lnew) || isinf(lnew)){

		lnew=-pow(10.0,20);
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);
// 		printf("Like: %g %g %g \n",lnew,Chisq,covdet);
		
	}

	delete[] DMVec;
	delete[] EFAC;
	delete[] EQUAD;
	delete[] powercoeff;
	delete[] freqs;
	delete[] Noise;
	delete[] Resvec;
	delete[] BATvec;
	delete[] likevals;
	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[]TNDM[j];
	}
	delete[]TNDM;
	delete[]TNDMVec;

	
	//printf("GPUChisq: %.8g %.8g %.8g %.8g %.8g %.8g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);

// 	if(isinf(lnew) || isinf(jointdet) || isinf(tdet) || isinf(freqdet) || isinf(timelike) || isinf(freqlike)){
 	//printf("Chisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
// 	}

}




