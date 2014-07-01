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


#include "tempo2.h"
#include "TempoNest.h"
#include "dgemm.h"
#include "dgesvd.h"
#include "dpotrf.h"
#include "dpotri.h"
#include "dgemv.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <iomanip>
#include "/usr/include/gsl/gsl_sf_gamma.h"

void writeCov(std::vector<double> Cube, int &ndim,void *context, std::string longname, int outputoption);

double iter_factorial(unsigned int n)
{
    double ret = 1;
    for(unsigned int i = 1; i <= n; ++i)
        ret *= i;
    return ret;
}


	

void readsummary(pulsar *psr, std::string longname, int ndim, void *context, long double *Tempo2Fit, int incRED, int ndims, int doTimeMargin, int doJumpMargin, int doLinear){

	int number_of_lines = 0;
	int outputoption=2;
	char *outname;
	int pcount=0;
	int fitcount=0;

	std::ifstream checkfile;
	std::string checkname = longname+"summary.txt";
	checkfile.open(checkname.c_str());
	std::string line;
	while (getline(checkfile, line))
		++number_of_lines;

	printf("number of lines %i \n",number_of_lines);
	checkfile.close();
	
	if(number_of_lines < 4){
		std::ifstream summaryfile;
		std::string fname = longname+"summary.txt";
		summaryfile.open(fname.c_str());
	
		for(int i=0;i<1;i++){
			
			std::string line;
			getline(summaryfile,line);
			//printf("line %s \n", line);
			std::istringstream myStream( line );
			std::istream_iterator< double > begin(myStream),eof;
			std::vector<double> paramlist(begin,eof);

			double **paramarray = new double*[ndims];
			for(int p =0;p < ndims; p++){
				paramarray[p]=new double[4];
				for(int pp=0; pp<4; pp++){
					paramarray[p][pp] = paramlist[pp*ndims+p];
				}
			}

			//getmaxlikeDM(psr,longname, ndims, context, paramarray);

			int numlongparams=((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps;
			long double *LDP = new long double[numlongparams];
			pcount=0;
			fitcount=0;
			for(int j=0;j<((MNStruct *)context)->numFitTiming;j++){
				if(((MNStruct *)context)->Dpriors[pcount][0] != ((MNStruct *)context)->Dpriors[pcount][1]){	
					LDP[j]=paramlist[fitcount+outputoption*ndim]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
					fitcount++;
				}
				else if(((MNStruct *)context)->Dpriors[pcount][0] == ((MNStruct *)context)->Dpriors[pcount][1]){  
					LDP[j]=((MNStruct *)context)->Dpriors[pcount][0]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
				}
			//	printf("LD: %.25Lg %g %g \n",LDP[j], ((MNStruct *)context)->Dpriors[pcount][0],((MNStruct *)context)->Dpriors[pcount][1]); 
				pcount++;
			}

			for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){
				if(((MNStruct *)context)->Dpriors[pcount][0] != ((MNStruct *)context)->Dpriors[pcount][1]){ 
					LDP[pcount]=paramlist[fitcount+outputoption*ndim]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
					fitcount++;
				}
				else if(((MNStruct *)context)->Dpriors[pcount][0] == ((MNStruct *)context)->Dpriors[pcount][1]){
					LDP[pcount]=((MNStruct *)context)->Dpriors[pcount][0]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
				}
				pcount++;
			}
			

				
			pcount=0;
			int jj=0;
			if(doLinear==1){
			
				double *linearParams=new double[ndim];
				double *nonLinearParams=new double[ndim];
				double *errorvec=new double[ndim];
				double *nonLinearerror=new double[ndim];
				fitcount=0;
				pcount=0;
	                        for(int j=0;j<((MNStruct *)context)->numFitTiming;j++){
	                                if(((MNStruct *)context)->Dpriors[pcount][0] != ((MNStruct *)context)->Dpriors[pcount][1]){
               		                        linearParams[j]=paramlist[fitcount+outputoption*ndim];
                                       		errorvec[j]=paramlist[fitcount+ndim];
						fitcount++;
                     			 }
                               		else if(((MNStruct *)context)->Dpriors[pcount][0] == ((MNStruct *)context)->Dpriors[pcount][1]){
                             			 linearParams[j]=((MNStruct *)context)->Dpriors[pcount][0];
						errorvec[j]=0;
					}
                               	 	
                        		//      printf("LD: %g %g \n",((MNStruct *)context)->Dpriors[pcount][0],((MNStruct *)context)->Dpriors[pcount][1]); 
                                
					nonLinearerror[pcount]=0;
					nonLinearParams[pcount]=0;
					pcount++;
				}
                       		for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){
                                	if(((MNStruct *)context)->Dpriors[pcount][0] != ((MNStruct *)context)->Dpriors[pcount][1]){
                                        	linearParams[pcount]=paramlist[fitcount+outputoption*ndim];
						errorvec[j]=paramlist[fitcount+ndim];
                  	                	fitcount++;
                        	      	 }
                                	else if(((MNStruct *)context)->Dpriors[pcount][0] == ((MNStruct *)context)->Dpriors[pcount][1]){
                                        	linearParams[pcount]=((MNStruct *)context)->Dpriors[pcount][0];
						errorvec[j]=0;
                               		 }
					nonLinearerror[pcount]=0;
                                        nonLinearParams[pcount]=0;

                                	pcount++;
                       		 }

				
				TNupdateParameters(psr,0,linearParams,errorvec, nonLinearParams);
				TNupdateParameters(psr,0,errorvec,errorvec, nonLinearerror);
                
				for(int j=0;j<((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps;j++){
                        		if(((MNStruct *)context)->LDpriors[j][1] != 0){
					//	 printf("before %i %g %g %Lg\n",j,nonLinearParams[j],nonLinearerror[j],((MNStruct *)context)->LDpriors[j][1]);
						nonLinearParams[j]=nonLinearParams[j]/((MNStruct *)context)->LDpriors[j][1];
						nonLinearerror[j]=nonLinearerror[j]/((MNStruct *)context)->LDpriors[j][1];
						//printf("after %i %g %g %Lg\n",j,nonLinearParams[j],nonLinearerror[j],((MNStruct *)context)->LDpriors[j][1]);
					}
					else{
						nonLinearParams[j]=0;
						nonLinearerror[j]=0;
					}
                		}

				pcount=1;
               			for(int j=1;j<((MNStruct *)context)->numFitTiming;j++){

	 	                       long double value=nonLinearParams[j]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
        		                long double error=nonLinearerror[j]*(((MNStruct *)context)->LDpriors[pcount][1]);
                		        ((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].val[((MNStruct *)context)->TempoFitNums[pcount][1]] = value;
                        		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].err[((MNStruct *)context)->TempoFitNums[pcount][1]] = error;
  	            		  //        printf("%i %Lg %Lg %Lg %Lg \n",j,(((MNStruct *)context)->LDpriors[pcount][0]),(((MNStruct *)context)->LDpriors[pcount][1]),value,error);
        	        	  //      printf("%i %g \n", j, paramlist[j]);
                	    	    pcount++;
                		}

                       		for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){

        	                        long double value=nonLinearParams[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
                	                long double error=nonLinearerror[pcount]*(((MNStruct *)context)->LDpriors[pcount][1]);
                        	        ((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[j]] = value;
                                	((MNStruct *)context)->pulse->jumpValErr[((MNStruct *)context)->TempoJumpNums[j]] = error;
                                	pcount++;
                        	}
			}
			else{ 
			pcount=1;
			fitcount=0;
			if(((MNStruct *)context)->LDpriors[0][2]==0)fitcount++;
			for(int j=1;j<((MNStruct *)context)->numFitTiming;j++){
			
				long double value;
				long double error;
				
				value=LDP[pcount];
	
				if(((MNStruct *)context)->LDpriors[j][2]==0){
					error=paramlist[fitcount+ndim]*(((MNStruct *)context)->LDpriors[pcount][1]);
					fitcount++;
				}
				else if(((MNStruct *)context)->LDpriors[j][2]==1){
					error=0;
				}
	
				
				((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].val[((MNStruct *)context)->TempoFitNums[pcount][1]] = value;
				((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].err[((MNStruct *)context)->TempoFitNums[pcount][1]] = error;
			//	printf("%i %Lg %Lg %Lg %Lg \n",j,(((MNStruct *)context)->LDpriors[pcount][0]),(((MNStruct *)context)->LDpriors[pcount][1]),value,error);
				//printf("%i %.25Lg %.25Lg \n", j, value, error);
				pcount++;
			}
			
			for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){
				
				long double value;
				long double error;
				
				value=LDP[pcount];
				

				if(((MNStruct *)context)->LDpriors[pcount][2]==0){
                                        error=paramlist[fitcount+ndim]*(((MNStruct *)context)->LDpriors[pcount][1]);
                                        fitcount++;
                                }
                                else if(((MNStruct *)context)->LDpriors[pcount][2]==1){
                                        error=0;
				}
				
				((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[j]] = value;
				((MNStruct *)context)->pulse->jumpValErr[((MNStruct *)context)->TempoJumpNums[j]] = error;
			//	printf("%i %Lg %Lg %Lg %Lg \n",pcount,(((MNStruct *)context)->LDpriors[pcount][0]),(((MNStruct *)context)->LDpriors[pcount][1]),value,error);
				pcount++;
			}	
			
			}
			
			//printf("here \n");
				
		 
			  if(((MNStruct *)context)->doLinear==0){
        		        if(((MNStruct *)context)->pulse->param[param_dmmodel].fitFlag[0] == 1){
					int DMnum=((MNStruct *)context)->pulse[0].dmoffsDMnum;
                        		for(int i =0; i < DMnum; i++){
						printf("update DMoff: %i %g \n", i, paramlist[outputoption*ndim + ndim-DMnum+i]);
                                		((MNStruct *)context)->pulse[0].dmoffsDM[i]=paramlist[outputoption*ndim + ndim-DMnum+i];
                       			}
                		}
			}
			
				
	
			formBatsAll(((MNStruct *)context)->pulse,1);           // Form Barycentric arrival times 
			//printf("formed bats \n");
			formResiduals(((MNStruct *)context)->pulse,1,1);       //Form residuals 
			//printf("done bats and stuff \n");	
			std::ofstream designfile;
			std::string dname = longname+"T2scaling.txt";
			printf("Writing Timing Model Design Matrix and Residuals\n");
			/*
			designfile.open(dname.c_str());
			double pdParamDeriv[MAX_PARAMS];
			int numtofit=((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps;
			for(int i =1;i<((MNStruct *)context)->numFitTiming; i++){
			designfile << psr->param[((MNStruct *)context)->TempoFitNums[i][0]].label[((MNStruct *)context)->TempoFitNums[pcount][1]];
			designfile << " ";
			std::stringstream ss;
                        ss.precision(std::numeric_limits<long double>::digits);//override the default
                        ss << ((MNStruct *)context)->LDpriors[i][0];
                        designfile << ss.str();
			designfile << " ";
			ss << ((MNStruct *)context)->LDpriors[i][1];
			designfile << ss.str();
			designfile << "\n";		
			} 
			*/
			designfile.close();


			std::ofstream planetfile;
			std::string planetname = longname+"planetdesignMatrix.dat";
			planetfile.open(planetname.c_str());
		//	printf("wWriting planets \n");
		/*	for(int i =0; i <((MNStruct *)context)->pulse->nobs; i++){
			
                	for(int j=0; j<9; j++) {
				if(j!=2){
                        		std::stringstream ss;
                        		ss.precision(std::numeric_limits<double>::digits10);//override the default
					printf("%i %i %g \n",i,j,dotproduct(((MNStruct *)context)->pulse->posPulsar,((MNStruct *)context)->pulse->obsn[i].planet_ssb[j]));
                        		ss << dotproduct(((MNStruct *)context)->pulse->posPulsar,((MNStruct *)context)->pulse->obsn[i].planet_ssb[j]);
                        		planetfile << ss.str();
                        		planetfile << " ";
                	}
			}

			planetfile << "\n";
			}
			planetfile.close();
		*/
			//writeCov(paramlist, ndim,context, longname,outputoption);
			
			double Evidence=paramlist[4*ndim];
			//printf("Ev: %i %g \n", ndim, Evidence);
			TNtextOutput(((MNStruct *)context)->pulse, 1, 0, Tempo2Fit,  context,incRED,ndims,paramlist, Evidence, doTimeMargin, doJumpMargin, doLinear, longname, paramarray);

			
			
			
			printf("finished output \n");
		}
		summaryfile.close();

	}
	else{
		printf("More than one mode has been detected, you'll have to handle that on your own for now.");
	}

}


void getmaxlikeDM(pulsar *pulse,std::string longname, int ndim, void *context, double **paramsarray){

	formBatsAll(pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
	formResiduals(pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
//	printf("Jump1: %g \n",pulse[0].jumpVal[((MNStruct *)context)->TempoJumpNums[0]]);


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Form the Design Matrix////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////  

	
	int TimetoFit=((MNStruct *)context)->numFitTiming;
	int JumpstoFit=((MNStruct *)context)->numFitJumps;
	int TimetoMargin=TimetoFit+JumpstoFit;
	int numtofit= TimetoFit+JumpstoFit;
	
	printf("num params, %i %i \n", TimetoFit, JumpstoFit);

	double **TNDM=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		TNDM[i]=new double[TimetoMargin];
	}
	


	//getDMatrix(((MNStruct *)context)->pulse, TimetoFit, JumpstoFit, numtofit, ((MNStruct *)context)->TempoFitNums,((MNStruct *)context)->TempoJumpNums, ((MNStruct *)context)->Dpriors,1, 2, TNDM);
	getCustomDMatrixLike(context, TNDM);
	
	
	std::ofstream designfile;
	std::string dname = longname+"DesignMatrix.txt";

	designfile.open(dname.c_str());

	designfile << ((MNStruct *)context)->pulse->nobs;
	designfile << " ";
	designfile << numtofit;
	designfile << "\n";
	designfile << ((MNStruct *)context)->pulse->param[param_raj].val[0];
	designfile << " ";
	designfile << ((MNStruct *)context)->pulse->param[param_decj].val[0];
	designfile << "\n";		
	for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {

		for(int j=0; j<numtofit; j++) {
			std::stringstream ss;
			ss.precision(std::numeric_limits<double>::digits10);//override the default
			ss << TNDM[i][j];
			designfile << ss.str();
			designfile << " ";
		//	if(j==6)printf("%i %g \n", i, TNDM[i][j]);
		} 
		designfile << "\n";

	} 

	designfile.close();


	double* S = new double[TimetoMargin];
	double** U = new double*[((MNStruct *)context)->pulse->nobs];
	for(int k=0; k < ((MNStruct *)context)->pulse->nobs; k++){
		U[k] = new double[((MNStruct *)context)->pulse->nobs];
	}
	double** VT = new double*[TimetoMargin]; 
	for (int k=0; k<TimetoMargin; k++) VT[k] = new double[TimetoMargin];

	dgesvd(TNDM,((MNStruct *)context)->pulse->nobs, TimetoMargin, S, U, VT);

	double **V=new double*[numtofit];

	for(int i=0;i<numtofit;i++){
		V[i]=new double[numtofit];
	}


	for(int j=0;j < numtofit;j++){
		for(int k=0;k < numtofit;k++){
			V[j][k]=VT[k][j];
		}
	}
		
	

	for(int j=0;j<((MNStruct *)context)->pulse->nobs;j++){
		for(int k=0;k < TimetoMargin;k++){
				TNDM[j][k]=U[j][k];
		}
	}



	///////////////Form the F Matrix////////////////////////////////////////



	int RedAmpParam=((MNStruct *)context)->numFitEFAC + ((MNStruct *)context)->numFitEQUAD;
	int RedSpecParam=RedAmpParam+1;
	int DMAmpParam=RedSpecParam+1;
	int DMSpecParam=DMAmpParam+1;


	double RedAmp=paramsarray[RedAmpParam][2];
	double RedIndex=paramsarray[RedSpecParam][2];
	double DMAmp=paramsarray[DMAmpParam][2];
	double DMIndex=paramsarray[DMSpecParam][2];

	printf("params %i %i %i %i\n", RedAmpParam, RedSpecParam, DMAmpParam, DMSpecParam);
	printf("params %g %g %g %g\n", RedAmp, RedIndex, DMAmp, DMIndex);

	int FitRedCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
	int FitDMCoeff=2*(((MNStruct *)context)->numFitDMCoeff);

    int totCoeff=0;
    if(((MNStruct *)context)->incRED != 0)totCoeff+=FitRedCoeff;
    if(((MNStruct *)context)->incDM != 0)totCoeff+=FitDMCoeff;

	double **FMatrix=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		FMatrix[i]=new double[totCoeff];
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

	double maxtspan=1*(end-start);


	double *freqs = new double[totCoeff];

	double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
	double DMKappa = 2.410*pow(10.0,-16);
	int startpos=0;
	double freqdet=0;

	double *powercoeff=new double[totCoeff];
	for(int o=0;o<totCoeff; o++){
		powercoeff[o]=0;
	}

	double Tspan = maxtspan;
	double f1yr = 1.0/3.16e7;

	RedAmp=pow(10.0, RedAmp);

	for (int i=0; i<FitRedCoeff/2; i++){

		freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
		freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
		
		double rho = (RedAmp*RedAmp/12.0/(M_PI*M_PI))*pow(f1yr,(-3)) * pow(freqs[i]*365.25,(-RedIndex))/(maxtspan*24*60*60);
		powercoeff[i]+= rho;
		powercoeff[i+FitRedCoeff/2]+= rho;
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



	DMAmp=pow(10.0, DMAmp);

	for (int i=0; i<FitDMCoeff/2; i++){

		freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[startpos/2 - ((MNStruct *)context)->incFloatRed +i]/maxtspan;
		freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
		
		double rho = (DMAmp*DMAmp)*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-DMIndex))/(maxtspan*24*60*60);	
		powercoeff[startpos+i]+=rho;
		powercoeff[startpos+i+FitDMCoeff/2]+=rho;
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


	///////////////////////Form Total Matrices////////////////////////////////////////////////

	int totalsize=numtofit+totCoeff;
	double **TotalMatrix=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i =0;i<((MNStruct *)context)->pulse->nobs;i++){
		TotalMatrix[i]=new double[totalsize];
		for(int j =0;j<totalsize; j++){
			TotalMatrix[i][j]=0;
		}
	}
	

	for(int i =0;i<((MNStruct *)context)->pulse->nobs;i++){
		for(int j =0;j<numtofit; j++){
			TotalMatrix[i][j]=TNDM[i][j];
		}
		
		for(int j =0;j<totCoeff; j++){
			TotalMatrix[i][j+numtofit]=FMatrix[i][j];
		}
	}

// 	printf("made TMatrix\n");

	double *Noise=new double[((MNStruct *)context)->pulse->nobs];
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2);
	}

// 	printf("made NMatrix\n");
	
	double **NG = new double*[((MNStruct *)context)->pulse->nobs]; for (int k=0; k<((MNStruct *)context)->pulse->nobs; k++) NG[k] = new double[totalsize];
	double **GNG = new double*[totalsize]; for (int k=0; k<totalsize; k++) GNG[k] = new double[totalsize];



	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		for(int j=0;j<totalsize; j++){

			NG[i][j]=TotalMatrix[i][j]/Noise[i];

		}
	}


	dgemm(TotalMatrix, NG,GNG,((MNStruct *)context)->pulse->nobs, totalsize,((MNStruct *)context)->pulse->nobs, totalsize, 'T','N');


	for(int j =0;j<totCoeff; j++){
		GNG[numtofit+j][numtofit+j]+=PPFM[j][j];
	}


	double tdet=0;
	dpotrf(GNG, totalsize, tdet);
	dpotri(GNG,totalsize);


	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Resvec[o]=(double)pulse->obsn[o].residual;
	}


	double *dG=new double[totalsize];
	dgemv(NG,Resvec,dG,((MNStruct *)context)->pulse->nobs,totalsize,'T');

	double *maxcoeff=new double[totalsize];
	dgemv(GNG,dG,maxcoeff,totalsize,totalsize,'N');

	double *Errorvec=new double[totalsize];

	for(int i =0; i < totalsize; i++){
		Errorvec[i]=pow(GNG[i][i], 0.5);
	}

        double *Scoeff=new double[numtofit];
        double *Serr=new double[numtofit];
        for(int i =0; i < numtofit; i++){
                Scoeff[i]=maxcoeff[i]/S[i];
                Serr[i]=Errorvec[i]/S[i];
        }

        double *TempoCoeff = new double[numtofit];
        double *TempoErr =  new double[numtofit];
        dgemv(V,Scoeff,TempoCoeff,numtofit,numtofit, 'N');
        dgemv(V,Serr,TempoErr,numtofit,numtofit, 'N');

        for(int i=0;i<numtofit; i++){
                double errsum=0;
          for(int j=0;j<numtofit; j++){
                                errsum += pow(V[i][j]*Serr[j],2);
                }
          TempoErr[i]=pow(errsum,0.5);
        }
        //updateParameters(((MNStruct *)context)->pulse,0,TempoCoeff,TempoErr);
	 for(int p=0;p<((MNStruct *)context)->numFitTiming;p++){
                        //printf("%i %.25Lg \n", p, ((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]]);
                }
                for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
                        //printf("%i %g \n", p,((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]);
                }




 	char name[1000];

        std::ofstream Finalfile;
        std::string Finalfilename = longname+"Final.dat";
        Finalfile.open(Finalfilename.c_str());


        for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
                double dsum=0;
                for(int j=0;j<totalsize; j++){
                        dsum=dsum+TotalMatrix[i][j]*maxcoeff[j];
                }
		Finalfile << std::fixed  << std::setprecision(8)  << (double)((MNStruct *)context)->pulse->obsn[i].bat << " ";
        Finalfile << std::scientific << Resvec[i] << " " << dsum <<  " " << pow(Noise[i],0.5) << " ";
	Finalfile << std::fixed  << std::setprecision(1)  <<(double)((MNStruct *)context)->pulse->obsn[i].freqSSB << " ";
		Finalfile << ((MNStruct *)context)->sysFlags[i]  << "\n";


        }
	Finalfile.close();

        std::ofstream Resfile;
        std::string Resfilename = longname+"Res.dat";
        Resfile.open(Resfilename.c_str());
        Resfile << ((MNStruct *)context)->pulse->name << "\n";
        Resfile << "obsid    Freq (MHz)   SAT    Res (s)    Error (s)\n";

        for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
                double dsum=0;
                for(int j=0;j<numtofit; j++){
                        dsum=dsum+TotalMatrix[i][j]*maxcoeff[j];
                }
                sscanf(((MNStruct *)context)->pulse->obsn[i].fname,"%s",name);
                Resfile <<  name << " ";
                Resfile << std::fixed  << std::setprecision(8)  << (double)((MNStruct *)context)->pulse->obsn[i].freq << " " << (double)((MNStruct *)context)->pulse->obsn[i].sat << " ";
                Resfile << std::scientific << Resvec[i]-dsum <<  " " << pow(Noise[i],0.5) << "\n";
        }

        Resfile.close();




        std::ofstream DMfile;
        std::string DMfilename = longname+"DM.dat";
        DMfile.open(DMfilename.c_str());
	DMfile << ((MNStruct *)context)->pulse->name << "\n";
	DMfile << "obsid    Freq (MHz)   SAT    DMRes (s)    Error (s)  DMVal    Error\n";


	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		double dsum=0;
		double scaledsum=0;
		double errsum=0;
		double scaleerrsum=0;
		for(int j=0;j<FitDMCoeff; j++){
			dsum=dsum+TotalMatrix[i][j+numtofit+FitRedCoeff]*maxcoeff[j+numtofit+FitRedCoeff];
			scaledsum=scaledsum+TotalMatrix[i][j+numtofit+FitRedCoeff]*maxcoeff[j+numtofit+FitRedCoeff]/DMVec[i];
			errsum=errsum+pow(Errorvec[j+numtofit+FitRedCoeff]*TotalMatrix[i][j+numtofit+FitRedCoeff],2);
			scaleerrsum=scaleerrsum+pow(Errorvec[j+numtofit+FitRedCoeff]*TotalMatrix[i][j+numtofit+FitRedCoeff]/DMVec[i],2);
		}
		
		double freq=(double)((MNStruct *)context)->pulse->obsn[i].freqSSB;
	        long double yrs = (((MNStruct *)context)->pulse->obsn[i].sat - ((MNStruct *)context)->pulse->param[param_dmepoch].val[0])/365.25;
 	        long double arg = 1.0;
		double dmDot=0;
		double dmDotErr=0;
                for (int d=1;d<9;d++){	
          		arg *= yrs;
          		if (((MNStruct *)context)->pulse->param[param_dm].paramSet[d]==1){
            			dmDot+=(double)(((MNStruct *)context)->pulse->param[param_dm].val[d]*arg);
				dmDotErr+=pow((double)(((MNStruct *)context)->pulse->param[param_dm].err[d]*arg),2);
			}
        	}
		dsum=dsum+dmDot*DMVec[i];
		scaledsum=scaledsum+dmDot;
		errsum=errsum+dmDotErr*pow(DMVec[i],2);
		scaleerrsum=scaleerrsum+dmDotErr;
		double restsum=0;
		for(int j=0;j<numtofit+FitRedCoeff; j++){
			restsum=restsum+TotalMatrix[i][j]*maxcoeff[j];
		}
		sscanf(((MNStruct *)context)->pulse->obsn[i].fname,"%s",name);
		DMfile <<  name << " ";
		DMfile << std::fixed  << std::setprecision(8)  << (double)((MNStruct *)context)->pulse->obsn[i].freq << " " << (double)((MNStruct *)context)->pulse->obsn[i].sat << " "; 
		DMfile << std::scientific << dsum <<  " " << pow(errsum,0.5) << " " << scaledsum << " " << pow(scaleerrsum,0.5) << "\n";
	}

	DMfile.close();



        std::ofstream Redfile;
        std::string Redfilename = longname+"Red.dat";
        Redfile.open(Redfilename.c_str());
        Redfile << ((MNStruct *)context)->pulse->name << "\n";
        Redfile << "obsid    Freq (MHz)   SAT    RedRes (s)    Error (s) \n";

        for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
                double dsum=0;
                double errsum=0;
                for(int j=0;j<FitRedCoeff; j++){
                        dsum=dsum+TotalMatrix[i][j+numtofit]*maxcoeff[j+numtofit];
                        errsum=errsum+pow(Errorvec[j+numtofit]*TotalMatrix[i][j+numtofit],2);
                }
                sscanf(((MNStruct *)context)->pulse->obsn[i].fname,"%s",name);
                Redfile <<  name << " ";
                Redfile << std::fixed  << std::setprecision(8)  << (double)((MNStruct *)context)->pulse->obsn[i].freq << " " << (double)((MNStruct *)context)->pulse->obsn[i].sat << " ";
                Redfile << std::scientific << dsum <<  " " << pow(errsum,0.5) << "\n";
        }

        Redfile.close();


	delete[]S;	

	for (int j = 0; j < TimetoMargin; j++){
		delete[]VT[j];
	}
	
	delete[]VT;

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[]U[j];
	}
	delete[]U;

}

void getCustomDMatrix(pulsar *pulse, int *MarginList, double **TNDM, int **TempoFitNums, int *TempoJumpNums, double **Dpriors, int incDM, int TimetoFit, int JumptoFit){
	
	double pdParamDeriv[MAX_PARAMS], dMultiplication;


	//Unset all fit flags for parameters we arn't marginalising over so they arn't in the design Matrix


		int pcount=1;
		int numToMargin=1;
		
		Dpriors[0][0]=0;
		Dpriors[0][1]=0;
		
		for (int p=1;p<TimetoFit;p++) {
			if(MarginList[pcount]!=1){
				pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 0;
			}
			else if(MarginList[pcount]==1){
				pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 1;
				Dpriors[pcount][0]=0;
				Dpriors[pcount][1]=0;
				numToMargin++;
			}
			pcount++;
		}
	
		for(int i=0; i < JumptoFit; i++){
			if(MarginList[pcount]!=1){
				pulse[0].fitJump[TempoJumpNums[i]]=0;
			}
			else if(MarginList[pcount]==1){
				pulse[0].fitJump[TempoJumpNums[i]]=1;
				Dpriors[pcount][0]=0;
				Dpriors[pcount][1]=0;
				numToMargin++;
			}
			pcount++;
		}
	
	
		for(int i=0; i < pulse->nobs; i++) {
			FITfuncs(pulse[0].obsn[i].bat - pulse[0].param[param_pepoch].val[0], pdParamDeriv, numToMargin, pulse, i,0);
			for(int j=0; j<numToMargin; j++) {
				TNDM[i][j]=pdParamDeriv[j];
					//printf("Dmatrix: %i %i %22.20e \n", i,j,pdParamDeriv[j]);
			} 
		} 

		//Now set fit flags back to how they were
	
		for (int p=1;p<TimetoFit;p++) {
			pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 1;
		}
	
		for(int i=0; i < JumptoFit; i++){
			pulse[0].fitJump[TempoJumpNums[i]]=1;
		}
}


void getCustomDMatrixLike(void *context, double **TNDM){
	
	double pdParamDeriv[MAX_PARAMS], dMultiplication;


	//Unset all fit flags for parameters we arn't marginalising over so they arn't in the design Matrix


		int pcount=1;
		int numToMargin=1;
		
		for (int p=1;p<((MNStruct *)context)->numFitTiming;p++) {
			if(((MNStruct *)context)->LDpriors[pcount][2]!=1){
				((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].fitFlag[((MNStruct *)context)->TempoFitNums[p][1]] = 0;
			}
			else if(((MNStruct *)context)->LDpriors[pcount][2]==1){
				((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].fitFlag[((MNStruct *)context)->TempoFitNums[p][1]] = 1;
				numToMargin++;
			}
			pcount++;
		}
	
		for(int i=0; i < ((MNStruct *)context)->numFitJumps; i++){
			if(((MNStruct *)context)->LDpriors[pcount][2]!=1){
				((MNStruct *)context)->pulse->fitJump[((MNStruct *)context)->TempoJumpNums[i]]=0;
			}
			else if(((MNStruct *)context)->LDpriors[pcount][2]==1){
				((MNStruct *)context)->pulse->fitJump[((MNStruct *)context)->TempoJumpNums[i]]=1;
				numToMargin++;
			}
			pcount++;
		}
	
	
		for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {
			FITfuncs(((MNStruct *)context)->pulse->obsn[i].bat - ((MNStruct *)context)->pulse->param[param_pepoch].val[0], pdParamDeriv, numToMargin, ((MNStruct *)context)->pulse, i,0);
			for(int j=0; j<numToMargin; j++) {
//				printf("CDML: %i %i %i %g\n", i,j,numToMargin,pdParamDeriv[j]);
				TNDM[i][j]=pdParamDeriv[j];
			} 
		} 

		//Now set fit flags back to how they were
	

		for (int p=1;p<((MNStruct *)context)->numFitTiming;p++) {
				((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].fitFlag[((MNStruct *)context)->TempoFitNums[p][1]] = 1;

		}
	
		for(int i=0; i < ((MNStruct *)context)->numFitJumps; i++){
			((MNStruct *)context)->pulse->fitJump[((MNStruct *)context)->TempoJumpNums[i]]=1;
		}
	
}




void getDMatrix(pulsar *pulse, int TimetoFit, int JumptoFit, int numToMargin, int **TempoFitNums, int *TempoJumpNums, double **Dpriors, int doJumpMargin, int doTimeMargin, double **TNDM){
	
	double pdParamDeriv[MAX_PARAMS], dMultiplication;


	//Unset all fit flags for parameters we arn't marginalising over so they arn't in the design Matrix
	//If we are marginalising, change the prior so you dont sample it
	if(JumptoFit>0){
		if(doJumpMargin==0){
			for(int i=0; i < JumptoFit; i++){
				pulse[0].fitJump[TempoJumpNums[i]]=0;
			}
		} 
		else if(doJumpMargin==1){
			for(int i=0; i < JumptoFit; i++){
// 				Dpriors[TimetoFit+i][0]=0;
// 				Dpriors[TimetoFit+i][1]=0;
				
			}
		}
	}


	if(doTimeMargin==0){
		for (int p=1;p<TimetoFit;p++) {
			pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 0;
		}
	}
	//If only marginalising over QSD (i.e. F0 and F1) then set all other flags to 0
	else if(doTimeMargin==1){	
		int pcount=0;
		for (int p=0;p<MAX_PARAMS;p++) {
			if(strcasecmp(pulse[0].param[p].shortlabel[0],"F0")!=0){
				for (int k=0;k<pulse[0].param[p].aSize;k++){
					if(pulse[0].param[p].fitFlag[k] == 1){
						//printf("in getD %s \n",pulse[0].param[p].shortlabel[k]);
						pulse[0].param[p].fitFlag[k]=0;
						pcount++;
					}
				}
			}
			else if(strcasecmp(pulse[0].param[p].shortlabel[0],"F0")==0){
				for (int k=0;k<pulse[0].param[p].aSize;k++){
					if(pulse[0].param[p].fitFlag[k] == 1){
						//printf("in getD %s setting prior %i to zero\n",pulse[0].param[p].shortlabel[k], pcount);
// 						Dpriors[pcount][0]=0;
// 						Dpriors[pcount][1]=0;
						pcount++;
					}
				}
			}
					
		}
	}
	else if(doTimeMargin==2){	
		for(int i=0; i < TimetoFit; i++){
// 			Dpriors[i][0]=0;
// 			Dpriors[i][1]=0;
		}
	}


	for(int i=0; i < pulse->nobs; i++) {
		FITfuncs(pulse[0].obsn[i].bat - pulse[0].param[param_pepoch].val[0], pdParamDeriv, numToMargin-1, pulse, i,0);
		for(int j=0; j<numToMargin; j++) {
			TNDM[i][j]=pdParamDeriv[j];
 			printf("%i %i %22.20e \n", i,j,pdParamDeriv[j]);

		} 

	} 


	//Now set fit flags back to how they were

	for (int p=1;p<TimetoFit;p++) {
		pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 1;
	}

	for(int i=0; i < JumptoFit; i++){
		pulse[0].fitJump[TempoJumpNums[i]]=1;
	}

}


void getMarginDMatrix(pulsar *pulse, int TimetoFit, int JumptoFit, int numToMargin, int **TempoFitNums, int *TempoJumpNums, double **Dpriors, int doJumpMargin, int doTimeMargin, double **TNDM, int linearFit){
	
	double pdParamDeriv[MAX_PARAMS], dMultiplication;


	//Unset all fit flags for parameters we arn't marginalising over so they arn't in the design Matrix
	//If we are marginalising, change the prior so you dont sample it


		if(JumptoFit>0){
		
				Dpriors[0][0]=0;  //Prior on Phase is zero
				Dpriors[0][1]=0;
				
				
			if(doJumpMargin==0){
				for(int i=0; i < JumptoFit; i++){
					pulse[0].fitJump[TempoJumpNums[i]]=0;
				}
			} 
			else if(doJumpMargin==1){
				for(int i=0; i < JumptoFit; i++){
					Dpriors[TimetoFit+i][0]=0;
					Dpriors[TimetoFit+i][1]=0;
					
				}
			}
		}
	
	
		if(doTimeMargin==0){
			for (int p=1;p<TimetoFit;p++) {
				pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 0;
			}
		}
		//If only marginalising over QSD (i.e. F0 and F1) then set all other flags to 0
		else if(doTimeMargin==1){	
			int pcount=0;
			Dpriors[pcount][0]=0;
			Dpriors[pcount][1]=0;
			pcount++;
			for (int p=0;p<MAX_PARAMS;p++) {
				if(strcasecmp(pulse[0].param[p].shortlabel[0],"F0")!=0){
					for (int k=0;k<pulse[0].param[p].aSize;k++){
						if(pulse[0].param[p].fitFlag[k] == 1){
							//printf("in getDNM %s \n",pulse[0].param[p].shortlabel[k]);
							pulse[0].param[p].fitFlag[k]=0;
							pcount++;
						}
					}
				}
				else if(strcasecmp(pulse[0].param[p].shortlabel[0],"F0")==0){
					for (int k=0;k<pulse[0].param[p].aSize;k++){
						if(pulse[0].param[p].fitFlag[k] == 1){
							//printf("in getDNM %s setting prior %i to zero\n",pulse[0].param[p].shortlabel[k], pcount);
							Dpriors[pcount][0]=0;
							Dpriors[pcount][1]=0;
							pcount++;
						}
					}
				}
						
			}
		}
		else if(doTimeMargin==2){	
			for(int i=0; i < TimetoFit; i++){
				Dpriors[i][0]=0;
				Dpriors[i][1]=0;
			}
		}
	
	
		for(int i=0; i < pulse->nobs; i++) {
			FITfuncs(pulse[0].obsn[i].bat - pulse[0].param[param_pepoch].val[0], pdParamDeriv, numToMargin, pulse, i,0);
			for(int j=0; j<numToMargin; j++) {
				TNDM[i][j]=pdParamDeriv[j];
 			//	printf("%i %i %22.20e \n", i,j,pdParamDeriv[j]);
	
			} 
	
		} 
	
	
		//Now set fit flags back to how they were
	
		for (int p=1;p<TimetoFit;p++) {
			pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 1;
		}
	
		for(int i=0; i < JumptoFit; i++){
			pulse[0].fitJump[TempoJumpNums[i]]=1;
		}
	




}



void makeGDesign(pulsar *pulse, int &Gsize, int numtofit, double** staticGMatrix, double **oneDesign){

	
	int numpulsars=1;
	double*** Umatrices=new double**[numpulsars];
	for(int i=0;i<numpulsars;i++){
		Umatrices[i]=new double*[pulse->nobs];
		for(int j=0;j<pulse->nobs;j++){
			Umatrices[i][j]=new double[pulse->nobs];
		}
	}

	for(int i=0;i<numpulsars;i++){

		double* S = new double[numtofit];
		double** U = new double*[pulse->nobs]; for (int k=0; k<pulse->nobs; k++) U[k] = new double[pulse->nobs];
		double** VT = new double*[numtofit]; for (int k=0; k<numtofit; k++) VT[k] = new double[numtofit];

		dgesvd(oneDesign,pulse->nobs, numtofit, S, U, VT);

		for(int j=0;j<pulse->nobs;j++){
// 			printf("h %i %i\n",i,j);
			for(int k=0;k<pulse->nobs;k++){
 				//printf("U %i %i %i %i %g \n",pulse->nobs, numtofit, j,k,U[j][k]);
				Umatrices[i][j][k]=U[j][k];
				
			}
		}

		delete[]S;	
		for (int j = 0; j < pulse->nobs; j++){
			delete[]U[j];
		}
		for (int j = 0; j < numtofit; j++){
			delete[]VT[j];
		}
		delete[]U;
		delete[]VT;
	
		
	}

// 	printf("Done SVD ALL");
	int Dsize=0;
	for(int i=0;i<numpulsars;i++){
		Dsize += numtofit;
	}




	int Osum=0;
	int Gsum=0;
	for(int i=0;i<numpulsars;i++){
		for(int j=0;j<pulse->nobs;j++){
			int nfsum=0;
			for(int k=0;k<pulse->nobs;k++){
				if(k>=numtofit){
					staticGMatrix[Osum+j][Osum-Gsum+nfsum]=Umatrices[i][j][k];
 					//printf("GBack %i %i %g \n",Osum+j,Osum-Gsum+nfsum,staticGMatrix[Osum+j][Osum-Gsum+nfsum]);
					nfsum++;

				}
			}
		}
		Gsum+=numtofit;
		Osum+=pulse->nobs;
	}


	Gsize=Osum-Gsum;
/*	for(int i=0;i<pulse->nobs;i++){
		for(int j=0;j<Gsize;j++){
			printf("%i %i %g \n ",i,j,staticGMatrix[i][j]);
		}
	}
	printf("%i %i \n",pulse->nobs,Gsize)*/;

// 	printf("Done SVD ALL");
}

void makeStaticGMatrix(pulsar *pulse, int Gsize, double **GMatrix, double** staticGMatrix, double &tdet){


	double *Noise=new double[pulse->nobs];
	for(int o=0;o < pulse->nobs; o++){
		Noise[o]=pow(pulse->obsn[o].toaErr*pow(10.0,-6),2);
	}

	double **NG = new double*[pulse->nobs]; for (int k=0; k< pulse->nobs; k++) NG[k] = new double[Gsize];
	for(int i=0;i< pulse->nobs;i++){
		for(int j=0;j< Gsize; j++){
			NG[i][j] = GMatrix[i][j]*Noise[i];
		}
	}

	double **GG = new double*[Gsize]; for (int k=0; k< Gsize; k++) GG[k] = new double[Gsize];

	dgemm(GMatrix, NG,GG,pulse->nobs, Gsize,pulse->nobs, Gsize, 'T','N');
	
	tdet=0;
	dpotrf(GG, Gsize, tdet);
	dpotri(GG,Gsize);
	


	dgemm(GMatrix, GG,NG, pulse->nobs, Gsize, Gsize, Gsize, 'N','N');

	dgemm(NG, GMatrix, staticGMatrix, pulse->nobs, Gsize, pulse->nobs, Gsize, 'N','T');
	
	delete[] Noise;
	
	for (int j = 0; j < pulse->nobs; j++){
		delete[] NG[j];
	}
	delete[] NG;

	for (int j = 0; j < Gsize; j++){
		delete[]GG[j];
	}
	delete[] GG;	
	
	
}

void makeStaticDiagGMatrix(pulsar *pulse, int Gsize, double **GMatrix, double** GNMatrix, double *SVec){


		




	double *Noise=new double[pulse->nobs];
	for(int o=0;o < pulse->nobs; o++){
		Noise[o]=pow(pulse->obsn[o].toaErr*pow(10.0,-6),2);
	}

	double **NG = new double*[pulse->nobs]; for (int k=0; k< pulse->nobs; k++) NG[k] = new double[Gsize];
	for(int i=0;i< pulse->nobs;i++){
		for(int j=0;j< Gsize; j++){
			NG[i][j] = GMatrix[i][j]*Noise[i];

		}
	}

	double **GG = new double*[Gsize]; for (int k=0; k< Gsize; k++) GG[k] = new double[Gsize];

	dgemm(GMatrix, NG,GG,pulse->nobs, Gsize,pulse->nobs, Gsize, 'T','N');
	

	double** U = new double*[Gsize]; for (int k=0; k<Gsize; k++) U[k] = new double[Gsize];
	double** VT = new double*[Gsize]; for (int k=0; k<Gsize; k++) VT[k] = new double[Gsize];

	dgesvd(GG,Gsize, Gsize, SVec, U, VT);
	
	double **GT = new double*[Gsize]; for (int k=0; k< Gsize; k++) GT[k] = new double[pulse->nobs];
	
	for(int i=0;i< Gsize;i++){
			for(int j=0;j< pulse->nobs; j++){
				GT[i][j]=GMatrix[j][i];
				
			}
		}
		
		
		
	dgemm(VT, GT,GNMatrix,Gsize,Gsize, Gsize, pulse->nobs, 'N','N');
	

	delete[] Noise;
	
	for (int j = 0; j < pulse->nobs; j++){
		delete[] NG[j];
		
	}
	delete[] NG;

	for (int j = 0; j < Gsize; j++){
		delete[]GG[j];
		delete[]VT[j];
		delete[] U[j];
		delete[] GT[j];

	}
	delete[] GG;	
	delete[]VT;
	delete[] U;
	delete[] GT;
	
}


void convertFromLinear(pulsar *psr, std::string longname, int ndim, void *context){

	int number_of_lines = 0;
	char *outname;

	std::ifstream checkfile;
	std::string checkname = longname+".txt";
	checkfile.open(checkname.c_str());
	std::string line;
	while (getline(checkfile, line))
		++number_of_lines;

	printf("CFL number of lines %i \n",number_of_lines);
	checkfile.close();
	
	std::ifstream summaryfile;
	std::string fname = longname+".txt";
	summaryfile.open(fname.c_str());

	std::ofstream newtxtfile;
	std::string newtxtname = longname+"NL.txt";
	newtxtfile.open(newtxtname.c_str());

	double *linearParams = new double[((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps+1];
	double *nonLinearParams = new double[((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps+1];
	double *errorvec = new double[((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps+1];
	for(int i=0;i<number_of_lines;i++){
		
		std::string line;
		getline(summaryfile,line);
		std::istringstream myStream( line );
		std::istream_iterator< double > begin(myStream),eof;
		std::vector<double> paramlist(begin,eof);

		
		for(int j=0;j<((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps;j++){
			linearParams[j]=paramlist[2+j];
			nonLinearParams[j]=0;
		}
		TNupdateParameters(psr,0,linearParams,errorvec, nonLinearParams);
		
		std::stringstream txtstream;
		txtstream.precision(std::numeric_limits<double>::digits10);//override the default


		txtstream << paramlist[0];
		txtstream << "\t";
		txtstream << paramlist[1];
		txtstream << "\t";
		for(int j=0;j<((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps;j++){
			if(j==0){
				txtstream << nonLinearParams[j];
				txtstream << "\t";
			}
			else{
				if(((MNStruct *)context)->LDpriors[j][1] != 0){
					txtstream << nonLinearParams[j]/((MNStruct *)context)->LDpriors[j][1];
					txtstream << "\t";
				}
				else{
					txtstream << 0;
                    txtstream << "\t";
				}
// 				printf("%i %i %g %g %g \n",i,j,linearParams[j],nonLinearParams[j],(double)((MNStruct *)context)->LDpriors[j-1][1]);
			}
		}
		for(int j=((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps;j<ndim;j++){

				txtstream << paramlist[j+2];
				txtstream << "\t";


		}
		newtxtfile << txtstream.str();
		newtxtfile << "\n";

	}
	summaryfile.close();
	newtxtfile.close();



}


void writeCov(std::vector<double> Cube, int &ndim,void *context, std::string longname, int outputoption)
{

	clock_t startClock,endClock;
	long double LDparams[ndim];
	double *EFAC;
	double EQUAD, redamp, redalpha;
	
			  for(int i =0; i< 4*ndim; i++){
		  	printf("bparm: %i %g \n ", i, Cube[i]);
			}	
	
	
	int pcount=outputoption*ndim;


	for(int p=0;p<((MNStruct *)context)->numFitTiming;p++){
		pcount++;
	}
	for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
		pcount++;
	}

	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double[1];
		EFAC[0]=1;
// 		
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double[1];
		EFAC[0]=Cube[pcount];
		pcount++;
		
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)context)->numFitEFAC];
		for(int p=0;p< ((MNStruct *)context)->numFitEFAC; p++){
			EFAC[p]=Cube[pcount];
			printf("EF: %i %g \n",p,EFAC[p]);
			pcount++;
			
		}
	}				

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=0;
// 		printf("EQUAD: %g \n",EQUAD);
	}
	else{
		
		EQUAD=pow(10.0,2*Cube[pcount]);
		printf("EQ: %g \n",Cube[pcount]);
		pcount++;
		

	}


	redamp=Cube[pcount];
	pcount++;
	redalpha=Cube[pcount];
	pcount++;

	printf("red: %g %g \n",redamp, redalpha);
	double secday=24*60*60;
	double LongestPeriod=1.0/pow(10.0,-5);
	double flo=1.0/LongestPeriod;

	double modelalpha=redalpha;
	double gwamp=pow(10.0,redamp);
	double gwampsquared=gwamp*gwamp*(pow((365.25*secday),2)/(12*M_PI*M_PI))*(pow(365.25,(1-modelalpha)))/(pow(flo,(modelalpha-1)));

	double timdiff=0;

	double covconst=gsl_sf_gamma(1-modelalpha)*sin(0.5*M_PI*modelalpha);

	
	double **CovMatrix = new double*[((MNStruct *)context)->pulse->nobs]; for(int o1=0;o1<((MNStruct *)context)->pulse->nobs;o1++)CovMatrix[o1]=new double[((MNStruct *)context)->pulse->nobs];

	for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){

		for(int o2=0;o2<((MNStruct *)context)->pulse->nobs; o2++){
			timdiff=((MNStruct *)context)->pulse->obsn[o1].bat-((MNStruct *)context)->pulse->obsn[o2].bat;	
			double tau=2.0*M_PI*fabs(timdiff);
			double covsum=0;

			for(int k=0; k <=4; k++){
				covsum=covsum+pow(-1.0,k)*(pow(flo*tau,2*k))/(iter_factorial(2*k)*(2*k+1-modelalpha));

			}

			CovMatrix[o1][o2]=gwampsquared*(covconst*pow((flo*tau),(modelalpha-1)) - covsum);
// 			printf("%i %i %g %g %g\n",o1,o2,CovMatrix[o1][o2],fabs(timdiff),covsum);

			if(o1==o2){
				CovMatrix[o1][o2] += pow(((((MNStruct *)context)->pulse->obsn[o1].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o1]],2) + EQUAD;
			}

		}
	}
	
	
	FILE *pFile;
	std::string covname = longname+"covMatrix.txt";
	pFile = fopen(covname.c_str(), "w+");
	
	for(int i=0; i<((MNStruct *)context)->pulse->nobs; i++) {

    	for(int j=0; j<((MNStruct *)context)->pulse->nobs; j++) {
     		 fprintf(pFile, "%22.20e", CovMatrix[i][j]);
     		 fprintf(pFile, " ");
     	}
        fprintf(pFile, "\n");
     }
     	
     fclose(pFile);	
     

	std::string maxstoc = longname+"maxStoc.txt";
	pFile = fopen(maxstoc.c_str(), "w+");
	
	for(int p=0;p< ((MNStruct *)context)->numFitEFAC; p++){
			fprintf(pFile, "Max EFAC %i %22.20e \n",p, EFAC[p]);
		}

	fprintf(pFile, "Max EQUAD 0 %22.20e \n",0.50*log10(EQUAD));
	fprintf(pFile, "Max Red Amp %22.20e \n",redamp);
	fprintf(pFile, "Max Red Index %22.20e \n",redalpha);

     fclose(pFile);	
}	 
     		 
