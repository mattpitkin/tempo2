#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <T2toolkit.h>
#include <fftw3.h>
#include "tempo2.h"
#include "makeRedNoise.h"

float CubicInterpolate( float y0,float y1, float y2,float y3, float mu);
float CatmullRomInterpolate( float y0,float y1, float y2,float y3, float mu);



rednoisemodel_t* setupRedNoiseModel(float start,float end, int npt, int nreal, float pwr_1yr, float index){
	if (nreal < 100)nreal=100;

	rednoisemodel_t* model = (rednoisemodel_t*) malloc(sizeof(rednoisemodel_t));
	model->start=start;
	model->end=end;
	model->npt=npt;
	model->nreal=nreal;
	model->pwr_1yr=pwr_1yr;
	model->index=index;
	model->flatten=0;
	model->cutoff=0;
	model->data=NULL;
	model->tres=0;
	model->mode=MODE_SIMPLE;
	return model;
}

void freeRedNoiseModel(rednoisemodel_t* model){
	if(model->data != NULL)free(model->data);
	free(model);
}

void populateRedNoiseModel(rednoisemodel_t* model,long seed){
	int t_npts;
	int i;
	double freq,A,index;
	double t_samp,f_bin,t_span;
	float *data;
	fftwf_plan plan;
	fftwf_complex *spectrum;
	double secperyear=365*86400.0;

	t_samp=(model->end - model->start)/(float)model->npt;

	model->start-=t_samp;
	model->end+=t_samp*2.0;

	if (seed == 0){
		seed=TKsetSeed();
	}else if (seed > 0){
		seed=-seed;
	}

	t_npts=model->npt*model->nreal;

	spectrum = (fftwf_complex*) fftwf_malloc((t_npts/2+1)*sizeof(fftwf_complex));
	data = (float*) fftwf_malloc(t_npts*sizeof(float));
	
	t_span=(model->end - model->start)/365.25; // years

	t_samp=t_span/(float)model->npt;
	model->tres=t_samp*365.25;
	f_bin=1.0/(t_span*model->nreal); // frequency in yr^-1


    // To convert PSD to FFT power, need to multiply by N^2/T.
	// BUT FFTW does not divide by N^2 on inverse transform
	// so we need to divide by N^2.
	// The factor of 4, converts one-sided P to 2-sided P as
	// we input one-sided, and FFTW expected 2-sided.
	//
	// I am now 99% sure this is correct. George looked at it too!
	// M. Keith 2013.
	A=model->pwr_1yr /(4*t_span*model->nreal);


	// we are forming "amplitudes" not powers
	// so square root.
	index=model->index/2.0;
	A = sqrt(A);

	// form a complex spectrum, then do a c2r transform.
	spectrum[0]=0;
	for (i=1; i < t_npts/2+1; i++){
	   freq=(double)i*f_bin;
	   double scale=0;
	   if(model->mode==MODE_T2CHOL){
		  // same definition as the Cholesky code (except index is negative)
		  scale=A*pow(1.0+pow(fabs(freq)/model->flatten,2),index/2.0);
	   }
	   if(model->mode==MODE_SIMPLE){
		  if (freq < model->flatten)
			 freq=model->flatten;

		  scale = A*pow(freq,index);
		  if (freq < model->cutoff)scale=0;
	   }
	   // complex spectrum
	   spectrum[i]=(scale*TKgaussDev(&seed) + I*scale*TKgaussDev(&seed));
	}
	
	plan=fftwf_plan_dft_c2r_1d(t_npts,spectrum,data,FFTW_ESTIMATE);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	fftwf_free(spectrum);

	model->data=data;
}



float getRedNoiseValue(rednoisemodel_t* model, float mjd,int real){
   int i = (int)((mjd-model->start)/model->tres);
   float mu = ((mjd-model->start)/model->tres) - (float)i;

   i+=real*model->npt;
   return CatmullRomInterpolate(model->data[i-1],model->data[i],model->data[i+1],model->data[i+2],mu);
}



float CubicInterpolate( float y0,float y1, float y2,float y3, float mu)
{
   float a0,a1,a2,a3,mu2;

   mu2 = mu*mu;
   a0 = y3 - y2 - y0 + y1;
   a1 = y0 - y1 - a0;
   a2 = y2 - y0;
   a3 = y1;

   return(a0*mu*mu2+a1*mu2+a2*mu+a3);
}


float CatmullRomInterpolate( float y0,float y1, float y2,float y3, float mu)
{
   float a0,a1,a2,a3,mu2;

   mu2 = mu*mu;
   a0 = -0.5*y0 + 1.5*y1 - 1.5*y2 + 0.5*y3;
   a1 = y0 - 2.5*y1 + 2*y2 - 0.5*y3;
   a2 = -0.5*y0 + 0.5*y2;
   a3 = y1;
   return(a0*mu*mu2+a1*mu2+a2*mu+a3);
}

float* getPowerSpectrum(rednoisemodel_t* model){
   fftwf_complex *spectrum;
   fftwf_plan plan;
   float *data;
   float *power_spectrum;
   spectrum = (fftwf_complex*) fftwf_malloc((model->npt/2+1)*sizeof(fftwf_complex));
   data = (float*) fftwf_malloc((model->npt)*sizeof(fftwf_complex));
   power_spectrum = (float*) malloc((model->npt/2+1)*sizeof(float));

   float t_span=(model->end - model->start)/365.25; // years
   float f_bin=1.0/(t_span); // frequency in yr^-1
   for (int p =0 ; p < (model->npt/2+1); p++){
	  power_spectrum[p]=0;
   }
   plan=fftwf_plan_dft_r2c_1d(model->npt,data,spectrum,FFTW_ESTIMATE);
   for (int r =0 ; r < model->nreal; r++){
	  int off=r*model->npt;
	  for (int p =0 ; p < model->npt; p++){
		 data[p]=model->data[p+off];
	  }

	  float dx=t_span/(float)(model->npt); // yr
	  int p;
	  for (p =0 ; p < model->npt-1; p++){
		 float dy=data[p+1]-data[p];
		 data[p] = dy/dx;
	  }
	  float dy=data[0]-data[p];
	  data[p] = 0;
	  fftwf_execute(plan);
	  for (int p =1 ; p < (model->npt/2+1); p++){
		 float factor=1.0/((float)p*f_bin);
		 spectrum[p]/=(float)(model->npt/2);
		 spectrum[p]*=factor;
		 power_spectrum[p]+=crealf(spectrum[p]*conjf(spectrum[p]));
	  }
   }
   fftwf_destroy_plan(plan);
   fftwf_free(spectrum);
   fftwf_free(data);
   for (int p =0 ; p < (model->npt/2+1); p++){
	  power_spectrum[p]/=(float)(model->nreal);
   }

   return power_spectrum;
}

