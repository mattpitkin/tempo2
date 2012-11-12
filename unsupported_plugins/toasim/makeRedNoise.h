#define MODE_T2CHOL 1
#define MODE_SIMPLE 0


typedef struct rednoisemodel {
	float start;   // start MJD
	float end;     // end MJD
	int npt;       // points per realisation
	int nreal;     // number of realisations
	float pwr_1yr; // power at 1 year
	float index;   // index in power spectrum
	float cutoff;  // model is zero below 'cutoff' (yr^-1)
	float flatten; // model is flat below 'flatten'(yr^-1)
	float* data;   // data
	float tres;
	char mode;
} rednoisemodel_t;


rednoisemodel_t* setupRedNoiseModel(float start,float end, int npt, int nreal, float amp_1yr, float index);
void populateRedNoiseModel(rednoisemodel_t* model,long seed);
float getRedNoiseValue(rednoisemodel_t* model, float mjd,int real);
void freeRedNoiseModel(rednoisemodel_t* model);
float* getPowerSpectrum(rednoisemodel_t* model);
