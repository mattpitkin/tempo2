#define MODE_T2CHOL 1
#define MODE_SIMPLE 0


typedef struct rednoisemodel {
	double start;   // start MJD
	double end;     // end MJD
	int npt;       // points per realisation
	int nreal;     // number of realisations
	double pwr_1yr; // power at 1 year
	double index;   // index in power spectrum
	double cutoff;  // model is zero below 'cutoff' (yr^-1)
	double flatten; // model is flat below 'flatten'(yr^-1)
	double* data;   // data
	double tres;
	char mode;
} rednoisemodel_t;


rednoisemodel_t* setupRedNoiseModel(double start,double end, int npt, int nreal, double amp_1yr, double index);
void populateRedNoiseModel(rednoisemodel_t* model,long seed);
double getRedNoiseValue(rednoisemodel_t* model, double mjd,int real);
void freeRedNoiseModel(rednoisemodel_t* model);
double* getPowerSpectrum(rednoisemodel_t* model);
