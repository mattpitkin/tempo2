
void T2Predictor_GetPhase_array_ld(const T2Predictor *t2p, long double *mjd, int nmjd, long double freq, double* out);
void T2Predictor_GetFrequency_array_ld(const T2Predictor *t2p, long double *mjd, int nmjd, long double freq, double* out);


void T2Predictor_GetPhase_array(const T2Predictor *t2p, double refmjd, double *mjd, const int nmjd, double freq, double* out);
void T2Predictor_GetFrequency_array(const T2Predictor *t2p, double refmjd, double *mjd, const int nmjd, double freq, double* out);
