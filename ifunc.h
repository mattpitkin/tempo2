
/**
 * Compute an ifunc gradient for a given 'k'
 */
double ifunc(const double* mjd, const double t,const int N, const int k);

/**
 * Compute an ifunc summed over all elements.
 */
double ifunc(const double* mjd, const double* yoffs, const double t,const int N);
double sinfunc(const double *T, const double t, const int k);


