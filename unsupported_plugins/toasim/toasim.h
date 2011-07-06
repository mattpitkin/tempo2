#ifndef TOASIM_H
#define TOASIM_H
#define TOASIM_STRLEN 1024
#define TOASIM_WRITER "libtoasim"
#define TOASIM_VERSION 2

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct toasim_header {
	uint32_t version;		// Arbitrary version number
	char writer[TOASIM_STRLEN];	// Name of the writing library
	char timfile_name[TOASIM_STRLEN];// Name of the idealised tim file used
	char parfile_name[TOASIM_STRLEN];// Name of the par file used
	char invocation[TOASIM_STRLEN];	// command used to generate this toasim file
	char short_desc[TOASIM_STRLEN];	// A short description of the correction
	char *description;		// A long description of the correction & parameters
	char *idealised_toas;		// The origianal .tim file
	char *orig_parfile;	// The origianal .par file
	char *gparam_desc;		// Description of global parameters
	char *gparam_vals;		// global parameter values
	char *rparam_desc;		// Description of per-realisation parameters
	uint32_t rparam_len;
	int64_t seed;			// The random seed used.
	uint32_t ntoa;			// The number of toas per realisation
	uint32_t nrealisations;		// The number of realisations generated
	uint32_t d_start;		// The byte offset that data begins
	uint32_t d_offset;		// The byte offset between realisations
} toasim_header_t;


typedef struct toasim_corrections {
	double *offsets;		// Absolute offsets in units of 10^-scale s.
	double a0;
	double a1;			// values of quadratic removed.
	double a2; 			// f(x) = a0 + a1*x + a2*x^2
	char *params;			// per-realisation parameter values
} toasim_corrections_t;

FILE *toasim_write_header(toasim_header_t *toasim_header, char* filename);
void *toasim_write_corrections(toasim_corrections_t* corr, toasim_header_t* header, FILE* file);
void *toasim_write_corrections_array(double* offsets,double a0, double a1, double a2, char* param, toasim_header_t* header, FILE* file);
void toasim_free_corrections(toasim_corrections_t* corr);


toasim_header_t *toasim_init_header();
toasim_header_t *toasim_read_header(FILE *file);
toasim_corrections_t *toasim_read_corrections(toasim_header_t *header, int nreal, FILE *file);

#ifdef __cplusplus
}
#endif

#endif
