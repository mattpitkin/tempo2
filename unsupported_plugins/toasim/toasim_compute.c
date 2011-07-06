#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "toasim.h"





double* toasim_apply_corrections(toasim_header_t* header, toasim_corrections_t* corr, double* vals){
	uint32_t i;
	if(vals==NULL){
		vals=(double*)calloc(header->ntoa,sizeof(double));
	}
	for (i=0; i < header->ntoa; i++){
		vals[i] += corr->offsets[i];
	}
	return vals;
}


