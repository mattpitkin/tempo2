#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "toasim.h"


toasim_header_t *toasim_init_header(){
	toasim_header_t *toasim_header=(toasim_header_t*)malloc(sizeof(toasim_header_t));
	toasim_header->version=TOASIM_VERSION;
	strcpy(toasim_header->writer,TOASIM_WRITER);
	strcpy(toasim_header->timfile_name,"");
	strcpy(toasim_header->parfile_name,"");
	strcpy(toasim_header->invocation,"");
	strcpy(toasim_header->short_desc,"");
	toasim_header->description=toasim_header->short_desc;
	toasim_header->idealised_toas=toasim_header->short_desc;
	toasim_header->orig_parfile=toasim_header->short_desc;
	toasim_header->gparam_desc=toasim_header->short_desc;
	toasim_header->gparam_vals=toasim_header->short_desc;
	toasim_header->rparam_desc=toasim_header->short_desc;
	toasim_header->seed=0;
	toasim_header->ntoa=0;
	toasim_header->nrealisations=0;
	toasim_header->d_start=0;
	toasim_header->d_offset=0;
	return toasim_header;
}


void toasim_free_corrections(toasim_corrections_t* corr){
	if(corr->offsets!=NULL)free(corr->offsets);
	if(corr->params!=NULL)free(corr->params);
	free(corr);
}
