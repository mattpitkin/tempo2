#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "toasim.h"


int main(int argc, char** argv){
	toasim_header_t* header;
	FILE* file;
	int i;
	char printtim=0;
	char printpar=0;
	char dumpit=0;

	if (argc < 2){
		printf("Need a toasim file to work on!\n");
		printf("Usage: %s toasim_file [--par|--tim]\n",argv[0]);
		exit(1);
	}

	// read the file
	file=fopen(argv[1],"r");
	header = toasim_read_header(file);

	if (header == NULL){
		printf("\nFailed to parse header\n");
		exit(1);
	}

	for(i=2;i<argc;i++){
		if(strcmp(argv[i],"--par")==0)printpar=1;
		if(strcmp(argv[i],"--tim")==0)printtim=1;
		if(strcmp(argv[i],"--dump")==0)dumpit=1;

	}

	if(printtim && header->idealised_toas != NULL){
		printf("%s\n",header->idealised_toas);
		exit(0);
	}
	if(printpar && header->orig_parfile != NULL){
		printf("%s\n",header->orig_parfile);
		exit(0);
	}




	printf("Read TOASim (libtoasim version %d)\n",TOASIM_VERSION);
	printf("\n");
	printf("========================\n");
	printf("Header values as read\n");
	printf("========================\n");
	printf("version             : %d\n",header->version);
	printf("writer              : %s\n",header->writer);
	printf("timfile             : %s\n",header->timfile_name);
	printf("parfile             : %s\n",header->parfile_name);
	printf("invocation          : %s\n",header->invocation);
	printf("description         : %s\n",header->short_desc);
	printf("seed                : %d\n",header->seed);
	printf("ntoa                : %d\n",header->ntoa);
	printf("nrealisations       : %d\n",header->nrealisations);
	printf("========================\n");
	printf("data_start          : %d\n",header->d_start);
	printf("data_offset         : %d\n",header->d_offset);
	printf("========================\n\n\n");

	FILE * dump;
	if(dumpit)dump=fopen("toasim.txt","w");
	int good=0;
	for (i = 0; i < header->nrealisations; i++){
		toasim_corrections_t *read_corr= toasim_read_corrections(header,i,file);
		if (dumpit){
		   int j;
		   for(j=0;j<header->ntoa; j++){
			fprintf(dump,"%d %lg\n",i,read_corr->offsets[j]);
		   }
		fprintf(dump,"\n",read_corr->offsets[j]);
		}
		if (read_corr == NULL){
			printf("Failed to parse corrections %d\n",i);
		} else {
			toasim_free_corrections(read_corr);
			good++;
		}
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("Read %03d/%03d/%03d realisations OK.",good,i,header->nrealisations);
		fflush(stdout);
	}
	if(dumpit)fclose(dump);

	printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	printf("Read %d/%d realisations OK.         \n",good,header->nrealisations);

	fclose(file);
}

