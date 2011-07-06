#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "toasim.h"

int fwrite_err(const void *ptr, size_t size, size_t count, FILE *stream){
	int c=fwrite(ptr,size,count,stream);
	if(c!=count){
		fprintf(stderr,"count %d != %d %d\n",count,c,size);
		perror("toasim: write error");
	}
}

int fread_err(void *ptr, size_t size, size_t count, FILE *stream){
	if(fread(ptr,size,count,stream)!=count){
		if (ferror(stream)){
			perror("toasim: read error");
		}
		if (feof(stream)){
			fprintf(stderr,"toasim: EOF reached!\n");
		}
	}
}

void toasim_read_strx(char** val, FILE* file){
	char key[8];
	uint32_t len;
	key[4]='\0';
	fread_err(key,4,1,file);
	fread_err(&len,4,1,file);
	(*val)=(char*)malloc(len+1);
	fread_err(*val,1,len,file);
	(*val)[len]='\0';
}

void toasim_read_str(char* val, FILE* file){
	char key[8];
	uint32_t len;
	key[4]='\0';
	fread_err(key,4,1,file);
	fread_err(&len,4,1,file);
	if(len > TOASIM_STRLEN){
		fprintf(stderr,"ERROR: TOASIM_STRLEN not large enough for this data file\n");
		len=TOASIM_STRLEN;
	}
	fread_err(val,1,len,file);
}
void toasim_read_32(void* val, FILE* file){
	char key[8];
	uint32_t len;
	key[4]='\0';
	fread_err(key,4,1,file);
	fread_err(&len,4,1,file);
	if (len!=4){
		fprintf(stderr,"Error reading key '%s' expected 64 bit value, got 32 bit value\n",key);
	}

	fread_err(val,len,1,file);
}
void toasim_read_64(void* val, FILE* file){
	char key[8];
	uint32_t len;
	key[4]='\0';
	fread_err(key,4,1,file);
	fread_err(&len,4,1,file);
	if (len!=8){
		fprintf(stderr,"Error reading key '%s' expected 64 bit value, got 32 bit value\n",key);
	}
	fread_err(val,len,1,file);
}





void toasim_write_str(char* fname, char* str,FILE *file){
	uint32_t sz=(uint32_t)strlen(str);
	fwrite_err(fname,4,1,file);
	fwrite_err(&sz,sizeof(uint32_t),1,file);
	fwrite_err(str,1,strlen(str),file);
}

void toasim_write_32(char* fname, void* val, FILE *file){
	uint32_t sz=4;
	fwrite_err(fname,4,1,file);
	fwrite_err(&sz,4,1,file);
	fwrite_err(val,sz,1,file);
}
void toasim_write_64(char* fname, void* val, FILE *file){
	uint32_t sz=8;
	fwrite_err(fname,4,1,file);
	fwrite_err(&sz,4,1,file);
	fwrite_err(val,sz,1,file);
}
/**
 *
 * Takes a toasim header, and a filename.
 * Returns the newly created file as a filepointer.
 * Does not close the file
 */
FILE *toasim_write_header(toasim_header_t *toasim_header, char* filename){
	FILE *file;
	file = fopen(filename,"w");
	// If the file cannot be opened return null.
	// Caller should check error status.
	if (file==NULL) return NULL;

	toasim_header->version=TOASIM_VERSION;
	strcpy(toasim_header->writer,TOASIM_WRITER);

	fwrite_err("TOASIM",6,1,file);
	toasim_write_32("VERS",&toasim_header->version,file);
	toasim_write_str("WRTR",toasim_header->writer,file);
	toasim_write_str("T_NM",toasim_header->timfile_name,file);
	toasim_write_str("P_NM",toasim_header->parfile_name,file);
	toasim_write_str("INVK",toasim_header->invocation,file);
	toasim_write_str("SHRT",toasim_header->short_desc,file);
	toasim_write_str("DESC",toasim_header->description,file);
	toasim_write_str("TOAS",toasim_header->idealised_toas,file);
	toasim_write_str("OPAR",toasim_header->orig_parfile,file);
	toasim_write_str("GP_D",toasim_header->gparam_desc,file);
	toasim_write_str("GP_V",toasim_header->gparam_vals,file);
	toasim_write_str("RP_D",toasim_header->rparam_desc,file);
	toasim_write_32("RP_L",&toasim_header->rparam_len,file);
	toasim_write_64("SEED",&toasim_header->seed,file);
	toasim_write_32("NTOA",&toasim_header->ntoa,file);
	toasim_write_32("NREA",&toasim_header->nrealisations,file);

	toasim_header->d_start=ftell(file)+4*6;
	toasim_header->d_offset=
		4 + 		// Marker
		toasim_header->ntoa*8 +	// ntoas * double
		3*8 +		// three quadratic params
		toasim_header->rparam_len;	// parameters.
	toasim_write_32("DSTT",&toasim_header->d_start,file);
	toasim_write_32("DOFF",&toasim_header->d_offset,file);
	return file;
}

toasim_header_t *toasim_read_header(FILE *file){
	toasim_header_t *header;
	char key[8];
	int32_t dmy_32;
	uint32_t len=0;
	fread_err(key,6,1,file);
	key[6]='\0';
	if(strcmp(key,"TOASIM")){
		fprintf(stderr,"ERROR: not a TOASIM file (%s)\n",key);
		return NULL;
	}
	header = toasim_init_header();
	toasim_read_32(&header->version,file);
	if (header->version > TOASIM_VERSION){
		fprintf(stderr,"ERROR: toasim file version (%d) > library version (%d)\n",header->version,TOASIM_VERSION);
	}
	toasim_read_str(header->writer,file);
	toasim_read_str(header->timfile_name,file);
	toasim_read_str(header->parfile_name,file);
	toasim_read_str(header->invocation,file);
	toasim_read_str(header->short_desc,file);
	toasim_read_strx(&header->description,file);
	toasim_read_strx(&header->idealised_toas,file);
	toasim_read_strx(&header->orig_parfile,file);
	toasim_read_strx(&header->gparam_desc,file);
	toasim_read_strx(&header->gparam_vals,file);
	toasim_read_strx(&header->rparam_desc,file);
	toasim_read_32(&header->rparam_len,file);
	toasim_read_64(&header->seed,file);
	toasim_read_32(&header->ntoa,file);
	toasim_read_32(&header->nrealisations,file);
	if(header->version==1)toasim_read_32(&dmy_32,file); // old scale factor
	toasim_read_32(&header->d_start,file);
	toasim_read_32(&header->d_offset,file);

	return header;
}


void *toasim_write_corrections_array(double* offsets,double a0, double a1, double a2, char* param, toasim_header_t* header, FILE* file){
	toasim_corrections_t *corr = (toasim_corrections_t*) malloc(sizeof(toasim_corrections_t));
	corr->offsets = offsets;
	corr->a0=a0;
	corr->a1=a1;
	corr->a2=a2;
	if(header->rparam_len > 0){
		corr->params=(char*)malloc(header->rparam_len+1);
		if(strlen(param) > header->rparam_len){
			fprintf(stderr,"Warn: param string too long, truncating\n");
			memcpy(corr->params,param,header->rparam_len);
		} else{
			strcpy(corr->params,param);
		}
	} else{
		corr->params=NULL;
	}
	toasim_write_corrections(corr,header,file);
	if(corr->params!=NULL)free(corr->params);
	free(corr);
}


void *toasim_write_corrections(toasim_corrections_t* corr, toasim_header_t* header, FILE* file){
	char key[8];
	strcpy(key,"CORR");
	fwrite_err(key,4,1,file);
	fwrite_err(&corr->a0,8,1,file);
	fwrite_err(&corr->a1,8,1,file);
	fwrite_err(&corr->a2,8,1,file);
	fwrite_err(corr->offsets,sizeof(double), header->ntoa, file);
	if(header->rparam_len){
		fwrite_err(corr->params,1,header->rparam_len,file);
	}

}

toasim_corrections_t *toasim_read_corrections(toasim_header_t *header, int nreal, FILE *file){
	//seek to the correction requested
	char key[8];
	uint32_t offset=header->d_start+header->d_offset*nreal;

	if(fseek(file,offset,SEEK_SET)!=0) perror("fseek");
	key[4]='\0';
	fread_err(key,4,1,file);

	if(strcmp(key,"CORR")){
		fprintf(stderr,"ERROR: Could not locate CORR keyword at offset for realisation %d\n",nreal);
		return NULL;
	}

	toasim_corrections_t *corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));
	corr->offsets=(double*)malloc(sizeof(double)*header->ntoa);
	fread_err(&corr->a0,8,1,file);
	fread_err(&corr->a1,8,1,file);
	fread_err(&corr->a2,8,1,file);
	fread_err(corr->offsets,sizeof(double),header->ntoa,file);
	if(header->rparam_len > 0){
		corr->params=(char*)malloc(header->rparam_len);
		fread_err(corr->params,1,header->rparam_len+1,file);
		corr->params[header->rparam_len]='\0';
	}
	return corr;
}

