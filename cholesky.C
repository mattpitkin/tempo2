#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sstream>
#include <fstream>
#include <ctype.h>
#include <math.h>
#include "config.h"
#include "TKlongdouble.h"
#include "TKmatrix.h"
#include "choleskyRoutines.h"
#include "T2accel.h"
#include "cholesky.h"

#define LINE_LENGTH 2048

void cholesky_readT2Model2(double **m, std::stringstream &ss ,double *resx,double *resy,double *rese,int np, int nc,int *ip, pulsar *psr);
void cholesky_readT2Model1(double **m, std::stringstream &ss,double *resx,double *resy,double *rese,int np, int nc,int *ip, pulsar *psr);

void cholesky_sqexp_kernel(double **m, double modelA, double modelSigma, double *resx,double *resy,double *rese,int np, int nc);
void cholesky_matern15_kernel(double **m, double modelA, double modelSigma, double *resx,double *resy,double *rese,int np, int nc);
void cholesky_matern25_kernel(double **m, double modelA, double modelSigma, double *resx,double *resy,double *rese,int np, int nc);
void cholesky_cov2dm(double **mm, pulsar *psr, int *ip,int np,int nc);

void cholesky_powerlawQPModel(double **m, double modelAlpha, double modelFc, double modelA,double log_qp_ratio, double qp_f0, double sig, double lam, double *resx,double *resy,double *rese,int np, int nc);

void cholesky_powerlawRedPinkModel(double **m, double modelAlpha, double modelFc, double modelA,double alpha2, double fknee, double *resx,double *resy,double *rese,int np, int nc);

/*
 * Derive a covariance matrix from the given filename.
 *
 * input:
 * univ		(double[np][np]	empty matrix that will be filled with the "Cholesky matrix"
 * fname	(char[*])		Input filename to read
 * psr		(*pulsar)		pointer to pulsar structure that univ is used for
 * resx		(double[np])	X values of residuals (+ 0 for constraints)
 * resy		(double[np])	Y values of residuals (+ 0 for constraints)
 * rese		(double[np])	Error on residuals (+ eta for constraints)
 * np		(int)			Number of data points in fit (nres+nc)
 * nc		(int)			Number of constraints in fit
 * ip		(int)			Mapping from fit point to observation number in pulsar struct.
 */
void getCholeskyMatrix(double **cholesky_L, const char* fname, pulsar *psr, double *resx,double *resy,double *rese, int np, int nc,int* ip){
// UNUSED VARIABLE //     FILE* modelDescriptionFile;
    char modelFileName[1024];
// UNUSED VARIABLE //     char tmp[1024];
    double **m;
    int i,j;

    int nrows=get_blas_rows(cholesky_L);
    int ncols=get_blas_cols(cholesky_L);

    if (ncols!=np){
        logmsg("np=%d ncols=%d",ncols,np);
        logerr("cholesky_L error. Either you did not use malloc_uinv() to create cholesky_L or np!=ncols");
        exit(1);
    }

    if (nrows==1){
        // we are just adding the errors to a diagonal matrix.
        getCholeskyDiagonals(cholesky_L,psr,resx,resy, rese, np, nc,ip);
        return;
    }
    if (nrows!=np){
        logmsg("np=%d nrows=%d",ncols,np);
        logerr("cholesky_L error. Either you did not use malloc_cholesky_L() to create cholesky_L or np!=nrows");
        exit(1);
    }

    if(cholesky_L[np-1] != cholesky_L[0]+(np-1)*np){
        logerr("cholesky_L matrix not declared as consecutive memory.");
        logmsg("Please use malloc_cholesky_L() and free_cholesky_L() from cholesky.C");
#ifdef ACCEL_UINV
        // if we are using the accelerated code then it will crash... otherwise it's just a warning
        exit(1);
#endif
    }

    logmsg("Reading Cholesky model '%s' for pulsar %s",fname,psr->name);
    m=(double**)malloc(sizeof(double*)*(np+1));
    if(!m)logerr("Could not allocate enough memory");
    for(i=0;i<np+1;i++){
        m[i]=(double*)malloc(sizeof(double)*(np+1));
        if(!m[i])logerr("Could not allocate enough memory");
    }

    if (strcmp(fname,"PSRJ")==0){
        // this is the old method where the covariance function is read based on pulsar name.
        // There is no model description file.
        // This method is depricated
        logdbg("Reading covariance function using 'PSRJ' is now depricated");
        sprintf(modelFileName,"covarFunc.dat_%s",psr->name);
        cholesky_readFromCovarianceFunction(m,modelFileName,resx,resy,rese,np,nc);
    } else {
        cholesky_readT2CholModel(m,fname,resx,resy,rese,np,nc,ip,psr);
    }

    if(psr->ToAextraCovar!=NULL){
        logmsg("adding extra covar function to m");
        for (i=0;i<np;i++)
        {
            for (j=0;j<np;j++){
                m[i][j]+=psr->ToAextraCovar[i][j];
            }
        }
    }
    logdbg("mbefore = ");
    if(debugFlag){
    for (i=0;i<5;i++)
    { 
        for (j=0;j<5;j++) fprintf(LOG_OUTFILE,"%10g ",m[i][j]); 
        fprintf(LOG_OUTFILE,"\n");
    }
    fprintf(LOG_OUTFILE,"\n");
    }

    // make sure constraints are not covariant with anything.
    logdbg("Ensuring constraints have zero co-variance, np = %d, nc = %d",np,nc);
    for (i=np-nc; i < np; i++){
        for (j=0; j < np; j++){
            m[i][j]=0;
            m[j][i]=0;
        }
    }
    logdbg("Adding errors");
    for(i=0;i<np;i++){
        m[i][i]+=rese[i]*rese[i];
    }

    if (debugFlag)
    {
        logdbg("m = ");
        for (i=0;i<5;i++)
        { 
            for (j=0;j<5;j++) fprintf(LOG_OUTFILE,"%10g ",m[i][j]); 
            fprintf(LOG_OUTFILE,"\n");
        }
        fprintf(LOG_OUTFILE,"\n");
        FILE* mFile = fopen("chol.covarMatrix","w");
        for (i=0;i<np;i++)
        { 
            for (j=0;j<np;j++) fprintf(mFile,"%d %d %.28lg\n",i,j,m[i][j]); 
            fprintf(mFile,"\n");
        }
        fclose(mFile);
    }



    logdbg("Form cholesky_L from cholesky matrix 'm'");

    int ret = cholesky_formL(cholesky_L,m,np);
    if (ret!=0) {
        logwarn("Error with formL... adding 0.0001%% to diagonal (%.1lg + %.1lg)",m[0][0],m[0][0]*0.000001);
        for(i=0;i<np;i++){
            m[i][i] *= 1.000001;
        }
        ret = cholesky_formL(cholesky_L,m,np);
        if (ret!=0) {
            logerr("Error with formL");
            exit(ret);
        }
    }
    logdbg("compute determinant of cholesky_L");
    double det = 1;
    for (i=0;i<np;i++){
        det += log(cholesky_L[i][i]);
    }
    psr->detL=det;
    logdbg("det(cholesky_L)=%lg",det);


    for(i=0;i<np+1;i++)free(m[i]);
    free(m);
}


void cholesky_readT2CholModel_R(double **m, double **mm, const char* fname,double *resx,double *resy,double *rese,int np, int nc,int *ip, pulsar *psr,char *_psrJ,double _mjd_start,double _mjd_end,int recursion){
// UNUSED VARIABLE //     int flag;
    char val[LINE_LENGTH];
    char dmy[LINE_LENGTH];
    char key[LINE_LENGTH];
    char psrJ[LINE_LENGTH];
    double mjd_start=_mjd_start;
    double mjd_end=_mjd_end;


    std::stringstream content;

    bool first=true;

    while(fname[0]==' ' || fname[0]=='\t')fname++;

    if (strncmp(fname,"MODEL",5)==0){
        content << std::string(fname);
    } else {
        std::ifstream file(fname);
        content << file.rdbuf();
        file.close();
    }

    recursion+=1;


    strcpy(psrJ,_psrJ);

    std::string line;
    while(std::getline(content,line)){
        const char* tmp = line.c_str();
        if (first){
            //if the first character in the file is a numeric digit, assume that we have a DCF file.
            if(isdigit(tmp[0])){
                // assume that we have 
                cholesky_readFromCovarianceFunction(m,fname,resx,resy,rese,np,nc);
                return;
            } else{
                logmsg("Reading a Tempo2 MODEL file: '%s'",fname);
            }
            first=false;
        }
        while(tmp[0]==' ' || tmp[0]=='\t')tmp++;
        if (tmp[0]=='#') continue;
        int ok=sscanf(tmp,"%s",key);
        if (ok==-1)continue;
        logdbg("key[%d]: %s ('%s' '%s')\n",recursion,key,psrJ,psr->name);

        // now test what the keyword is...
        //
        if (strcmp(key,"CLEAR")==0 || strcmp(key,"END")==0){
            psrJ[0]='\0';
            mjd_start=0.0;
            mjd_end=99999.0;
            continue;
        }

        if (strcmp(key,"PSR")==0){
            sscanf(tmp,"%s %s",dmy,psrJ);
        }

        if(psrJ[0]!='\0' && strcmp(psrJ,psr->name))continue;
        if (strcmp(key,"MJD")==0){
            sscanf(tmp,"%s %lf %lf",dmy,&mjd_start,&mjd_end);
        }
        if (strcmp(key,"INCLUDE")==0){
            sscanf(tmp,"%s %s",dmy,val);
            cholesky_readT2CholModel_R(m,mm,val,resx,resy,rese,np,nc,ip,psr,psrJ,mjd_start,mjd_end,recursion);
        }
        if (strcmp(key,"DCF")==0){
            sscanf(tmp,"%s %s",dmy,val);
            cholesky_readFromCovarianceFunction(mm,val,resx,resy,rese,np,nc);
            addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
        }
        else if (strcmp(key,"ECM")==0){ // Extra covariance matrix
            sscanf(tmp,"%s %s",dmy,val);
            cholesky_ecm(mm,val,resx,resy,rese,np,nc);
            addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
        }

        if (strcmp(key,"MODEL")==0){
            sscanf(tmp,"%s %s",dmy,val);
            if (strcmp(val,"T2PowerLaw")==0){
                double alpha,amp,fc;
                sscanf(tmp,"%s %s %lf %lf %lf",dmy,val,&alpha,&amp,&fc);
                cholesky_powerlawModel(mm,alpha,fc,amp, resx, resy,rese,np, nc);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if (strcmp(val,"T2PowerLaw_QPc")==0){
                double alpha,amp,fc,log_qp_ratio,qp_f0,sig,lam;
                sscanf(tmp,"%s %s %lf %lf %lf %lf %lf %lf %lf",dmy,val,&alpha,&amp,&fc,&log_qp_ratio,&qp_f0,&sig,&lam);
                cholesky_powerlawQPModel(mm,alpha,fc,amp, log_qp_ratio,qp_f0,sig,lam, resx, resy,rese,np, nc);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if (strcmp(val,"T2PowerLaw_RedPink")==0){
                double alpha,amp,fc,alpha2,fknee;
                sscanf(tmp,"%s %s %lf %lf %lf %lf %lf",dmy,val,&alpha,&amp,&fc,&alpha2,&fknee);
                cholesky_powerlawRedPinkModel(mm,alpha,fc,amp, alpha2,fknee, resx, resy,rese,np, nc);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if (strcmp(val,"T2YrPowerLaw_X")==0){
                double alpha,amp,fc;
                sscanf(tmp,"%s %s %lf %lf",dmy,val,&alpha,&amp);
                fc=0.05;
                cholesky_powerlawModel(mm,alpha,fc,amp*pow(fc,-alpha), resx, resy,rese,np, nc);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
           } else if (strcmp(val,"T2YrPowerLaw")==0){
                double alpha,amp,fc;
                sscanf(tmp,"%s %s %lf %lf %lf",dmy,val,&alpha,&amp,&fc);
                cholesky_powerlawModel(mm,alpha,fc,amp*pow(fc,-alpha), resx, resy,rese,np, nc);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if (strcmp(val,"T2YrPowerLawDM")==0){
                double alpha,amp,fc;
                sscanf(tmp,"%s %s %lf %lf %lf",dmy,val,&alpha,&amp,&fc);
                cholesky_powerlawModel(mm,alpha,fc,amp*pow(fc,-alpha), resx, resy,rese,np, nc);
                cholesky_cov2dm(mm,psr,ip,np,nc);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if (strcmp(val,"T2PowerLawDM")==0){
                double alpha,amp,fc;
                sscanf(tmp,"%s %s %lf %lf %lf",dmy,val,&alpha,&amp,&fc);
                cholesky_powerlawModel(mm,alpha,fc,amp, resx, resy,rese,np, nc);
                cholesky_cov2dm(mm,psr,ip,np,nc);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if (strcmp(val,"T2SqExpDM")==0){
                double amp,sigma;
                sscanf(tmp,"%s %s %lf %lf",dmy,val,&amp,&sigma);
                cholesky_sqexp_kernel(mm,amp,sigma, resx, resy,rese,np, nc);
                cholesky_cov2dm(mm,psr,ip,np,nc);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if (strcmp(val,"T2SqExp")==0){
                double amp,sigma;
                sscanf(tmp,"%s %s %lf %lf",dmy,val,&amp,&sigma);
                cholesky_sqexp_kernel(mm,amp,sigma, resx, resy,rese,np, nc);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if (strcmp(val,"T2Matern1.5")==0){
                double amp,sigma;
                sscanf(tmp,"%s %s %lf %lf",dmy,val,&amp,&sigma);
                cholesky_matern15_kernel(mm,amp,sigma, resx, resy,rese,np, nc);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if (strcmp(val,"T2Matern1.5DM")==0){
                double amp,sigma;
                sscanf(tmp,"%s %s %lf %lf",dmy,val,&amp,&sigma);
                cholesky_matern15_kernel(mm,amp,sigma, resx, resy,rese,np, nc);
                cholesky_cov2dm(mm,psr,ip,np,nc);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if (strcmp(val,"T2Matern2.5")==0){
                double amp,sigma;
                sscanf(tmp,"%s %s %lf %lf",dmy,val,&amp,&sigma);
                cholesky_matern25_kernel(mm,amp,sigma, resx, resy,rese,np, nc);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if (strcmp(val,"T2Matern2.5DM")==0){
                double amp,sigma;
                sscanf(tmp,"%s %s %lf %lf",dmy,val,&amp,&sigma);
                cholesky_matern25_kernel(mm,amp,sigma, resx, resy,rese,np, nc);
                cholesky_cov2dm(mm,psr,ip,np,nc);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if(strcmp(val,"DM")==0){

                double D,d,ref_freq;
                sscanf(tmp,"%s %s %lf %lf %lf",dmy,val,&D,&d,&ref_freq);
                cholesky_dmModel(mm,D,d,ref_freq,resx, resy,rese,np, nc);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if(strcmp(val,"DMCovarParam")==0){ // Added by Daniel Reardon and George Hobbs
                double alpha,a,b;
                sscanf(tmp,"%s %s %lf %lf %lf",dmy,val,&alpha,&a,&b);
                cholesky_dmModelCovarParam(mm,alpha,a,b,resx, resy,rese,np, nc);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if (strcmp(val,"1")==0){
                cholesky_readT2Model1(mm,content,resx,resy,rese,np,nc,ip,psr);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if (strcmp(val,"2")==0){
                cholesky_readT2Model2(mm,content,resx,resy,rese,np,nc,ip,psr);
                addCovar(m,mm,resx,resy,rese,np,nc,ip,psr,mjd_start,mjd_end);
                continue;
            } else if (strcmp(val,"T2")==0){
                continue;
            }
            logerr("Model type '%s' not understood",val);
            exit(1);
        }

    }

}

void cholesky_readFromCovarianceFunction(double **m, const char* fname,double *resx,double *resy,double *rese,int np, int nc){
    int ndays = ceil((resx[np-1-nc]-resx[0])+1e-10);
    double covarFunc[ndays+1];
    double escaleFactor = 1.0;
    int i;
    FILE* fin;
    logmsg("Parsing '%s' as covariance function",fname);

    logtchk("reading covariance function from disk, ndays = %d",ndays);
    if (!(fin = fopen(fname,"r")))
    {
        logerr("Unable to open covariance function file: %s",fname);
        exit(1);
    }
    logdbg("ndays = %d",ndays);
    fscanf(fin,"%lf",&escaleFactor);
    for (i=0;i<=ndays;i++)
    {
        if (fscanf(fin,"%lf",&covarFunc[i])!=1)
        {
            logerr("Incorrect number of days in the Cholesky matrix: %s, trying to read %d days",fname,ndays);
            exit(1);
        }
    }
    fclose(fin);
    if(escaleFactor!=1.0){
        logerr("Ignoring 'error scale factor': %g",escaleFactor);
    }
    logtchk("complete reading covariance function from disk");

    cholesky_covarFunc2matrix(m,covarFunc,ndays,resx,resy,rese,np,nc);

}

void cholesky_covarFunc2matrix(double** m, double* covarFunc, int ndays,double *resx,double *resy,double *rese,int np, int nc){
// UNUSED VARIABLE //     double escaleFactor = 1.0;
    int i,j;
    int ix,iy;
    double t0,cint,t;
    int t1,t2;


    logtchk("forming Cholesky matrix ... determing m[ix][iy] = fabs(resx[ix]-resx[iy])");
    for (ix=0;ix<(np);ix++)
    {
        for (iy=0;iy<(np);iy++)
            m[ix][iy] = fabs(resx[ix]-resx[iy]);
    }
    logtchk("forming Cholesky matrix ... complete determing m[ix][iy] = fabs(resx[ix]-resx[iy])");
    if (debugFlag==1)
    {
        logdbg("First m = ");
        for (i=0;i<5;i++)
        { 
            for (j=0;j<5;j++) fprintf(LOG_OUTFILE,"%10g ",m[i][j]); 
            fprintf(LOG_OUTFILE,"\n");
        }
        fprintf(LOG_OUTFILE,"\n");
        logdbg("CovarFunc = ");
        for (i=0;i<10;i++)
        { 
            fprintf(LOG_OUTFILE,"%10g\n",covarFunc[i]); 
        }

    }

    // Insert the covariance which depends only on the time difference.
    // Linearly interpolate between elements on the covariance function because
    // valid covariance matrix must have decreasing off diagonal elements.
    // logdbg("Inserting into the covariance matrix");
    logtchk("forming Cholesky matrix ... determing covariance based on time difference");
    for (ix=0;ix<(np);ix++)
    {
        for (iy=0;iy<(np);iy++)
        {
            if (ix >= np-nc || iy >= np-nc)
            {
                m[ix][iy] = 0;
            }
            else
            {
                t0 = m[ix][iy];
                t1 = (int)floor(t0);
                t2 = t1+1;
                t  = t0-t1;
                if (t1 > ndays || t2 > ndays)
                {
                    logerr("Problem that t1 or t2 > ndays: t1 = %d, t2 = %d, ndays = %d, ix = %d, iy = %d, np = %d",t1,t2,ndays,ix,iy,np);
                    exit(1);
                }
                cint = covarFunc[t1]*(1-t)+covarFunc[t2]*t; // Linear interpolation
                m[ix][iy] = cint;
            }
        }
    }
    logtchk("forming Cholesky matrix ... complete determing covariance based on time difference");
    // add the values for the constraints
    // Constraints are not covariant with anything so it's all zero!
    for (i=np-nc; i < np; i++){
        for (j=0; j < np; j++){
            m[i][j]=0;
            m[j][i]=0;
        }
    }
}



void getCholeskyDiagonals(double **L, pulsar *psr, double *resx,double *resy,double *rese, int np, int nc,int* ip){
    int i;
    double det = 0;
    for(i=0;i<np;i++){
        L[0][i]=rese[i];
        det += log(L[0][i]);
    }
    psr->detL=det;
    logdbg("det(L)=%lg (C diagonal)",det);

}


/**
 * L is a lower triangluar matrix.
 * Matricies are row-major order, i.e. cholesky_L[r][c].
 * returns 0 if ok.
 *
 * Note that this matrix used to be called "uinv", and was
 * the inverse of the lower-triangular matrix.
 * 
 * Since 2021 we no longer invert the matrix and instead use it
 * with forward substitution as this is much faster. Rather than call
 * the non-inverted lower triangular matrix "uinv" I (MJK) decided to
 * re-name it "L" everywhere!. Hope this is not too confusing.
 */
int cholesky_formL(double **cholesky_L,double** m,int np){
    int i,j,k;
    logtchk("forming Cholesky matrix ... do Cholesky decomposition");
#ifdef ACCEL_UINV
    if(useT2accel){
        logdbg("Doing ACCELERATED Chol Decomp (M.Keith/LAPACK method)");
        for(i =0;i<np;i++){
            memcpy(cholesky_L[i],m[i],np*sizeof(double));
        }
        logdbg("NO_UINV!! (M.Keith/LAPACK method)");
        int ret = accel_cholfac(cholesky_L[0],np);
        logtchk("forming Cholesky matrix ... complete calculate U");
        if (ret != 0) return ret;
    } else {
#endif
        double sum;

        double *cholp  = (double *)malloc(sizeof(double)*(np+1));
        logmsg("Doing Cholesky decomp and inverting matrix (SLOW method)");
        if(!cholp)logerr("Could not allocate enough memory");

        T2cholDecomposition(m,np,cholp);
        logtchk("forming Cholesky matrix ... complete do Cholesky decomposition");
        // We have to construct L from the output of the above.
        for (i=0;i<np;i++)
        {
            for (j=i;j>=0;j--) { // just get the lower triangular part...
                cholesky_L[i][j] = m[i][j];
            }
            cholesky_L[i][i]=cholp[i]; // and the diagonal terms
        }
        logtchk("forming Cholesky matrix ... complete calculate cholesky L");
        if (debugFlag)
        {
            logdbg("L = ");
            for (i=0;i<5;i++)
            { 
                for (j=0;j<5;j++) fprintf(LOG_OUTFILE,"%10g ",cholesky_L[i][j]); 
                fprintf(LOG_OUTFILE,"\n");
            }
            fprintf(LOG_OUTFILE,"\n");
        }


        logtchk("forming Cholesky matrix ... free memory");
        free(cholp);
        logtchk("forming Cholesky matrix ... complete free memory");

#ifdef ACCEL_UINV
    } // end the if clause when we have the option of accelerated cholesky.
#endif

    if(debugFlag){
        logdbg("Write cholesky_L");
        FILE* file=fopen("chol.L","w");
        for(i =0;i<np;i++){
            for(j =0;j<np;j++){
                fprintf(file,"%d %d %lg\n",i,j,cholesky_L[i][j]);
            }
            fprintf(file,"\n");
        }
        fclose(file);
    }
    return 0;

}


int cholesky_formLinv(double **Linv,double** m,int np){
    int i,j,k;
    double sum;
    double **L = malloc_uinv(np);
    double *cholp = (double*)malloc(sizeof(double)*np);
    cholesky_formL(L,m,np);
    for (i=0;i<np;i++){
        cholp[i]=L[i][i];
        L[i][i]=0;
    }
    logtchk("forming Cholesky matrix ... calculate uinv");
    for (i=0;i<np;i++)
    {
        L[i][i] = 1.0/cholp[i];
        Linv[i][i] = L[i][i];
        for (j=0;j<i;j++)
            Linv[j][i] = 0.0;
        for (j=i+1;j<np;j++)
        {
            sum=0.0;
            for (k=i;k<j;k++) sum-=L[j][k]*L[k][i];
            L[j][i]=sum/cholp[j];
            Linv[j][i] = L[j][i];
        }
    }
    free(cholp);
    free_uinv(L);
    return 0;
}


void cholesky_readT2Model1(double **m, std::stringstream &ss,double *resx,double *resy,double *rese,int np, int nc,int *ip, pulsar *psr){
    char dmy[LINE_LENGTH];
    double alpha;
    double fc;
    double amp;
    std::string s;
    std::getline(ss,s);
    sscanf(s.c_str(),"%s %lg\n",dmy,&alpha);
    std::getline(ss,s);
    sscanf(s.c_str(),"%s %lg\n",dmy,&fc);
    std::getline(ss,s);
    sscanf(s.c_str(),"%s %lg\n",dmy,&amp);
    cholesky_powerlawModel(m,alpha,fc,amp, resx, resy,rese,np, nc);
}

void cholesky_readT2Model2(double **m, std::stringstream &ss ,double *resx,double *resy,double *rese,int np, int nc,int *ip, pulsar *psr){
    char dmy[LINE_LENGTH];
    double alpha;
    double beta;
    double fc;
    double amp;

    std::string s;
    std::getline(ss,s);
    sscanf(s.c_str(),"%s %lg\n",dmy,&alpha);
    std::getline(ss,s);
    sscanf(s.c_str(),"%s %lg\n",dmy,&beta);
    std::getline(ss,s);
    sscanf(s.c_str(),"%s %lg\n",dmy,&fc);
    std::getline(ss,s);
    sscanf(s.c_str(),"%s %lg\n",dmy,&amp);
    cholesky_powerlawModel_withBeta(m,alpha,beta,fc,amp, resx, resy,rese,np, nc);
}

void cholesky_ecm(double **m, char* fileName,double *resx,double *resy,double *rese,int np, int nc){
// UNUSED VARIABLE //     char dmy[LINE_LENGTH];
    int i,j;
    FILE *fin;
    double dummy;

    if (!(fin = fopen(fileName,"r")))
    {
        printf("ERROR: Unable to open ECM file: %s\n",fileName);
        exit(1);
    }

    //   cholesky_powerlawModel(m,alpha,fc,amp, resx, resy,rese,np, nc);
    logmsg("Using extra covariance matrix from: %s\n",fileName);
    for (i=0;i<np+nc;i++)
    {
        for (j=0;j<np+nc;j++)
        {
            if (fscanf(fin,"%lf",&m[i][j])!=1)
            {
                logwarn("Error reading element (%d,%d) from %s\n",i,j,fileName);
                exit(1);
            }
        }
    }
    if (fscanf(fin,"%lf",&dummy)==1)
    {
        logwarn("WARNING: %s seems to have too many elements\n",fileName);
    }
    fclose(fin);
}



void cholesky_readT2CholModel(double **m, const char* fname,double *resx,double *resy,double *rese,int np, int nc,int *ip, pulsar *psr){
    char psrJ[LINE_LENGTH];
    int i,j;
    psrJ[0]='\0';
    // mm is a duplicate covariance matrix that is added on top of m each time.
    double** mm=(double**)malloc(sizeof(double*)*(np+1));
    if(!mm){
        logerr("Could not allocate enough memory");
        exit(1);
    }
    for(i=0;i<np+1;i++){
        mm[i]=(double*)malloc(sizeof(double)*(np+1));
        if(!mm[i]){
            logerr("Could not allocate enough memory");
            exit(1);
        }
        for(j=0;j<np+1;j++){
            m[i][j]=0;
            mm[i][j]=0;
        }
    }


    cholesky_readT2CholModel_R(m,mm,fname,resx,resy,rese,np,nc,ip,psr,psrJ,0.0,99999.0,0);
    for(i=0;i<np+1;i++)free(mm[i]);
    free(mm);


}



void addCovar(double **m,double **mm,double *resx,double *resy,double *rese,int np, int nc,int *ip, pulsar *psr,double mjd_start,double mjd_end){
// UNUSED VARIABLE //     int *goodvals;
    int i,j;

    int istart=0;
    int iend=0;

    // if we aren't really doing a pulsar, e.g. model for ifunc or dmmodel
    // Then we don't have ip, and don't deal with start and end
    // Maybe we can fix this better in future, but for now this should be fine.
    if (ip==NULL){
    logdbg("Adding matrix m to mm (MJD range disabled as ip=NULL)");
        for(i=0;i<np-nc;i++){
            for(j=0;j<np-nc;j++){
                m[i][j]+=mm[i][j];
            }
         }
    } else {

    logdbg("Adding matrix m to mm (MJD range: %lf -> %lf) %lf %lf",mjd_start,mjd_end,(double)psr->obsn[ip[0]].sat,(double)psr->obsn[ip[np-nc-1]].sat);
    for(i=0;i<np-nc;i++){
        if(psr->obsn[ip[i]].sat>mjd_start){
            istart=i;
            break;
        }
    }
    for(i=np-nc-1; i>=0;i--){
        if(psr->obsn[ip[i]].sat < mjd_end){
            iend=i;
            break;
        }
    }

    logdbg("istart=%d iend=%d (%d)",istart,iend,np-nc);


    /****
     *
     * Needs to have point closest to start/end as reference point for each way.
     *
     */
    for(i=0;i<np-nc;i++){
        for(j=0;j<np-nc;j++){
            if( (psr->obsn[ip[j]].sat > mjd_start && psr->obsn[ip[j]].sat < mjd_end)
                    && (psr->obsn[ip[i]].sat > mjd_start && psr->obsn[ip[i]].sat < mjd_end))
                m[i][j]+=mm[i][j]; // 1
            else{
                if((psr->obsn[ip[j]].sat < mjd_start) && (psr->obsn[ip[i]].sat < mjd_start))
                    m[i][j]+=mm[istart][istart]; // 2
                else if((psr->obsn[ip[j]].sat > mjd_end) && (psr->obsn[ip[i]].sat > mjd_end))
                    m[i][j]+=mm[iend][iend]; // 3
                else if((psr->obsn[ip[j]].sat < mjd_start) && (psr->obsn[ip[i]].sat > mjd_end))
                    m[i][j]+=mm[istart][iend]; //4
                else if((psr->obsn[ip[j]].sat > mjd_end) && (psr->obsn[ip[i]].sat < mjd_start))
                    m[i][j]+=mm[istart][iend]; //5
                else if((psr->obsn[ip[j]].sat > mjd_end)
                        && ((psr->obsn[ip[i]].sat > mjd_start ) && (psr->obsn[ip[i]].sat < mjd_end)))
                    m[i][j]+=mm[i][iend]; // 6
                else if((psr->obsn[ip[j]].sat < mjd_start)
                        && ((psr->obsn[ip[i]].sat > mjd_start ) && (psr->obsn[ip[i]].sat < mjd_end)))
                    m[i][j]+=mm[i][istart]; // 7
                else if((psr->obsn[ip[i]].sat > mjd_end)
                        && ((psr->obsn[ip[j]].sat > mjd_start ) && (psr->obsn[ip[j]].sat < mjd_end)))
                    m[i][j]+=mm[iend][j]; // 8
                else if((psr->obsn[ip[i]].sat < mjd_start)
                        && ((psr->obsn[ip[j]].sat > mjd_start ) && (psr->obsn[ip[j]].sat < mjd_end)))
                    m[i][j]+=mm[istart][j]; // 9
            }

        }
    }
    }
}


// a in seconds^2
// b in days
//
void cholesky_dmModelCovarParam(double **m, double alpha, double a, double b,double *resx,double *resy,double *rese,int np, int nc){
// UNUSED VARIABLE //     double secperyear=365.25*86400.0;
// UNUSED VARIABLE //     double tobs=ceil((resx[np-1])-(resx[0]))/365.25;
    double *covarFunc;
    double x;
    int i;
    int ndays=ceil((resx[np-1])-(resx[0])+1e-10);
    covarFunc=(double*)malloc(sizeof(double)*(ndays+1));

    logmsg("In DMCovarParam function\n");
    for (i=0; i <= ndays; i++){
        x = (i+1e-10);
        covarFunc[i]=a*exp(-pow(x/b,alpha));
        logmsg("DMCovarParam: %g\n",covarFunc[i]);
    }
    cholesky_covarFunc2matrix(m,covarFunc,ndays,resx,resy,rese,np,nc);


    free(covarFunc);
}


void cholesky_dmModel(double **m, double D_d, double d, double ref_freq,double *resx,double *resy,double *rese,int np, int nc){
    double secperyear=365.25*86400.0;
    double tobs=ceil((resx[np-1])-(resx[0]))/365.25;
    double alpha=8.0/3.0;
    D_d *=1e-12; // convert us to seconds
    d*=86400.0;  // convert days to seconds.

    logdbg("D_d(%f)  = %f (us^2)\n",d,D_d    );
    logdbg("ref_freq = %f (MHz)\n",ref_freq    );
    logmsg("Warning: cholesky DM model currently ignores ref_freq");


    //power at 1 year
    //Equation 12 from Keith et al. (2012)
    // 0.0112 scales a structure function into a power spectrum
    double pism = 0.0112 * D_d * pow(d,(-5.0/3.0)) * pow(secperyear,-1.0/3.0);
    logdbg("pism(1yr^-1)  = %g (yr^3)",pism);

    double fc=1/tobs;


    // power at fc
    pism*=pow(fc,-alpha);

    logdbg("pism(%lfyr^-1)  = %g (yr^3)",fc,pism);

    cholesky_powerlawModel(m,alpha,fc,pism,resx,resy,rese,np,nc);

}

void cholesky_matern15_kernel(double **m, double modelA, double modelSigma, double *resx,double *resy,double *rese,int np, int nc){

    int ndays;
    double *covarFunc;
    logmsg("Generate covar matrix from Matern Kernel model (A=%lg sigma=%lf nu=1.5) (np=%d, nc=%d)",modelA,modelSigma,np,nc);
    ndays=ceil((resx[np-1])-(resx[0])+1e-10);
    covarFunc=(double*)malloc(sizeof(double)*(ndays+1));

    double x;
    double k1 = sqrt(3.0)/modelSigma;
    for (int i=0; i <= ndays ; ++i){
        x=i;
        covarFunc[i] = modelA*(1.0 + k1*x) * exp(-x*k1);
    }

    if(debugFlag){
        FILE* outFile = fopen("newDCF","w");
        for(int i=0;i<ndays;i++){
            fprintf(outFile,"%lg\n",covarFunc[i]);
        }
        fclose(outFile);
    }

    cholesky_covarFunc2matrix(m,covarFunc,ndays,resx,resy,rese,np,nc);
    free(covarFunc);
}
void cholesky_matern25_kernel(double **m, double modelA, double modelSigma, double *resx,double *resy,double *rese,int np, int nc){

    int ndays;
    double *covarFunc;
    logmsg("Generate covar matrix from Matern Kernel model (A=%lg sigma=%lf nu=2.5) (np=%d, nc=%d)",modelA,modelSigma,np,nc);
    ndays=ceil((resx[np-1])-(resx[0])+1e-10);
    covarFunc=(double*)malloc(sizeof(double)*(ndays+1));

    double x;
    double k1 = sqrt(5.0)/modelSigma;
    double k2 = 5.0/modelSigma/modelSigma/3.0;
    for (int i=0; i <= ndays ; ++i){
        x=i;
        covarFunc[i] = modelA*(1.0 + k1*x + k2*x*x) * exp(-x*k1);
    }

    if(debugFlag){
        FILE* outFile = fopen("newDCF","w");
        for(int i=0;i<ndays;i++){
            fprintf(outFile,"%lg\n",covarFunc[i]);
        }
        fclose(outFile);
    }

    cholesky_covarFunc2matrix(m,covarFunc,ndays,resx,resy,rese,np,nc);
    free(covarFunc);
}



void cholesky_sqexp_kernel(double **m, double modelA, double modelSigma, double *resx,double *resy,double *rese,int np, int nc){

    int ndays;
    double *covarFunc;
    logmsg("Generate covar matrix from Square Exponential model (A=%lg sigma=%lf) (np=%d, nc=%d)",modelA,modelSigma,np,nc);
    ndays=ceil((resx[np-1])-(resx[0])+1e-10);
    covarFunc=(double*)malloc(sizeof(double)*(ndays+1));

    double x;
    double K = 1.0/(2.*modelSigma*modelSigma);
    for (int i=0; i <= ndays ; ++i){
        x=i;
        covarFunc[i] = modelA*exp(-x*x*K);
    }

    if(debugFlag){
        FILE* outFile = fopen("newDCF","w");
        for(int i=0;i<ndays;i++){
            fprintf(outFile,"%lg\n",covarFunc[i]);
        }
        fclose(outFile);
    }

    cholesky_covarFunc2matrix(m,covarFunc,ndays,resx,resy,rese,np,nc);
    free(covarFunc);
}



void cholesky_powerlawModel(double **m, double modelAlpha, double modelFc, double modelA,double *resx,double *resy,double *rese,int np, int nc){

    int ndays,i;
    double *covarFunc;
    logmsg("Generate covar matrix from powerlaw model (a=%lf fc=%lf A=%lg) (np=%d, nc=%d)",modelAlpha,modelFc,modelA,np,nc);
    ndays=ceil((resx[np-1])-(resx[0])+1e-10);
    covarFunc=(double*)malloc(sizeof(double)*(ndays+1));
    int ndays_out = T2calculateCovarFunc(modelAlpha,modelFc,modelA,0,0,covarFunc,resx,resy,rese,np);
    if(ndays!=ndays_out){
        logerr("Ndays in != Ndays out!");
    }
    if(debugFlag){
        FILE* outFile = fopen("newDCF","w");
        for(i=0;i<ndays;i++){
            fprintf(outFile,"%lg\n",covarFunc[i]);
        }
        fclose(outFile);
    }

    cholesky_covarFunc2matrix(m,covarFunc,ndays,resx,resy,rese,np,nc);
    free(covarFunc);
}

void cholesky_powerlawRedPinkModel(double **m, double modelAlpha, double modelFc, double modelA,double alpha2, double fknee, double *resx,double *resy,double *rese,int np, int nc){

    int ndays,i;
    double *covarFunc;
    logmsg("Generate covar matrix from 'RedPink' powerlaw model (a=%lf fc=%lf A=%lg a2=%lg, fknee=%lg) (np=%d, nc=%d)",modelAlpha,modelFc,modelA,alpha2,fknee,np,nc);
    ndays=ceil((resx[np-1])-(resx[0])+1e-10);
    covarFunc=(double*)malloc(sizeof(double)*(ndays+1));
    int ndays_out = T2calculateCovarFunc_redpink(modelAlpha,modelFc,modelA,alpha2,fknee,covarFunc,resx,resy,rese,np);
    if(ndays!=ndays_out){
        logerr("Ndays in != Ndays out!");
    }
    if(debugFlag){
        FILE* outFile = fopen("newDCF","w");
        for(i=0;i<ndays;i++){
            fprintf(outFile,"%lg\n",covarFunc[i]);
        }
        fclose(outFile);
    }

    cholesky_covarFunc2matrix(m,covarFunc,ndays,resx,resy,rese,np,nc);
    free(covarFunc);
}
void cholesky_powerlawQPModel(double **m, double modelAlpha, double modelFc, double modelA,double log_qp_ratio, double qp_f0, double sig, double lam, double *resx,double *resy,double *rese,int np, int nc){

    int ndays,i;
    double *covarFunc;
    logmsg("Generate covar matrix from QP + powerlaw model (a=%lf fc=%lf A=%lg qpr=%lg, f0=%lg sig=%lg lam=%lg) (np=%d, nc=%d)",modelAlpha,modelFc,modelA,log_qp_ratio,qp_f0,sig,lam,np,nc);
    ndays=ceil((resx[np-1])-(resx[0])+1e-10);
    covarFunc=(double*)malloc(sizeof(double)*(ndays+1));
    int ndays_out = T2calculateCovarFunc_with_QP(modelAlpha,modelFc,modelA,log_qp_ratio,qp_f0,sig,lam,covarFunc,resx,resy,rese,np);
    if(ndays!=ndays_out){
        logerr("Ndays in != Ndays out!");
    }
    if(debugFlag){
        FILE* outFile = fopen("newDCF","w");
        for(i=0;i<ndays;i++){
            fprintf(outFile,"%lg\n",covarFunc[i]);
        }
        fclose(outFile);
    }

    cholesky_covarFunc2matrix(m,covarFunc,ndays,resx,resy,rese,np,nc);
    free(covarFunc);
}

void cholesky_powerlawModel_withBeta(double **m, double modelAlpha, double modelBeta,double modelFc, double modelA,double *resx,double *resy,double *rese,int np, int nc){

    int ndays,i;
    double *covarFunc;
    logmsg("Generate covar matrix from powerlaw model (a=%lf b=%lf fc=%lf A=%lg)",modelAlpha,modelBeta,modelFc,modelA);
    ndays=ceil((resx[np-1])-(resx[0])+1e-10);
    covarFunc=(double*)malloc(sizeof(double)*(ndays+1));
    int ndays_out = T2calculateCovarFunc(modelAlpha,modelFc,modelA,1,modelBeta,covarFunc,resx,resy,rese,np);
    if(ndays!=ndays_out){
        logerr("Ndays in != Ndays out!");
    }
    if(debugFlag){
        FILE* outFile = fopen("newDCF","w");
        for(i=0;i<ndays;i++){
            fprintf(outFile,"%lg\n",covarFunc[i]);
        }
        fclose(outFile);
    }

    cholesky_covarFunc2matrix(m,covarFunc,ndays,resx,resy,rese,np,nc);
    free(covarFunc);
}


// this could be made more efficient!
void cholesky_cov2dm(double **mm, pulsar *psr, int *ip,int np,int nc){
    for (int i = 0; i < np; ++i){
        double f1 = psr->obsn[ip[i]].freq;
        double g1 = 1.0/(DM_CONST*f1*f1);
        for (int j = 0; j < np; ++j){
            double f2 = psr->obsn[ip[j]].freq;
            double g2 = 1.0/(DM_CONST*f2*f2);
            mm[i][j] *= g1*g2;
        }
    }
}
