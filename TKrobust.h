
double TKrobust(double* data, double* white_data,
        double** designMatrix, double** white_designMatrix,
        double** constraintsMatrix, int ndata,int nparams, int nconstraints, double tol, char rescale_errors,
        double* outP, double* e, double** Ocvm, char robust);



