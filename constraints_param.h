struct constraint_param_info {
    int param;
    int param_k;
    double val;
    double err;
};

double constraint_param_function(pulsar *psr,int ipsr, int iconstraint,int iparam,int constraintk,int k,void* special);

