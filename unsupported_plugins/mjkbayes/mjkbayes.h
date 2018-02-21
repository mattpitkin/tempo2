#include "tempo2.h"
#include <vector>
#include <string>

#define FITTYPE_PARAM 0
#define FITTYPE_CHOL 1
#define FITTYPE_EFAC 2
#define FITTYPE_EQUAD 3
#define FITTYPE_BIN 4

#define FITTYPE_CHOL_K_AMP 0
#define FITTYPE_CHOL_K_ALPHA 1
#define FITTYPE_CHOL_K_FC 2

#define FITTYPE_BIN_K_ASINI 11
#define FITTYPE_BIN_K_M2    12
#define FITTYPE_BIN_K_ECC   13
#define FITTYPE_BIN_K_PB    14
#define FITTYPE_BIN_K_INC  15
#define FITTYPE_BIN_K_M1  16


class mjkparam {
    public:
        mjkparam(int fittype, label fitlabel, int fitk,
                double centre, double halfrange) : 
            fittype(fittype),
               fitlabel(fitlabel),
               fitk(fitk),
               fitscale(halfrange),
               fitoffset(centre) {
               }
        mjkparam(int fittype, label fitlabel, int fitk) : 
            fittype(fittype),
               fitlabel(fitlabel),
               fitk(fitk),
               fitscale(0),
               fitoffset(0) {
               }

        std::string shortlabel(const pulsar* psr) const {
            switch(fittype) {
                case FITTYPE_PARAM:
                    return psr->param[fitlabel].shortlabel[fitk];
                case FITTYPE_BIN:
                    switch(fitk){
                        case FITTYPE_BIN_K_ASINI:
                            return "bin_ASINI";
                        case FITTYPE_BIN_K_M2:
                            return "bin_M2";
                        case FITTYPE_BIN_K_ECC:
                            return "bin_ECC";
                        case FITTYPE_BIN_K_PB:
                            return "bin_PB";
                        case FITTYPE_BIN_K_INC:
                            return "bin_inc";
                        case FITTYPE_BIN_K_M1:
                            return "bin_M1";
                    }
                    break;
                case FITTYPE_CHOL:
                    switch(fitk){
                        case FITTYPE_CHOL_K_ALPHA:
                            return "CHOL_alpha";
                        case FITTYPE_CHOL_K_AMP:
                            return "CHOL_log(amp)";
                        case FITTYPE_CHOL_K_FC:
                            return "CHOL_fc";
                    }

                    break;
                case FITTYPE_EFAC:
                    return "EFAC";
                case FITTYPE_EQUAD:
                    return "EQUAD";
            }
            return "????";
        }

    public:
        int fittype;
        label fitlabel;
        int fitk;
        double fitscale;
        double fitoffset;
        char* flagmask;
        std::string flagid;
        std::string flagval;
};

struct mjkcontext {
    pulsar* psr;
    FILE* debugfile;
    char root[1024];
    std::vector<mjkparam> params;
    std::vector<mjkparam> xtra;
};

void mjkbayes_analyse(pulsar* psr, struct mjkcontext *context);


void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &o, void *context);

void loadmjkbayescfg(const char* cfg, pulsar* psr, mjkcontext *context);

char* mjkbayesflagmask(pulsar* psr, const char* flag, const char* flagval);

