#include "tempo2.h"
#include <vector>
#include <string>

#define FITTYPE_PARAM 0
#define FITTYPE_CVM 1
#define FITTYPE_EFAC 2
#define FITTYPE_EQUAD 3
#define FITTYPE_BIN 4
#define FITTYPE_CONCVM 5

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
               fitoffset(centre),
               exp(false),fixed(false) {
               }
        mjkparam(int fittype, label fitlabel, int fitk) : 
            fittype(fittype),
               fitlabel(fitlabel),
               fitk(fitk),
               fitscale(0),
               fitoffset(0),
               exp(false),fixed(false) {
               }

        char* parseScaleoffset(char* string){
            while(string[0] != '\0' && isspace(string[0])){
                ++string;
            }
            if (string[0]=='\0') return 0;
            if (string[0]=='^') {
                exp=true;
                ++string;
            }
            if (string[0]=='\0') return 0;
            if (string[0]=='!') {
                // Fixed value
                fixed=true;
                ++string;
                int end=strlen(string);
                for (int i=0; i < end; ++i) {
                    if (string[i]=='!'){
                        string[i]='\0';
                        end=i+1;
                        break;
                    }
                }
                sscanf(string,"%lg",&fitoffset);
                fitscale=0;
                string += end;
                return string;

            }
            if (string[0]=='[') {
                ++string;
                // start and finish values
                int end=strlen(string);
                for (int i=0; i < end; ++i) {
                    if (string[i]==']'){
                        string[i]='\0';
                        end=i+1;
                        break;
                    }
                }
                double start,fin;
                sscanf(string,"%lg %lg",&start,&fin);
                fitscale = (fin-start)/2.0;
                fitoffset = (fin+start)/2.0;
                string += end;
                return string;
            }
            if (string[0]=='(') {
                ++string;
                // centre and half-range
                int end=strlen(string);
                for (int i=0; i < end; ++i) {
                    if (string[i]==')'){
                        end=i+1;
                        string[i]='\0';
                        break;
                    }
                }
                sscanf(string,"%lg %lg",&fitoffset,&fitscale);
                string += end;
                return string;
            }
            if (isdigit(string[0])) {
                // simple centre and half-range
                sscanf(string,"%lg %lg",&fitoffset,&fitscale);
                return 0;
            }

        }

        std::string fitdesc(const pulsar* psr) const {
            std::stringstream ss;
            if (exp) {
                ss << "log("<< shortlabel(psr) << ")";
            }else{
                ss << shortlabel(psr);
            }
            ss << " from " << fitoffset-fitscale << " to " << fitoffset+fitscale << "\n";
            return ss.str();

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
                case FITTYPE_CVM:
                    return txt;
                case FITTYPE_CONCVM:
                    return txt;
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
        bool exp;
        std::string flagid;
        std::string flagval;
        std::string txt;
        bool fixed;
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

