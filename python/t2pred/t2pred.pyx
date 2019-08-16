from libc.stdlib cimport malloc, free
cimport numpy as np
import numpy as np


cdef extern from "tempo2pred.h":
    ctypedef struct T2Predictor:
        pass


    ctypedef enum T2PredictorKind:
        NonePredType, Cheby, T1

    void T2Predictor_Init(T2Predictor *t2p)
    int T2Predictor_Read(T2Predictor *t2p, char *fname)

    char * T2Predictor_GetPSRName(T2Predictor *t2p)
    char * T2Predictor_GetSiteName(T2Predictor *t2p)
    long double T2Predictor_GetStartMJD(T2Predictor *t2p)
    long double T2Predictor_GetEndMJD(T2Predictor *t2p)
    long double T2Predictor_GetStartFreq(T2Predictor *t2p)
    long double T2Predictor_GetEndFreq(T2Predictor *t2p)
    T2PredictorKind T2Predictor_Kind(T2Predictor *t2p)

    long double T2Predictor_GetPhase(const T2Predictor *t2p, long double mjd, long double freq)
    long double T2Predictor_GetFrequency(const T2Predictor *t2p, long double mjd, long double freq)

cdef extern from "t2pred_loops.h":
    void T2Predictor_GetPhase_array_ld(const T2Predictor *t2p, long double *mjd, int nmjd, long double freq, double* out)
    void T2Predictor_GetFrequency_array_ld(const T2Predictor *t2p, long double *mjd, int nmjd, long double freq, double* out)

    void T2Predictor_GetPhase_array(const T2Predictor *t2p, double refmjd, double *mjd, const int nmjd, double freq, double* out)
    void T2Predictor_GetFrequency_array(const T2Predictor *t2p, double refmjd, double *mjd, const int nmjd, double freq, double* out)

cdef class phase_predictor:
    cdef T2Predictor *thisptr
    cdef int count
    def __init__(self,filename=None):
        self.count=0
        self.thisptr = <T2Predictor*> malloc(sizeof(T2Predictor));
        T2Predictor_Init(self.thisptr)
        if not filename is None:
            self.read(filename)

    def __dealloc___(self):
        free(self.thisptr)


    def read(self,filename):
        bytestring = filename.encode('UTF-8')
        cdef char* c_string = bytestring
        return T2Predictor_Read(self.thisptr,c_string)

    def getPhase_array(self, mjd,freq):
        cdef int n = np.product(mjd.shape)
        phase=np.zeros(n)
        times=mjd.flatten(order='C')
        cdef long double[::1] times_view = times
        cdef double[::1] phase_view = phase
        T2Predictor_GetPhase_array_ld(self.thisptr,&times_view[0],n,freq,&phase_view[0])
        
        return phase.reshape(mjd.shape)

    def getFrequency_array(self, mjd,freq):
        cdef int n = np.product(mjd.shape)
        out=np.zeros(n)
        times=mjd.flatten(order='C')
        cdef long double[::1] times_view = times
        cdef double[::1] out_view = out
        T2Predictor_GetFrequency_array_ld(self.thisptr,&times_view[0],n,freq,&out_view[0])
        return out.reshape(mjd.shape)



    def getPhase(self,mjd,freq):
        return T2Predictor_GetPhase(self.thisptr,mjd,freq)

    def getFrequency(self,mjd,freq):
        return T2Predictor_GetFrequency(self.thisptr,mjd,freq)

    def kind(self):
        cdef T2PredictorKind kind
        kind = T2Predictor_Kind(self.thisptr)
        if kind == Cheby:
            return "Cheby"
        elif kind==T1:
            return "T1"
        else:
            return "NonePredType"


    def name(self):
        cdef char* name = T2Predictor_GetPSRName(self.thisptr)
        return name.decode('UTF-8')

    def site(self):
        cdef char* name = T2Predictor_GetSiteName(self.thisptr)
        return name.decode('UTF-8')

    def startMJD(self):
        return T2Predictor_GetStartMJD(self.thisptr)

    def endMJD(self):
        return T2Predictor_GetEndMJD(self.thisptr)

    def startFreq(self):
        return T2Predictor_GetStartFreq(self.thisptr)

    def endFreq(self):
        return T2Predictor_GetEndFreq(self.thisptr)

