#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <map>
#include "TKmatrix.h"
#include "TKlog.h"
#include "T2accel.h"


template<typename DataType>
std::map<DataType*,TKmatrix<DataType>* > TKmatrix<DataType>::_map;

template <typename T>
void TKmatrix_register(TKmatrix<T> *m){
    //logmsg("register %lp",m->getRaw());
    m->_map[m->getRaw()] = m;
}

template <typename T>
void TKmatrix_deregister(TKmatrix<T> *m){
    //logmsg("deregister %lp",m->getRaw());
    m->_map.erase(m->getRaw());
}

template <typename T>
TKmatrix<T>* TKmatrix_getFromPtr(T* p){
    return TKmatrix<T>::_map[p];
}

template<typename DataType>
DataType** TKmatrix<DataType>::getIdx(){
    if (_idx[0]==NULL){
        const size_t fast = _rowmajor ? _cols : _rows;
        const size_t slow = _rowmajor ? _rows : _cols;
        for (size_t islow = 0; islow < slow; islow++){
            _idx[islow] = &_raw[islow*fast];
        }
    }
    return &_idx[0];

}


template <typename DataType>
TKmatrix<DataType> *TKmatrix<DataType>::T(bool noswap){
    TKmatrix<DataType> *ret;
    if (noswap){
        ret = new TKmatrix<DataType>(_cols,_rows,!_rowmajor);
        ret->_raw += _raw; 
    } else {
        ret = new TKmatrix<DataType>(_cols,_rows,_rowmajor);
        const size_t fast = _rowmajor ? _cols : _rows;
        const size_t slow = _rowmajor ? _rows : _cols;
        for (size_t ifast = 0; ifast < fast; ifast++){
            ret->_raw[std::slice(ifast*slow,slow,1)] += _raw[std::slice(ifast,slow,fast)];
        }
    }
    return ret;
}



template void TKmatrix_register(TKmatrix<float>*);
template void TKmatrix_register(TKmatrix<double>*);
template void TKmatrix_register(TKmatrix<longdouble>*);

template void TKmatrix_deregister(TKmatrix<float>*);
template void TKmatrix_deregister(TKmatrix<double>*);
template void TKmatrix_deregister(TKmatrix<longdouble>*);

template TKmatrix<float>* TKmatrix_getFromPtr(float*);
template TKmatrix<double>* TKmatrix_getFromPtr(double*);
template TKmatrix<longdouble>* TKmatrix_getFromPtr(longdouble*);

template class TKmatrix<float>;
template class TKmatrix<double>;
template class TKmatrix<longdouble>;



/// C functions


TKmatrix_d malloc_matrix_sq_d(size_t rc){
    return malloc_matrix_d(rc,rc);
}

TKmatrix_d malloc_matrix_d(size_t r, size_t c){
    return static_cast<TKmatrix_d>((new TKmatrix<double>(r,c))->getIdx());
}

TKvector_d malloc_vector_d(size_t r) {
    return static_cast<TKvector_d>((new TKvector<double>(r))->getRaw());
}

void free_vector_d(TKvector_d A){
    free_matrix_d(&A);
}
void free_matrix_d(TKmatrix_d A){
    TKmatrix<double> *p =TKmatrix_getFromPtr(A[0]);
    if (p){
        delete p;
    } else {
        logerr("Tried to free non-existant TKmatrix");
    }
}


