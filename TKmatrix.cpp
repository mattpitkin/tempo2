
template<typename DataType>
std::map<DataType*,TK::matrix<DataType>* > TK::matrix<DataType>::_map;

template <typename T>
void TK::matrix_register(TK::matrix<T> *m){
    //logdbg("register %lp",m->getRaw());
    m->_map[m->getRaw()] = m;
}

template <typename T>
void TK::matrix_deregister(TK::matrix<T> *m){
    //logdbg("deregister %lp",m->getRaw());
    m->_map.erase(m->getRaw());
}

template<typename D>
TK::matrix<D> *TK::matrix<D>::getFromPtr(D* p){
    TK::matrix<D> *ret = TK::matrix<D>::_map[p];
    assert(ret);
    return static_cast<TK::matrix<D>*>(ret);
}

template<typename D>
TK::vector<D> *TK::vector<D>::getFromPtr(D* p){
    TK::matrix<D> *ret = TK::matrix<D>::_map[p];
    assert(ret);
    return static_cast<TK::vector<D>*>(ret);
}


template<typename D>
TK::diagonal<D> *TK::diagonal<D>::getFromPtr(D* p){
    TK::matrix<D> *ret = TK::matrix<D>::_map[p];
    assert(ret);
    return static_cast<TK::diagonal<D>*>(ret);
}


template<typename DataType>
DataType** TK::matrix<DataType>::getIdx(){
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
TK::matrix<DataType> TK::matrix<DataType>::T(bool noswap) const{
    if (noswap){
        TK::matrix<DataType> ret(_cols,_rows,!_rowmajor);
        ret._raw += _raw; 
        return ret;
    } else {
        TK::matrix<DataType> ret(_cols,_rows,_rowmajor);
        const size_t fast = _rowmajor ? _cols : _rows;
        const size_t slow = _rowmajor ? _rows : _cols;
        for (size_t ifast = 0; ifast < fast; ifast++){
            ret._raw[std::slice(ifast*slow,slow,1)] += _raw[std::slice(ifast,slow,fast)];
        }
        return ret;
    }
}


template class TK::matrix<float>;
template class TK::matrix<double>;
template class TK::matrix<longdouble>;


// template specialisations to use BLAS multiply when avaliable.


#if defined(HAVE_BLAS) && !defined(DEFINED_DBL_MM)
#define DEFINED_DBL_MM

#define F77_dgemm F77_FUNC(dgemm,DGEMM)
extern "C" {
    extern void F77_dgemm(const char* ta, const char* tb, int* m, int* n, int* k, double* alpha, 
            double* a, int* lda, double* b, int* ldb, double* beta, double* c, int* ldc);
#define F77_dgemv F77_FUNC(dgemv,DGEMV)
    extern void F77_dgemv(const char* trans, int* m, int* n, double* alpha, 
            double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);
}

namespace TK {

template<>
TK::vector<double> operator*(const TK::matrix<double> &lhs, const TK::vector<double> &rhs){
    assert(lhs._cols == rhs._rows);
    double alpha = 1;
    double beta = 0;
    int m = lhs._rowmajor ? lhs._cols : lhs._rows;
    int n = lhs._rowmajor ? lhs._rows : lhs._cols;
    int lda = m;
    int ldb = 1;

    TK::vector<double> out(lhs._rows);
    int ldc = 1;

    logdbg("dgemv");

    F77_dgemv(lhs._rowmajor ? "T":"N", // transpose A, B
            &m,&n,
            &alpha, const_cast<double*>(lhs.getRaw()), &lda,
            const_cast<double*>(rhs.getRaw()), &ldb, &beta,
            out.getRaw(), &ldc);

    return out;
}

// specialisation of multiply for doubles - use BLAS if avaliable
template<>
TK::matrix<double> operator*(const TK::matrix<double> &lhs, const TK::matrix<double> &rhs){
    if(lhs.triangular()=='D' || rhs.triangular()=='D'){
        logdbg("fall");
        return TK::fallbackMultiply(lhs,rhs);
    }
    logdbg("BLAS");
    double alpha = 1;
    double beta = 0;
    int m = lhs._rows;
    int n = rhs._cols;
    int k = lhs._cols;
    int lda = lhs._rowmajor ? k : m;
    int ldb = rhs._rowmajor ? n : k;
    int ldc = m;

    TK::matrix<double> out(m,n,false); // col_major

    logdbg("dgemm");
    F77_dgemm(lhs._rowmajor ? "T":"N", rhs._rowmajor ? "T":"N", // transpose A, B
            &m,&n,&k,
            &alpha, const_cast<double*>(lhs.getRaw()), &lda,
            const_cast<double*>(rhs.getRaw()), &ldb, &beta,
            out.getRaw(), &ldc);

    return out;
}
} // namespace TK

#endif // BLAS


#ifdef HAVE_MPACK
#ifndef DEFINED_DBL_MM
#define DEFINED_DBL_MM

#include <mblas_double.h>
namespace TK {

template<>
TK::vector<double> operator*(const TK::matrix<double> &lhs, const TK::vector<double> &rhs){
    assert(lhs._cols == rhs._rows);
    double alpha = 1;
    double beta = 0;
    int m = lhs._rowmajor ? lhs._cols : lhs._rows;
    int n = lhs._rowmajor ? lhs._rows : lhs._cols;
    int lda = m;
    int ldb = 1;

    TK::vector<double> out(lhs._rows);
    int ldc = 1;


    Rgemv(lhs._rowmajor ? "T":"N", // transpose A, B
            m,n,
            alpha, const_cast<double*>(lhs.getRaw()), lda,
            const_cast<double*>(rhs.getRaw()), ldb, beta,
            out.getRaw(), ldc);

    return out;
}

template<>
TK::matrix<double> operator*(const TK::matrix<double> &lhs, const TK::matrix<double> &rhs){
    if(lhs.triangular()=='D' || rhs.triangular()=='D'){
        logdbg("fall");
        return TK::fallbackMultiply(lhs,rhs);
    }
    double alpha = 1;
    double beta = 0;
    int m = lhs._rows;
    int n = rhs._cols;
    int k = lhs._cols;
    int lda = lhs._rowmajor ? k : m;
    int ldb = rhs._rowmajor ? n : k;
    int ldc = m;

    TK::matrix<double> out(lhs._cols,rhs._rows,false); // col_major

    logdbg("Rgemm [%d,%d] x [%d,%d] => [%d,%d]",lhs._rows,lhs._cols, rhs._rows, rhs._cols, out._rows,out._cols);
    Rgemm(lhs._rowmajor ? "T":"N", rhs._rowmajor ? "T":"N", // transpose A, B
            m,n,k,
            alpha, const_cast<double*>(lhs.getRaw()), lda,
            const_cast<double*>(rhs.getRaw()), ldb, beta,
            out.getRaw(), ldc);

    return out;
}
} // namespace TK

#endif // ifndef DEFINED_DBL_MM

#if defined(LONGDOUBLE_IS_FLOAT128) && ! defined(DEFINED_LD_MM)
#define DEFINED_LD_MM
#include <mblas___float128.h>
namespace TK {
// specialisation of multiply for float128 - use MBLAS if avaliable
template<>
TK::matrix<longdouble> operator*(const TK::matrix<longdouble> &lhs, const TK::matrix<longdouble> &rhs){
    if(lhs.triangular()=='D' || rhs.triangular()=='D'){
        logdbg("fall");
        return TK::fallbackMultiply(lhs,rhs);
    }
    longdouble alpha = 1;
    longdouble beta = 0;
    int m = lhs._rows;
    int n = rhs._cols;
    int k = lhs._cols;
    int lda = lhs._rowmajor ? lhs._cols : lhs._rows;
    int ldb = rhs._rowmajor ? rhs._cols : rhs._rows;

    TK::matrix<longdouble> out(lhs._cols,rhs._rows,false); // col_major
    int ldc = out._cols;

    logdbg("Rgemm");
    Rgemm(lhs._rowmajor ? "T":"N", rhs._rowmajor ? "T":"N", // transpose A, B
            m,n,k,
            alpha, const_cast<longdouble*>(lhs.getRaw()), lda,
            const_cast<longdouble*>(rhs.getRaw()), ldb, beta,
            out.getRaw(), ldc);

    return out;
}

template<>
TK::vector<longdouble> operator*(const TK::matrix<longdouble> &lhs, const TK::vector<longdouble> &rhs){
    assert(lhs._cols == rhs._rows);
    longdouble alpha = 1;
    longdouble beta = 0;
    int m = lhs._rowmajor ? lhs._cols : lhs._rows;
    int n = lhs._rowmajor ? lhs._rows : lhs._cols;
    int lda = m;
    int ldb = 1;

    TK::vector<longdouble> out(lhs._rows);
    int ldc = 1;

    logdbg("Rgemv");

    Rgemv(lhs._rowmajor ? "T":"N", // transpose A, B
            m,n,
            alpha, const_cast<longdouble*>(lhs.getRaw()), lda,
            const_cast<longdouble*>(rhs.getRaw()), ldb, beta,
            out.getRaw(), ldc);

    return out;
}

} // namespace TK
#endif // LONGDOUBLE_IS_FLOAT128

#endif // HAVE_MPACK


namespace TK {
#ifndef DEFINED_DBL_MM
template<>
TK::vector<double> operator*(const TK::matrix<double> &lhs, const TK::vector<double> &rhs){
    return TK::fallbackMultiply(lhs,rhs);
}

template<>
TK::matrix<double> operator*(const TK::matrix<double> &lhs, const TK::matrix<double> &rhs){
        return TK::fallbackMultiply(lhs,rhs);
}
#endif

#ifndef DEFINED_LD_MM
template<>
TK::vector<longdouble> operator*(const TK::matrix<longdouble> &lhs, const TK::vector<longdouble> &rhs){
    return TK::fallbackMultiply(lhs,rhs);
}

template<>
TK::matrix<longdouble> operator*(const TK::matrix<longdouble> &lhs, const TK::matrix<longdouble> &rhs){
        return TK::fallbackMultiply(lhs,rhs);
}
#endif

} // namespace TK
