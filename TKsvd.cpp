
namespace TK {

template<typename D>
void _SVD(matrix<D> &in, matrix<D> &U, diagonal<D> &S, matrix<D> &V){
    assert(in._rows >= in._cols);
    TKsingularValueDecomposition_lsq(in.getIdx(),in._rows,in._cols,V.getIdx(), S.getRaw(), U.getIdx());
    U.is(in);
}



template<typename D>
void _SVD(matrix<D> &in, matrix<D> &U, diagonal<D> &S, matrix<D> &V);

template<typename D>
TK::SVD<D>::SVD(TK::matrix<D> &in) : in(in,true), U(in._rows,in._cols), V(in._cols,in._cols), S(in._cols) {
    TK::_SVD(this->in, U,S,V);
}


template<typename D>
TK::SVD<D>::SVD(const TK::matrix<D> &in) : in(in,true), U(in._rows,in._cols), V(in._cols,in._cols), S(in._cols) { 
    TK::_SVD(this->in, U,S,V);
}
} 



#ifdef HAVE_MPACK

#define DEFINED_SVD_d
#define DEFINED_SVD_ld
#include <mblas_double.h>
#include <mblas___float128.h>
#include <mlapack_double.h>
#include <mlapack___float128.h>
#include <algorithm>

namespace TK {

template<>
void _SVD(matrix<double> &in, matrix<double> &U, diagonal<double> &S, matrix<double> &V){
    assert(!in._rowmajor);
    assert(!U._rowmajor);
    assert(V._rowmajor);

    mpackint M = in._rows;
    mpackint N = in._cols;
    assert(M >= N );
    mpackint info;
    mpackint iwork[8*M];
    mpackint nwork = 3*std::min(M,N) + std::max(std::max(M,N),5*std::min(M,N)*std::min(M,N)+4*std::min(M,N));
    nwork*=4;
    double work[nwork];
    logmsg("Rgesdd");

    Rgesdd("O", M, N, in.getRaw(), M, S.getRaw(), U.getRaw(), M, V.getRaw(), N, work, nwork, &iwork[0], &info );
    logmsg("info=%d",info);

    U.is(in);
    logmsg("BBBB");
    logmsg("BBBB");

}


template<>
TK::SVD<double>::SVD(TK::matrix<double> &in) : in(in,false), U(in._rows,in._cols,false), V(in._cols,in._cols,true), S(in._cols) {
    TK::_SVD(this->in, U,S,V);
}


template<>
TK::SVD<double>::SVD(const TK::matrix<double> &in) : in(in,false), U(in._rows,in._cols,false), V(in._cols,in._cols,true), S(in._cols) {
    TK::_SVD(this->in, U,S,V);
}


}

#endif


namespace TK{

template class SVD<float>;
template class SVD<double>;
template class SVD<longdouble>;

}
