

namespace TK {
template<typename D>
void _SVD(matrix<D> &in, matrix<D> &U, diagonal<D> &S, matrix<D> &V){
    TKsingularValueDecomposition_lsq(in.getIdx(),in._rows,in._cols,V.getIdx(), S.getRaw(), U.getIdx());
}

template<typename D>
TK::SVD<D>::SVD(TK::matrix<D> &in) : in(in), U(in._rows,in._rows), V(in._cols,in._cols), S((in._rows>in._cols) ? in._cols : in._rows) {
    TK::_SVD(this->in, U,S,V);
}


template<typename D>
TK::SVD<D>::SVD(const TK::matrix<D> &in) : in(in), U(in._rows,in._rows), V(in._cols,in._cols), S((in._rows>in._cols) ? in._cols : in._rows) { 
    TK::_SVD(this->in, U,S,V);
}

template class SVD<float>;
template class SVD<double>;
template class SVD<longdouble>;
}
