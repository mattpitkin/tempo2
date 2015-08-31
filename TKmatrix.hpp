
// Template method definitions
template<typename D, typename F>
TK::vector<D> TK::fallbackMultiply(const TK::matrix<D> &lhs, const TK::vector<F> &rhs){
    assert(lhs._cols == rhs._rows);
    TK::vector<D> ret(lhs._rows,false);
    size_t i,k;
    for (i=0;i<lhs._rows;i++){
        for (k=0;k<lhs._cols;k++) {
            assert(k < rhs._rows);
            assert(i < ret._rows);
            ret[i] += lhs[i][k]*rhs[k];
        }
    }
    return ret;
}
template<typename D, typename F>
TK::matrix<D> TK::fallbackMultiply(const TK::matrix<D> &lhs, const TK::matrix<F> &rhs){
    assert(lhs._cols == rhs._rows);
    TK::matrix<D> ret(lhs._rows,rhs._cols,lhs._rowmajor);
    size_t i,j,k;
    for (k=0;k<lhs._cols;k++) {
        for (i=0;i<lhs._rows;i++){
            for (j=0;j<rhs._cols;j++){
                ret[i][j] += lhs[i][k]*rhs[k][j];
            }
        }
    }
    return ret;
}



namespace TK {
// operator definitions

template<typename D, typename F>
TK::diagonal<D> operator*(const TK::diagonal<D> &lhs, const TK::diagonal<F> &rhs){
    assert(lhs._cols == rhs._cols);
    TK::diagonal<D> ret(lhs._cols);
    ret._raw = lhs._raw;
    ret._raw *= rhs._raw;
    return ret;
}

template<typename D, typename F>
TK::matrix<D> operator*(const TK::matrix<D> &lhs, const TK::diagonal<F> &rhs){
    assert(lhs._cols == rhs._cols);
    TK::diagonal<D> ret(lhs._cols);
    ret._raw = lhs._raw;
    ret._raw *= rhs._raw;
    return ret;
}


template<typename D, typename F>
TK::matrix<D> operator*(const TK::matrix<D> &lhs, const TK::matrix<F> &rhs){
    return TK::fallbackMultiply(lhs,rhs);
}
template<typename D, typename F>
TK::vector<D> operator*(const TK::matrix<D> &lhs, const TK::vector<F> &rhs){
    return TK::fallbackMultiply(lhs,rhs);
}




template<typename D, typename F>
TK::matrix<D> operator+(const TK::matrix<D> &lhs, const TK::matrix<F> &rhs){
    assert(lhs._rows == rhs._rows);
    assert(lhs._cols == rhs._cols);
    TK::matrix<D> ret(lhs._rows,lhs._cols);
    if (lhs._rowmajor == rhs._rowmajor){
        ret._raw = lhs._raw + rhs._raw;
    } else {
        ret._raw = lhs._raw;
        const size_t fast = lhs._rowmajor ? lhs._cols : lhs._rows;
        const size_t slow = lhs._rowmajor ? lhs._rows : lhs._cols;
        for (size_t ifast = 0; ifast < fast; ifast++){
            ret._raw[std::slice(ifast,slow,fast)] +=
                rhs._raw[std::slice(ifast*slow,slow,1)];
        }

    }
    return ret;
}

}

