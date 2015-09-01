#include <ostream>

// Template method definitions
template<typename D, typename F>
TK::vector<D> TK::fallbackMultiply(const TK::matrix<D> &lhs, const TK::vector<F> &rhs){
    assert(lhs._cols == rhs._rows);
    logmsg("Fallback MV");
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
    logmsg("Fallback MM");
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
template<typename charT, typename traits, typename D>
std::basic_ostream<charT, traits> &
operator<< (std::basic_ostream<charT, traits> &lhs, matrix<D> const &rhs){
    for(size_t r = 0; r < rhs._rows; r++){
        for(size_t c = 0; c < rhs._cols; c++){
            lhs << rhs[r][c]<< " ";
        }
        lhs << std::endl;
    }
    return lhs;
}





template<typename D, typename F>
TK::diagonal<D> operator*(const TK::diagonal<D> &lhs, const TK::diagonal<F> &rhs){
    assert(lhs._cols == rhs._rows);
    TK::diagonal<D> ret(lhs._rows,rhs._cols);
    ret._raw = lhs._raw;
    if (lhs._rows == rhs._cols){
    ret._raw *= rhs._raw;
    } else {
        const size_t smallest = ret._rows < ret._cols ? ret._rows : ret._cols;
        for(size_t i=0; i < smallest; i++){
            ret[i][i] = lhs[i][i]*rhs[i][i];
        }
    }
    return ret;
}

template<typename D, typename F>
TK::matrix<D> operator*(const TK::matrix<D> &lhs, const TK::matrix<F> &rhs){
    logmsg("tmp mm*");
    return TK::fallbackMultiply(lhs,rhs);
}
template<typename D, typename F>
TK::vector<D> operator*(const TK::matrix<D> &lhs, const TK::vector<F> &rhs){
    logmsg("tmp mv*");
    return TK::fallbackMultiply(lhs,rhs);
}


// template specialisations
// double M*V
template<>
TK::vector<double> operator*(const TK::matrix<double> &lhs, const TK::vector<double> &rhs);

// double M*M
template<>
TK::matrix<double> operator*(const TK::matrix<double> &lhs, const TK::matrix<double> &rhs);

// ld M*V
template<>
TK::vector<longdouble> operator*(const TK::matrix<longdouble> &lhs, const TK::vector<longdouble> &rhs);

// ld M*M
template<>
TK::matrix<longdouble> operator*(const TK::matrix<longdouble> &lhs, const TK::matrix<longdouble> &rhs);

// end of specialisations

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

