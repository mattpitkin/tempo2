//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards
#ifndef __TKmatrix_h
#define __TKmatrix_h
#include "TKlongdouble.h"
#include "TKlog.h"

/*
 *    This file is part of TEMPO2. 
 * 
 *    TEMPO2 is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    TEMPO2 is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
 */

/*
 *    If you use TEMPO2 then please acknowledge it by citing 
 *    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
 *    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
 *    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
 *    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
 *    timing model.
 */
#ifdef __cplusplus
extern "C" {
#endif
    // THE C INTERFACE
    typedef double** TKmatrix_d;
    typedef double* TKvector_d;

    typedef float** TKmatrix_f;
    typedef float*  TKvector_f;

    typedef longdouble** TKmatrix_L;
    typedef longdouble*  TKvector_L;


#define TK_POSTFIX _d
#include "TKmatrix.detail.h"
#undef TK_POSTFIX

#define TK_POSTFIX _f
#include "TKmatrix.detail.h"
#undef TK_POSTFIX

#define TK_POSTFIX _L
#include "TKmatrix.detail.h"
#undef TK_POSTFIX

    // END OF C INTERFACE
#ifdef __cplusplus
}
// THE GLORIOUS C++ INTERFACE!
#include <valarray>
#include <vector>
#include <map>
#include <assert.h>

namespace TK {

// class declarations
template<typename DataType=double>
class vector;

template<typename DataType=double>
class diagonal;

template<typename DataType=double>
class matrix;


// function declarations
template<typename D, typename F>
TK::vector<D> fallbackMultiply(const TK::matrix<D> &lhs, const TK::vector<F> &rhs);

template<typename D, typename F>
TK::matrix<D> fallbackMultiply(const TK::matrix<D> &lhs, const TK::matrix<F> &rhs);

template<typename D>
void matrix_register(TK::matrix<D> *m);

template<typename D>
void matrix_deregister(TK::matrix<D> *m);

}
// ** End of function declarations

// Class definitions.
template<typename DataType>
class TK::matrix{
    class Row {
        private:
            TK::matrix<DataType> *_src;
            const size_t _row;
        public:
            Row(TK::matrix<DataType> *src,size_t row) : _src(src), _row(row) {
            }

            DataType& operator[](size_t col) {
                return _src->getRef(_row,col);
            }

            const DataType& operator[](size_t col) const {
                return _src->getRef(_row,col);
            }

    };
    private:
    matrix (size_t rows) : 
        _cols(rows), _rows(rows), _rowmajor(false),
        _slice(0,rows,1),
        _raw((DataType)0, rows), _idx(1,NULL){
            TK::matrix_register(this);
            _triangular='D';
        }

    public:
    matrix (size_t rows, size_t cols, bool rowmajor=true) : 
        _cols(cols), _rows(rows), _rowmajor(rowmajor),
        _slice(0,rowmajor ? _rows : _cols,rowmajor ? _cols : _rows),
        _raw((DataType)0, rows*cols), _idx(rowmajor ? _rows : _cols,NULL){
            TK::matrix_register(this);
            _triangular='N';
        }

    virtual ~matrix() {
        TK::matrix_deregister(this);
    }

    DataType** getIdx();

    TK::matrix<DataType> T(bool noswap=false) const;

    template<typename F>
TK::matrix<DataType> &operator*=(const TK::matrix<F> rhs);

    template<typename F>
TK::matrix<DataType> &operator+=(const TK::matrix<F> &rhs);

    DataType *packTri() const {
        assert(_triangular=='L' || _triangular=='U');
        DataType *_t = static_cast<DataType*>(malloc(sizeof(DataType)*(_rows*(_rows+1))/2));

        if(_triangular=='U'){
            for (size_t j=0;j < _cols;j++){ // cols
                for (size_t i=0; i <=j; i++) { // rows
                    // i=0 j=1 => row0, col1
                    int fi = i+1;
                    int fj = j+1;
                    _t[fi + (fj-1)*fj/2 - 1] = this->get(i,j);
                }
            }
        } else {
            int jc=0;
            for (size_t j=0;j<_cols;j++){ // cols
                for (size_t i=j;i<_rows;i++){ //rows
                    _t[jc+i-j] = this->get(i,j);
                }
                jc=jc+_rows-j;
            }
        }
        return _t;
    }

    void unPackTri(const char UPLO, const DataType* _t) {
        assert(UPLO=='L' || UPLO=='U');
        this->triangular(UPLO);

        if(_triangular=='U'){
            for (size_t j=0;j < _cols;j++){ // cols
                for (size_t i=0; i <=j; i++) { // rows
                    // i=0 j=1 => row0, col1
                    int fi = i+1;
                    int fj = j+1;
                    this->set(i,j,_t[fi + (fj-1)*fj/2 - 1]);
                }
            }
        } else {
            int jc=0;
            for (size_t j=0;j<_cols;j++){ // cols
                for (size_t i=j;i<_rows;i++){ //rows
                    this->set(i,j,_t[jc+i-j]);
                }
                jc=jc+_rows-j;
            }
        }
    }

    TK::matrix<DataType> inv() const;


    Row operator[](size_t row){
        assert(row < _rows);
        return Row(this,row);
    }

    const Row operator[](size_t row) const{
        assert(row < _rows);
        return Row(const_cast<TK::matrix<DataType>*>(this),row);
    }


    virtual DataType &getRef(size_t row, size_t col){
        assert(row < _rows);
        assert(col< _cols);
        return _raw[_rowmajor ? (row*_cols + col) : (col*_rows + row)];
    }
    virtual const DataType &getRef(size_t row, size_t col) const {
        assert(row < _rows);
        assert(col< _cols);
        return _raw[_rowmajor ? (row*_cols + col) : (col*_rows + row)];
    }

    char triangular() const{
        return _triangular;
    }

    void triangular(const char UPLO){
        assert(UPLO=='L' || UPLO=='N' || UPLO=='U');
        assert(UPLO=='N' || (_rows==_cols));
        _triangular = UPLO;
    }

    template<typename F> void is(const TK::matrix<F> &rhs){
        assert(rhs._rows == _rows);
        assert(rhs._cols == _cols);
        if(rhs._rowmajor == _rowmajor){
            this->_raw = rhs._raw;
        } else {
            this->_raw = rhs.T()._raw;
        }
    }

    inline const DataType* getRaw() const {
        return &_raw[0];
    }

    inline DataType* getRaw() {
        return &_raw[0];
    }

    inline void set(size_t row,size_t col,DataType d){
        this->getRef(row,col) = d;
    }

    inline const DataType get(size_t row,size_t col) const{
        return this->getRef(row,col);
    }

    template<typename D, typename F>
friend TK::matrix<D> operator+(const TK::matrix<D> &lhs, const TK::matrix<F> &rhs);
    template<typename D, typename F>
friend TK::matrix<D> operator*(const TK::matrix<D> &lhs, const TK::matrix<F> &rhs);

    template<typename D, typename F>
friend TK::diagonal<D> operator*(const TK::diagonal<D> &lhs, const TK::diagonal<F> &rhs);

    static TK::matrix<DataType> *getFromPtr(DataType* p);

    public:
    const size_t _cols;
    const size_t _rows;
    const bool _rowmajor;
    static std::map<DataType*,TK::matrix<DataType>* > _map;
    private:
    char _triangular;
    std::slice _slice;
    std::valarray<DataType> _raw;
    std::vector<DataType*> _idx;

    friend class TK::vector<DataType>;
    friend class TK::diagonal<DataType>;

};

template<typename DataType>
class TK::vector : public TK::matrix<DataType> {
    public:
        vector(size_t n,bool colvector=false) : TK::matrix<DataType>(colvector ? 1:n, colvector ? n:1, colvector)
            , _colvector(colvector) {
            }

        DataType &operator[](size_t i) {
            return this->getRef(_colvector ? 0:i, _colvector ? i:0);
        }

        const DataType &operator[](size_t i) const {
            return this->getRef(_colvector ? 0:i, _colvector ? i:0);
        }


        static TK::vector<DataType> *getFromPtr(DataType* p);

    private:
        const bool _colvector;
};

template<typename DataType>
class TK::diagonal : public TK::matrix<DataType> {

    public:
        diagonal(size_t n) : TK::matrix<DataType>(n), czero(0) {
        }

        DataType &getRef(size_t row, size_t col) {
            assert(row < this->_rows);
            assert(col < this->_cols);
            if (row==col){
                return this->_raw[row];
            } else {
                return zero;
            }
        }

        const DataType &getRef(size_t row, size_t col) const {
            assert(row < this->_rows);
            assert(col < this->_cols);
            if (row==col){
                return this->_raw[row];
            } else {
                return czero;
            }
        }

        TK::diagonal<DataType> inv() const{
            TK::diagonal<DataType> ret(this->_rows);
            ret->_raw = 1.0 / this->_raw;
        }

        static TK::diagonal<DataType> *getFromPtr(DataType* p);
    private:
        const DataType czero;
        DataType zero;
};

#include "TKmatrix.hpp"

// end of C++ 
#endif

// end of TK::matrix.h
#endif
