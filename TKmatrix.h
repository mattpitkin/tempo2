//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards
#ifndef __TKmatrix_h
#define __TKmatrix_h
#include "TKlongdouble.h"

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

    void TKmultMatrix_sq_d   (TKmatrix_d A, TKmatrix_d B, size_t rcA, size_t rowsB, TKmatrix_d O);
    void TKmultMatrixVec_sq_d(TKmatrix_d A, TKvector_d b, size_t rcA, TKvector_d o);
    void TKmultMatrix_d      (TKmatrix_d A, TKmatrix_d B, size_t rowsA, size_t colsA, size_t rowsB, TKmatrix_d O);
    void TKmultMatrixVec_d   (TKmatrix_d A, TKvector_d b, size_t rowsA, size_t colsA, TKvector_d o);

    TKmatrix_d malloc_matrix_sq_d(size_t rc);
    TKmatrix_d malloc_matrix_d(size_t r, size_t c);
    TKvector_d malloc_vector_d(size_t r);

    void free_matrix_d(TKmatrix_d A);
    void free_vector_d(TKvector_d A);

    size_t getRows_TKmatrix_d(TKmatrix_d matrix);
    size_t getCols_TKmatrix_d(TKmatrix_d matrix);
    size_t getRows_TKvector_d(TKvector_d matrix);

    typedef float** TKmatrix_f;
    typedef float*  TKvector_f;

    typedef longdouble** TKmatrix_L;
    typedef longdouble*  TKvector_L;


#ifdef __cplusplus
}
// THE C++ INTERFACE!
#include <valarray>
#include <vector>
#include <map>

template<typename DataType=double>
class TKmatrix {
    public:
        TKmatrix (size_t rows, size_t cols, bool rowmajor=true) : 
            _cols(cols), _rows(rows), _rowmajor(rowmajor),
            _slice(0,rowmajor ? _rows : _cols,rowmajor ? _cols : _rows),
            _raw((DataType)0, rows*cols), _idx(rowmajor ? _rows : _cols,NULL){
                TKmatrix_register(this);
            }

        ~TKmatrix() {
            TKmatrix_deregister(this);
        }

        DataType** getIdx();

        TKmatrix<DataType> *T(bool noswap=false);

        inline const DataType* getRaw() const {
            return &_raw[0];
        }

        inline DataType* getRaw() {
            return &_raw[0];
        }

        inline void set(size_t row,size_t col,DataType d){
            assert(row < _rows);
            assert(col< _cols);
            _raw[_rowmajor ? (row*_cols + col) : (col*_rows + row)] = d;
        }

        inline DataType get(size_t row,size_t col){
            assert(row < _rows);
            assert(col< _cols);
            return _raw[_rowmajor ? (row*_cols + col) : (col*_rows + row)];
        }

    public:
        const size_t _cols;
        const size_t _rows;
        const bool _rowmajor;
        static std::map<DataType*,TKmatrix<DataType>* > _map;
    private:
        std::slice _slice;
        std::valarray<DataType> _raw;
        std::vector<DataType*> _idx;

};

template<typename DataType=double>
class TKvector : public TKmatrix<DataType> {
    public:
        TKvector(size_t n) : TKmatrix<DataType>(n,1,false){
    }
};


template <typename T>
void TKmatrix_register(TKmatrix<T> *m);

template <typename T>
void TKmatrix_deregister(TKmatrix<T> *m);

template <typename T>
TKmatrix<T>* TKmatrix_getFromPtr(T* p);


// end of C++ 
#endif

// end of TKmatrix.h
#endif
