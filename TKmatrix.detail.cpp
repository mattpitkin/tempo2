#define _cat(A,B) A ## B
#define cat(A,B) _cat(A,B)
#define _TKmatrix_function(fname) cat(fname, TK_POSTFIX)
#define _TKmatrix_type cat(TKmatrix, TK_POSTFIX )
#define _TKvector_type cat(TKvector,TK_POSTFIX)

void _TKmatrix_function(TKmultMatrix)       (_TKmatrix_type A, _TKmatrix_type B, _TKmatrix_type O){
    TK::matrix<TKmatrix_D> *a = TK::matrix<TKmatrix_D>::getFromPtr(*A);
    TK::matrix<TKmatrix_D> *b = TK::matrix<TKmatrix_D>::getFromPtr(*B);
    TK::matrix<TKmatrix_D> *o = TK::matrix<TKmatrix_D>::getFromPtr(*O);
    assert(a);
    assert(b);
    assert(o);

    o->is((*a)*(*b));
}

void _TKmatrix_function(TKmultMatrixVec)    (_TKmatrix_type A, _TKvector_type b, _TKvector_type o){
    TK::matrix<TKmatrix_D> *m = TK::matrix<TKmatrix_D>::getFromPtr(*A);
    TK::vector<TKmatrix_D> *u = TK::vector<TKmatrix_D>::getFromPtr(b);
    TK::vector<TKmatrix_D> *v = TK::vector<TKmatrix_D>::getFromPtr(o);
    assert(m);
    assert(u);
    assert(v);

    TK::vector<TKmatrix_D> v2 = *m * (*u);
    v->is(v2);
}

_TKmatrix_type _TKmatrix_function(malloc_matrix_sq) (size_t rc){
    return _TKmatrix_function(malloc_matrix)(rc,rc);
}
_TKmatrix_type _TKmatrix_function(malloc_matrix)  (size_t r, size_t c){
    return static_cast<_TKmatrix_type>((new TK::matrix<TKmatrix_D>(r,c))->getIdx());
}

_TKvector_type _TKmatrix_function(malloc_vector)    (size_t r){
    return static_cast<_TKvector_type>((new TK::vector<TKmatrix_D>(r))->getRaw());
}

void _TKmatrix_function(free_matrix) (_TKmatrix_type A){
    TK::matrix<TKmatrix_D> *p =TK::matrix<TKmatrix_D>::getFromPtr(A[0]);
    if (p){
        delete p;
    } else {
        logerr("Tried to free non-existant TKmatrix");
    }
}


void _TKmatrix_function(free_vector) (_TKvector_type A){
    _TKmatrix_function(free_matrix) (&A);
}


size_t _TKmatrix_function(getRows_TKmatrix) (_TKmatrix_type matrix);
size_t _TKmatrix_function(getCols_TKmatrix) (_TKmatrix_type matrix);
size_t _TKmatrix_function(getRows_TKvector) (_TKvector_type matrix);

#undef _TKmatrix_type
#undef _TKvector_type
#undef _TKmatrix_function
#undef _cat
#undef cat
