#define _cat(A,B) A ## B
#define cat(A,B) _cat(A,B)
#define _TKmatrix_function(fname) cat(fname, TK_POSTFIX)
#define _TKmatrix_type cat(TKmatrix, TK_POSTFIX )
#define _TKvector_type cat(TKvector,TK_POSTFIX)

void _TKmatrix_function(TKmultMatrix)       (_TKmatrix_type A, _TKmatrix_type B, _TKmatrix_type O);
void _TKmatrix_function(TKmultMatrixVec)    (_TKmatrix_type A, _TKvector_type b, _TKvector_type o);

_TKmatrix_type _TKmatrix_function(malloc_matrix_sq) (size_t rc);
_TKmatrix_type _TKmatrix_function(malloc_matrix)  (size_t r, size_t c);
_TKvector_type _TKmatrix_function(malloc_vector)    (size_t r);

void _TKmatrix_function(free_matrix) (_TKmatrix_type A);
void _TKmatrix_function(free_vector) (_TKvector_type A);

size_t _TKmatrix_function(getRows_TKmatrix) (_TKmatrix_type matrix);
size_t _TKmatrix_function(getCols_TKmatrix) (_TKmatrix_type matrix);
size_t _TKmatrix_function(getRows_TKvector) (_TKvector_type matrix);

#undef _TKmatrix_type
#undef _TKvector_type
#undef _TKmatrix_function
#undef _cat
#undef cat
