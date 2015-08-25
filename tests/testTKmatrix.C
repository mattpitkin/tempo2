#include <gtest/gtest.h>
#include <src/gtest_main.cc>

#include"TKmatrix.h"


TEST(testTKmatrix_d, construct){
    TKmatrix<double> M(1,2);
    ASSERT_EQ(M._rows,1);
    ASSERT_EQ(M._cols,2);
}

TEST(testTKmatrix_d, getset){
    TKmatrix<double> M(2,2,true);
    M.set(0,0,1.0);
    M.set(0,1,-2.0);
    M.set(1,0,2.0);

    ASSERT_EQ(M.get(0,0),1.0);
    ASSERT_EQ(M.get(0,1),-2.0);
    ASSERT_EQ(M.get(1,0),2.0);
}

TEST(testTKmatrix_d, transpose){
    TKmatrix<double> M(3,3,true);
    M.set(0,0,1.0);
    M.set(0,1,-2.0);
    M.set(1,0,2.0);

    M.set(1,2,-3.0);
    M.set(2,1,3.0);

    M.set(0,2,-4.0);
    M.set(2,0,4.0);

    {
    TKmatrix<double> *MT = M.T(false);

    ASSERT_EQ(MT->get(0,0), 1.0);
    ASSERT_EQ(MT->get(1,0),-2.0);
    ASSERT_EQ(MT->get(0,1), 2.0);
    ASSERT_EQ(MT->get(2,1),-3.0);
    ASSERT_EQ(MT->get(1,2), 3.0);
    ASSERT_EQ(MT->get(2,0),-4.0);
    ASSERT_EQ(MT->get(0,2), 4.0);

    // check that the raw data was transposed
    for(size_t j =0 ; j < 3; j++){
        for(size_t k =0 ; k < 3; k++){
            size_t i1 = j+3*k;
            size_t i2 = k+3*j;
            EXPECT_EQ(MT->getRaw()[i1],M.getRaw()[i2]);
        }
    }
    delete MT;
    }

    {
        TKmatrix<double> *MT = M.T(true);

        ASSERT_EQ(MT->get(0,0), 1.0);
        ASSERT_EQ(MT->get(1,0),-2.0);
        ASSERT_EQ(MT->get(0,1), 2.0);
        ASSERT_EQ(MT->get(2,1),-3.0);
        ASSERT_EQ(MT->get(1,2), 3.0);
        ASSERT_EQ(MT->get(2,0),-4.0);
        ASSERT_EQ(MT->get(0,2), 4.0);

        for(size_t i =0 ; i < 9; i++){
            EXPECT_EQ(MT->getRaw()[i],M.getRaw()[i]);
        }
        delete MT;
    }

}

TEST(testTKmatrix_d, transposeBig){
    TKmatrix<double> M(1024,2048,true);
    for (int i=0;i<1024*2048;i++){
        M.getRaw()[i]=(double)i;
    }
    TKmatrix<double> *MT = M.T(false);
    ASSERT_EQ(MT->_rows,M._cols);
    ASSERT_EQ(MT->_cols,M._rows);

    for (int i=0;i<1024;i++){
        for (int j=0;j<2048;j++){
            ASSERT_EQ(M.get(i,j),MT->get(j,i));
        }
    }
    delete MT;
}

TEST(testTKmatrix_d, transposeBig2){
    TKmatrix<double> M(1024,2048,true);
    for (int i=0;i<1024*2048;i++){
        M.getRaw()[i]=(double)i;
    }
    TKmatrix<double> *MT = M.T(true);
    ASSERT_EQ(MT->_rows,M._cols);
    ASSERT_EQ(MT->_cols,M._rows);

    for (int i=0;i<1024;i++){
        for (int j=0;j<2048;j++){
            ASSERT_EQ(M.get(i,j),MT->get(j,i));
        }
    }
    delete MT;
}

TEST(testTKmatrix_d, getIdx){
    TKmatrix<double> M(2,2,true);
    M.set(0,0,1.0);
    M.set(0,1,-2.0);
    M.set(1,0,2.0);
    TKmatrix_d m = M.getIdx();

    ASSERT_EQ(m[0][0],1.0);
    ASSERT_EQ(m[0][1],-2.0);
    ASSERT_EQ(m[1][0],2.0);

}

TEST(testTKmatrix_d, toFromPtr){
    TKmatrix<double> M(2,2,true);
    double* p = M.getRaw();

    TKmatrix<double>* Mp = TKmatrix_getFromPtr(p);
    ASSERT_EQ(&M,Mp);
}

TEST(testTKmatrix_d, mallocfree){
    TKmatrix_d m = malloc_matrix_d(128,32);
    free_matrix_d(m);
}

TEST(testTKvector_d, mallocfree){
    TKvector_d m = malloc_vector_d(128);
    free_vector_d(m);
}



