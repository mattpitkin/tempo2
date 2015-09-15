#include <gtest/gtest.h>
#include <src/gtest_main.cc>
#include <iostream>
#include <iomanip>
#include <cfloat>

#include"TKmatrix.h"
#include"TKsvd.h"
#include"TKlog.h"



TEST(testTKmatrix_d, construct){
    TK::matrix<double> M(1,2);
    ASSERT_EQ(M._rows,1);
    ASSERT_EQ(M._cols,2);
}

TEST(testTKmatrix_d, getset){
    TK::matrix<double> M(2,2,true);
    M.set(0,0,1.0);
    M.set(0,1,-2.0);
    M.set(1,0,2.0);

    ASSERT_EQ(M.get(0,0),1.0);
    ASSERT_EQ(M.get(0,1),-2.0);
    ASSERT_EQ(M.get(1,0),2.0);
}

TEST(testTKmatrix_d, transpose){
    TK::matrix<double> M(3,3,true);
    M.set(0,0,1.0);
    M.set(0,1,-2.0);
    M.set(1,0,2.0);

    M.set(1,2,-3.0);
    M.set(2,1,3.0);

    M.set(0,2,-4.0);
    M.set(2,0,4.0);

    {
        TK::matrix<double> MT = M.T(false);

        ASSERT_EQ(MT.get(0,0), 1.0);
        ASSERT_EQ(MT.get(1,0),-2.0);
        ASSERT_EQ(MT.get(0,1), 2.0);
        ASSERT_EQ(MT.get(2,1),-3.0);
        ASSERT_EQ(MT.get(1,2), 3.0);
        ASSERT_EQ(MT.get(2,0),-4.0);
        ASSERT_EQ(MT.get(0,2), 4.0);

        // check that the raw data was transposed
        for(size_t j =0 ; j < 3; j++){
            for(size_t k =0 ; k < 3; k++){
                size_t i1 = j+3*k;
                size_t i2 = k+3*j;
                EXPECT_EQ(MT.getRaw()[i1],M.getRaw()[i2]);
            }
        }
    }

    {
        TK::matrix<double> MT = M.T(true);

        ASSERT_EQ(MT.get(0,0), 1.0);
        ASSERT_EQ(MT.get(1,0),-2.0);
        ASSERT_EQ(MT.get(0,1), 2.0);
        ASSERT_EQ(MT.get(2,1),-3.0);
        ASSERT_EQ(MT.get(1,2), 3.0);
        ASSERT_EQ(MT.get(2,0),-4.0);
        ASSERT_EQ(MT.get(0,2), 4.0);

        for(size_t i =0 ; i < 9; i++){
            EXPECT_EQ(MT.getRaw()[i],M.getRaw()[i]);
        }
    }

}

TEST(testTKmatrix_d, transposeBig){
    TK::matrix<double> M(1024,2048,true);
    for (int i=0;i<1024*2048;i++){
        M.getRaw()[i]=(double)i;
    }
    TK::matrix<double> MT = M.T(false);
    ASSERT_EQ(MT._rows,M._cols);
    ASSERT_EQ(MT._cols,M._rows);

    for (int i=0;i<1024;i++){
        for (int j=0;j<2048;j++){
            ASSERT_EQ(M.get(i,j),MT.get(j,i));
        }
    }
}

TEST(testTKmatrix_d, transposeBig2){
    TK::matrix<double> M(1024,2048,true);
    for (int i=0;i<1024*2048;i++){
        M.getRaw()[i]=(double)i;
    }
    TK::matrix<double> MT = M.T(true);
    ASSERT_EQ(MT._rows,M._cols);
    ASSERT_EQ(MT._cols,M._rows);

    for (int i=0;i<1024;i++){
        for (int j=0;j<2048;j++){
            ASSERT_EQ(M.get(i,j),MT.get(j,i));
        }
    }
}

TEST(testTKmatrix_d, getIdx){
    TK::matrix<double> M(2,2,true);
    M.set(0,0,1.0);
    M.set(0,1,-2.0);
    M.set(1,0,2.0);
    TKmatrix_d m = M.getIdx();

    ASSERT_EQ(m[0][0],1.0);
    ASSERT_EQ(m[0][1],-2.0);
    ASSERT_EQ(m[1][0],2.0);

}

TEST(testTKmatrix_d, opPlus){
    TK::matrix<double> m1(16,32);
    TK::matrix<double> m2(16,32);
    TK::matrix<double> m2b(16,32,false);
    for (int i=0;i<16*32;i++){
        m1.getRaw()[i]=i;
        m2.getRaw()[i]=16.5-i;
        m2b.getRaw()[i]=16.5-i;
    }
    TK::matrix<double> m3 = m1+m2;
    TK::matrix<double> m4 = m1+m2b;

    for (int i=0;i<16;i++){
        for (int j=0;j<32;j++){
            ASSERT_EQ(m3.get(i,j), m1.get(i,j)+m2.get(i,j));
            ASSERT_EQ(m4.get(i,j), m1.get(i,j)+m2b.get(i,j));
        }
    }
}


TEST(testTKmatrix_d, opMult){
    TK::matrix<double> m1(5,3);
    TK::matrix<double> m2(3,2);
    TK::matrix<double> m2b(3,2,false);

    for (int i=0;i<5;i++){
        for (int j=0;j<3;j++){
            double d = i*100.0+j;
            m1.set(i,j,d);
        }
    }
    for (int i=0;i<2;i++){
        for (int j=0;j<3;j++){
            double d = i*100.0+j;
            m2.set(j,i,d*M_PI);
            m2b.set(j,i,d*M_PI);
        }
    }
    TK::matrix<double> m3 = m1*m2;
    TK::matrix<double> m4 = m1*m2b;

    for (int i=0;i<5;i++){
        for (int j=0;j<2;j++){
            double element=0;
            for (int k=0; k < 3; k++){
                element += m1.get(i,k)*m2.get(k,j);
            }
            ASSERT_EQ(element,m3.get(i,j));
            ASSERT_EQ(element,m3.get(i,j));
        }
    }
}


TEST(testTKmatrix_d, toFromPtr){
    TK::matrix<double> M(2,2,true);
    double* p = M.getRaw();

    TK::matrix<double>* Mp = TK::matrix<double>::getFromPtr(p);
    ASSERT_EQ(Mp,&M);
}

TEST(testTKmatrix_d, mallocfree){
    TKmatrix_d m = malloc_matrix_d(128,32);
    free_matrix_d(m);
}

TEST(testTKvector_d, mallocfree){
    TKvector_d m = malloc_vector_d(128);
    free_vector_d(m);
}



TEST(testTKvector_d, opMVmult){
    TK::matrix<double> m1(2,3,true);
    TK::matrix<double> m2(2,3,false);
    TK::vector<double> v(3);
    for (int i=0;i<2;i++){
        for (int j=0;j<3;j++){
            double d = i*100.0+j;
            m1[i][j] = d;
            m2[i][j] = d;
        }
    }

    for (int j=0;j<3;j++){
        v[j] = j+1;
    }
    TK::vector<double> u = m1*v;
    TK::vector<double> u2 = m2*v;

    for (int i=0;i<2;i++){
        double element=0;
        for (int k=0; k < 3; k++){
            element += m1.get(i,k)*v[k];
        }
        ASSERT_EQ(element,u[i]);
        ASSERT_EQ(element,u2[i]);
    }
}


TEST(testTKvector_d, multMatrixVec){
    TKmatrix_d m1 = malloc_matrix_d(2,3);
    TKvector_d v  = malloc_vector_d(3);
    for (int i=0;i<2;i++){
        for (int j=0;j<3;j++){
            double d = i*100.0+j;
            m1[i][j] = d;
        }
    }

    for (int j=0;j<3;j++){
        v[j] = j+1;
    }
    TKvector_d u = malloc_vector_d(2);
    TKmultMatrixVec_d(m1,v,u);

    for (int i=0;i<2;i++){
        double element=0;
        for (int k=0; k < 3; k++){
            element += m1[i][k]*v[k];
        }
        ASSERT_EQ(element,u[i]);
    }
}


TEST(testTKdiagonal_d, alloc){
    TK::diagonal<double> M(8);
    ASSERT_EQ(0,M[0][0]);
    M[2][2] = 1;
    ASSERT_EQ(1,M[2][2]);
}

TEST(testTKdiagonal_d, opMult){
    TK::diagonal<double> M1(8);
    TK::diagonal<double> M2(8);
    TK::matrix<double> M3(8,8);

    for (int i=0;i<8;i++){
        M1[i][i] = i;
        M2[i][i] = 1.0/(double)(i+1);
        M3[i][i] = M2[i][i];
    }

    TK::diagonal<double> O1 = M1*M2;
    TK::matrix<double> O2 = M1*M3;

    for (int i=0;i<8;i++){
        ASSERT_EQ(M1[i][i]*M2[i][i], O1[i][i]);
        ASSERT_EQ(M1[i][i]*M3[i][i], O2[i][i]);
    }
}


TEST(testTKmatrix_ld, construct){
    TK::matrix<longdouble> M(1,2);
    ASSERT_EQ(M._rows,1);
    ASSERT_EQ(M._cols,2);
}

TEST(testTKmatrix_ld, getset){
    TK::matrix<longdouble> M(2,2,true);
    M.set(0,0,1.0);
    M.set(0,1,-2.0);
    M.set(1,0,2.0);

    ASSERT_EQ(M.get(0,0),longdouble(1.0));
    ASSERT_EQ(M.get(0,1),longdouble(-2.0));
    ASSERT_EQ(M.get(1,0),longdouble(2.0));
}


TEST(testTKmatrix_ld, opMult){
    debugFlag=1;
    TK::matrix<longdouble> m1(2,3);
    TK::matrix<longdouble> m2(3,2);
    TK::matrix<longdouble> m2b(3,2,false);
    for (int i=0;i<2;i++){
        for (int j=0;j<3;j++){
            longdouble d = i*longdouble(100.0)+j;
            m1.set(i,j,d);
            m2.set(j,i,d*LD_PI);
            m2b.set(j,i,d*LD_PI);
        }
    }
    TK::matrix<longdouble> m3 = m1*m2;
    TK::matrix<longdouble> m4 = m1*m2b;

    for (int i=0;i<2;i++){
        for (int j=0;j<2;j++){
            longdouble element=0;
            for (int k=0; k < 3; k++){
                element += m1.get(i,k)*m2.get(k,j);
            }
            ASSERT_EQ(element,m3.get(i,j));
            ASSERT_EQ(element,m3.get(i,j));
        }
    }
}



TEST(svd,trivial_d){
    std::cout << std::scientific << std::setprecision(3);
    TK::matrix<double> m1(6,3);
    TK::matrix<double> m2(6,3,false);
    for (int i=0;i<6;i++){
        for (int j=0;j<3;j++){
            if(i==j){
                m1[i][j] = i+1;
                m2[i][j] = i+1;
            }
        }
    }
    m1[0][1]=-0.4;
    m1[3][1]=4.2;
    m2[0][1]=-0.4;
    m2[3][1]=4.2;
    {
        TK::SVD<double> svd(m1);
        ASSERT_EQ(svd.U._rows,m1._rows);
        TK::matrix<double> o = svd.U*(svd.S*svd.V.T(true));
        for (int i=0;i<6;i++){
            for (int j=0;j<3;j++){
                EXPECT_NEAR(m1[i][j],o[i][j],5*svd.S[0][0]*DBL_EPSILON) << "Error in rowmajor SVD";
            }
        }
    }

    {
        TK::SVD<double> svd(m2);
        std::cout << svd.S << std::endl;
        ASSERT_EQ(svd.U._rows,m2._rows);
        TK::matrix<double> o = svd.U*(svd.S*svd.V.T(true));
        for (int i=0;i<6;i++){
            for (int j=0;j<3;j++){
                EXPECT_NEAR(m2[i][j],o[i][j],5*svd.S[0][0]*DBL_EPSILON) << "Error in colmajor SVD";
            }
        }
    }



}
