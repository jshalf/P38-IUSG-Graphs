
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include "PositBase.hh"
#include "Posit32.hh"
#include "Matrix.hh"
#include "SuperMatrix.hh"
#include "half.hpp"
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstdlib>
using namespace std;

void RegressionTest1(){
    Matrix<double> A(10,10), B(10,10),C(10,10),D(3,3),E(2,2);
    using vec = vector<double>;
    vec b(2);
    
    D.m[0][0]= 3; D.m[0][1]=-2; D.m[0][2]=0;
    D.m[1][0]= 1; D.m[1][1]= 0; D.m[1][2]=1;
    D.m[2][0]=-2; D.m[2][1]= 3; D.m[2][2]=0;
    
    //D.m[0][0]= 1; D.m[0][1]= 2; D.m[0][2]=3;
    //D.m[1][0]= 4; D.m[1][1]= 5; D.m[1][2]=6;
    //D.m[2][0]= 7; D.m[2][1]= 8; D.m[2][2]=9;
    
    // cout << "DetD()=" << D.det() << '\n';
    // D.transpose();
    // cout << "DetD(Transpose)" << D.det() << '\n';
    
    // now implement a matrix inversion using the cofactor matrix
    
    A=multiply(D,D);
    puts("Print A");
    D.print();
    A.print();
}

void regressionTest2(){
    Matrix<double> A(10,10), B(10,10),C(10,10),D(3,3),E(2,2);
    using vec = vector<double>;
    vec b(2);
    
    b[0]=1; b[1]=2;
    A.setIdentity();
    A.m[2][5]=5;
    A.print();
    printf("IsLowerTriangular  isUpperTriangular %u,%u\n",A.isLowerTriangular(),A.isUpperTriangular());
    A.transpose();
    A.print();
    printf("IsLowerTriangular  isUpperTriangular %u,%u\n",A.isLowerTriangular(),A.isUpperTriangular());
    cout << "DetA()=" << A.det() << '\n';
    
    B.setIdentity();
    B.setZero();
    A.print();
    A.print();
    // Matrix<double>::multiply(A,B,C);
    // multiply(A,B,C);
    C.print();
}


int testInnerProduct() {
    cout << "Testing inner product..." << endl;
    bool success = 0;
    Matrix<double> U;
    vector<double> a(3), b(3);
    
    a[0] = 0; a[1] = 2; a[2] = 3;
    b[0] = 5; b[1] = -2; b[2] = 0;
    
    if (U.innerProduct(a, a) != 13 || U.innerProduct(a, a) != sqrt (U.vectorNorm(a))) success = 0;
    
    if (U.innerProduct(b, a) != -4) success = 0;
    
    cout << "Good" << endl << endl;
    return success;
} 

int testMatVec() {
    bool success = 0;
    Matrix<double> U(2, 2);
    vector<double> a(2), b(2);
    U.setIdentity();
    a[0] = 1; a[1] = -1; 
    b[0] = 3; b[1] = 0; 
    if (matVec(U, a) != a) success = 0;
    if (matVec(U, b) != b) success = 0;
    if (success == 0) cout << "Identity matrix vector multiplication good." << endl;
    cout << "Testing mat vec multiply. Should see A, x1, Ax1, x2, Ax2, ..." << endl;

    U.m[0][0] = 1; U.m[0][1] = 2;
    U.m[1][0] = 2.5; U.m[1][1] = -5;

    //Doing some more matVec operations.
    U.print();
    double g = U.m[1][0];
    printVector(a);
    printVector(matVec(U, a));
    cout << endl;
    printVector(b);
    printVector(matVec(U, b));
    return 1;
} 

int testMultiply() {
    Matrix<double> I(3,3), A(3, 3), B(3, 3), C(3, 4);
    I.setIdentity();

    A.m[0][0] = 1; A.m[0][1] = 2; A.m[0][2] = 3;
    A.m[1][0] = 3; A.m[1][1] = 4; A.m[1][2] = 1;
    A.m[2][0] = 2; A.m[2][1] = 1; A.m[2][2] = 0;

    B.m[0][0] = 2; B.m[0][1] = -3; B.m[0][2] = 5;
    B.m[1][0] = -4; B.m[1][1] = 2; B.m[1][2] = 0;
    B.m[2][0] = 2; B.m[2][1] = -3; B.m[2][2] = 1;

    C.m[0][0] = 1; C.m[0][1] = 1; C.m[0][2] = 1; C.m[0][3] = 1;
    C.m[1][0] = 1; C.m[1][1] = 1; C.m[1][2] = 1; C.m[1][3] = 1;
    C.m[2][0] = 1; C.m[2][1] = 1; C.m[2][2] = 1; C.m[2][3] = 1;
    
    cout << "Multiplying by identity, should see two of each matrix..." << endl;
    A.print();
    multiply(I, A).print();
    B.print();
    multiply(I, B).print();
    C.print();
    multiply(I, C).print();

    cout << "More multiplications of form op1, op2, op1*op2 ..." << endl;
    A.print();
    B.print();
    multiply(A, B).print();
    B.print();
    C.print();
    multiply(B, C).print();

    return 1;
} 

int testTriangleCheckers() {
    cout << "Checking isUpperTriangular and isLowerTriangular" << endl;
    bool success = 1;
    Matrix<double> I(3,3), A(3, 3), B(3, 3);
    
    I.setIdentity();
    A.m[0][0] = 1; A.m[0][1] = 2; A.m[0][2] = 3;
    A.m[1][0] = 0; A.m[1][1] = 4; A.m[1][2] = 1;
    A.m[2][0] = 0; A.m[2][1] = 1; A.m[2][2] = 0;
    
    if (!(I.isUpperTriangular())) success = 0;
    if (A.isUpperTriangular()) success = 0;
    A.m[2][1] = 0; //Now A is upper triangular.
    if (!A.isUpperTriangular()) success = 0;

    B.m[0][0] = 2; B.m[0][1] = 0; B.m[0][2] = 5;
    B.m[1][0] = -4; B.m[1][1] = 0; B.m[1][2] = 0;
    B.m[2][0] = 2; B.m[2][1] = -3; B.m[2][2] = 2;
    if (!(I.isLowerTriangular())) success = 0;
    if (B.isLowerTriangular()) success = 0;
    B.m[0][2] = 0;
    if (!(B.isLowerTriangular())) success = 0;
    if (success) cout << "Passed tests" << endl;
    else cout << "failed" << endl;
    return success;
}

void testUpperLowerTriSolves() {
    cout << "Checking tri solves" << endl;
    bool success = 1;
    Matrix<double> A(3, 3), B(3, 3);
    vector<double> b(3), b2(3), x(3); //Ax = b
    b[0] = 1; b[1] = 1; b[2] = 1;

    cout << "Form A, x, b, Ax. Should see Ax ~= b." << endl;
    
    A.m[0][0] = 1; A.m[0][1] = 1; A.m[0][2] = 1;
    A.m[1][0] = 0; A.m[1][1] = 1; A.m[1][2] = 1;
    A.m[2][0] = 0; A.m[2][1] = 0; A.m[2][2] = 1;
    A.solveUpperTriangularSystem(x, b);
    A.print();
    printVector<double>(x);
    printVector<double>(b);
    printVector<double>(matVec(A, x));

    A.m[0][0] = 6.5; A.m[0][1] = -1; A.m[0][2] = -33;
    A.m[1][0] = 0; A.m[1][1] = 2; A.m[1][2] = 6;
    A.m[2][0] = 0; A.m[2][1] = 0; A.m[2][2] = 3;
    b[0] = 1.5; b[1] = -2.4; b[2] = 1.05;
    A.solveUpperTriangularSystem(x, b);
    A.print();
    printVector<double>(x);
    printVector<double>(b);
    printVector<double>(matVec(A, x));

    A.setSize(5, 5);  
    A.m[0][0] = 6.75; A.m[0][1] = -9.37; A.m[0][2] = -33; A.m[0][3] = -1; A.m[0][4] = 3;
    A.m[1][0] = 0; A.m[1][1] = 2.314; A.m[1][2] = 6; A.m[1][3] = 80; A.m[1][4] = -56.74;
    A.m[2][0] = 0; A.m[2][1] = 0; A.m[2][2] = 3.9998; A.m[2][3] = 0; A.m[2][4] = -70.4;
    A.m[3][0] = 0; A.m[3][1] = 0; A.m[3][2] = 0; A.m[3][3] = -78.9; A.m[3][4] = 45;
    A.m[4][0] = 0; A.m[4][1] = 0; A.m[4][2] = 0; A.m[4][3] = 0; A.m[4][4] = -19;
    b.resize(5);
    x.resize(5);
    b[0] = 1.5; b[1] = -2.4; b[2] = 1.05; b[3] = 19; b[4] = 0.01;
    A.solveUpperTriangularSystem(x, b);
    A.print();
    printVector<double>(x);
    printVector<double>(b);
    printVector<double>(matVec(A, x));

    
    x.resize(3); b.resize(3);
    B.m[0][0] = 6.5; B.m[0][1] = 0; B.m[0][2] = 0;
    B.m[1][0] = -1; B.m[1][1] = 2; B.m[1][2] = 0;
    B.m[2][0] = -33; B.m[2][1] = 6; B.m[2][2] = 3;
    b[0] = 1.508; b[1] = -2.4; b[2] =  -19.07;
    B.solveLowerTriangularSystem(x, b);
    B.print();
    printVector<double>(x);
    printVector<double>(b);
    printVector<double>(matVec(B, x));
}

void testLU() {
    cout << "Testing LU decomp. Of form A, L, U, L*U. A should match L*U." << endl;
    Matrix<double> A(3, 3), L(3, 3), U(3, 3);
    A.m[0][0] = 2; A.m[0][1] = -1; A.m[0][2] = -2;
    A.m[1][0] = -4; A.m[1][1] = 6; A.m[1][2] = 3;
    A.m[2][0] = -4; A.m[2][1] = -2; A.m[2][2] = 8;
    A.LUdecomposition(L, U);
    A.print();
    L.print();
    U.print();
    multiply(L, U).print();
}

void testTriSolve() {
    cout << "Testing TriSolve" << endl;
    cout << "Outputs of form A, x, b, A*x. Should see A*x ~= b." << endl;
    Matrix<double> A(3, 3);
    vector<double> x(3), b(3);
    A.m[0][0] = .908; A.m[0][1] = -1.65732; A.m[0][2] = -2.045678;
    A.m[1][0] = -4243; A.m[1][1] = 6; A.m[1][2] = 3789;
    A.m[2][0] = .00032452; A.m[2][1] = -27.8905; A.m[2][2] = -89;
    b[0] = 1.5; b[1] = -9.7; b[2] = 19.09;
    
    A.triSolve(x, b);
    A.print();
    printVector<double>(x);
    printVector<double>(b);
    printVector<double>(matVec(A, x)); 
}

void testHalfPrecision() {
    cout << "Testing Half precision" << endl;
    Matrix<double> A(5, 5); 
    A.m[0][0] = 1000; A.m[0][1] = 0; A.m[0][2] = 0; A.m[0][3] = 0; A.m[0][4] = 0;
    A.m[1][0] = 0; A.m[1][1] = 3; A.m[1][2] = 0; A.m[1][3] = 0; A.m[1][4] = 0;
    A.m[2][0] = 0; A.m[2][1] = 0; A.m[2][2] = 3; A.m[2][3] = 0; A.m[2][4] = 0;
    A.m[3][0] = 0; A.m[3][1] = 0; A.m[3][2] = 0; A.m[3][3] = 3; A.m[3][4] = 0;
    A.m[4][0] = 0; A.m[4][1] = 0; A.m[4][2] = 0; A.m[4][3] = 0; A.m[4][4] = 7;

    //Matrix<half> AT = A;
    //AT.transpose();
    //Matrix<half> ATA = multiply(AT, A);
    //ATA.print();

    vector<double> x(5), b(5);
    for (int i=0;i<5;i++) x[i] = 0;
    for (int i=0;i<5;i++) b[i] = 1;

    A.conjugateGradientSolver(0, A, b, x);
    cout << A.getMin() << endl;
    cout << A.getMax() << endl;

    printVector(matVec(A, x));
}

void testPosit32s() {
    //Posit32s a(1), b(2);
    //Posit32s r = a+b;
    //cout << r << endl;
}

void testLUHalf(char* matrix, mpf_class xmax, mpf_class theta) {
    //Matrix<mpf_class> systemM;
    //systemM.loadMPF(matrix);
    //int n = systemM.nCols();
    //Matrix<mpf_class> rescaled = systemM.rescale(0);
//
    //Matrix<half> systemH(n,n);
    //systemH.set(rescaled);
    //Matrix<half> L(n,n), U(n,n);
    //systemH.LUdecomposition(L, U);
}

void testPositLoading() {
    Matrix<Posit32gmp> M;
    M.loadMPF("allMatrices/nos4.mtx");
    M.print();
}

int main(int argc, char*argv[]){
    Matrix<double> A, L, U;
    A.load(argv[1]);
    A.LUdecomposition(L, U);

    vector<vector<double>> x(L.nCols());
    vector<vector<double>> b(L.nCols());
    for (int i=0;i<L.nCols();i++) b[i] = vector<double> {1};

    SuperMatrix<double> S(L);
    S.solveLowerTriangularRandom(x, b);

    vector<double> x1(A.nCols());
    for (int i=0;i<A.nCols();i++) x1[i]=x[i][0];

    vector<double> result = matVec(L, x1);
    printVector(result);

    return 0;
}







