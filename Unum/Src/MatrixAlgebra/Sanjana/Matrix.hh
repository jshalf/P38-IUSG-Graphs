#ifndef __Matrix_hh__
#define __Matrix_hh__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "mmio.hh"
#include "Addition.cpp"
#include "Multiplication.cpp"

using namespace std;

const double NEARZERO = 1.0e-10;
template<class T>
class Matrix {
    using vec = vector<T>;
    using mat = vector<vec>;
public:
    mat m;
public:
    inline int nRows() { return m.size();}
    inline int nCols() { return m[0].size();}
    void setSize(int nrows,int ncols);
    // bool isSymmetric(){
    // }
    bool setZero();
    bool setIdentity();
    
    void printMatrix();
    void print();
    
    Matrix(size_t rows, size_t cols);
    
    inline T &index(int i,int j){ return m[i][j];}

    friend Matrix multiply(Matrix &A, Matrix &B);
    
    void transpose();
    
    friend T innerProduct(vec &a, vec &b);
    
    friend vec matVec(const Matrix &A,const vec &v);
    
    friend vec vectorCombination(T a,vec &U,T b,vec &V);
    
    friend T vectorNorm(const vec & V);
    
    bool isSymmetric();
    
    bool isUpperTriangular();
    
    bool isLowerTriangular();
    
    bool solveLowerTriangularSystem(vec &x,const vec b);
    
    bool solveLowerTriangularSystem(const Matrix L, vec &x,const vec b);
    
    bool solveUpperTriangularSystem(vec &x,const vec b);
    bool solveUpperTriangularSystem(const Matrix U,vec &x,const vec b);
    
    bool LUdecomposition(Matrix &L, Matrix &U);
    
    bool conjugateGradientSolver(double tolerance, const Matrix &A, const vec &B ,vec & X);
    
    //T det(Matrix<T> &A){
    //    // do right diagonals - left diagonals with wrap-around (modulo)
    //   int i,j,k;
    //   T accum=0;
    //   if(A.nCols()!=A.nRows()){fprintf(stderr,"Fail: det() nRows not equal to ncols\n"); return 0;}
    //   else if(A.nCols()==1) return A.m[0][0];
    //   for(k=0;k<A.nCols();k++){
    //       T a=1.0,b=1.0;
    //       for(i=0;i<A.nCols();i++){
    //           j=k; // start row/column
    //           a*=A.m[i][j%(A.nCols())]; // positive diagonal contribution
    //           b*=A.m[A.nCols()-1][j%(A.nCols())]; // negative diagonal contribution
    //       }
    //       printf("A=%lf B=%lf\n",a,b);
    //       accum+=(a-b); // sume the positive and negative diagonal contributions together.
    //  }
    //  return accum;
    //}
    
    bool getCofactorSubmatrix(int r,int c,Matrix<T> &A, Matrix<T> &S);
    
    bool getCofactorMatrix(Matrix <T> &A);
    
    bool getInvertedMatrix(Matrix <T> &A);
    
    /* Recursive function for finding determinant of matrix.
     n is current dimension of mat[][]. */
    T det();
    
    // compute the condition number for a matrix
    
    bool load(char *filename);
};

int main(){
	/*AddControl a;
	a.add();
	MultiplyControl m;
	m.mult();*/
	
    Matrix<double> A(10,10), B(10,10),C(10,10),D(3,3),E(2,2);
    using vec = vector<double>;
    vec b(2);
    
  //  A.printMatrix();
  //  B.printMatrix();
    
    
    //C.load("can24.mtx");
   // cout << "C[" << C.nRows() << ',' << C.nCols() << ']' << '\n';
    // C.printMatrix();
   // C.det();
   // cout << "DetC()=" << A.det() << '\n';
    
    D.m[0][0]= 3; D.m[0][1]=2; D.m[0][2]=0;
    D.m[1][0]= 1; D.m[1][1]= 0; D.m[1][2]=1;
    D.m[2][0]=2; D.m[2][1]= 3; D.m[2][2]=0;
    
    D.printMatrix();
    D.printMatrix();
    A = multiply(D, D);
    
    //D.m[0][0]= 1; D.m[0][1]= 2; D.m[0][2]=3;
    //D.m[1][0]= 4; D.m[1][1]= 5; D.m[1][2]=6;
    //D.m[2][0]= 7; D.m[2][1]= 8; D.m[2][2]=9;
    
   // cout << "DetD()=" << D.det() << '\n';
   // D.transpose();
   // cout << "DetD(Transpose)" << D.det() << '\n';
    
    // now implement a matrix inversion using the cofactor matrix
    /*
    A=multiply(D,D);
    puts("Print A");
    D.printMatrix();
    A.printMatrix();*/
    
    return 1;
    
    b[0]=1; b[1]=2;
    A.setIdentity();
    A.m[2][5]=5;
    A.printMatrix();
    printf("IsLowerTriangular  isUpperTriangular %u,%u\n",A.isLowerTriangular(),A.isUpperTriangular());
    A.transpose();
    A.printMatrix();
    printf("IsLowerTriangular  isUpperTriangular %u,%u\n",A.isLowerTriangular(),A.isUpperTriangular());
    cout << "DetA()=" << A.det() << '\n';
    return 0;
    B.setIdentity();
    B.setZero();
    A.print();
    A.printMatrix();
    // Matrix<double>::multiply(A,B,C);
    // multiply(A,B,C);
    C.printMatrix();
    return 0;
}

#endif

