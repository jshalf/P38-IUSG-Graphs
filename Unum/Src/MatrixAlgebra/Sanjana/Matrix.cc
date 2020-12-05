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
    int nRows() { return m.size();}
    int nCols() { return m[0].size();}
    void setSize(int nrows,int ncols){
        if(nrows<0 || nrows<0){fprintf(stderr,"error: setSize() to negative value\n"); return; }
        m.resize(nrows);
        for(int i=0;i<m.size();i++){
            m[i].resize(ncols);
        }
        this->setZero(); // straight resize will set to zero.  Submatrix is different
    }
    // bool isSymmetric(){
    // }
    bool setZero(){
        for(int i=0;i< nRows();i++){
            for(int j=0;j< nCols();j++){
                m[i][j]=0;
            }
        }
    }
    bool setIdentity(){
        if(m.size()<1 || this->nRows()!=this->nCols()) fprintf(stderr,"setIdentity failed on non-square matrix [%u,%u]\n",this->nRows(),this->nCols());
        this->setZero();
        for(int i=0;i<m.size();i++) m[i][i]=1;
        return true;
    }
    
    void printMatrix(){
        unsigned i,j;
        for(j=0;j < this->nCols();j++) printf("\tc[%u]",j);
        printf("\n");
        for(i=0;i < this->nRows();i++){
            printf("r[%u]",i);
            for(j=0;j < this->nCols();j++){
                std::cout << '\t' << m[i][j];
            }
            printf("\n");
        }
        printf("\n");
    }
    void print(){
        unsigned i,j;
        for(i=0;i < this->nRows();i++){
            for(j=0;j < this->nCols();j++){
                std::cout << "[" << i << "," << j << "]=" << m[i][j] << '\t';
            }
            printf("\n");
        }
    }
    
    Matrix(size_t rows, size_t cols){
        setSize(rows,cols);
        setZero();
    }
    
    T &index(int i,int j){
        return m[i][j];
    }

    friend Matrix multiply(Matrix &A, Matrix &B){
		GivenTwoPosit gtp;
		MultiplyControl mc;
		AddPositToDec ptd;
		AddDecToPosit dtp;
		
		Matrix out(A.nRows(),B.nCols()); // ensure output matrix is correct size
        // we should probably check to make sure A.nCols()==B.nRows()....
        out.setZero();
        bitset<32> result = dtp.control(0);
        bitset<32> newResult = dtp.control(0);
        bitset<32> eachResult;
        
        for(int i=0;i<A.nRows();i++){
            for(int j=0;j<B.nCols();j++){
                for(int k=0;k<out.nCols();k++){
				   newResult = 0;
				   if(A.m[i][k] != 0 && B.m[k][j] != 0){
						newResult = mc.mult(A.m[i][k], B.m[k][j]);
						
				   }
				   if(ptd.control(result) == 0){
					   result = newResult;
				   } else {
					   result = gtp.addPosits(result, newResult);
				   }
                }
                out.m[i][j] = ptd.control(result);
                result = 0;
            }
        }
        
        out.printMatrix();
        
        return out;
    }
    
    void transpose(){
        // matrix transpose (row to column
        Matrix t(this->nCols(),this->nRows());
        for(int i=0;i<this->nCols();i++){
            for(int j=0;j<this->nRows();j++){
                t.m[j][i]=this->m[i][j];
            }
        }
        this->m=t.m; // will it deep copy?
    }
    
    friend T innerProduct(vec &a, vec &b){
        T accum=0;
        if(a.size()!=b.size()) {fprintf(stderr,"Error: inner product vectors not same length %u,%u)",a.size(),b.size()); return 0;}
        for(int i=0;i<a.size();i++) accum+=(a[i]*b[i]);
        return accum;
    }
    
    friend vec matVec(const Matrix &A,const vec &v){
        vec out;
        if(A.nCols()!=v.size()){fprintf(stderr,"Error MatVec is mismatched A[%u,%u] v[%u]\n",A.nRows(),A.nCols(),v.size()); return out;}
        for(int i=0;i<A.nRows();i++){
            out[i]=0;
            for(int j=0;j<A.nCols;j++){
                out[j]+=A.m[i][j]*v[j];
            }
        }
        return out;
    }
    
    friend vec vectorCombination(T a,vec &U,T b,vec &V){
        vec W(U.size());
        if(U.size() != V.size()) {fprintf(stderr,"vectorCombo: sizes don't match %u,%u\n",U.size(),V.size()); return W;}
        for(int i=0;i<W.size();i++) W[i] = a * U[i] + b * V[i];
        return W;
    }
    
    friend T vectorNorm(const vec & V){
        return sqrt( innerProduct(V,V) );
    }
    
    bool isSymmetric(){
        if(nRows()==nCols()) return 1;
        else return 0;
    }
    
    bool isUpperTriangular(){
        if(!this->isSymmetric()) return 0;
        if(nRows()<=1) return 0;
        for(int i=0;i<nRows();i++)
            for(int j=0;j<i;j++)
                if(m[i][j]>0) return 0;
        return 1;
    }
    
    bool isLowerTriangular(){
        if(!this->isSymmetric()) return 0;
        if(nRows()<=1) return 0;
        for(int i=0;i<nRows();i++)
            for(int j=i+1;j<nCols();j++)
                if(m[i][j]>0) return 0;
        return 1;
    }
    
    bool solveLowerTriangularSystem(vec &x,const vec b){
        if(!this->isLowerTrangular()) return false;
        // solve sequentially by forward-substitution from top to bottom
        for(int i=0;i<this->nRows();i++){
            T sum=0;
            int j;
            for(j=0;j<(i-1);j++) sum+=m[i][j]*x[j];
            x[i]=(b[i]-sum)/m[i][i];
        }
        return true; // success
    }
    
    bool solveLowerTriangularSystem(const Matrix L, vec &x,const vec b){
        if(!L.isLowerTrangular()) return false;
        // solve sequentially by forward-substitution from top to bottom
        for(int i=0;i<L.nRows();i++){
            T sum=0;
            int j;
            for(j=0;j<(i-1);j++) sum+=L.m[i][j]*x[j];
            x[i]=(b[i]-sum)/L.m[i][i];
        }
        return true; // success
    }
    
    bool solveUpperTriangularSystem(vec &x,const vec b){
        if(!this->isUpperTriangular()) return false;
        // solve sequentially by back-substitution from bottom to top
        for(int i=0;i<nRows();i++){
            T sum=0;
            int j;
            for(j=i;j<nRows();j++) sum+=m[i][j]*x[j];
            x[i]=(b[i]-sum)/m[i][i];
        }
    }
    bool solveUpperTriangularSystem(const Matrix U,vec &x,const vec b){
        if(!U.isUpperTriangular()) return false;
        // solve sequentially by back-substitution from bottom to top
        for(int i=0;i<U.nRows();i++){
            T sum=0;
            int j;
            for(j=i;j<U.nRows();j++) sum+=U.m[i][j]*x[j];
            x[i]=(b[i]-sum)/U.m[i][i];
        }
    }
    
    bool LUdecomposition(Matrix &L, Matrix &U){
        int n=this->nRows();
        //for(int k=0;k<n-1;k++){
        //    int j;
        //    for (j=k+1;j<n;j++) A.m[j][k]=A.m[j][k]/A.m[k][k];
        //    for(j=k;j<n;j++){
        //        for(int i=k;i<n;i++){
        //            A.m[i][j]=A.m[i][j]-A.m[i][j]*A.m[k][j];
        //        }
        //    }
        //}
        // now separate out the pieces of the resulting matrix into upper and lower triangular parts
        
    
    // Decomposing matrix into Upper and Lower
    // triangular matrix
    for (int i = 0; i < n; i++) {
        
        // Upper Triangular
        for (int k = i; k < n; k++) {
            
            // Summation of L(i, j) * U(j, k)
            T sum = 0;
            for (int j = 0; j < i; j++)
                sum += (L.m[i][j] * L.m[j][k]);
            
            // Evaluating U(i, k)
            U.m[i][k] = m[i][k] - sum;
        }
        
        // Lower Triangular
        for (int k = i; k < n; k++) {
            if (i == k)
                L.m[i][i] = 1; // Diagonal as 1
            else {
                
                // Summation of L(k, j) * U(j, i)
                T sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (L.m[k][j] * U.m[j][i]);
                
                // Evaluating L(k, i)
                if(U.m[i][i]==0.0) return false; // fail... cannot divide by zero
                L.m[k][i] = (m[k][i] - sum) / U.m[i][i];
            }
        }
        return true;
    }
    }
    
    bool conjugateGradientSolver(double tolerance, const Matrix &A, const vec &B ,vec & X)
    {
       // double TOLERANCE = 1.0e-10;
        
        int n = A.nRows();
        // vec X( n, 0.0 );
        vec R = B;
        cout << R << endl;
        vec P = R;
        int k = 0;
        for(int i=0;i<n;i++) X[i]=0; // zero out the vector
        
        while ( k < n )
        {
            vec Rold = R;                                         // Store previous residual
            vec AP = matrixTimesVector( A, P );
            
            double alpha = innerProduct( R, R ) / max( innerProduct( P, AP ), NEARZERO );
            X = vectorCombination( 1.0, X, alpha, P );            // Next estimate of solution
            R = vectorCombination( 1.0, R, -alpha, AP );          // Residual
            
            if ( vectorNorm( R ) < tolerance ) break;             // Convergence test
            
            double beta = innerProduct( R, R ) / max( innerProduct( Rold, Rold ), NEARZERO );
            P = vectorCombination( 1.0, R, beta, P );             // Compute Next gradient
            k++;
        }
        
        return true;
    }
    
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
    
    bool getCofactorSubmatrix(int r,int c,Matrix<T> &A, Matrix<T> &S){
        // extract a submatrix (S) from matrix A.
        if(A.nCols()>9) fprintf(stderr,".%u.",A.nCols());
        if(A.nCols() != A.nRows()){fprintf(stderr,"Fail: det() nRows not equal to ncols\n"); return 0;}
        if(A.nCols()<=1){fprintf(stderr,"Fail: matrix only has one col & one rows\n"); return 0;}
        // S.setSize(A.nRows()-1,A.nCols()-1);
        for(int i=0,ii=0;i<S.nRows();i++,ii++){
            if(i==r) ii++; // skip row
            for(int j=0,jj=0;j<S.nRows();j++,jj++){
                if(j==c) jj++; // skip col
                S.m[i][j]=A.m[ii][jj];
            }
        }
        return true;
    }
    
    bool getCofactorMatrix(Matrix <T> &A){
        // not yet implemented
        return true;
    }
    
    bool getInvertedMatrix(Matrix <T> &A){
        // invert self... probably should just be a subroutine
        return true;
    }
    
    /* Recursive function for finding determinant of matrix.
     n is current dimension of mat[][]. */
    T det() {
        T accum = 0; // Initialize result
        if(this->nRows() != this->nCols()){fprintf(stderr,"Fail: det() nRows not equal to ncols\n"); return 0;}
        //  Base case : if matrix contains single element
        if (this->nRows() == 1)
            return m[0][0];
        else if (this->nRows()==2) return (m[0][0]*m[1][1] - m[0][1]*m[1][0]);
        Matrix<T> A(this->nRows()-1,this->nCols()-1); // To store cofactors  (should it be cols-1?)
        
        int sign = 1;  // To store sign multiplier
        
        // Iterate for each element of first row
        for (int j = 0; j < this->nCols(); j++)
        {
            // Getting Cofactor of mat[0][f]
            //getCofactor(mat, temp, 0, f, n);
            getCofactorSubmatrix(0,j,*this,A);   // construct a submatrix of cofactors
            accum += sign * m[0][j] * A.det();
            
            // terms are to be added with alternate sign
            sign = -sign;
        }
        return accum;
    }
    
    // compute the condition number for a matrix
    
    bool load(char *filename){
        // from MMIO example
        int ret_code;
        FILE *f;
        int M, N, nz,n;
        double *val;
        MM_typecode matcode;

        if ((f = fopen(filename, "r")) == NULL)
            exit(1);
        
        if (mm_read_banner(f, &matcode) != 0)
        {
            printf("Could not process Matrix Market banner.\n");
            exit(1);
        }
        if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
            mm_is_sparse(matcode) )
        {
            printf("Sorry, this application does not support ");
            printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
            exit(1);
        }
        
        if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
            exit(1);
        this->setSize(M,N);
        for (n=0; n<nz; n++)
        {
            int i,j;
            T val;
            fscanf(f, "%d %d %lg\n", &i, &j, &val);
            m[i-1][j-1] = val;
        }
        
        if (f !=stdin) fclose(f);
        
        return 1; // load a matrix from a matrixmarket file
    }
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
