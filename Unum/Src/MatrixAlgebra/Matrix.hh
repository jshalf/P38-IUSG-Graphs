#ifndef __MATRIX_HH_
#define __MATRIX_HH_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <bitset>
#include <gmpxx.h>
#include <iomanip>
#include <numeric> 
#include <fenv.h>
#include <stdio.h>
#include "mmio.hh"
#include "helpers.hh"
using namespace std;

const double NEARZERO = 1.0e-10;
template<class T>
class Matrix {
    using vec = vector<T>;
    using mat = vector<vec>;
public:
    mat m;
public:

    int nRows() const { return m.size();}
    int nCols() const { return m[0].size();}
    
    void setSize(int nrows,int ncols){
        if(nrows<0 || nrows<0){fprintf(stderr,"error: setSize() to negative value\n"); return; }
        m.resize(nrows);
        for(int i=0;i<m.size();i++){
            m[i].resize(ncols);
        }
        this->setZero(); // straight resize will set to zero.  Submatrix is different
    }
    // bool isSquare(){
    // }
    bool setZero(){
        for(int i=0;i< nRows();i++){
            for(int j=0;j< nCols();j++){
                m[i][j] = 0;
            }
        }
        return true;
    }
    bool setIdentity(){
        if(m.size()<1 || this->nRows()!=this->nCols()) fprintf(stderr,"setIdentity failed on non-square matrix [%u,%u]\n",this->nRows(),this->nCols());
        this->setZero();
        for(int i=0;i<m.size();i++) {
            m[i][i] = 1; 
        } 
        return true;
    }

    void print() const {
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
    
    Matrix(){
    }
 
    template <typename E>
    Matrix(const Matrix<E> &A){
        this->setSize(A.nRows(), A.nCols());
        for(int i=0;i<A.nRows();i++) {
            for(int j=0;j<A.nCols();j++) {
                cast(this->m[i][j], A.m[i][j]);
            }
        }
    }
    
    T &index(int i,int j){
        return m[i][j];
    }

    int nnz() const {
        int nnz=0;
        for(int i=0;i<nRows();i++) {
            for(int j=0;j<nCols();j++) {
                if (m[i][j] != 0) { nnz++; }
            }
        }
        return nnz;
    }
    
    /*friend Matrix multiply(Matrix &A, Matrix &B){
        Matrix out(A.nRows(),B.nCols()); // ensure output matrix is correct size
        // we should probably check to make sure A.nCols()==B.nRows()....
        out.setZero();
        for(int i=0;i<A.nRows();i++){
            for(int j=0;j<B.nCols();j++){
                for(int k=0;k<out.nCols();k++){
                    cout << i << j << k << endl;
                    out.m[i][j]+=(A.m[i][k] * B.m[k][j]);
                }
            }
        }
        return out;
    } */

    void set(const Matrix<mpf_class> &M) {
        this->setSize(M.nRows(), M.nCols());
        for (int i=0;i<nRows();i++) downcast(m[i], M.m[i]);
    }

    //Changed bound on k.
    friend Matrix multiply(const Matrix &A, const Matrix &B){
        Matrix out(A.nRows(),B.nCols()); // ensure output matrix is correct size
        // we should probably check to make sure A.nCols()==B.nRows()....
        out.setZero();
        for(int i=0;i<A.nRows();i++){
            for(int j=0;j<B.nCols();j++){
                for(int k=0;k<A.nCols();k++){
                    out.m[i][j]+=(A.m[i][k] * B.m[k][j]);
                }
            }
        }
        return out;
    } 

    //Retrieve rectangular submatrix
    Matrix getRectangularSubMatrix(int i, int j, int length, int width) {
        if (i+length >= nRows() || j+width >= nCols()) 
            { fprintf(stderr, "Submatrix out of bounds"); }

        Matrix subMatrix(length, width);
        
        for (int k; k<length; k++) {
            for (int l; l<width; l++) {
                subMatrix.m[k][l]=m[i+k][j+l];
            }
        }

        return subMatrix;
    }

    //Retrieve triangular submatrix from columns i...j
    Matrix getTriangularSubMatrix(int start, int end) {
        if (end<=start || end >= nCols() || !isLowerTriangular())
            { fprintf(stderr, "Invalid Arguments"); }

        Matrix subMatrix(end-start, end-start);
        
        for (int i=start;i<end;i++) {
            for (int j=start;j<i+1;j++) {
                subMatrix.m[i][j]=m[i][j];
            }
        }

        return subMatrix;
    }

    vec getSubVector(vec x, int start, int end) {
        return vector<T>(x.begin()+start, x.begin()+end);
    }

    friend Matrix multiplyDiagRight(const Matrix &A, const Matrix &D){
        Matrix out(A.nRows(),A.nCols()); // ensure output matrix is correct size
        
        out.setZero();
        for (int j=0;j<A.nCols();j++) {
            for (int i=0;i<A.nRows();i++) {
                out.m[i][j] = A.m[i][j]*D.m[j][j];
            }
        }
        return out;
    } 

    friend Matrix multiplyDiagLeft(const Matrix &D, const Matrix &A){
        Matrix out(A.nRows(),A.nCols()); // ensure output matrix is correct size
        
        out.setZero();
        for (int i=0;i<A.nRows();i++) {
            for (int j=0;j<A.nCols();j++) {
                out.m[i][j] = A.m[i][j]*D.m[i][i];
            }
        }
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
    
    T innerProduct(const vec &a, const vec &b) const {
        T accum=0;
        if(a.size()!=b.size()) {fprintf(stderr,"Error: inner product vectors not same length %zu,%zu)",a.size(),b.size()); return 0;}
        for(int i=0;i<a.size();i++) accum+=(a[i]*b[i]);
        return accum;
    }
    
    friend vec matVec (const Matrix &A,const vec &v) {
        vec out(v.size()); //vec out;
        if(A.nCols()!=v.size()){fprintf(stderr,"Error MatVec is mismatched A[%u,%u] v[%zu]\n",A.nRows(),A.nCols(),v.size()); return out;}
        for(int i=0;i<A.nRows();i++){
            out[i] = 0;
            for(int j=0;j<A.nCols();j++){
                out[i]+=A.m[i][j]*v[j]; //changed out[j] to out[i]
            }
        }
        return out;
    }
    
    vec vectorCombination(T a, const vec &U,T b, const vec &V) const { //was this->minally friend. 
        vec W(U.size());
        if(U.size() != V.size()) {fprintf(stderr,"vectorCombo: sizes don't match %zu,%zu\n",U.size(),V.size()); return W;}
        for(int i=0;i<W.size();i++) W[i] = a * U[i] + b * V[i];
        return W;
    }
    
    T vectorNorm(const vec & V) const {
        return sqrt( innerProduct(V,V) ); 
    }

    T vectorNormOne(const vec &v) const {
        T sum=0;
        for(int i=0;i<v.size();i++) sum += abs(v[i]);
        return sum;
    }

    T vectorNormInf(const vec &v) const {
        return abs(*max_element(v.begin(), v.end(), [](T a, T b) {return abs(a) < abs(b);}));
    }
    
    // only checks to see if the matrix is square (not Symmetric in the mathematical sense)
    bool isSquare() const {
        if(nRows()==nCols()) return 1;
        else return 0;
    }
    
    bool isSymmetric() const {
        if(!this->isSquare()) return 0;
        for(int i=0;i<this->nRows();i++){
            for(int j=0;j<this->nCols();j++){
                if(this->m[i][j]!=this->m[j][i]) return 0;
            }
        }
        return 1;
    }

    int numAsymmetricEntries() const {
        int nonSymmetrics = 0;
        if(!this->isSquare()) return -1;
        for(int i=0;i<this->nRows();i++){
            for(int j=0;j<this->nCols();j++){
                if(this->m[i][j]!=this->m[j][i]) {nonSymmetrics+=1; cout << m[i][j] << ", " << m[j][i] << endl;}
            }
        }
        return nonSymmetrics;
    }

    bool isSymPD() const {
        if (!isSymmetric()) return 0;
        size_t n = nCols();
        Matrix<T> L(n, n), U(n, n);
        this->LUdecomposition(L, U);

        for (int i=0;i<nCols();i++)
            if (U.m[i][i] < 0) return 0;
        return 1;
    }
    
    bool isUpperTriangular() const {
        if(!this->isSquare()) return 0;
        if(nRows()<=1) return 0;
        for(int i=0;i<nRows();i++)
            for(int j=0;j<i;j++)
                if(m[i][j]>0) return 0;
        return 1;
    }

    //Returns the smallest entry in the matrix
    T getMin() {
        vec mins(nRows()); 
        for (vec row : m) mins.push_back(*min_element(row.begin(), row.end(), [](T a, T b) { return abs(a) < abs(b); }));
        return *min_element(mins.begin(), mins.end(), [](T a, T b) { return abs(a) < abs(b); });
    }

    //Returns the largest entry in the matrix
    T getMax() {
        vec maxs(nRows()); 
        for (vec row : m) maxs.push_back(*max_element(row.begin(), row.end(), [](T a, T b) { return abs(a) < abs(b); }));
        return *max_element(maxs.begin(), maxs.end(), [](T a, T b) { return abs(a) < abs(b); });
    }

    friend Matrix add(const Matrix &A, const Matrix &B) {
        int r = A.nRows();
        int c = A.nCols();
        if (B.nRows() != r || B.nCols() != c) { fprintf(stderr, "Invalid Dimensions."); exit(1); }

        Matrix R(r,c);
            
        for (int i=0;i<r;i++) {
            for (int j=0;j<c;j++) {
                R.m[i][j] = A.m[i][j] + B.m[i][j];
            }
        }

        return R;
    }

    friend Matrix scale(T k, Matrix A) {
        Matrix R(A.nRows(),A.nCols());
        for (int i=0;i<A.nRows();i++) {
            for  (int j=0;j<A.nCols();j++) {
                R.m[i][j] = k*A.m[i][j];
            }
        }
        return R;
    }

    vec scaleVec(T k, vec x) const {
        vec out(x.size());
        for (int i=0;i<x.size();i++) out[i]=k*x[i];
        return out;
    }

    T frobenius() {
        T total=0;
        for (int i=0;i<nRows();i++) {
            for (int j=0;j<nCols();j++) {
                total += m[i][j]*m[i][j];
            }
        }
        total = sqrt(total);
        return total;
    }

    T infinityNorm() {
        T high=0;
        for (int i=0;i<nRows();i++) {
            T sum=0;
            for (int j=0;j<nCols();j++) sum += abs(m[i][j]);
            high = max(high, sum);
        }
        return high;
    }
    
    bool solveUpperTriangularSystem(vec &x,const vec b) const { 
        if(!this->isUpperTriangular()) return false;
        if (!(x.size() == nRows())) return false; 
        // solve sequentially by back-substitution from bottom to top
        for(int i=nRows()-1;i>=0;i--){
            T sum = 0;
            int j;
            for(j=i+1;j<nRows();j++) sum+=m[i][j]*x[j]; //Changed j to start at i+1 from i. 
            x[i]=(b[i]-sum)/m[i][i];
        }
        return true;
    }
    
    bool solveUpperTriangularSystem(const Matrix U,vec &x,const vec b) const {
        if(!U.isUpperTriangular()) return false;
        // solve sequentially by back-substitution from bottom to top
        for(int i=0;i<U.nRows();i++){
            T sum = 0;
            int j;
            for(j=i+1;j<U.nRows();j++) sum+=U.m[i][j]*x[j]; //Changed j to start at i+1 from i.  
            x[i]=(b[i]-sum)/U.m[i][i];
        }
    }
    
    bool isLowerTriangular() const {
        if(!this->isSquare()) return 0;
        if(nRows()<1) return 0;
        for(int i=0;i<nRows();i++)
            for(int j=i+1;j<nCols();j++)
                if(m[i][j]!=0) return 0;
        return 1;
    }
    
    bool solveLowerTriangularSystem(vec &x,const vec &b) const {
        if(!this->isLowerTriangular()) return false;
        if (!(x.size() == nRows())) return false;
        // solve sequentially by forward-substitution from top to bottom
        
        for(int i=0;i<this->nRows();i++){
            T sum = 0;
            int j;
            for(j=0;j<i;j++) sum += m[i][j] * x[j]; //changed loop bound from i-1 to i
            x[i] = (b[i]-sum) / m[i][i];

        }
        return true; // success
    }
    
    bool solveLowerTriangularSystem(const Matrix L, vec &x,const vec b) const {
        if(!L.isLowerTrangular()) return false;
        if (!(x.size() == nRows())) return false;
        // solve sequentially by forward-substitution from top to bottom
        for(int i=0;i<L.nRows();i++){
            T sum = 0;
            int j;
            for(j=0;j<i;j++) {         //changed loop bound from i-1 to i
                sum+=L.m[i][j] * x[j]; // check for dependencies here.
            }
            x[i] = (b[i]-sum) / L.m[i][i];
        }
        return true; // success
    }

    //Generalization of above method
    /*
    bool solveBlockLowerTriangularSystem(const Matrix<Matrix> L, vec &x,const vec b) const {
        if(!L.isLowerTrangular()) return false;
        if (!(x.size() == nRows())) return false;
        fill(x.begin(), x.end(), 0);
        
        trueIndex = 0;
        
        for(int i=0;i<L.nRows();i++){
            
            vector<T> subX(L[i][0].nCols());
            vector<T> subB = getSubVector(b, trueIndex, trueIndex+L[i][0].nCols());
            vector<T> sum(L[i][0].nCols());
            int j;
            for(j=0;j<i;j++) {         //changed loop bound from i-1 to i
                vector<T> product = matVec(L.m[i][j], x[j]); // check for dependencies here.
                transform (sum.begin(), sum.end(), product.begin(), product.end(), plus<T>());
            }
            
            transform (b.begin(), b.end(), sum.begin(), sum.end(), minus<T>());

            L.m[i][i].solveLowerTriangularSystem(subX, subB);

            for (int i=trueIndex;i<trueIndex+L[i][0].nCols();i++) x[i]=subX[i-trueIndex];

            trueIndex+=L[i][0].nCols();
        }
        return true; // success
    }
    */

    /*
    Solve system. Iterations is an optional variable that specifies number of iterations past the initial solve
    that we want to perform. (iterative refinement)
    */
    bool triSolve(vec &x, const vec b, int iterations=1) const {
        int n = this->nRows();

        if (!(isSquare())) { fprintf(stderr, "Matrix not square"); return false; }
        if (x.size() != n || b.size() != n) { fprintf(stderr, "Invalid dimensions"); return false; }
        
        for (int i=0;i<n;i++) x[i]=0;
        
        Matrix<T> L(n, n), U(n, n);
        LUdecomposition(L, U);
        vec r = b;
        vec y(n), d(n);
        
        for (int i=0;i<iterations;i++) {
            L.solveLowerTriangularSystem(y, r);
            U.solveUpperTriangularSystem(d, y);
            x = vectorCombination(1, x, 1, d);

            r = vectorCombination(1, b, 1, matVec(*this, x));
        }
        
        return true; 
    }

    void symmetricTriSolve(vec &x, const vec b) {
        int n = this->nRows();

        if (!isSymmetric()) { fprintf(stderr, "Matrix not symmetric"); return; }
        if (x.size() != n || b.size() != n) { fprintf(stderr, "Invalid dimensions"); return; }
        
        for (int i=0;i<n;i++) x[i]=0;
        
        Matrix<T> R(n, n), Rstar(n,n);

        cholesky(R);
        Rstar = R;
        Rstar.transpose();
        
        vec y(n);
        Rstar.solveLowerTriangularSystem(y, b);
        R.solveUpperTriangularSystem(x, y);
    }

    //solve Ax=b with A=LU done in 'E' precision. Use iterative refinement to reach tolerance threshold.
    //Tolerance metric is taken from Higham.  
    template<class E>
    int triSolveMixed(vec &x, const vec b, double tolerance) {
        int n = this->nRows();

        if (!isSquare()) { fprintf(stderr, "Matrix not square"); return -1; }
        if (x.size() != n || b.size() != n) { fprintf(stderr, "Invalid dimensions"); return -1;   }
        
        for (int i=0;i<n;i++) x[i]=0;
        
        Matrix<E> LH(n,n), UH(n,n);
        Matrix<E> A(*this);
        A.LUdecomposition(LH, UH);
        
        Matrix<T> L(LH), U(UH);
        
        vec y(n), r(n), d(n);
        r = vectorCombination(1, b, -1, matVec(*this, x));
        
        T error = vectorNormInf(r)/(infinityNorm()*vectorNormInf(x)+vectorNormInf(b)); //from Higham.
        
        int k=0;
        while (error > tolerance && k < 1000) {
            L.solveLowerTriangularSystem(y, r);
            U.solveUpperTriangularSystem(d, y);
            x = vectorCombination(1, x, 1, d);
            r = vectorCombination(1, b, -1, matVec(*this, x));

            error = vectorNormInf(r)/(infinityNorm()*vectorNormInf(x)+vectorNormInf(b));
            k++;
        }

        return k;
    }

    //Use three precisions. From Carson/Higham.
    //E is factorization/correction-equation solving precision. R is for calculating residual
    //T is working precision. 
    template<class E, class R>
    int triSolveThreePrecision(vec &x, const vec b, double tolerance) {
        int n = this->nRows();
        Matrix<R> B;

        if (!isSquare()) { fprintf(stderr, "Matrix not square"); return -1; }
        if (x.size() != n || b.size() != n) { fprintf(stderr, "Invalid dimensions"); return -1;   }
        
        for (int i=0;i<n;i++) x[i]=0;
        
        Matrix<E> L(n,n), U(n,n); 
        Matrix<E> A(*this);
        A.LUdecomposition(L, U);
        
        vector<E> y(n), d(n), rs(n);
        vector<R> r(n), br(n), axr(n);
        vec dw(n);
        
        convert(br, b);
        convert(axr, matVec(*this, x));
        
        r = B.vectorCombination(1, br, -1, axr);
        T error = B.vectorNormInf(r)/(infinityNorm()*vectorNormInf(x)+vectorNormInf(b)); //from Carson/Higham
        
        int k=0;
        while (error>tolerance && k<500) {
            convert(rs, r); //convert residual to factorization/correction-equation precision (probably half)
            L.solveLowerTriangularSystem(y, rs);
            U.solveUpperTriangularSystem(d, y);
            convert(dw, d); //convert correction to working precision. 

            x = vectorCombination(1, x, 1, dw); //update solution. 

            convert(axr, matVec(*this, x));
            r = B.vectorCombination(1, br, -1, axr); //update residual in (presumably) higher precision. 

            error = B.vectorNormInf(r)/(infinityNorm()*vectorNormInf(x)+vectorNormInf(b)); //update error
            k++;
        }

        return k;
    }

    //mixed precision for symmetric PD matrix.
    //file input is to capture backward error of factorization in the case that an error occurs when 
    //computing the solution. 
    template<class E>
    pair<int, T> symmetricTriSolveMixed(vec &x, const vec b, mpf_class tolerance, string filename="") {
        int n = this->nRows();

        if (x.size() != n || b.size() != n) { fprintf(stderr, "Invalid dimensions"); exit(1); }
        for (int i=0;i<n;i++) x[i]=0;
        
        Matrix<E> RH(n,n);
        Matrix<E> A(*this);
        A.cholesky(RH);

        Matrix<T> R(RH);
        Matrix<T> Rstar(R);  
        Rstar.transpose();

        T backwardError = (add(multiply(Rstar, R), scale(-1, *this)).frobenius())/frobenius();
        if (!filename.empty()) {
            ofstream f(filename, ofstream::app);
            f << "\tR*R backward error: " << backwardError << endl; //Useful information in case of failure.
            f.close();
        }
        
        vec y(n), d(n), r(n);
        r=b;
        T error = vectorNormInf(r)/(infinityNorm()*vectorNormInf(x)+vectorNormInf(b));

        cout << tolerance << endl;
        int k=0;
        while (error > tolerance && k<1000) {
            Rstar.solveLowerTriangularSystem(y, r);
            R.solveUpperTriangularSystem(d, y);
            x = vectorCombination(1, x, 1, d);
            r = vectorCombination(1, b, -1, matVec(*this, x));

            error = vectorNormInf(r)/(infinityNorm()*vectorNormInf(x)+vectorNormInf(b));
            k++;
        }

        return pair<int, T>(k, backwardError);
    }

    template<class E, class G>
    int symTriSolveThreePrecision(vec &x, const vec b, double tolerance) {
        int n = this->nRows();
        Matrix<G> B; //For accessing friend members !

        if (x.size() != n || b.size() != n) { fprintf(stderr, "Invalid dimensions"); return -1;   }
        
        for (int i=0;i<n;i++) x[i]=0;
        
        Matrix<E> R(n,n), Rstar(n,n);
        Matrix<E> A(*this);
        A.cholesky(R);
        Rstar=R;
        Rstar.transpose();
        
        vector<E> y(n), d(n), rs(n);
        vector<G> r(n), br(n), axr(n);
        vec dw(n);
        
        convert(br, b);
        convert(axr, matVec(*this, x));
        
        r = B.vectorCombination(1, br, -1, axr);
        T error = B.vectorNormInf(r)/(infinityNorm()*vectorNormInf(x)+vectorNormInf(b));
        
        int k=0;
        while (error>tolerance && k<500) {
            convert(rs, r);
            Rstar.solveLowerTriangularSystem(y, rs);
            R.solveUpperTriangularSystem(d, y);
            convert(dw, d);

            x = vectorCombination(1, x, 1, dw);

            convert(axr, matVec(*this, x));
            r = B.vectorCombination(1, br, -1, axr);

            error = B.vectorNormInf(r)/(infinityNorm()*vectorNormInf(x)+vectorNormInf(b));
            k++;
        }

        return k;
    }
    
    // void findIndependentOps(){ findIndependentOps(*this); }
    
    void findIndependentOps(Matrix<int> &backDeps){
        // or should we just do the TriSolve and compute costs for sending the independent components
        // this will be word granularity
        
        // how many sends required (and output the message map for that)
        // also count how many pre-sends would be required?
        
        // computing row depends (Fig 1 of Liu paper)
        // for each row, for each non-zero in row add back-depend to row based on column the nonzero is in
        // then you can sort them into level-sets
        // can fill out table of forward depends (make depends bi-directional) so that you can exec sync free
        // back-depends need to be converted into forward dependencies so that you can message targets (free row x)
        // can round-robin the rows?  (or blocks if done in block-mode)
        
        // any row with just a non-zero on the diagonal is independent (within that row then what?)
        // need to follow the CSR for the row compression
        
        // can pack into a CSR form using a row pointer array and a dependency array
        // or row pointers and then do a brute force search for the back-dependencies....
        // or do a T(r,c) for each point in the matrix?
        // can we brute-force compute the dependencies using the x[j] and m[i][j] and
        //     x[i] = (b[i]) and then the accumulated dependencies that are accumulated into the "sum"
        // the encoding of the matrices (CSR) should be row-oriented, so partitioning can be row-oriented
        // and then round-robin (modulo) to processor cores....
        //      we can repeat this analysis for block-oriented structures and do BLAS3-like operations
        
        // each row can be a vec (worse case depend is full matrix that is populated incrementally)
        
        Matrix<int> forwardDeps;
        vector<int> ndepsBack,ndepsForward,inDegree,inDegreeDone;  // for forward deps
        int maxpar=0; // stats
        int nlevels=0;
        // printf("Prepare to resize to %u\n",this->nRows());
        forwardDeps.setSize(this->nRows(),this->nCols());
        backDeps.setSize(this->nRows(),this->nCols());
        ndepsBack.resize(this->nRows());  ndepsForward.resize(this->nRows());  inDegree.resize(this->nRows()); inDegreeDone.resize(this->nRows());
        for(int i=0;i<this->nRows();i++) { ndepsForward[i]=ndepsBack[i]=inDegree[i]=0; }// zero out the deps array
        // first lets compute the ackward dependencies
        for(int r=0;r<this->nRows();r++){
            //fprintf(stderr,"for Row[%u], compute nDepsBack\n",r);
            ndepsBack[r]=0; // for current row to compute the back depends... accumulate the forward depends by walking through backDeps mat.
            // for each row, for each nonzero in row (up until diagonal) lets compute backDepends()
            for(int c=0;c<(r-1);c++) {
                //fprintf(stderr,"c[%u]=>%lf\t",c,this->m[r][c]);
                if(this->m[r][c]!=0) {
                    //fprintf(stderr,"\n\t Location[%u][%u]\n",r,c);
                    backDeps.m[r][ndepsBack[r]]=c+1; (ndepsBack[r])++;
                }
            }
            // backDeps.m[r][colptr]=-1;
        } // backdeps complete
        // now convert into forward dependencies
        for(int r=0;r<this->nRows();r++){
            // walk through the entries of the forward dependency matrix
            for(int c=0;c<ndepsBack[r];c++){
                // search for forward dependency and accumulate into the forward deps matrix
                // so my backward dependency should be recorded into the forward Deps matrix...
                int dep = backDeps.m[r][c]-1; // record the backward dep
                // now the forward dep will be forwardDeps[dependentRow][current number of deps for that row]
                forwardDeps.m[dep][ndepsForward[dep]] = r;  // set that location to point to the row that had the backDep
                (ndepsForward[dep])++; // increment the number of deps for that dependent row that was just updated
            }
        }
        printf("PrintBackwardDeps********************\n");
        backDeps.print();
        printf("PrintForward Deps********************\n");
        forwardDeps.print();
        // can I use the backdeps to count # of independent sets?
        // first compute indegree for each row (indegree =0 can fire)
        printf("\n NdepsBack--->");
        for(int r=0;r<this->nRows();r++){
            printf("[%u].",ndepsBack[r]);
            inDegree[r]=ndepsBack[r]; // just copy because we will destroy the info in ndeps in subsequent steps
            inDegreeDone[r]=0;
        }
        puts("\n");
        // use graph coloring?  start from first row?
        //  How about do a brute-force search for tasks that have deps satisfied for each forward step
        //  So can do "step 1" fire forward dependencies from first row and see how many forward deps there are
        //  But need to track vertex in-degree to know it is ready to fire
        //  so compute in-degree for each row to count the deps....  (or do back-dep to compute firing order?)
        for(int i=0;i<this->nRows();i++){
            int npar=0,notdone;
            printf("Firing for Level[%u]  .",i);
            for(int r=0;r<this->nRows();r++){
                // scan for rows with no deps
                if(inDegreeDone[r]==0 && inDegree[r]==0) {
                    printf("%03u.",r);
                    npar++;
                    if(npar>maxpar) maxpar=npar; // increment max parallelism of levelset (could compute average too)
                    inDegreeDone[r]=-1;
                }
            }
            puts(""); // end that levelset of independent rows for this level
            // for(int r=0;r<this->nRows();r++){ printf("\tLevelset[%d]: inDegree[%d]=%d done=%d\n",i,r,inDegree[r],inDegreeDone[r]);}
            for(int r=0;r<this->nRows();r++){
                // now we release forward deps for those rows
                if(inDegreeDone[r]<0 && inDegree[r]==0){
                    for(int c=0;c<ndepsForward[r];c++){
                        // release the deps
                        // printf("Release Forward Dep[%d][%d] value=%d\n",r,c,forwardDeps.m[r][c]);
                        (inDegree[(forwardDeps.m[r][c])])--; // decrement deps
                        // printf("\tNow inDegree[%d]==>%d\n",forwardDeps.m[r][c],inDegree[forwardDeps.m[r][c]]);
                    }
                }
                // printf(".....Updated Row[%u] inDegree=%d done=%d\n",r,inDegree[r],inDegreeDone[r]);
            }
            
            // Now fix up the "done" flags
            // printf("fixing up notdone:**********************\n\n");
            notdone=0;
            for(int r=0;r<this->nRows();r++) {
                //printf("r[%u] inDegree[%u] notDone[%d]\n",r,inDegree[r],inDegreeDone[r]);
                if(inDegreeDone[r]<0){
                    inDegreeDone[r]=i+1; // color by levelset
                }
                if(inDegreeDone[r]==0) notdone=1;
                //printf("\t....r[%u] inDegree[%u] notDone[%d]\n",r,inDegree[r],inDegreeDone[r]);
            } // mark to ignore for next cycle
            if(notdone==0) { nlevels=i+1; break; }  // exit loop because no more work to do
            // inDegreeDone[r]=1; // mark that row as completed and don't count it any more
            //for(int r=0;r<this->nRows();r++){ printf("\tRescan Levelset[%d]: inDegree[%d]=%d\n",i,r,inDegree[r]);}
            // scan quickly for level
            // notdone=0;
            // for(int r=0;r<U.nRows();r++) { if(inDegree[r]==0) notdone=1; }
            // if(notdone==0) break; // break out of the loop... we are done
            // this should print out the level sets for each row of the matrix
        }
        // how many parallel ops are available?
        printf("\n\n\tFinal Tally... maxlevels=%u nlevels=%u \t maxpar=%u\n",this->nRows(),nlevels,maxpar);
        // or output in DOT format?  (AT&T Graph Viz... what is it called for Brew?
        // how many levels in a levelset are required?
        
        printf("now dump in graphviz format\n*************\n\n");
        printf("use \n\tneato -Tps fdeps.gv -o fdeps.ps\n\t or dot\n");
        printf("Digraph ForwardDeps {\n");
        printf("\tsize=\"7.5,10\";\n \tratio = compress;\n");
        // node[style=filled, colorscheme=bugn9, color=7];
        for(int r=0;r<this->nRows();r++){
            // now we release forward deps for those rows
            // if(ndepsForward[r]==0){ printf("\trow%u -> row%u;\n",r,r);}
            for(int c=0;c<ndepsForward[r];c++){
                printf("\t\t node[style=filled, colorscheme=pastel13, color=%u];\n",inDegreeDone[forwardDeps.m[r][c]]);
                printf("\trow%u -> row%u;\n",r,forwardDeps.m[r][c]);
            }
        }
        puts("}");
        
        printf("Digraph BackDeps {\n");
        for(int r=0;r<this->nRows();r++){
            // now we release forward deps for those rows
            for(int c=0;c<ndepsBack[r];c++){
                printf("\trow%u -> row%u;\n",r,backDeps.m[r][c]-1);
            }
        }
        puts("}");
        
        // could create a proxy matrix with 1 for blocks that are nonzero and 0 for zero blocks
        // can recompute the dependencies for that
    }

    
    //Input blocksize
    template<int N>
    void findIndependentOps(Matrix<uint16_t> &backDeps){
        // or should we just do the TriSolve and compute costs for sending the independent components
        // this will be word granularity
        
        // how many sends required (and output the message map for that)
        // also count how many pre-sends would be required?
        
        // computing row depends (Fig 1 of Liu paper)
        // for each row, for each non-zero in row add back-depend to row based on column the nonzero is in
        // then you can sort them into level-sets
        // can fill out table of forward depends (make depends bi-directional) so that you can exec sync free
        // back-depends need to be converted into forward dependencies so that you can message targets (free row x)
        // can round-robin the rows?  (or blocks if done in block-mode)
        
        // any row with just a non-zero on the diagonal is independent (within that row then what?)
        // need to follow the CSR for the row compression
        
        // can pack into a CSR form using a row pointer array and a dependency array
        // or row pointers and then do a brute force search for the back-dependencies....
        // or do a T(r,c) for each point in the matrix?
        // can we brute-force compute the dependencies using the x[j] and m[i][j] and
        //     x[i] = (b[i]) and then the accumulated dependencies that are accumulated into the "sum"
        // the encoding of the matrices (CSR) should be row-oriented, so partitioning can be row-oriented
        // and then round-robin (modulo) to processor cores....
        //      we can repeat this analysis for block-oriented structures and do BLAS3-like operations
        
        // each row can be a vec (worse case depend is full matrix that is populated incrementally)
        //for (int i=0;i<N-(nCols()%N);i++) m.push_back(vec(nCols()));
        int numBlocks=((int) (nRows()/N));
        Matrix<uint16_t> forwardDeps;
        vector<uint16_t> ndepsBack,ndepsForward,inDegree,inDegreeDone;  // for forward deps
        int maxpar=0; // stats
        int nlevels=0;
        // printf("Prepare to resize to %u\n",this->nRows());
        forwardDeps.setSize(numBlocks,numBlocks);
        backDeps.setSize(numBlocks,numBlocks);
        ndepsBack.resize(numBlocks);  ndepsForward.resize(numBlocks);  inDegree.resize(numBlocks); inDegreeDone.resize(numBlocks);
        for(int i=0;i<numBlocks;i++) { ndepsForward[i]=ndepsBack[i]=inDegree[i]=0; }// zero out the deps array
        // first lets compute the ackward dependencies
        for(int r=0;r<numBlocks;r++){
            //fprintf(stderr,"for Row[%u], compute nDepsBack\n",r);
            ndepsBack[r]=0; // for current row to compute the back depends... accumulate the forward depends by walking through backDeps mat.
            // for each row, for each nonzero in row (up until diagonal) lets compute backDepends()
            for(int c=0;c<r;c++) { //changed from (r-1)
                //fprintf(stderr,"c[%u]=>%lf\t",c,this->m[r][c]);
                bool nonzero=0;
                for (int j=0;j<N;j++) {
                    for (int k=0;k<N;k++) {
                        if(this->m[r*N+j][c*N+k]!=0) nonzero=1;
                    }
                }
                if (nonzero) {
                    backDeps.m[r][ndepsBack[r]]=c+1; 
                    (ndepsBack[r])++;
                }
            }
        // backDeps.m[r][colptr]=-1;
        } // backdeps complete
        // now convert into forward dependencies
        for(int r=0;r<numBlocks;r++){
            // walk through the entries of the forward dependency matrix
            for(int c=0;c<ndepsBack[r];c++){
                // search for forward dependency and accumulate into the forward deps matrix
                // so my backward dependency should be recorded into the forward Deps matrix...
                int dep = backDeps.m[r][c]-1; // record the backward dep
                // now the forward dep will be forwardDeps[dependentRow][current number of deps for that row]
                forwardDeps.m[dep][ndepsForward[dep]] = r;  // set that location to point to the row that had the backDep
                (ndepsForward[dep])++; // increment the number of deps for that dependent row that was just updated
            }
        }
        printf("PrintBackwardDeps********************\n");
        //backDeps.print();
        printf("PrintForward Deps********************\n");
        //forwardDeps.print();

        // can I use the backdeps to count # of independent sets?
        // first compute indegree for each row (indegree =0 can fire)
        printf("\n NdepsBack--->");
        for(int r=0;r<numBlocks;r++){
            printf("[%u].",ndepsBack[r]);
            inDegree[r]=ndepsBack[r]; // just copy because we will destroy the info in ndeps in subsequent steps
            inDegreeDone[r]=0;
        }
        puts("\n");
        // use graph coloring?  start from first row?
        //  How about do a brute-force search for tasks that have deps satisfied for each forward step
        //  So can do "step 1" fire forward dependencies from first row and see how many forward deps there are
        //  But need to track vertex in-degree to know it is ready to fire
        //  so compute in-degree for each row to count the deps....  (or do back-dep to compute firing order?)
        for(int i=0;i<numBlocks;i++){
            int npar=0,notdone;
            printf("Firing for Level[%u]  .",i);
            for(int r=0;r<numBlocks;r++){
                // scan for rows with no deps
            
                if(inDegreeDone[r]==0 && inDegree[r]==0) {
                    printf("%03u.",r);
                    npar++;
                    if(npar>maxpar) maxpar=npar; // increment max parallelism of levelset (could compute average too)
                    inDegreeDone[r]=-1;
                }
            }
            puts(""); // end that levelset of independent rows for this level
            // for(int r=0;r<this->nRows();r++){ printf("\tLevelset[%d]: inDegree[%d]=%d done=%d\n",i,r,inDegree[r],inDegreeDone[r]);}
            for(int r=0;r<numBlocks;r++){
                // now we release forward deps for those rows
                if(inDegreeDone[r]<0 && inDegree[r]==0){
                    for(int c=0;c<ndepsForward[r];c++){
                        // release the deps
                        // printf("Release Forward Dep[%d][%d] value=%d\n",r,c,forwardDeps.m[r][c]);
                        (inDegree[(forwardDeps.m[r][c])])--; // decrement deps
                        // printf("\tNow inDegree[%d]==>%d\n",forwardDeps.m[r][c],inDegree[forwardDeps.m[r][c]]);
                    }
                }
                // printf(".....Updated Row[%u] inDegree=%d done=%d\n",r,inDegree[r],inDegreeDone[r]);
            }
            
            // Now fix up the "done" flags
            // printf("fixing up notdone:**********************\n\n");
            notdone=0;
            for(int r=0;r<numBlocks;r++) {
                //printf("r[%u] inDegree[%u] notDone[%d]\n",r,inDegree[r],inDegreeDone[r]);
                if(inDegreeDone[r]<0){
                    inDegreeDone[r]=i+1; // color by levelset
                }
                if(inDegreeDone[r]==0) notdone=1;
                //printf("\t....r[%u] inDegree[%u] notDone[%d]\n",r,inDegree[r],inDegreeDone[r]);
            } // mark to ignore for next cycle
            if(notdone==0) { nlevels=i+1; break; }  // exit loop because no more work to do
            // inDegreeDone[r]=1; // mark that row as completed and don't count it any more
            //for(int r=0;r<this->nRows();r++){ printf("\tRescan Levelset[%d]: inDegree[%d]=%d\n",i,r,inDegree[r]);}
            // scan quickly for level
            // notdone=0;
            // for(int r=0;r<U.nRows();r++) { if(inDegree[r]==0) notdone=1; }
            // if(notdone==0) break; // break out of the loop... we are done
            // this should print out the level sets for each row of the matrix
        }
        
        // how many parallel ops are available?
        printf("\n\n\tFinal Tally... maxlevels=%u nlevels=%u \t maxpar=%u\n",this->nRows(),nlevels,maxpar);
        // or output in DOT format?  (AT&T Graph Viz... what is it called for Brew?
        // how many levels in a levelset are required?
        
        //printf("now dump in graphviz format\n*************\n\n");
        //printf("use \n\tneato -Tps fdeps.gv -o fdeps.ps\n\t or dot\n");
        //printf("Digraph ForwardDeps {\n");
        //printf("\tsize=\"7.5,10\";\n \tratio = compress;\n");
        printf("now dump in graphviz format\n*************\n\n");
        printf("use \n\tneato -Tps fdeps.gv -o fdeps.ps\n\t or dot\n");
        
        FILE *file;
        file = fopen("graph.gv", "w");
        //if ((file = fopen("fdeps.txt", "w")) == NULL) {printf("asdf"); return;}
        
        fprintf(file, "Digraph ForwardDeps {\n");
        fprintf(file, "\tsize=\"2.5,5\";\n \tratio = compress;\n");
        // node[style=filled, colorscheme=bugn9, color=7];
        for(int r=0;r<numBlocks;r++){
        //    // now we release forward deps for those rows
        //    // if(ndepsForward[r]==0){ printf("\trow%u -> row%u;\n",r,r);}
            
            for(int c=0;c<ndepsForward[r];c++){
                fprintf(file,  "\t\t node[style=filled, colorscheme=set312, color=%u];\n",inDegreeDone[forwardDeps.m[r][c]] % 12 + 1);
                fprintf(file,  "\trow%u -> row%u;\n",r,forwardDeps.m[r][c]);
            }
        }
        ////puts("}");
        fprintf(file, "}");
        //
        //fprintf(file, "Digraph BackDeps {\n");
        //for(int r=0;r<numBlocks;r++){
        ////    // now we release forward deps for those rows
        //    for(int c=0;c<ndepsBack[r];c++){
        //        //fprintf(file, "\trow%u -> row%u;\n",r,backDeps.m[r][c]);
        //    }
        //}
        //////puts("}");
        //fprintf(file, "}");
        if (file !=stdin) fclose(file);
        printVector(inDegree);
        

        // could create a proxy matrix with 1 for blocks that are nonzero and 0 for zero blocks
        // can recompute the dependencies for that
    }
    
    void propagateDependencies(const Matrix U){
        // work dependencies backwards from X[row];
        // x[row+1] depends on prior rows
        // if compressed with CSR, then could just hand out row-pointers to the processors
        // and walk through row pointers to create the table
        // count up operations required for each (account for pipeline depth)
        // could gang-schedule stuff if at least 6 cycles for pipeline depth...
    }
    
    
    int findNumBlocks(const Matrix U,int rblock, int cblock){ // for given row and column blocking, how many blocks are there
        // next compute the dependencies between those blocks  (and then can compute how many cycles to
        // similar to previous one
    }
    // next compute what kind of compression can be achieved with delta-compression
    // what is the overall largest delta for a RLE or CSR encoding for the matrix
    // for the BCSR encoding of the matrix?
    // at what point is hypersparse the best encoding (can we do delta encoding for hypersparse?
    // This function would return how much space the encoding would take
    
    // performance density for FFT
    // what is FFT performance density and area density?

    //Matrix must be symmetric PD.
    
    void cholesky(Matrix &R) {
        if (!this->isSymmetric()) { fprintf(stderr, "Matrix should be symmetric and positive definite."); return; }
        R.setSize(this->nRows(), this->nRows());
        R.setZero();
        int n = R.nRows();

        for (int i = 0; i < n; i++) { 
            for (int j = 0; j <= i; j++) { 
                T sum = 0; 
    
                if (j == i) // summation for diagnols 
                { 
                    for (int k = 0; k < j; k++) 
                        sum += R.m[j][k]*R.m[j][k];
                    if ((m[j][j] - sum) < 0) { fprintf(stderr, "Matrix not Positive Definite.\n"); exit(1); }
                    R.m[j][j] = sqrt(m[j][j] - sum); 
                } else { 
                    // Evaluating L(i, j) using L(j, j) 
                    for (int k = 0; k < j; k++) 
                        sum += (R.m[i][k] * R.m[j][k]); 
                    R.m[i][j] = (m[i][j] - sum) / R.m[j][j]; 
                } 
            } 
        }
        R.transpose(); 
    }
    
    bool LUdecomposition(Matrix &L, Matrix &U) const {
        L.setSize(m.size(), m.size());
        U.setSize(m.size(), m.size());
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
                    sum += (L.m[i][j] * U.m[j][k]); 
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
                    if(U.m[i][i]==0) { fprintf(stderr, "divide by zero error."); return false; } // fail... cannot divide by zero
                    L.m[k][i] = (m[k][i] - sum) / U.m[i][i];
                }
            }
        }
        return true;
    }

    bool LUdecomposition(Matrix<Matrix> blocks) const {
        ///
    }

    
    // now we need to overload the greater than operator or the less than operator
    T max(T a,T b){
        if(a>b) return a;
        else return b;
    }

    //easy row, column normalization. Doesn't preserve symmetry. From Higham "Squeezing matrix..." paper
    //normalize according to infinity norm so max entry in each row/column is 1. 
    Matrix rescale(int type) {
        int m=nRows();
        int n=nCols();

        Matrix<T> D(m,n), S(m,n), A(m,n);
        for (int i=0;i<m;i++) D.m[i][i]=(1/vectorNormInf(this->m[i]));
        A = multiplyDiagLeft(D, *this);
        
        A.transpose(); //easy access to columns
        for (int i=0;i<m;i++) S.m[i][i]=(1/vectorNormInf(A.m[i]));
        A.transpose(); //flip back.

        A = multiplyDiagRight(A, S);

        T u;
        if (type==0 || type==1) u=4096; //take better advantage of half range. 
        if (type==2)            u=pow(2, pow(2, ES)); //try to entries into good posit interval (this is useed). 
        
        A = scale(u, A);

        return A;
    }

    //uses knight/ruiz approach for equilibrating matrix according to inf-norm. 
    //similar to above but preserves symmetry. 
    Matrix rescaleSymmetric(int type) {
        if (!this->isSymmetric()) { fprintf(stderr, "Matrix should be symmetric."); exit(1); }

        double theta = .1;
        double tolerance=0;

        int n = nCols();
        Matrix A(*this);

        Matrix r(n,n);
        while (abs(r.getMax()-1) > tolerance) {

            for (int i=0;i<A.nCols();i++) {
                r.m[i][i] = 1/sqrt(vectorNormInf(A.m[i]));
            }

            A=multiplyDiagLeft(r, A);
            A=multiplyDiagRight(A, r);
        }
        
        T u=1;
        if (type ==0 || type==1) {
            u=4096;
        }
        if (type==2) {
            u=pow(2, pow(2, ES));
        } 
        
        A = scale(u, A);
        return A;
    }

    //This uses knight, ruiz scaling applied to the 1-norm \
    //Convergence isn't guaranteed as it is when applied to inf-norm.
    Matrix rescaleSymmetricOne(int type) {
        if (!this->isSymmetric()) { fprintf(stderr, "Matrix should be symmetric."); exit(1); }
        double tolerance=0; 

        int n = nCols();
        Matrix A(*this);

        Matrix r(n,n);
        while (abs(r.maxEntry()-1) > tolerance) {

            for (int i=0;i<A.nCols();i++) {
                r.m[i][i] = 1/sqrt(vectorNormOne(A.m[i]));
            }

            A=multiplyDiagLeft(r, A);
            A=multiplyDiagRight(A, r);
        }
        
        return A;
    }

    //A simple loss-free scaling method for improving posit performance for cholesky factorization.
    //Should only use on diagonally dominant/positive-definite matrices. 
    //Works by normalizest largest diagonal entry to 1 and then scaling depending on the purpose.  
    void diagScale(int type, Matrix &A, vec &b) {
        T high=0;
        for (int i=0;i<nCols();i++) high=max(high, m[i][i]);
        T scaler = 1/high;
        double scalerD;
        cast(scalerD, scaler);
        scalerD=closestPowerTwo(scalerD);
        
        if (type==1) scalerD *= 4096;
        if (type==2) scalerD *= pow(2, pow(2, ES));
        A =    scale(scalerD, A);
        b = scaleVec(scalerD, b);
    }


    //This is an alternative approach to the above method which tries to push diagonals into posit interval.  
    void diagScaleAvg(int type, Matrix &A, vec &b) {
        T sum=0;
        for (int i=0;i<nCols();i++) sum+=abs(m[i][i]);
        T avg = sum/nCols();
        T scaler = 1/avg;
        cout << scaler << endl;
        
        double scalerD;
        cast(scalerD, scaler);
        scalerD = closestPowerTwo(scalerD);
        
        if(type==1) scalerD *= 4096; //If we're working with half then we should give it extra room. 

        A = scale(scalerD, *this);
        b = scaleVec(scalerD, b);//for (int i=0;i<b.size();i++) b[i]*=scalerD;
        //for (int i=0;i<nCols();i++) cout  << abs(A.m[i][i]) << endl;
    }

    //This is an alternative approach to the above method which tries to push diagonals into posit interval.  
    void diagScaleGeomAvg(int type, Matrix &A, vec &b) {
        T sum=0;
        for (int i=0;i<nCols();i++) sum+=log2(abs(m[i][i]));
        T avg = sum/nCols();
        T scaler = pow(2, -avg);
        
        double scalerD;
        cast(scalerD, scaler);
        //scalerD = closestPowerTwo(scalerD);
        scalerD = scalerD;
        
        if(type==1) scalerD *= 4096; //If we're working with half then we should give it extra room. 

        A = scale(scalerD, *this);
        b = scaleVec(scalerD, b);//for (int i=0;i<b.size();i++) b[i]*=scalerD;
    }

    //This is an alternative approach to the above method which tries to push diagonals into posit interval.  
    void scaleAvg(int type, Matrix &A, vec &b) {
        T sum=0;
        int nnz=0;
        for (int i=0;i<nRows();i++) {
            for (int j=0;j<nCols();j++) {
                if (m[i][j]!=0) { sum += log2(abs(m[i][j])); nnz++; }
            }
        }
        T avg = sum/nnz;
        T scaler = pow(2, -avg);
        
        double scalerD;
        cast(scalerD, scaler);
        scalerD = closestPowerTwo(scalerD);
        
        if(type==1) scalerD *= 4096; //If we're working with half then we should give it extra room. 

        A = scale(scalerD, *this);
        b = scaleVec(scalerD, b);//for (int i=0;i<b.size();i++) b[i]*=scalerD;
    }

    //In the case of CG, matrix-norm affects viability of posit. This scaling method tries to 
    //push matrix-norm to where posit tends to perform better. 
    void scaleNorm(Matrix &A, vec &b) {
        double optimal=(1<<10); //magic number

        double norm;
        cast(norm, infinityNorm());

        A = scale(optimal/closestPowerTwo(norm), *this);
        b = scaleVec(optimal/closestPowerTwo(norm),  b);
    }
    
    //returns iterations till convergence.
    int conjugateGradientSolver(double tolerance, const Matrix &A, vec B ,vec &X, string plotfile="", \
        string trafficfile="", bool clean=0) 
    {
        int n = A.nRows();
        
		int k = 0;
        for(int i=0;i<n;i++) X[i]=0; //zero out the vector

        vec R = B;
        vec P = R;
        
        ofstream plot;
        if (!plotfile.empty()) plot.open(plotfile, ofstream::app);

        ofstream traffic;
        if (!trafficfile.empty()) traffic.open(trafficfile, ofstream::trunc);

        T residual = vectorNorm( R );
        
        if (!trafficfile.empty()) Posit32::clearCounter();
        if (!plotfile.empty())    plot << residual;
        
        string delimiter="";
        T sum=0;
        T bNorm = vectorNorm(B);
        while ( residual > tolerance && k < 5000)
        {
            if (!trafficfile.empty()) Posit32::clearCounter();
            cout << k << " " << residual/bNorm << endl;
            
            if (k % 50 == 0 && clean) A.ConjugateGradientStep(A, P, R, X, B, 1); 
            else                      A.ConjugateGradientStep(A, P, R, X, B, 0); 

            if (!trafficfile.empty()) traffic << delimiter << Posit32::distillAdvantage(); 
            
            residual = vectorNorm( R );
            if (!plotfile.empty()) plot << "," << residual;
        
            k++;
            delimiter=",";
        }

        if (!trafficfile.empty()) traffic << endl;
        if (!plotfile.empty())    plot    << endl;
        plot.close();
        traffic.close();
        
        return k;
    }

    bool ConjugateGradientStep(const Matrix &A, vec &P, vec &R, vec &X, const vec &B, bool clean=0) const {
        vec Rold = R;
                                       
        vec AP = matVec( A, P );
        
        T alpha = innerProduct(R, R) / innerProduct( P, AP );

        X = vectorCombination( 1.0, X, alpha, P );  
        
        if (clean) R = vectorCombination( 1.0, B, -1.0,   matVec(A, X));
        else       R = vectorCombination( 1.0, R, -alpha, AP          );  

        T beta = innerProduct(R, R) / innerProduct( Rold, Rold ); 

        P = vectorCombination( 1.0, R, beta, P );

        return true;
    } 
    
    //T det(Matrix<T> &A){
    //    // do right diagonals - left diagonals with wrap-around (modulo)
    //   int i,j,k;
    //   T accum=0;
    //   if(A.nCols()!=A.nRows()){fprintf(stderr,"Fail: det() nRows not equal to ncols\n"); return 0;}
    //   else if(A.nCols()==1) return A.m[0][0];
    //   for(k=0;k<A.nCols();k++){
    //       T a=T(1.0),b=T(1.0);
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
        T accum=0; // Initialize result
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
            getCofactorSubmatrix(0,j,*this,A);  // construct a submatrix of cofactors
            accum += sign * m[0][j] * A.det();
            
            // terms are to be added with alternate sign
            sign = -sign;
        }
        return accum;
    }
    
    // compute the condition number for a matrix
    
    bool load(const char *filename){
        // from MMIO example
        int ret_code;
        FILE *f;
        int M, N, nz,n;
        T *val;
        MM_typecode matcode;
        string type = typeid(T).name();
        
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
            if (type == "f") {
                fscanf(f, "%d %d %f\n", &i, &j, &val); //means to extract from mtx depends on type. 
            } else if (type == "d") {
                fscanf(f, "%d %d %lg\n", &i, &j, &val);
            } else if (type == "N10half_float4halfE") {
                float fval;
                fscanf(f, "%d %d %f\n", &i, &j, &fval);
                val = fval; //extract as float. Cast down to half precision.
            } else { 
                printf("Type not supported. For Posits use loadMPF.");
                exit(1);
            }
            m[i-1][j-1] = val;
            if (mm_is_symmetric(matcode)) m[j-1][i-1] = val;
        }
        if (f !=stdin) fclose(f);
        
        return 1; // load a matrix from a matrixmarket file
    }

    //Load through maximum precision. Useful for Posits. Type must be
    //convertible from std::string e.g mpf_class or Posit. 
    bool loadMPF(const char *filename){
        // from MMIO example
        int ret_code;
        FILE *f;
        int M, N, nz,n;
        MM_typecode matcode;

        
        if ((f = fopen(filename, "r")) == NULL)
            return 1;
        
        if (mm_read_banner(f, &matcode) != 0)
        {
            printf("Could not process Matrix Market banner.\n");
            exit(1);
        }

        if (mm_is_pattern(matcode) || mm_is_complex(matcode)) {
            puts("Pattern or complex matrix.");
            return 1;
        }
        
        if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
            mm_is_sparse(matcode))
        {
            printf("Sorry, this application does not support ");
            printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
            exit(1);
        }

        
        if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
            return 1;
        this->setSize(M,N);

        char* entry = (char*) calloc(200, sizeof(char));
        for (n=0; n<nz; n++)
        {
            int i,j;
            fscanf(f, "%d %d %s\n", &i, &j, entry);
            m[i-1][j-1] = string(entry);
            if (mm_is_symmetric(matcode)) m[j-1][i-1] = string(entry);
            memset(entry, 0, 200*sizeof(char));
        }
        free(entry);
        
        if (f !=stdin) fclose(f);
        
        return 1; // load a matrix from a matrixmarket file
    }

    bool recordMatrix(const char *filename){
        // from MMIO example
        int ret_code;
        FILE *f;
        int M, N, nz,n;
        MM_typecode matcode;
        Posit32::initializeCounter();

        
        if ((f = fopen(filename, "r")) == NULL)
            return 0;
        
        if (mm_read_banner(f, &matcode) != 0)
        {
            return 0;
        }

        if (mm_is_pattern(matcode)) {
            puts("Pattern matrix.");
            return 0;
        }
        
        if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
            mm_is_sparse(matcode))
        {
            //printf("Sorry, this application does not support ");
            //printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
            //exit(1);
            return 0;
        }

        
        if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
            return 0;

        char* real = (char*) calloc(200, sizeof(char));
        char* imag = (char*) calloc(200, sizeof(char));

        if (!mm_is_complex(matcode)) {
            for (n=0; n<nz; n++)
            {
                int i,j;
                fscanf(f, "%d %d %s\n", &i, &j, real);
                Posit32gmp(string(real));
                if (mm_is_symmetric(matcode)) Posit32gmp(string(real));
                memset(real, 0, 200*sizeof(char));
            }
        } else{
            for (n=0; n<nz; n++)
            {
                int i,j;
                fscanf(f, "%d %d %s %s\n", &i, &j, real, imag);
                Posit32gmp(string(real));
                Posit32gmp(string(imag));
                if (mm_is_symmetric(matcode)) { Posit32gmp(string(real)); Posit32gmp(string(imag));}
                memset(real, 0, 200*sizeof(char));
                memset(imag, 0, 200*sizeof(char));
            }
        }
        free(real);
        free(imag);
        
        if (f !=stdin) fclose(f);
        return 1;
    }
};

#endif

