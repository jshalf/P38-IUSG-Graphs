#ifndef  __QUIRE_HH_
#define  __QUIRE_HH_

#include <gmpxx.h>
#include <cmath>
#include "Posit32.hh"
#include "Matrix.hh"
#include "bfloat16.hh"
#include "helpers.hh"
using namespace std;

/*
Class used to perform scratchpaper operations for Posit arithmetic. 
Abbreviates Gustaffson's idea of Quire/Just does the stuff we care about (overrides 
basic linalg operations with deferred rounding). 
*/
template<class T>
T innerProductQ(vector<T> x, vector<T> y) {
    if(x.size() != y.size()) { fprintf(stderr, "Invalid dimensions."); return 0; }

    int n = x.size();

    vector<mpf_class> xm(n), ym(n);
    mpf_class rm(0);
    T rp;

    upcast(xm, x);
    upcast(ym, y);

    for (int i=0;i<n;i++) rm += xm[i]*ym[i]; 

    cast(rp, rm);
    return rp;
}

template <class T> 
void matVecQ(vector<T> &out, Matrix<T> A, vector<T> x) {
    if(A.nCols() != x.size() || out.size() != A.nRows()) { fprintf(stderr, "Invalid dimensions."); return; }

    Matrix<mpf_class> AM(A.nRows(), A.nCols());
    vector<mpf_class> xm(x.size()), bm(A.nRows());

    for (int i=0;i<A.nRows();i++) upcast(AM.m[i], A.m[i]);

    upcast(xm, x);

    bm = matVec(AM, xm);
    downcast(out, bm);
}

template <class T>
bool ConjugateGradientStepQ(Matrix<T> &A, vector<T>
    &P, vector<T> &R, vector<T> &X, vector<T> &B, bool clean) {
    int n = X.size();

    vector<T> Rold = R;                                         
    vector<T> AP(n,n);
    matVecQ(AP, A, P);

    T alpha = innerProductQ(R, R) / innerProductQ( P, AP );
    X = A.vectorCombination( 1.0, X, alpha, P );            
    
    if (clean) {
        vector<T> currentOutput(n,n);
        matVecQ(currentOutput, A, X);
        R = A.vectorCombination( 1.0, B, -1.0, currentOutput ); 
    }
    else {
        R = A.vectorCombination( 1.0, R, -alpha, AP );
    }

    T beta = innerProductQ(R, R) / innerProductQ( Rold, Rold );
    P = A.vectorCombination( 1.0, R, beta, P );

    return true;
}

template <class T>
int conjugateGradientSolverQ(double tolerance, Matrix<T> &A, vector<T> \
    &B ,vector<T> & X, string plotfile, string trafficplot, bool clean) {
    
    int n = A.nRows();

    int k = 0;
    for(int i=0;i<n;i++) X[i]=0; 

    bool plotResidual = !(plotfile.empty());
    bool plotTraffic  = !(trafficplot.empty());

    ofstream file, tplot;
    if (plotResidual) file.open(plotfile, ofstream::app);
    if (plotTraffic)  tplot.open(trafficplot, ofstream::trunc);

    vector<T> R = B;
    vector<T> P = R;
    T r = sqrt( innerProductQ(R, R) );
    //vector<T> residual(n);
    
    if (plotResidual) file << r;

    string delimiter = "";
    while ( r > tolerance && k < 12000 )
    {
        if (plotTraffic) Posit32::clearCounter();
        if (k % 50 == 0 && clean) ConjugateGradientStepQ(A, P, R, X, B, 1); 
        else                      ConjugateGradientStepQ(A, P, R, X, B, 0); 
        if (plotTraffic) tplot << delimiter << Posit32::distillAdvantage();

        r = sqrt( innerProductQ(R, R) );

        //residual = A.vectorCombination(1.0, B, -1.0, matVec(A, X));
        //r = A.vectorNorm(residual);
        cout << k << ": " << r << endl;
        if (plotResidual) file << "," << r;

        delimiter=",";
        k++;
    }
    file << endl;
    tplot << endl;
    file.close();
    tplot.close();

    return k;
}       

        /*
        Posit32gmp innerProductQ(vector<Posit32gmp> x, vector<Posit32gmp> y);
        
        void matVecQ(vector<Posit32gmp> &out, Matrix<Posit32gmp> A, vector<Posit32gmp> x);

        bool ConjugateGradientStepQ(Matrix<Posit32gmp> &A, vector<Posit32gmp> &P, \
            vector<Posit32gmp> &R, vector<Posit32gmp> &X, vector<Posit32gmp> &B, bool clean);

        int conjugateGradientSolverQ(double tolerance, Matrix<Posit32gmp> &A, \
            vector<Posit32gmp> &B, vector<Posit32gmp> & X, string plotfile, string trafficPlot, bool clean);

        //float

        float innerProductQ(vector<float> x, vector<float> y);
        
        void matVecQ(vector<float> &out, Matrix<float> A, vector<float> x);

        bool ConjugateGradientStepQ(Matrix<float> &A, vector<float> &P, \
            vector<float> &R, vector<float> &X, vector<float> &B, bool clean);

        int conjugateGradientSolverQ(double tolerance, Matrix<float> &A, \
            vector<float> &B, vector<float> & X, string plotfile, bool clean);

        //bfloat

        bfloat16 innerProductQ(vector<bfloat16> x, vector<bfloat16> y);
        
        void matVecQ(vector<bfloat16> &out, Matrix<bfloat16> A, vector<bfloat16> x);

        bool ConjugateGradientStepQ(Matrix<bfloat16> &A, vector<bfloat16> &P, \
            vector<bfloat16> &R, vector<bfloat16> &X, vector<bfloat16> &B, bool clean);

        int conjugateGradientSolverQ(double tolerance, Matrix<bfloat16> &A, \
            vector<bfloat16> &B, vector<bfloat16> & X, string plotfile, bool clean);    
        */
//};



#endif
