#include "Quire.hh"
#include "Matrix.hh"




/*Perform mpf inner product and cast down to posit. Would correspond with using quire register to 
accumulate and then clearing quire at the end. Supposedly as long as the sum of the entries in each vector
does not exceed 2^30 this method should return the same result as doing arithmetic in the quire. */
/* 
Posit32gmp Quire::innerProductQ(vector<Posit32gmp> x, vector<Posit32gmp> y) {
    if(x.size() != y.size()) { fprintf(stderr, "Invalid dimensions."); return 0; }

    int n = x.size();

    vector<mpf_class> xm(n), ym(n);
    mpf_class rm(0);
    Posit32gmp rp;

    upcast(xm, x);
    upcast(ym, y);

    for (int i=0;i<n;i++) rm += xm[i]*ym[i]; 

    rp.set(rm);
    return rp;
}

//Basically a repeated version of inner product. We would clear quire after calculating
//each entry of the output vector. 
void Quire::matVecQ(vector<Posit32gmp> &out, Matrix<Posit32gmp> A, vector<Posit32gmp> x) {
    if(A.nCols() != x.size() || out.size() != A.nRows()) { fprintf(stderr,  "Invalid dimensions."); return; }
    
    Matrix<mpf_class> AM(A.nRows(), A.nCols());
    vector<mpf_class> xm(x.size()), bm(A.nRows());

    for (int i=0;i<A.nRows();i++) upcast(AM.m[i], A.m[i]);

    upcast(xm, x);

    bm = matVec(AM, xm);
    downcast(out, bm);
}

bool Quire::ConjugateGradientStepQ(Matrix<Posit32gmp> &A, vector<Posit32gmp>
  &P, vector<Posit32gmp> &R, vector<Posit32gmp> &X, vector<Posit32gmp> &B, bool clean) {
    int n = X.size();
        
    vector<Posit32gmp> Rold = R;                                         
    vector<Posit32gmp> AP(n,n);
    matVecQ(AP, A, P);
    
    Posit32gmp alpha = innerProductQ(R, R) / innerProductQ( P, AP );
    X = A.vectorCombination( 1.0, X, alpha, P );            
    
    vector<Posit32gmp> currentOutput(n,n);
    matVecQ(currentOutput, A, X);
    
    if (clean) R = A.vectorCombination( 1.0, B, -1.0, currentOutput ); 
    else R = A.vectorCombination( 1.0, R, -alpha, AP );

    Posit32gmp beta = innerProductQ(R, R) / innerProductQ( Rold, Rold );
    P = A.vectorCombination( 1.0, R, beta, P );

    return true;
}

int Quire::conjugateGradientSolverQ(double tolerance, Matrix<Posit32gmp> &A, vector<Posit32gmp> \
    &B ,vector<Posit32gmp> & X, string plotfile, string trafficplot, bool clean) {

    int n = A.nRows();
    
    int k = 0;
    for(int i=0;i<n;i++) X[i]=0; 

    ofstream file, tplot;
    if (!plotfile.empty())    file.open(plotfile, ofstream::app);
    if (!trafficplot.empty()) tplot.open(trafficplot, ofstream::trunc);

    vector<Posit32gmp> R = B;
    vector<Posit32gmp> P = R;
    Posit32gmp r = sqrt( innerProductQ(R, R) );
    if (!plotfile.empty()) file << r;

    string delimiter = "";
    while ( r > tolerance && k < 15*n )
    {
        if (!trafficplot.empty()) Posit32::clearCounter();
        if (k % 50 == 0 && clean) ConjugateGradientStepQ(A, P, R, X, B, 1); 
        else                      ConjugateGradientStepQ(A, P, R, X, B, 0); 
        if (!trafficplot.empty() && !trafficplot.empty()) tplot << delimiter << Posit32::distillAdvantage();
        
        r = sqrt( innerProductQ(R, R) );
        if (!plotfile.empty()) file << "," << r;
        cout << k << " " << r;

        delimiter=",";
        k++;
    }

    file << endl;
    tplot << endl;
    file.close();
    tplot.close();
    
    return k;
}


//For floats. 


float Quire::innerProductQ(vector<float> x, vector<float> y) {
    if(x.size() != y.size()) { fprintf(stderr, "Invalid dimensions."); return 0; }

    int n = x.size();

    vector<mpf_class> xm(n), ym(n);
    mpf_class rm(0);
    float rp;

    upcast(xm, x);
    upcast(ym, y);

    for (int i=0;i<n;i++) rm += xm[i]*ym[i]; 

    rp = rm.get_d();
    return rp;
}

//Basically a repeated version of inner product. We would clear quire after calculating
//each entry of the output vector. 
void Quire::matVecQ(vector<float> &out, Matrix<float> A, vector<float> x) {
    if(A.nCols() != x.size() || out.size() != A.nRows()) { fprintf(stderr,  "Invalid dimensions."); return; }
    
    Matrix<mpf_class> AM(A.nRows(), A.nCols());
    vector<mpf_class> xm(x.size()), bm(A.nRows());

    for (int i=0;i<A.nRows();i++) upcast(AM.m[i], A.m[i]);

    upcast(xm, x);

    bm = matVec(AM, xm);
    downcast(out, bm);
}

bool Quire::ConjugateGradientStepQ(Matrix<float> &A, vector<float>
  &P, vector<float> &R, vector<float> &X, vector<float> &B, bool clean) {
    int n = X.size();
        
    vector<float> Rold = R;                                         
    vector<float> AP(n,n);
    matVecQ(AP, A, P);
    
    float alpha = innerProductQ(R, R) / innerProductQ( P, AP );
    X = A.vectorCombination( 1.0, X, alpha, P );            
    
    vector<float> currentOutput(n,n);
    matVecQ(currentOutput, A, X);
    
    if (clean) R = A.vectorCombination( 1.0, B, -1.0, currentOutput ); 
    else       R = A.vectorCombination( 1.0, R, -alpha, AP );

    float beta = innerProductQ(R, R) / innerProductQ( Rold, Rold );
    P = A.vectorCombination( 1.0, R, beta, P );

    return true;
}

int Quire::conjugateGradientSolverQ(double tolerance, Matrix<float> &A, vector<float> \
    &B ,vector<float> & X, string plotfile, bool clean) {

    int n = A.nRows();
    
    int k = 0;
    for(int i=0;i<n;i++) X[i]=0; 

    ofstream file;
    if (!plotfile.empty()) file.open(plotfile, ofstream::app);

    vector<float> R = B;
    vector<float> P = R;
    float r = sqrt( innerProductQ(R, R) );
    if (!plotfile.empty()) file << r;
    
    while ( r > tolerance && k < 15*n )
    {
        if (k % 50 == 0   && clean) ConjugateGradientStepQ(A, P, R, X, B, 1); 
        else                        ConjugateGradientStepQ(A, P, R, X, B, 0); 
        r = sqrt( innerProductQ(R, R) );
        if (!plotfile.empty()) file << "," << r;
        cout << k << " " << r;

        k++;
    }

    file << endl;
    file.close();
    
    return k;
}

//for bfloat


bfloat16 Quire::innerProductQ(vector<bfloat16> x, vector<bfloat16> y) {
    if(x.size() != y.size()) { fprintf(stderr, "Invalid dimensions."); return 0; }

    int n = x.size();

    vector<mpf_class> xm(n), ym(n);
    mpf_class rm(0);
    bfloat16 rp;

    for (int i=0;i<n;i++) { xm[i]=x[i].f; ym[i]=y[i].f; }

    for (int i=0;i<n;i++) { rm += xm[i]*ym[i]; } 

    rp = rm.get_d();
    return rp;
}

//Basically a repeated version of inner product. We would clear quire after calculating
//each entry of the output vector. 
void Quire::matVecQ(vector<bfloat16> &out, Matrix<bfloat16> A, vector<bfloat16> x) {
    if(A.nCols() != x.size() || out.size() != A.nRows()) { fprintf(stderr,  "Invalid dimensions."); return; }
    
    Matrix<mpf_class> AM(A.nRows(), A.nCols());
    vector<mpf_class> xm(x.size()), bm(A.nRows());

    for (int i=0;i<A.nRows();i++) {
        for (int j=0;j<A.nCols();j++) {
            AM.m[i][j] = A.m[i][j].f;
        }
    }

    for (int i=0;i<x.size();  i++) { xm[i]=x[i].f;   }
    bm = matVec(AM, xm);
    for (int i=0;i<out.size();i++) { out[i]=bm[i].get_d(); }
}

bool Quire::ConjugateGradientStepQ(Matrix<bfloat16> &A, vector<bfloat16>
  &P, vector<bfloat16> &R, vector<bfloat16> &X, vector<bfloat16> &B, bool clean) {
    int n = X.size();
        
    vector<bfloat16> Rold = R;                                         
    vector<bfloat16> AP(n,n);
    matVecQ(AP, A, P);
    
    bfloat16 alpha = innerProductQ(R, R) / innerProductQ( P, AP );
    X = A.vectorCombination( 1.0, X, alpha, P );            
    
    vector<bfloat16> currentOutput(n,n);
    matVecQ(currentOutput, A, X);
    
    if (clean) R = A.vectorCombination( 1.0, B, -1.0, currentOutput ); 
    else       R = A.vectorCombination( 1.0, R, -alpha, AP );

    bfloat16 beta = innerProductQ(R, R) / innerProductQ( Rold, Rold );
    P = A.vectorCombination( 1.0, R, beta, P );

    return true;
}

int Quire::conjugateGradientSolverQ(double tolerance, Matrix<bfloat16> &A, vector<bfloat16> \
    &B ,vector<bfloat16> & X, string plotfile, bool clean) {

    int n = A.nRows();
    
    int k = 0;
    for(int i=0;i<n;i++) X[i]=0; 

    ofstream file;
    if (!plotfile.empty()) file.open(plotfile, ofstream::app);

    vector<bfloat16> R = B;
    vector<bfloat16> P = R;
    bfloat16 r = sqrt( innerProductQ(R, R) );
    if (!plotfile.empty()) file << r;

    while ( r > tolerance && k < 15*n )
    {
        if (k % 50 == 0 && clean) ConjugateGradientStepQ(A, P, R, X, B, 1); 
        else                      ConjugateGradientStepQ(A, P, R, X, B, 0); 
        
        r = sqrt( innerProductQ(R, R) );
        
        if (!plotfile.empty()) file << "," << r;
        cout << k << " " << r;

        k++;
    }

    file << endl;
    file.close();
    
    return k;
}
*/

