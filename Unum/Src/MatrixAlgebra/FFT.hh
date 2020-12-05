#ifndef __FFT_HH_
#define __FFT_HH_

#include <iostream>
#include <complex>
#include "Posit32.hh"
#include "bfloat16.hh"
#include "Matrix.hh"

using namespace std;

/*Convenience Functions*/

//mod that only returns positive values. 
int mod(int m, int n)  { return m % n >= 0 ? (m % n) : ((m % n) + n); }

int powerOfTwo(uint n) {return (n & (n-1)) == 0;}

template <class T>
void partition(vector<T> &out, vector<T> &in, int odd) {
	for (int i=0;i<out.size();++i) out[i]=in[2*i+odd];
}

//flatten matrix into 1D vector.
template <class T>
void flatten(vector<T> &f, Matrix<T> A) {
    if (f.size()!=A.nCols()*A.nRows()) {fprintf(stderr, "Invalid dimensions."); return;}

    for (int i=0;i<A.nRows();i++) {
        for (int j=0;j<A.nCols();j++) {
            f[i*A.nCols()+j] = A.m[i][j];
        }
    }
}

/*DFT/IDFT algorithms */

//1D computation of dft using recursive fft<T>. 
template<class T>
void fft(vector<complex<T> >& out, vector<complex<T> > &in)
{
    if (out.size() != in.size()) {fprintf(stderr, "Invalid dimensions."); return;}
    if (!powerOfTwo(in.size()))  {fprintf(stderr, "Signal should have length a power of two."); return;}

    int n = in.size();
    if (n <= 1) return;
    T pi = mpf_class(PI).get_d(); //Cast pi to double and then to T. 
 
    // divide
    vector<complex<T> > even(n/2); 
    vector<complex<T> >  odd(n/2); 
    partition<complex<T> >(even, in, 0);
    partition<complex<T> >(odd,  in, 1);
 
    // conquer
    fft<T>(even, even);
    fft<T>(odd, odd);
 
    // combine
    for (int k = 0; k < n/2; ++k)
    {
        T half_periods = ((T) -2) * ((T) k) / ((T) n); //make sure precision stays in T. 
        complex<T> t = complex<T>(cos(half_periods*pi), sin(half_periods*pi)) * odd[k];
        out[k    ] = even[k] + t;
        out[k+n/2] = even[k] - t;
    }
}

//inverse dft.
template <class T>
void ifft(vector<complex<T> >& out, vector<complex<T> > &in) {
    if (!powerOfTwo(in.size())) {fprintf(stderr, "Input signal should have length a power of two."); return;}

    int n = in.size();
    
    for (int i=0;i<n;i++) out[i]=conj(in[i]);
 
    // forward fft<T>
    fft<T>(out, out);
 
    // conjugate the complex numbers again
     
    for (int i=0;i<n;i++) out[i] = conj(out[i]);
    
    // scale the numbers
    for (int i=0;i<n;i++) out[i] /= n;
}

//Performs a circular convolution of x, y the slow and stable way.  
template <class T>
void convolve(vector<complex<T> > &out, vector<complex<T> > x, vector<complex<T> > y) {
    
    int n = x.size();
    out = vector<complex<T> >(n);
    for (int i=0;i<n;i++) out[i]=complex<T>(0,0);
    
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            out[i] += x[j]*y[mod(i-j,n)];
        }
    }
}

template <class T>
void convolveFFT(vector<complex<T> > &out, vector<complex<T> > x, vector<complex<T> > y) {
    if (x.size() != y.size()) 
        {fprintf(stderr, "Input signals should have the same length."); return;}
    
    fft<T>(x, x);
    fft<T>(y, y);

    for (int i=0;i<out.size();i++) out[i]=x[i]*y[i]; 
    ifft<T>(out, out);
}

//Performs a 2D circular convolution the straightforward way. in and kernel should be same size.
template <class T>
void convolve2D(Matrix<complex<T> > &out, Matrix<complex<T> > x,\
    Matrix<complex<T> > y) {
    
    if (x.nCols() != y.nCols() || x.nRows() != y.nRows()) 
        {fprintf(stderr, "in/out differing in size."); return; }
    
    out.setSize(x.nRows(),x.nCols());
    
    for (int i=0;i<x.nRows();i++) {
        for (int j=0;j<x.nCols();j++) {
            out.m[i][j] = complex<T>(0,0);
            for (int k=0;k<x.nRows();k++) {
                for (int l=0;l<x.nCols();l++) {
                    out.m[i][j] += x.m[k][l] * y.m[mod(i-k,y.nRows())][mod(j-l,y.nCols())];
                }
            }
        }
    }
}

//find dft using fft<T> algorithm.
template <class T>
void fft2D(Matrix<complex<T> > &out, Matrix<complex<T> > in) {

    if (!(in.isSquare()) || !powerOfTwo(in.nRows())) 
        {fprintf(stderr, "Please enter square input with power of two dimensions."); return; }

    int N = in.nRows();
    out.setSize(N,N);
    Matrix<complex<T> > partial(in);
    
    //Store matrix where rows are 1D fft<T>s of input rows. 
    for (int m=0;m<N;m++) {
        //vector<complex<T> > v = partial.m[m];
        //fft<T>(v, v);
        fft<T>(partial.m[m], partial.m[m]);
        //for(int i=0;i<N;i++) partial.m[m][i]=v[i]; //copy
    }

    partial.transpose(); //For easy access to columns. 

    //Take fft<T>s of columns of partial. 
    for (int n=0;n<N;n++) {
        //vector<complex<T> > v = partial.m[n];
        //fft<T>(v, v);
        //for(int i=0;i<N;i++) out.m[i][n]=v[i]; //copy
        fft<T>(partial.m[n], partial.m[n]);
        for(int i=0;i<N;i++) out.m[i][n]=partial.m[n][i];
    }
}

//Compute inverse DFT using FFT. (Same strategy as 1D case)
template <class T>
void ifft2D(Matrix<complex<T> > &out, Matrix<complex<T> > in) {
    if (!(in.isSquare()) || !powerOfTwo(in.nRows())) 
        {fprintf(stderr, "Please enter square input with power of two dimensions."); return; }

    int N = in.nRows();
    out.setSize(N, N);

    //Conjugate all entries.  
    for (int i=0;i<N;i++) {
        for (int j=0;j<N;j++) {
            out.m[i][j] = conj(in.m[i][j]);
        }
    }

    fft2D<T>(out, out);

    for (int i=0;i<N;i++) {
        for (int j=0;j<N;j++) {
            out.m[i][j] = conj(out.m[i][j]);
            out.m[i][j] /= (N*N);
        }
    }
}

template <class T>
void convolve2DFFT(Matrix<complex<T> > &out, Matrix<complex<T> > x,\
    Matrix<complex<T> > y) {
    if (x.nCols() != y.nCols() || x.nRows() != y.nRows()) 
        {fprintf(stderr, "in/out differing in size."); return; }
    if (!(x.isSquare()) || !(y.isSquare()) || !powerOfTwo(x.nRows())) 
        {fprintf(stderr, "Please enter square inputs with power of two dimensions."); return; }
    
    int n = x.nRows();
    out.setSize(n, n);
    
    fft2D<T>(x, x);
    fft2D<T>(y, y);

    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            out.m[i][j] = x.m[i][j]*y.m[i][j];
        }
    }

    ifft2D<T>(out, out);
}

template <class T>
void correlate2DFFT(Matrix<complex<T> > &out, Matrix<complex<T> > x,\
    Matrix<complex<T> > y) {
    if (x.nCols() != y.nCols() || x.nRows() != y.nRows()) 
        {fprintf(stderr, "in/out differing in size."); return; }
    if (!(x.isSquare()) || !(y.isSquare()) || !powerOfTwo(x.nRows())) 
        {fprintf(stderr, "Please enter square inputs with power of two dimensions."); return; }
    
    int n = x.nRows();
    out.setSize(n, n);
    
    fft2D<T>(x, x);
    fft2D<T>(y, y);

    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            y.m[i][j] = conj(y.m[i][j]);
            complex<T> hprod = x.m[i][j]*y.m[i][j];
            hprod /= abs(hprod);
            out.m[i][j] = hprod;
        }
    }

    ifft2D<T>(out, out);
}

#endif
