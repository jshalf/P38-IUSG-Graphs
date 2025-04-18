#ifndef __HELPERS_HH_
#define __HELPERS_HH_

#include <vector>
#include <complex>
//#include <experimental/filesystem>
#include <gmpxx.h>
#include "mmio.hh"
#include "../Posit32.hh"
#include "../half.hpp"
#include "../bfloat16.hh"

using namespace std;
//using namespace experimental::filesystem;
using half_float::half;

//Get random mpf value in [0, 1]
mpf_class rand(gmp_randstate_t r);

//Get vector of random values in [0, itemSize]
void getVector(vector<mpf_class> &v, int itemSize);

//Same as getVector but complex entries with imaginary components set to zero. 
void getVectorComplex(vector<complex<mpf_class> > &v, int itemSize);

//mpf_class mse(vector<complex<mpf_class> >  modified, vector<complex<mpf_class> > original);

mpf_class mse(vector<vector<mpf_class> > modified, vector<vector<mpf_class> > original);
// mpf_class mse(vector<vector<complex<mpf_class> > > modified, vector<vector<complex<mpf_class> > > original);


int closestPowerTwo(int n);
double closestPowerTwo(double d);

double closestPowerFour(double d);

void cast(float &a, half       b);
void cast(float &a, bfloat16   b);
void cast(float &a, Posit32gmp b);
void cast(half  &a,      float b);
void cast(bfloat16 &a,   float b);
void cast(Posit32gmp &a, float b);

void cast(double &a, half       b);
void cast(double &a, bfloat16   b);
void cast(double &a, Posit32gmp b);
void cast(double &a, double     b);
void cast(double &a, mpf_class  b);

void cast(half  &a,      double b);
void cast(bfloat16 &a,   double b);
void cast(Posit32gmp &a, double b);

void cast(double &a,     mpf_class b);
void cast(float &a,      mpf_class b);
void cast(bfloat16 &a,   mpf_class b);
void cast(half &a,       mpf_class b);
void cast(Posit32gmp &a, mpf_class b);

void cast(mpf_class &a, double     b);
void cast(mpf_class &a, float      b);
void cast(mpf_class &a, bfloat16   b);   
void cast(mpf_class &a, half       b);
void cast(mpf_class &a, Posit32gmp b);

void downcast(vector<double    > &out, vector<mpf_class> in);
void downcast(vector<float     > &out, vector<mpf_class> in);
void downcast(vector<bfloat16  > &out, vector<mpf_class> in);
void downcast(vector<half      > &out, vector<mpf_class> in);
void downcast(vector<Posit32gmp> &out, vector<mpf_class> in);
void downcast(vector<mpf_class > &out, vector<mpf_class> in);

void upcast(vector<mpf_class> &out, vector<double    > in);
void upcast(vector<mpf_class> &out, vector<float     > in);
void upcast(vector<mpf_class> &out, vector<bfloat16  > in);
void upcast(vector<mpf_class> &out, vector<half      > in);
void upcast(vector<mpf_class> &out, vector<Posit32gmp> in);
void upcast(vector<mpf_class> &out, vector<mpf_class > in);

template<class T, class E>
void convert(vector<T> &out, vector<E> in) {
    vector<mpf_class> m(in.size());
    upcast(m, in);
    downcast(out, m);
}

template<class T>
void downcastComplex(vector<complex<T> > &out, const vector<complex<mpf_class> > &in) {
    int n = in.size();
    vector<mpf_class> reM(n), imM(n);
    for (int i=0;i<n;i++) { reM[i]=in[i].real(); imM[i]=in[i].imag(); } 

    vector<T> re(n), im(n);
    downcast(re, reM);
    downcast(im, imM);

    for (int i=0;i<n;i++) out[i]=complex<T>(re[i], im[i]);
}

template<class T>
void upcastComplex(vector<complex<mpf_class> > &out, const vector<complex<T> > &in) {
    int n = in.size();
    vector<T> re(n), im(n);
    for (int i=0;i<n;i++) { re[i]=in[i].real(); im[i]=in[i].imag(); } 

    vector<mpf_class> reM(n), imM(n);
    upcast(reM, re);
    upcast(imM, im);

    for (int i=0;i<n;i++) out[i]=complex<mpf_class>(reM[i], imM[i]);
}

template <typename E>
void printVector(vector<E> v) { 
    cout << '{';
    for (auto i = v.begin(); i < v.end() - 1; i++) cout << *i << ",";
    cout << *(v.end()-1) << '}' << endl << endl;
}

//bool recordToPosit(const char *filename);

//void writeAllMatrixInfo(string directory, string filename);

#endif
