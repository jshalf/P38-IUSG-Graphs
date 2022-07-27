#ifndef MISC_HPP
#define MISC_HPP

#define CACHE_LINE_SIZE 64

#include "Main.hpp"
//#include "Matrix.hpp"

enum class FileType {
   bin,
   text
};

enum class ParallelCommType {
   atomic,
   MsgQ
};

struct alignas(CACHE_LINE_SIZE) cache {
   int c;
};

template<typename T>
struct alignas(CACHE_LINE_SIZE) Cache {
   T data;
};

double RandDouble(double low, double high);

int SumInt(int *x, int n);

int RandInt(int low, int high, double seed);

double InnerProd(double *x, double *y, int n);

//double Residual2Norm(Matrix A, double *x, double *b);
//
//double Residual2Norm_CSC(Matrix A, /* sparse matrix data (input) */
//                         double *x, /* solution (input) */
//                         double *b /* right-hand side (input) */
//                         );

double SumDouble(double *x, int n);

uint64_t rdtsc();

void PrintDummy(double dummy);

void PrintTraces(char *filename, TraceData trace_loc);

#endif
