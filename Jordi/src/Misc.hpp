#ifndef MISC_HPP
#define MISC_HPP

#include "Main.hpp"
#include "Matrix.hpp"

double RandDouble(double low, double high);

int SumInt(int *x, int n);

int RandInt(int low, int high, double seed);

double InnerProd(double *x, double *y, int n);

double Residual2Norm(CSR A, double *x, double *b);

void LevelSets(CSR A, LevelSetData *lvl_set);

double SumDouble(double *x, int n);

uint64_t rdtsc();

#endif
