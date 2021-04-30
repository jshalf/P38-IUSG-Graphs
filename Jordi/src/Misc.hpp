#ifndef MISC_HPP
#define MISC_HPP

#include "Main.hpp"
#include "Matrix.hpp"

double RandDouble(double low, double high);

int SumInt(int *x, int n);

int RandInt(int low, int high, double seed);

double InnerProd(double *x, double *y, int n);

double Residual2Norm(Matrix A, double *x, double *b);

double Residual2Norm_CSC(Matrix A, /* sparse matrix data (input) */
                         double *x, /* solution (input) */
                         double *b /* right-hand side (input) */
                         );

void LevelSets(Matrix A, LevelSetData *lvl_set, int L_flag, int csc_flag);

void LevelSetsDestroy(LevelSetData *lvl_set);

double SumDouble(double *x, int n);

uint64_t rdtsc();

#endif
