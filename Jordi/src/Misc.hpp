#ifndef MISC_HPP
#define MISC_HPP

#include "Main.hpp"
#include "Matrix.hpp"

double RandDouble(double low, double high);

double RandInt(int low, int high, double seed);

double InnerProd(double *x, double *y, int n);

double Residual2Norm(CSR A, double *x, double *b);

#endif
