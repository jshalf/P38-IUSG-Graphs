#ifndef __BENCHMARKS_HH_
#define __BENCHMARKS_HH_

#include <string.h>
#include <gmpxx.h>
#include <fstream>
#include "helpers.hh"
#include "Matrix.hh"

void CGTest(Matrix<mpf_class> M, string matrixname, ofstream::openmode mode);

void trisolveTest(Matrix<mpf_class> M);

void fftTest();

void fftTest2D();

void convolveImage(string imgFilename, string kernelFilename, string output);

#endif