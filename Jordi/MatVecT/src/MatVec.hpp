#ifndef JACOBI_HPP
#define JACOBI_HPP

#include "../../src/Main.hpp"
#include "../../src/Matrix.hpp"

void MatVec_CSR(MatVecData *mv,
                CSR A,
                double *x,
                double *y);

void MatVecT_CSR(MatVecData *mv,
                 CSR A,
                 double *x,
                 double *y1,
                 double *y2);

void MatVecT_COO(MatVecData *mv,
                 CSR A,
                 double *x,
                 double *y1,
                 double *y2);

#endif
