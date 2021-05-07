#ifndef JACOBI_HPP
#define JACOBI_HPP

#include "../../src/Main.hpp"
#include "../../src/Matrix.hpp"

void MatVecT_CSR_Seq(MatVecData *mv,
                     Matrix A,
                     double *x,
                     double *y);

void MatVecT_CSR(MatVecData *mv,
                 Matrix A,
                 double *x,
                 double *y);

#endif
