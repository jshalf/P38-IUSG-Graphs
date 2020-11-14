#ifndef JACOBI_HPP
#define JACOBI_HPP

#include "Main.hpp"
#include "Matrix.hpp"

void MatVec(MatVecData *mv,
            CSR A,
            double *x,
            double *y);

void MatVecT(MatVecData *mv,
             CSR A,
             double *x,
             double *y1,
             double *y2);

#endif
