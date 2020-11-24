#ifndef JACOBI_HPP
#define JACOBI_HPP

#include "../../src/Main.hpp"
#include "../../src/Matrix.hpp"

void TriSolve_CSR(TriSolveData *ts,
                  CSR L,
                  CSR U,
                  int *L_perm,
                  int *U_perm,
                  double *x,
                  double *y,
                  double *b);

void TriSolve_LevelSets_CSR(TriSolveData *ts,
                            CSR L,
                            CSR U,
                            double *x,
                            double *y,
                            double *b);

void TriSolve_FineGrained_COO(TriSolveData *ts,
                              CSR L,
                              CSR U,
                              double *x,
                              double *y,
                              double *b);


#endif
