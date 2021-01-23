#ifndef TRISOLVE_HPP
#define TRISOLVE_HPP

#include "../../src/Main.hpp"
#include "../../src/Matrix.hpp"

#define TRISOLVE_OMPFOR_SCHED static

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
