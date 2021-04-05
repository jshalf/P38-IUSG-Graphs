#ifndef TRISOLVE_HPP
#define TRISOLVE_HPP

#include "../../src/Main.hpp"
#include "../../src/Matrix.hpp"

#define TRISOLVE_OMPFOR_SCHED static

void TriSolve_CSR(TriSolveData *ts,
                  CSR T,
                  int *T_perm,
                  double *x,
                  double *b);

void TriSolve_LevelSets_CSR(TriSolveData *ts,
                            LevelSetData lvl_set,
                            CSR T,
                            double *x,
                            double *b);

void TriSolve_Async_COO(TriSolveData *ts,
                        CSR T,
                        double *x,
                        double *b);

void TriSolve_Async_CSR(TriSolveData *ts,
                        CSR T,
                        double *x,
                        double *b);


#endif
