#ifndef TRISOLVE_HPP
#define TRISOLVE_HPP

#include "../../src/Main.hpp"
#include "../../src/Matrix.hpp"

#define TRISOLVE_OMPFOR_SCHED static

void TriSolve_Seq(TriSolveData *ts,
                  Matrix T,
                  int *T_perm,
                  double *x,
                  double *b);

void TriSolve_LevelSchedule(TriSolveData *ts,
                            LevelSetData lvl_set,
                            Matrix T,
                            double *x,
                            double *b);

void TriSolve_Async(TriSolveData *ts,
                    Matrix T,
                    double *x,
                    double *b);

void TriSolve_AtomicCounter(TriSolveData *ts,
                            Matrix T,
                            double *x,
                            double *b);


#endif
