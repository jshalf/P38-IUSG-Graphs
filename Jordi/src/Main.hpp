#ifndef MAIN_HPP
#define MAIN_HPP

#define SYNC_JACOBI  0
#define ASYNC_JACOBI 1
#define SYNC_BLOCK_JACOBI  2
#define ASYNC_BLOCK_JACOBI 3

#define TRISOLVE_ASYNC 0
#define TRISOLVE_LEVEL_SCHEDULED 1

#define ILU_ASYNC 0
#define ILU_LEVEL_SCHEDULED 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream> 
#include <algorithm>
#include <omp.h>
#include <vector>
#include <random>
#include <queue>

using namespace std;

typedef struct{
   int i;
   int j;
   double val;
}Triplet_AOS;

typedef struct{
   int solver_type;
   int num_threads;
   int num_iters;
   int async_flag;
   int atomic_flag;
   int AAT_flag;
   int expand_flag;
   int coo_flag;
   int omp_for_flag;
   int MsgQ_flag;
}InputData;

typedef struct{
   double solve_wtime;
   double setup_wtime;
   double solve_wtime_thread;
   double *atomic_wtime_vec;
   double *solve_wtime_vec;
   int *num_relax;
   int **relax_hist;
}OutputData;

typedef struct{
   InputData input;
   OutputData output;
   double *y1_expand;
   double *y2_expand;
}MatVecData;

typedef struct{
   InputData input;
   OutputData output;
}SolverData;

typedef struct{
   int *perm;
   vector<int> level_size;
   vector<int> level_start;
   int num_levels;
}LevelSetData;

typedef struct{
   InputData input;
   OutputData output;
   LevelSetData L_lvl_set;
   LevelSetData U_lvl_set;
}TriSolveData;

typedef struct{
   InputData input;
   OutputData output;
   LevelSetData L_lvl_set;
   LevelSetData U_lvl_set;
}ILUData;

#endif
