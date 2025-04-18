/* *
 * header file that defines structs to be used by all benchmarks
 * */

#ifndef MAIN_HPP
#define MAIN_HPP

/* types of solvers used in the AsyncJacobi benchmark */
#define SYNC_JACOBI  0
#define ASYNC_JACOBI 1
#define SYNC_BLOCK_JACOBI  2
#define ASYNC_BLOCK_JACOBI 3

#define MATRIX_STORAGE_CSR 0
#define MATRIX_STORAGE_CSC 1
#define MATRIX_STORAGE_COO 2
#define MATRIX_STORAGE_DENSE 3

#define LEVEL_SETS_SEQ_SETUP 0
#define LEVEL_SETS_ASYNC_SETUP 1

/* types of solvers used in the ILU benchmark */
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
#include <inttypes.h>
#include <numeric>
#include <smmintrin.h>
#include <pthread.h>
#include <unordered_map>
#include <atomic>
#include <functional>
#include <thread>
#include <mutex>
//#include <barrier>

using namespace std;

/* triplet struct used for reading binary matrix files */
typedef struct{
   int i; /* row index */
   int j; /* column index */
   double val; /* value of (i,j) element */
}Triplet_AOS;

/* input data (primarily taken in at command line) */
typedef struct{
   int solver_type; /* solver type */
   int num_threads; /* number of threads */
   int num_iters; /* number of iterations */
   int async_flag; /* is solver asynchronous? */
   int atomic_flag; /* is solver using atomics? */
   int AAT_flag; /* will both the sparse matrix-vector and sparse matrix-transpose-vector products both be computed? (only used in MatVecT benchmark) */
   int expand_flag; /* will the ``expand''scheme be used? (only used in MatVecT benchmark) */
   int omp_for_flag; /* are we using OpenMP for loops? */
   int MsgQ_flag; /* Are we using message queues instead of atomics? */
   int block_size;
   int hybrid_async_flag;
   int fine_grained_flag;
   int comp_wtime_flag;
   int MsgQ_wtime_flag;
   int MsgQ_cycles_flag;
   int comp_noop_flag;
   int MsgQ_noop_flag;
   int setup_type;
   int reduce_flag;
   int print_traces_flag;
}InputData;

typedef struct{
   double wtime;
   int num_spins;
   int num_atomics;
}RowOutputData;

/* Output data */
typedef struct{
   double solve_wtime; /* solve wall-clock time */
   double setup_wtime; /* setup wall-clock time */
   double *atomic_wtime_vec; /* atomic wall-clock time per thread */
   double *solve_wtime_vec; /* solve wall-clock time per thread */
   double *setup_wtime_vec;
   double *MsgQ_get_wtime_vec;
   double *MsgQ_put_wtime_vec;
   double *comp_wtime_vec;
   uint64_t *MsgQ_get_cycles_vec;
   uint64_t *MsgQ_put_cycles_vec;
   int *num_iters; /* number of iterations */
   int *num_relax; /* number of relaxations */
   int *num_qPuts_vec;
   int *num_qGets_vec;
   double *row_wtime;
   RowOutputData *row_output;
}OutputData;

/* Struct used by MatVecT benchmark */
typedef struct{
   InputData input;
   OutputData output;
   double *y_expand;
}MatVecData;

/* Struct used by AsyncJacobi benchmark */
typedef struct{
   InputData input;
   OutputData output;
}SolverData;

///* Level set data used in level scheduling algorithms */
//typedef struct{
//   vector<int> order; /* row ordering (rows are ordered by level)*/
//   vector<int> level_size; /* number of rows in each level */
//   vector<int> level_start; /* starting point in ``perm'' for each level */
//   int num_levels; /* number of level sets */
//}LevelSetData;

///* Struct used by TriSolve benchmark */
//typedef struct{
//   InputData input;
//   OutputData output;
//   LevelSetData L_lvl_set; /* level set data for lower triangular part */
//   LevelSetData U_lvl_set; /* level set data for upper triangular part */
//   PthreadData pt;
//}TriSolveData;

///* Struct used by ILU benchmark */
//typedef struct{
//   InputData input;
//   OutputData output;
//   LevelSetData L_lvl_set; /* level set data for lower triangular part */
//   LevelSetData U_lvl_set; /* level set data for upper triangular part */
//}ILUData;

/* Struct used by Histogram benchmark */
typedef struct{
   InputData input;
   OutputData output;
   int *Tally_expand;
}HistogramData;

/* Struct used by Histogram benchmark */
typedef struct{
   vector<int> phase;
   vector<int> source;
   vector<int> dest;
   vector<size_t> msg_size;
}TraceData;

#endif
