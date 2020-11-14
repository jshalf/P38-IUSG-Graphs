#ifndef MAIN_HPP
#define MAIN_HPP

#define SYNC_JACOBI  0
#define ASYNC_JACOBI 1
#define SYNC_BLOCK_JACOBI  2
#define ASYNC_BLOCK_JACOBI 3

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include<iostream> 
#include<algorithm>
#include <omp.h>

typedef struct{
   int num_threads;
   int num_iters;
   int async_flag;
   int atomic_flag;
   int AAT_flag;
   int expand_flag;
   int coo_flag;
}InputData;

typedef struct{
   double solve_wtime;
}OutputData;

typedef struct{
   InputData input;
   OutputData output;
   double *y1_expand;
   double *y2_expand;
}MatVecData;

using namespace std;

#endif
