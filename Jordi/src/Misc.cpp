#include "Misc.hpp"

using namespace std;

/* random double in the range [low,high) (naive)*/
double RandDouble(double low, double high)
{
   double x;

   //uniform_real_distribution<double> unif(low, high);
   //default_random_engine re;
   //x = unif(re);

   x = low + (high - low) * ((double)rand() / RAND_MAX);

   return x;
}

/* Random integer in the range [low,high] (naive)*/
int RandInt(int low, int high, double seed)
{
   int x;

   //random_device rd;
   //mt19937 gen(rd()); // seed the generator
   //uniform_int_distribution<> distr(low, high); // define the range
   //return distr(gen);

   x = (int)round(RandDouble((double)low, (double)high));

   return x;
}

/* Inner product x'*y of two vectors of size ``n'' */
double InnerProd(double *x, double *y, int n)
{
   double inner_prod = 0;
   #pragma omp parallel for reduction(+:inner_prod)
   for (int i = 0; i < n; i++){
      inner_prod += x[i] * y[i];
   }

   return inner_prod;
}

/* Sum of array ``x'' of ints of length ``n'' */
int SumInt(int *x, int n)
{
   int sum = 0;
   for (int i = 0; i < n; i++){
      sum += x[i];
   }

   return sum;
}

/* Sum of array ``x'' of doubles of length ``n'' */
double SumDouble(double *x, int n)
{
   double sum = 0;
   for (int i = 0; i < n; i++){
      sum += x[i];
   }

   return sum;
}

/* Residual L2-norm using OpenMP reduction */
double Residual2Norm(Matrix A, /* sparse matrix data (input) */
                     double *x, /* solution (input) */
                     double *b /* right-hand side (input) */
                     )
{
   int n = A.n;
   double r_2norm = 0, b_2norm = 0;
   #pragma omp parallel for reduction(+:r_2norm,b_2norm)
   for (int i = 0; i < n; i++){
      /* compute residual inner product */
      double res = b[i];
      for (int jj = A.start[i]; jj < A.start[i+1]; jj++)
      {
         int ii = A.j[jj];
         res -= A.data[jj] * x[ii];
      }
      r_2norm += res*res;
      /* compute right-hand side inner product */
      b_2norm += b[i]*b[i];
   }

   return sqrt(r_2norm)/sqrt(b_2norm);
}

/* Compute level sets for level scheduling algorithms */
void LevelSets(Matrix A, /* matrix data (input) */
              LevelSetData *lvl_set, /* level set data (output) */
              int L_flag, /* is the matrix lower or upper triangular ? */
              int csc_flag /* is the matrix in CSC format? */
              )
{
   int n = A.n;
   vector<int> depth(n, -1), level_size(n, 0);
   lvl_set->perm = (int *)calloc(n, sizeof(int));
   lvl_set->num_levels = 0;
   /* loop over all nodes to compute depths which are levels */
   for (int I = 0; I < n; I++){
      int i;
      if (L_flag == 1){ 
         i = I;
      }
      else {
         i = n-I-1;
      }
      int max_depth;
      
      if (csc_flag == 1){
         depth[i]++;
         for (int kk = A.start[i]; kk < A.start[i+1]; kk++){
            int j = A.i[kk];
            depth[j] = max(depth[j], depth[i]);
         }
      }
      else {
         int max_depth = -1;
         for (int kk = A.start[i]; kk < A.start[i+1]; kk++){
            int j = A.j[kk];
            if (max_depth < depth[j]){
               max_depth = depth[j];
            }
         }
         depth[i] = 1 + max_depth;
      }

      level_size[depth[i]]++;
      if (lvl_set->num_levels < depth[i]){
         lvl_set->num_levels = depth[i];
      }
   }
   lvl_set->num_levels++;

   lvl_set->level_size.resize(lvl_set->num_levels); 
   lvl_set->level_start.resize(lvl_set->num_levels+1);
   lvl_set->level_start[0] = 0;
   /* compute the starting points in ``perm'' of each level set */
   for (int i = 0; i < lvl_set->num_levels; i++){
      lvl_set->level_size[i] = level_size[i];
      lvl_set->level_start[i+1] = lvl_set->level_start[i] + lvl_set->level_size[i];
   }

   vector<int> level_counts(lvl_set->num_levels, 0);
   for (int i = 0; i < n; i++){
      int d = depth[i];
      lvl_set->perm[lvl_set->level_start[d]+level_counts[d]] = i;
      level_counts[d]++;
   }
}

void LevelSetsDestroy(LevelSetData *lvl_set)
{
   free(lvl_set->perm);
}

/* Functions for computing number of clock cycles */
#ifdef WINDOWS
#include <intrin.h>
uint64_t rdtsc(){
    return __rdtsc();
}
#else
#include <x86intrin.h>
uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}
#endif
