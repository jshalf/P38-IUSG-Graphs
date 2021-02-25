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
double Residual2Norm(CSR A, /* sparse matrix data (input) */
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
      for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++)
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
void LevelSets(CSR A, /* matrix data (input) */
               LevelSetData *lvl_set /* level set data (output) */
               )
{
   int n = A.n;
   vector<int> nodes(n);
   vector<int>::iterator it;
   lvl_set->perm = (int *)calloc(n, sizeof(int));
   int *nonzero_flags = (int *)calloc(n, sizeof(int));
   int *nonzero_flags_prev = (int *)calloc(n, sizeof(int));

   for (int i = 0; i < n; i++){
      nonzero_flags[i] = 1;
      nodes[i] = i;
   }

   int k = 0;
   /* loop until all nodes have been eliminated, i.e., until all nodes have been placed in a level set */
   while(nodes.size() > 0){
      int level_size = 0;
      for (int i = 0; i < n; i++){
         nonzero_flags_prev[i] = nonzero_flags[i];
      }
      /* loop over nodes to find roots */
      for(it = nodes.begin(); it != nodes.end(); it++){
         int i = *it;
         int root_flag = 1;
         /* if all nodes connected to node ``i'' have been eliminated, ``i'' is a root */
         for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
            if (nonzero_flags_prev[A.j[jj]] == 1){
               root_flag = 0;
               break;
            }
         }
         /* if ``i'' is a root, include it in current level set and remove it from list of nodes */ 
         if (root_flag == 1){
            nonzero_flags[i] = 0;
            lvl_set->perm[k] = i;
            k++;
            level_size++;
            nodes.erase(it--);
         }
      }
      lvl_set->level_size.push_back(level_size);
   }

   lvl_set->num_levels = lvl_set->level_size.size();
   lvl_set->level_start.resize(lvl_set->num_levels+1);
   lvl_set->level_start[0] = 0;
   /* compute the starting points in ``perm'' of each level set */
   for (int i = 0; i < lvl_set->num_levels; i++){
      lvl_set->level_start[i+1] = lvl_set->level_start[i] + lvl_set->level_size[i];
   } 

   free(nonzero_flags);
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
