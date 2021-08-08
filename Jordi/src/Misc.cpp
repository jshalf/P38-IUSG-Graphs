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

/* Inner product x'*y of two arrays of size n */
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
      for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
         int ii = A.j[jj];
         res -= A.data[jj] * x[ii];
      }
      r_2norm += res*res;
      /* compute right-hand side inner product */
      b_2norm += b[i]*b[i];
   }

   return sqrt(r_2norm)/sqrt(b_2norm);
}

/* Residual L2-norm using OpenMP reduction */
double Residual2Norm_CSC(Matrix A, /* sparse matrix data (input) */
                         double *x, /* solution (input) */
                         double *b /* right-hand side (input) */
                         )
{
   int n = A.n;
   double r_2norm = 0, b_2norm = 0;
   double *r = (double *)calloc(n, sizeof(double));
   #pragma omp parallel
   {
      #pragma omp for
      for (int i = 0; i < n; i++){
         r[i] = b[i];
      }
      #pragma omp for
      for (int i = 0; i < n; i++){
         double xi = x[i];
         for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
            int ii = A.i[jj];
            #pragma omp atomic
            r[ii] -= A.data[jj] * xi;
         }
      }

      #pragma omp for reduction(+:r_2norm,b_2norm)
      for (int i = 0; i < n; i++){
         r_2norm += r[i]*r[i];
         b_2norm += b[i]*b[i];
      }
   }
   free(r);

   return sqrt(r_2norm)/sqrt(b_2norm);
}

/* Compute level sets for level scheduling algorithms */
void LevelSets(InputData input,
               Matrix A, /* matrix data (input) */
               LevelSetData *lvl_set, /* level set data (output) */
               int L_flag /* is the matrix lower or upper triangular ? */
               )
{
   int n = A.n;
   int nnz = A.nnz;
   vector<int> depth(n, -1), level_size(n, 0);
   lvl_set->perm = (int *)calloc(n, sizeof(int));
   int lump = 1;

   /* loop over all nodes to compute depths which are levels */
   if (input.setup_type == LEVEL_SETS_ASYNC_SETUP){

      int *row_done_flags, *row_counts;
      if (input.mat_storage_type == MATRIX_STORAGE_CSC){
         row_counts = (int *)calloc(n, sizeof(int));
      }
      else {
         row_done_flags = (int *)calloc(n, sizeof(int));
      }

      #pragma omp parallel
      {
         int i_loc, j_loc, jj_loc, n_loc, nnz_loc;
         int *nz_done_flags_loc;

         n_loc = 0;
         nnz_loc = 0;
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < n; i++){
            n_loc++;
            for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
               nnz_loc++;
            }
         }

         if (input.mat_storage_type == MATRIX_STORAGE_CSC){
            #pragma omp for schedule(static, lump) nowait
            for (int k = 0; k < nnz; k++){
               #pragma omp atomic
               row_counts[A.i[k]]++;
            }
         }
         else {
            nz_done_flags_loc = (int *)calloc(nnz_loc, sizeof(int));
         }

         if (input.mat_storage_type == MATRIX_STORAGE_CSC){
            #pragma omp for schedule(static, lump) nowait
            for (int J = 0; J < n; J++){
               int j;
               if (L_flag == 1){
                  j = J;
               }
               else {
                  j = n-J-1;
               }

               int row_counts_j;
               /* stay idle until element j of x is ready to be used */
               do {
                  #pragma omp atomic read
                  row_counts_j = row_counts[j];
               } while (row_counts_j > 0);
               depth[j]++;
               for (int kk = A.start[j]; kk < A.start[j+1]; kk++){
                  int i = A.i[kk];
                  depth[i] = max(depth[i], depth[j]);
                  #pragma omp atomic
                  row_counts[i]--;
               }
            }
         }
         else {
            jj_loc = 0;
            #pragma omp for schedule(static, lump) nowait
            for (int I = 0; I < n; I++){
               int i;
               if (L_flag == 1){
                  i = I;
               }
               else {
                  i = n-I-1;
               }

               int max_depth = -1;
               int jj_start = A.start[i];
               int jj_end = A.start[i+1];
               int jj_diff = jj_end - jj_start;
               int row_count_i = jj_diff;
               while (row_count_i > 0){
                  int jj_loc_temp = jj_loc;
                  for (int jj = jj_start; jj < jj_end; jj++){
                     if (nz_done_flags_loc[jj_loc_temp] == 0){
                        int j = A.j[jj];
                        int row_done_flag_j;
                        #pragma omp atomic read
                        row_done_flag_j = row_done_flags[j];
                        if (row_done_flag_j == 1){
                           max_depth = max(max_depth, depth[j]);
                           row_count_i--;
                           nz_done_flags_loc[jj_loc_temp] = 1;
                        }
                     }
                     jj_loc_temp++;
                  }
               }
               depth[i] = 1 + max_depth;
               #pragma omp atomic write
               row_done_flags[i] = 1;
               jj_loc += jj_diff;
            }
         }

         if (input.mat_storage_type == MATRIX_STORAGE_CSC){

         }
         else {
            free(nz_done_flags_loc);
         }
      }

      if (input.mat_storage_type == MATRIX_STORAGE_CSC){
         free(row_counts);
      }
      else {
         free(row_done_flags);
      }
   }
   else {
      for (int I = 0; I < n; I++){
         int i;
         if (L_flag == 1){ 
            i = I;
         }
         else {
            i = n-I-1;
         }
         
         if (input.mat_storage_type == MATRIX_STORAGE_CSC){
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
               max_depth = max(max_depth, depth[j]);
            }
            depth[i] = 1 + max_depth;
         }
      }
   }


   lvl_set->num_levels = 0;
   for (int i = 0; i < n; i++){
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

void PrintDummy(double dummy)
{
   int tid = omp_get_thread_num();
   char dummy_filename[100];
   sprintf(dummy_filename, "dummy_%d_%jd.txt", tid, (intmax_t)time(NULL));
   FILE *dummy_file = fopen(dummy_filename, "w");
   fprintf(dummy_file, "%e\n", dummy);
   fclose(dummy_file);
   remove(dummy_filename);
}
