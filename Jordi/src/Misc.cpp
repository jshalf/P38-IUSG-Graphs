#include "Misc.hpp"
#include <sys/stat.h>
#include <sys/types.h>

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

///* Residual L2-norm using OpenMP reduction */
//double Residual2Norm(Matrix A, /* sparse matrix data (input) */
//                     double *x, /* solution (input) */
//                     double *b /* right-hand side (input) */
//                     )
//{
//   int n = A.n;
//   double r_2norm = 0, b_2norm = 0;
//   #pragma omp parallel for reduction(+:r_2norm,b_2norm)
//   for (int i = 0; i < n; i++){
//      /* compute residual inner product */
//      double res = b[i];
//      for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
//         int ii = A.j[jj];
//         res -= A.data[jj] * x[ii];
//      }
//      r_2norm += res*res;
//      /* compute right-hand side inner product */
//      b_2norm += b[i]*b[i];
//   }
//
//   return sqrt(r_2norm)/sqrt(b_2norm);
//}
//
///* Residual L2-norm using OpenMP reduction */
//double Residual2Norm_CSC(Matrix A, /* sparse matrix data (input) */
//                         double *x, /* solution (input) */
//                         double *b /* right-hand side (input) */
//                         )
//{
//   int n = A.n;
//   double r_2norm = 0, b_2norm = 0;
//   double *r = (double *)calloc(n, sizeof(double));
//   #pragma omp parallel
//   {
//      #pragma omp for
//      for (int i = 0; i < n; i++){
//         r[i] = b[i];
//      }
//      #pragma omp for
//      for (int i = 0; i < n; i++){
//         double xi = x[i];
//         for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
//            int ii = A.i[jj];
//            #pragma omp atomic
//            r[ii] -= A.data[jj] * xi;
//         }
//      }
//
//      #pragma omp for reduction(+:r_2norm,b_2norm)
//      for (int i = 0; i < n; i++){
//         r_2norm += r[i]*r[i];
//         b_2norm += b[i]*b[i];
//      }
//   }
//   free(r);
//
//   return sqrt(r_2norm)/sqrt(b_2norm);
//}

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
   mkdir("dummy_files", 0777);
   int tid = omp_get_thread_num();
   char dummy_filename[100];
   sprintf(dummy_filename, "dummy_files/dummy_%d_%jd.txt", tid, (intmax_t)time(NULL));
   FILE *dummy_file = fopen(dummy_filename, "w");
   fprintf(dummy_file, "%e\n", dummy);
   fclose(dummy_file);
   //remove(dummy_filename);
}

void PrintTraces(char *filename, TraceData trace_loc)
{
   int tid = omp_get_thread_num();
   int num_threads = omp_get_num_threads();
   for (int t = 0; t < num_threads; t++){
      if (t == tid){
         if (tid == 0) remove(filename);
         FILE *file = fopen(filename, "a");
         int num_traces = trace_loc.phase.size();
         for (int i = 0; i < num_traces; i++){
            fprintf(file, "%d %d %d %zu\n", 
                           trace_loc.phase[i],
                           trace_loc.source[i],
                           trace_loc.dest[i],
                           trace_loc.msg_size[i]);
         }
         fclose(file);
      }
      #pragma omp barrier
   }
}
