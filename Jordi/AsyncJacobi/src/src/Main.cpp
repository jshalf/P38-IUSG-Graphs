#include "Main.hpp"
#include "Misc.hpp"
#include "Matrix.hpp"
#include "Jacobi.hpp"

int main (int argc, char *argv[])
{
   SolverData solver;
   solver.input.solver_type = SYNC_JACOBI;
   solver.input.num_threads = 1;
   solver.input.num_iters = 50;
   solver.input.atomic_flag = 1;
   int verbose_output = 0;
   int num_runs = 1;
   int m = 10; 
   double w = 1.0;

   int arg_index = 0;
   while (arg_index < argc){
      if (strcmp(argv[arg_index], "-n") == 0){
         arg_index++;
         m = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-max_iters") == 0){
         arg_index++;
         solver.input.num_iters = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_threads") == 0){
         arg_index++;
         solver.input.num_threads = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-w") == 0){
         arg_index++;
         w = atof(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-solver") == 0){
         arg_index++;
         if (strcmp(argv[arg_index], "sync_jacobi") == 0){
            solver.input.solver_type = SYNC_JACOBI;
         }
         else if (strcmp(argv[arg_index], "async_jacobi") == 0){
            solver.input.solver_type = ASYNC_JACOBI;
         }
         else if (strcmp(argv[arg_index], "async_block_jacobi") == 0){
            solver.input.solver_type = ASYNC_BLOCK_JACOBI;
         }
      }
      else if (strcmp(argv[arg_index], "-atomic") == 0){
         arg_index++;
         solver.input.atomic_flag = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_runs") == 0){
         arg_index++;
         num_runs = max(1, atoi(argv[arg_index]));
      }
      else if (strcmp(argv[arg_index], "-verb_out") == 0){
         verbose_output = 1;
      }
      arg_index++;
   }
   
   omp_set_num_threads(solver.input.num_threads);

   CSR A;
   Laplace_2D_5pt(&A, m);
   int n = A.n;
   double *x = (double *)calloc(n, sizeof(double));
   double *b = (double *)calloc(n, sizeof(double));
 
   for (int run = 1; run <= num_runs; run++){
      srand(0);
      for (int i = 0; i < n; i++){
         b[i] = RandDouble(-1.0, 1.0);
         x[i] = 0;
      }
      Jacobi(&solver, A, b, &x);
      if (verbose_output){
         printf("Rel res. 2-norm %e, Solve wall-clock time %e\n", Residual2Norm(A, x, b), solver.output.solve_wtime);
      }
      else {
         printf("%e %e\n", Residual2Norm(A, x, b), solver.output.solve_wtime);
      }
   }

   free(x);
   free(b);

   return 0;
}