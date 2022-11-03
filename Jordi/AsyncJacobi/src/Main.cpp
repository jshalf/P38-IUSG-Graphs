#include "Main.hpp"
#include "Misc.hpp"
#include "Matrix.hpp"
#include "Jacobi.hpp"

double Residual2Norm(SparseMatrix A, double *x, double *b);
double Residual2Norm_CSC(SparseMatrix A, /* sparse matrix data (input) */
                         double *x, /* solution (input) */
                         double *b /* right-hand side (input) */
                         );

int main (int argc, char *argv[])
{
   /* matrix defaults */
   SparseMatrixInput spm_input;
   spm_input.mat_type = MatrixType::symm;
   spm_input.file_type = FileType::bin;
   spm_input.num_rows = 100;
   spm_input.num_cols = 100;
   spm_input.grid_len.resize(1);
   spm_input.grid_len[0] = 10;
   spm_input.store_type = SparseMatrixStorageType::CSR;
   spm_input.store_diag_in_vec = true;

   /* set defaults */
   SolverData solver;
   solver.input.solver_type = ASYNC_JACOBI;
   solver.input.num_threads = 1;
   solver.input.num_iters = 200;
   solver.input.atomic_flag = 1;
   solver.input.MsgQ_flag = 0;
   solver.input.comp_wtime_flag = 0;
   solver.input.MsgQ_wtime_flag = 0;
   solver.input.MsgQ_cycles_flag = 0;
   solver.input.comp_noop_flag = 0;
   solver.input.MsgQ_noop_flag = 0;
   int verbose_output = 0;
   int num_runs = 1;
   int m = 10; 
   double w = 1.0;
   int problem_type = PROBLEM_5PT_POISSON;
   char mat_file_str[128];

   int arg_index = 0;
   int print_usage = 0;
   while (arg_index < argc){
      if (strcmp(argv[arg_index], "-problem") == 0){ /* test problem */
         arg_index++;
         if (strcmp(argv[arg_index], "5pt") == 0){ /* five-point centered-difference Laplace problem */
            problem_type = PROBLEM_5PT_POISSON;
         }
         else if (strcmp(argv[arg_index], "file") == 0){ /* read matrix from binary file */
            arg_index++;
            problem_type = PROBLEM_FILE;
            strcpy(mat_file_str, argv[arg_index]);
         }
      }
      else if (strcmp(argv[arg_index], "-solver") == 0){ /* solver name */
         arg_index++;
         if (strcmp(argv[arg_index], "sj") == 0){ /* synchronous Jacobi */
            solver.input.solver_type = SYNC_JACOBI;
         }
         else if (strcmp(argv[arg_index], "aj") == 0){ /* asynchronous Jacobi */
            solver.input.solver_type = ASYNC_JACOBI;
         }
      }
      else if (strcmp(argv[arg_index], "-n") == 0){ /* size of matrix. n*n rows for Laplace, n rows otherwise */
         arg_index++;
         spm_input.grid_len[0] = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_iters") == 0){ /* max number of iterations */
         arg_index++;
         solver.input.num_iters = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_threads") == 0){ /* number of threads */
         arg_index++;
         solver.input.num_threads = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-w") == 0){ /* weight for weighsolver Jacobi (currently not used) */
         arg_index++;
         w = atof(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-no_atomic") == 0){ /* use atomics for async Jacobi */
         solver.input.atomic_flag = 0;
      }
      else if (strcmp(argv[arg_index], "-num_runs") == 0){ /* number of separate Jacobi runs */
         arg_index++;
         num_runs = std::max(1, atoi(argv[arg_index]));
      }
      else if (strcmp(argv[arg_index], "-verb_out") == 0){ /* verbose output */
         verbose_output = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ") == 0){ /* use message queues in async solvers */
         solver.input.MsgQ_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-help") == 0){ /* print command line options */
         print_usage = 1;
      }
      //else if (strcmp(argv[arg_index], "-sp_store_type") == 0){ /* matrix storage type */
      //   arg_index++;
      //   if (strcmp(argv[arg_index], "csr") == 0){ /* compressed sparse row */
      //      solver.input.mat_storage_type = MATRIX_STORAGE_CSR;
      //   }
      //   else if (strcmp(argv[arg_index], "csc") == 0){ /* compressed sparse column */
      //      solver.input.mat_storage_type = MATRIX_STORAGE_CSC;
      //   }
      //}
      else if (strcmp(argv[arg_index], "-MsgQ_wtime") == 0){
         solver.input.MsgQ_wtime_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-comp_wtime") == 0){
         solver.input.comp_wtime_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ_cycles") == 0){
         solver.input.MsgQ_cycles_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ_noop") == 0){
         solver.input.MsgQ_noop_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-comp_noop") == 0){
         solver.input.comp_noop_flag = 1;
      }
      arg_index++;
   }

   if (print_usage == 1){
      printf("\n");
      printf("-problem <problem_name>:  test problem.\n");
      printf("      5pt:                five-point centered difference discretization of the Poisson equation.\n");
      printf("      file:               read matrix from binary file.\n"
             "                          Must be in (i,j,val) format starting with (num rows, num cols, num nnz) on first line.\n");
      printf("-solver <solver_name>:    method for solving Ax=b.\n");
      printf("      sj:                 synchronous Jacobi.\n");
      printf("      aj:                 asynchronous Jacobi.\n");
      printf("-n <int value>:           size of test problem.  For 5pt, this is the length of the 2D grid, i.e., the matrix has n^2 rows.\n");
      printf("-num_threads <int value>: number of OpenMP threads.\n");
      printf("-no_atomic:               turn off atomics.  Only meant for performance measuremensolver and will not produce a correct result.\n");
      printf("-num_runs <int value>:    number of independent runs.  Used for data collection.\n");
      printf("-verb_out:                verbose output.\n");
      printf("-MsgQ:                    use message queues instead of atomics.\n");
      printf("-num_iters <int value>:   number of Jacobi iterations.\n");
      printf("\n");
      return 0;
   }

   if (solver.input.solver_type == SYNC_JACOBI){
      solver.input.MsgQ_flag = 0;
   }
   
   omp_set_num_threads(solver.input.num_threads);

   /* set up problem */
   SparseMatrix A(spm_input);

   if (problem_type == PROBLEM_FILE){
      A.ConstructMatrixFromFile();
   }
   else {
      A.ConstructLaplace2D5pt();
   }

   int n = A.GetNumRows();
   double *x = (double *)calloc(n, sizeof(double));
   double *b = (double *)calloc(n, sizeof(double));

   solver.output.solve_wtime_vec = (double *)calloc(solver.input.num_threads, sizeof(double));
   solver.output.num_qGets_vec = (int *)calloc(solver.input.num_threads, sizeof(int));
   solver.output.num_qPuts_vec = (int *)calloc(solver.input.num_threads, sizeof(int));

   if (solver.input.MsgQ_flag == 1){
      //if (solver.input.MsgQ_wtime_flag == 1){
         solver.output.MsgQ_put_wtime_vec = (double *)calloc(solver.input.num_threads, sizeof(double));
         solver.output.MsgQ_get_wtime_vec = (double *)calloc(solver.input.num_threads, sizeof(double));
      //}
      //else if (solver.input.MsgQ_cycles_flag == 1){
         solver.output.MsgQ_put_cycles_vec = (uint64_t *)calloc(solver.input.num_threads, sizeof(uint64_t));
         solver.output.MsgQ_get_cycles_vec = (uint64_t *)calloc(solver.input.num_threads, sizeof(uint64_t));
      //}
      //else if (solver.input.comp_wtime_flag == 1){
         solver.output.comp_wtime_vec = (double *)calloc(solver.input.num_threads, sizeof(double));
      //}
   }

   //int lump = 1;
   //#pragma omp parallel
   //{
   //   double dummy = 0.0;
   //   #pragma omp for schedule(static, lump) nowait
   //   for (int i = 0; i < n; i++){
   //      dummy += A.diag[i];
   //      for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
   //         dummy += A.j[jj] + A.data[jj];
   //      }
   //   }
   //   PrintDummy(dummy);
   //}

   for (int run = 1; run <= num_runs; run++){
      for (int t = 0; t < solver.input.num_threads; t++){
         solver.output.solve_wtime_vec[t] = 0.0;
      }
      srand(0);
      for (int i = 0; i < n; i++){
         b[i] = 1.0;//RandDouble(-1.0, 1.0);
         x[i] = 0;
      }
      /* run Jacobi solver */
      Jacobi(&solver, A, b, &x);

      for (int i = 0; i < n; i++){
         //printf("%f\n", x[i]);
      }

      double solve_wtime_sum = accumulate(solver.output.solve_wtime_vec, solver.output.solve_wtime_vec+solver.input.num_threads, (double)0.0);
      solver.output.solve_wtime = solve_wtime_sum / (double)solver.input.num_threads;

      double res_norm;
      res_norm = Residual2Norm(A, x, b);

      double MsgQ_put_wtime_sum = 0.0, MsgQ_put_wtime_mean = 0.0;
      double MsgQ_get_wtime_sum = 0.0, MsgQ_get_wtime_mean = 0.0;
      uint64_t uint64_t_accum_start = 0;
      uint64_t MsgQ_put_cycles_sum = 0;
      uint64_t MsgQ_get_cycles_sum = 0;
      double MsgQ_put_cycles_mean = 0.0;
      double MsgQ_get_cycles_mean = 0.0;
      double comp_wtime_sum = 0.0, comp_wtime_mean = 0.0;
      int num_qGets_sum = 0, num_qPuts_sum = 0;
      if (solver.input.MsgQ_flag == 1){
         //if (solver.input.MsgQ_wtime_flag == 1){
            MsgQ_put_wtime_sum = accumulate(solver.output.MsgQ_put_wtime_vec, solver.output.MsgQ_put_wtime_vec+solver.input.num_threads, (double)0.0);
            MsgQ_put_wtime_mean = MsgQ_put_wtime_sum / (double)solver.input.num_threads;

            MsgQ_get_wtime_sum = accumulate(solver.output.MsgQ_get_wtime_vec, solver.output.MsgQ_get_wtime_vec+solver.input.num_threads, (double)0.0);
            MsgQ_get_wtime_mean = MsgQ_get_wtime_sum / (double)solver.input.num_threads;
         //}
         //else if (solver.input.MsgQ_cycles_flag == 1){
            MsgQ_put_cycles_sum = accumulate(solver.output.MsgQ_put_cycles_vec, solver.output.MsgQ_put_cycles_vec+solver.input.num_threads, (uint64_t)0);
            MsgQ_put_cycles_mean = MsgQ_put_cycles_sum / (double)solver.input.num_threads;

            MsgQ_get_cycles_sum = accumulate(solver.output.MsgQ_get_cycles_vec, solver.output.MsgQ_get_cycles_vec+solver.input.num_threads, (uint64_t)0);
            MsgQ_get_cycles_mean = MsgQ_get_cycles_sum / (double)solver.input.num_threads;
         //}
         //else if (solver.input.comp_wtime_flag == 1){
            comp_wtime_sum = accumulate(solver.output.comp_wtime_vec, solver.output.comp_wtime_vec+solver.input.num_threads, (double)0.0);
            comp_wtime_mean = comp_wtime_sum / (double)solver.input.num_threads;
         //}

         num_qGets_sum = accumulate(solver.output.num_qGets_vec, solver.output.num_qGets_vec+solver.input.num_threads, (int)0);
         num_qPuts_sum = accumulate(solver.output.num_qPuts_vec, solver.output.num_qPuts_vec+solver.input.num_threads, (int)0);
      }

      /* print solver stats */
      if (verbose_output){
         printf("Rel res. 2-norm = %e\n"
                "Solve wall-clock time = %e\n"
                "Mean comp wtime = %e\n"
                "MsgQ put wtime = %e\n"
                "MsgQ get wtime = %e\n"
                "MsgQ put cycles = %" PRIu64 "\n"
                "MsgQ get cycles = %" PRIu64 "\n"
                "Num qPuts = %d\n"
                "Num qGets = %d\n",
                res_norm,
                solver.output.solve_wtime,
                comp_wtime_mean,
                MsgQ_put_wtime_mean,
                MsgQ_get_wtime_mean,
                MsgQ_put_cycles_sum,
                MsgQ_get_cycles_sum,
                num_qPuts_sum,
                num_qGets_sum);
      }
      else {
         printf("%e %e %e %e %e %" PRIu64 " %" PRIu64 " %d %d\n",
                res_norm,
                solver.output.solve_wtime,
                comp_wtime_mean,
                MsgQ_put_wtime_mean,
                MsgQ_get_wtime_mean,
                MsgQ_put_cycles_sum,
                MsgQ_get_cycles_sum,
                num_qPuts_sum,
                num_qGets_sum);
      }
   }

   free(x);
   free(b);

   free(solver.output.solve_wtime_vec);
   free(solver.output.num_qGets_vec);
   free(solver.output.num_qPuts_vec);

   if (solver.input.MsgQ_flag == 1){
      //if (solver.input.MsgQ_wtime_flag == 1){
         free(solver.output.MsgQ_put_wtime_vec);
         free(solver.output.MsgQ_get_wtime_vec);
      //}
      //else if (solver.input.MsgQ_cycles_flag == 1){
         free(solver.output.MsgQ_put_cycles_vec);
         free(solver.output.MsgQ_get_cycles_vec);
      //}
      //else if (solver.input.comp_wtime_flag == 1){
         free(solver.output.comp_wtime_vec);
      //}
   }

   return 0;
}


/* Residual L2-norm, i.e., ||b - Ax||_2, using OpenMP reduction */
double Residual2Norm(SparseMatrix A, /* sparse matrix data (input) */
                     double *x, /* solution (input) */
                     double *b /* right-hand side (input) */
                     )
{
   int n = A.GetNumRows();
   vector<int> start = A.GetIndexStarts();
   vector<int> col_idx = A.GetColIndices();
   vector<int> row_idx = A.GetRowIndices();
   vector<double> mat_values = A.GetValues();
   vector<double> diag = A.GetDiagValues();
   double r_2norm = 0, b_2norm = 0;
   if (A.GetStorageType() == SparseMatrixStorageType::CSC){
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
            for (int jj = start[i]; jj < start[i+1]; jj++){
               int ii = row_idx[jj];
               #pragma omp atomic
               r[ii] -= mat_values[jj] * xi;
            }
         }

         #pragma omp for reduction(+:r_2norm,b_2norm)
         for (int i = 0; i < n; i++){
            r_2norm += r[i]*r[i];
            b_2norm += b[i]*b[i];
         }
      }
      free(r);
   }
   else {
      #pragma omp parallel for reduction(+:r_2norm,b_2norm)
      for (int i = 0; i < n; i++){
         /* compute residual inner product */
         double res = b[i];
         for (int jj = start[i]; jj < start[i+1]; jj++){
            int ii = col_idx[jj];
            res -= mat_values[jj] * x[ii];
         }
         r_2norm += res*res;
         /* compute right-hand side inner product */
         b_2norm += b[i]*b[i];
      }
   }

   return sqrt(r_2norm)/sqrt(b_2norm);
}
