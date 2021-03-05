#include "../../src/Main.hpp"
#include "TriSolve.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/Misc.hpp"

#ifdef USE_SUPERLU
#include "slu_mt_ddefs.h"
#endif

using namespace std;

//TODO: fix SuperLU stuff
#ifdef USE_SUPERLU
void SuperLU_MT_Setup(CSR A, SuperMatrix *A_slu, SuperMatrix *B_slu, SuperMatrix *L_slu, SuperMatrix *U_slu, SuperMatrix *AA_slu, SuperMatrix *AC_slu, double *b, int num_threads);
void Print_SuperMatrix(SuperMatrix A_slu, CSR *A);
#endif

int main (int argc, char *argv[])
{
   /* set defaults */
   TriSolveData ts;
   ts.input.num_threads = 1;
   ts.input.atomic_flag = 1;
   ts.input.coo_flag = 0;
   ts.input.async_flag = 1;
   ts.input.omp_for_flag = 0;
   ts.input.MsgQ_flag = 0;
   int verbose_output = 0;
   int num_runs = 1;
   int m = 10; 
   int max_row_nnz = 3;
   int solver_type = TRISOLVE_ASYNC;
   int problem_type = PROBLEM_5PT_POISSON;
   double start;
   char mat_file_str[128];

   /* command line arguments */
   int arg_index = 0;
   int print_usage = 0;
   while (arg_index < argc){
      if (strcmp(argv[arg_index], "-problem") == 0){ /* test problem */
         arg_index++;
         if (strcmp(argv[arg_index], "5pt") == 0){ /* five-point centered-difference Laplace problem*/
            problem_type = PROBLEM_5PT_POISSON;
         }
         else if (strcmp(argv[arg_index], "rand") == 0){ /* random matrix */
            problem_type = PROBLEM_RANDOM;
         }
         else if (strcmp(argv[arg_index], "file") == 0){ /* read matrix from binary file */
            arg_index++;
            problem_type = PROBLEM_FILE;
            strcpy(mat_file_str, argv[arg_index]);
         }
      }
      else if (strcmp(argv[arg_index], "-solver") == 0){ /* solver type */
         arg_index++;
         if (strcmp(argv[arg_index], "async") == 0){ /* asynchronous iterative fine-grained */
            solver_type = TRISOLVE_ASYNC;
         }
         else if (strcmp(argv[arg_index], "lev_sched") == 0){ /* classical level-scheduled */
            solver_type = TRISOLVE_LEVEL_SCHEDULED;
         }
      }
      else if (strcmp(argv[arg_index], "-n") == 0){ /* ``size'' of matrix. n*n rows for Laplace, n rows otherwise. */
         arg_index++;
         m = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_threads") == 0){ /* number of threads */
         arg_index++;
         ts.input.num_threads = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-no_atomic") == 0){ /* use no-atomics implementations */
         ts.input.atomic_flag = 0;
      }
      else if (strcmp(argv[arg_index], "-num_runs") == 0){ /* number of TriSolve runs */
         arg_index++;
         num_runs = std::max(1, atoi(argv[arg_index]));
      }
      else if (strcmp(argv[arg_index], "-verb_out") == 0){ /* verbose output */
         verbose_output = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ") == 0){ /* use message queues in async solvers */
         ts.input.MsgQ_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-mxr_nnz") == 0){ /* max number of non-zeros per row (only used for generating random amtrices) */
         arg_index++;
         max_row_nnz = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-sync") == 0){ /* make the async solver synchronous */
         ts.input.async_flag = 0;
      }
      else if (strcmp(argv[arg_index], "-omp_for") == 0){ /* use OpenMP for loops in the async solver */
         ts.input.omp_for_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-help") == 0){ /* print command line options */
         print_usage = 1;
      }
      arg_index++;
   }

   if (print_usage == 1){
      printf("\n");
      printf("-problem <problem_name>:  test problem.\n");
      printf("      5pt:                five-point centered difference discretization of the Poisson equation.\n");
      printf("      rand:               random sparse matrix.\n");
      printf("      file:               read matrix from binary file.\n"
             "                          Must be in (i,j,val) format starting with (num rows, num cols, num nnz) on first line.\n");
      printf("-solver <solver_name>:    method for solving LUx=b.\n");
      printf("      lev_sched:          level scheduled method.\n");
      printf("      async:              asynchronous iterative fine-grained method.\n");
      printf("-n <int value>:           size of test problem.  For 5pt, this is the length of the 2D grid, i.e., the matrix has n^2 rows.\n");
      printf("                          for the random matrix, this is the number of rows.\n");
      printf("-num_threads <int value>: number of OpenMP threads.\n");
      printf("-no_atomic:               turn off atomics.  Only meant for performance measurements and will not produce a correct result.\n");
      printf("-num_runs <int value>:    number of independent runs.  Used for data collection.\n");
      printf("-verb_out:                verbose output.\n");
      printf("-MsgQ:                    use message queues instead of atomics.\n");
      printf("-mxr_nnz <int value>:     maximum number of non-zero values per row for the random matrix.\n");
      printf("-sync:                    use synchronization in the asynchronous method.\n");
      printf("-omp_for:                 use OpenMP parallel for loops in the asynchronous method.\n");
      printf("\n");
      return 0;
   }

   int num_rows;
   
   omp_set_num_threads(ts.input.num_threads);

   /* set up problem (TODO: add 5pt problem) */
   CSR A, L, U;
   //Laplace_2D_5pt(ts.input, &A, m);
   //num_rows = A.n;
   //double *b = (double *)calloc(num_rows, sizeof(double));
   //srand(0);
   //for (int i = 0; i < num_rows; i++){
   //   b[i] = 1.0;//RandDouble(-1.0, 1.0);
   //}   

   //SuperMatrix *A_slu, B_slu, L_slu, U_slu, *AA_slu, AC_slu;
   //SuperLU_MT_Setup(A, A_slu, &B_slu, &L_slu, &U_slu, AA_slu, &AC_slu, b, ts.input.num_threads);

   //Print_SuperMatrix(L_slu, &L);
   //char A_filename[100];
   //sprintf(A_filename, "./matlab/A.txt");
   //PrintCOO(A, A_filename);
   if (problem_type == PROBLEM_FILE){
      char L_mat_file_str[128];
      char U_mat_file_str[128];
      sprintf(L_mat_file_str, "%s_L.txt.bin", mat_file_str);
      sprintf(U_mat_file_str, "%s_U.txt.bin", mat_file_str);
      freadBinaryMatrix(L_mat_file_str, &L, 0);
      freadBinaryMatrix(U_mat_file_str, &U, 0);

      //char L_outfile[128];
      //char U_outfile[128];
      //sprintf(L_outfile, "./matlab/L.txt");
      //sprintf(U_outfile, "./matlab/U.txt");
      //PrintCOO(L, L_outfile, 0);
      //PrintCOO(U, U_outfile, 0);
   }
   else { 
      RandomMatrix(ts.input, &L, m, max_row_nnz, MATRIX_LOWER);
      RandomMatrix(ts.input, &U, m, max_row_nnz, MATRIX_UPPER);
   }
   
   ts.output.setup_wtime = 0.0;
   /* set up for level scheduling method */
   if (solver_type == TRISOLVE_LEVEL_SCHEDULED){
      start = omp_get_wtime();
      LevelSets(L, &(ts.L_lvl_set));
      LevelSets(U, &(ts.U_lvl_set));
      ts.output.setup_wtime = omp_get_wtime() - start;
   }
   num_rows = L.n;
   /* set up stuff for serial solver */
   int *L_perm = (int *)calloc(num_rows, sizeof(int));
   int *U_perm = (int *)calloc(num_rows, sizeof(int));
   for (int i = 0; i < num_rows; i++){
      L_perm[i] = i;
      U_perm[i] = num_rows - (i+1);
   }

//   char L_filename[100], U_filename[100], A_filename[100];
//   sprintf(L_filename, "./matlab/L.txt");
//   sprintf(U_filename, "./matlab/U.txt");
//   //sprintf(A_filename, "./matlab/A.txt");
//   PrintCOO(L, L_filename);
//   PrintCOO(U, U_filename);
//   //PrintCOO(A, A_filename);

   double *x = (double *)calloc(num_rows, sizeof(double));
   double *x_exact = (double *)calloc(num_rows, sizeof(double));
   double *y = (double *)calloc(num_rows, sizeof(double));
   double *y_exact = (double *)calloc(num_rows, sizeof(double));
   double *e_x = (double *)calloc(num_rows, sizeof(double));
   double *e_y = (double *)calloc(num_rows, sizeof(double));
   double *b = (double *)calloc(num_rows, sizeof(double));
   srand(0);
   for (int i = 0; i < num_rows; i++){
      b[i] = 1.0;//RandDouble(-1.0, 1.0);
   }

   ts.output.atomic_wtime_vec = (double *)calloc(ts.input.num_threads, sizeof(double));
   ts.output.solve_wtime_vec = (double *)calloc(ts.input.num_threads, sizeof(double));
   ts.output.num_relax = (int *)calloc(ts.input.num_threads, sizeof(int));

   for (int run = 1; run <= num_runs; run++){
      for (int i = 0; i < num_rows; i++){
         x[i] = 0;
         y[i] = 0;
         x_exact[i] = 0;
         y_exact[i] = 0;
      } 
      for (int t = 0; t < ts.input.num_threads; t++){
         ts.output.atomic_wtime_vec[t] = 0.0;
         ts.output.solve_wtime_vec[t] = 0.0;
         ts.output.num_relax[t] = 0;
      }
      /* serial solver first */
      TriSolve_CSR(&ts, L, U, L_perm, U_perm, x_exact, y_exact, b);

      /* parallel solver */
      if (solver_type == TRISOLVE_LEVEL_SCHEDULED){
         start = omp_get_wtime();
         TriSolve_LevelSets_CSR(&ts, L, U, x, y, b);
         ts.output.solve_wtime = omp_get_wtime() - start;
         for (int t = 0; t < ts.input.num_threads; t++){
            ts.output.num_relax[t] = ts.L_lvl_set.num_levels + ts.U_lvl_set.num_levels;
         }
      }
      else {
         TriSolve_FineGrained_COO(&ts, L, U, x, y, b);
         double solve_wtime_sum = SumDouble(ts.output.solve_wtime_vec, ts.input.num_threads);
         ts.output.solve_wtime = solve_wtime_sum / (double)ts.input.num_threads;
      }

      /* compute the error between the serial and parallel solvers */
      for (int i = 0; i < num_rows; i++){
         e_x[i] = x_exact[i] - x[i];
         e_y[i] = y_exact[i] - y[i];
         //printf("%e %e\n", y_exact[i], y[i]);
      }
      double error_x = sqrt(InnerProd(e_x, e_x, num_rows))/sqrt(InnerProd(x_exact, x_exact, num_rows));
      double error_y = sqrt(InnerProd(e_y, e_y, num_rows))/sqrt(InnerProd(y_exact, y_exact, num_rows));
      double atomic_wtime_sum = SumDouble(ts.output.atomic_wtime_vec, ts.input.num_threads);
      int num_relax_sum = SumInt(ts.output.num_relax, ts.input.num_threads);
      /* print output stats */
      if (verbose_output){
         printf("Solve wall-clock time = %e\nSetup wall-clock time = %e\nAtomics wall-clock time = %e\nL solve forward-error L2-norm = %e\nU solve forward-error L2-norm = %e\nmean relaxations = %f\n",
                 ts.output.solve_wtime, ts.output.setup_wtime, atomic_wtime_sum/(double)ts.input.num_threads, error_x, error_y, (double)num_relax_sum/(double)ts.input.num_threads);
      }
      else {
         printf("%e %e %e %e %e %f\n",
                ts.output.solve_wtime, ts.output.setup_wtime, atomic_wtime_sum/(double)ts.input.num_threads, error_x, error_y, (double)num_relax_sum/(double)ts.input.num_threads);
      }
   }

   free(x);
   free(x_exact);
   free(y);
   free(y_exact);
   free(b);
   free(e_x);
   free(e_y);

   return 0;
}

#ifdef USE_SUPERLU
void SuperLU_MT_Setup(CSR A, SuperMatrix *A_slu, SuperMatrix *B_slu, SuperMatrix *L_slu, SuperMatrix *U_slu, SuperMatrix *AA_slu, SuperMatrix *AC_slu, double *b, int num_threads)
{
   A_slu = (SuperMatrix *)SUPERLU_MALLOC(sizeof(SuperMatrix));
   //double *at;
   //int *rowind, *colptr;
   //dCompRow_to_CompCol(A.n, A.n, A.nnz, A.data, A.j, A.i_ptr, &at, &rowind, &colptr);
   dCreate_CompCol_Matrix(A_slu, A.n, A.n, A.nnz, A.data, A.j, A.i_ptr, SLU_NC, SLU_D, SLU_GE);
   int nrhs = 1;
   dCreate_Dense_Matrix(B_slu, A.n, nrhs, b, A.n, SLU_DN, SLU_S, SLU_GE);
   int_t *perm_r, *perm_c;
   perm_r = intMalloc(A.n);
   perm_c = intMalloc(A.n);
   int_t info;
   int_t permc_spec = 0;
   get_perm_c(permc_spec, A_slu, perm_c);
   superlumt_options_t superlumt_options;
   Gstat_t Gstat;
   int nprocs = num_threads;
   fact_t fact = DOFACT;
   yes_no_t refact = NO;
   trans_t trans = NOTRANS;
   int_t panel_size = 1;//sp_ienv(1);
   int_t relax = sp_ienv(2);
   double diag_pivot_thresh = 1.0;
   yes_no_t usepr = NO;
   double drop_tol = 0.0;
   void *work = NULL;
   int_t lwork = 0;

   //pdgssv(1, &A_slu, perm_c, perm_r, &L_slu, &U_slu, &B_slu, &info);
   StatAlloc(A.n, nprocs, panel_size, relax, &Gstat);
   StatInit(A.n, nprocs, &Gstat);
   if (A_slu->Stype == SLU_NR){
      NRformat *Astore = (NRformat *)A_slu->Store;
      AA_slu = (SuperMatrix *)SUPERLU_MALLOC(sizeof(SuperMatrix));
      dCreate_CompCol_Matrix(AA_slu, A_slu->ncol, A_slu->nrow, Astore->nnz,
                             (double *)Astore->nzval, Astore->colind, Astore->rowptr,
                             SLU_NC, A_slu->Dtype, A_slu->Mtype);
   }
   else if (A_slu->Stype == SLU_NC){
      AA_slu = A_slu;
   }
   pdgstrf_init(nprocs, fact, trans, refact, panel_size, relax,
                diag_pivot_thresh, usepr, drop_tol, perm_c, perm_r,
                work, lwork, AA_slu, AC_slu, &superlumt_options, &Gstat);
   pdgstrf(&superlumt_options, AC_slu, perm_r, L_slu, U_slu, &Gstat, &info);
   dgstrs(trans, L_slu, U_slu, perm_r, perm_c, B_slu, &Gstat, &info);
}

void Print_SuperMatrix(SuperMatrix A_slu, CSR *A)
{
   if (A_slu.Stype == SLU_SCP){
      SCPformat *Astore = (SCPformat *)A_slu.Store;
      vector<vector<int>> colind_vec(A_slu.nrow, vector<int>());
      vector<vector<double>> nzval_vec(A_slu.nrow, vector<double>());
      double *nzval = (double *)Astore->nzval;
      for (int k = 0; k <= Astore->nsuper; ++k) {
         int c = Astore->sup_to_colbeg[k];
         int nsup = Astore->sup_to_colend[k] - c;
         for (int j = c; j < c + nsup; j++) {
           int d = Astore->nzval_colbeg[j];
           for (int i = Astore->rowind_colbeg[c]; i < Astore->rowind_colend[c+1]; i++) {
              printf("%d %d %e\n", Astore->rowind[i], j, nzval[d++]);
           }
         }
      } 
   }
   else if (A_slu.Stype == SLU_NR){
   }
}
#endif
