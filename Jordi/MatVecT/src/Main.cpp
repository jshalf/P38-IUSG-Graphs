#include "../../src/Main.hpp"
#include "MatVec.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/Misc.hpp"

int main (int argc, char *argv[])
{
   MatVecData mv;
   mv.input.num_threads = 1;
   mv.input.atomic_flag = 1;
   mv.input.AAT_flag = 0;
   mv.input.expand_flag = 0;
   mv.input.coo_flag = 0;
   mv.input.MsgQ_flag = 0;
   int verbose_output = 0;
   int num_runs = 1;
   int m = 10; 
   int problem_type = PROBLEM_5PT_POISSON;
   char mat_file_str[128];

   int print_usage = 0;
   int arg_index = 0;
   while (arg_index < argc){
      if (strcmp(argv[arg_index], "-problem") == 0){ /* test problem */
         arg_index++;
         if (strcmp(argv[arg_index], "5pt") == 0){ /* five-point centered-difference Laplace problem*/
            problem_type = PROBLEM_5PT_POISSON;
         }
         else if (strcmp(argv[arg_index], "file") == 0){ /* read matrix from binary file */
            arg_index++;
            problem_type = PROBLEM_FILE;
            strcpy(mat_file_str, argv[arg_index]);
         }
      }
      else if (strcmp(argv[arg_index], "-method") == 0){ /* method */
         arg_index++;
         if (strcmp(argv[arg_index], "atomic") == 0){ /* atomic method that treats Matrix as CSC */
            mv.input.expand_flag = 0;
         }
         else if (strcmp(argv[arg_index], "expand") == 0){ /* each thread uses a local vector which are summed at the end */
            mv.input.expand_flag = 1;
         }
      }
      else if (strcmp(argv[arg_index], "-n") == 0){ /* ``size'' of matrix. n*n rows for Laplace, n rows otherwise. */
         arg_index++;
         m = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_threads") == 0){ /* number of threads */
         arg_index++;
         mv.input.num_threads = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-no_atomic") == 0){ /* use no-atomics implementations */
         mv.input.atomic_flag = 0;
      }
      else if (strcmp(argv[arg_index], "-num_runs") == 0){ /* number of TriSolve runs */
         arg_index++;
         num_runs = std::max(1, atoi(argv[arg_index]));
      }
      else if (strcmp(argv[arg_index], "-verb_out") == 0){ /* verbose output */
         verbose_output = 1;
      }
      else if (strcmp(argv[arg_index], "-AAT") == 0){ /* Compute Ax and A^Tx together */
         mv.input.AAT_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-coo") == 0){ /* use coordinate format */
         mv.input.coo_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ") == 0){ /* use message queues */
         mv.input.MsgQ_flag = 1;
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
      printf("      file:               read matrix from binary file.\n"
             "                          Must be in (i,j,val) format starting with (num rows, num cols, num nnz) on first line.\n");
      printf("-method <method_name>:    method for computing y=A^Tx.\n");
      printf("      atomic:             method that uses atomics.\n");
      printf("      expand:             baseline method that does not require atomics.\n");
      printf("-n <int value>:           size of test problem.  For 5pt, this is the length of the 2D grid, i.e., the matrix has n^2 rows.\n");
      printf("-num_threads <int value>: number of OpenMP threads.\n");
      printf("-no_atomic:               turn off atomics.  Only meant for performance measurements and will not produce a correct result.\n");
      printf("-num_runs <int value>:    number of independent runs.  Used for data collection.\n");
      printf("-verb_out:                verbose output.\n");
      printf("-AAT:                     compute both A^Tx and Ax together.\n");
      printf("-coo:                     use coordinate format instead of Matrix.\n");
      printf("-MsgQ:                    use message queues instead of atomics.\n");
      printf("\n");
      return 0;
   }
   
   omp_set_num_threads(mv.input.num_threads);

   /* set up problem */
   Matrix A;
   if (problem_type == PROBLEM_FILE){
      char A_mat_file_str[128];
      sprintf(A_mat_file_str, "%s_A.txt.bin", mat_file_str);
      freadBinaryMatrix(A_mat_file_str, &A, 1, 0, 1);

      //char A_outfile[128];
      //sprintf(A_outfile, "./matlab/A.txt");
      //PrintCOO(A, A_outfile, 0);
   }
   else {
      Laplace_2D_5pt(mv.input, &A, m);
   }
   int num_rows = A.n;
   int num_cols = A.m;

   if (num_rows != num_cols){
      mv.input.expand_flag = 0;
   }

   double *x = (double *)calloc(num_rows, sizeof(double));
   double *y1 = (double *)calloc(num_cols, sizeof(double));
   double *y1_exact = (double *)calloc(num_cols, sizeof(double));
   double *e1 = (double *)calloc(num_cols, sizeof(double));

   double *y2, *y2_exact, *e2;
   if (mv.input.AAT_flag == 1){
      y2 = (double *)calloc(num_rows, sizeof(double));
      y2_exact = (double *)calloc(num_rows, sizeof(double));
      e2 = (double *)calloc(num_rows, sizeof(double));
   }

   if (mv.input.expand_flag == 1){
      mv.y1_expand = (double *)calloc(mv.input.num_threads * num_cols, sizeof(double));
      mv.y2_expand = (double *)calloc(mv.input.num_threads * num_rows, sizeof(double));
   }

   srand(0);
   for (int i = 0; i < num_rows; i++){
      x[i] = RandDouble(-1.0, 1.0);
   } 
   for (int run = 1; run <= num_runs; run++){
      for (int i = 0; i < num_cols; i++){
         y1[i] = 0;
         if (mv.input.AAT_flag == 1){
            y2[i] = 0;
         }
      }
      /* serial MatVec (since we are only considering aymmetric matrices right now, Ax = A^Tx).
       * TODO: handle non-symmetric case.  */ 
      MatVec_CSR(&mv, A, x, y1_exact);

      /* parallel MatVecT */
      double start = omp_get_wtime();
      if (mv.input.coo_flag == 1){
         MatVecT_COO(&mv, A, x, y1, y2);
      }
      else {
         MatVecT_CSR(&mv, A, x, y1, y2);
      }
      mv.output.solve_wtime = omp_get_wtime() - start;

      /* compute error */
      for (int i = 0; i < num_cols; i++){
         e1[i] = y1_exact[i] - y1[i];
         //printf("%e %e\n", y1_exact[i], y1[i]);
      }
      if (mv.input.AAT_flag == 1){
         for (int i = 0; i < num_rows; i++){
            y2_exact[i] = y1_exact[i];
            e2[i] = y2_exact[i] - y2[i];
         }
      }
      double error1 = sqrt(InnerProd(e1, e1, num_cols))/sqrt(InnerProd(y1_exact, y1_exact, num_rows));
      double error2 = 0.0;
      if (mv.input.AAT_flag == 1) error2 = sqrt(InnerProd(e2, e2, num_rows))/sqrt(InnerProd(y2_exact, y2_exact, num_rows));
      /* print stats */
      if (verbose_output){
         printf("MatVec wall-clock time %e, AT error L2-norm %e, A error L2-norm = %e\n", mv.output.solve_wtime, error1, error2);
      }
      else {
         printf("%e %e %e\n", mv.output.solve_wtime, error1, error2);
      }
   }

   free(x);
   free(y1);
   free(y1_exact);
   free(e1);

   if (mv.input.AAT_flag == 1){
      free(y2);
      free(y2_exact);
      free(e2);
   }

   if (mv.input.expand_flag == 1){
      free(mv.y1_expand);
      free(mv.y2_expand);
   }

   return 0;
}
