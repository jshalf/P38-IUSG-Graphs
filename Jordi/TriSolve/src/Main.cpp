#include "../../src/Main.hpp"
#include "TriSolve.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/Misc.hpp"
#include "slu_mt_ddefs.h"

using namespace std;

void SuperLU_MT_Setup(CSR A, SuperMatrix *A_slu, SuperMatrix *B_slu, SuperMatrix *L_slu, SuperMatrix *U_slu, SuperMatrix *AA_slu, SuperMatrix *AC_slu, double *b, int num_threads);
void SuperMatrix_to_CSR(SuperMatrix A_slu, CSR *A);

int main (int argc, char *argv[])
{
   TriSolveData ts;
   ts.input.num_threads = 1;
   ts.input.num_iters = 1;
   ts.input.atomic_flag = 0;
   ts.input.coo_flag = 0;
   int verbose_output = 0;
   int num_runs = 1;
   int m = 10; 
   int max_row_nnz = 3;

   int arg_index = 0;
   while (arg_index < argc){
      if (strcmp(argv[arg_index], "-n") == 0){
         arg_index++;
         m = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-max_iters") == 0){
         arg_index++;
         ts.input.num_iters = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_threads") == 0){
         arg_index++;
         ts.input.num_threads = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-atomic") == 0){
         ts.input.atomic_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-num_runs") == 0){
         arg_index++;
         num_runs = std::max(1, atoi(argv[arg_index]));
      }
      else if (strcmp(argv[arg_index], "-verb_out") == 0){
         verbose_output = 1;
      }
      else if (strcmp(argv[arg_index], "-coo") == 0){
         ts.input.coo_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-mxr_nnz") == 0){
         arg_index++;
         max_row_nnz = atoi(argv[arg_index]);
      }
      arg_index++;
   }
   
   omp_set_num_threads(ts.input.num_threads);
   srand(0);

   CSR A, L, U;
   //Laplace_2D_5pt(ts.input, &A, m);
   RandomMatrix(ts.input, &L, m, max_row_nnz, MATRIX_LOWER);
   RandomMatrix(ts.input, &U, m, max_row_nnz, MATRIX_UPPER);
   RandomMatrix(ts.input, &A, m, max_row_nnz, MATRIX_NONSYMMETRIC);
   LevelSets(L, &(ts.L_lvl_set));
   LevelSets(U, &(ts.U_lvl_set));
   int num_rows = A.n;
   int *L_perm = (int *)calloc(num_rows, sizeof(int));
   int *U_perm = (int *)calloc(num_rows, sizeof(int));
   for (int i = 0; i < num_rows; i++){
      L_perm[i] = i;
      U_perm[i] = num_rows - (i+1);
   }

   //char L_filename[100], U_filename[100], A_filename[100];
   //sprintf(L_filename, "./matlab/L.txt");
   //sprintf(U_filename, "./matlab/U.txt");
   //sprintf(A_filename, "./matlab/A.txt");
   //PrintCOO(L, L_filename);
   //PrintCOO(U, U_filename);
   //PrintCOO(A, A_filename);

   double *x = (double *)calloc(num_rows, sizeof(double));
   double *x_exact = (double *)calloc(num_rows, sizeof(double));
   double *y = (double *)calloc(num_rows, sizeof(double));
   double *y_exact = (double *)calloc(num_rows, sizeof(double));
   double *b = (double *)calloc(num_rows, sizeof(double));
   double *e_x = (double *)calloc(num_rows, sizeof(double));
   double *e_y = (double *)calloc(num_rows, sizeof(double));

   //SuperMatrix *A_slu, B_slu, L_slu, U_slu, *AA_slu, AC_slu;
   //SuperLU_MT_Setup(A, A_slu, &B_slu, &L_slu, &U_slu, AA_slu, &AC_slu, b, ts.input.num_threads);

   //SuperMatrix_to_CSR(L_slu, &L);

   srand(0);
   for (int i = 0; i < num_rows; i++){
      b[i] = 1.0;//RandDouble(-1.0, 1.0);
   } 
   for (int run = 1; run <= num_runs; run++){
      int iter;
      double start = omp_get_wtime();
      for (iter = 0; iter < ts.input.num_iters; iter++){
         for (int i = 0; i < num_rows; i++){
            x[i] = 0;
            y[i] = 0;
         } 
         TriSolve_CSR(&ts, L, U, L_perm, U_perm, x_exact, y_exact, b);

         //if (ts.input.coo_flag == 1){
         //   TriSolve_FineGrained_COO(&ts, L, U, x, y, b);
         //}
         //else {
            TriSolve_LevelSets_CSR(&ts, L, U, x, y, b);
         //}

         for (int i = 0; i < num_rows; i++){
            //printf("%e %e\n", y_exact[i], y[i]); 
            e_x[i] = x_exact[i] - x[i];
            e_y[i] = y_exact[i] - y[i];
         }
      }
      ts.output.solve_wtime = omp_get_wtime() - start;
      double error_x = sqrt(InnerProd(e_x, e_x, num_rows));
      double error_y = sqrt(InnerProd(e_y, e_y, num_rows));
      if (verbose_output){
         printf("TriSolve wall-clock time %e, L solve error L2-norm = %e, iterations = %d\n", ts.output.solve_wtime, error_x, iter);
      }
      else {
         printf("%e %e %e %d\n", ts.output.solve_wtime, error_x, error_y, iter);
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
   int_t panel_size = sp_ienv(1);
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
   printf("info %d\n", info);
   dgstrs(trans, L_slu, U_slu, perm_r, perm_c, B_slu, &Gstat, &info);
}

void SuperMatrix_to_CSR(SuperMatrix A_slu, CSR *A)
{
   if (A_slu.Stype == SLU_SCP){
      SCPformat *SCP = (SCPformat *)A_slu.Store;
      vector<vector<int>> colind_vec(A_slu.nrow, vector<int>());
      vector<vector<double>> nzval_vec(A_slu.nrow, vector<double>());
      double *nzval = (double *)SCP->nzval;
      for (int j = 0; j < A_slu.ncol; j++){
         printf("%d %d %d %d\n",
                SCP->rowind_colbeg[j], SCP->rowind_colend[j], SCP->nzval_colbeg[j], SCP->nzval_colend[j]);
         //for (int ii = SCP->colptr[j]; ii < SCP->colptr[j+1]; ii++){
         //   int row = SCP->rowind[ii];
         //   int col = j;
         //   double v = nzval[ii];
         //   colind_vec[row].push_back(col);
         //   nzval_vec[row].push_back(v);
         //   printf ("%d %d %e\n", row, col, v);
         //}
      }
      for (int j = 0; j < SCP->nnz; j++){
         printf("%d %e\n", SCP->rowind[j], nzval[j]);
      }
   }
   else if (A_slu.Stype == SLU_NR){
      ;
   }
}
