/*
 * functions related to generating test matrices
 */

#include "Matrix.hpp"
#include "Misc.hpp"

using namespace std;

/* five-point centered difference discretization of the Laplace equation */
void Laplace_2D_5pt(InputData input, /* input data */
                    Matrix *A, /* sparse matrix data (output) */
                    int n /* size of grid (n*n rows) */
                    )
{
   int col;
   int N = n*n;

   A->n = N;
   A->m = N;
   A->nnz = 5*n*n - 4*n;
   A->start = (int *)calloc(N+1, sizeof(int));
   if (input.coo_flag == 1){
      A->i = (int *)calloc(A->nnz, sizeof(int));
      A->j = (int *)calloc(A->nnz, sizeof(int));
   }
   else {
      if (input.mat_storage_type == MATRIX_STORAGE_CSC){   
         A->i = (int *)calloc(A->nnz, sizeof(int));
      }
      else {
         A->j = (int *)calloc(A->nnz, sizeof(int));
      }
   }
   A->data = (double *)calloc(A->nnz, sizeof(double));
   A->diag = (double *)calloc(A->n, sizeof(double));

   int *idx = (int *)calloc(A->nnz, sizeof(int));

   int block_end = n-1;
   int block_start = 0;
   int k = 0;
   A->start[0] = 0;
   for(int i = 0; i < N; i++){
      A->diag[i] = 4.0;
      col = i - n;
      if (col >= 0){
         A->data[k] = -1.0;
         idx[k] = col;
         k++;
      }
      if (i > block_start){
         col = i - 1;
         A->data[k] = -1.0;
         idx[k] = col;
         k++;
      }
      A->data[k] = 4.0;
      idx[k] = i;
      k++;
      if (i < block_end){
         col = i + 1;
         A->data[k] = -1.0;
         idx[k] = col;
         k++;
      }
      col = i + n;
      if (col < N){
         A->data[k] = -1.0;
         idx[k] = col;
         k++;
      }
      A->start[i+1] = k;

      if (i == block_end){
         block_end += n;
         block_start += n;
      }
   }

   for (int k = 0; k < A->nnz; k++){
      if (input.mat_storage_type == MATRIX_STORAGE_CSC){
         A->i[k] = idx[k];
      }
      else {
         A->j[k] = idx[k];
      }
   }

   if (input.coo_flag == 1){
      k = 0;
      for (int i = 0; i < A->n; i++){
         for (int jj = A->start[i]; jj < A->start[i+1]; jj++){
            if (input.mat_storage_type == MATRIX_STORAGE_CSC){
               A->j[k] = i;
            }
            else {
               A->i[k] = i;
            }
            k++;
         }
      }
   }
}

/* Generate a random sparse matrix */
void RandomMatrix(InputData input, /* input data */
                  Matrix *A, /* matrix data (output) */
                  int n, /* number of rows */
                  int max_row_nnz, /* maximum number of non-zeros per row */ 
                  int mat_type, /* matrix type can be symmetric, lower triangular, or upper triangular */
                  int csc_flag,
                  int coo_flag
                  )
{
   A->diag = (double *)calloc(n, sizeof(double));
   double low = -1.0/(double)(max_row_nnz);
   double high = 1.0/(double)(max_row_nnz);

   //srand(time(NULL));
   srand(0);
   if (mat_type == MATRIX_LOWER || 
       mat_type == MATRIX_UPPER ||
       mat_type == MATRIX_NONSYMMETRIC){
      vector<vector<int>> rows(n, vector<int>()); 
      vector<vector<double>> nzval(n, vector<double>());
      int nnz = 0;
      for (int i = 0; i < n; i++){
         int i_max_row_nnz;
         if (mat_type == MATRIX_LOWER){
            i_max_row_nnz = min(i+1, max_row_nnz);
         }
         else if (mat_type == MATRIX_UPPER){
            i_max_row_nnz = min(n-i, max_row_nnz);
         }
         else {
            i_max_row_nnz = max_row_nnz;
         }
         int row_nnz = (int)RandDouble(1, i_max_row_nnz);
         int count = 1;
         A->diag[i] = RandDouble(0.0, 1.0);
         //rows[i].push_back(i);
         //nzval[i].push_back(RandDouble(-1.0, 1.0));
         //nnz++;
         while(count < row_nnz){
            int j;
            if (mat_type == MATRIX_LOWER){
               j = (int)RandInt(0, i, time(NULL));
            }
            else if (mat_type == MATRIX_UPPER){
               j = (int)RandInt(i, n-1, time(NULL));
            }
            else {
               j = (int)RandInt(0, n-1, time(NULL));
            }
            vector<int>::iterator it;

            it = find(rows[i].begin(), rows[i].end(), j);
            if (it == rows[i].end() && i != j){
               rows[i].push_back(j);
               nzval[i].push_back(RandDouble(0.0, high));
               count++;
               nnz++;
            }
         }
      }

      A->nnz = nnz;
      A->n = n;
      A->m = n;
      if (csc_flag == 1 || coo_flag == 1){
         A->i = (int *)calloc(A->nnz, sizeof(int));
      }
      if (csc_flag == 0 || coo_flag == 1){
         A->j = (int *)calloc(A->nnz, sizeof(int));
      }
      A->data = (double *)calloc(A->nnz, sizeof(double));
      A->start = (int *)calloc(A->n+1, sizeof(int));

      vector<int> counts(A->n, 0);

      for (int i = 0; i < n; i++){
         for (int k = 0; k < rows[i].size(); k++){
            int row = i;
            int col = rows[i][k];
            int idx;
            if (csc_flag == 1){
               idx = col;
            }
            else {
               idx = row;
            }
            counts[idx]++;
         }
      }

      for (int idx = 0; idx < A->n; idx++){
         A->start[idx+1] = A->start[idx] + counts[idx];
         counts[idx] = 0;
      }

      for (int i = 0; i < n; i++){
         for (int k = 0; k < rows[i].size(); k++){
            int row = i;
            int col = rows[i][k];
            double elem = nzval[i][k];

            int idx;
            if (csc_flag == 1){
               idx = col;
            }
            else {
               idx = row;
            }
            int kk = A->start[idx] + counts[idx];
            counts[idx]++;

            if (coo_flag == 1){
               A->i[kk] = row;
               A->j[kk] = col;
            }
            else {
               if (csc_flag == 1){
                  A->i[kk] = row;
               }
               else {
                  A->j[kk] = col;
               }
            }
            A->data[kk] = elem;
         }
      }
   }
}

/* print coordinate format of sparse matrix */
void PrintMatrix(Matrix A, char *filename, int print_diag_flag, int csc_flag)
{
   FILE *file = fopen(filename, "w");
   
   if (print_diag_flag == 1){
      for (int i = 0; i < A.n; i++){
         fprintf(file, "%d %d %.15e\n", i+1, i+1, A.diag[i]);
      }
   }
   if (csc_flag == 1){
      for (int j = 0; j < A.n; j++){
         for (int kk = A.start[j]; kk < A.start[j+1]; kk++){
            int row = A.i[kk]+1;
            int col = j+1;
            fprintf(file, "%d %d %.15e\n", row, col, A.data[kk]);
         }
      }
   }
   else {
      for (int i = 0; i < A.n; i++){
         for (int kk = A.start[i]; kk < A.start[i+1]; kk++){
            int row = i+1;
            int col = A.j[kk]+1;
            fprintf(file, "%d %d %.15e\n", row, col, A.data[kk]);
         }
      }
   }
   
   fclose(file);
}

/* read matrix from binary file. matrix entries must be ordered by increasing row index then increasing column index */
void freadBinaryMatrix(char *mat_file_str,
                       Matrix *A,
                       int include_diag_flag,
                       int csc_flag,
                       int coo_flag)
{
   size_t size;
   int temp_size;
   Triplet_AOS *buffer;

   FILE *mat_file = fopen(mat_file_str, "rb");

   fseek(mat_file, 0, SEEK_END);
   size = ftell(mat_file);
   rewind(mat_file);
   buffer = (Triplet_AOS *)malloc(sizeof(Triplet_AOS) * size);
   fread(buffer, sizeof(Triplet_AOS), size, mat_file);

   int file_lines = size/sizeof(Triplet_AOS);
   int num_rows = (int)buffer[0].i;
   int nnz = (int)buffer[0].j;

   A->n = num_rows;
   A->m = A->n;
   A->nnz = nnz;
   if (include_diag_flag == 0){
      A->nnz -= A->n;
   }

   vector<int> counts(A->n, 0);

   if (csc_flag == 1 || coo_flag == 1){
      A->i = (int *)calloc(A->nnz, sizeof(int));
   }
   if (csc_flag == 0 || coo_flag == 1){
      A->j = (int *)calloc(A->nnz, sizeof(int));
   }
   A->data = (double *)calloc(A->nnz, sizeof(double));
   A->start = (int *)calloc(A->n+1, sizeof(int));
   A->diag = (double *)calloc(A->n, sizeof(double));

   for (int k = 1; k < file_lines; k++){
      int row = buffer[k].i-1;
      int col = buffer[k].j-1;
      
      int include_elem_flag = 1;
      if (row == col){
         if (include_diag_flag == 0){
            include_elem_flag = 0;
         }
      }

      if (include_elem_flag == 1){
         int idx;
         if (csc_flag == 1){
            idx = col;
         }
         else {
            idx = row;
         }
         counts[idx]++;
      }
   }

   for (int idx = 0; idx < A->n; idx++){
      A->start[idx+1] = A->start[idx] + counts[idx];
      counts[idx] = 0;
   }

   for (int k = 1; k < file_lines; k++){
      int row = buffer[k].i-1;
      int col = buffer[k].j-1;
      double elem = buffer[k].val;
    
 
      int include_elem_flag = 1;
      if (row == col){
         A->diag[row] = elem;
         if (include_diag_flag == 0){
            include_elem_flag = 0;
         }
      }
     
      if (include_elem_flag == 1) {   
         int kk, idx;
         if (csc_flag == 1){
            idx = col;
         }
         else {
            idx = row;
         }
         kk = A->start[idx] + counts[idx];
         counts[idx]++;
         if (coo_flag == 1){
            A->i[kk] = row;
            A->j[kk] = col; 
         }
         else {
            if (csc_flag == 1){
               A->i[kk] = row;
            }
            else {
               A->j[kk] = col;
            }
         }
         A->data[kk] = elem;
      }
   }

   fclose(mat_file);
}
