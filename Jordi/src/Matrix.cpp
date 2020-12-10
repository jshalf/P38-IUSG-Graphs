#include "Matrix.hpp"
#include "Misc.hpp"

using namespace std;

void Laplace_2D_5pt(InputData input, CSR *A, int n)
{
   int col;
   int N = n*n;

   A->n = N;
   A->m = N;
   A->nnz = 5*n*n - 4*n;
   A->i_ptr = (int *)calloc(N+1, sizeof(int));
   A->j = (int *)calloc(A->nnz, sizeof(int));
   A->data = (double *)calloc(A->nnz, sizeof(double));
   A->diag = (double *)calloc(A->n, sizeof(double));

   int block_end = n-1;
   int block_start = 0;
   int k = 0;
   A->i_ptr[0] = 0;
   for(int i = 0; i < N; i++){
      A->diag[i] = 4.0;
      col = i - n;
      if (col >= 0){
         A->data[k] = -1.0;
         A->j[k] = col;
         k++;
      }
      if (i > block_start){
         col = i - 1;
         A->data[k] = -1.0;
         A->j[k] = col;
         k++;
      }
      A->data[k] = 4.0;
      A->j[k] = i;
      k++;
      if (i < block_end){
         col = i + 1;
         A->data[k] = -1.0;
         A->j[k] = col;
         k++;
      }
      col = i + n;
      if (col < N){
         A->data[k] = -1.0;
         A->j[k] = col;
         k++;
      }
      A->i_ptr[i+1] = k;

      if (i == block_end){
         block_end += n;
         block_start += n;
      }
   }

   //if (input.coo_flag == 1){
      A->i = (int *)calloc(A->nnz, sizeof(int));
      k = 0;
      for (int i = 0; i < A->n; i++){
         for (int jj = A->i_ptr[i]; jj < A->i_ptr[i+1]; jj++){
            A->i[k] = i;
            k++;
         }
      }
   //}
}

void RandomMatrix(InputData input, CSR *A, int n, int max_row_nnz, int mat_type)
{
   A->diag = (double *)calloc(n, sizeof(double));
   double low = -1.0/(double)(max_row_nnz);
   double high = 1.0/(double)(max_row_nnz);

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
            int col;
            if (mat_type == MATRIX_LOWER){
               col = (int)RandInt(0, i, time(NULL));
            }
            else if (mat_type == MATRIX_UPPER){
               col = (int)RandInt(i, n-1, time(NULL));
            }
            else {
               col = (int)RandInt(0, n-1, time(NULL));
            }
            vector<int>::iterator it;

            it = find(rows[i].begin(), rows[i].end(), col);
            if (it == rows[i].end() && i != col){
               rows[i].push_back(col);
               nzval[i].push_back(RandDouble(0.0, high));
               count++;
               nnz++;
            }
         }
      }

      A->nnz = nnz;
      A->n = n;
      A->m = n;
      A->i = (int *)calloc(nnz, sizeof(int));
      A->j = (int *)calloc(nnz, sizeof(int));
      A->data = (double *)calloc(nnz, sizeof(double));
      A->i_ptr = (int *)calloc(n+1, sizeof(int));
      int k = 0;
      A->i_ptr[0] = 0;
      for (int i = 0; i < n; i++){
         A->i_ptr[i+1] = A->i_ptr[i] + rows[i].size();
         for (int jj = 0; jj < rows[i].size(); jj++){
            A->i[k] = i;
            A->j[k] = rows[i][jj];
            A->data[k] = nzval[i][jj];
            k++;
         }
      }
   }
}

void PrintCOO(CSR A, char *filename)
{
   FILE *file = fopen(filename, "w");
   for (int i = 0; i < A.n; i++){
      fprintf(file, "%d %d %.15e\n", i+1, i+1, A.diag[i]);
   }
   for (int k = 0; k < A.nnz; k++){
      fprintf(file, "%d %d %.15e\n", A.i[k]+1, A.j[k]+1, A.data[k]);
   }
   fclose(file);
}
