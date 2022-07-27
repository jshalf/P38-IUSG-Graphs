#include "Matrix.hpp"
#include "Misc.hpp"

/************************************
 * Member functions for Matrix class
 ************************************/

int Matrix::GetNumRows(void)
{
   return num_rows;
}

int Matrix::GetNumCols(void)
{
   return num_cols;
}




/******************************************
 * Member functions for SparseMatrix class
 ******************************************/

SparseMatrix::SparseMatrix(SparseMatrixInput input)
{
   num_rows = input.num_rows;
   num_cols = input.num_cols;
   grid_len = input.grid_len;
   max_row_nnz = input.max_row_nnz;
   input_file_name = input.file_name;
   mat_type = input.mat_type;
   store_type = input.store_type;
   input_file_type = input.file_type;
   store_diag_in_vec = input.store_diag_in_vec;
}

SparseMatrix::~SparseMatrix(void)
{

}

vector<double> SparseMatrix::GetValues(void)
{
   return data;
}

vector<double> SparseMatrix::GetDiagValues(void)
{
   return diag;
}

vector<int> SparseMatrix::GetColIndices(void)
{
   return col_idx;
}

vector<int> SparseMatrix::GetRowIndices(void)
{
   return row_idx;
}

vector<int> SparseMatrix::GetIndexStarts(void)
{
   return start_idx;
}

int SparseMatrix::GetNNZ(void)
{
   return nnz;
}

SparseMatrixStorageType SparseMatrix::GetStorageType(void)
{
   return store_type;
}

void SparseMatrix::ConstructMatrixFromFile(void)
{
   FreadBinaryMatrix();
}

/* five-point centered difference discretization of the Laplace equation */
void SparseMatrix::ConstructLaplace2D5pt(void)
{
   int col;
   int nx = grid_len[0];
   int N = nx * nx;

   num_rows = N;
   num_rows = N;
   nnz = 5*nx*nx - 4*nx;
   start_idx.resize(N);

   switch(store_type){
      case SparseMatrixStorageType::COO:
         row_idx.resize(nnz);
         col_idx.resize(nnz);
         break;
      case SparseMatrixStorageType::CSC:
         row_idx.resize(nnz);
         break;
      case SparseMatrixStorageType::CSR:
         col_idx.resize(nnz);
         break;
      default:
         ;
   }

   data.resize(nnz);
   diag.resize(num_rows);

   vector<int> idx;
   idx.resize(nnz);

   int block_end = nx-1;
   int block_start = 0;
   int k = 0;
   start_idx[0] = 0;
   for(int i = 0; i < N; i++){
      diag[i] = 4.0;
      col = i - nx;
      if (col >= 0){
         data[k] = -1.0;
         idx[k] = col;
         k++;
      }
      if (i > block_start){
         col = i - 1;
         data[k] = -1.0;
         idx[k] = col;
         k++;
      }
      data[k] = 4.0;
      idx[k] = i;
      k++;
      if (i < block_end){
         col = i + 1;
         data[k] = -1.0;
         idx[k] = col;
         k++;
      }
      col = i + nx;
      if (col < N){
         data[k] = -1.0;
         idx[k] = col;
         k++;
      }
      start_idx[i+1] = k;

      if (i == block_end){
         block_end += nx;
         block_start += nx;
      }
   }

   for (int k = 0; k < nnz; k++){
      if (store_type == SparseMatrixStorageType::CSC){
         row_idx[k] = idx[k];
      }
      else {
         col_idx[k] = idx[k];
      }
   }

   if (store_type == SparseMatrixStorageType::COO){
      k = 0;
      for (int i = 0; i < num_rows; i++){
         for (int jj = start_idx[i]; jj < start_idx[i+1]; jj++){
            if (store_type == SparseMatrixStorageType::CSC){
               col_idx[k] = i;
            }
            else {
               row_idx[k] = i;
            }
            k++;
         }
      }
   }
}

/* Generate a random sparse matrix */
void SparseMatrix::ConstructRandomMatrix(void)
{
   diag.resize(num_rows);
   double low = -1.0/(double)(max_row_nnz);
   double high = 1.0/(double)(max_row_nnz);

   //srand(time(NULL));
   srand(0);
   if (mat_type == MatrixType::lower || 
       mat_type == MatrixType::upper ||
       mat_type == MatrixType::non_symm){
      vector<vector<int>> rows(num_rows, vector<int>()); 
      vector<vector<double>> nzval(num_rows, vector<double>());
      nnz = 0;
      for (int i = 0; i < num_rows; i++){
         int i_max_row_nnz;
         if (mat_type == MatrixType::lower){
            i_max_row_nnz = std::min(i+1, max_row_nnz);
         }
         else if (mat_type == MatrixType::upper){
            i_max_row_nnz = std::min(num_rows-i, max_row_nnz);
         }
         else {
            i_max_row_nnz = max_row_nnz;
         }
         int row_nnz = (int)RandDouble(1, i_max_row_nnz);
         int count = 1;
         diag[i] = RandDouble(0.0, 1.0);
         //rows[i].push_back(i);
         //nzval[i].push_back(RandDouble(-1.0, 1.0));
         //nnz++;

         //while(count < row_nnz){
         //   int j;
         //   if (mat_type == MatrixType::lower){
         //      j = (int)RandInt(0, i, time(NULL));
         //   }
         //   else if (mat_type == MatrixType::upper){
         //      j = (int)RandInt(i, num_rows-1, time(NULL));
         //   }
         //   else {
         //      j = (int)RandInt(0, num_rows-1, time(NULL));
         //   }
         //   vector<int>::iterator it;

         //   it = find(rows[i].begin(), rows[i].end(), j);
         //   if (it == rows[i].end() && i != j){
         //      rows[i].push_back(j);
         //      nzval[i].push_back(RandDouble(0.0, high));
         //      count++;
         //      nnz++;
         //   }
         //}

         for (int k = 0; k < row_nnz; k++){
            int j;
            if (mat_type == MatrixType::lower){
               j = (int)RandInt(0, i, time(NULL));
            }
            else if (mat_type == MatrixType::upper){
               j = (int)RandInt(i, num_rows-1, time(NULL));
            }
            else {
               j = (int)RandInt(0, num_rows-1, time(NULL));
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

      num_cols = num_rows;
      if (store_type == SparseMatrixStorageType::CSC ||
          store_type == SparseMatrixStorageType::COO){
         row_idx.resize(nnz);
      }
      if (store_type == SparseMatrixStorageType::CSR ||
          store_type == SparseMatrixStorageType::COO){
         col_idx.resize(nnz);
      }
      data.resize(nnz);
      start_idx.resize(num_rows+1);

      vector<int> counts(num_rows, 0);

      for (int i = 0; i < num_rows; i++){
         for (int k = 0; k < rows[i].size(); k++){
            int row = i;
            int col = rows[i][k];
            int idx;
            if (store_type == SparseMatrixStorageType::CSC){
               idx = col;
            }
            else {
               idx = row;
            }
            counts[idx]++;
         }
      }

      for (int idx = 0; idx < num_rows; idx++){
         start_idx[idx+1] = start_idx[idx] + counts[idx];
         counts[idx] = 0;
      }

      for (int i = 0; i < num_rows; i++){
         for (int k = 0; k < rows[i].size(); k++){
            int row = i;
            int col = rows[i][k];
            double elem = nzval[i][k];

            int idx;
            if (store_type == SparseMatrixStorageType::CSC){
               idx = col;
            }
            else {
               idx = row;
            }
            int kk = start_idx[idx] + counts[idx];
            counts[idx]++;

            if (store_type == SparseMatrixStorageType::COO){
               row_idx[kk] = row;
               col_idx[kk] = col;
            }
            else {
               if (store_type == SparseMatrixStorageType::CSC){
                  row_idx[kk] = row;
               }
               else {
                  col_idx[kk] = col;
               }
            }
            data[kk] = elem;
         }
      }
   }
}

/* print coordinate format of sparse matrix */
void SparseMatrix::PrintMatrix(char *filename)
{
   FILE *file = fopen(filename, "w");
   
   //if (print_diag_flag == 1){
   //   for (int i = 0; i < n; i++){
   //      fprintf(file, "%d %d %.15e\n", i+1, i+1, diag[i]);
   //   }
   //}
   if (store_type == SparseMatrixStorageType::CSC){
      for (int j = 0; j < num_rows; j++){
         for (int kk = start_idx[j]; kk < start_idx[j+1]; kk++){
            int row = row_idx[kk]+1;
            int col = j+1;
            fprintf(file, "%d %d %.15e\n", row, col, data[kk]);
         }
      }
   }
   else {
      for (int i = 0; i < num_rows; i++){
         for (int kk = start_idx[i]; kk < start_idx[i+1]; kk++){
            int row = i+1;
            int col = col_idx[kk]+1;
            fprintf(file, "%d %d %.15e\n", row, col, data[kk]);
         }
      }
   }
   
   fclose(file);
}

/* read matrix from binary file. matrix entries must be ordered by increasing row index then increasing column index */
void SparseMatrix::FreadBinaryMatrix(void)
{
   size_t size;
   int temp_size;
   Triplet_AOS *buffer;

   FILE *mat_file = fopen(input_file_name, "rb");

   fseek(mat_file, 0, SEEK_END);
   size = ftell(mat_file);
   rewind(mat_file);
   buffer = (Triplet_AOS *)malloc(sizeof(Triplet_AOS) * size);
   fread(buffer, sizeof(Triplet_AOS), size, mat_file);

   int file_lines = size/sizeof(Triplet_AOS);
   num_rows = (int)buffer[0].i;
   num_rows = num_rows;
   nnz = (int)buffer[0].j;

   if (!store_diag_in_vec){
      nnz -= num_rows;
   }

   vector<int> counts(num_rows, 0);
   vector<int> include_elem_flag(file_lines, 1);

   for (int k = 1; k < file_lines; k++){
      int row = buffer[k].i-1;
      int col = buffer[k].j-1;
      
      include_elem_flag[k] = 1;
      if (row == col){
         if (store_diag_in_vec){
            include_elem_flag[k] = 0;
         }
      }
      else if (row < col){
         if (mat_type == MatrixType::lower){
            include_elem_flag[k] = 0;
            nnz--;
         }
      }
      else if (row > col){
         if (mat_type == MatrixType::upper){
            include_elem_flag[k] = 0;
            nnz--;
         }
      }

      if (include_elem_flag[k] == 1){
         int idx;
         if (store_type == SparseMatrixStorageType::CSC){
            idx = col;
         }
         else {
            idx = row;
         }
         counts[idx]++;
      }
   }

   if (store_type == SparseMatrixStorageType::CSC ||
       store_type == SparseMatrixStorageType::COO){
      row_idx.resize(nnz);
   }
   if (store_type == SparseMatrixStorageType::CSR ||
       store_type == SparseMatrixStorageType::COO){
      col_idx.resize(nnz);
   }
   data.resize(nnz);
   start_idx.resize(num_rows+1);
   diag.resize(num_rows);
   for (int idx = 0; idx < num_rows; idx++){
      start_idx[idx+1] = start_idx[idx] + counts[idx];
      counts[idx] = 0;
   }

   for (int k = 1; k < file_lines; k++){
      int row = buffer[k].i-1;
      int col = buffer[k].j-1;
      double elem = buffer[k].val;

      if (elem == 0.0) elem = 1.0;

      if (row == col){
         diag[row] = elem;
      }
     
      if (include_elem_flag[k] == 1) {   
         int kk, idx;
         if (store_type == SparseMatrixStorageType::CSC){
            idx = col;
         }
         else {
            idx = row;
         }
         kk = start_idx[idx] + counts[idx];
         counts[idx]++;
         if (store_type == SparseMatrixStorageType::COO){
            row_idx[kk] = row;
            col_idx[kk] = col; 
         }
         else {
            if (store_type == SparseMatrixStorageType::CSC){
               row_idx[kk] = row;
            }
            else {
               col_idx[kk] = col;
            }
         }
         data[kk] = elem;
      }
   }

   fclose(mat_file);
}
