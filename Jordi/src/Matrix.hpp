#ifndef MATRIX_H
#define MATRIX_H

/* matrix pattern types */
#define MATRIX_NONSYMMETRIC 0
#define MATRIX_LOWER 1
#define MATRIX_UPPER 2

/* test matrix types */
#define PROBLEM_RANDOM 0
#define PROBLEM_5PT_POISSON 1
#define PROBLEM_FILE 2

#include "Main.hpp"
#include "Misc.hpp"

///* CSR struct */
//typedef struct{
//   int *start; /* pointer to row starts */
//   int *i; /* columns indices (used in COO and CSC) */
//   int *j; /* row indices (used in COO and CSR) */
//   double *data; /* matrix values */
//   double *diag; /* diagonal elements */
//   int n; /* number of rows */
//   int m; /* number of columns */
//   int nnz; /* number of non-zero values */
//}Matrix;
//
//void Laplace_2D_5pt(InputData input, Matrix *A, int n);
//
//void RandomMatrix(InputData input,
//                  Matrix *A,
//                  int n,
//                  int max_row_nnz,
//                  int mat_type,
//                  int csc_flag,
//                  int coo_flag);
//
//void PrintMatrix(Matrix A,
//                 char *filename,
//                 int print_diag_flag,
//                 int csc_flag);
//
//void freadBinaryMatrix(char *mat_file_str,
//                       Matrix *A,
//                       int include_diag_flag,
//                       int csc_flag,
//                       int coo_flag,
//                       int mat_type);

enum class MatrixType {
   lower,
   upper,
   symm,
   non_symm,
   spd
};

enum class SparseMatrixStorageType {
   CSR,
   CSC,
   COO
};

typedef struct {
   MatrixType mat_type;
   FileType file_type = FileType::bin;
   int num_rows; /* number of rows */
   int num_cols; /* number of columns */
   vector<int> grid_len; /* for structured grids, length of each grid dimension */
   int max_row_nnz;
   char file_name[256];
}MatrixInput;

typedef struct SparseMatrixInputStruct : MatrixInput {
   SparseMatrixStorageType store_type;
   bool store_diag_in_vec;
}SparseMatrixInput;

class Matrix
{
public:
   //Matrix(MatrixInput input) {};
   //~Matrix(void) {};
   virtual void ConstructMatrixFromFile(void) = 0;
   virtual void ConstructRandomMatrix(void) = 0;
   virtual void PrintMatrix(char *out_filename) = 0;
   int GetNumRows(void);
   int GetNumCols(void);

protected:
   int num_rows, num_cols;
   int max_row_nnz;
   MatrixType mat_type; 
   char* input_file_name;
   FileType input_file_type;
};

class SparseMatrix : public Matrix
{
public:
   SparseMatrix(SparseMatrixInput input);
   ~SparseMatrix(void);
   void ConstructLaplace2D5pt(void);
   void ConstructMatrixFromFile(void) override;
   void ConstructRandomMatrix(void) override;
   void PrintMatrix(char *out_filename) override;

   vector<double> GetValues(void);
   vector<double> GetDiagValues(void);
   vector<int> GetColIndices(void);
   vector<int> GetRowIndices(void);
   vector<int> GetIndexStarts(void);
   int GetNNZ(void);
   SparseMatrixStorageType GetStorageType(void);

protected:
   vector<int> start_idx, row_idx, col_idx;
   vector<double> data, diag;
   vector<int> grid_len;
   int nnz;
   SparseMatrixStorageType store_type;
   bool store_diag_in_vec;
   
private:
   void FreadBinaryMatrix(void);
};

#endif
