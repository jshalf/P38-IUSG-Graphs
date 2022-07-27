#include "LevelSets.hpp"

/* Compute level sets for level scheduling algorithms */
void LevelSets(SparseMatrix A, /* matrix data (input) */
               LevelSetData *lvl_set, /* level set data (output) */
               int L_flag) /* is the matrix lower or upper triangular ? */
{
   int num_rows = A.GetNumRows();
   int nnz = A.GetNNZ();
   vector<int> start = A.GetIndexStarts();
   vector<int> row_idx = A.GetRowIndices();
   vector<int> col_idx = A.GetColIndices();
   vector<double> mat_values = A.GetValues();
   SparseMatrixStorageType spmat_store_type = A.GetStorageType();

   vector<int> depth(num_rows, -1), level_size(num_rows, 0);
   lvl_set->order.resize(num_rows);
   int lump = 1;

   ///* loop over all nodes to compute depths which are levels */
   //if (input.setup_type == LEVEL_SETS_ASYNC_SETUP){

   //   vector<int> row_done_flags, deps_counts;
   //   if (spmat_store_type == SparseMatrixStorageType::CSC){
   //      deps_counts.resize(num_rows);
   //   }
   //   else {
   //      row_done_flags.resize(num_rows);
   //   }

   //   #pragma omp parallel
   //   {
   //      int i_loc, j_loc, jj_loc, n_loc, nnz_loc;
   //      int *nz_done_flags_loc;

   //      n_loc = 0;
   //      nnz_loc = 0;
   //      #pragma omp for schedule(static, lump) nowait
   //      for (int i = 0; i < num_rows; i++){
   //         n_loc++;
   //         for (int jj = start[i]; jj < start[i+1]; jj++){
   //            nnz_loc++;
   //         }
   //      }

   //      if (spmat_store_type == SparseMatrixStorageType::CSC){
   //         #pragma omp for schedule(static, lump) nowait
   //         for (int k = 0; k < nnz; k++){
   //            #pragma omp atomic
   //            deps_counts[row_idx[k]]++;
   //         }
   //      }
   //      else {
   //         nz_done_flags_loc = (int *)calloc(nnz_loc, sizeof(int));
   //      }

   //      if (spmat_store_type == SparseMatrixStorageType::CSC){
   //         #pragma omp for schedule(static, lump) nowait
   //         for (int J = 0; J < num_rows; J++){
   //            int j;
   //            if (L_flag == 1){
   //               j = J;
   //            }
   //            else {
   //               j = num_rows-J-1;
   //            }

   //            int deps_counts_j;
   //            /* stay idle until element j of x is ready to be used */
   //            do {
   //               #pragma omp atomic read
   //               deps_counts_j = deps_counts[j];
   //            } while (deps_counts_j > 0);
   //            depth[j]++;
   //            for (int kk = start[j]; kk < start[j+1]; kk++){
   //               int i = row_idx[kk];
   //               depth[i] = max(depth[i], depth[j]);
   //               #pragma omp atomic
   //               deps_counts[i]--;
   //            }
   //         }
   //      }
   //      else {
   //         jj_loc = 0;
   //         #pragma omp for schedule(static, lump) nowait
   //         for (int I = 0; I < num_rows; I++){
   //            int i;
   //            if (L_flag == 1){
   //               i = I;
   //            }
   //            else {
   //               i = num_rows-I-1;
   //            }

   //            int max_depth = -1;
   //            int jj_start = start[i];
   //            int jj_end = start[i+1];
   //            int jj_diff = jj_end - jj_start;
   //            int row_count_i = jj_diff;
   //            while (row_count_i > 0){
   //               int jj_loc_temp = jj_loc;
   //               for (int jj = jj_start; jj < jj_end; jj++){
   //                  if (nz_done_flags_loc[jj_loc_temp] == 0){
   //                     int j = col_idx[jj];
   //                     int row_done_flag_j;
   //                     #pragma omp atomic read
   //                     row_done_flag_j = row_done_flags[j];
   //                     if (row_done_flag_j == 1){
   //                        max_depth = max(max_depth, depth[j]);
   //                        row_count_i--;
   //                        nz_done_flags_loc[jj_loc_temp] = 1;
   //                     }
   //                  }
   //                  jj_loc_temp++;
   //               }
   //            }
   //            depth[i] = 1 + max_depth;
   //            #pragma omp atomic write
   //            row_done_flags[i] = 1;
   //            jj_loc += jj_diff;
   //         }
   //      }

   //      if (spmat_store_type == SparseMatrixStorageType::CSC){

   //      }
   //      else {
   //         free(nz_done_flags_loc);
   //      }
   //   }
   //}
   //else {
      for (int I = 0; I < num_rows; I++){
         int i;
         if (L_flag == 1){
            i = I;
         }
         else {
            i = num_rows-I-1;
         }

         if (spmat_store_type == SparseMatrixStorageType::CSC){
            depth[i]++;
            for (int kk = start[i]; kk < start[i+1]; kk++){
               int j = row_idx[kk];
               depth[j] = max(depth[j], depth[i]);
            }
         }
         else {
            int max_depth = -1;
            for (int kk = start[i]; kk < start[i+1]; kk++){
               int j = col_idx[kk];
               max_depth = max(max_depth, depth[j]);
            }
            depth[i] = 1 + max_depth;
         }
      }
   //}


   lvl_set->num_levels = 0;
   for (int i = 0; i < num_rows; i++){
      level_size[depth[i]]++;
      if (lvl_set->num_levels < depth[i]){
         lvl_set->num_levels = depth[i];
      }
   }
   lvl_set->num_levels++;

   lvl_set->level_size.resize(lvl_set->num_levels);
   lvl_set->level_start.resize(lvl_set->num_levels+1);
   lvl_set->level_start[0] = 0;
   /* compute the starting points for the ordering for each level set */
   for (int i = 0; i < lvl_set->num_levels; i++){
      lvl_set->level_size[i] = level_size[i];
      lvl_set->level_start[i+1] = lvl_set->level_start[i] + lvl_set->level_size[i];
   }

   vector<int> level_counts(lvl_set->num_levels, 0);
   for (int i = 0; i < num_rows; i++){
      int d = depth[i];
      lvl_set->order[lvl_set->level_start[d]+level_counts[d]] = i;
      level_counts[d]++;
   }
}

void LevelSetsDestroy(LevelSetData *lvl_set)
{
}
