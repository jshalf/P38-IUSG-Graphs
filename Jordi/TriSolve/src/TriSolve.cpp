#include "TriSolve.hpp"
#include "Matrix.hpp"
#include "MsgQ.hpp"
#include "Misc.hpp"
#include "Solver.hpp"

//#ifdef USE_DEVA
//#include "deva_includes.hpp"
//#endif


int qstride = 1;

/*****************************************************
 * Sparse triangular solver member functions
 *****************************************************/
void SparseTriSolver::Setup(OutputData &output_data_input)
{
   double start;
   output_data_ptr = &output_data_input;

   int num_rows = A->GetNumRows();
   start = omp_get_wtime();

   idx_solve_order.resize(num_rows);
   for (int i = 0; i < num_rows; i++){
      idx_solve_order[i] = i;
   }

   output_data_ptr->setup_wtime += omp_get_wtime() - start;
}

void SparseTriSolver::InitSolveData(void)
{

}

void SparseTriSolver::Solve(std::vector<double> b_input,
                            std::vector<double> &x_input,
                            OutputData &output_data_input)
{
   double wtime_start;

   output_data_ptr = &output_data_input;
   b_ptr = &b_input;
   x_ptr = &x_input;

   int num_threads = num_procs;

   wtime_start = omp_get_wtime();
   //if (num_threads == 1){
   //   SequentialSolve();
   //}
   //else {
#ifdef USE_STDTHREAD
      thread_info.th.resize(num_threads);
      vector<TriSolveParArg> arg_vec(num_threads);
      //thread_info.bar = new std::barrier(num_threads);
      pthread_barrier_init(&(pthread_info.barrier), NULL, num_threads);
      for (int t = 0; t < num_threads; t++){
         arg_vec[t].proc_id = t;
         arg_vec[t].b = b_ptr;
         arg_vec[t].x = x_ptr;
         //pthread_create(&(pthread_info.threads[t]),
         //               NULL,
         //               (THREADFUNCPTR) &SparseTriSolver::ParallelSolveVoidStar,
         //               &(arg_vec[t]));
         thread_info.th[t] = new std::thread(&SparseTriSolver::ParallelSolveFunc,
                                             this,
                                             &(arg_vec[t]));
      }
      for (int t = 0; t < num_threads; t++) {
         thread_info.th[t]->join();
         //pthread_join(pthread_info.threads[t], NULL);
      }
      pthread_barrier_destroy(&(pthread_info.barrier));
#else
      TriSolveParArg arg;
      arg.b = b_ptr;
      arg.x = x_ptr;
      ParallelSolveFunc(&arg);
#endif
   //}
}

void SparseTriSolver::SequentialSolve(void)
{
   double wtime_start = omp_get_wtime();

   int num_rows = A->GetNumRows();
   int nnz = A->GetNNZ();
   vector<int> start = A->GetIndexStarts();
   vector<int> row_idx = A->GetRowIndices();
   vector<int> col_idx = A->GetColIndices();
   vector<double> mat_values = A->GetValues();
   vector<double> diag = A->GetDiagValues();
   SparseMatrixStorageType spmat_store_type = A->GetStorageType();

   std::vector<double> &x = *x_ptr;
   std::vector<double> b =  *b_ptr;

   if (spmat_store_type == SparseMatrixStorageType::CSC){
      for (int i = 0; i < num_rows; i++) x[i] = b[i];
      for (int ii = 0; ii < num_rows; ii++){
         int i = idx_solve_order[ii];
         x[i] /= diag[i];
         for (int kk = start[i]; kk < start[i+1]; kk++){
            x[row_idx[kk]] -= mat_values[kk] * x[i];
         }
      }
   }
   else {
      for (int ii = 0; ii < num_rows; ii++){
         int i = idx_solve_order[ii];
         x[i] = b[i];
         for (int kk = start[i]; kk < start[i+1]; kk++){
            x[i] -= mat_values[kk] * x[col_idx[kk]];
         }
         x[i] /= diag[i];
      }
   }

   output_data_ptr->solve_wtime_vec[0] = omp_get_wtime() - wtime_start;
}













/*****************************************************
 * Level scheduled triangular solver member functions
 *****************************************************/
void LevelSchedTriSolver::Setup(OutputData &output_data_input)
{
   SparseTriSolver::Setup(output_data_input);
   output_data_ptr = &output_data_input;

   double start;
   ParallelInput para_input;
   int num_rows = A->GetNumRows();

   start = omp_get_wtime();

   LevelSets(*A, &(level_set_info), 1);

   level_para_info.resize(level_set_info.num_levels);
   for (int l = 0; l < level_set_info.num_levels; l++){
      /* communication parameters */
      para_input.comm_input.comm_type = CommunicationType::atomic;
      /* partition parameters */
      para_input.part_input.idx_start = level_set_info.level_start[l];
      para_input.part_input.size_glob = level_set_info.level_start[l+1] - level_set_info.level_start[l];
      para_input.part_input.idx_glob = level_set_info.order;
      level_para_info[l] = new ParallelInfo(num_procs, para_input);
      level_para_info[l]->Part()->ConstructPartition();
   }

   for (int i = 0; i < num_rows; i++){
      idx_solve_order[i] = level_set_info.order[i];
   }

   output_data_ptr->setup_wtime += omp_get_wtime() - start;
}

void LevelSchedTriSolver::InitSolveData(void)
{
   SparseTriSolver::InitSolveData();
}


void LevelSchedTriSolver::ParallelSolveFunc(TriSolveParArg *arg)
{
   std::vector<double>& b = *(arg->b);
   std::vector<double>& x = *(arg->x);

   int num_rows = A->GetNumRows();
   vector<int> start = A->GetIndexStarts();
   vector<int> row_idx = A->GetRowIndices();
   vector<int> col_idx = A->GetColIndices();
   vector<double> mat_values = A->GetValues();
   vector<double> diag = A->GetDiagValues();
   SparseMatrixStorageType spmat_store_type = A->GetStorageType();

   int num_threads = num_procs;

   if (spmat_store_type == SparseMatrixStorageType::CSC){
      x = b;
   }

#ifdef USE_STDTHREAD
   int tid = arg->proc_id;
#else
   int lump = 1;
   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
#endif
      double solve_start, solve_stop;
      int num_relax = 0, num_iters = 0;

#ifdef USE_STDTHREAD
      solve_start = omp_get_wtime();
#else
      solve_start = omp_get_wtime();
#endif

      if (spmat_store_type == SparseMatrixStorageType::CSC){
         for (int l = 0; l < level_set_info.num_levels; l++){ /* loop over level sets */
            int n_loc = level_para_info[l]->Part()->GetPartitionSize()[tid];
            vector<unsigned int> my_rows = level_para_info[l]->Part()->GetPartition()[tid];
            for (int i_loc = 0; i_loc < n_loc; i_loc++) {
               int i = my_rows[i_loc];
               x[i] /= diag[i];
               for (int ii = start[i]; ii < start[i+1]; ii++){
#ifdef USE_STDTHREAD
#else
                  #pragma omp atomic
#endif
                  x[row_idx[ii]] -= mat_values[ii] * x[i];
                  num_relax++;
               }
            }
            num_iters++;

#ifdef USE_STDTHREAD
            //printf("%d %d\n", l, tid);
            pthread_barrier_wait(&(pthread_info.barrier));
#else
            #pragma omp barrier
#endif
         }
      }
      else {
         for (int l = 0; l < level_set_info.num_levels; l++){ /* loop over level sets */
            int n_loc = level_para_info[l]->Part()->GetPartitionSize()[tid];
            vector<unsigned int> my_rows = level_para_info[l]->Part()->GetPartition()[tid];
            for (int i_loc = 0; i_loc < n_loc; i_loc++) {
               int i = my_rows[i_loc];
               x[i] = b[i];
               for (int jj = start[i]; jj < start[i+1]; jj++){
                  x[i] -= mat_values[jj] * x[col_idx[jj]];
                  num_relax++;
               }
               x[i] /= diag[i];
            }
            num_iters++;
#ifdef USE_STDTHREAD
            pthread_barrier_wait(&(pthread_info.barrier));
#else
            #pragma omp barrier
#endif
            num_iters++;
         }
      }

#ifdef USE_STDTHREAD
      solve_stop = omp_get_wtime();
#else
      solve_stop = omp_get_wtime();
#endif

      output_data_ptr->solve_wtime_vec[tid] = solve_stop - solve_start;
      output_data_ptr->num_relax[tid] = num_relax;
      output_data_ptr->num_iters[tid] = num_iters;

#ifdef USE_STDTHREAD
#else
   }
#endif

   
}


















/*****************************************************
 * Asynchronous triangular solver member functions
 *****************************************************/
void AsyncTriSolver::Setup(OutputData &output_data_input,
                           TriSolverSolveOrder solve_order_input,
                           CommunicationType comm_type_input)
{
   SparseTriSolver::Setup(output_data_input);
   output_data_ptr = &output_data_input;
   solve_order = solve_order_input;
   
   int num_rows = A->GetNumRows();
   int nnz = A->GetNNZ();
   vector<int> start = A->GetIndexStarts();
   vector<int> row_idx = A->GetRowIndices();
   vector<int> col_idx = A->GetColIndices();
   vector<double> mat_values = A->GetValues();
   SparseMatrixStorageType spmat_store_type = A->GetStorageType();

   double wtime_start;
   ParallelInput para_input;

   wtime_start = omp_get_wtime();

   if (solve_order == TriSolverSolveOrder::level_sched){
      LevelSets(*A, &(level_set_info), 1);
      para_input.part_input.idx_glob = level_set_info.order;
      idx_solve_order = level_set_info.order;
   }
   else {
      para_input.part_input.idx_glob = idx_solve_order;
   }

   /* communication parameters */
   para_input.comm_input.comm_type = comm_type_input;
   /* partition parameters */
   para_input.part_input.idx_start = 0;
   para_input.part_input.size_glob = num_rows;
   para_info = new ParallelInfo(num_procs, para_input);
   para_info->Part()->ConstructPartition();
   para_info->Part()->ConstructIndexToProcMap();
   para_info->Part()->ConstructGlobToLocIndexMap();

   if (para_info->Comm()->GetCommType() == CommunicationType::MsgQ){
#if USE_DEVA
      Q = new MessageQueue<TriSolveCSRMessage>(num_procs);
#else
      Q.resize(num_procs);
      for (int p = 0; p < num_procs; p++){
         Q[p] = new MessageQueue<TriSolveCSRMessage>(qstride * num_procs);
      }
#endif
      qPut_dst_rows.resize(num_rows, vector<unsigned int>());
      for (int i = 0; i < num_rows; i++){
         for (int kk = start[i]; kk < start[i+1]; kk++){
            qPut_dst_rows[col_idx[kk]].push_back(i);
         }
      }
   }
   else {
#ifdef USE_STDTHREAD
      row_solved = (std::atomic<unsigned short> *)calloc
                           (sizeof(std::atomic<unsigned short>), num_rows);
#else
      row_solved.resize(num_rows, 0);
#endif
   }

   output_data_ptr->setup_wtime += omp_get_wtime() - wtime_start;
}

void AsyncTriSolver::InitSolveData(void)
{
   SparseTriSolver::InitSolveData();

   int num_rows = A->GetNumRows();
   int nnz = A->GetNNZ();

   if (para_info->Comm()->GetCommType() == CommunicationType::MsgQ){
   }
   else {
#ifdef USE_STDTHREAD
      for (int i = 0; i < num_rows; i++){
         row_solved[i] = 0;
      }
#else
      std::fill(row_solved.begin(), row_solved.end(), 0);
#endif
   }
}

void AsyncTriSolver::ParallelSolveFunc(TriSolveParArg *arg)
{
   if (para_info->Comm()->GetCommType() == CommunicationType::MsgQ){
      return ParallelSolveFunc_MsgQ(arg);
   }
   else {
      return ParallelSolveFunc_Atomic(arg);
   }
}


//void AsyncTriSolver::FindUpdateRow(unsigned int src_row,
//                                   unsigned int dst_row,
//                                   std::vector<unsigned int>& deps_counts,
//                                   std::queue<unsigned int>& ready_to_solve,
//                                   std::vector<double>& x)
//{
//
//}


void AsyncTriSolver::ParallelSolveFunc_MsgQ(TriSolveParArg *arg)
{
   int num_rows = A->GetNumRows();
   vector<int> start = A->GetIndexStarts();
   vector<int> row_idx = A->GetRowIndices();
   vector<int> col_idx = A->GetColIndices();
   vector<double> mat_values = A->GetValues();
   vector<double> diag = A->GetDiagValues();
   SparseMatrixStorageType spmat_store_type = A->GetStorageType();
   vector<unsigned int> row_to_proc_map = para_info->Part()->GetIndexToProcMap();

   std::vector<double>& b = *(arg->b);
   std::vector<double>& x = *(arg->x);


   int num_threads = num_procs;

   x = b;


#ifdef USE_STDTHREAD
   int tid = arg->proc_id;
#elif USE_DEVA
   deva::run( [this, &x, b] () {
   int tid = deva::rank_me();   
#else
   int lump = 1;
   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
#endif

      double solve_start, solve_stop;
      int num_relax = 0, num_iters = 0;
      int num_qPuts = 0, num_qGets = 0;
      int num_spins = 0, nothing_spins = 0;

      vector<unsigned int> my_rows = para_info->Part()->GetPartition()[tid];
      unordered_map<unsigned int, unsigned int> glob_to_loc_map = para_info->Part()->GetGlobToLocIndexMap()[tid];
      int n_loc = para_info->Part()->GetPartitionSize()[tid];
      vector<double> x_loc(n_loc);
      int num_unsolved = n_loc;
      std::vector<unsigned int> deps_counts(n_loc, 0);
      std::queue<unsigned int> ready_to_solve;
      //std::queue<TriSolveCSRMessage> recvd_msgs;
      for (int i_loc = 0; i_loc < n_loc; i_loc++){
         int i = my_rows[i_loc];
         int i_nnz = start[i+1] - start[i];
 
         if (i_nnz == 0){
            ready_to_solve.push(i_loc);
         }
         deps_counts[i_loc] = i_nnz;

         x_loc[i_loc] = b[i];
      }

      solve_start = omp_get_wtime();

      while (1){
         //int nothing_flag = 1;
         /* solve equations that are ready to be solved */
         while (ready_to_solve.size()){
            unsigned int i_loc = ready_to_solve.front();
            unsigned int i = my_rows[i_loc];
            x_loc[i_loc] /= diag[i];
            double xi = x_loc[i_loc];

            for (int j = 0; j < qPut_dst_rows[i].size(); j++){
               int dst_row = qPut_dst_rows[i][j];
               int p = row_to_proc_map[dst_row];
               if (p == tid){
                  unsigned int dst_row_loc = glob_to_loc_map[dst_row];
                  for (int jj = start[dst_row]; jj < start[dst_row+1]; jj++){
                     if (col_idx[jj] == i){
                        x_loc[dst_row_loc] -= mat_values[jj] * xi;
                        break;
                     }
                  }
                  deps_counts[dst_row_loc]--;
                  if (deps_counts[dst_row_loc] == 0){
                     ready_to_solve.push(dst_row_loc);
                  }
               }
               else {
                  TriSolveCSRMessage msg;
                  msg.src_row = i;
                  msg.dst_row = dst_row;
                  msg.data = xi;
#if USE_DEVA
                  Q->qPut(p, msg);
#else
                  Q[p]->qPut(qstride * tid, msg);
#endif
                  num_qPuts++;
               }
            } 
            ready_to_solve.pop();
            num_unsolved--;
            num_iters++;

            //nothing_flag = 0;
         }

         if (num_unsolved == 0) break;

         /* receive messages */
         while (ready_to_solve.empty()){
            TriSolveCSRMessage msg;
#if USE_DEVA
            deva::progress();
            int get_flag = Q->qGet(tid, &msg);
#else
            for (int t = 0; t < num_threads; t++){
               int get_flag = Q[tid]->qGet(qstride * t, &msg);
#endif
               if (get_flag){
                  num_qGets++;
   
                  int src_row = msg.src_row;
                  int i = msg.dst_row;
                  double xj = msg.data;
                  unsigned int i_loc = glob_to_loc_map[i];
   
                  for (int jj = start[i]; jj < start[i+1]; jj++){
                     if (col_idx[jj] == src_row){
                        x_loc[i_loc] -= mat_values[jj] * xj;
                        break;
                     }
                  }
                  deps_counts[i_loc]--;
                  if (deps_counts[i_loc] == 0){
                     ready_to_solve.push(i_loc);
                  }

                  //nothing_flag = 0;
               }
#if USE_DEVA
#else
            }
#endif
         }

//#if USE_DEVA
//         deva::progress();
//#endif

         num_spins++;
         //if (nothing_flag == 1){
         //   nothing_spins++;
         //}

         ///* receive messages */
         //while (1){
         //   TriSolveCSRMessage msg;
         //   int get_flag = Q->qGet(tid, &msg);
         //   if (get_flag == 0) break;
         //   num_qGets++;
         //   recvd_msgs.push(msg);
         //}
         //while (recvd_msgs.size()){
         //   TriSolveCSRMessage msg = recvd_msgs.front();
         //   recvd_msgs.pop();
         //   int src_row = msg.src_row;
         //   int i = msg.dst_row;
         //   double xj = msg.data;
         //   unsigned int i_loc = glob_to_loc_map[i];
         //   //printf("%d %d %d %d %e\n", tid, i, i_loc, n_loc, xj);
         //   for (int jj = start[i]; jj < start[i+1]; jj++){
         //      if (col_idx[jj] == src_row){ 
         //         x_loc[i_loc] -= mat_values[jj] * xj;
         //         break;
         //      }
         //   }
         //   deps_counts[i_loc]--;
         //   if (deps_counts[i_loc] == 0){
         //      ready_to_solve.push(i_loc);
         //   }
         //}
      }

      solve_stop = omp_get_wtime();

      //printf("%d %d %d %d\n", tid, n_loc, num_spins, nothing_spins);

      for (int i_loc = 0; i_loc < n_loc; i_loc++){
         int i = my_rows[i_loc];
         x[i] = x_loc[i_loc];
      }

      output_data_ptr->solve_wtime_vec[tid] = solve_stop - solve_start;
      output_data_ptr->num_relax[tid] = num_relax;
      output_data_ptr->num_iters[tid] = num_iters;
#ifdef USE_STDTHREAD
#elif USE_DEVA
   });
#else
   }
#endif

   
}


//TODO: fix STD atomic load and stores (not getting correct answer)
void AsyncTriSolver::ParallelSolveFunc_Atomic(TriSolveParArg *arg)
{
   std::vector<double>& b = *(arg->b);
   std::vector<double>& x = *(arg->x);

   int num_rows = A->GetNumRows();
   vector<int> start = A->GetIndexStarts();
   vector<int> row_idx = A->GetRowIndices();
   vector<int> col_idx = A->GetColIndices();
   vector<double> mat_values = A->GetValues();
   vector<double> diag = A->GetDiagValues();
   SparseMatrixStorageType spmat_store_type = A->GetStorageType();

   int num_threads = num_procs;

   x = b;


#ifdef USE_STDTHREAD
   int tid = arg->proc_id;
#else
   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
#endif
      double solve_start, solve_stop;
      int num_relax = 0, num_iters = 0;

      solve_start = omp_get_wtime();

      vector<unsigned int> my_rows = para_info->Part()->GetPartition()[tid];
      int n_loc = para_info->Part()->GetPartitionSize()[tid];

      for (int i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
         int num_spins = 0;
         int i = my_rows[i_loc];
         int jj_start = start[i];
         int jj_end = start[i+1];
         int jj_diff = jj_end - jj_start;
         int deps_counter_i = jj_diff;
         vector<bool> nz_used_flags(jj_diff, false);
         while (deps_counter_i > 0){ /* loop until x[i] has been completed */
            int jj_loc = 0;
            for (int jj = jj_start; jj < jj_end; jj++){
               if (nz_used_flags[jj_loc] == false){ /* has this non-zero been used? */
                  int j = col_idx[jj];
                  unsigned short int row_solved_flag_j;
                  /* check if x[j] is available (deps_counters[j] must be zero) */
#ifdef USE_STDTHREAD
                  row_solved_flag_j = row_solved[j].load();
#else
                  #pragma omp atomic read
                  row_solved_flag_j = row_solved[j];//row_solved_flag_j = row_solved_flags[j].c;
#endif

                  num_spins++;

                  /* if x[j] is available, update x[i] and deps_counters[i] */
                  if (row_solved_flag_j){
                     x[i] -= mat_values[jj] * x[j];

                     num_relax++;

                     deps_counter_i--;

                     nz_used_flags[jj_loc] = true;
                  }
               }
               jj_loc++;
            }
         }
         x[i] /= diag[i];
#ifdef USE_STDTHREAD
         row_solved[i].store(1);
#else
         #pragma omp atomic write
         row_solved[i] = 1;//row_solved_flags[i].c = 1;
#endif
      }

      solve_stop = omp_get_wtime();

      output_data_ptr->solve_wtime_vec[tid] = solve_stop - solve_start;
      output_data_ptr->num_relax[tid] = num_relax;
      output_data_ptr->num_iters[tid] = num_iters;
#ifdef USE_STDTHREAD
#else
   }
#endif

   
}


///******************
// * Serial TriSolve 
// ******************/
//void TriSolve_Seq(TriSolveData *ts,
//                  SparseMatrix A, /* triangular matrix */
//                  int *A_perm, /* ordering for computing elements of x */
//                  double *x, /* solution (output) */
//                  double *b /* right-hand side */
//                  )
//{
//   int num_rows = A.GetNumRows();
//   int nnz = A.GetNNZ();
//   vector<int> start = A.GetIndexStarts();
//   vector<int> row_idx = A.GetRowIndices();
//   vector<int> col_idx = A.GetColIndices();
//   vector<double> mat_values = A.GetValues();
//   vector<double> diag = A.GetDiagValues();
//   SparseMatrixStorageType spmat_store_type = A.GetStorageType();
//
//   if (spmat_store_type == SparseMatrixStorageType::CSC){
//      for (int i = 0; i < num_rows; i++) x[i] = b[i];
//      for (int I = 0; I < num_rows; I++){
//         int i = A_perm[I];
//         x[i] /= diag[i];
//         for (int kk = start[i]; kk < start[i+1]; kk++){
//            x[row_idx[kk]] -= mat_values[kk] * x[i];
//         }
//      }
//   }
//   else {
//      for (int I = 0; I < num_rows; I++){
//         int i = A_perm[I];
//         x[i] = b[i];
//         for (int kk = start[i]; kk < start[i+1]; kk++){
//            x[i] -= mat_values[kk] * x[col_idx[kk]];
//         }
//         x[i] /= diag[i];
//      }
//   }
//}
//
//
//
//
//
//
//
//
//
//
//
//
//typedef struct {
//   TriSolveData *ts;
//   LevelSetData *level_set_info;
//   SparseMatrix *A;
//   double *x;
//   double *b;
//} TriSolveArg;
//
//TriSolveArg *ts_arg;
//
//typedef struct {
//   int src_row;
//   double data;
//} TriSolveCSRMessage;
//
//void *TriSolve_LevelSchedule(void *arg);
//void *TriSolve_Async(void *arg);
//
///****************************
// * Level-scheduled TriSolve
// ****************************/
//void TriSolve(TriSolveData *ts,
//              LevelSetData level_set_info, /* level set data */
//              SparseMatrix A, /* triangular matrix */
//              double *x, /* solution (output) */
//              double *b) /* right-hand side */
//{
//   //TriSolveArg *ts_arg = (TriSolveArg *)malloc(sizeof(TriSolveArg));
//   ts_arg = (TriSolveArg *)malloc(sizeof(TriSolveArg));
//   ts_arg->ts = ts;
//   ts_arg->level_set_info = &level_set_info;
//   ts_arg->A = &A;
//   ts_arg->x = x;
//   ts_arg->b = b;
//
//   int num_threads = ts->input.num_threads;
//
//#ifdef USE_STDTHREAD
//   ts->pthread_info.threads.resize(num_threads);
//   vector<int> tid_vec(num_threads);
//   for (int t = 0; t < num_threads; t++) tid_vec[t] = t;
//   pthread_barrier_init(&(ts->pthread_info.barrier), NULL, num_threads);
//   for (int t = 0; t < num_threads; t++){
//      if (ts->input.solver_type == TRISOLVE_LEVEL_SCHEDULED){
//         pthread_create(&(ts->pthread_info.threads[t]), NULL, TriSolve_LevelSchedule, &(tid_vec[t]));
//      }
//      else {
//         pthread_create(&(ts->pthread_info.threads[t]), NULL, TriSolve_Async, &(tid_vec[t]));
//      }
//   }
//   for (int t = 0; t < num_threads; t++) pthread_join(ts->pthread_info.threads[t], NULL);
//   pthread_barrier_destroy(&(ts->pthread_info.barrier));
//#else
//   if (ts->input.solver_type == TRISOLVE_LEVEL_SCHEDULED){
//      TriSolve_LevelSchedule(ts_arg);
//   }
//   else {
//      TriSolve_Async(ts_arg);
//   }
//#endif
//
//   //free(ts_arg);
//}
//
//void *TriSolve_LevelSchedule(void *arg)
//{
//   //TriSolveArg *ts_arg = (TriSolveArg *)arg;
//   TriSolveData *ts = ts_arg->ts;
//   LevelSetData level_set_info = *(ts_arg->level_set_info);
//   SparseMatrix A = *(ts_arg->A);
//   double *x = ts_arg->x;
//   double *b = ts_arg->b;
//
//   int num_rows = A.GetNumRows();
//   vector<int> start = A.GetIndexStarts();
//   vector<int> row_idx = A.GetRowIndices();
//   vector<int> col_idx = A.GetColIndices();
//   vector<double> mat_values = A.GetValues();
//   vector<double> diag = A.GetDiagValues();
//   SparseMatrixStorageType spmat_store_type = A.GetStorageType();
//
//   int num_threads = ts->input.num_threads;
//
//#ifdef USE_STDTHREAD
//   int *tid_ptr = (int *)arg;
//   int tid = *tid_ptr;
//#elif USE_DEVA
//#else
//   int lump = 1;
//   #pragma omp parallel
//   {
//      int tid = omp_get_thread_num();
//#endif
//      double solve_start, solve_stop;
//      vector<int> *my_rows, n_loc;
//
//      n_loc.resize(level_set_info.num_levels);
//      my_rows = (vector<int> *)calloc(level_set_info.num_levels, sizeof(vector<int>));
//      for (int l = 0; l < level_set_info.num_levels; l++){
//         int m = level_set_info.level_start[l+1] - level_set_info.level_start[l];
//         n_loc[l] = partition::GetThreadPartitionSize(tid, num_threads, m);
//         my_rows[l] = partition::GetThreadPartition(tid, num_threads, n_loc[l], level_set_info.level_start[l], level_set_info.order);
//         for (int ii = 0; ii < n_loc[l]; ii++) {
//            int i = my_rows[l][ii];
//            x[i] = b[i];
//         }
//      }
//
//      int num_relax = 0, num_iters = 0;
//
//#ifdef USE_STDTHREAD
//#elif USE_DEVA
//      solve_start = omp_get_wtime();
//#else
//      solve_start = omp_get_wtime();
//#endif
//
//      if (spmat_store_type == SparseMatrixStorageType::CSC){
//         for (int l = 0; l < level_set_info.num_levels; l++){ /* loop over level sets */
//            for (int i_loc = 0; i_loc < n_loc[l]; i_loc++) {
//               int i = my_rows[l][i_loc];
//               x[i] /= diag[i];
//               for (int ii = start[i]; ii < start[i+1]; ii++){
//#ifdef USE_STDTHREAD
//#else
//                  #pragma omp atomic
//#endif
//                  x[row_idx[ii]] -= mat_values[ii] * x[i];
//                  num_relax++;
//               }
//            }
//            num_iters++;
//         }
//      }
//      else {
//         for (int l = 0; l < level_set_info.num_levels; l++){ /* loop over level sets */
//            for (int i_loc = 0; i_loc < n_loc[l]; i_loc++){
//               int i = my_rows[l][i_loc];
//               x[i] = b[i];
//               for (int jj = start[i]; jj < start[i+1]; jj++){
//                  x[i] -= mat_values[jj] * x[col_idx[jj]];
//                  num_relax++;
//               }
//               x[i] /= diag[i];
//            }
//            num_iters++;
//#ifdef USE_STDTHREAD
//            //printf("%d %d\n", l, tid);
//            pthread_barrier_wait(&(ts->pthread_info.barrier));
//#else
//            #pragma omp barrier
//#endif
//            num_iters++;
//         }
//      }
//
//#ifdef USE_STDTHREAD
//      solve_stop = omp_get_wtime();
//#else
//      solve_stop = omp_get_wtime();
//#endif
//
//      ts->output.solve_wtime_vec[tid] = solve_stop - solve_start;
//
//      ts->output.num_relax[tid] = num_relax;
//      ts->output.num_iters[tid] = num_iters;
//
//      free(my_rows);
//#ifdef USE_STDTHREAD
//#else
//   }
//#endif
//
//   
//}
//
//
///******************
// * Async TriSolve  
// ******************/
//void *TriSolve_Async(void *arg)
//{
//   TriSolveData *ts = ts_arg->ts;
//   LevelSetData level_set_info = *(ts_arg->level_set_info);
//   SparseMatrix A = *(ts_arg->A);
//   double *x = ts_arg->x;
//   double *b = ts_arg->b;
//
//   int num_threads = ts->input.num_threads;
//
//   int num_rows = A.GetNumRows();
//   int nnz = A.GetNNZ();
//   vector<int> start = A.GetIndexStarts();
//   vector<int> row_idx = A.GetRowIndices();
//   vector<int> col_idx = A.GetColIndices();
//   vector<double> mat_values = A.GetValues();
//   vector<double> diag = A.GetDiagValues();
//   SparseMatrixStorageType spmat_store_type = A.GetStorageType();
//
//   int lump = 1; 
//   vector<int> deps_counters, row_solved_flags;
//   vector<vector<int>> deps_counters_thread;
//   //cache *row_solved_flags;
//   
//   MessageQueue<TriSolveCSRMessage> *Q, *Q_comm;
//   vector<int> row_to_thread;
//   vector<int> t_sum;
//   vector<vector<int>> t_sum_thread;
//   if (ts->input.MsgQ_flag == 1){ /* set up message queue data */
//      if (spmat_store_type == SparseMatrixStorageType::CSC){
//         deps_counters.resize(num_rows);
//         deps_counters_thread.resize(num_threads, vector<int>());
//         for (int t = 0; t < num_threads; t++){
//            deps_counters_thread[t].resize(num_rows, 0);
//         }
//
//         Q = new MessageQueue<TriSolveCSRMessage>(num_rows);
//      }
//      else {
//         row_to_thread.resize(num_rows);
//         t_sum.resize(num_threads);
//         t_sum_thread.resize(num_threads, vector<int>());
//         for (int t = 0; t < ts->input.num_threads; t++){
//            t_sum_thread[t].resize(num_threads, 0);
//         }
//
//         Q_comm = new MessageQueue<TriSolveCSRMessage>(num_rows);
//         Q = new MessageQueue<TriSolveCSRMessage>(nnz);
//      }
//   }
//   else {
//      if (spmat_store_type == SparseMatrixStorageType::CSC){
//         deps_counters.resize(num_rows);
//      }
//      else {
//         if (ts->input.fine_grained_flag == 1){
//            deps_counters.resize(num_rows);
//         }
//         else {
//            row_solved_flags.resize(num_rows);
//         }
//      }
//   }
//
//   vector<int> row_order, nz_order;
//   if (ts->input.fine_grained_flag == 1){
//      nz_order.resize(nnz);
//      if (ts->input.solver_type == TRISOLVE_ASYNC_LEVEL_SCHEDULED){
//         int k = 0; 
//         for (int ii = 0; ii < num_rows; ii++){
//            int i = level_set_info.order[ii];
//            for (int jj = start[i]; jj < start[i+1]; jj++){
//               nz_order[k] = jj;
//               k++;
//            }
//         }
//      }
//      else {
//         for (int i = 0; i < nnz; i++){
//            nz_order[i] = i;
//         }
//      }
//   }
//   row_order.resize(num_rows);
//   if (ts->input.solver_type == TRISOLVE_ASYNC_LEVEL_SCHEDULED){
//      for (int i = 0; i < num_rows; i++){
//         row_order[i] = level_set_info.order[i];
//      }
//   }
//   else {
//      for (int i = 0; i < num_rows; i++){
//         row_order[i] = i;
//      }
//   }
//
//   vector<int> row_to_tid_map(num_rows);
//
//#ifdef USE_STDTHREAD
//   int *tid_ptr = (int *)arg;
//   int tid = *tid_ptr;
//#else
//   #pragma omp parallel
//   {
//      int tid = omp_get_thread_num();
//#endif
//      int i_loc, j_loc, jj_loc, n_loc, nnz_loc, i_prev, i_loc_prev;
//      int kk, k;
//      double solve_wtime, setup_wtime, wtime_start, wtime_stop;
//      double comp_wtime_start, comp_wtime_stop, comp_wtime = 0.0;
//      double MsgQ_wtime_start, MsgQ_wtime_stop, MsgQ_wtime = 0.0;
//      double MsgQ_put_wtime_start, MsgQ_put_wtime_stop, MsgQ_put_wtime = 0.0;
//      double MsgQ_get_wtime_start, MsgQ_get_wtime_stop, MsgQ_get_wtime = 0.0;
//      uint64_t MsgQ_cycles_start, MsgQ_cycles_stop, MsgQ_cycles = 0;
//      uint64_t MsgQ_put_cycles_start, MsgQ_put_cycles_stop, MsgQ_put_cycles = 0;
//      uint64_t MsgQ_get_cycles_start, MsgQ_get_cycles_stop, MsgQ_get_cycles = 0;
//      int num_qGets = 0, num_qPuts = 0;
//      int dummy_num_qGets = 0, dummy_num_qPuts = 0;
//      int dummy_num_spins = 0;
//      int temp_num_qPuts = 0, temp_num_qGets = 0;
//      vector<int> nz_used_flags_loc, row_solved_flags_loc;
//      int num_relax, num_iters;
//      int atomic_flag = ts->input.atomic_flag;
//      vector<double> x_loc;
//      double dummy;
//      int num_gets;
//      vector<vector<int>> put_targets;
//
//
//
//
///********************************************************************
// *
// *                            SETUP
// *
// ********************************************************************/
//      wtime_start = omp_get_wtime();
//      vector<int> my_rows, my_nzs;
//      n_loc = partition::GetThreadPartitionSize(tid, num_threads, num_rows);
//      my_rows = partition::GetThreadPartition(tid, num_threads, n_loc, 0, row_order);
//      if (ts->input.fine_grained_flag == 1){
//         nnz_loc = partition::GetThreadPartitionSize(tid, num_threads, nnz);
//         my_nzs = partition::GetThreadPartition(tid, num_threads, nnz_loc, 0, nz_order);
//      }
//      else {
//         nnz_loc = 0;
//         for (int i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
//            int i = my_rows[i_loc];
//            for (int jj = start[i]; jj < start[i+1]; jj++){
//               nnz_loc++;
//            }
//         }
//
//         x_loc.resize(n_loc);
//      }
//
//      /* compute row counts (number of non-zeros per row) and other data.
//       * atomics are required for computing row counts in the csc case. */
//      i_loc = 0;
//      jj_loc = 0;
//      if (spmat_store_type == SparseMatrixStorageType::CSC){
//         if (ts->input.MsgQ_flag == 1){
//            #pragma omp for schedule(static, lump)
//            for (int k = 0; k < nnz; k++){
//               deps_counters_thread[tid][row_idx[k]]++;
//            }
//            #pragma omp for schedule(static, lump) nowait
//            for (int i = 0; i < num_rows; i++){
//               for (int t = 0; t < ts->input.num_threads; t++){
//                  deps_counters[i] += deps_counters_thread[t][i];
//               }
//            }
//         }
//         else {
//            #pragma omp for schedule(static, lump) nowait
//            for (int k = 0; k < nnz; k++){
//               #pragma omp atomic
//               deps_counters[row_idx[k]]++;
//            }
//         }
//      }
//      else {
//         if (ts->input.fine_grained_flag == 1){
//            for (int i_loc = 0; i_loc < n_loc; i_loc++){
//               int i = my_rows[i_loc];
//
//               int jj_low = start[i];
//               int jj_high = start[i+1];
//
//               deps_counters[i] = jj_high - jj_low;
//            }
//         }
//         else {
//
//            if (ts->input.MsgQ_flag == 1){
//               put_targets.resize(n_loc);
//
//               for (i_loc = 0; i_loc < n_loc; i_loc++){
//                  int i = my_rows[i_loc];
//                  row_to_thread[i] = tid;
//               }
//#ifdef USE_STDTHREAD
//               pthread_barrier_wait(&(ts->pthread_info.barrier));
//#else
//               #pragma omp barrier
//#endif
//
//               for (i_loc = 0; i_loc < n_loc; i_loc++){
//                  int i = my_rows[i_loc];
//                  for (int jj = start[i]; jj < start[i+1]; jj++){
//                     int ii = col_idx[jj];
//                     t_sum_thread[tid][row_to_thread[ii]]++;
//                  }
//               }
//#ifdef USE_STDTHREAD
//               pthread_barrier_wait(&(ts->pthread_info.barrier));
//#else
//               #pragma omp barrier
//#endif
//
//               for (int t = 0; t < ts->input.num_threads; t++){
//                  t_sum[tid] += t_sum_thread[t][tid];
//               }
//#ifdef USE_STDTHREAD
//               pthread_barrier_wait(&(ts->pthread_info.barrier));
//#else
//               #pragma omp barrier
//#endif
//
//               int num_get = t_sum[tid];
//               int temp_count = 0;
//               for (int i_loc = 0; i_loc < n_loc; i_loc++){
//                  int i = my_rows[i_loc];
//                  for (int jj = start[i]; jj < start[i+1]; jj++){
//                     int ii = col_idx[jj];
//                     TriSolveCSRMessage msg{jj, 0.0};
//                     Q_comm->qPut(ii, msg);
//                     temp_count++;
//                  }
//               }
//               temp_count = 0;
//               while (temp_count < num_get) {
//                  for (int i_loc = 0; i_loc < n_loc; i_loc++){
//                     int i = my_rows[i_loc];
//                     TriSolveCSRMessage msg;
//                     if (Q_comm->qGet(i, &msg)){
//                        put_targets[i_loc].push_back(msg.src_row);
//                        temp_count++;
//                     }
//                  }
//               }
//               num_gets = 0;
//               for (int i_loc = 0; i_loc < n_loc; i_loc++){
//                  int i = my_rows[i_loc];
//                  int jj_start = start[i];
//                  int jj_end = start[i+1];
//                  num_gets += jj_end - jj_start;
//               }
//            }
//            else {
//            }
//         }
//         nz_used_flags_loc.resize(nnz_loc);
//      }
//
//      /* initialize solution */
//      if (ts->input.MsgQ_flag == 1){
//         if (ts->input.fine_grained_flag == 1){
//            for (int i_loc = 0; i_loc < n_loc; i_loc++){
//               int i = my_rows[i_loc];
//               x_loc[i_loc] = b[i] / diag[i];
//            }
//         }
//         else {
//            for (int i_loc = 0; i_loc < n_loc; i_loc++){
//               int i = my_rows[i_loc];
//               x_loc[i_loc] = b[i];
//            }
//#ifdef USE_STDTHREAD
//            pthread_barrier_wait(&(ts->pthread_info.barrier));
//#else
//            #pragma omp barrier
//#endif
//         }
//      }
//      else {
//         if (ts->input.fine_grained_flag == 1){
//            for (int i_loc = 0; i_loc < n_loc; i_loc++){
//               int i = my_rows[i_loc];
//               x[i] = b[i] / diag[i];
//            }
//         }
//         else {
//            for (int i_loc = 0; i_loc < n_loc; i_loc++){
//               int i = my_rows[i_loc];
//               x[i] = b[i];
//            }
//         }
//
//#ifdef USE_STDTHREAD
//         pthread_barrier_wait(&(ts->pthread_info.barrier));
//#else
//         #pragma omp barrier
//#endif
//      }
//
//      num_relax = 0, num_iters = 0;
//      jj_loc = 0;
//      i_loc = 0;
//
//      wtime_stop = omp_get_wtime();
//      ts->output.setup_wtime_vec[tid] = wtime_stop - wtime_start;
//     /*********************************
//      *          end SETUP
//      ********************************/
//
//
//
//
//
//
//
//
//      vector<double> row_wtime_loc(n_loc);
//      vector<int> row_num_spins_loc(n_loc);
//
//
//
//
//
//
//
//
//      wtime_start = omp_get_wtime();
///********************************************************************
// *            
// *                       MESSAGE QUEUES
// *
// ********************************************************************/
//      if (ts->input.MsgQ_flag == 1){
//
//         /**********************
//          *
//          *      MsgQ CSC
//          *
//          **********************/         
//         if (spmat_store_type == SparseMatrixStorageType::CSC){
//            #pragma omp for schedule(static, lump) nowait
//            for (int j = 0; j < num_rows; j++){ /* loop over rows */
//               int deps_counters_j = deps_counters[j];
//               double z;
//               /* stay idle until element j of x is ready to be used */
//               while (deps_counters_j > 0){
//                  //int get_flag = qGet(&Q, j, &z);
//                  int get_flag = 1;
//                  if (get_flag == 1){
//                     x_loc[j_loc] -= z;
//                     deps_counters_j--;
//                  }
//               }
//               /* for row j, update elements of x and deps_counters */
//               x_loc[j_loc] /= diag[j];
//               double xj = x_loc[j_loc];
//               for (int kk = start[j]; kk < start[j+1]; kk++){
//                  int i = col_idx[kk];
//                  z = mat_values[kk] * xj;
//                  //qPut(&Q, i, z);
//                  num_relax++;
//               }
//               j_loc++;
//            }
//         }
//
//
//
//
//
//
//
//
//
//
//         /**********************
//          *
//          *      MsgQ CSR
//          *
//          **********************/
//         else {
//            jj_loc = 0;
//            wtime_start = omp_get_wtime();
//            for (i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
//               int i = my_rows[i_loc];
//               int jj_start = start[i];
//               int jj_end = start[i+1];
//               int jj_diff = jj_end - jj_start;
//               int deps_counter_i = jj_diff;
//               while (deps_counter_i > 0){ /* loop until x[i] has been completed */
//                  int jj_loc_temp = jj_loc;
//                  for (int jj = jj_start; jj < jj_end; jj++){
//                     int this_nz_used = nz_used_flags_loc[jj_loc_temp];
//                     if (this_nz_used == 0){ /* has this non-zero been used? */
//                        TriSolveCSRMessage msg;
//                        int get_flag = Q->qGet(jj, &msg);
//                        /* if x[j] is available, update x[i] and other data */
//                        if (get_flag == 1){
//                           x_loc[i_loc] -= mat_values[jj] * msg.data;
//                           deps_counter_i--;
//                           nz_used_flags_loc[jj_loc_temp] = 1;
//                           num_qGets++;
//                        }
//                     }
//                     jj_loc_temp++;
//                  }
//               }
//               x_loc[i_loc] /= diag[i];
//               double xi = x_loc[i_loc];
//               int num_puts_i = put_targets[i_loc].size();
//               for (int j = 0; j < num_puts_i; j++){
//                  int put_target = put_targets[i_loc][j];
//                  TriSolveCSRMessage msg{0, xi};
//                  Q->qPut(put_target, msg);
//                  num_qPuts++;
//               }
//
//               jj_loc += jj_diff;
//            }
//         }
//      }
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
///*****************************************************************
// * 
// *                         ATOMIC
// *
// *****************************************************************/
//      else {
//         if (spmat_store_type == SparseMatrixStorageType::CSC){
//            /**************************
//             * 
//             * Atomic CSC fine-grained
//             *
//             **************************/ 
//            if (ts->input.fine_grained_flag == 1){
//               kk = 0;
//               /************
//                * atomic
//                ************/
//               if (atomic_flag == 1){
//                  for (int J = 0; J < num_rows; J++){ /* loop over rows */
//                     int j = row_order[J];
//                     k = my_nzs[kk];
//                     int k_start = start[j]; 
//                     int k_end = start[j+1];
//                     if (k < k_end && k >= k_start){ /* does this thread own non-zeros for this row? */
//                        int deps_counters_j;
//                        /* stay idle until element j of x is ready to be used */
//                        do {
//                           #pragma omp atomic read
//                           deps_counters_j = deps_counters[j];
//                        } while (deps_counters_j > 0);
//                        double xj = x[j];
//                        /* for row j, update elements of x and deps_counters */
//                        while (k < k_end && k >= k_start){
//                           int i = col_idx[k];
//                           double Tij = mat_values[k];
//                           double Tii = diag[i];
//
//                           #pragma omp atomic
//                           x[i] -= xj * (Tij / Tii);
//
//                           #pragma omp atomic
//                           deps_counters[i]--;
//
//                           kk++;
//                           k = my_nzs[kk];
//                          
//                           num_relax++;
//                        }
//                     }
//                  }
//               }
//               /************
//                * no atomic
//                ************/
//               else {
//                  for (int J = 0; J < num_rows; J++){ /* loop over rows */
//                     int j = row_order[J];
//                     k = my_nzs[kk];
//                     int k_start = start[j];
//                     int k_end = start[j+1];
//                     if (k < k_end && k >= k_start){ /* does this thread own non-zeros for this row? */
//                        int deps_counters_j;
//                        /* stay idle until element j of x is ready to be used */
//                        do {
//                           //#pragma omp atomic read
//                           deps_counters_j = deps_counters[j];
//                        } while (deps_counters_j > 0);
//                        double xj = x[j];
//                        /* for row j, update elements of x and deps_counters */
//                        while (k < k_end && k >= k_start){
//                           int i = col_idx[k];
//                           double Tij = mat_values[k];
//                           double Tii = diag[i];
//
//                           x[i] -= xj * (Tij / Tii);
//
//                           #pragma omp atomic
//                           deps_counters[i]--;
//
//                           kk++;
//                           k = my_nzs[kk];
//
//                           num_relax++;
//                        }
//                     }
//                  }
//               }
//            }
//
//
//
//
//
//
//
//
//
//
//            /**********************
//             * 
//             *     Atomic CSC
//             *
//             **********************/
//            else {
//               /************
//                * atomic
//                ************/
//               if (atomic_flag == 1){
//                  #pragma omp for schedule(static, lump) nowait
//                  for (int j = 0; j < num_rows; j++){ /* loop over rows */
//                     int deps_counters_j;
//                     /* stay idle until element j of x is ready to be used */
//                     do {
//                        #pragma omp atomic read
//                        deps_counters_j = deps_counters[j];
//                     } while (deps_counters_j > 0);
//                     /* for row j, update elements of x and deps_counters */
//                     x[j] /= diag[j];
//                     double xj = x[j];
//                     for (int kk = start[j]; kk < start[j+1]; kk++){
//                        int i = row_idx[kk];
//
//                        #pragma omp atomic
//                        x[i] -= mat_values[kk] * xj;
//
//                        #pragma omp atomic
//                        deps_counters[i]--;
//
//                        num_relax++;
//                     }
//                  }
//               }
//               /************ 
//                * no atomic 
//                ************/
//               else {
//                  #pragma omp for schedule(static, lump) nowait
//                  for (int j = 0; j < num_rows; j++){ /* loop over rows */
//                     int deps_counters_j;
//                     /* stay idle until element j of x is ready to be used */
//                     do {
//                        //#pragma omp atomic read
//                        deps_counters_j = deps_counters[j];
//                     } while (deps_counters_j > 0);
//                     /* for row j, update elements of x and deps_counters */
//                     x[j] /= diag[j];
//                     double xj = x[j];
//                     for (int kk = start[j]; kk < start[j+1]; kk++){
//                        int i = row_idx[kk];
//
//                        x[i] -= mat_values[kk] * xj;
//
//                        #pragma omp atomic
//                        deps_counters[i]--;
//
//                        num_relax++;
//                     }
//                  }
//               }
//            }
//
//
//         }
//         else {
//
//
//
//
//
//
//
//
//
//
//            /**************************
//             * 
//             * Atomic CSR fine-grained
//             *
//             **************************/
//            if (ts->input.fine_grained_flag == 1){
//               kk = 0;
//               /************
//                * atomic
//                ************/
//               if (atomic_flag == 1){
//                  for (int I = 0; I < num_rows; I++){ /* loop over rows */
//                     int i = row_order[I];
//                     k = my_nzs[kk];
//                     int kk_start = kk;
//                     int done_flag;
//                     int k_start = start[i];
//                     int k_end = start[i+1];
//                     if (k < k_end && k >= k_start){ /* does this thread own non-zeros for this row? */
//                        double Tii = diag[i];
//                        do{
//                           kk = kk_start;
//                           k = my_nzs[kk];
//                           done_flag = 1;
//                           while (k < k_end && k >= k_start){
//                              if (nz_used_flags_loc[kk] == 0){ /* has this non-zero been used? */
//                                 int j = col_idx[k];
//                                 int deps_counters_j;
//                                 /* check if x[j] is available (deps_counters[j] must be zero) */
//                                 #pragma omp atomic read
//                                 deps_counters_j = deps_counters[j];
//                                 
//                                 /* if x[j] is available, update x[i] and deps_counters[i] */
//                                 if (deps_counters_j == 0){
//                                    double Tij = mat_values[k];
//                                    double xj = x[j];
//                                  
//                                    #pragma omp atomic
//                                    x[i] -= xj * Tij / Tii;
//
//                                    num_relax++;
//
//                                    #pragma omp atomic
//                                    deps_counters[i]--;
//
//                                    nz_used_flags_loc[kk] = 1;
//                                 }
//                                 else {
//                                    done_flag = 0;
//                                 }
//                              }
//                              kk++;
//                              k = my_nzs[kk];
//                           }
//                        } while (done_flag == 0); /* loop until x[i] has been completed */
//                     }
//                  }
//               }
//               /************
//                * no atomic
//                ************/
//               else {
//                  for (int I = 0; I < num_rows; I++){ /* loop over rows */
//                     int i = row_order[I];
//                     k = my_nzs[kk];
//                     int kk_start = kk;
//                     int done_flag;
//                     int k_start = start[i];
//                     int k_end = start[i+1];
//                     if (k < k_end){ /* does this thread own non-zeros for this row? */
//                        double Tii = diag[i];
//                        do{
//                           kk = kk_start;
//                           k = my_nzs[kk];
//                           done_flag = 1;
//                           while (k < k_end && k >= k_start){
//                              if (nz_used_flags_loc[kk] == 0){ /* has this non-zero been used? */
//                                 int j = col_idx[k];
//                                 int deps_counters_j;
//                                 /* check if x[j] is available (deps_counters[j] must be zero) */
//                                 deps_counters_j = deps_counters[j];
//
//                                 /* if x[j] is available, update x[i] and deps_counters[i] */
//                                 if (deps_counters_j == 0){
//                                    double Tij = mat_values[k];
//                                    double xj = x[j];
//
//                                    x[i] -= xj * Tij / Tii;
//
//                                    num_relax++;
//
//                                    #pragma omp atomic
//                                    deps_counters[i]--;
//
//                                    nz_used_flags_loc[kk] = 1;
//                                 }
//                                 else {
//                                    done_flag = 0;
//                                 }
//                              }
//                              kk++;
//                              k = my_nzs[kk];
//                           }
//                        } while (done_flag == 0); /* loop until x[i] has been completed */
//                     }
//                  }
//               }
//            } 
//
//
//
//
//
//
//
//
//
//
//            /**********************
//             *        
//             *    Atomic CSR
//             *
//             **********************/
//            else {
//               /************
//                * atomic
//                ************/
//               if (atomic_flag == 1){
//                  for (int i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
//                     int num_spins = 0;
//                     double row_wtime_start = omp_get_wtime();
//                     int i = my_rows[i_loc];
//                     int jj_start = start[i];
//                     int jj_end = start[i+1];
//                     int jj_diff = jj_end - jj_start;
//                     int deps_counter_i = jj_diff;
//                     while (deps_counter_i > 0){ /* loop until x[i] has been completed */
//                        int jj_loc_temp = jj_loc;
//                        for (int jj = jj_start; jj < jj_end; jj++){
//                           if (nz_used_flags_loc[jj_loc_temp] == 0){ /* has this non-zero been used? */
//                              int j = col_idx[jj];
//                              int row_solved_flag_j;
//                              /* check if x[j] is available (deps_counters[j] must be zero) */
//#ifdef USE_STDTHREAD
//#else
//                              #pragma omp atomic read
//#endif
//                              row_solved_flag_j = row_solved_flags[j];//row_solved_flag_j = row_solved_flags[j].c;
//
//                              num_spins++;
//
//                              /* if x[j] is available, update x[i] and deps_counters[i] */
//                              if (row_solved_flag_j == 1){
//                                 x[i] -= mat_values[jj] * x[j];
//
//                                 num_relax++;
//
//                                 deps_counter_i--;
//
//                                 nz_used_flags_loc[jj_loc_temp] = 1;
//                              }
//                           }
//                           jj_loc_temp++;
//                        }
//                     }
//                     x[i] /= diag[i];
//#ifdef USE_STDTHREAD
//#else
//                     #pragma omp atomic write
//#endif
//                     row_solved_flags[i] = 1;//row_solved_flags[i].c = 1;
//
//                     jj_loc += jj_diff;
//
//                     double row_wtime_stop = omp_get_wtime();
//                     row_wtime_loc[i_loc] = row_wtime_stop - row_wtime_start;
//                     row_num_spins_loc[i_loc] = num_spins;
//                  }
//               }
//               /************
//                * no atomic
//                ************/
//               else {
//                  for (int i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
//                     int i = my_rows[i_loc];
//                     int jj_start = start[i];
//                     int jj_end = start[i+1];
//                     int jj_diff = jj_end - jj_start;
//                     int deps_counter_i = jj_diff;
//                     while (deps_counter_i > 0){ /* loop until x[i] has been completed */
//                        int jj_loc_temp = jj_loc;
//                        for (int jj = jj_start; jj < jj_end; jj++){
//                           if (nz_used_flags_loc[jj_loc_temp] == 0){ /* has this non-zero been used? */
//                              int j = col_idx[jj];
//                              int row_solved_flag_j;
//                              /* check if x[j] is available (deps_counters[j] must be zero) */
//                              row_solved_flag_j = row_solved_flags[j];//row_solved_flag_j = row_solved_flags[j].c;
//
//                              /* if x[j] is available, update x[i] and deps_counters[i] */
//                              if (row_solved_flag_j == 1){
//                                 x[i] -= mat_values[jj] * x[j];
//
//                                 num_relax++;
//
//                                 deps_counter_i--;
//
//                                 nz_used_flags_loc[jj_loc_temp] = 1;
//                              }
//                           }
//                           jj_loc_temp++;
//                        }
//                     }
//                     x[i] /= diag[i];
//                     #pragma omp atomic write
//                     row_solved_flags[i] = 1;//row_solved_flags[i].c = 1;
//
//                     jj_loc += jj_diff;
//                  }
//               }
//            }
//         }
//      }
//
//      wtime_stop = omp_get_wtime();
//      ts->output.solve_wtime_vec[tid] = wtime_stop - wtime_start;
//
//      ts->output.num_relax[tid] = num_relax;
//      ts->output.num_iters[tid] = num_iters;
//
//      for (i_loc = 0; i_loc < n_loc; i_loc++){
//         int i = my_rows[i_loc];
//         ts->output.row_output[i].wtime     += row_wtime_loc[i_loc];
//         ts->output.row_output[i].num_spins += row_num_spins_loc[i_loc];
//         int jj_start = start[i];
//         int jj_end = start[i+1];
//         for (int jj = jj_start; jj < jj_end; jj++){
//#ifdef USE_STDTHREAD
//#else
//            #pragma omp atomic
//#endif
//            ts->output.row_output[col_idx[jj]].num_atomics++;
//         }
//      }
//
//      if (ts->input.MsgQ_flag == 1){
//         for (i_loc = 0; i_loc < n_loc; i_loc++){
//            int i = my_rows[i_loc];
//            x[i] = x_loc[i_loc];
//         }
//
//         if (ts->input.MsgQ_wtime_flag == 1){
//            ts->output.MsgQ_put_wtime_vec[tid] = MsgQ_put_wtime;
//            ts->output.MsgQ_get_wtime_vec[tid] = MsgQ_get_wtime;
//         }
//         else if (ts->input.MsgQ_cycles_flag == 1){
//            ts->output.MsgQ_put_cycles_vec[tid] = MsgQ_put_cycles;
//            ts->output.MsgQ_get_cycles_vec[tid] = MsgQ_get_cycles;
//         }
//         else if (ts->input.comp_wtime_flag == 1){
//            ts->output.comp_wtime_vec[tid] = comp_wtime;
//         }
//
//         ts->output.num_qGets_vec[tid] = num_qGets;
//         ts->output.num_qPuts_vec[tid] = num_qPuts;
//
//         PrintDummy(dummy);
//      }
//   }
//
//   
//}
