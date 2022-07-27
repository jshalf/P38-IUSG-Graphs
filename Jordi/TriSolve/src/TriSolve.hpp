#ifndef TRISOLVE_HPP
#define TRISOLVE_HPP

#include "Main.hpp"
#include "Matrix.hpp"
#include "Parallel.hpp"
#include "MsgQ.hpp"
#include "Solver.hpp"
#include "LevelSets.hpp"


#define TRISOLVE_OMPFOR_SCHED static

/* types of trisolve algorithms */
#define TRISOLVE_ASYNC 0
#define TRISOLVE_LEVEL_SCHEDULED 1
#define TRISOLVE_ASYNC_LEVEL_SCHEDULED 2

/* Struct used by TriSolve benchmark */
typedef struct{
   InputData input;
   OutputData output;
   LevelSetData L_lvl_set; /* level set data for lower triangular part */
   LevelSetData U_lvl_set; /* level set data for upper triangular part */
#ifdef USE_STDTHREAD
   PthreadInfo pt;
#endif
}TriSolveData;

//typedef struct {
//   Cache<int> src_row;
//   Cache<int> dst_row;
//   Cache<double> data;
//} TriSolveCSRMessage;

typedef struct {
   int src_row;
   int dst_row;
   double data;
} TriSolveCSRMessage;

enum class TriSolverMethod {
   sequential,
   sync,
   async,
   async_fine_grained
};

enum class TriSolverSolveOrder {
   natural,
   level_sched
};

typedef struct TriSolveParArgStruct : SolverParArg {
   ;
} TriSolveParArg;

class SparseTriSolver : public SparseSolver
{
public:
   //SparseTriSolver(void);
   using SparseSolver::SparseSolver;
   ~SparseTriSolver(void) {};

   void Setup(OutputData &output_data_ptr);
   void InitSolveData(void);
   void Solve(std::vector<double> b,
              std::vector<double> &x,
              OutputData &output_data_input);

protected:
   void SequentialSolve(void);
   virtual void ParallelSolveFunc(TriSolveParArg *arg) = 0;

   /* order in which equations are solved */
   std::vector<unsigned int> idx_solve_order;
};

class LevelSchedTriSolver : public SparseTriSolver
{
public:
   //LevelSchedTriSolver();
   using SparseTriSolver::SparseTriSolver;
   ~LevelSchedTriSolver() {};

   void Setup(OutputData &output_data_ptr);
   void InitSolveData(void);
   using SparseTriSolver::Solve;

   vector<ParallelInfo *> GetParInfo(void) {
      return level_para_info;
   }

   LevelSetData GetLevelSetInfo(void) {
      return level_set_info;
   }

protected:
   using SparseTriSolver::SequentialSolve;
   //void *ParallelSolveVoidStar(void *arg) override;
   void ParallelSolveFunc(TriSolveParArg *arg);

   LevelSetData level_set_info;
   vector<ParallelInfo *> level_para_info;
};

class AsyncTriSolver : public SparseTriSolver
{
public:
   using SparseTriSolver::SparseTriSolver;
   ~AsyncTriSolver() {};

   void Setup(OutputData &output_data_ptr,
              TriSolverSolveOrder solve_order_input = TriSolverSolveOrder::natural,
              CommunicationType comm_type_input = CommunicationType::atomic);
   void InitSolveData(void);
   using SparseTriSolver::Solve;

   ParallelInfo* GetParInfo(void) {
      return para_info;
   }

   LevelSetData GetLevelSetInfo(void) {
      return level_set_info;
   }

protected:
   using SparseTriSolver::SequentialSolve;

   //void *ParallelSolveVoidStar(void *arg) override;
   //void *ParallelSolveVoidStar_Atomic(void *arg);
   //void *ParallelSolveVoidStar_MsgQ(void *arg);

   void ParallelSolveFunc(TriSolveParArg *arg);
   void ParallelSolveFunc_Atomic(TriSolveParArg *arg);
   void ParallelSolveFunc_MsgQ(TriSolveParArg *arg);

   LevelSetData level_set_info;
   ParallelInfo *para_info;
   TriSolverSolveOrder solve_order;
#ifdef USE_STDTHREAD
   //vector<AtomicType<unsigned short>> row_solved;
   std::atomic<unsigned short> *row_solved;
#else
   vector<unsigned short int> row_solved;
#endif

   vector<vector<unsigned int>> qPut_dst_rows;
#if USE_DEVA
   MessageQueue<TriSolveCSRMessage> *Q;
#else
   vector<MessageQueue<TriSolveCSRMessage>*> Q;
#endif
};

//void TriSolve_Seq(TriSolveData *ts,
//                  SparseMatrix A,
//                  int *T_perm,
//                  double *x,
//                  double *b);
//
//void TriSolve(TriSolveData *ts,
//              LevelSetData lvl_set,
//              SparseMatrix A,
//              double *x,
//              double *b);

#endif
