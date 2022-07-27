#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "Parallel.hpp"

typedef struct {
   std::vector<double> *b;
   std::vector<double> *x;
   int proc_id;
} SolverParArg;

class Solver
{

public:
   Solver(unsigned int num_procs_input = 1):
      num_procs(num_procs_input)
   {};
   ~Solver(void){};

   virtual void Setup(OutputData &output_data_ptr) = 0;
   virtual void InitSolveData(void) = 0;
   virtual void Solve(std::vector<double> b,
                      std::vector<double> &x,
                      OutputData &output_data_input) = 0;

protected:
   virtual void SequentialSolve(void) = 0;
   //virtual void ParallelSolve(void) = 0;
   //virtual void *ParallelSolveVoidStar(void *arg) = 0;

   OutputData *output_data_ptr;
   std::vector<double> *x_ptr, *b_ptr;
   unsigned int num_procs;
#ifdef USE_STDTHREAD
   //PthreadInfo pthread_info;
   STDthreadInfo thread_info;
#else
#endif
};



class SparseSolver : public Solver
{

public:
   SparseSolver(SparseMatrix &A_input,
                unsigned int num_procs_input):
                A(&A_input),
                Solver(num_procs_input)
                {};
   virtual ~SparseSolver() = 0;

protected:
   SparseMatrix *A;
};

inline SparseSolver::~SparseSolver(void) {};



class IterSparseSolver : public SparseSolver 
{

public:
   using SparseSolver::SparseSolver;
   virtual ~IterSparseSolver() = 0;

   void SetTol(double tol_new) {tol = tol_new;};
   void SetNumIters(unsigned int num_iters_new) {num_iters = num_iters_new;};

protected:
   double tol;
   unsigned int num_iters;
};

#endif
