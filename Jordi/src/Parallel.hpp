#ifndef PARALLEL_HPP
#define PARALLEL_HPP

#include "Main.hpp"

template<typename T>
struct AtomicType {
   std::atomic<T> val;
};

enum class PartitionType {
   round_robin,
   metis
};

enum class CommunicationType {
   atomic,
   MsgQ
};

enum class ParallelProgram {
   OpenMP,
   pthreads,
   devastator,
   MPI
};

typedef struct {
   unsigned int idx_start;
   unsigned int size_glob;
   vector<unsigned int> idx_glob;
   unsigned int lump = 1;
} PartitionerInput;

typedef struct {
   CommunicationType comm_type;
} CommunicatorInput;

typedef struct {
   CommunicatorInput comm_input;
   PartitionerInput part_input;
} ParallelInput;

#ifdef USE_PTHREADS
typedef struct{
   vector<pthread_t> threads;
   pthread_barrier_t barrier;
}PthreadInfo;

typedef struct{
   vector<std::thread *> th;
   //std::barrier *bar;
}STDthreadInfo;
#else

#endif

class Communicator
{
public:
   Communicator(unsigned int num_procs,
                CommunicatorInput input);   
   ~Communicator(void);

   CommunicationType GetCommType(void);

protected:
   CommunicationType comm_type;
   unsigned int num_procs;
};

class Partitioner
{
public:
   Partitioner(unsigned int num_parts,
               PartitionerInput input);
   ~Partitioner(void);

   vector<unsigned int> GetPartitionSize(void);
   vector<vector<unsigned int>> GetPartition(void);
   vector<unsigned int> GetIndexToProcMap(void);
   vector<unordered_map<unsigned int, unsigned int>> GetGlobToLocIndexMap(void);

   void ConstructPartition(void);
   void ConstructIndexToProcMap(void);
   void ConstructGlobToLocIndexMap(void);

protected:
   /* input values */
   unsigned int size_glob;
   unsigned int num_parts;
   vector<unsigned int> idx_glob;
   unsigned int idx_start;
   unsigned int lump = 1;

   /* computed by partitioner */
   vector<unsigned int> size;
   vector<vector<unsigned int>> idx;
   vector<unsigned int> idx_to_proc_map;
   vector<unordered_map<unsigned int, unsigned int>> glob_to_loc_idx_map;
};

class ParallelInfo 
{
public:
   ParallelInfo(unsigned int num_procs,
                ParallelInput para_input);
   ~ParallelInfo(void);

   Partitioner *Part(void);
   Communicator *Comm(void);
   unsigned int GetNumProcs(void);
#ifdef USE_PTHREADS
   PthreadInfo GetPthreadInfo(void)
   {
      return pthread_info;
   }
#else
#endif
   
protected:
   Partitioner *part;
   Communicator *comm;
   unsigned int num_procs;
   ParallelProgram para_program;
#ifdef USE_PTHREADS
   PthreadInfo pthread_info;
#else
#endif
};

#endif
