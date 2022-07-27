#include "Parallel.hpp"

using namespace std;

/********************************
 * Communicator member functions
 ********************************/
Communicator::Communicator(unsigned int num_procs,
                           CommunicatorInput comm_input):
   num_procs(num_procs),
   comm_type(comm_input.comm_type)
{

}

Communicator::~Communicator(void)
{

}

CommunicationType Communicator::GetCommType(void)
{
   return comm_type;
}

/*******************************
 * Partitioner member functions
 *******************************/

Partitioner::Partitioner(unsigned int num_parts, 
                         PartitionerInput part_input):
   num_parts(num_parts),
   idx_start(part_input.idx_start),
   size_glob(part_input.size_glob),
   idx_glob(part_input.idx_glob)
{

}

Partitioner::~Partitioner(void)
{

}

vector<unsigned int> Partitioner::GetPartitionSize(void)
{
   return size;
}

vector<vector<unsigned int>> Partitioner::GetPartition(void)
{
   return idx;
}

vector<unsigned int> Partitioner::GetIndexToProcMap(void)
{
   return idx_to_proc_map;
}

vector<unordered_map<unsigned int, unsigned int>> Partitioner::GetGlobToLocIndexMap(void)
{
   return glob_to_loc_idx_map;
}

void Partitioner::ConstructPartition(void)
{
   size.resize(num_parts, 0);
   idx.resize(num_parts, vector<unsigned int>());
   for (int p = 0; p < num_parts; p++){
      if (p < size_glob % num_parts){
         size[p] = size_glob / num_parts + 1;
      }
      else {
         size[p] = size_glob / num_parts;
      }
      /* round-robin */
      idx[p].resize(size[p], 0);
      for (int i = 0; i < size[p]; i++){
         int ii = idx_start + num_parts * i + p;
         idx[p][i] = idx_glob[ii];
      }
   }
}

void Partitioner::ConstructIndexToProcMap(void)
{
   /* round-robin */
   idx_to_proc_map.resize(size_glob);
   unsigned int p = 0;
   for (int ii = 0; ii < size_glob; ii++){
      unsigned int i = idx_glob[ii];
      idx_to_proc_map[i] = p;
      if (p == num_parts-1){
         p = 0;
      }
      else {
         p++;
      }
   }
}

void Partitioner::ConstructGlobToLocIndexMap(void)
{
   glob_to_loc_idx_map.resize(num_parts, unordered_map<unsigned int, unsigned int>());
   for (int p = 0; p < num_parts; p++){
      for (int i_loc = 0; i_loc < size[p]; i_loc++){
         int i = idx[p][i_loc];
         glob_to_loc_idx_map[p][i] = i_loc;
      }
   }
}

/*******************************
 * ParallelInfo member functions
 *******************************/
ParallelInfo::ParallelInfo(unsigned int num_procs,
                           ParallelInput para_input): 
   num_procs(num_procs),
   part(new Partitioner(num_procs, para_input.part_input)),
   comm(new Communicator(num_procs, para_input.comm_input))
{
   //part->ConstructPartition();
   //part->ConstructIndexToProcMap();
   //part->ConstructGlobToLocIndexMap();
//#ifdef USE_STDTHREAD
//   pthread_info.threads.resize(num_procs);
//#else
//
//#endif
}

ParallelInfo::~ParallelInfo(void)
{

}

Partitioner *ParallelInfo::Part(void)
{
   return part;
}

Communicator *ParallelInfo::Comm(void)
{
   return comm;
}

unsigned int ParallelInfo::GetNumProcs(void)
{
   return num_procs;
}
