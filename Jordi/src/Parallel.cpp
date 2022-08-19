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
   idx_glob(part_input.idx_glob),
   lump(part_input.lump)
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

   /* round-robin */
   unsigned int p = 0;
   unsigned int j = 0;
   for (int ii = idx_start; ii < size_glob; ii++){
      idx[p].push_back(idx_glob[ii]);

      if (j == lump-1){
         if (p == num_parts-1){
            p = 0;
         }
         else{
            p++;
         }
      }
      if (j == lump-1){
         j = 0;
      }
      else {
         j++;
      }
   }

   for (int p = 0; p < num_parts; p++){
      size[p] = idx[p].size();
   }

   //for (int p = 0; p < num_parts; p++){
   //   printf("%d:", p);
   //   for (int j = 0; j < size[p]; j++){
   //      printf(" %d", idx[p][j]);
   //   }
   //   printf("\n");
   //}
}

void Partitioner::ConstructIndexToProcMap(void)
{
   /* round-robin */
   idx_to_proc_map.resize(size_glob);
   unsigned int p = 0;
   unsigned int j = 0;
   for (int ii = idx_start; ii < size_glob; ii++){
      unsigned int i = idx_glob[ii];
      idx_to_proc_map[i] = p;

      if (j == lump-1){
         if (p == num_parts-1){
            p = 0;
         }
         else {
            p++;
         }
      }
      if (j == lump-1){
         j = 0;
      }
      else {
         j++;
      }
   }

   //for (int ii = 0; ii < size_glob; ii++){
   //   unsigned int i = idx_glob[ii];
   //   printf("%d\n", idx_to_proc_map[i]);
   //}
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
