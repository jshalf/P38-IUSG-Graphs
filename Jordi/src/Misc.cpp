#include "Misc.hpp"

using namespace std;

double RandDouble(double low, double high)
{
   return low + (high - low) * ((double)rand() / RAND_MAX);
}

double RandInt(int low, int high, double seed)
{
   //std::random_device rd;
   //std::mt19937 gen(rd()); // seed the generator
   //std::uniform_int_distribution<> distr(low, high); // define the range
   //return distr(gen);

   return (int)round((rand() % (high - low + 1)) + low);
}

double InnerProd(double *x, double *y, int n)
{
   double inner_prod = 0;
   #pragma omp parallel for reduction(+:inner_prod)
   for (int i = 0; i < n; i++){
      inner_prod += x[i] * y[i];
   }

   return inner_prod;
}

double Residual2Norm(CSR A, double *x, double *b)
{
   int n = A.n;
   double r_2norm = 0, b_2norm = 0;
   #pragma omp parallel for reduction(+:r_2norm,b_2norm)
   for (int i = 0; i < n; i++){
      double res = b[i];
      for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++)
      {
         int ii = A.j[jj];
         res -= A.data[jj] * x[ii];
      }
      r_2norm += res*res;
      b_2norm += b[i]*b[i];
   }

   return sqrt(r_2norm)/sqrt(b_2norm);
}

void LevelSets(CSR A, LevelSetData *lvl_set)
{
   int n = A.n;
   vector<int> nodes(n);
   vector<int>::iterator it;
   lvl_set->perm = (int *)calloc(n, sizeof(int));
   int *nonzero_flags = (int *)calloc(n, sizeof(int));
   int *nonzero_flags_prev = (int *)calloc(n, sizeof(int));

   for (int i = 0; i < n; i++){
      nonzero_flags[i] = 1;
      nodes[i] = i;
   }

   int k = 0;
   while(nodes.size() > 0){
      int level_size = 0;
      for (int i = 0; i < n; i++){
         nonzero_flags_prev[i] = nonzero_flags[i];
      }
      for(it = nodes.begin(); it != nodes.end(); it++){
         int i = *it;
         int root_flag = 1;
         for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
            if (nonzero_flags_prev[A.j[jj]] == 1){
               root_flag = 0;
               break;
            }
         }
         if (root_flag == 1){
            nonzero_flags[i] = 0;
            lvl_set->perm[k] = i;
            k++;
            level_size++;
            nodes.erase(it--);
         }
      }
      lvl_set->level_size.push_back(level_size);
   }

   lvl_set->num_levels = lvl_set->level_size.size();
   lvl_set->level_start.resize(lvl_set->num_levels+1);
   lvl_set->level_start[0] = 0;
   for (int i = 0; i < lvl_set->num_levels; i++){
      lvl_set->level_start[i+1] = lvl_set->level_start[i] + lvl_set->level_size[i];
   } 

   free(nonzero_flags);
}
