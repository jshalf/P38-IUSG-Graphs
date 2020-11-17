#include "Misc.hpp"

double Residual2Norm(CSR A, double *x, double *b)
{
   int n = A.n;
   double r_2norm = 0, b_2norm = 0;
   #pragma omp parallel for reduction(+:r_2norm,b_2norm)
   for (int i = 0; i < n; i++){
      double res = b[i];
      for (int jj = A.i[i]; jj < A.i[i+1]; jj++)
      {
         int ii = A.j[jj];
         res -= A.data[jj] * x[ii];
      }
      r_2norm += res*res;
      b_2norm += b[i]*b[i];
   }
   
   return sqrt(r_2norm)/sqrt(b_2norm);
}

double RandDouble(double low, double high)
{
   return low + (high - low) * ((double)rand() / RAND_MAX);
}
