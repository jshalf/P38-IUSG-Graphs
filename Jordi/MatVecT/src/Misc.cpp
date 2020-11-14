#include "Misc.hpp"

double RandDouble(double low, double high)
{
   return low + (high - low) * ((double)rand() / RAND_MAX);
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
