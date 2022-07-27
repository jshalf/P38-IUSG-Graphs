#ifndef LEVEL_SETS_HPP
#define LEVEL_SETS_HPP

#include "Main.hpp"
#include "Matrix.hpp"

/* Level set data used in level scheduling algorithms */
typedef struct{
   vector<unsigned int> order; /* row ordering (rows are ordered by level)*/
   vector<unsigned int> level_size; /* number of rows in each level */
   vector<unsigned int> level_start; /* starting point in ``perm'' for each level */
   int num_levels; /* number of level sets */
}LevelSetData;

void LevelSets(SparseMatrix A, /* matrix data (input) */
               LevelSetData *lvl_set, /* level set data (output) */
               int L_flag); /* is the matrix lower or upper triangular ? */

void LevelSetsDestroy(LevelSetData *lvl_set);


#endif
