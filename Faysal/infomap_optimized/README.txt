This is an optimized version of the Infomap code that in the other repository

Optimization Include:
==========================

	1. Fully parallelized (MPI + OpenMP) version of the community optimization kernel (prioritizeMove, prioritizeMoveSpNode)

	2. ConvertModulestoSuperNode parallelized

	3. Timing codes are changed with <chrono> library instead of high overhead MPI_Wtime() function.

	4. Function parameters are cleaned up with global extern variables across multiple files.

	5. Which function takes how much time and called how many times in a particular run is reported.

	6. Further optimization may be possible


Known Issue:
==========================

	1. There seems to be some accuracy issue (modularity and conductance) for running code in very large number of processing cores (e.g. more than thousands core). The reason is unknown. An investigation needs to be made after the synchronization operation to see every MPI process has the uniform community assignment for each vertex.


How to compile in cori KNL:
===========================

	1. navigate to the directory where the code is saved

	2. issue the command "module load metis"

	3. issue the command "module swap craype-haswell craype-mic-knl"

	4. issue the command "make"
