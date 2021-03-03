	----------
	 Overview
	----------
This code contains a collection of benchmarks that test the scalability
of atomic operations.  While there are different options that can be set
within each benchamrk, each benchmark includes three important ingredients:
	1. A baseline method that may not require atomics.
	2. A method that requires atomics.
	3. The ability to turn off the atomics in the atomics method in order to measure the performance difference.
In the third case, note that removing atomics may not give the correct result and
is merely there to measure the potential performance gain if the atomics 
had no cost.

For the sparse linear algebra bechmarks, there are three available test problems:
	1. The five-point centered difference discretization of the Poisson equation.
	2. Random symmetric matrix.
	3. Matrix read from file.
The matrices are stored in either compressed sparse row (CSR) or coordinate format (COO).


	-------
	 Build
	-------
To build a benchmark, navigate to the benchmark's directory 
and use make.



	-------------------
	 MatVecT Benchmark
	-------------------
In this benchmark, a sparse matrix-transpose vector product is computed, i.e., y = A^Tx,
without explicitely forming the matrix transpose.  Multiple OpenMP threads are used for parallelism.
The baseline method has each thread write to a local output vector all of which are 
summed to obtain the correct result.
The baseline is does not require atomics.
The atomics method can be thought of as computing a standard matrix-vector product
(no transpose) using a matrix in compressed sparse column (CSC) format.

Example run:
	
	./main -num_threads 8 -problem 5pt -n 100

This example uses the atomics method, 8 threads, and a five-point matrix with 100x100 rows.



	-----------------------
	 AsyncJacobi Benchmark
	-----------------------
In this benchmark, the solution of a linear system is computed using the Jacobi method with multiple OpenMP threads.
The baseline method is the standard synchronous Jacobi method and does not require atomic operations.
The atomics method is the asynchronous Jacobi method.

Example run:
	
	./main -num_threads 8 -solver aj -problem 5pt -n 1000

This example uses the asynchronous Jacobi method, 8 threads, and a five-point matrix with 100x100 rows.


	--------------------
	 TriSolve Benchmark
	--------------------
In this benchmark, the solution to a triangular linear system is computed using multiple OpenMP threads.
The baseline method is the standard level-scheduled algorithm which is synchronous and does not require atomic operations.
The atomics method is an asynchronous iterative fine-grained method.

Example run:
	
	./main -num_threads 8 -solver async -problem rand -n 1000 -mxr_nnz 10

This example uses the asynchronous atomics method, 8 threads, and a random matrix with 1000 rows and a maximum of
10 non-zero values per row.