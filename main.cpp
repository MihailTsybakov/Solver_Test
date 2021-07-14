#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "mpi.h"
#include "mkl.h"
#include "mkl_cluster_sparse_solver.h"

void MPI_test_invoke(int argc = 0, char* argv[] = nullptr)
{
	int rank, size, name_len;
	char machine_name[256];

	std::cout << " <logs> Invoking MPI..." << std::endl;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Get_processor_name(machine_name, &name_len);

	std::cout << " <logs> MPI initialized." << std::endl;
	std::cout << " <logs> Current process rank: " << rank << std::endl;
	std::cout << " <logs> Current process group size: " << size << std::endl;
	std::cout << " <logs> Current machine name: " << machine_name << std::endl;

	MPI_Finalize();
}

void check_solver_errors(int error_code)
{
	std::map<int, std::string> errors;

	errors[0] = " <logs> No solver errors occured.";
	errors[-1] = " <logs> Solver error: inconsistent input.";
	errors[-2] = " <logs> Solver error: not enough memory.";
	errors[-3] = " <logs> Solver error: reordering failed.";
	errors[-4] = " <logs> Solver error: zero pivot: try to change pivoting pertuberation.";
	errors[-5] = " <logs> Solver error: unclassified internal error.";
	errors[-6] = " <logs> Solver error: reordering failed (matrix types 11 & 13).";
	errors[-7] = " <logs> Solver error: diagonal matrix is singular.";
	errors[-8] = " <logs> Solver error: 32-bit integer overflow encountered.";
	errors[-9] = " <logs> Solver error: not enough memory for OOC.";
	errors[-10] = " <logs> Solver error: OOC opening failed.";
	errors[-11] = " <logs> Solver error: OOC reading/writing failed.";

	if (error_code > 0 || error_code < -11)
	{
		std::cout << " <logs> Error: unknown error code encountered." << std::endl;
	}
	else
	{
		std::cout << errors[error_code] << std::endl;
	}
}

void pack_params(int* solver_params, bool default_settings)
{
	std::vector<int> setzero_indexes = { 2, 3, 8, 18, 19, 23, 24, 25, 28, 31, 32, 33, 37, 38,
	42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 60, 61, 63 };

	for (auto index : setzero_indexes) solver_params[index] = 0;

	if (default_settings)
	{
		solver_params[0] = 0;
		return;
	}

	solver_params[0] = 1;    // Non-Zero value ==> Non-Default settings
	solver_params[1] = 2;    // MPI nested dissection
	solver_params[4] = 0;    // Ignoring user permutations
	solver_params[5] = 0;    // Writing result to x vector
	solver_params[6] = 0;    // Output value
	solver_params[7] = 2;    /// 2 refinement steps  ~!!
	solver_params[9] = 13;   // Pivoting perturbations for non-symmetric matricies, 1e-13
	solver_params[10] = 1;   // Vector scaling enabled - default for non-symmetric matricies
	solver_params[11] = 0;   // No conjugation/transposition
	solver_params[12] = 1;   // Matching enabled - default for non-symmetric matricies
	solver_params[13] = 0;   // Output value
	solver_params[14] = 0;   // Output value
	solver_params[15] = 0;   // Output value
	solver_params[16] = 0;   // Output value
	solver_params[17] = 0;   // Disable non-zero factors reporting
	solver_params[20] = 0;   // 1x1 Pivoting for symmetric matricies
	solver_params[21] = 0;   // Output value
	solver_params[22] = 0;   // Output value
	solver_params[26] = 1;   // Check inputted matrix correctness enabled
	solver_params[27] = 1;   /// Single precision calculations ~!
	solver_params[29] = 0;   // Output value
	solver_params[30] = 0;   /// Partial solve disabled ~!!
	solver_params[34] = 1;   // C-Style indexing
	solver_params[35] = 0;   /// Schur complementation is not calculated ~!
	solver_params[36] = 0;   // CSR Matrix storage format
	solver_params[39] = 0;   /// a-Matrix is provided from 0-hostprocess ~!!
	solver_params[40] = -1;  // Not actual for current format settings
	solver_params[41] = -1;  // Not actual for current format settings
	solver_params[59] = 2;   /// Cluster sparse solver mode set to OOC MODE ~!!!
	solver_params[62] = 0;   // Output value
}

int main()
{
	bool logs = true;
	size_t matrix_size = 4;			 // NxN Matrix size = equations number

	void* pt[64] = {0};			 // Dont modify after first cpardiso call
	MKL_INT maxfct = 1;			 // Factorizations
	MKL_INT mnum = 1;			 // Ignored
	MKL_INT mtype = 11;			 // Real non-symmetric matrix
	MKL_INT phase = 13;			 // Full solving process
	MKL_INT  n = 4;				 /// Equations number ~!!
	float* a;                        	 // A Matrix
	MKL_INT* ia;				 // ia index array
	MKL_INT* ja;				 // ja index array
	MKL_INT service_int;			 // Ignored?
	MKL_INT nrhs = 1;		         // Number of right-hand sides
	MKL_INT* iparm;				 // Solver parameters
	MKL_INT msglvl = 1;			 // Report messages ON
	float* b;				 // b vector
	float* x;				 // x vector
	MKL_INT error = 0;			 // Error Status

	a = new float[5];  
	b = new float[4]; 
	x = new float[4]; 
	iparm = new MKL_INT[64]; 
	ia = new MKL_INT[5]; 
	ja = new MKL_INT[5]; 

	//    -======= Matrix Description ========-
	/// =========================================
	x[0] = 0; x[1] = 0; x[2] = 0; x[3] = 0; 

	//			Some random matrix
	a[0] = 1; a[1] = 5; a[2] = 7; a[3] = 6; a[4] = 9;
	b[0] = 1; b[1] = 1; b[2] = 1; b[3] = 1;
	ja[0] = 0; ja[1] = 2; ja[2] = 1; ja[3] = 1; ja[4] = 3;
	ia[0] = 0; ia[1] = 1; ia[2] = 2; ia[3] = 3; ia[4] = 5;

	//       Sample matrix from intel example
	
	//MKL_INT n = 5;
	//MKL_INT* ia = new MKL_INT[6];
	//ia[0] = 0; ia[1] = 3; ia[2] = 5; ia[3] = 8; ia[4] = 11; ia[5] = 13;
	//b[0] = 1; b[1] = 1; b[2] = 1; b[3] = 1; b[4] = 1;

	////MKL_INT* ja = new MKL_INT[13];
	//ja[0] = 0; ja[1] = 1; ja[2] = 3;
	//ja[3] = 0; ja[4] = 1;
	//ja[5] = 1; ja[6] = 2; ja[7] = 3;
	//ja[8] = 0; ja[9] = 2; ja[10] = 3;
	//ja[11] = 1; ja[12] = 4;

	////double* a = new double[13];
	//a[0] = 1; a[1] = -1; a[2] = -3;
	//a[3] = -2; a[4] = 5;
	//a[5] = 4; a[6] = 6; a[7] = 4;
	//a[8] = -4; a[9] = 2; a[10] = 7;
	//a[11] = 8; a[12] = -5;

	//			   Unit Matrix
	/*a[0] = 1; a[1] = 1; a[2] = 1; a[3] = 1;
	b[0] = 1; b[1] = 1; b[2] = 1; b[3] = 1;
	ja[0] = 0; ja[1] = 1; ja[2] = 2; ja[3] = 3;
	ia[0] = 0; ia[1] = 1; ia[2] = 2; ia[3] = 3; ia[4] = 4;*/
	/// =========================================

	pack_params(iparm, false);

	int rank, size, comm, argc = 0;
	char** argv;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	comm = MPI_Comm_c2f(MPI_COMM_WORLD);

	if (!rank && logs) std::cout << " <logs> MPI Initialized" << std::endl;

	if (!rank && logs) std::cout << " <logs> Invoking cluster solver... " << std::endl;
	
	/// Symbolic factorization & reordering
	phase = 11;
	if (!rank && logs) std::cout << " <logs> Phase 1: Symbolic factorization" << std::endl;
	cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &service_int, &nrhs, iparm, &msglvl, b, x, &comm, &error);
	if (!rank) check_solver_errors(error);

	if (error != 0)
	{
		if (!rank) std::cout << " <logs> Shutting down" << std::endl;
		phase = -1;
		// Releasing memory
		cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &service_int, &nrhs, iparm, &msglvl, b, x, &comm, &error);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD, -1);
		delete[] a; delete[] b; delete[] x; delete ia; delete[] ja; delete[] iparm;
		return -1;
	}

	/// Numerical Factorization
	phase = 22;
	if (!rank && logs) std::cout << " <logs> Phase 2: Numerical factorization" << std::endl;
	cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &service_int, &nrhs, iparm, &msglvl, b, x, &comm, &error);
	if (!rank) check_solver_errors(error);

	if (error != 0)
	{
		if (!rank) std::cout << " <logs> Shutting down" << std::endl;
		phase = -1;
		// Releasing memory
		cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &service_int, &nrhs, iparm, &msglvl, b, x, &comm, &error);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD, -1);
		delete[] a; delete[] b; delete[] x; delete ia; delete[] ja; delete[] iparm;
		return -1;
	}

	/// Solving
	phase = 33;
	if (!rank && logs) std::cout << " <logs> Phase 3: Solving sparse system" << std::endl;
	cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &service_int, &nrhs, iparm, &msglvl, b, x, &comm, &error);
	if (!rank) check_solver_errors(error);

	if (error != 0)
	{
		if (!rank) std::cout << " <logs> Shutting down" << std::endl;
		phase = -1;
		// Releasing memory
		cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &service_int, &nrhs, iparm, &msglvl, b, x, &comm, &error);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD, -1);
		delete[] a; delete[] b; delete[] x; delete ia; delete[] ja; delete[] iparm;
		return -1;
	}

	// Checking result:
	if (!rank && logs)
	{
		std::cout << " SOLUTION CHECK: " << std::endl;
		for (int i = 0; i < n; ++i)
		{
			std::cout << " x[" << i << "] = " << x[i] << std::endl;
		}
	}
	// Releasing solver memory
	phase = -1;
	cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &service_int, &nrhs, iparm, &msglvl, b, x, &comm, &error);

	delete[] ia;
	delete[] ja;
	delete[] x;
	delete[] b;
	delete[] a;
	delete[] iparm;

	MPI_Finalize();
	if (rank == 0 && logs) std::cout << " <logs> MPI shut down." << std::endl;

	return 0;
}
