#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "mpi.h"
#include "mkl.h"
#include "mkl_cluster_sparse_solver.h"

#include "matrixio.h"

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
	solver_params[27] = 2;   /// Double precision calculations ~!
	solver_params[29] = 0;   // Output value
	solver_params[30] = 0;   /// Partial solve disabled ~!!
	solver_params[34] = 1;   // C-Style indexing
	solver_params[35] = 0;   /// Schur complementation is not calculated ~!
	solver_params[36] = 0;   // CSR Matrix storage format
	solver_params[39] = 0;   /// a-Matrix is provided from 0-hostprocess ~!!
	solver_params[40] = -1;  // Not actual for current format settings
	solver_params[41] = -1;  // Not actual for current format settings
	solver_params[59] = 1;   /// Cluster sparse solver mode set to OOC MODE ~!!!
	solver_params[62] = 0;   // Output value
}

int main()
{
	bool logs = true;
	size_t matrix_size = 4;		 // NxN Matrix size = equations number

	void* pt[64] = {0};			 // Dont modify after first cpardiso call
	MKL_INT maxfct = 1;			 // Factorizations
	MKL_INT mnum = 1;			 // Ignored
	MKL_INT mtype = 11;			 // Real non-symmetric matrix
	MKL_INT phase = 13;			 // Full solving process
	MKL_INT  n = 4;				 /// Equations number ~!!
	double* a;                   // A Matrix
	MKL_INT* ia;				 // ia index array
	MKL_INT* ja;				 // ja index array
	MKL_INT service_int;		 // Ignored?
	MKL_INT nrhs = 1;	         // Number of right-hand sides
	MKL_INT* iparm;				 // Solver parameters
	MKL_INT msglvl = 1;			 // Report messages ON
	double* b, *b_test;			 // b vector
	double* x;					 // x vector
	MKL_INT error = 0;			 // Error Status

	iparm = new MKL_INT[64];

	//    -======= Matrix Description ========-
	/// =========================================
	std::string source_path = "repository\\Testing\\Examples\\fidesys_static_3d";

	std::ifstream B_vector;
	std::ifstream X_vector;
	std::ifstream A_matrix;

	B_vector.open(source_path + "\\B.vec", std::ios::in);
	X_vector.open(source_path + "\\X.vec", std::ios::in);
	A_matrix.open(source_path + "\\A.txt", std::ios::in);

	if (!B_vector.is_open() || !X_vector.is_open() || !A_matrix.is_open())
	{
		std::cout << " <logs> Error: cannot open one or several files in " << source_path << std::endl;
		return -1;
	}

	MKL_INT A_size, A_nonzero;

	std::vector<double> A;
	std::vector<double> B;
	std::vector<double> X_true;
	std::vector<MKL_INT> IA;
	std::vector<MKL_INT> JA;

	read_csr_matrix_file<double>(A_matrix, A_size, A_nonzero, IA, JA, A);

	MKL_INT B_size, X_size;
	B_vector >> B_size;
	X_vector >> X_size;
	
	read_vector<double>(B_vector, B_size, B);
	read_vector<double>(X_vector, X_size, X_true);

	B_vector.close();
	X_vector.close();
	A_matrix.close();
	
	a = new double[A_nonzero];
	b = new double[B_size];
	b_test = new double[B_size];
	x = new double[X_size];
	ia = new MKL_INT[IA.size()];
	ja = new MKL_INT[JA.size()];

	// Optimization needed
	for (int i = 0; i < A_nonzero; ++i) a[i] = A[i];
	for (int i = 0; i < B_size; ++i) b[i] = B[i];
	for (size_t i = 0; i < IA.size(); ++i) ia[i] = IA[i];
	for (size_t i = 0; i < JA.size(); ++i) ja[i] = JA[i];
	for (int i = 0; i < X_size; ++i) x[i] = X_true[i];

	n = B_size;

	/// =========================================


	pack_params(iparm, false);

	int rank, size, comm, argc = 0;
	char** argv;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	comm = MPI_Comm_c2f(MPI_COMM_WORLD);

	// Validating Given solution: computing residual
	if (rank == 0)
	{
		struct matrix_descr descr_A;
		sparse_matrix_t csr_A;
		sparse_status_t csr_status;
		csr_status = mkl_sparse_d_create_csr(&csr_A, SPARSE_INDEX_BASE_ZERO, n, n, ia, (ia + 1), ja, a);

		if (csr_status != SPARSE_STATUS_SUCCESS)
		{
			std::cout << " <logs> Error occured while tried to create sparse struct." << std::endl;
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Abort(MPI_COMM_WORLD, -1);
			delete[] a; delete[] b; delete[] x; delete ia; delete[] ja; delete[] iparm;
			mkl_sparse_destroy(csr_A);
			return -1;
		}

		descr_A.type = SPARSE_MATRIX_TYPE_GENERAL;
		csr_status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csr_A, descr_A, x, 0.0, b_test);

		if (csr_status != SPARSE_STATUS_SUCCESS)
		{
			std::cout << " <logs> Error occured while tried to create sparse struct." << std::endl;
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Abort(MPI_COMM_WORLD, -1);
			delete[] a; delete[] b; delete[] x; delete ia; delete[] ja; delete[] iparm;
			mkl_sparse_destroy(csr_A);
			return -1;
		}

		mkl_sparse_destroy(csr_A);
		double resid = 0.0, resid0 = 0.0;
		for (int i = 0; i < n; ++i)
		{
			resid += (b_test[i] - b[i]) * (b_test[i] - b[i]);
			resid0 += b[i] * b[i];
		}

		resid = sqrt(resid) / sqrt(resid0);
		std::cout << " <logs> Given system residual: " << resid << std::endl;
	}

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

	// Calculating mean deviation
	
	if (!rank)
	{
		double mean_dev = 0.0;
		for (int i = 0; i < X_size; ++i) mean_dev += std::abs(X_true[i] - x[i]);
		mean_dev /= X_size;
		std::cout << " <logs> Mean deviation from Fidesys vector: " << mean_dev << std::endl;
	}

	if (rank == 0)
	{
		struct matrix_descr descr_A;
		sparse_matrix_t csr_A;
		sparse_status_t csr_status;
		csr_status = mkl_sparse_d_create_csr(&csr_A, SPARSE_INDEX_BASE_ZERO, n, n, ia, ia + 1, ja, a);

		if (csr_status != SPARSE_STATUS_SUCCESS)
		{
			std::cout << " <logs> Error occured while tried to create sparse struct." << std::endl;
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Abort(MPI_COMM_WORLD, -1);
			delete[] a; delete[] b; delete[] x; delete ia; delete[] ja; delete[] iparm;
			mkl_sparse_destroy(csr_A);
			return -1;
		}

		descr_A.type = SPARSE_MATRIX_TYPE_GENERAL;
		csr_status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csr_A, descr_A, x, 0.0, b_test);

		if (csr_status != SPARSE_STATUS_SUCCESS)
		{
			std::cout << " <logs> Error occured while tried to create sparse struct." << std::endl;
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Abort(MPI_COMM_WORLD, -1);
			delete[] a; delete[] b; delete[] x; delete ia; delete[] ja; delete[] iparm;
			mkl_sparse_destroy(csr_A);
			return -1;
		}

		mkl_sparse_destroy(csr_A);
		double resid = 0.0, resid0 = 0.0;
		for (int i = 0; i < n; ++i)
		{
			resid += (b_test[i] - b[i]) * (b_test[i] - b[i]);
			resid0 += b[i] * b[i];
		}

		resid = sqrt(resid) / sqrt(resid0);
		std::cout << " <logs> Solution residual: " << resid << std::endl;
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
