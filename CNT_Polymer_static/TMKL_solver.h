#ifndef TMKL_SOLVER
#define TMKL_SOLVER
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include <iostream>
using namespace std;
class TMKL_solver
{
public:
	TMKL_solver();
	void Numerical_factorization();
	//Factorize the matrix and store the results in pt 
	void Solving(double* b, double* dis);
	//Solving the equation with a RHS b
	void Release_Memory();
	//Release the Memory
	void S_freedom_degree(int num_freedom) { n = num_freedom; };
	//Set the freedom of degree of the equation
	void S_matrix(int* ik, int* jk, double*k) { ia = ik, ja = jk, a = k; };
private:
	MKL_INT n, mtype, nrhs, iparm[64], maxfct, mnum, phase, error, msglvl, idum;
	double ddum; 
	void *pt[64];
	//The sparse matrix
	int* ia;
	int* ja;
	double* a;
};
#endif
