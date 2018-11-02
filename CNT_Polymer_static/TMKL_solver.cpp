#include "TMKL_solver.h"
TMKL_solver::TMKL_solver()
{
	for (int i = 0; i < 64; i++)
	{
		iparm[i] = 0;
	}
	iparm[0] = 1;			/* No solver default */
	iparm[1] = 2;			/* Fill-in reordering from METIS */
	iparm[2] = 1;           /* Numbers of processors, value of OMP_NUM_THREADS */
	iparm[3] = 0;			/* No iterative-direct algorithm */
	iparm[4] = 0;			/* No user fill-in reducing permutation */
	iparm[5] = 0;			/* Write solution into x */
	iparm[6] = 0;			/* Not in use */
	iparm[7] = 2;			/* Max numbers of iterative refinement steps */
	iparm[8] = 0;			/* Not in use */
	iparm[9] = 13;		    /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1;		    /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0;		    /* Not in use */
	iparm[12] = 0;		    /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm[13] = 0;		    /* Output: Number of perturbed pivots */
	iparm[14] = 0;		    /* Not in use */
	iparm[15] = 0;		    /* Not in use */
	iparm[16] = 0;		    /* Not in use */
	iparm[17] = -1;		    /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;		    /* Output: Mflops for LU factorization */
	iparm[19] = 0;		    /* Output: Numbers of CG Iterations */
	maxfct = 1;			    /* Maximum number of numerical factorizations. */
	mnum = 1;			    /* Which factorization to use. */
	msglvl = 0;			    /* Print statistical information in file */
	mtype = -2;             /*Real symmetric matrix */  
	nrhs = 1;               /*Number of right hand sides */
	error = 0;			    /* Initialize error flag */
}
void TMKL_solver::Numerical_factorization()
{
	//-------------------------------------------------------------------- 
	 // Initialize the internal solver memory pointer. This is only 
	 // necessary for the FIRST call of the PARDISO solver.
	 //--------------------------------------------------------------------
	for (int i = 0; i < 64; i++)
	{
		pt[i] = 0;
	}
	//--------------------------------------------------------------------
	// Reordering and Symbolic Factorization. This step also allocates
	// all memory that is necessary for the factorization
	//--------------------------------------------------------------------
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during symbolic factorization: %d", error);
		system("Pause");
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
	//--------------------------------------------------------------------
	// Numerical factorization
	//--------------------------------------------------------------------
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during numerical factorization: %d", error);
		system("Pause");
		exit(2);
	}
	printf("\nFactorization completed ... ");
	return;
}
void TMKL_solver::Solving(double* b, double* dis)
{
	phase = 33;
	iparm[7] = 2;			// Max numbers of iterative refinement steps.
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, dis, &error);
	if (error != 0)
	{
		printf("\nERROR during solution: %d", error);
		system("Pause");
		exit(3);
	}
	//printf("\nSolve completed ... ");
}
void TMKL_solver::Release_Memory()
{
	phase = -1;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	return;
}