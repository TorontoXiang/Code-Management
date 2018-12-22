#include "TResistance_net.h"
#include <algorithm>
#include "mkl_pardiso.h"
#include "mkl_types.h"
using namespace std;
double Snode::Calculate_diagonal_element()
{
	double r = 0;
	for (int i = 0; i < _connected_node.size(); i++)
	{
		r = r + 1.0 / _connected_resistance[i];
	}
	return r;
}
double Snode::Calculate_non_diagonal_element(int j)
{
	double r = 0;
	for (int i = 0; i < _connected_node.size(); i++)
	{
		if (_connected_node[i]==j)
		{
			r = r - 1.0 / _connected_resistance[i];
		}
	}
	return r;
}
void TResisitance_net::input_resistance_net(ifstream &input)
{
	int num_resistance;
	input >> num_resistance;
	int n1, n2, n_max;
	double R;
	for (int i = 0; i < num_resistance; i++)
	{
		input >> n1 >> n2 >> R;
		n_max = _Max_value(n1, n2);
		if (n_max>=_node_list.size())
		{
			_node_list.resize(n_max+1);
		}
		_node_list[n1]._connected_node.push_back(n2);
		_node_list[n1]._connected_resistance.push_back(R);
		_node_list[n2]._connected_node.push_back(n1);
		_node_list[n2]._connected_resistance.push_back(R);
	}
	input >> _nv1 >> _nv2 >> _v;
	_num_freedom = _num_node=_node_list.size();
	return;
}
void TResisitance_net::Generate_Sparse_Matrix()
{
	_iK = new int[_num_freedom+1];
	_iK[0] = 1;
	vector<int> connected_node;
	vector<int> jK;
	vector<double> K;
	for (int i = 1; i < _num_node; i++)
	{
		connected_node = _node_list[i]._connected_node;
		sort(connected_node.begin(), connected_node.end());
		connected_node.erase(unique(connected_node.begin(), connected_node.end()),connected_node.end());
		int num_non_zero = 1;
		jK.push_back(i);
		K.push_back(_node_list[i].Calculate_diagonal_element());
		for (int j = 0; j < connected_node.size(); j++)
		{
			if (connected_node[j]>=i)
			{
				num_non_zero = num_non_zero + 1;
				jK.push_back(connected_node[j]);
				K.push_back(_node_list[i].Calculate_non_diagonal_element(connected_node[j]));
			}
		}
		if (_nv1==i)
		{
			num_non_zero = num_non_zero + 1;
			jK.push_back(_num_freedom);
			K.push_back(1.0);
		}
		if (_nv2==i)
		{
			num_non_zero = num_non_zero + 1;
			jK.push_back(_num_freedom);
			K.push_back(-1.0);
		}
		_iK[i] = _iK[i - 1] + num_non_zero;
	}
	_iK[_num_freedom] = _iK[_num_freedom-1];
	int non_zero = jK.size();
	_jK = new int[non_zero];
	_K = new double[non_zero];
	for (int i = 0; i < non_zero; i++)
	{
		_jK[i] = jK[i];
		_K[i] = K[i];
	}
	return;
}
double TResisitance_net::Calculate_effective_resistance()
{
	double* b;
	double* x;
	b = new double[_num_freedom];
	x = new double[_num_freedom];
	for (int i = 0; i < _num_freedom; i++)
	{
		b[i] = x[i] = 0;
	}
	b[_num_freedom - 1] = _v;
	MKL_INT n, mtype, nrhs, iparm[64], maxfct, mnum, phase, error, msglvl, idum;
	n = _num_freedom;
	double ddum;
	void *pt[64];
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
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, _K, _iK, _jK, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during symbolic factorization: %d", error);
		system("Pause");
		exit(1);
	}
	//printf("\nReordering completed ... ");
	//printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	//printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
	//--------------------------------------------------------------------
	// Numerical factorization
	//--------------------------------------------------------------------
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, _K, _iK, _jK, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during numerical factorization: %d", error);
		system("Pause");
		exit(2);
	}
	phase = 33;
	iparm[7] = 2;			// Max numbers of iterative refinement steps.
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, _K, _iK, _jK, &idum, &nrhs, iparm, &msglvl, b, x, &error);
	if (error != 0)
	{
		printf("\nERROR during solution: %d", error);
		system("Pause");
		exit(3);
	}
	phase = -1;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, _iK, _jK, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	return _v / x[_num_freedom - 1];
}