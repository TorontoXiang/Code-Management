#include "Tgrid_smooth_continum.h"
#include "TMATIsotropic_elastic.h"
#include "Sparse_matrix.h"
//#include "mkl_pardiso.h"
//#include "mkl_types.h"
Tgrid_smooth_continum::Tgrid_smooth_continum(Tbody_ALE* the_body):Tgrid_smooth_base(the_body)
{
	_smooth_name="Continum_analogy";
	_nume=_ALE_body->_nume;_nump=_ALE_body->_nump;
	//Create the pseudo cell list
	create_pseudo_cell_list();
	//Calculate the ID vector
	calculate_ID();
	//Calculate the structure of the sparse matrix
	calculate_iK_jK();
};
void Tgrid_smooth_continum::generate_new_grid()
{
	double* dis;
	dis=new double[_num_freedom_degree];
	set_initial_and_final_position();
	int num_move=10;
	for (int i = 0; i < num_move; i++)
	{
		calculate_current_position_at_conode(i,num_move);
		assemble_equations();
		solve_equations(dis);
		for (int i = 0; i < _nump; i++)
		{
			vec3D pos0=_ALE_body->_grid_old->_node_list[i]->_position0;
			vec3D pos_old=_ALE_body->_grid_old->_node_list[i]->_position;
			for (int j = 0; j < 3; j++)
			{
				if (_ID[3*i+j]>0)
				{
					_ALE_body->_grid_new->_node_list[i]->_position.value(j+1,pos0.access(j+1)+dis[_ID[3*i+j]-1]);
				}
				else if (_ID[3*i+j]==-2)
				{
					_ALE_body->_grid_new->_node_list[i]->_position.value(j+1,pos_old.access(j+1));
				}
				else
				{
					_ALE_body->_grid_new->_node_list[i]->_position.value(j+1,pos0.access(j+1));
				}
			}
		}
		solidify_node();
		//_ALE_body->plot_grid_new();
		//system("Pause");
	}
	return;	
}
void Tgrid_smooth_continum::assemble_equations()
{
	double Ke[24][24],load_e[24];
	int LM[24];
	//Initilize global stiffness matrix
	for (int i = 0; i < _num_non_zero; i++)
	{
		_K[i]=0;
	}
	//Initialize the load vector
	for (int i = 0; i < _num_freedom_degree; i++)
	{
		_load[i]=0;
	}
	//Assemble global matrix
	for (int n = 0; n < _nume; n++)
	{
		//Calculate the LM vector
		for (int i = 0; i < 8; i++)
		{
			int node_id=_pseudo_cell_list[n]._node_ptr[i]->_id-1;
			for (int j = 0; j < 3; j++)
			{
				LM[3*i+j]=_ID[3*node_id+j];
			}
		}
		//Calculate the element stiffness matrix
		_pseudo_cell_list[n].modify_Youngs();
		_pseudo_cell_list[n].calculate_element_stiffness(Ke);
		//Assemble the global stiffness
		for (int i = 0; i < 24; i++)
		{
			int freedom_i=LM[i];
			if (freedom_i>0)
			{
				//freedom_i is freedom
				for (int j = 0; j < 24; j++)
				{
					int freedom_j=LM[j];		
					if (freedom_j>0)
					{
						//freedom_j is freedom
						if (freedom_i<=freedom_j)
						{
							//The element locates in the upper triangle
							assmeble_stifness(freedom_i,freedom_j,_iK,_jK,_K,Ke[i][j]);
						}
					}
				}
			}
			else if (freedom_i==-2)
			{
				//The displacement of freedom_i is given
				int node_id=_pseudo_cell_list[n]._node_ptr[i/3]->_id-1;
				int direction=i%3;
				double dis=_ALE_body->_grid_old->_node_list[node_id]->calculate_displacement().access(direction+1);
				for (int j = 0; j < 24; j++)
				{
					load_e[j]=dis*Ke[j][i];
				}
				//Assemble element load to global load
				for (int j = 0; j < 24; j++)
				{
					int freedom_j=LM[j];
					if (freedom_j>0)
					{
						_load[freedom_j-1]=_load[freedom_j-1]-load_e[j];
					}
				}
			}
		}
	}
	return;
}
void Tgrid_smooth_continum::create_pseudo_cell_list()
{
	_pseudo_cell_list.resize(_nume);
	int nGauss=2;
	Tnode* node_ptr[8];
	TMAT_base** mat;
	mat=new TMAT_base*[nGauss*nGauss*nGauss];
	for (int i = 0; i < _nume; i++)
	{
		for (int j = 0; j < nGauss*nGauss*nGauss; j++)
		{
			TMATIsotropic_elastic* pseudo_material=new TMATIsotropic_elastic(1,100,0.1);
			mat[j]=pseudo_material;
		}
		int id=_ALE_body->_grid_old->_cell_list[i]->_id;
		for (int j = 0; j < 8; j++)
		{
			node_ptr[j]=_ALE_body->_grid_old->_cell_list[i]->_node_ptr[j];
		}
		Tcell_8H new_cell(id,nGauss,node_ptr,mat);
		_pseudo_cell_list[i]=new_cell;
	}
	return;
}
void Tgrid_smooth_continum::calculate_ID()
{
	_num_freedom_degree=_num_fixed=0;
	_ID.resize(3*_nump);
	for (int i = 0; i < _nump; i++)
	{
		if (_ALE_body->_grid_old->_node_list[i]->_is_conode==1)
		{
			_ALE_body->_grid_new->_node_list[i]->_is_conode=1;
			//Set the co-node flag for the other grid
			for (int j = 0; j < 3; j++)
			{
				if (_ALE_body->_grid_old->_node_list[i]->_bc_type_position[j]==1)
				{
					_ID[3*i+j]=-1;
				}
				else
				{
					_ID[3*i+j]=-2;
					_num_fixed=_num_fixed+1;
				}
			}
		}
		else
		{
			for (int j = 0; j < 3; j++)
			{
				if (_ALE_body->_grid_old->_node_list[i]->_bc_type_position[j]==1)
				{
					_ID[3*i+j]=-1;
				}
				else
				{
					_num_freedom_degree=_num_freedom_degree+1;
					_ID[3*i+j]=_num_freedom_degree;
				}
			}
		}
	}
	//_ID[6]=_ID[21]=_ID[39]=_ID[54]=-2;
	_load=new double[_num_freedom_degree];
	for (int i = 0; i < _num_freedom_degree; i++)
	{
		_load[i]=0;
	}
	return;
}
void Tgrid_smooth_continum::calculate_iK_jK()
{
	vector<vector<int>> IEN;
	IEN.resize(_nume);
	for (int i = 0; i < _nume; i++)
	{
		IEN[i].resize(8);
		for (int j = 0; j < 8; j++)
		{
			IEN[i][j]=_pseudo_cell_list[i]._node_ptr[j]->_id;
		}
	}
	calculate_matrix_structure(_ID,IEN,_num_freedom_degree,3,_iK,_jK,_num_non_zero);
	_K=new double[_num_non_zero];
	return;
}
void Tgrid_smooth_continum::set_initial_and_final_position()
{
	Tnode_fluid* node_ptr;
	for (int i = 0; i < _nump; i++)
	{
		if (_ALE_body->_grid_old->_node_list[i]->_is_conode==1)
		{
			node_ptr=_ALE_body->_grid_old->_node_list[i];
			node_ptr->_position_beginning=node_ptr->_position0;//Store the initial position in _position_beginning
			node_ptr->_position_temp=node_ptr->_position;      //Store the final position in _position_temp;
		}
	}
	return;
}
void Tgrid_smooth_continum::calculate_current_position_at_conode(int n,int num_move)
{
	Tnode_fluid* node_ptr;
	for (int i = 0; i < _nump; i++)
	{
		if (_ALE_body->_grid_old->_node_list[i]->_is_conode==1)
		{
			node_ptr=_ALE_body->_grid_old->_node_list[i];
			vec3D d_displacement;
			d_displacement=(node_ptr->_position_temp-node_ptr->_position_beginning)/num_move;
			node_ptr->_position=node_ptr->_position_beginning+d_displacement*(n+1);
		}
	}
	return;
}
void Tgrid_smooth_continum::solidify_node()
{
	Tnode_fluid* node_ptr_old;
	Tnode_fluid* node_ptr_new;
	for (int i = 0; i < _nump; i++)
	{
		node_ptr_old=_ALE_body->_grid_old->_node_list[i];
		node_ptr_new=_ALE_body->_grid_new->_node_list[i];
		node_ptr_old->_position0=node_ptr_new->_position;
	}
	return;
}
void Tgrid_smooth_continum::solve_equations(double* dis)
{
//	MKL_INT n=_num_freedom_degree;          //Number of cows of the matrix    
//	MKL_INT mtype=-2;                       //Real symmetric matrix  
//	MKL_INT nrhs=1;                         //Number of right hand sides.
//	void *pt[64];
//	MKL_INT iparm[64];
//	MKL_INT maxfct, mnum, phase, error, msglvl;
//	// Auxiliary variables.
//	MKL_INT i;
//	double ddum;			//Double dummy 
//	MKL_INT idum;			//Integer dummy
////-------------------------------------------------------------------- 
//// Setup Pardiso control parameters. */
////--------------------------------------------------------------------
//	for (i = 0; i < 64; i++)
//    {
//      iparm[i] = 0;
//    }
//  iparm[0] = 1;			/* No solver default */
//  iparm[1] = 2;			/* Fill-in reordering from METIS */
//  /* Numbers of processors, value of OMP_NUM_THREADS */
//  iparm[2] = 1;
//  iparm[3] = 0;			/* No iterative-direct algorithm */
//  iparm[4] = 0;			/* No user fill-in reducing permutation */
//  iparm[5] = 0;			/* Write solution into x */
//  iparm[6] = 0;			/* Not in use */
//  iparm[7] = 2;			/* Max numbers of iterative refinement steps */
//  iparm[8] = 0;			/* Not in use */
//  iparm[9] = 13;		/* Perturb the pivot elements with 1E-13 */
//  iparm[10] = 1;		/* Use nonsymmetric permutation and scaling MPS */
//  iparm[11] = 0;		/* Not in use */
//  iparm[12] = 0;		/* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
//  iparm[13] = 0;		/* Output: Number of perturbed pivots */
//  iparm[14] = 0;		/* Not in use */
//  iparm[15] = 0;		/* Not in use */
//  iparm[16] = 0;		/* Not in use */
//  iparm[17] = -1;		/* Output: Number of nonzeros in the factor LU */
//  iparm[18] = -1;		/* Output: Mflops for LU factorization */
//  iparm[19] = 0;		/* Output: Numbers of CG Iterations */
//  maxfct = 1;			/* Maximum number of numerical factorizations. */
//  mnum = 1;			/* Which factorization to use. */
//  msglvl = 1;			/* Print statistical information in file */
//  error = 0;			/* Initialize error flag */
////-------------------------------------------------------------------- 
//// Initialize the internal solver memory pointer. This is only 
//// necessary for the FIRST call of the PARDISO solver.
////--------------------------------------------------------------------
//	for (i = 0; i < 64; i++)
//	{
//		pt[i] = 0;
//	}
////--------------------------------------------------------------------
//// Reordering and Symbolic Factorization. This step also allocates
//// all memory that is necessary for the factorization
////--------------------------------------------------------------------
//	phase = 11;
//	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,&n, _K, _iK, _jK, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
//	if (error != 0)
//    {
//		printf ("\nERROR during symbolic factorization: %d", error);
//		system("Pause");
//		exit (1);
//    }
//	printf ("\nReordering completed ... ");
//	printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
//	printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
////--------------------------------------------------------------------
//// Numerical factorization
////--------------------------------------------------------------------
//	phase = 22;
//	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,&n, _K, _iK, _jK, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
//	if (error != 0)
//    {
//		printf ("\nERROR during numerical factorization: %d", error);
//		system("Pause");
//		exit (2);
//    }
//	printf ("\nFactorization completed ... ");
////--------------------------------------------------------------------
//// Back substitution and iterative refinement
////--------------------------------------------------------------------
//	phase = 33;
//	iparm[7] = 2;			// Max numbers of iterative refinement steps.
//	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,&n, _K, _iK, _jK, &idum, &nrhs, iparm, &msglvl, _load, dis, &error);
//	if (error != 0)
//    {
//		printf ("\nERROR during solution: %d", error);
//		system("Pause");
//		exit (3);
//	}
//	printf ("\nSolve completed ... ");
//	printf ("\nThe solution of the system is: ");
//	//for (i = 0; i < n; i++)
// //   {
//	//	printf ("\n x [%d] = % f", i, dis[i]);
// //   }
//	//printf ("\n");
////--------------------------------------------------------------------
//// Termination and release of memory.
//// -------------------------------------------------------------------
//  phase = -1;
//  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,&n, &ddum, _iK, _jK, &idum, &nrhs,iparm, &msglvl, &ddum, &ddum, &error);
//  return;
}