#include "TFE_grid.h"
#include <fstream>
#include "Sparse_matrix.h"
#include "Shape_function.h"
#include "mkl.h"
#include "mkl_blas.h"
#include <cmath>
using namespace std;
void Tgrid::calculate_ID()
{
	_num_freedom_degree = _num_fixed = 0;
	_ID.resize(3 * _nump);
	for (int i = 0; i < _nump; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (_node_list[i]._bc[j] == 1)
			{
				_num_fixed = _num_fixed + 1;
				_ID[3 * i + j] = -_num_fixed;

			}
			else
			{
				_num_freedom_degree = _num_freedom_degree + 1;
				_ID[3 * i + j] = _num_freedom_degree;
			}
		}
	}
	_p = new double[_num_freedom_degree];
	_dis = new double[_num_freedom_degree];
	for (int i = 0; i < _num_freedom_degree; i++)
	{
		_p[i] = _dis[i] = 0;
	}
	return;
}
void Tgrid_Polymer::assemble_equations()
{
	double Ke[24][24];
	int LM[24];
	//Assemble global matrix
	for (int n = 0; n < _nume; n++)
	{
		//Calculate the LM vector
		for (int i = 0; i < 8; i++)
		{
			int node_id = _cell_list[n]._node_ptr[i]->_id - 1;
			for (int j = 0; j < 3; j++)
			{
				LM[3 * i + j] = _ID[3 * node_id + j];
			}
		}
		//Calculate the element stiffness matrix
		Smat mat = _mat_list[_cell_list[n]._mid];
		_cell_list[n].calculate_element_stiffness(Ke,mat);
		//Assemble the global stiffness
		for (int i = 0; i < 24; i++)
		{
			int freedom_i = LM[i];
			if (freedom_i > 0)
			{
				//freedom_i is freedom
				for (int j = 0; j < 24; j++)
				{
					int freedom_j = LM[j];
					if (freedom_j > 0)
					{
						//freedom_j is freedom
						if (freedom_i <= freedom_j)
						{
							//The element locates in the upper triangle
							assmeble_stifness(freedom_i, freedom_j, _iK, _jK, _K, Ke[i][j]);
						}
					}
				}
			}
		}
	}
	return;
}
void Tgrid_CNT::assemble_equations()
{
	double Ke[24][24];
	int LM[24];
	int dis_ID[24];
	//Assemble global matrix
	for (int n = 0; n < _nume; n++)
	{
		//Calculate the LM vector
		for (int i = 0; i < 8; i++)
		{
			int node_id = _cell_list[n]._node_ptr[i]->_id - 1;
			for (int j = 0; j < 3; j++)
			{
				LM[3 * i + j] = _ID[3 * node_id + j];
				dis_ID[3 * i + j] = 3 * node_id + j+1;
			}
		}
		//Calculate the element stiffness matrix
		Smat mat = _mat_list[_cell_list[n]._mid];
		_cell_list[n].calculate_element_stiffness(Ke, mat);
		//Assemble the global stiffness
		for (int i = 0; i < 24; i++)
		{
			int freedom_i = LM[i];
			if (freedom_i > 0)
			{
				//freedom_i is freedom
				for (int j = 0; j < 24; j++)
				{
					int freedom_j = LM[j];
					if (freedom_j > 0)
					{
						//freedom_j is freedom
						if (freedom_i <= freedom_j)
						{
							//The element locates in the upper triangle
							assmeble_stifness(freedom_i, freedom_j, _iK, _jK, _K, Ke[i][j]);
						}
					}
				}
			}
			else if (freedom_i<0)
			{
				for (int j = 0; j < 24; j++)
				{
					int dis_ID_j = dis_ID[j];
					//Assemble this element to _KB
					assmeble_stifness(-freedom_i, dis_ID_j, _iKB, _jKB, _KB, Ke[i][j]);
				}
			}
		}
	}
	return;
}
void Tgrid_Polymer::calculate_iK_jK()
{
	vector<vector<int>> IEN;
	IEN.resize(_nume);
	for (int i = 0; i < _nume; i++)
	{
		IEN[i].resize(8);
		for (int j = 0; j < 8; j++)
		{
			IEN[i][j] = _cell_list[i]._node_ptr[j]->_id;
		}
	}
	int num_non_zero;
	calculate_matrix_structure_KII(_ID, IEN, _num_freedom_degree, 3, _iK, _jK, num_non_zero);
	_K = new double[num_non_zero];
	for (int i = 0; i < num_non_zero; i++)
	{
		_K[i] = 0;
	}
	//Allocate memory for reacting froce
	_reacting_force = new double[_num_fixed];
	for (int i = 0; i < _num_fixed; i++)
	{
		_reacting_force[i] = 0;
	}
	return;
}
void Tgrid_CNT::calculate_iK_jK()
{
	vector<vector<int>> IEN;
	IEN.resize(_nume);
	for (int i = 0; i < _nume; i++)
	{
		IEN[i].resize(8);
		for (int j = 0; j < 8; j++)
		{
			IEN[i][j] = _cell_list[i]._node_ptr[j]->_id;
		}
	}
	//Calculate the sparse matrix structure for solving the equilibrium equation (KII)
	int num_non_zero;
	calculate_matrix_structure_KII(_ID, IEN, _num_freedom_degree, 3, _iK, _jK, num_non_zero);
	_K = new double[num_non_zero];
	for (int i = 0; i < num_non_zero; i++)
	{
		_K[i] = 0;
	}
	//Calculate the sparse matrix structure for calculating the racting force (KB)
	calculate_matrix_structure_KB(_ID, IEN, 3 * _nump, 3, _iKB, _jKB, num_non_zero);
	_KB = new double[num_non_zero];
	for (int i = 0; i < num_non_zero; i++)
	{
		_KB[i] = 0;
	}
	//Allocate memory for reacting froce
	_reacting_force = new double[3*_nump];
	_dis_whole = new double[3 * _nump];
	for (int i = 0; i < 3 * _nump; i++)
	{
		_reacting_force[i] = _dis_whole[i] = 0;
	}
	return;
}
void Tgrid_CNT::Create_MKL_solver()
{
	calculate_ID();
	calculate_iK_jK();
	assemble_equations();
	MKL_solver.S_freedom_degree(_num_freedom_degree);
	MKL_solver.S_matrix(_iK, _jK, _K);
	MKL_solver.Numerical_factorization();
	return;
}
void Tgrid_CNT::Solving_equilibrium_equation()
{
	MKL_solver.Solving(_p, _dis);
}
Sstress** Tgrid::Calculate_cell_stress()
{
	Sstress** cell_stress;
	cell_stress = new Sstress*[_nume];
	for (int i = 0; i < _nume; i++)
	{
		cell_stress[i] = new Sstress[8];
	}
	double d_cell[8][3];
	int mid;
	for (int i = 0; i < _nume; i++)
	{
		calculate_cell_displacement(i, d_cell);
		mid = _cell_list[i]._mid;
		Smat mat = _mat_list[mid];
		_cell_list[i].calculate_cell_stress(cell_stress[i], d_cell, mat);
	}
	return cell_stress;
}
void Tgrid::calculate_cell_displacement(int cell_id,double(&d_cell)[8][3])
{
	int node_id;
	for (int i = 0; i < 8; i++)
	{
		node_id = _cell_list[cell_id]._node_ptr[i]->_id-1;
		for (int j = 0; j < 3; j++)
		{
			if (_ID[3 * node_id + j] > 0)
			{
				d_cell[i][j] = _dis[_ID[3 * node_id + j] - 1];
			}
			else if (_ID[3 * node_id + j] < 0)
			{
				d_cell[i][j] = _node_list[node_id]._bc_value[j];
			}
		}
	}
	return;
}
void Tgrid::Output_Tecplot(string file_name,int mid)
{
	//Calculate cell stress
	Sstress** cell_stress;
	cell_stress = Calculate_cell_stress();
	Snode_stress* node_stress;
	node_stress = new Snode_stress[_nump];
	int node_id;
	for (int i = 0; i < _nume; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			if (mid == -1 || mid == _cell_list[i]._mid)
			{
				node_id = _cell_list[i]._node_ptr[j]->_id - 1;
				node_stress[node_id].num = node_stress[node_id].num + 1;
				node_stress[node_id].add(cell_stress[i][j]);
			}
		}
	}
	for (int i = 0; i < _nump; i++)
	{
		node_stress[i].Average();
	}
	for (int i = 0; i < _nume; i++)
	{
		delete cell_stress[i];
	}
	delete cell_stress;
	//Update the grid position
	double** pos;
	pos = new double*[_nump];
	for (int i = 0; i < _nump; i++)
	{
		pos[i] = new double[3];
	}
	for (int i = 0; i < _nump; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (_ID[3 * i + j] > 0)
			{
				pos[i][j] = _node_list[i]._pos[j] + _dis[_ID[3 * i + j] - 1];
			}
			else if (_ID[3 * i + j] < 0)
			{
				pos[i][j] = _node_list[i]._pos[j] + _node_list[i]._bc_value[j];
			}
		}
	}
	//Output the Tecplot result
	int nume_plot = 0;
	if (mid==-1)
	{
		nume_plot = _nume;
	}
	else
	{
		for (int i = 0; i < _nume; i++)
		{
			if (_cell_list[i]._mid==mid)
			{
				nume_plot = nume_plot + 1;
			}
		}
	}
	ofstream output;
	output.open(file_name);
	output << "TITLE = \"Tecplot Grid\"" << endl;
	output << "VARIABLES = \"X\",\"Y\",\"Z\",\"SXX\",\"SYY\",\"SZZ\",\"SXY\",\"SXZ\",\"SYZ\",\"EQS\"" << endl;
	output << "ZONE F=FEPOINT,N=" << _nump << "," << "E=" << nume_plot << "," << "ET=BRICK" << endl;
	for (int i = 0; i < _nump; i++)
	{
		output << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << " ";
		output << node_stress[i].sxx << " " << node_stress[i].syy << " " << node_stress[i].szz << " ";
		output << node_stress[i].sxy << " " << node_stress[i].sxz << " " << node_stress[i].syz << " ";
		output << node_stress[i].equivalent_stress() << endl;
	}
	for (int i = 0; i < _nume; i++)
	{
		if (mid==-1 || mid==_cell_list[i]._mid)
		{
			for (int j = 0; j < 8; j++)
			{
				output << _cell_list[i]._node_ptr[j]->_id << " ";
			}
		}
		output << endl;
	}
	delete node_stress;
	for (int i = 0; i < _nump; i++)
	{
		delete pos[i];
	}
	delete pos;
	return;
}
void Tgrid_Polymer::Input_Polymer()
{
	//Read the input information
	ifstream input;
	input.open("RegularGridTest.k");
	Skeyword keyword;
	read_in_keyword_file(input, keyword);

	//Generate the regular polymer grid
	_nx = keyword.regular_grid.nx; _ny = keyword.regular_grid.ny; _nz = keyword.regular_grid.nz;
	_nump = (_nx + 1)*(_ny + 1)*(_nz + 1);
	_nume = _nx * _ny*_nz;
	_mat_list.resize(1);
	for (int i = 0; i < 3; i++)
	{
		_x_min[i] = keyword.regular_grid.x_min[i];
		_x_max[i] = keyword.regular_grid.x_max[i];
	}
	_interval[0] = (_x_max[0] - _x_min[0]) / _nx;
	_interval[1] = (_x_max[1] - _x_min[1]) / _ny;
	_interval[2] = (_x_max[2] - _x_min[2]) / _nz;
	_node_list.resize(_nump);
	_cell_list.resize(_nume);
	double dx[3], coor[3];
	int index;
	for (int i = 0; i < _nx + 1; i++)
	{
		for (int j = 0; j < _ny + 1; j++)
		{
			for (int k = 0; k < _nz + 1; k++)
			{
				index = (_ny + 1)*(_nz + 1)*i + (_nz + 1)*j + k;
				dx[0] = i * _interval[0]; dx[1] = j * _interval[1]; dx[2] = k * _interval[2];
				coor[0] = _x_min[0] + dx[0]; coor[1] = _x_min[1] + dx[1]; coor[2] = _x_min[2] + dx[2];
				Snode new_node(index + 1, coor);
				_node_list[index] = new_node;
			}
		}
	}
	int index_node[8];
	Snode* node_ptr[8];
	for (int i = 0; i < _nx; i++)
	{
		for (int j = 0; j < _ny; j++)
		{
			for (int k = 0; k < _nz; k++)
			{
				index = _ny * _nz*i + _nz * j + k; //The index convertion for cell
				index_node[0] = (_ny + 1)*(_nz + 1)*i + (_nz + 1)*j + k; index_node[1] = (_ny + 1)*(_nz + 1)*(i + 1) + (_nz + 1)*j + k;
				index_node[2] = (_ny + 1)*(_nz + 1)*(i + 1) + (_nz + 1)*(j + 1) + k; index_node[3] = (_ny + 1)*(_nz + 1)*i + (_nz + 1)*(j + 1) + k;
				index_node[4] = (_ny + 1)*(_nz + 1)*i + (_nz + 1)*j + k + 1; index_node[5] = (_ny + 1)*(_nz + 1)*(i + 1) + (_nz + 1)*j + k + 1;
				index_node[6] = (_ny + 1)*(_nz + 1)*(i + 1) + (_nz + 1)*(j + 1) + k + 1; index_node[7] = (_ny + 1)*(_nz + 1)*i + (_nz + 1)*(j + 1) + k + 1;
				for (int n = 0; n < 8; n++)
				{
					node_ptr[n] = &_node_list[index_node[n]];
				}
				Tcell new_cell(index, node_ptr,0);
				_cell_list[index] = new_cell;
			}
		}
	}
	int num_bc = keyword.regular_bc_list.size();
	for (int i = 0; i < num_bc; i++)
	{
		Sregular_grid_bc bc = keyword.regular_bc_list[i];
		apply_regular_bc(bc);
	}
	_mat_list[0]._E = keyword.material_list[0].Youngs;
	_mat_list[0]._mu = keyword.material_list[0].Possion;
	detect_boundary_cells();
	return;
}
void Tgrid_Polymer::apply_regular_bc(Sregular_grid_bc& bc)
{
	int nx_min = 0, ny_min = 0, nz_min = 0;
	int nx_max = _nx + 1, ny_max = _ny + 1, nz_max = _nz + 1;
	if (bc.FaceName == "x_min")
	{
		nx_min = 0; nx_max = 1;
	}
	else if (bc.FaceName == "x_max")
	{
		nx_min = _nx; nx_max = _nx + 1;
	}
	else if (bc.FaceName == "y_min")
	{
		ny_min = 0; ny_max = 1;
	}
	else if (bc.FaceName == "y_max")
	{
		ny_min = _ny; ny_max = _ny + 1;
	}
	else if (bc.FaceName == "z_min")
	{
		nz_min = 0; nz_max = 1;
	}
	else if (bc.FaceName == "z_max")
	{
		nz_min = _nz; nz_max = _nz + 1;
	}
	else
	{
		cout << "Invalid boundary contition for regular grid" << endl;
		system("Pause");
		exit(0);
	}
	for (int i = nx_min; i < nx_max; i++)
	{
		for (int j = ny_min; j < ny_max; j++)
		{
			for (int k = nz_min; k < nz_max; k++)
			{
				int index = (_ny + 1)*(_nz + 1)*i + (_nz + 1)*j + k;
				Snode* nodeptr = &_node_list[index];
				nodeptr->_bc[bc.direction] = 1;
				nodeptr->_bc_value[bc.direction] = bc.value;
			}
		}
	}
	return;
}
void Tgrid_CNT::Input_CNT(Skeyword& keyword)
{
	_nume = keyword.cell_8_list.size();
	_nump = keyword.node_list.size();
	_node_list.resize(_nump);
	_cell_list.resize(_nume);
	for (int i = 0; i < _nump; i++)
	{
		double coor[3];
		coor[0] = keyword.node_list[i].x; coor[1] = keyword.node_list[i].y; coor[2] = keyword.node_list[i].z;
		int id = keyword.node_list[i].id;
		Snode new_node(id, coor);
		_node_list[id-1] = new_node;
	}
	Snode* node_ptr[8];
	int mid, IENj, cell_id;
	for (int i = 0; i < _nume; i++)
	{
		mid = keyword.cell_8_list[i].mid-1;
		cell_id = keyword.cell_8_list[i].cell_id;
		for (int j = 0; j < 8; j++)
		{
			IENj = keyword.cell_8_list[i].IEN[j];
			node_ptr[j] = &_node_list[IENj - 1];
		}
		Tcell new_cell(cell_id, node_ptr, mid);
		_cell_list[cell_id-1] = new_cell;
	}
	int num_material = keyword.material_list.size();
	_mat_list.resize(num_material);
	for (int i = 0; i < num_material; i++)
	{
		_mat_list[i]._E = keyword.material_list[i].Youngs;
		_mat_list[i]._mu = keyword.material_list[i].Possion;
	}
	detect_boundary_nodes();
	detect_boundary_cells();
	return;
}
void Tgrid_CNT::detect_boundary_nodes()
{
	int* node_flag;
	int* cell_flag;
	node_flag = new int[_nump];
	cell_flag = new int[_nume];
	for (int i = 0; i < _nump; i++)
	{
		node_flag[i] = 0;
	}
	for (int i = 0; i < _nume; i++)
	{
		cell_flag[i] = 0;
	}
	//Get cells connected to every node
	vector<vector<int>> node2cell;
	node2cell.resize(_nump);
	for (int i = 0; i < _nume; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			int node_id = _cell_list[i]._node_ptr[j]->_id - 1;
			node2cell[node_id].push_back(i);
		}
	}
	//Calculate boundary nodes
	int cell_face[6][4] = { {0,3,2,1},{0,1,5,4},{1,2,6,5},{2,3,7,6},{0,4,7,3},{4,5,6,7} };
	vector<int> cell2cell;
	for (int i = 0; i < _nume; i++)
	{
		//Calculate the cells connected to cell i
		cell_flag[i] = 1;
		for (int j = 0; j < 8; j++)
		{
			int node_id = _cell_list[i]._node_ptr[j]->_id - 1;
			int num_node2cell = node2cell[node_id].size();
			for (int k = 0; k < num_node2cell; k++)
			{
				int cell_id = node2cell[node_id][k];
				if (cell_flag[cell_id] == 0)
				{
					cell2cell.push_back(cell_id);
					cell_flag[cell_id] = 1;
				}
			}
		}
		//Rest the cell_flag
		cell_flag[i] = 0;
		for (int j = 0; j < 8; j++)
		{
			int node_id = _cell_list[i]._node_ptr[j]->_id - 1;
			int num_node2cell = node2cell[node_id].size();
			for (int k = 0; k < num_node2cell; k++)
			{
				int cell_id = node2cell[node_id][k];
				cell_flag[cell_id] = 0;
			}
		}
		//Calculate the adjacent cells of cell i
		int adjacent[6];
		for (int j = 0; j < 6; j++)
		{
			adjacent[j] = -1;
		}
		for (int j = 0; j < 6; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				int local_node_id = cell_face[j][k];
				int node_id = _cell_list[i]._node_ptr[local_node_id]->_id - 1;
				node_flag[node_id] = 1;
			}
			int num_around = cell2cell.size();
			for (int k = 0; k < num_around; k++)
			{
				int cell_id = cell2cell[k];
				int temp = 0;
				for (int l = 0; l < 8; l++)
				{
					int node_id = _cell_list[cell_id]._node_ptr[l]->_id - 1;
					temp = temp + node_flag[node_id];
				}
				if (temp == 4)
				{
					adjacent[j] = cell_id;
					break;
				}
			}
			for (int k = 0; k < 4; k++)
			{
				int local_node_id = cell_face[j][k];
				int node_id = _cell_list[i]._node_ptr[local_node_id]->_id - 1;
				node_flag[node_id] = 0;
			}
		}
		for (int j = 0; j < 6; j++)
		{
			if (adjacent[j]==-1)
			{
				for (int k = 0; k < 4; k++)
				{
					int local_node_id = cell_face[j][k];
					int node_id = _cell_list[i]._node_ptr[local_node_id]->_id - 1;
					_node_list[node_id]._is_surface = true;
					_node_list[node_id]._bc[0] = _node_list[node_id]._bc[1] = _node_list[node_id]._bc[2] = 1;
				}
			}
		}
		cell2cell.resize(0);
	}
	delete cell_flag;
	delete node_flag;
	return;
}
void Tgrid::detect_boundary_cells()
{
	for (int i = 0; i < _nume; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			int node_id = _cell_list[i]._node_ptr[j]->_id - 1;
			int sum = _node_list[node_id]._bc[0] + _node_list[node_id]._bc[1] + _node_list[node_id]._bc[2];
			if (sum>0 || _node_list[node_id]._is_surface)
			{
				_cell_list[i]._is_boundary = true;
				break;
			}
		}
	}
	return;
}
void Tgrid_Polymer::calculate_reacting_force()
{
	for (int i = 0; i < _num_fixed; i++)
	{
		_reacting_force[i] = 0;
	}
	double d_cell[24];
	int LM[24];
	double Ke[24][24];
	for (int i = 0; i < _nume; i++)
	{
		if (_cell_list[i]._is_boundary)
		{
			for (int j = 0; j < 8; j++)
			{
				int node_id = _cell_list[i]._node_ptr[j]->_id - 1;
				for (int k = 0; k < 3; k++)
				{
					LM[3 * j + k] = _ID[3 * node_id + k];
					if (LM[3*j+k]>0)
					{
						d_cell[3 * j + k] = _dis[LM[3 * j + k] - 1];
					}
					else
					{
						d_cell[3 * j + k] = _node_list[node_id]._bc_value[k];
					}
				}
			}
			Smat mat = _mat_list[_cell_list[i]._mid];
			_cell_list[i].calculate_element_stiffness(Ke, mat);
			for (int j = 0; j < 24; j++)
			{
				int freedom_j = LM[j];
				if (freedom_j<0)
				{
					int fixed_id = -freedom_j-1;
					double re_force = 0;
					for (int k = 0; k < 24; k++)
					{
						re_force = re_force + Ke[j][k] * d_cell[k];
					}
					_reacting_force[fixed_id] = _reacting_force[fixed_id] + re_force;
				}
			}
		}
	}
	return;
}
void Tgrid_CNT::calculate_reacting_force()
{
	for (int i = 0; i < 3*_nump; i++)
	{
		_reacting_force[i] = 0;
	}
	for (int i = 0; i < _nump; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (_ID[3*i+j]>0)
			{
				_dis_whole[3 * i + j] = _dis[_ID[3 * i + j] - 1];
			}
			else
			{
				_dis_whole[3 * i + j] = _node_list[i]._bc_value[j];
			}
		}
	}
	char uplo = 'n';
	int m = _nump * 3;
	mkl_dcsrgemv(&uplo, &m, _KB, _iKB, _jKB, _dis_whole, _reacting_force);
	return;
}
vec3D Tgrid_Polymer::calculate_external_load(string face_name)
{
	int nx_min = 0, ny_min = 0, nz_min = 0;
	int nx_max = _nx + 1, ny_max = _ny + 1, nz_max = _nz + 1;
	if (face_name == "x_min")
	{
		nx_min = 0; nx_max = 1;
	}
	else if (face_name == "x_max")
	{
		nx_min = _nx; nx_max = _nx + 1;
	}
	else if (face_name == "y_min")
	{
		ny_min = 0; ny_max = 1;
	}
	else if (face_name == "y_max")
	{
		ny_min = _ny; ny_max = _ny + 1;
	}
	else if (face_name == "z_min")
	{
		nz_min = 0; nz_max = 1;
	}
	else if (face_name == "z_max")
	{
		nz_min = _nz; nz_max = _nz + 1;
	}
	else
	{
		cout << "Invalid boundary contition for regular grid" << endl;
		system("Pause");
		exit(0);
	}
	vec3D external_force;
	for (int i = nx_min; i < nx_max; i++)
	{
		for (int j = ny_min; j < ny_max; j++)
		{
			for (int k = nz_min; k < nz_max; k++)
			{
				int index = (_ny + 1)*(_nz + 1)*i + (_nz + 1)*j + k;
				vec3D nodal_re_froce;
				for (int l = 0; l < 3; l++)
				{
					int freedom_id = _ID[3 * index + l];
					if (freedom_id<0)
					{
						nodal_re_froce.value(l + 1, _reacting_force[-freedom_id-1]);
					}
				}
				external_force = external_force + nodal_re_froce;
			}
		}
	}
	return external_force;
}
void Tgrid_Polymer::calculate_load_from_constraint()
{
	int LM[24];
	double Ke[24][24], load_e[24];
	for (int i = 0; i < _num_freedom_degree; i++)
	{
		//_load_constraint[i] = 0;
		_p[i] = 0;
	}
	for (int n = 0; n < _nume; n++)
	{
		if (_cell_list[n]._is_boundary)
		{
			//Calculate the LM vector
			for (int i = 0; i < 8; i++)
			{
				int node_id = _cell_list[n]._node_ptr[i]->_id - 1;
				for (int j = 0; j < 3; j++)
				{
					LM[3 * i + j] = _ID[3 * node_id + j];
				}
			}
			//Calculate the element stiffness matrix
			Smat mat = _mat_list[_cell_list[n]._mid];
			_cell_list[n].calculate_element_stiffness(Ke, mat);
			for (int i = 0; i < 24; i++)
			{
				int freedom_i = LM[i];
				if (freedom_i < 0)
				{
					//The displacement of freedom_i is given
					int node_id = _cell_list[n]._node_ptr[i / 3]->_id - 1;
					int direction = i % 3;
					double dis = _node_list[node_id]._bc_value[direction];
					for (int j = 0; j < 24; j++)
					{
						load_e[j] = dis * Ke[j][i];
					}
					//Assemble element load to global load
					for (int j = 0; j < 24; j++)
					{
						int freedom_j = LM[j];
						if (freedom_j > 0)
						{
							_p[freedom_j - 1] = _p[freedom_j - 1] - load_e[j];
						}
					}
				}
			}
		}
	}
	return;
}
void Tgrid_CNT::calculate_load_from_constraint()
{
	for (int i = 0; i < _num_freedom_degree; i++)
	{
		_p[i] = 0;
	}
	for (int i = 0; i < _nump; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int fixed_id = _ID[3 * i + j];
			if (fixed_id <0)
			{
				_dis_whole[-fixed_id - 1] = _node_list[i]._bc_value[j];
			}
		}
	}
	char uplo = 't';
	int m = _nump * 3;
	//for (int i = 0; i < m; i++)
	//{
	//	if (isnan(_dis_whole[i]))
	//	{
	//		cout << "nan occur in _dis_whole" << endl;
	//		system("Pause");
	//	}
	//}
	//for (int i = 0; i < m; i++)
	//{
	//	if (isnan(_reacting_force[i]))
	//	{
	//		cout << "nan occur in _dis_whole" << endl;
	//		system("Pause");
	//	}
	//}
	mkl_dcsrgemv(&uplo, &m, _KB, _iKB, _jKB, _dis_whole, _reacting_force);
	for (int i = 0; i < _nump; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int freedom_id= _ID[3 * i + j];
			if (freedom_id>0)
			{
				_p[freedom_id - 1] = -_reacting_force[3 * i + j];
				//if (isnan(_p[freedom_id - 1]))
				//{
				//	cout << "nan occur in calculate_load_from_constraint()" << endl;
				//	ofstream output;
				//	output.open("_KB.k");
				//	int num_non_zero = _iKB[m] - 1;
				//	for (int i = 0; i < num_non_zero; i++)
				//	{
				//		output << _KB[i] << endl;
				//	}
				//	system("Pause");
				//}
			}
		}
	}
	return;
}
void Tgrid_CNT::calculate_CNT_location(Tgrid_Polymer* grid_polymer)
{
	double x_min[3], x_max[3], interval[3];
	int nx_max, ny_max, nz_max;
	for (int i = 0; i < 3; i++)
	{
		x_min[i] = grid_polymer->_x_min[i];
		x_max[i] = grid_polymer->_x_max[i];
		interval[i] = grid_polymer->_interval[i];
	}
	nx_max = grid_polymer->_nx; ny_max = grid_polymer->_ny; nz_max = grid_polymer->_nz;
	for (int i = 0; i < _nump; i++)
	{
		_node_list[i].calculate_location(x_min, x_max, interval, nx_max, ny_max, nz_max);
	}
	return;
}
void Tgrid_CNT::calculate_CNT_boundary_displacement(Tgrid_Polymer* grid_polymer)
{
	for (int i = 0; i < _nump; i++)
	{
		if (_node_list[i]._is_surface)
		{
			_node_list[i]._bc_value[0] = _node_list[i]._bc_value[1] = _node_list[i]._bc_value[2] = 0;
			int cell_id = _node_list[i]._location->cell_id;
			double xi = _node_list[i]._location->iso_coor[0];
			double eta = _node_list[i]._location->iso_coor[1];
			double zeta = _node_list[i]._location->iso_coor[2];
			double cell_d[8][3];
			grid_polymer->calculate_cell_displacement(cell_id, cell_d);
			for (int j = 0; j < 8; j++)
			{
				double Nj = N_brick_8(j, xi, eta, zeta);
				_node_list[i]._bc_value[0] = _node_list[i]._bc_value[0] + Nj *cell_d[j][0];
				_node_list[i]._bc_value[1] = _node_list[i]._bc_value[1] + Nj *cell_d[j][1];
				_node_list[i]._bc_value[2] = _node_list[i]._bc_value[2] + Nj *cell_d[j][2];
			}
		}
	}
}
void Tgrid_Polymer::assemble_Fp(Tgrid_CNT* grid_CNT)
{
	int nump_CNT = grid_CNT->_nump;
	for (int i = 0; i < nump_CNT; i++)
	{
		if (grid_CNT->_node_list[i]._is_surface)
		{
			int cell_id = grid_CNT->_node_list[i]._location->cell_id;
			double xi = grid_CNT->_node_list[i]._location->iso_coor[0];
			double eta = grid_CNT->_node_list[i]._location->iso_coor[1];
			double zeta = grid_CNT->_node_list[i]._location->iso_coor[2];
			double force_P2C[3];
			for (int j = 0; j < 3; j++)
			{
				int fixed_id = -grid_CNT-> _ID[3 * i + j] - 1;
				force_P2C[j] = grid_CNT->_reacting_force[fixed_id];
			}
			for (int j = 0; j < 8; j++)
			{
				int node_id = _cell_list[cell_id]._node_ptr[j]->_id - 1;
				double Nj = N_brick_8(j, xi, eta, zeta);
				for (int k = 0; k < 3; k++)
				{
					int freedom_id = _ID[3 * node_id + k];
					if (freedom_id>0)
					{
						_Fp[freedom_id - 1] = _Fp[freedom_id - 1] - Nj * force_P2C[k];
					}
				}
			}
		}
	}
	return;
}
void Tgrid_Polymer::assemble_reacting_froce_from_CNT(Tgrid_CNT* grid_CNT)
{
	int nump_CNT = grid_CNT->_nump;
	for (int i = 0; i < nump_CNT; i++)
	{
		if (grid_CNT->_node_list[i]._is_surface)
		{
			int cell_id = grid_CNT->_node_list[i]._location->cell_id;
			double xi = grid_CNT->_node_list[i]._location->iso_coor[0];
			double eta = grid_CNT->_node_list[i]._location->iso_coor[1];
			double zeta = grid_CNT->_node_list[i]._location->iso_coor[2];
			double force_P2C[3];
			for (int j = 0; j < 3; j++)
			{
				int fixed_id = -grid_CNT->_ID[3 * i + j] - 1;
				force_P2C[j] = grid_CNT->_reacting_force[fixed_id];
			}
			for (int j = 0; j < 8; j++)
			{
				int node_id = _cell_list[cell_id]._node_ptr[j]->_id - 1;
				double Nj = N_brick_8(j, xi, eta, zeta);
				for (int k = 0; k < 3; k++)
				{
					int freedom_id = _ID[3 * node_id + k];
					if (freedom_id<0)
					{
						int fixed_id = -freedom_id - 1;
						_reacting_force[fixed_id] = _reacting_force[fixed_id] + Nj * force_P2C[k];
					}
				}
			}
		}
	}
	return;
}
void Tgrid_Polymer::Calculate_stiffness_matrix()
{
	calculate_ID();
	calculate_iK_jK();
	assemble_equations();
	_r = new double[_num_freedom_degree];
	_Fp = new double[_num_freedom_degree];
	_F0 = new double[_num_freedom_degree];
	_Kp = new double[_num_freedom_degree];
	for (int i = 0; i < _num_freedom_degree; i++)
	{
		_r[i] = _Fp[i] = _F0[i] = _Kp[i] = 0;
	}
	return;
}
void Tgrid_Polymer::initialize_Fp()
{
	for (int i = 0; i < _num_freedom_degree; i++)
	{
		_Fp[i] = 0;
	}
}
void Tgrid_Polymer::Set_F0()
{
	for (int i = 0; i < _num_freedom_degree; i++)
	{
		_F0[i] = _Fp[i];
	}
}
void Tgrid_Polymer::Calculate_Kp()
{
	char uplo = 'u';
	mkl_dcsrsymv(&uplo, &_num_freedom_degree, _K, _iK, _jK, _p, _Kp);
	for (int i = 0; i < _num_freedom_degree; i++)
	{
		_Kp[i] = _Kp[i] - (_Fp[i] - _F0[i]);
	}
	return;
}
void Tgrid_Polymer::CG_update()
{
	double rr, pKp;
	rr = cblas_ddot(_num_freedom_degree, _r, 1, _r, 1);
	pKp = cblas_ddot(_num_freedom_degree, _p, 1, _Kp, 1);
	double alpha = rr / pKp;
	for (int i = 0; i < _num_freedom_degree; i++)
	{
		_dis[i] = _dis[i] + alpha * _p[i];
		_r[i] = _r[i] - alpha * _Kp[i];
	}
	double rr1= cblas_ddot(_num_freedom_degree, _r, 1, _r, 1);
	double b = rr1 / rr;
	for (int i = 0; i < _num_freedom_degree; i++)
	{
		_p[i] = _r[i] + b * _p[i];
	}
	return;
}
void Tgrid_CNT::calculate_CNT_boundary_p(Tgrid_Polymer* grid_polymer)
{
	for (int i = 0; i < _nump; i++)
	{
		if (_node_list[i]._is_surface)
		{
			_node_list[i]._bc_value[0] = _node_list[i]._bc_value[1] = _node_list[i]._bc_value[2] = 0;
			int cell_id = _node_list[i]._location->cell_id;
			double xi = _node_list[i]._location->iso_coor[0];
			double eta = _node_list[i]._location->iso_coor[1];
			double zeta = _node_list[i]._location->iso_coor[2];
			double cell_d[8][3];
			grid_polymer->calculate_cell_p(cell_id, cell_d);
			for (int j = 0; j < 8; j++)
			{
				double Nj = N_brick_8(j, xi, eta, zeta);
				_node_list[i]._bc_value[0] = _node_list[i]._bc_value[0] + Nj * cell_d[j][0];
				_node_list[i]._bc_value[1] = _node_list[i]._bc_value[1] + Nj * cell_d[j][1];
				_node_list[i]._bc_value[2] = _node_list[i]._bc_value[2] + Nj * cell_d[j][2];
			}
		}
	}
}
void Tgrid_Polymer::calculate_cell_p(int cell_id, double(&p_cell)[8][3])
{
	int node_id;
	for (int i = 0; i < 8; i++)
	{
		node_id = _cell_list[cell_id]._node_ptr[i]->_id - 1;
		for (int j = 0; j < 3; j++)
		{
			if (_ID[3 * node_id + j] > 0)
			{
				p_cell[i][j] = _p[_ID[3 * node_id + j] - 1];
			}
			else if (_ID[3 * node_id + j] < 0)
			{
				p_cell[i][j] = _node_list[node_id]._bc_value[j];
			}
		}
	}
	return;
}
double Tgrid_Polymer::calculate_diff()
{
	return abs(_r[cblas_idamax(_num_freedom_degree, _r, 1)]);
}
void Tgrid_Polymer::CG_initialization()
{
	//ofstream output1, output2;
	//output1.open("_p.k");
	//output2.open("_f.k");
	for (int i = 0; i < _num_freedom_degree; i++)
	{
		_r[i] = _p[i]+_F0[i];
		//if (isnan(_p[i]))
		//{
		//	output1 << i << " " << _p[i] << endl;
		//}
		//if (isnan(_F0[i]))
		//{
		//	output2 << i << " " << _F0[i] << endl;
		//}
		_p[i] = _r[i];
	}
}
void Tgrid_CNT::output_tecplot(ofstream& output)
{
	//Calculate cell stress
	Sstress** cell_stress;
	cell_stress = Calculate_cell_stress();
	Snode_stress* node_stress;
	node_stress = new Snode_stress[_nump];
	int node_id;
	for (int i = 0; i < _nume; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			node_id = _cell_list[i]._node_ptr[j]->_id - 1;
			node_stress[node_id].num = node_stress[node_id].num + 1;
			node_stress[node_id].add(cell_stress[i][j]);
		}
	}
	for (int i = 0; i < _nump; i++)
	{
		node_stress[i].Average();
	}
	for (int i = 0; i < _nume; i++)
	{
		delete cell_stress[i];
	}
	delete cell_stress;
	//Update the grid position
	double** pos;
	pos = new double*[_nump];
	for (int i = 0; i < _nump; i++)
	{
		pos[i] = new double[3];
	}
	for (int i = 0; i < _nump; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (_ID[3 * i + j] > 0)
			{
				pos[i][j] = _node_list[i]._pos[j] + _dis[_ID[3 * i + j] - 1];
			}
			else if (_ID[3 * i + j] < 0)
			{
				pos[i][j] = _node_list[i]._pos[j] + _node_list[i]._bc_value[j];
			}
		}
	}
	//Output the Tecplot result
	output << "TITLE = \"Tecplot Grid\"" << endl;
	output << "VARIABLES = \"X\",\"Y\",\"Z\",\"SXX\",\"SYY\",\"SZZ\",\"SXY\",\"SXZ\",\"SYZ\",\"EQS\"" << endl;
	output << "ZONE F=FEPOINT,N=" << _nump << "," << "E=" << _nume << "," << "ET=BRICK" << endl;
	for (int i = 0; i < _nump; i++)
	{
		output << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << " ";
		output << node_stress[i].sxx << " " << node_stress[i].syy << " " << node_stress[i].szz << " ";
		output << node_stress[i].sxy << " " << node_stress[i].sxz << " " << node_stress[i].syz << " ";
		output << node_stress[i].equivalent_stress() << endl;
	}
	for (int i = 0; i < _nume; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			output << _cell_list[i]._node_ptr[j]->_id << " ";
		}
		output << endl;
	}
	delete node_stress;
	for (int i = 0; i < _nump; i++)
	{
		delete pos[i];
	}
	delete pos;
	return;
}