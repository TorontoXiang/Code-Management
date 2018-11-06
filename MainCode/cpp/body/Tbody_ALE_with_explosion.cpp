#include "Tbody_ALE_with_explosion.h"
#include "Tcell_mixed_fluid.h"
#include "Tcell_intersect.h"
#include <iomanip>
Tbody_ALE_with_explosion::Tbody_ALE_with_explosion(int i, int num_material, string body_name) :Tbody_ALE(i, num_material, body_name) {};
void Tbody_ALE_with_explosion::input_body(ifstream& input)
{
	Skeyword keyword;
	read_in_keyword_file(input, keyword);
	//Read the body global information from keyword
	_nume = keyword.cell_8_list.size();
	_nump = keyword.node_list.size();
	_endtime = keyword.time_control.endtime;
	_remapping_start_time = -1;
	_CFL = keyword.time_control.CFL;
	_hourglass_option = keyword.hourglass_option;
	//----------------------------------------------------
	//Read node list
	//----------------------------------------------------
	int id;
	double x, y, z;
	_grid1._node_list.resize(_nump);
	_grid2._node_list.resize(_nump);
	_nodeptr_list.resize(_nump);
	for (int i = 0; i < _nump; i++)
	{
		id = keyword.node_list[i].id;
		x = keyword.node_list[i].x;
		y = keyword.node_list[i].y;
		z = keyword.node_list[i].z;
		Tnode_fluid new_node(id, x, y, z);
		_grid1._node_list[i] = new Tnode_fluid(id, x, y, z);
		_grid2._node_list[i] = new Tnode_fluid(id, x, y, z);
	}
	//-----------------------------------------------------
	//Read cell list
	//-----------------------------------------------------
	_grid1._cell_list.resize(_nume);
	_grid2._cell_list.resize(_nume);
	_cellptr_list.resize(_nume);
	int part_id, material_id, EOS_id;
	Tnode_fluid* node_ptr1[8];
	Tnode_fluid* node_ptr2[8];
	int mid_TNT;
	for (int i = 0; i < _nume; i++)
	{
		//Read cell id and IEN
		id = keyword.cell_8_list[i].cell_id;
		for (int j = 0; j < 8; j++)
		{
			int IENj = keyword.cell_8_list[i].IEN[j];
			node_ptr1[j] = _grid1._node_list[IENj - 1];
			node_ptr2[j] = _grid2._node_list[IENj - 1];
		}
		//Allocate the materials in this problem
		TMAT_base** mat1;
		TMAT_base** mat2;
		mat1 = new TMAT_base*[_num_material];
		mat2 = new TMAT_base*[_num_material];
		for (int j = 0; j < keyword.part_list.size(); j++)
		{
			EOS_id = keyword.part_list[j].EOS_id;
			if (EOS_id == 0 || keyword.EOS_list.size() == 0)
			{
				cout << "Error:There must be an EOS in body_explosion" << endl;
				system("Pause");
				exit(0);
			}
			material_id = keyword.part_list[j].material_id;
			Smaterial this_mat = keyword.material_list[material_id - 1];
			SEOS this_EOS = keyword.EOS_list[EOS_id - 1];
			if (this_EOS.EOS_type=="*EOS_JWL")
			{
				mid_TNT = material_id;
			}
			mat1[j] = generate_material(this_mat, this_EOS, 1);
			mat2[j] = generate_material(this_mat, this_EOS, 1);
		}
		//Create a new cell
		_grid1._cell_list[i] = new Tcell_mixed_fluid(id, _num_material, node_ptr1, mat1);
		_grid2._cell_list[i] = new Tcell_mixed_fluid(id, _num_material, node_ptr2, mat2);
		//Initialize the volume fraction and material centroid
		part_id = keyword.cell_8_list[i].part_id;
		material_id = keyword.part_list[part_id - 1].material_id;
		_grid1._cell_list[i]->S_material_fraction(1, material_id - 1);
		_grid1._cell_list[i]->S_material_centroid(_grid1._cell_list[i]->G_cell_polyhedron()->G_centroid(), material_id - 1);
		_grid2._cell_list[i]->S_material_fraction(1, material_id - 1);
		_grid2._cell_list[i]->S_material_centroid(_grid2._cell_list[i]->G_cell_polyhedron()->G_centroid(), material_id - 1);
		if (material_id==mid_TNT)
		{
			_grid1._cell_list[i]->_flag = 1;
		}
		else
		{
			_grid1._cell_list[i]->_flag = 0;
		}
	}
	//-------------------------------------------------------
	//Read boundary condition
	//-------------------------------------------------------
	int node_group_id;   //Boundary condition (load) is applied on node_group_id

	int num_bd_group;    //Number of boundary condition group
	int num_node;        //Number of nodes in the node group
	num_bd_group = keyword.boundary_list.size();
	int bc_type;
	for (int i = 0; i < num_bd_group; i++)
	{
		node_group_id = keyword.boundary_list[i].id;
		num_node = keyword.node_group_list[node_group_id - 1].node_id.size();
		bc_type = keyword.boundary_list[i].type;
		if (bc_type == 0)
		{
			for (int j = 0; j < 3; j++)
			{
				if (keyword.boundary_list[i].pos[j] == 1)
				{
					for (int k = 0; k < num_node; k++)
					{
						_grid1._node_list[keyword.node_group_list[node_group_id - 1].node_id[k] - 1]->_bc_type_position[j] = 1;
						_grid2._node_list[keyword.node_group_list[node_group_id - 1].node_id[k] - 1]->_bc_type_position[j] = 1;
					}
				}
			}
		}
		else if (bc_type == 1)
		{
			for (int j = 0; j < 3; j++)
			{
				if (keyword.boundary_list[i].dofvel[j] == 1)
				{
					for (int k = 0; k < num_node; k++)
					{
						_grid1._node_list[keyword.node_group_list[node_group_id - 1].node_id[k] - 1]->_bc_type_position[j] = 2;
						_grid1._node_list[keyword.node_group_list[node_group_id - 1].node_id[k] - 1]->_bc_velocity[j] = keyword.boundary_list[i].vel[j];
						_grid2._node_list[keyword.node_group_list[node_group_id - 1].node_id[k] - 1]->_bc_type_position[j] = 2;
						_grid2._node_list[keyword.node_group_list[node_group_id - 1].node_id[k] - 1]->_bc_velocity[j] = keyword.boundary_list[i].vel[j];
					}
				}
			}
		}
		else if (bc_type == 2 || bc_type == 3)
		{
			vec2D center(keyword.boundary_list[i].x, keyword.boundary_list[i].y);
			double v = keyword.boundary_list[i].v;
			for (int k = 0; k < num_node; k++)
			{
				for (int j = 0; j < 3; j++)
				{
					_grid1._node_list[keyword.node_group_list[node_group_id - 1].node_id[k] - 1]->_bc_type_position[j] = 3;
					_grid2._node_list[keyword.node_group_list[node_group_id - 1].node_id[k] - 1]->_bc_type_position[j] = 3;
				}
			}
		}
	}
	//---------------------------------------------------------------------------------
	//Read external load
	//---------------------------------------------------------------------------------
	int num_load_group;          //Number of load group
	int curve_id;                //The curve id of the load
	int direction;               //The load direction
	double magnitude;            //The load magnitude
	Sexternal_load new_load;
	num_load_group = keyword.load_list.size();
	for (int i = 0; i < num_load_group; i++)
	{
		node_group_id = keyword.load_list[i].id;
		new_load.node_id = keyword.node_group_list[node_group_id - 1].node_id;
		new_load.num_node = new_load.node_id.size();
		direction = keyword.load_list[i].dof;
		new_load.direction = direction;
		curve_id = keyword.load_list[i].curve_id;
		magnitude = keyword.curve_list[curve_id - 1].v2[0];      //Only for constant load
		new_load.magnitude = magnitude;
		_external_load_list.push_back(new_load);
	}
	overlapping_with_explosion_region();

	//_grid_old point to grid1 and _grid_new point to grid2 initially
	_grid_old = &_grid1; _grid_new = &_grid2;
	bool negative;
	for (int i = 0; i < _nume; i++)
	{
		_grid1._cell_list[i]->calculate_cell_volume("Full", negative);
	}

	//Get the cell connection for grid2
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i] = _grid2._cell_list[i];
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i] = _grid2._node_list[i];
	}
	set_grid_topology();

	//Get the cell connection for grid1
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i] = _grid1._cell_list[i];
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i] = _grid1._node_list[i];
	}
	set_grid_topology();

	//Initialize the nodal variables
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->calculate_nodal_variables();
	}
	//Clear the variables on grid2
	for (int i = 0; i < _nume; i++)
	{
		_grid2._cell_list[i]->clear_cell_variable();
	}
	for (int i = 0; i < _nump; i++)
	{
		_grid2._node_list[i]->clear_nodal_variable();
	}

	//_grid_old point to grid1 and _grid_new point to grid2 initially
	_grid_old = &_grid1; _grid_new = &_grid2;

	//Open the file for remapping information
	_iremapping.open("Remapping_info.txt");
	return;
}
void Tbody_ALE_with_explosion::overlapping_with_explosion_region()
{
	//Clear the information in the explosion region
	for (int i = 0; i < _nume; i++)
	{
		if (_grid1._cell_list[i]->_flag==1)
		{
			_grid1._cell_list[i]->clear_cell_variable();
		}
	}
	//Read in explosion grid
	Sfluid_grid* grid_exp=input_explosion_region();
	int nume = grid_exp->_cell_list.size();

	//----------------------------------------prepare for overlapping---------------------------------------------
	//Reconstruct the variable gradient and material surface on explosion grid
	old_grid_reconstruction();
	//Generate a bucket search for explosion grid
	vec3D coor_min, coor_max;
	double cell_edge_max;
	calculate_explosion_region_size(grid_exp,coor_min, coor_max, cell_edge_max);
	Tbucket_searching buckets(coor_min, coor_max, cell_edge_max);
	bool access;
	for (int i = 0; i < nume; i++)
	{
		int ix, iy, iz;
		vec3D centroid = grid_exp->_cell_list[i]->_cell_polyhedron.G_centroid();
		buckets.search_bucket_ijk(centroid, ix, iy, iz, access);
		buckets.add_element(i, ix, iy, iz);
		if (access)
		{
			cout << "Error: The cell centroid is out of body range" << endl;
			system("Pause");
			exit(0);
		}
	}

	//------------------------------remapping begin-------------------------------------------
	Tcell_fluid_base* new_cell; Tcell_fluid_base* exp_cell;
	int nx_min, ny_min, nz_min, nx_max, ny_max, nz_max;
	double error = 0, error_max = 0; int n_error = 0;
	vec3D vertex_list[8];
	vector<int> error_id;
	Tcell_intersect intersection;
	cout << "remapping beginning!" << endl;
	cout << "Completing ratio: ";
	int cell_interval = maxval(_nume / 1000, 1);
	bool first = true;
	double num_intersection = 0;
	#pragma omp parallel for schedule(dynamic) private(intersection,nx_min,ny_min,nz_min,nx_max,ny_max,nz_max,exp_cell,new_cell,vertex_list)
	for (int i = 0; i < _nume; i++)
	{
		vec3D coor_min_new, coor_max_new;
		new_cell = _grid1._cell_list[i];
		if (new_cell->_flag==1)
		{
			if (i%cell_interval == 0)
			{
				if (first)
				{
					cout << setw(4) << 100 * i / _nume << "%";
					first = false;
				}
				else
				{
					cout << '\b' << '\b' << '\b' << '\b' << '\b';
					cout << setw(4) << 100 * i / _nume << "%";
				}
			}
			for (int j = 0; j < 8; j++)
			{
				vertex_list[j] = new_cell->_node_ptr[j]->_position;
			}
			buckets.search_bucket_range(vertex_list, 8, nx_min, ny_min, nz_min, nx_max, ny_max, nz_max, access);
			if (access)
			{
				cout << "Error: The new grid is out of the range of old grid" << endl;
				system("Pause");
				exit(0);
			}
			intersection.S_new_cell(new_cell);
			intersection.create_tet_list();
			for (int i1 = nx_min - 1; i1 <= nx_max + 1; i1++)
			{
				for (int i2 = ny_min - 1; i2 <= ny_max + 1; i2++)
				{
					for (int i3 = nz_min - 1; i3 <= nz_max + 1; i3++)
					{
						Sbucket* the_bucket = &buckets._bucket_list[i1][i2][i3];
						for (int j = 0; j < the_bucket->_nume; j++)
						{
							exp_cell = grid_exp->_cell_list[the_bucket->eleid[j]];
							vec3D coor_min_old, coor_max_old;
							if (!comparison_cell_coordinate(exp_cell->_coor_min, exp_cell->_coor_max, new_cell->_coor_min, new_cell->_coor_max))
							{
								intersection.S_old_cell(exp_cell);
								intersection.intersect();
								num_intersection = num_intersection + 1;
							}
						}
					}
				}
			}
		}
	}
	cout << endl;
	//------------------Update the variables on new grid and out put the remapping result------------
	double volume_error, av_volume_error = 0, max_volume_error = 0;
	Tcell_fluid_base* newcell_ptr = NULL;
	bool overlapping_failed = false;
	int nume_act = 0;
	for (int i = 0; i < _nume; i++)
	{
		newcell_ptr = _grid1._cell_list[i];
		if (newcell_ptr->_flag == 1)
		{
			nume_act = nume_act + 1;
			volume_error = abs(newcell_ptr->_volume_intersection - newcell_ptr->_cell_volume);
			if (volume_error > 1e-10)
			{
				cout << "Error:Overlapping fails in cell  " << i << " with volume error= " << volume_error << endl;
				overlapping_failed = true;
			}
			if (volume_error > max_volume_error)
			{
				max_volume_error = volume_error;
			}
			av_volume_error = av_volume_error + volume_error;
			//Update the variables in new grid
			newcell_ptr->reset_state_after_overlapping_exp();
		}
	}
	cout << "------Overlapping result---------" << endl;
	cout << "The average volume error in this remapping phase is: " << av_volume_error / nume_act << endl;
	if (overlapping_failed)
	{
		system("Pause");
	}
	return;
}
void Tbody_ALE_with_explosion::calculate_explosion_region_size(Sfluid_grid* grid_exp,vec3D& coor_min, vec3D& coor_max, double cell_edge_max)
{
	Tcell_fluid_base* cell_ptr;
	cell_edge_max = 0;
	double l_ele;
	coor_min.value(1e10, 1e10, 1e10); coor_max.value(-1e10, -1e10, -1e10);
	int nume = grid_exp->_cell_list.size();
	for (int i = 0; i < nume; i++)
	{
		cell_ptr = grid_exp->_cell_list[i];
		l_ele = cell_ptr->calculate_maximal_edge();
		if (l_ele > cell_edge_max)
		{
			cell_edge_max = l_ele;
		}
	}
	int nump = grid_exp->_node_list.size();
	for (int i = 0; i < nump; i++)
	{
		for (int j = 1; j < 4; j++)
		{
			if (grid_exp->_node_list[i]->_position.access(j) < coor_min.access(j)) coor_min.value(j, grid_exp->_node_list[i]->_position.access(j));
			if (grid_exp->_node_list[i]->_position.access(j) > coor_max.access(j)) coor_max.value(j, grid_exp->_node_list[i]->_position.access(j));
		}
	}
	return;
}
Tbody_ALE_with_explosion::Sfluid_grid* Tbody_ALE_with_explosion::input_explosion_region()
{
	//Read in the explosion information
	ifstream input;
	input.open("Explosion_data.k");
	if (input.fail())
	{
		cout << "Error:File Explosion_data.k does not exist for Tbody_ALE_with_explosion!" << endl;
		system("Pause");
		exit(0);
	}
	Skeyword keyword;
	read_in_keyword_file(input, keyword);
	//Generate the explosion region grid
	Sfluid_grid grid_exp;
	int nume_exp = keyword.variable_cell_8_list.size();
	int nump_exp = keyword.vel_node_list.size();
	grid_exp._node_list.resize(nump_exp);
	grid_exp._cell_list.resize(nume_exp);
	//Generate explosion grid node
	int id;
	for (int i = 0; i < nump_exp; i++)
	{
		id = keyword.vel_node_list[i].id;
		double x = keyword.vel_node_list[i].x;
		double y = keyword.vel_node_list[i].y;
		double z = keyword.vel_node_list[i].z;
		double vx = keyword.vel_node_list[i].vx;
		double vy = keyword.vel_node_list[i].vy;
		double vz = keyword.vel_node_list[i].vz;
		Tnode_fluid new_node(id, x, y, z);
		grid_exp._node_list[i] = new Tnode_fluid(id, x, y, z);
		grid_exp._node_list[i]->_velocity.value(vx, vy, vz);
	}
	//Generate explosion grid cell
	int part_id, material_id, EOS_id;
	Tnode_fluid* node_ptr[8];
	double density, internal_energy;
	for (int i = 0; i < nume_exp; i++)
	{
		//Read cell id and IEN
		id = keyword.variable_cell_8_list[i].cell_id;
		density = keyword.variable_cell_8_list[i].density;
		internal_energy = keyword.variable_cell_8_list[i].internal_energy;
		for (int j = 0; j < 8; j++)
		{
			int IENj = keyword.variable_cell_8_list[i].IEN[j];
			node_ptr[j] = grid_exp._node_list[IENj - 1];
		}
		//Allocate the materials in this problem
		TMAT_base** mat;
		mat = new TMAT_base*[_num_material];
		for (int j = 0; j < keyword.part_list.size(); j++)
		{
			EOS_id = keyword.part_list[j].EOS_id;
			if (EOS_id == 0 || keyword.EOS_list.size() == 0)
			{
				cout << "Error:There must be an EOS in body_explosion" << endl;
				system("Pause");
				exit(0);
			}
			material_id = keyword.part_list[j].material_id;
			Smaterial this_mat = keyword.material_list[material_id - 1];
			SEOS this_EOS = keyword.EOS_list[EOS_id - 1];
			mat[j] = generate_material(this_mat, this_EOS, 1);
		}
		//Create a new cell
		grid_exp._cell_list[i] = new Tcell_mixed_fluid(id, _num_material, node_ptr, mat);
		//Initialize the volume fraction and material centroid
		part_id = keyword.cell_8_list[i].part_id;
		material_id = keyword.part_list[part_id - 1].material_id;
		grid_exp._cell_list[i]->S_material_fraction(1, material_id - 1);
		grid_exp._cell_list[i]->S_material_centroid(grid_exp._cell_list[i]->G_cell_polyhedron()->G_centroid(), material_id - 1);
		grid_exp._cell_list[i]->_gausspoint[material_id - 1]->S_density(density);
		grid_exp._cell_list[i]->_gausspoint[material_id - 1]->S_density(internal_energy);
	}
	return &grid_exp;
}