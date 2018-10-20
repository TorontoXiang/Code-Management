#include "readin.h"
#include <string>
#include <iostream>
#include "public_function.h"
#include "TMATIsotropic_elastic.h"
#include "TMATIsotropic_plastic.h"
using namespace std;
void Skeyword::clear_keyword()
{
	node_list.resize(0);
	cell_8_list.resize(0);
	boundary_list.resize(0);
	node_group_list.resize(0);
	load_list.resize(0);
	part_list.resize(0);
	material_list.resize(0);
	//time_control;
	//background;
	initial_velocity_list.resize(0);
}
void new_line(ifstream& input)
{
	char a;
	do
	{
		input.get(a);
	} while (a!='\n');
}
bool exam_keyword(ifstream& input)
{
	char a;
	new_line(input);
	do
	{
		input.get(a);
		if (a=='*')
		{
			input.putback(a);
			return true;
		}
		else if (a=='$')
		{
			new_line(input);
		}
		else
		{
			input.putback(a);
			return false;
		}
	} while (a=='$');
	return 0;
}
void next_keyword(ifstream& input)
{
	char a;
	do
	{
		input.get(a);
		if (a=='*')
		{
			input.putback(a);
			return;
		}
		else  if (a!='\n')
		{
			//If a=='\n' the input stream will automatically come to the new line
			new_line(input);
		}
	} while (true);
	return;
}
void next_data(ifstream& input)
{
	new_line(input);
	char a;
	do
	{
		input.get(a);
		if (a=='$')
		{
			new_line(input);
		}
		else
		{
			input.putback(a);
			return;
		}
	} while (true);
}
TMAT_base* generate_material(Smaterial mat,int type)
{
	TMAT_base* mat_ptr=NULL;
	if (mat.material_type=="*MAT_ELASTIC")
	{
		TMATIsotropic_elastic* new_mat;
		new_mat=new TMATIsotropic_elastic(mat.density,mat.Youngs,mat.Possion);
		mat_ptr=new_mat;
	}
	else if (mat.material_type == "*MAT_PLASTIC_KINEMATIC")
	{
		TMATIsotropic_plastic* new_mat;
		new_mat = new TMATIsotropic_plastic(mat.density, mat.Youngs, mat.Possion, mat.Sigma_y, mat.ET);
		mat_ptr = new_mat;
	}
	return mat_ptr;
}
//-------------------------------------------------------------------
//Read in keyword file
//-------------------------------------------------------------------
void read_in_keyword_file(ifstream& input,Skeyword& keyword)
{
	string a;
	string temp;
	while (true)
	{
		next_keyword(input);
		input>>a;
		cout<<a<<endl;
		if (a=="*NODE")
		{
			while (!exam_keyword(input))
			{
				Snode_temp new_node;
				input>>new_node.id>>new_node.x>>new_node.y>>new_node.z;
				keyword.node_list.push_back(new_node);
			}
		}
		else if (a=="*CONTROL_TERMINATION")
		{
			next_data(input);
			input>>keyword.time_control.endtime;
		}
		else if (a=="*CONTROL_TIMESTEP")
		{
			next_data(input);
			input>>temp>>keyword.time_control.CFL;
		}
		else if (a=="*BOUNDARY_SPC_SET")
		{
			next_data(input);
			Sboundary_type new_bd;
			input>>new_bd.id>>temp>>new_bd.pos[0]>>new_bd.pos[1]>>new_bd.pos[2];
			input>>new_bd.rot[0]>>new_bd.rot[1]>>new_bd.rot[2];
			new_bd.type=0;
			keyword.boundary_list.push_back(new_bd);
		}
		else if (a=="*VELOCITY_BOUNDARY")
		{
			next_data(input);
			Sboundary_type new_bd;
			input>>new_bd.id>>new_bd.dofvel[0]>>new_bd.dofvel[1]>>new_bd.dofvel[2]>>new_bd.vel[0]>>new_bd.vel[1]>>new_bd.vel[2];
			new_bd.type=1;
			keyword.boundary_list.push_back(new_bd);
		}
		else if (a=="*SET_NODE_LIST")
		{
			next_data(input);
			Snode_group new_group;
			input>>new_group.group_id;
			while (!exam_keyword(input))
			{
				int node_id;
				for (int i = 0; i < 8; i++)
				{
					input>>node_id;
					if (node_id>0)
					{
						new_group.node_id.push_back(node_id);
					}
				}
			}
			//Add new_group to node_group_list and keep in order
			unsigned int group_id=new_group.group_id;
			if (group_id>keyword.node_group_list.size())
			{
				keyword.node_group_list.resize(group_id);
			}
			keyword.node_group_list[group_id-1]=new_group;
		}
		else if (a=="*LOAD_NODE_SET")
		{
			while (!exam_keyword(input))
			{
				Sload_on_node new_load;
				input>>new_load.id>>new_load.dof>>new_load.curve_id;
				keyword.load_list.push_back(new_load);
			}
		}
		else if (a=="*PART")
		{
			Spart new_part;
			next_data(input);
			next_data(input);
			input>>new_part.part_id>>new_part.section_id>>new_part.material_id>>new_part.EOS_id;
			//Add new_part to part_list and keep in order
			unsigned int part_id=new_part.part_id;
			if (part_id>keyword.part_list.size())
			{
				keyword.part_list.resize(part_id);
			}
			keyword.part_list[part_id-1]=new_part;
		}
		else if (a=="*MAT_ELASTIC")
		{
			Smaterial new_material;
			next_data(input);
			input>>new_material.material_id;
			input>>new_material.density;
			input>>new_material.Youngs>>new_material.Possion;
			new_material.material_type=a;
			//Add new_material into material_list and keep in order
			unsigned int material_id=new_material.material_id;
			if (material_id>keyword.material_list.size())
			{
				keyword.material_list.resize(material_id);
			}
			keyword.material_list[material_id-1]=new_material;
		}
		else if (a=="*MAT_PLASTIC_KINEMATIC")
		{
			Smaterial new_material;
			next_data(input);
			input>>new_material.material_id;
			input>>new_material.density>>new_material.Youngs>>new_material.Possion;
			input>>new_material.Sigma_y>>new_material.ET;
			new_material.material_type=a;
			//Add new_material into material_list and keep in order
			unsigned int material_id=new_material.material_id;
			if (material_id>keyword.material_list.size())
			{
				keyword.material_list.resize(material_id);
			}
			keyword.material_list[material_id-1]=new_material;
		}
		else if (a=="*ELEMENT_SOLID")
		{
			while (!exam_keyword(input))
			{
				Scell_8 new_cell;
				input>> new_cell.cell_id>> new_cell.part_id;
				input>> new_cell.IEN[0]>> new_cell.IEN[1]>> new_cell.IEN[2]>> new_cell.IEN[3];
				input>> new_cell.IEN[4]>> new_cell.IEN[5]>> new_cell.IEN[6]>> new_cell.IEN[7];
				keyword.cell_8_list.push_back(new_cell);
			}
		}
		else if (a=="*REGULAR_GRID")
		{
			keyword.is_regular_grid = true;
			next_data(input);
			input>>keyword.regular_grid.x_min.x>>keyword.regular_grid.x_min.y>>keyword.regular_grid.x_min.z;
			next_data(input);
			input>>keyword.regular_grid.x_max.x>>keyword.regular_grid.x_max.y>>keyword.regular_grid.x_max.z;
			next_data(input);
			input>>keyword.regular_grid.nx>>keyword.regular_grid.ny>>keyword.regular_grid.nz;
		}
		else if (a=="*REGULAR_GRID_BD")
		{
			Sregular_grid_bc new_bc;
			next_data(input);
			input >> new_bc.FaceName >> new_bc.bc_type >> new_bc.direction >> new_bc.velocity;
			keyword.regular_bc_list.push_back(new_bc);
		}
		else if (a=="*INITIAL_VELOCITY")
		{
			next_data(input);
			Sinitial_velocity_temp new_initial_velocity;
			input>>new_initial_velocity.part_id>>new_initial_velocity.velocity.x>>new_initial_velocity.velocity.y;
			input>>new_initial_velocity.velocity.z;
			keyword.initial_velocity_list.push_back(new_initial_velocity);
		}
		else if (a=="*VELOCITY_INITIAL")
		{
			next_data(input);
			vec3D v;
			input >> v.x >> v.y >> v.z;
			keyword.is_initial_vel = true;
			keyword.initial_vel = v;
		}
		else if (a == "*DEFINE_CURVE")
		{
			double v1, v2;
			Scurve new_curve;
			next_data(input);
			input >> new_curve.curve_id;
			//new_line(input);
			while (!exam_keyword(input))
			{
				input >> v1 >> v2;
				//new_line(input);
				new_curve.v1.push_back(v1);
				new_curve.v2.push_back(v2);
			}
			//Add new_curve to curve_list and keep in order
			unsigned int curve_id = new_curve.curve_id;
			if (curve_id > keyword.curve_list.size())
			{
				keyword.curve_list.resize(curve_id);
			}
			keyword.curve_list[curve_id - 1] = new_curve;
		}
		else if (a=="*END")
		{
			return;
		}
	}
	return;
}