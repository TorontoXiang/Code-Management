#include "readin.h"
#include <string>
#include <iostream>
using namespace std;
void Skeyword::clear_keyword()
{
	node_list.resize(0);
	cell_8_list.resize(0);
	node_group_list.resize(0);
	part_list.resize(0);
	material_list.resize(0);
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
			double temp;
			next_data(input);
			input>>new_material.material_id;
			input>>temp;
			input>>new_material.Youngs>>new_material.Possion;
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
				input>> new_cell.cell_id>> new_cell.mid;
				input>> new_cell.IEN[0]>> new_cell.IEN[1]>> new_cell.IEN[2]>> new_cell.IEN[3];
				input>> new_cell.IEN[4]>> new_cell.IEN[5]>> new_cell.IEN[6]>> new_cell.IEN[7];
				keyword.cell_8_list.push_back(new_cell);
			}
		}
		else if (a=="*REGULAR_GRID")
		{
			next_data(input);
			input>>keyword.regular_grid.x_min[0]>>keyword.regular_grid.x_min[1]>>keyword.regular_grid.x_min[2];
			next_data(input);
			input>>keyword.regular_grid.x_max[0]>>keyword.regular_grid.x_max[1]>>keyword.regular_grid.x_max[2];
			next_data(input);
			input>>keyword.regular_grid.nx>>keyword.regular_grid.ny>>keyword.regular_grid.nz;
		}
		else if (a=="*REGULAR_GRID_BD")
		{
			Sregular_grid_bc new_bc;
			next_data(input);
			input >> new_bc.FaceName >> new_bc.direction >> new_bc.value;
			keyword.regular_bc_list.push_back(new_bc);
		}
		else if (a=="*END")
		{
			return;
		}
	}
	return;
}