#include "Tmaterial_grid.h"
void Tmaterial_grid::input_material_grid(ifstream &input)
{
	Skeyword keyword;
	read_in_keyword_file(input,keyword);
	_nume=keyword.cell_4_list.size();
	_nump=keyword.node_list.size();
	//Read in the material state information of this material
	int material_id=keyword.part_list[0].material_id;
	int EOS_id=keyword.part_list[0].EOS_id;
	if (EOS_id==0 || keyword.EOS_list.size()==0)
	{
		cout<<"Error:There must be an EOS in material_grid"<<endl;
		system("Pause");
		exit(0);
	}
	_material=keyword.material_list[material_id-1];
	_EOS=keyword.EOS_list[EOS_id-1];
	//Read in the material tetrahedron
	vec3D p[4];
	for (int i = 0; i < _nump; i++)
	{
		vec3D new_node;
		new_node.value(keyword.node_list[i].x,keyword.node_list[i].y,keyword.node_list[i].z);
		_node_list.push_back(new_node);
	}
	for (int i = 0; i < _nume; i++)
	{
		SIEN new_IEN;
		for (int j = 0; j < 4; j++)
		{
			int nid=keyword.cell_4_list[i].IEN[j]-1;	
			new_IEN.IEN[j]=nid+1;
			p[j].value(keyword.node_list[nid].x,keyword.node_list[nid].y,keyword.node_list[nid].z);
		}
		Ttetrahedron new_tet(p[0],p[1],p[2],p[3]);
		_material_tet_list.push_back(new_tet);
		_IEN_list.push_back(new_IEN);
	}
	calculate_geometry_property();
	return;
}
void Tmaterial_grid::calculate_geometry_property()
{
	for (int i = 0; i < _nume; i++)
	{
		double cell_volume=_material_tet_list[i].G_volume();
		_total_volume=_total_volume+cell_volume;
		_total_centroid=_total_centroid+_material_tet_list[i].G_centroid()*cell_volume;
	}
	_total_centroid=_total_centroid/_total_volume;
	return;
}
void Tmaterial_grid::plot_material_grid(ofstream &output)
{
	output<<"TITLE=\"material mesh\""<<endl;
	output<<"Variables=\"X\",\"Y\",\"Z\",\"mid\""<<endl;
	output<<"ZONE,F=FEPOINT N="<<_nump<<" E="<<_nume<<" ET=TETRAHEDRON"<<endl;
	for (int i = 0; i < _nump; i++)
	{
		vec3D node=_node_list[i];
		output<<node.x<<" "<<node.y<<" "<<node.z<<" "<<_material_id<<endl;
	}
	for (int i = 0; i < _nume; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			output<<_IEN_list[i].IEN[j]<<" ";
		}
		output<<endl;
	}
	return;
}