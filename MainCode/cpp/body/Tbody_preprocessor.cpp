#include "Tbody_preprocessor.h"
#include "public_function.h"
#include <iomanip>
using namespace std;
Tbody_preprocessor::Tbody_preprocessor(int i,string body_name):Tbody_base(i,body_name)
{
	_output.open("modified_part.k");
}
void Tbody_preprocessor::modify_body(ifstream &input)
{
	input_body(input);
	modify_k_file();
	output_modified_part(_output);
	return;
}
void Tbody_preprocessor::input_body(ifstream &input)
{
	read_in_keyword_file(input,_keyword);
	return;
}
void Tbody_preprocessor::modify_k_file()
{
	int num_modified_node=_keyword.node_group_list[0].node_id.size();
	vec3D coor,coor_sphere;
	double eps=0.0005,dis;
	int node_id;
	for (int i = 0; i < num_modified_node; i++)
	{
		node_id=_keyword.node_group_list[0].node_id[i];
		coor.value(_keyword.node_list[node_id-1].x,_keyword.node_list[node_id-1].y,_keyword.node_list[node_id-1].z);
		//coor_sphere=convert_to_sphere_coor(coor);
		double r=coor.get_length();
		double theta=acos(coor.x/r);
		//dis=eps*sin(8*coor_sphere.z);
		coor.x=coor.x*(1+eps*cos(8.0*theta));
		coor.y=coor.y*(1+eps*cos(8.0*theta));
		_keyword.node_list[node_id-1].x=coor.x;
		_keyword.node_list[node_id-1].y=coor.y;
		_keyword.node_list[node_id-1].z=coor.z;
	}
}
void Tbody_preprocessor::output_modified_part(ofstream &output)
{
	int nump=_keyword.node_list.size();
	output.precision(8);
	for (int i = 0; i < nump; i++)
	{
		output<<setw(8)<<_keyword.node_list[i].id<<setw(16)<<_keyword.node_list[i].x;
		output<<setw(16)<<_keyword.node_list[i].y<<setw(16)<<_keyword.node_list[i].z;
		output<<setw(8)<<0<<setw(8)<<0<<endl;
	}
	return;
}