#include "Tbody_base.h"
#include <iomanip>
using namespace std;
Tbody_base::Tbody_base(int i,string body_name)
{
	_body_name=body_name;
	_id=i;
	_dt=_endtime=_current_time=_CFL=0;
	_current_step=_nume=_nump=0;
	_mesh_out=NULL;
	_Arti_Vis=NULL;
	_global_output=new Sbody_global_variable(i,false);
	_test_name="no";
	_dimension=3;
	_hourglass_option=true;
}
void Tbody_base::advance_in_time()
{
	_current_time=_current_time+_dt;
	_current_step=_current_step+1;
	return;
}
Tbody_base::Scurve_output::Scurve_output(int body_id,int type,int id) : _body_id(body_id),_type(type),_id(id)
{
	char a_id[32],a_bodyid[32];
	string s_id,s_bodyid;
	_itoa_s(_id,a_id,10);_itoa_s(_body_id,a_bodyid,10);
	s_id=a_id;s_bodyid=a_bodyid;
	if (type==0)
	{
		file_name="Curve_Body"+s_bodyid+"_Node"+s_id+".dat";
		output.open(file_name);
		output<<"Time-history curve for node "<<_id<<endl;
		output<<setw(20)<<"Time"<<setw(20)<<"x-displacement"<<setw(20)<<"y-displacement"<<setw(20)<<"z-displacement"<<endl;
	}
	else if (type==1)
	{
		file_name="Curve_Body"+s_bodyid+"_Cell"+s_id+".dat";
		output.open(file_name);
		output<<"Time-history curve for cell "<<_id<<endl;
	}
	return;
}
Tbody_base::Smesh_output::Smesh_output(int body_id,double amplify_ratio,int total_frame,double endtime) : 
	                      _body_id(body_id),_amplify_ratio(amplify_ratio),_total_frame(total_frame)
{
	char a_bodyid[32];
	string s_bodyid;
	_itoa_s(_body_id,a_bodyid,10);
	 s_bodyid=a_bodyid;
	file_name="Mesh_Body"+s_bodyid+"_initial.dat";
	output_initial.open(file_name);
	file_name="Mesh_Body"+s_bodyid+"_process.dat";
	output_process.open(file_name);
	file_name="Mesh_Body"+s_bodyid+"_final.dat";
	output_final.open(file_name);
	distribution.open("variable_distribution.dat");
	_output_point.resize(_total_frame-1);
	double dt=endtime/_total_frame;
	for (int i = 0; i < _total_frame-1; i++)
	{
		_output_point[i]=dt*(i+1);
	}	
}
bool Tbody_base::Smesh_output::is_output(double current_time,double dt)
{
	for (int i = 0; i < _total_frame-1; i++)
	{
		if (current_time<=_output_point[i])
		{
			if (current_time+dt>_output_point[i])
			{
				_output_point[i]=-1;
				return true;
			}
			else
			{
				return false;
			}
			
		}
	}
	return false;
}
Tbody_base::Scell_monitor::Scell_monitor(int body_id,int cell_id)
{
	char a_id[32],a_bodyid[32];
	string s_id,s_bodyid;
	_body_id=body_id;_cell_id=cell_id;
	_itoa_s(_cell_id,a_id,10);_itoa_s(_body_id,a_bodyid,10);
	s_id=a_id;s_bodyid=a_bodyid;
	file_name="Cell_"+s_id+"deformation_in_body"+s_bodyid+".dat";
	output.open(file_name);
	return;
}
Tbody_base::Sbody_global_variable::Sbody_global_variable(int body_id,bool is_output)
{
	_body_id=body_id;_is_output=is_output;
}
void Tbody_base::Sbody_global_variable::create_ofstream()
{
	_is_output=true;
	char a_bodyid[32];
	string s_bodyid;
	_itoa_s(_body_id,a_bodyid,10);
	s_bodyid=a_bodyid;
	file_name="Global variables for body "+s_bodyid+".dat";
	output.open(file_name);
	return;
}
void Tbody_base::create_curve(int body_id,int type,int id)
{
	Scurve_output* new_curve;
	new_curve=new Scurve_output(body_id,type,id);
	_curve_out_list.push_back(new_curve);
	return;
}
void Tbody_base::create_mesh(int body_id,double amplify_ratio,int total_frame)
{
	Smesh_output* new_mesh=NULL;
	new_mesh= new Smesh_output(body_id,amplify_ratio,total_frame,_endtime);
	_mesh_out=new_mesh;
	return;
}
void Tbody_base::create_cell_monitor(int body_id,int cell_id)
{
	Scell_monitor* new_monitor=NULL;
	new_monitor= new Scell_monitor(body_id,cell_id);
	_cell_monitor_list.push_back(new_monitor);
	return;
}
void Tbody_base::output_mesh(string type)
{
	if (_mesh_out==NULL)
	{
		return;
	}
	if (type=="Initial")
	{
		cout<<"Output grid at initial for body: "<<_id<<endl;
		output_tecplot(_mesh_out->output_initial,_mesh_out->_amplify_ratio);
		//cout<<"after initial plot"<<endl;
		output_distribution(_mesh_out->distribution);
		//cout<<"after initial output distribution"<<endl;
	}
	else if (type=="Process")
	{
		
		if (_mesh_out->is_output(_current_time,_dt))
		{
			cout<<"Output grid at time="<<setw(8)<<_current_time<<" for body: "<<_id<<endl;
			output_tecplot(_mesh_out->output_process,_mesh_out->_amplify_ratio);
			output_distribution(_mesh_out->distribution);
		}
	}
	else if (type=="Final")
	{
		cout<<"Output the final grid for body: "<<_id<<endl;
		output_tecplot(_mesh_out->output_final,_mesh_out->_amplify_ratio);
	}
	else
	{
		cout<<"Error: Invalid type in output mesh"<<endl;
	}
	return;
}
void Tbody_base::output_global_curve()
{
	if (_global_output->_is_output)
	{
		_global_output->output<<setw(20)<<_current_time<<" ";
		if (_global_output->_internal_energy)
		{
			_global_output->output<<setw(20)<<calculate_total_internal_energy()<<" ";
		}
		if (_global_output->_epeqv)
		{
			_global_output->output<<setw(20)<<calculate_maximal_effective_plastic_strain()<<" ";
		}
		if (_global_output->_max_p)
		{
			_global_output->output<<setw(20)<<calculate_maximal_pressure()<<" ";
		}
		_global_output->output<<endl;
	}
	return;
}
void Tbody_base::exam_av_type()
{
	if (_Arti_Vis==NULL)
	{
		cout<<"Error: No artificial viscosity in body "<<_id<<endl;
		system("Pause");
		exit(0);
	}
	if (_body_name=="body_shell")
	{
		if (_Arti_Vis->G_name()!="Solid_viscosity")
		{
			cout<<"Error: Body shell can not use type "<<_Arti_Vis->G_name()<<endl;
			system("Pause");
			exit(0);
		}
	}
	else if (_body_name=="body_pure_fluid")
	{
		if (_Arti_Vis->G_name()=="Solid_viscosity")
		{
			cout<<"Error: Body pure fluid can not use type "<<_Arti_Vis->G_name()<<endl;
			system("Pause");
			exit(0);
		}
	}
	return;
}
string Tbody_base::calculate_body_type_id()
{
	//Return the body type id,ordered by programming time
	if (_body_name=="body_shell")
	{
		return "00";
	}
	else if (_body_name=="body_Lagrangian")
	{
		return "01";
	}
	else if (_body_name=="body_ALE")
	{
		return "02";
	}
	else
	{
		cout<<"Error:Invalid body name in calling calculate_body_type_id()"<<endl;
		system("Pause");
		exit(0);
	}
}
//--------------------------------------------------------------------------------
//Test functions
//--------------------------------------------------------------------------------
void Tbody_base::test_for_MoF_derivative()
{
	TMoF_solver MoF_solver(2,0,3);
	vec3D p1(0,0,0),p2(0,1,0),p3(0,0,1),p4(1,0,0);
	Tpolyhedron poly;
	poly.add_vertex(p1);poly.add_vertex(p2);poly.add_vertex(p3);poly.add_vertex(p4);
	poly.add_piece(0,1,3);poly.add_piece(1,2,3);poly.add_piece(0,2,1);poly.add_piece(0,3,2);
	double fraction[2];
	vec3D centroid[2];
	Tpolyhedron* sub_poly[2];
	fraction[0]=0.3;fraction[1]=1-fraction[0];
	centroid[0].value(0.1,0.2,0.3);
	sub_poly[0]=sub_poly[1]=NULL;
	MoF_solver.MoF_for_pairs(&poly,fraction,centroid,sub_poly);
	return;
}
void Tbody_base::test_for_MoF_patch_test()
{
	TMoF_solver MoF_solver(2,0,3);
	vec3D p1(0,0,0),p2(0,1,0),p3(0,0,1),p4(1,0,0);
	Tpolyhedron poly;
	poly.add_vertex(p1);poly.add_vertex(p2);poly.add_vertex(p3);poly.add_vertex(p4);
	poly.add_piece(0,1,3);poly.add_piece(1,2,3);poly.add_piece(0,2,1);poly.add_piece(0,3,2);
	double fraction=0.3;
	MoF_solver.patch_test(&poly,fraction);
	return;
}
void Tbody_base::test_for_MoF_patch_test_3_material()
{
	TMoF_solver MoF_solver(3,0,3);
	vec3D p1(0,0,0),p2(0,1,0),p3(0,0,1),p4(1,0,0);
	Tpolyhedron poly;
	poly.add_vertex(p1);poly.add_vertex(p2);poly.add_vertex(p3);poly.add_vertex(p4);
	poly.add_piece(0,1,3);poly.add_piece(1,2,3);poly.add_piece(0,2,1);poly.add_piece(0,3,2);
	double fraction[2];
	vec2D angle[2];
	fraction[0]=fraction[1]=0.3;
	double pi=3.141592654;
	angle[0].value(0,0);angle[1].value(0.3*pi,0.7*pi);
	MoF_solver.three_material_patch_test(&poly,fraction,angle);
	return;
}