#include "Tbody_Lagrangian.h"
#include "public_function.h"
#include <iomanip>
#include <Windows.h>
using namespace std;
Tbody_Lagrangian::Smaterial_plot::Smaterial_plot(int body_id,int mid,string output_type):_mid(mid),_output_type(output_type)
{
	char a_mid[32],a_body_id[32];
	string s_mid,s_body_id;
	_itoa_s(_mid,a_mid,10);
	_itoa_s(body_id,a_body_id,10);
	s_mid=a_mid;s_body_id=a_body_id;
	if (_output_type=="surface")
	{
		_file_name="Material surface of material "+s_mid+" in body "+s_body_id+".dat";
	}
	else if (_output_type=="polyhedron")
	{
		_file_name="Material polyhedron of material "+s_mid+" in body "+s_body_id+".dat";
	}
	else if (_output_type=="all")
	{
		_file_name="All material polyhedron in body "+s_body_id+".dat";
	}
	else
	{
		cout<<"Error: Invalid material output type"<<endl;
		system("Pause");
		exit(0);
	}
	_output.open(_file_name);
}
Tbody_Lagrangian::Tbody_Lagrangian(int i,int num_material,string body_name):Tbody_base(i,body_name),_num_material(num_material){}
void Tbody_Lagrangian::calculate_force_in_body()
{
	//Calcualte the corner force for every cell
	calculate_corner_force(_dt);
	//Calculate the final nodal force
	assemble_nodal_force();
	return;
}
void Tbody_Lagrangian::S_material_plot(int mid,string type)
{
	Smaterial_plot* new_plot=NULL;
	new_plot= new Smaterial_plot(_id,mid,type);
	_material_output.push_back(new_plot);
	return;
}
double Tbody_Lagrangian::calculate_time_step()
{
	double min_time_step=100000;
	int cell_id;
	double cell_time_step;
	for (int i = 0; i < _nume; i++)
	{
		cell_time_step=_cellptr_list[i]->calculate_time_step();
		if (cell_time_step<min_time_step)
		{
			min_time_step=cell_time_step;
			cell_id=i;
		}
		//min_time_step=minval(min_time_step,_cellptr_list[i]->calculate_time_step());
	}
	return _CFL*min_time_step;
	if (_current_step==0)
	{
		return _CFL*min_time_step;
	}
	else
	{
		return minval(_CFL*min_time_step,1.05*_dt);
	}
}
bool Tbody_Lagrangian::predictor_step()
{
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->set_velocity_bc(_current_time,_dimension);
	}
	bool negative;
	//negative=volume_check("Predict");
	//if (negative)
	//{
	//	return negative;
	//}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->predict_position(_dt);
	}
	calculate_force_in_body();
	return false;
}
void Tbody_Lagrangian::calculate_final_acceleration()
{
	for (int i = 0; i < _nump; i++)
	{
		//Calculate the final nodal acceleration
		_nodeptr_list[i]->calculate_acceleration();
	}
	return;
}
bool Tbody_Lagrangian::volume_check(double dt)
{
	#pragma omp parallel for
	for (int i = 0; i < _nump; i++)
	{
		//Calculate the trial velocity at the next time step
		_nodeptr_list[i]->trial_velocity(_dt);
		//Calculate the trial position at the next time point
		_nodeptr_list[i]->update_pos(_dt);
	}
	//Check whether the volume is validated
	bool negative;
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]->calculate_cell_volume("Full",negative);
		if (negative)
		{
			for (int i = 0; i < _nump; i++)
			{
				_nodeptr_list[i]->back_to_origin_pos(_dt);
			}
			return negative;
		}
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->back_to_origin_pos(_dt);
	}
	return false;
}
bool Tbody_Lagrangian::corrector_step()
{
	#pragma omp parallel for
	for (int i = 0; i < _nump; i++)
	{
		//Calculate the velocity at the next time step
		_nodeptr_list[i]->update_vel(_dt);
		//Calculate the position at the next time point
		_nodeptr_list[i]->update_pos(_dt);
	}
	#pragma omp parallel for
	//Update cell geometry and material state
	for (int i = 0; i < _nume; i++)
	{
		if (_cellptr_list[i]->update_state(_dt))
		{
			cout<<"Error:Negative volume in the corrector step"<<endl;
			cout<<"body id= "<<_id<<" "<<"cell id= "<<_cellptr_list[i]->_id<<endl;
			system("Pause");
			exit(0);
		}
	}
	#pragma omp parallel for
	//Update the nodal variables
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->calculate_nodal_variables();
	}
	return false;
}
void Tbody_Lagrangian::calculate_nodal_inertance()
{
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->_mass=_nodeptr_list[i]->_mass0=0;
	}
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]->assemble_corner_inertance();
	}
	return;
}
void Tbody_Lagrangian::calculate_corner_force(double dt)
{
	Tviscosity_base* Vis_ptr=_Arti_Vis;
	if (_hourglass_option)
	{
		#pragma omp parallel for firstprivate(Vis_ptr)
		for (int i = 0; i < _nume; i++)
		{
			_cellptr_list[i]->calculate_corner_force(Vis_ptr,_dt);
		}
	}
	else
	{
		#pragma omp parallel for firstprivate(Vis_ptr)
		for (int i = 0; i < _nume; i++)
		{
			_cellptr_list[i]->calculate_corner_force_without_hourglass(Vis_ptr,_dt);
		}
	}

}
void Tbody_Lagrangian::assemble_nodal_force()
{
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]->assemble_corner_force();
	}
	//Apply the gravity
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->apply_gravity(_gravity);
	}
	return;
}
void Tbody_Lagrangian::set_grid_topology()
{
	//Calculate cells connected to a node
	for (int i=0;i<_nume;i++)
	{
		_cellptr_list[i]->add_to_node();
	}
	//Calculate _adjacent[] and _around
	for (int i=0;i<_nume;i++)
	{
		_cellptr_list[i]->set_cell_connection();
	}
	return;
}
double Tbody_Lagrangian::calculate_total_internal_energy()
{
	double total_internal_energy=0;
	for (int i = 0; i < _nume; i++)
	{
		for (int j = 0; j < _num_material; j++)
		{
			total_internal_energy=total_internal_energy+_cellptr_list[i]->calculate_cell_internal_energy();
		}
	}
	return total_internal_energy;
}
double Tbody_Lagrangian::calculate_maximal_pressure()
{
	double max_pressure=0;
	for (int i = 0; i < _nume; i++)
	{
		double cell_pressuer=_cellptr_list[i]->_gausspoint[0]->G_pressure();
		max_pressure=maxval(max_pressure,cell_pressuer);
	}
	return max_pressure;
}
void Tbody_Lagrangian::calculate_body_size(vec3D &coor_min,vec3D &coor_max,double &cell_edge_max,double &cell_edge_min)
{
	cell_edge_max=0,cell_edge_min=1000000;
	double l_ele;
	vec3D eps(1e-5,1e-5,1e-5);
	for (int i=0;i<_nume;i++)
	{
		l_ele=_cellptr_list[i]->calculate_characteristic_length();
		if (l_ele>cell_edge_max)
		{
			cell_edge_max=l_ele;
		}
		if (l_ele<cell_edge_min)
		{
			cell_edge_min=l_ele;
		}
	}
	coor_min.value(1e10,1e10,1e10);coor_max.value(-1e10,-1e10,-1e10);
	for (int i=0;i<_nump;i++)
	{
		for (int j = 1; j < 4; j++)
		{
			if (_nodeptr_list[i]->_position.access(j)<coor_min.access(j)) coor_min.value(j,_nodeptr_list[i]->_position.access(j));
			if (_nodeptr_list[i]->_position.access(j)>coor_max.access(j)) coor_max.value(j,_nodeptr_list[i]->_position.access(j));
		}
	}
	//To avoid a node is exactly located on the boundary
	coor_min=coor_min-eps;coor_max=coor_max+eps;
	return;
}
void Tbody_Lagrangian::surface_reconstrcution()
{
	if (_num_material==1)
	{
		return;
	}
	//Surface reconstruction only applied on multi-material problem
	TMoF_solver MoF_solver(_num_material,0,_dimension);
	Tcell_fluid_base* cell_ptr;
	double* fraction_reference;vec3D* centroid_reference;
	fraction_reference=new double[_num_material];
	centroid_reference=new vec3D[_num_material];
	double t_begin=GetTickCount();
	for (int i = 0; i < _nume; i++)
	{
		int pre=MoF_solver.G_nun_non_convergence();
		cell_ptr=_cellptr_list[i];
		cell_ptr->update_cell_polyhedron();
		cell_ptr->update_material_centroid();
		for (int j = 0; j < _num_material; j++)
		{
			fraction_reference[j]=cell_ptr->G_material_fraction(j);
			centroid_reference[j]=cell_ptr->G_material_centroid(j);
		}
		MoF_solver.set_initial_condition(&cell_ptr->_cell_polyhedron,fraction_reference,centroid_reference,_num_material);
		MoF_solver.MoF_surface_reconstrcution();
		for (int j = 0; j < _num_material; j++)
		{
			cell_ptr->S_material_polyhedron(j,MoF_solver.G_result(j));
		}
		MoF_solver.Finalize();
	}
	double t_end=GetTickCount();
	cout<<"Time consuming in surface reconstruction is "<<t_end-t_begin<<endl;
	cout<<"Number of non-convergence BFGS interation is "<<MoF_solver.G_nun_non_convergence()<<endl;
	return;
}
void Tbody_Lagrangian::reset_final_initial_condition()
{
	if (_num_material==1)
	{
		return;
	}
	else
	{
		for (int i = 0; i < _nume; i++)
		{
			_cellptr_list[i]->reset_geometry_initial_condition();
		}
		if (_test_name=="no")
		{
			for (int i = 0; i < _nume; i++)
			{
				_cellptr_list[i]->reset_conservative_variable_initial_condition();
			}
		}
		else if (_test_name=="2D_Taylor")
		{
			for (int i = 0; i < _nume; i++)
			{
				_cellptr_list[i]->reset_TR_conditions("2D");
				_cellptr_list[i]->reset_conservative_variable_initial_condition();
			}
		}
		else if (_test_name=="3D_Taylor")
		{
			for (int i = 0; i < _nume; i++)
			{
				_cellptr_list[i]->reset_TR_conditions("3D");
				_cellptr_list[i]->reset_conservative_variable_initial_condition();
			}
		}
		else
		{
			cout<<"Invalid standard test name!"<<endl;
			system("Pause");
			exit(0);
		}
	}
	return;
}
void Tbody_Lagrangian::standard_test_initial_conditions()
{
	if (_test_name=="no")
	{
		return;
	}
	else if (_test_name=="2D_Taylor")
	{
		for (int i = 0; i < _nume; i++)
		{
			_cellptr_list[i]->reset_TR_conditions("2D");
		}
		for (int i = 0; i < _nump; i++)
		{
			_nodeptr_list[i]->calculate_nodal_variables();
		}
		_dimension=2;
	}
	else if (_test_name=="3D_Taylor")
	{
		for (int i = 0; i < _nume; i++)
		{
			_cellptr_list[i]->reset_TR_conditions("3D");
		}
		for (int i = 0; i < _nump; i++)
		{
			_nodeptr_list[i]->calculate_nodal_variables();
		}
	}
	else if (_test_name=="dam_break")
	{
		for (int i = 0; i < _nume; i++)
		{
			_cellptr_list[i]->reset_dam_breaking_conditions(0.02,30,1,_num_material);
		}
		calculate_nodal_inertance();
		for (int i = 0; i < _nump; i++)
		{
			_nodeptr_list[i]->calculate_nodal_variables();
		}
	}
	else if (_test_name=="Sedov")
	{
		TMAT_base* mat_ptr=_cellptr_list[0]->G_gausspoint(0);
		mat_ptr->S_internal_energy(6808);
		double pressure=mat_ptr->calculate_pressure();
		mat_ptr->S_pressure(pressure);
		bool fail;
		double c=mat_ptr->calculate_sound_speed(fail);
		mat_ptr->S_soundspeed(c);
		cout<<c<<endl;
	}
	else
	{
		cout<<"Invalid standard test name!"<<endl;
		system("Pause");
		exit(0);
	}
}
bool Tbody_Lagrangian::volume_check(string type)
{
	bool negative;
	if (type=="Predict")
	{
		for (int i = 0; i < _nump; i++)
		{
			_nodeptr_list[i]->predict_position(_dt);
		}
		//If any temporary volume is negative,perform the remapping
		for (int i = 0; i < _nume; i++)
		{
			_cellptr_list[i]->calculate_cell_volume("Temporary",negative);
			if (negative)
			{
				return true;
			}
		}
	}
	else if (type=="Correct")
	{
		bool negative;
		#pragma omp parallel for
		for (int i = 0; i < _nump; i++)
		{
			//Calculate the velocity at the next time step
			_nodeptr_list[i]->trial_velocity(_dt);
			//Calculate the position at the next time point
			_nodeptr_list[i]->trial_pos(_dt);
		}
		//Check the volume before the update
		for (int i = 0; i < _nume; i++)
		{
			_cellptr_list[i]->calculate_cell_volume("Trial",negative);
			if (negative)
			{
				cout<<i<<endl;
				return true;
			}
		}
	}
	else
	{
		cout<<"Invalid type"<<endl;
		system("Pause");
	}
	return false;
}
//--------------------------------------------------------------------------
//Output functions
//--------------------------------------------------------------------------
void Tbody_Lagrangian::plot_reconstructed_material()
{
	string type;
	int mid;
	int num_material_plot=_material_output.size();
	for (int i = 0; i < num_material_plot; i++)
	{
		mid=_material_output[i]->_mid;
		type=_material_output[i]->_output_type;
		if (type=="surface")
		{
			plot_material_surface(mid,_material_output[i]->_output);
		}
		else if (type=="polyhedron")
		{
			plot_material_polyhedron(mid,_material_output[i]->_output);
		}
		else if (type=="all")
		{
			plot_all_material_polyhedron(_material_output[i]->_output);
		}
	}
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]->delete_material_polyhedron();
	}
	return;
}
void Tbody_Lagrangian::output_tecplot(ofstream &output,double ratio)
{
	output<<"TITLE = \"Mesh for Brick fluid Element\""<<endl;
	output<<"VARIABLES = \"X\",\"Y\",\"Z\",\"density\",\"pressure\",\"vx\",\"vy\",\"vz\",\"vmag\"£¬\"soundspeed\""<<endl;	
	//output<<"VARIABLES = \"X\",\"Y\",\"Z\""<<endl;	
	output<<"ZONE T=\""<<_current_time<<"\" F=FEPOINT,N="<<_nump<<","<<"E="<<_nume<<","<<"ET=BRICK"<<endl;
	vec3D dis;
	for (int i = 0; i < _nump; i++)
	{
		ratio=1;
		dis=_nodeptr_list[i]->calculate_displacement_amplify(ratio);
		output<<dis.x<<" "<<dis.y<<" "<<dis.z<<" "<<_nodeptr_list[i]->_density<<" "<<_nodeptr_list[i]->_pressure
			  <<" "<<_nodeptr_list[i]->_velocity.x<<" "<<_nodeptr_list[i]->_velocity.y<<" "<<_nodeptr_list[i]->_velocity.z
			  <<" "<<_nodeptr_list[i]->_velocity.get_length()<<" "<<_nodeptr_list[i]->_mass<<endl;
		//output<<dis.x<<" "<<dis.y<<" "<<dis.z<<endl;
	}
	for (int i = 0; i < _nume; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			output<<_cellptr_list[i]->_node_ptr[j]->_id<<" ";
		}
		output<<endl;
	}
	ofstream curve;
	curve.open("curve.txt");
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]->update_cell_polyhedron();
		double r=_cellptr_list[i]->G_cell_polyhedron()->G_centroid().get_length();
		//double den=0;
		//for (int j = 0; j < _num_material; j++)
		//{
			//den=den+_cellptr_list[i]->_gausspoint[j]->G_density()*_cellptr_list[i]->G_material_fraction(j);
		//}
		//curve<<r<<" "<<den<<endl;
		curve<<r<<" "<<_cellptr_list[i]->_gausspoint[0]->G_density()<<endl;
	}
	curve.close();
	surface_reconstrcution();
	plot_reconstructed_material();
	
	return;
}
void Tbody_Lagrangian::output_distribution(ofstream &output)
{
	vec3D dis;
	output<<"VARIABLES = \"X\",\"p\""<<endl;
	for (int i = 0; i < _nump; i++)
	{
		dis=_nodeptr_list[i]->_position;
		if (abs(dis.y)<1e-10 && abs(dis.z)<1e-10)
		{
			output<<dis.x<<" "<<_nodeptr_list[i]->_pressure<<endl;
		}
	}
	return;
}
void Tbody_Lagrangian::output_curve()
{
	Tbody_base::output_curve();
	int num_curve=_curve_out_list.size();
	int type,id;
	Tnode_fluid* node_ptr;
	Tcell_fluid_base* cell_ptr;
	for (int i = 0; i < num_curve; i++)
	{
		type=_curve_out_list[i]->_type;
		id=_curve_out_list[i]->_id-1;
		_curve_out_list[i]->output.precision(10);
		if (type==0)
		{
			node_ptr=_nodeptr_list[id];
			vec3D dis=node_ptr->calculate_displacement();
			vec3D vel=node_ptr->_velocity;
			_curve_out_list[i]->output<<setw(20)<<_current_time<<" ";
			_curve_out_list[i]->output<<setw(20)<<vel.x<<setw(20)<<vel.y<<setw(20)<<vel.z<<endl;
		}
		else if (type==1)
		{
			cell_ptr=_cellptr_list[id];
			_curve_out_list[i]->output<<setw(20)<<_current_time<<" ";
			_curve_out_list[i]->output<<setw(20)<<cell_ptr->_gausspoint[0]->G_pressure()<<endl;
		}
	}
	return;
}
void Tbody_Lagrangian::output_cell_monitor()
{
	int num_monitor=_cell_monitor_list.size();
	int cell_id;
	for (int i = 0; i < num_monitor; i++)
	{
		cell_id=_cell_monitor_list[i]->_cell_id-1;
		_cell_monitor_list[i]->output.precision(10);
		_cellptr_list[cell_id]->output_cell_monitor(_cell_monitor_list[i]->output,_current_time);
	}
	return;
}
void Tbody_Lagrangian::plot_material_surface(int mid,ofstream &output)
{
	int numv=0,nump=0;                 //Total vertexes number and total surface triangle pieces number
	Tcell_fluid_base* cell_ptr=NULL;
	vec3D vertex;
	int* IEN;
	int vertex_start;                 //Start vertex in every polyhedron
	int** piece_list;                 //Total pieces list
	vec3D* vertex_list;               //Total vertexes list
	int* global2used;                 //The id transformation from total vertex to used vertex
	int* used2global;                 //The id transformation from used vertex to total vertex
	
	//Calculate the total vertexes number and total surface pieces number before compression
	for (int i=0;i<_nume;i++)
	{
		cell_ptr=_cellptr_list[i];
		if (cell_ptr->G_material_fraction(mid)>0 && cell_ptr->G_material_fraction(mid)!=1)
		{
			//For mixed cell,the surface piece is the pieces on cutting plane
			numv=numv+cell_ptr->G_material_polyhedron(mid)->G_numv();
			nump=nump+cell_ptr->G_material_polyhedron(mid)->G_num_surface_piece();
		}
		else if (cell_ptr->G_material_fraction(mid)==1)
		{
			//For pure cell,the surface is the cell face which is adjacent to a cell which does not contain material mid
			numv=numv+cell_ptr->G_cell_polyhedron()->G_numv();
			nump=nump+cell_ptr->is_surface_cell(mid);
		}
	}
	//Allocate the memory for total vertex and piece and their topology
	piece_list=new int*[nump];
	for (int i=0;i<nump;i++)
	{
		piece_list[i]=new int[3];
	}
	vertex_list=new vec3D[numv];
	global2used=new int[numv];
	for (int i=0;i<numv;i++)
	{
		global2used[i]=-1;
	}
	//Get the total vertexes' coordinate and pieces' IEN before vertex compression
	numv=0;
	for (int i=0;i<_nume;i++)
	{
		cell_ptr=_cellptr_list[i];
		if (cell_ptr->G_material_fraction(mid)>0)
		{
			for (int j=0;j<cell_ptr->G_material_polyhedron(mid)->G_numv();j++)
			{
				vertex=cell_ptr->G_material_polyhedron(mid)->G_vertex(j);
				vertex_list[numv]=vertex;
				numv=numv+1;
			}	
		}
	}
	vertex_start=1;
	nump=0;
	for (int i=0;i<_nume;i++)
	{		
		cell_ptr=_cellptr_list[i];
		if (cell_ptr->G_material_fraction(mid)>0 && cell_ptr->G_material_fraction(mid)!=1)
		{
			//The cutting plane in mixed cell
			for (int j=0;j<cell_ptr->G_material_polyhedron(mid)->G_num_surface_piece();j++)
			{
				int k=cell_ptr->G_material_polyhedron(mid)->G_surface_piece_id(j);
				IEN=cell_ptr->G_material_polyhedron(mid)->G_IEN(k);
				for (int i=0;i<3;i++)
				{
					piece_list[nump][i]=IEN[i]+vertex_start;
				}
				nump=nump+1;
			}
			vertex_start=vertex_start+cell_ptr->G_material_polyhedron(mid)->G_numv();
		}
		else if (cell_ptr->G_material_fraction(mid)==1)
		{
			cell_ptr->add_surface_cell_surface(mid,piece_list,vertex_start,nump);
			vertex_start=vertex_start+cell_ptr->G_cell_polyhedron()->G_numv();
		}
	}
	//Compress the vertex:only plot the used vertexes
	int global_id;
	int num_used=0;
	for (int i=0;i<nump;i++)
	{
		for (int j=0;j<3;j++)
		{
			global_id=piece_list[i][j]-1;
			if (global2used[global_id]==-1)
			{
				num_used=num_used+1;
				global2used[global_id]=num_used;
			}
		}
	}
	used2global=new int[num_used+1];
	for (int i=0;i<numv;i++)
	{
		if (global2used[i]>0)
		{
			used2global[global2used[i]]=i;
		}
	}
	output<<"TITLE=\"material polyhedron\""<<endl;
	output<<"Variables=\"X\",\"Y\",\"Z\""<<endl;
	//output<<"ZONE F=FEPOINT N= "<<num_used<<" E= "<<nump<<" ET=TRIANGLE"<<endl;
	output<<"ZONE T=\""<<_current_time<<"\" F=FEPOINT N= "<<num_used<<" E= "<<nump<<" ET=TRIANGLE"<<endl;
	for (int i=1;i<num_used+1;i++)
	{
		output<<vertex_list[used2global[i]].x<<" "<<vertex_list[used2global[i]].y<<" "<<vertex_list[used2global[i]].z<<endl;
	}
	for (int i=0;i<nump;i++)
	{
		output<<global2used[piece_list[i][0]-1]<<" "<<global2used[piece_list[i][1]-1]<<" "<<global2used[piece_list[i][2]-1]<<endl;
	}
	//Recover the memory
	delete used2global;
	delete global2used;
	delete vertex_list;
	for (int i=0;i<nump;i++)
	{
		delete piece_list[i];
	}
	delete piece_list;
	return;
}
void Tbody_Lagrangian::plot_material_polyhedron(int mid,ofstream &output)
{
	int numv=0,nump=0;          //Total point and total triangle piece
	Tcell_fluid_base* cell_ptr=NULL;
	vec3D vertex;
	int* IEN;
	int vertex_start;          //Start vertex in every polyhedron
	for (int i=0;i<_nume;i++)
	{
		cell_ptr=_cellptr_list[i];
		if (cell_ptr->G_material_fraction(mid)>0)
		{
			numv=numv+cell_ptr->G_material_polyhedron(mid)->G_numv();
			nump=nump+cell_ptr->G_material_polyhedron(mid)->G_nump();
		}
	}
	output<<"TITLE=\"material polyhedron\""<<endl;
	output<<"Variables=\"X\",\"Y\",\"Z\",\"mid\""<<endl;
	//output<<"ZONE F=FEPOINT N= "<<numv<<" E= "<<nump<<" ET=TRIANGLE"<<endl;
	output<<"ZONE T=\""<<_current_time<<"\" F=FEPOINT N= "<<numv<<" E= "<<nump<<" ET=TRIANGLE"<<endl;
	//Output the vertexes
	for (int i=0;i<_nume;i++)
	{
		cell_ptr=_cellptr_list[i];
		if (cell_ptr->G_material_fraction(mid)>0)
		{
			//output<<i<<endl;
			for (int j=0;j<cell_ptr->G_material_polyhedron(mid)->G_numv();j++)
			{
				vertex=cell_ptr->G_material_polyhedron(mid)->G_vertex(j);
				output<<vertex.x<<" "<<vertex.y<<" "<<vertex.z<<" "<<mid<<endl;
			}	
		}
	}
	//Output the pieces
	vertex_start=1;
	for (int i=0;i<_nume;i++)
	{
		cell_ptr=_cellptr_list[i];
		if (cell_ptr->G_material_fraction(mid)>0)
		{
			for (int j=0;j<cell_ptr->G_material_polyhedron(mid)->G_nump();j++)
			{
				IEN=cell_ptr->G_material_polyhedron(mid)->G_IEN(j);
				output<<IEN[0]+vertex_start<<" "<<IEN[1]+vertex_start<<" "<<IEN[2]+vertex_start<<" "<<endl;
			}
			vertex_start=vertex_start+cell_ptr->G_material_polyhedron(mid)->G_numv();
		}
	}
	return;
}
void Tbody_Lagrangian::plot_all_material_polyhedron(ofstream &output)
{
	int numv=0,nump=0;          //Total point and total triangle piece
	Tcell_fluid_base* cell_ptr=NULL;
	vec3D vertex;
	int* IEN;
	int vertex_start;          //Start vertex in every polyhedron
	for (int i=0;i<_nume;i++)
	{
		cell_ptr=_cellptr_list[i];
		for (int mid = 0; mid < _num_material; mid++)
		{
			if (cell_ptr->G_material_fraction(mid)>0)
			{
				numv=numv+cell_ptr->G_material_polyhedron(mid)->G_numv();
				nump=nump+cell_ptr->G_material_polyhedron(mid)->G_nump();
			}
		}
	}
	output<<"TITLE=\"material polyhedron\""<<endl;
	output<<"Variables=\"X\",\"Y\",\"Z\",\"mid\""<<endl;
	//output<<"ZONE F=FEPOINT N= "<<numv<<" E= "<<nump<<" ET=TRIANGLE"<<endl;
	output<<"ZONE T=\""<<_current_time<<"\" F=FEPOINT N= "<<numv<<" E= "<<nump<<" ET=TRIANGLE"<<endl;
	//Output the vertexes
	for (int i=0;i<_nume;i++)
	{
		cell_ptr=_cellptr_list[i];
		for (int mid = 0; mid < _num_material; mid++)
		{
			if (cell_ptr->G_material_fraction(mid)>0)
			{
				for (int j=0;j<cell_ptr->G_material_polyhedron(mid)->G_numv();j++)
				{
					vertex=cell_ptr->G_material_polyhedron(mid)->G_vertex(j);
					output<<vertex.x<<" "<<vertex.y<<" "<<vertex.z<<" "<<mid<<endl;
				}	
			}
		}
	}
	//Output the pieces
	vertex_start=1;
	for (int i=0;i<_nume;i++)
	{
		cell_ptr=_cellptr_list[i];
		for (int mid = 0; mid < _num_material; mid++)
		{
			if (cell_ptr->G_material_fraction(mid)>0)
			{
				for (int j=0;j<cell_ptr->G_material_polyhedron(mid)->G_nump();j++)
				{
					IEN=cell_ptr->G_material_polyhedron(mid)->G_IEN(j);
					output<<IEN[0]+vertex_start<<" "<<IEN[1]+vertex_start<<" "<<IEN[2]+vertex_start<<" "<<endl;
				}
				vertex_start=vertex_start+cell_ptr->G_material_polyhedron(mid)->G_numv();
			}
		}
	}
	return;
}

void Tbody_Lagrangian::input_body(ifstream& input)
{
	Skeyword keyword;
	read_in_keyword_file(input,keyword);
	//Read the body global information from keyword
	_nume=keyword.cell_8_list.size();
	_nump=keyword.node_list.size();
	_endtime=keyword.time_control.endtime;
	_CFL=keyword.time_control.CFL;
	_hourglass_option=keyword.hourglass_option;
	//----------------------------------------------------
	//Read node list
	//----------------------------------------------------
	int id;
	double x,y,z;
	if (_num_material==1)
	{
		_grid_pure_fluid._node_list.resize(_nump);
	}
	else
	{
		_grid_mix_fluid._node_list.resize(_nump);
	}
	_nodeptr_list.resize(_nump);
	for (int i = 0; i < _nump; i++)
	{
		id=keyword.node_list[i].id;
		x=keyword.node_list[i].x;
		y=keyword.node_list[i].y;
		z=keyword.node_list[i].z;
		Tnode_fluid new_node(id,x,y,z);
		if (_num_material==1)
		{
			_grid_pure_fluid._node_list[i]=new_node;
		}
		else
		{
			_grid_mix_fluid._node_list[i]=new_node;
		}
		_nodeptr_list[i]=new Tnode_fluid(id,x,y,z);
	}
	//-----------------------------------------------------
	//Read cell list
	//-----------------------------------------------------
	if (_num_material==1)
	{
		_grid_pure_fluid._cell_list.resize(_nume);
	}
	else
	{
		_grid_mix_fluid._cell_list.resize(_nume);
	}
	_cellptr_list.resize(_nume);
	int part_id,material_id,EOS_id;
	Tnode_fluid* node_ptr[8];
	for (int i = 0; i < _nume; i++)
	{
		//Read cell id and IEN
		id=keyword.cell_8_list[i].cell_id;
		for (int j = 0; j < 8; j++)
		{
			int IENj=keyword.cell_8_list[i].IEN[j];
			node_ptr[j]=_nodeptr_list[IENj-1];
		}
		//Read cell properties (nGauss and thickness)
		part_id=keyword.cell_8_list[i].part_id;
		//Read cell material
		material_id=keyword.part_list[part_id-1].material_id;
		EOS_id=keyword.part_list[part_id-1].EOS_id;
		if (EOS_id==0)
		{
			cout<<"Error:There must be an EOS in body_Lagrangian"<<endl;
			system("Pause");
			exit(0);
		}
		TMAT_base* mat;
		Smaterial this_mat=keyword.material_list[material_id-1];
		SEOS this_EOS=keyword.EOS_list[EOS_id-1];
		mat=generate_material(this_mat,this_EOS,1);
		if (_num_material==1)
		{
			_cellptr_list[i]=new Tcell_pure_fluid(id,node_ptr,mat);
		}
	}
	//-------------------------------------------------------
	//Read boundary condition
	//-------------------------------------------------------
	int node_group_id;   //Boundary condition (load) is applied on node_group_id

	int num_bd_group;    //Number of boundary condition group
	int num_node;        //Number of nodes in the node group
	num_bd_group=keyword.boundary_list.size();
	int bc_type;
	for (int i = 0; i < num_bd_group; i++)
	{
		node_group_id=keyword.boundary_list[i].id;
		num_node=keyword.node_group_list[node_group_id-1].node_id.size();
		bc_type=keyword.boundary_list[i].type;
		if (bc_type==0)
		{
			for (int j = 0; j < 3; j++)
			{
				if (keyword.boundary_list[i].pos[j]==1)
				{
					for (int k = 0; k < num_node; k++)
					{
						_nodeptr_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_type_position[j]=1;
					}
				}
			}
		}
		else if (bc_type==1)
		{
			for (int j = 0; j < 3; j++)
			{
				if (keyword.boundary_list[i].dofvel[j]==1)
				{
					for (int k = 0; k < num_node; k++)
					{
						_nodeptr_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_type_position[j]=2;
						_nodeptr_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_velocity[j]=keyword.boundary_list[i].vel[j];
					}
				}
			}
		}
		else if (bc_type==2)
		{
			vec2D center(keyword.boundary_list[i].x,keyword.boundary_list[i].y);
			double v=keyword.boundary_list[i].v;
			for (int k = 0; k < num_node; k++)
			{
				_nodeptr_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_type_position[0]=2;
				_nodeptr_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_type_position[1]=2;
				_nodeptr_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_type_position[2]=2;
				vec3D coor=_nodeptr_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_position;
				vec2D r(coor.x-center.x,coor.y-center.y);
				r=r/r.get_length();
				_nodeptr_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_velocity[0]=v*r.x;
				_nodeptr_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_velocity[1]=v*r.y;
				_nodeptr_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_velocity[2]=0;
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
	num_load_group=keyword.load_list.size();
	for (int i = 0; i < num_load_group; i++)
	{
		node_group_id=keyword.load_list[i].id;
		new_load.node_id=keyword.node_group_list[node_group_id-1].node_id;
		new_load.num_node=new_load.node_id.size();
		direction=keyword.load_list[i].dof;
		new_load.direction=direction;
		curve_id=keyword.load_list[i].curve_id;
		magnitude=keyword.curve_list[curve_id-1].v2[0];      //Only for constant load
		new_load.magnitude=magnitude;
		_external_load_list.push_back(new_load);
	}
	
	set_grid_topology();
	//Get the cell connection

	//Initialize the nodal variables
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->calculate_nodal_variables();
	}
	return;
}