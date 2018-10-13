#include "Tregion.h"
#include "Tbody_shell.h"
#include "Tbody_ALE.h"
#include "Tbody_Lagrangian.h"
#include "Tbody_explosion.h"
#include "Tbody_preprocessor.h"
#include "public_function.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "readin.h"
#include <iomanip>
#include "Tbody_ALEMPM.h"
#include "Tbody_MPM.h"
//#include "Tmaterial_grid.h"
using namespace std;
Tregion::Tregion()
{
	_num_body=_current_step=_num_interaction=0;_time_step=_current_time=_endtime=0;
}
int Tregion::calculate_time_step()
{
	double min_step=100000;
	int id=-1;
	//Calculate the minimal time step for all bodies
	for (int i = 0; i < _num_body; i++)
	{
		double body_step=_body_list[i]->calculate_time_step();
		if (body_step<min_step)
		{
			min_step=body_step;
			id=i;
		}
	}
	if (_current_time+min_step>_endtime)
	{
		min_step=_endtime-_current_time;
	}
	//Set the minval step for all bodies
	for (int i = 0; i < _num_body; i++)
	{
		_body_list[i]->S_dt(min_step);
	}
	_time_step=min_step;
	return id;
}
void Tregion::time_integration()
{

	//Start up for solid body
	for (int i = 0; i < _num_body; i++)
	{
		_body_list[i]->start_up_step();
	}
	//Output the mesh in initial phase
	for (int i = 0; i < _num_body; i++)
	{
		_body_list[i]->output_mesh("Initial");
	}
	//Advance in time
	while (_current_time<_endtime)
	{
		//cout<<"new cycle"<<endl;
		bool is_remapping=false;
		int id;
		//Calculate the time step
		id=calculate_time_step();
		//cout<<"after calculate time step"<<endl;
		//Advance time
		_current_time=_current_time+_time_step;
		_current_step=_current_step+1;
		cout<<setw(6)<<_current_step<<" cycles at time="<<setw(8)<<_current_time<<"  dt="<<setw(8)<<_time_step<<" controlled by body "<<id<<endl;
		//Predictor step:calculate the nodal force in bodies
		for (int i = 0; i < _num_body; i++)
		{
			if (_body_list[i]->predictor_step())
			{
				cout<<"The temperory volume is negative in the predict step"<<endl;
				system("Pause");
			}
		}
		//cout<<"after predictor step"<<endl;
		//Modify the nodal force by calculating the interaction between bodies
		for (int i = 0; i < _num_interaction; i++)
		{
			_interaction_list[i]->apply_interaction_force();
		}
		//Calculate the final acceleraction
		for (int i = 0; i < _num_body; i++)
		{
			_body_list[i]->calculate_final_acceleration();
		}
		//cout<<"after calculate acceleration"<<endl;
		//Corrector step:update position,velocity and material state
		for (int i = 0; i < _num_body; i++)
		{
			bool is_remapping_in_body;
			is_remapping_in_body=_body_list[i]->corrector_step();
			if (is_remapping_in_body)
			{
				is_remapping=true;
			}
		}
		//cout<<"after corrector step"<<endl;
		if (is_remapping)
		{
			for (int i = 0; i < _num_interaction; i++)
			{
				_interaction_list[i]->reset_nodal_inertance();
			}
			for (int i = 0; i < _num_interaction; i++)
			{
				_interaction_list[i]->modify_after_remapping();
			}
			for (int i = 0; i < _num_interaction; i++)
			{
				_interaction_list[i]->exam_continuity();
			}
		}
		for (int i = 0; i < _num_interaction; i++)
		{
			_interaction_list[i]->exam_continuity();
		}
		//Advance in time for bodies
		for (int i = 0; i < _num_body; i++)
		{
			_body_list[i]->advance_in_time();
		}
		//Output the curve
		for (int i = 0; i < _num_body; i++)
		{
			_body_list[i]->output_curve();
		}
		//Output global veriables
		for (int i = 0; i < _num_body; i++)
		{
			_body_list[i]->output_global_curve();
		}
		//Output the mesh in process
		for (int i = 0; i < _num_body; i++)
		{
			_body_list[i]->output_mesh("Process");
		}
		//if (_current_step==16000 || _current_step==17200 || _current_step==20000)
		//{
		//	_body_list[0]->output_mesh("Final");
		//}
		//Output the cell monitor
		for (int i = 0; i < _num_body; i++)
		{
			_body_list[i]->output_cell_monitor();
		}
	}
	//Out put the mesh in final phase
	for (int i = 0; i < _num_body; i++)
	{
		_body_list[i]->output_mesh("Final");
	}
	return;
}
//-------------------------------------------------------------------------------
//Function for input
//-------------------------------------------------------------------------------
bool Tregion::input_region()
{
	bool state=true;
	ifstream input;
	string a;
	input.open("body_list.dat");
	if (input.fail())
	{
		cout<<"Error in openning file body_list.dat!!!"<<endl;
		system("Pause");
		exit(0);
	}
	//Read the name and type of every body file
	vector<string> name_list;
	vector<string> body_type_list;
	vector<int> materials_in_body;
	int num_exist_body=0;
	_num_body=0;
	while (true)
	{
		next_keyword(input);
		input>>a;
		if (a=="*BODY_LIST")
		{
			unsigned int body_id;
			while (!exam_keyword(input))
			{
				input>>body_id;
				if (body_id+1>name_list.size())
				{
					name_list.resize(body_id+1);
					body_type_list.resize(body_id+1);
					materials_in_body.resize(body_id+1);
					input>>name_list[body_id]>>body_type_list[body_id]>>materials_in_body[body_id];
				}
				else
				{
					input>>name_list[body_id]>>body_type_list[body_id]>>materials_in_body[body_id];
				}
				num_exist_body=num_exist_body+1;
			}
			_num_body=name_list.size();
			if (_num_body<=0)
			{
				cout<<"Error:The number of body is 0"<<endl;
				system("Pause");
				exit(0);
			}
			if (_num_body>num_exist_body)
			{
				cout<<"Error:The body id is not continuous"<<endl;
				system("Pause");
				exit(0);
			}
			else if (_num_body<num_exist_body)
			{
				cout<<"Error:Coincident body id exists"<<endl;
				system("Pause");
				exit(0);
			}
			_body_list.resize(_num_body);
		}
		else if (a=="*END")
		{
			input.close();
			break;
		}
	}
	//Exam whether all the input files exist
	for (int i = 0; i < _num_body; i++)
	{
		input.open(name_list[i]);
		if (input.fail())
		{
			cout<<"Error:File "<<name_list[i]<<" is not existed"<<endl;
			system("Pause");
			state=false;
		}
		input.close();
	}
	//If all the input files exist,valid,read every input file
	if (state)
	{
		int num_fluid_cell=0;
		int num_solid_cell=0;
		for (int i = 0; i < _num_body; i++)
		{
			input.open(name_list[i]);
			if (body_type_list[i]=="body_shell")
			{
				Tbody_shell* body_shell;
				body_shell=new Tbody_shell(i,body_type_list[i]);
				_body_list[i]=body_shell;
				_body_list[i]->input_body(input);
				num_solid_cell=num_solid_cell+_body_list[i]->_nume;
			}
			else if (body_type_list[i]=="body_Lagrangian")
			{
				if (materials_in_body[i]==1)
				{
					Tbody_Lagrangian* body_Lagrangian;
					body_Lagrangian=new Tbody_Lagrangian(i,materials_in_body[i],body_type_list[i]);
					_body_list[i]=body_Lagrangian;
					_body_list[i]->input_body(input);
					num_fluid_cell=num_fluid_cell+_body_list[i]->_nume;
				}
			}
			else if (body_type_list[i]=="body_ALE")
			{
				Tbody_ALE* body_ALE;
				body_ALE=new Tbody_ALE(i,materials_in_body[i],body_type_list[i]);
				_body_list[i]=body_ALE;
				_body_list[i]->input_body(input);
				num_fluid_cell=num_fluid_cell+_body_list[i]->_nume;
			}
			else if (body_type_list[i]=="body_MPM")
			{
				Tbody_MPM* body_MPM;
				body_MPM=new Tbody_MPM(i,materials_in_body[i],body_type_list[i]);
				_body_list[i]=body_MPM;
				_body_list[i]->input_body(input);
			}
			else if (body_type_list[i]=="body_ALEMPM")
			{
				Tbody_ALEMPM* body_ALEMPM;
				body_ALEMPM=new Tbody_ALEMPM(i,materials_in_body[i],body_type_list[i]);
				_body_list[i]=body_ALEMPM;
				_body_list[i]->input_body(input);
			}
			else if (body_type_list[i]=="body_explosion")
			{
				Tbody_explosion* body_explosion;
				body_explosion=new Tbody_explosion(i,materials_in_body[i],body_type_list[i]);
				_body_list[i]=body_explosion;
				_body_list[i]->input_body(input);
			}
			else if (body_type_list[i]=="body_preprocessor")
			{
				Tbody_preprocessor* body_preprocessor;
				body_preprocessor=new Tbody_preprocessor(i,body_type_list[i]);
				_body_list[i]=body_preprocessor;
				_body_list[i]->modify_body(input);
				cout<<"Finish the body modification"<<endl;
				system("Pause");
				exit(0);
			}
			else 
			{
				cout<<"Error:Invalid body name in function input_region()"<<endl;
				system("Pause");
				exit(0);
			}
			input.close();
		}
	//system("Pause");
	}

	//Set the endtime for the region
	for (int i = 0; i < _num_body; i++)
	{
		_endtime=maxval(_body_list[i]->_endtime,_endtime);
	}
	//Calculate the nodal inertance
	for (int i = 0; i < _num_body; i++)
	{
		_body_list[i]->calculate_nodal_inertance();
	}
	//Input the output information
	input_extra_info();
	for (int i = 0; i < _num_body; i++)
	{
		_body_list[i]->standard_test_initial_conditions();
	}
	return state;
}
void Tregion::input_extra_info()
{
	ifstream input;
	string a;
	input.open("extra_information.dat");
	if (input.fail())
	{
		cout<<"Error: No extra information input"<<endl;
		system("Pause");
		exit(0);
	}
	while (true)
	{
		next_keyword(input);
		input>>a;
		if (a=="*CURVE_NODE")
		{
			int body_id,id;
			while (!exam_keyword(input))
			{
				input>>body_id>>id;
				if (body_id>=_num_body)
				{
					cout<<"Error: Body "<<body_id<<" in *CURVE_NODE excess"<<endl;
					cout<<"_num_body="<<_num_body<<endl;
					system("Pause");
					exit(0);
				}
				if (id>=_body_list[body_id]->_nump)
				{
					cout<<"Error: Node "<<id<<" in body "<<body_id<<" excess"<<endl;
					cout<<"Number of node in body "<<body_id<<" is "<<_body_list[body_id]->_nump<<endl;
					system("Pause");
					exit(0);
				}
				_body_list[body_id]->create_curve(body_id,0,id);
			}
		}
		else if (a=="*CURVE_CELL")
		{
			int body_id,id;
			while (!exam_keyword(input))
			{
				input>>body_id>>id;
				if (body_id>=_num_body)
				{
					cout<<"Error: Body "<<body_id<<" in *CURVE_CELL excess"<<endl;
					cout<<"_num_body="<<_num_body<<endl;
					system("Pause");
					exit(0);
				}
				if (id>=_body_list[body_id]->_nume)
				{
					cout<<"Error: Cell "<<id<<" in body "<<body_id<<" excess"<<endl;
					cout<<"Number of cell in body "<<body_id<<" is "<<_body_list[body_id]->_nume<<endl;
					system("Pause");
					exit(0);
				}
				_body_list[body_id]->create_curve(body_id,1,id);
			}
		}
		else if (a=="*MESH")
		{
			int body_id,total_frame;
			double amplify_ratio;
			while (!exam_keyword(input))
			{
				input>>body_id>>amplify_ratio>>total_frame;
				if (body_id>=_num_body)
				{
					cout<<"Error: Body "<<body_id<<" in *Mesh excess"<<endl;
					cout<<"_num_body="<<_num_body<<endl;
					system("Pause");
					exit(0);
				}
				_body_list[body_id]->create_mesh(body_id,amplify_ratio,total_frame);
			}
		}
		else if (a=="*ARTIFICIAL_VISCOSITY")
		{
			int body_id;
			while (!exam_keyword(input))
			{
				Sarti_vis new_av;
				input>>body_id>>new_av.av_type>>new_av.k1>>new_av.k2;
				if (body_id>=_num_body)
				{
					cout<<"Error: Body"<<body_id<<" in *ARTIFICIAL_VISCOSITY excess"<<endl;
					cout<<"_num_body="<<_num_body<<endl;
					system("Pause");
					exit(0);
				}
				_body_list[body_id]->_Arti_Vis=generate_AV(new_av);
			}
			//Exam whether the artificial viscosity type is valide
			for (int i = 0; i < _num_body; i++)
			{
				_body_list[i]->exam_av_type();
			}
		}
		else if (a=="*INTERACTION")
		{
			int num_exist_interaction=0;
			while (!exam_keyword(input))
			{
				unsigned int interaction_id;
				int id1,id2;
				Sinteraction new_interaction;
				input>>interaction_id>>new_interaction.interaction_name>>id1>>id2;
				if (id1>=_num_body || id2>=_num_body)
				{
					cout<<"Error: Body "<<id1<<"("<<id2<<")"<<" in *INTERACTION excess"<<endl;
					system("Pause");
					exit(0);
				}
				//Ensure the number of nodes in body1 is larger than that in body2
				if (_body_list[id1]->_nump>_body_list[id2]->_nump)
				{
					new_interaction.body1=_body_list[id1];
					new_interaction.body2=_body_list[id2];
				}
				else
				{
					new_interaction.body1=_body_list[id2];
					new_interaction.body2=_body_list[id1];
				}
				Tinteraction_base* interaction_ptr=generate_interaction(new_interaction);
				if (interaction_id+1>_interaction_list.size())
				{
					_interaction_list.resize(interaction_id+1);
					_interaction_list[interaction_id]=interaction_ptr;
				}
				else
				{
					_interaction_list[interaction_id]=interaction_ptr;
				}
				num_exist_interaction=num_exist_interaction+1;
			}
			_num_interaction=_interaction_list.size();
			if (_num_body>1 && _num_interaction==0)
			{
				cout<<"Warning:No interaction between bodies"<<endl;
				system("Pause");
			}
			if (_num_interaction>num_exist_interaction)
			{
				cout<<"Error:The interaction id is not continuous"<<endl;
				system("Pause");
				exit(0);
			}
			else if (_num_interaction<num_exist_interaction)
			{
				cout<<"Error:Coincident interaction id exists"<<endl;
				system("Pause");
				exit(0);
			}
		}
		else if (a=="*REMAPPING_START_TIME")
		{
			int body_id;
			double remapping_start_time;
			while (!exam_keyword(input))
			{
				input>>body_id>>remapping_start_time;
				if (body_id>=_num_body)
				{
					cout<<"Error: Body "<<body_id<<" in *REMAPPING_SCHEME excess"<<endl;
					system("Pause");
					exit(0);
				}
				string body_name=_body_list[body_id]->_body_name;
				if (body_name!="body_ALE" && body_name!="body_explosion")
				{
					cout<<"Error:Detonation time must be applied body_explosion or body_ALE"<<endl;
					system("Pause");
					exit(0);
				}
				else
				{
					_body_list[body_id]->S_remapping_start_time(remapping_start_time);
				}
			}
		}
		else if (a=="*DETONATION_TIME")
		{
			int body_id;
			double detonation_time;
			while (!exam_keyword(input))
			{
				input>>body_id>>detonation_time;
				if (body_id>=_num_body)
				{
					cout<<"Error: Body "<<body_id<<" in *REMAPPING_SCHEME excess"<<endl;
					system("Pause");
					exit(0);
				}
				string body_name=_body_list[body_id]->_body_name;
				if (body_name!="body_explosion")
				{
					cout<<"Error:Detonation time must be applied body_explosion"<<endl;
					system("Pause");
					exit(0);
				}
				else
				{
					_body_list[body_id]->S_detonation_time(detonation_time);
				}
			}
		}
		else if (a=="*REMAPPING_SCHEME_DTREDUCTION" || a=="*REMAPPING_SCHEME_TIMES" || a=="*REMAPPING_SCHEME_ADAPTION" || a=="*REMAPPING_SCHEME_NONE")
		{
			while (!exam_keyword(input))
			{
				int body_id;
				Sremapping_scheme scheme;
				if (a=="*REMAPPING_SCHEME_DTREDUCTION")
				{
					input>>body_id>>scheme.name>>scheme.DtReductionRatio;
				}
				else if (a=="*REMAPPING_SCHEME_ADAPTION")
				{
					input>>body_id>>scheme.name>>scheme.tolerance_angle>>scheme.tolerance_length;
				}
				else if (a=="*REMAPPING_SCHEME_TIMES")
				{
					input>>body_id>>scheme.name>>scheme.num_remapping;
				}
				else if (a=="*REMAPPING_SCHEME_NONE")
				{
					input>>body_id>>scheme.name;
				}
				scheme.exam_remapping_scheme();
				if (body_id>=_num_body)
				{
					cout<<"Error: Body "<<body_id<<" in *REMAPPING_SCHEME excess"<<endl;
					system("Pause");
					exit(0);
				}
				string body_name=_body_list[body_id]->_body_name;
				if (body_name!="body_ALE" && body_name!="body_ALEMPM" && body_name!="body_explosion")
				{
					cout<<"Error:Remapping scheme must be applied on body_ALE, body_ALEMPM or body_explosion"<<endl;
					system("Pause");
					exit(0);
				}
				else
				{
					_body_list[body_id]->S_remapping_scheme(scheme);
				}
			}
		}
		else if (a=="*REMESH_SCHEME")
		{
			int body_id;
			string name;
			while (!exam_keyword(input))
			{
				input>>body_id>>name;
				if (body_id>=_num_body)
				{
					cout<<"Error: Body "<<body_id<<" in *REMESH_SCHEME excess"<<endl;
					system("Pause");
					exit(0);
				}
				string body_name=_body_list[body_id]->_body_name;
				if (body_name!="body_ALE" && body_name!="body_ALEMPM" && body_name!="body_explosion")
				{
					cout<<"Error:Remesh scheme must be applied on body_ALE, body_ALEMPM or body_explosion"<<endl;
					system("Pause");
					exit(0);
				}
				_body_list[body_id]->S_remesh_scheme(name);
			}
		}
		else if (a=="*CELL_MONITOR")
		{
			int body_id,cell_id;
			while (!exam_keyword(input))
			{
				input>>body_id>>cell_id;
				if (body_id>=_num_body)
				{
					cout<<"Error: Body"<<body_id<<" in *CELL_MONITOR excess"<<endl;
					cout<<"_num_body="<<_num_body<<endl;
					system("Pause");
					exit(0);
				}
				if (cell_id>=_body_list[body_id]->_nume)
				{
					cout<<"Error: Cell "<<cell_id<<" in body "<<body_id<<" excess"<<endl;
					cout<<"Number of cell in body "<<body_id<<" is "<<_body_list[body_id]->_nume<<endl;
					system("Pause");
					exit(0);
				}
				_body_list[body_id]->create_cell_monitor(body_id,cell_id);
			}
		}
		else if (a=="*GLOBAL_VARIABLES")
		{
			int body_id;
			Tbody_base* body_ptr;
			while (!exam_keyword(input))
			{
				input>>body_id;
				if (body_id>=_num_body)
				{
					cout<<"Error: Body"<<body_id<<" in *GLOBAL_VARIABLES excess"<<endl;
					cout<<"_num_body="<<_num_body<<endl;
					system("Pause");
					exit(0);
				}
				body_ptr=_body_list[body_id];
				body_ptr->_global_output->create_ofstream();
				input>>body_ptr->_global_output->_internal_energy>>body_ptr->_global_output->_max_p;
				input>>body_ptr->_global_output->_epeqv;
			}
		}
		else if (a=="*MATERIAL_OUTPUT_CONTROL")
		{
			int body_id,mid;
			string type;
			while (!exam_keyword(input))
			{
				input>>body_id>>mid>>type;
				if (body_id>=_num_body)
				{
					cout<<"Error: Body"<<body_id<<" in *MATERIAL_OUTPUT_CONTROL excess"<<endl;
					cout<<"_num_body="<<_num_body<<endl;
					system("Pause");
					exit(0);
				}
				_body_list[body_id]->S_material_plot(mid,type);
			}
		}
		else if (a=="*GRAVITY")
		{
			int body_id;
			while (!exam_keyword(input))
			{
				vec3D gravity;
				input>>body_id>>gravity.x>>gravity.y>>gravity.z;
				if (body_id>=_num_body)
				{
					cout<<"Error: Body"<<body_id<<" in *GRAVITY excess"<<endl;
					cout<<"_num_body="<<_num_body<<endl;
					system("Pause");
					exit(0);
				}
				_body_list[body_id]->_gravity=gravity;
			}
		}
		else if (a=="*TEST_EXAMPLE")
		{
			int body_id;
			while (!exam_keyword(input))
			{
				string name;
				input>>body_id>>name;
				if (body_id>=_num_body)
				{
					cout<<"Error: Body"<<body_id<<" in *TEST_EXAMPLE excess"<<endl;
					cout<<"_num_body="<<_num_body<<endl;
					system("Pause");
					exit(0);
				}
				_body_list[body_id]->_test_name=name;
			}
		}
		else if (a=="*DIMENSION")
		{
			int body_id;
			while (!exam_keyword(input))
			{
				int dimension;
				input>>body_id>>dimension;
				if (body_id>=_num_body)
				{
					cout<<"Error: Body"<<body_id<<" in *DIMENSION excess"<<endl;
					cout<<"_num_body="<<_num_body<<endl;
					system("Pause");
					exit(0);
				}
				_body_list[body_id]->_dimension=dimension;
			}
		}
		else if (a=="*END")
		{
			input.close();	
			return;
		}
	}
	return;
}
//----------------------------------------------------------------
//Functions for test
//----------------------------------------------------------------
void Tregion::test()
{
	//_body_list[0]->test_for_MM_remapping();
	_body_list[0]->test_for_MoF();
	//_body_list[0]->test_for_MoF_patch_test_3_material();
	//_body_list[0]->test_for_MoF_patch_test();
	//_body_list[0]->test_for_MoF_derivative();
	//_body_list[0]->test_for_volume_solver();
	//_body_list[0]->set_displacemant_at_conode(0.4);
	//_body_list[0]->test_for_fraction_below();
	return;
}