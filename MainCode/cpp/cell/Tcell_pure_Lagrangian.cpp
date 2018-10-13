#include "Tcell_pure_Lagrangian.h"
#include "Tnode_fluid.h"
#include <iostream>
#include <fstream>
#include "public_function.h"
using namespace std;
//--------------------------------------------------------------------------------
//Initialize the static variables
//-------------------------------------------------------------------------------
int Tcell_pure_Lagrangian::cell_face[6][4]={{0,3,2,1},{0,1,5,4},{1,2,6,5},{2,3,7,6},{0,4,7,3},{4,5,6,7}};
int Tcell_pure_Lagrangian::cell_edge_f[12][2]={{4,2},{1,3},{2,4},{3,1},{4,2},{1,3},{2,4},{3,1},{0,5},{0,5},{0,5},{0,5}};
int Tcell_pure_Lagrangian::cell_edge_f1[12][2]={{1,0},{2,0},{3,0},{4,0},{5,1},{5,2},{5,3},{5,4},{4,1},{1,2},{2,3},{3,4}};
int Tcell_pure_Lagrangian::cell_edge_v[12][2]={{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
int Tcell_pure_Lagrangian::cell_corner_v[8][3]={{1,3,4},{2,0,5},{3,1,6},{0,2,7},{7,5,0},{4,6,1},{5,7,2},{6,4,3}};
int Tcell_pure_Lagrangian::cell_corner_f[8][3]={{4,1,0},{1,2,0},{2,3,0},{3,4,0},{1,4,5},{2,1,5},{3,2,5},{4,3,5}};
//----------------------------------------------------------------------------------
//Member Functions
//----------------------------------------------------------------------------------
Tcell_pure_Lagrangian::Tcell_pure_Lagrangian(int cell_id,int nGauss,Tnode_fluid* (&node_ptr)[8],TMAT_base** mat) : Tcell_base(cell_id,nGauss)
{
	//Allocate memory
	_corner_stress_force.resize(8);
	//Set cell node
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]=node_ptr[i];
	}
	//Set Gauss point material
	for (int i = 0; i < nGauss; i++)
	{
		_gausspoint[i]=mat[i];
	}
	//Set geometry information
	bool negative;
	_cell_volume=calculate_cell_volume("Full",negative);
	if (negative)
	{
		cout<<"Error: Negative volume of the initial input"<<endl;
		cout<<"cell id="<<_id;
		system("Pause");
		exit(0);
	}
	_cell_length=0;
	_cell_mass=_cell_volume*_gausspoint[0]->G_density();
	calculate_corner_mass();
	for (int i = 0; i < 6; i++)
	{
		_adjacent[i]=NULL;
	}
	_cell_energy=0;
	_flag=0;
}
void Tcell_pure_Lagrangian::calculate_corner_force(Tviscosity_base* Arti_Vis,double dt)
{
	vec3D assistant[8][3];      //The 24 assistant direction area pointing from i to cell_corner_v[i][j]
	vec3D surfacial[8][3];      //The 24 surfacial direction area at cell_corner_f[i][j]
	calculate_auxiliary_area(surfacial,assistant);
	calculate_corner_stress_force(surfacial,assistant);
	calculate_corner_hourglass_force(surfacial,assistant);
	Arti_Vis->calculate_viscosity_corner_force(this,assistant,dt);
	return;
}
double Tcell_pure_Lagrangian::calculate_time_step()
{
	double v,v_max;                    //Calculate the velocity of node to get the maximal velocity
	v_max=0;
	for (int i=0;i<8;i++)
	{
		v=_node_ptr[i]->_velocity*_node_ptr[i]->_velocity;
		if (v>v_max)
		{
			v_max=v;
		}
	}
	//Calcualte characteristic length when calculating time step
	_cell_length=calculate_characteristic_length();
	return _cell_length/(_gausspoint[0]->G_soundspeed()+sqrt(v_max));
}
bool Tcell_pure_Lagrangian::calculate_corner_stress_force(vec3D (&surfacial)[8][3],vec3D (&assistant)[8][3])
{
	double temp_volume;                     //The volume of cell at predict step
	double temp_pressure;                   //The pressure at predict step
	bool negative;
	temp_volume=calculate_cell_volume("Temporary",negative);
	double ro,e;
	ro=_gausspoint[0]->G_density();e=_gausspoint[0]->G_internal_energy();
	temp_pressure=_gausspoint[0]->G_EOS()->predict_p(ro,e,temp_volume,_cell_volume);
	for (int i=0;i<8;i++)
	{
		for (int j=0;j<3;j++)
		{
			_corner_stress_force[i]=_corner_stress_force[i]+surfacial[i][j]*temp_pressure;
		}
	}
	return negative;
}
void Tcell_pure_Lagrangian::calculate_corner_hourglass_force(vec3D (&surfacial)[8][3],vec3D (&assistant)[8][3])
{
	double dp[8];                    //The difference of pressure due to hourglass
	for (int i = 0; i < 8; i++)
	{
		dp[i]=0;
	}
	double ro=_gausspoint[0]->G_density();
	double p=_gausspoint[0]->G_pressure();
	double e=_gausspoint[0]->G_internal_energy();
	double c=_gausspoint[0]->G_soundspeed();
	for (int i=0;i<8;i++)
	{
		dp[i]=_gausspoint[0]->G_EOS()->ro_e_p(_corner_mass[i]/_corner_volume[i],e)-p;
	}
	for (int i=0;i<8;i++)
	{
		for (int j=0;j<3;j++)
		{
 			_corner_hourglass_force[i]=_corner_hourglass_force[i]+surfacial[i][j]*dp[i]
 			+assistant[i][j]*(dp[i]-dp[Tcell_pure_Lagrangian::cell_corner_v[i][j]])*0.5;
		}
	}
	return;
}
double Tcell_pure_Lagrangian::calculate_characteristic_length()
{
	double l,l_min=100000;
	for (int i = 0; i < 4; i++)
	{
		l=(_node_ptr[i]->_position-_node_ptr[(i+1)%4]->_position).self_multuply();
		l_min=minval(l,l_min);
		l=(_node_ptr[i+4]->_position-_node_ptr[(i+1)%4+4]->_position).self_multuply();
		l_min=minval(l,l_min);
		l=(_node_ptr[i]->_position-_node_ptr[i+4]->_position).self_multuply();
		l_min=minval(l,l_min);
	}
	return sqrt(l_min);
}
void Tcell_pure_Lagrangian::calculate_segment_center(vec3D &c_cell,vec3D (&c_face)[6],vec3D (&c_edge)[12],string type)
{
	vec3D node_pos[8];
	if (type=="Temporary")
	{
		for (int i = 0; i < 8; i++)
		{
			node_pos[i]=_node_ptr[i]->_position_temp;
		}
	}
	else if (type=="Full")
	{
		for (int i = 0; i < 8; i++)
		{
			node_pos[i]=_node_ptr[i]->_position;
		}
	}
	else
	{
		cout<<"Error: the type in Function \"calculate_segment_center\" must be \"Full\" or \"Temporary\"!"<<endl;
		system("Pause");
		exit(0);
	}
	c_cell.value(0,0,0);
	 for (int i=0;i<8;i++)
	 {
		 c_cell=c_cell+node_pos[i];	
	 }
	 c_cell=c_cell*0.125;
	 for (int i=0;i<6;i++)
	 {
		 c_face[i].value(0,0,0);
		 for (int j=0;j<4;j++)
		 {
			 c_face[i]=c_face[i]+node_pos[cell_face[i][j]];
		 }
		 c_face[i]=c_face[i]*0.25;
	 }
	 for (int i=0;i<12;i++)
	 {
		 c_edge[i]=(node_pos[cell_edge_v[i][0]]+node_pos[cell_edge_v[i][1]])*0.5;
	 }
	 return;
}
double Tcell_pure_Lagrangian::calculate_cell_volume(string type,bool &negative)
{
	vec3D c_cell,c_face[6],c_edge[12];     //The cell center,face center and edge center
	int n1,n2,s1,s2;
	double temp;
	calculate_segment_center(c_cell,c_face,c_edge,type);
	negative=false;
	if (type=="Full")
	{
		for (int i=0;i<8;i++)
		{
			_corner_volume[i]=0;
		}
		for (int i=0;i<12;i++)
		{
			n1=cell_edge_v[i][0];n2=cell_edge_v[i][1];
			s1=cell_edge_f1[i][0];s2=cell_edge_f1[i][1];
			temp=((c_face[s1]-c_cell)%(c_edge[i]-c_cell)+(c_edge[i]-c_cell)%(c_face[s2]-c_cell))*(_node_ptr[n2]->_position-_node_ptr[n1]->_position);
			temp=temp/12.0;
			_corner_volume[n1]=_corner_volume[n1]+temp;
			_corner_volume[n2]=_corner_volume[n2]+temp;
		}
		_cell_volume=0;
		for (int i=0;i<8;i++)
		{
			if (_corner_volume[i]<0)
			{
				negative=true;
			}
			_cell_volume=_cell_volume+_corner_volume[i];
		}
		return _cell_volume;
	}
	else if (type=="Temporary")
	{
		double temp_volume=0;
		for (int i=0;i<12;i++)
		{
			n1=cell_edge_v[i][0];n2=cell_edge_v[i][1];
			s1=cell_edge_f1[i][0];s2=cell_edge_f1[i][1];
			temp=((c_face[s1]-c_cell)%(c_edge[i]-c_cell)+(c_edge[i]-c_cell)%(c_face[s2]-c_cell))*(_node_ptr[n2]->_position_temp-_node_ptr[n1]->_position_temp);
			temp=temp/6.0;
			temp_volume=temp_volume+temp;
		}
	 if (temp_volume<=0)
	 {
		negative=true;
	 }
	 return temp_volume;
	}
	else
	{
		cout<<"Error:Invalid type in calling calculate_cell_volume()"<<endl;
		system("Pause");
		exit(0);
	}
	return 0;
}
void Tcell_pure_Lagrangian::calculate_auxiliary_area(vec3D (&surfacial)[8][3],vec3D (&assistant)[8][3])
{
	vec3D c_cell,c_face[6],c_edge[12];                    //The center of the cell
	int n0,n1,n2,s1,s2;             
	vec3D b0,b1,b2;
	vec3D temp;

	calculate_segment_center(c_cell,c_face,c_edge,"Full");
	//Calculate the surfacial direction
	for (int i=0;i<8;i++)
	{
		for (int j=0;j<3;j++)
		{
			n1=Tcell_pure_Lagrangian::cell_corner_v[i][(j+1)%3];
			n2=Tcell_pure_Lagrangian::cell_corner_v[i][(j+2)%3];
			s1=Tcell_pure_Lagrangian::cell_corner_f[i][j];
			b0=_node_ptr[i]->_position;
			b1=(_node_ptr[n1]->_position+b0)*0.5;
			b2=(_node_ptr[n2]->_position+b0)*0.5;
			surfacial[i][j]=((b2-b0)%(c_face[s1]-b0)+(c_face[s1]-b0)%(b1-b0))*0.5;
		}
	}
	//Calculate the assistant direction
	for (int i=0;i<12;i++)
	{
		n1=Tcell_pure_Lagrangian::cell_edge_v[i][0];n2=Tcell_pure_Lagrangian::cell_edge_v[i][1];
		s1=Tcell_pure_Lagrangian::cell_edge_f1[i][0];s2=Tcell_pure_Lagrangian::cell_edge_f1[i][1];
		temp=((c_face[s1]-c_cell)%(c_edge[i]-c_cell)+(c_edge[i]-c_cell)%(c_face[s2]-c_cell))*0.5;
		n0=cell_corner_v_reverse(n1,n2);
		assistant[n1][n0]=temp;
		n0=cell_corner_v_reverse(n2,n1);
		assistant[n2][n0]=temp*(-1);
	}
	return;
}
int Tcell_pure_Lagrangian::cell_corner_v_reverse(int p0,int p1)
{
	for (int i=0;i<3;i++)
	{
		if (p1==cell_corner_v[p0][i])
		{
			return i;
		}
	}
	cout<<"Error:This node isn't belong to this corner!!!"<<endl;
	system("Pause");
	exit(0);
	return 0;
}
int Tcell_pure_Lagrangian::vertex_reverse(int nid)
{
	for (int i=0;i<8;i++)
	{
		if (_node_ptr[i]->_id==nid)
		{
			return i;
		}
	}
	cout<<"Error: this node isn't belong to this cell!!"<<endl;
	system("Pause");
	exit(0);
	return 0;
}
Tnode_fluid* Tcell_pure_Lagrangian::access_extended_node(int local_face_id,int local_node_id)
{
	int global_id_on_face[4];            //The global number of nodes on this face
	int global_id_target;                //The global id of _node_ptr[local_node_id]
	int global_id_trial;                 //The global id of node on the adjacent cell
	Tcell_pure_Lagrangian* adjacent_cell;    //The adjacent cell to this face 
	int local_id;                        //The local node id of _node_ptr[local_node_id] on the adjacent cell
	int temp;
	for (int i=0;i<4;i++)
	{
		global_id_on_face[i]=_node_ptr[cell_face[local_face_id][i]]->_id;
	}
	global_id_target=_node_ptr[local_node_id]->_id;
	//If the node is located on boundary, return null
	if (_adjacent[local_face_id]==NULL)
	{
		return NULL;
	}
	else
	{
		adjacent_cell=_adjacent[local_face_id];
	}
	local_id=adjacent_cell->vertex_reverse(global_id_target);
	for (int i=0;i<3;i++)
	{
		temp=0;
		global_id_trial=adjacent_cell->_node_ptr[cell_corner_v[local_id][i]]->_id;
		for (int j=0;j<4;j++)
		{
			if (global_id_trial==global_id_on_face[j])
			{
				temp=temp+1;
			}
		}
		if (temp==0)
		{
			return adjacent_cell->_node_ptr[cell_corner_v[local_id][i]];
		}
	}
	cout<<"Error: the adjacent node isn't found"<<endl;
	system("Pause");
	exit(0);
	return 0;
}
void Tcell_pure_Lagrangian::set_cell_connection()
{
	//Calculate the cells which has at least one common node to this cell
	_flag=1;
	int num_surround;
	for (int i=0;i<8;i++)
	{
		num_surround=_node_ptr[i]->_surrounding.size();
		for (int j=0;j<num_surround;j++)
		{
			if (_node_ptr[i]->_surrounding[j]->_flag==0)
			{
				_around.push_back(_node_ptr[i]->_surrounding[j]);
				//Set flag=1 will ensure this cell only being pushed back once
				_node_ptr[i]->_surrounding[j]->_flag=1;
			}
		}
	}
	_flag=0;
	int num_around=_around.size();
	for (int i=0;i<8;i++)
	{
		num_surround=_node_ptr[i]->_surrounding.size();
		for (int j=0;j<num_surround;j++)
		{
			//Recover the flags
			_node_ptr[i]->_surrounding[j]->_flag=0;
		}
	}
	//Calculate the cells which is adjacent to this cell
	for (int i=0;i<6;i++)
	{
		for (int j=0;j<4;j++)
		{
			_node_ptr[cell_face[i][j]]->_flag=1;
		}
		for (int j=0;j<num_around;j++)
		{
			int temp=0;
			for (int k=0;k<8;k++)
			{
				temp=temp+_around[j]->_node_ptr[k]->_flag;
			}
			if (temp==4)
			{
				_adjacent[i]=_around[j];
				break;
			}
		}
		for (int j=0;j<4;j++)
		{
			_node_ptr[cell_face[i][j]]->_flag=0;
		}
	}
	return;
}
bool Tcell_pure_Lagrangian::update_state(double dt)
{
	//Calculate the work by pressure and viscosity force
	double de_p=0,de_v=0,de_g=0;	
	for (int i=0;i<8;i++)
	{
		de_p=de_p+_node_ptr[i]->_velocity_temp*_corner_stress_force[i];
		de_v=de_v+_node_ptr[i]->_velocity_temp*_corner_viscosity_force[i];
		de_g=de_g+_node_ptr[i]->_velocity_temp*_corner_hourglass_force[i];
	}
	if (de_v>1e-5)
	{
		//cout<<"Error in Aritficitial force work!!!"<<endl;
		//cout<<de_v<<endl;
	}
	//Calculate internal energy before cell volume is updated
	double total_internal_energy=_gausspoint[0]->G_density()*_gausspoint[0]->G_internal_energy()*_cell_volume;
	//Update geometry variables
	bool negative;
	calculate_cell_volume("Full",negative);
	//Update density after cell volume is updated
	_gausspoint[0]->S_density(_cell_mass/_cell_volume);
	//Update energy by compatible scheme	
	total_internal_energy=total_internal_energy-(de_p+de_v+de_g)*dt;
	bool failed;
	_gausspoint[0]->update_state(_cell_mass,_cell_volume,total_internal_energy,failed);
	if (failed)
	{
		cout<<"Error:soundspeed in cell "<<_id<<" is negative"<<endl;
		system("Pause");
		exit(0);
	}

	//The corner force must be set to 0 when the energy is calculated
	for (int i=0;i<8;i++)
	{
		_corner_stress_force[i].value(0,0,0);
		_corner_viscosity_force[i].value(0,0,0);
		_corner_hourglass_force[i].value(0,0,0);
	}
	return negative;
}
void Tcell_pure_Lagrangian::assemble_corner_inertance()
{
	double mass=_cell_volume*_gausspoint[0]->G_density();
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]->accumulate_inertance(mass/8.0);
	}
	return;
}
void Tcell_pure_Lagrangian::add_to_node()
{
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]->add_connected_cell(this);
	}
	return;
}
void Tcell_pure_Lagrangian::assemble_corner_force()
{
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]->_force=_node_ptr[i]->_force+_corner_stress_force[i]+_corner_hourglass_force[i]+_corner_viscosity_force[i];
		_node_ptr[i]->_force0=_node_ptr[i]->_force0+_corner_stress_force[i]+_corner_hourglass_force[i]+_corner_viscosity_force[i];
	}
	return;
}
void Tcell_pure_Lagrangian::calculate_average_variable(double  &av_density,double &av_pressure,double &av_soundspeed)
{
	av_density=_gausspoint[0]->G_density();
	av_pressure=_gausspoint[0]->G_pressure();
	av_soundspeed=_gausspoint[0]->G_soundspeed();
	return;
}
void Tcell_pure_Lagrangian::calculate_distorted_ratio(double &ratio_angle,double &ratio_length)
{
	double l_max=0;             //The maximum length of a cell
	double l_min=1000000;       //The minimum length of a cell
	ratio_angle=0;              //The maximum cos of angle
	vec3D p1,p2,p3;
	int i1,i2,i3;
	//Calculate ratio_length
	for (int i=0;i<12;i++)
	{
		i1=cell_edge_v[i][0];i2=cell_edge_v[i][1];
		p1=_node_ptr[i1]->_position;p2=_node_ptr[i2]->_position;
		double l=sqrt((p1-p2)*(p1-p2));
		if (l<l_min) l_min=l;
		if (l>l_max) l_max=l;
	}
	ratio_length=l_max/l_min;
	//Calculate ratio angle
	for (int i=0;i<8;i++)
	{
		i1=i;
		p1=_node_ptr[i1]->_position;
		for (int j=0;j<3;j++)
		{
			i2=cell_corner_v[i][j];i3=cell_corner_v[i][(j+1)%3];
			p2=_node_ptr[i2]->_position;p3=_node_ptr[i3]->_position;
			double l1=sqrt((p2-p1)*(p2-p1));
			double l2=sqrt((p3-p1)*(p3-p1));
			double dot=abs((p3-p1)*(p2-p1));
			double angle=dot/(l1*l2);
			if (angle>ratio_angle) ratio_angle=angle;
		}
	}
	return;
}
void Tcell_pure_Lagrangian::calculate_corner_mass()
{
	for (int i = 0; i < 8; i++)
	{
		_corner_mass[i]=_corner_volume[i]*_gausspoint[0]->G_density();
	}
	return;
}
//---------------------------------------------------------------------------
//test functions
//---------------------------------------------------------------------------
void Tcell_pure_Lagrangian::output_cell_monitor(ofstream &output)
{
	output<<"TITLE = \"Mesh for Brick fluid Element\""<<endl;
	output<<"VARIABLES = \"X\",\"Y\",\"Z\""<<endl;//,\"density\",\"pressure\",\"vx\",\"vy\",\"vz\",\"mass\""<<endl;
	output<<"ZONE F=FEPOINT,N="<<8<<","<<"E="<<1<<","<<"ET=BRICK"<<endl;
	vec3D dis;
	for (int i = 0; i < 8; i++)
	{
		dis=_node_ptr[i]->_position;
		output<<dis.x<<" "<<dis.y<<" "<<dis.z<<endl;
		   //   <<" "<<_nodeptr_list[i]->_density<<" "<<_nodeptr_list[i]->_pressure
			  //<<" "<<_nodeptr_list[i]->_velocity.x<<" "<<_nodeptr_list[i]->_velocity.y<<" "<<_nodeptr_list[i]->_velocity.z
			  //<<" "<<_nodeptr_list[i]->_mass<<endl;
	}
	output<<"1 2 3 4 5 6 7 8"<<endl;
	return;
}