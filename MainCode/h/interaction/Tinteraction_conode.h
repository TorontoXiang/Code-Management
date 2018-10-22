#ifndef TINTERACTION_CONODE
#define TINTERACTION_CONODE
//Define the class for conode interaction
#include "Tinteraction_base.h"
#include "Tbody_shell.h"
#include "Tbody_ALE.h"
#include "Tbody_Lagrangian.h"
#include <vector>
#include "Tbucket_searching.h"
#include "Tbody_brick.h"
using namespace std;
template<class T1,class T2>
class Tinteraction_conode : public Tinteraction_base
{
public:
	Tinteraction_conode(string name,T1 body1,T2 body2);//:Tinteraction_base(name),_body1(body1),_body2(body2){};

	virtual void apply_interaction_force();
	//Calculate the interaction force by conode interaction
	virtual void reset_nodal_inertance();
	//Reset the nodal inertance after remapping
	virtual void modify_after_remapping();
	//Modify the variables after fluid body remapping
	virtual void exam_continuity();
	//Exam the continuity of position,velocity and mass
protected:
	T1 _body1;        //Input body1
	T2 _body2;        //Input body2
	int _num_conode;  //number of conode

	vector<int> _conode_in_body1;     //The global node id of conode in body1
	vector<int> _conode_in_body2;     //The global node id of conode in body2
	                                  //_conode_in_body1[i] and _conode_in_body2[i] are the conode nodes

	void create_conode_list();
	//Create _conode_in_body1 and _conode_in_body2
	void modify_nodal_mass();
	//Modify the mass of conodes
	bool exam_boundary_condition();
	//Exam whether the boundary condition on co-nodes are the same (Only exam the translation freedom)
	
};

template<class T1,class T2> 
Tinteraction_conode<T1,T2>::Tinteraction_conode(string name,T1 body1,T2 body2):Tinteraction_base(name),_body1(body1),_body2(body2)
{
	create_conode_list();
	modify_nodal_mass();
}
template<class T1,class T2> 
void Tinteraction_conode<T1 ,T2>::apply_interaction_force()
{
	int id1,id2;
	vec3D force10,moment10,force20,moment20;
	for (int i = 0; i < _num_conode; i++)
	{
		id1=_conode_in_body1[i];id2=_conode_in_body2[i];
		_body1->_nodeptr_list[id1]->access_nodal_initial_force(force10,moment10);
		_body2->_nodeptr_list[id2]->access_nodal_initial_force(force20,moment20);
		_body1->_nodeptr_list[id1]->conode_force_addition(force20,moment20);
		_body2->_nodeptr_list[id2]->conode_force_addition(force10,moment10);
	}
	return;
	id1=_conode_in_body1[0];id2=_conode_in_body2[0];
	return;
}
template<class T1,class T2>
void Tinteraction_conode<T1,T2>::create_conode_list()
{
	vec3D coor_min,coor_max;
	double cell_edge_max,cell_edge_min;
	_body1->calculate_body_size(coor_min,coor_max,cell_edge_max,cell_edge_min);
	Tbucket_searching bucket(coor_min,coor_max,cell_edge_min);
	//Put the node in body1 into the bucket
	int ix,iy,iz;
	bool access;
	for (int i = 0; i < _body1->_nump; i++)
	{
		bucket.search_bucket_ijk(_body1->_nodeptr_list[i]->_position,ix,iy,iz,access);
		bucket.add_element(i,ix,iy,iz);
	}
	for (int i = 0; i < _body2->_nump; i++)
	{
		bucket.search_bucket_ijk(_body2->_nodeptr_list[i]->_position,ix,iy,iz,access);
		if (access)
		{
			goto breakLoop;
		}
		for (int nx = ix-1; nx < ix+1; nx++)
		{
			for (int ny = iy-1; ny < iy+1; ny++)
			{
				for (int nz = iz-1; nz < iz+1; nz++)
				{
					Sbucket* the_bucket=&bucket._bucket_list[nx][ny][nz];	
					for (int j=0;j<the_bucket->_nume;j++)
					{
						int id_in_body1=the_bucket->eleid[j];
						if (_body1->_nodeptr_list[id_in_body1]->_position==_body2->_nodeptr_list[i]->_position)
						{
							_body1->_nodeptr_list[id_in_body1]->_is_conode=1;
							_body2->_nodeptr_list[i]->_is_conode=1;
							_conode_in_body1.push_back(id_in_body1);
							_conode_in_body2.push_back(i);
							goto breakLoop;
						}
					}
				}
			}
		}
		breakLoop:;
	}
	_num_conode=_conode_in_body1.size();
	if (!exam_boundary_condition())
	{
		cout<<"Error: The boundary condition in co-node conflict"<<endl;
		cout<<"body_id1= "<<_body1->_id<<" body_id2= "<<_body2->_id<<endl;
		system("Pause");
		exit(0);
	}
	return;
}
template<class T1,class T2> 
void Tinteraction_conode<T1,T2>::modify_nodal_mass()
{
	int id1,id2;
	double mass10,mass_rotation10,mass20,mass_rotation20;
	//double total_mass,total_mass_rotation;
	for (int i = 0; i < _num_conode; i++)
	{
		id1=_conode_in_body1[i];id2=_conode_in_body2[i];
		_body1->_nodeptr_list[id1]->access_nodal_initial_inertance(mass10,mass_rotation10);
		_body2->_nodeptr_list[id2]->access_nodal_initial_inertance(mass20,mass_rotation20);
		//total_mass=mass1+mass2;total_mass_rotation=mass_rotation1+mass_rotation2;
		_body1->_nodeptr_list[id1]->conode_mass_addition(mass20,mass_rotation20);
		_body2->_nodeptr_list[id2]->conode_mass_addition(mass10,mass_rotation10);
	}
	return;
}
template<class T1,class T2>
bool Tinteraction_conode<T1,T2>::exam_boundary_condition()
{
	int bc_type_pos1[3],bc_type_pos2[3];
	double bc_vel1[3],bc_vel2[3];
	int id1,id2;
	for (int i = 0; i < _num_conode; i++)
	{
		id1=_conode_in_body1[i];id2=_conode_in_body2[i];
		_body1->_nodeptr_list[id1]->access_nodal_bc(bc_type_pos1,bc_vel1);
		_body2->_nodeptr_list[id2]->access_nodal_bc(bc_type_pos2,bc_vel2);
		for (int j = 0; j < 3; j++)
		{
			if ((bc_type_pos1[j]!=bc_type_pos2[j])||(bc_vel1[j]!=bc_vel2[j]))
			{
				cout<<bc_type_pos1[0]<<" "<<bc_type_pos1[1]<<" "<<bc_type_pos1[2]<<endl;
				cout<<bc_type_pos2[0]<<" "<<bc_type_pos2[1]<<" "<<bc_type_pos2[2]<<endl;
				return false;
			}
		}
	}
	return true;
}
template<class T1,class T2>
void Tinteraction_conode<T1,T2>::reset_nodal_inertance()
{
	//Initilize the nodal inertance
	for (int i = 0; i < _body1->_nump; i++)
	{
		_body1->_nodeptr_list[i]->reset_nodal_inertance();
	}
	for (int i = 0; i < _body2->_nump; i++)
	{
		_body2->_nodeptr_list[i]->reset_nodal_inertance();
	}
}
template<class T1,class T2>
void Tinteraction_conode<T1,T2>::modify_after_remapping()
{
	modify_nodal_mass();
	return;
}
template<class T1,class T2>
void Tinteraction_conode<T1,T2>::exam_continuity()
{
	int id1,id2;
	double pos_error,vel_error,mass_error;
	vec3D pos1,pos2,vel1,vel2;
	double mass1,mass2;
	for (int i = 0; i < _num_conode; i++)
	{
		id1=_conode_in_body1[i];id2=_conode_in_body2[i];
		pos1=_body1->_nodeptr_list[id1]->_position;pos2=_body2->_nodeptr_list[id2]->_position;
		vel1=_body1->_nodeptr_list[id1]->_velocity;vel2=_body2->_nodeptr_list[id2]->_velocity;
		mass1=_body1->_nodeptr_list[id1]->_mass;mass2=_body2->_nodeptr_list[id2]->_mass;
		pos_error=(pos1-pos2).get_length();
		vel_error=(vel1-vel2).get_length();
		mass_error=abs(mass1-mass2);
		if (pos_error>1e-10)
		{
			cout<<"pos_error: "<<pos_error<<" in"<<id1<<" and "<<id2<<endl;
			cout<<"pos of id1: "<<pos1.x<<" "<<pos1.y<<" "<<pos1.z<<endl;
			cout<<"pos of id2: "<<pos2.x<<" "<<pos2.y<<" "<<pos2.z<<endl;
		}
		if (vel_error>1e-10)
		{
			cout<<"vel_error: "<<vel_error<<" in"<<id1<<" and "<<id2<<endl;
			cout<<"vel of id1: "<<vel1.x<<" "<<vel1.y<<" "<<vel1.z<<endl;
			cout<<"vel of id2: "<<vel2.x<<" "<<vel2.y<<" "<<vel2.z<<endl;
		}
		if (mass_error>1e-10)
		{
			cout<<"mass_error: "<<mass_error<<" in "<<id1<<" and "<<id2<<endl;
			cout<<"mass of id1: "<<mass1<<endl;
			cout<<"mass of id2: "<<mass2<<endl;
		}
	}
}
#endif