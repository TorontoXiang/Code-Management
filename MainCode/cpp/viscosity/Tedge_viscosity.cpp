#include "Tedge_viscosity.h"
#include "Tnode_fluid.h"
#include <iostream>
#include "public_function.h"
using namespace std;
void Tedge_viscosity::calculate_viscosity_corner_force(Tcell_fluid_base* this_cell,vec3D (&assistant)[8][3],double dt)
{
	int n1,n2;                        //The local id of the vertexes at edge
	int s1,s2;                        //The face id for n1 and n2
	Tnode_fluid *p1,*p2;              //The nodal information of edge endpoint
	Tnode_fluid *q1,*q2;              //The nodal information of the adjacent node of p1 and p2
	vec3D dx0,udx0,dv0,udv0;
	vec3D dxl,dvl,dxr,dvr;
	double lx,lv;
	double rl,rr;
	vec3D arti_corner[8];
	for (int i=0;i<12;i++)
	{
		n1=Tcell_fluid_base::cell_edge_v[i][0];n2=Tcell_fluid_base::cell_edge_v[i][1];
		p1=this_cell->_node_ptr[n1];p2=this_cell->_node_ptr[n2];
		s1=Tcell_fluid_base::cell_edge_f[i][0];s2=Tcell_fluid_base::cell_edge_f[i][1];
		q1=this_cell->access_extended_node(s1,n1);q2=this_cell->access_extended_node(s2,n2);
		dx0=p2->_position-p1->_position;dv0=p2->_velocity-p1->_velocity;
		lx=sqrt(dx0*dx0);lv=sqrt(dv0*dv0);
		if (lx==0)
		{
			cout<<"Error:Two node are coincident!!!"<<endl;
			system("Pause");
			exit(0);
		}
		//Calculate rl and rr
		if (lv*dt/lx<1e-14)
		{
			//Set the rl and rr as 1 for constant velocity field
			rl=rr=1;
		}
		else
		{
			udx0=dx0/lx;udv0=dv0/lv;
			//Calculate rl
			if (q1==NULL)
			{
				//Set the rl as 1 when p1 is boundary node
				rl=1;
			}
			else
			{
				dxl=p1->_position-q1->_position;dvl=p1->_velocity-q1->_velocity;
				double temp=dxl*udx0;
				if (temp==0)
				{
					temp=1;
				}
				rl=dvl*udv0*lx/(temp*lv);
			}
			//Calculate rr
			if (q2==NULL)
			{
				//Set rr=1 when p2 is boundary
				rr=1;
			}
			else
			{
				dxr=p2->_position-q2->_position;dvr=p2->_velocity-q2->_velocity;
				double temp=dxr*udx0;
				if (temp==0)
				{
					temp=1;
				}
				rr=dvr*udv0*lx/(temp*lv);
			}
		}
		//Calculate limiter
		double phy=maxval(0.0,minval(minval(0.5*rl+0.5*rr,2*rl),minval(2*rr,1.0)));
		int n0=this_cell->cell_corner_v_reverse(n1,n2);
		double temp=assistant[n1][n0]*dv0;
		vec3D edge_force;
		if (temp<0)
		{
			double aro,ac,artificial;
			aro=2*p1->_density*p2->_density/(p1->_density+p2->_density);
			ac=minval(p1->_soundspeed,p2->_soundspeed);
			artificial=aro*(_k1*lv+_k2*ac)*(1-phy);
			edge_force=udv0*artificial*temp;
			this_cell->_corner_viscosity_force[n1]=this_cell->_corner_viscosity_force[n1]-edge_force;
			this_cell->_corner_viscosity_force[n2]=this_cell->_corner_viscosity_force[n2]+edge_force;
		}
	}
	return;
}