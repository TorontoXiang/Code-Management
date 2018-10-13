#include "TFEM_viscosity.h"
#include "data_structure.h"
#include "public_function.h"
#include "Shape_function.h"
#include "Gauss_weight.h"
TFEM_viscosity::TFEM_viscosity(string name,double k1):Tviscosity_base(name),_k1(k1){}
void TFEM_viscosity::initialize_variables(double (&phy)[8][8],vec3D (&dphy)[8][8])
{
	int id;
	double gp[2];
	gp[0]=Gauss_coordinate[1][0],gp[1]=Gauss_coordinate[1][1];
	for (int i=0;i<8;i++)
	{
		for (int j1=0;j1<2;j1++)
		{
			for (int j2=0;j2<2;j2++)
			{
				for (int j3=0;j3<2;j3++)
				{
					id=4*j1+2*j2+j3;
					phy[i][id]=N_brick_8(i,gp[j3],gp[j2],gp[j1]);
					dphy[i][id]=GN_brick_8(i,gp[j3],gp[j2],gp[j1]);
				}
			}
		}
	}
}
void TFEM_viscosity::calculate_viscosity_corner_force(Tcell_fluid_base* this_cell,vec3D (&assistant)[8][3],double dt)
{
	double phy[8][8];                //The value of 8 shape function at Gauss point
	vec3D dphy[8][8];                //The gradient of 8 shape function at Gauss point
	vec3D key[8][8];                 //The key vector at Gauss point
	double det[8];                   //The determinant of Jacobi at Gauss point
	double S[8][8];                  //The "stiffness" matrix
	int id;                          //The id for every Gauss point
	double qs;					     //The magnitude of artificial viscosity
	vec3D temp;                      //The temporary corner force
	initialize_variables(phy,dphy);
	get_key_variable(this_cell,dphy,det,key);
	//Calculate the "stiffness" matrix
	qs=get_qs(this_cell);
  	for (int i=0;i<8;i++)
  	{
  		for (int j=i;j<8;j++)
 		{
  			S[i][j]=0;
  			for (int i1=0;i1<2;i1++)
 			{
 				for (int i2=0;i2<2;i2++)
 				{
					for (int i3=0;i3<2;i3++)
					{
						id=4*i1+2*i2+i3;
 						S[i][j]=S[i][j]+key[i][id]*key[j][id]*abs(det[id]);  //The weight at each Gausspoint is 1,so they are ignored
					}
				}
			}
			S[j][i]=S[i][j];
 		}
 	}
	//Calculate the corner force
 	for (int i=0;i<8;i++)
 	{
 		temp.x=temp.y=temp.z=0;
  		for (int j=0;j<8;j++)
  		{
  			temp=temp-this_cell->_node_ptr[j]->_velocity_temp*S[i][j]*qs;
  		}
		this_cell->_corner_viscosity_force[i]=temp*_k1;
 	}
	//Calculate the corner limiter
	/*for (int i=0;i<8;i++)
	{
		temp.x=temp.y=temp.z=0;
		for (int j=0;j<8;j++)
		{
			temp=temp-target->G_node(j)->G_velocity_temp()*S[i][j];
		}
		corner_limiter[i]=temp;
	}
	assemble();
	assemble_limiter();*/
	return;
}

void TFEM_viscosity::get_key_variable(Tcell_fluid_base* this_cell,vec3D (dphy)[8][8],double (&det)[8],vec3D (&key)[8][8])
{
	double jacobi[3][3];
	bool degenerate;
	for (int i=0;i<8;i++)
	{
		vec3D temp[3];
		for (int j=0;j<8;j++)
		{
			for (int k = 0; k < 3; k++)
			{
				temp[k]=temp[k]+dphy[j][i]*this_cell->_node_ptr[j]->_position.access(k+1);
			}
		}
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				jacobi[k][j]=temp[j].access(k+1);
			}
		}
		inverse(jacobi,det[i],degenerate);
		if (degenerate)
		{
			cout<<"Warning: Jacobi matrix is degenerated!"<<endl;
			for (int i=0;i<3;i++)
			{
				for (int j=0;j<3;j++)
				{
					jacobi[i][j]=0;
				}
			}
		}
		//Calculate the 8 key vector at ith Gauss point
		for (int j=0;j<8;j++)
		{
			key[j][i]=dphy[j][i].multiply_by_matrix(jacobi);
		}
	}
	return;
}

double TFEM_viscosity::get_qs(Tcell_fluid_base* this_cell)
{
	double phy0[8];           //The shape function at (0,0,0)
	vec3D dphy0[8];            //The gradient of shape function at (0,0,0)
	vec3D key0[8];             //The key vector at (0,0,0)
	double qs=0;
	double jacobi[3][3];
	bool degenerate;
	double det0;
	for (int i=0;i<8;i++)
	{
		phy0[i]=N_brick_8(i,0,0,0);
		dphy0[i]=GN_brick_8(i,0,0,0);
	}
	//Calculate the key vector at (0,0,0)
	vec3D temp[3];
	for (int j=0;j<8;j++)
	{
		for (int k = 0; k < 3; k++)
		{
			temp[k]=temp[k]+dphy0[j]*this_cell->_node_ptr[j]->_position.access(k+1);
		}
	}
	for (int j = 0; j < 3; j++)
	{
		for (int k = 0; k < 3; k++)
		{
			jacobi[k][j]=temp[j].access(k+1);
		}
	}
	inverse(jacobi,det0,degenerate);
	if (degenerate)
	{
		cout<<"Warning: Jacobi matrix is degenerated!"<<endl;
		for (int i=0;i<3;i++)
		{
			for (int j=0;j<3;j++)
			{
				jacobi[i][j]=0;
			}
		}
	}
	for (int j=0;j<8;j++)
	{
		key0[j]=dphy0[j].multiply_by_matrix(jacobi);
	}	
	//Calculate the velocity gradient
	vec3D gvx,gvy,gvz;
	for (int i=0;i<8;i++)
	{
		gvx=gvx+key0[i]*this_cell->_node_ptr[i]->_velocity.x;
		gvy=gvy+key0[i]*this_cell->_node_ptr[i]->_velocity.y;
		gvz=gvz+key0[i]*this_cell->_node_ptr[i]->_velocity.z;
	}
	//Calculate the divergence and the magnitude of curl of velocity
	double dive,curl;
	dive=gvx.x+gvy.y+gvz.z;
	curl=sqrt((gvz.y-gvy.z)*(gvz.y-gvy.z)+(gvz.x-gvx.z)*(gvz.x-gvx.z)+(gvy.x-gvx.y)*(gvy.x-gvx.y));
	//Calculate the magnitude of artificial viscosity
	double k1,k2,h;
	if (dive<0)
	{
		k1=1;
	}
	else
	{
		k1=0;
	}
	k2=1;
	if (curl>1e-6)
	{
		double num=-abs(dive)/curl;
		k2=1-exp(num);
	}
	h=pow(this_cell->_cell_volume,1.0/3.0);
	double av_density,av_soundspeed,av_pressure;
	this_cell->calculate_average_variable(av_density,av_pressure,av_soundspeed);
	return k1*av_density*h*(h*abs(dive)+k2*av_soundspeed);
}