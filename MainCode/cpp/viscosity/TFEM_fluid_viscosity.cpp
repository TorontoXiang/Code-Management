#include "TFEM_fluid_viscosity.h"
#include "public_function.h"
#include "Shape_function.h"
#include "Gauss_weight.h"
void TFEM_fluid_viscosity::initialize_variables(double (&phy)[8][8],vec3D (&dphy)[8][8])
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
void TFEM_fluid_viscosity::calculate_fluid_viscosity_corner_force(Tcell_fluid_base* this_cell)
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
	_mu=1.25e-5;
	//Calculate the "stiffness" matrix
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
  			temp=temp-this_cell->_node_ptr[j]->_velocity_temp*S[i][j]*_mu;
  		}
		this_cell->_corner_stress_force[i]=this_cell->_corner_stress_force[i]+temp;
 	}
	return;
}
void TFEM_fluid_viscosity::get_key_variable(Tcell_fluid_base* this_cell,vec3D (dphy)[8][8],double (&det)[8],vec3D (&key)[8][8])
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