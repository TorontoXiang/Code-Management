#ifndef TCELL_SHELL
#define TCELL_SHELL
#include "Tcell_base.h"
#include "Tnode.h"
#include "data_structure.h"
#include <vector>
#include "Tviscosity_base.h"
//Define cell shell
class Tcell_shell : public Tcell_base
{
public:
	Tcell_shell(){};
	Tcell_shell(int cell_id,int nGauss,Tnode_rotation* (&node_ptr)[4],double h,TMAT_base** mat);

	
	void calculate_corner_force(Tviscosity_base* Arti_vis,double dt);
	//Calculate the corner force
	void assemble_corner_force();
	//Assemble the corner force into nodal force and clear the corner force
	double calculate_time_step();
	//Calculate the cell time step
	void assemble_corner_inertance();
	//Assemble the corner inertance into node
	void update_strain_energy(double dt);
protected:
	Tnode_rotation* _node_ptr[4];
	double _h;     //The thickness of the shell
	double _area;  //The area of the shell
	double _clength; //The characteristic length of the shell
	double _soundspeed_max;    //The maximal sound speed at Gausspoints
	vec3D _corner_stress_moment[4];   //The corner moment 
	double _HGstress[5];              //The hourglass viscosity
	double _strain_energy;            //The strain energy in this cell

	struct auxiliary_solver
	{
		vec3D r21,r31,r41,r42,r34;
		vec3D e1,e2,e3;          //The local coordinate 
		double transfer_matrix[3][3];   //The matrix from local system to global_system
		double y24,y31,x42,x13;  //The component for B materix
		vec3D local_position[4]; //The node coordinate in local coordinate system
		vec3D local_velocity[4]; //The node velocity in local coordinate system
		vec3D local_Omega[4];    //The node Omega in local coordinate system
		double Bvelocity[2][3],BOmega[2][2],NOemga[2];
		//The component for velocity gradient
		//Bvelocity[i][j]=B[i][I]*local_velocity[j][I] (i=1,2;j represents the direction for velocity)
		//BOmega[i][j]=B[i][I]*local_Omega[j][I] (z-direction is zero)
		//NOmega[i]=N[I]*local_Omega[i][I] 
		double invJ;             //The inverse of Jacobi materix
		double integral_sigma[6],integral_z_sigma[3];   //The integral of stress in h direction

		auxiliary_solver();
		//Calculate the local coordinate system
		void calculate_local_system(Tcell_shell* this_cell);
		//Calculate the local variables
		void calculate_local_coordinate(Tcell_shell* this_cell);
		//Calculate the characterstic length and area of cell
		void calculate_geometry_property(Tcell_shell* this_cell);
		//Calculate the component for velocity gradient
		void calculate_component_velocity_gratient();
		//Calculate the integral on sigma and z*sigma
		void calculate_integral_at_Gausspoint(double weight, double z, Tcell_shell* this_cell);
		//Calculate the sum of 0.5*B[I][J]*B[I][J]*area*area
		double calculate_sumb(){return y24*y24+y31*y31+x13*x13+x42*x42;};
	};
	void calculate_strain_rate(double z, const auxiliary_solver &tool, double (&de)[6]);
	//Calculate the strain rate at Gauss point
	void calculate_local_corner_force(auxiliary_solver &tool,vec3D (&local_force)[4],vec3D (&local_moment)[4]);
	//Calculate the corner force and moment in local coordinate system
	void calculate_hourglass_viscosity(auxiliary_solver &tool,double dt,vec3D (&local_force)[4],vec3D (&local_moment)[4]);
	//Calculate the hourglass viscosity term
	void calculate_global_corner_force(auxiliary_solver &tool,vec3D (&local_force)[4],vec3D (&local_moment)[4]);
	//Calculate the corner force in global system
	double calculate_area();
	//Calculate the cell area
	double calculate_characteristic_length();
	//Calculate the characteristic length

	friend class Tbody_shell;

};
#endif