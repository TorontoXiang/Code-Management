#ifndef TBODY_MPM
#define TBODY_MPM
//Define the body for MPM method
#include <vector>
#include <iostream>
#include "Tbody_base.h"
#include "Tcell_MPM.h"
#include "TMPM_particle.h"
#include "TMPM_background.h"
using namespace std;
class Tbody_MPM : public Tbody_base
{
public:
	Tbody_MPM(int i,int num_part,string body_name):Tbody_base(i,body_name),_num_part(num_part){};

	virtual void input_body(ifstream& input);
	//Input the MPM body
	virtual void start_up_step();
	//Initialize the background and calculate the nodal force in the first step
	virtual double calculate_time_step();
	//Calculate the time step 
	virtual bool predictor_step(){return false;};
	//Calculate the nodal force at n+1/2 step in the body
	//For shell element, the force at n+1/2 equals to force at n
	virtual bool corrector_step();
	//Update the DOF in the corrector step and calculate the in-body force for next time step
	virtual void calculate_force_in_body();
	//Calculate the nodal force in this body
	virtual void calculate_final_acceleration();
	//Calculte the final acceleration
	virtual void output_tecplot(ofstream& output,double ratio);
	//Output the result in tecplot format
	virtual void output_curve();
	//Output curve variables

	void map_particle_to_background();
	//Map the variable on particles to the background grid
	void calculate_corner_force();
	//Calculate the corner of cells
	void apply_external_force();
	//Apply the external force
private:
	int _num_part;                                 //The number of particle parts
	vector<TMPM_particle>* _particle_group_list;   //The particle group list
	TMPM_background _background;                   //The background


};
#endif