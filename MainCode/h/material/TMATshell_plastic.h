#ifndef TMATSHELL_PLASTIC
#define TMATSHELL_PLASTIC
//Define the plane stress elastic-plastic material for shell
#include "TMATIsotropic_elastic.h"
class TMATshell_plastic : public TMATIsotropic_elastic
{
public:
	TMATshell_plastic(){};
	TMATshell_plastic(double density,double Young,double Possion,double sigma_y,double ET);

	//Update the state at Gauss point including:
	//stress, volume, density and sound speed
	virtual void update_state(double (&de)[6],double (&vort)[3],double dt,double v=0);
	virtual double G_epeff(){return _Epeqv;};
protected:
	double _sigma_y;        //The yield stress
	double _ET;             //The tangent model
	double _Epeqv;          //The equivalent plastic strain

	void elastic_trial(double (&de)[6],double (&s_old)[6],double &sm_old,double (&s_new)[6],double &sm_new,double dt);
	//The elastic trial solution
	double calculate_trial_equivalent_stress(double (&sigma)[6]);
	//Calculate the equivalent stress of the trial solution
};
#endif