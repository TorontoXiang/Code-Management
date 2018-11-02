#ifndef TMAT_BASE
#define TMAT_BASE
//Define the base class of material
#include <iostream>
class TMAT_base
{
public:
	TMAT_base(){};
	
	TMAT_base(double density)    //Constructor for solid material
	{
		_density=density;_internal_energy=_soundspeed=_pressure=0;
	}	

	virtual void update_state(double (&de)[6],double (&vort)[3],double dt,double v=0){};
	//Update the state on Gausspoint for solid material

	virtual void update_state(double mass,double volume,double total_internal_energy,bool &failed){};
	//Update the state on Gausspoint for fluid material

	virtual double calculate_strain_energy_increment_for_shell(double (&de)[6],double q,double v,double dt){return 0;};
	//Calculate the strain energy increment in solid material

	virtual double calculate_sound_speed(){return 0;};
	virtual double calculate_sound_speed(bool &failed){return 0;}
	//Calculate the sound speed on Gausspoint

	virtual void calculate_stress_integral_for_shell(double (&sigma)[6], double (&z_sigma)[3], double weight, double z,double q){};
	//Accumulate the integral for stress on each Gausspoint (for shell element) 

	virtual double calculate_pressure(){return 0;};
	//Calculate the pressure on Gausspoint

	//Access functions
	double G_soundspeed(){return _soundspeed;};
	double G_density(){return _density;};
	double G_internal_energy(){return _internal_energy;};
	double G_roe(){return _density*_internal_energy;};
	double G_pressure(){return _pressure;};
	virtual double G_Youngs(){return 0;};
	virtual double G_Poisson(){return 0;};
	virtual double G_G(){return 0;};
	virtual double G_epeff(){return 0;};
	virtual void G_stress(double (&stress)[6]){};
	virtual double G_seqv(){return 0;};
	virtual double calculate_equivalent_stress() { return 0; };

	//Value functions
	void S_soundspeed(double soundspeed){_soundspeed=soundspeed;};
	void S_density(double density){_density=density;};
	void S_internal_energy(double internal_energy){_internal_energy=internal_energy;};
	void S_pressure(double pressure){_pressure=pressure;};
	virtual void S_Youngs(double Young){};
protected:
	//The material state
	double _pressure;
	double _density;
	double _soundspeed;
	double _internal_energy;

};
#endif