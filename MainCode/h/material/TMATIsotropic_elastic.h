#pragma once
#ifndef TMATISOTROPIC_ELASTIC
#define TMATIS0TROPIC_ELASTIC
//Define the isotropic material
#include "TMAT_base.h"

class TMATIsotropic_elastic : public TMAT_base
{
public:
	TMATIsotropic_elastic(){};
	TMATIsotropic_elastic(double density,double Young,double Possion);

	double calculate_Lambda();      //Calculate Lambda
	double calculate_G();           //Calculate G
	double calculate_K();           //Calculate K
	virtual double calculate_equivalent_stress();   //calculate equivalent stress

	//Initialize the material state
	virtual double calculate_sound_speed();
	//Accumulate the integral for stress on each Gausspoint (for shell element) 
	virtual void calculate_stress_integral_for_shell(double (&sigma)[6],double (&z_sigma)[3],double weight,double z, double q);
	//Calculate the strain energy increment 
	virtual double calculate_strain_energy_increment_for_shell(double (&de)[6],double q,double v,double dt);
	//Update material state
	virtual void update_state(double (&de)[6],double (&vort)[3],double dt,double v=0);
	//Roate stress
	void roate_stress(double (&sR)[6],double (&de)[6],double (&vort)[3],double dt);
	//sR:The roated deviate stress;smR:The roated main stress
	void elastic_trial(double (&sR_trial)[6],double (&de)[6],double dt);
	//Elastic tiral for the deviate stress
	void elastice_p(double (&de)[6],double dt);
	//Update the main stress

	//Access functions
	virtual double G_Youngs(){return _Young;};
	virtual double G_G(){return calculate_G();};
	virtual double G_Poisson(){return _Poisson;};
	virtual void G_stress(double (&stress)[6]){stress[0]=_Sxx+_Sm;stress[1]=_Syy+_Sm;stress[2]=_Szz+_Sm;stress[3]=_Sxy;stress[4]=_Sxz;stress[5]=_Syz;}
	virtual double G_seqv(){return calculate_equivalent_stress();};
	//Value functions
	virtual void S_Youngs(double Young){_Young=Young;};
protected:
	//Material properties
	double _Young;         //Young's modulus
	double _Poisson;       //Poisson's ratio
	//Material stress
	double _Sxx,_Syy,_Szz,_Sxy,_Sxz,_Syz;     //The deviatoric stress
	double _Sm;                               //The spheric stress
	double _Seqv;                             //The equivalent stress
};
#endif