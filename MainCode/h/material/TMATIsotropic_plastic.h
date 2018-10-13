#pragma once
#ifndef TMATISOTROPIC_PLASTIC
#define TMATIS0TROPIC_PLASTIC
//Define the isotropic material
#include "TMAT_base.h"
#include "TMATIsotropic_elastic.h"

class TMATIsotropic_plastic : public TMATIsotropic_elastic
{
public:
	TMATIsotropic_plastic(){};
	TMATIsotropic_plastic(double density,double Young,double Possion,double sigma_y,double ET);

	virtual void update_state(double (&de)[6],double (&vort)[3],double dt,double v=0);
	//Update the material state
protected:
	//Material properties
	double _sigma_y;        //The yield stress
	double _ET;             //The tangent model
	double _Epeqv;          //The equivalent plastic strain
};
#endif