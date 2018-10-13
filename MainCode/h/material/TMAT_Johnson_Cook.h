#pragma once
#ifndef TMAT_JOHNSON_COOK
#define TMAT_JOHNSON_COOK
//Define the isotropic material
#include "TMAT_base.h"
#include "TMATIsotropic_elastic.h"

class TMAT_Johnson_Cook : public TMATIsotropic_elastic
{
public:
	TMAT_Johnson_Cook(){};
	TMAT_Johnson_Cook(double density,double Young,double Possion,double A,double B,double C,double n,double epso,TEOS_base* EOS);

	virtual void update_state(double (&de)[6],double (&vort)[3],double dt,double v);
	//Update the material state
protected:
	//Material properties
	double _A,_B,_C,_n,_epso;      //The parameters for the Johnson Cook model
	double _sigma_y;        //The yield stress
	double _Epeqv;          //The equivalent plastic strain
};
#endif