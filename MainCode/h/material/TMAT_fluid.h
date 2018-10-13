#ifndef TMAT_FLUID
#define TMAT_FLUID
//Define the material for fluid
#include "TMAT_base.h"

class TMAT_fluid : public TMAT_base
{
public:
	TMAT_fluid(){};
	TMAT_fluid(double density,double internal_energy,TEOS_base* EOS);

	virtual void update_state(double mass,double volume,double internal_energy,bool &failed);
	//Update the state on Gausspoint for fluid material

	virtual double calculate_pressure();
	//Calculate the pressure on Gausspoint

	virtual double calculate_sound_speed(bool &failed);
	//Calculate the sound speed on Gausspoint

};
#endif