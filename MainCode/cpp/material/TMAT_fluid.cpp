#include "TMAT_fluid.h"
TMAT_fluid::TMAT_fluid(double density,double internal_energy,TEOS_base* EOS) : TMAT_base(density,internal_energy,EOS)
{
	bool failed;
	_pressure=calculate_pressure();
	_soundspeed=calculate_sound_speed(failed);
}
void TMAT_fluid::update_state(double mass,double volume,double total_internal_energy,bool &failed)
{
	_density=mass/volume;
	_internal_energy=total_internal_energy/mass;
	_pressure=calculate_pressure();
	_soundspeed=calculate_sound_speed(failed);
	return;
}
double TMAT_fluid::calculate_pressure()
{
	return _EOS->ro_e_p(_density,_internal_energy);
}
double TMAT_fluid::calculate_sound_speed(bool &failed)
{
	return _EOS->ro_e_c(_density,_internal_energy,failed);
}