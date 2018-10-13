#include "TMATIsotropic_plastic.h"
TMATIsotropic_plastic::TMATIsotropic_plastic(double density,double Young,double Possion,double sigma_y,double ET):TMATIsotropic_elastic(density,Young,Possion)
{
	_ET=ET;_sigma_y=sigma_y;
}
void TMATIsotropic_plastic::update_state(double (&de)[6],double (&vort)[3],double dt,double v)
{
	double sR[6]={0};
	roate_stress(sR,de,vort,dt);
	//Update the deviate stress
	elastic_trial(sR,de,dt);
	_Sxx=sR[0];_Syy=sR[1];_Szz=sR[2];_Syz=sR[3];_Sxz=sR[4];_Sxy=sR[5];
	double seqv_trial=calculate_equivalent_stress();
	double G=calculate_G();
	double Plamod=_Young*_ET/(_Young-_ET);
	if (seqv_trial>_sigma_y)
	{
		double dEpsilon=(seqv_trial-_sigma_y)/(3*G+Plamod);
		_sigma_y=_sigma_y+Plamod*dEpsilon;
		double ratio=_sigma_y/seqv_trial;
		_Sxx=_Sxx*ratio;_Syy=_Syy*ratio;_Szz=_Szz*ratio;
		_Syz=_Syz*ratio;_Sxz=_Sxz*ratio;_Sxy=_Sxy*ratio;
	}
	//Update the main stress
	elastice_p(de,dt);
	//Update the particle density,soundspeed
	_density=_density/(1+(de[0]+de[1]+de[2])*dt);
	//cout<<"Stress: "<<_Sxx<<" "<<_Syy<<" "<<_Szz<<" "<<endl;
	_soundspeed=calculate_sound_speed();
	return;
}