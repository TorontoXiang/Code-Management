#include "TMATshell_elastic.h"
#include <cmath>
TMATshell_elastic::TMATshell_elastic(double density,double Young,double Possion) : TMATIsotropic_elastic(density,Young,Possion){}
void TMATshell_elastic::update_state(double (&de)[6],double (&vort)[3],double dt,double v)
{
	double Lambda=calculate_Lambda();
	double G=calculate_G();
	double K=calculate_K();
	//Calculate epsilon[z] explicitly to ensure sigma[z] equals 0
	de[2]=-Lambda*(de[0]+de[1])/(Lambda+2*G);
	double Epsilon_m=(de[0]+de[1]+de[2])/3.0;
	//Update main stress
	_Sm=_Sm+3*K*Epsilon_m*dt;
	//Update stress
	_Sxx=_Sxx+2*G*(de[0]-Epsilon_m)*dt;_Syy=_Syy+2*G*(de[1]-Epsilon_m)*dt;_Szz=_Szz+2*G*(de[2]-Epsilon_m)*dt;
	_Syz=_Syz+G*de[3]*dt;_Sxz=_Sxz+G*de[4]*dt;_Sxy=_Sxy+G*de[5]*dt;
	//Update equivalent stress
	_Seqv=calculate_equivalent_stress();
	//Update density
	_density=_density*(1-1.5*Epsilon_m*dt)/(1+1.5*Epsilon_m*dt);
	//Update sound speed (do not use the EOS)
	_soundspeed=calculate_sound_speed();
	return;
}
