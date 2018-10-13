#include "TMATIsotropic_elastic.h"
#include <iostream>
using namespace std;
#include <cmath>
TMATIsotropic_elastic::TMATIsotropic_elastic(double density,double Young,double Possion) : TMAT_base(density)
{	
	_Young=Young;_Poisson=Possion;
	_Sxx=_Sxy=_Sxz=_Syy=_Syz=_Szz=_Seqv=_Sm=0;
	_soundspeed=calculate_sound_speed();
}
double TMATIsotropic_elastic::calculate_Lambda()
{
	return _Young*_Poisson/((1+_Poisson)*(1-2*_Poisson));
}
double TMATIsotropic_elastic::calculate_G()
{
	return 0.5*_Young/(1+_Poisson);
}
double TMATIsotropic_elastic::calculate_K()
{
	return _Young/(3*(1-2*_Poisson));
}
double TMATIsotropic_elastic::calculate_equivalent_stress()
{
	double J2=0.5*(_Sxx*_Sxx+_Syy*_Syy+_Szz*_Szz)+_Syz*_Syz+_Sxz*_Sxz+_Sxy*_Sxy;
	return sqrt(J2*3);
}
double TMATIsotropic_elastic::calculate_sound_speed()
{
	double temp=_Young*(1-_Poisson)/((1+_Poisson)*(1-2*_Poisson));
	double c=sqrt(temp/_density);
	return c;
}
void TMATIsotropic_elastic::roate_stress(double (&sR)[6],double (&de)[6],double (&vort)[3],double dt)
{
	double q[3],rot[6];
	sR[0]=_Sxx;sR[1]=_Syy;sR[2]=_Szz;
	sR[3]=_Syz;sR[4]=_Sxz;sR[5]=_Sxy;
    q[0]=2*sR[5]*vort[2];q[1]=2*sR[4]*vort[1];q[2]=2*sR[3]*vort[0];
	rot[0]=-q[0]+q[1];rot[1]=q[0]-q[2];rot[2]=-q[1]+q[2];
	rot[3]=vort[0]*(sR[1]-sR[2])+vort[2]*sR[4]-vort[1]*sR[5];
	rot[4]=vort[1]*(sR[2]-sR[0])+vort[0]*sR[5]-vort[2]*sR[3];
	rot[5]=vort[2]*(sR[0]-sR[1])+vort[1]*sR[3]-vort[0]*sR[4];
	for (int i = 0; i < 6; i++)
	{
		sR[i]=sR[i]+rot[i]*dt;
	}
	sR[0]=sR[0];sR[1]=sR[1];sR[2]=sR[2];
	return;
}
void TMATIsotropic_elastic::elastic_trial(double (&sR_trial)[6],double (&de)[6],double dt)
{
	double dem=(de[0]+de[1]+de[2])/3.0;
	double G=calculate_G();
	sR_trial[0]=sR_trial[0]+2*G*(de[0]-dem)*dt;
	sR_trial[1]=sR_trial[1]+2*G*(de[1]-dem)*dt;
	sR_trial[2]=sR_trial[2]+2*G*(de[2]-dem)*dt;
	sR_trial[3]=sR_trial[3]+G*de[3]*dt;
	sR_trial[4]=sR_trial[4]+G*de[4]*dt;
	sR_trial[5]=sR_trial[5]+G*de[5]*dt;
	return;
}
void TMATIsotropic_elastic::elastice_p(double (&de)[6],double dt)
{
	double dem=(de[0]+de[1]+de[2])*dt;
	double K=calculate_K();
	double dsm=K*dem;
	_Sm=_Sm+dsm;
	return;
}
void TMATIsotropic_elastic::update_state(double (&de)[6],double (&vort)[3],double dt,double v)
{
	double sR[6]={0};
	roate_stress(sR,de,vort,dt);
	//Update the deviate stress
	elastic_trial(sR,de,dt);
	_Sxx=sR[0];_Syy=sR[1];_Szz=sR[2];_Syz=sR[3];_Sxz=sR[4];_Sxy=sR[5];
	//Update the main stress
	elastice_p(de,dt);
	//Update the particle density,soundspeed
	_density=_density/(1+(de[0]+de[1]+de[2])*dt);
	//cout<<"Stress: "<<_Sxx<<" "<<_Syy<<" "<<_Szz<<" "<<endl;
	_soundspeed=calculate_sound_speed();
	return;
}
void TMATIsotropic_elastic::calculate_stress_integral_for_shell(double (&sigma)[6], double (&z_sigma)[3],double weight, double z,double q)
{
	sigma[0]=sigma[0]+weight*(_Sxx+_Sm-q);
	sigma[1]=sigma[1]+weight*(_Syy+_Sm-q);
	sigma[3]=sigma[3]+weight*_Sxy;
	sigma[4]=sigma[4]+weight*_Syz;
	sigma[5]=sigma[5]+weight*_Sxz;
	z_sigma[0]=z_sigma[0]-weight*z*(_Sxx+_Sm);
	z_sigma[1]=z_sigma[1]-weight*z*(_Syy+_Sm);
	z_sigma[2]=z_sigma[2]-weight*z*_Sxy;
	return;
}
double TMATIsotropic_elastic::calculate_strain_energy_increment_for_shell(double (&de)[6],double q,double v,double dt)
{
	double energy_increment;
	energy_increment=(_Sxx+_Sm+q)*de[0]+(_Syy+_Sm+q)*de[1]+(_Szz+_Sm+q)*de[2]+_Syz*de[3]+_Sxz*de[4]+_Sxy*de[5];
	energy_increment=energy_increment*v*dt;
	return energy_increment;
}