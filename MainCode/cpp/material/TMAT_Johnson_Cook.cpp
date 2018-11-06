#include "TMAT_Johnson_Cook.h"
#include <iostream>
using namespace std;
TMAT_Johnson_Cook::TMAT_Johnson_Cook(double density,double Young,double Possion,double A,double B,double C,double n,double epso,TEOS_base* EOS):TMATIsotropic_elastic(density,Young,Possion)
{
	_A=A;_B=B;_C=C;_n=n;_epso=epso;
	_sigma_y=A;_Epeqv=0;
	_EOS=EOS;
	_internal_energy=0;
}
void TMAT_Johnson_Cook::update_state(double (&de)[6],double (&vort)[3],double dt,double v)
{
	double sR[6]={0};
	double sxx_old,syy_old,szz_old,syz_old,sxz_old,sxy_old;
	sxx_old=_Sxx;syy_old=_Syy;szz_old=_Szz;syz_old=_Syz;sxz_old=_Sxz;sxy_old=_Sxy;
	roate_stress(sR,de,vort,dt);
	//Update the deviate stress
	elastic_trial(sR,de,dt);
	_Sxx=sR[0];_Syy=sR[1];_Szz=sR[2];_Syz=sR[3];_Sxz=sR[4];_Sxy=sR[5];
	double seqv_trial=calculate_equivalent_stress();
	double G=calculate_G();
	_Epeqv=_Epeqv+0.0001;
	double Plamod=_B*_n*pow(_Epeqv,_n-1);
    _Epeqv=_Epeqv-0.0001;
	if (seqv_trial>_sigma_y)
	{
		double dEpsilon=(seqv_trial-_sigma_y)/(3*G+Plamod);
		_Epeqv=_Epeqv+dEpsilon;
		_sigma_y=(_A+_B*pow(_Epeqv,_n))*(1+_C*log(dEpsilon/_epso/dt));
		double ratio=_sigma_y/seqv_trial;
		_Sxx=_Sxx*ratio;_Syy=_Syy*ratio;_Szz=_Szz*ratio;
		_Syz=_Syz*ratio;_Sxz=_Sxz*ratio;_Sxy=_Sxy*ratio;
	}
	if (_EOS==NULL)
	{
		//Update the main stress
		elastice_p(de,dt);
		_density=_density/(1+(de[0]+de[1]+de[2])*dt);
	}
	else
	{
		//Update the volume
		double v_new=v*(1+(de[0]+de[1]+de[2])*dt);
		double dv=v_new-v;
		//Update the internal energy
		double ie=((_Sxx+sxx_old)*de[0]+(_Syy+syy_old)*de[1]+(_Szz+szz_old)*de[2]+(_Syz+syz_old)*de[3]+(_Sxz+sxz_old)*de[4]+(_Sxy+sxy_old)*de[5])*dt;
		ie=ie*(v+v_new)*0.25+0.5*dv*_Sm;
		_internal_energy=_internal_energy+ie;
		//Calculate v0
		double v0=_density*v/_EOS->G_ro0();
		double e=_internal_energy/v0;
		_density=_density/(1+(de[0]+de[1]+de[2])*dt);
		double p=_EOS->ro_e_p(_density,e,dv,v0);
		_internal_energy=_internal_energy-0.5*dv*p;
		_Sm=-p;
	}
	//Update the main stress
	//elastice_p(de,dt);
	//Update the particle density,soundspeed
	//_density=_density/(1+(de[0]+de[1]+de[2])*dt);
	//cout<<"Stress: "<<_Sxx<<" "<<_Syy<<" "<<_Szz<<" "<<endl;
	_soundspeed=calculate_sound_speed();
	return;
}