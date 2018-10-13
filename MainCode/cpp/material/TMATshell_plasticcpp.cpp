#include "TMATshell_plastic.h"
#include <iostream>
using namespace std;
TMATshell_plastic::TMATshell_plastic(double density,double Young,double Possion,double sigma_y,double ET):TMATIsotropic_elastic(density,Young,Possion),
	_sigma_y(sigma_y),_ET(ET),_Epeqv(0){}

void TMATshell_plastic::update_state(double (&de)[6],double (&vort)[3],double dt,double v)
{
	double s_old[6],sm_old;               //The stress at previous iteration
	double s_new[6],sm_new;               //The stress at this iteration
	double Lambda=calculate_Lambda();
	double G=calculate_G();
	double K=calculate_K();
	//Initialize the iterative stress
	s_old[0]=_Sxx;s_old[1]=_Syy;s_old[2]=_Szz;s_old[3]=_Syz;s_old[4]=_Sxz;s_old[5]=_Sxy;sm_old=_Sm;
	//The elastic trial solution
	de[2]=-Lambda*(de[0]+de[1])/(Lambda+2*G);
	elastic_trial(de,s_old,sm_old,s_new,sm_new,dt);
	double seqv=calculate_trial_equivalent_stress(s_new);
	double dEpsilon;
	bool yield=false;
	double sigma_y0=_sigma_y;
	if (seqv>sigma_y0)
	{
		dEpsilon=(seqv-sigma_y0)/(3*G+_ET);
		double sigma_y1=sigma_y0+dEpsilon*_ET;
		double ratio=sigma_y1/seqv;
		//Itration to let sigma_zz equals zero
		int max_itearation=20;
		double eps1=1e-6,eps2=1e-12;
		double Szz_old=0;
		double Szz_new=s_new[2]*ratio+sm_new;       //The sigma_zz at the beginning of the iteration
		double dezz_old=0;
		double dezz_new=de[2];                      //The epsilon_zz at the beginning of the iteration
		de[2]=-(de[0]+de[1]);                       //de[2] is the iterative variable to make sigma_zz equals zero
		for (int i = 0; i < max_itearation; i++)
		{
			elastic_trial(de,s_old,sm_old,s_new,sm_new,dt);
			seqv=calculate_trial_equivalent_stress(s_new);
			dEpsilon=(seqv-sigma_y0)/(3*G+_ET);
			sigma_y1=sigma_y0+dEpsilon*_ET;
			ratio=sigma_y1/seqv;
			//It is a secant iterative method.
			//Szz=f(dezz) so (dezz_old,Szz_old) and (dezz_new,Szz_new) are the known iterative point
			//Using linear approximation, we have dezz_new=dezz_old-Szz_old*ddezz/dSzz
			Szz_old=Szz_new;Szz_new=s_new[2]*ratio+sm_new;
			dezz_old=dezz_new;dezz_new=de[2];
			double dSzz=Szz_new-Szz_old+eps2;
			double ddezz=dezz_new-dezz_old;
			de[2]=dezz_old-Szz_old*ddezz/dSzz;
			if (i>=max_itearation-1)
			{
				cout<<"Warning: The iteration in calculating de[2] in TMATshell_plastic is not convergent"<<endl;
			}
			if (abs(de[2]-dezz_new)/(abs(de[2])+eps2)<eps1)
			{
				//If the relative increment of de[2] is smaller than eps1, break the iteration
				break;
			}
		}
		//Finish the iteration and update the stress by return mapping
		_Sxx=s_new[0]*ratio;_Syy=s_new[1]*ratio;_Szz=s_new[2]*ratio;
		_Syz=s_new[3]*ratio;_Sxz=s_new[4]*ratio;_Sxy=s_new[5]*ratio;
		_Epeqv=_Epeqv+dEpsilon;
		_sigma_y=_sigma_y+dEpsilon*_ET;
		_Seqv=_Seqv*ratio;
		_Sm=sm_new;
		yield=true;
	}
	else
	{
		//Update the stress in elastic
		_Sxx=s_new[0];_Syy=s_new[1];_Szz=s_new[2];
		_Syz=s_new[3];_Sxz=s_new[4];_Sxy=s_new[5];
		_Seqv=seqv;_Sm=sm_new;
	}
	double Epsilon_m=(de[0]+de[1]+de[2])/3.0;
	//Update density
	_density=_density*(1-1.5*Epsilon_m*dt)/(1+1.5*Epsilon_m*dt);
	//Update sound speed (do not use the EOS)
	_soundspeed=calculate_sound_speed();
	return;
}
void TMATshell_plastic::elastic_trial(double (&de)[6],double (&s_old)[6],double &sm_old,double (&s_new)[6],double &sm_new,double dt)
{
	double G=calculate_G();
	double K=calculate_K();
	double Epsilon_m=(de[0]+de[1]+de[2])/3.0;
	//Update main stress
	sm_new=sm_old+3*K*Epsilon_m*dt;
	//Update stress
	s_new[0]=s_old[0]+2*G*(de[0]-Epsilon_m)*dt;s_new[1]=s_old[1]+2*G*(de[1]-Epsilon_m)*dt;s_new[2]=s_old[2]+2*G*(de[2]-Epsilon_m)*dt;
	s_new[3]=s_old[3]+G*de[3]*dt;s_new[4]=s_old[4]+G*de[4]*dt;s_new[5]=s_old[5]+G*de[5]*dt;
	return;
}
double TMATshell_plastic::calculate_trial_equivalent_stress(double (&sigma)[6])
{
	double J2=0.5*(sigma[0]*sigma[0]+sigma[1]*sigma[1]+sigma[2]*sigma[2])+sigma[3]*sigma[3]+sigma[4]*sigma[4]+sigma[5]*sigma[5];
	return sqrt(J2*3);
}