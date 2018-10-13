#include "TEOS.h"
#include <iostream>
#include "public_function.h"
using namespace std;
//------------------------------------------------------------------
//EOS of linear polynomial 
//------------------------------------------------------------------
TEOS_linear_polynomial::TEOS_linear_polynomial(double c0,double c1,double c2,double c3,double c4,double c5,double c6,double ro0) : 
	_c0(c0),_c1(c1),_c2(c2),_c3(c3),_c4(c4),_c5(c5),_c6(c6),_ro0(ro0){}

double TEOS_linear_polynomial::ro_e_p(double ro,double e)
{
	double u;
	u=ro/_ro0-1;
	if (u>0)
	{
		return _c0+_c1*u+_c2*u*u+_c3*u*u*u+(_c4+_c5*u+_c6*u*u)*_ro0*e;
		//return maxval(_c0+_c1*u+_c2*u*u+_c3*u*u*u+(_c4+_c5*u+_c6*u*u)*_ro0*e,1e-5);
	}
	else
	{
		return _c0+_c1*u+_c3*u*u*u+(_c4+_c5*u)*_ro0*e;
		//return maxval(_c0+_c1*u+_c3*u*u*u+(_c4+_c5*u)*_ro0*e,1e-5);
	}
}
double TEOS_linear_polynomial::ro_e_c(double ro,double e,bool &failed)
{
	double u=ro/_ro0-1;
	double v=_ro0/ro;
	double p=ro_e_p(ro,e);
	double temp;
	failed=false;
	if (u>0)
	{
		temp=(_c1+2*_c2*u+3*_c3*u*u+(_c5+2*_c6*u)*_ro0*e)/_ro0+p*v*v*(_c4+_c5*u+_c6*u*u)/_ro0;
	}
	else
	{
		temp=(_c1+3*_c3*u*u+_c5*_ro0*e)/_ro0+p*v*v*(_c4+_c5*u)/_ro0;
	}
	if (temp<0)
	{
		//cout<<"Nagitive sound speed in linear polynomial EOS"<<endl;
		//cout<<"density= "<<ro<<" pressure= "<<p<<" internal energy= "<<e<<" c2= "<<temp<<endl;
		//failed=true;
		temp=1e-10;
	}
	return sqrt(temp);
}
double TEOS_linear_polynomial::p_ro_e(double p,double ro)
{
	return p/(_c4*ro);
}
double TEOS_linear_polynomial::predict_p(double ro,double e,double temp_v,double v)
{
	bool failed;
	double p=ro_e_p(ro,e);
	double c=ro_e_c(ro,e,failed);
	return p-ro*c*c*(temp_v-v)/v;
}
//-------------------------------------------------------------------------
//EOS of polytropic
//-------------------------------------------------------------------------
TEOS_polytropic::TEOS_polytropic(double p0,double ro0,double gama,double b):_p0(p0),_ro0(ro0),_gama(gama),_b(b){};

double TEOS_polytropic::ro_e_p(double ro,double e)
{
	return _p0*(pow(ro/_ro0,_gama)-_b);
}
double TEOS_polytropic::ro_e_c(double ro,double e,bool &failed)
{
	failed=false;
	double p=ro_e_p(ro,e);
	double temp=_gama*p/ro;
	if (temp<=0)
	{
		//cout<<"Neagtive sound speed in polytropic EOS"<<endl;
		//cout<<"density= "<<ro<<" pressure= "<<p<<" internal energy= "<<e<<" c2= "<<temp<<endl;
		//failed=true;
		temp=1e-10;
	}
	return sqrt(temp);
}
double TEOS_polytropic::predict_p(double ro,double e,double temp_v,double v)
{
	double temp_ro=ro*v/temp_v;
	return ro_e_p(temp_ro,e);
}
//-----------------------------------------------------------------------------
//EOS of weak impressible water
//-----------------------------------------------------------------------------
TEOS_weak_impressible::TEOS_weak_impressible(double c,double ro0,double p0):_c(c),_ro0(ro0),_p0(p0){};

double TEOS_weak_impressible::ro_e_p(double ro,double e)
{
	return _c*_c*(ro-_ro0)+_p0;
}
double TEOS_weak_impressible::ro_e_c(double ro,double e,bool &failed)
{
	failed=false;
	return _c;
}
double TEOS_weak_impressible::predict_p(double ro,double e,double temp_v,double v)
{
	double temp_ro=ro*v/temp_v;
	return ro_e_p(temp_ro,e);
}
//----------------------------------------------------------------------
//EOS for JWL
//----------------------------------------------------------------------
TEOS_JWL::TEOS_JWL(double A,double B,double R1,double R2,double w,double ro0):_A(A),_B(B),_R1(R1),_R2(R2),_w(w),_ro0(ro0){};

double TEOS_JWL::ro_e_p(double ro,double e)
{
	double V=_ro0/ro;
	return _A*(1-_w/(_R1*V))*exp(-_R1*V)+_B*(1-_w/(_R2*V))*exp(-_R2*V)+_w*ro*e;
}
double TEOS_JWL::ro_e_c(double ro,double e,bool &failed)
{
	failed=false;
	double V=_ro0/ro;
	double p=ro_e_p(ro,e);
	double temp=(_A*_R1*(1-_w/(_R1*V))-_A*_w/(_R1*V*V))*exp(-_R1*V);
	temp=temp+(_B*_R2*(1-_w/(_R2*V))-_B*_w/(_R2*V*V))*exp(-_R2*V);
	temp=temp*V*V/_ro0+_w*e+p*_w/ro;
	if (temp<=0)
	{
		failed=true;
	}
	return sqrt(temp);
}
double TEOS_JWL::predict_p(double ro,double e,double temp_v,double v)
{
	bool failed;
	double p=ro_e_p(ro,e);
	double c=ro_e_c(ro,e,failed);
	return p-ro*c*c*(temp_v-v)/v;
}
//----------------------------------------------------------------------
//EOS for Mie-Gruneisen //only used in mental material to update sm
//----------------------------------------------------------------------
TEOS_MieGruneisen::TEOS_MieGruneisen(double c0,double s,double gama0,double ro0)
{
	_c1=ro0*c0*c0;_c2=_c1*(2*s-1);_c3=_c1*(s-1)*(3*s-1);_c4=gama0*ro0;_c5=gama0;
	_ro0=ro0;
}
double TEOS_MieGruneisen::ro_e_p(double ro,double e,double dv,double v0)
{
	double u=ro/_ro0-1;
	double A,B=_c5;
	if (u>0)
	{
		A=(_c1*u+_c2*u*u+_c3*u*u*u)*(1-0.5*_c4*u/ro);
	}
	else
	{
		A=_c1*u;
	}
	return (A+B*e)/(1+B*0.5*dv/v0);
}
//-----------------------------------------------------------------------
//Functions for pseudo fluid
//-----------------------------------------------------------------------
TEOS_PSEUDO::TEOS_PSEUDO(double c0,double c1,double c2,double c3,double c4,double c5,double c6,double ro0) : 
	_c0(c0),_c1(c1),_c2(c2),_c3(c3),_c4(c4),_c5(c5),_c6(c6),_ro0(ro0){}