#ifndef EOS_BASE
#define EOS_BASE
//Define the base class of EOS

class TEOS_base
{
public:
	TEOS_base(){};
	virtual double ro_e_p(double ro,double e){return 0;};
	virtual double ro_e_p(double ro,double e,double dv,double v0){return 0;};
	//Calculate pressure by density and internal energy
	virtual double ro_e_c(double ro,double e,bool &failed){return 0;};
	//Calculate sound speed by density and internal energy
	virtual	double p_ro_e(double p,double ro){return 0;};
	//Calculate internal energy by pressure and density
	virtual double predict_p(double ro,double e,double temp_v,double v){return 0;};
	//Predict the pressure by isentropic process
	virtual double G_ro0(){return 0;};
protected:

};
class TEOS_linear_polynomial :public TEOS_base
{
public:
	TEOS_linear_polynomial(){};
	TEOS_linear_polynomial(double c0,double c1,double c2,double c3,double c4,double c5,double c6,double ro0);

	virtual double ro_e_p(double ro,double e);
	//Calculate pressure by density and internal energy
	virtual double ro_e_c(double ro,double e,bool &failed);
	//Calculate sound speed by density and internal energy
	virtual	double p_ro_e(double p,double ro);
	//Calculate internal energy by pressure and density
	virtual double predict_p(double ro,double e,double temp_v,double v);
	//Predict the pressure by isentropic process
protected:
	double _c0,_c1,_c2,_c3,_c4,_c5,_c6;
	double _ro0;
};
class TEOS_polytropic : public TEOS_base
{
public:
	TEOS_polytropic(){};
	TEOS_polytropic(double p0,double ro0,double gama,double b);

	virtual double ro_e_p(double ro,double e);
	//Calculate pressure by density and internal energy
	virtual double ro_e_c(double ro,double e,bool &failed);
	//Calculate sound speed by density and internal energy
	virtual double predict_p(double ro,double e,double temp_v,double v);
	//Calculate the pressure in predict step

protected:
	double _p0;
	double _ro0;
	double _gama;
	double _b;
};
class TEOS_weak_impressible : public TEOS_base
{
public:
	TEOS_weak_impressible(){};
	TEOS_weak_impressible(double c,double ro0,double p0);

	virtual double ro_e_p(double ro,double e);
	//Calculate pressure by density and internal energy
	virtual double ro_e_c(double ro,double e,bool &failed);
	//Calculate sound speed by density and internal energy
	virtual double predict_p(double ro,double e,double temp_v,double v);
	//Predict the pressure by isentropic process
protected:
	double _ro0;                //The initial density
	double _c;                  //The artificial soundspeed
	double _p0;                 //The initial pressure
};
class TEOS_JWL : public TEOS_base
{
public:
	TEOS_JWL(){};
	TEOS_JWL(double A,double B,double R1,double R2,double w,double ro0);

	virtual double ro_e_p(double ro,double e);
	//Calculate pressure by density and internal energy
	virtual double ro_e_c(double ro,double e,bool &failed);
	//Calculate sound speed by density and internal energy
	virtual double predict_p(double ro,double e,double temp_v,double v);
	//Predict the pressure by isentropic process
protected:
	double _ro0;
	double _A,_B,_R1,_R2,_w;
};
class TEOS_MieGruneisen : public TEOS_base
{
public:
	TEOS_MieGruneisen(){};
	TEOS_MieGruneisen(double c0,double s,double gama0,double ro0);

	virtual double ro_e_p(double ro,double e,double dv,double v0);
	virtual double G_ro0(){return _ro0;};
protected:
	double _ro0;
	double _c1,_c2,_c3,_c4,_c5;
};
class TEOS_PSEUDO : public TEOS_base
{
public:
	TEOS_PSEUDO(){};
	TEOS_PSEUDO(double c0,double c1,double c2,double c3,double c4,double c5,double c6,double ro0);
protected:
	double _c0,_c1,_c2,_c3,_c4,_c5,_c6;
	double _ro0;

protected:
	virtual double ro_e_p(double ro,double e){return 0;};
	//Calculate pressure by density and internal energy
	virtual double ro_e_c(double ro,double e,bool &failed){return 1e-5;};
	//Calculate sound speed by density and internal energy
	virtual	double p_ro_e(double p,double ro){return 0;};
	//Calculate internal energy by pressure and density
	virtual double predict_p(double ro,double e,double temp_v,double v){return 0;};
	//Predict the pressure by isentropic process
};
#endif