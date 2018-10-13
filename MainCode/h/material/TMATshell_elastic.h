#ifndef TMATSHELL_ELASTIC
#define TMATSHELL_ELASTIC
//Define the plane stress elastic material
#include "TMATIsotropic_elastic.h"

class TMATshell_elastic : public TMATIsotropic_elastic
{
public:
	TMATshell_elastic(){};
	TMATshell_elastic(double density,double Young,double Possion);

	//Update the state at Gauss point including:
	//stress, volume, density and sound speed
	virtual void update_state(double (&de)[6],double (&vort)[3],double dt,double v=0);
protected:
};
#endif