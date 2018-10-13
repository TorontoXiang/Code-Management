#ifndef INTERACTION_BASE
#define INTERACTION_BASE
//Define the base class for body interaction
#include <string>
using namespace std;
class Tinteraction_base
{
public:
	Tinteraction_base(string name):_interaction_name(name){};

	virtual void apply_interaction_force()=0;
	//Apply the interaction force between bodies
	virtual void reset_nodal_inertance(){};
	//Reset the nodal inertance after remapping
	virtual void modify_after_remapping(){};
	//Modify the interaction object after remapping
	virtual void exam_continuity(){};
	//Exam the continuity of position,velocity and mass in co-node interaction 
protected:
	string _interaction_name;
};
#endif