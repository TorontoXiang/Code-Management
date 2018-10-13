#ifndef TREGION
#define TREGION
#include "Tbody_base.h"
#include "Tinteraction_base.h"
#include <vector>
using namespace std;
//Define the region class

class Tregion
{
public:

	Tregion();

	bool input_region();
	//Input the computational information from input files
	void input_extra_info();
	//Input the output information
	int calculate_time_step();
	//Calculate the minimal time step for the all bodies
	//Return the controlled body
	void time_integration();
	//Advance in time for all bodies

//--------------------------------------------------------------
//Functions for test
//--------------------------------------------------------------
	void test();
private:
	int _num_body;                        //The number of body in the region
	int _num_interaction;                 //The number of interaction pairs in the region
	//Computation control
	double _endtime;                      //The terminal time
	double _current_time;                 //The current time
	double _time_step;                    //The time step
	int _current_step;                     //The current step

	vector<Tbody_base*> _body_list;      //The bodies in the region
	vector<Tinteraction_base*> _interaction_list;
};
#endif