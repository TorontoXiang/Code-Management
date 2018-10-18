#pragma once
#ifndef OUTPUT_CONTROL
#define OUTPUT_CONTROL
//Define the struct for output control
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
struct Smesh_output
{
	double _amplify_ratio;   //The amplify ratio in mesh plot
	int _total_frame;        //The step interval for mesh plot
	ofstream output_initial;
	ofstream output_process;
	ofstream output_final;
	vector<double> _output_point;

	Smesh_output(double amplify_ratio, int total_frame, double endtime);
	Smesh_output() {};
	bool is_output(double current_time, double dt);
};
struct Scurve_output
{
	int _type;               //0-node curve;1-cell curve; 
	int _id;                 //Node or cell id
	ofstream output;

	Scurve_output(int type, int id);
};
#endif