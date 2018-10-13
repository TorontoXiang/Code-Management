#ifndef TBODY_PREPROCESSOR
#define TBODY_PREPROCESSOR
//Define the Lagrangian body
#include "Tbody_base.h"
#include "Tcell_fluid_base.h"
#include "Tcell_pure_fluid.h"
#include "Tnode_fluid.h"
#include <vector>
#include <iostream>
#include "readin.h"
class Tbody_preprocessor : public Tbody_base
{
public:
	Tbody_preprocessor(int i,string body_name);
	//Create the body
	virtual void modify_body(ifstream &input);
	//Modify the input body
protected:
	void input_body(ifstream &input);
	//Input the body to be modified
	void output_modified_part(ofstream &output);
	//Out put the modified part
	void modify_k_file();
	//Modify the k file

	Skeyword _keyword;
	ofstream _output;
};
#endif