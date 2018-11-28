#include "Tsubstrate.h"
#include <iostream>
#include <fstream>
using namespace std;
int main()
{
	srand((int)time(0));
	ofstream output_tec;
	output_tec.open("subtrate.dat");
	Tsubstrate coupled_substrate;
	bool is_success = false;
	//Tsubstrate straight_substrate(300, 300, 15, 0, 300, 40, 0);
	//straight_substrate.output_substrate(output_tec,0);
	//straight_substrate.generate_coupled_straight_substrate(0, coupled_substrate, is_success);
	Tsubstrate wavy_substrate(300, 300, 15, 0, 300, 25, 0.5, 20);
	wavy_substrate.output_substrate(output_tec,0);
	wavy_substrate.generate_coupled_wavy_substrate(0.5, coupled_substrate, is_success);
	if (is_success)
	{
		cout << "Generate the dual substrate!" << endl;
		coupled_substrate.output_substrate(output_tec,1);
	}
	else
	{
		cout << "Fail to generate the dual substrate..." << endl;
	}
	system("Pause");
	return 0;
}