#include <iostream>
#include <fstream>
#include <string>
#include "TCNT_net.h"
using namespace std;
int main()
{
	double dE = 1 * 1.60219e-19;    //The height of the barrier 1eV
	int M = 460;                    //The total number of conduction channels
	double sigma_CNT = 1e3;         //The sigma of the CNT
	TCNT_net CNT_net(dE,M,sigma_CNT);
	CNT_net.input_CNT_net();
	CNT_net.Calculate_effective_resistance("x_min", "x_max");
	system("Pause");
	return 0;
}