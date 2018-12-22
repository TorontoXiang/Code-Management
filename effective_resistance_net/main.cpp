#include <iostream>
#include <fstream>
#include <string>
#include "TCNT_net.h"
using namespace std;
int main()
{
	TCNT_net CNT_net(10000000,400,0.35);
	CNT_net.input_CNT_net();
	CNT_net.Create_tunneling_resistance();
	CNT_net.Create_intrinsic_resistance();
	system("Pause");
	return 0;
}