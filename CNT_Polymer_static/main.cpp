#include "Tcomposite.h"
#include <time.h>
#include <fstream>
using namespace std;
int main()
{
	ofstream output;
	output.open("Time.txt");
	Tcomposite composite;
	composite.Input();
	time_t start = time(NULL);
	composite.CG_iteration();
	time_t end = time(NULL);
	output << end - start << endl;
	composite.output_result();
	return 0;
}
