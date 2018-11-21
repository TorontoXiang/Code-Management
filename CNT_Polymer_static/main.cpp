#include "TFE_grid.h"
#include <iostream>
#include "readin.h"
#include <fstream>
#include "Tcomposite.h"
#include <Windows.h>
#include <omp.h>
using namespace std;

int main(int argc, char* argv[])
{
	int num_thread = atoi(argv[1]);
	omp_set_num_threads(num_thread);
	Tcomposite composite;
	composite.Input();
	double t_begin = GetTickCount();
	composite.CG_iteration();
	double t_end = GetTickCount();
	int total_seconds = int((t_end - t_begin) / 1000);
	int hour = total_seconds / 3600;
	int min = (total_seconds - 3600 * hour) / 60;
	int second = total_seconds - 3600 * hour - 60 * min;
	cout << "Time consuming is " << hour << " hours " << min << " mins " << second << " seconds" << endl;
	composite.output_result();
	system("Pause");
	return 0;
}
