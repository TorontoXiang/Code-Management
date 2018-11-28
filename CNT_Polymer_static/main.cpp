#include "TFE_grid.h"
#include <iostream>
#include "readin.h"
#include <fstream>
#include "Tcomposite.h"
#include <Windows.h>
#include <omp.h>
using namespace std;

int main()
{
	omp_set_num_threads(16);
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
	//system("Pause");
	return 0;
}
