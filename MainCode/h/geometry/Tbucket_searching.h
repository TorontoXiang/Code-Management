#ifndef TBUCKET_SEARCHING
#define TBUCKET_SEARCHING
//Define the algorithm for bucket searching
#include <vector>
#include "data_structure.h"
using namespace std;
struct Sbucket 
{
	Sbucket():_nume(0){};
	int _nume;           //Number of element(cell or node) in the bucket
	vector<int> eleid;   //Element(cell or node) id in the bucket
};
class Tbucket_searching
{
public:
	Tbucket_searching(vec3D coor_min,vec3D coor_max,double cell_edge_max);
	//Input date to construct a bucket searching object
	~Tbucket_searching();

	void search_bucket_ijk(vec3D the_point,int &i,int &j,int &k,bool &asccess);
	//Search the bucket _bucket_list[i][j][k] where the point located
	//access:whether the point is out of range
	void add_element(int id,int i,int j,int k);
	//Add id into _bucket_list[i][j][k]
	void search_bucket_range(vec3D point_list[],int num_point,int &nx_min,int &ny_min,int &nz_min,int &nx_max,int &ny_max,int &nz_max,bool &access);
	//Search the bucket id range for serval points
	double _x_min,_y_min,_z_min;
	double _x_max,_y_max,_z_max;     //Range of the grid
	double _l_max;                   //Maximum length of cells in the grid
	int _nx,_ny,_nz;                 //The number of blocks in three directions
	double _dx,_dy,_dz;              //Size of the blocks in three directions
	Sbucket*** _bucket_list;	     //Information of the back-grid	
};
#endif