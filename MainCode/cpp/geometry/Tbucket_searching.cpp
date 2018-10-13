#include "Tbucket_searching.h"
#include "public_function.h"
#include <iostream>
using namespace std;
Tbucket_searching::Tbucket_searching(vec3D coor_min,vec3D coor_max,double cell_edge_max)
	              :_x_min(coor_min.x),_y_min(coor_min.y),_z_min(coor_min.z),
				   _x_max(coor_max.x),_y_max(coor_max.y),_z_max(coor_max.z),
				   _l_max(cell_edge_max)
{
	//Calculate the number of bucket in each dimension
	double ratio=0.8;
	_nx=maxval(1,int(ratio*(_x_max-_x_min)/_l_max));
	_ny=maxval(1,int(ratio*(_y_max-_y_min)/_l_max));
	_nz=maxval(1,int(ratio*(_z_max-_z_min)/_l_max));
	cout<<"dx="<<_x_max-_x_min<<" "<<"dy="<<_y_max-_y_min<<" "<<"dz="<<_z_max-_z_min<<" "<<"l_max="<<_l_max<<" "<<endl;
	cout<<ratio*(_x_max-_x_min)/_l_max<<endl;
	cout<<int(ratio*(_x_max-_x_min)/_l_max)<<endl;
	cout<<"The divisions of the bucket is: "<<endl;
	cout<<"nx="<<_nx<<" ny="<<_ny<<" nz="<<_nz<<endl;
	//Calculate the interval in each dimension
	_dx=(_x_max-_x_min)/_nx;_dy=(_y_max-_y_min)/_ny;_dz=(_z_max-_z_min)/_nz;

	//Allocate memory for bucket_list
	//The first and last bucket in each dimension of bucket_list is empty for convenience in cycle
	_bucket_list=new Sbucket**[_nx+2];
	for (int i=0;i<_nx+2;i++)
	{
		_bucket_list[i]=new Sbucket*[_ny+2];
		for (int j=0;j<_ny+2;j++)	
		{
			_bucket_list[i][j]=new Sbucket[_nz+2];
		}
	}
}
Tbucket_searching::~Tbucket_searching()
{
	for (int i=0;i<_nx+2;i++)
	{
		for (int j=0;j<_ny+2;j++)
		{
			delete[] _bucket_list[i][j];
		}
	}
	for (int i=0;i<_nx+2;i++)
	{
		delete[] _bucket_list[i];
	}
	delete[] _bucket_list;
}
void Tbucket_searching::search_bucket_ijk(vec3D the_point,int &i,int &j,int &k,bool &access)
{
	access=false;
	i=int((the_point.x-_x_min)/_dx)+1;
	j=int((the_point.y-_y_min)/_dy)+1;
	k=int((the_point.z-_z_min)/_dz)+1;
	if (i<1 || i>_nx+1 || j<1 || j>_ny+1 ||k<1 || k>_nz+1)
	{
		access=true;
	}
	return;
}	
void Tbucket_searching::add_element(int id,int i,int j,int k)
{
	_bucket_list[i][j][k]._nume=_bucket_list[i][j][k]._nume+1;
	_bucket_list[i][j][k].eleid.push_back(id);
	return;
}
void Tbucket_searching::search_bucket_range(vec3D point_list[],int num_point,int &nx_min,int &ny_min,int &nz_min,int &nx_max,int &ny_max,int &nz_max,bool &access)
{
	access=false;
	nx_min=ny_min=nz_min=1000000;
	nx_max=ny_max=nz_max=-1000000;
	int nx,ny,nz;
	vec3D temp;
	bool access_local;
	for (int i=0;i<num_point;i++)
	{
		temp=point_list[i];
		search_bucket_ijk(temp,nx,ny,nz,access_local);
		if (nx<nx_min){nx_min=nx;}
		if (nx>nx_max){nx_max=nx;}
		if (ny<ny_min){ny_min=ny;}
		if (ny>ny_max){ny_max=ny;}
		if (nz<nz_min){nz_min=nz;}
		if (nz>nz_max){nz_max=nz;}
	}
	if (access_local)
	{
		access=true;
	}
	//If a node is exactly in the boundary of the region,it locates in the last variable of block_list
	//and n_max+1 in cycle will excess, so we enforce it to locate in n_max-1
	if (nx_max==_nx+1){nx_max=_nx;}
	if (ny_max==_ny+1){ny_max=_ny;}
	if (nz_max==_nz+1){nz_max=_nz;}
	return;
}
