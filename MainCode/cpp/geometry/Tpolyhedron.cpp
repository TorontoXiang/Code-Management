#include "Tpolyhedron.h"
#include <iostream>
using namespace std;
void Scutting_ring::add_point(vec3D &new_point)
{
	if (_length==_point_list.size())
	{
		_point_list.push_back(new_point);
		_length=_length+1;
	}
	else
	{
		_point_list[_length]=new_point;
		_length=_length+1;
	}
	return;
}
void Scutting_ring::final_one_ring()
{
	_num_ring=_num_ring+1;
	_strat_point_id[_num_ring]=_length;
	return;
}
void Scutting_ring::calcualte_volume_integral(vec3D &n,double &cutting_area,double &int_1,vec3D &int_x,string type)
{
	double direction;
	vec3D center;
	int ring_strat_id,ring_end_id,ring_length;
	double int_1_in_ring=0;
	vec3D int_x_in_ring;
	cutting_area=0;int_1=0;int_x.value(0,0,0);
	for (int i=0;i<_num_ring;i++)
	{
		//Calculate the center of each ring for triangle decomposition
		//The cutting must be decomposed into triangle pieces in case it is non-convex
		ring_strat_id=_strat_point_id[i];ring_end_id=_strat_point_id[i+1];ring_length=ring_end_id-ring_strat_id;
		center.value(0,0,0);
		for (int j=ring_strat_id;j<ring_end_id;j++)
		{
			center=center+_point_list[j];
		}
		center=center/ring_length;
		double s,nz;
		direction=0;
		int_1_in_ring=0;
		int_x_in_ring.value(0,0,0);
		for (int j=ring_strat_id;j<ring_end_id;j++)		
		{
			//Either the cutting ring is convex or not, this process could give the right answer
			area(center,_point_list[j],_point_list[ring_strat_id+(j+1-ring_strat_id)%ring_length],s,nz);
			int_1_in_ring=int_1_in_ring+gauss(center,_point_list[j],_point_list[ring_strat_id+(j+1-ring_strat_id)%ring_length],"z")*s*nz;
			cutting_area=cutting_area+s;
			if (type=="centroid")
			{
				int_x_in_ring=int_x_in_ring+gauss_for_centroid(center,_point_list[j],_point_list[ring_strat_id+(j+1-ring_strat_id)%ring_length])*s*nz;
			}
			direction=direction+s*nz;
		}
		//Direction is a representative of the order of the cutting ring
		if (direction*n.z<0)
		{
			int_1_in_ring=-int_1_in_ring;
			int_x_in_ring=-int_x_in_ring;
		}
		int_1=int_1+int_1_in_ring;
		int_x=int_x+int_x_in_ring;
	}
	return;
}
void Scutting_ring::determine_ring_order(vec3D &n,int &order,vec3D_int &ring_centroid)
{
	ring_centroid._coor.value(0,0,0);
	double direction=0,s_total=0;	
	double s;
	int max_direction;           //The direction of maximal component of n
	double max_component;        //The value of maximal component of n
	vec3D n_triangle;            //The outward normal of triangle piece in such order
	absmax_component(n,max_direction,max_component);   //Compare with the maximal component to avoid numerical error
	for (int i=1;i<_length-1;i++)
	{
		area(_point_list[0],_point_list[i],_point_list[i+1],s,n_triangle);
		ring_centroid._coor.x=ring_centroid._coor.x+gauss(_point_list[0],_point_list[i],_point_list[i+1],"x")*s;
		ring_centroid._coor.y=ring_centroid._coor.y+gauss(_point_list[0],_point_list[i],_point_list[i+1],"y")*s;
		ring_centroid._coor.z=ring_centroid._coor.z+gauss(_point_list[0],_point_list[i],_point_list[i+1],"z")*s;
		direction=direction+n_triangle.access(max_direction)*s;
		s_total=s_total+s;
	}
	//Direction is a representative of the order of the cutting ring
	if (direction*max_component<0)
	{
		order=1;                //The order of cutting_id and normal is opposite
	}
	else
	{
		order=0;
	}
	if (s_total<1e-21)
	{
		cout<<"s_total<1e-21"<<endl;
		cout<<s_total<<endl;
	}
	ring_centroid._coor=ring_centroid._coor/s_total;
	ring_centroid._type=-1;    //The centroid of this ring must be a new point
	return;
}
vec3D Scutting_ring::calculate_surface_integral(string type)
{
	vec3D result;
	double direction;
	vec3D center;
	int ring_strat_id,ring_end_id,ring_length;
	double cutting_area=0;
	for (int i=0;i<_num_ring;i++)
	{
		//Calculate the center of each ring for triangle decomposition
		//The cutting must be decomposed into triangle pieces in case it is non-convex
		ring_strat_id=_strat_point_id[i];ring_end_id=_strat_point_id[i+1];ring_length=ring_end_id-ring_strat_id;
		center.value(0,0,0);
		for (int j=ring_strat_id;j<ring_end_id;j++)
		{
			center=center+_point_list[j];
		}
		center=center/ring_length;
		double s,nz;
		for (int j=ring_strat_id;j<ring_end_id;j++)		
		{
			//Either the cutting ring is convex or not, this process could give the right answer
			area(center,_point_list[j],_point_list[ring_strat_id+(j+1-ring_strat_id)%ring_length],s,nz);
			if (type=="centroid")
			{
				cutting_area=cutting_area+s;
				result.x=result.x+gauss(center,_point_list[j],_point_list[ring_strat_id+(j+1-ring_strat_id)%ring_length],"x")*s;
				result.y=result.y+gauss(center,_point_list[j],_point_list[ring_strat_id+(j+1-ring_strat_id)%ring_length],"y")*s;
				result.z=result.z+gauss(center,_point_list[j],_point_list[ring_strat_id+(j+1-ring_strat_id)%ring_length],"z")*s;
			}
			else if (type=="derivative")
			{
				//result=result+gauss_for_MoF_derivative(center,_point_list[j],_point_list[ring_strat_id+(j+1-ring_strat_id)%ring_length])*s;
				result.x=result.x+gauss(center,_point_list[j],_point_list[ring_strat_id+(j+1-ring_strat_id)%ring_length],"xy")*s;
				result.y=result.y+gauss(center,_point_list[j],_point_list[ring_strat_id+(j+1-ring_strat_id)%ring_length],"yy")*s;
				result.z=result.z+gauss(center,_point_list[j],_point_list[ring_strat_id+(j+1-ring_strat_id)%ring_length],"xx")*s;
			}
		}
	}
	if (type=="centroid")
	{
		if (cutting_area<1e-21)
		{
			cout<<"cutting_area<1e-21"<<endl;
		}
		result=result/cutting_area;
	}
	return result;
}
void Scutting_ring::calculate_local_coordinate_system(double theta,double phy,vec3D &l1,vec3D &l2,vec3D &l3)
{
	double stheta,ctheta,sphy,cphy;
	ctheta=cos(theta);stheta=sin(theta);
	sphy=sin(phy);cphy=cos(phy);
	l1.x=sphy;l1.y=-cphy;l1.z=0;
	l2.x=ctheta*cphy;l2.y=ctheta*sphy;l2.z=-stheta;
	l3.x=stheta*cphy;l3.y=stheta*sphy;l3.z=ctheta;
	return;
}
void Scutting_ring::clear_ring()
{
	_length=_num_ring=0;
	for (int i = 0; i < 10; i++)
	{
		_strat_point_id[i]=0;
	}
	_point_list.resize(0);
	return;
}
void Scutting_ring::print_points()
{
	for (int i = 0; i < _length; i++)
	{
		cout<<_point_list[i].x<<" "<<_point_list[i].y<<" "<<_point_list[i].z<<endl;
	}
	return;
}
Tpolyhedron::Spiece::Spiece(int n1,int n2,int n3)
{
	_IEN[0]=n1;_IEN[1]=n2;_IEN[2]=n3;
	_adjacent[0]=_adjacent[1]=_adjacent[2]=-1; 
	_type[0]=_type[1]=-1; //Initialize as an invalid variable
	_is_visit=0;
}
void Tpolyhedron::add_vertex(vec3D_double &new_vertex)
{
	_vertex_list.push_back(new_vertex);
	_numv=_numv+1;
	return;
}
void Tpolyhedron::add_vertex(vec3D &new_coordinate)
{
	vec3D_double new_vertex(new_coordinate);
	_vertex_list.push_back(new_vertex);
	_numv=_numv+1;
	return;
}
void Tpolyhedron::add_piece(int n1,int n2,int n3)
{
	Spiece new_piece(n1,n2,n3);
	_piece_list.push_back(new_piece);
	_nump=_nump+1;
	return;
}
double Tpolyhedron::calculate_polyhedron_volume()
{
	double volume=0;
	double s,nz;
	int i1,i2,i3;
	for (int i=0;i<_nump;i++)
	{
		i1=_piece_list[i]._IEN[0];i2=_piece_list[i]._IEN[1];i3=_piece_list[i]._IEN[2];
		area(_vertex_list[i1]._coor,_vertex_list[i2]._coor,_vertex_list[i3]._coor,s,nz);
		volume=volume+gauss(_vertex_list[i1]._coor,_vertex_list[i2]._coor,_vertex_list[i3]._coor,"z")*s*nz;
	}
	_volume=volume;
	return volume;
}
vec3D Tpolyhedron::calculate_polyhedron_centroid()
{
	vec3D centroid;
	double volume=0;
	double s,nz;
	int i1,i2,i3;
	for (int i=0;i<_nump;i++)
	{
		i1=_piece_list[i]._IEN[0];i2=_piece_list[i]._IEN[1];i3=_piece_list[i]._IEN[2];
		area(_vertex_list[i1]._coor,_vertex_list[i2]._coor,_vertex_list[i3]._coor,s,nz);
		volume=volume+gauss(_vertex_list[i1]._coor,_vertex_list[i2]._coor,_vertex_list[i3]._coor,"z")*s*nz;
		centroid.x=centroid.x+gauss(_vertex_list[i1]._coor,_vertex_list[i2]._coor,_vertex_list[i3]._coor,"zx")*s*nz;
		centroid.y=centroid.y+gauss(_vertex_list[i1]._coor,_vertex_list[i2]._coor,_vertex_list[i3]._coor,"zy")*s*nz;
		centroid.z=centroid.z+gauss(_vertex_list[i1]._coor,_vertex_list[i2]._coor,_vertex_list[i3]._coor,"zz")*s*nz*0.5;
	}
	_volume=volume;
	if (volume<1e-21)
	{
		cout<<"volume<1e-21"<<endl;
	}
	centroid=centroid/volume;
	_centroid=centroid;
	return centroid;
}
vec3D_double Tpolyhedron::calculate_polyhedron_center()
{
	vec3D_double p_center;
	//Calcualte the center of the polyhedron
	for (int i=0;i<_numv;i++)
	{
		p_center._coor=p_center._coor+_vertex_list[i]._coor;
	}
	p_center._coor=p_center._coor/_numv;
	return p_center;
}
void Tpolyhedron::clear_the_polyhedron()
{
	_vertex_list.resize(0);
	_piece_list.resize(0);
	_numv=_nump=0;
	_volume=0;
	_centroid.value(0,0,0);
	_cutting_ring.clear_ring();
	_surface_piece_list.resize(0);
}
void Tpolyhedron::calculate_vertex_range(vec3D &coor_min,vec3D &coor_max)
{
	coor_max.value(-1e10,-1e10,-1e10);
	coor_min.value(1e10,1e10,1e10);
	for (int i=0;i<_numv;i++)
	{
		for (int j = 1; j < 4; j++)
		{
			if (_vertex_list[i]._coor.access(j)>coor_max.access(j)) coor_max.value(j,_vertex_list[i]._coor.access(j));
			if (_vertex_list[i]._coor.access(j)<coor_min.access(j)) coor_min.value(j,_vertex_list[i]._coor.access(j));
		}
	}
	return;
}
void Tpolyhedron::plot_this_polyehdron(ofstream &output,double value)
{
	output<<"TITLE=\"material polyhedron\""<<endl;
	output<<"Variables=\"X\",\"Y\",\"Z\",\"V\""<<endl;
	output<<"ZONE F=FEPOINT N= "<<_numv<<" E= "<<_nump<<" ET=TRIANGLE"<<endl;
	for (int i=0;i<_numv;i++)
	{
		output<<_vertex_list[i]._coor.x<<" "<<_vertex_list[i]._coor.y<<" "<<_vertex_list[i]._coor.z<<" "<<value<<endl;
	}
	for (int i=0;i<_nump;i++)
	{
		output<<_piece_list[i]._IEN[0]+1<<" "<<_piece_list[i]._IEN[1]+1<<" "<<_piece_list[i]._IEN[2]+1<<endl;
	}
	return;
}
void Tpolyhedron::move_polyhedron(vec3D &n,double d)
{
	for (int i = 0; i < _numv; i++)
	{
		_vertex_list[i]._coor=_vertex_list[i]._coor+n*d;
	}
	return;
}
vec3D Tpolyhedron::calculate_fraction_below(vec3D &n,double d)
{
	vec3D result=calculate_integral_below(n,d,"volume");
	if (_volume<1e-21)
	{
		cout<<" Tpolyhedron::calculate_fraction_below(vec3D &n,double d)"<<endl;
	}
	return result/_volume;
}
vec3D Tpolyhedron::calcualte_centroid_below(vec3D &n,double d)
{
	return calculate_integral_below(n,d,"centroid");
}
vec3D Tpolyhedron::calculate_integral_below(vec3D &n,double d,string type)
{
	vec3D result;
	classify_pieces(n,d);
	_cutting_ring._length=_cutting_ring._num_ring=0;
	_cutting_ring._strat_point_id[0]=0;
	double volume_below=0;      //The volume below the cutting plane
	double int_1=0;             //The result of the integral on 1 for volume
	vec3D int_x;                //The result of the integral on x for centroid
	double sub_int_1;           //The integral on sub-piece on 1 for volume
	vec3D sub_int_x;            //The integral on sub-piece on x for centroid
	double cutting_plane_int_1; //The integral on cutting-plane on 1 for volume
	vec3D cutting_plane_int_x;  //The integral on cutting-plane on x for centroid
	double cutting_area;        //The area of the cutting plane
	for (int i=0;i<_nump;i++)
	{
		double s,nz;
		if (_piece_list[i]._is_visit==0)
		{
			if (_piece_list[i]._type[0]==2)
			{
				//This piece is below the cutting plane so the integral is calculated in the whole piece
				area(_vertex_list[_piece_list[i]._IEN[0]]._coor,_vertex_list[_piece_list[i]._IEN[1]]._coor,_vertex_list[_piece_list[i]._IEN[2]]._coor,s,nz);
				int_1=int_1+gauss(_vertex_list[_piece_list[i]._IEN[0]]._coor,_vertex_list[_piece_list[i]._IEN[1]]._coor,_vertex_list[_piece_list[i]._IEN[2]]._coor,"z")*s*nz;
				if (type=="centroid")
				{
					//int_x=int_x+gauss_for_centroid(_vertex_list[_piece_list[i]._IEN[0]]._coor,_vertex_list[_piece_list[i]._IEN[1]]._coor,_vertex_list[_piece_list[i]._IEN[2]]._coor)*s*nz;
					int_x.x=int_x.x+gauss(_vertex_list[_piece_list[i]._IEN[0]]._coor,_vertex_list[_piece_list[i]._IEN[1]]._coor,_vertex_list[_piece_list[i]._IEN[2]]._coor,"zx")*s*nz;
					int_x.y=int_x.y+gauss(_vertex_list[_piece_list[i]._IEN[0]]._coor,_vertex_list[_piece_list[i]._IEN[1]]._coor,_vertex_list[_piece_list[i]._IEN[2]]._coor,"zy")*s*nz;
					int_x.z=int_x.z+gauss(_vertex_list[_piece_list[i]._IEN[0]]._coor,_vertex_list[_piece_list[i]._IEN[1]]._coor,_vertex_list[_piece_list[i]._IEN[2]]._coor,"zz")*s*nz*0.5;
				}
				_piece_list[i]._is_visit=1;
			} 
			else if (_piece_list[i]._type[0]==1)
			{
				//This piece is above the cutting plane and doing nothing
				_piece_list[i]._is_visit=1;
			}
			else if (_piece_list[i]._type[0]==0)
			{
				//The cutting plane is across this piece and one ring will be generated
				vec3D_int pre_point,next_point;
				//-----------------------------------------------------------------------------------------------------
				//vec3D_int::_coor represents the coordiante of the point in the cutting ring
				//vec3D_int:_type represents the state of this point i.e.
				//_type>=0:this point is degenerate and it is coincident with _vertex_list[_type]
				//_type=-1,-2,-3:this point is not degenerate and locates on IEN[abs(type+1)] edge of the next piece
				//------------------------------------------------------------------------------------------------------
				int next_piece_id,current_piece_id,start_piece_id;
				calculate_initial_cutting_point(i,pre_point,next_piece_id);
				_cutting_ring.add_point(pre_point._coor);
				start_piece_id=next_piece_id;current_piece_id=start_piece_id;next_piece_id=-1;
				while (next_piece_id!=start_piece_id)   //If next_piece_id==start_piece_id this ring is closed
				{
					calculate_next_cutting_point(current_piece_id,pre_point,next_piece_id,next_point);			
					_piece_list[current_piece_id]._is_visit=1;
					if (next_piece_id!=start_piece_id)
					{
						//next_point is not the initial point, so it will be added to _cutting_ring
						_cutting_ring.add_point(next_point._coor);
					}
					calculate_integral_on_subpiece(pre_point,next_point,current_piece_id,sub_int_1,sub_int_x,type);
					int_1=int_1+sub_int_1;
					if (type=="centroid")
					{
						int_x=int_x+sub_int_x;
					}
					current_piece_id=next_piece_id;pre_point=next_point;
				}
				_cutting_ring.final_one_ring();
			}
		}

	}
	//The contribution of the surface integral on cutting plane will be added at last
	_cutting_ring.calcualte_volume_integral(n,cutting_area,cutting_plane_int_1,cutting_plane_int_x,type);
	int_1=int_1+cutting_plane_int_1;
	if (type=="centroid")
	{
		int_x=int_x+cutting_plane_int_x;
	}
	for (int i=0;i<_nump;i++)
	{
		_piece_list[i]._is_visit=0;
	}
	if (type=="volume")
	{
		result.value(int_1,cutting_area,0);
	}
	else if (type=="centroid")
	{
		if (int_1<1e-21)
		{
			cout<<"int_1<1e-21"<<endl;
		}
		result=int_x/int_1;
	}
	return result;
}
void Tpolyhedron::cut_polyhedron_by_plane(vec3D &n,double d,Tpolyhedron* poly_below,Tpolyhedron* poly_above)
{
	//The cutting process is below:
	//1.Calculate all vertexes including origin vertexes and generated vertexes on the cutting plane
	//2.Allocate the IEN of pieces into poly_below and poly_above (IEN is the global id of origin polyhedron)
	//3.Allocate the vertexes into poly_below and poly_above
	//4.Compress the vertexes' id in poly_below and poly_above
	//The process is similar to calculate the integral below the cutting plane
	classify_pieces(n,d);
	vector<vec3D> new_point_list;       //Newly generated points in the cutting process (in-degenerated points in the cutting ring)
	vector<int>   id_in_ring;           //Global vertex id of points in the cutting ring
	             
	for (int i=0;i<_nump;i++)
	{
		if (_piece_list[i]._is_visit==0)
		{
			if (_piece_list[i]._type[0]==2)
			{
				//The piece is below the cutting plane so add this piece into below
				//The adding vertexes are in global id of the origin polyhedron
				poly_below->add_piece(_piece_list[i]._IEN[0],_piece_list[i]._IEN[1],_piece_list[i]._IEN[2]);
				_piece_list[i]._is_visit=1;
			} 
			else if (_piece_list[i]._type[0]==1)
			{
				//The piece is above the cutting plane so add this piece into above
				poly_above->add_piece(_piece_list[i]._IEN[0],_piece_list[i]._IEN[1],_piece_list[i]._IEN[2]);
				_piece_list[i]._is_visit=1;
			}
			else if (_piece_list[i]._type[0]==0)
			{
				int ring_start_position=id_in_ring.size();          //The start point position id of this ring in cutting_id
				vec3D_int pre_point,next_point;int next_piece_id,start_piece_id,current_piece_id;
				_cutting_ring._length=0;            //The data of previous rings need not to be saved, 
				                                    //so _cutting_ring in this process only contains one ring
				calculate_initial_cutting_point(i,pre_point,next_piece_id);
				classify_new_point(id_in_ring,new_point_list,pre_point);
				_cutting_ring.add_point(pre_point._coor);
				start_piece_id=next_piece_id;current_piece_id=start_piece_id;next_piece_id=-1;
				int id_pre,id_next;                 //The global id of pre_point and next_point
				int order;                          //The order of the cutting ring
				while (next_piece_id!=start_piece_id)                
				{
					calculate_next_cutting_point(current_piece_id,pre_point,next_piece_id,next_point);				
					_piece_list[current_piece_id]._is_visit=1;
					if (next_piece_id!=start_piece_id)
					{
						classify_new_point(id_in_ring,new_point_list,next_point);    
						_cutting_ring.add_point(next_point._coor);
						//id_pre and id_next is the last two variables in id_in_ring
						id_pre=id_in_ring[id_in_ring.size()-2];
						id_next=id_in_ring[id_in_ring.size()-1];
					}
					else
					{
						//The cutting ring is closed so no new point is added
						//The id_pre is the last element in cutting_id and id_next is the first point in this ring
						id_pre=id_in_ring[id_in_ring.size()-1];
						id_next=id_in_ring[ring_start_position];
					}
					allocate_subpieces(pre_point,next_point,current_piece_id,id_pre,id_next,poly_below,poly_above);
					current_piece_id=next_piece_id;pre_point=next_point;
				}
				vec3D_int ring_centroid;                                    //The centroid of this ring
				_cutting_ring.determine_ring_order(n,order,ring_centroid);  //Determine the order of this ring and calcualte the ring centroid
				classify_new_point(id_in_ring,new_point_list,ring_centroid);     //Add the centroid of cutting ring into cutting_id and new_point
				allocate_pieces_in_ring(ring_start_position,id_in_ring,poly_below,poly_above,order);      //Allocate the cutting plane into below and above polyhedron
			}
			//The cycle will be continued because there might be more than one cutting rings if the polygon is not convex
		}
	}
	//The new polyhedrons' pieces with id in ordinary polyhedron point is complete
	//The final task is allocating the vertexes and compressing the vertex id
	compress_vertex_id(id_in_ring,new_point_list,poly_below);
	compress_vertex_id(id_in_ring,new_point_list,poly_above);
	for (int i=0;i<_nump;i++)
	{
		_piece_list[i]._is_visit=0;
	}
	return;
}
void Tpolyhedron::calculate_piece_topology()
{
	int sumi[3],subi[3],sumj[3],subj[3];
	for (int i=0;i<_nump;i++)
	{
		for (int k=0;k<3;k++)
		{
			sumi[k]=_piece_list[i]._IEN[k%3]+_piece_list[i]._IEN[(k+1)%3];
			subi[k]=abs(_piece_list[i]._IEN[k%3]-_piece_list[i]._IEN[(k+1)%3]);
		}
		for (int j=0;j<_nump;j++)
		{
			if (j!=i)
			{
				for (int k=0;k<3;k++)
				{
					sumj[k]=_piece_list[j]._IEN[k%3]+_piece_list[j]._IEN[(k+1)%3];
					subj[k]=abs(_piece_list[j]._IEN[k%3]-_piece_list[j]._IEN[(k+1)%3]);
				}
				for (int k=0;k<3;k++)
				{
					for (int l=0;l<3;l++)
					{
						if (sumi[k]==sumj[l] && subi[k]==subj[l])
						{
							_piece_list[i]._adjacent[k]=j;
							//If the sum and abs of difference of global number on two edge are identical,
							//they are coincident and the two pieces are adjacent.
							goto breakLoop;
						}
					}
				}
				breakLoop: ;
			}
		}
	}
}
void Tpolyhedron::calculate_sign_distance(vec3D &n,double d)
{
	for (int i=0;i<_numv;i++)
	{
		_vertex_list[i]._variable=_vertex_list[i]._coor*n-d;
	}
	return;
}
void Tpolyhedron::calculate_altitude_range(vec3D &n,double &d_min,double &d_max)
{
	d_max=-100000;d_min=100000;
	double temp;
	for (int i=0;i<_numv;i++)
	{
		temp=_vertex_list[i]._coor*n;
		if (temp>d_max)
		{
			d_max=temp;
		}
		if (temp<d_min)
		{
			d_min=temp;
		}
	}
	return;
}
void Tpolyhedron::classify_pieces(vec3D &n,double d)
{
	calculate_sign_distance(n,d);
	double sign_distance[3];   //The sign distance of the vertexes of pieces
	for (int i=0;i<_nump;i++)
	{
		for (int j=0;j<3;j++)
		{
			sign_distance[j]=_vertex_list[_piece_list[i]._IEN[j]]._variable;
		}
		//Determine the first component of piece type: whether the cutting plane is across the piece
		if (sign_distance[0]>=0 && sign_distance[1]>=0 && sign_distance[2]>=0)
		{
			//This piece is above the cutting plane
			_piece_list[i]._type[0]=1;
		} 
		else if (sign_distance[0]<=0 && sign_distance[1]<=0 && sign_distance[2]<=0)
		{
			//This piece is below the cutting plane
			_piece_list[i]._type[0]=2;
		} 
		else
		{
			//The cutting plane is across this piece
			_piece_list[i]._type[0]=0;
		}
		//Determine the second component of piece type: whether the piece is degenerate
		if (sign_distance[0]*sign_distance[1]*sign_distance[2]!=0)
		{
			//This piece is not degenerated
			_piece_list[i]._type[1]=0;
		} 
		else if (sign_distance[0]==0 && sign_distance[1]*sign_distance[2]!=0)
		{
			//This piece is degenerate and the cutting only cross IEN[0]
			_piece_list[i]._type[1]=1;
		} 
		else if (sign_distance[1]==0 && sign_distance[0]*sign_distance[2]!=0)
		{
			//This piece is degenerate and the cutting only cross IEN[1]
			_piece_list[i]._type[1]=2;
		} 
		else if (sign_distance[2]==0 && sign_distance[0]*sign_distance[1]!=0)
		{
			//This piece is degenerate and the cutting only cross IEN[2]
			_piece_list[i]._type[1]=3;
		} 
		else if (sign_distance[2]!=0)
		{
			//The cutting plane is across IEN[0] and IEN[1] and we assume this situation also as degenerate case
			//although this piece is above or below the cutting plane indeed
			_piece_list[i]._type[1]=4;
			//0-1 is degenerated
			//In order to distingish this kind of above/below, the type[0] of this piece should be modified to 0
			_piece_list[i]._type[0]=0;      
			                                
		} 
		else if (sign_distance[0]!=0)
		{
			//The cutting plane is across IEN[1] and IEN[2]
			_piece_list[i]._type[1]=5;
			//1-2 is degenerated
			_piece_list[i]._type[0]=0;
		} 
		else if (sign_distance[1]!=0)
		{
			//The cutting plane is across IEN[0] and IEN[2]
			_piece_list[i]._type[1]=6;
			//2-0 is degenerate
			_piece_list[i]._type[0]=0;
		} 
		else
		{
			cout<<"Error: The cutting plane is coincident with a piece in calling classify_pieces()"<<endl;
			system("Pause");
			exit(0);
		}
	}
	return;
}
void Tpolyhedron::calculate_initial_cutting_point(int this_piece_id,vec3D_int &initial_point,int &next_piece_id)
{
	double sign_distance[3]; 
	Spiece* this_piece=&_piece_list[this_piece_id];
	for (int j=0;j<3;j++)
	{
		sign_distance[j]=_vertex_list[this_piece->_IEN[j]]._variable;
	}
	int type=this_piece->_type[1];
	if (type==0)
	{
		//This initial piece is not degenerate
		for (int j=0;j<3;j++)
		{
			if (sign_distance[j%3]*sign_distance[(j+1)%3]<0)
			{
				initial_point._coor=interpolate(_vertex_list[this_piece->_IEN[j%3]],_vertex_list[this_piece->_IEN[(j+1)%3]]);
				next_piece_id=this_piece->_adjacent[j];
				initial_point._type=-calculate_common_edge_id(next_piece_id,this_piece_id,j)-1;
				return;
			}
		}
		cout<<"Error:conflict in calling calculate_initial_cutting_point()"<<endl;
		cout<<sign_distance[0]<<" "<<sign_distance[1]<<" "<<sign_distance[2]<<endl;
		for (int j=0;j<3;j++)
		{
			cout<<_vertex_list[this_piece->_IEN[j]]._coor.x<<" "<<_vertex_list[this_piece->_IEN[j]]._coor.y<<" "<<_vertex_list[this_piece->_IEN[j]]._coor.z<<endl;;
		}
		system("Pause");
		exit(0);
	} 
	else if (type==1||type==2||type==3)	
	{
		//The initial piece is one point degenerate and coincident with _vertex_list[_piece_list[_IEN[type-1]]]
		initial_point._type=this_piece->_IEN[type-1];
		initial_point._coor=_vertex_list[initial_point._type]._coor;
		next_piece_id=find_certain_piece_id(this_piece_id,this_piece_id,initial_point._type);
		//The next_piece_id of this cutting point must fulfill
		//1.not equal to this_piece_id
		//2.is across by the cutting plane i.e. _type[0]=0
		//3.has a vertex id equals to initial_point._type
	} 
	else if (type==4||type==5||type==6)
	{
		//The initial piece is two point degenerate
		initial_point._type=this_piece->_IEN[type-4];
		initial_point._coor=_vertex_list[initial_point._type]._coor;
		next_piece_id=find_certain_piece_id(this_piece_id,this_piece->_adjacent[type-4],initial_point._type);
		//The next_piece_id of this cutting point must fulfill
		//1.not equal to this_piece_id and the adjacent piece of the degenerated edge i.e. this_piece->_adjacent[type-4]
		//2.is across by the cutting plane i.e. _type[0]=0
		//3.has a vertex id equals to initial_point._type
	}
	return;
}
void Tpolyhedron::calculate_next_cutting_point(int current_piece_id,vec3D_int &pre_point,int &next_piece_id,vec3D_int &next_point)
{
	double sign_distance[3]; 
	Spiece* this_piece=&_piece_list[current_piece_id];
	for (int j=0;j<3;j++)
	{
		sign_distance[j]=_vertex_list[this_piece->_IEN[j]]._variable;
	}
	int type=this_piece->_type[1];
	if (type==0)
	{
		//If The current piece is not degenerate, the pre and next points must be not degenerate
		for (int j=0;j<3;j++)
		{
			if (j!=abs(pre_point._type+1))   //The _type of the pre_point is the its located edge id of current_piece
			{
				if (sign_distance[j%3]*sign_distance[(j+1)%3]<0)
				{
					next_point._coor=interpolate(_vertex_list[this_piece->_IEN[j%3]],_vertex_list[this_piece->_IEN[(j+1)%3]]);
					next_piece_id=this_piece->_adjacent[j];
					next_point._type=-calculate_common_edge_id(next_piece_id,current_piece_id,j)-1;
					return;
				}
			}
		}
		cout<<"Error:conflict in calling calcualte_next_cutting_point()"<<endl;
		system("Pause");
		exit(0);
	} 
	else if (type==1||type==2||type==3)	
	{
		//The current piece is one point degenerate
		if (pre_point._type<0)		
		{
			//Pre_point is not degenerate so the next point must be degenerate
			next_point._type=this_piece->_IEN[type-1];
			next_point._coor=_vertex_list[next_point._type]._coor;
			next_piece_id=find_certain_piece_id(current_piece_id,current_piece_id,next_point._type);
		} 
		else	
		{
			//Pre_point is degenerate so the next point must be not degenerate
			//As the defination of _type[1], vertex IEN[type-1] is degenerated so edge IEN[type%3]-IEN[(type+1)%3] cross the cutting plane 
			next_point._coor=interpolate(_vertex_list[this_piece->_IEN[type%3]],_vertex_list[this_piece->_IEN[(type+1)%3]]);
			int edge_id=type%3; //The local edge id of the across edge in current_piece
			next_piece_id=this_piece->_adjacent[edge_id];
			next_point._type=-calculate_common_edge_id(next_piece_id,current_piece_id,edge_id)-1;
		}
	} 
	else if (type==4||type==5||type==6)
	{
		//Current piece is two point degenerate so pre_point and next_point must be degenerated
		//Compare the global id of pre_point with _IEN[type-4] and _IEN[(type-3)%3 to obtain the next point
		if (this_piece->_IEN[type-4]==pre_point._type)
		{
			next_point._type=this_piece->_IEN[(type-3)%3];
		} 
		else
		{
			next_point._type=this_piece->_IEN[type-4];
		}
		next_point._coor=_vertex_list[next_point._type]._coor;
		next_piece_id=find_certain_piece_id(current_piece_id,this_piece->_adjacent[type-4],next_point._type);
	}
	return;
}
void Tpolyhedron::calculate_integral_on_subpiece(vec3D_int &pre_point,vec3D_int &next_point,int &current_piece_id,double &int_1,vec3D &int_x,string type)
{
	int_1=0;int_x.value(0,0,0);
	vec3D point_list[6];
	int num_point;
	calculate_subpiece_point_below(pre_point,next_point,current_piece_id,point_list,num_point);
	double s,nz;
	//subpiece point_list[0],point_list[1],point_list[2] must exist
	area(point_list[0],point_list[1],point_list[2],s,nz);
	int_1=gauss(point_list[0],point_list[1],point_list[2],"z")*s*nz;
	if (type=="centroid")
	{
		int_x=gauss_for_centroid(point_list[0],point_list[1],point_list[2])*s*nz;
	}
	if (num_point==4)
	{
		area(point_list[0],point_list[2],point_list[3],s,nz);
		int_1=int_1+gauss(point_list[0],point_list[2],point_list[3])*s*nz;
		if (type=="centroid")
		{
			int_x=int_x+gauss_for_centroid(point_list[0],point_list[2],point_list[3])*s*nz;
		}
	}
	else if (num_point==6)
	{
		area(point_list[3],point_list[4],point_list[5],s,nz);
		int_1=int_1+gauss(point_list[0],point_list[2],point_list[3])*s*nz;
		if (type=="centroid")
		{
			int_x=int_x+gauss_for_centroid(point_list[0],point_list[2],point_list[3])*s*nz;
		}
	}
	return;
}
void Tpolyhedron::calculate_subpiece_point_below(vec3D_int &pre_point,vec3D_int &next_point,int &current_piece_id,vec3D (&point_list)[6],int &num_point)
{
	double sign_distance[3]; 
	Spiece* this_piece=&_piece_list[current_piece_id];
	for (int j=0;j<3;j++)
	{
		sign_distance[j]=_vertex_list[this_piece->_IEN[j]]._variable;
	}	
	int piece_type=this_piece->_type[1];
	if (piece_type==0)
	{
		//Current piece is not degenerate
		int next=(-pre_point._type)%3;   //pre_point->IEN[next] is the along outward order   
		if (sign_distance[0]*sign_distance[1]*sign_distance[2]>0)
		{
			//Two points are below the cutting plane in current piece, the integral on the sub-piece will be divided into 
			//the integral on two triangles
			if (sign_distance[next]<0)		
			{
				point_list[0]=pre_point._coor;point_list[1]=_vertex_list[this_piece->_IEN[next]]._coor;
				point_list[2]=_vertex_list[this_piece->_IEN[(next+1)%3]]._coor;
				point_list[3]=next_point._coor;
				num_point=4;
				return;
			} 
			else	
			{
				//Next point is not below:
				point_list[0]=pre_point._coor;point_list[1]=next_point._coor;
				point_list[2]=_vertex_list[this_piece->_IEN[(next+1)%3]]._coor;
				point_list[3]=_vertex_list[this_piece->_IEN[(next+2)%3]]._coor;
				num_point=4;
				return;
			}
		} 
		else	
		{
			//One point is below the cutting plane in current piece
			if (sign_distance[next]<0)
			{
				//Next point is below: represent one order
				point_list[0]=pre_point._coor;
				point_list[1]=_vertex_list[this_piece->_IEN[next]]._coor;
				point_list[2]=next_point._coor;
				num_point=3;
				return;
			} 
			else	
			{
				//Next point is not below
				point_list[0]=pre_point._coor;
				point_list[1]=next_point._coor;
				point_list[2]=_vertex_list[this_piece->_IEN[(next+2)%3]]._coor;
				num_point=3;
				return;
			}
		}
	} 
	else if (piece_type==1||piece_type==2||piece_type==3)	
	{
		//Current piece is one point degenerate
		if (pre_point._type>=0)		
		{
			//Pre_point is degenerate and next_point is not degenerate
			int l=find_vertex_local_id(current_piece_id,pre_point._type);    //The local id of pre_point
			                                                                 //pre_point is coincident with a vertex so it has a local id in piece_id)
			if (sign_distance[(l+1)%3]>0)
			{
				//The next point is not below:
				point_list[0]=pre_point._coor;
				point_list[1]=next_point._coor;
				point_list[2]=_vertex_list[this_piece->_IEN[(l+2)%3]]._coor;
				num_point=3;
				return;
			} 
			else	
			{
				//The next point is below
				point_list[0]=pre_point._coor;
				point_list[1]=_vertex_list[this_piece->_IEN[(l+1)%3]]._coor;
				point_list[2]=next_point._coor;
				num_point=3;
				return;
			}
		}
		else	
		{

			//Pre_point is not degenerate and next_point is degenerate
			int next=(-pre_point._type)%3;     //The next point local id of pre_point with outward order
			if (sign_distance[next]<0)
			{
				//Next point is below
				point_list[0]=pre_point._coor;
				point_list[1]=_vertex_list[this_piece->_IEN[next]]._coor;
				point_list[2]=next_point._coor;
				num_point=3;
				return;
			} 
			else	
			{
				//Next point is not below: represent the other order
				point_list[0]=pre_point._coor;
				point_list[1]=next_point._coor;
				point_list[2]=_vertex_list[this_piece->_IEN[(next+2)%3]]._coor;
				num_point=3;
				return;
			}
		}
	}
	else
	{
		int num=0;   //The number of point from current piece (equals 0 or 3)
		//Current piece is two-point degenerate so the piece must entirely below or above the cutting plane
		if (sign_distance[0]+sign_distance[1]+sign_distance[2]<0)	
		{
			//The piece is below the cutting plane
			point_list[0]=_vertex_list[this_piece->_IEN[0]]._coor;
			point_list[1]=_vertex_list[this_piece->_IEN[1]]._coor;
			point_list[2]=_vertex_list[this_piece->_IEN[2]]._coor;
			num=3;
		} 
		//The adjacent piece must be taken into consideration because the cutting plane also cross this piece
		//and it will not be visited in this ring
		int adjacent=this_piece->_adjacent[piece_type-4];
		_piece_list[adjacent]._is_visit=1;    //Mark the adjacent piece visited
		for (int j=0;j<3;j++)
		{
			sign_distance[j]=_vertex_list[_piece_list[adjacent]._IEN[j]]._variable;
		}
		if (sign_distance[0]+sign_distance[1]+sign_distance[2]<0)
		{
			//The adjacent piece is entirely below the cutting plane
			point_list[num]=_vertex_list[_piece_list[adjacent]._IEN[0]]._coor;
			point_list[num+1]=_vertex_list[_piece_list[adjacent]._IEN[1]]._coor;
			point_list[num+2]=_vertex_list[_piece_list[adjacent]._IEN[2]]._coor;
			num_point=num+3;
		}
		return;
	}
}
void Tpolyhedron::allocate_subpieces(vec3D_int &pre_point,vec3D_int &next_point,int current_piece_id,int id_pre,int id_next,Tpolyhedron* poly_below,Tpolyhedron* poly_above)
{
	//Allcoate the below and above piece IEN to poly_below and poly_above
	//This process is similar to calculate_subpiece_point_below()
	Spiece* this_piece=&_piece_list[current_piece_id];
	int piece_type=this_piece->_type[1];
	double sign_distance[3];
	for (int j=0;j<3;j++)
	{
		sign_distance[j]=_vertex_list[this_piece->_IEN[j]]._variable;
	}
	if (piece_type==0)
	{
		int next=(-pre_point._type)%3;
		if (sign_distance[0]*sign_distance[1]*sign_distance[2]>0)
		{
			if (sign_distance[next]<0)
			{
				poly_below->add_piece(id_pre,this_piece->_IEN[next],this_piece->_IEN[(next+1)%3]);
				poly_below->add_piece(id_pre,this_piece->_IEN[(next+1)%3],id_next);
				poly_above->add_piece(id_pre,id_next,this_piece->_IEN[(next+2)%3]);
			} 
			else
			{
				poly_below->add_piece(id_pre,id_next,this_piece->_IEN[(next+1)%3]);
				poly_below->add_piece(id_pre,this_piece->_IEN[(next+1)%3],this_piece->_IEN[(next+2)%3]);
				poly_above->add_piece(id_pre,this_piece->_IEN[next],id_next);
			}
		} 
		else
		{
			if (sign_distance[next]<0)
			{
				poly_below->add_piece(id_pre,this_piece->_IEN[next],id_next);
				poly_above->add_piece(id_pre,id_next,this_piece->_IEN[(next+1)%3]);
				poly_above->add_piece(id_pre,this_piece->_IEN[(next+1)%3],this_piece->_IEN[(next+2)%3]);
			} 
			else
			{
				poly_below->add_piece(id_pre,id_next,this_piece->_IEN[(next+2)%3]);
				poly_above->add_piece(id_pre,this_piece->_IEN[next],this_piece->_IEN[(next+1)%3]);
				poly_above->add_piece(id_pre,this_piece->_IEN[(next+1)%3],id_next);
			}
		}
	} 
	else if (piece_type==1||piece_type==2||piece_type==3)
	{
		if (pre_point._type>=0)
		{
			int l=find_vertex_local_id(current_piece_id,pre_point._type);
			if (sign_distance[(l+1)%3]<0)
			{
				poly_below->add_piece(id_pre,this_piece->_IEN[(l+1)%3],id_next);
				poly_above->add_piece(id_pre,id_next,this_piece->_IEN[(l+2)%3]);
			} 
			else
			{
				poly_below->add_piece(id_pre,id_next,this_piece->_IEN[(l+2)%3]);
				poly_above->add_piece(id_pre,this_piece->_IEN[(l+1)%3],id_next);
			}
		}
		else
		{
			int next=(-pre_point._type)%3;
			if (sign_distance[next]<0)
			{
				poly_below->add_piece(id_pre,this_piece->_IEN[next],id_next);
				poly_above->add_piece(id_pre,id_next,this_piece->_IEN[(next+2)%3]);
			} 
			else
			{
				poly_below->add_piece(id_pre,id_next,this_piece->_IEN[(next+2)%3]);
				poly_above->add_piece(id_pre,this_piece->_IEN[next],id_next);
			}
		}
	}
	else
	{
		if (sign_distance[0]+sign_distance[1]+sign_distance[2]<0)
		{
			poly_below->add_piece(this_piece->_IEN[0],this_piece->_IEN[1],this_piece->_IEN[2]);
		} 
		else
		{
			poly_above->add_piece(this_piece->_IEN[0],this_piece->_IEN[1],this_piece->_IEN[2]);
		}
		int adjacent=this_piece->_adjacent[piece_type-4];
		_piece_list[adjacent]._is_visit=1;
		for (int j=0;j<3;j++)
		{
			sign_distance[j]=_vertex_list[_piece_list[adjacent]._IEN[j]]._variable;
		}
		if (sign_distance[0]+sign_distance[1]+sign_distance[2]<0)
		{
			poly_below->add_piece(_piece_list[adjacent]._IEN[0],_piece_list[adjacent]._IEN[1],_piece_list[adjacent]._IEN[2]);
	
		}
		else
		{
			poly_above->add_piece(_piece_list[adjacent]._IEN[0],_piece_list[adjacent]._IEN[1],_piece_list[adjacent]._IEN[2]);
		}
	}
	return;
}
void Tpolyhedron::allocate_pieces_in_ring(int ring_start_position,vector<int> &id_in_ring,Tpolyhedron* poly_below,Tpolyhedron* poly_above,int order)
{
	int final_position=id_in_ring.size();     //The final position of this ring is the last element in id_in_ring
	                                          //The last element is the ring centroid (position final_position-1)
	                                          //The position of the last element in the ring is final_position-2
	if (order==0)
	{
		for (int i=ring_start_position;i<final_position-2;i++)
		{
			poly_below->add_piece(id_in_ring[final_position-1],id_in_ring[i],id_in_ring[i+1]);
			poly_below->_surface_piece_list.push_back(poly_below->_piece_list.size()-1);   //The pieces on cutting ring is the material surface
			poly_above->add_piece(id_in_ring[final_position-1],id_in_ring[i+1],id_in_ring[i]);
			poly_above->_surface_piece_list.push_back(poly_above->_piece_list.size()-1);
		}
		poly_below->add_piece(id_in_ring[final_position-1],id_in_ring[final_position-2],id_in_ring[ring_start_position]);
		poly_below->_surface_piece_list.push_back(poly_below->_piece_list.size()-1);
		poly_above->add_piece(id_in_ring[final_position-1],id_in_ring[ring_start_position],id_in_ring[final_position-2]);
		poly_above->_surface_piece_list.push_back(poly_above->_piece_list.size()-1);
	}
	else
	{
		for (int i=ring_start_position;i<final_position-2;i++)
		{
			poly_below->add_piece(id_in_ring[final_position-1],id_in_ring[i+1],id_in_ring[i]);
			poly_below->_surface_piece_list.push_back(poly_below->_piece_list.size()-1);   //The pieces on cutting ring is the material surface
			poly_above->add_piece(id_in_ring[final_position-1],id_in_ring[i],id_in_ring[i+1]);
			poly_above->_surface_piece_list.push_back(poly_above->_piece_list.size()-1);
		}
		poly_below->add_piece(id_in_ring[final_position-1],id_in_ring[ring_start_position],id_in_ring[final_position-2]);
		poly_below->_surface_piece_list.push_back(poly_below->_piece_list.size()-1);
		poly_above->add_piece(id_in_ring[final_position-1],id_in_ring[final_position-2],id_in_ring[ring_start_position]);
		poly_above->_surface_piece_list.push_back(poly_above->_piece_list.size()-1);
	}
	return;
}
void Tpolyhedron::classify_new_point(vector<int> &id_in_ring,vector<vec3D> &new_point_list,vec3D_int &new_point_in_ring)
{
	if (new_point_in_ring._type>=0)
	{
		id_in_ring.push_back(new_point_in_ring._type);      //This new cutting point is degenerate so the vertex id is its global id      
	}
	else
	{
		//This new cutting point is not degenerated so it's a new point 
		int num=new_point_list.size();
		new_point_list.push_back(new_point_in_ring._coor);    //The cutting point is a new point
		id_in_ring.push_back(_numv+num);                 //The new created point's id is begin with nump
	}
	return;
}
void Tpolyhedron::compress_vertex_id(vector<int> &id_in_ring,vector<vec3D> &new_point_list,Tpolyhedron* poly_ptr)
{
	vector<int> new2old;    //i th point in created polyhedron p_ploy is new2old[i] in origin polyhedron
	                        //Point in in created polyhedron p_ploy is new2old.size()
	
	int old_id;
	int max_vertex_id_in_old=0;           //maximum vertex id in origin created polyhedron
	//Calculate all vertexes in generated polyhedron and its id in the origin polyhedron
	for (int i=0;i<poly_ptr->_nump;i++)
	{
		for (int j=0;j<3;j++)
		{
			old_id=poly_ptr->_piece_list[i]._IEN[j];
			if (old_id>max_vertex_id_in_old)
			{
				max_vertex_id_in_old=old_id;
			}
			if (!is_exist(old_id,new2old))
			{
				new2old.push_back(old_id);
			}
		}
	}
	//Calculate the IEN in generated polyehdron 
	int* old2new;                             //id in origin polyhedron is old2new[i] in created polyhedron
	old2new=new int[max_vertex_id_in_old+1];  //old2new[max_vertex_id_in_old_id] is exist,so the length of old2new is max_vertex_id_in_old_id+1
	for (int i=0;i<new2old.size();i++)
	{
		old2new[new2old[i]]=i;
	}	
	for (int i=0;i<poly_ptr->_nump;i++)
	{
		for (int j=0;j<3;j++)
		{
			old_id=poly_ptr->_piece_list[i]._IEN[j];
			poly_ptr->_piece_list[i]._IEN[j]=old2new[old_id];
		}
	}
	delete old2new;
	//Finally obtain the vertexes of the generated polyhedron
	for (int i=0;i<new2old.size();i++)
	{
		old_id=new2old[i];
		if (old_id<_numv)           
		{
			//This point is origin point
			poly_ptr->add_vertex(_vertex_list[old_id]._coor);
		}
		else
		{
			//This point is not the origin point,it will comes from new_point_list
			poly_ptr->add_vertex(new_point_list[old_id-_numv]);
		}
	}
	return;
}
int Tpolyhedron::calculate_common_edge_id(int &piece_id1,int &piece_id2,int &j)
{
	int sum0=_piece_list[piece_id2]._IEN[j%3]+_piece_list[piece_id2]._IEN[(j+1)%3];
	int sub0=abs(_piece_list[piece_id2]._IEN[j%3]-_piece_list[piece_id2]._IEN[(j+1)%3]);
	for (int i=0;i<3;i++)
	{
		int sumi=_piece_list[piece_id1]._IEN[i%3]+_piece_list[piece_id1]._IEN[(i+1)%3];
		int subi=abs(_piece_list[piece_id1]._IEN[i%3]-_piece_list[piece_id1]._IEN[(i+1)%3]);
		if (sumi==sum0 && subi==sub0)
		{
			return i;
		}
	}
	cout<<"Error: Failed to find the approperate edge id in calling calculate_common_edge_id()"<<endl;
	system("Pasue");
	exit(0);
}
int Tpolyhedron::find_certain_piece_id(int &piece_id1,int &piece_id2,int &_vertex_id)
{
	for (int i=0;i<_nump;i++)
	{
		if (i!=piece_id1 && i!=piece_id2 && _piece_list[i]._type[0]==0)
		{
			for (int j=0;j<3;j++)
			{
				if (_piece_list[i]._IEN[j]==_vertex_id)
				{
					return i;
				}
			}
		}
	}
	cout<<"Error: Failed to find the next piece of a degenerated cutting point in calling find_certain_piece_id()"<<endl;
	system("Pause");
	exit(0);
}
int Tpolyhedron::find_vertex_local_id(int &piece_id,int &vertex_id)
{
	for (int i=0;i<3;i++)
	{
		if (_piece_list[piece_id]._IEN[i]==vertex_id)
		{
			return i;
		}
	}
	cout<<"Error: there is no vertex "<<vertex_id<<" in piece "<<piece_id<<" in calling find_vertex_local_id()"<<endl;
	system("Pause");
	exit(0);
}