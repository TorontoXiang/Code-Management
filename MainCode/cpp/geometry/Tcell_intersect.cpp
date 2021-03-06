#include "Tcell_intersect.h"
#include "Shape_function.h"
void Tcell_intersect::intersect()
{
	double volume,ix,iy,iz;                      //Intersected volume and its integrals for every tetrahedron sub-cells in new cell and material polyhedron in old cell
	vec3D int_x;
	double ro_old,roe_old;                       
	vec3D gradient;
	double mass,energy;                          
	double av_density,av_pressure,av_soundspeed;
	vec3D material_centroid;
	vec3D intersect_centroid,coordinate_std;     //The centroid of the intersection portion and its corresponding standard coordinate
	vec3D vel,moment;                            //The interpolated velocity at the centroid of the intersection portion and its moment
	bool discontinue=whether_dis(120);
	if (discontinue)
	{
		create_corner_polyhedron();
		_old_cell->calculate_average_variable(av_density,av_pressure,av_soundspeed);
	}
	int num_material=_old_cell->_gausspoint.size();
	for (int i=0;i<48;i++)
	{
		for (int j = 0; j < num_material; j++)
		{
			if (_old_cell->G_material_polyhedron(j)!=NULL)
			{
				//Calculate the integral on intersection portion
				_tet_list[i].intersect_with_polyhedron(_old_cell->G_material_polyhedron(j),volume,ix,iy,iz);
				int_x.value(ix,iy,iz);
				//Density,energy and cell volume accumulation
				material_centroid=_old_cell->G_material_centroid(j);
				ro_old=_old_cell->_gausspoint[j]->G_density();gradient=_old_cell->G_g_ro(j);
				mass=(ro_old-material_centroid*gradient)*volume+gradient.x*ix+gradient.y*iy+gradient.z*iz;
				roe_old=_old_cell->_gausspoint[j]->G_roe();gradient=_old_cell->G_g_roe(j);
				energy=(roe_old-material_centroid*gradient)*volume+gradient.x*ix+gradient.y*iy+gradient.z*iz;
				_new_cell->accumulate_variables(j,volume,int_x,mass,energy);
				_new_cell->_volume_intersection=_new_cell->_volume_intersection+volume;

				//Nodal mass and moment accumulation
				//Stored in corner_mass and corner_moment for the convenience of the OpenMP parallel
				if (volume>1e-15)                                     
				{
					//Accumulate the moment when intersection portion exists
					if (!discontinue)
					{
						//If the velocity in the old cell is continuous
						//The intersection portion moment is generated by shape function interpolation
						_new_cell->_corner_mass[_ass_node[i]]=_new_cell->_corner_mass[_ass_node[i]]+mass;
						intersect_centroid.value(ix/volume,iy/volume,iz/volume);
						coordinate_std=_old_cell->calculate_std_coordinate(intersect_centroid);
						vel=_old_cell->shape_function_interpolate(coordinate_std,"Velocity");
						moment=vel*mass;
						_new_cell->_corner_moment[_ass_node[i]]=_new_cell->_corner_moment[_ass_node[i]]+moment;
					}
				}
			}
			if (discontinue)
			{
				//If the velocity in the old cell is discontinuous
				//The intersection portion of sub-cell in cell base and corner-cell in cell change have to be calculated		
				double ivolume;
				for (int k=0;k<8;k++)
				{
					_tet_list[i].intersect_with_polyhedron(&_corner_polyhedron[k],ivolume,ix,iy,iz);
					vel=_old_cell->_node_ptr[k]->_velocity;
					mass=ivolume*av_density;
					moment=vel*mass;
					_new_cell->_corner_mass[_ass_node[i]]=_new_cell->_corner_mass[_ass_node[i]]+mass;
					_new_cell->_corner_moment[_ass_node[i]]=_new_cell->_corner_moment[_ass_node[i]]+moment;
				}
			}
		}
	}
	//Clear the corner polyhedron on the old cell
	for (int i=0;i<8;i++)
	{
		_corner_polyhedron[i].clear_the_polyhedron();
	}
	return;
}
void Tcell_intersect::overlapping(int mid)
{
	double volume,ix,iy,iz;                      //Intersected volume and its integrals for every tetrahedron sub-cells in new cell and material polyhedron in old cell
	vec3D int_x;
    //Calculate the integral on intersection portion
	_material_tet->intersect_with_polyhedron(_old_cell->G_cell_polyhedron(),volume,ix,iy,iz);
	//Density,energy and cell volume accumulation
	int_x.value(ix,iy,iz);
	_old_cell->accumulate_variables(mid,volume,int_x,0,0);
	return;
}
void Tcell_intersect::create_tet_list()
{
	vec3D v[8],c_face[6],c_edge[12],c_cell;   //points on node,surface,edge and center
	//Node point
	for (int i=0;i<8;i++)
	{
		v[i]=_new_cell->_node_ptr[i]->_position;
	}
	//center
	c_cell=(v[0]+v[1]+v[2]+v[3]+v[4]+v[5]+v[6]+v[7])/8.0;
	//Surface center
	c_face[0]=(v[0]+v[4]+v[7]+v[3])/4.0;c_face[1]=(v[0]+v[1]+v[5]+v[4])/4.0;
	c_face[2]=(v[1]+v[5]+v[6]+v[2])/4.0;c_face[3]=(v[2]+v[6]+v[7]+v[3])/4.0;
	c_face[4]=(v[0]+v[1]+v[2]+v[3])/4.0;c_face[5]=(v[4]+v[5]+v[6]+v[7])/4.0;
	//Edge center
	for (int i=0;i<4;i++)
	{
		c_edge[i]=(v[(i+3)%4]+v[(i+4)%4])/2.0;
		c_edge[i+4]=(v[(i+3)%4+4]+v[(i+4)%4+4])/2.0;
		c_edge[i+8]=(v[i]+v[i+4])/2.0;
	}
	//tet_list[6*i~6*i+5] are associated to p[i]
	for (int i=0;i<4;i++)
	{
		_tet_list[6*i+0]=Ttetrahedron(c_cell,v[i],c_face[i],c_edge[i]);
		_tet_list[6*i+1]=Ttetrahedron(c_cell,v[i],c_face[i],c_edge[i+8]);
		_tet_list[6*i+2]=Ttetrahedron(c_cell,v[i],c_face[(i+1)%4],c_edge[i+8]);
		_tet_list[6*i+3]=Ttetrahedron(c_cell,v[i],c_face[(i+1)%4],c_edge[(i+1)%4]);
		_tet_list[6*i+4]=Ttetrahedron(c_cell,v[i],c_face[4],c_edge[(i+1)%4]);
		_tet_list[6*i+5]=Ttetrahedron(c_cell,v[i],c_face[4],c_edge[i]);
	}
	for (int i=4;i<8;i++)
	{
		_tet_list[6*i+0]=Ttetrahedron(c_cell,v[i],c_face[i-4],c_edge[i]);
		_tet_list[6*i+1]=Ttetrahedron(c_cell,v[i],c_face[i-4],c_edge[i+4]);
		_tet_list[6*i+2]=Ttetrahedron(c_cell,v[i],c_face[(i-3)%4],c_edge[i+4]);
		_tet_list[6*i+3]=Ttetrahedron(c_cell,v[i],c_face[(i-3)%4],c_edge[(i+1)%4+4]);
		_tet_list[6*i+4]=Ttetrahedron(c_cell,v[i],c_face[5],c_edge[(i+1)%4+4]);
		_tet_list[6*i+5]=Ttetrahedron(c_cell,v[i],c_face[5],c_edge[i]);
	}
	//Set the associated node of every tetrahedrons
	for (int i=0;i<8;i++)
	{
		for (int j=0;j<6;j++)
		{
			_ass_node[6*i+j]=i;
		}
	}
	return ;
}
void Tcell_intersect::create_corner_polyhedron()
{
	vec3D v[8],c_face[6],c_edge[12],c_cell;   //points on node,surface,edge and center
	//Node point
	for (int i=0;i<8;i++)
	{
		v[i]=_old_cell->_node_ptr[i]->_position;
	}
	//center
	c_cell=(v[0]+v[1]+v[2]+v[3]+v[4]+v[5]+v[6]+v[7])/8.0;
	//Surface center
	c_face[0]=(v[0]+v[4]+v[7]+v[3])/4.0;c_face[1]=(v[0]+v[1]+v[5]+v[4])/4.0;
	c_face[2]=(v[1]+v[5]+v[6]+v[2])/4.0;c_face[3]=(v[2]+v[6]+v[7]+v[3])/4.0;
	c_face[4]=(v[0]+v[1]+v[2]+v[3])/4.0;c_face[5]=(v[4]+v[5]+v[6]+v[7])/4.0;
	//Edge center
	for (int i=0;i<4;i++)
	{
		c_edge[i]=(v[(i+3)%4]+v[(i+4)%4])/2.0;
		c_edge[i+4]=(v[(i+3)%4+4]+v[(i+4)%4+4])/2.0;
		c_edge[i+8]=(v[i]+v[i+4])/2.0;
	}
	//Generate the corner polyhedrons
	for (int i=0;i<4;i++)
	{
		//The 4 corner polyhedron below
		_corner_polyhedron[i].add_vertex(c_cell);_corner_polyhedron[i].add_vertex(v[i]);
		_corner_polyhedron[i].add_vertex(c_edge[i]);_corner_polyhedron[i].add_vertex(c_edge[(i+1)%4]);
		_corner_polyhedron[i].add_vertex(c_face[i]);_corner_polyhedron[i].add_vertex(c_face[(i+1)%4]);
		_corner_polyhedron[i].add_vertex(c_face[4]);_corner_polyhedron[i].add_vertex(c_edge[i+8]);
		_corner_polyhedron[i].add_piece(6,3,1);_corner_polyhedron[i].add_piece(6,1,2);
		_corner_polyhedron[i].add_piece(0,4,7);_corner_polyhedron[i].add_piece(0,7,5);
		_corner_polyhedron[i].add_piece(4,2,1);_corner_polyhedron[i].add_piece(4,1,7);
		_corner_polyhedron[i].add_piece(5,7,1);_corner_polyhedron[i].add_piece(5,1,3);
		_corner_polyhedron[i].add_piece(0,5,3);_corner_polyhedron[i].add_piece(0,3,6);
		_corner_polyhedron[i].add_piece(0,6,2);_corner_polyhedron[i].add_piece(0,2,4);
		//The 4 corner polyhedron above
		_corner_polyhedron[i+4].add_vertex(c_cell);_corner_polyhedron[i+4].add_vertex(v[i+4]);
		_corner_polyhedron[i+4].add_vertex(c_edge[i+4]);_corner_polyhedron[i+4].add_vertex(c_edge[(i+1)%4+4]);
		_corner_polyhedron[i+4].add_vertex(c_face[i]);_corner_polyhedron[i+4].add_vertex(c_face[(i+1)%4]);
		_corner_polyhedron[i+4].add_vertex(c_face[5]);_corner_polyhedron[i+4].add_vertex(c_edge[i+8]);
		_corner_polyhedron[i+4].add_piece(6,1,3);_corner_polyhedron[i+4].add_piece(6,2,1);
		_corner_polyhedron[i+4].add_piece(0,5,7);_corner_polyhedron[i+4].add_piece(0,7,4);
		_corner_polyhedron[i+4].add_piece(4,7,1);_corner_polyhedron[i+4].add_piece(4,1,2);
		_corner_polyhedron[i+4].add_piece(5,3,1);_corner_polyhedron[i+4].add_piece(5,1,7);
		_corner_polyhedron[i+4].add_piece(0,6,3);_corner_polyhedron[i+4].add_piece(0,3,5);
		_corner_polyhedron[i+4].add_piece(0,4,2);_corner_polyhedron[i+4].add_piece(0,2,6);
	}
	return;
}
bool Tcell_intersect::whether_dis(double tolerance)
{
	vec3D dux,duy,duz;
	double du;
	_old_cell->calculate_velocity_gradient(dux,duy,duz);
	du=sqrt(dux*dux+duy*duy+duz*duz);
	if (du>tolerance) 
	{
		return true;
	}
	else 
	{
		return false;
	}
}