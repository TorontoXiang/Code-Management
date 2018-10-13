#include "Tgrid_smooth.h"
Tgrid_smooth_OST::Tgrid_smooth_OST(Tbody_ALE* the_body):Tgrid_smooth_base(the_body)
{
	_smooth_name="orthogonality similarity transformation";
	_x_min0=_y_min0=_z_min0=1e10;
	_x_max0=_y_max0=_z_max0=-1e10;
	vec3D coor;
	for (int i=0;i<_ALE_body->_nump;i++)
	{
		coor=_ALE_body->_grid1._node_list[i]->_position0;
		if (coor.x<_x_min0) _x_min0=coor.x;
		if (coor.y<_y_min0) _y_min0=coor.y;
		if (coor.z<_z_min0) _z_min0=coor.y;
		if (coor.x>_x_max0) _x_max0=coor.x;
		if (coor.y>_y_max0) _y_max0=coor.y;
		if (coor.z>_z_max0) _z_max0=coor.y;
	}
	return;
}
void Tgrid_smooth_OST::generate_new_grid()
{
	vec3D coor,coor_new;
	double x_min,y_min,z_min,x_max,y_max,z_max;
	x_min=y_min=z_min=1e10;
	x_max=y_max=z_max=-1e10;
	for (int i=0;i<_ALE_body->_nump;i++)
	{
		coor=_ALE_body->_grid_old->_node_list[i]->_position;
		if (coor.x<x_min) x_min=coor.x;
		if (coor.y<y_min) y_min=coor.y;
		if (coor.z<z_min) z_min=coor.y;
		if (coor.x>x_max) x_max=coor.x;
		if (coor.y>y_max) y_max=coor.y;
		if (coor.z>z_max) z_max=coor.y;
	}
	for (int i = 0; i < _ALE_body->_nump; i++)
	{
		coor=_ALE_body->_grid_old->_node_list[i]->_position0;
		coor_new.x=x_min+(coor.x-_x_min0)*(x_max-x_min)/(_x_max0-_x_min0);
		coor_new.y=y_min+(coor.y-_y_min0)*(y_max-y_min)/(_y_max0-_y_min0);
		coor_new.z=z_min+(coor.z-_z_min0)*(z_max-z_min)/(_z_max0-_z_min0);
		_ALE_body->_grid_new->_node_list[i]->_position=coor_new;
	}
}
Tgrid_smooth_CST::Tgrid_smooth_CST(Tbody_ALE* the_body,vec2D center):Tgrid_smooth_base(the_body)
{
	_center=center;
	_smooth_name="cylinder similarity transformation";
	vec3D coor;
	double avr=0;
	int num=0;
	for (int i=0;i<_ALE_body->_nump;i++)
	{
		if (_ALE_body->_grid1._node_list[i]->_bc_type_position[0]==2)
		{
			coor=_ALE_body->_grid1._node_list[i]->_position0;
			vec2D r(coor.x-center.x,coor.y-center.y);
			avr=avr+r.get_length();
			num=num+1;
		}
	}
	avr=avr/num;
	_r0=avr;
	for (int i = 0; i < _ALE_body->_nump; i++)
	{
		if (_ALE_body->_grid1._node_list[i]->_bc_type_position[0]==2)
		{
			coor=_ALE_body->_grid1._node_list[i]->_position0;
			vec2D r(coor.x-center.x,coor.y-center.y);
			vec2D r_new=r*avr/r.get_length();
			vec3D coor_new(r_new.x,r_new.y,coor.z);
			_ALE_body->_grid1._node_list[i]->_position0=coor_new;
			_ALE_body->_grid2._node_list[i]->_position0=coor_new;
		}
	}
	return;
}
void Tgrid_smooth_CST::generate_new_grid()
{
	vec3D coor,coor_new;
	double avr=0;
	int num=0;
	for (int i=0;i<_ALE_body->_nump;i++)
	{
		if (_ALE_body->_grid_old->_node_list[i]->_bc_type_position[0]==2)
		{
			coor=_ALE_body->_grid_old->_node_list[i]->_position;
			vec2D r(coor.x-_center.x,coor.y-_center.y);
			avr=avr+r.get_length();
			num=num+1;
		}
	}
	avr=avr/num;
	for (int i=0;i<_ALE_body->_nump;i++)
	{
		Tnode_fluid* node_ptr=_ALE_body->_grid_old->_node_list[i];
		coor=node_ptr->_position0;
		if (node_ptr->_bc_type_position[0]==0)
		{
			coor_new.x=(coor.x-_center.x)*avr/_r0+_center.x;
			//cout<<coor_new.x<<" "<<coor.x<<endl;
		}
		else if (node_ptr->_bc_type_position[0]==1)
		{
			coor_new.x=coor.x;
		}
		if (node_ptr->_bc_type_position[1]==0)
		{
			coor_new.y=(coor.y-_center.y)*avr/_r0+_center.y;
		}
		else if (node_ptr->_bc_type_position[1]==1)
		{
			coor_new.y=coor.y;
		}
		if (node_ptr->_bc_type_position[0]==2)
		{
			coor_new=node_ptr->_position;
			//cout<<sqrt(coor.x*coor.x+coor.y*coor.y)<<endl;
		}
		coor_new.z=coor.z;
		_ALE_body->_grid_new->_node_list[i]->_position=coor_new;
	}

	return;
}
Tgrid_smooth_PST::Tgrid_smooth_PST(Tbody_ALE* the_body,vec3D center):Tgrid_smooth_base(the_body)
{
	_center=center;
	_smooth_name="cylinder similarity transformation";
	vec3D coor;
	double avr=0;
	int num=0;
	for (int i=0;i<_ALE_body->_nump;i++)
	{
		if (_ALE_body->_grid1._node_list[i]->_bc_type_position[0]==2)
		{
			coor=_ALE_body->_grid1._node_list[i]->_position0;
			vec3D r=coor-_center;
			avr=avr+r.get_length();
			num=num+1;
		}
	}
	avr=avr/num;
	_r0=avr;
	for (int i = 0; i < _ALE_body->_nump; i++)
	{
		if (_ALE_body->_grid1._node_list[i]->_bc_type_position[0]==2)
		{
			coor=_ALE_body->_grid1._node_list[i]->_position0;
			vec3D r=coor-_center;
			vec3D r_new=r*avr/r.get_length();
			vec3D coor_new(r_new.x,r_new.y,r_new.z);
			_ALE_body->_grid1._node_list[i]->_position0=coor_new;
			_ALE_body->_grid2._node_list[i]->_position0=coor_new;
		}
	}
	return;
}
void Tgrid_smooth_PST::generate_new_grid()
{
	vec3D coor,coor_new;
	double avr=0;
	int num=0;
	for (int i=0;i<_ALE_body->_nump;i++)
	{
		if (_ALE_body->_grid_old->_node_list[i]->_bc_type_position[0]==2)
		{
			coor=_ALE_body->_grid_old->_node_list[i]->_position;
			vec3D r=coor-_center;
			avr=avr+r.get_length();
			num=num+1;
		}
	}
	avr=avr/num;
	for (int i=0;i<_ALE_body->_nump;i++)
	{
		Tnode_fluid* node_ptr=_ALE_body->_grid_old->_node_list[i];
		coor=node_ptr->_position0;
		if (node_ptr->_bc_type_position[0]==0)
		{
			coor_new.x=(coor.x-_center.x)*avr/_r0+_center.x;
			//cout<<coor_new.x<<" "<<coor.x<<endl;
		}
		else if (node_ptr->_bc_type_position[0]==1)
		{
			coor_new.x=coor.x;
		}
		if (node_ptr->_bc_type_position[1]==0)
		{
			coor_new.y=(coor.y-_center.y)*avr/_r0+_center.y;
		}
		else if (node_ptr->_bc_type_position[1]==1)
		{
			coor_new.y=coor.y;
		}
		if (node_ptr->_bc_type_position[2]==0)
		{
			coor_new.z=(coor.z-_center.z)*avr/_r0+_center.z;
		}
		else if (node_ptr->_bc_type_position[2]==1)
		{
			coor_new.z=coor.z;
		}
		if (node_ptr->_bc_type_position[0]==2)
		{
			coor_new=node_ptr->_position;
			//cout<<sqrt(coor.x*coor.x+coor.y*coor.y)<<endl;
		}
		_ALE_body->_grid_new->_node_list[i]->_position=coor_new;
	}

	return;
}