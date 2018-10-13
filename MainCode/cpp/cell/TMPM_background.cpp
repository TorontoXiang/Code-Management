#include "TMPM_background.h"
void TMPM_background::generate_background(vec3D& x_min,vec3D& x_max,int nx,int ny,int nz)
{
	_x_min=x_min;_x_max=x_max;
	_nx=nx;_ny=ny;_nz=nz;
	//Allocate the node and cell information in sutrcuture mesh
	vec3D*** coordinate_list;
	coordinate_list=new vec3D**[_nx+1];
	for (int i = 0; i < _nx+1; i++)
	{
		coordinate_list[i]=new vec3D*[_ny+1];
	}
	for (int i = 0; i < _nx+1; i++)
	{
		for (int j = 0; j < _ny+1; j++)
		{
			coordinate_list[i][j]=new vec3D[_nz+1];
		}
	}
	Tcell_MPM*** cell_list;
	cell_list=new Tcell_MPM**[_nx];
	for (int i = 0; i < _nx; i++)
	{
		cell_list[i]=new Tcell_MPM*[_ny];
	}
	for (int i = 0; i < _nx; i++)
	{
		for (int j = 0; j < _ny; j++)
		{
			cell_list[i][j]=new Tcell_MPM[_nz];
		}
	}
	_interval.x=(_x_max-_x_min).x/_nx;
	_interval.y=(_x_max-_x_min).x/_ny;
	_interval.z=(_x_max-_x_min).z/_nz;
	Tcell_MPM::_length=_interval;
	vec3D dx;
	for (int i = 0; i < _nx+1; i++)
	{
		for (int j = 0; j < _ny+1; j++)
		{
			for (int k = 0; k < _nz+1; k++)
			{
				dx.value(i*_interval.x,j*_interval.y,k*_interval.z);
				coordinate_list[i][j][k]=_x_min+dx;
			}
		}
	}
	_total_node=(_nx+1)*(_ny+1)*(_nz+1);
	_total_cell=_nx*_ny*_nz;
	_node_list.resize(_total_node);
	_cell_list.resize(_total_cell);
	//Convert the structure mesh into the unstructure mesh
	int index;
	for (int i = 0; i < _nx+1; i++)
	{
		for (int j = 0; j < _ny+1; j++)
		{
			for (int k = 0; k < _nz+1; k++)
			{
				index=(_ny+1)*(_nz+1)*i+(_nz+1)*j+k;
				_node_list[index]._position0=coordinate_list[i][j][k];
				_node_list[index]._id=index;
			}
		}
	}
	int index_node[8];
	Tnode* node_ptr[8];
	for (int i = 0; i < _nx; i++)
	{
		for (int j = 0; j < _ny; j++)
		{
			for (int k = 0; k < _nz; k++)
			{
				index=_ny*_nz*i+_nz*j+k; //The index convertion for cell
				_cell_list[index]=cell_list[i][j][k];
				index_node[0]=(_ny+1)*(_nz+1)*i+(_nz+1)*j+k;index_node[1]=(_ny+1)*(_nz+1)*(i+1)+(_nz+1)*j+k;
				index_node[2]=(_ny+1)*(_nz+1)*(i+1)+(_nz+1)*(j+1)+k;index_node[3]=(_ny+1)*(_nz+1)*i+(_nz+1)*(j+1)+k;
				index_node[4]=(_ny+1)*(_nz+1)*i+(_nz+1)*j+k+1;index_node[5]=(_ny+1)*(_nz+1)*(i+1)+(_nz+1)*j+k+1;
				index_node[6]=(_ny+1)*(_nz+1)*(i+1)+(_nz+1)*(j+1)+k+1;index_node[7]=(_ny+1)*(_nz+1)*i+(_nz+1)*(j+1)+k+1;
				for (int n = 0; n < 8; n++)
				{
					node_ptr[n]=&_node_list[index_node[n]];
				}
				_cell_list[index].calculate_node_ptr(node_ptr);
				_cell_list[index]._cell_id=index;
			}
		}
	}

	//Delete the structure mesh
	for (int i = 0; i < _nx; i++)
	{
		for (int j = 0; j < _ny; j++)
		{
			delete cell_list[i][j];
		}
	}
	for (int i = 0; i < _nx; i++)
	{
		delete cell_list[i];
	}
	delete cell_list;
	for (int i = 0; i < _nx+1; i++)
	{
		for (int j = 0; j < _ny+1; j++)
		{
			delete coordinate_list[i][j];
		}
	}
	for (int i = 0; i < _nx+1; i++)
	{
		delete coordinate_list[i];
	}
	delete coordinate_list;
	return;
}
void TMPM_background::initialize()
{
	for (int i = 0; i < _total_node; i++)
	{
		_node_list[i].initialize();
	}
	for (int i = 0; i < _total_cell; i++)
	{
		_cell_list[i].initialize();
	}
	return;
}
void TMPM_background::distribute_particle_force(TMPM_particle* particle_ptr)
{
	int cell_id=particle_ptr->_cell_index;
	_cell_list[cell_id].distribute_particle_force(particle_ptr);
	return;
}
void TMPM_background::map_particle_to_cell(TMPM_particle* particle_ptr)
{
	int cell_id=calculate_paritcle_loaction(particle_ptr);
	particle_ptr->_cell_index=cell_id;
	_cell_list[cell_id].distribute_particle_to_node(particle_ptr);
	return;
}
void TMPM_background::calculate_character_length()
{
	for (int i = 0; i < _total_cell; i++)
	{
		_cell_list[i].calculate_cell_length();
	}
	return;
}
int TMPM_background::calculate_paritcle_loaction(TMPM_particle* particle_ptr)
{
	vec3D particle_coordinate=particle_ptr->_position;
	for (int i = 1; i <= 3; i++)
	{
		if (particle_coordinate.access(i)<_x_min.access(i) || particle_coordinate.access(i)>_x_max.access(i))
		{
			cout<<"The particle is out of the background ground"<<endl;
			cout<<"Particleis located in "<<particle_ptr->_part_id<<"th part "<<particle_ptr->_id_in_part<<"th particle"<<endl;
			cout<<"The coordinate is: "<<particle_coordinate.x<<" "<<particle_coordinate.y<<" "<<particle_coordinate.z<<endl;
			system("Pause");
			exit(0);
		}
	}
	int ix,iy,iz;
	ix=int((particle_coordinate-_x_min).x/_interval.x);
	iy=int((particle_coordinate-_x_min).y/_interval.y);
	iz=int((particle_coordinate-_x_min).z/_interval.z);
	return (_ny*_nz)*ix+_nz*iy+iz;
}
void TMPM_background::update_background_DOF(double dt)
{
	for (int i = 0; i < _total_node; i++)
	{
		//Calculate the velocity at the next time step
		_node_list[i].update_vel(dt);
		//Calculate the position at the next time point
		_node_list[i].update_pos(dt);
	}
	return;
}
void TMPM_background::update_particle_DOF(TMPM_particle* particle_ptr,double dt)
{
	int cell_id=particle_ptr->_cell_index;
	_cell_list[cell_id].update_particle_DOF(particle_ptr,dt);
}
void TMPM_background::assemble_nodal_force()
{
	for (int i = 0; i < _total_cell; i++)
	{
		_cell_list[i].assemble_corner_force();
	}
	return;
}
void TMPM_background::update_particle_state(TMPM_particle* particle_ptr,double dt)
{
	int cell_id=particle_ptr->_cell_index;
	_cell_list[cell_id].update_particle_state(particle_ptr,dt);
	return;
}
void TMPM_background::reset_nodal_velocity()
{
	for (int i = 0; i < _total_node; i++)
	{
		if (_node_list[i]._mass>0)
		{
			_node_list[i]._velocity=_node_list[i]._moment/_node_list[i]._mass;
			for (int j = 0; j < 3; j++)
			{
				if (_node_list[i]._bc_type_position[j]==1)
				{
					_node_list[i]._velocity.value(j+1,0);
				}
			}
			
		}
	}
	return;
}
void TMPM_background::apply_background_load()
{
	for (int i = 0; i < _total_node; i++)
	{
		if (_node_list[i]._mass>0)
		{
			_node_list[i]._force=_node_list[i]._force+_node_list[i]._external_load;	
		}
	}
	return;
}
void TMPM_background::plot_background(ofstream& output)
{
	output<<"TITLE = \"Background grid\""<<endl;
	output<<"VARIABLES = \"X\",\"Y\",\"Z\""<<endl;
	output<<"ZONE F=FEPOINT,N="<<_total_node<<","<<"E="<<_total_cell<<","<<"ET=BRICK"<<endl;
	for (int i = 0; i < _total_node; i++)
	{
		vec3D coordiante=_node_list[i]._position;
		output<<coordiante.x<<" "<<coordiante.y<<" "<<coordiante.z<<endl;
	}
	for (int i = 0; i < _total_cell; i++)
	{
		Tcell_MPM* cell_ptr;
		cell_ptr=&_cell_list[i];
		for (int j = 0; j < 8; j++)
		{
			output<<cell_ptr->_node_ptr[j]->_id+1<<" ";
		}
		output<<endl;
	}
	return;
}
void TMPM_background::plot_background_cell(ofstream& output,int cell_id)
{
	output<<"TITLE = \"Background grid\""<<endl;
	output<<"VARIABLES = \"X\",\"Y\",\"Z\""<<endl;
	output<<"ZONE F=FEPOINT,N="<<8<<","<<"E="<<1<<","<<"ET=BRICK"<<endl;
	Tcell_MPM* cell_ptr;
	cell_ptr=&_cell_list[cell_id];
	for (int j = 0; j < 8; j++)
	{
		vec3D coordinate=cell_ptr->_node_ptr[j]->_position;;
		output<<coordinate.x<<" "<<coordinate.y<<" "<<coordinate.z<<endl;
	}
	for (int j = 0; j < 8; j++)
	{
		output<<j+1<<" ";
	}
	output<<endl;
	return;
}
