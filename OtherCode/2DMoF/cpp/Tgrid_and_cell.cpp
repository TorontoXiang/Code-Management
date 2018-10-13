#include "Tgrid_and_cell.h"
#include "Ttriangle.h"
void Tgrid::plot_material_polygon(ofstream& output,int id)
{
	int num_vertex=0,num_polygon=0;
	for (int i = 0; i < nume; i++)
	{
		if (cell_list[i].fraction[id]==1)
		{
			num_vertex=num_vertex+5;num_polygon=num_polygon+4;
		}
		else if (cell_list[i].fraction[id]>0)
		{
			num_vertex=num_vertex+cell_list[i].material_polygon[id]->G_num_vertex()+1;
			num_polygon=num_polygon+cell_list[i].material_polygon[id]->G_num_vertex();
		}
	}
	output<<"TITLE=\"material polygon\""<<endl;
	output<<"Variables=\"X\",\"Y\",\"V\""<<endl;
	output<<"ZONE F=FEPOINT N= "<<num_vertex<<" E= "<<num_polygon<<" ET=TRIANGLE"<<endl;
	for (int i = 0; i < nume; i++)
	{
		if (cell_list[i].fraction[id]==1)
		{
			for (int j = 0; j < cell_list[i].polygon.G_num_vertex(); j++)
			{
				output<<cell_list[i].polygon.G_vertex(j).x<<" "<<cell_list[i].polygon.G_vertex(j).y<<" "<<id<<endl;
			}
			output<<cell_list[i].polygon.G_centroid().x<<" "<<cell_list[i].polygon.G_centroid().y<<" "<<id<<endl;
		}
		else if (cell_list[i].fraction[id]>0)
		{
			for (int j = 0; j < cell_list[i].material_polygon[id]->G_num_vertex(); j++)
			{
				output<<cell_list[i].material_polygon[id]->G_vertex(j).x<<" "<<cell_list[i].material_polygon[id]->G_vertex(j).y<<" "<<id<<endl;
			}
			cell_list[i].material_polygon[id]->calculate_properities();
			output<<cell_list[i].material_polygon[id]->G_centroid().x<<" "<<cell_list[i].material_polygon[id]->G_centroid().y<<" "<<id<<endl;
		}
	}
	int id_begin=1;
	for (int i = 0; i < nume; i++)
	{
		int vertex_this;
		if (cell_list[i].fraction[id]==1)
		{
			vertex_this=cell_list[i].polygon.G_num_vertex();
			for (int i = 0; i < vertex_this; i++)
			{
				output<<id_begin+vertex_this<<" "<<id_begin+i<<" "<<id_begin+(i+1)%vertex_this<<endl;
			}
			id_begin=id_begin+vertex_this+1;
		}
		else if (cell_list[i].fraction[id]>0)
		{
			vertex_this=cell_list[i].material_polygon[id]->G_num_vertex();
			for (int i = 0; i < vertex_this; i++)
			{
				output<<id_begin+vertex_this<<" "<<id_begin+i<<" "<<id_begin+(i+1)%vertex_this<<endl;
			}
			id_begin=id_begin+vertex_this+1;
		}
	}
	return;
}
void Tgrid::plot_grid(ofstream& output)
{
	int num_vertex=0,num_polygon=0;
	for (int i = 0; i < nume; i++)
	{
		num_vertex=num_vertex+5;num_polygon=num_polygon+4;
	}
	output<<"TITLE=\"material polygon\""<<endl;
	output<<"Variables=\"X\",\"Y\",\"V\""<<endl;
	output<<"ZONE F=FEPOINT N= "<<num_vertex<<" E= "<<num_polygon<<" ET=TRIANGLE"<<endl;
	for (int i = 0; i < nume; i++)
	{
		for (int j = 0; j < cell_list[i].polygon.G_num_vertex(); j++)
		{
			output<<cell_list[i].polygon.G_vertex(j).x<<" "<<cell_list[i].polygon.G_vertex(j).y<<" "<<0<<endl;
		}
		output<<cell_list[i].polygon.G_centroid().x<<" "<<cell_list[i].polygon.G_centroid().y<<" "<<0<<endl;
	}
	int id_begin=1;
	for (int i = 0; i < nume; i++)
	{
		int vertex_this;
		vertex_this=cell_list[i].polygon.G_num_vertex();
		for (int i = 0; i < vertex_this; i++)
		{
			output<<id_begin+vertex_this<<" "<<id_begin+i<<" "<<id_begin+(i+1)%vertex_this<<endl;
		}
		id_begin=id_begin+vertex_this+1;
	}
	return;
}
void Tgrid_polygon::generate_node_list()
{
	int nx=660,ny=200;
	_nx_node=nx;_ny_node=ny;
	double dx=3.2/nx,dy=1.0/ny;
	_node_list=new vec2D*[nx+1];
	for (int i = 0; i < nx+1; i++)
	{
		_node_list[i]=new vec2D[ny+1];
	}
	for (int i = 0; i < nx+1; i++)
	{
		for (int j = 0; j < ny+1; j++)
		{
			_node_list[i][j].value(i*dx,j*dy);
		}
	}
	return;
}
void Tgrid_polygon::generate_cell_list()
{
	int nx=220;
	int i=0,ny=51;
	for (int j = 0; j < ny; j++)
	{
		Tcell_polygon new_cell;
		new_cell._fraction[0]=new_cell._fraction[1]=0;
		if (j==0)
		{
			new_cell._node_list.push_back(&_node_list[0][0]);
			new_cell._node_list.push_back(&_node_list[2][0]);
			new_cell._node_list.push_back(&_node_list[1][2]);
			new_cell._node_list.push_back(&_node_list[0][2]);
		}
		else if (j==ny-1)
		{
			new_cell._node_list.push_back(&_node_list[0][_ny_node]);
			new_cell._node_list.push_back(&_node_list[0][_ny_node-2]);
			new_cell._node_list.push_back(&_node_list[1][_ny_node-2]);
			new_cell._node_list.push_back(&_node_list[2][_ny_node]);
		}
		else
		{
			new_cell._node_list.push_back(&_node_list[0][4*j-2]);
			new_cell._node_list.push_back(&_node_list[1][4*j-2]);
			new_cell._node_list.push_back(&_node_list[2][4*j]);
			new_cell._node_list.push_back(&_node_list[1][4*j+2]);
			new_cell._node_list.push_back(&_node_list[0][4*j+2]);
		}
		_cell_list.push_back(new_cell);
	}
	i=110;
	for (int j = 0; j < ny; j++)
	{
		Tcell_polygon new_cell;
		new_cell._fraction[0]=new_cell._fraction[1]=0;
		if (j==0)
		{
			new_cell._node_list.push_back(&_node_list[_nx_node-2][0]);
			new_cell._node_list.push_back(&_node_list[_nx_node][0]);
			new_cell._node_list.push_back(&_node_list[_nx_node][2]);
			new_cell._node_list.push_back(&_node_list[_nx_node-1][2]);
		}
		else if (j==ny-1)
		{
			new_cell._node_list.push_back(&_node_list[_nx_node][_ny_node]);
			new_cell._node_list.push_back(&_node_list[_nx_node-2][_ny_node]);
			new_cell._node_list.push_back(&_node_list[_nx_node-1][_ny_node-2]);
			new_cell._node_list.push_back(&_node_list[_nx_node][_ny_node-2]);
		}
		else
		{
			new_cell._node_list.push_back(&_node_list[_nx_node-2][4*j]);
			new_cell._node_list.push_back(&_node_list[_nx_node-1][4*j-2]);
			new_cell._node_list.push_back(&_node_list[_nx_node][4*j-2]);
			new_cell._node_list.push_back(&_node_list[_nx_node][4*j+2]);
			new_cell._node_list.push_back(&_node_list[_nx_node-1][4*j+2]);
		}
		_cell_list.push_back(new_cell);
	}
	for (int i = 1; i < 110; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			Tcell_polygon new_cell;
			new_cell._fraction[0]=new_cell._fraction[1]=0;
			if (j==0)
			{
				new_cell._node_list.push_back(&_node_list[6*i-2][0]);
				new_cell._node_list.push_back(&_node_list[6*i-1][2]);
				new_cell._node_list.push_back(&_node_list[6*i+1][2]);
				new_cell._node_list.push_back(&_node_list[6*i+2][0]);
			}
			else if (j==ny-1)
			{
				new_cell._node_list.push_back(&_node_list[6*i-2][_ny_node]);
				new_cell._node_list.push_back(&_node_list[6*i-1][_ny_node-2]);
				new_cell._node_list.push_back(&_node_list[6*i+1][_ny_node-2]);
				new_cell._node_list.push_back(&_node_list[6*i+2][_ny_node]);
			}
			else
			{
				int ic=6*i,jc=4*j;
				new_cell._node_list.push_back(&_node_list[ic-2][jc]);
				new_cell._node_list.push_back(&_node_list[ic-1][jc-2]);
				new_cell._node_list.push_back(&_node_list[ic+1][jc-2]);
				new_cell._node_list.push_back(&_node_list[ic+2][jc]);
				new_cell._node_list.push_back(&_node_list[ic+1][jc+2]);
				new_cell._node_list.push_back(&_node_list[ic-1][jc+2]);
			}
			_cell_list.push_back(new_cell);
		}
	}
	nx=219;ny=50;
	for (int i = 1; i < 111; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			Tcell_polygon new_cell;
			new_cell._fraction[0]=new_cell._fraction[1]=0;
			int ic=6*i-3,jc=2+4*(j);
			new_cell._node_list.push_back(&_node_list[ic-2][jc]);
			new_cell._node_list.push_back(&_node_list[ic-1][jc-2]);
			new_cell._node_list.push_back(&_node_list[ic+1][jc-2]);
			new_cell._node_list.push_back(&_node_list[ic+2][jc]);
			new_cell._node_list.push_back(&_node_list[ic+1][jc+2]);
			new_cell._node_list.push_back(&_node_list[ic-1][jc+2]);
			_cell_list.push_back(new_cell);
		}
	}
	for (int i = 0; i < _cell_list.size(); i++)
	{
		_cell_list[i].generate_cell_polygon();
	}
	return;
}
void Tcell_polygon::generate_cell_polygon()
{
	int edge_number=_node_list.size();
	_cell_polygon.S_num_edge(edge_number);
	for (int i = 0; i < edge_number; i++)
	{
		_cell_polygon.Add_vertex(*_node_list[i],i);
	}
	_cell_polygon.calculate_properities();
	return;
}
void Tgrid_polygon::plot_grid(ofstream& output,int id)
{
	int nume=_cell_list.size();
	int num_vertex=0,num_polygon=0;
	for (int i = 0; i < nume; i++)
	{
		num_vertex=num_vertex+_cell_list[i]._cell_polygon.G_num_vertex()+1;
		num_polygon=num_polygon+_cell_list[i]._cell_polygon.G_num_vertex();
	}
	output<<"TITLE=\"material polygon\""<<endl;
	output<<"Variables=\"X\",\"Y\",\"V\""<<endl;
	output<<"ZONE F=FEPOINT N= "<<num_vertex<<" E= "<<num_polygon<<" ET=TRIANGLE"<<endl;
	for (int i = 0; i < nume; i++)
	{
		for (int j = 0; j < _cell_list[i]._cell_polygon.G_num_vertex(); j++)
		{
			output<<_cell_list[i]._cell_polygon.G_vertex(j).x<<" "<<_cell_list[i]._cell_polygon.G_vertex(j).y<<" "<<id<<endl;
		}
			output<<_cell_list[i]._cell_polygon.G_centroid().x<<" "<<_cell_list[i]._cell_polygon.G_centroid().y<<" "<<id<<endl;
	}
	int id_begin=1;
	for (int i = 0; i < nume; i++)
	{
		int vertex_this;
		vertex_this=_cell_list[i]._cell_polygon.G_num_vertex();
		for (int i = 0; i < vertex_this; i++)
		{
			output<<id_begin+vertex_this<<" "<<id_begin+i<<" "<<id_begin+(i+1)%vertex_this<<endl;
		}
		id_begin=id_begin+vertex_this+1;
	}
	return;
}
void Tgrid_polygon::remapping_surface(Tgrid* grid_old)
{
	ofstream output;
	output.open("dam_breaking_MOF.dat");
	int nume=_cell_list.size();
	vec2D x_range_old,y_range_old,x_range_new,y_range_new;
	double error_max=0;
	for (int i = 0; i < nume; i++)
	{
		cout<<i<<" ";
		_cell_list[i]._cell_polygon.calculate_range(x_range_new,y_range_new);
		for (int j = 0; j < grid_old->nume; j++)
		{
			grid_old->cell_list[j].polygon.calculate_range(x_range_old,y_range_old);
			if (is_intersect(x_range_old,y_range_old,x_range_new,y_range_new))
			{
				_cell_list[i].intersecting_with_other(&grid_old->cell_list[j]);
			}
		}
		vec2D centroid=(_cell_list[i]._centroid[0]+_cell_list[i]._centroid[1])/_cell_list[i]._cell_polygon.G_area();
		double centroid_error=sqrt((centroid-_cell_list[i]._cell_polygon.G_centroid())*(centroid-_cell_list[i]._cell_polygon.G_centroid()));
		//cout<<centroid_error<<endl;
		for (int j = 0; j < 2; j++)
		{
			if (_cell_list[i]._fraction[j]>0)
			{
				_cell_list[i]._centroid[j]=_cell_list[i]._centroid[j]/_cell_list[i]._fraction[j];
				_cell_list[i]._fraction[j]=_cell_list[i]._fraction[j]/_cell_list[i]._cell_polygon.G_area();
			}
		}
		for (int j = 0; j < 2; j++)
		{
			if (abs(_cell_list[i]._fraction[j])<1e-6)
			{
				_cell_list[i]._fraction[j]=0;
			}
			if (abs(_cell_list[i]._fraction[j]-1)<1e-6)
			{
				_cell_list[i]._fraction[j]=1;
			}
		}
		error_max=maxval(error_max,centroid_error);
	}
	for (int i = 0; i < nume; i++)
	{
		output<<_cell_list[i]._fraction[0]<<" "<<_cell_list[i]._centroid[0].x<<" "<<_cell_list[i]._centroid[0].y<<endl;
		output<<_cell_list[i]._fraction[1]<<" "<<_cell_list[i]._centroid[1].x<<" "<<_cell_list[i]._centroid[1].y<<endl;
	}
	cout<<_cell_list[3435]._fraction[0]<<" "<<_cell_list[3435]._fraction[1]<<endl;
	cout<<error_max<<endl;
	return;
}
void Tcell_polygon::intersecting_with_other(Tcell* other)
{
	int num_edge=_cell_polygon.G_num_vertex();
	double area;
	vec2D moment;
	for (int i = 0; i < num_edge; i++)
	{
		triangle triangle_piece(_cell_polygon.G_centroid(),_cell_polygon.G_vertex(i),_cell_polygon.G_vertex((i+1)%num_edge));
		for (int j = 0; j < 2; j++)
		{
			if (other->material_polygon[j]!=NULL)
			{
				triangle_piece.intersecting_with_polygon(other->material_polygon[j],area,moment);
				_fraction[j]=_fraction[j]+area;
				_centroid[j]=_centroid[j]+moment;
			}
		}
	}
	return;
}
void Tgrid_polygon::plot_material(ofstream& output,int id)
{
	int nume=_cell_list.size();
	int num_vertex=0,num_polygon=0;
	for (int i = 0; i < nume; i++)
	{
		//cout<<i<<endl;
		if (_cell_list[i]._fraction[id]==1)
		{
			int vertex_cell=_cell_list[i]._cell_polygon.G_num_vertex();
			num_vertex=num_vertex+vertex_cell+1;num_polygon=num_polygon+vertex_cell;
		}
		else if (_cell_list[i]._fraction[id]>0)
		{
			num_vertex=num_vertex+_cell_list[i]._material_polygon[id]->G_num_vertex()+1;
			num_polygon=num_polygon+_cell_list[i]._material_polygon[id]->G_num_vertex();
		}
	}
	output<<"TITLE=\"material polygon\""<<endl;
	output<<"Variables=\"X\",\"Y\",\"V\""<<endl;
	output<<"ZONE F=FEPOINT N= "<<num_vertex<<" E= "<<num_polygon<<" ET=TRIANGLE"<<endl;
	for (int i = 0; i < nume; i++)
	{
		if (_cell_list[i]._fraction[id]==1)
		{
			for (int j = 0; j < _cell_list[i]._cell_polygon.G_num_vertex(); j++)
			{
				output<<_cell_list[i]._cell_polygon.G_vertex(j).x<<" "<<_cell_list[i]._cell_polygon.G_vertex(j).y<<" "<<id<<endl;
			}
			output<<_cell_list[i]._cell_polygon.G_centroid().x<<" "<<_cell_list[i]._cell_polygon.G_centroid().y<<" "<<id<<endl;
		}
		else if (_cell_list[i]._fraction[id]>0)
		{
			for (int j = 0; j < _cell_list[i]._material_polygon[id]->G_num_vertex(); j++)
			{
				output<<_cell_list[i]._material_polygon[id]->G_vertex(j).x<<" "<<_cell_list[i]._material_polygon[id]->G_vertex(j).y<<" "<<id<<endl;
			}
			_cell_list[i]._material_polygon[id]->calculate_properities();
			output<<_cell_list[i]._material_polygon[id]->G_centroid().x<<" "<<_cell_list[i]._material_polygon[id]->G_centroid().y<<" "<<id<<endl;
		}
	}
	int id_begin=1;
	for (int i = 0; i < nume; i++)
	{
		int vertex_this;
		if (_cell_list[i]._fraction[id]==1)
		{
			vertex_this=_cell_list[i]._cell_polygon.G_num_vertex();
			for (int i = 0; i < vertex_this; i++)
			{
				output<<id_begin+vertex_this<<" "<<id_begin+i<<" "<<id_begin+(i+1)%vertex_this<<endl;
			}
			id_begin=id_begin+vertex_this+1;
		}
		else if (_cell_list[i]._fraction[id]>0)
		{
			vertex_this=_cell_list[i]._material_polygon[id]->G_num_vertex();
			for (int i = 0; i < vertex_this; i++)
			{
				output<<id_begin+vertex_this<<" "<<id_begin+i<<" "<<id_begin+(i+1)%vertex_this<<endl;
			}
			id_begin=id_begin+vertex_this+1;
		}
	}
	return;
}