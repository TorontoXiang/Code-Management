#include "Tcell_pure_ALE.h"
#include "data_structure.h"
#include "public_function.h"
#include "Shape_function.h"
Tcell_pure_ALE::Tcell_pure_ALE(int cell_id,int nGauss,Tnode_fluid* (&node_ptr)[8],TMAT_base** mat)
	           :Tcell_pure_Lagrangian(cell_id,nGauss,node_ptr,mat)
{
	create_cell_polyhedron();
	_g_ro.resize(1);_g_roe.resize(1);
	_cell_centroid=_cell_polyhedron.G_centroid();
	_volume_intersection=0;
}
void Tcell_pure_ALE::create_cell_polyhedron()
{
	//Create the polyhedron vertexes
	vec3D temp;
	for (int i=0;i<8;i++)
	{
		vec3D_double new_vertex(_node_ptr[i]->_position);
		_cell_polyhedron.add_vertex(_node_ptr[i]->_position);
	}
	for (int i=0;i<4;i++)
	{
		temp=(_node_ptr[i]->_position+_node_ptr[(i+1)%4]->_position+_node_ptr[i+4]->_position+_node_ptr[(i+1)%4+4]->_position)*0.25;
		vec3D_double new_vertex(temp);
		_cell_polyhedron.add_vertex(temp);
	}
	for (int i=0;i<2;i++)
	{
		temp=(_node_ptr[4*i]->_position+_node_ptr[4*i+1]->_position+_node_ptr[4*i+2]->_position+_node_ptr[4*i+3]->_position)*0.25;
		vec3D_double new_vertex(temp);
		_cell_polyhedron.add_vertex(temp);
	}
	//Create the triangle pieces
	for (int i=0;i<4;i++)
	{
		_cell_polyhedron.add_piece(i+8,i,(i+1)%4);
		_cell_polyhedron.add_piece(i+8,(i+1)%4,(i+5)%4+4);
		_cell_polyhedron.add_piece(i+8,(i+5)%4+4,i+4);
		_cell_polyhedron.add_piece(i+8,i+4,i);
	}
	for (int i=0;i<4;i++)
	{
		_cell_polyhedron.add_piece(12,(4-i)%4,(7-i)%4);
		_cell_polyhedron.add_piece(13,i+4,(i+5)%4+4);
	}
	//Calculate the volume and centroid of polyhedron after generating it
	_cell_polyhedron.calculate_polyhedron_centroid();
	return;
}
void Tcell_pure_ALE::update_cell_polyhedron()
{
	vec3D temp;
	for (int i=0;i<8;i++)
	{
		_cell_polyhedron.S_vertex(i,_node_ptr[i]->_position);
	}
	for (int i=0;i<4;i++)
	{
		temp=(_node_ptr[i]->_position+_node_ptr[(i+1)%4]->_position+_node_ptr[i+4]->_position+_node_ptr[(i+1)%4+4]->_position)*0.25;
		_cell_polyhedron.S_vertex(i+8,temp);
	}
	for (int i=0;i<2;i++)
	{
		temp=(_node_ptr[4*i]->_position+_node_ptr[4*i+1]->_position+_node_ptr[4*i+2]->_position+_node_ptr[4*i+3]->_position)*0.25;
		_cell_polyhedron.S_vertex(i+12,temp);
	}
	//Calculate the volume and centroid of polyhedron after updating it
	_cell_polyhedron.calculate_polyhedron_centroid();
	return;
}
void Tcell_pure_ALE::clear_cell_variable()
{
	_cell_mass=_cell_volume=_cell_length=_volume_intersection=_cell_energy=0;                
	for (int i = 0; i < 8; i++)
	{
		_corner_mass[i]=_corner_volume[i]=0;
		_corner_hourglass_force[i].value(0,0,0); 
		_corner_viscosity_force[i].value(0,0,0);
		_corner_moment[i].value(0,0,0);
	}      
	_gausspoint[0]->S_density(0);
	_gausspoint[0]->S_internal_energy(0);
	_gausspoint[0]->S_pressure(0);
	_gausspoint[0]->S_soundspeed(0);
	return;     
}
void Tcell_pure_ALE::calculate_gradient_in_cell()
{
	double k[3][3];
	vec3D b,c;
	bool degenerate;
	double det;
	for (int j=0;j<3;j++)
	{
		for (int l=0;l<3;l++)
		{
			k[j][l]=0;
		}
	}
	for (unsigned int j=0;j<_around.size();j++)
	{
		double ro,roe,ro0,roe0;
		ro=_around[j]->G_gausspoint(0)->G_density();ro0=_gausspoint[0]->G_density();
		roe=_around[j]->G_gausspoint(0)->G_roe();roe0=_gausspoint[0]->G_roe();
		vec3D centroid_around=_around[j]->G_cell_centroid();
		b=b+(centroid_around-_cell_centroid)*(ro-ro0);
		c=c+(centroid_around-_cell_centroid)*(roe-roe0);
		for (int s = 0; s < 3; s++)
		{
			for (int t = s; t < 3; t++)
			{
				k[s][t]=k[s][t]+(centroid_around-_cell_centroid).access(s+1)*(centroid_around-_cell_centroid).access(t+1);
			}
		}
	}
	k[1][0]=k[0][1];k[2][0]=k[0][2];k[2][1]=k[1][2];
	inverse(k,det,degenerate);
	if (degenerate)
	{
		//If k is degenerated, set the gradient to be zero
		_g_ro[0].x=_g_ro[0].y=_g_ro[0].z=0;
		_g_roe[0].x=_g_roe[0].y=_g_roe[0].z=0;
	}
	else
	{
		_g_ro[0]=b.multiply_by_matrix(k);_g_roe[0]=c.multiply_by_matrix(k);
	}
	return;
}
void Tcell_pure_ALE::limit_gradient_in_cell()
{
	double ro_min,ro_max,roe_min,roe_max;           //The minimum and maximum density and energy in surrounding cells
	double ro_cell,roe_cell;
	double ro_around,roe_around;                                                                
	double ro_vertex,roe_vertex;                     
	double limiter_ro,limiter_roe;                  //limiter factor at each vertexes
	double limiter_ro_min,limiter_roe_min;          //The minimum limiter factor which is used in gradient limitation
	vec3D vertex;
	ro_max=ro_min=_gausspoint[0]->G_density();roe_max=roe_min=_gausspoint[0]->G_roe();

	//Calculate the minimum and maximum density and energy in surrounding cells
	for (unsigned int j=0;j<_around.size();j++)
	{
		ro_around=_around[j]->G_gausspoint(0)->G_density();
		roe_around=_around[j]->G_gausspoint(0)->G_roe();
		if (ro_around>ro_max) ro_max=ro_around;
		if (ro_around<ro_min) ro_min=ro_around;
		if (roe_around>roe_max) roe_max=roe_around;
		if (roe_around<roe_min) roe_min=roe_around;
	}
	//Calculate the minimum limiter factor
	limiter_ro_min=limiter_roe_min=1e50;		
	ro_cell=_gausspoint[0]->G_density();
	roe_cell=_gausspoint[0]->G_roe();
	for (int i=0;i<_cell_polyhedron.G_numv();i++)
	{
		vertex=_cell_polyhedron.G_vertex(i);

		//Calculate the limiter for ro
		ro_vertex=ro_cell+_g_ro[0]*(vertex-_cell_centroid);
		if (ro_vertex>ro_cell)
		{
			limiter_ro=minval(1.0,(ro_max-ro_cell)/(ro_vertex-ro_cell));
		}
		else if (ro_vertex<ro_cell)
		{
			limiter_ro=minval(1.0,(ro_min-ro_cell)/(ro_vertex-ro_cell));
		}
		else
		{
			limiter_ro=1;
		}
		if (limiter_ro<limiter_ro_min) limiter_ro_min=limiter_ro;

		//Calculate the limiter for roe
		roe_vertex=roe_cell+_g_roe[0]*(vertex-_cell_centroid);
		if (roe_vertex>roe_cell)
		{
			limiter_roe=minval(1.0,(roe_max-roe_cell)/(roe_vertex-roe_cell));
		}
		else if (roe_vertex<roe_cell)
		{
			limiter_roe=minval(1.0,(roe_min-roe_cell)/(roe_vertex-roe_cell));
		}
		else
		{
			limiter_roe=1;
		}
		if (limiter_roe<limiter_roe_min) limiter_roe_min=limiter_roe;
	}
	_g_ro[0]=_g_ro[0]*limiter_ro_min;
	_g_roe[0]=_g_roe[0]*limiter_roe_min;
	return;
}
void Tcell_pure_ALE::calculate_coordinate_range()
{
	_coor_max.value(-1e10,-1e10,-1e10);
	_coor_min.value(1e10,1e10,1e10);
	for (int i=0;i<8;i++)
	{
		for (int j = 1; j < 4; j++)
		{
			if (_node_ptr[i]->_position.access(j)>_coor_max.access(j)) _coor_max.value(j,_node_ptr[i]->_position.access(j));
			if (_node_ptr[i]->_position.access(j)<_coor_min.access(j)) _coor_min.value(j,_node_ptr[i]->_position.access(j));
		}
	}
	return;
}
double Tcell_pure_ALE::calculate_maximal_edge()
{
	double edge_max=0;
	double edge[4];
	edge[0]=(_node_ptr[0]->_position-_node_ptr[6]->_position).self_multuply();
	edge[1]=(_node_ptr[1]->_position-_node_ptr[7]->_position).self_multuply();
	edge[2]=(_node_ptr[2]->_position-_node_ptr[4]->_position).self_multuply();
	edge[3]=(_node_ptr[3]->_position-_node_ptr[5]->_position).self_multuply();
	for (int i=0;i<4;i++)
	{
		if (edge[i]>edge_max)
		{
			edge_max=edge[i];
		}
	}
	return sqrt(edge_max);
}
vec3D Tcell_pure_ALE::calculate_std_coordinate(vec3D coordinate_physical)
{
	//Calculate the initial intertion point
	double v1,v2,v3,v4,v5,v6,v;
	vec3D p[14];
	vec3D result;
	for (int i=0;i<14;i++)
	{
		p[i]=_cell_polyhedron.G_vertex(i);
	}
	v4=(coordinate_physical-p[10])*((p[6]-p[3])%(p[7]-p[2]))/6.0;
	v3=(coordinate_physical-p[8])*((p[5]-p[0])%(p[1]-p[4]))/6.0;
	v1=(coordinate_physical-p[11])*((p[4]-p[3])%(p[0]-p[7]))/6.0;
	v2=(coordinate_physical-p[9])*((p[6]-p[1])%(p[2]-p[5]))/6.0;
	v5=(coordinate_physical-p[12])*((p[1]-p[3])%(p[2]-p[0]))/6.0;
	v6=(coordinate_physical-p[13])*((p[5]-p[7])%(p[4]-p[6]))/6.0;
	v=(v1+v2)*(v3+v4)*(v5+v6);
	double f0,f1,f2,f3,f4,f5,f6,f7;
	//f3=v2*v4*v6/v;f0=v1*v4*v6/v;f2=v2*v3*v6/v;f1=v1*v3*v6/v;
	//f7=v2*v4*v5/v;f4=v1*v4*v5/v;f6=v2*v3*v5/v;f5=v1*v3*v5/v;
	f0=v2*v4*v6/v;f1=v1*v4*v6/v;f2=v1*v3*v6/v;f3=v2*v3*v6/v;
	f4=v2*v4*v5/v;f5=v1*v4*v5/v;f6=v1*v3*v5/v;f7=v2*v3*v5/v;
	result.x=f1+f2+f5+f6-f0-f3-f4-f7;
	result.y=f2+f3+f6+f7-f0-f1-f4-f5;
	result.z=f4+f5+f6+f7-f0-f1-f2-f3;
	//Begin the iteration
	iteration_for_std_coordinate(result,coordinate_physical);
	return result;
}
void Tcell_pure_ALE::iteration_for_std_coordinate(vec3D &iterator_coordinate,vec3D &coordinate_physical)
{
	double a[3][3],error,det;
	bool degenerate;
	vec3D trial_result;
	int num_step=0;
	trial_result=shape_function_interpolate(iterator_coordinate,"Position");
	error=sqrt((trial_result-coordinate_physical).self_multuply());
	while (error>1e-12)
	{
		for (int i=0;i<3;i++)
		{
			for (int j=0;j<3;j++)
			{
				a[i][j]=get_iteration_matrix(i,j,iterator_coordinate);
			}
		}
		inverse(a,det,degenerate);
		if (degenerate)
		{
			cout<<"The gradient matrix for shape function is degenerate!"<<endl;
			exit(0);
		}
		iterator_coordinate=iterator_coordinate-(trial_result-coordinate_physical).multiply_by_matrix(a);
		trial_result=shape_function_interpolate(iterator_coordinate,"Position");
		error=sqrt((trial_result-coordinate_physical).self_multuply());
		num_step=num_step+1;
		if (num_step>=20)
		{
			iterator_coordinate.value(0,0,0);
			cout<<"Warning:The iteration for standard coordinate is not convergence in calling iteration_for_std_coordinate()"<<endl;
			return;
		}
	}
	return;
}
double Tcell_pure_ALE::get_iteration_matrix(int i,int j,vec3D &iterator_coordinate)
{
	double v[8],f[8],aij=0;
	for (int n=0;n<8;n++)
	{
		v[n]=_node_ptr[n]->_position.access(i+1);
	}
	double xi=iterator_coordinate.x,eta=iterator_coordinate.y,zeta=iterator_coordinate.z;
	for (int n = 0; n < 8; n++)
	{
		f[n]=GN_brick_8(n,xi,eta,zeta).access(j+1);
	}
	for (int n=0;n<8;n++)
	{
		aij=aij+f[n]*v[n];
	}
	return aij;
}
vec3D Tcell_pure_ALE::shape_function_interpolate(vec3D std_coordinate,string name)
{
	vec3D v[8],result;
	if (name=="Position")
	{
		for (int j=0;j<8;j++)
		{
			v[j]=_node_ptr[j]->_position;
		}
	}
	else if (name=="Velocity")
	{
		for (int j=0;j<8;j++)
		{
			v[j]=_node_ptr[j]->_velocity;
		}
	}
	else
	{
		cout<<"Error:Invalid variable name in calling shape_function_interpolate()"<<endl;
		system("Pause");
		exit(0);
	}
	double xi=std_coordinate.x,eta=std_coordinate.y,zeta=std_coordinate.z;
	for (int j=0;j<8;j++)
	{
		result=result+v[j]*N_brick_8(j,xi,eta,zeta);
	}
	return result;
}
void Tcell_pure_ALE::calculate_velocity_gradient(vec3D &dvx,vec3D &dvy,vec3D &dvz)
{
	vec3D nf[6],vf[6];
	double sf[6];
	double s1,s2,s3,s4;
	vec3D n1,n2,n3,n4;
	for (int i=0;i<4;i++)
	{
		area(_cell_polyhedron.G_vertex(i+8),_cell_polyhedron.G_vertex(i),_cell_polyhedron.G_vertex((i+1)%4),s1,n1);
		area(_cell_polyhedron.G_vertex(i+8),_cell_polyhedron.G_vertex((i+1)%4),_cell_polyhedron.G_vertex((i+5)%4+4),s2,n2);
		area(_cell_polyhedron.G_vertex(i+8),_cell_polyhedron.G_vertex((i+5)%4+4),_cell_polyhedron.G_vertex(i+4),s3,n3);
		area(_cell_polyhedron.G_vertex(i+8),_cell_polyhedron.G_vertex(i+4),_cell_polyhedron.G_vertex(i),s4,n4);
		sf[i]=s1+s2+s3+s4;
		nf[i]=(n1*s1+n2*s2+n3*s3+n4*s4)/sf[i];
		vf[i]=(_node_ptr[i]->_velocity+_node_ptr[(i+1)%4]->_velocity+_node_ptr[(i+5)%4+4]->_velocity+_node_ptr[i+4]->_velocity)*0.25;
	}
	area(_cell_polyhedron.G_vertex(12),_cell_polyhedron.G_vertex(0),_cell_polyhedron.G_vertex(3),s1,n1);
	area(_cell_polyhedron.G_vertex(12),_cell_polyhedron.G_vertex(3),_cell_polyhedron.G_vertex(2),s2,n2);
	area(_cell_polyhedron.G_vertex(12),_cell_polyhedron.G_vertex(2),_cell_polyhedron.G_vertex(1),s3,n3);
	area(_cell_polyhedron.G_vertex(12),_cell_polyhedron.G_vertex(1),_cell_polyhedron.G_vertex(0),s4,n4);
	sf[4]=s1+s2+s3+s4;
	nf[4]=(n1*s1+n2*s2+n3*s3+n4*s4)/sf[4];
	vf[4]=(_node_ptr[0]->_velocity+_node_ptr[1]->_velocity+_node_ptr[2]->_velocity+_node_ptr[3]->_velocity)*0.25;
	area(_cell_polyhedron.G_vertex(13),_cell_polyhedron.G_vertex(4),_cell_polyhedron.G_vertex(5),s1,n1);
	area(_cell_polyhedron.G_vertex(13),_cell_polyhedron.G_vertex(5),_cell_polyhedron.G_vertex(6),s2,n2);
	area(_cell_polyhedron.G_vertex(13),_cell_polyhedron.G_vertex(6),_cell_polyhedron.G_vertex(7),s3,n3);
	area(_cell_polyhedron.G_vertex(13),_cell_polyhedron.G_vertex(7),_cell_polyhedron.G_vertex(4),s4,n4);
	sf[5]=s1+s2+s3+s4;
	nf[5]=(n1*s1+n2*s2+n3*s3+n4*s4)/sf[5];
	vf[5]=(_node_ptr[4]->_velocity+_node_ptr[5]->_velocity+_node_ptr[6]->_velocity+_node_ptr[7]->_velocity)*0.25;
	dvx.value(0,0,0);
	dvz=dvy=dvx;
	for (int i=0;i<6;i++)
	{
		dvx=dvx+nf[i]*sf[i]*vf[i].x;
		dvy=dvy+nf[i]*sf[i]*vf[i].y;
		dvz=dvz+nf[i]*sf[i]*vf[i].z;
	}
	dvx=dvx/_cell_volume;
	dvy=dvy/_cell_volume;
	dvz=dvz/_cell_volume;
	return;
}
void Tcell_pure_ALE::assemble_nodal_mass_moment()
{
	for (int i=0;i<8;i++)
	{
		_node_ptr[i]->_mass0=_node_ptr[i]->_mass0+_corner_mass[i];
		_node_ptr[i]->_mass=_node_ptr[i]->_mass+_corner_mass[i];
		_node_ptr[i]->_moment=_node_ptr[i]->_moment+_corner_moment[i];
	}
	return;
}