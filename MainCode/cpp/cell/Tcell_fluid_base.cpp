#include "Tcell_fluid_base.h"
#include "data_structure.h"
#include "public_function.h"
#include "Shape_function.h"
//--------------------------------------------------------------------------------
//Initialize the static variables
//-------------------------------------------------------------------------------
int Tcell_fluid_base::cell_face[6][4]={{0,3,2,1},{0,1,5,4},{1,2,6,5},{2,3,7,6},{0,4,7,3},{4,5,6,7}};
int Tcell_fluid_base::cell_edge_f[12][2]={{4,2},{1,3},{2,4},{3,1},{4,2},{1,3},{2,4},{3,1},{0,5},{0,5},{0,5},{0,5}};
int Tcell_fluid_base::cell_edge_f1[12][2]={{1,0},{2,0},{3,0},{4,0},{5,1},{5,2},{5,3},{5,4},{4,1},{1,2},{2,3},{3,4}};
int Tcell_fluid_base::cell_edge_v[12][2]={{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
int Tcell_fluid_base::cell_corner_v[8][3]={{1,3,4},{2,0,5},{3,1,6},{0,2,7},{7,5,0},{4,6,1},{5,7,2},{6,4,3}};
int Tcell_fluid_base::cell_corner_f[8][3]={{4,1,0},{1,2,0},{2,3,0},{3,4,0},{1,4,5},{2,1,5},{3,2,5},{4,3,5}};
Tcell_fluid_base::Tcell_fluid_base(int cell_id,int num_material,Tnode_fluid* (&node_ptr)[8],TMAT_base* mat)
	:Tcell_base(cell_id,num_material)
{
	//Allocate memory and initialize
	_corner_stress_force.resize(8);
	_cell_length=0;_flag=0;_volume_intersection=0;
	//Set cell node
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]=node_ptr[i];
	}
	//Initial adjacent
	for (int i = 0; i < 6; i++)
	{
		_adjacent[i]=NULL;
	}
	//Set the background material EOS
	if (mat!=NULL)
	{
		_gausspoint[0]=mat;
	}
	//Set geometry information
	bool negative;
	_cell_volume=calculate_cell_volume("Full",negative);
	if (negative)
	{
		cout<<"Error: Negative volume of the initial input"<<endl;
		cout<<"cell id="<<_id;
		system("Pause");
		exit(0);
	}
	create_cell_polyhedron();
	_flag=0;
}
Tcell_fluid_base::Tcell_fluid_base(int cell_id,int num_material,Tnode_fluid* (&node_ptr)[8],TMAT_base* mat[])
	:Tcell_base(cell_id,num_material)
{
	//Allocate memory and initialize
	_corner_stress_force.resize(8);
	_cell_length=0;_flag=0;_volume_intersection=0;
	//Set cell node
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]=node_ptr[i];
	}
	//Initial adjacent
	for (int i = 0; i < 6; i++)
	{
		_adjacent[i]=NULL;
	}
	//Set all the materials
	for (int i = 0; i < num_material; i++)
	{
		_gausspoint[i]=mat[i];
	}
	//Set geometry information
	bool negative;
	_cell_volume=calculate_cell_volume("Full",negative);
	if (negative)
	{
		cout<<"Error: Negative volume of the initial input"<<endl;
		cout<<"cell id="<<_id;
		system("Pause");
		exit(0);
	}
	create_cell_polyhedron();
	_flag=0;
}
void Tcell_fluid_base::calculate_corner_force(Tviscosity_base* Arti_Vis,double dt)
{
	vec3D assistant[8][3];      //The 24 assistant direction area pointing from i to cell_corner_v[i][j]
	vec3D surfacial[8][3];      //The 24 surfacial direction area at cell_corner_f[i][j]
	calculate_auxiliary_area(surfacial,assistant);
	calculate_corner_stress_force(surfacial,assistant,dt);
	calculate_corner_hourglass_force(surfacial,assistant);
	Arti_Vis->calculate_viscosity_corner_force(this,assistant,dt);
	return;
}
void Tcell_fluid_base::calculate_corner_force_without_hourglass(Tviscosity_base* Arti_Vis,double dt)
{
	vec3D assistant[8][3];      //The 24 assistant direction area pointing from i to cell_corner_v[i][j]
	vec3D surfacial[8][3];      //The 24 surfacial direction area at cell_corner_f[i][j]
	calculate_auxiliary_area(surfacial,assistant);
	calculate_corner_stress_force(surfacial,assistant,dt);
	Arti_Vis->calculate_viscosity_corner_force(this,assistant,dt);
	return;
}
void Tcell_fluid_base::calculate_viscosity_force()
{
	//Calculate the gradient of the velocity
	vec3D GNshape[8];
	double detJ;
	calculate_DN_matrix(GNshape,detJ,0,0,0);
	double de[6]={0};
	//de[0]:ux,x;de[1]:uy,y;de[2]:uz,z;de[3]:uy,z+uz,y;de[4]:uz,x+ux,z;de[5]:ux,y+uy,x;
	for (int i = 0; i < 8; i++)
	{
		de[0]=de[0]+_node_ptr[i]->_velocity.x*GNshape[i].x;
		de[1]=de[1]+_node_ptr[i]->_velocity.y*GNshape[i].y;
		de[2]=de[2]+_node_ptr[i]->_velocity.z*GNshape[i].z;
		de[3]=de[3]+_node_ptr[i]->_velocity.y*GNshape[i].z+_node_ptr[i]->_velocity.z*GNshape[i].y;
		de[4]=de[4]+_node_ptr[i]->_velocity.z*GNshape[i].x+_node_ptr[i]->_velocity.x*GNshape[i].z;
		de[5]=de[5]+_node_ptr[i]->_velocity.x*GNshape[i].y+_node_ptr[i]->_velocity.y*GNshape[i].x;
	}
	//Calculate the stress from fluid viscosity
	double mu=1.25e-4,k=-2.0/3.0;
	double sxx,syy,szz,sxy,syz,sxz;
	double divu=de[0]+de[1]+de[2];
	sxx=2*mu*de[0]+k*mu*divu;syy=2*mu*de[1]+k*mu*divu;szz=2*mu*de[2]+k*mu*divu;
	syz=mu*de[3];sxz=mu*de[4];sxy=mu*de[5];
	//Calculate the corner force from viscosity force
	for (int i = 0; i < 8; i++)
	{
		_corner_stress_force[i].x=_corner_stress_force[i].x-(sxx*GNshape[i].x+sxy*GNshape[i].y+sxz*GNshape[i].z)*_cell_volume;
		_corner_stress_force[i].y=_corner_stress_force[i].y-(sxy*GNshape[i].x+syy*GNshape[i].y+syz*GNshape[i].z)*_cell_volume;
		_corner_stress_force[i].z=_corner_stress_force[i].z-(sxz*GNshape[i].x+syz*GNshape[i].y+szz*GNshape[i].z)*_cell_volume;
	}
	return;
}
double Tcell_fluid_base::calculate_time_step()
{
	double v,v_max;                    //Calculate the velocity of node to get the maximal velocity
	v_max=0;
	for (int i=0;i<8;i++)
	{
		v=_node_ptr[i]->_velocity*_node_ptr[i]->_velocity;
		if (v>v_max)
		{
			v_max=v;
		}
	}
	//Calcualte characteristic length when calculating time step
	_cell_length=calculate_characteristic_length();
	double max_sound_speed=0;
	int num_material=_gausspoint.size();
	if (num_material==1)
	{
		max_sound_speed=_gausspoint[0]->G_soundspeed();
	}
	else
	{
		for (int i = 0; i < num_material; i++)
		{
			max_sound_speed=maxval(max_sound_speed,_gausspoint[i]->G_soundspeed());
			//max_sound_speed=max_sound_speed+_gausspoint[i]->G_soundspeed()*G_material_fraction(i);
		}
	}
	return _cell_length/(max_sound_speed+sqrt(v_max));
}
void Tcell_fluid_base::assemble_corner_force()
{
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]->_force=_node_ptr[i]->_force+_corner_stress_force[i]+_corner_hourglass_force[i]+_corner_viscosity_force[i];
		_node_ptr[i]->_force0=_node_ptr[i]->_force0+_corner_stress_force[i]+_corner_hourglass_force[i]+_corner_viscosity_force[i];
	}
	return;
}
void Tcell_fluid_base::assemble_nodal_mass_moment()
{
	for (int i=0;i<8;i++)
	{
		_node_ptr[i]->_mass0=_node_ptr[i]->_mass0+_corner_mass[i];
		_node_ptr[i]->_mass=_node_ptr[i]->_mass+_corner_mass[i];
		_node_ptr[i]->_moment=_node_ptr[i]->_moment+_corner_moment[i];
	}
	return;
}
double Tcell_fluid_base::calculate_cell_volume(string type,bool &negative)
{
	vec3D c_cell,c_face[6],c_edge[12];     //The cell center,face center and edge center
	int n1,n2,s1,s2;
	double temp;
	calculate_segment_center(c_cell,c_face,c_edge,type);
	negative=false;
	if (type=="Full" || type=="Trial")
	{
		for (int i=0;i<8;i++)
		{
			_corner_volume[i]=0;
		}
		for (int i=0;i<12;i++)
		{
			n1=cell_edge_v[i][0];n2=cell_edge_v[i][1];
			s1=cell_edge_f1[i][0];s2=cell_edge_f1[i][1];
			temp=((c_face[s1]-c_cell)%(c_edge[i]-c_cell)+(c_edge[i]-c_cell)%(c_face[s2]-c_cell))*(_node_ptr[n2]->_position-_node_ptr[n1]->_position);
			temp=temp/12.0;
			_corner_volume[n1]=_corner_volume[n1]+temp;
			_corner_volume[n2]=_corner_volume[n2]+temp;
		}
		if (type=="Full")
		{
			_cell_volume=0;
		}
		for (int i=0;i<8;i++)
		{
			if (_corner_volume[i]<0)
			{
				negative=true;
			}
			if (type=="Full")
			{
				_cell_volume=_cell_volume+_corner_volume[i];
			}
		}
		return _cell_volume;
	}
	else if (type=="Temporary")
	{
		double temp_volume=0;
		for (int i=0;i<12;i++)
		{
			n1=cell_edge_v[i][0];n2=cell_edge_v[i][1];
			s1=cell_edge_f1[i][0];s2=cell_edge_f1[i][1];
			temp=((c_face[s1]-c_cell)%(c_edge[i]-c_cell)+(c_edge[i]-c_cell)%(c_face[s2]-c_cell))*(_node_ptr[n2]->_position_temp-_node_ptr[n1]->_position_temp);
			temp=temp/6.0;
			temp_volume=temp_volume+temp;
		}
		if (temp_volume<=0)
		{
			negative=true;
		}
		return temp_volume;
	}
	else
	{
		cout<<"Error:Invalid type in calling calculate_cell_volume()"<<endl;
		system("Pause");
		exit(0);
	}
	return 0;
}
void Tcell_fluid_base::create_cell_polyhedron()
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
void Tcell_fluid_base::update_cell_polyhedron()
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
void Tcell_fluid_base::clear_cell_variable()
{
	int num_material=_gausspoint.size();
	for (int i = 0; i < num_material; i++)
	{
		_gausspoint[i]->S_density(0);
		_gausspoint[i]->S_internal_energy(0);
		_gausspoint[i]->S_pressure(0);
		_gausspoint[i]->S_soundspeed(0);
	}
	_cell_volume=_cell_length=_volume_intersection=0;                
	for (int i = 0; i < 8; i++)
	{
		_corner_mass[i]=_corner_volume[i]=0;
		_corner_hourglass_force[i].value(0,0,0); 
		_corner_viscosity_force[i].value(0,0,0);
		_corner_moment[i].value(0,0,0);
	}     
	return;     
}
void Tcell_fluid_base::calculate_coordinate_range()
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
double Tcell_fluid_base::calculate_maximal_edge()
{
	//double edge_max=0;
	//double edge[4];
	//double l_max=(_coor_max-_coor_min)*(_coor_max-_coor_min);
	//edge[0]=(_node_ptr[0]->_position-_node_ptr[6]->_position).self_multuply();
	//edge[1]=(_node_ptr[1]->_position-_node_ptr[7]->_position).self_multuply();
	//edge[2]=(_node_ptr[2]->_position-_node_ptr[4]->_position).self_multuply();
	//edge[3]=(_node_ptr[3]->_position-_node_ptr[5]->_position).self_multuply();
	//for (int i=0;i<4;i++)
	//{
	//	if (edge[i]>edge_max)
	//	{
	//		edge_max=edge[i];
	//	}
	//}
	return (_coor_max-_coor_min).get_length();
}
vec3D Tcell_fluid_base::calculate_std_coordinate(vec3D coordinate_physical)
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
	f0=v2*v4*v6/v;f1=v1*v4*v6/v;f2=v1*v3*v6/v;f3=v2*v3*v6/v;
	f4=v2*v4*v5/v;f5=v1*v4*v5/v;f6=v1*v3*v5/v;f7=v2*v3*v5/v;
	result.x=f1+f2+f5+f6-f0-f3-f4-f7;
	result.y=f2+f3+f6+f7-f0-f1-f4-f5;
	result.z=f4+f5+f6+f7-f0-f1-f2-f3;
	//Begin the iteration
	iteration_for_std_coordinate(result,coordinate_physical);
	return result;
}
void Tcell_fluid_base::iteration_for_std_coordinate(vec3D &iterator_coordinate,vec3D &coordinate_physical)
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
double Tcell_fluid_base::get_iteration_matrix(int i,int j,vec3D &iterator_coordinate)
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
vec3D Tcell_fluid_base::shape_function_interpolate(vec3D std_coordinate,string name)
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
void Tcell_fluid_base::calculate_velocity_gradient(vec3D &dvx,vec3D &dvy,vec3D &dvz)
{
	//vec3D nf[6],vf[6];
	//double sf[6];
	//double s1,s2,s3,s4;
	//vec3D n1,n2,n3,n4;
	//for (int i=0;i<4;i++)
	//{
	//	area(_cell_polyhedron.G_vertex(i+8),_cell_polyhedron.G_vertex(i),_cell_polyhedron.G_vertex((i+1)%4),s1,n1);
	//	area(_cell_polyhedron.G_vertex(i+8),_cell_polyhedron.G_vertex((i+1)%4),_cell_polyhedron.G_vertex((i+5)%4+4),s2,n2);
	//	area(_cell_polyhedron.G_vertex(i+8),_cell_polyhedron.G_vertex((i+5)%4+4),_cell_polyhedron.G_vertex(i+4),s3,n3);
	//	area(_cell_polyhedron.G_vertex(i+8),_cell_polyhedron.G_vertex(i+4),_cell_polyhedron.G_vertex(i),s4,n4);
	//	sf[i]=s1+s2+s3+s4;
	//	nf[i]=(n1*s1+n2*s2+n3*s3+n4*s4)/sf[i];
	//	vf[i]=(_node_ptr[i]->_velocity+_node_ptr[(i+1)%4]->_velocity+_node_ptr[(i+5)%4+4]->_velocity+_node_ptr[i+4]->_velocity)*0.25;
	//}
	//area(_cell_polyhedron.G_vertex(12),_cell_polyhedron.G_vertex(0),_cell_polyhedron.G_vertex(3),s1,n1);
	//area(_cell_polyhedron.G_vertex(12),_cell_polyhedron.G_vertex(3),_cell_polyhedron.G_vertex(2),s2,n2);
	//area(_cell_polyhedron.G_vertex(12),_cell_polyhedron.G_vertex(2),_cell_polyhedron.G_vertex(1),s3,n3);
	//area(_cell_polyhedron.G_vertex(12),_cell_polyhedron.G_vertex(1),_cell_polyhedron.G_vertex(0),s4,n4);
	//sf[4]=s1+s2+s3+s4;
	//nf[4]=(n1*s1+n2*s2+n3*s3+n4*s4)/sf[4];
	//vf[4]=(_node_ptr[0]->_velocity+_node_ptr[1]->_velocity+_node_ptr[2]->_velocity+_node_ptr[3]->_velocity)*0.25;
	//area(_cell_polyhedron.G_vertex(13),_cell_polyhedron.G_vertex(4),_cell_polyhedron.G_vertex(5),s1,n1);
	//area(_cell_polyhedron.G_vertex(13),_cell_polyhedron.G_vertex(5),_cell_polyhedron.G_vertex(6),s2,n2);
	//area(_cell_polyhedron.G_vertex(13),_cell_polyhedron.G_vertex(6),_cell_polyhedron.G_vertex(7),s3,n3);
	//area(_cell_polyhedron.G_vertex(13),_cell_polyhedron.G_vertex(7),_cell_polyhedron.G_vertex(4),s4,n4);
	//sf[5]=s1+s2+s3+s4;
	//nf[5]=(n1*s1+n2*s2+n3*s3+n4*s4)/sf[5];
	//vf[5]=(_node_ptr[4]->_velocity+_node_ptr[5]->_velocity+_node_ptr[6]->_velocity+_node_ptr[7]->_velocity)*0.25;
	//dvx.value(0,0,0);
	//dvz=dvy=dvx;
	//for (int i=0;i<6;i++)
	//{
	//	dvx=dvx+nf[i]*sf[i]*vf[i].x;
	//	dvy=dvy+nf[i]*sf[i]*vf[i].y;
	//	dvz=dvz+nf[i]*sf[i]*vf[i].z;
	//}
	//dvx=dvx/_cell_volume;
	//dvy=dvy/_cell_volume;
	//dvz=dvz/_cell_volume;
	//return;
}
int Tcell_fluid_base::vertex_reverse(int nid)
{
	for (int i=0;i<8;i++)
	{
		if (_node_ptr[i]->_id==nid)
		{
			return i;
		}
	}
	cout<<"Error: this node isn't belong to this cell!!"<<endl;
	system("Pause");
	exit(0);
	return 0;
}
int Tcell_fluid_base::cell_corner_v_reverse(int p0,int p1)
{
	for (int i=0;i<3;i++)
	{
		if (p1==cell_corner_v[p0][i])
		{
			return i;
		}
	}
	cout<<"Error:This node isn't belong to this corner!!!"<<endl;
	system("Pause");
	exit(0);
	return 0;
}
Tnode_fluid* Tcell_fluid_base::access_extended_node(int local_face_id,int local_node_id)
{
	int global_id_on_face[4];            //The global number of nodes on this face
	int global_id_target;                //The global id of _node_ptr[local_node_id]
	int global_id_trial;                 //The global id of node on the adjacent cell
	Tcell_fluid_base* adjacent_cell;    //The adjacent cell to this face 
	int local_id;                        //The local node id of _node_ptr[local_node_id] on the adjacent cell
	int temp;
	for (int i=0;i<4;i++)
	{
		global_id_on_face[i]=_node_ptr[cell_face[local_face_id][i]]->_id;
	}
	global_id_target=_node_ptr[local_node_id]->_id;
	//If the node is located on boundary, return null
	if (_adjacent[local_face_id]==NULL)
	{
		return NULL;
	}
	else
	{
		adjacent_cell=_adjacent[local_face_id];
	}
	local_id=adjacent_cell->vertex_reverse(global_id_target);
	for (int i=0;i<3;i++)
	{
		temp=0;
		global_id_trial=adjacent_cell->_node_ptr[cell_corner_v[local_id][i]]->_id;
		for (int j=0;j<4;j++)
		{
			if (global_id_trial==global_id_on_face[j])
			{
				temp=temp+1;
			}
		}
		if (temp==0)
		{
			return adjacent_cell->_node_ptr[cell_corner_v[local_id][i]];
		}
	}
	cout<<"Error: the adjacent node isn't found"<<endl;
	system("Pause");
	exit(0);
	return 0;
}
void Tcell_fluid_base::set_cell_connection()
{
	//Calculate the cells which has at least one common node to this cell
	_flag=1;
	int num_surround;
	for (int i=0;i<8;i++)
	{
		num_surround=_node_ptr[i]->_surroundings.size();
		for (int j=0;j<num_surround;j++)
		{
			if (_node_ptr[i]->_surroundings[j]->_flag==0)
			{
				_around.push_back(_node_ptr[i]->_surroundings[j]);
				//Set flag=1 will ensure this cell only being pushed back once
				_node_ptr[i]->_surroundings[j]->_flag=1;
			}
		}
	}
	_flag=0;
	int num_around=_around.size();
	for (int i=0;i<8;i++)
	{
		num_surround=_node_ptr[i]->_surroundings.size();
		for (int j=0;j<num_surround;j++)
		{
			//Recover the flags
			_node_ptr[i]->_surroundings[j]->_flag=0;
		}
	}
	//Calculate the cells which is adjacent to this cell
	for (int i=0;i<6;i++)
	{
		for (int j=0;j<4;j++)
		{
			_node_ptr[cell_face[i][j]]->_flag=1;
		}
		for (int j=0;j<num_around;j++)
		{
			int temp=0;
			for (int k=0;k<8;k++)
			{
				temp=temp+_around[j]->_node_ptr[k]->_flag;
			}
			if (temp==4)
			{
				_adjacent[i]=_around[j];
				break;
			}
		}
		for (int j=0;j<4;j++)
		{
			_node_ptr[cell_face[i][j]]->_flag=0;
		}
	}
	return;
}
void Tcell_fluid_base::add_to_node()
{
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]->add_connected_cell(this);
	}
	return;
}
void Tcell_fluid_base::calculate_distorted_ratio(double &ratio_angle,double &ratio_length)
{
	double l_max=0;             //The maximum length of a cell
	double l_min=1000000;       //The minimum length of a cell
	ratio_angle=0;              //The maximum cos of angle
	vec3D p1,p2,p3;
	int i1,i2,i3;
	//Calculate ratio_length
	for (int i=0;i<12;i++)
	{
		i1=cell_edge_v[i][0];i2=cell_edge_v[i][1];
		p1=_node_ptr[i1]->_position;p2=_node_ptr[i2]->_position;
		double l=sqrt((p1-p2)*(p1-p2));
		if (l<l_min) l_min=l;
		if (l>l_max) l_max=l;
	}
	ratio_length=l_max/l_min;
	//Calculate ratio angle
	for (int i=0;i<8;i++)
	{
		i1=i;
		p1=_node_ptr[i1]->_position;
		for (int j=0;j<3;j++)
		{
			i2=cell_corner_v[i][j];i3=cell_corner_v[i][(j+1)%3];
			p2=_node_ptr[i2]->_position;p3=_node_ptr[i3]->_position;
			double l1=sqrt((p2-p1)*(p2-p1));
			double l2=sqrt((p3-p1)*(p3-p1));
			double dot=abs((p3-p1)*(p2-p1));
			double angle=dot/(l1*l2);
			if (angle>ratio_angle) ratio_angle=angle;
		}
	}
	return;
}
void Tcell_fluid_base::calculate_auxiliary_area(vec3D (&surfacial)[8][3],vec3D (&assistant)[8][3])
{
	vec3D c_cell,c_face[6],c_edge[12];                    //The center of the cell
	int n0,n1,n2,s1,s2;             
	vec3D b0,b1,b2;
	vec3D temp;

	calculate_segment_center(c_cell,c_face,c_edge,"Full");
	//Calculate the surfacial direction
	for (int i=0;i<8;i++)
	{
		for (int j=0;j<3;j++)
		{
			n1=Tcell_fluid_base::cell_corner_v[i][(j+1)%3];
			n2=Tcell_fluid_base::cell_corner_v[i][(j+2)%3];
			s1=Tcell_fluid_base::cell_corner_f[i][j];
			b0=_node_ptr[i]->_position;
			b1=(_node_ptr[n1]->_position+b0)*0.5;
			b2=(_node_ptr[n2]->_position+b0)*0.5;
			surfacial[i][j]=((b2-b0)%(c_face[s1]-b0)+(c_face[s1]-b0)%(b1-b0))*0.5;
		}
	}
	//Calculate the assistant direction
	for (int i=0;i<12;i++)
	{
		n1=Tcell_fluid_base::cell_edge_v[i][0];n2=Tcell_fluid_base::cell_edge_v[i][1];
		s1=Tcell_fluid_base::cell_edge_f1[i][0];s2=Tcell_fluid_base::cell_edge_f1[i][1];
		temp=((c_face[s1]-c_cell)%(c_edge[i]-c_cell)+(c_edge[i]-c_cell)%(c_face[s2]-c_cell))*0.5;
		n0=cell_corner_v_reverse(n1,n2);
		assistant[n1][n0]=temp;
		n0=cell_corner_v_reverse(n2,n1);
		assistant[n2][n0]=temp*(-1);
	}
	return;
}
double Tcell_fluid_base::calculate_characteristic_length()
{
	double l,l_min=100000;
	for (int i = 0; i < 4; i++)
	{
		l=(_node_ptr[i]->_position-_node_ptr[(i+1)%4]->_position).self_multuply();
		l_min=minval(l,l_min);
		l=(_node_ptr[i+4]->_position-_node_ptr[(i+1)%4+4]->_position).self_multuply();
		l_min=minval(l,l_min);
		l=(_node_ptr[i]->_position-_node_ptr[i+4]->_position).self_multuply();
		l_min=minval(l,l_min);
	}
	return sqrt(l_min);
}
void Tcell_fluid_base::calculate_segment_center(vec3D &c_cell,vec3D (&c_face)[6],vec3D (&c_edge)[12],string type)
{
	vec3D node_pos[8];
	if (type=="Temporary")
	{
		for (int i = 0; i < 8; i++)
		{
			node_pos[i]=_node_ptr[i]->_position_temp;
		}
	}
	else if (type=="Full")
	{
		for (int i = 0; i < 8; i++)
		{
			node_pos[i]=_node_ptr[i]->_position;
		}
	}
	else if (type=="Trial")
	{
		for (int i = 0; i < 8; i++)
		{
			node_pos[i]=_node_ptr[i]->_position_trial;
			//cout<<_node_ptr[i]->_position_trial.x<<" "<<_node_ptr[i]->_position_trial.y<<" "<<_node_ptr[i]->_position_trial.z<<endl;
			//cout<<_node_ptr[i]->_position_trial.x<<" "<<_node_ptr[i]->_position_trial.y<<" "<<_node_ptr[i]->_position_trial.z<<endl;
		}
	
	}
	else
	{
		cout<<"Error: the type in Function \"calculate_segment_center\" must be \"Full\" or \"Temporary\"!"<<endl;
		system("Pause");
		exit(0);
	}
	c_cell.value(0,0,0);
	 for (int i=0;i<8;i++)
	 {
		 c_cell=c_cell+node_pos[i];	
	 }
	 c_cell=c_cell*0.125;
	 for (int i=0;i<6;i++)
	 {
		 c_face[i].value(0,0,0);
		 for (int j=0;j<4;j++)
		 {
			 c_face[i]=c_face[i]+node_pos[cell_face[i][j]];
		 }
		 c_face[i]=c_face[i]*0.25;
	 }
	 for (int i=0;i<12;i++)
	 {
		 c_edge[i]=(node_pos[cell_edge_v[i][0]]+node_pos[cell_edge_v[i][1]])*0.5;
	 }
	 return;
}
void Tcell_fluid_base::reset_TR_conditions(string type)
{
	double pi=3.141592654;
	if (type=="2D")
	{
		double h,x,y;
		double density,pressure,internal_energy,soundspeed;
		bool failed;
		for (int i = 0; i < 2; i++)
		{
			if (G_material_fraction(i)>0)
			{
				x=G_material_centroid(i).x;y=G_material_centroid(i).y;
				h=0.5+0.01*cos(6*pi*x);
				density=_gausspoint[i]->G_density();
				if (i==0)
				{
					pressure=1-0.2*(h-1)-0.1*(y-h);
				}
				else
				{
					pressure=pressure=1-0.2*(y-1);
				}
				internal_energy=_gausspoint[i]->G_EOS()->p_ro_e(pressure,density);
				_gausspoint[i]->S_pressure(pressure);
				_gausspoint[i]->S_internal_energy(internal_energy);
				soundspeed=_gausspoint[i]->calculate_sound_speed(failed);
				_gausspoint[i]->S_soundspeed(soundspeed);
			}
		}
	}
	else
	{
		double h,x,y,z;
		double density,pressure,internal_energy,soundspeed;
		bool failed;
		for (int i = 0; i < 2; i++)
		{
			if (G_material_fraction(i)>0)
			{
				x=G_material_centroid(i).x;y=G_material_centroid(i).y;z=G_material_centroid(i).z;
				h=0.5+0.05*(cos(8*pi*x)+cos(8*pi*y))/4.0;
				density=_gausspoint[i]->G_density();
				if (i==0)
				{
					pressure=1-0.2*(h-1)-0.1*(z-h);
				}
				else
				{
					pressure=pressure=1-0.2*(h-1);
				}
				internal_energy=_gausspoint[i]->G_EOS()->p_ro_e(pressure,density);
				_gausspoint[i]->S_pressure(pressure);
				_gausspoint[i]->S_internal_energy(internal_energy);
				soundspeed=_gausspoint[i]->calculate_sound_speed(failed);
				_gausspoint[i]->S_soundspeed(soundspeed);
			}
		}
	}
	return;
}
void Tcell_fluid_base::reset_dam_breaking_conditions(double dz,int nz,int mat_id,int num_material)
{
	if (num_material==1 || G_material_fraction(mat_id)==1)
	{
		double ro0=_gausspoint[mat_id]->G_density();
		double ro1=2500*ro0/(2500-4.9*dz);
		double z=_cell_polyhedron.G_centroid().y;
		z=z+0.5*dz;
		int n=find_int(z/dz);
		//int n=1;
		n=nz-n;
		double temp=(2500+4.9*dz)/(2500-4.9*dz);
		double density=pow(temp,n)*ro1;
		_gausspoint[mat_id]->S_density(density);
		double pressure=_gausspoint[mat_id]->calculate_pressure();
		bool fail;
		double soundspeed=_gausspoint[mat_id]->calculate_sound_speed(fail);
		_gausspoint[mat_id]->S_pressure(pressure);
		_gausspoint[mat_id]->S_soundspeed(soundspeed);
		double mass=density*_cell_volume;
		S_mass(mat_id,mass);
		for (int i = 0; i < 8; i++)
		{
			_corner_mass[i]=_corner_volume[i]*density;
		}
	}
	return;
}
void Tcell_fluid_base::calculate_DN_matrix(vec3D (&DN)[8],double &detJ,double xi,double eta,double zeta)
{
	double Jacobi[3][3];
	vec3D GN[8];
	bool degenerate;
	calcualte_Jacobi(Jacobi,GN,xi,eta,zeta);
	inverse(Jacobi,detJ,degenerate);
	for (int i = 0; i < 8; i++)
	{
		DN[i]=GN[i].multiply_by_matrix(Jacobi);
	}
	return;
}
void Tcell_fluid_base::calcualte_Jacobi(double (&Jacobi)[3][3],vec3D (&GN)[8],double xi,double eta,double zeta)
{
	for (int i = 0; i < 8; i++)
	{
		GN[i]=GN_brick_8(i,xi,eta,zeta);
	}
	vec3D temp[3];
	//Calculate the Jacobi matrix with respected to the current position
	for (int i=0;i<8;i++)
	{
		for (int j = 0; j < 3; j++)
		{
			temp[j]=temp[j]+GN[i]*_node_ptr[i]->_position.access(j+1);
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Jacobi[i][j]=temp[j].access(i+1);
		}
	}
	return;
}
