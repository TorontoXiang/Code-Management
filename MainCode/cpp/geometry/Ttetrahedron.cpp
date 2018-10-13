#include "Ttetrahedron.h"
#include "Ttriangle.h"
#include <cmath>
#include <iostream>
#include "public_function.h"
Ttetrahedron::Ttetrahedron(vec3D p1,vec3D p2,vec3D p3,vec3D p4):_v1(p1),_v2(p2),_v3(p3),_v4(p4),_distort_ratio(0)
{
	//Calculate the transfer matrix
	double a11=_v2.x-_v1.x,a12=_v3.x-_v1.x,a13=_v4.x-_v1.x;
	double a21=_v2.y-_v1.y,a22=_v3.y-_v1.y,a23=_v4.y-_v1.y;
	double a31=_v2.z-_v1.z,a32=_v3.z-_v1.z,a33=_v4.z-_v1.z;
	double det=(a11*a22*a33-a11*a23*a32-a12*a21*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31);
	_transfer_matrix[0][0]=(a22*a33-a23*a32)/det,_transfer_matrix[0][1]=(a13*a32-a12*a33)/det,_transfer_matrix[0][2]=(a12*a23-a13*a22)/det;
	_transfer_matrix[1][0]=(a23*a31-a21*a33)/det,_transfer_matrix[1][1]=(a11*a33-a13*a31)/det,_transfer_matrix[1][2]=(a13*a21-a11*a23)/det;
	_transfer_matrix[2][0]=(a21*a32-a22*a31)/det,_transfer_matrix[2][1]=(a12*a31-a11*a32)/det,_transfer_matrix[2][2]=(a11*a22-a12*a21)/det;
	//Calculate the volume
	double x1=_v1.x,x2=_v2.x,x3=_v3.x,x4=_v4.x;
	double y1=_v1.y,y2=_v2.y,y3=_v3.y,y4=_v4.y;
	double z1=_v1.z,z2=_v2.z,z3=_v3.z,z4=_v4.z;
	_volume=abs(x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4
		      -x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1
			  -x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2)/6.0;
}
Ttetrahedron::Ttetrahedron(vec3D p1,vec3D p2,vec3D p3,vec3D p4,int temp):_v1(p1),_v2(p2),_v3(p3),_v4(p4)
{
	//Calculate the volume
	double x1=_v1.x,x2=_v2.x,x3=_v3.x,x4=_v4.x;
	double y1=_v1.y,y2=_v2.y,y3=_v3.y,y4=_v4.y;
	double z1=_v1.z,z2=_v2.z,z3=_v3.z,z4=_v4.z;
	_volume=abs(x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4
		      -x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1
			  -x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2)/6.0;
	//Calculate surfical area
	double s[4],nz;
	area(_v1,_v2,_v3,s[0],nz);area(_v1,_v2,_v4,s[1],nz);
	area(_v2,_v3,_v4,s[2],nz);area(_v1,_v3,_v4,s[3],nz);
	double rt=3*_volume/(s[0]+s[1]+s[2]+s[3]);
	double l_max=0;
	l_max=maxval(l_max,(_v1-_v2).self_multuply());l_max=maxval(l_max,(_v1-_v3).self_multuply());
	l_max=maxval(l_max,(_v1-_v4).self_multuply());l_max=maxval(l_max,(_v4-_v2).self_multuply());
	l_max=maxval(l_max,(_v3-_v2).self_multuply());l_max=maxval(l_max,(_v3-_v4).self_multuply());
	_distort_ratio=sqrt(l_max)/rt;
}
void Ttetrahedron::transfer_vertex(vec3D &p)
{
	double x=p.x-_v1.x,y=p.y-_v1.y,z=p.z-_v1.z;
	p.x=_transfer_matrix[0][0]*x+_transfer_matrix[0][1]*y+_transfer_matrix[0][2]*z;
	p.y=_transfer_matrix[1][0]*x+_transfer_matrix[1][1]*y+_transfer_matrix[1][2]*z;
	p.z=_transfer_matrix[2][0]*x+_transfer_matrix[2][1]*y+_transfer_matrix[2][2]*z;
}
void Ttetrahedron::intersect_with_polyhedron(Tpolyhedron* mp,double& ivolume,double& ix,double& iy,double& iz )
{
	ivolume=ix=iy=iz=0;
	vec3D p1,p2,p3;
	vec3D integration_standard,integration_physical;
	double temp_v,temp_ix,temp_iy,temp_iz;
	for (int i=0;i<mp->G_nump();i++)
	{
		p1=mp->G_vertex_on_piece(i,0);p2=mp->G_vertex_on_piece(i,1);p3=mp->G_vertex_on_piece(i,2);
		transfer_vertex(p1);transfer_vertex(p2);transfer_vertex(p3);
		Ttriangle face(p1,p2,p3);
		face.column_integration(temp_v,temp_ix,temp_iy,temp_iz);
		ivolume=ivolume+temp_v;ix=ix+temp_ix;iy=iy+temp_iy;iz=iz+temp_iz;
	}
	integration_standard.x=ix;integration_standard.y=iy;integration_standard.z=iz;
	//In order to avoid divide operator because ivolume may equal to zero
	integration_physical=transfer_result(integration_standard,ivolume);
	ivolume=6*abs(ivolume)*_volume;
	ix=integration_physical.x*_volume;iy=integration_physical.y*_volume;iz=integration_physical.z*_volume;
	return;
}
void Ttetrahedron::calculate_vertex_range(vec3D &coor_min,vec3D &coor_max)
{
	for (int i = 1; i < 4; i++)
	{
		coor_min.value(i,minval(minval(_v1.access(i),_v2.access(i)),minval(_v3.access(i),_v4.access(i))));
		coor_max.value(i,maxval(maxval(_v1.access(i),_v2.access(i)),maxval(_v3.access(i),_v4.access(i))));
	}
	return;
}
vec3D Ttetrahedron::transfer_result(vec3D &p,double& ivolume)
{
	vec3D temp;
	double a11=_v2.x-_v1.x,a12=_v3.x-_v1.x,a13=_v4.x-_v1.x;
	double a21=_v2.y-_v1.y,a22=_v3.y-_v1.y,a23=_v4.y-_v1.y;
	double a31=_v2.z-_v1.z,a32=_v3.z-_v1.z,a33=_v4.z-_v1.z;
	temp.x=a11*p.x+a12*p.y+a13*p.z;
	temp.y=a21*p.x+a22*p.y+a23*p.z;
	temp.z=a31*p.x+a32*p.y+a33*p.z;
	if (ivolume<0)
	{
		temp=temp*(-1);
	}
	temp=temp*6+_v1*6*abs(ivolume);
	return temp;
}
vec3D Ttetrahedron::G_vertex(int i)
{
	if (i==0)
	{
		return _v1;
	}
	else if (i==1)
	{
		return _v2;
	}
	else if (i==2)
	{
		return _v3;
	}
	else if (i==3)
	{
		return _v4;
	}
	else
	{
		cout<<"Error: the vertex id of a tetrahedron is larger than 4"<<endl;
		system("Pause");
		exit(0);
	}
}

