#include "public_function.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
vec3D gauss_for_centroid(vec3D &p1,vec3D &p2,vec3D &p3)
{
	double x_m1=(p1.x+p2.x)*0.5,y_m1=(p1.y+p2.y)*0.5,z_m1=(p1.z+p2.z)*0.5;
	double x_m2=(p1.x+p3.x)*0.5,y_m2=(p1.y+p3.y)*0.5,z_m2=(p1.z+p3.z)*0.5;
	double x_m3=(p3.x+p2.x)*0.5,y_m3=(p3.y+p2.y)*0.5,z_m3=(p3.z+p2.z)*0.5;
	vec3D result;
	result.x=(x_m1*z_m1+x_m2*z_m2+x_m3*z_m3)/3.0;
	result.y=(y_m1*z_m1+y_m2*z_m2+y_m3*z_m3)/3.0;
	result.z=(z_m1*z_m1+z_m2*z_m2+z_m3*z_m3)/6.0;
	return result;
}
vec3D gauss_for_MoF_derivative(vec3D &p1,vec3D &p2,vec3D &p3)
{
	double x_m1=(p1.x+p2.x)*0.5,y_m1=(p1.y+p2.y)*0.5,z_m1=(p1.z+p2.z)*0.5;
	double x_m2=(p1.x+p3.x)*0.5,y_m2=(p1.y+p3.y)*0.5,z_m2=(p1.z+p3.z)*0.5;
	double x_m3=(p3.x+p2.x)*0.5,y_m3=(p3.y+p2.y)*0.5,z_m3=(p3.z+p2.z)*0.5;
	vec3D result;
	result.x=(x_m1*y_m1+x_m2*y_m2+x_m3*y_m3)/3.0;
	result.y=(y_m1*y_m1+y_m2*y_m2+y_m3*y_m3)/3.0;
	result.z=(x_m1*x_m1+x_m2*x_m2+x_m3*x_m3)/3.0;
	return result;
}
double gauss(vec3D &p1,vec3D &p2,vec3D &p3,string type)
{
	//Integral functions for calculating volume and centroid
	if (type=="z")
	{
		return (p1.z+p2.z+p3.z)/3.0;
	} 
	else if (type=="xz" || type=="zx")
	{
		//double x1=p1.x,x2=p2.x,x3=p3.x;
		//double z1=p1.z,z2=p2.z,z3=p3.z;
		return 0.25*((p1.x+p2.x)*(p1.z+p2.z)+(p1.x+p3.x)*(p1.z+p3.z)+(p2.x+p3.x)*(p2.z+p3.z))/3.0;
	} 
	else if (type=="yz" || type=="zy")
	{
		double y1=p1.y,y2=p2.y,y3=p3.y;
		double z1=p1.z,z2=p2.z,z3=p3.z;
		return 0.25*((y1+y2)*(z1+z2)+(y1+y3)*(z1+z3)+(y2+y3)*(z2+z3))/3.0;
	} 
	else if (type=="zz")
	{
		double z1=p1.z,z2=p2.z,z3=p3.z;
		return 0.25*((z1+z2)*(z1+z2)+(z1+z3)*(z1+z3)+(z2+z3)*(z2+z3))/3.0;
	}
	//Integral fucntion for calculating the anaytical dervatives
	else if (type=="x")
	{
		return (p1.x+p2.x+p3.x)/3.0;
	}
	else if (type=="y")
	{
		return (p1.y+p2.y+p3.y)/3.0;
	}
	else if (type=="xy" || type=="yx")
	{
		double x1=p1.x,x2=p2.x,x3=p3.x;
		double y1=p1.y,y2=p2.y,y3=p3.y;
		return 0.25*((x1+x2)*(y1+y2)+(x1+x3)*(y1+y3)+(x2+x3)*(y2+y3))/3.0;
	}
	else if (type=="yy")
	{
		double y1=p1.y,y2=p2.y,y3=p3.y;
		return 0.25*((y1+y2)*(y1+y2)+(y1+y3)*(y1+y3)+(y2+y3)*(y2+y3))/3.0;
	}
	else if (type=="xx")
	{
		double x1=p1.x,x2=p2.x,x3=p3.x;
		return 0.25*((x1+x2)*(x1+x2)+(x1+x3)*(x1+x3)+(x2+x3)*(x2+x3))/3.0;
	}
	else
	{
		cout<<"Error: Invalid integral type in calling gauss()"<<endl;
		system("Pause");
		exit(0);
	}
}
void area(vec3D &p1,vec3D &p2,vec3D &p3,double &s,double &nz)
{
	double nx=(p2.y-p1.y)*(p3.z-p1.z)-(p2.z-p1.z)*(p3.y-p1.y);
	double ny=(p2.z-p1.z)*(p3.x-p1.x)-(p2.x-p1.x)*(p3.z-p1.z);
	nz=(p2.x-p1.x)*(p3.y-p1.y)-(p2.y-p1.y)*(p3.x-p1.x);
	s=sqrt(nx*nx+ny*ny+nz*nz);
	if (s<1e-20)
	{
		s=nz=0;
		return;
	}
	nz=nz/s;
	s=0.5*s;
	return;
}

void area( vec3D &p1,vec3D &p2,vec3D &p3,double &s,vec3D &n )
{
	double nx=(p2.y-p1.y)*(p3.z-p1.z)-(p2.z-p1.z)*(p3.y-p1.y);
	double ny=(p2.z-p1.z)*(p3.x-p1.x)-(p2.x-p1.x)*(p3.z-p1.z);
	double nz=(p2.x-p1.x)*(p3.y-p1.y)-(p2.y-p1.y)*(p3.x-p1.x);
	s=sqrt(nx*nx+ny*ny+nz*nz);
	if (s<1e-20)
	{
		s=nz=nx=ny=0;
		return;
	}
	nx=nx/s;ny=ny/s;nz=nz/s;
	n.x=nx;n.y=ny;n.z=nz;
	s=0.5*s;
	return;
}
double random( double a,double b )
{
	return a+(b-a)*rand()/(RAND_MAX+1.0);
}

void inverse( double a[][3],double& det,bool& degenerate )
{
	double inva[3][3];
	degenerate=false;
	det=a[0][0]*a[1][1]*a[2][2]-a[0][0]*a[1][2]*a[2][1]-a[0][1]*a[1][0]*a[2][2]+a[0][1]*a[1][2]*a[2][0]+a[0][2]*a[1][0]*a[2][1]-a[0][2]*a[1][1]*a[2][0];
	if (abs(det)<1e-20)
	{
		degenerate=true;
		return;
	}
	inva[0][0]=a[1][1]*a[2][2]-a[1][2]*a[2][1];inva[0][1]=a[0][2]*a[2][1]-a[0][1]*a[2][2];inva[0][2]=a[0][1]*a[1][2]-a[0][2]*a[1][1];
	inva[1][0]=a[1][2]*a[2][0]-a[1][0]*a[2][2];inva[1][1]=a[0][0]*a[2][2]-a[0][2]*a[2][0];inva[1][2]=a[0][2]*a[1][0]-a[0][0]*a[1][2];
	inva[2][0]=a[1][0]*a[2][1]-a[1][1]*a[2][0];inva[2][1]=a[0][1]*a[2][0]-a[0][0]*a[2][1];inva[2][2]=a[0][0]*a[1][1]-a[0][1]*a[1][0];
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			a[i][j]=inva[i][j]/det;
		}
	}
	return;
}

double tetredron_volume( vec3D& v1,vec3D& v2,vec3D& v3,vec3D& v4 )
{
	double x1=v1.x,x2=v2.x,x3=v3.x,x4=v4.x;
	double y1=v1.y,y2=v2.y,y3=v3.y,y4=v4.y;
	double z1=v1.z,z2=v2.z,z3=v3.z,z4=v4.z;
	return abs(x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4
		-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1
		-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2)/6.0;
}

void absmax_component(vec3D& n,int&id,double& absmax)
{
	if (abs(n.x)>=abs(n.y) && abs(n.x)>=abs(n.z))
	{
		id=1,absmax=n.x;
	}
	else if (abs(n.y)>=abs(n.x) && abs(n.y)>=abs(n.z))
	{
		id=2,absmax=n.y;
	}
	else
	{
		id=3,absmax=n.z;
	}
	return;
}
vec3D interpolate(vec3D_double &p1,vec3D_double &p2)
{
	double ratio=p1._variable/(p1._variable-p2._variable);
	return p1._coor+(p2._coor-p1._coor)*ratio;
}
int sign(double a)
{
	if (a>0)	{
		return 1;
	}
	else if (a<0)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}
bool comparison_cell_coordinate(vec3D const &coor_min_old,vec3D const &coor_max_old,vec3D const &coor_min_new,vec3D const &coor_max_new)
{
	//if (coor_min_old.x>coor_max_new.x) return true;
	//if (coor_max_old.x<coor_min_new.x) return true;
	//if (coor_min_old.y>coor_max_new.y) return true;
	//if (coor_max_old.y<coor_min_new.y) return true;
	//if (coor_min_old.z>coor_max_new.z) return true;
	//if (coor_max_old.z<coor_min_new.z) return true;
	//return false;
	return (coor_min_old.x>coor_max_new.x || coor_max_old.x<coor_min_new.x || 
		    coor_min_old.y>coor_max_new.y || coor_max_old.y<coor_min_new.y || 
			coor_min_old.z>coor_max_new.z || coor_max_old.z<coor_min_new.z);
}
void bisect(double& a,double& b,double& f)
{
	if (f<0)
	{
		a=(a+b)/2.0;
	} 
	else
	{
		b=(a+b)/2.0;
	}
	return;
}
void secant(double& a,double& b,double& fa,double& fb,double& f)
{
	if (f<0)
	{
		a=a-fa*(b-a)/(fb-fa);
		fa=f;
	}
	else
	{
		b=a-fa*(b-a)/(fb-fa);
		fb=f;
	}
	return;
}

double newton( double& x,double& f,double& df )
{
	return x-f/df;
}
double quad(double a1,double a2,double f1,double df1,double f2)
{
	double temp=1+(f1-f2)/((a2-a1)*df1);
	return a1+0.5*(a2-a1)/temp;
}
double cubic1(double a1,double a2,double f1,double df1,double f2,double df2)
{
	//double s,z,w,temp;
	//s=3*(f2-f1)/(a2-a1);
	//z=s-df1-df2;
	//w=sqrt(z*z-df1*df2);
	//temp=1-(df2+w+z)/(df2-df1+2*w);
	//return a1+(a2-a1)*temp;
	double temp=2*df1+df2-3*(f2-f1)/(a2-a1);
	temp=sqrt(abs((temp-df1)*(temp-df1)-df1*df2))-temp;
	return a1-df1*(a2-a1)/temp;
}
double cubic2(double a1,double a,double f1,double df1,double f,double df)
{
	double temp=2*df+df1-3*(f1-f)/(a1-a);
	if ((temp-df)*(temp-df)-df1*df<0)
	{
		//The cubic interpolation does not have a local minimal point, the step should be amplified
		return 2*a;
	}
	else
	{
		//The local minimal point is the next iterative point
		temp=sqrt(((temp-df)*(temp-df)-df1*df))+temp;
	}
	return a-df*(a-a1)/temp;
}
bool is_exist(int i,vector<int> &list)
{
	for (int j=0;j<list.size();j++)
	{
		if (i==list[j])
		{
			return true;
		}
	}
	return false;
}
int find_int(double number)
{
	for (int i = 0; i < 10000000; i++)
	{
		if (abs(number-i)<1e-10)
		{
			return i;
		}
	}
	return -1;
}
double root_of_approximate_volume_function(vec3D &min_info,vec3D &max_info,double &reference_fraction)
{
	double d_min=min_info.x,f_min=min_info.y,df_min=min_info.z;
	double d_max=max_info.x,f_max=max_info.y,df_max=max_info.z;
	//Calculate the coefficients of the cubic function
	double d1=d_min-d_max,d2=d1*d1,d3=d2*d1;
	double a=(df_min+df_max)-(2*f_min-2*f_max)/d1;
	double b=-(d_min*d_min*(df_min+2*df_max)-d_max*d_max*(2*df_min+df_max)-d_min*(3*f_min-3*f_max)+d_max*(3*f_max-3*f_min+d_min*(df_min-df_max)))/d1;
	double c=df_max*d2+(2*d_max*(df_min+2*df_max))*d1-(6*d_max*d_max*(f_min-f_max))/d1 + (3*d_max*(2*f_max-2*f_min+d_max*df_min+d_max*df_max));
	double d=(d_max*d_max*f_min*(3*d_min-d_max))/d1 - (d_min*d_max*d_max*df_min) - (d_min*d_min*d_max*df_max) + (d_min*d_min*f_max*(d_min-3*d_max))/d1;
	double iterator=0.5*(d_min+d_max);
	double iterator2,iterator3;
	double f,df;
	iterator2=iterator*iterator;iterator3=iterator2*iterator;
	f=a*iterator3+b*iterator2+c*iterator+d;
	df=3*a*iterator2+2*b*iterator+c;
	int num_Newton=0;
	while (abs(f)>1e-14 && num_Newton<20)
	{
		iterator=newton(iterator,f,df);
		if (iterator<d_min || iterator>d_max)
		{
			iterator=0.5*(d_max+d_min);
		}
		num_Newton=num_Newton+1;
		iterator2=iterator*iterator;iterator3=iterator2*iterator;
		f=a*iterator3+b*iterator2+c*iterator+d;
		df=3*a*iterator2+2*b*iterator+c;
		if (abs(f)<1e-14)
		{
			return iterator;
		}
		if (f<0)
		{
			d_min=iterator;
		}
		else if (f>0)
		{
			d_max=iterator;
		}
	}
	//while (abs(f)>1e-14)
	//{
	//	iterator=0.5*(d_min+d_max);
	//	iterator2=iterator*iterator;iterator3=iterator2*iterator;
	//	f=a*iterator3+b*iterator2+c*iterator+d;
	//	if (abs(f)<1e-14)
	//	{
	//		return iterator;
	//	}
	//	if (f<0)
	//	{
	//		d_min=iterator;
	//	}
	//	else if (f>0)
	//	{
	//		d_max=iterator;
	//	}
	//}
	return iterator;
}
vec3D convert_to_sphere_coor(vec3D coor)
{
	double r,theta,phy;
	r=coor.get_length();
	phy=asin(coor.z/r);
	double rxy=sqrt(coor.x*coor.x+coor.y*coor.y);
	theta=asin(coor.y/rxy);
	vec3D temp(r,theta,phy);
	return temp;
}
double calcualte_alpha(double d0,int N,double r_in,double r_out)
{
	double temp=(r_out-r_in)/d0;
	double alpha_min=1,alpha_max=N;
	double trial=(alpha_min+alpha_max)*0.5;
	while (alpha_max-alpha_min>1e-15)
	{
		double f=f_calculate_alpha(N,temp,trial);
		if (f>0)
		{
			alpha_max=trial;
		}
		else
		{
			alpha_min=trial;
		}
		trial=0.5*(alpha_min+alpha_max);
	}
	return (alpha_min+alpha_max)*0.5;
}
double f_calculate_alpha(int N,double k,double x)
{
	return pow(x,N)-k*x+k-1;
}