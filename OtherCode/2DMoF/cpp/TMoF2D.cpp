#include "TMoF2D.h"
#include <cmath>
TMoF2D::TMoF2D(int highest_order,int num_vertex):_hightest_order(highest_order),_num_vertex(num_vertex)
{
	for (int i = 0; i < 2; i++)
	{
		_inclination[i]=new double[_num_vertex];
	}
}
void TMoF2D::Access_to_polygon(Tpolygon* polygon_ptr)
{
	_polygon_ptr=polygon_ptr;
	_num_vertex=_polygon_ptr->_num_vertex;
	//Calculate the incliantion of each edge
	vec2D dp;
	double theta;
	for (int i = 0; i < _num_vertex; i++)
	{
		dp=_polygon_ptr->_vertex_list[(i+1)%_num_vertex]-_polygon_ptr->_vertex_list[i];
		dp.normalize();
		if (dp.y>0)
		{
			_inclination[0][i]=dp.y;_inclination[1][i]=dp.x;
		}
		else
		{
			_inclination[0][i]=-dp.y;_inclination[1][i]=-dp.x;
		}
	}
}
void TMoF2D::Calculate_interation_information(double theta)
{
	double offset_solution,length;
	_polygon_ptr->calculate_properities();
	vec2D n(cos(theta),sin(theta));
	_iteration_point.theta=theta;
	_polygon_ptr->calculate_offset_and_span(n);
	_polygon_ptr->flood_algorithm(_fraction,offset_solution,length);
	double area_below;
	vec2D centroid_below;
	_polygon_ptr->calculate_geometry_for_derivatives(n,offset_solution,area_below,centroid_below);
	_iteration_point.value=(_centroid_ref-centroid_below)*(_centroid_ref-centroid_below);
	//Calculate the derivatives
	Calculate_derivatives(n,length,centroid_below);

	return;
}
double TMoF2D::G_derivative(int i)
{
	if (i==1)
	{
		return _iteration_point.d1;
	}
	else if (i==2)
	{
		return _iteration_point.d2;
	}
	else if (i==3)
	{
		return _iteration_point.d3;
	}
	return 0;
}
void TMoF2D::Calculate_derivatives(vec2D& n,double length,vec2D& centroid_below)
{
	double temp=_fraction*_polygon_ptr->G_area();
	double s=n.y,c=n.x;
	temp=length*length*length/(12*temp);
	double s1,s2,c1,c2;
	double a,b;
	double dx=centroid_below.x-_centroid_ref.x;
	double dy=centroid_below.y-_centroid_ref.y;
	if (_hightest_order>=1)
	{
		_iteration_point.d1=(dx*s*temp-dy*c*temp)*2;
	}
	if (_hightest_order>=2)
	{
		s1=_s1*c-_c1*s;s2=_s2*c-_c2*s;
		c1=_c1*c+_s1*s;c2=_c2*c+_s2*s;
		a=s1/c1+s2/c2;
		_iteration_point.d2=2*temp*temp+(dx*(-3*a*s+2*c)+dy*(3*a*c+2*s))*temp;
	}
	if (_hightest_order>=3)
	{
		b=-(1/(c1*c1)+1/(c2*c2));
		_iteration_point.d3=-3*temp*temp*(-s*(-3*a*s+2*c)+c*(3*a*c+2*s))+(dx*(-1.5*a*(-3*a*s+2*c)-3*(b*s+a*c)-2*s)+dy*(-1.5*a*(3*a*c+2*s)+3*(b*c-a*s)+2*c))*temp;
	}
	return;
}
double TMoF2D::initial_guess()
{
	vec2D centroid_polygon=_polygon_ptr->G_centroid();
	vec2D d=centroid_polygon-_centroid_ref;
	double theta_initial=acos(d.x/sqrt(d*d));
	if (d.y<0)
	{
		theta_initial=2*pi-theta_initial;
	}
	return theta_initial;
}
void TMoF2D::Calculate_single_summit_area(double& theta_min,double& theta_max,double& theta_initial)
{
	int n_divided=2*_num_vertex;double dtheta=2*pi/n_divided;
	double theta0=theta_initial;
	int direction;
	Calculate_interation_information(theta0);
	if (abs(_iteration_point.d1)<1e-10)
	{
		theta_min=theta_max=theta_initial;
		return;
	}
	else if (_iteration_point.d1>0)
	{
		direction=-1;
		theta_max=theta_initial;
	}
	else
	{
		direction=1;
		theta_min=theta_initial;
	}
	for (int i = 0; i < n_divided; i++)
	{
		theta0=theta0+direction*dtheta;
		Calculate_interation_information(theta0);
		if (abs(_iteration_point.d1)<1e-10)
		{
			theta_min=theta_max=theta0;
			return;
		}
		if (direction==1)
		{
			if (_iteration_point.d1<0)
			{
				theta_min=theta0;
			}
			else if (_iteration_point.d1>0)
			{
				theta_max=theta0;
				return;
			}
		}
		else if (direction==-1)
		{
			if (_iteration_point.d1>0)
			{
				theta_max=theta0;
			}
			else if (_iteration_point.d1<0)
			{
				theta_min=theta0;
				return;
			}
		}
	}
	return;
}
double TMoF2D::Zoom_algorithm(double theta_min,double theta_max)
{
	double eps=1e-7;
	double theta1=theta_min,theta2=theta_max;
	iteration_point iterator1,iterator2;
	Calculate_interation_information(theta1);
	iterator1=_iteration_point;
	Calculate_interation_information(theta2);
	iterator2=_iteration_point;
	double theta,value,value1=iterator1.value,value2=iterator2.value,derivative,derivative1=iterator1.d1,derivative2=iterator2.d1;
	if (iterator1.value>iterator2.value)
	{
		theta=theta1;value=iterator1.value;derivative=iterator1.d1;
		theta1=theta2;value1=iterator2.value;derivative1=iterator2.d1;
		theta2=theta;value2=value;derivative2=derivative;
	}
	double d1,d2;
	while (abs(theta1-theta2)>eps)
	{
		d1=derivative1+derivative2-3*(value1-value2)/(theta1-theta2);
		if (theta2-theta1>0)
		{
			d2=sqrt(d1*d1-derivative1*derivative2);
		}
		else if (theta2-theta1<0)
		{
			d2=-sqrt(d1*d1-derivative1*derivative2);
		}
		else
		{
			d2=0;
		}
		theta=theta2-(theta2-theta1)*(derivative2+d2-d1)/(derivative2-derivative1+2*d2);
		if ((theta-theta1)/(theta2-theta1)<0.05)
		{
			theta=theta1+(theta2-theta1)*0.05;
		}
		Calculate_interation_information(theta);
		value=_iteration_point.value;derivative=_iteration_point.d1;
		if (abs(derivative)<1e-7)
		{
			return theta;
		}
		if (value1<value)
		{
			theta2=theta;value2=value;derivative2=derivative;
		}
		else if (derivative*(theta2-theta1)<0)
		{
			theta1=theta;value1=value;derivative1=derivative;
		}
		else
		{
			theta2=theta1;value2=value1;derivative2=derivative1;
			theta1=theta;value1=value;derivative1=derivative;
		}
	}
	return 0.5*(theta1+theta2);
}
vec2D TMoF2D::calculate_centroid_below(double theta)
{
	double offset_solution,length;
	vec2D n(cos(theta),sin(theta));
	_polygon_ptr->calculate_offset_and_span(n);
	_polygon_ptr->flood_algorithm(_fraction,offset_solution,length);
	double area_below;
	vec2D centroid_below;
	_polygon_ptr->calculate_geometry_for_derivatives(n,offset_solution,area_below,centroid_below);
	return centroid_below;
}
double TMoF2D::MoF_solver_new()
{
	//Determine the initial id_base
	vec2D n_trial=_polygon_ptr->_centroid-_centroid_ref;
	n_trial.normalize();
	int id_base,id_base0;
	int order;
	id_base0=id_base=_polygon_ptr->trial_id_base(n_trial,_fraction,order);
	int id_below;
	int time=0;
	vec2D *begin,*end,*temp_ptr;
	vec2D pm,pm0,pb_trial,pre_begin;
	vec2D n,d1,d2,d_centroid;
	double index,s1,s2,sx,l1,l2,l;
	discontinuous_point dp_pre,dp;
	//Find the initial discontinuous point
	double area_below=_fraction*_polygon_ptr->_area;
	vec2D pb,p_below,p_high;
	//The discontinuity theta for anti-clockwise order
	//dp._id=id_base;
	pb=_polygon_ptr->_vertex_list[id_base];
	double lack=area_below;
	double xc=0,yc=0;
	double s,ratio;
	for (int i = 1; i < _num_vertex-1; i++)
	{
		if (order==0)
		{
			p_below=_polygon_ptr->_vertex_list[(id_base+i)%_num_vertex];
			p_high=_polygon_ptr->_vertex_list[(id_base+i+1)%_num_vertex];
		}
		else
		{
			p_high=_polygon_ptr->_vertex_list[(id_base-i+_num_vertex)%_num_vertex];
			p_below=_polygon_ptr->_vertex_list[(id_base-i-1+_num_vertex)%_num_vertex];
		}
		s=area(pb,p_below,p_high);
		if (s<lack)
		{
			lack=lack-s;
			xc+=(pb.x+p_below.x+p_high.x)*s/3.0;
			yc+=(pb.y+p_below.y+p_high.y)*s/3.0;
		}
		else
		{
			ratio=lack/s;			
			if (order==0)
			{
				id_below=(id_base+i)%_num_vertex;
				pm.value(p_below.x+ratio*(p_high.x-p_below.x),p_below.y+ratio*(p_high.y-p_below.y));
				begin=&pb;end=&pm;
				xc+=(pb.x+p_below.x+pm.x)*lack/3.0;
				yc+=(pb.y+p_below.y+pm.y)*lack/3.0;
			}
			else
			{
				id_below=(id_base-i-1+_num_vertex)%_num_vertex;
				pm.value(p_high.x+ratio*(p_below.x-p_high.x),p_high.y+ratio*(p_below.y-p_high.y));
				begin=&pm;end=&pb;
				xc+=(pb.x+p_high.x+pm.x)*lack/3.0;
				yc+=(pb.y+p_high.y+pm.y)*lack/3.0;
			}
			xc=xc/area_below;yc=yc/area_below;
			dp._xc=xc;
			dp._yc=yc;
			dp._length=ratio=sqrt((pm-pb)*(pm-pb));
			n.value(begin->y-end->y,end->x-begin->x);
			n.normalize();
			double angle=acos(n.x);
			if (n.y>=0)
			{
				dp._theta=angle;
			}
			else
			{
				dp._theta=2*pi-angle;
			}
			ratio=ratio*ratio*ratio/(12*area_below);
			vec2D x1(n.y*ratio,-n.x*ratio);
			dp._d1=((xc-_centroid_ref.x)*x1.x+(yc-_centroid_ref.y)*x1.y)*2;
			if (dp._d1<0)
			{
				_s1=_inclination[0][id_base];_c1=_inclination[1][id_base];
			}
			else
			{
				_s2=_inclination[0][(id_base-1+_num_vertex)%_num_vertex];
				_c2=_inclination[1][(id_base-1+_num_vertex)%_num_vertex];
				_theta0=dp._theta;_centroid0.value(dp._xc,dp._yc);
				_pb1=_polygon_ptr->_vertex_list[id_below];_d1=_polygon_ptr->_vertex_list[(id_below+1)%_num_vertex]-_pb1;
				_pb2=_polygon_ptr->_vertex_list[(id_base-1+_num_vertex)%_num_vertex];_d2=_polygon_ptr->_vertex_list[id_base]-_pb2;
				_i1=pm;_i2=_polygon_ptr->_vertex_list[id_base];
			}
			break;
		}
	}
	int revise_time1=0,time_revise2=0;
	double theta_optimal,value=1000000;
	if (dp._d1>0)
	{
		order=1;
	}
	else
	{
		order=0;
	}
	while (time!=2)
	{
		dp_pre=dp;
		pre_begin=*begin;
		//Find the next discomtious point
		//Suppose the id_base+1 is the next id_base;
		if (order==0)
		{
			pb_trial=_polygon_ptr->_vertex_list[(id_base+1)%_num_vertex];
		}
		else
		{
			pb_trial=_polygon_ptr->_vertex_list[(id_base-1+_num_vertex)%_num_vertex];
		}
		d1=pm-pb_trial;d2=p_high-p_below;
		index=d2.x*d1.y-d2.y*d1.x;
		ratio=-1;
		if (abs(index)>1e-10)
		{
			ratio=(pb.x*d1.y-pb.y*d1.x+p_below.y*d1.x-p_below.x*d1.y)/index;
		}
		if (ratio>0 && ratio<1)
		{
			time_revise2=time_revise2+1;
			if (order==0)
			{
				id_base=(id_base+1)%_num_vertex;
			}
			else
			{
				id_base=(id_base-1+_num_vertex)%_num_vertex;
			}
			s1=area(pb,pb_trial,pm);
			l1=sqrt((pm-pb_trial)*(pm-pb_trial));
			pm=p_below+d2*ratio;
			s2=area(pb,pb_trial,pm);
			l2=sqrt((pm-pb)*(pm-pb));
			sx=s1*s2/(s1+s2);
			l=(l1+l2)/3.0;
			pb=pb_trial;
			n.value(begin->y-end->y,end->x-begin->x);
			dp._length=ratio=sqrt(n*n);
			n.normalize();
			double angle=acos(n.x);
			if (n.y>=0)
			{
				dp._theta=angle;
			}
			else
			{
				dp._theta=2*pi-angle;
			}
			if (order==0)
			{
				d1=*end-pre_begin;
			}
			else
			{
				d1=pre_begin-*end;
			}
			d1.normalize();
			dp._xc=xc+d1.x*(l*sx/area_below);
			dp._yc=yc+d1.y*(l*sx/area_below);
			xc=dp._xc;yc=dp._yc;
			ratio=ratio*ratio*ratio/(12*area_below);
			dp._d1=((xc-_centroid_ref.x)*n.y*ratio-(yc-_centroid_ref.y)*n.x*ratio)*2;
			//dp._id=id_base;
			pb=pb_trial;
			if (dp._d1<0)
			{
				_s1=_inclination[0][id_base];_c1=_inclination[1][id_base];
			}
			else
			{
				_s2=_inclination[0][(id_base-1+_num_vertex)%_num_vertex];
				_c2=_inclination[1][(id_base-1+_num_vertex)%_num_vertex];
				_theta0=dp._theta;_centroid0.value(dp._xc,dp._yc);
				_pb1=_polygon_ptr->_vertex_list[id_below];_d1=_polygon_ptr->_vertex_list[(id_below+1)%_num_vertex]-_pb1;
				_pb2=_polygon_ptr->_vertex_list[(id_base-1+_num_vertex)%_num_vertex];_d2=_polygon_ptr->_vertex_list[id_base]-_pb2;
				_i1=pm;_i2=_polygon_ptr->_vertex_list[id_base];
			}
		}
		else
		{
			revise_time1=revise_time1+1;
			if (order==0)
			{
				pb_trial=p_high;
			}
			else
			{
				pb_trial=p_below;
			}
			if (order==0)
			{
				p_below=pb;p_high=_polygon_ptr->_vertex_list[(id_base+1)%_num_vertex];
			}
			else
			{
				p_below=_polygon_ptr->_vertex_list[(id_base-1+_num_vertex)%_num_vertex];;p_high=pb;
			}
			d1=pb_trial-pb;d2=p_high-p_below;
			ratio=(pm.x*d1.y-pm.y*d1.x+pb.y*d1.x-pb.x*d1.y)/(d2.x*d1.y-d2.y*d1.x);
			pm0=pm;pm=pb+d2*ratio;
			s1=area(pb,pm,pb_trial);
			s2=area(pb,pm,pm0);
			sx=s1*s2/(s1+s2);
			l1=sqrt((pb-pb_trial)*(pb-pb_trial));
			l2=sqrt((pm-pm0)*(pm-pm0));
			l=(l1+l2)/3.0;
			int temp=id_base;
			if (order==0)
			{
				id_base=(id_below+1)%_num_vertex;
			}
			else
			{
				id_base=id_below;
			}
			if (order==0)
			{
				id_below=temp;
			}
			else
			{
				id_below=(temp-1+_num_vertex)%_num_vertex;
			}
			pb=pb_trial;
			temp_ptr=begin;begin=end;end=temp_ptr;
			n.value(begin->y-end->y,end->x-begin->x);
			dp._length=ratio=sqrt(n*n);
			n.normalize();
			double angle=acos(n.x);
			if (n.y>=0)
			{
				dp._theta=angle;
			}
			else
			{
				dp._theta=2*pi-angle;
			}
			if (order==0)
			{
				d1=*end-pre_begin;
			}
			else
			{
				d1=pre_begin-*end;
			}
			d1.normalize();
			dp._xc=xc+d1.x*(l*sx/area_below);
			dp._yc=yc+d1.y*(l*sx/area_below);
			xc=dp._xc;yc=dp._yc;
			ratio=ratio*ratio*ratio/(12*area_below);
			dp._d1=((xc-_centroid_ref.x)*n.y*ratio-(yc-_centroid_ref.y)*n.x*ratio)*2;
			//dp._id=id_base;
			pb=pb_trial;
			if (dp._d1<0)
			{
				_s1=_inclination[0][id_base];_c1=_inclination[1][id_base];
			}
			else
			{
				_s2=_inclination[0][(id_base-1+_num_vertex)%_num_vertex];
				_c2=_inclination[1][(id_base-1+_num_vertex)%_num_vertex];
				_theta0=dp._theta;_centroid0.value(dp._xc,dp._yc);
				_pb1=_polygon_ptr->_vertex_list[id_below];_d1=_polygon_ptr->_vertex_list[(id_below+1)%_num_vertex]-_pb1;
				_pb2=_polygon_ptr->_vertex_list[(id_base-1+_num_vertex)%_num_vertex];_d2=_polygon_ptr->_vertex_list[id_base]-_pb2;
				_i1=pm;_i2=_polygon_ptr->_vertex_list[id_base];
			}
		}
		if (id_base==id_base0)
		{
			time=time+1;
		}
		if ((order==0 && dp_pre._d1<0 && dp._d1>0) || (order==1 && dp_pre._d1>0 && dp._d1<0))
		{
			iteration_point theta_min,theta_max;
			if (order==0)
			{
				theta_min.calculate_d2(dp_pre,_s1,_s2,_c1,_c2,_centroid_ref,area_below);
				theta_max.calculate_d2(dp,_s1,_s2,_c1,_c2,_centroid_ref,area_below);
			}
			else
			{
				theta_min.calculate_d2(dp,_s1,_s2,_c1,_c2,_centroid_ref,area_below);
				theta_max.calculate_d2(dp_pre,_s1,_s2,_c1,_c2,_centroid_ref,area_below);
			}
			if (theta_max.theta<theta_min.theta)
			{
				theta_max.theta=theta_max.theta+2*pi;
			}
			double theta_trial=root_of_cubic_interpolation(theta_min,theta_max);
			double eps=1e-7;
			Calculate_iteration_information_new(theta_trial-_theta0);
			while (abs(_iteration_point.d1)>eps)
			{
				if (_iteration_point.d1<0)
				{
					theta_min=_iteration_point;
				}
				else
				{
					theta_max=_iteration_point;
				}
				theta_trial=theta_trial-(2*_iteration_point.d1*_iteration_point.d2)/(2*_iteration_point.d2*_iteration_point.d2-_iteration_point.d3*_iteration_point.d1);
				if (theta_trial<theta_min.theta || theta_trial>theta_max.theta)
				{
					theta_trial=root_of_cubic_interpolation(theta_min,theta_max);
				}
				Calculate_iteration_information_new(theta_trial-_theta0);
			}
			return theta_trial;
			if (_iteration_point.value<value)
			{
				value=_iteration_point.value;
				theta_optimal=theta_trial;
			}
		}
	}
	return 0;
}
void TMoF2D::Calculate_iteration_information_new(double dtheta)
{
	vec2D n(cos(_theta0+dtheta),sin(_theta0+dtheta));
	double nx=n.x,ny=n.y;
	double d1x=_d1.x,d1y=_d1.y,pb1x=_pb1.x,pb1y=_pb1.y;
	double d2x=_d2.x,d2y=_d2.y,pb2x=_pb2.x,pb2y=_pb2.y;
	double i1x=_i1.x,i1y=_i1.y,i2x=_i2.x,i2y=_i2.y;
	double A,B,C,D,k,b;
	double a1,a2,a3;
	double r1,r2;
	k=(d2x*nx+d2y*ny)/(d1x*nx+d1y*ny);b=(-nx*(pb1x-pb2x)-ny*(pb1y-pb2y))/(d1x*nx+d1y*ny);
	A=d1y*d2x-d1x*d2y;
	B=d1x*(i1y-pb2y)-d1y*(i1x-pb2x);
	C=d2y*(i2x-pb1x)-d2x*(i2y-pb1y);
	D=(i1x-pb2x)*(i2y-pb1y)-(i1y-pb2y)*(i2x-pb1x);
	a1=A*k;a2=(C+B*k+A*b);a3=B*b+D;
	if (abs(a1)<1e-15)
	{
		r2=-a3/a2;
	}
	else
	{
		double delta=sqrt(a2*a2-4*a1*a3);
		r2=(-a2+delta)/(2*a1);
		if (r2>1 || r2<0)
		{
			r2=(-a2-delta)/(2*a1);
		}
	}
	r1=k*r2+b;
	vec2D i1=_pb1+_d1*r1,i2=_pb2+_d2*r2;
	double length=sqrt((i1-i2)*(i1-i2));
	double s1=area(_i1,i1,_i2),s2=area(_i1,i1,i2);
	double sx=s1*s2/(s1+s2);
	double l1=sqrt((i1-_i2)*(i1-_i2)),l2=sqrt((i2-_i1)*(i2-_i1));
	double l=(l1+l2)/3.0;
	vec2D d=i2-_i1;
	vec2D dr(-d.y,d.x);
	if (dr*n>0)
	{
		d.x=-d.x;d.y=-d.y;
	}
	d.normalize();
	vec2D centroid_below=_centroid0+d*(l*sx/(_fraction*_polygon_ptr->_area));
	_iteration_point.theta=_theta0+dtheta;
	_iteration_point.value=(centroid_below-_centroid_ref)*(centroid_below-_centroid_ref);
	Calculate_derivatives(n,length,centroid_below);
	return;
}