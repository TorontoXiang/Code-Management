#include "Tpolygon.h"
#include <algorithm>
#include <functional>
void Tpolygon::Add_vertex(vec2D& new_vertex,int id)
{
	_vertex_list[id]=new_vertex;
	return;
}
void Tpolygon::calculate_offset_and_span(vec2D& n)
{
	//vec2D n(cos(theta),sin(theta));
	//Calculate the offset
	for (int i = 0; i < _num_vertex; i++)
	{
		_offset_list[i]._offset=_vertex_list[i]*n;
		_offset_list[i]._id=i;
	}
	//Calculate the span
	for (int i = 0; i < _num_vertex; i++)
	{
		_offset_list[i]._span=calculate_span(i);
	}
	return;
}
double Tpolygon::calculate_span(int id)
{
	double offset0=_offset_list[id]._offset;
	double d1,d2;
	for (int i = 1; i < _num_vertex-1; i++)
	{
		d1=_offset_list[(id+i)%_num_vertex]._offset-offset0;
		d2=_offset_list[((id+i)%_num_vertex+1)%_num_vertex]._offset-offset0;
		if (d1*d2<=0)
		{
			vec2D p1=_vertex_list[(id+i)%_num_vertex];
			vec2D p2=_vertex_list[((id+i)%_num_vertex+1)%_num_vertex];
			vec2D temp=p1+(p2-p1)*(d1/(d1-d2));
			vec2D p0=_vertex_list[id];
			return sqrt((p0-temp)*(p0-temp));
		}
	}
	return 0;
}
void Tpolygon::flood_algorithm(double fraction,double& offset_solution,double& length)
{
	double eps=1e-15;

	//The flood algorithm
	//First step: arange the offset list in increasing order
	std::sort(_offset_list.begin(),_offset_list.end(),greater<offset>());

	//Second step: solve the volume equation
	double f1=0,f2;
	double doffset;
	for (int i = 1; i < _num_vertex; i++)
	{
		doffset=_offset_list[i]._offset-_offset_list[i-1]._offset;
		f2=f1+(0.5*doffset*(_offset_list[i]._span+_offset_list[i-1]._span))/_area;
		if (f2<fraction)
		{
			//The offset of the solution is higher than _offset_list[i]._offset
			f1=f2;
		}
		else
		{
			//The solution is between _offset_list[i-1]._offset and _offset_list[i]._offset
			double dspan=_offset_list[i]._span-_offset_list[i-1]._span;
			double offset0=_offset_list[i-1]._offset;
			if (abs(dspan)<eps)
			{
				length=(_offset_list[i]._span+_offset_list[i-1]._span)/2.0;
				offset_solution=offset0+(fraction-f1)*doffset/(f2-f1);
			}
			else
			{
				double s1=_offset_list[i-1]._span;
				double s2=_offset_list[i]._span;
				length=s1*s1+(fraction-f1)*(s2*s2-s1*s1)/(f2-f1);
				length=sqrt(length);
				offset_solution=offset0+(length-_offset_list[i-1]._span)*doffset/dspan;
			}
			break;
		}
	}
}
void Tpolygon::divide_polygon(vec2D& n,double fraction,Tpolygon* &below,Tpolygon* &above)
{
	double offset_solution,length;
	calculate_offset_and_span(n);
	flood_algorithm(fraction,offset_solution,length);
	int polygon_below[20],polygon_above[20];
	vec2D endpoint[2];
	int counter_below=0,counter_above=0;
	int counter_endpoint=0;
	double d1,d2;
	int begin=_offset_list[0]._id;
	for (int i = 0; i < _num_vertex; i++)
	{
		d1=_vertex_list[(begin+i)%_num_vertex]*n-offset_solution;
		d2=_vertex_list[((begin+i)%_num_vertex+1)%_num_vertex]*n-offset_solution;
		if (d1<=0)
		{
			polygon_below[counter_below]=(begin+i)%_num_vertex;
			counter_below=counter_below+1;
		}
		if (d1>=0)
		{
			polygon_above[counter_above]=(begin+i)%_num_vertex;
			counter_above=counter_above+1;
		}
		if (d1*d2<0)
		{
			vec2D p1=_vertex_list[(begin+i)%_num_vertex];
			vec2D p2=_vertex_list[((begin+i)%_num_vertex+1)%_num_vertex];
			vec2D temp=p1+(p2-p1)*(d1/(d1-d2));
			vec2D dp=p2-p1;
			endpoint[counter_endpoint]=temp;
			polygon_below[counter_below]=-counter_endpoint-1;
			polygon_above[counter_above]=-counter_endpoint-1;
			counter_endpoint=counter_endpoint+1;	
			counter_below=counter_below+1;
			counter_above=counter_above+1;
		}
	}
	vec2D vertex;
	below=new Tpolygon(counter_below);
	for (int i = 0; i < counter_below; i++)
	{
		if (polygon_below[i]>=0)
		{
			vertex=_vertex_list[polygon_below[i]];
		}
		else
		{
			vertex=endpoint[-polygon_below[i]-1];
		}
		below->Add_vertex(vertex,i);
	}
	above=new Tpolygon(counter_above);
	for (int i = 0; i < counter_above; i++)
	{
		if (polygon_above[i]>=0)
		{
			vertex=_vertex_list[polygon_above[i]];
		}
		else
		{
			vertex=endpoint[-polygon_above[i]-1];
		}
		above->Add_vertex(vertex,i);
	}
	return;
}
void Tpolygon::plot_polygon(ofstream& output,double value)
{
	output<<"TITLE=\"material polygon\""<<endl;
	output<<"Variables=\"X\",\"Y\",\"V\""<<endl;
	output<<"ZONE F=FEPOINT N= "<<_num_vertex+1<<" E= "<<_num_vertex<<" ET=TRIANGLE"<<endl;
	for (int i = 0; i < _num_vertex; i++)
	{
		output<<_vertex_list[i].x<<" "<<_vertex_list[i].y<<" "<<value<<endl;
	}
	output<<_centroid.x<<" "<<_centroid.y<<" "<<value<<endl;
	for (int i = 0; i < _num_vertex; i++)
	{
		output<<_num_vertex+1<<" "<<i+1<<" "<<(i+1)%_num_vertex+1<<endl;
	}
	return;
}
void Tpolygon::calculate_properities()
{
	_area=0;_centroid.value(0,0);
	double s;
	vec2D centroid_triangle;
	for (int i = 1; i < _num_vertex-1; i++)
	{
		s=area(_vertex_list[0],_vertex_list[i],_vertex_list[i+1]);
		centroid_triangle=(_vertex_list[0]+_vertex_list[i]+_vertex_list[i+1])/3.0;
		_area=_area+s;
		_centroid=_centroid+centroid_triangle*s;
	}
	_centroid=_centroid/_area;
	return;
}
void Tpolygon::calculate_geometry_for_derivatives(vec2D& n,double offset_solution,double& area_below,vec2D& centroid_below)
{
	//vec2D n(cos(theta),sin(theta));
	//calculate the endpoints of the surface and the corresponding inclinations
	int polygon_below[20];
	vec2D endpoint[2];
	int counter1=0;
	int counter2=0;
	double d1,d2;
	int begin=_offset_list[0]._id;
	for (int i = 0; i < _num_vertex; i++)
	{
		d1=_vertex_list[(begin+i)%_num_vertex]*n-offset_solution;
		d2=_vertex_list[((begin+i)%_num_vertex+1)%_num_vertex]*n-offset_solution;
		if (d1<=0)
		{
			polygon_below[counter1]=(begin+i)%_num_vertex;
			counter1=counter1+1;
		}
		if (d1*d2<0)
		{
			vec2D p1=_vertex_list[(begin+i)%_num_vertex];
			vec2D p2=_vertex_list[((begin+i)%_num_vertex+1)%_num_vertex];
			vec2D temp=p1+(p2-p1)*(d1/(d1-d2));
			vec2D dp=p2-p1;
			endpoint[counter2]=temp;
			polygon_below[counter1]=-counter2-1;
			counter2=counter2+1;	
			counter1=counter1+1;
		}
	}
	//calculate the centroid below the surface
	area_below=0;
	centroid_below.value(0,0);
	double s;
	vec2D p1=_vertex_list[begin],p2,p3,centroid_triangle;
	for (int i = 1; i < counter1-1; i++)
	{
		if (polygon_below[i]>=0)
		{
			p2=_vertex_list[polygon_below[i]];
		}
		else
		{
			p2=endpoint[-polygon_below[i]-1];
		}
		if (polygon_below[i+1]>=0)
		{
			p3=_vertex_list[polygon_below[i+1]];
		}
		else
		{
			p3=endpoint[-polygon_below[i+1]-1];
		}
		s=area(p1,p2,p3);
		centroid_triangle=(p1+p2+p3)/3.0;
		area_below=area_below+s;
		centroid_below=centroid_below+centroid_triangle*s;
	}
	centroid_below=centroid_below/area_below;
	return;
}
//void Tpolygon::rotation_algorithm(double fraction,vec2D& centroid_ref,int id,discontinuous_point (&dp)[2])
//{
//	double area_below=fraction*_area;
//	vec2D p0,p1,p2;
//	//The discontinuity theta for anti-clockwise order
//	dp[0]._id=id;
//	p0=_vertex_list[id];
//	double lack=area_below;
//	vec2D centroid(0,0);
//	double s,ratio;
//	for (int i = 1; i < _num_vertex-1; i++)
//	{
//		p1=_vertex_list[(id+i)%_num_vertex];
//		p2=_vertex_list[(id+i+1)%_num_vertex];
//		s=area(p0,p1,p2);
//		if (s<lack)
//		{
//			lack=lack-s;
//			centroid.x+=(p0.x+p1.x+p2.x)*s/3.0;
//			centroid.y+=(p0.y+p1.y+p2.y)*s/3.0;
//		}
//		else
//		{
//			ratio=lack/s;
//			vec2D pp(p1.x+ratio*(p2.x-p1.x),p1.y+ratio*(p2.y-p1.y));
//			dp[0]._pm=pp;
//			centroid.x+=(p0.x+p1.x+pp.x)*lack/3.0;
//			centroid.y+=(p0.y+p1.y+pp.y)*lack/3.0;
//			centroid=centroid/area_below;
//			dp[0]._centroid_below=centroid;
//			dp[0]._length=ratio=sqrt((pp-p0)*(pp-p0));
//			vec2D d(-(pp-p0).y,(pp-p0).x);
//			d.normalize();
//			double angle=acos(d.x);
//			if (d.y>=0)
//			{
//				dp[0]._theta=angle;
//			}
//			else
//			{
//				dp[0]._theta=2*pi-angle;
//			}
//			ratio=ratio*ratio*ratio/(12*area_below);
//			vec2D x1(d.y*ratio,-d.x*ratio);
//			dp[0]._d1=(centroid-centroid_ref)*x1*2;
//			break;
//		}
//	}
//
//	//The discontinuity theta for clockwise order
//	dp[1]._id=id;
//	p0=_vertex_list[id];
//	lack=area_below;
//	centroid.value(0,0);
//	for (int i = 1; i < _num_vertex-1; i++)
//	{
//		p1=_vertex_list[(id-i+_num_vertex)%_num_vertex];
//		p2=_vertex_list[(id-i-1+_num_vertex)%_num_vertex];
//		s=area(p0,p1,p2);
//		if (s<lack)
//		{
//			lack=lack-s;
//			centroid.x+=(p0.x+p1.x+p2.x)*s/3.0;
//			centroid.y+=(p0.y+p1.y+p2.y)*s/3.0;
//		}
//		else
//		{
//			double ratio=lack/s;
//			vec2D pp(p1.x+ratio*(p2.x-p1.x),p1.y+ratio*(p2.y-p1.y));
//			dp[1]._pm=pp;
//			centroid.y+=(p0.y+p1.y+pp.y)*lack/3.0;
//			centroid.x+=(p0.x+p1.x+pp.x)*lack/3.0;
//			centroid=centroid/area_below;
//			dp[1]._centroid_below=centroid;
//			dp[1]._length=ratio=sqrt((pp-p0)*(pp-p0));
//			vec2D d((pp-p0).y,-(pp-p0).x);
//			d.normalize();
//			double angle=acos(d.x);
//			if (d.y>=0)
//			{
//				dp[1]._theta=angle;
//			}
//			else
//			{
//				dp[1]._theta=2*pi-angle;
//			}
//			ratio=ratio*ratio*ratio/(12*area_below);
//			vec2D x1(d.y*ratio,-d.x*ratio);
//			dp[1]._d1=(centroid-centroid_ref)*x1*2;
//			break;
//		}
//	}
//	return;
//}
int Tpolygon::trial_id_base(vec2D& n,double fraction,int& order)
{
	calculate_offset_and_span(n);
	std::sort(_offset_list.begin(),_offset_list.end(),greater<offset>());

	//Second step: solve the volume equation
	double f1=0,f2;
	double doffset;
	int id0,id1;
	for (int i = 1; i < _num_vertex; i++)
	{
		doffset=_offset_list[i]._offset-_offset_list[i-1]._offset;
		f2=f1+(0.5*doffset*(_offset_list[i]._span+_offset_list[i-1]._span))/_area;
		if (f2<fraction)
		{
			//The offset of the solution is higher than _offset_list[i]._offset
			f1=f2;
		}
		else
		{
			id0=_offset_list[i]._id;
			id1=(id0+1)%_num_vertex;
			vec2D d=_vertex_list[id1]-_vertex_list[id0];
			if (d*n<0)
			{
				order=0;
				return id0;
			}
			else
			{
				order=1;
				return (id0-1+_num_vertex)%_num_vertex;
			}
		}
	}
	return 0;
}
void Tpolygon::calculate_range(vec2D& x_range,vec2D& y_range)
{
	x_range.value(1000000,-1000000);y_range=x_range;
	for (int i = 0; i < _num_vertex; i++)
	{
		if (_vertex_list[i].x<x_range.x)
		{
			x_range.x=_vertex_list[i].x;
		}
		if (_vertex_list[i].x>x_range.y)
		{
			x_range.y=_vertex_list[i].x;
		}
		if (_vertex_list[i].y<y_range.x)
		{
			y_range.x=_vertex_list[i].y;
		}
		if (_vertex_list[i].y>y_range.y)
		{
			y_range.y=_vertex_list[i].y;
		}
	}
	return;
}