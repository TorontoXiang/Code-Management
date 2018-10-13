#include "Ttriangle.h"

void triangle::intersecting_with_polygon(Tpolygon* polygon,double& s,vec2D& moment)
{
	s=0;moment.value(0,0);
	//Calculate the transfer matrix
	double transfer_matrix[2][2];
	double x1=p1.x,x2=p2.x,x3=p3.x,y1=p1.y,y2=p2.y,y3=p3.y;
    double det=(y3-y1)*(x2-x1)-(y2-y1)*(x3-x1);
    transfer_matrix[0][0]=(y3-y1)/det;transfer_matrix[0][1]=(x1-x3)/det;
    transfer_matrix[1][0]=(y1-y2)/det;transfer_matrix[1][1]=(x2-x1)/det;

	//Calculate the coordinate of the polygon on the new system
	//for (int i = 0; i < polygon->G_num_vertex(); i++)
	//{
	//	vec2D x0=polygon->G_vertex(i);
	//	vec2D x1;
	//	x1.x=transfer_matrix[0][0]*(x0.x-p1.x)+transfer_matrix[0][1]*(x0.y-p1.y);
	//	x1.y=transfer_matrix[1][0]*(x0.x-p1.x)+transfer_matrix[1][1]*(x0.y-p1.y);
	//	polygon->S_vertex(x1,i);
	//}

	//Calculate the area and centroid below each line
	double sl;
	vec2D momentl;
	for (int i = 0; i < polygon->G_num_vertex(); i++)
	{
		vec2D q10=polygon->G_vertex(i),q20=polygon->G_vertex((i+1)%polygon->G_num_vertex());
		vec2D q1,q2;
		q1.x=transfer_matrix[0][0]*(q10.x-p1.x)+transfer_matrix[0][1]*(q10.y-p1.y);
		q1.y=transfer_matrix[1][0]*(q10.x-p1.x)+transfer_matrix[1][1]*(q10.y-p1.y);
		q2.x=transfer_matrix[0][0]*(q20.x-p1.x)+transfer_matrix[0][1]*(q20.y-p1.y);
		q2.y=transfer_matrix[1][0]*(q20.x-p1.x)+transfer_matrix[1][1]*(q20.y-p1.y);
		triangle_segment_calculation(q1,q2,sl,momentl);
		if (q1.x<q2.x)
		{
			s=s-sl;
			moment=moment-momentl;
		}
		else
		{
			s=s+sl;
			moment=moment+momentl;
		}
	}
	vec2D centroid0;
	if (s>1e-20)
	{
		centroid0=moment/s;
	}
	s=s*det;
	vec2D centroid;
	centroid.x=centroid0.x*(x2-x1)+centroid0.y*(x3-x1)+x1;
	centroid.y=centroid0.x*(y2-y1)+centroid0.y*(y3-y1)+y1;
	moment=centroid*s;
	return;
}
void triangle::triangle_segment_calculation(vec2D q1,vec2D q2,double& area,vec2D& moment)
{
	double p1[2],p2[2];
	p1[0]=p2[1]=0;p1[1]=p2[0]=1;
	area=0;moment.value(0,0);
	double xl,yl,xr,yr;
	if (q1.x<=q2.x)
	{
		xl=q1.x;yl=q1.y;xr=q2.x;yr=q2.y;
	}
	else
	{
		xl=q2.x;yl=q2.y;xr=q1.x;yr=q1.y;
	}
	if (xr==xl || xr<0 ||xl>1)
	{
		return;
	}
	if (xl<=0)
	{
		yl=yr-xr*(yr-yl)/(xr-xl);xl=0;
	}
	if (xr>=1)
	{
		yr=yr+(1-xr)*(yr-yl)/(xr-xl);xr=1;
	}
	if (yl<=0 && yr<=0)
	{
		return;
	}
	if (yl<=0 && yr>=0)
	{
		xl=xr-yr*(xr-xl)/(yr-yl);yl=0;
	}
	if (yl>=0 && yr<=0)
	{
		xr=xr-yr*(xr-xl)/(yr-yl);yr=0;
	}
	vec2D ql(xl,yl),qr(xr,yr);
	if (xl+yl>=1 && xr+yr>=1)
	{
		area=(2-xl-xr)*(xr-xl)/2.0;
		moment=calculate_moment(xl,1-xl,xr,1-xr);
	}
	else if (xl+yl<=1 && xr+yr<=1)
	{
		area=(yl+yr)*(xr-xl)/2.0;
		moment=calculate_moment(xl,yl,xr,yr);
	}
	else if (xl+yl>1 && xr+yr<=1)
	{
		double dl=xl+yl-1,dr=xr+yr-1;
		vec2D qm=ql+(qr-ql)*(dl/(dl-dr));
		area=(1-xl+qm.y)*(qm.x-xl)/2.0+(qm.y+yr)*(xr-qm.x)/2.0;
		moment=calculate_moment(xl,1-xl,qm.x,qm.y)+calculate_moment(qm.x,qm.y,xr,yr);
	}
	else if (xl+yl<=1 && xr+yr>1)
	{
		double dl=xl+yl-1,dr=xr+yr-1;
		vec2D qm=ql+(qr-ql)*(dl/(dl-dr));
		area=(yl+qm.y)*(qm.x-xl)/2.0+(qm.y+1-xr)*(xr-qm.x)/2.0;
		moment=calculate_moment(xl,yl,qm.x,qm.y)+calculate_moment(qm.x,qm.y,xr,1-xr);
	}
	return;
}
