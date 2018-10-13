#include "Ttriangle.h"
#include <cstddef>
#include <cmath>
Ttriangle::Ttriangle(vec3D v1,vec3D v2,vec3D v3):_v1(v1),_v2(v2),_v3(v3)
{
	//Create the initial ring
	_ring[1]._v=_v1;_ring[2]._v=_v2;_ring[3]._v=_v3;
	_ring[1]._discard=_ring[2]._discard=_ring[3]._discard=0;
	_ring[1]._next=&_ring[2];_ring[2]._next=&_ring[3];_ring[3]._next=&_ring[1];
	//The end node should point to the start node
	_head=&_ring[1];
	_used=3;
}
void Ttriangle::column_integration(double& volume,double& ix,double& iy,double& iz)
{
 	double s,nz;
	area(_v1,_v2,_v3,s,nz);// nz is used to determine the sign of the integration of this triangle piece
	for (int i=1;i<=4;i++)
	{
		value_ring(i);
		get_head();
		if (_head==NULL)
		{
			//If every node in a ring is out of any cutting plane, the integration of this piece equals zero
			volume=ix=iy=iz=0;
			return;
		}
		insert_new_node();
		cut_outside_node();//The first 4 plane will cut the outward node while the last plane only project the outward node
	}
	value_ring(5);
	insert_new_node();
	project_onto_plane();//Projecting the outward nodes rather than cut them
	final_integration(volume,ix,iy,iz);
	if (nz<0)
	{
		volume=-volume;ix=-ix;iy=-iy;iz=-iz;
		//If nz<0, the integrations of this piece should be negative
	} 
}
void Ttriangle::value_ring(int i)
{
	if (_head==NULL)
	{
		return;
	}
	Snode* ptr=_head;
	switch(i)
	{
	case 1:
		do 
		{
			ptr->_variable=ptr->_v.x;
			ptr=ptr->_next;
		}while(ptr!=_head);
		break;
	case 2:
		do
		{
			ptr->_variable=ptr->_v.y;
			ptr=ptr->_next;
		}while (ptr!=_head);
		break;
	case 3:
		do
		{
			ptr->_variable=ptr->_v.z;
			ptr=ptr->_next;
		}while (ptr!=_head);
		break;
	case 4:
		do
		{
			ptr->_variable=1-ptr->_v.x-ptr->_v.y;
			ptr=ptr->_next;
		}while (ptr!=_head);
		break;
	case 5:
		do
		{
			ptr->_variable=1-ptr->_v.x-ptr->_v.y-ptr->_v.z;
			ptr=ptr->_next;
		}while (ptr!=_head);
		break;
	}
	
}
void Ttriangle::get_head()
{
	Snode* ptr=_head;
	do
	{
		if (ptr->_variable>0)
		{
			_head=ptr;
			return;
		}
		ptr=ptr->_next;
	}while (ptr!=_head);
	//If all the node is out of the plane, head is none and the integration goes to zero
	_head=NULL;
}
void Ttriangle::insert_new_node()
{
	if (_head==NULL)
	{
		return;
	}
	Snode* ptr=_head;Snode* temp=NULL;double ratio;
	do
	{
		if (ptr->_variable<0)
		{
			ptr->_discard=1;
			//Discard the nodes which are outside of the cutting plane
		}
		if ((ptr->_variable)*(ptr->_next->_variable)<0)
		{
			//If the variable sign of any adjacent nodes is different, insert a new node lies on this 
			//cutting plane and change the node connection
			_used=_used+1;
			temp=&_ring[_used];
			//---------------------------------------------------
			temp->_discard=0;temp->_variable=0;//可以不要
			//----------------------------------------------------
			ratio=ptr->_variable/(ptr->_variable-ptr->_next->_variable);
			temp->_v.x=ptr->_v.x+ratio*(ptr->_next->_v.x-ptr->_v.x);
			temp->_v.y=ptr->_v.y+ratio*(ptr->_next->_v.y-ptr->_v.y);
			temp->_v.z=ptr->_v.z+ratio*(ptr->_next->_v.z-ptr->_v.z);
			temp->_next=ptr->_next;
			ptr->_next=temp;
			ptr=ptr->_next->_next;
		}
		else
		{
			ptr=ptr->_next;
		}
	}while (ptr!=_head);
}
void Ttriangle::cut_outside_node()
{
	Snode* ptr=_head;Snode* start=NULL;Snode* end=NULL;
	while (ptr->_next!=_head)
	{
		if (ptr->_discard==0 && ptr->_next->_discard==1)
		{
			start=ptr;
		}
		if (ptr->_discard==1 && ptr->_next->_discard==0)
		{
			end=ptr->_next;
			if (start!=NULL)
			{
				start->_next=end;
				start=end=NULL;
			}
			//If there is no numerical error, the marked nodes are continuous because the ring must be a convex polygon
			//The numerical error is random and the ring might become an non-convex polygon and the marked node will discontinuous in this circumstance
			//Any time we find an end node, the cutting operation should be performed to avoid this problem
		}
		ptr=ptr->_next;
	}
	return;
}
void Ttriangle::project_onto_plane()
{
	Snode* ptr=_head;Snode* temp=_head;
	do 
	{
		if (ptr->_variable<0)
		{
			ptr->_v.z=1-ptr->_v.x-ptr->_v.y;
		}
		if (ptr->_variable==0)
		{
			//After the projection,the nodes of the final ring may not in a same plane
			//In order to calculate the integration in a same way, we use the node which lies on the x+y+z=1 plane as the head
			temp=ptr;
		}
		ptr=ptr->_next;
	}while(ptr!=_head);
	
	_head=temp;
	return;
}

void Ttriangle::final_integration(double& volume,double& ix,double& iy,double& iz)
{
	Snode* ptr=_head->_next;
	vec3D p1=_head->_v,p2,p3;
	double s,nz;
	volume=ix=iy=iz=0;
	while (ptr->_next!=_head)
	{
		//Convert the volume integration into a surface integration and utilize the Gauss scheme to calculate
		p2=ptr->_v;p3=ptr->_next->_v;
		area(p1,p2,p3,s,nz);
		volume=volume+gauss(p1,p2,p3,"z")*s*abs(nz);
		ix=ix+gauss(p1,p2,p3,"zx")*s*abs(nz);
		iy=iy+gauss(p1,p2,p3,"zy")*s*abs(nz);
		iz=iz+gauss(p1,p2,p3,"zz")*s*abs(nz)*0.5;
		ptr=ptr->_next;
	}
	volume=abs(volume);ix=abs(ix);iy=abs(iy);iz=abs(iz);
	//The integration here must be positive and the final sign is controlled by the nz of the triangle piece
	return;
}
