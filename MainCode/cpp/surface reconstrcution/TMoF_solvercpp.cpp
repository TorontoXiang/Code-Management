#include "TMoF_solver.h"
#include <Windows.h>
TMoF_solver::SMoF_input::SMoF_input(Tpolyhedron* poly):_poly(poly),_fraction_reference(0)
{
	_centroid_reference.value(0,0,0);
	_poly->calculate_polyhedron_centroid();
}
TMoF_solver::TMoF_solver(int maximal_num_material,int type,int dimension):_maximal_num_material(maximal_num_material),_derivative_type(type)
	         ,_num_material(0),_num_non_convergence(0),_dimension(dimension)
{
	_fraction_list.resize(_maximal_num_material);
	_centroid_list.resize(_maximal_num_material);
	_result.resize(_maximal_num_material);
	_material_id.resize(_maximal_num_material);
	_cell_poly=NULL;
	for (int i = 0; i < _maximal_num_material; i++)
	{
		_result[i]=NULL;
	}
}
void TMoF_solver::set_initial_condition(Tpolyhedron* cell_poly,double fraction_list[],vec3D centroid_list[],int num_material)
{
	if (num_material!=_maximal_num_material)
	{
		cout<<"Error:The number of material is in invalid"<<endl;
		system("Pause");
		exit(0);
	}
	for (int i = 0; i < _maximal_num_material; i++)
	{
		_fraction_list[i]=fraction_list[i];
		_centroid_list[i]=centroid_list[i];
	}
	_cell_poly=cell_poly;
	_cell_poly->calculate_polyhedron_centroid();
}
void TMoF_solver::Finalize()
{
	_cell_poly=NULL;
	for (int i = 0; i < _num_material; i++)
	{
		_fraction_list[i]=0;
		_centroid_list[i].value(0,0,0);
		_result[i]=NULL;
	}
	return;
}
double TMoF_solver::MoF_surface_reconstrcution()
{
	double result;
	classify_mix_cell();
	for (int i = 0; i < _maximal_num_material; i++)
	{
		_result[i]=NULL;
	}
	if (_num_material==1)
	{
		//This is a pure material
		int material_id=_material_id[0];
		_result[material_id]=new Tpolyhedron;
		*_result[material_id]=*_cell_poly;
		result=0;
	}
	else if (_num_material==2)
	{
		//Two materials in this cell
		double fraction[2];
		vec3D centroid[2];
		Tpolyhedron* sub_poly[2];
		for (int i = 0; i < 2; i++)
		{
			fraction[i]=_fraction_list[_material_id[i]];centroid[i]=_centroid_list[_material_id[i]];
			sub_poly[i]=NULL;
		}
		MoF_for_pairs(_cell_poly,fraction,centroid,sub_poly);
		for (int i = 0; i < 2; i++)
		{
			_result[_material_id[i]]=sub_poly[i];
		}
	}
	else if (_num_material==3)
	{
		//Three materials in this cell
		double min_error=1e10;
		Tpolyhedron* trial_sub_poly[3][3];
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				trial_sub_poly[i][j]=NULL;
			}
		}
		for (int i = 0; i < 3; i++)
		{
			int pre=_num_non_convergence;
			double error_i=three_material_trial(i,trial_sub_poly[i]);
			if (error_i<min_error)
			{
				//This trial solution is temporary minimal,update min_error and _result
				min_error=error_i;
				for (int j = 0; j < 3; j++)
				{
					if (_result[j]==NULL)
					{
						_result[j]=trial_sub_poly[i][j];
					}
					else
					{
						//Delete the previous polyhedron to avoid memory leak
						delete _result[j];
						_result[j]=trial_sub_poly[i][j];
					}
				}
			}
			else
			{
				//This trial solution is not the optimal solution,delete the generated sub-polyhedrons
				for (int j = 0; j < 3; j++)
				{
					delete trial_sub_poly[i][j];
				}
			}
		}
		result=min_error;
	}
	return result;
}
void TMoF_solver::classify_mix_cell()
{
	_num_material=0;
	for (int i = 0; i < _maximal_num_material; i++)
	{
		if (_fraction_list[i]>0)
		{
			_material_id[_num_material]=i;
			_num_material=_num_material+1;
		}
	}
	return;
}
void TMoF_solver::MoF_for_pairs(Tpolyhedron* poly_in,double (&fraction)[2],vec3D (&centroid)[2],Tpolyhedron* (&sub_poly)[2])
{
	vec2D result;
	if (fraction[0]==0)
	{
		sub_poly[0]=NULL;
		sub_poly[1]=new Tpolyhedron;
		*sub_poly[1]=*poly_in;
	}
	else if (fraction[0]==1)
	{
		sub_poly[0]=new Tpolyhedron;
		*sub_poly[0]=*poly_in;
		sub_poly[1]=NULL;
	} 
	else
	{
		bool is_alter=false;
		SMoF_input MoF_input(poly_in);
		vec2D initial_geuss,optimal_direction;
		vec3D n;
		bool is_failed;
		if (fraction[0]<=0.5)
		{
			MoF_input._fraction_reference=fraction[0];
			MoF_input._centroid_reference=centroid[0];
		}
		else
		{
			MoF_input._fraction_reference=fraction[1];
			MoF_input._centroid_reference=centroid[1];
			is_alter=true;
		}
		_volume_equation_solver.set_target_polyhedron(poly_in);
		initial_geuss=calculate_initial_guess(MoF_input);
		optimal_direction=BFGS_iteration(MoF_input,initial_geuss,is_failed);
		double theta=optimal_direction.x,phy=optimal_direction.y;
		n.value(sin(theta)*cos(phy),sin(theta)*sin(phy),cos(theta));
		double d=_volume_equation_solver.calculate_plane_constant_advanced(n,MoF_input._fraction_reference);
		Tpolyhedron *poly_below,*poly_above;
		poly_below=new Tpolyhedron;poly_above=new Tpolyhedron;
		poly_in->cut_polyhedron_by_plane(n,d,poly_below,poly_above);
		if (!is_alter)
		{
			sub_poly[0]=poly_below;sub_poly[1]=poly_above;
		}
		else
		{
			sub_poly[1]=poly_below;sub_poly[0]=poly_above;
		}
	}
	return;
}
double TMoF_solver::three_material_trial(int i,Tpolyhedron* (&sub_poly)[3])
{
	double result=0;
	double fraction[2];
	vec3D centroid[2];
	Tpolyhedron* temp_sub_poly1[2];
	Tpolyhedron* temp_sub_poly2[2];
	//----------------------------------------------------------------------------
	//The first division
	//----------------------------------------------------------------------------
	fraction[0]=_fraction_list[i];fraction[1]=1-fraction[0];
	centroid[0]=_centroid_list[i];centroid[1]=(_cell_poly->G_centroid()-centroid[0]*fraction[0])/fraction[1];
	temp_sub_poly1[0]=temp_sub_poly1[1]=NULL;
	MoF_for_pairs(_cell_poly,fraction,centroid,temp_sub_poly1);
	//After MoF_for_pairs process, the material polyhedron i is temp_sub_poly1[0]
	//material (i+1)%3 and (i+2)%3 are in temp_sub_poly1[1]
	sub_poly[i]=temp_sub_poly1[0];
	//-------------------------------------------------------------------------------
	//The second division
	//--------------------------------------------------------------------------------
	fraction[0]=_fraction_list[(i+1)%3]/(_fraction_list[(i+1)%3]+_fraction_list[(i+2)%3]);fraction[1]=1-fraction[0];
	centroid[0]=_centroid_list[(i+1)%3];centroid[1]=_centroid_list[(i+2)%3];
	temp_sub_poly2[0]=temp_sub_poly2[1]=NULL;
	MoF_for_pairs(temp_sub_poly1[1],fraction,centroid,temp_sub_poly2);
	sub_poly[(i+1)%3]=temp_sub_poly2[0];sub_poly[(i+2)%3]=temp_sub_poly2[1];
	//Delete the temporary polyhedron
	delete temp_sub_poly1[1];
	//Finally calculate the discrepance of approxiamte centroids and the reference centroids
	for (int j = 0; j < 3; j++)
	{
		result=result+(_centroid_list[j]-sub_poly[j]->calculate_polyhedron_centroid()).self_multuply();
	}
	return result;
}
vec2D TMoF_solver::BFGS_iteration(SMoF_input &MoF_input,vec2D initial_geuss,bool &is_failed)
{
	Siterator_point iterator_pre,iterator_next;   //The iterator point
	double theta,phy;
	vec3D n;                                 //The coressponding normal at each iterative point
	vec2D d_derivative;                      //The plane constant derivative
	vec2D direction;                         //The searching direction
	double eps=1e-9;                         //The termination conditions
	double L2;                               //The L2 normal of gradiant
	int num=0;                               //Number of iterative steps in BGFS process
	//Intermediary variables
	vec2D ss,yy;       
	matrix2D H={1,0,0,1},temp1,temp2;
	double temp3;
	//double predict_d;
	vec3D iterator_info;

	double fraction_reference=MoF_input._fraction_reference;
	vec3D centroid_reference=MoF_input._centroid_reference;
	theta=initial_geuss.x;phy=initial_geuss.y;
	n.value(sin(theta)*cos(phy),sin(theta)*sin(phy),cos(theta));
	iterator_next._d=_volume_equation_solver.calculate_plane_constant_advanced(n,fraction_reference);
	iterator_info=calculate_iterator_info(MoF_input,theta,phy,iterator_next._d,d_derivative);
	iterator_next._point=initial_geuss;
	iterator_next._value=iterator_info.x;iterator_next._gradient.value(iterator_info.y,iterator_info.z);
	iterator_next._d_derivative=d_derivative;
	L2=sqrt(iterator_next._gradient*iterator_next._gradient);
	//Iteration begin
	//convergence=true;
	while (L2>eps && num<=100)
	{
		iterator_pre=iterator_next;
		//der_d_pre=der_d_next;
		direction=H*iterator_pre._gradient*(-1);
		//double df=direction*grad_pre;
		iterator_next=line_searching(MoF_input,iterator_pre,direction,is_failed);
		num=num+1;
		if (num==100)
		{
			//If the iterative steps is lager than 100, the iteration fails
			//convergence=false;
			_num_non_convergence=_num_non_convergence+1;
		}
		L2=sqrt(iterator_next._gradient*iterator_next._gradient);
		if (is_failed)
		{
			//If the line searching is failed, the quasi-Newton matrix will set to be the unit matrix
			H.a11=H.a22=1;H.a12=H.a21=0;
		}
		else
		{
			//Update the quasi Newton matrix
			ss=iterator_next._point-iterator_pre._point;
			yy=iterator_next._gradient-iterator_pre._gradient;
			if (sqrt(ss*ss)<1e-20)
			{
				//If iterator_pre and iterator_next is very close,the linear searching is failed
				H.a11=H.a22=1;H.a12=H.a21=0;
			}
			else
			{
				temp1.a11=yy.x*ss.x;temp1.a12=yy.x*ss.y;temp1.a21=yy.y*ss.x;temp1.a22=yy.y*ss.y;
				temp2.a11=ss.x*ss.x;temp2.a12=temp2.a21=ss.x*ss.y;temp2.a22=ss.y*ss.y;
				temp3=ss*yy;
				H=H-(H*temp1+(H*temp1).transfer())/temp3+temp2*(1+H%yy/temp3)/temp3;
			}
		}
	}
	return iterator_next._point;
}
TMoF_solver::Siterator_point TMoF_solver::line_searching(SMoF_input &MoF_input,Siterator_point &iterator,vec2D direction,bool &is_failed)
{
	double fracrion_reference=MoF_input._fraction_reference;
	vec3D centroid_reference=MoF_input._centroid_reference;
	vec3D n;
	double b1=0.1;                            //The parameter for objective function decrement
	double b2=1;                              //The parameter for derivative decrement
	double index=0.1;                         //The parameter for searching area constraint
	double d_min=0,d_max=2e10;                //The searching range,from 0-infinite at first
	double d_iterator=1;                      //The iterator distance in line searching,iterator=1 at first
	double d_pre;                             //The previous iterator
	double f_min,f_max,f_iterator;            //The values of objective function at specific points
	double df_min,df_max,df_iterator;         //The derivative alone this line at specific points 
	double theta0,phy0;                       //The start point
	double dtheta0,dphy0;                     //The unit increment of theta and phy
	double f0;                                //The initial value of objective function
	double theta_iterator,phy_iterator;       //The theta and phy at d_iterator
	double df_theta,df_phy;                   //The partial derivative with respect to theta and phy at d_iterator
	double condition;                         //The termination condition
	double d;                                 //The plane constant at each iterator point
	vec2D d_derivative;                       //The plane constant derivative at each iterator point
	double dtheta,dphy;                       //The increment of theta and phy in each step 
	int num=0;                                //The mumber of iterative steps in this 1D searching
	vec3D iterator_info;                      //The information at each iterator point
	double predict_d;

	theta0=iterator._point.x;phy0=iterator._point.y;
	dtheta0=direction.x;dphy0=direction.y;
	f0=iterator._value;
	d=iterator._d;
	d_derivative=iterator._d_derivative;
	condition=direction*iterator._gradient;
	f_min=iterator._value;df_min=direction*iterator._gradient;df_max=-1;    //df2 is unknown at first, set to be negative
	//Decrement examination
	if (condition>=0)
	{
		//If the searching direction is not the decrement direction,output warning
		cout<<"Warning: The searching direction is not the decrement direction"<<endl;
		is_failed=true;
		return iterator;
	}
	//Information at the initial iterative point
	theta_iterator=theta0+d_iterator*dtheta0;phy_iterator=phy0+d_iterator*dphy0;   
	n.value(sin(theta_iterator)*cos(phy_iterator),sin(theta_iterator)*sin(phy_iterator),cos(theta_iterator));
	dtheta=d_iterator*dtheta0;dphy=d_iterator*dphy0;
	d_pre=d_iterator;
	predict_d=d+d_iterator*dtheta0*d_derivative.x+d_iterator*dphy0*d_derivative.y;
	//predict_d=predict_altitude(der_current,d_current,ith,iphy);
	d=_volume_equation_solver.calculate_plane_constant_advanced(n,fracrion_reference,predict_d);
	iterator_info=calculate_iterator_info(MoF_input,theta_iterator,phy_iterator,d,d_derivative);
	f_iterator=iterator_info.x;df_theta=iterator_info.y;df_phy=iterator_info.z;
	//Iteration begin
	is_failed=false;
	while (true)
	{
		while (f0-f_iterator<-d_iterator*b1*condition)       //Condition for objective function decrement
		{
			//If the objective function is not sufficiently decreased,
			//using cubic or quad function to calculate the trial optimal point
			d_max=d_iterator;f_max=f_iterator;df_max=direction.x*df_theta+direction.y*df_phy;
			if (abs(d_max-d_min)<1e-10)
			{
				//The sufficient decrement can not be attained, this 1D searching is failed
				is_failed=true;
				return iterator;
			}
			if (df_max>0)
			{
				//The derivative at d_max is positive, using cubic interpolation to get the trial point
				d_iterator=cubic1(d_min,d_max,f_min,df_min,f_max,df_max);
			}
			else
			{
				//The derivative at d_max is negative, using quad interpolation to get the trial point
				d_iterator=quad(d_min,d_max,f_min,df_min,f_max);
			}
			//Constraint a into a specific area
			//The length of [d_min,d_max] decrease min(index,1-index) at least
			d_iterator=minval(maxval(d_iterator,d_min+index*(d_max-d_min)),d_max-index*(d_max-d_min));
			//The information at the next iterator point
			theta_iterator=theta0+d_iterator*dtheta0;phy_iterator=phy0+d_iterator*dphy0;
			n.value(sin(theta_iterator)*cos(phy_iterator),sin(theta_iterator)*sin(phy_iterator),cos(theta_iterator));
			dtheta=(d_iterator-d_pre)*dtheta0;dphy=(d_iterator-d_pre)*dphy0;
			d_pre=d_iterator;
			//predict_d=predict_altitude(der_current,d_current,ith,iphy);
			predict_d=d+dtheta*d_derivative.x+dphy*d_derivative.y;
			d=_volume_equation_solver.calculate_plane_constant_advanced(n,fracrion_reference,predict_d);
			iterator_info=calculate_iterator_info(MoF_input,theta_iterator,phy_iterator,d,d_derivative);
			f_iterator=iterator_info.x;df_theta=iterator_info.y;df_phy=iterator_info.z;
			num=num+1;
			//Output warning message if the iterative step is too much
			if (num>1000)
			{
				//cout<<"Warning: The iterative step in 1D searching is too much"<<endl;
				is_failed=true;
				return iterator;
			}
		}
		//The objective decrement requirement is fulfilled and begin to fulfill the derivative decrement requirement
		df_iterator=direction.x*df_theta+direction.y*df_phy;
		if (df_iterator>b2*condition)
		{
			//The derivative requirement is fulfilled, output the information
			Siterator_point temp;
			temp._point.value(theta_iterator,phy_iterator);
			temp._gradient.value(df_theta,df_phy);
			temp._value=f_iterator;
			temp._d=d;
			temp._d_derivative=d_derivative;
			return temp;
		}
		else if (df_iterator>0)
		{
			//The derivative at current point is positive, using cubic interpolation to get the next point
			d_max=d_iterator;f_max=f_iterator;df_max=df_iterator;
			d_iterator=cubic1(d_min,d_max,f_min,df_min,f_max,df_max);
			d_iterator=minval(maxval(d_iterator,d_min+index*(d_max-d_min)),d_max-index*(d_max-d_min));
		}
		else if (df_iterator<0)
		{
			//The derivative at current point is negative, a1=a
			if (d_max>1e10)
			{
				//The searching range is infinite, the iterative step should be amplified 
				double temp=cubic2(d_min,d_iterator,f_min,df_min,f_iterator,df_iterator);
				//The d_iterator will at least going ahead (d_iterator-d_min),at most 9*(d_iterator-d_min)
				temp=minval(maxval(temp,d_iterator+(d_iterator-d_min)),d_iterator+9*(d_iterator-d_min));
				d_min=d_iterator;f_min=f_iterator;df_min=df_iterator;
				d_iterator=temp;
			}
			else
			{
				//The searching rang is finite, interpolation is used
				d_min=d_iterator;f_min=f_iterator;df_min=df_iterator;
				if (df_max>0)
				{
					//Cubic interpolation for df_max>0
					d_iterator=cubic1(d_min,d_max,f_min,df_min,f_max,df_max);
				}
				else
				{
					//Quad interpolation for df_max<0
					d_iterator=quad(d_min,d_max,f_min,df_min,f_max);
				}
				d_iterator=minval(maxval(d_iterator,d_min+index*(d_max-d_min)),d_max-index*(d_max-d_min));
			}
		}
		if (abs(d_max-d_min)<1e-10)
		{
			//The sufficient decrement can not be attained, this 1D searching is failed
			//The input information should not be changed
			is_failed=true;
			return iterator;
		}
		//The information at the next iterator point
		theta_iterator=theta0+d_iterator*dtheta0;phy_iterator=phy0+d_iterator*dphy0;
		n.value(sin(theta_iterator)*cos(phy_iterator),sin(theta_iterator)*sin(phy_iterator),cos(theta_iterator));
		dtheta=(d_iterator-d_pre)*dtheta0;dphy=(d_iterator-d_pre)*dphy0;
		d_pre=d_iterator;
		predict_d=d+dtheta*d_derivative.x+dphy*d_derivative.y;
		//predict_d=predict_altitude(der_current,d_current,ith,iphy);
		d=_volume_equation_solver.calculate_plane_constant_advanced(n,fracrion_reference,predict_d);
		iterator_info=calculate_iterator_info(MoF_input,theta_iterator,phy_iterator,d,d_derivative);
		f_iterator=iterator_info.x;df_theta=iterator_info.y;df_phy=iterator_info.z;
		num=num+1;
		//Output warning message if the iterative step is too much
		if (num>1000)
		{
			//cout<<"Warning: The iterative step in 1D searching is too much"<<endl;
			is_failed=true;
			return iterator;
		}
	}
	return iterator;
}
vec3D TMoF_solver::calculate_iterator_info(SMoF_input &MoF_input,double theta,double phy,double d,vec2D &d_derivative)
{
	vec3D result;
	double fraction_reference=MoF_input._fraction_reference;
	vec3D centroid_reference=MoF_input._centroid_reference;
	Tpolyhedron* poly=MoF_input._poly;
	vec3D centroid_below;           //The centroid below the cutting plane
	vec3D centroid_ring;            //The centroid of the cutting ring
	vec3D n;                        //The normal of the approximate surface
	Scutting_ring* cutting_ring;    
	cutting_ring=poly->G_cutting_ring();
	n.x=sin(theta)*cos(phy),n.y=sin(theta)*sin(phy),n.z=cos(theta);
	centroid_below=poly->calcualte_centroid_below(n,d);
	centroid_ring=cutting_ring->calculate_surface_integral("centroid");
	d_derivative.x=cos(theta)*cos(phy)*centroid_ring.x+cos(theta)*sin(phy)*centroid_ring.y-sin(theta)*centroid_ring.z;
	d_derivative.y=-sin(theta)*sin(phy)*centroid_ring.x+sin(theta)*cos(phy)*centroid_ring.y;
	result.x=(centroid_below-centroid_reference).self_multuply();
	double poly_volume=poly->G_volume(),sth=sin(theta);
	if (_derivative_type==0)
	{
		//Semi-analytical derivatives
		vec3D l1,l2,l3;              //The new coordinate system: the 3 axises direction and the origin point
		vec3D temp_th,temp_phy;      //Vector x-x_ref at new coordinate system
		vec3D integration;           //Integration on the cutting plane
		cutting_ring->calculate_local_coordinate_system(theta,phy,l1,l2,l3);
		temp_th=centroid_below.transfer_into_new_system(l1,l2,l3,centroid_ring)-centroid_reference.transfer_into_new_system(l1,l2,l3,centroid_ring);
		temp_phy.x=temp_th.y;temp_phy.y=-temp_th.x;temp_phy.z=temp_th.z;
		for (int i=0;i<cutting_ring->_length;i++)
		{
			cutting_ring->_point_list[i]=cutting_ring->_point_list[i].transfer_into_new_system(l1,l2,l3,centroid_ring);
		}
		integration=cutting_ring->calculate_surface_integral("derivative");
		result.y=-(temp_th.x*integration.x+temp_th.y*integration.y)*2/(fraction_reference*poly_volume);
		result.z=(temp_phy.x*integration.x-temp_phy.y*integration.z)*2/(fraction_reference*poly_volume)*sth;
	}
	else if (_derivative_type==1)
	{
		//differential derivatives
		double dd=1e-6;                //The differential step
		double f;                      //Objective function at (th,phy)
		//Intermediary variables
		double d,predict_d;

		predict_d=100000;
		//Derivative for th
		theta=theta+dd;
		n.x=sin(theta)*cos(phy),n.y=sin(theta)*sin(phy),n.z=cos(theta);
		d=_volume_equation_solver.calculate_plane_constant_advanced(n,fraction_reference);
		centroid_below=poly->calcualte_centroid_below(n,d);
		f=(centroid_below-centroid_reference).self_multuply();
		result.y=(f-result.x)/dd;
		//Derivative for phy
		theta=theta-dd;phy=phy+dd;
		n.x=sin(theta)*cos(phy),n.y=sin(theta)*sin(phy),n.z=cos(theta);
		d=_volume_equation_solver.calculate_plane_constant_advanced(n,fraction_reference);
		centroid_below=poly->calcualte_centroid_below(n,d);
		f=(centroid_below-centroid_reference).self_multuply();
		result.z=(f-result.x)/dd;
		phy=phy-dd;
	}
	if (_dimension==2)
	{
		result.y=0;
	}
	return result;
}
vec2D TMoF_solver::calculate_initial_guess(SMoF_input &MoF_input)
{
	vec2D initial_geuss;
	vec3D poly_centroid=MoF_input._poly->G_centroid();
	vec3D centroid_reference=MoF_input._centroid_reference;
	double pi=3.141592653;
	vec3D temp;               //The unit vector from centroid_reference to polyhedron centroid
	double theta,phy;
	double s;
	temp=poly_centroid-centroid_reference;
	temp.normalize();
	//n=(sin(theta)cos(phy),sin(theta)sin(phy),cos(theta))
	if (temp.z==1)
	{
		//theta=0 and phy can be any variable
		theta=0;phy=0;
		initial_geuss.x=theta;initial_geuss.y=phy;
		return initial_geuss;
	}
	theta=acos(temp.z);  //acos() is in [0,pi],so theta=acos(temp.z)
	if (temp.y>=0)
	{
		if (temp.x>=0)
		{
			s=temp.y/sin(theta);
			if (s>=1)
			{
				phy=pi/2.0;
			}
			else
			{
				phy=asin(temp.y/sin(theta));
			}

		}
		else
		{
			s=temp.y/sin(theta);
			if (s>=1)
			{
				phy=pi/2.0;
			}
			else
			{
				phy=pi-asin(temp.y/sin(theta));
			}
		}
	}
	else
	{
		if (temp.x>=0)
		{
			s=temp.y/sin(theta);
			if (s<=-1)
			{
				phy=1.5*pi;
			}
			else
			{
				phy=2*pi+asin(temp.y/sin(theta));
			}	
		}
		else
		{
			s=temp.y/sin(theta);
			if (s<=-1)
			{
				phy=1.5*pi;
			}
			else
			{
				phy=pi-asin(temp.y/sin(theta));
			}	
		}
	}
	if (_dimension==2)
	{
		theta=pi/2;
	}
	initial_geuss.x=theta;initial_geuss.y=phy;
	return initial_geuss;
}
vec2D TMoF_solver::limit_angle(vec2D direction)
{
	double theta=direction.x,phy=direction.y;
	double pi=3.141592654;
	while (theta<0 || theta>pi)
	{
		if (theta<0)
		{
			theta=-theta;
			phy=phy+pi;
		}
		else if (theta>pi)
		{
			theta=2*pi-theta;
			phy=phy+pi;
		}
	}
	while (phy<0 || phy>2*pi)
	{
		if (phy<0)
		{
			phy=phy+2*pi;
		}
		else if (phy>2*pi)
		{
			phy=phy-2*pi;
		}
	}
	vec2D temp;
	temp.value(theta,phy);
	return temp;
}
//--------------------------------------------------------------------------------------------
//Test functions
//--------------------------------------------------------------------------------------------
vec2D TMoF_solver::derivative_comparison(Tpolyhedron* poly_in,double fraction,vec3D centroid)
{
	vec2D result;
	SMoF_input MoF_input(poly_in);
	MoF_input._fraction_reference=fraction;
	MoF_input._centroid_reference=centroid;
	double theta,phy;
	double pi=3.141592653;
	vec3D n;
	_volume_equation_solver.set_target_polyhedron(poly_in);
	int i_max_error_theta,j_max_error_theta;
	int i_max_error_phy,j_max_error_phy;
	int n_divide=100;
//	for (int i = 37; i < 38; i++)
//	{
//		for (int j = 7; j < 8; j++)
//		{
			vec2D initial_guess=calculate_initial_guess(MoF_input);
			theta=initial_guess.x;phy=initial_guess.y;
			//theta=3;phy=5.5;
			n.x=sin(theta)*cos(phy),n.y=sin(theta)*sin(phy),n.z=cos(theta);
			double d=_volume_equation_solver.calculate_plane_constant_advanced(n,MoF_input._fraction_reference);
			vec3D result_numerical,result_analytical;
			vec2D d_derivative;
			_derivative_type=0;
			result_analytical=calculate_iterator_info(MoF_input,theta,phy,d,d_derivative);
			_derivative_type=1;
			result_numerical=calculate_iterator_info(MoF_input,theta,phy,d,d_derivative);
			_derivative_type=0;
			//if (abs((result_analytical-result_numerical).y)>result.x)
			//{
			//	result.x=(result_analytical-result_numerical).y;
			//	i_max_error_theta=i;j_max_error_theta=j;
			//}
			//if (abs((result_analytical-result_numerical).z)>result.y)
			//{
			//	result.y=(result_analytical-result_numerical).z;
			//	i_max_error_phy=i;j_max_error_phy=j;
			//}
		//}
	//}
	
	cout<<"i,j for d_theta max: "<<i_max_error_theta<<" "<<j_max_error_theta<<endl;
	cout<<"i,j for d_phy max: "<<i_max_error_phy<<" "<<j_max_error_phy<<endl;
	return result;
}
vec2D TMoF_solver::patch_test(Tpolyhedron* poly_in,double fraction)
{
	vec2D result;
	double pi=3.141592654;
	_volume_equation_solver.set_target_polyhedron(poly_in);
	vec3D n;
	double theta,phy;
	int n_divide=100;
	double error_max_theta=0,error_max_phy=0;
	int i_max_error_theta,j_max_error_theta;
	int i_max_error_phy,j_max_error_phy;
	double error_phy;
	for (int i = 1; i < n_divide-1; i++)
	{
		for (int j = 0; j < n_divide+1; j++)
		{
			theta=i*pi/n_divide;phy=2*pi*j/n_divide;
			n.value(sin(theta)*cos(phy),sin(theta)*sin(phy),cos(theta));
			double d=_volume_equation_solver.calculate_plane_constant_advanced(n,fraction);
			vec3D centroid_reference=poly_in->calcualte_centroid_below(n,d);
			SMoF_input MoF_input(poly_in);
			MoF_input._fraction_reference=fraction;
			MoF_input._centroid_reference=centroid_reference;
			bool is_failed;
			vec2D initial_geuss=calculate_initial_guess(MoF_input);
			vec2D optimal_direction=BFGS_iteration(MoF_input,initial_geuss,is_failed);
			optimal_direction=limit_angle(optimal_direction);
			if (abs(theta-optimal_direction.x)>error_max_theta)
			{
				error_max_theta=abs(theta-optimal_direction.x);
				i_max_error_theta=i;j_max_error_theta=j;
			}
			error_phy=minval(abs(phy-optimal_direction.y),2*pi-abs(phy-optimal_direction.y));
			if (error_phy>error_max_phy)
			{
				error_max_phy=error_phy;
				i_max_error_phy=i;j_max_error_phy=j;
			}
		}
	}
	cout<<"maximal error of theta is "<<error_max_theta<<" at "<<i_max_error_theta<<","<<j_max_error_theta<<endl;
	cout<<"maximal error of phy is "<<error_max_phy<<" at "<<i_max_error_phy<<","<<j_max_error_phy<<endl;
	result.value(error_max_theta,error_max_phy);
	return result;
}
double TMoF_solver::three_material_patch_test(Tpolyhedron* poly_in,double (&fraction)[2],vec2D (&angle)[2])
{
	//ofstream output;
	//output.open("MoF_test.dat");
	//poly_in->plot_this_polyehdron(output);
	//poly_in->calculate_polyhedron_centroid();
	double pi=3.141592654;
	vec3D n1,n2,centroid1,centroid2,centroid3;
	double theta1=angle[0].x,phy1=angle[0].y;
	double theta2=angle[1].x,phy2=angle[1].y;
	double fraction1=fraction[0],fraction2=fraction[1]/(1-fraction1);
	n1.value(sin(theta1)*cos(phy1),sin(theta1)*sin(phy1),cos(theta1));
	n2.value(sin(theta2)*cos(phy2),sin(theta2)*sin(phy2),cos(theta2));
	_volume_equation_solver.set_target_polyhedron(poly_in);
	double d1=_volume_equation_solver.calculate_plane_constant_advanced(n1,fraction1);
	centroid1=poly_in->calcualte_centroid_below(n1,d1);
	Tpolyhedron *poly_below,*poly_above;
	poly_below=new Tpolyhedron;poly_above=new Tpolyhedron;
	poly_in->cut_polyhedron_by_plane(n1,d1,poly_below,poly_above);
	_volume_equation_solver.set_target_polyhedron(poly_above);
	double d2=_volume_equation_solver.calculate_plane_constant_advanced(n2,fraction2);
	centroid2=poly_above->calcualte_centroid_below(n2,d2);
	delete poly_below;delete poly_above;
	centroid3=(poly_in->G_centroid()-centroid1*fraction[0]-centroid2*fraction[1])/(1-fraction[0]-fraction[1]);
	//Set the initial condition of 3 material MoF surface reconstruction 
	_fraction_list[0]=fraction[0];_fraction_list[1]=fraction[1];_fraction_list[2]=1-fraction[0]-fraction[1];
	_centroid_list[0]=centroid1;_centroid_list[1]=centroid2;_centroid_list[2]=centroid3;
	_cell_poly=poly_in;
	double centroid_error=MoF_surface_reconstrcution();
	//Release the memory of approximate polyhedron
	for (int i = 0; i < 3; i++)
	{
		delete _result[i];
		_result[i]=NULL;
	}
	return centroid_error;
}
