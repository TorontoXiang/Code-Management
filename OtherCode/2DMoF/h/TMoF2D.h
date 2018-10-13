#include "Tpolygon.h"
#include "public.h"
#include <vector>
#include <algorithm>
#include <functional>
using namespace std;

class TMoF2D
{
public:
	TMoF2D(int highest_order,int num_vertex);

	void Access_to_polygon(Tpolygon* polygon_ptr);
	//Access the MoF solver to the polygon pointer
	void Calculate_interation_information(double theta);
	//Calculate the value and derivatives of the objective function at theta
	void Calculate_iteration_information_new(double dtheta);

	double MoF_solver();
	double MoF_solver_new();
	//Calculate the feasible region of the optimal solution and the initial guess
	double iteration(iteration_point theta_min,iteration_point theta_max,double theta_initial);
	//Solve the iteration problem and return the optimal theta
	void Calculate_derivatives(vec2D& n,double length,vec2D& centroid_below);
	//Calculate the derivatives of the objective function

	void S_fration(double fraction){_fraction=fraction;};
	void S_highest_order(int order){_hightest_order=order;};
	void S_centroid(vec2D centroid){_centroid_ref=centroid;};
	void S_centroid(double x1,double x2){_centroid_ref.value(x1,x2);};
	double G_value(){return _iteration_point.value;};
	double G_derivative(int i);

	//Comparison with previous method
	double initial_guess();
	void Calculate_single_summit_area(double& theta_min,double& theta_max,double& theta_initial);
	double Zoom_algorithm(double theta_min,double theta_max);

	//Patch test function
	vec2D calculate_centroid_below(double theta);
public:
//private:
	Tpolygon* _polygon_ptr;     //The cell polygon to be reconstructed
	int _num_vertex;            //The number of vertex of the polygon
	double* _inclination[2];    //The inclination of each edge
	int _hightest_order;        //The highest order derivatives
	iteration_point _iteration_point;   //The current iteration point
	double _fraction;           //The reference fraction
	vec2D _centroid_ref;        //The reference centroid

	double _alpha1,_alpha2,_s1,_c1,_s2,_c2,_theta0;
	vec2D _pb1,_pb2,_d1,_d2,_i1,_i2,_centroid0;
};