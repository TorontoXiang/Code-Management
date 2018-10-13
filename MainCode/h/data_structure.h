//Define the data structure in program
#ifndef DATA_STRUCTURE
#define DATA_STRUCTURE
//3D vector
struct vec3D
{
	double x,y,z;

	vec3D();
	vec3D(double xx,double yy,double zz){x=xx,y=yy,z=zz;};
	//Normalize a vec3D
	void normalize();
	//Value a vec3D
	void value(double const xx,double const yy,double const zz);
	//Value the ith variable 
	void value(int i,double num);
	//Self multiply
	double self_multuply() const;
	//Get length
	double get_length() const;
	//Multiply by a matrix
	vec3D multiply_by_matrix(double (&A)[3][3]);
	//Transfer the point into a new coordinate system
	vec3D transfer_into_new_system(vec3D &n1,vec3D &n2,vec3D &n3,vec3D &x0);
	//Get ith variable
	double access(int i);

	//Opreations
	bool operator ==(vec3D const &other);
	vec3D operator +(vec3D const &other);
	vec3D operator -(vec3D const &other);
	vec3D operator -();
	vec3D operator *(double const a);
	vec3D operator /(double const a);
	vec3D operator %(vec3D const &other);    //vector product
	double operator *(vec3D const &other);  //dot product
	bool operator <(vec3D const &other);
	bool operator >(vec3D const &other);
};
//A coordiante with a double variable
struct vec3D_double
{
	vec3D _coor;
	double _variable;
	vec3D_double():_variable(0){};
	vec3D_double(vec3D &coor):_coor(coor),_variable(0){};
};
//A coordinate with a int type
struct vec3D_int
{
	vec3D _coor;
	int _type;
	vec3D_int():_type(0){};
	vec3D_int(vec3D &coor):_coor(coor),_type(0){};
};
//2D vector
struct vec2D
{
	double x;
	double y;

	vec2D():x(0),y(0){};
	vec2D(double xx,double yy):x(xx),y(yy){};
	void value(double xx,double yy){x=xx;y=yy;};
	double operator *(vec2D other);     //Dot product
	vec2D operator *(double a);       
	vec2D operator -(vec2D other);      
	vec2D operator /(double a); 
	double get_length();
};
//2*2 matrix
struct matrix2D
{
	double a11,a12,a21,a22;
	matrix2D operator +(matrix2D other);   
	matrix2D operator -(matrix2D other);   
	matrix2D operator *(matrix2D other);   //Matrix multiply Matrix
	matrix2D operator *(double a);         //Matrix multiply Matrix
	double operator %(vec2D v);            //Multiply twice into a number
	vec2D operator *(vec2D v);             //Multiply once into a vec
	matrix2D operator /(double a);         //Divided by a number
	matrix2D transfer();                   //Transfer the matrix
};
//The topology of a cell
struct cell_topology8H
{
	int id;
	int IEN[8];
};
#endif