#ifndef SPARSE_MATRIX
#define SPARSE_MATRIX
//Define the matrix assembly of sparse matrix
#include <vector>
using namespace std;

void calculate_matrix_structure_KII(vector<int> &ID,vector<vector<int>> &IEN,int num_freedom,int num_dimension,int* &ia,int* &ja,int &num_non_zero);
void calculate_matrix_structure_KB(vector<int> &ID, vector<vector<int>> &IEN, int num_ID, int num_dimension, int* &ia, int* &ja, int &num_non_zero);
//Calculte the structure of a sparse matrix for KII (for solving equilibrium equation) 
//and KB (for calculating the reacting force) when a grid topology is given
//Input & Output:
//ID -The number of freedom degree and fixed ID of node i in j direction is ID[3*i+j]
//IEN-The node id connected to a cell
//NE -The cell id connected to a node
//num_freedom-The number of freedom degree of the equation
//num_ID-3*_nump
//num_dimension-The dimension of the problem
//ia,ja-The structure of the sparse matrix

void assmeble_stifness(int i,int j,int* ia,int* ja,double* a,double aij);
//Assemble the aij into sparse matrix K

bool have_freedom(vector<int> &ID,int id,int num_dimension);
//Whether there is a freedom at id
bool have_fixed(vector<int> &ID, int id, int num_dimension);
//Whether there is a constraint at id
#endif