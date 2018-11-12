#include "Sparse_matrix.h"
#include <algorithm>
#include <iostream>
using namespace std;
void calculate_matrix_structure_KII(vector<int> &ID,vector<vector<int>> &IEN,int num_freedom,int num_dimension,int* &ia,int* &ja,int &num_non_zero)
{
	int nump=ID.size()/num_dimension,nume=IEN.size();
	vector<int> node_flag;            //The node flag
	vector<int> node_c_node;   
	vector<int> node_c_freedom;
	vector<vector<int>> freedom_c_freedom;
	ia=new int[num_freedom+1];
	ia[0]=1;
	freedom_c_freedom.resize(num_freedom);

	//Initialize the node_flag
	node_flag.resize(nump);
	for (int i = 0; i < nump; i++)
	{
		node_flag[i]=0;
	}

	//Calculate the cells connected to a node
	vector<vector<int>> NE;
	NE.resize(nump);
	for (int i = 0; i < nume; i++)
	{
		int num_node_in_cell=IEN[i].size();
		for (int j = 0; j < num_node_in_cell; j++)
		{
			int node_id=IEN[i][j]-1;
			NE[node_id].push_back(i);
		}
	}

	//Calculate freedom_c_freedom
	for (int i = 0; i < nump; i++)
	{
		if (have_freedom(ID,i,num_dimension))
		{
			//If there is a freedom at the node, calculate the freedom id connected to this node

			//Calculate the node_c_node for node i
			int num_connected_cell=NE[i].size();
			for (int j = 0; j < num_connected_cell; j++)
			{
				int cell_id=NE[i][j];
				int num_node_in_cell=IEN[cell_id].size();
				for (int k = 0; k < num_node_in_cell; k++)
				{
					int node_id=IEN[cell_id][k]-1;
					if (node_flag[node_id]==0)
					{
						node_c_node.push_back(node_id+1);
						node_flag[node_id]=1;
					}
				}
			}
			//Clear the node_flag after obtain the node_c_node
			for (int j = 0; j < num_connected_cell; j++)
			{
				int cell_id=NE[i][j];
				int num_node_in_cell=IEN[cell_id].size();
				for (int k = 0; k < num_node_in_cell; k++)
				{
					int node_id=IEN[cell_id][k]-1;
					node_flag[node_id]=0;
				}
			}
			//Order the node_c_node
			std::sort(node_c_node.begin(),node_c_node.end());
			//Calulate the id of freedom degree connected to this node
			int num_connected_node=node_c_node.size();
			for (int j = 0; j < num_connected_node; j++)
			{
				int node_id=node_c_node[j]-1;
				for (int k = 0; k < num_dimension; k++)
				{
					if (ID[num_dimension*node_id+k]>0)
					{
						node_c_freedom.push_back(ID[num_dimension*node_id+k]);
					}
				}
			}
			//Calculate the freedom_c_freedom of the freedom on this node
			for (int j = 0; j < num_dimension; j++)
			{
				int num_connected_freedom=node_c_freedom.size();
				if (ID[num_dimension*i+j]>0)
				{
					for (int k = 0; k < num_connected_freedom; k++)
					{
						if (node_c_freedom[k]>=ID[num_dimension*i+j])
						{
							freedom_c_freedom[ID[num_dimension*i+j]-1].push_back(node_c_freedom[k]);
						}
					}
				}
			}
		}
		node_c_node.clear();
		node_c_freedom.clear();
	}
	//Calcualte ia and number of non-zero element
	int non_zero=0;
	for (int i = 0; i < num_freedom; i++)
	{
		ia[i+1]=ia[i]+freedom_c_freedom[i].size();
		non_zero=non_zero+freedom_c_freedom[i].size();
	}
	//Calculate ja
	vector<int> test;
	test.resize(non_zero);
	ja=new int[non_zero];
	non_zero=0;
	for (int i = 0; i < num_freedom; i++)
	{
		int num_freedom_connected_freedom=freedom_c_freedom[i].size();
		for (int j = 0; j < num_freedom_connected_freedom; j++)
		{
			ja[non_zero+j]=freedom_c_freedom[i][j];
			test[non_zero+j]=freedom_c_freedom[i][j];
		}
		non_zero=non_zero+num_freedom_connected_freedom;
	}
	//Return the non-zero elements in the matrix
	num_non_zero=non_zero;
	return;
}
void calculate_matrix_structure_KB(vector<int> &ID, vector<vector<int>> &IEN, int num_ID, int num_dimension, int* &ia, int* &ja, int &num_non_zero)
{
	int nump = ID.size() / num_dimension, nume = IEN.size();
	vector<int> node_flag;            //The node flag
	vector<int> node_c_node;
	//vector<int> node_c_freedom;
	vector<vector<int>> constraint_c_disID;
	ia = new int[num_ID + 1];
	ia[0] = 1;
	constraint_c_disID.resize(num_ID);

	//Initialize the node_flag
	node_flag.resize(nump);
	for (int i = 0; i < nump; i++)
	{
		node_flag[i] = 0;
	}

	//Calculate the cells connected to a node
	vector<vector<int>> NE;
	NE.resize(nump);
	for (int i = 0; i < nume; i++)
	{
		int num_node_in_cell = IEN[i].size();
		for (int j = 0; j < num_node_in_cell; j++)
		{
			int node_id = IEN[i][j] - 1;
			NE[node_id].push_back(i);
		}
	}

	//Calculate constraint_c_freedom
	for (int i = 0; i < nump; i++)
	{
		if (have_fixed(ID, i, num_dimension))
		{
			//If there is a constraint at the node, calculate the displacement ID connected to this constraint

			//Calculate the node_c_node for node i
			int num_connected_cell = NE[i].size();
			for (int j = 0; j < num_connected_cell; j++)
			{
				int cell_id = NE[i][j];
				int num_node_in_cell = IEN[cell_id].size();
				for (int k = 0; k < num_node_in_cell; k++)
				{
					int node_id = IEN[cell_id][k] - 1;
					if (node_flag[node_id] == 0)
					{
						node_c_node.push_back(node_id + 1);
						node_flag[node_id] = 1;
					}
				}
			}
			//Clear the node_flag after obtain the node_c_node
			for (int j = 0; j < num_connected_cell; j++)
			{
				int cell_id = NE[i][j];
				int num_node_in_cell = IEN[cell_id].size();
				for (int k = 0; k < num_node_in_cell; k++)
				{
					int node_id = IEN[cell_id][k] - 1;
					node_flag[node_id] = 0;
				}
			}
			//Order the node_c_node
			std::sort(node_c_node.begin(), node_c_node.end());
			int num_connected_node = node_c_node.size();
			for (int j = 0; j < num_dimension; j++)
			{
				//int num_connected_freedom = node_c_freedom.size();
				if (ID[num_dimension*i + j] < 0)
				{
					for (int k = 0; k < num_connected_node; k++)
					{
						int node_id = node_c_node[k] - 1;
						constraint_c_disID[-ID[num_dimension*i + j] - 1].push_back(num_dimension*node_id+1);
						constraint_c_disID[-ID[num_dimension*i + j] - 1].push_back(num_dimension*node_id+2);
						constraint_c_disID[-ID[num_dimension*i + j] - 1].push_back(num_dimension*node_id+3);
					}
				}
			}
			////Calculate the displacement id connected to the fixed id
			//int num_connected_node = node_c_node.size();
			////Calulate the id of freedom degree connected to this node
			//int num_connected_node = node_c_node.size();
			//for (int j = 0; j < num_connected_node; j++)
			//{
			//	int node_id = node_c_node[j] - 1;
			//	for (int k = 0; k < num_dimension; k++)
			//	{
			//		if (ID[num_dimension*node_id + k] > 0)
			//		{
			//			node_c_freedom.push_back(ID[num_dimension*node_id + k]);
			//		}
			//	}
			//}
			////Calculate the constraint_c_freedom of the freedom on this node
			//for (int j = 0; j < num_dimension; j++)
			//{
			//	int num_connected_freedom = node_c_freedom.size();
			//	if (ID[num_dimension*i + j] < 0)
			//	{
			//		for (int k = 0; k < num_connected_freedom; k++)
			//		{
			//			if (node_c_freedom[k] >= ID[num_dimension*i + j])
			//			{
			//				freedom_c_freedom[ID[num_dimension*i + j] - 1].push_back(node_c_freedom[k]);
			//			}
			//		}
			//	}
			//}
		}
		node_c_node.clear();
		//node_c_freedom.clear();
	}
	//Calcualte ia and number of non-zero element
	int non_zero = 0;
	for (int i = 0; i < num_ID; i++)
	{
		ia[i + 1] = ia[i] + constraint_c_disID[i].size();
		non_zero = non_zero + constraint_c_disID[i].size();
	}
	//Calculate ja
	//vector<int> test;
	//test.resize(non_zero);
	ja = new int[non_zero];
	non_zero = 0;
	for (int i = 0; i < num_ID; i++)
	{
		int num_constraint_connected_disID = constraint_c_disID[i].size();
		for (int j = 0; j < num_constraint_connected_disID; j++)
		{
			ja[non_zero + j] = constraint_c_disID[i][j];
			//test[non_zero + j] = freedom_c_freedom[i][j];
		}
		non_zero = non_zero + num_constraint_connected_disID;
	}
	//Return the non-zero elements in the matrix
	num_non_zero = non_zero;
	return;
}
void assmeble_stifness(int i,int j,int* ia,int* ja,double* a,double aij)
{
	//if (i>j)
	//{
	//	cout<<"Error:only aij with i>=j need to assemble in calling assmeble_stifness()"<<endl;
	//	system("Pause");
	//	exit(0);
	//}
	//Calcualte the position of the first element in i row
	int row_begin=ia[i-1];
	//Calculate the position of the first in i+1 row
	int row_begin_next=ia[i];
	//Calculate the position of aij in ja;
	int pos_aij=-1;
	for (int n = row_begin-1; n < row_begin_next-1; n++)
	{
		if (ja[n]==j)
		{
			pos_aij=n;
			break;
		}
	}
	if (pos_aij==-1)
	{
		cout<<"Error: aij is not in ja in calling assmeble_stifness"<<endl;
		system("Pause");
		exit(0);
	}
	//Assemble aij into a
	a[pos_aij]=a[pos_aij]+aij;
	return;
}
bool have_freedom(vector<int> &ID,int id,int num_dimension)
{
	if (num_dimension==3)
	{
		if (ID[3*id]>0 || ID[3*id+1]>0 || ID[3*id+2]>0)
		{
			return true;
		}
	}
	else if (num_dimension==2)
	{
		if (ID[2*id]>0 || ID[2*id+1]>0)
		{
			return true;
		}
	}
	else if (num_dimension==1)
	{
		if (ID[id]>0)
		{
			return true;
		}
	}
	else
	{
		cout<<"Error:Invalid dimension in calling have_freedom()"<<endl;
	}
	return false;
}
bool have_fixed(vector<int> &ID, int id, int num_dimension)
{
	if (num_dimension == 3)
	{
		if (ID[3 * id] < 0 || ID[3 * id + 1] < 0 || ID[3 * id + 2] < 0)
		{
			return true;
		}
	}
	else if (num_dimension == 2)
	{
		if (ID[2 * id] < 0 || ID[2 * id + 1] < 0)
		{
			return true;
		}
	}
	else if (num_dimension == 1)
	{
		if (ID[id] < 0)
		{
			return true;
		}
	}
	else
	{
		cout << "Error:Invalid dimension in calling have_freedom()" << endl;
	}
	return false;
}