#include "Tcell_CNT.h"
Tcell_CNT::Tcell_CNT(int cell_id, Tnode_CNT* (&node_ptr)[2], double area, TMAT_base* mat)
{
	_id = cell_id;
	_node_ptr[0] = node_ptr[0]; _node_ptr[1] = node_ptr[1];
	_area = area;
	_mat_ptr = mat;
	_l0 = (node_ptr[0]->_position - node_ptr[1]->_position).get_length();
}
void Tcell_CNT::calculate_corner_force()
{
	double l= (_node_ptr[0]->_position - _node_ptr[1]->_position).get_length();
	double force = (l - _l0)*_mat_ptr->G_Youngs()*_area/_l0;
	vec3D direction = (_node_ptr[0]->_position - _node_ptr[1]->_position) / l;
	_corner_force[0] = -direction * force; _corner_force[1] = direction * force;
	return;
}
void Tcell_CNT::assemble_corner_force()
{
	_node_ptr[0]->_force = _node_ptr[0]->_force + _corner_force[0];
	_node_ptr[1]->_force = _node_ptr[1]->_force + _corner_force[1];
}
