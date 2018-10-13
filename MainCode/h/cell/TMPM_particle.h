#ifndef TMPM_PARTICLE
#define TMPM_PARTICLE
//Define the class of particle in MPM method
#include "data_structure.h"
#include "TMAT_base.h"
class TMPM_particle
{
public:
	TMPM_particle(vec3D position,double volume,int group_id,int part_id,int id_in_part):
		_position(position),_volume(volume),_mass(0),_internal_energy(0),_part_id(part_id),
		_cell_index(-1),_gausspoint(NULL),_group_id(group_id),_id_in_part(id_in_part){};
protected:
	int _group_id;
	int _part_id;
	int _id_in_part;
	vec3D _position;
	vec3D _velocity;
	vec3D _moment;
	double _volume;
	double _mass;
	double _internal_energy;
	TMAT_base* _gausspoint;
	int _cell_index;             //The background cell that this particle in
	vec3D _standard_coordinate;  //The standard coordinate in the background cell
	vec3D _external_force;       //The external force applied on the particle
	double _de[6];
	double _detJ;


	friend class Tbody_MPM;
	friend class TMPM_background;
	friend class Tcell_MPM;
	friend class Tbody_ALEMPM;
	friend class Tcell_pure_fluid_MPM;
	friend class Tcell_mixed_fluid_MPM;
	friend class Tbody_explosion;
};
#endif