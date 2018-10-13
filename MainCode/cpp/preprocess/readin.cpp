#include "readin.h"
#include <string>
#include <iostream>
#include "public_function.h"
#include "Tedge_viscosity.h"
#include "TFEM_viscosity.h"
#include "TMATshell_elastic.h"
#include "TMATshell_plastic.h"
#include "TMATIsotropic_plastic.h"
#include "TMAT_fluid.h"
#include "TMAT_Johnson_Cook.h"
using namespace std;
void Skeyword::clear_keyword()
{
	node_list.resize(0);
	particle_list.resize(0);
	cell_4_list.resize(0);
	cell_8_list.resize(0);
	boundary_list.resize(0);
	node_group_list.resize(0);
	load_list.resize(0);
	part_list.resize(0);
	section_shell_list.resize(0);
	material_list.resize(0);
	EOS_list.resize(0);
	curve_list.resize(0);
	//time_control;
	//background;
	initial_velocity_list.resize(0);
	MPM_boundary_condition_list.resize(0);
	MPM_background_load_list.resize(0);
	MPM_particle_load_list.resize(0);
	ALEMPM_background_property.resize(0);
}
string Sinteraction::calculate_interaction_type_id()
{
	//Return the interaction type id,order by programming time
	if (interaction_name=="conode")
	{
		return "00";
	}
	else
	{
		cout<<"Error: Invalid interaction type in calling calculate_interaction_type_id()"<<endl;
		system("Pause");
		exit(0);
	}
}
void Sremapping_scheme::exam_remapping_scheme()
{
	if (name!="Times" && name!="Adaption" && name!="None" && name!="DtReduction")
	{
		cout<<"Error:Invalid remapping scheme type"<<endl;
		system("Pause");
		exit(0);
		return;
	}
	else if (name=="Times")
	{
		if (num_remapping<=0)
		{
			cout<<"Error:The number of remapping times is "<<num_remapping<<endl;
			system("Pause");
			exit(0);
			return;
		}
	}
	else if (name=="Adaption")
	{
		if (tolerance_angle<=0)
		{
			cout<<"Error:The tolerance_angle in remapping scheme is "<<tolerance_angle<<endl;
			system("Pause");
			exit(0);
			return;
		}
		if (tolerance_length<=0)
		{
			cout<<"Error:The tolerance_length in remapping scheme is "<<tolerance_length<<endl;
			system("Pause");
			exit(0);
			return;
		}
	}
	return;
}
void new_line(ifstream& input)
{
	char a;
	do
	{
		input.get(a);
	} while (a!='\n');
}
bool exam_keyword(ifstream& input)
{
	char a;
	new_line(input);
	do
	{
		input.get(a);
		if (a=='*')
		{
			input.putback(a);
			return true;
		}
		else if (a=='$')
		{
			new_line(input);
		}
		else
		{
			input.putback(a);
			return false;
		}
	} while (a=='$');
	return 0;
}
void next_keyword(ifstream& input)
{
	char a;
	do
	{
		input.get(a);
		if (a=='*')
		{
			input.putback(a);
			return;
		}
		else  if (a!='\n')
		{
			//If a=='\n' the input stream will automatically come to the new line
			new_line(input);
		}
	} while (true);
	return;
}
void next_data(ifstream& input)
{
	new_line(input);
	char a;
	do
	{
		input.get(a);
		if (a=='$')
		{
			new_line(input);
		}
		else
		{
			input.putback(a);
			return;
		}
	} while (true);
}
TEOS_base* generate_EOS(SEOS EOS,Smaterial mat)
{
	TEOS_base* EOS_ptr;
	if (EOS.EOS_type=="*EOS_LINEAR_POLYNOMIAL")
	{
		TEOS_linear_polynomial* new_EOS;
		new_EOS=new TEOS_linear_polynomial(EOS.c0,EOS.c1,EOS.c2,EOS.c3,EOS.c4,EOS.c5,EOS.c6,mat.density);
		EOS_ptr=new_EOS;
	}
	else if (EOS.EOS_type=="*EOS_POLYTROPIC")
	{
		TEOS_polytropic* new_EOS;
		new_EOS=new TEOS_polytropic(EOS.p0,mat.density,EOS.gama,EOS.b);
		EOS_ptr=new_EOS;
	}
	else if (EOS.EOS_type=="*EOS_WEAK_IMPRESSIBLE")
	{
		TEOS_weak_impressible* new_EOS;
		new_EOS=new TEOS_weak_impressible(EOS.c0,mat.density,EOS.c1);
		EOS_ptr=new_EOS;
	}
	else if (EOS.EOS_type=="*EOS_JWL")
	{
		TEOS_JWL* new_EOS;
		new_EOS=new TEOS_JWL(EOS.A,EOS.B,EOS.R1,EOS.R2,EOS.w,mat.density);
		EOS_ptr=new_EOS;
	}
	else if (EOS.EOS_type=="*EOS_MieGruneisen")
	{
		TEOS_MieGruneisen* new_EOS;
		new_EOS=new TEOS_MieGruneisen(EOS.c0,EOS.s,EOS.gama0,mat.density);
		EOS_ptr=new_EOS;
	}
	else if (EOS.EOS_type=="*EOS_PSUEDO")
	{
		TEOS_PSEUDO* new_EOS;
		new_EOS=new TEOS_PSEUDO(EOS.c0,EOS.c1,EOS.c2,EOS.c3,EOS.c4,EOS.c5,EOS.c6,mat.density);
		EOS_ptr=new_EOS;
	}
	else
	{
		EOS_ptr=NULL;
	}
	return EOS_ptr;
}
TMAT_base* generate_material(Smaterial mat,SEOS EOS,int type)
{
	TMAT_base* mat_ptr;
	if (mat.material_type=="*MAT_ELASTIC")
	{
		if (type==1)
		{
			TMATshell_elastic* new_mat;
			new_mat=new TMATshell_elastic(mat.density,mat.Youngs,mat.Possion);
			mat_ptr=new_mat;
		}
		else if (type==0)
		{
			TMATIsotropic_elastic* new_mat;
			new_mat=new TMATIsotropic_elastic(mat.density,mat.Youngs,mat.Possion);
			mat_ptr=new_mat;
		}
	}
	else if (mat.material_type=="*MAT_PLASTIC_KINEMATIC")
	{
		if (type==1)
		{
			TMATshell_plastic* new_mat;
			new_mat=new TMATshell_plastic(mat.density,mat.Youngs,mat.Possion,mat.Sigma_y,mat.ET);
			mat_ptr=new_mat;
		}
		else if (type==0)
		{
			TMATIsotropic_plastic* new_mat;
			new_mat=new TMATIsotropic_plastic(mat.density,mat.Youngs,mat.Possion,mat.Sigma_y,mat.ET);
			mat_ptr=new_mat;
		}
	}
	else if (mat.material_type=="*MAT_Johnson_Cook")
	{
		TEOS_base* EOS_ptr;
		EOS_ptr=generate_EOS(EOS,mat);
		TMAT_Johnson_Cook* new_mat;
		new_mat=new TMAT_Johnson_Cook(mat.density,mat.Youngs,mat.Possion,mat.A,mat.B,mat.C,mat.n,mat.epso,EOS_ptr);
		mat_ptr=new_mat;
	}
	else if (mat.material_type=="*MAT_NULL")
	{
		TEOS_base* EOS_ptr;
		EOS_ptr=generate_EOS(EOS,mat);
		TMAT_fluid* new_mat;
		double internal_energy=maxval(EOS.internal_energy_per_volume/mat.density,1e-5);
		new_mat=new TMAT_fluid(mat.density,internal_energy,EOS_ptr);
		mat_ptr=new_mat;
	}
	return mat_ptr;
}
Tviscosity_base* generate_AV(Sarti_vis arti_vis)
{
	Tviscosity_base* av_ptr;
	if (arti_vis.av_type==0)
	{
		Tsolid_viscosity* new_av;
		new_av=new Tsolid_viscosity("Solid_viscosity",arti_vis.k1,arti_vis.k2);
		av_ptr=new_av;
	}
	else if (arti_vis.av_type==1)
	{
		Tedge_viscosity* new_av;
		new_av=new Tedge_viscosity("Edge_viscosity",arti_vis.k1,arti_vis.k2);
		av_ptr=new_av;
	}
	else if (arti_vis.av_type==2)
	{
		TFEM_viscosity* new_av;
		new_av=new TFEM_viscosity("FEM_viscosity",arti_vis.k1);
		av_ptr=new_av;
	}
	return av_ptr;
}

//-------------------------------------------------------------------
//Read in keyword file
//-------------------------------------------------------------------
void read_in_keyword_file(ifstream& input,Skeyword& keyword)
{
	string a;
	string temp;
	while (true)
	{
		next_keyword(input);
		input>>a;
		cout<<a<<endl;
		if (a=="*NODE")
		{
			while (!exam_keyword(input))
			{
				Snode_temp new_node;
				input>>new_node.id>>new_node.x>>new_node.y>>new_node.z;
				keyword.node_list.push_back(new_node);
			}
		}
		else if (a=="*MPM_PARTICLE")
		{
			next_data(input);
			int group_id,part_id;
			input>>group_id>>part_id;
			while (!exam_keyword(input))
			{
				Sparticle_temp new_particle;
				new_particle.group_id=group_id;
				new_particle.part_id=part_id;
				input>>new_particle.position.x>>new_particle.position.y;
				input>>new_particle.position.z>>new_particle.volume;
				keyword.particle_list.push_back(new_particle);
			}
		}
		else if (a=="*CONTROL_TERMINATION")
		{
			next_data(input);
			input>>keyword.time_control.endtime;
		}
		else if (a=="*CONTROL_TIMESTEP")
		{
			next_data(input);
			input>>temp>>keyword.time_control.CFL;
		}
		else if (a=="*BOUNDARY_SPC_SET")
		{
			next_data(input);
			Sboundary_type new_bd;
			input>>new_bd.id>>temp>>new_bd.pos[0]>>new_bd.pos[1]>>new_bd.pos[2];
			input>>new_bd.rot[0]>>new_bd.rot[1]>>new_bd.rot[2];
			new_bd.type=0;
			keyword.boundary_list.push_back(new_bd);
		}
		else if (a=="*VELOCITY_BOUNDARY")
		{
			next_data(input);
			Sboundary_type new_bd;
			input>>new_bd.id>>new_bd.dofvel[0]>>new_bd.dofvel[1]>>new_bd.dofvel[2]>>new_bd.vel[0]>>new_bd.vel[1]>>new_bd.vel[2];
			new_bd.type=1;
			keyword.boundary_list.push_back(new_bd);
		}
		else if (a=="*VELOCITY_BOUNDARY_CYLINDER")
		{
			next_data(input);
			Sboundary_type new_bd;
			input>>new_bd.id>>new_bd.x>>new_bd.y>>new_bd.v; 
			new_bd.type=2;
			keyword.boundary_list.push_back(new_bd);
		}
		else if (a=="*VELOCITY_BOUNDARY_SPHERE")
		{
			next_data(input);
			Sboundary_type new_bd;
			input>>new_bd.id>>new_bd.x>>new_bd.y>>new_bd.z>>new_bd.v; 
			new_bd.type=3;
			keyword.boundary_list.push_back(new_bd);
		}
		else if (a=="*NODE_PERTURBATION")
		{
			next_data(input);
			SNodePerturbation new_per;
			input>>new_per.id>>new_per.dimension>>new_per.x>>new_per.y>>new_per.z>>new_per.ratio; 
			keyword.im_node_perturbation=new_per;
		}
		else if (a=="*IMPLOSION")
		{
			next_data(input);
			keyword.is_implosion=true;
			input>>keyword.im.d0>>keyword.im.N>>keyword.im.r_in>>keyword.im.r_out;
		}
		else if (a=="*SET_NODE_LIST")
		{
			next_data(input);
			Snode_group new_group;
			input>>new_group.group_id;
			while (!exam_keyword(input))
			{
				int node_id;
				for (int i = 0; i < 8; i++)
				{
					input>>node_id;
					if (node_id>0)
					{
						new_group.node_id.push_back(node_id);
					}
				}
			}
			//Add new_group to node_group_list and keep in order
			unsigned int group_id=new_group.group_id;
			if (group_id>keyword.node_group_list.size())
			{
				keyword.node_group_list.resize(group_id);
			}
			keyword.node_group_list[group_id-1]=new_group;
		}
		else if (a=="*LOAD_NODE_SET")
		{
			while (!exam_keyword(input))
			{
				Sload_on_node new_load;
				input>>new_load.id>>new_load.dof>>new_load.curve_id;
				keyword.load_list.push_back(new_load);
			}
		}
		else if (a=="*PART" || a=="*PART_PARTICLE")
		{
			Spart new_part;
			next_data(input);
			next_data(input);
			input>>new_part.part_id>>new_part.section_id>>new_part.material_id>>new_part.EOS_id;
			//Add new_part to part_list and keep in order
			unsigned int part_id=new_part.part_id;
			if (part_id>keyword.part_list.size())
			{
				keyword.part_list.resize(part_id);
			}
			keyword.part_list[part_id-1]=new_part;
			if (a=="*PART_PARTICLE")
			{
				keyword.num_particle_part=keyword.num_particle_part+1;
			}
		}
		else if (a=="*SECTION_SHELL")
		{
			Ssection_shell new_section;
			next_data(input);
			input>>new_section.section_id>>temp>>temp>>new_section.nGauss;
			next_data(input);
			input>>new_section.thickness;
			//Add new_section to section_shell_list and keep in order
			unsigned int section_id=new_section.section_id;
			if (section_id>keyword.section_shell_list.size())
			{
				keyword.section_shell_list.resize(section_id);
			}
			keyword.section_shell_list[section_id-1]=new_section;
		}
		else if (a=="*MAT_ELASTIC")
		{
			Smaterial new_material;
			next_data(input);
			input>>new_material.material_id;
			input>>new_material.density;
			input>>new_material.Youngs>>new_material.Possion;
			new_material.material_type=a;
			//Add new_material into material_list and keep in order
			unsigned int material_id=new_material.material_id;
			if (material_id>keyword.material_list.size())
			{
				keyword.material_list.resize(material_id);
			}
			keyword.material_list[material_id-1]=new_material;
		}
		else if (a=="*MAT_PLASTIC_KINEMATIC")
		{
			Smaterial new_material;
			next_data(input);
			input>>new_material.material_id;
			input>>new_material.density>>new_material.Youngs>>new_material.Possion;
			input>>new_material.Sigma_y>>new_material.ET;
			new_material.material_type=a;
			//Add new_material into material_list and keep in order
			unsigned int material_id=new_material.material_id;
			if (material_id>keyword.material_list.size())
			{
				keyword.material_list.resize(material_id);
			}
			keyword.material_list[material_id-1]=new_material;
		}
		else if (a=="*MAT_Johnson_Cook")
		{
			Smaterial new_material;
			next_data(input);
			input>>new_material.material_id;
			input>>new_material.density>>new_material.Youngs>>new_material.Possion;
			input>>new_material.A>>new_material.B>>new_material.n>>new_material.C>>new_material.epso;
			new_material.material_type=a;
			//Add new_material into material_list and keep in order
			unsigned int material_id=new_material.material_id;
			if (material_id>keyword.material_list.size())
			{
				keyword.material_list.resize(material_id);
			}
			keyword.material_list[material_id-1]=new_material;
		}
		else if (a=="*MAT_NULL")
		{
			Smaterial new_material;
			next_data(input);
			input>>new_material.material_id;
			input>>new_material.density;
			new_material.material_type=a;
			//Add new_material into material_list and keep in order
			unsigned int material_id=new_material.material_id;
			if (material_id>keyword.material_list.size())
			{
				keyword.material_list.resize(material_id);
			}
			keyword.material_list[material_id-1]=new_material;
		}
		else if (a=="*EOS_LINEAR_POLYNOMIAL")
		{
			SEOS new_EOS;
			next_data(input);
			input>>new_EOS.EOS_id;
			input>>new_EOS.c0>>new_EOS.c1>>new_EOS.c2>>new_EOS.c3;
			input>>new_EOS.c4>>new_EOS.c5>>new_EOS.c6;
			next_data(input);
			input>>new_EOS.internal_energy_per_volume;
			new_EOS.EOS_type=a;
			//Add new_EOS into EOS_list and keep in order
			unsigned int EOS_id=new_EOS.EOS_id;
			if (EOS_id>keyword.EOS_list.size())
			{
				keyword.EOS_list.resize(EOS_id);
			}
			keyword.EOS_list[EOS_id-1]=new_EOS;
		}
		else if (a=="*EOS_POLYTROPIC")
		{
			SEOS new_EOS;
			next_data(input);
			input>>new_EOS.EOS_id;
			input>>new_EOS.p0>>new_EOS.gama>>new_EOS.b>>new_EOS.internal_energy_per_volume;
			next_data(input);
			new_EOS.EOS_type=a;
			//Add new_EOS into EOS_list and keep in order
			unsigned int EOS_id=new_EOS.EOS_id;
			if (EOS_id>keyword.EOS_list.size())
			{
				keyword.EOS_list.resize(EOS_id);
			}
			keyword.EOS_list[EOS_id-1]=new_EOS;
		}
		else if (a=="*EOS_WEAK_IMPRESSIBLE")
		{
			SEOS new_EOS;
			next_data(input);
			input>>new_EOS.EOS_id>>new_EOS.c0>>new_EOS.c1>>new_EOS.internal_energy_per_volume;
			next_data(input);
			new_EOS.EOS_type=a;
			//Add new_EOS into EOS_list and keep in order
			unsigned int EOS_id=new_EOS.EOS_id;
			if (EOS_id>keyword.EOS_list.size())
			{
				keyword.EOS_list.resize(EOS_id);
			}
			keyword.EOS_list[EOS_id-1]=new_EOS;
		}
		else if (a=="*EOS_JWL")
		{
			SEOS new_EOS;
			next_data(input);
			input>>new_EOS.EOS_id>>new_EOS.A>>new_EOS.B>>new_EOS.R1>>new_EOS.R2>>new_EOS.w>>new_EOS.internal_energy_per_volume;
			next_data(input);
			new_EOS.EOS_type=a;
			//Add new_EOS into EOS_list and keep in order
			unsigned int EOS_id=new_EOS.EOS_id;
			if (EOS_id>keyword.EOS_list.size())
			{
				keyword.EOS_list.resize(EOS_id);
			}
			keyword.EOS_list[EOS_id-1]=new_EOS;
		}
		else if (a=="*EOS_MieGruneisen")
		{
			SEOS new_EOS;
			next_data(input);
			input>>new_EOS.EOS_id>>new_EOS.c0>>new_EOS.s>>new_EOS.gama0;
			next_data(input);
			new_EOS.EOS_type=a;
			unsigned int EOS_id=new_EOS.EOS_id;
			if (EOS_id>keyword.EOS_list.size())
			{
				keyword.EOS_list.resize(EOS_id);
			}
			keyword.EOS_list[EOS_id-1]=new_EOS;
		}
		else if (a=="*EOS_PESUDO")
		{
			SEOS new_EOS;
			next_data(input);
			input>>new_EOS.EOS_id;
			input>>new_EOS.c0>>new_EOS.c1>>new_EOS.c2>>new_EOS.c3;
			input>>new_EOS.c4>>new_EOS.c5>>new_EOS.c6;
			next_data(input);
			input>>new_EOS.internal_energy_per_volume;
			new_EOS.EOS_type=a;
			//Add new_EOS into EOS_list and keep in order
			unsigned int EOS_id=new_EOS.EOS_id;
			if (EOS_id>keyword.EOS_list.size())
			{
				keyword.EOS_list.resize(EOS_id);
			}
			keyword.EOS_list[EOS_id-1]=new_EOS;
		}
		else if (a=="*DEFINE_CURVE")
		{
			double v1,v2;
			Scurve new_curve;
			next_data(input);
			input>>new_curve.curve_id;
			//new_line(input);
			while (!exam_keyword(input))
			{
				input>>v1>>v2;
				//new_line(input);
				new_curve.v1.push_back(v1);
				new_curve.v2.push_back(v2);
			}
			//Add new_curve to curve_list and keep in order
			unsigned int curve_id=new_curve.curve_id;
			if (curve_id>keyword.curve_list.size())
			{
				keyword.curve_list.resize(curve_id);
			}
			keyword.curve_list[curve_id-1]=new_curve;
		}
		else if (a=="*ELEMENT_SHELL")
		{
			while (!exam_keyword(input))
			{
				Scell_4 new_shell;
				input>>new_shell.cell_id>>new_shell.part_id;
				input>>new_shell.IEN[0]>>new_shell.IEN[1];
				input>>new_shell.IEN[2]>>new_shell.IEN[3];
				keyword.cell_4_list.push_back(new_shell);
			}
		}
		else if (a=="*ELEMENT_SOLID")
		{
			while (!exam_keyword(input))
			{
				Scell_8 new_shell;
				input>>new_shell.cell_id>>new_shell.part_id;
				input>>new_shell.IEN[0]>>new_shell.IEN[1]>>new_shell.IEN[2]>>new_shell.IEN[3];
				input>>new_shell.IEN[4]>>new_shell.IEN[5]>>new_shell.IEN[6]>>new_shell.IEN[7];
				keyword.cell_8_list.push_back(new_shell);
			}
		}
		else if (a=="*ELEMENT_SOLID_TET")
		{
			bool is_first=true;
			while (!exam_keyword(input))
			{
				Scell_4 new_tet;
				input>>new_tet.cell_id>>new_tet.part_id;
				if (is_first)
				{
					next_data(input);
					is_first=false;
				}
				input>>new_tet.IEN[0]>>new_tet.IEN[1]>>new_tet.IEN[2]>>new_tet.IEN[3];
				keyword.cell_4_list.push_back(new_tet);
			}
		}
		else if (a=="*BACKGROUND")
		{
			next_data(input);
			input>>keyword.background.x_min.x>>keyword.background.x_min.y>>keyword.background.x_min.z;
			next_data(input);
			input>>keyword.background.x_max.x>>keyword.background.x_max.y>>keyword.background.x_max.z;
			next_data(input);
			input>>keyword.background.nx>>keyword.background.ny>>keyword.background.nz;
		}
		else if (a=="*INITIAL_VELOCITY")
		{
			next_data(input);
			Sinitial_velocity_temp new_initial_velocity;
			input>>new_initial_velocity.part_id>>new_initial_velocity.velocity.x>>new_initial_velocity.velocity.y;
			input>>new_initial_velocity.velocity.z;
			keyword.initial_velocity_list.push_back(new_initial_velocity);
		}
		else if (a=="*MPM_BOUNDARY_CONDITION")
		{
			next_data(input);
			SMPM_boundary_condrition new_boundary_condition;
			input>>new_boundary_condition.strat_subscript[0]>>new_boundary_condition.strat_subscript[1]>>new_boundary_condition.strat_subscript[2];
			next_data(input);
			input>>new_boundary_condition.terminate_subscript[0]>>new_boundary_condition.terminate_subscript[1]>>new_boundary_condition.terminate_subscript[2];
			next_data(input);
			input>>new_boundary_condition.constraint[0]>>new_boundary_condition.constraint[1]>>new_boundary_condition.constraint[2];
			keyword.MPM_boundary_condition_list.push_back(new_boundary_condition);
		}
		else if (a=="*MPM_BOUNDARY_LOAD_BACKGROUND")
		{
			next_data(input);
			SMPM_background_load new_load;
			input>>new_load.strat_subscript[0]>>new_load.strat_subscript[1]>>new_load.strat_subscript[2];
			next_data(input);
			input>>new_load.terminate_subscript[0]>>new_load.terminate_subscript[1]>>new_load.terminate_subscript[2];
			next_data(input);
			input>>new_load.load.x>>new_load.load.y>>new_load.load.z;
			keyword.MPM_background_load_list.push_back(new_load);
		}
		else if (a=="*MPM_PARTICLE_LOAD")
		{
			next_data(input);
			SMPM_particle_load new_load;
			input>>new_load.part_id>>new_load.particle_id;
			next_data(input);
			input>>new_load.load.x>>new_load.load.y>>new_load.load.z;
			next_data(input);
			keyword.MPM_particle_load_list.push_back(new_load);
		}
		else if (a=="*ALEMPM_BACKGROUND_PROPERTY")
		{
			next_data(input);
			Sbackground_property new_property;
			input>>new_property.part_id;
			next_data(input);
			input>>new_property.nx_begin>>new_property.ny_begin>>new_property.nz_begin;
			next_data(input);
			input>>new_property.nx_end>>new_property.ny_end>>new_property.nz_end;
			keyword.ALEMPM_background_property.push_back(new_property);
		}
		else if (a=="*HOURGLASS_OPTION")
		{
			next_data(input);
			input>>keyword.hourglass_option;
		}
		else if (a=="*END")
		{
			return;
		}
	}
	return;
}