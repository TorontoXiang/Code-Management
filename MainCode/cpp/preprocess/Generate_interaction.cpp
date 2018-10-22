#include "readin.h"
#include "Tinteraction_conode.h"

template<class T1,class T2> Tinteraction_base* generate_interaction_with_type(Sinteraction interaction)
{
	Tinteraction_base* interaction_ptr=NULL;
	T1 body1_ptr=dynamic_cast<T1>(interaction.body1);
	T2 body2_ptr=dynamic_cast<T2>(interaction.body2);
	if (interaction.interaction_name=="conode")
	{
		Tinteraction_conode<T1,T2>* new_interaction;
		new_interaction=new Tinteraction_conode<T1,T2>(interaction.interaction_name,body1_ptr,body2_ptr);
		interaction_ptr=new_interaction;
	}
	else
	{
		cout<<"Error:Invalid interaction_name in calling generate_interaction_with_type"<<endl;
		system("Pause");
		exit(0);
	}
	return interaction_ptr;
}
Tinteraction_base* generate_interaction(Sinteraction interaction)
{
	Tinteraction_base* interaction_ptr=NULL;

	string interaction_type_id=interaction.calculate_interaction_type_id();
	string body_type_id1=interaction.body1->calculate_body_type_id();
	string body_type_id2=interaction.body2->calculate_body_type_id();
	string index=body_type_id1+"-"+body_type_id2+"-"+interaction_type_id;
	if (index=="00-00-00")         //00-00-00:shell-shell-conode
	{
		typedef Tbody_shell* T1;typedef Tbody_shell* T2;
		return generate_interaction_with_type<T1,T2>(interaction);
	}
	else if (index=="01-01-00")    //01-01-00:Lagrangian-Lagrangian-conode
	{
		typedef Tbody_Lagrangian* T1;typedef Tbody_Lagrangian* T2;
		return generate_interaction_with_type<T1,T2>(interaction);
	}
	else if (index=="00-02-00")    //00-02-00:shell-ALE-conode
	{
		typedef Tbody_shell* T1;typedef Tbody_ALE* T2;
		return generate_interaction_with_type<T1,T2>(interaction);
	}
	else if (index=="02-00-00")    //02-00-00:ALE-shell-conode
	{
		typedef Tbody_ALE* T1;typedef Tbody_shell* T2;
		return generate_interaction_with_type<T1,T2>(interaction);
	}
	else if (index=="02-02-00")    //04-00-00:ALE-ALE-conode
	{
		typedef Tbody_ALE* T1;typedef Tbody_ALE* T2;
		return generate_interaction_with_type<T1,T2>(interaction);
	}
	else if (index=="02-03-00")
	{
		typedef Tbody_ALE* T1; typedef Tbody_brick* T2;
		return generate_interaction_with_type<T1, T2>(interaction);
	}
	else
	{
		cout<<"Error: Invalid interaction group: "<<index<<endl;
		cout<<"Body type id:"<<endl;
		cout<<"00-shell;01-pure Lagrangian;02-pure ALE"<<endl;
		cout<<"Interaction type id:"<<endl;
		cout<<"00-conode interaction"<<endl;
		system("Pause");
		exit(0);
	}
}
