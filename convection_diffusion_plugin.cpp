/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/**
 * File for registration of ConvectionDiffusion routines.
 */

#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "convection_diffusion_base.h"
#include "convection_diffusion_sss.h"
#include "fv1/convection_diffusion_fv1.h"
#include "mp/convection_diffusion_mp.h"
#include "fe/convection_diffusion_fe.h"
#include "fe/convection_diffusion_stab_fe.h"
#include "fvcr/convection_diffusion_fvcr.h"
#include "fv/convection_diffusion_fv.h"
#include "fractfv1/convection_diffusion_fractfv1.h"

#include "convection_diffusion_plugin.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace ConvectionDiffusionPlugin{

/** 
 *  \defgroup convection_diffusion Convection Diffusion
 *  \ingroup plugins
 *  This plugin provides the discretization of convection and diffusion problems.
 *  \{
 */

/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain dependent parts
 * of the plugin. All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TRegistry=ug::bridge::Registry>
static void Domain(TRegistry& reg, string grp)
{
	static const int dim = TDomain::dim;
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

//	Convection Diffusion Base
	{
		typedef ConvectionDiffusionBase<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("ConvectionDiffusionBase").append(suffix);
		reg.template add_class_<T, TBase >(name, grp)
			.add_method("set_diffusion", static_cast<void (T::*)(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> >)>(&T::set_diffusion), "", "Diffusion")
			.add_method("set_diffusion", static_cast<void (T::*)(number)>(&T::set_diffusion), "", "Diagonal Diffusion")
#ifdef UG_FOR_LUA
			.add_method("set_diffusion", static_cast<void (T::*)(const char*)>(&T::set_diffusion), "", "Diffusion")
			.add_method("set_diffusion", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_diffusion), "", "Diffusion")
#endif

			.add_method("set_diffusion_Sw", static_cast<void (T::*)(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> >)>(&T::set_diffusion_Sw), "", "Diffusion_Sw")
			.add_method("set_diffusion_Sw", static_cast<void (T::*)(number)>(&T::set_diffusion_Sw), "", "Diagonal Diffusion_Sw")
#ifdef UG_FOR_LUA
			.add_method("set_diffusion_Sw", static_cast<void (T::*)(const char*)>(&T::set_diffusion_Sw), "", "Diffusion_Sw")
			.add_method("set_diffusion_Sw", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_diffusion_Sw), "", "Diffusion_Sw")
#endif

			.add_method("set_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_velocity), "", "Velocity Field")
			.add_method("set_velocity", static_cast<void (T::*)(const std::vector<number>&)>(&T::set_velocity), "", "Velocity Field")
#ifdef UG_FOR_LUA
			.add_method("set_velocity", static_cast<void (T::*)(const char*)>(&T::set_velocity), "", "Velocity Field")
			.add_method("set_velocity", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_velocity), "", "Velocity Field")
#endif

			.add_method("set_darcyW", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_darcyW), "", "DarcyW Field")
			.add_method("set_darcyW", static_cast<void (T::*)(const std::vector<number>&)>(&T::set_darcyW), "", "DarcyW Field")
#ifdef UG_FOR_LUA
			.add_method("set_darcyW", static_cast<void (T::*)(const char*)>(&T::set_darcyW), "", "DarcyW Field")
			.add_method("set_darcyW", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_darcyW), "", "DarcyW Field")
#endif

			.add_method("set_darcyN", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_darcyN), "", "DarcyN Field")
			.add_method("set_darcyN", static_cast<void (T::*)(const std::vector<number>&)>(&T::set_darcyN), "", "DarcyN Field")
#ifdef UG_FOR_LUA
			.add_method("set_darcyN", static_cast<void (T::*)(const char*)>(&T::set_darcyN), "", "DarcyN Field")
			.add_method("set_darcyN", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_darcyN), "", "DarcyN Field")
#endif


			.add_method("set_saturationW", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_saturationW), "", "SaturationW")
			.add_method("set_saturationW", static_cast<void (T::*)(number)>(&T::set_saturationW), "", "SaturationW")
#ifdef UG_FOR_LUA
			.add_method("set_saturationW", static_cast<void (T::*)(const char*)>(&T::set_saturationW), "", "SaturationW")
			.add_method("set_saturationW", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_saturationW), "", "SaturationW")
#endif

			.add_method("set_upwindsaturationW", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_upwindsaturationW), "", "UpwindSaturationW")
			.add_method("set_upwindsaturationW", static_cast<void (T::*)(number)>(&T::set_upwindsaturationW), "", "UpwindSaturationW")
#ifdef UG_FOR_LUA
			.add_method("set_upwindsaturationW", static_cast<void (T::*)(const char*)>(&T::set_upwindsaturationW), "", "UpwindSaturationW")
			.add_method("set_upwindsaturationW", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_upwindsaturationW), "", "UpwindSaturationW")
#endif

			.add_method("set_modifiedsaturationW", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_modifiedsaturationW), "", "ModifiedSaturationW")
			.add_method("set_modifiedsaturationW", static_cast<void (T::*)(number)>(&T::set_modifiedsaturationW), "", "ModifiedSaturationW")
#ifdef UG_FOR_LUA
			.add_method("set_modifiedsaturationW", static_cast<void (T::*)(const char*)>(&T::set_modifiedsaturationW), "", "ModifiedSaturationW")
			.add_method("set_modifiedsaturationW", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_modifiedsaturationW), "", "ModifiedSaturationW")
#endif

			.add_method("set_upwindmodifiedsaturationW", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_upwindmodifiedsaturationW), "", "UpwindModifiedSaturationW")
			.add_method("set_upwindmodifiedsaturationW", static_cast<void (T::*)(number)>(&T::set_upwindmodifiedsaturationW), "", "UpwindModifiedSaturationW")
#ifdef UG_FOR_LUA
			.add_method("set_upwindmodifiedsaturationW", static_cast<void (T::*)(const char*)>(&T::set_upwindmodifiedsaturationW), "", "UpwindModifiedSaturationW")
			.add_method("set_upwindmodifiedsaturationW", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_upwindmodifiedsaturationW), "", "UpwindModifiedSaturationW")
#endif

			.add_method("set_MassFractionWc", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_MassFractionWc), "", "MassFractionWc")
			.add_method("set_MassFractionWc", static_cast<void (T::*)(number)>(&T::set_MassFractionWc), "", "MassFractionWc")
#ifdef UG_FOR_LUA
			.add_method("set_MassFractionWc", static_cast<void (T::*)(const char*)>(&T::set_MassFractionWc), "", "MassFractionWc")
			.add_method("set_MassFractionWc", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_MassFractionWc), "", "MassFractionWc")
#endif

			.add_method("set_PressurePn", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_PressurePn), "", "PressurePn")
			.add_method("set_PressurePn", static_cast<void (T::*)(number)>(&T::set_PressurePn), "", "PressurePn")
#ifdef UG_FOR_LUA
			.add_method("set_PressurePn", static_cast<void (T::*)(const char*)>(&T::set_PressurePn), "", "PressurePn")
			.add_method("set_PressurePn", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_PressurePn), "", "PressurePn")
#endif

			.add_method("set_permeability", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_permeability), "", "Permeability")
			.add_method("set_permeability", static_cast<void (T::*)(number)>(&T::set_permeability), "", "Permeability")
#ifdef UG_FOR_LUA
			.add_method("set_permeability", static_cast<void (T::*)(const char*)>(&T::set_permeability), "", "Permeability")
			.add_method("set_permeability", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_permeability), "", "Permeability")
#endif
			.add_method("set_porosity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_porosity), "", "Porosity")
			.add_method("set_porosity", static_cast<void (T::*)(number)>(&T::set_porosity), "", "Porosity")
#ifdef UG_FOR_LUA
			.add_method("set_porosity", static_cast<void (T::*)(const char*)>(&T::set_porosity), "", "Porosity")
			.add_method("set_porosity", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_porosity), "", "Porosity")
#endif
			.add_method("set_pd", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_pd), "", "EntryPressure")
			.add_method("set_pd", static_cast<void (T::*)(number)>(&T::set_pd), "", "EntryPressure")
#ifdef UG_FOR_LUA
			.add_method("set_pd", static_cast<void (T::*)(const char*)>(&T::set_pd), "", "EntryPressure")
			.add_method("set_pd", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_pd), "", "EntryPressure")
#endif

			.add_method("set_swr", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_swr), "", "ResidualAqueous")
			.add_method("set_swr", static_cast<void (T::*)(number)>(&T::set_swr), "", "ResidualAqueous")
#ifdef UG_FOR_LUA
			.add_method("set_swr", static_cast<void (T::*)(const char*)>(&T::set_swr), "", "ResidualAqueous")
			.add_method("set_swr", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_swr), "", "ResidualAqueous")
#endif
			.add_method("set_snr", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_snr), "", "ResidualCarbonic")
			.add_method("set_snr", static_cast<void (T::*)(number)>(&T::set_snr), "", "ResidualCarbonic")
#ifdef UG_FOR_LUA
			.add_method("set_snr", static_cast<void (T::*)(const char*)>(&T::set_snr), "", "ResidualCarbonic")
			.add_method("set_snr", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_snr), "", "ResidualCarbonic")
#endif
			.add_method("set_lambda", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_lambda), "", "BrooksCoreyNumber")
			.add_method("set_lambda", static_cast<void (T::*)(number)>(&T::set_lambda), "", "BrooksCoreyNumber")
#ifdef UG_FOR_LUA
			.add_method("set_lambda", static_cast<void (T::*)(const char*)>(&T::set_lambda), "", "BrooksCoreyNumber")
			.add_method("set_lambda", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_lambda), "", "BrooksCoreyNumber")
#endif

			.add_method("set_minPd", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_minPd), "", "MinPd")
			.add_method("set_minPd", static_cast<void (T::*)(number)>(&T::set_minPd), "", "MinPd")
#ifdef UG_FOR_LUA
			.add_method("set_minPd", static_cast<void (T::*)(const char*)>(&T::set_minPd), "", "MinPd")
			.add_method("set_minPd", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_minPd), "", "MinPd")
#endif

			.add_method("set_minSwr", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_minSwr), "", "MinSwr")
			.add_method("set_minSwr", static_cast<void (T::*)(number)>(&T::set_minSwr), "", "MinSwr")
#ifdef UG_FOR_LUA
			.add_method("set_minSwr", static_cast<void (T::*)(const char*)>(&T::set_minSwr), "", "MinSwr")
			.add_method("set_minSwr", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_minSwr), "", "MinSwr")
#endif

			.add_method("set_minSnr", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_minSnr), "", "MinSnr")
			.add_method("set_minSnr", static_cast<void (T::*)(number)>(&T::set_minSnr), "", "MinSnr")
#ifdef UG_FOR_LUA
			.add_method("set_minSnr", static_cast<void (T::*)(const char*)>(&T::set_minSnr), "", "MinSnr")
			.add_method("set_minSnr", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_minSnr), "", "MinSnr")
#endif

			.add_method("set_minLambda", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_minLambda), "", "MinLambda")
			.add_method("set_minLambda", static_cast<void (T::*)(number)>(&T::set_minLambda), "", "MinLambda")
#ifdef UG_FOR_LUA
			.add_method("set_minLambda", static_cast<void (T::*)(const char*)>(&T::set_minLambda), "", "MinLambda")
			.add_method("set_minLambda", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_minLambda), "", "MinLambda")
#endif


			.add_method("set_flux", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_flux), "", "Flux")
			.add_method("set_flux", static_cast<void (T::*)(const std::vector<number>&)>(&T::set_flux), "", "Flux")
#ifdef UG_FOR_LUA
			.add_method("set_flux", static_cast<void (T::*)(const char*)>(&T::set_flux), "", "Flux")
			.add_method("set_flux", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_flux), "", "Flux")
#endif

			.add_method("set_reaction_rate", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_reaction_rate), "", "Reaction Rate")
			.add_method("set_reaction_rate", static_cast<void (T::*)(number)>(&T::set_reaction_rate), "", "Reaction Rate")
#ifdef UG_FOR_LUA
			.add_method("set_reaction_rate", static_cast<void (T::*)(const char*)>(&T::set_reaction_rate), "", "Reaction Rate")
			.add_method("set_reaction_rate", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_reaction_rate), "", "Reaction Rate")
#endif

			.add_method("set_reaction", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_reaction), "", "Reaction")
			.add_method("set_reaction", static_cast<void (T::*)(number)>(&T::set_reaction), "", "Reaction")
#ifdef UG_FOR_LUA
			.add_method("set_reaction", static_cast<void (T::*)(const char*)>(&T::set_reaction), "", "Reaction")
			.add_method("set_reaction", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_reaction), "", "Reaction")
#endif

			.add_method("set_reaction_rate_explicit", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_reaction_rate_explicit), "", "Reaction Rate Explicit")
			.add_method("set_reaction_rate_explicit", static_cast<void (T::*)(number)>(&T::set_reaction_rate_explicit), "", "Reaction Rate Explicit")
#ifdef UG_FOR_LUA
			.add_method("set_reaction_rate_explicit", static_cast<void (T::*)(const char*)>(&T::set_reaction_rate_explicit), "", "Reaction Rate Explicit")
#endif

			.add_method("set_reaction_explicit", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_reaction_explicit), "", "Reaction Explicit")
			.add_method("set_reaction_explicit", static_cast<void (T::*)(number)>(&T::set_reaction_explicit), "", "Reaction Explicit")
#ifdef UG_FOR_LUA
			.add_method("set_reaction_explicit", static_cast<void (T::*)(const char*)>(&T::set_reaction_explicit), "", "Reaction Explicit")
#endif

			.add_method("set_source_explicit", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_source_explicit), "", "Source Explicit")
			.add_method("set_source_explicit", static_cast<void (T::*)(number)>(&T::set_source_explicit), "", "Source Explicit")
			#ifdef UG_FOR_LUA
			.add_method("set_source_explicit", static_cast<void (T::*)(const char*)>(&T::set_source_explicit), "", "Source Explicit")
			#endif

			.add_method("set_source", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_source), "", "Source")
			.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source), "", "Source")
#ifdef UG_FOR_LUA
			.add_method("set_source", static_cast<void (T::*)(const char*)>(&T::set_source), "", "Source")
			.add_method("set_source", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_source), "", "Source")
#endif

			.add_method("set_vector_source", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_vector_source), "", "Vector Source")
#ifdef UG_FOR_LUA
			.add_method("set_vector_source", static_cast<void (T::*)(const char*)>(&T::set_vector_source), "", "Vector Source")
			.add_method("set_vector_source", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_vector_source), "", "Vector Source")
#endif

			.add_method("set_mass_scale", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_mass_scale), "", "Mass Scale")
			.add_method("set_mass_scale", static_cast<void (T::*)(number)>(&T::set_mass_scale), "", "Mass Scale")
#ifdef UG_FOR_LUA
			.add_method("set_mass_scale", static_cast<void (T::*)(const char*)>(&T::set_mass_scale), "", "Mass Scale")
			.add_method("set_mass_scale", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_mass_scale), "", "Mass Scale")
#endif

			.add_method("set_mass", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_mass), "", "Mass")
			.add_method("set_mass", static_cast<void (T::*)(number)>(&T::set_mass), "", "Mass")
#ifdef UG_FOR_LUA
			.add_method("set_mass", static_cast<void (T::*)(const char*)>(&T::set_mass), "", "Mass")
			.add_method("set_mass", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_mass), "", "Mass")
#endif

			.add_method("upwindmodifiedvalue", &T::upwindmodifiedvalue)
			.add_method("upwindvalue", &T::upwindvalue)
			.add_method("modifiedvalue", &T::modifiedvalue)
			.add_method("value", &T::value)
		  .add_method("gradient", &T::gradient)
		  .add_method("gradient_pd", &T::gradient_pd)
		  .add_method("gradient_pd2", &T::gradient_pd2);
		  /*
			.add_method("set_partial_velocity", &T::set_partial_velocity)
			.add_method("set_partial_flux", &T::set_partial_flux)
			.add_method("set_partial_mass", &T::set_partial_mass);
		  */
		reg.add_class_to_group(name, "ConvectionDiffusionBase", tag);
	}

//	Convection Diffusion FV1
	{
		typedef ConvectionDiffusionFV1<TDomain> T;
		typedef ConvectionDiffusionBase<TDomain> TBase;
		string name = string("ConvectionDiffusionFV1").append(suffix);
		reg.template add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_condensed_FV", &T::set_condensed_FV, "", "[De-]Activates the condensed FV scvf ip's")
			.add_method("set_upwind", &T::set_upwind, "", "Sets the upwind type for the convective terms")
			.add_method("set_singular_sources_and_sinks", &T::set_sss_manager, "", "Sets the singular sources and sinks manager")
			.add_method("singular_sources_and_sinks", &T::sss_manager, "", "Returns the singular sources and sinks manager")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusionFV1", tag);
	}
	
//	Convection Diffusion MP
	{
		typedef ConvectionDiffusionMP<TDomain> T;
		typedef ConvectionDiffusionBase<TDomain> TBase;
		string name = string("ConvectionDiffusionMP").append(suffix);
		reg.template add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*,const int)>("Function(s)#Subset(s)#FormulationIndex")
			.add_method("set_condensed_FV", &T::set_condensed_FV, "", "[De-]Activates the condensed FV scvf ip's")
			.add_method("set_upwind", &T::set_upwind, "", "Sets the upwind type for the convective terms")
			.add_method("set_singular_sources_and_sinks", &T::set_sss_manager, "", "Sets the singular sources and sinks manager")
			.add_method("singular_sources_and_sinks", &T::sss_manager, "", "Returns the singular sources and sinks manager")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusionMP", tag);
	}

//	Convection Diffusion FE
	{
		typedef ConvectionDiffusionFE<TDomain> T;
		typedef ConvectionDiffusionBase<TDomain> TBase;
		string name = string("ConvectionDiffusionFE").append(suffix);
		reg.template add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_quad_order", &T::set_quad_order)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusionFE", tag);
	}

//	Convection Diffusion (FE) stabilization
	{
		typedef ConvectionDiffusionStabFE<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("ConvectionDiffusionStabFE").append(suffix);
		reg.template add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.template add_constructor<void (*)(const char*,const char*,number)>("Function(s)#Subset(s)#stabilization")
			.template add_constructor<void (*)(const char*,const char*,number,number)>("Function(s)#Subset(s)#stabilization")
			.add_method("set_quad_order", &T::set_quad_order)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusionStabFE", tag);
	}


//	Convection Diffusion FVCR
	{
		typedef ConvectionDiffusionFVCR<TDomain> T;
		typedef ConvectionDiffusionBase<TDomain> TBase;
		string name = string("ConvectionDiffusionFVCR").append(suffix);
		reg.template add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_upwind", &T::set_upwind)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusionFVCR", tag);
	}

//	Convection Diffusion FV
	{
		typedef ConvectionDiffusionFV<TDomain> T;
		typedef ConvectionDiffusionBase<TDomain> TBase;
		string name = string("ConvectionDiffusionFV").append(suffix);
		reg.template add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_quad_order", &T::set_quad_order)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusionFV", tag);
	}
}

template <int dim, typename TRegistry=ug::bridge::Registry>
static void Dimension(TRegistry& reg, string grp)
{
	string dimSuffix = GetDimensionSuffix<dim>();
	string dimTag = GetDimensionTag<dim>();

	//	singular sources and sinks 
	{
		typedef CDSingularSourcesAndSinks<dim> T;
		typedef typename T::point_sss_type TPointSSS;
		typedef typename T::line_sss_type TLineSSS;

		string point_name = string("CDPointSourcesSink").append(dimSuffix);
		reg.template add_class_<TPointSSS>(point_name, grp)
			.template add_constructor<void (*) (const std::vector<number>&)> ()
			.add_method ("set", static_cast<void (TPointSSS::*) (number)> (&TPointSSS::set))
			.add_method ("set", static_cast<void (TPointSSS::*) (typename TPointSSS::user_data_type)> (&TPointSSS::set))
#ifdef UG_FOR_LUA
			.add_method ("set", static_cast<void (TPointSSS::*) (LuaFunctionHandle)> (&TPointSSS::set))
#endif
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(point_name, "CDPointSourcesSink", dimTag);

		string line_name = string("CDLineSourcesSink").append(dimSuffix);
		reg.template add_class_<TLineSSS>(line_name, grp)
			.template add_constructor<void (*) (const std::vector<number>&, const std::vector<number>&)> ()
			.add_method ("set", static_cast<void (TLineSSS::*) (number)> (&TLineSSS::set))
			.add_method ("set", static_cast<void (TLineSSS::*) (LuaFunctionHandle)> (&TLineSSS::set))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(line_name, "CDLineSourcesSink", dimTag);

		string name = string("CDSingularSourcesAndSinks").append(dimSuffix);
		reg.template add_class_<T>(name, grp)
			.add_constructor()
			.add_method ("add_point", static_cast<void (T::*) (SmartPtr<TPointSSS>)> (&T::add_point))
			.add_method ("add_line", static_cast<void (T::*) (SmartPtr<TLineSSS>)> (&T::add_line))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CDSingularSourcesAndSinks", dimTag);
	}
}

}; // end Functionality

/**
 * Class exporting the functionality of the plugin restricted to 2 and 3 spatial
 * dimensions. All functionality that is to be used in scripts or visualization
 * only in 2d and 3d must be registered here.
 */
struct Functionality2d3d
{

/**
 * Function called for the registration of Domain dependent parts
 * of the plugin. All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TRegistry=ug::bridge::Registry>
static void Domain(TRegistry& reg, string grp)
{
	static const int dim = TDomain::dim;
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

//	Convection Diffusion FV1 for the low-dimensional fractures
	{
		typedef ConvectionDiffusionFractFV1<TDomain> T;
		typedef ConvectionDiffusionBase<TDomain> TBase;
		string name = string("ConvectionDiffusionFractFV1").append(suffix);
		reg.template add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_fract_manager", static_cast<void (T::*)(SmartPtr<DegeneratedLayerManager<dim> >)>(&T::set_fract_manager), "Sets the fracture manager", "Deg. fracture manager")
			.add_method("set_upwind", &T::set_upwind)
			.add_method("set_aperture", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_aperture), "", "Fract. aperture")
			.add_method("set_aperture", static_cast<void (T::*)(number)>(&T::set_aperture), "", "Fract. aperture")
#ifdef UG_FOR_LUA
			.add_method("set_aperture", static_cast<void (T::*)(const char*)>(&T::set_aperture), "", "Fract. aperture")
			.add_method("set_aperture", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_aperture), "", "Fract. aperture")
#endif
			.add_method("set_ortho_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_ortho_velocity), "", "Orthogonal Velocity Field")
			.add_method("set_ortho_velocity", static_cast<void (T::*)(number)>(&T::set_ortho_velocity), "", "Orthogonal Velocity Field")
#ifdef UG_FOR_LUA
			.add_method("set_ortho_velocity", static_cast<void (T::*)(const char*)>(&T::set_ortho_velocity), "", "Orthogonal Velocity Field")
			.add_method("set_ortho_velocity", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_ortho_velocity), "", "Orthogonal Velocity Field")
#endif
			.add_method("set_ortho_diffusion", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_ortho_diffusion), "", "Orthogonal Diffusion")
			.add_method("set_ortho_diffusion", static_cast<void (T::*)(number)>(&T::set_ortho_diffusion), "", "Orthogonal Diffusion")
#ifdef UG_FOR_LUA
			.add_method("set_ortho_diffusion", static_cast<void (T::*)(const char*)>(&T::set_ortho_diffusion), "", "Orthogonal Diffusion")
			.add_method("set_ortho_diffusion", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_ortho_diffusion), "", "Orthogonal Diffusion")
#endif
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusionFractFV1", tag);
	}
}

}; // end Functionality2d3d

// end group convection_diffusion
/// \}

} // end namespace ConvectionDiffusionPlugin



#ifndef UG_USE_PYBIND11


//! This function is called when the plugin is loaded.
extern "C"
void InitUGPlugin_ConvectionDiffusion(ug::bridge::Registry* reg, string grp)
{
	grp.append("/SpatialDisc/ElemDisc");
	typedef ConvectionDiffusionPlugin::Functionality Functionality;
	typedef ConvectionDiffusionPlugin::Functionality2d3d Functionality2d3d;

	try{
		RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dDependent<Functionality2d3d>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

#else // UG_USE_PYBIND11

template <typename TRegistry=ug::bridge::Registry>
void InitUGPlugin_ConvectionDiffusion_(TRegistry* reg, string grp)
{
	grp.append("/SpatialDisc/ElemDisc");
	typedef ConvectionDiffusionPlugin::Functionality Functionality;
	typedef ConvectionDiffusionPlugin::Functionality2d3d Functionality2d3d;

	try{
		RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality, TRegistry>(*reg,grp);
		RegisterDomain2d3dDependent<Functionality2d3d, TRegistry>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}


//! This function is called when the plugin is loaded.
extern "C" void
InitUGPlugin_ConvectionDiffusion(ug::bridge::Registry* reg, string grp)
{ InitUGPlugin_ConvectionDiffusion_(reg, grp); }

#ifdef UG_USE_PYBIND11
namespace ConvectionDiffusionPlugin{

	void Init(ug::pybind::Registry* reg, string grp)
	{ InitUGPlugin_ConvectionDiffusion_<ug::pybind::Registry>(reg, grp); }
}
#endif

#endif // UG_USE_PYBIND11

}// namespace ug
