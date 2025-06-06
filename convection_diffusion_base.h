/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_BASE__
#define __H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_BASE__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

namespace ug{
namespace ConvectionDiffusionPlugin{

/// \addtogroup convection_diffusion
/// \{

/// Discretization for the Convection-Diffusion Equation
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the convection diffusion equation.
 * The Equation has the form
 * \f[
 * 	\partial_t (m_1 c + m_2) - \nabla \cdot \left ( D \nabla c - \vec{v} c - \vec{F} \right )
 * 		+ r_1 \cdot c + r_2 = f + \nabla \cdot \vec{f}_2
 * \f]
 * with
 * <ul>
 * <li>	\f$ c \f$ is the unknown solution
 * <li>	\f$ m_1 \equiv m_1(\vec{x},t) \f$ is the Mass Scaling Term
 * <li>	\f$ m_2 \equiv m_1(\vec{x},t) \f$ is the Mass Term
 * <li>	\f$ D \equiv D(\vec{x},t) \f$ is the Diffusion Tensor
 * <li>	\f$ v \equiv \vec{v}(\vec{x},t) \f$ is the Velocity Field
 * <li>	\f$ F \equiv \vec{F}(\vec{x},t) \f$ is the Flux
 * <li>	\f$ r_1 \equiv r_1(\vec{x},t) \f$ is the Reaction Rate
 * <li>	\f$ r_2 \equiv r_2(\vec{x},t) \f$ is a Reaction Term
 * <li>	\f$ f \equiv f(\vec{x},t) \f$ is a Source Term
 * <li> \f$ \vec{f}_2 \equiv \vec{f}_2(\vec{x},t) \f$ is a Vector Source Term
 * </ul>
 *
 * \tparam	TDomain		Domain
 */
template<	typename TDomain>
class ConvectionDiffusionBase
: public IElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IElemDisc<TDomain> base_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor
		ConvectionDiffusionBase(const char* functions, const char* subsets);
	protected:
		void init_imports();
	public:
	///	sets the diffusion tensor
	/**
	 * This method sets the Diffusion tensor used in computations. If no
	 * Tensor is set, a zero value is assumed.
	 */
	///	\{
		void set_diffusion(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> > user);
		void set_diffusion(number val);
#ifdef UG_FOR_LUA
		void set_diffusion(const char* fctName);
		void set_diffusion(LuaFunctionHandle fct);
#endif
	///	\}

	///	sets the diffusion_Sw tensor
	/**
	 * This method sets the Diffusion tensor used in computations. If no
	 * Tensor is set, a zero value is assumed.
	 */
	///	\{
		void set_diffusion_Sw(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> > user);
		void set_diffusion_Sw(number val);
#ifdef UG_FOR_LUA
		void set_diffusion_Sw(const char* fctName);
		void set_diffusion_Sw(LuaFunctionHandle fct);
#endif
	///	\}
	
	///	sets the velocity field
	/**
	 * This method sets the Velocity field. If no field is provided a zero
	 * value is assumed.
	 */
	/// \{
		void set_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
		void set_velocity(const std::vector<number>& vVel);
#ifdef UG_FOR_LUA
		void set_velocity(const char* fctName);
		void set_velocity(LuaFunctionHandle fct);
#endif
	/// \}

	///	sets the set_darcyW field
	/**
	 * This method sets the Velocity field. If no field is provided a zero
	 * value is assumed.
	 */
	/// \{
		void set_darcyW(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
		void set_darcyW(const std::vector<number>& vVel);
#ifdef UG_FOR_LUA
		void set_darcyW(const char* fctName);
		void set_darcyW(LuaFunctionHandle fct);
#endif
	/// \}

	///	sets the set_darcyN field
	/**
	 * This method sets the Velocity field. If no field is provided a zero
	 * value is assumed.
	 */
	/// \{
		void set_darcyN(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
		void set_darcyN(const std::vector<number>& vVel);
#ifdef UG_FOR_LUA
		void set_darcyN(const char* fctName);
		void set_darcyN(LuaFunctionHandle fct);
#endif
	/// \}

	///	sets saturation
	/**
	 * This method sets the saturation value. The default value is 0.0.
	 */
	///	\{
		void set_saturationW(SmartPtr<CplUserData<number, dim> > user);
		void set_saturationW(number val);
#ifdef UG_FOR_LUA
		void set_saturationW(const char* fctName);
		void set_saturationW(LuaFunctionHandle fct);
#endif
	///	\}
	
	
	///	sets upwind saturation
	/**
	 * This method sets the saturation value. The default value is 0.0.
	 */
	///	\{
		void set_upwindsaturationW(SmartPtr<CplUserData<number, dim> > user);
		void set_upwindsaturationW(number val);
#ifdef UG_FOR_LUA
		void set_upwindsaturationW(const char* fctName);
		void set_upwindsaturationW(LuaFunctionHandle fct);
#endif
	///	\}
	
	
		///	sets modified saturation
	/**
	 * This method sets the saturation value. The default value is 0.0.
	 */
	///	\{
		void set_modifiedsaturationW(SmartPtr<CplUserData<number, dim> > user);
		void set_modifiedsaturationW(number val);
#ifdef UG_FOR_LUA
		void set_modifiedsaturationW(const char* fctName);
		void set_modifiedsaturationW(LuaFunctionHandle fct);
#endif
	///	\}


	///	sets upwind modified saturation
	/**
	 * This method sets the saturation value. The default value is 0.0.
	 */
	///	\{
		void set_upwindmodifiedsaturationW(SmartPtr<CplUserData<number, dim> > user);
		void set_upwindmodifiedsaturationW(number val);
#ifdef UG_FOR_LUA
		void set_upwindmodifiedsaturationW(const char* fctName);
		void set_upwindmodifiedsaturationW(LuaFunctionHandle fct);
#endif
	///	\}
	
	
		///	sets mass fraction Wc
	/**
	 * This method sets the mass fraction Wc value. The default value is 0.0.
	 */
	///	\{
		void set_MassFractionWc(SmartPtr<CplUserData<number, dim> > user);
		void set_MassFractionWc(number val);
#ifdef UG_FOR_LUA
		void set_MassFractionWc(const char* fctName);
		void set_MassFractionWc(LuaFunctionHandle fct);
#endif
	///	\}
	
		///	sets pressure of nowetting phase Pn
	/**
	 * This method sets the pressure of nowetting phase Pn value. The default value is 0.0.
	 */
	///	\{
		void set_PressurePn(SmartPtr<CplUserData<number, dim> > user);
		void set_PressurePn(number val);
#ifdef UG_FOR_LUA
		void set_PressurePn(const char* fctName);
		void set_PressurePn(LuaFunctionHandle fct);
#endif
	///	\}
	
	
		///	sets permeability
	/**
	 * This method sets the permeability value. The default value is 0.0.
	 */
	///	\{
		void set_permeability(SmartPtr<CplUserData<number, dim> > user);
		void set_permeability(number val);
#ifdef UG_FOR_LUA
		void set_permeability(const char* fctName);
		void set_permeability(LuaFunctionHandle fct);
#endif
	///	\}
	
	
		///	sets porosity
	/**
	 * This method sets the porosity value. The default value is 0.0.
	 */
	///	\{
		void set_porosity(SmartPtr<CplUserData<number, dim> > user);
		void set_porosity(number val);
#ifdef UG_FOR_LUA
		void set_porosity(const char* fctName);
		void set_porosity(LuaFunctionHandle fct);
#endif
	///	\}
	
	
		///	sets pd
	/**
	 * This method sets the pd value. The default value is 0.0.
	 */
	///	\{
		void set_pd(SmartPtr<CplUserData<number, dim> > user);
		void set_pd(number val);
#ifdef UG_FOR_LUA
		void set_pd(const char* fctName);
		void set_pd(LuaFunctionHandle fct);
#endif
	///	\}
	
	
		///	sets swr
	/**
	 * This method sets the pd value. The default value is 0.0.
	 */
	///	\{
		void set_swr(SmartPtr<CplUserData<number, dim> > user);
		void set_swr(number val);
#ifdef UG_FOR_LUA
		void set_swr(const char* fctName);
		void set_swr(LuaFunctionHandle fct);
#endif
	///	\}
	
		///	sets snr
	/**
	 * This method sets the pd value. The default value is 0.0.
	 */
	///	\{
		void set_snr(SmartPtr<CplUserData<number, dim> > user);
		void set_snr(number val);
#ifdef UG_FOR_LUA
		void set_snr(const char* fctName);
		void set_snr(LuaFunctionHandle fct);
#endif
	///	\}
	
	
		///	sets lambda
	/**
	 * This method sets the pd value. The default value is 0.0.
	 */
	///	\{
		void set_lambda(SmartPtr<CplUserData<number, dim> > user);
		void set_lambda(number val);
#ifdef UG_FOR_LUA
		void set_lambda(const char* fctName);
		void set_lambda(LuaFunctionHandle fct);
#endif
	///	\}
	
	
	
		///	sets minPd
	/**
	 * This method sets the minPd value. The default value is 0.0.
	 */
	///	\{
		void set_minPd(SmartPtr<CplUserData<number, dim> > user);
		void set_minPd(number val);
#ifdef UG_FOR_LUA
		void set_minPd(const char* fctName);
		void set_minPd(LuaFunctionHandle fct);
#endif
	///	\}
	
	
		///	sets minSwr
	/**
	 * This method sets the minSwr value. The default value is 0.0.
	 */
	///	\{
		void set_minSwr(SmartPtr<CplUserData<number, dim> > user);
		void set_minSwr(number val);
#ifdef UG_FOR_LUA
		void set_minSwr(const char* fctName);
		void set_minSwr(LuaFunctionHandle fct);
#endif
	///	\}
	
		///	sets minSnr
	/**
	 * This method sets the minSnr value. The default value is 0.0.
	 */
	///	\{
		void set_minSnr(SmartPtr<CplUserData<number, dim> > user);
		void set_minSnr(number val);
#ifdef UG_FOR_LUA
		void set_minSnr(const char* fctName);
		void set_minSnr(LuaFunctionHandle fct);
#endif
	///	\}
	
		///	sets minLambda
	/**
	 * This method sets the minLambda value. The default value is 0.0.
	 */
	///	\{
		void set_minLambda(SmartPtr<CplUserData<number, dim> > user);
		void set_minLambda(number val);
#ifdef UG_FOR_LUA
		void set_minLambda(const char* fctName);
		void set_minLambda(LuaFunctionHandle fct);
#endif
	///	\}
	
	///	sets the flux
	/**
	 * This method sets the Flux. If no field is provided a zero
	 * value is assumed.
	 */
	/// \{
		void set_flux(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
		void set_flux(const std::vector<number>& vVel);
#ifdef UG_FOR_LUA
		void set_flux(const char* fctName);
		void set_flux(LuaFunctionHandle fct);
#endif
	/// \}

	///	sets the reaction rate
	/**
	 * This method sets the Reaction Rate. A zero value is assumed as default.
	 */
	///	\{
		void set_reaction_rate(SmartPtr<CplUserData<number, dim> > user);
		void set_reaction_rate(number val);
#ifdef UG_FOR_LUA
		void set_reaction_rate(const char* fctName);
		void set_reaction_rate(LuaFunctionHandle fct);
#endif
	///	\}

	///	sets the reaction
	/**
	 * This method sets the Reaction. A zero value is assumed as default.
	 */
	///	\{
		void set_reaction(SmartPtr<CplUserData<number, dim> > user);
		void set_reaction(number val);
#ifdef UG_FOR_LUA
		void set_reaction(const char* fctName);
		void set_reaction(LuaFunctionHandle fct);
#endif
	///	\}

		void set_reaction_rate_explicit(SmartPtr<CplUserData<number, dim> > user);
		void set_reaction_rate_explicit(number val);
#ifdef UG_FOR_LUA
		void set_reaction_rate_explicit(const char* fctName);
#endif

		void set_reaction_explicit(SmartPtr<CplUserData<number, dim> > user);
		void set_reaction_explicit(number val);
#ifdef UG_FOR_LUA
		void set_reaction_explicit(const char* fctName);
#endif

		void set_source_explicit(SmartPtr<CplUserData<number, dim> > user);
		void set_source_explicit(number val);
#ifdef UG_FOR_LUA
		void set_source_explicit(const char* fctName);
#endif

	///	sets the source / sink term
	/**
	 * This method sets the source/sink value. A zero value is assumed as
	 * default.
	 */
	///	\{
		void set_source(SmartPtr<CplUserData<number, dim> > user);
		void set_source(number val);
#ifdef UG_FOR_LUA
		void set_source(const char* fctName);
		void set_source(LuaFunctionHandle fct);
#endif
	///	\}

	///	sets the vector source term
	/**
	 * This method sets the divergence of the source as an effect of an
	 * external field. A zero value is assumed as default, thus this term is
	 * ignored then.
	 */
	///	\{
		void set_vector_source(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
		void set_vector_source(const std::vector<number>& vVel);
#ifdef UG_FOR_LUA
		void set_vector_source(const char* fctName);
		void set_vector_source(LuaFunctionHandle fct);
#endif
	///	\}

	///	sets mass scale
	/**
	 * This method sets the mass scale value. The default value is 1.0.
	 */
	///	\{
		void set_mass_scale(SmartPtr<CplUserData<number, dim> > user);
		void set_mass_scale(number val);
#ifdef UG_FOR_LUA
		void set_mass_scale(const char* fctName);
		void set_mass_scale(LuaFunctionHandle fct);
#endif
	///	\}

	///	sets mass
	/**
	 * This method sets the mass value. The default value is 0.0.
	 */
	///	\{
		void set_mass(SmartPtr<CplUserData<number, dim> > user);
		void set_mass(number val);
#ifdef UG_FOR_LUA
		void set_mass(const char* fctName);
		void set_mass(LuaFunctionHandle fct);
#endif
	///	\}

	protected:
	///	Data import for Diffusion
		DataImport<MathMatrix<dim,dim>, dim> m_imDiffusion;

	///	Data import for the Velocity field
		DataImport<MathVector<dim>, dim > m_imVelocity;
		
		
	///	Data import for the Darcy Velocity field of wetting phase
		DataImport<MathVector<dim>, dim > m_imDarcyW;
		
	///	Data import for the Darcy Velocity field of non-wetting phase
		DataImport<MathVector<dim>, dim > m_imDarcyN;
		
	///	Data import for the saturationW
		DataImport<number, dim> m_imSaturationW;
		
	///	Data import for the upwind saturationW
		DataImport<number, dim> m_imUpwindSaturationW;
		
	///	Data import for the upwind saturationW
		DataImport<number, dim> m_imModifiedSaturationW;
		
	///	Data import for the upwind saturationW
		DataImport<number, dim> m_imUpwindModifiedSaturationW;
	
	///	Data import for Diffusion_Sw
		DataImport<MathMatrix<dim,dim>, dim> m_imDiffusion_Sw;
	
	///	Data import for the mass fraction Wc
		DataImport<number, dim> m_imMassFractionWc;
		
	///	Data import for the pressure of nowetting phase Pn
		DataImport<number, dim> m_imPressurePn;
		
		
	///	Data import for the Permeability
		DataImport<number, dim> m_imPermeability;

	///	Data import for the Porosity
		DataImport<number, dim> m_imPorosity;

	///	Data import for the EntryPressure
		DataImport<number, dim> m_imEntryPressure;		
		
	///	Data import for the ResidualAqueous
		DataImport<number, dim> m_imResidualAqueous;
		
	///	Data import for the ResidualCarbonic
		DataImport<number, dim> m_imResidualCarbonic;
		
	///	Data import for the BrooksCoreyNumber
		DataImport<number, dim> m_imBrooksCoreyNumber;
		
		
	///	Data import for the MinPd
		DataImport<number, dim> m_imMinPd;
	///	Data import for the MinSwr
		DataImport<number, dim> m_imMinSwr;
	///	Data import for the MinSnr
		DataImport<number, dim> m_imMinSnr;
	///	Data import for the MinLambda
		DataImport<number, dim> m_imMinLambda;
		
	///	Data import for the Flux
		DataImport<MathVector<dim>, dim > m_imFlux;

	///	Data import for the reaction term
		DataImport<number, dim> m_imReactionRate;

	///	Data import for the reaction term
		DataImport<number, dim> m_imReaction;

	///	Data import for the reaction_rate term explicit
		DataImport<number, dim> m_imReactionRateExpl;

	///	Data import for the reaction term explicit
		DataImport<number, dim> m_imReactionExpl;

	///	Data import for the source term explicit
		DataImport<number, dim> m_imSourceExpl;

	///	Data import for the right-hand side (volume)
		DataImport<number, dim> m_imSource;

	///	Data import for the right-hand side (vector)
		DataImport<MathVector<dim>, dim > m_imVectorSource;

	///	Data import for the mass scale
		DataImport<number, dim> m_imMassScale;

	///	Data import for the mass
		DataImport<number, dim> m_imMass;

	private:
	///	returns if local time series is needed
		virtual bool requests_local_time_series() {return false;}

	public:
		typedef SmartPtr<CplUserData<number, dim> > NumberExport;
		typedef SmartPtr<CplUserData<MathVector<dim>, dim> > GradExport;
		
	///	returns the export of the value of associated unknown function
		virtual SmartPtr<CplUserData<number, dim> > upwindmodifiedvalue();
		
	///	returns the export of the value of associated unknown function
		virtual SmartPtr<CplUserData<number, dim> > upwindvalue();	
		
	///	returns the export of the value of associated unknown function
		virtual SmartPtr<CplUserData<number, dim> > modifiedvalue();
		
	///	returns the export of the value of associated unknown function
		virtual SmartPtr<CplUserData<number, dim> > value();

	///	returns the export of the gradient of associated unknown function
		virtual SmartPtr<CplUserData<MathVector<dim>, dim> > gradient();
		
	///	returns the export of the gradient of associated unknown function
		virtual SmartPtr<CplUserData<MathVector<dim>, dim> > gradient_pd();
		
	///	returns the export of the gradient of associated unknown function
		virtual SmartPtr<CplUserData<MathVector<dim>, dim> > gradient_pd2();

	protected:
	///	Export for the concentration
		SmartPtr<DataExport<number, dim> > m_exUpwindModifiedValue;

	///	Export for the concentration
		SmartPtr<DataExport<number, dim> > m_exUpwindValue;	
	
	///	Export for the concentration
		SmartPtr<DataExport<number, dim> > m_exModifiedValue;

	///	Export for the concentration
		SmartPtr<DataExport<number, dim> > m_exValue;
		
	///	Export for the gradient of concentration
		SmartPtr<DataExport<MathVector<dim>, dim> > m_exGrad;
		
	///	Export for the gradient of concentration * pd
		SmartPtr<DataExport<MathVector<dim>, dim> > m_exGrad_pd;
		
	///	Export for the gradient of concentration * pd
		SmartPtr<DataExport<MathVector<dim>, dim> > m_exGrad_pd2;
};

// end group convection_diffusion
/// \}

} // end ConvectionDiffusionPlugin
} // end namespace ug


#endif /*__H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_BASE__*/
