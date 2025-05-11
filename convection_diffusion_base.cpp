/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#include "convection_diffusion_base.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace ConvectionDiffusionPlugin{

////////////////////////////////////////////////////////////////////////////////
//	user data
////////////////////////////////////////////////////////////////////////////////

//////// Diffusion

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_diffusion(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> > user)
{
	m_imDiffusion.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_diffusion(number val)
{
	if(val == 0.0) set_diffusion(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> >());
	else set_diffusion(make_sp(new ConstUserMatrix<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_diffusion(const char* fctName)
{
	set_diffusion(LuaUserDataFactory<MathMatrix<dim,dim>, dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_diffusion(LuaFunctionHandle fct)
{
	set_diffusion(make_sp(new LuaUserData<MathMatrix<dim,dim>, dim>(fct)));
}
#endif

//////// Diffusion_Sw

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_diffusion_Sw(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> > user)
{
	m_imDiffusion_Sw.set_data(user);
	m_imDiffusion_Sw.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_diffusion_Sw(number val)
{
	if(val == 0.0) set_diffusion_Sw(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> >());
	else set_diffusion_Sw(make_sp(new ConstUserMatrix<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_diffusion_Sw(const char* fctName)
{
	set_diffusion_Sw(LuaUserDataFactory<MathMatrix<dim,dim>, dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_diffusion_Sw(LuaFunctionHandle fct)
{
	set_diffusion_Sw(make_sp(new LuaUserData<MathMatrix<dim,dim>, dim>(fct)));
}
#endif

//////// Velocity

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> > user)
{
	m_imVelocity.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_velocity(const std::vector<number>& vVel)
{
	bool bZero = true;
	for(size_t i = 0; i < vVel.size(); ++i){
		if(vVel[i] != 0.0) bZero = false;
	}

	if(bZero) set_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> >());
	else set_velocity(SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_velocity(const char* fctName)
{
	set_velocity(LuaUserDataFactory<MathVector<dim>,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_velocity(LuaFunctionHandle fct)
{
	set_velocity(make_sp(new LuaUserData<MathVector<dim>,dim>(fct)));
}
#endif


//////// darcyW

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_darcyW(SmartPtr<CplUserData<MathVector<dim>, dim> > user)
{
	m_imDarcyW.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_darcyW(const std::vector<number>& vVel)
{
	bool bZero = true;
	for(size_t i = 0; i < vVel.size(); ++i){
		if(vVel[i] != 0.0) bZero = false;
	}

	if(bZero) set_darcyW(SmartPtr<CplUserData<MathVector<dim>, dim> >());
	else set_darcyW(SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_darcyW(const char* fctName)
{
	set_darcyW(LuaUserDataFactory<MathVector<dim>,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_darcyW(LuaFunctionHandle fct)
{
	set_darcyW(make_sp(new LuaUserData<MathVector<dim>,dim>(fct)));
}
#endif

//////// darcyN
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_darcyN(SmartPtr<CplUserData<MathVector<dim>, dim> > user)
{
	m_imDarcyN.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_darcyN(const std::vector<number>& vVel)
{
	bool bZero = true;
	for(size_t i = 0; i < vVel.size(); ++i){
		if(vVel[i] != 0.0) bZero = false;
	}

	if(bZero) set_darcyN(SmartPtr<CplUserData<MathVector<dim>, dim> >());
	else set_darcyN(SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_darcyN(const char* fctName)
{
	set_darcyN(LuaUserDataFactory<MathVector<dim>,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_darcyN(LuaFunctionHandle fct)
{
	set_darcyN(make_sp(new LuaUserData<MathVector<dim>,dim>(fct)));
}
#endif

//////// SaturationW

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_saturationW(SmartPtr<CplUserData<number, dim> > user)
{
	m_imSaturationW.set_data(user);
	m_imSaturationW.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_saturationW(number val)
{
	if(val == 0.0) set_saturationW(SmartPtr<CplUserData<number, dim> >());
	else set_saturationW(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_saturationW(const char* fctName)
{
	set_saturationW(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_saturationW(LuaFunctionHandle fct)
{
	set_saturationW(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif


//////// Upwind SaturationW

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_upwindsaturationW(SmartPtr<CplUserData<number, dim> > user)
{
	m_imUpwindSaturationW.set_data(user);
	m_imUpwindSaturationW.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_upwindsaturationW(number val)
{
	if(val == 0.0) set_upwindsaturationW(SmartPtr<CplUserData<number, dim> >());
	else set_upwindsaturationW(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_upwindsaturationW(const char* fctName)
{
	set_upwindsaturationW(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_upwindsaturationW(LuaFunctionHandle fct)
{
	set_upwindsaturationW(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif


//////// Modified SaturationW

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_modifiedsaturationW(SmartPtr<CplUserData<number, dim> > user)
{
	m_imModifiedSaturationW.set_data(user);
	m_imModifiedSaturationW.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_modifiedsaturationW(number val)
{
	if(val == 0.0) set_modifiedsaturationW(SmartPtr<CplUserData<number, dim> >());
	else set_modifiedsaturationW(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_modifiedsaturationW(const char* fctName)
{
	set_modifiedsaturationW(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_modifiedsaturationW(LuaFunctionHandle fct)
{
	set_modifiedsaturationW(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif

//////// Upwind Modified SaturationW

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_upwindmodifiedsaturationW(SmartPtr<CplUserData<number, dim> > user)
{
	m_imUpwindModifiedSaturationW.set_data(user);
	m_imUpwindModifiedSaturationW.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_upwindmodifiedsaturationW(number val)
{
	if(val == 0.0) set_upwindmodifiedsaturationW(SmartPtr<CplUserData<number, dim> >());
	else set_upwindmodifiedsaturationW(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_upwindmodifiedsaturationW(const char* fctName)
{
	set_upwindmodifiedsaturationW(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_upwindmodifiedsaturationW(LuaFunctionHandle fct)
{
	set_upwindmodifiedsaturationW(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif


//////// MassFractionWc

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_MassFractionWc(SmartPtr<CplUserData<number, dim> > user)
{
	m_imMassFractionWc.set_data(user);
	m_imMassFractionWc.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_MassFractionWc(number val)
{
	if(val == 0.0) set_MassFractionWc(SmartPtr<CplUserData<number, dim> >());
	else set_MassFractionWc(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_MassFractionWc(const char* fctName)
{
	set_MassFractionWc(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_MassFractionWc(LuaFunctionHandle fct)
{
	set_MassFractionWc(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif

//////// PressurePn

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_PressurePn(SmartPtr<CplUserData<number, dim> > user)
{
	m_imPressurePn.set_data(user);
	m_imPressurePn.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_PressurePn(number val)
{
	if(val == 0.0) set_PressurePn(SmartPtr<CplUserData<number, dim> >());
	else set_PressurePn(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_PressurePn(const char* fctName)
{
	set_PressurePn(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_PressurePn(LuaFunctionHandle fct)
{
	set_PressurePn(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif


//////// Permeability

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_permeability(SmartPtr<CplUserData<number, dim> > user)
{
	m_imPermeability.set_data(user);
	m_imPermeability.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_permeability(number val)
{
	if(val == 0.0) set_permeability(SmartPtr<CplUserData<number, dim> >());
	else set_permeability(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_permeability(const char* fctName)
{
	set_permeability(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_permeability(LuaFunctionHandle fct)
{
	set_permeability(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif

//////// Porosity

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_porosity(SmartPtr<CplUserData<number, dim> > user)
{
	m_imPorosity.set_data(user);
	m_imPorosity.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_porosity(number val)
{
	if(val == 0.0) set_porosity(SmartPtr<CplUserData<number, dim> >());
	else set_porosity(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_porosity(const char* fctName)
{
	set_porosity(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_porosity(LuaFunctionHandle fct)
{
	set_porosity(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif


//////// EntryPressure
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_pd(SmartPtr<CplUserData<number, dim> > user)
{
	m_imEntryPressure.set_data(user);
	m_imEntryPressure.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_pd(number val)
{
	if(val == 0.0) set_pd(SmartPtr<CplUserData<number, dim> >());
	else set_pd(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_pd(const char* fctName)
{
	set_pd(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_pd(LuaFunctionHandle fct)
{
	set_pd(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif


//////// ResidualAqueous
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_swr(SmartPtr<CplUserData<number, dim> > user)
{
	m_imResidualAqueous.set_data(user);
	m_imResidualAqueous.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_swr(number val)
{
	if(val == 0.0) set_swr(SmartPtr<CplUserData<number, dim> >());
	else set_swr(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_swr(const char* fctName)
{
	set_swr(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_swr(LuaFunctionHandle fct)
{
	set_swr(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif

//////// ResidualCarbonic
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_snr(SmartPtr<CplUserData<number, dim> > user)
{
	m_imResidualCarbonic.set_data(user);
	m_imResidualCarbonic.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_snr(number val)
{
	if(val == 0.0) set_snr(SmartPtr<CplUserData<number, dim> >());
	else set_snr(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_snr(const char* fctName)
{
	set_snr(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_snr(LuaFunctionHandle fct)
{
	set_snr(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif


//////// BrooksCoreyNumber
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_lambda(SmartPtr<CplUserData<number, dim> > user)
{
	m_imBrooksCoreyNumber.set_data(user);
	m_imBrooksCoreyNumber.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_lambda(number val)
{
	if(val == 0.0) set_lambda(SmartPtr<CplUserData<number, dim> >());
	else set_lambda(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_lambda(const char* fctName)
{
	set_lambda(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_lambda(LuaFunctionHandle fct)
{
	set_lambda(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif


//////// MinPd

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minPd(SmartPtr<CplUserData<number, dim> > user)
{
	m_imMinPd.set_data(user);
	m_imMinPd.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minPd(number val)
{
	if(val == 0.0) set_minPd(SmartPtr<CplUserData<number, dim> >());
	else set_minPd(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minPd(const char* fctName)
{
	set_minPd(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minPd(LuaFunctionHandle fct)
{
	set_minPd(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif


//////// MinSwr

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minSwr(SmartPtr<CplUserData<number, dim> > user)
{
	m_imMinSwr.set_data(user);
	m_imMinSwr.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minSwr(number val)
{
	if(val == 0.0) set_minSwr(SmartPtr<CplUserData<number, dim> >());
	else set_minSwr(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minSwr(const char* fctName)
{
	set_minSwr(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minSwr(LuaFunctionHandle fct)
{
	set_minSwr(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif


//////// MinSnr

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minSnr(SmartPtr<CplUserData<number, dim> > user)
{
	m_imMinSnr.set_data(user);
	m_imMinSnr.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minSnr(number val)
{
	if(val == 0.0) set_minSnr(SmartPtr<CplUserData<number, dim> >());
	else set_minSnr(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minSnr(const char* fctName)
{
	set_minSnr(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minSnr(LuaFunctionHandle fct)
{
	set_minSnr(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif


//////// MinLambda

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minLambda(SmartPtr<CplUserData<number, dim> > user)
{
	m_imMinLambda.set_data(user);
	m_imMinLambda.set_comp_lin_defect(false);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minLambda(number val)
{
	if(val == 0.0) set_minLambda(SmartPtr<CplUserData<number, dim> >());
	else set_minLambda(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minLambda(const char* fctName)
{
	set_minLambda(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_minLambda(LuaFunctionHandle fct)
{
	set_minLambda(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif


//////// Flux

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_flux(SmartPtr<CplUserData<MathVector<dim>, dim> > user)
{
	m_imFlux.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_flux(const std::vector<number>& vVel)
{
	bool bZero = true;
	for(size_t i = 0; i < vVel.size(); ++i){
		if(vVel[i] != 0.0) bZero = false;
	}

	if(bZero) set_flux(SmartPtr<CplUserData<MathVector<dim>, dim> >());
	else set_flux(SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_flux(const char* fctName)
{
	set_flux(LuaUserDataFactory<MathVector<dim>,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_flux(LuaFunctionHandle fct)
{
	set_flux(make_sp(new LuaUserData<MathVector<dim>,dim>(fct)));
}
#endif

//////// Reaction Rate

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_rate(SmartPtr<CplUserData<number, dim> > user)
{
	m_imReactionRate.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_rate(number val)
{
	if(val == 0.0) set_reaction_rate(SmartPtr<CplUserData<number, dim> >());
	else set_reaction_rate(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_rate(const char* fctName)
{
	set_reaction_rate(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_rate(LuaFunctionHandle fct)
{
	set_reaction_rate(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif

//////// Reaction Rate Explicit

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_rate_explicit(SmartPtr<CplUserData<number, dim> > user)
{
	m_imReactionRateExpl.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_rate_explicit(number val)
{
	if(val == 0.0) set_reaction_rate_explicit(SmartPtr<CplUserData<number, dim> >());
	else set_reaction_rate_explicit(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_rate_explicit(const char* fctName)
{
	set_reaction_rate_explicit(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif

//////// Reaction

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction(SmartPtr<CplUserData<number, dim> > user)
{
	m_imReaction.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction(number val)
{
	if(val == 0.0) set_reaction(SmartPtr<CplUserData<number, dim> >());
	else set_reaction(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction(const char* fctName)
{
	set_reaction(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction(LuaFunctionHandle fct)
{
	set_reaction(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif

//////// Reaction Explicit

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_explicit(SmartPtr<CplUserData<number, dim> > user)
{
	m_imReactionExpl.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_explicit(number val)
{
	if(val == 0.0) set_reaction_explicit(SmartPtr<CplUserData<number, dim> >());
	else set_reaction_explicit(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_explicit(const char* fctName)
{
	set_reaction_explicit(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif


//////// Source

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_source(SmartPtr<CplUserData<number, dim> > user)
{
	m_imSource.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_source(number val)
{
	if(val == 0.0) set_source(SmartPtr<CplUserData<number, dim> >());
	else set_source(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_source(const char* fctName)
{
	set_source(LuaUserDataFactory<number,dim>::create(fctName));
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_source(LuaFunctionHandle fct)
{
	set_source(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif

//////// Source explicit

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_source_explicit(SmartPtr<CplUserData<number, dim> > user)
{
	m_imSourceExpl.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_source_explicit(number val)
{
	if(val == 0.0) set_source_explicit(SmartPtr<CplUserData<number, dim> >());
	else set_source_explicit(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_source_explicit(const char* fctName)
{
	set_source_explicit(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif


//////// Vector Source

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_vector_source(SmartPtr<CplUserData<MathVector<dim>, dim> > user)
{
	m_imVectorSource.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_vector_source(const std::vector<number>& vVel)
{
	bool bZero = true;
	for(size_t i = 0; i < vVel.size(); ++i){
		if(vVel[i] != 0.0) bZero = false;
	}

	if(bZero) set_vector_source(SmartPtr<CplUserData<MathVector<dim>, dim> >());
	else set_vector_source(SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_vector_source(const char* fctName)
{
	set_vector_source(LuaUserDataFactory<MathVector<dim>,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_vector_source(LuaFunctionHandle fct)
{
	set_vector_source(make_sp(new LuaUserData<MathVector<dim>,dim>(fct)));
}
#endif

//////// Mass Scale

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_mass_scale(SmartPtr<CplUserData<number, dim> > user)
{
	m_imMassScale.set_data(user);
	//Shuai Debug
	//m_imMassScale.set_comp_lin_defect(false);
	//Shuai Debug
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_mass_scale(number val)
{
	if(val == 0.0) set_mass_scale(SmartPtr<CplUserData<number, dim> >());
	else set_mass_scale(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_mass_scale(const char* fctName)
{
	set_mass_scale(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_mass_scale(LuaFunctionHandle fct)
{
	set_mass_scale(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif

//////// Mass

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_mass(SmartPtr<CplUserData<number, dim> > user)
{
	m_imMass.set_data(user);
	//Shuai Debug
	//m_imMass.set_comp_lin_defect(false);
	//Shuai Debug
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_mass(number val)
{
	if(val == 0.0) set_mass(SmartPtr<CplUserData<number, dim> >());
	else set_mass(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_mass(const char* fctName)
{
	set_mass(LuaUserDataFactory<number,dim>::create(fctName));
}
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_mass(LuaFunctionHandle fct)
{
	set_mass(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif

////////////////////////////////////////////////////////////////////////////////
//	Exports
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
typename ConvectionDiffusionBase<TDomain>::NumberExport
ConvectionDiffusionBase<TDomain>::
upwindmodifiedvalue() {return m_exUpwindModifiedValue;}

template <typename TDomain>
typename ConvectionDiffusionBase<TDomain>::NumberExport
ConvectionDiffusionBase<TDomain>::
upwindvalue() {return m_exUpwindValue;}

template <typename TDomain>
typename ConvectionDiffusionBase<TDomain>::NumberExport
ConvectionDiffusionBase<TDomain>::
modifiedvalue() {return m_exModifiedValue;}

template <typename TDomain>
typename ConvectionDiffusionBase<TDomain>::NumberExport
ConvectionDiffusionBase<TDomain>::
value() {return m_exValue;}

template <typename TDomain>
typename ConvectionDiffusionBase<TDomain>::GradExport
ConvectionDiffusionBase<TDomain>::
gradient() {return m_exGrad;}

template <typename TDomain>
typename ConvectionDiffusionBase<TDomain>::GradExport
ConvectionDiffusionBase<TDomain>::
gradient_pd() {return m_exGrad_pd;}

template <typename TDomain>
typename ConvectionDiffusionBase<TDomain>::GradExport
ConvectionDiffusionBase<TDomain>::
gradient_pd2() {return m_exGrad_pd2;}

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
init_imports()
{
	//	register imports
		this->register_import(m_imDiffusion);
		this->register_import(m_imVelocity);
		
		this->register_import(m_imDarcyW);
		this->register_import(m_imDarcyN);
		this->register_import(m_imSaturationW);
		this->register_import(m_imUpwindSaturationW);
		this->register_import(m_imModifiedSaturationW);
		this->register_import(m_imUpwindModifiedSaturationW);
		this->register_import(m_imDiffusion_Sw);
		this->register_import(m_imMassFractionWc);
		this->register_import(m_imPressurePn);
		
		this->register_import(m_imPermeability);
		this->register_import(m_imPorosity);
		this->register_import(m_imEntryPressure);
		this->register_import(m_imResidualAqueous);
		this->register_import(m_imResidualCarbonic);
		this->register_import(m_imBrooksCoreyNumber);
		
		this->register_import(m_imMinPd);
		this->register_import(m_imMinSwr);
		this->register_import(m_imMinSnr);
		this->register_import(m_imMinLambda);
		
		this->register_import(m_imFlux);
		this->register_import(m_imReactionRate);
		this->register_import(m_imReaction);
		this->register_import(m_imReactionRateExpl);
		this->register_import(m_imReactionExpl);
		this->register_import(m_imSourceExpl);
		this->register_import(m_imSource);
		this->register_import(m_imVectorSource);
		this->register_import(m_imMassScale);
		this->register_import(m_imMass);

		m_imMassScale.set_mass_part();
		m_imMass.set_mass_part();
		m_imSource.set_rhs_part();
		m_imVectorSource.set_rhs_part();
		m_imSourceExpl.set_expl_part();
		m_imReactionExpl.set_expl_part();
		m_imReactionRateExpl.set_expl_part();
}

template<typename TDomain>
ConvectionDiffusionBase<TDomain>::
ConvectionDiffusionBase(const char* functions, const char* subsets)
 : IElemDisc<TDomain>(functions,subsets),
   m_exUpwindModifiedValue(new DataExport<number, dim>(functions)),
   m_exUpwindValue(new DataExport<number, dim>(functions)),
   m_exModifiedValue(new DataExport<number, dim>(functions)),
   m_exValue(new DataExport<number, dim>(functions)),
   m_exGrad(new DataExport<MathVector<dim>, dim>(functions)),
   m_exGrad_pd(new DataExport<MathVector<dim>, dim>(functions)),
   m_exGrad_pd2(new DataExport<MathVector<dim>, dim>(functions))
{
//	check number of functions
	if(this->num_fct() != 1)
		UG_THROW("Wrong number of functions: The ElemDisc 'ConvectionDiffusion'"
					   " needs exactly "<<1<<" symbolic function.");
// init all imports
	init_imports();

//	default value for mass scale
	set_mass_scale(1.0);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class ConvectionDiffusionBase<Domain1d>;
#endif
#ifdef UG_DIM_2
template class ConvectionDiffusionBase<Domain2d>;
#endif
#ifdef UG_DIM_3
template class ConvectionDiffusionBase<Domain3d>;
#endif

} // end namespace ConvectionDiffusionPlugin
} // namespace ug
