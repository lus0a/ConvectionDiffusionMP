/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Shuai Lu
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

#include "convection_diffusion_mp.h"

#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape.h"

namespace ug{
namespace ConvectionDiffusionPlugin{

DebugID DID_CONV_DIFF_MP("CONV_DIFF_MP");

////////////////////////////////////////////////////////////////////////////////
//	general
////////////////////////////////////////////////////////////////////////////////

static number krw_max = 0.9;
static number krn_max = 0.5;
static number m = 2;
static number n = 2;
//static number n = 1;

//static number Permeability = 4.845e-13;
static number ViscosityW = 1e-3;

static number swr = 0.0;
static number snr = 0.0;
//static number swr = 0.2;
//static number snr = 0.1;

//static number lambda = 4.2;

static number lambda = 2.0;

/*
static number krw(number sw)
{
	number sw_eff = (sw - swr) / (1-snr-swr);
	if (sw_eff>0 && sw_eff<1)
		return pow(sw_eff, 2/lambda+3);
	else if (sw_eff<=0)
		return 0;
	else
		return 1;
}

static number D_krw(number sw)
{
	number sw_eff = (sw - swr) / (1-snr-swr);
	if (sw_eff>0 && sw_eff<1)
		return (2/lambda+3) * pow(sw_eff, 2/lambda+2) / (1-snr-swr);
	else if (sw_eff<=0)
		return 0;
	else
		return (2/lambda+3) / (1-snr-swr);
}

static number krn(number sw)
{
    number sw_eff = (sw - swr) / (1-snr-swr);
	if (sw_eff>0 && sw_eff<1)
		return (1-sw_eff)*(1-sw_eff) * (1-pow(sw_eff, 2/lambda+1));
	else if (sw_eff>=1)
		return 0;
	else
		return 1;
}
*/

/*
// For Extended Buckley Leverett
static number krw1(number sw)
{
	number sw_eff = (sw - swr) / (1-snr-swr);
	if (sw_eff>0 && sw_eff<1)
		return 1.831*pow(sw_eff, 4);
	else if (sw_eff<=0)
		return 0;
	else
		return 1.831;
}

static number D_krw1(number sw)
{
	number sw_eff = (sw - swr) / (1-snr-swr);
	if (sw_eff>0 && sw_eff<1)
		return 1.831*4*pow(sw_eff, 3)/(1-snr-swr);
	else if (sw_eff<=0)
		return 0;
	else
		return 1.831*4;
}

static number krn1(number sw)
{
    number sw_eff = (sw - swr) / (1-snr-swr);
	if (sw_eff>0 && sw_eff<1)
		return 0.75*pow(1-1.25*sw_eff, 2)* (1-1.652*pow(sw_eff, 2));
	else if (sw_eff>=1)
		return 0.75*pow(-0.25, 2)* (-0.652);
	else
		return 0.75;
}

static number krw2(number sw)
{
	number sw_eff = (sw - swr) / (1-snr-swr);
	if (sw_eff>0 && sw_eff<1)
		return 0.4687*pow(sw_eff, 2);
	else if (sw_eff<=0)
		return 0;
	else
		return 0.4687;
}

static number D_krw2(number sw)
{
	number sw_eff = (sw - swr) / (1-snr-swr);
	if (sw_eff>0 && sw_eff<1)
		return 0.4687*2*sw_eff;
	else if (sw_eff<=0)
		return 0;
	else
		return 0.4687*2;
}

static number krn2(number sw)
{
    number sw_eff = (sw - swr) / (1-snr-swr);
	if (sw_eff>0 && sw_eff<1)
		return 0.25 * pow(1-1.25*sw_eff, 2);
	else if (sw_eff>=1)
		return 0.25 * pow(-0.25, 2);
	else
		return 0.25;
}
*/




///*
static number krw(number sw)
{
	number sw_eff = (sw - swr) / (1-snr-swr);
	if (sw_eff>0 && sw_eff<1)
		return krw_max * pow(sw_eff, m);
	else if (sw_eff<=0)
		return 0;
	else
		return krw_max;
}

static number D_krw(number sw)
{
	number sw_eff = (sw - swr) / (1-snr-swr);
	if (sw_eff>0 && sw_eff<1)
		return krw_max * m * pow(sw_eff, m-1) / (1-snr-swr);
	else if (sw_eff<=0)
		return 0;
	else
		return krw_max * m / (1-snr-swr);
}

static number krn(number sw)
{	
	number sw_eff = (sw - swr) / (1-snr-swr);
	if (sw_eff>0 && sw_eff<1)
		return krn_max * pow(1-sw_eff, n);
	else if (sw_eff>=1)
		return 0;
	else
		return krn_max;
}
//*/


///*
static number DensityW(number wC, number pN)
{
	pN=pN/10000000;
	//int pN=5;
	number result = (1005+307.5*wC-6.98*3
		-589.8*wC*wC+52.53*wC*pN); //For Pn in [25,35]Mpa
	//number result = (1000+259.3*wC
	//	-446.8*wC*wC+30.83*wC*pN);
	//number result = (1006+259.3*wC-2.853*pN
	//	-446.8*wC*wC+30.83*wC*pN);
	/*
	number result = (988.1+473.2*wC-6.628*pN
		-797.3*wC*wC+17.55*wC*pN+1.672*pN*pN
		-82.44*wC*wC*wC+98.26*wC*wC*pN-8.403*wC*pN*pN);
	*/
	return result;
}
//*/

// Leverett J
static number LeverJ(number sw_eff) {
	if (sw_eff>0 && sw_eff<1)
		return pow(sw_eff, -1/lambda);
	else if (sw_eff>=1)
		return 1;
	else
		return pow(10,10);
}
// Inverse Leverett J
static number InverseLeverJ(number J) {
	return pow(J, -lambda);
}

///*
static number Pd(number perm, number poro)
{
	return 7.37 * pow(poro/perm,0.43);
}
static number D_diffusion_Sw(number sw, number perm, number poro)
{
	number sw_eff = (sw - swr) / (1-snr-swr);
	if (sw_eff>0 && sw_eff<1)
		return -perm/ViscosityW *Pd(perm, poro)*(-1/lambda)*krw_max* (m-1/lambda-1) * pow(sw_eff, (m-1/lambda-2)) / pow(1-snr-swr, 2);
	else if (sw_eff>=1)
		return -perm/ViscosityW *Pd(perm, poro)*(-1/lambda)*krw_max* (m-1/lambda-1) / pow(1-snr-swr, 2);
	else
		return 0;
}
//*/

/*
static number Pd(number perm, number poro)
{
	return 1000;
	//return 0.10;
	//return 7.37 * pow(1/perm,0.43);
}
static number D_diffusion_Sw(number sw, number perm, number poro)
{
	number sw_eff = (sw - swr) / (1-snr-swr);
	if (sw_eff>0 && sw_eff<1)
		return -perm/ViscosityW *Pd(perm, poro)*(-1/lambda)* (1/lambda+2) * pow(sw_eff, (1/lambda+1)) / pow(1-snr-swr, 2);
	else if (sw_eff>=1)
		return -perm/ViscosityW *Pd(perm, poro)*(-1/lambda)* (1/lambda+2) / pow(1-snr-swr, 2);
	else
		return 0;
		
}
*/

static number Modify_sw(number sw, number minPd, number perm, number poro)
{	
	number elePd = Pd(perm, poro);
	if ( abs(minPd - elePd) < 1 ){
		return sw;
	}
	number sw_eff = (sw - swr) / (1-snr-swr);
	if ( elePd >= minPd*LeverJ(sw_eff) ){
		return 1-snr;
	}
	sw_eff = InverseLeverJ( minPd*LeverJ(sw_eff)/elePd );	
	return sw_eff*(1-snr-swr)+swr;
}

template<typename TDomain>
ConvectionDiffusionMP<TDomain>::
ConvectionDiffusionMP(const char* functions, const char* subsets)
 : ConvectionDiffusionBase<TDomain>(functions,subsets),
   m_spConvShape(new ConvectionShapesNoUpwind<dim>),
   m_bNonRegularGrid(false), m_bCondensedFV(false)
{
	register_all_funcs(m_bNonRegularGrid);
}

template<typename TDomain>
void ConvectionDiffusionMP<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
//	check number
	if(vLfeID.size() != 1)
		UG_THROW("ConvectionDiffusion: Wrong number of functions given. "
				"Need exactly "<<1);

	if(vLfeID[0].order() != 1 || vLfeID[0].type() != LFEID::LAGRANGE)
		UG_THROW("ConvectionDiffusion FV Scheme only implemented for 1st order.");

//	remember
	m_bNonRegularGrid = bNonRegularGrid;

//	update assemble functions
	register_all_funcs(m_bNonRegularGrid);
}

template<typename TDomain>
bool ConvectionDiffusionMP<TDomain>::
use_hanging() const
{
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// Assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void ConvectionDiffusionMP<TDomain>::
prep_assemble_loop()
{
	if (m_sss_mngr.valid ())
	//	reset the markers of the point sources and sinks (as there are no marks for the lines)
		m_sss_mngr->init_all_point_sss ();
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{

//	check, that upwind has been set
//TODO: It is really used only if we have the convection
	if(m_spConvShape.invalid())
		UG_THROW("ConvectionDiffusionMP::prep_elem_loop:"
						" Upwind has not been set.");

//	set local positions
	if(!TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;
		TFVGeom& geo = GeomProvider<TFVGeom>::get();
		const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
		const size_t numSCVFip = geo.num_scvf_ips();
		const MathVector<refDim>* vSCVip = geo.scv_local_ips();
		const size_t numSCVip = geo.num_scv_ips();
		m_imDiffusion.template 		set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imVelocity.template 		set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		
		m_imDarcyW.template 		set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imDarcyN.template 		set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imSaturationW.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imDiffusion_Sw.template 	set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imMassFractionWc.template set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imPressurePn.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		
		m_imPermeability.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imPorosity.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imMinPd.template 			set_local_ips<refDim>(vSCVip,numSCVip, false);
		
		m_imFlux.template 			set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imSource.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imVectorSource.template 	set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imReactionRate.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imReaction.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imReactionRateExpl.template set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imReactionExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imSourceExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imMassScale.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imMass.template 			set_local_ips<refDim>(vSCVip,numSCVip, false);

		//	init upwind for element type
		if(!m_spConvShape->template set_geometry_type<TFVGeom>(geo))
			UG_THROW("ConvectionDiffusionMP::prep_elem_loop:"
						" Cannot init upwind for element type.");
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
fsh_elem_loop()
{}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
// 	Update Geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

	try{
		UG_DLOG(DID_CONV_DIFF_MP, 2, ">>OCT_DISC_DEBUG: " << "convection_diffusion_fv1.cpp: " << "prep_elem(): update(): "<< roid << std::endl);
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}UG_CATCH_THROW("ConvectionDiffusionMP::prep_elem:"
						" Cannot update Finite Volume Geometry.");

//	set local positions
	if(TFVGeom::usesHangingNodes)
	{
		const int refDim = TElem::dim;
		const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
		const size_t numSCVFip = geo.num_scvf_ips();
		const MathVector<refDim>* vSCVip = geo.scv_local_ips();
		const size_t numSCVip = geo.num_scv_ips();
		m_imDiffusion.template 		set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imVelocity.template 		set_local_ips<refDim>(vSCVFip,numSCVFip);
		
		m_imDarcyW.template 		set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imDarcyN.template 		set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imSaturationW.template 	set_local_ips<refDim>(vSCVip,numSCVip);
		
		m_imFlux.template 			set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imSource.template 		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imVectorSource.template 	set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imReactionRate.template 	set_local_ips<refDim>(vSCVip,numSCVip);
		m_imReaction.template 		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imReactionRateExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip);
		m_imReactionExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip);
		m_imSourceExpl.template		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imMassScale.template 		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imMass.template 			set_local_ips<refDim>(vSCVip,numSCVip);

		if(m_spConvShape.valid())
			if(!m_spConvShape->template set_geometry_type<TFVGeom>(geo))
				UG_THROW("ConvectionDiffusionMP::prep_elem:"
								" Cannot init upwind for element type.");
	}

	//	set global positions
	const MathVector<dim>* vSCVFip = geo.scvf_global_ips();
	const size_t numSCVFip = geo.num_scvf_ips();
	const MathVector<dim>* vSCVip = geo.scv_global_ips();
	const size_t numSCVip = geo.num_scv_ips();
	m_imDiffusion.			set_global_ips(vSCVFip, numSCVFip);
	m_imVelocity.			set_global_ips(vSCVFip, numSCVFip);
	
	m_imDarcyW.				set_global_ips(vSCVFip, numSCVFip);
	m_imDarcyN.				set_global_ips(vSCVFip, numSCVFip);
	m_imSaturationW.		set_global_ips(vSCVip, numSCVip);
	m_imDiffusion_Sw.		set_global_ips(vSCVFip, numSCVFip);
	m_imMassFractionWc.		set_global_ips(vSCVip, numSCVip);
	m_imPressurePn.			set_global_ips(vSCVip, numSCVip);
	
	m_imPermeability.		set_global_ips(vSCVip, numSCVip);
	m_imPorosity.			set_global_ips(vSCVip, numSCVip);
	m_imMinPd.				set_global_ips(vSCVip, numSCVip);
	
	m_imFlux.				set_global_ips(vSCVFip, numSCVFip);
	m_imSource.				set_global_ips(vSCVip, numSCVip);
	m_imVectorSource.		set_global_ips(vSCVFip, numSCVFip);
	m_imReactionRate.		set_global_ips(vSCVip, numSCVip);
	m_imReactionRateExpl.	set_global_ips(vSCVip, numSCVip);
	m_imReactionExpl.		set_global_ips(vSCVip, numSCVip);
	m_imSourceExpl.			set_global_ips(vSCVip, numSCVip);
	m_imReaction.			set_global_ips(vSCVip, numSCVip);
	m_imMassScale.			set_global_ips(vSCVip, numSCVip);
	m_imMass.				set_global_ips(vSCVip, numSCVip);
}

template <class TVector>
static TVector CalculateCenter(GridObject* o, const TVector* coords)
{
	TVector v;
	VecSet(v, 0);

	size_t numCoords = 0;
	switch(o->base_object_id()){
		case VERTEX: numCoords = 1; break;
		case EDGE: numCoords = static_cast<Edge*>(o)->num_vertices(); break;
		case FACE: numCoords = static_cast<Face*>(o)->num_vertices(); break;
		case VOLUME: numCoords = static_cast<Volume*>(o)->num_vertices(); break;
		default: UG_THROW("Unknown element type."); break;
	}

	for(size_t i = 0; i < numCoords; ++i)
		VecAdd(v, v, coords[i]);

	if(numCoords > 0)
		VecScale(v, v, 1. / (number)numCoords);

	return v;
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
add_sss_jac_elem
(
	LocalMatrix& J, ///< the matrix to update
	const LocalVector& u, ///< current solution
	GridObject* elem, ///< the element
	const TFVGeom& geo, ///< the FV geometry for that element
	size_t i, ///< index of the SCV
	number flux ///< flux through source/sink (premultiplied by the length for lines)
)
{
	size_t co = geo.scv(i).node_id ();
	
	if (flux < 0.0)
		// sink
		J(_C_, co, _C_, co) -= flux;
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	Diff. Tensor times Gradient
	MathVector<dim> Dgrad;
	
//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo, false);
	
//	Consider capillary trapping	or not
	number u_modified[convShape.num_sh()];
	if ( m_imPermeability.data_given() && m_imMinPd.data_given() ){
		for (size_t i = 0; i < convShape.num_sh(); ++i){
			u_modified[i] = Modify_sw(u(_C_, i), m_imMinPd[i], m_imPermeability[i], m_imPorosity[i]);
		}
	}
	else{
		//UG_THROW ("Permeability or MinPd is missing in some elements! "<<"Permeability:"<<m_imPermeability.data_given()<<" MinPd:"<<m_imMinPd.data_given());
		for (size_t i = 0; i < convShape.num_sh(); ++i)
			u_modified[i] = u(_C_, i);
	}

//	Diffusion and Velocity Term
	if(m_imDiffusion.data_given() || m_imDiffusion_Sw.data_given() || m_imVelocity.data_given() || m_imDarcyW.data_given())
	{
	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
		////////////////////////////////////////////////////
		// Diffusive Term
		////////////////////////////////////////////////////
			if(m_imDiffusion.data_given())
			{
				#ifdef UG_ENABLE_DEBUG_LOGS
				//	DID_CONV_DIFF_MP
				number D_diff_flux_sum = 0.0;
				#endif

			// 	loop shape functions
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
				// 	Compute Diffusion Tensor times Gradient
					MatVecMult(Dgrad, m_imDiffusion[ip], scvf.global_grad(sh));

				//	Compute flux at IP
					const number D_diff_flux = VecDot(Dgrad, scvf.normal());
					UG_DLOG(DID_CONV_DIFF_MP, 2, ">>OCT_DISC_DEBUG: " << "convection_diffusion_mp.cpp: " << "add_jac_A_elem(): " << "sh # "  << sh << " ; normalSize scvf # " << ip << ": " << VecLength(scvf.normal()) << "; \t from "<< scvf.from() << "; to " << scvf.to() << "; D_diff_flux: " << D_diff_flux << "; scvf.global_grad(sh): " << scvf.global_grad(sh) << std::endl);

				// 	Add flux term to local matrix // HIER MATRIXINDIZES!!!
					UG_ASSERT((scvf.from() < J.num_row_dof(_C_)) && (scvf.to() < J.num_col_dof(_C_)),
							  "Bad local dof-index on element with object-id " << elem->base_object_id()
							  << " with center: " << CalculateCenter(elem, vCornerCoords));

					J(_C_, scvf.from(), _C_, sh) -= D_diff_flux;
					J(_C_, scvf.to()  , _C_, sh) += D_diff_flux;

					#ifdef UG_ENABLE_DEBUG_LOGS
					//	DID_CONV_DIFF_MP
					D_diff_flux_sum += D_diff_flux;
					#endif
				}

				UG_DLOG(DID_CONV_DIFF_MP, 2, "D_diff_flux_sum = " << D_diff_flux_sum << std::endl << std::endl);
			}
		////////////////////////////////////////////////////
		// Diffusive_Sw Term
		////////////////////////////////////////////////////
			if(m_imDiffusion_Sw.data_given() && m_imPermeability.data_given())
			{
				#ifdef UG_ENABLE_DEBUG_LOGS
				//	DID_CONV_DIFF_MP
				number D_diff_flux_Sw_sum = 0.0;
				#endif

			// 	loop shape functions
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
					MathMatrix<dim,dim> Real_Diffusion_Sw;
					MatScale(Real_Diffusion_Sw, m_imPermeability[0]*Pd(m_imPermeability[0], m_imPorosity[0]), m_imDiffusion_Sw[ip]);
				// 	Compute Diffusion Tensor times Gradient
					MatVecMult(Dgrad, Real_Diffusion_Sw, scvf.global_grad(sh));
					
				//  Compute D_diffusion * lambda_c * sum( sw_i*grad(lambda_i))
					MathVector<dim> sum; VecSet(sum, 0.0);
					number sw_integ = 0.0;
					MathVector<dim> Dgrad2;
					
					for (size_t i = 0; i < scvf.num_sh(); ++i)
					{
						VecScaleAppend(sum, u_modified[i], scvf.global_grad(i));
						sw_integ += u_modified[i]*scvf.shape(i);
					}
					VecScale(sum, sum, scvf.shape(sh));
					MathMatrix<dim,dim> D_diffusion_Sw_Mat;
					for (size_t j = 0; j < dim; ++j)
						for (size_t k = 0; k < dim; ++k)
						{
							if (j == k)
							{
								if (m_imMassFractionWc.data_given() && m_imPressurePn.data_given())
									D_diffusion_Sw_Mat(j,k) = DensityW(m_imMassFractionWc[ip], m_imPressurePn[ip]) * D_diffusion_Sw(sw_integ, m_imPermeability[0], m_imPorosity[0]);
									//D_diffusion_Sw_Mat(j,k) = DensityW(m_imMassFractionWc[ip]) * D_diffusion_Sw(sw_integ, m_imPermeability[0], m_imPorosity[0]);
								else if (m_imPressurePn.data_given())
									D_diffusion_Sw_Mat(j,k) = DensityW(0, m_imPressurePn[ip]) * D_diffusion_Sw(sw_integ, m_imPermeability[0], m_imPorosity[0]);
									//D_diffusion_Sw_Mat(j,k) = DensityW(0) * D_diffusion_Sw(sw_integ, m_imPermeability[0], m_imPorosity[0]);
								else //immiscible & incompressible, densityW=1000
									D_diffusion_Sw_Mat(j,k) = 1000 * D_diffusion_Sw(sw_integ, m_imPermeability[0], m_imPorosity[0]);
								/*
								//Shuai, Debug
								if (D_diffusion_Sw_Mat(j,k)<0)
								{
									number xxxx= 0;
									xxxx =D_diffusion_Sw_Mat(j,k);
								}
								*/
							}
							else
								D_diffusion_Sw_Mat(j,k) = 0.0;
						}
					MatVecMult(Dgrad2, D_diffusion_Sw_Mat, sum);
					
					VecAdd(Dgrad, Dgrad, Dgrad2);
					
				//	Compute flux at IP
					const number D_diff_flux_Sw = VecDot(Dgrad, scvf.normal());
					UG_DLOG(DID_CONV_DIFF_MP, 2, ">>OCT_DISC_DEBUG: " << "convection_diffusion_mp.cpp: " << "add_jac_A_elem(): " << "sh # "  << sh << " ; normalSize scvf # " << ip << ": " << VecLength(scvf.normal()) << "; \t from "<< scvf.from() << "; to " << scvf.to() << "; D_diff_flux_Sw: " << D_diff_flux_Sw << "; scvf.global_grad(sh): " << scvf.global_grad(sh) << std::endl);

				// 	Add flux term to local matrix // HIER MATRIXINDIZES!!!
					UG_ASSERT((scvf.from() < J.num_row_dof(_C_)) && (scvf.to() < J.num_col_dof(_C_)),
							  "Bad local dof-index on element with object-id " << elem->base_object_id()
							  << " with center: " << CalculateCenter(elem, vCornerCoords));
					
					
					J(_C_, scvf.from(), _C_, sh) -= D_diff_flux_Sw;
					J(_C_, scvf.to()  , _C_, sh) += D_diff_flux_Sw;

					#ifdef UG_ENABLE_DEBUG_LOGS
					//	DID_CONV_DIFF_MP
					D_diff_flux_Sw_sum += D_diff_flux_Sw;
					#endif
				}

				UG_DLOG(DID_CONV_DIFF_MP, 2, "D_diff_flux_Sw_sum = " << D_diff_flux_Sw_sum << std::endl << std::endl);
			}

		////////////////////////////////////////////////////
		// Convective Term
		////////////////////////////////////////////////////
			if(m_imVelocity.data_given())
			{
			//	Add Flux contribution
				for(size_t sh = 0; sh < convShape.num_sh(); ++sh)
				{
					const number D_conv_flux = convShape(ip, sh);

				//	Add flux term to local matrix
					J(_C_, scvf.from(), _C_, sh) += D_conv_flux;
					J(_C_, scvf.to(),   _C_, sh) -= D_conv_flux;
				}
			}
			
			if(m_imDarcyW.data_given())
			{
			//	Add Flux contribution
			
			/* This is the old implementation of upwind
				number D_krw_up = 0.0;
				number factor = 0.0;
				// 	loop shape functions
				for(size_t sh = 0; sh < convShape.num_sh(); ++sh)
				{
					D_krw_up += D_krw( u_modified[sh] )* convShape(ip, sh);
					factor += convShape(ip, sh);
				}
				//Shuai, debugging
				sw = sw / factor;
			*/		
				
				for(size_t sh = 0; sh < convShape.num_sh(); ++sh)
				{
					/*//For Extended Buckley Leverett
					number D_krw_value = 0.0;
					if (m_imPermeability.data_given())
					{
						if (m_imPermeability[0]>5*pow(10,-14))
							D_krw_value = D_krw1( u_modified[sh] );
						else
							D_krw_value = D_krw2( u_modified[sh] );
					}
					const number D_conv_flux = D_krw_value * convShape(ip, sh);
					*/
					
					const number D_conv_flux = D_krw( u_modified[sh] ) * convShape(ip, sh);
					
	
				//	Add flux term to local matrix
					J(_C_, scvf.from(), _C_, sh) += D_conv_flux;
					J(_C_, scvf.to(),   _C_, sh) -= D_conv_flux;
				}
			}

			// no explicit dependency on flux import
		}
	}

	//UG_LOG("Local Matrix is: \n"<<J<<"\n");

////////////////////////////////////////////////////
// Reaction Term (using lumping)
////////////////////////////////////////////////////

	if(m_imReactionRate.data_given())
	{
	// 	loop Sub Control Volume (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);
			
		// 	get associated node
			const int co = scv.node_id();
			
		// 	Add to local matrix
			J(_C_, co, _C_, co) += m_imReactionRate[ip] * scv.volume();
		}
	}
	
//	reaction term does not explicitly depend on the associated unknown function

////////////////////////////////
// Singular sources and sinks
////////////////////////////////

	if (m_sss_mngr.valid () && (m_sss_mngr->num_points () != 0 || m_sss_mngr->num_lines () != 0))
    {
    	typedef typename TDomain::position_accessor_type t_pos_accessor;
    	typedef typename CDSingularSourcesAndSinks<dim>::template
    		point_iterator<TElem,t_pos_accessor,TFVGeom> t_pnt_sss_iter;
    	typedef typename CDSingularSourcesAndSinks<dim>::template
    		line_iterator<TElem,t_pos_accessor,TFVGeom> t_lin_sss_iter;

		t_pos_accessor& aaPos = this->domain()->position_accessor();
		Grid& grid = (Grid&) *this->domain()->grid();
		
		for(size_t i = 0; i < geo.num_scv(); i++)
		{
			size_t co = geo.scv(i).node_id ();
			
		//	point sources
			for (t_pnt_sss_iter pnt (m_sss_mngr.get(), (TElem *) elem, grid, aaPos, geo, co);
				! pnt.is_over (); ++pnt)
			{
				FVPointSourceOrSink<dim, cd_point_sss_data<dim> > * pnt_sss = *pnt;
				if (! pnt_sss->marked_for (elem, co))
					continue;
				pnt_sss->compute (pnt_sss->position (), this->time (), -1); //TODO: set the subset id instead of -1
				this->template add_sss_jac_elem<TElem, TFVGeom> (J, u, elem, geo, i,
					pnt_sss->flux ());
			}
			
		//	line sources
			for (t_lin_sss_iter line (m_sss_mngr.get(), (TElem *) elem, grid, aaPos, geo, co);
				! line.is_over (); ++line)
			{
				FVLineSourceOrSink<dim, cd_line_sss_data<dim> > * line_sss = *line;
				number len = VecDistance (line.seg_start (), line.seg_end ());
				line_sss->compute (line.seg_start (), this->time (), -1); //TODO: set the subset id instead of -1
				this->template add_sss_jac_elem<TElem, TFVGeom> (J, u, elem, geo, i,
					line_sss->flux () * len);
			}
		}
    }
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	if(!m_imMassScale.data_given()) return;

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local matrix
		J(_C_, co, _C_, co) += scv.volume() * m_imMassScale[ip];
	}

//	m_imMass part does not explicitly depend on associated unknown function
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
add_sss_def_elem
(
	LocalVector& d, ///< the defect to update
	const LocalVector& u, ///< current solution
	GridObject* elem, ///< the element
	const TFVGeom& geo, ///< the FV geometry for that element
	size_t i, ///< index of the SCV
	number flux ///< flux through source/sink (premultiplied by the length for lines)
)
{
	size_t co = geo.scv(i).node_id ();

	//	Consider capillary trapping	or not
	number u_modified;
	if ( m_imPermeability.data_given() && m_imMinPd.data_given() ){
		u_modified = Modify_sw(u(_C_, co), m_imMinPd[co], m_imPermeability[co], m_imPorosity[co]);
	}
	else{
		u_modified = u(_C_, co);
	}
	
	if (flux > 0.0)
		// source
		d(_C_, co) -= flux;
	else
		// sink
		d(_C_, co) -= flux * u_modified;
		
	if (i != co)
	{
		UG_THROW("More than one shape function per scv ? Please check!");
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo, false);
	
	//	Consider capillary trapping	or not
	number u_modified[convShape.num_sh()];
	if ( m_imPermeability.data_given() && m_imMinPd.data_given() ){
		for (size_t i = 0; i < convShape.num_sh(); ++i)
			u_modified[i] = Modify_sw(u(_C_, i), m_imMinPd[i], m_imPermeability[i], m_imPorosity[i]);
	}
	else{
		for (size_t i = 0; i < convShape.num_sh(); ++i)
			u_modified[i] = u(_C_, i);
	}
	

	if(m_imDiffusion.data_given() || m_imDiffusion_Sw.data_given() || m_imVelocity.data_given() || m_imFlux.data_given() || m_imDarcyW.data_given() || (m_imDarcyN.data_given() && m_imSaturationW.data_given()))
	{
	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		/////////////////////////////////////////////////////
		// Diffusive Term
		/////////////////////////////////////////////////////
			if(m_imDiffusion.data_given())
			{
			//	to compute D \nabla c
				MathVector<dim> Dgrad_c, grad_c;

			// 	compute gradient and shape at ip
				VecSet(grad_c, 0.0);
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecScaleAppend(grad_c, u_modified[sh], scvf.global_grad(sh));

			//	scale by diffusion tensor
				MatVecMult(Dgrad_c, m_imDiffusion[ip], grad_c);

			// 	Compute flux
				const number diff_flux = VecDot(Dgrad_c, scvf.normal());

			// 	Add to local defect
				d(_C_, scvf.from()) -= diff_flux;
				d(_C_, scvf.to()  ) += diff_flux;
			}
			
		/////////////////////////////////////////////////////
		// Diffusive_Sw Term
		/////////////////////////////////////////////////////
			if(m_imDiffusion_Sw.data_given())
			{
			//	to compute D \nabla c
				MathVector<dim> Dgrad_c, grad_c;

			// 	compute gradient and shape at ip
				VecSet(grad_c, 0.0);
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecScaleAppend(grad_c, u_modified[sh], scvf.global_grad(sh));

				MathMatrix<dim,dim> Real_Diffusion_Sw;
				if (m_imPermeability.data_given())
					MatScale(Real_Diffusion_Sw, m_imPermeability[0]*Pd(m_imPermeability[0], m_imPorosity[0]), m_imDiffusion_Sw[ip]);
				else
					MatScale(Real_Diffusion_Sw, 1, m_imDiffusion_Sw[ip]);
			//	scale by diffusion tensor
				MatVecMult(Dgrad_c, Real_Diffusion_Sw, grad_c);

			// 	Compute flux
				const number diff_flux_Sw = VecDot(Dgrad_c, scvf.normal());

			// 	Add to local defect
				d(_C_, scvf.from()) -= diff_flux_Sw;
				d(_C_, scvf.to()  ) += diff_flux_Sw;
			}

		/////////////////////////////////////////////////////
		// Convective Term
		/////////////////////////////////////////////////////
			if(m_imVelocity.data_given())
			{
			//	sum up convective flux using convection shapes
				number conv_flux = 0.0;
				for(size_t sh = 0; sh < convShape.num_sh(); ++sh)
				{
					conv_flux += u_modified[sh] * convShape(ip, sh);

				}
					

			//  add to local defect
				d(_C_, scvf.from()) += conv_flux;
				d(_C_, scvf.to()  ) -= conv_flux;
			}
			
			if(m_imDarcyW.data_given())
			{
				
			//	sum up convective flux using convection shapes
				number conv_flux = 0.0;
				for(size_t sh = 0; sh < convShape.num_sh(); ++sh)
				{
					/*//For Extended Buckley Leverett
					if (m_imPermeability.data_given())
					{
						if (m_imPermeability[0]>5*pow(10,-14))
							conv_flux += krw1( u_modified[sh] ) * convShape(ip, sh);
						else
							conv_flux += krw2( u_modified[sh] ) * convShape(ip, sh);
					}
					*/
					
					conv_flux += krw( u_modified[sh] ) * convShape(ip, sh);
						
				}
					
			//  add to local defect
				d(_C_, scvf.from()) += conv_flux;
				d(_C_, scvf.to()  ) -= conv_flux;
			}


		//Flux Term - krn(sW)*DarcyN
			if(m_imDarcyN.data_given() && m_imSaturationW.data_given())
			{
			//	sum up convective flux using convection shapes
				number conv_flux = 0.0;
				for(size_t sh = 0; sh < convShape.num_sh(); ++sh)
				{
					/*//For Extended Buckley Leverett
					if (m_imPermeability.data_given())
					{
						if (m_imPermeability[0]>5*pow(10,-14))
							conv_flux += krn1( m_imSaturationW[sh] ) * convShape(ip, sh);
						else
							conv_flux += krn2( m_imSaturationW[sh] ) * convShape(ip, sh);
					}
					*/
					
					// Nothing to change, since m_imSaturationW is modified before import
					conv_flux += krn( m_imSaturationW[sh] ) * convShape(ip, sh);
					
				}
					
				//  add to local defect
				d(_C_, scvf.from()) += conv_flux;
				d(_C_, scvf.to()  ) -= conv_flux;
				
			}
		/////////////////////////////////////////////////////
		// Flux Term
		/////////////////////////////////////////////////////
			if(m_imFlux.data_given())
			{
			//	sum up flux
				const number flux = VecDot(m_imFlux[ip], scvf.normal());

			//  add to local defect
				d(_C_, scvf.from()) += flux;
				d(_C_, scvf.to()  ) -= flux;
				// For debug
				//flux_Flux = flux;
			}
			
		}
	}

			
//	reaction rate
	if(m_imReactionRate.data_given())
	{
	// 	loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

		// 	get associated node
			const int co = scv.node_id();

		// 	Add to local defect
			d(_C_, co) += u_modified[co] * m_imReactionRate[ip] * scv.volume();
		}
	}

//	reaction term
	if(m_imReaction.data_given())
	{
	// 	loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

		// 	get associated node
			const int co = scv.node_id();

		// 	Add to local defect
			d(_C_, co) += m_imReaction[ip] * scv.volume();
		}
	}

////////////////////////////////
// Singular sources and sinks
////////////////////////////////

    if (m_sss_mngr.valid () && (m_sss_mngr->num_points () != 0 || m_sss_mngr->num_lines () != 0))
    {
    	typedef typename TDomain::position_accessor_type t_pos_accessor;
    	typedef typename CDSingularSourcesAndSinks<dim>::template
    		point_iterator<TElem,t_pos_accessor,TFVGeom> t_pnt_sss_iter;
    	typedef typename CDSingularSourcesAndSinks<dim>::template
    		line_iterator<TElem,t_pos_accessor,TFVGeom> t_lin_sss_iter;
    	
		t_pos_accessor& aaPos = this->domain()->position_accessor();
		Grid& grid = (Grid&) *this->domain()->grid();
		
		for(size_t i = 0; i < geo.num_scv(); i++)
		{
			size_t co = geo.scv(i).node_id ();
			
		//	point sources
			for (t_pnt_sss_iter pnt (m_sss_mngr.get(), (TElem *) elem, grid, aaPos, geo, co);
				! pnt.is_over (); ++pnt)
			{
				FVPointSourceOrSink<dim, cd_point_sss_data<dim> > * pnt_sss = *pnt;
				if (! pnt_sss->marked_for (elem, co))
					continue;
				pnt_sss->compute (pnt_sss->position (), this->time (), -1); //TODO: set the subset id instead of -1
				this->template add_sss_def_elem<TElem, TFVGeom> (d, u, elem, geo, i,
					pnt_sss->flux ());
			}
			
		//	line sources
			for (t_lin_sss_iter line (m_sss_mngr.get(), (TElem *) elem, grid, aaPos, geo, co);
				! line.is_over (); ++line)
			{
				FVLineSourceOrSink<dim, cd_line_sss_data<dim> > * line_sss = *line;
				number len = VecDistance (line.seg_start (), line.seg_end ());
				line_sss->compute (line.seg_start (), this->time (), -1); //TODO: set the subset id instead of -1
				this->template add_sss_def_elem<TElem, TFVGeom> (d, u, elem, geo, i,
					line_sss->flux () * len);
			}
		}
    }
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
add_def_A_expl_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reaction rate
	if(m_imReactionRateExpl.data_given())
	{
	// 	loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

		// 	get associated node
			const int co = scv.node_id();

		//	Consider capillary trapping	or not
			number u_modified;
			if ( m_imPermeability.data_given() && m_imMinPd.data_given() ){
				u_modified = Modify_sw(u(_C_, co), m_imMinPd[co], m_imPermeability[co], m_imPorosity[co]);
			}
			else{
				u_modified = u(_C_, co);
			}

		// 	Add to local defect
			d(_C_, co) += u_modified * m_imReactionRateExpl[ip] * scv.volume();
		}
	}

//	reaction
	if(m_imReactionExpl.data_given())
	{
	// 	loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

		// 	get associated node
			const int co = scv.node_id();

		// 	Add to local defect
			d(_C_, co) += m_imReactionExpl[ip] * scv.volume();
		}
	}

	if(m_imSourceExpl.data_given())
	{
		// 	loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
			// 	get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

			// 	get associated node
			const int co = scv.node_id();

			// 	Add to local rhs
			d(_C_, co) -= m_imSourceExpl[ip] * scv.volume();
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	if(!m_imMassScale.data_given() && !m_imMass.data_given()) return;
	
// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	//	Consider capillary trapping	or not
		number u_modified;
		if ( m_imPermeability.data_given() && m_imMinPd.data_given() ){
			u_modified = Modify_sw(u(_C_, co), m_imMinPd[co], m_imPermeability[co], m_imPorosity[co]);
		}
		else{
			u_modified = u(_C_, co);
		}
		
		//	mass value
		number val = 0.0;

	//	multiply by scaling
		if(m_imMassScale.data_given())
			val += m_imMassScale[ip] * u_modified;

	//	add mass
		if(m_imMass.data_given())
			val += m_imMass[ip];

	// 	Add to local defect
		d(_C_, co) += val * scv.volume();
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	// loop Sub Control Volumes (SCV)
	if ( m_imSource.data_given() ) {
		for ( size_t ip = 0; ip < geo.num_scv(); ++ip ) {
			// get current SCV
			const typename TFVGeom::SCV& scv = geo.scv( ip );

			// get associated node
			const int co = scv.node_id();

			// Add to local rhs
			d(_C_, co) += m_imSource[ip] * scv.volume();
			//UG_LOG("d(_C_, co) = " << d(_C_, co) << "; \t ip " << ip << "; \t co " << co << "; \t scv_vol " << scv.volume() << "; \t m_imSource[ip] " << m_imSource[ip] << std::endl);
		}
	}

	// loop Sub Control Volumes (SCVF)
	if ( m_imVectorSource.data_given() ) {
		for ( size_t ip = 0; ip < geo.num_scvf(); ++ip ) {
			// get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf( ip );

			// Add to local rhs
			number flux = VecDot(m_imVectorSource[ip], scvf.normal());
			d(_C_, scvf.from()) -= flux;
			d(_C_, scvf.to()  ) += flux;
		}
	}
}


// ////////////////////////////////
//   error estimation (begin)   ///

//	prepares the loop over all elements of one type for the computation of the error estimator
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
{
	//	get the error estimator data object and check that it is of the right type
	//	we check this at this point in order to be able to dispense with this check later on
	//	(i.e. in prep_err_est_elem and compute_err_est_A_elem())
	if (this->m_spErrEstData.get() == NULL)
	{
		UG_THROW("No ErrEstData object has been given to this ElemDisc!");
	}

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (!err_est_data)
	{
		UG_THROW("Dynamic cast to SideAndElemErrEstData failed."
				<< std::endl << "Make sure you handed the correct type of ErrEstData to this discretization.");
	}


//	check that upwind has been set
	if (m_spConvShape.invalid())
		UG_THROW("ConvectionDiffusionMP::prep_err_est_elem_loop: "
				 "Upwind has not been set.");

//	set local positions
	if (!TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;

		// get local IPs
		size_t numSideIPs, numElemIPs;
		const MathVector<refDim>* sideIPs;
		const MathVector<refDim>* elemIPs;
		try
		{
			numSideIPs = err_est_data->num_all_side_ips(roid);
			numElemIPs = err_est_data->num_elem_ips(roid);
			sideIPs = err_est_data->template side_local_ips<refDim>(roid);
			elemIPs = err_est_data->template elem_local_ips<refDim>(roid);

			if (!sideIPs || !elemIPs) return;	// are NULL if TElem is not of the same dim as TDomain
		}
		UG_CATCH_THROW("Integration points for error estimator cannot be set.");

		// set local IPs in imports
		m_imDiffusion.template 		set_local_ips<refDim>(sideIPs, numSideIPs, false);
		m_imVelocity.template 		set_local_ips<refDim>(sideIPs, numSideIPs, false);
		
		m_imDarcyW.template 		set_local_ips<refDim>(sideIPs, numSideIPs, false);
		m_imDarcyN.template 		set_local_ips<refDim>(sideIPs, numSideIPs, false);
		
		m_imFlux.template 			set_local_ips<refDim>(sideIPs, numSideIPs, false);
		m_imSource.template 		set_local_ips<refDim>(elemIPs, numElemIPs, false);
		m_imVectorSource.template 	set_local_ips<refDim>(sideIPs, numSideIPs, false);
		m_imReactionRate.template 	set_local_ips<refDim>(elemIPs, numElemIPs, false);
		m_imReaction.template 		set_local_ips<refDim>(elemIPs, numElemIPs, false);
		m_imMassScale.template 		set_local_ips<refDim>(elemIPs, numElemIPs, false);
		m_imMass.template 			set_local_ips<refDim>(elemIPs, numElemIPs, false);

		//	init upwind for element type
		TFVGeom& geo = GeomProvider<TFVGeom>::get();
		if (!m_spConvShape->template set_geometry_type<TFVGeom>(geo))
			UG_THROW("ConvectionDiffusionMP::prep_err_est_elem_loop: "
					 "Cannot init upwind for element type.");

		// store values of shape functions in local IPs
		LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> trialSpace
					= Provider<LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> >::get();

		m_shapeValues.resize(numElemIPs, numSideIPs, trialSpace.num_sh());
		for (size_t ip = 0; ip < numElemIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtElemIP(ip), elemIPs[ip]);
		for (size_t ip = 0; ip < numSideIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtSideIP(ip), sideIPs[ip]);
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

// 	update geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try
	{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("ConvectionDiffusionMP::prep_err_est_elem: Cannot update Finite Volume Geometry.");

//	roid
	ReferenceObjectID roid = elem->reference_object_id();

//	set local positions
	if (TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;

		size_t numSideIPs, numElemIPs;
		const MathVector<refDim>* sideIPs;
		const MathVector<refDim>* elemIPs;
		try
		{
			numSideIPs = err_est_data->num_all_side_ips(roid);
			numElemIPs = err_est_data->num_elem_ips(roid);
			sideIPs = err_est_data->template side_local_ips<refDim>(roid);
			elemIPs = err_est_data->template elem_local_ips<refDim>(roid);

			if (!sideIPs || !elemIPs) return;	// are NULL if TElem is not of the same dim as TDomain
		}
		UG_CATCH_THROW("Integration points for error estimator cannot be set.");

		m_imDiffusion.template 		set_local_ips<refDim>(sideIPs, numSideIPs);
		m_imVelocity.template 		set_local_ips<refDim>(sideIPs, numSideIPs);
		
		m_imDarcyW.template 		set_local_ips<refDim>(sideIPs, numSideIPs);
		m_imDarcyN.template 		set_local_ips<refDim>(sideIPs, numSideIPs);
		
		m_imFlux.template 			set_local_ips<refDim>(sideIPs, numSideIPs);
		m_imSource.template 		set_local_ips<refDim>(elemIPs, numElemIPs);
		m_imVectorSource.template 	set_local_ips<refDim>(sideIPs, numSideIPs);
		m_imReactionRate.template 	set_local_ips<refDim>(elemIPs, numElemIPs);
		m_imReaction.template 		set_local_ips<refDim>(elemIPs, numElemIPs);
		m_imMassScale.template 		set_local_ips<refDim>(elemIPs, numElemIPs);
		m_imMass.template 			set_local_ips<refDim>(elemIPs, numElemIPs);

		//	init upwind for element type
		TFVGeom& geo = GeomProvider<TFVGeom>::get();
		if (!m_spConvShape->template set_geometry_type<TFVGeom>(geo))
			UG_THROW("ConvectionDiffusionMP::prep_err_est_elem_loop: "
					 "Cannot init upwind for element type.");

		// store values of shape functions in local IPs
		LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> trialSpace
					= Provider<LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> >::get();

		m_shapeValues.resize(numElemIPs, numSideIPs, trialSpace.num_sh());
		for (size_t ip = 0; ip < numElemIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtElemIP(ip), elemIPs[ip]);
		for (size_t ip = 0; ip < numSideIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtSideIP(ip), sideIPs[ip]);
	}

//	set global positions
	size_t numSideIPs, numElemIPs;
	const MathVector<dim>* sideIPs;
	const MathVector<dim>* elemIPs;

	try
	{
		numSideIPs = err_est_data->num_all_side_ips(roid);
		numElemIPs = err_est_data->num_elem_ips(roid);

		sideIPs = err_est_data->all_side_global_ips(elem, vCornerCoords);
		elemIPs = err_est_data->elem_global_ips(elem, vCornerCoords);
	}
	UG_CATCH_THROW("Global integration points for error estimator cannot be set.");

	m_imDiffusion.			set_global_ips(&sideIPs[0], numSideIPs);
	m_imVelocity.			set_global_ips(&sideIPs[0], numSideIPs);
	
	m_imDarcyW.				set_global_ips(&sideIPs[0], numSideIPs);
	m_imDarcyN.				set_global_ips(&sideIPs[0], numSideIPs);
	
	m_imFlux.				set_global_ips(&sideIPs[0], numSideIPs);
	m_imSource.				set_global_ips(&elemIPs[0], numElemIPs);
	m_imVectorSource.		set_global_ips(&sideIPs[0], numSideIPs);
	m_imReactionRate.		set_global_ips(&elemIPs[0], numElemIPs);
	m_imReaction.			set_global_ips(&elemIPs[0], numElemIPs);
	m_imMassScale.			set_global_ips(&elemIPs[0], numElemIPs);
	m_imMass.				set_global_ips(&elemIPs[0], numElemIPs);
}

//	computes the error estimator contribution (stiffness part) for one element
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (err_est_data->surface_view().get() == NULL) {UG_THROW("Error estimator has NULL surface view.");}
	MultiGrid* pErrEstGrid = (MultiGrid*) (err_est_data->surface_view()->subset_handler()->multi_grid());

//	request geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();


// SIDE TERMS //

//	get the sides of the element
	//	We have to cast elem to a pointer of type SideAndElemErrEstData::elem_type
	//	for the SideAndElemErrEstData::operator() to work properly.
	//	This cannot generally be achieved by casting to TElem*, since this method is also registered for
	//	lower-dimensional types TElem, and must therefore be compilable, even if it is never EVER to be executed.
	//	The way we achieve this here, is by calling associated_elements_sorted() which has an implementation for
	//	all possible types. Whatever comes out of it is of course complete nonsense if (and only if)
	//	SideAndElemErrEstData::elem_type != TElem. To be on the safe side, we throw an error if the number of
	//	entries in the list is not as it should be.

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::side_type>::secure_container side_list;
	pErrEstGrid->associated_elements_sorted(side_list, (TElem*) elem);
	if (side_list.size() != (size_t) ref_elem_type::numSides)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionMP::compute_err_est_elem'");

// 	some help variables
	MathVector<dim> fluxDensity, gradC, normal;

	// FIXME: The computation of the gradient has to be reworked.
	// In the case of P1 shape functions, it is valid. For Q1 shape functions, however,
	// the gradient is not constant (but bilinear) on the element - and along the sides.
	// We cannot use the FVGeom here. Instead, we need to calculate the gradient in each IP!

	// calculate grad u (take grad from first scvf ip (grad u is constant on the entire element))
/*	if (geo.num_scvf() < 1) {UG_THROW("Element has no SCVFs!");}
	const typename TFVGeom::SCVF& scvf = geo.scvf(0);

	VecSet(gradC, 0.0);
	for (size_t j=0; j<m_shapeValues.num_sh(); j++)
		VecScaleAppend(gradC, u(_C_,j), scvf.global_grad(j));*/

	// calculate grad u as average (over scvf)
	VecSet(gradC, 0.0);
	for(size_t ii = 0; ii < geo.num_scvf(); ++ii)
	{
		const typename TFVGeom::SCVF& scvf = geo.scvf(ii);
		for (size_t j=0; j<m_shapeValues.num_sh(); j++)
				VecScaleAppend(gradC, u(_C_,j), scvf.global_grad(j));
	}
	VecScale(gradC, gradC, (1.0/geo.num_scvf()));

	/*VecSet(gradC, 0.0);
	for(size_t ii = 0; ii < geo.num_scv(); ++ii)
	{
		const typename TFVGeom::SCV& scv = geo.scv(ii);
			for (size_t j=0; j<m_shapeValues.num_sh(); j++)
					VecScaleAppend(gradC, u(_C_,j), scv.global_grad(j));
	}
	VecScale(gradC, gradC, (1.0/geo.num_scvf()));
*/

// calculate flux through the sides
	size_t passedIPs = 0;
	for (size_t side=0; side < (size_t) ref_elem_type::numSides; side++)
	{
		// normal on side
		SideNormal<ref_elem_type,dim>(normal, side, vCornerCoords);
		VecNormalize(normal, normal);

		try
		{
			for (size_t sip = 0; sip < err_est_data->num_side_ips(side_list[side]); sip++)
			{
				size_t ip = passedIPs + sip;

				VecSet(fluxDensity, 0.0);

			// diffusion //
				if (m_imDiffusion.data_given())
					MatVecScaleMultAppend(fluxDensity, -1.0, m_imDiffusion[ip], gradC);

			// convection //
				if (m_imVelocity.data_given())
				{
					number val = 0.0;
					for (size_t sh = 0; sh < m_shapeValues.num_sh(); sh++)
						val += u(_C_,sh) * m_shapeValues.shapeAtSideIP(sh,sip);

					VecScaleAppend(fluxDensity, val, m_imVelocity[ip]);
				}

			// general flux //
				if (m_imFlux.data_given())
					VecAppend(fluxDensity, m_imFlux[ip]);

				(*err_est_data)(side_list[side],sip) += scale * VecDot(fluxDensity, normal);
			}

			passedIPs += err_est_data->num_side_ips(side_list[side]);
		}
		UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
				<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
	}

// VOLUME TERMS //

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::elem_type>::secure_container elem_list;
	pErrEstGrid->associated_elements_sorted(elem_list, (TElem*) elem);
	if (elem_list.size() != 1)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionMP::compute_err_est_elem'");

	try
	{
		for (size_t ip = 0; ip < err_est_data->num_elem_ips(elem->reference_object_id()); ip++)
		{
			number total = 0.0;

		// diffusion //	TODO ONLY FOR (PIECEWISE) CONSTANT DIFFUSION TENSOR SO FAR!
		// div(D*grad(c))
		// nothing to do, as c is piecewise linear and div(D*grad(c)) disappears
		// if D is diagonal and c bilinear, this should also vanish (confirm this!)

		// convection // TODO ONLY FOR (PIECEWISE) CONSTANT OR DIVERGENCE-FREE
					  //      VELOCITY FIELDS SO FAR!
		// div(v*c) = div(v)*c + v*grad(c) -- gradC has been calculated above
			if (m_imVelocity.data_given())
				total += VecDot(m_imVelocity[ip], gradC);

		// general flux // TODO ONLY FOR DIVERGENCE-FREE FLUX FIELD SO FAR!
		// nothing to do

		// reaction //
			if (m_imReactionRate.data_given())
			{
				number val = 0.0;
				for (size_t sh = 0; sh < geo.num_sh(); sh++)
					val += u(_C_,sh) * m_shapeValues.shapeAtElemIP(sh,ip);

				total += m_imReactionRate[ip] * val;
			}

			if (m_imReaction.data_given())
			{
				total += m_imReaction[ip];
			}

			(*err_est_data)(elem_list[0],ip) += scale * total;
		}
	}
	UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
			<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
}

//	computes the error estimator contribution (mass part) for one element
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
compute_err_est_M_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
// note: mass parts only enter volume term

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (err_est_data->surface_view().get() == NULL) {UG_THROW("Error estimator has NULL surface view.");}
	MultiGrid* pErrEstGrid = (MultiGrid*) (err_est_data->surface_view()->subset_handler()->multi_grid());

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::elem_type>::secure_container elem_list;
	pErrEstGrid->associated_elements_sorted(elem_list, (TElem*) elem);
	if (elem_list.size() != 1)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionMP::compute_err_est_elem'");

//	request geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop integration points
	try
	{
		for (size_t ip = 0; ip < err_est_data->num_elem_ips(elem->reference_object_id()); ip++)
		{
			number total = 0.0;

		// mass scale //
			if (m_imMassScale.data_given())
			{
				number val = 0.0;
				for (size_t sh = 0; sh < geo.num_sh(); sh++)
					val += u(_C_,sh) * m_shapeValues.shapeAtElemIP(sh,ip);

				total += m_imMassScale[ip] * val;
			}

		// mass //
			if (m_imMass.data_given())
			{
				total += m_imMass[ip];
			}

			(*err_est_data)(elem_list[0],ip) += scale * total;
		}
	}
	UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
			<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
}

//	computes the error estimator contribution (rhs part) for one element
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
compute_err_est_rhs_elem(GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (err_est_data->surface_view().get() == NULL) {UG_THROW("Error estimator has NULL surface view.");}
	MultiGrid* pErrEstGrid = (MultiGrid*) (err_est_data->surface_view()->subset_handler()->multi_grid());

// SIDE TERMS //

//	get the sides of the element
	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::side_type>::secure_container side_list;
	pErrEstGrid->associated_elements_sorted(side_list, (TElem*) elem);
	if (side_list.size() != (size_t) ref_elem_type::numSides)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionMP::compute_err_est_elem'");

// loop sides
	size_t passedIPs = 0;
	for (size_t side = 0; side < (size_t) ref_elem_type::numSides; side++)
	{
		// normal on side
		MathVector<dim> normal;
		SideNormal<ref_elem_type,dim>(normal, side, vCornerCoords);
		VecNormalize(normal, normal);

		try
		{
			for (size_t sip = 0; sip < err_est_data->num_side_ips(side_list[side]); sip++)
			{
				size_t ip = passedIPs + sip;

			// vector source //
				if (m_imVectorSource.data_given())
					(*err_est_data)(side_list[side],sip) += scale * VecDot(m_imVectorSource[ip], normal);
			}

			passedIPs += err_est_data->num_side_ips(side_list[side]);
		}
		UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
				<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
	}

// VOLUME TERMS //

	if (!m_imSource.data_given()) return;

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::elem_type>::secure_container elem_list;
	pErrEstGrid->associated_elements_sorted(elem_list, (TElem*) elem);
	if (elem_list.size() != 1)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionMP::compute_err_est_elem'");

// source //
	try
	{
		for (size_t ip = 0; ip < err_est_data->num_elem_ips(elem->reference_object_id()); ip++)
			(*err_est_data)(elem_list[0],ip) += scale * m_imSource[ip];
	}
	UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
			<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
}

//	postprocesses the loop over all elements of one type in the computation of the error estimator
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
fsh_err_est_elem_loop()
{
//	finish the element loop in the same way as the actual discretization
	this->template fsh_elem_loop<TElem, TFVGeom> ();
};

//    error estimation (end)     ///
// /////////////////////////////////

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
lin_def_velocity(const LocalVector& u,
                 std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                 const size_t nip)
{
// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo, true);
	
	//	Consider capillary trapping	or not
	number u_modified[convShape.num_sh()];
	if ( m_imPermeability.data_given() && m_imMinPd.data_given() ){
		for (size_t i = 0; i < convShape.num_sh(); ++i)
			u_modified[i] = Modify_sw(u(_C_, i), m_imMinPd[i], m_imPermeability[i], m_imPorosity[i]);
	}
	else{
		for (size_t i = 0; i < convShape.num_sh(); ++i)
			u_modified[i] = u(_C_, i);
	}
	
//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

	//	sum up contributions of convection shapes
		MathVector<dim> linDefect;
		VecSet(linDefect, 0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(linDefect, u_modified[sh], convShape.D_vel(ip, sh));

	//	add parts for both sides of scvf
		vvvLinDef[ip][_C_][scvf.from()] += linDefect;
		vvvLinDef[ip][_C_][scvf.to()] -= linDefect;
	}
}

//	computes the linearized defect w.r.t to the darcyW
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
lin_def_darcyW(const LocalVector& u,
                 std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                 const size_t nip)
{
// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo, true);

	//	Consider capillary trapping	or not
	number u_modified[convShape.num_sh()];
	if ( m_imPermeability.data_given() && m_imMinPd.data_given() ){
		for (size_t i = 0; i < convShape.num_sh(); ++i)
			u_modified[i] = Modify_sw(u(_C_, i), m_imMinPd[i], m_imPermeability[i], m_imPorosity[i]);
	}
	else{
		for (size_t i = 0; i < convShape.num_sh(); ++i)
			u_modified[i] = u(_C_, i);
	}
	
//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
		
	//	sum up contributions of convection shapes
		MathVector<dim> linDefect;
		VecSet(linDefect, 0.0);
		
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			VecScaleAppend(linDefect, krw( u_modified[sh] ), convShape.D_vel(ip, sh));
		/*//For Extended Buckley Leverett
		if (m_imPermeability.data_given())
		{
			if (m_imPermeability[0]>5*pow(10,-14))
				VecScaleAppend(linDefect, krw1( u_modified[sh] ), convShape.D_vel(ip, sh));
			else
				VecScaleAppend(linDefect, krw2( u_modified[sh] ), convShape.D_vel(ip, sh));
		}
		*/
		}
		
	//	add parts for both sides of scvf
		vvvLinDef[ip][_C_][scvf.from()] += linDefect;
		vvvLinDef[ip][_C_][scvf.to()] -= linDefect;
	}
}

//	computes the linearized defect w.r.t to the darcyN
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
lin_def_darcyN(const LocalVector& u,
                 std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                 const size_t nip)
{
// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo, true);

//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

	//	sum up contributions of convection shapes
		MathVector<dim> linDefect;
		VecSet(linDefect, 0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			// 	m_imSaturationW is modified before import
			VecScaleAppend(linDefect, krn( m_imSaturationW[sh] ), convShape.D_vel(ip, sh));
		
		/*//For Extended Buckley Leverett
		if (m_imPermeability.data_given())
		{
			if (m_imPermeability[0]>5*pow(10,-14))
				VecScaleAppend(linDefect, krn1( m_imSaturationW[sh] ), convShape.D_vel(ip, sh));
			else
				VecScaleAppend(linDefect, krn2( m_imSaturationW[sh] ), convShape.D_vel(ip, sh));
		}
		*/
		}
	//	add parts for both sides of scvf
		vvvLinDef[ip][_C_][scvf.from()] += linDefect;
		vvvLinDef[ip][_C_][scvf.to()] -= linDefect;
		
	}
}


//	computes the linearized defect w.r.t to the diffusion
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
lin_def_diffusion(const LocalVector& u,
                  std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
                  const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo, true);

	//	Consider capillary trapping	or not
	number u_modified[convShape.num_sh()];
	if ( m_imPermeability.data_given() && m_imMinPd.data_given() ){
		for (size_t i = 0; i < convShape.num_sh(); ++i)
			u_modified[i] = Modify_sw(u(_C_, i), m_imMinPd[i], m_imPermeability[i], m_imPorosity[i]);
	}
	else{
		for (size_t i = 0; i < convShape.num_sh(); ++i)
			u_modified[i] = u(_C_, i);
	}
	
//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

	// 	compute gradient at ip
		MathVector<dim> grad_u;	VecSet(grad_u, 0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(grad_u, u_modified[sh], scvf.global_grad(sh));

	//	compute the lin defect at this ip
		MathMatrix<dim,dim> linDefect;

	//	part coming from -\nabla u * \vec{n}
		for(size_t k=0; k < (size_t)dim; ++k)
			for(size_t j = 0; j < (size_t)dim; ++j)
				linDefect(j,k) = (scvf.normal())[j] * grad_u[k];

	//	add contribution from convection shapes
		if(convShape.non_zero_deriv_diffusion())
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				MatAdd(linDefect, convShape.D_diffusion(ip, sh), u_modified[sh]);

	//	add contributions
		vvvLinDef[ip][_C_][scvf.from()] -= linDefect;
		vvvLinDef[ip][_C_][scvf.to()  ] += linDefect;
	}
}

template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
lin_def_flux(const LocalVector& u,
             std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
             const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

	//	add parts for both sides of scvf
		vvvLinDef[ip][_C_][scvf.from()] += scvf.normal();
		vvvLinDef[ip][_C_][scvf.to()] -= scvf.normal();
	}
}

//	computes the linearized defect w.r.t to the reaction rate
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
lin_def_reaction_rate(const LocalVector& u,
                      std::vector<std::vector<number> > vvvLinDef[],
                      const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	//	Consider capillary trapping	or not
		number u_modified;
		if ( m_imPermeability.data_given() && m_imMinPd.data_given() ){
			u_modified = Modify_sw(u(_C_, co), m_imMinPd[co], m_imPermeability[co], m_imPorosity[co]);
		}
		else{
			u_modified = u(_C_, co);
		}
		
	// 	set lin defect
		vvvLinDef[ip][_C_][co] = u_modified * scv.volume();
	}
}

//	computes the linearized defect w.r.t to the reaction
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
lin_def_reaction(const LocalVector& u,
                 std::vector<std::vector<number> > vvvLinDef[],
                 const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	set lin defect
		vvvLinDef[ip][_C_][co] = scv.volume();
	}
}

//	computes the linearized defect w.r.t to the source
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
lin_def_source(const LocalVector& u,
               std::vector<std::vector<number> > vvvLinDef[],
               const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	set lin defect
		vvvLinDef[ip][_C_][co] = scv.volume();
	}
}

//	computes the linearized defect w.r.t to the vector source
//	(in analogy to velocity)
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
lin_def_vector_source(const LocalVector& u,
                      std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                      const size_t nip)
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

	// loop Sub Control Volumes Faces (SCVF)
	for ( size_t ip = 0; ip < geo.num_scvf(); ++ip ) {
		// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf( ip );

		// add parts for both sides of scvf
		vvvLinDef[ip][_C_][scvf.from()] -= scvf.normal();
		vvvLinDef[ip][_C_][scvf.to()] += scvf.normal();
	}
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
lin_def_mass_scale(const LocalVector& u,
                   std::vector<std::vector<number> > vvvLinDef[],
                   const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
	
	number u_modified[geo.num_scv()];
	if ( m_imPermeability.data_given() && m_imMinPd.data_given() ){
		for (size_t i = 0; i < geo.num_scv(); ++i)
			u_modified[i] = Modify_sw(u(_C_, i), m_imMinPd[i], m_imPermeability[i], m_imPorosity[i]);
	}
	else{
		for (size_t i = 0; i < geo.num_scv(); ++i)
			u_modified[i] = u(_C_, i);
	}
	
// 	loop Sub Control Volumes (SCV)
	for(size_t co = 0; co < geo.num_scv(); ++co)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(co);

	// 	Check associated node
		UG_ASSERT(co == scv.node_id(), "Only one shape per SCV");

	// 	set lin defect
		vvvLinDef[co][_C_][co] = u_modified[co] * scv.volume();
	}
}

//	computes the linearized defect w.r.t to the mass
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
lin_def_mass(const LocalVector& u,
             std::vector<std::vector<number> > vvvLinDef[],
             const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t co = 0; co < geo.num_scv(); ++co)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(co);

	// 	Check associated node
		UG_ASSERT(co == scv.node_id(), "Only one shape per SCV");

	// 	set lin defect
		vvvLinDef[co][_C_][co] = scv.volume();
	}
}

template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
ex_modified_value(number vValue[],
         const MathVector<dim> vGlobIP[],
         number time, int si,
         const LocalVector& u,
         GridObject* elem,
         const MathVector<dim> vCornerCoords[],
         const MathVector<TFVGeom::dim> vLocIP[],
         const size_t nip,
         bool bDeriv,
         std::vector<std::vector<number> > vvvDeriv[])
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::numCorners;
	
	//consider about it !!!!!!! it requires the order of u should be the same as the order of MinPd !!!
	// calculate modified sw
	//	Consider capillary trapping	or not
	number u_modified[numSH];
	if ( m_imPermeability.data_given() && m_imMinPd.data_given() && (nip != 1) ){
		for (size_t i = 0; i < numSH; ++i)
			u_modified[i] = Modify_sw(u(_C_, i), m_imMinPd[i], m_imPermeability[i], m_imPorosity[i]);
	}
	else{
		for (size_t i = 0; i < numSH; ++i)
			u_modified[i] = u(_C_, i);
	}
	
//	MP SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//	Loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
		
		//consider about it !!!!!!! it means the order of vValue is the order of geo.scvf(), the order of u is the order of shape function
		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					vValue[ip] += u_modified[sh] * scvf.shape(sh);
			
		//	compute derivative w.r.t. to unknowns iff needed
			if(bDeriv)
			{
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					vvvDeriv[ip][_C_][sh] = scvf.shape(sh);

				// do not forget that number of DoFs (== vvvDeriv[ip][_C_])
				// might be > scvf.num_sh() in case of hanging nodes!
				size_t ndof = vvvDeriv[ip][_C_].size();
				for (size_t sh = scvf.num_sh(); sh < ndof; ++sh)
					vvvDeriv[ip][_C_][sh] = 0.0;
			}
		}
	}
//	MP SCV ip
	else if(vLocIP == geo.scv_local_ips())
	{
	//	Loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
		// 	Get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

		//	get corner of SCV
			const size_t co = scv.node_id();
			
		//consider about it !!!!!!! it means the order of vValue is the order of geo.scv(), but the order of u is not the same as theirs
		//	solution at ip
			vValue[ip] = u_modified[co];

		//	set derivatives if needed
			if(bDeriv)
			{
				size_t ndof = vvvDeriv[ip][_C_].size();
				for(size_t sh = 0; sh < ndof; ++sh)
					vvvDeriv[ip][_C_][sh] = (sh==co) ? 1.0 : 0.0;
			}
		}
	}
// 	general case
	else
	{
	//	get trial space
		LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	//	storage for shape function at ip
		number vShape[numSH];

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.shapes(vShape, vLocIP[ip]);

		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < numSH; ++sh)
				vValue[ip] += u_modified[sh] * vShape[sh];

		//	compute derivative w.r.t. to unknowns iff needed
		//	\todo: maybe store shapes directly in vvvDeriv
			if(bDeriv)
			{
				for(size_t sh = 0; sh < numSH; ++sh)
					vvvDeriv[ip][_C_][sh] = vShape[sh];

				// beware of hanging nodes!
				size_t ndof = vvvDeriv[ip][_C_].size();
				for (size_t sh = numSH; sh < ndof; ++sh)
					vvvDeriv[ip][_C_][sh] = 0.0;
			}
		}
	}
}

template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
ex_value(number vValue[],
         const MathVector<dim> vGlobIP[],
         number time, int si,
         const LocalVector& u,
         GridObject* elem,
         const MathVector<dim> vCornerCoords[],
         const MathVector<TFVGeom::dim> vLocIP[],
         const size_t nip,
         bool bDeriv,
         std::vector<std::vector<number> > vvvDeriv[])
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::numCorners;
	
//	MP SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//	Loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
		
		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					vValue[ip] += u(_C_, sh) * scvf.shape(sh);
			
		//	compute derivative w.r.t. to unknowns iff needed
			if(bDeriv)
			{
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					vvvDeriv[ip][_C_][sh] = scvf.shape(sh);

				// do not forget that number of DoFs (== vvvDeriv[ip][_C_])
				// might be > scvf.num_sh() in case of hanging nodes!
				size_t ndof = vvvDeriv[ip][_C_].size();
				for (size_t sh = scvf.num_sh(); sh < ndof; ++sh)
					vvvDeriv[ip][_C_][sh] = 0.0;
			}
		}
	}
//	MP SCV ip
	else if(vLocIP == geo.scv_local_ips())
	{
	//	Loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
		// 	Get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

		//	get corner of SCV
			const size_t co = scv.node_id();
			
		//consider about it !!!!!!! it means the order of vValue is the order of geo.scv(), but the order of u is not the same as theirs
		//	solution at ip
			vValue[ip] = u(_C_, co);

		//	set derivatives if needed
			if(bDeriv)
			{
				size_t ndof = vvvDeriv[ip][_C_].size();
				for(size_t sh = 0; sh < ndof; ++sh)
					vvvDeriv[ip][_C_][sh] = (sh==co) ? 1.0 : 0.0;
			}
		}
	}
// 	general case
	else
	{
	//	get trial space
		LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	//	storage for shape function at ip
		number vShape[numSH];

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.shapes(vShape, vLocIP[ip]);

		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < numSH; ++sh)
				vValue[ip] += u(_C_, sh) * vShape[sh];

		//	compute derivative w.r.t. to unknowns iff needed
		//	\todo: maybe store shapes directly in vvvDeriv
			if(bDeriv)
			{
				for(size_t sh = 0; sh < numSH; ++sh)
					vvvDeriv[ip][_C_][sh] = vShape[sh];

				// beware of hanging nodes!
				size_t ndof = vvvDeriv[ip][_C_].size();
				for (size_t sh = numSH; sh < ndof; ++sh)
					vvvDeriv[ip][_C_][sh] = 0.0;
			}
		}
	}
}

// grad_Sw output is needed in Wc equation
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
ex_grad(MathVector<dim> vValue[],
        const MathVector<dim> vGlobIP[],
        number time, int si,
        const LocalVector& u,
        GridObject* elem,
        const MathVector<dim> vCornerCoords[],
        const MathVector<TFVGeom::dim> vLocIP[],
        const size_t nip,
        bool bDeriv,
        std::vector<std::vector<MathVector<dim> > > vvvDeriv[])
{
// 	Get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	reference dimension
	static const int refDim = ref_elem_type::dim;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::numCorners;

	number u_modified[numSH];
	if ( m_imPermeability.data_given() && m_imMinPd.data_given() && (nip != 1) ){
		for (size_t i = 0; i < numSH; ++i)
			u_modified[i] = Modify_sw(u(_C_, i), m_imMinPd[i], m_imPermeability[i], m_imPorosity[i]);
	}
	else{
		for (size_t i = 0; i < numSH; ++i)
			u_modified[i] = u(_C_, i);
	}

//	FV1 SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//	Loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

			VecSet(vValue[ip], 0.0);

			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecScaleAppend(vValue[ip], u_modified[sh], scvf.global_grad(sh));
				//VecScaleAppend(vValue[ip], u(_C_, sh), scvf.global_grad(sh));

			if(bDeriv)
			{
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					vvvDeriv[ip][_C_][sh] = scvf.global_grad(sh);

				// beware of hanging nodes!
				size_t ndof = vvvDeriv[ip][_C_].size();
				for (size_t sh = scvf.num_sh(); sh < ndof; ++sh)
					vvvDeriv[ip][_C_][sh] = 0.0;
			}
		}
	}
// 	general case
	else
	{
	//	get trial space
		LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	//	storage for shape function at ip
		MathVector<refDim> vLocGrad[numSH];
		MathVector<refDim> locGrad;

	//	Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		ReferenceMapping<ref_elem_type, dim> mapping(vCornerCoords);

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.grads(vLocGrad, vLocIP[ip]);

		//	compute grad at ip
			VecSet(locGrad, 0.0);
			for(size_t sh = 0; sh < numSH; ++sh)
				VecScaleAppend(locGrad, u_modified[sh], vLocGrad[sh]);
				//VecScaleAppend(locGrad, u(_C_, sh), vLocGrad[sh]);

		//	compute global grad
			mapping.jacobian_transposed_inverse(JTInv, vLocIP[ip]);
			MatVecMult(vValue[ip], JTInv, locGrad);

		//	compute derivative w.r.t. to unknowns iff needed
			if(bDeriv)
			{
				for(size_t sh = 0; sh < numSH; ++sh)
					MatVecMult(vvvDeriv[ip][_C_][sh], JTInv, vLocGrad[sh]);

				// beware of hanging nodes!
				size_t ndof = vvvDeriv[ip][_C_].size();
				for (size_t sh = numSH; sh < ndof; ++sh)
					vvvDeriv[ip][_C_][sh] = 0.0;
			}
		}
	}
};


// grad_Sw output is needed in Wc equation
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
ex_grad_pd(MathVector<dim> vValue[],
        const MathVector<dim> vGlobIP[],
        number time, int si,
        const LocalVector& u,
        GridObject* elem,
        const MathVector<dim> vCornerCoords[],
        const MathVector<TFVGeom::dim> vLocIP[],
        const size_t nip,
        bool bDeriv,
        std::vector<std::vector<MathVector<dim> > > vvvDeriv[])
{
// 	Get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	reference dimension
	static const int refDim = ref_elem_type::dim;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::numCorners;

	number u_modified[numSH];
	if ( m_imPermeability.data_given() && m_imMinPd.data_given() && (nip != 1) ){
		for (size_t i = 0; i < numSH; ++i)
			u_modified[i] = Pd(m_imPermeability[i], m_imPorosity[i]) * Modify_sw(u(_C_, i), m_imMinPd[i], m_imPermeability[i], m_imPorosity[i]);
	}
	else{
		for (size_t i = 0; i < numSH; ++i)
			u_modified[i] = Pd(m_imPermeability[i], m_imPorosity[i]) * u(_C_, i);
	}

//	FV1 SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//	Loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

			VecSet(vValue[ip], 0.0);

			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecScaleAppend(vValue[ip], u_modified[sh], scvf.global_grad(sh));
				//VecScaleAppend(vValue[ip], u(_C_, sh), scvf.global_grad(sh));

			if(bDeriv)
			{
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					vvvDeriv[ip][_C_][sh] = scvf.global_grad(sh);

				// beware of hanging nodes!
				size_t ndof = vvvDeriv[ip][_C_].size();
				for (size_t sh = scvf.num_sh(); sh < ndof; ++sh)
					vvvDeriv[ip][_C_][sh] = 0.0;
			}
		}
	}
// 	general case
	else
	{
	//	get trial space
		LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	//	storage for shape function at ip
		MathVector<refDim> vLocGrad[numSH];
		MathVector<refDim> locGrad;

	//	Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		ReferenceMapping<ref_elem_type, dim> mapping(vCornerCoords);

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.grads(vLocGrad, vLocIP[ip]);

		//	compute grad at ip
			VecSet(locGrad, 0.0);
			for(size_t sh = 0; sh < numSH; ++sh)
				VecScaleAppend(locGrad, u_modified[sh], vLocGrad[sh]);
				//VecScaleAppend(locGrad, u(_C_, sh), vLocGrad[sh]);

		//	compute global grad
			mapping.jacobian_transposed_inverse(JTInv, vLocIP[ip]);
			MatVecMult(vValue[ip], JTInv, locGrad);

		//	compute derivative w.r.t. to unknowns iff needed
			if(bDeriv)
			{
				for(size_t sh = 0; sh < numSH; ++sh)
					MatVecMult(vvvDeriv[ip][_C_][sh], JTInv, vLocGrad[sh]);

				// beware of hanging nodes!
				size_t ndof = vvvDeriv[ip][_C_].size();
				for (size_t sh = numSH; sh < ndof; ++sh)
					vvvDeriv[ip][_C_][sh] = 0.0;
			}
		}
	}
};


////////////////////////////////////////////////////////////////////////////////
//	upwind
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void ConvectionDiffusionMP<TDomain>::
set_upwind(SmartPtr<IConvectionShapes<dim> > shapes) {m_spConvShape = shapes;}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
const typename ConvectionDiffusionMP<TDomain>::conv_shape_type&
ConvectionDiffusionMP<TDomain>::
get_updated_conv_shapes(const FVGeometryBase& geo, bool compute_deriv)
{
//	compute upwind shapes for transport equation
//	\todo: we should move this computation into the preparation part of the
//			disc, to only compute the shapes once, reusing them several times.

/*
	if(m_imVelocity.data_given() && m_imDarcyW.data_given())//for upwind Wc
	{
		//	get diffusion at ips
		const MathMatrix<dim, dim>* vDiffusion = NULL;
		if(m_imDiffusion.data_given()) vDiffusion = m_imDiffusion.values();
		
		MathVector<dim> Uw = m_imDarcyW[ip];
		VecScaleAppend(Uw, m_imPermeability[0]*Pd(m_imPermeability[0]), m_imVelocity[ip]);

		if(!m_spConvShape->update(&geo, Uw, vDiffusion, compute_deriv))
		{
			UG_LOG("ERROR in 'ConvectionDiffusionMP::get_updated_conv_shapes': "
					"Cannot compute convection shapes.\n");
		}
	}
*/

	if(m_imVelocity.data_given())
	{
	//	get diffusion at ips
		const MathMatrix<dim, dim>* vDiffusion = NULL;
		if(m_imDiffusion.data_given()) vDiffusion = m_imDiffusion.values();

	//	update convection shapes
		if(!m_spConvShape->update(&geo, m_imVelocity.values(), vDiffusion, compute_deriv))
		{
			UG_LOG("ERROR in 'ConvectionDiffusionMP::get_updated_conv_shapes': "
					"Cannot compute convection shapes.\n");
		}
	}
	
	if(m_imDarcyW.data_given())
	{
	//	get diffusion at ips
		const MathMatrix<dim, dim>* vDiffusion = NULL;
		if(m_imDiffusion_Sw.data_given()) vDiffusion = m_imDiffusion_Sw.values();

	//	update convection shapes
		if(!m_spConvShape->update(&geo, m_imDarcyW.values(), vDiffusion, compute_deriv))
		{
			UG_LOG("ERROR in 'ConvectionDiffusionMP::get_updated_conv_shapes': "
					"Cannot compute convection shapes.\n");
		}
	}
	
	if(m_imDarcyN.data_given())
	{
	//	get diffusion at ips
		const MathMatrix<dim, dim>* vDiffusion = NULL;
		if(m_imDiffusion.data_given()) vDiffusion = m_imDiffusion.values();

	//	update convection shapes
		if(!m_spConvShape->update(&geo, m_imDarcyN.values(), vDiffusion, compute_deriv))
		{
			UG_LOG("ERROR in 'ConvectionDiffusionMP::get_updated_conv_shapes': "
					"Cannot compute convection shapes.\n");
		}
	}

//	return a const (!!) reference to the upwind
	return *const_cast<const IConvectionShapes<dim>*>(m_spConvShape.get());
}


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

/**
 * Registers the assembling functions for all the element types
 */
template<typename TDomain>
void ConvectionDiffusionMP<TDomain>::
register_all_funcs(bool bHang)
{
//	list of element types for assembling
	typedef typename domain_traits<dim>::AllElemList AssembleElemList;
	
//	register the local assembling functions
	boost::mpl::for_each<AssembleElemList> (RegisterLocalDiscr (this, bHang));
}

/**
 * Registers the assembling functions for an element type.
 */
template<typename TDomain>
template<typename TElem>
void
ConvectionDiffusionMP<TDomain>::
register_func_for_(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		if(!m_bCondensedFV)
			register_func<TElem, FV1Geometry<TElem, dim> >();
		else
			register_func<TElem, FV1CondensedGeometry<TElem, dim> >();
	}
	else
	{
		if(m_bCondensedFV)
			UG_THROW("ConvectionDiffusionMP: Condensed FV not supported for hanging nodes.");
		register_func<TElem, HFV1Geometry<TElem, dim> >();
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionMP<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;
	static const int refDim = reference_element_traits<TElem>::dim;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_A_expl_elem_fct(id, &T::template add_def_A_expl_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(  id, &T::template add_rhs_elem<TElem, TFVGeom>);

// error estimator parts
	this->set_prep_err_est_elem_loop(id, &T::template prep_err_est_elem_loop<TElem, TFVGeom>);
	this->set_prep_err_est_elem(id, &T::template prep_err_est_elem<TElem, TFVGeom>);
	this->set_compute_err_est_A_elem(id, &T::template compute_err_est_A_elem<TElem, TFVGeom>);
	this->set_compute_err_est_M_elem(id, &T::template compute_err_est_M_elem<TElem, TFVGeom>);
	this->set_compute_err_est_rhs_elem(id, &T::template compute_err_est_rhs_elem<TElem, TFVGeom>);
	this->set_fsh_err_est_elem_loop(id, &T::template fsh_err_est_elem_loop<TElem, TFVGeom>);

//	set computation of linearized defect w.r.t velocity
	m_imDiffusion.set_fct(id, this, &T::template lin_def_diffusion<TElem, TFVGeom>);
	m_imVelocity. set_fct(id, this, &T::template lin_def_velocity<TElem, TFVGeom>);
	
	m_imDarcyW. set_fct(id, this, &T::template lin_def_darcyW<TElem, TFVGeom>);
	m_imDarcyN. set_fct(id, this, &T::template lin_def_darcyN<TElem, TFVGeom>);
	//m_imSaturationW.	set_fct(id, this, &T::template lin_def_saturationW<TElem, TFVGeom>);
	
	m_imFlux.set_fct(id, this, &T::template lin_def_flux<TElem, TFVGeom>);
	m_imReactionRate. set_fct(id, this, &T::template lin_def_reaction_rate<TElem, TFVGeom>);
	m_imReaction. set_fct(id, this, &T::template lin_def_reaction<TElem, TFVGeom>);
	m_imSource.	  set_fct(id, this, &T::template lin_def_source<TElem, TFVGeom>);
	m_imVectorSource.set_fct(id, this, &T::template lin_def_vector_source<TElem, TFVGeom>);
	m_imMassScale.set_fct(id, this, &T::template lin_def_mass_scale<TElem, TFVGeom>);
	m_imMass.	set_fct(id, this, &T::template lin_def_mass<TElem, TFVGeom>);

//	exports
	m_exModifiedValue->		template set_fct<T,refDim>(id, this, &T::template ex_modified_value<TElem, TFVGeom>);
	m_exValue->		template set_fct<T,refDim>(id, this, &T::template ex_value<TElem, TFVGeom>);
	m_exGrad->		template set_fct<T,refDim>(id, this, &T::template ex_grad<TElem, TFVGeom>);
	m_exGrad_pd->		template set_fct<T,refDim>(id, this, &T::template ex_grad_pd<TElem, TFVGeom>);

}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class ConvectionDiffusionMP<Domain1d>;
#endif
#ifdef UG_DIM_2
template class ConvectionDiffusionMP<Domain2d>;
#endif
#ifdef UG_DIM_3
template class ConvectionDiffusionMP<Domain3d>;
#endif

} // end namespace ConvectionDiffusionPlugin
} // namespace ug

