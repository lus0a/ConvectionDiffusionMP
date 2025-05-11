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

static size_t saturated_Node = 439; 

//static number krw_max = 0.9;
//static number krn_max = 0.5;

static number krw_max = 1.0;
static number krn_max = 1.0;

//static number m = 1;
static number m = 1;
static number n = 1;
//static number n = 1;

//static number m = 4;
//static number n = 4;

static number DensityH = 1000;
//static number DensityC = 1000;
//static number DensityC = 1.225;
static number DensityC = 900;

//static number Permeability = 4.845e-13;
//static number ViscosityW = 1e-3;

static number ViscosityW = 0.00047100526797888;
static number ViscosityN = 0.00006019;

//static number ViscosityW = 0.0004734;
//static number ViscosityN = 0.00007755;

//static number ViscosityW = 0.001;
//static number ViscosityN = 0.001;
//static number ViscosityN = 0.0000181;
//static number ViscosityN = 0.00009;

//static number swr = 0.0;
//static number snr = 0.0;
//static number swr = 0.2;
//static number snr = 0.1;

//static number lambda = 4.2;

//static number lambda = 2.0;
static number pc_max = 2000000;
///*
static number krw(number sw, number swr, number lambda)
{
	number sw_eff = (sw - swr) / (1-swr);
	if (sw_eff>0 && sw_eff<1)
		return pow(sw_eff, 2/lambda+3);
	else if (sw_eff<=0)
		return 0;
	else
		return 1;
}

static number D_krw(number sw, number swr, number lambda)
{
	number sw_eff = (sw - swr) / (1-swr);
	if (sw_eff>0 && sw_eff<1)
		return (2/lambda+3) * pow(sw_eff, 2/lambda+2) / (1-swr);
	else if (sw_eff<=0)
		return 0;
	else
		return (2/lambda+3) / (1-swr);
}

static number krn(number sw, number snr, number lambda)
{
    number sw_eff = sw / (1-snr);
	if (sw_eff>0 && sw_eff<1)
		return (1-sw_eff)*(1-sw_eff) * (1-pow(sw_eff, 2/lambda+1));
	else if (sw_eff>=1)
		return 0;
	else
		return 1;
}

static number D_krn(number sw, number snr, number lambda)
{	
	number sw_eff = sw / (1-snr);
	if (sw_eff>0 && sw_eff<1)
		return ( -2*(1-sw_eff)*(1-pow(sw_eff, 2/lambda+1))-(1-sw_eff)*(1-sw_eff)*(2/lambda+1)*pow(sw_eff, 2/lambda) )/ (1-snr);
	else if (sw_eff>=1)
		return 0;
	else
		return -2/(1-snr);
}
//*/

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

/*
static number krw(number sw, number swr, number lambda)
{
	number sw_eff = (sw - swr) / (1-swr);
	if (sw_eff>0 && sw_eff<1)
		return krw_max * pow(sw_eff, m);
	else if (sw_eff<=0)
		return 0;
	else
		return krw_max;
}

static number D_krw(number sw, number swr, number lambda)
{
	number sw_eff = (sw - swr) / (1-swr);
	if (sw_eff>0 && sw_eff<1)
		return krw_max * m * pow(sw_eff, m-1) / (1-swr);
	else if (sw_eff<=0)
		return krw_max * m / (1-swr);
		//return 0;
	else
		return krw_max * m / (1-swr);
}

static number krn(number sw, number snr, number lambda)
{
    number sw_eff = sw / (1-snr);
	if (sw_eff>0 && sw_eff<1)
		return krn_max * pow(1-sw_eff, n);
	else if (sw_eff>=1)
		return 0;
	else
		return krn_max;
}

static number D_krn(number sw, number snr, number lambda)
{	
	number sw_eff = sw / (1-snr);
	if (sw_eff>0 && sw_eff<1)
		return -krn_max * n * pow(1-sw_eff, n-1) / (1-snr);
	else if (sw_eff>=1)
		return -krn_max * n / (1-snr);
		//return 0;
	else
		return -krn_max * n / (1-snr);
}
*/

///*
static number DensityW(number wC, number pN)
{
	pN=pN/10000000;
	//pN = 4.5;
	//int pN=5;
	//number result = (939.3+25.9*2-7.406*2*2);
	//number result = 1000;
	
	number result = (956.9+408.3*wC+1.873*pN-723.1*wC*wC+6.93*wC*pN-0.02482*pN*pN+303*wC*wC*wC+13.51*wC*wC*pN-0.1472*wC*pN*pN-424.2*wC*wC*wC*wC+17.63*wC*wC*wC*pN-0.5314*wC*wC*pN*pN);
	
	//number result = (939.3+497.9*wC+25.9*pN-288.5*wC*wC-214.1*wC*pN-7.406*pN*pN-741.2*wC*wC*wC+94.47*wC*wC*pN+84.17*wC*pN*pN-345.2*wC*wC*wC*wC+510.9*wC*wC*wC*pN-137.3*wC*wC*pN*pN);

	//number result = (939.3+497.9*wC+25.9*2-288.5*wC*wC-214.1*wC*2-7.406*2*2-741.2*wC*wC*wC+94.47*wC*wC*2+84.17*wC*2*2-345.2*wC*wC*wC*wC+510.9*wC*wC*wC*2-137.3*wC*wC*2*2);

	//number result = (1005+307.5*wC-6.98*3
	//	-589.8*wC*wC+52.53*wC*pN); //For Pn in [25,35]Mpa
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

static number D_densityW_X(number wC, number pN)
{
	pN=pN/10000000;
	//pN = 4.5;
	//int pN=5;
	//number result = (939.3+25.9*2-7.406*2*2);
	//number result = 1000;
	
	//number result = 0;
	number result = (408.3-723.1*2*wC+6.93*pN+303*3*wC*wC+13.51*2*wC*pN-0.1472*pN*pN-424.2*4*wC*wC*wC+17.63*3*wC*wC*pN-0.5314*2*wC*pN*pN);
	//number result = (497.9-288.5*2*wC-214.1*pN-741.2*3*wC*wC+94.47*2*wC*pN+84.17*pN*pN-345.2*4*wC*wC*wC+510.9*3*wC*wC*pN-137.3*2*wC*pN*pN);

	//number result = (497.9-288.5*2*wC-214.1*2-741.2*3*wC*wC+94.47*2*wC*2+84.17*2*2-345.2*4*wC*wC*wC+510.9*3*wC*wC*2-137.3*2*wC*2*2);

	//number result = (307.5
	//	-589.8*2*wC+52.53*pN); //For Pn in [25,35]Mpa
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

static number D_densityW_P(number wC, number pN)
{
	pN=pN/10000000;
	
	number result = (1.873+6.93*wC-0.02482*2*pN+13.51*wC*wC-0.1472*wC*2*pN+17.63*wC*wC*wC-0.5314*wC*wC*2*pN)/10000000;
	
	return result;
}


static number DensityN(number pN)
{
	pN=pN/10000000;
	//pN = 4.5;
	//int pN=5;
	//number result = (939.3+497.9+25.9*2-288.5-214.1*2-7.406*2*2-741.2+94.47*2+84.17*2*2-345.2+510.9*2-137.3*2*2);
	//number result = 1000;
	//number result = 1.225;
	number result = (939.3+497.9+25.9*pN-288.5-214.1*pN-7.406*pN*pN-741.2+94.47*pN+84.17*pN*pN-345.2+510.9*pN-137.3*pN*pN);
	
	//number result = (1005+307.5*1-6.98*3
	//	-589.8*1*1+52.53*1*pN); //For Pn in [25,35]Mpa
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

static number D_densityN(number pN)
{
	pN=pN/10000000;
	//pN = 4.5;
	//int pN=5;
	//number result = (939.3+497.9+25.9*2-288.5-214.1*2-7.406*2*2-741.2+94.47*2+84.17*2*2-345.2+510.9*2-137.3*2*2);
	//number result = 1000;
	//number result = 0;
	number result = (25.9-214.1-7.406*2*pN+94.47+84.17*2*pN+510.9-137.3*2*pN)/10000000;
	//number result = 52.53/10000000; //For Pn in [25,35]Mpa
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
static number LeverJ(number sw_eff, number lambda) {
	if (sw_eff>0 && sw_eff<1)
		return pow(sw_eff, -1/lambda);
	else if (sw_eff>=1)
		return 1;
	else
		return 100;
}
// Inverse Leverett J
static number InverseLeverJ(number J, number lambda) {
	return pow(J, -lambda);
}

///*
static number Pd(number perm, number poro)
{
	return 0.01*7.37 * pow(poro/perm,0.43);
	
	//return 0;
	//return 7.37 * pow(poro/perm,0.43);
	
	//return 7.37 * pow(1/perm,0.43);
}

static number Dpc(number sw, number swr, number snr, number lambda, number perm, number poro)
{
	number sw_eff = (sw - swr) / (1-snr-swr);
	number pd = Pd(perm, poro);
	// get sw_max crossponding to pc_max
	number pd0 = 10000;
	number sw_max = pow(pc_max/pd0, -lambda) * (1-snr-swr) + swr;
	number sw_star_eff = (pd < 1/pd0) ? sw_max : pow( (pc_max - pd)/sw_max * (1-snr-swr) * lambda /pd, -1/(1+1/lambda) );
	number sw_star = sw_star_eff * (1-snr-swr) + swr;
	
	if (sw <= sw_star)
		return perm * pd *(-1/lambda) * pow(sw_star_eff, (-1/lambda-1)) / (1-snr-swr);
	else if (sw_eff <= 1)
		return perm * pd *(-1/lambda) * pow(sw_eff, (-1/lambda-1)) / (1-snr-swr);
	else
		return perm * pd *(-1/lambda) / (1-snr-swr);
}

static number D_Dpc(number sw, number swr, number snr, number lambda, number perm, number poro)
{
	number sw_eff = (sw - swr) / (1-snr-swr);
	number pd = Pd(perm, poro);
	// get sw_max crossponding to pc_max
	number pd0 = 10000;
	number sw_max = pow(pc_max/pd0, -lambda) * (1-snr-swr) + swr;
	number sw_star_eff = (pd < 1/pd0) ? sw_max : pow( (pc_max - pd)/sw_max * (1-snr-swr) * lambda /pd, -1/(1+1/lambda) );
	number sw_star = sw_star_eff * (1-snr-swr) + swr;
	
	if (sw <= sw_star)
		return 0;
	else if (sw_eff <= 1)
		return perm * pd *(-1/lambda) * (-1/lambda-1) * pow(sw_eff, (-1/lambda-2)) / pow(1-snr-swr, 2);
	else
		return perm * pd *(-1/lambda) * (-1/lambda-1) / pow(1-snr-swr, 2);
}

static number Dpc(number sw, number swr, number snr, number lambda, number perm, number poro, number pd)
{
	number sw_eff = (sw - swr) / (1-snr-swr);
	// get sw_max crossponding to pc_max
	number pd0 = 10000;
	number sw_max = pow(pc_max/pd0, -lambda) * (1-snr-swr) + swr;
	number sw_star_eff = (pd < 1/pd0) ? sw_max : pow( (pc_max - pd)/sw_max * (1-snr-swr) * lambda /pd, -1/(1+1/lambda) );
	number sw_star = sw_star_eff * (1-snr-swr) + swr;
	
	//pd = 100*pd;
	/*
	sw_star = 0.1;
	for (size_t i = 0; i < 4; ++i)
	{
		sw_star = swr / (lambda*log10(pc_max/pd)+log10(x)-log10(1-Swr)-1);
	}
	*/

	if (sw <= sw_star)
		return perm * pd *(-1/lambda) * pow(sw_star_eff, (-1/lambda-1)) / (1-snr-swr);
	else if (sw_eff <= 1)
		return perm * pd *(-1/lambda) * pow(sw_eff, (-1/lambda-1)) / (1-snr-swr);
	else
		return perm * pd *(-1/lambda) / (1-snr-swr);
}

		
static number D_Dpc(number sw, number swr, number snr, number lambda, number perm, number poro, number pd)
{
	number sw_eff = (sw - swr) / (1-snr-swr);
	// get sw_max crossponding to pc_max
	number pd0 = 10000;
	number sw_max = pow(pc_max/pd0, -lambda) * (1-snr-swr) + swr;
	number sw_star_eff = (pd < 1/pd0) ? sw_max : pow( (pc_max - pd)/sw_max * (1-snr-swr) * lambda /pd, -1/(1+1/lambda) );
	number sw_star = sw_star_eff * (1-snr-swr) + swr;
	

	//pd = 100*pd;
	
	/*
	sw_star = 0.1;
	for (size_t i = 0; i < 4; ++i)
	{
		sw_star = swr / (lambda*log10(pc_max/pd)+log10(x)-log10(1-Swr)-1);
	}
	*/
	
	if (sw <= sw_star)
		return 0;
	else if (sw_eff <= 1)
		return perm * pd *(-1/lambda) * (-1/lambda-1) * pow(sw_eff, (-1/lambda-2)) / pow(1-snr-swr, 2);
	else
		return perm * pd *(-1/lambda) * (-1/lambda-1) / pow(1-snr-swr, 2);
}
		

static number Diffusion_S(size_t FormulationIndex, number sw, number swr, number snr, number lambda, number perm, number poro)
{	//Diffusion_S = kr*Dpc;
	if (FormulationIndex == 1){
		return -krw(sw, swr, lambda)*Dpc(sw, swr, snr, lambda, perm, poro)/ViscosityW;
	}
	else if (FormulationIndex == 2){
		sw = 1-sw;
		return -krn(sw, snr, lambda)*Dpc(sw, swr, snr, lambda, perm, poro)/ViscosityN;
	}
	
}
static number Diffusion_S(size_t FormulationIndex, number sw, number swr, number snr, number lambda, number perm, number poro, number pd)
{	//Diffusion_S = kr*Dpc;
	if (FormulationIndex == 1){
		return -krw(sw, swr, lambda)*Dpc(sw, swr, snr, lambda, perm, poro, pd)/ViscosityW;
	}
	else if (FormulationIndex == 2){
		sw = 1-sw;
		return -krn(sw, snr, lambda)*Dpc(sw, swr, snr, lambda, perm, poro, pd)/ViscosityN;
	}
	
}



static number D_diffusion_S(size_t FormulationIndex, number sw, number swr, number snr, number lambda, number perm, number poro)
{	//D_diffusion_S = dkr*Dpc+kr*D_Dpc;
	if (FormulationIndex == 1){
		return -(D_krw(sw, swr, lambda)*Dpc(sw, swr, snr, lambda, perm, poro)+krw(sw, swr, lambda)*D_Dpc(sw, swr, snr, lambda, perm, poro))/ViscosityW;
	}
	else if (FormulationIndex == 2){
		sw = 1-sw;
		return (D_krn(sw, snr, lambda)*Dpc(sw, swr, snr, lambda, perm, poro)+krn(sw, snr, lambda)*D_Dpc(sw, swr, snr, lambda, perm, poro))/ViscosityN;
	}
	
}
static number D_diffusion_S(size_t FormulationIndex, number sw, number swr, number snr, number lambda, number perm, number poro, number pd)
{	//D_diffusion_S = dkr*Dpc+kr*D_Dpc;
	if (FormulationIndex == 1){
		return -(D_krw(sw, swr, lambda)*Dpc(sw, swr, snr, lambda, perm, poro, pd)+krw(sw, swr, lambda)*D_Dpc(sw, swr, snr, lambda, perm, poro, pd))/ViscosityW;
	}
	else if (FormulationIndex == 2){
		sw = 1-sw;
		return (D_krn(sw, snr, lambda)*Dpc(sw, swr, snr, lambda, perm, poro, pd)+krn(sw, snr, lambda)*D_Dpc(sw, swr, snr, lambda, perm, poro, pd))/ViscosityN;
	}
	
}
//*/

/*
static number Pd(number perm, number poro)
{
	//return 1000;
	return 5000;
	//return 0.10;
	//return 7.37 * pow(1/perm,0.43);
}
static number D_diffusion_Sw(number sw, number perm, number poro)
{
	number sw_eff = (sw - swr) / (1-snr-swr);
	if (sw_eff>0 && sw_eff<1)
		return perm/ViscosityW *Pd(perm, poro)*(-1/lambda)* (1/lambda+2) * pow(sw_eff, (1/lambda+1)) / pow(1-snr-swr, 2);
	else if (sw_eff>=1)
		return perm/ViscosityW *Pd(perm, poro)*(-1/lambda)* (1/lambda+2) / pow(1-snr-swr, 2);
	else
		return 0;
		
}
*/

static number Modify_s(size_t FormulationIndex, number sw, number swr, number snr, number lambda, number minPd, number minSwr, number minSnr, number minLambda, number perm, number poro)
{	
	if (FormulationIndex == 2){
		sw = 1-sw;
		number elePd = Pd(perm, poro);

		number refsw_eff = (sw - minSwr) / (1-minSnr-minSwr);
		if ( elePd >= minPd*LeverJ(refsw_eff, minLambda) ){
			return snr;
		}
		refsw_eff = InverseLeverJ( minPd*LeverJ(refsw_eff, minLambda)/elePd, lambda );	
		return 1-(refsw_eff*(1-snr-swr)+swr);
	}
	
	number elePd = Pd(perm, poro);

	number refsw_eff = (sw - minSwr) / (1-minSnr-minSwr);
	if ( elePd >= minPd*LeverJ(refsw_eff, minLambda) ){
		return 1-snr;
	}
	refsw_eff = InverseLeverJ( minPd*LeverJ(refsw_eff, minLambda)/elePd, lambda );	
	return refsw_eff*(1-snr-swr)+swr;
}

static number Modify_s(size_t FormulationIndex, number sw, number swr, number snr, number lambda, number minPd, number minSwr, number minSnr, number minLambda, number pd)
{	
	if (FormulationIndex == 2){
		sw = 1-sw;
		number elePd = pd;

		number refsw_eff = (sw - minSwr) / (1-minSnr-minSwr);
		
		//test pd
		//if (elePd >10*minPd) return snr;		
		
		if ( elePd >= minPd*LeverJ(refsw_eff, minLambda) ){
			return snr;
		}

		refsw_eff = pow( minPd*pow(refsw_eff, -1/minLambda)/elePd, -lambda );	
		return 1-(refsw_eff*(1-snr-swr)+swr);
	}
	
	number elePd = pd;

	number refsw_eff = (sw - minSwr) / (1-minSnr-minSwr);
	if ( elePd >= minPd*LeverJ(refsw_eff, minLambda) ){
		return 1-snr;
	}
	refsw_eff = InverseLeverJ( minPd*LeverJ(refsw_eff, minLambda)/elePd, lambda );	
	return refsw_eff*(1-snr-swr)+swr;
}

template<typename TDomain>
ConvectionDiffusionMP<TDomain>::
ConvectionDiffusionMP(const char* functions, const char* subsets, const size_t FormulationIndex)
 : ConvectionDiffusionBase<TDomain>(functions,subsets),
   m_spConvShape(new ConvectionShapesNoUpwind<dim>),
   m_bNonRegularGrid(false), m_bCondensedFV(false), iFormulationIndex(FormulationIndex)
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
	
	if(iFormulationIndex == 0)
		UG_THROW("Please specify the FormulationIndex, ConvectionDiffusionMP(functions, subsets, FormulationIndex), 1:(Sw,Pn); 2:(Sn,Pw)");
		
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
		m_imModifiedSaturationW.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imUpwindSaturationW.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imUpwindModifiedSaturationW.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imDiffusion_Sw.template 	set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imMassFractionWc.template set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imPressurePn.template 	set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		
		m_imPermeability.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imPorosity.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imEntryPressure.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imResidualAqueous.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imResidualCarbonic.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imBrooksCoreyNumber.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imMinPd.template 			set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imMinSwr.template 			set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imMinSnr.template 			set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imMinLambda.template 			set_local_ips<refDim>(vSCVip,numSCVip, false);
		
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
	m_imModifiedSaturationW.		set_global_ips(vSCVip, numSCVip);
	m_imUpwindSaturationW.	set_global_ips(vSCVip, numSCVip);
	m_imUpwindModifiedSaturationW.	set_global_ips(vSCVip, numSCVip);
	m_imDiffusion_Sw.		set_global_ips(vSCVFip, numSCVFip);
	m_imMassFractionWc.		set_global_ips(vSCVFip, numSCVFip);
	m_imPressurePn.			set_global_ips(vSCVFip, numSCVFip);
	
	m_imPermeability.		set_global_ips(vSCVip, numSCVip);
	m_imPorosity.			set_global_ips(vSCVip, numSCVip);
	m_imEntryPressure.		set_global_ips(vSCVip, numSCVip);
	m_imResidualAqueous.	set_global_ips(vSCVip, numSCVip);
	m_imResidualCarbonic.	set_global_ips(vSCVip, numSCVip);
	m_imBrooksCoreyNumber.	set_global_ips(vSCVip, numSCVip);
	
	
	
	m_imMinPd.				set_global_ips(vSCVip, numSCVip);
	m_imMinSwr.				set_global_ips(vSCVip, numSCVip);
	m_imMinSnr.				set_global_ips(vSCVip, numSCVip);
	m_imMinLambda.			set_global_ips(vSCVip, numSCVip);
	
	
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
		if (m_imEntryPressure.data_given())
			for (size_t i = 0; i < convShape.num_sh(); ++i){
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imEntryPressure[i]);
			}
		else
			for (size_t i = 0; i < convShape.num_sh(); ++i){
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imPermeability[i], m_imPorosity[i]);
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
		// Diffusive_Sw Term--There are two case: 1. in Sw equation or 2. in Wc equation. For case2, there is no contribution for jac matrix.
		////////////////////////////////////////////////////
		
			// For case1, since Sw is given in case2
			if(m_imDiffusion_Sw.data_given() && !m_imModifiedSaturationW.data_given() && m_imPermeability.data_given())
			{
				#ifdef UG_ENABLE_DEBUG_LOGS
				//	DID_CONV_DIFF_MP
				number D_diff_flux_Sw_sum = 0.0;
				#endif

			// 	loop shape functions
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
					MathVector<dim> sum; VecSet(sum, 0.0);
					number s_integ = 0.0;
					MathVector<dim> Dgrad2;
					
					for (size_t i = 0; i < scvf.num_sh(); ++i)
					{
						VecScaleAppend(sum, u_modified[i], scvf.global_grad(i));
						s_integ += u_modified[i]*scvf.shape(i);
					}
					VecScale(sum, sum, scvf.shape(sh));
					
					MathMatrix<dim,dim> Real_Diffusion_Sw;
					if (m_imEntryPressure.data_given())
						MatScale(Real_Diffusion_Sw, Diffusion_S(iFormulationIndex, s_integ, m_imResidualAqueous[0], m_imResidualCarbonic[0], m_imBrooksCoreyNumber[0], m_imPermeability[0], m_imPorosity[0], m_imEntryPressure[0])
, m_imDiffusion_Sw[ip]);
					else
						MatScale(Real_Diffusion_Sw, Diffusion_S(iFormulationIndex, s_integ, m_imResidualAqueous[0], m_imResidualCarbonic[0], m_imBrooksCoreyNumber[0], m_imPermeability[0], m_imPorosity[0])
, m_imDiffusion_Sw[ip]);
				
				// 	Compute Diffusion Tensor times Gradient
					MatVecMult(Dgrad, Real_Diffusion_Sw, scvf.global_grad(sh));
					
				//  Compute D_diffusion * sum( s_i*grad(lambda_i))

					MathMatrix<dim,dim> D_diffusion_S_Mat;
					for (size_t j = 0; j < dim; ++j)
						for (size_t k = 0; k < dim; ++k)
						{
							if (j == k)
							{	
								number D_diffusion_S_value;
								if (m_imEntryPressure.data_given())
									D_diffusion_S_value = D_diffusion_S(iFormulationIndex, s_integ, m_imResidualAqueous[0], m_imResidualCarbonic[0], m_imBrooksCoreyNumber[0], m_imPermeability[0], m_imPorosity[0], m_imEntryPressure[0]);
								else
									D_diffusion_S_value = D_diffusion_S(iFormulationIndex, s_integ, m_imResidualAqueous[0], m_imResidualCarbonic[0], m_imBrooksCoreyNumber[0], m_imPermeability[0], m_imPorosity[0]);
								
								if (m_imMassFractionWc.data_given() && m_imPressurePn.data_given())
								{
									if (iFormulationIndex == 1)
										D_diffusion_S_Mat(j,k) = DensityW(m_imMassFractionWc[ip], m_imPressurePn[ip]) * D_diffusion_S_value;
									else if (iFormulationIndex == 2)
										D_diffusion_S_Mat(j,k) = DensityN(m_imPressurePn[ip]) * D_diffusion_S_value;
									//D_diffusion_Sw_Mat(j,k) = DensityW(m_imMassFractionWc[ip]) * D_diffusion_Sw(sw_integ, m_imPermeability[0], m_imPorosity[0]);
								}
								else if (m_imPressurePn.data_given())
								{
									if (iFormulationIndex == 1)
										D_diffusion_S_Mat(j,k) = DensityW(0, m_imPressurePn[ip]) * D_diffusion_S_value;
									else if (iFormulationIndex == 2)
										D_diffusion_S_Mat(j,k) = DensityN(m_imPressurePn[ip]) * D_diffusion_S_value;
									//D_diffusion_Sw_Mat(j,k) = DensityW(0) * D_diffusion_Sw(sw_integ, m_imPermeability[0], m_imPorosity[0]);
								
								}
								else //immiscible & incompressible, densityW=1000
								{
									if (iFormulationIndex == 1)
										D_diffusion_S_Mat(j,k) = DensityH * D_diffusion_S_value;
									else if (iFormulationIndex == 2)
										D_diffusion_S_Mat(j,k) = DensityC * D_diffusion_S_value;
								}
							}
							else
								D_diffusion_S_Mat(j,k) = 0.0;
							//Shuai debug , here the D_diffusion_S_Mat is set to 0 in the SnPw formulation
							//if (iFormulationIndex == 2)
							//	D_diffusion_S_Mat(j,k) = 0.0;
							//debug end
						}
					MatVecMult(Dgrad2, D_diffusion_S_Mat, sum);
					
					VecAdd(Dgrad, Dgrad, Dgrad2);
					
					//Shuai debug , here the D_diffusion_S_Mat is set to 0 in the SnPw formulation
					//	if (iFormulationIndex == 2)
					//		VecSet(Dgrad, 0.0);
					//debug end
					
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
					//const number D_conv_flux = convShape(ip, sh);
					const number D_conv_flux = (D_densityW_X(u_modified[sh], m_imPressurePn[ip])/DensityW(u_modified[sh], m_imPressurePn[ip])*u_modified[sh]+1)*convShape(ip, sh);
					
					//const number D_conv_flux = krw( m_imSaturationW[sh] ) * convShape(ip, sh);
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
					number D_kr_value = 0.0;
					if (iFormulationIndex == 1)
						D_kr_value = D_krw( u_modified[sh], m_imResidualAqueous[0], m_imBrooksCoreyNumber[0] );
					else if (iFormulationIndex == 2)
						D_kr_value = -D_krn( 1-u_modified[sh], m_imResidualCarbonic[0], m_imBrooksCoreyNumber[0] );
					
					/*
					//Shuai debug
					if ( u_modified[sh]>0.999 )//For each ip, find the saturated node
					{ 	
						//check the node is related to the subcontrol faces(from or to)
						//the upwind node must be from or to
						// 1. avoid inflow to the (downwind && saturated) node
						// 2. avoid outflow from the (upwind && saturated) node
						
						//the saturated node is inside the subcontrolvolume //flow in, the upwind node is scvf.to()
						if ( ( sh==scvf.from() ) && ( convShape(ip, scvf.to())<0 ) ) 
						{
								D_kr_value = 0;
						}
						//the saturated node is outside the subcontrolvolume //flow out, the upwind node is scvf.from()
						else if ( ( sh==scvf.to() ) && ( convShape(ip, scvf.from())>0 ) )
						{
								D_kr_value = 0;	
						}
					}
					//Shuai debug end
					*/
					
					const number D_conv_flux = D_kr_value * convShape(ip, sh);
					
				//	Add flux term to local matrix
					J(_C_, scvf.from(), _C_, sh) += D_conv_flux;
					J(_C_, scvf.to(),   _C_, sh) -= D_conv_flux;
				}
			}
			
			if( m_imDarcyN.data_given() && m_imMassFractionWc.data_given() )
			{	
				
				for(size_t sh = 0; sh < convShape.num_sh(); ++sh)
				{
					number D_kr_value = 0.0;
					if (iFormulationIndex == 1)
						D_kr_value = D_densityN(u_modified[sh])/DensityN(u_modified[sh]);
					else if (iFormulationIndex == 2)
						D_kr_value = D_densityW_P(m_imMassFractionWc[ip], u_modified[sh])/DensityW(m_imMassFractionWc[ip], u_modified[sh]);
					
					/*
					//Shuai debug
					if ( u_modified[sh]>0.999 )//For each ip, find the saturated node
					{ 	
						//check the node is related to the subcontrol faces(from or to)
						//the upwind node must be from or to
						// 1. avoid inflow to the (downwind && saturated) node
						// 2. avoid outflow from the (upwind && saturated) node
						
						//the saturated node is inside the subcontrolvolume //flow in, the upwind node is scvf.to()
						if ( ( sh==scvf.from() ) && ( convShape(ip, scvf.to())<0 ) ) 
						{
								D_kr_value = 0;
						}
						//the saturated node is outside the subcontrolvolume //flow out, the upwind node is scvf.from()
						else if ( ( sh==scvf.to() ) && ( convShape(ip, scvf.from())>0 ) )
						{
								D_kr_value = 0;	
						}
					}
					//Shuai debug end
					*/
					
					const number D_conv_flux = D_kr_value * convShape(ip, sh);
					
	
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
		if (m_imEntryPressure.data_given())
			u_modified = Modify_s(iFormulationIndex, u(_C_, co), m_imResidualAqueous[co], m_imResidualCarbonic[co], m_imBrooksCoreyNumber[co], m_imMinPd[co], m_imMinSwr[co], m_imMinSnr[co], m_imMinLambda[co], m_imEntryPressure[co]);
		else
			u_modified = Modify_s(iFormulationIndex, u(_C_, co), m_imResidualAqueous[co], m_imResidualCarbonic[co], m_imBrooksCoreyNumber[co], m_imMinPd[co], m_imMinSwr[co], m_imMinSnr[co], m_imMinLambda[co], m_imPermeability[co], m_imPorosity[co]);
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
	if ( m_imPermeability.data_given() && m_imMinPd.data_given() && !m_imModifiedSaturationW.data_given() ){
		if (m_imEntryPressure.data_given())
			for (size_t i = 0; i < convShape.num_sh(); ++i){
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imEntryPressure[i]);
			}
		else
			for (size_t i = 0; i < convShape.num_sh(); ++i){
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imPermeability[i], m_imPorosity[i]);
			}
	}
	else{
		for (size_t i = 0; i < convShape.num_sh(); ++i)
			u_modified[i] = u(_C_, i);
	}
	

	if(m_imDiffusion.data_given() || m_imDiffusion_Sw.data_given() || m_imVelocity.data_given() || m_imFlux.data_given() || m_imDarcyW.data_given() || (m_imDarcyN.data_given() && m_imModifiedSaturationW.data_given()))
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
		
		/*
		    // In Wc equation
			if(m_imDiffusion_Sw.data_given() && m_imSaturationW.data_given())
			{
			//	to compute D \nabla c
				MathVector<dim> Dgrad_c, grad_c;

			// 	compute Sw gradient and shape at ip
				VecSet(grad_c, 0.0);
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecScaleAppend(grad_c, m_imSaturationW[sh], scvf.global_grad(sh));
				
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
				
				number conv_flux = 0.0;
				for(size_t sh = 0; sh < convShape.num_sh(); ++sh)
				{
					////For Extended Buckley Leverett
					//if (m_imPermeability.data_given())
					//{
					//	if (m_imPermeability[0]>5*pow(10,-14))
					//		conv_flux += krw1( u_modified[sh] ) * convShape(ip, sh);
					//	else
					//		conv_flux += krw2( u_modified[sh] ) * convShape(ip, sh);
					//}
					//
					
					conv_flux += krw( u_modified[sh] ) * convShape(ip, sh);
						
				}
				
				d(_C_, scvf.from()) += conv_flux;
				d(_C_, scvf.to()  ) -= conv_flux;
			}
			
		*/	
			// In Sw equation
			if(m_imDiffusion_Sw.data_given() && !m_imModifiedSaturationW.data_given())
			{
			//	to compute D \nabla c
				MathVector<dim> Dgrad_c, grad_c;

			// 	compute gradient and shape at ip
				VecSet(grad_c, 0.0);
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecScaleAppend(grad_c, u_modified[sh], scvf.global_grad(sh));
				
				number s_integ = 0.0;
				for (size_t i = 0; i < scvf.num_sh(); ++i)
				{
					s_integ += u_modified[i]*scvf.shape(i);
				}
				
				MathMatrix<dim,dim> Real_Diffusion_Sw;
				if (m_imPermeability.data_given())
				{
					if (m_imEntryPressure.data_given())
						MatScale(Real_Diffusion_Sw, Diffusion_S(iFormulationIndex, s_integ, m_imResidualAqueous[0], m_imResidualCarbonic[0], m_imBrooksCoreyNumber[0], m_imPermeability[0], m_imPorosity[0], m_imEntryPressure[0])
, m_imDiffusion_Sw[ip]);
					else
						MatScale(Real_Diffusion_Sw, Diffusion_S(iFormulationIndex, s_integ, m_imResidualAqueous[0], m_imResidualCarbonic[0], m_imBrooksCoreyNumber[0], m_imPermeability[0], m_imPorosity[0])
, m_imDiffusion_Sw[ip]);
				}
				else
					MatScale(Real_Diffusion_Sw, 1, m_imDiffusion_Sw[ip]);
			//	scale by diffusion tensor
				MatVecMult(Dgrad_c, Real_Diffusion_Sw, grad_c);

			//Shuai debug
			//	if (iFormulationIndex == 2)
			//		VecSet(Dgrad_c, 0.0);
			//debug end


			// 	Compute flux
				const number diff_flux_Sw = VecDot(Dgrad_c, scvf.normal());
				
			// 	Add to local defect
				d(_C_, scvf.from()) -= diff_flux_Sw;
				d(_C_, scvf.to()  ) += diff_flux_Sw;
			}

		/////////////////////////////////////////////////////
		// Convective Term
		/////////////////////////////////////////////////////
			if(m_imVelocity.data_given() && m_imModifiedSaturationW.data_given())
			{
			//	sum up convective flux using convection shapes
				number conv_flux = 0.0;
				for(size_t sh = 0; sh < convShape.num_sh(); ++sh)
				{
					//conv_flux += krw( m_imSaturationW[sh] ) * u_modified[sh] * convShape(ip, sh);
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
				bool saturated = false;
				
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
					
					///*
					//Shuai debug
					size_t xx = ((TElem*) elem)->vertex(sh)->grid_data_index();
					// find node
					//if ( u(_C_, sh)>0.65 &&  m_imEntryPressure[sh]>200)
					//{
					//	int sss;
					//	sss=xx;
					//}
					if ( u(_C_, sh)>0.755 && xx==saturated_Node )
					//if ( (u(_C_, sh)>0.69) && (m_imEntryPressure[0]>1999) )
					{
						saturated = true;
						//saturated_Node = xx;
						//UG_LOG("BadNode"<<saturated_Node);
					}
						//UG_THROW("DataImport::set_roid: Setting unknown ReferenceObjectId.");
						//
					//UG_ASSERT( u(_C_, sh)<0.01 || m_imEntryPressure[0]<9, "BadNode "<<xx<<' '<<u(_C_, sh)<<' '<<m_imEntryPressure[0]);
						
					
						
					/*
					if ( u(_C_, sh)>0.999 || xx==saturated_Node )//For each ip, find the saturated node
					{ 	
						//check the node is related to the subcontrol faces(from or to)
						//the upwind node must be from or to
						// 1. avoid inflow to the (downwind && saturated) node
						// 2. avoid outflow from the (upwind && saturated) node
						
						if (saturated_Node==-1)
						{
							saturated_Node = xx;
						}
						if ( ( (sh==scvf.from()) || (sh==scvf.to()) ) && (xx==saturated_Node) )
						{
								saturated = true;
						}
						
						
						//the saturated node is inside the subcontrolvolume //flow in, the upwind node is scvf.to()
						if ( ( sh==scvf.from() ) && ( convShape(ip, scvf.to())<0 ) ) 
						{
								saturated = true;
						}
						//the saturated node is outside the subcontrolvolume //flow out, the upwind node is scvf.from()
						else if ( ( sh==scvf.to() ) && ( convShape(ip, scvf.from())>0 ) )
						{
							saturated = true;
						}
						
					}
					//Shuai debug end
					*/
					
					if (iFormulationIndex == 1)
						conv_flux += krw( u(_C_, sh), m_imMinSwr[sh], m_imMinLambda[sh] ) * convShape(ip, sh);
					else if (iFormulationIndex == 2)
						conv_flux += krn( 1-u_modified[sh], m_imResidualCarbonic[0], m_imBrooksCoreyNumber[0] ) * convShape(ip, sh);

	
				}
				///*
				if (saturated){
					
					size_t node_from = ((TElem*) elem)->vertex(scvf.from())->grid_data_index();
					size_t node_to = ((TElem*) elem)->vertex(scvf.to())->grid_data_index();
					size_t ele_index = elem->grid_data_index();
					UG_LOGN("             "<<"at Ip# "<<ip<<" "<<"Normal: node# "<<node_from<<"--> node# "<<node_to);
					
					double Sum_convShape=0.0;
					size_t upwind_node=scvf.from();
					for (size_t i=0; i<convShape.num_sh(); i++)
					{
						Sum_convShape += convShape(ip, i);
						if (fabs(convShape(ip, i))>0.000000000001)
							upwind_node = i;
					}
					Sum_convShape = (fabs(Sum_convShape)>0)? Sum_convShape : 1;
					size_t node_upwind = ((TElem*) elem)->vertex(upwind_node)->grid_data_index();
					
					UG_LOGN("  Flux of Sn-Eq. : "<<conv_flux<<"  krn : "<<conv_flux/Sum_convShape<<"   Modified_Sn : "<<u_modified[upwind_node]<<"   Sn : "<<u(_C_, upwind_node)<<"       Upwind_node : "<<node_upwind<<'\n');

					//conv_flux =0;
				}
				//*/
				
			//  add to local defect
				d(_C_, scvf.from()) += conv_flux;
				d(_C_, scvf.to()  ) -= conv_flux;
			}


		//Flux Term - krn(sW)*DarcyN
			if(m_imDarcyN.data_given() && m_imModifiedSaturationW.data_given())
			{
			//	sum up convective flux using convection shapes
				number conv_flux = 0.0;
				
				bool saturated = false;
				
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
					
					///*
					//Shuai debug
					size_t xx = ((TElem*) elem)->vertex(sh)->grid_data_index();
					
					if ( (1-m_imSaturationW[sh])>0.755 && xx==saturated_Node )
					//if ( (1-m_imSaturationW[sh])>0.65 && m_imEntryPressure[0]>99999 )
						saturated = true;
					/*
					if ( (1-m_imUpwindModifiedSaturationW[sh])>0.999 || xx==saturated_Node )//For each ip, find the saturated node
					{ 	
						//check the node is related to the subcontrol faces(from or to)
						//the upwind node must be from or to
						// 1. avoid inflow to the (downwind && saturated) node
						// 2. avoid outflow from the (upwind && saturated) node
						if (saturated_Node==-1)
						{
							saturated_Node = xx;
						}
						if ( ( (sh==scvf.from()) || (sh==scvf.to()) ) && (xx==saturated_Node) )
						{
								saturated = true;
						}
					*/	
					/*
						//the saturated node is inside the subcontrolvolume //flow in, the upwind node is scvf.to()
						if ( ( sh==scvf.from() ) && ( convShape(ip, scvf.to())<0 ) ) 
						{
								saturated = true;
						}
						//the saturated node is outside the subcontrolvolume //flow out, the upwind node is scvf.from()
						else if ( ( sh==scvf.to() ) && ( convShape(ip, scvf.from())>0 ) )
						{
							saturated = true;
						}
						
					}
					*/
					
					
					// Nothing to change, since m_imSaturationW is modified before import
					
					//Shuai debug
					/*
					if (iFormulationIndex == 1)
						conv_flux += krn( m_imUpwindModifiedSaturationW[sh], m_imResidualCarbonic[0], m_imBrooksCoreyNumber[0] ) * convShape(ip, sh);
					else if (iFormulationIndex == 2)
						conv_flux += krw( m_imUpwindModifiedSaturationW[sh], m_imResidualAqueous[0], m_imBrooksCoreyNumber[0] ) * convShape(ip, sh);
					*/
					
					if (iFormulationIndex == 1)
						conv_flux += krn( m_imModifiedSaturationW[sh], m_imResidualCarbonic[0], m_imBrooksCoreyNumber[0] ) * convShape(ip, sh);
					else if (iFormulationIndex == 2)
						conv_flux += krw( m_imSaturationW[sh], m_imMinSwr[sh], m_imMinLambda[sh] ) * convShape(ip, sh);
						//conv_flux += krw( m_imUpwindSaturationW[sh], m_imResidualAqueous[0], m_imBrooksCoreyNumber[0] ) * convShape(ip, sh);
					//Shuai debug end 
					
					
					
					//conv_flux += krn( m_imSaturationW[sh] ) * convShape(ip, sh);
					
					/*
					//Shuai debug
					if ( (1-m_imSaturationW[sh])>0.999 )//For each ip, find the saturated node
					{ 	
						//check the node is related to the subcontrol faces(from or to)
						//the upwind node must be from or to
						// 1. avoid inflow to the (downwind && saturated) node
						// 2. avoid outflow from the (upwind && saturated) node
						
						//the saturated node is inside the subcontrolvolume //flow in, the upwind node is scvf.to()
						if ( ( sh==scvf.from() ) && ( convShape(ip, scvf.to())<pow(10,-10) ) ) 
						{
								conv_flux = 0;
						}
						//the saturated node is outside the subcontrolvolume //flow out, the upwind node is scvf.from()
						else if ( ( sh==scvf.to() ) && ( convShape(ip, scvf.from())>-pow(10,-10) ) )
						{
								conv_flux = 0;	
						}
					}
					//Shuai debug end
					*/
				}
				
				///*
				if (saturated){
					
					size_t node_from = ((TElem*) elem)->vertex(scvf.from())->grid_data_index();
					size_t node_to = ((TElem*) elem)->vertex(scvf.to())->grid_data_index();
					size_t ele_index = elem->grid_data_index();
					UG_LOGN("             "<<"at Ip# "<<ip<<" "<<"Normal: node# "<<node_from<<"--> node# "<<node_to);
					
					double Sum_convShape=0.0;
					size_t upwind_node=scvf.from();
					for (size_t i=0; i<convShape.num_sh(); i++)
					{
						Sum_convShape += convShape(ip, i);
						if (fabs(convShape(ip, i))>0.000000000001)
							upwind_node = i;
					}
					Sum_convShape = (fabs(Sum_convShape)>0)? Sum_convShape : 1;
					
					size_t node_upwind = ((TElem*) elem)->vertex(upwind_node)->grid_data_index();
					UG_LOGN("  Flux of Pw-Eq. : "<<conv_flux<<"  krw : "<<conv_flux/Sum_convShape<<"   UpwindModified_Sn : "<<1-m_imModifiedSaturationW[upwind_node]<<"   UpwindSn : "<<1-m_imSaturationW[upwind_node]<<"       Upwind_node : "<<node_upwind<<"     Swr_upwindnode : "<<m_imMinSwr[upwind_node]<<"     Swr_ele : "<<m_imResidualAqueous[0]<<'\n');

					//conv_flux =0;
				}	
				//*/
					
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
				if (m_imEntryPressure.data_given())
					u_modified = Modify_s(iFormulationIndex, u(_C_, co), m_imResidualAqueous[co], m_imResidualCarbonic[co], m_imBrooksCoreyNumber[co], m_imMinPd[co], m_imMinSwr[co], m_imMinSnr[co], m_imMinLambda[co], m_imEntryPressure[co]);
				else
					u_modified = Modify_s(iFormulationIndex, u(_C_, co), m_imResidualAqueous[co], m_imResidualCarbonic[co], m_imBrooksCoreyNumber[co], m_imMinPd[co], m_imMinSwr[co], m_imMinSnr[co], m_imMinLambda[co], m_imPermeability[co], m_imPorosity[co]);
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
			if (m_imEntryPressure.data_given())
				u_modified = Modify_s(iFormulationIndex, u(_C_, co), m_imResidualAqueous[co], m_imResidualCarbonic[co], m_imBrooksCoreyNumber[co], m_imMinPd[co], m_imMinSwr[co], m_imMinSnr[co], m_imMinLambda[co], m_imEntryPressure[co]);
			else
				u_modified = Modify_s(iFormulationIndex, u(_C_, co), m_imResidualAqueous[co], m_imResidualCarbonic[co], m_imBrooksCoreyNumber[co], m_imMinPd[co], m_imMinSwr[co], m_imMinSnr[co], m_imMinLambda[co], m_imPermeability[co], m_imPorosity[co]);	
		}
		else{
			u_modified = u(_C_, co);
		}
		
		//	mass value
		number val = 0.0;

	//	multiply by scaling
		if(m_imMassScale.data_given())
			val += m_imMassScale[ip] * u_modified;
		
		/*
		
		bool saturated = false;		
								
		size_t xx = ((TElem*) elem)->vertex(ip)->grid_data_index();
		
		if ( ( m_imUpwindModifiedSaturationW.data_given() && (1-m_imUpwindModifiedSaturationW[ip])>0.999 ) || xx==saturated_Node )//For each ip, find the saturated node
		{ 	
			//check the node is related to the subcontrol faces(from or to)
			//the upwind node must be from or to
			// 1. avoid inflow to the (downwind && saturated) node
			// 2. avoid outflow from the (upwind && saturated) node
			if (saturated_Node==-1)
			{
				saturated_Node = xx;
			}			
			
			if (xx==saturated_Node)
			{
				saturated = true;
			}

		}
		
		
		if (saturated && m_imMass.data_given()){
			size_t node = ((TElem*) elem)->vertex(ip)->grid_data_index();
			size_t ele_index = elem->grid_data_index();
			UG_LOGN("In ele# "<<ele_index<<"  "<<"at node# "<<node);
			UG_LOGN("  Mass of Pw-Eq. : "<<m_imMass[ip]<<"   UpwindModified_Sn = "<<1-m_imUpwindModifiedSaturationW[ip]<<"   UpwindSn = "<<1-m_imUpwindSaturationW[ip]<<'\n');
			//conv_flux =0;
		}
		
		*/
		

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
	if ( m_imPermeability.data_given() && m_imMinPd.data_given() && !m_imModifiedSaturationW.data_given() ){
		if (m_imEntryPressure.data_given())
			for (size_t i = 0; i < convShape.num_sh(); ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imEntryPressure[i]);
		else
			for (size_t i = 0; i < convShape.num_sh(); ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imPermeability[i], m_imPorosity[i]);
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
			//VecScaleAppend(linDefect, krw( m_imSaturationW[sh] ) * u_modified[sh], convShape.D_vel(ip, sh));
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
		if (m_imEntryPressure.data_given())
			for (size_t i = 0; i < convShape.num_sh(); ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imEntryPressure[i]);
		else
			for (size_t i = 0; i < convShape.num_sh(); ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imPermeability[i], m_imPorosity[i]);
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
		
		//bool saturated = false;
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			
		/*//For Extended Buckley Leverett
		if (m_imPermeability.data_given())
		{
			if (m_imPermeability[0]>5*pow(10,-14))
				VecScaleAppend(linDefect, krw1( u_modified[sh] ), convShape.D_vel(ip, sh));
			else
				VecScaleAppend(linDefect, krw2( u_modified[sh] ), convShape.D_vel(ip, sh));
		}
		*/
		
			/*		
			//Shuai debug
			if ( u_modified[sh]>0.999 )//For each ip, find the saturated node
			{ 	
				//check the node is related to the subcontrol faces(from or to)
				//the upwind node must be from or to
				// 1. avoid inflow to the (downwind && saturated) node
				// 2. avoid outflow from the (upwind && saturated) node
						
				//the saturated node is inside the subcontrolvolume //flow in, the upwind node is scvf.to()
				if ( ( sh==scvf.from() ) && ( convShape(ip, scvf.to())<0 ) ) 
				{
					saturated = true;
					//continue;
				}
				//the saturated node is outside the subcontrolvolume //flow out, the upwind node is scvf.from()
				else if ( ( sh==scvf.to() ) && ( convShape(ip, scvf.from())>0 ) )
				{
					saturated = true;
					//continue;
				}
			}
			//Shuai debug end
			*/
		
			if (iFormulationIndex == 1)
				VecScaleAppend(linDefect, krw( u(_C_, sh), m_imMinSwr[sh], m_imMinLambda[sh] ), convShape.D_vel(ip, sh));
			else if (iFormulationIndex == 2)
				VecScaleAppend(linDefect, krn( 1-u_modified[sh], m_imResidualCarbonic[0], m_imBrooksCoreyNumber[0] ), convShape.D_vel(ip, sh));
						
			//VecScaleAppend(linDefect, krw( u_modified[sh] ), convShape.D_vel(ip, sh));
			
			
		}
		/*
		if (saturated){
			VecSet(linDefect, 0.0);
		}
		*/
		
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
			//Shuai debug
			/*
			if (iFormulationIndex == 1)
				VecScaleAppend(linDefect, krn( m_imUpwindModifiedSaturationW[sh], m_imResidualCarbonic[0], m_imBrooksCoreyNumber[0] ), convShape.D_vel(ip, sh));
			else if (iFormulationIndex == 2)
				VecScaleAppend(linDefect, krw( m_imUpwindModifiedSaturationW[sh], m_imResidualAqueous[0], m_imBrooksCoreyNumber[0] ), convShape.D_vel(ip, sh));
			*/
			if (iFormulationIndex == 1)
				VecScaleAppend(linDefect, krn( m_imModifiedSaturationW[sh], m_imResidualCarbonic[0], m_imBrooksCoreyNumber[0] ), convShape.D_vel(ip, sh));
			else if (iFormulationIndex == 2)
				VecScaleAppend(linDefect, krw( m_imSaturationW[sh], m_imMinSwr[sh], m_imMinLambda[sh] ), convShape.D_vel(ip, sh));
				//VecScaleAppend(linDefect, krw( m_imUpwindSaturationW[sh], m_imResidualAqueous[0], m_imBrooksCoreyNumber[0] ), convShape.D_vel(ip, sh));
			
			/*
			if (sh>0)
			{
				if ( abs(m_imMinSwr[sh]-m_imMinSwr[sh-1])>0.00000001 )
				{
					int xx;
					xx = 0;
					
				}
			}
			*/
			
			//Shuai debug end
			
			
			//VecScaleAppend(linDefect, krn(  ), convShape.D_vel(ip, sh));
		
		/*//For Extended Buckley Leverett
		if (m_imPermeability.data_given())
		{
			if (m_imPermeability[0]>5*pow(10,-14))
				VecScaleAppend(linDefect, krn1( m_imSaturationW[sh] ), convShape.D_vel(ip, sh));
			else
				VecScaleAppend(linDefect, krn2( m_imSaturationW[sh] ), convShape.D_vel(ip, sh));
		}
		*/
		
			/*
			//Shuai debug
			if ( (1-m_imSaturationW[sh])>0.999 )//For each ip, find the saturated node
			{ 	
				//check the node is related to the subcontrol faces(from or to)
				//the upwind node must be from or to
				// 1. avoid inflow to the (downwind && saturated) node
				// 2. avoid outflow from the (upwind && saturated) node
						
				//the saturated node is inside the subcontrolvolume //flow in, the upwind node is scvf.to()
				if ( ( sh==scvf.from() ) && ( convShape(ip, scvf.to())<pow(10,-10) ) ) 
				{
					VecSet(linDefect, 0.0);
				}
				//the saturated node is outside the subcontrolvolume //flow out, the upwind node is scvf.from()
				else if ( ( sh==scvf.to() ) && ( convShape(ip, scvf.from())>-pow(10,-10) ) )
				{
					VecSet(linDefect, 0.0);
				}
			}
			//Shuai debug end
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
		if (m_imEntryPressure.data_given())
			for (size_t i = 0; i < convShape.num_sh(); ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imEntryPressure[i]);
		else
			for (size_t i = 0; i < convShape.num_sh(); ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imPermeability[i], m_imPorosity[i]);
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
			if (m_imEntryPressure.data_given())
				u_modified = Modify_s(iFormulationIndex, u(_C_, co), m_imResidualAqueous[co], m_imResidualCarbonic[co], m_imBrooksCoreyNumber[co], m_imMinPd[co], m_imMinSwr[co], m_imMinSnr[co], m_imMinLambda[co], m_imEntryPressure[co]);
			else
				u_modified = Modify_s(iFormulationIndex, u(_C_, co), m_imResidualAqueous[co], m_imResidualCarbonic[co], m_imBrooksCoreyNumber[co], m_imMinPd[co], m_imMinSwr[co], m_imMinSnr[co], m_imMinLambda[co], m_imPermeability[co], m_imPorosity[co]);
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
		if (m_imEntryPressure.data_given())
			for (size_t i = 0; i < geo.num_scv(); ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imEntryPressure[i]);
		else
			for (size_t i = 0; i < geo.num_scv(); ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imPermeability[i], m_imPorosity[i]);
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
ex_upwind_modified_value(number vValue[],
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
	
	//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo, false);
	
	// calculate modified sw
	//	Consider capillary trapping	or not
	number u_modified[numSH];
	if ( m_imPermeability.data_given() && m_imMinPd.data_given() && (nip != 1) ){
		if (m_imEntryPressure.data_given())
			for (size_t i = 0; i < numSH; ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imEntryPressure[i]);
		else
			for (size_t i = 0; i < numSH; ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imPermeability[i], m_imPorosity[i]);
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
		//	compute upwind value at ip
			vValue[ip] = 0.0;
			number sum_convShape = 0.0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
				vValue[ip] += u_modified[sh] * convShape(ip, sh);
				//sum_convShape += fabs(convShape(ip, sh));
				sum_convShape += convShape(ip, sh);
			}
			// Sw will be 0 when velocity is 0
			//if (sum_convShape > pow(10,-12))
			if (fabs(sum_convShape)>0)
				vValue[ip] = vValue[ip]/sum_convShape;
			else
				vValue[ip] = u_modified[0];
			/*
			else
			{
				vValue[ip] = 0.0;
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
					vValue[ip] += u_modified[sh] * scvf.shape(sh);
				}
			}	
			*/
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
			
		//	compute upwind value at ip
			vValue[ip] = 0.0;
			number sum_convShape = 0.0;
			for(size_t sh = 0; sh < numSH; ++sh)
			{
				vValue[ip] += u_modified[sh] * convShape(ip, sh);
				//sum_convShape += fabs(convShape(ip, sh));
				sum_convShape += convShape(ip, sh);
			}
			// Sw will be 0 when velocity is 0
			//if (sum_convShape > pow(10,-12))
			if (fabs(sum_convShape)>0)
				vValue[ip] = vValue[ip]/sum_convShape;
			else
				vValue[ip] = u_modified[0];
			/*
			else
			{
				vValue[ip] = 0.0;
				for(size_t sh = 0; sh < numSH; ++sh)
				{
					vValue[ip] += u_modified[sh] * vShape[sh];
				}
			}
			*/

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
ex_upwind_value(number vValue[],
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
	
	//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo, false);
	
//	MP SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//	Loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
		
		//consider about it !!!!!!! it means the order of vValue is the order of geo.scvf(), the order of u is the order of shape function
		//	compute upwind value at ip
			vValue[ip] = 0.0;
			number sum_convShape = 0.0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
				vValue[ip] += u(_C_, sh) * convShape(ip, sh);
				//sum_convShape += fabs(convShape(ip, sh));
				sum_convShape += convShape(ip, sh);
			}
			// Sw will be 0 when velocity is 0
			//if (sum_convShape > pow(10,-12))
			if (fabs(sum_convShape)>0)
				vValue[ip] = vValue[ip]/sum_convShape;
			else
				vValue[ip] = u(_C_, 0);
			/*
			else
			{
				vValue[ip] = 0.0;
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
					vValue[ip] += u_modified[sh] * scvf.shape(sh);
				}
			}	
			*/
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
			
		//	compute upwind value at ip
			vValue[ip] = 0.0;
			number sum_convShape = 0.0;
			for(size_t sh = 0; sh < numSH; ++sh)
			{
				vValue[ip] += u(_C_, sh) * convShape(ip, sh);
				//sum_convShape += fabs(convShape(ip, sh));
				sum_convShape += convShape(ip, sh);
			}
			// Sw will be 0 when velocity is 0
			//if (sum_convShape > pow(10,-12))
			if (fabs(sum_convShape)>0)
				vValue[ip] = vValue[ip]/sum_convShape;
			else
				vValue[ip] = u(_C_, 0);
			/*
			else
			{
				vValue[ip] = 0.0;
				for(size_t sh = 0; sh < numSH; ++sh)
				{
					vValue[ip] += u_modified[sh] * vShape[sh];
				}
			}
			*/

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
		if (m_imEntryPressure.data_given())
			for (size_t i = 0; i < numSH; ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imEntryPressure[i]);
		else
			for (size_t i = 0; i < numSH; ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imPermeability[i], m_imPorosity[i]);
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
		if (m_imEntryPressure.data_given())
			for (size_t i = 0; i < numSH; ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imEntryPressure[i]);
		else
			for (size_t i = 0; i < numSH; ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imPermeability[i], m_imPorosity[i]);
	}
	else{
		for (size_t i = 0; i < numSH; ++i)
		{
			u_modified[i] = u(_C_, i);
				
		}
			
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

	//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo, false);
	
	
	number u_modified[numSH];
	if ( m_imPermeability.data_given() && m_imMinPd.data_given() && (nip != 1) ){
		if (m_imEntryPressure.data_given())
			for (size_t i = 0; i < numSH; ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imEntryPressure[i]);
		else
			for (size_t i = 0; i < numSH; ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imPermeability[i], m_imPorosity[i]);
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
			
			//	to compute D \nabla Sw
				MathVector<dim> grad_c;
				VecSet(grad_c, 0.0);
				
			//gradient Sw
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecScaleAppend(grad_c, u_modified[sh], scvf.global_grad(sh));
			

			number s_integ = 0.0;
			for (size_t i = 0; i < scvf.num_sh(); ++i)
			{
				s_integ += u_modified[i]*scvf.shape(i);
			}
				
			MathMatrix<dim,dim> Real_Diffusion_Sw;
			if (m_imPermeability.data_given())
				if (m_imEntryPressure.data_given())
					MatScale(Real_Diffusion_Sw, Diffusion_S(iFormulationIndex, s_integ, m_imResidualAqueous[0], m_imResidualCarbonic[0], m_imBrooksCoreyNumber[0], m_imPermeability[0], m_imPorosity[0], m_imEntryPressure[0])
, m_imDiffusion_Sw[ip]);
				else
					MatScale(Real_Diffusion_Sw, Diffusion_S(iFormulationIndex, s_integ, m_imResidualAqueous[0], m_imResidualCarbonic[0], m_imBrooksCoreyNumber[0], m_imPermeability[0], m_imPorosity[0])
, m_imDiffusion_Sw[ip]);
			else
				MatScale(Real_Diffusion_Sw, 1, m_imDiffusion_Sw[ip]);

			
			//diffusion_Sw part
			MatVecMult(vValue[ip], Real_Diffusion_Sw, grad_c);
					
			
			//	compute upwind relative permeability value at ip 
			number krw_up = 0.0;
			number sum_convShape = 0.0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
				krw_up += krw( u(_C_, sh), m_imMinSwr[sh], m_imMinLambda[sh] ) * convShape(ip, sh);
				//sum_convShape += fabs(convShape(ip, sh));
				sum_convShape += convShape(ip, sh);
			}
			// Sw will be 0 when velocity is 0
			//if (sum_convShape > pow(10,-12))
			if (fabs(sum_convShape) > 0)
				krw_up = krw_up/sum_convShape;
			else
				krw_up = krw( u(_C_, 0), m_imMinSwr[0], m_imMinLambda[0] );
				
			VecScaleAppend(vValue[ip], -krw_up, m_imDarcyW[ip]);
			//VecScaleAppend(vValue[ip], -krw(Sw), m_imDarcyW[ip]);
			

			
			if(bDeriv)
			{
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
					vvvDeriv[ip][_C_][sh] = scvf.global_grad(sh);
					MatVecMult(vvvDeriv[ip][_C_][sh], Real_Diffusion_Sw, vvvDeriv[ip][_C_][sh]);
				}
					

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
ex_grad_pd2(MathVector<dim> vValue[],
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

	//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo, false);
	
	
	number u_modified[numSH];
	if ( m_imPermeability.data_given() && m_imMinPd.data_given() && (nip != 1) ){
		if (m_imEntryPressure.data_given())
			for (size_t i = 0; i < numSH; ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imEntryPressure[i]);
		else
			for (size_t i = 0; i < numSH; ++i)
				u_modified[i] = Modify_s(iFormulationIndex, u(_C_, i), m_imResidualAqueous[i], m_imResidualCarbonic[i], m_imBrooksCoreyNumber[i], m_imMinPd[i], m_imMinSwr[i], m_imMinSnr[i], m_imMinLambda[i], m_imPermeability[i], m_imPorosity[i]);
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
			//+ convection_Sw part
			//Shuai debug
			//	compute upwind relative permeability value at ip 
			number krw_up = 0.0;
			number sum_convShape = 0.0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
				krw_up += krw( m_imSaturationW[sh], m_imMinSwr[sh], m_imMinLambda[sh] ) * convShape(ip, sh);
				//sum_convShape += fabs(convShape(ip, sh));
				sum_convShape += convShape(ip, sh);
			}
			// Sw will be 0 when velocity is 0
			//if (sum_convShape > pow(10,-12))
			if (fabs(sum_convShape) > 0)
				krw_up = krw_up/sum_convShape;
			else
				krw_up = krw( m_imSaturationW[0], m_imMinSwr[0], m_imMinLambda[0] );
				
			//VecScaleAppend(vValue[ip], -krw(m_imUpwindModifiedSaturationW[ip], m_imResidualAqueous[0], m_imBrooksCoreyNumber[0]), m_imDarcyN[ip]);
			VecScaleAppend(vValue[ip], -krw_up, m_imDarcyN[ip]);
			
			//Shuai debug end
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
		///*
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
		//*/
		//	VecScaleAppend(vValue[ip], -krw(m_imSaturationW[ip], m_imResidualAqueous[0]), m_imDarcyN[ip]);
			
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
		if(m_imDiffusion.data_given()) vDiffusion = m_imDiffusion.values();
		//if(m_imDiffusion_Sw.data_given()) vDiffusion = m_imDiffusion_Sw.values();

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

/*
// error estimator parts
	this->set_prep_err_est_elem_loop(id, &T::template prep_err_est_elem_loop<TElem, TFVGeom>);
	this->set_prep_err_est_elem(id, &T::template prep_err_est_elem<TElem, TFVGeom>);
	this->set_compute_err_est_A_elem(id, &T::template compute_err_est_A_elem<TElem, TFVGeom>);
	this->set_compute_err_est_M_elem(id, &T::template compute_err_est_M_elem<TElem, TFVGeom>);
	this->set_compute_err_est_rhs_elem(id, &T::template compute_err_est_rhs_elem<TElem, TFVGeom>);
	this->set_fsh_err_est_elem_loop(id, &T::template fsh_err_est_elem_loop<TElem, TFVGeom>);
*/

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
	m_exUpwindModifiedValue->		template set_fct<T,refDim>(id, this, &T::template ex_upwind_modified_value<TElem, TFVGeom>);
	m_exUpwindValue->		template set_fct<T,refDim>(id, this, &T::template ex_upwind_value<TElem, TFVGeom>);
	m_exModifiedValue->		template set_fct<T,refDim>(id, this, &T::template ex_modified_value<TElem, TFVGeom>);
	m_exValue->		template set_fct<T,refDim>(id, this, &T::template ex_value<TElem, TFVGeom>);
	m_exGrad->		template set_fct<T,refDim>(id, this, &T::template ex_grad<TElem, TFVGeom>);
	m_exGrad_pd->		template set_fct<T,refDim>(id, this, &T::template ex_grad_pd<TElem, TFVGeom>);
	m_exGrad_pd2->		template set_fct<T,refDim>(id, this, &T::template ex_grad_pd2<TElem, TFVGeom>);

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

