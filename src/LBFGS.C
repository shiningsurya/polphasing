#include "LBFGS.hpp"
/*
 * There is no absolute phase in an interferometer. 
 * Arbitrarily set the phase of complex gain of some hand of some antenna to zero.
 * In this code, we set the complex gain phase of first hand of the first antenna to zero.
 * This means the imaginary part of RR::C00 is zero and the real part of RR::C00 is positive.
 * If the sign is not enforced, we would have 0<->pi phase jumps.
 *
 * Our LBFGS implementation only solves unconstrained optimization problem. 
 * Since we are adding only one constraint, it would be excessive to implement and use the bounded version of LBFGS. 
 * Instead, we use parameter transformation. 
 * Gain C00:RR = (exp(parameter), 0.0 )
 *
 * Encapsulate the reference part in GREF
 *
 * In hindsight, this could be done after solving by dividing solved complex gains by that of any reference antenna.3C138_unphased_rr.gains.
 *
 * why do we not need to reference when solving for full jones?
 *
 * COMBINED SOLVING
 * - We fit for unpolarized I. 
 *   Which should be positive. 
 *   So we use the same trick we applied for GREF
 *   I = exp(parameter)
 *   derivative_parameter = derivative_I * I
*/

//#define DPRINT

#ifdef DPRINT
#include <iostream>
#endif

using std::conj;

LBFGS::real_type LBFGS::combined_jones (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	c_data_t *pkg = static_cast<c_data_t*>(instance);

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    ngains ( 4 * pkg->nantennas );

	/* populate complex gains vector */
	vc_type  cgains ( ngains, complex_type(0.0f, 0.0f) );
	for ( int igain = 0; igain < ngains; igain++ ) {
		cgains[igain]   = complex_type ( rgains[2*igain+0], rgains[2*igain+1] );
	}

	/* create complex grad vector */
	vc_type       grad ( ngains, complex_type(0.0f, 0.0f) );
	/* for the unpol StokesI fitting, parameterization to ensure positiveness */
	const real_type  Iup ( std::exp(rgains[2*ngains]) );
	/* gradient wrt I is purely real, but we keep complex because intermediates are complex */
	complex_type  igrad (0.0f, 0.0f);

	/* iterate over the polar baselines */
	for ( int ibl = 0; ibl < pkg->npolarbaselines; ibl++ ) {

		// fetch the antenna index
		const int iant1 ( pkg->iant1[ibl] );
		const int iant2 ( pkg->iant2[ibl] );

		// fetch the par angles of unpolarized scan
		const real_type __t1 ( pkg->uol_par[iant1] );
		const real_type __t2 ( pkg->uol_par[iant2] );

		// complex exponentials of parang for unpolarized source
		const complex_type z1 ( cos (__t1), sin(__t1) );
		const complex_type z2 ( cos (__t2), sin(__t2) );

		// fetch the pb2corr
		const int pb2corr    = pkg->pb2corr [ ibl ];

		// fetch complex data
		const complex_type pol_data ( pkg->pol_data[ibl] );
		const complex_type pol_cata ( conj(pol_data) );

		const complex_type uol_data ( pkg->uol_data[ibl] );
		const complex_type uol_cata ( conj(uol_data) );

		// fetch the par corrected model
		const complex_type mrr ( pkg->polpar_model_rr[ibl] );
		const complex_type mrl ( pkg->polpar_model_rl[ibl] );
		const complex_type mlr ( pkg->polpar_model_lr[ibl] );
		const complex_type mll ( pkg->polpar_model_ll[ibl] );

		// set the gain indices
		const int iprr ( 4*iant1 + 0 );
		const int iprl ( 4*iant1 + 1 );
		const int iplr ( 4*iant1 + 2 );
		const int ipll ( 4*iant1 + 3 );

		const int iqrr ( 4*iant2 + 0 );
		const int iqrl ( 4*iant2 + 1 );
		const int iqlr ( 4*iant2 + 2 );
		const int iqll ( 4*iant2 + 3 );

		// fetch full gains for both antennas
		const complex_type gprr ( cgains[iprr] );
		const complex_type gprl ( cgains[iprl] );
		const complex_type gplr ( cgains[iplr] );
		const complex_type gpll ( cgains[ipll] );

		const complex_type gqrr ( cgains[iqrr] );
		const complex_type gqrl ( cgains[iqrl] );
		const complex_type gqlr ( cgains[iqlr] );
		const complex_type gqll ( cgains[iqll] );

		// the following long expressions come from sympy
		// see :math_gradient.py:
		// see :math_gradient.code:

		// model forward depends on pb2corr
		complex_type pol_model;
		complex_type uol_model;

		/*
		We do this over polarbaseline loop, so that we keep track of all the baselines
		This is for polarized data
		*/
#if 1
		if ( pb2corr == 0 ) {

const complex_type Dgprr_coeff_dpqrr = -gqrl*conj(mrl) - gqrr*conj(mrr) ; 
const complex_type Dgprl_coeff_dpqrr = -gqrl*conj(mll) - gqrr*conj(mlr) ; 
const complex_type Dgqrr_coeff_dqprr = -gprl*mlr - gprr*mrr ; 
const complex_type Dgqrl_coeff_dqprr = -gprl*mll - gprr*mrl ; 
const complex_type Dgprr_coeff_gprr = gqll*mrl*conj(gqll)*conj(mrl) + gqll*mrr*conj(gqlr)*conj(mrl) + gqlr*mrl*conj(gqll)*conj(mrr) + gqlr*mrr*conj(gqlr)*conj(mrr) + gqrl*mrl*conj(gqrl)*conj(mrl) + gqrl*mrr*conj(gqrr)*conj(mrl) + gqrr*mrl*conj(gqrl)*conj(mrr) + gqrr*mrr*conj(gqrr)*conj(mrr) ; 
const complex_type Dgprl_coeff_gprr = gqll*mrl*conj(gqll)*conj(mll) + gqll*mrr*conj(gqlr)*conj(mll) + gqlr*mrl*conj(gqll)*conj(mlr) + gqlr*mrr*conj(gqlr)*conj(mlr) + gqrl*mrl*conj(gqrl)*conj(mll) + gqrl*mrr*conj(gqrr)*conj(mll) + gqrr*mrl*conj(gqrl)*conj(mlr) + gqrr*mrr*conj(gqrr)*conj(mlr) ; 
const complex_type Dgqrr_coeff_gqrr = gpll*mlr*conj(gpll)*conj(mlr) + gpll*mlr*conj(gplr)*conj(mrr) + gplr*mrr*conj(gpll)*conj(mlr) + gplr*mrr*conj(gplr)*conj(mrr) + gprl*mlr*conj(gprl)*conj(mlr) + gprl*mlr*conj(gprr)*conj(mrr) + gprr*mrr*conj(gprl)*conj(mlr) + gprr*mrr*conj(gprr)*conj(mrr) ; 
const complex_type Dgqrl_coeff_gqrr = gpll*mll*conj(gpll)*conj(mlr) + gpll*mll*conj(gplr)*conj(mrr) + gplr*mrl*conj(gpll)*conj(mlr) + gplr*mrl*conj(gplr)*conj(mrr) + gprl*mll*conj(gprl)*conj(mlr) + gprl*mll*conj(gprr)*conj(mrr) + gprr*mrl*conj(gprl)*conj(mlr) + gprr*mrl*conj(gprr)*conj(mrr) ; 

			grad[iprr] += Dgprr_coeff_dpqrr*pol_data + Dgprr_coeff_gprr*gprr;
			grad[iprl] += Dgprl_coeff_dpqrr*pol_data + Dgprl_coeff_gprr*gprr;

			grad[iqrr] += Dgqrr_coeff_dqprr*pol_cata + Dgqrr_coeff_gqrr*gqrr;
			grad[iqrl] += Dgqrl_coeff_dqprr*pol_cata + Dgqrl_coeff_gqrr*gqrr;

pol_model = gprl*mll*conj(gqrl) + gprl*mlr*conj(gqrr) + gprr*mrl*conj(gqrl) + gprr*mrr*conj(gqrr) ;

		} // rr
		else if ( pb2corr == 1 ) {

const complex_type Dgprr_coeff_dpqrl = -gqll*conj(mrl) - gqlr*conj(mrr) ; 
const complex_type Dgprl_coeff_dpqrl = -gqll*conj(mll) - gqlr*conj(mlr) ; 
const complex_type Dgqlr_coeff_dqprl = -gprl*mlr - gprr*mrr ; 
const complex_type Dgqll_coeff_dqprl = -gprl*mll - gprr*mrl ; 
const complex_type Dgprr_coeff_gprl = gqll*mll*conj(gqll)*conj(mrl) + gqll*mlr*conj(gqlr)*conj(mrl) + gqlr*mll*conj(gqll)*conj(mrr) + gqlr*mlr*conj(gqlr)*conj(mrr) + gqrl*mll*conj(gqrl)*conj(mrl) + gqrl*mlr*conj(gqrr)*conj(mrl) + gqrr*mll*conj(gqrl)*conj(mrr) + gqrr*mlr*conj(gqrr)*conj(mrr) ; 
const complex_type Dgprl_coeff_gprl = gqll*mll*conj(gqll)*conj(mll) + gqll*mlr*conj(gqlr)*conj(mll) + gqlr*mll*conj(gqll)*conj(mlr) + gqlr*mlr*conj(gqlr)*conj(mlr) + gqrl*mll*conj(gqrl)*conj(mll) + gqrl*mlr*conj(gqrr)*conj(mll) + gqrr*mll*conj(gqrl)*conj(mlr) + gqrr*mlr*conj(gqrr)*conj(mlr) ; 
const complex_type Dgqrr_coeff_gqrl = gpll*mlr*conj(gpll)*conj(mll) + gpll*mlr*conj(gplr)*conj(mrl) + gplr*mrr*conj(gpll)*conj(mll) + gplr*mrr*conj(gplr)*conj(mrl) + gprl*mlr*conj(gprl)*conj(mll) + gprl*mlr*conj(gprr)*conj(mrl) + gprr*mrr*conj(gprl)*conj(mll) + gprr*mrr*conj(gprr)*conj(mrl) ; 
const complex_type Dgqrl_coeff_gqrl = gpll*mll*conj(gpll)*conj(mll) + gpll*mll*conj(gplr)*conj(mrl) + gplr*mrl*conj(gpll)*conj(mll) + gplr*mrl*conj(gplr)*conj(mrl) + gprl*mll*conj(gprl)*conj(mll) + gprl*mll*conj(gprr)*conj(mrl) + gprr*mrl*conj(gprl)*conj(mll) + gprr*mrl*conj(gprr)*conj(mrl) ; 

			grad[iprr] += Dgprr_coeff_dpqrl*pol_data + Dgprr_coeff_gprl*gprl;
			grad[iprl] += Dgprl_coeff_dpqrl*pol_data + Dgprl_coeff_gprl*gprl;

			grad[iqlr] += Dgqlr_coeff_dqprl*pol_cata;
			grad[iqll] += Dgqll_coeff_dqprl*pol_cata;

			grad[iqrr] += Dgqrr_coeff_gqrl*gqrl;
			grad[iqrl] += Dgqrl_coeff_gqrl*gqrl;

pol_model = gprl*mll*conj(gqll) + gprl*mlr*conj(gqlr) + gprr*mrl*conj(gqll) + gprr*mrr*conj(gqlr) ;

		} // rl
		else if ( pb2corr == 2 ) {

const complex_type Dgplr_coeff_dpqlr = -gqrl*conj(mrl) - gqrr*conj(mrr) ; 
const complex_type Dgpll_coeff_dpqlr = -gqrl*conj(mll) - gqrr*conj(mlr) ; 
const complex_type Dgqrr_coeff_dqplr = -gpll*mlr - gplr*mrr ; 
const complex_type Dgqrl_coeff_dqplr = -gpll*mll - gplr*mrl ; 
const complex_type Dgplr_coeff_gplr = gqll*mrl*conj(gqll)*conj(mrl) + gqll*mrr*conj(gqlr)*conj(mrl) + gqlr*mrl*conj(gqll)*conj(mrr) + gqlr*mrr*conj(gqlr)*conj(mrr) + gqrl*mrl*conj(gqrl)*conj(mrl) + gqrl*mrr*conj(gqrr)*conj(mrl) + gqrr*mrl*conj(gqrl)*conj(mrr) + gqrr*mrr*conj(gqrr)*conj(mrr) ; 
const complex_type Dgpll_coeff_gplr = gqll*mrl*conj(gqll)*conj(mll) + gqll*mrr*conj(gqlr)*conj(mll) + gqlr*mrl*conj(gqll)*conj(mlr) + gqlr*mrr*conj(gqlr)*conj(mlr) + gqrl*mrl*conj(gqrl)*conj(mll) + gqrl*mrr*conj(gqrr)*conj(mll) + gqrr*mrl*conj(gqrl)*conj(mlr) + gqrr*mrr*conj(gqrr)*conj(mlr) ; 
const complex_type Dgqlr_coeff_gqlr = gpll*mlr*conj(gpll)*conj(mlr) + gpll*mlr*conj(gplr)*conj(mrr) + gplr*mrr*conj(gpll)*conj(mlr) + gplr*mrr*conj(gplr)*conj(mrr) + gprl*mlr*conj(gprl)*conj(mlr) + gprl*mlr*conj(gprr)*conj(mrr) + gprr*mrr*conj(gprl)*conj(mlr) + gprr*mrr*conj(gprr)*conj(mrr) ; 
const complex_type Dgqll_coeff_gqlr = gpll*mll*conj(gpll)*conj(mlr) + gpll*mll*conj(gplr)*conj(mrr) + gplr*mrl*conj(gpll)*conj(mlr) + gplr*mrl*conj(gplr)*conj(mrr) + gprl*mll*conj(gprl)*conj(mlr) + gprl*mll*conj(gprr)*conj(mrr) + gprr*mrl*conj(gprl)*conj(mlr) + gprr*mrl*conj(gprr)*conj(mrr) ; 

			grad[iplr] += Dgplr_coeff_dpqlr*pol_data + Dgplr_coeff_gplr*gplr;
			grad[ipll] += Dgpll_coeff_dpqlr*pol_data + Dgpll_coeff_gplr*gplr;

			grad[iqrr] += Dgqrr_coeff_dqplr*pol_cata;
			grad[iqrl] += Dgqrl_coeff_dqplr*pol_cata;

			grad[iqlr] += Dgqlr_coeff_gqlr*gqlr;
			grad[iqll] += Dgqll_coeff_gqlr*gqlr;

pol_model = gpll*mll*conj(gqrl) + gpll*mlr*conj(gqrr) + gplr*mrl*conj(gqrl) + gplr*mrr*conj(gqrr) ;

		} // lr
		else if ( pb2corr == 3 ) {

const complex_type Dgplr_coeff_dpqll = -gqll*conj(mrl) - gqlr*conj(mrr) ; 
const complex_type Dgpll_coeff_dpqll = -gqll*conj(mll) - gqlr*conj(mlr) ; 
const complex_type Dgqlr_coeff_dqpll = -gpll*mlr - gplr*mrr ; 
const complex_type Dgqll_coeff_dqpll = -gpll*mll - gplr*mrl ; 
const complex_type Dgplr_coeff_gpll = gqll*mll*conj(gqll)*conj(mrl) + gqll*mlr*conj(gqlr)*conj(mrl) + gqlr*mll*conj(gqll)*conj(mrr) + gqlr*mlr*conj(gqlr)*conj(mrr) + gqrl*mll*conj(gqrl)*conj(mrl) + gqrl*mlr*conj(gqrr)*conj(mrl) + gqrr*mll*conj(gqrl)*conj(mrr) + gqrr*mlr*conj(gqrr)*conj(mrr) ; 
const complex_type Dgpll_coeff_gpll = gqll*mll*conj(gqll)*conj(mll) + gqll*mlr*conj(gqlr)*conj(mll) + gqlr*mll*conj(gqll)*conj(mlr) + gqlr*mlr*conj(gqlr)*conj(mlr) + gqrl*mll*conj(gqrl)*conj(mll) + gqrl*mlr*conj(gqrr)*conj(mll) + gqrr*mll*conj(gqrl)*conj(mlr) + gqrr*mlr*conj(gqrr)*conj(mlr) ; 
const complex_type Dgqlr_coeff_gqll = gpll*mlr*conj(gpll)*conj(mll) + gpll*mlr*conj(gplr)*conj(mrl) + gplr*mrr*conj(gpll)*conj(mll) + gplr*mrr*conj(gplr)*conj(mrl) + gprl*mlr*conj(gprl)*conj(mll) + gprl*mlr*conj(gprr)*conj(mrl) + gprr*mrr*conj(gprl)*conj(mll) + gprr*mrr*conj(gprr)*conj(mrl) ; 
const complex_type Dgqll_coeff_gqll = gpll*mll*conj(gpll)*conj(mll) + gpll*mll*conj(gplr)*conj(mrl) + gplr*mrl*conj(gpll)*conj(mll) + gplr*mrl*conj(gplr)*conj(mrl) + gprl*mll*conj(gprl)*conj(mll) + gprl*mll*conj(gprr)*conj(mrl) + gprr*mrl*conj(gprl)*conj(mll) + gprr*mrl*conj(gprr)*conj(mrl) ; 

			grad[iplr] += Dgplr_coeff_dpqll*pol_data + Dgplr_coeff_gpll*gpll;
			grad[ipll] += Dgpll_coeff_dpqll*pol_data + Dgpll_coeff_gpll*gpll;

			grad[iqlr] += Dgqlr_coeff_dqpll*pol_cata + Dgqlr_coeff_gqll*gqll;
			grad[iqll] += Dgqll_coeff_dqpll*pol_cata + Dgqll_coeff_gqll*gqll;

pol_model = gpll*mll*conj(gqll) + gpll*mlr*conj(gqlr) + gplr*mrl*conj(gqll) + gplr*mrr*conj(gqlr) ;

		} // ll
#endif
		/*
		We do this over polarbaseline loop, so that we keep track of all the baselines
		This is for unpolarized data
		*/
		if ( pb2corr == 0 ) {

const complex_type Dgprr_coeff_dpqrr = -Iup*gqrr*z1*conj(z2) ; 
const complex_type Dgprl_coeff_dpqrr = -Iup*gqrl*z2*conj(z1) ; 

const complex_type Dgprr_coeff_gprr = (Iup * Iup)*gqlr*z1*z2*conj(gqlr)*conj(z1)*conj(z2) + (Iup * Iup)*gqrr*z1*z2*conj(gqrr)*conj(z1)*conj(z2) ; 
const complex_type Dgprl_coeff_gprr = (Iup * Iup)*gqll*(z2 * z2)*conj(gqlr)*(conj(z1) * conj(z1)) + (Iup * Iup)*gqrl*(z2 * z2)*conj(gqrr)*(conj(z1) * conj(z1)) ; 

const complex_type Dgqrr_coeff_dqprr = -Iup*gprr*z2*conj(z1) ; 
const complex_type Dgqrl_coeff_dqprr = -Iup*gprl*z1*conj(z2) ; 

const complex_type Dgqrr_coeff_gqrr = (Iup * Iup)*gplr*z1*z2*conj(gplr)*conj(z1)*conj(z2) + (Iup * Iup)*gprr*z1*z2*conj(gprr)*conj(z1)*conj(z2) ; 
const complex_type Dgqrl_coeff_gqrr = (Iup * Iup)*gpll*(z1 * z1)*conj(gplr)*(conj(z2) * conj(z2)) + (Iup * Iup)*gprl*(z1 * z1)*conj(gprr)*(conj(z2) * conj(z2)) ; 

const complex_type DI_coeff_dpqrr = -gqrl*z2*conj(gprl)*conj(z1) - gqrr*z1*conj(gprr)*conj(z2) ; 
const complex_type DI_coeff_dqprr = -gprl*z1*conj(gqrl)*conj(z2) - gprr*z2*conj(gqrr)*conj(z1) ; 

/* this is the constant that must be added once*/
const complex_type DI_coeff_constant = 2*Iup*gpll*gqll*z1*z2*conj(gpll)*conj(gqll)*conj(z1)*conj(z2) + 2*Iup*gpll*gqlr*(z1 * z1)*conj(gplr)*conj(gqll)*(conj(z2) * conj(z2)) + 2*Iup*gpll*gqrl*z1*z2*conj(gpll)*conj(gqrl)*conj(z1)*conj(z2) + 2*Iup*gpll*gqrr*(z1 * z1)*conj(gplr)*conj(gqrl)*(conj(z2) * conj(z2)) + 2*Iup*gplr*gqll*(z2 * z2)*conj(gpll)*conj(gqlr)*(conj(z1) * conj(z1)) + 2*Iup*gplr*gqlr*z1*z2*conj(gplr)*conj(gqlr)*conj(z1)*conj(z2) + 2*Iup*gplr*gqrl*(z2 * z2)*conj(gpll)*conj(gqrr)*(conj(z1) * conj(z1)) + 2*Iup*gplr*gqrr*z1*z2*conj(gplr)*conj(gqrr)*conj(z1)*conj(z2) + 2*Iup*gprl*gqll*z1*z2*conj(gprl)*conj(gqll)*conj(z1)*conj(z2) + 2*Iup*gprl*gqlr*(z1 * z1)*conj(gprr)*conj(gqll)*(conj(z2) * conj(z2)) + 2*Iup*gprl*gqrl*z1*z2*conj(gprl)*conj(gqrl)*conj(z1)*conj(z2) + 2*Iup*gprl*gqrr*(z1 * z1)*conj(gprr)*conj(gqrl)*(conj(z2) * conj(z2)) + 2*Iup*gprr*gqll*(z2 * z2)*conj(gprl)*conj(gqlr)*(conj(z1) * conj(z1)) + 2*Iup*gprr*gqlr*z1*z2*conj(gprr)*conj(gqlr)*conj(z1)*conj(z2) + 2*Iup*gprr*gqrl*(z2 * z2)*conj(gprl)*conj(gqrr)*(conj(z1) * conj(z1)) + 2*Iup*gprr*gqrr*z1*z2*conj(gprr)*conj(gqrr)*conj(z1)*conj(z2) ; 

			grad[iprr] += Dgprr_coeff_dpqrr*uol_data + Dgprr_coeff_gprr*gprr;
			grad[iprl] += Dgprl_coeff_dpqrr*uol_data + Dgprl_coeff_gprr*gprr;

			grad[iqrr] += Dgqrr_coeff_dqprr*uol_cata + Dgqrr_coeff_gqrr*gqrr;
			grad[iqrl] += Dgqrl_coeff_dqprr*uol_cata + Dgqrl_coeff_gqrr*gqrr;

			igrad += DI_coeff_dpqrr*uol_data + DI_coeff_dqprr*uol_cata;
			igrad += DI_coeff_constant;

			uol_model = Iup*gprl*z1*conj(gqrl)*conj(z2) + Iup*gprr*z2*conj(gqrr)*conj(z1);

		} else if ( pb2corr == 1 ) {

const complex_type Dgprr_coeff_dpqrl = -Iup*gqlr*z1*conj(z2) ; 
const complex_type Dgprl_coeff_dpqrl = -Iup*gqll*z2*conj(z1) ; 

const complex_type Dgprr_coeff_gprl = (Iup * Iup)*gqlr*(z1 * z1)*conj(gqll)*(conj(z2) * conj(z2)) + (Iup * Iup)*gqrr*(z1 * z1)*conj(gqrl)*(conj(z2) * conj(z2)) ; 
const complex_type Dgprl_coeff_gprl = (Iup * Iup)*gqll*z1*z2*conj(gqll)*conj(z1)*conj(z2) + (Iup * Iup)*gqrl*z1*z2*conj(gqrl)*conj(z1)*conj(z2) ; 

const complex_type Dgqrr_coeff_gqrl = (Iup * Iup)*gplr*(z2 * z2)*conj(gpll)*(conj(z1) * conj(z1)) + (Iup * Iup)*gprr*(z2 * z2)*conj(gprl)*(conj(z1) * conj(z1)) ; 
const complex_type Dgqrl_coeff_gqrl = (Iup * Iup)*gpll*z1*z2*conj(gpll)*conj(z1)*conj(z2) + (Iup * Iup)*gprl*z1*z2*conj(gprl)*conj(z1)*conj(z2) ; 

const complex_type Dgqlr_coeff_dqprl = -Iup*gprr*z2*conj(z1) ; 
const complex_type Dgqll_coeff_dqprl = -Iup*gprl*z1*conj(z2) ; 

const complex_type DI_coeff_dpqrl = -gqll*z2*conj(gprl)*conj(z1) - gqlr*z1*conj(gprr)*conj(z2) ; 
const complex_type DI_coeff_dqprl = -gprl*z1*conj(gqll)*conj(z2) - gprr*z2*conj(gqlr)*conj(z1) ; 
			
			igrad += DI_coeff_dpqrl*uol_data + DI_coeff_dqprl*uol_cata;

			uol_model = Iup*gprl*z1*conj(gqll)*conj(z2) + Iup*gprr*z2*conj(gqlr)*conj(z1);

			grad[iprr] += Dgprr_coeff_dpqrl*uol_data + Dgprr_coeff_gprl*gprl;
			grad[iprl] += Dgprl_coeff_dpqrl*uol_data + Dgprl_coeff_gprl*gprl;

			grad[iqrr] += Dgqrr_coeff_gqrl*gqrl;
			grad[iqrl] += Dgqrl_coeff_gqrl*gqrl;

			grad[iqlr] += Dgqlr_coeff_dqprl*uol_cata;
			grad[iqll] += Dgqll_coeff_dqprl*uol_cata;

		} else if ( pb2corr == 2 ) {

const complex_type Dgqrr_coeff_dqplr = -Iup*gplr*z2*conj(z1) ; 
const complex_type Dgqrl_coeff_dqplr = -Iup*gpll*z1*conj(z2) ; 

const complex_type Dgplr_coeff_dpqlr = -Iup*gqrr*z1*conj(z2) ; 
const complex_type Dgpll_coeff_dpqlr = -Iup*gqrl*z2*conj(z1) ; 

const complex_type Dgplr_coeff_gplr = (Iup * Iup)*gqlr*z1*z2*conj(gqlr)*conj(z1)*conj(z2) + (Iup * Iup)*gqrr*z1*z2*conj(gqrr)*conj(z1)*conj(z2) ; 
const complex_type Dgpll_coeff_gplr = (Iup * Iup)*gqll*(z2 * z2)*conj(gqlr)*(conj(z1) * conj(z1)) + (Iup * Iup)*gqrl*(z2 * z2)*conj(gqrr)*(conj(z1) * conj(z1)) ; 

const complex_type Dgqlr_coeff_gqlr = (Iup * Iup)*gplr*z1*z2*conj(gplr)*conj(z1)*conj(z2) + (Iup * Iup)*gprr*z1*z2*conj(gprr)*conj(z1)*conj(z2) ; 
const complex_type Dgqll_coeff_gqlr = (Iup * Iup)*gpll*(z1 * z1)*conj(gplr)*(conj(z2) * conj(z2)) + (Iup * Iup)*gprl*(z1 * z1)*conj(gprr)*(conj(z2) * conj(z2)) ; 

const complex_type DI_coeff_dpqlr = -gqrl*z2*conj(gpll)*conj(z1) - gqrr*z1*conj(gplr)*conj(z2) ; 
const complex_type DI_coeff_dqplr = -gpll*z1*conj(gqrl)*conj(z2) - gplr*z2*conj(gqrr)*conj(z1) ; 
			
			igrad += DI_coeff_dpqlr*uol_data + DI_coeff_dqplr*uol_cata;

			uol_model = Iup*gpll*z1*conj(gqrl)*conj(z2) + Iup*gplr*z2*conj(gqrr)*conj(z1);

			grad[iplr] += Dgplr_coeff_dpqlr*uol_data + Dgplr_coeff_gplr*gplr;
			grad[ipll] += Dgpll_coeff_dpqlr*uol_data + Dgpll_coeff_gplr*gplr;

			grad[iqrr] += Dgqrr_coeff_dqplr*uol_cata;
			grad[iqrl] += Dgqrl_coeff_dqplr*uol_cata;

			grad[iqlr] += Dgqlr_coeff_gqlr*gqlr;
			grad[iqll] += Dgqll_coeff_gqlr*gqlr;

		} else if ( pb2corr == 3 ) {

const complex_type Dgplr_coeff_dpqll = -Iup*gqlr*z1*conj(z2) ; 
const complex_type Dgpll_coeff_dpqll = -Iup*gqll*z2*conj(z1) ; 

const complex_type Dgplr_coeff_gpll = (Iup * Iup)*gqlr*(z1 * z1)*conj(gqll)*(conj(z2) * conj(z2)) + (Iup * Iup)*gqrr*(z1 * z1)*conj(gqrl)*(conj(z2) * conj(z2)) ; 
const complex_type Dgpll_coeff_gpll = (Iup * Iup)*gqll*z1*z2*conj(gqll)*conj(z1)*conj(z2) + (Iup * Iup)*gqrl*z1*z2*conj(gqrl)*conj(z1)*conj(z2) ; 

const complex_type Dgqlr_coeff_dqpll = -Iup*gplr*z2*conj(z1) ; 
const complex_type Dgqll_coeff_dqpll = -Iup*gpll*z1*conj(z2) ; 

const complex_type Dgqlr_coeff_gqll = (Iup * Iup)*gplr*(z2 * z2)*conj(gpll)*(conj(z1) * conj(z1)) + (Iup * Iup)*gprr*(z2 * z2)*conj(gprl)*(conj(z1) * conj(z1)) ; 
const complex_type Dgqll_coeff_gqll = (Iup * Iup)*gpll*z1*z2*conj(gpll)*conj(z1)*conj(z2) + (Iup * Iup)*gprl*z1*z2*conj(gprl)*conj(z1)*conj(z2) ; 

const complex_type DI_coeff_dpqll = -gqll*z2*conj(gpll)*conj(z1) - gqlr*z1*conj(gplr)*conj(z2) ; 
const complex_type DI_coeff_dqpll = -gpll*z1*conj(gqll)*conj(z2) - gplr*z2*conj(gqlr)*conj(z1) ; 

			igrad += DI_coeff_dpqll*uol_data + DI_coeff_dqpll*uol_cata;

			uol_model = Iup*gpll*z1*conj(gqll)*conj(z2) + Iup*gplr*z2*conj(gqlr)*conj(z1);

			grad[iplr] += Dgplr_coeff_dpqll*uol_data + Dgplr_coeff_gpll*gpll;
			grad[ipll] += Dgpll_coeff_dpqll*uol_data + Dgpll_coeff_gpll*gpll;

			grad[iqlr] += Dgqlr_coeff_dqpll*uol_cata + Dgqlr_coeff_gqll*gqll;
			grad[iqll] += Dgqll_coeff_dqpll*uol_cata + Dgqll_coeff_gqll*gqll;

		}
			
		// update cost
		cost += std::norm ( pol_data - pol_model );
		cost += std::norm ( uol_data - uol_model );
		//std::cout << " iterationcost=" << cost << " ";
		//std::cout <<  " cost uol=" << std::norm(uol_data - uol_model) << " pol=" << std::norm(pol_data - pol_model) << std::endl;

	} // iterate over polar baselines
	
	//std::cout << " rgrads=";

	/* load complex grad into real and imaginary parts */
	for ( int igain = 0; igain < ngains; igain++ ) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain + 0] = gg.real();
		rgrad[2*igain + 1] = gg.imag();
	}

	rgrad[2*ngains] = igrad.real() * Iup;
	/* igrad is purely real mathematically, but we use complex because intermediates are complex */

#ifdef DPRINT
	std::cout << " full_jones_cost=" << cost << " Iup=" << Iup << std::endl; 
#endif

	return cost;
}


LBFGS::real_type LBFGS::full_jones (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const data_t *pkg = reinterpret_cast<const data_t*>(instance);

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    ngains ( 4 * pkg->nantennas );

	/* populate complex gains vector */
	vc_type  cgains ( ngains, complex_type(0.0f, 0.0f) );
	for ( int igain = 0; igain < ngains; igain++ ) {
		cgains[igain]   = complex_type ( rgains[2*igain+0], rgains[2*igain+1] );
	}

	/* create complex grad vector */
	vc_type  grad ( ngains, complex_type(0.0f, 0.0f) );

	/* iterate over the polar baselines */
	for ( int ibl = 0; ibl < pkg->npolarbaselines; ibl++ ) {

		// fetch the antenna index
		const int iant1 ( pkg->iant1[ibl] );
		const int iant2 ( pkg->iant2[ibl] );

		// fetch the pb2corr
		const int pb2corr    = pkg->pb2corr [ ibl ];

		// fetch complex data
		const complex_type data ( pkg->data[ibl] );
		const complex_type cata ( conj(data) );

		// fetch the par corrected model
		const complex_type mrr ( pkg->par_model_rr[ibl] );
		const complex_type mrl ( pkg->par_model_rl[ibl] );
		const complex_type mlr ( pkg->par_model_lr[ibl] );
		const complex_type mll ( pkg->par_model_ll[ibl] );

		// set the gain indices
		const int iprr ( 4*iant1 + 0 );
		const int iprl ( 4*iant1 + 1 );
		const int iplr ( 4*iant1 + 2 );
		const int ipll ( 4*iant1 + 3 );

		const int iqrr ( 4*iant2 + 0 );
		const int iqrl ( 4*iant2 + 1 );
		const int iqlr ( 4*iant2 + 2 );
		const int iqll ( 4*iant2 + 3 );

		// fetch full gains for both antennas
		const complex_type gprr ( cgains[iprr] );
		const complex_type gprl ( cgains[iprl] );
		const complex_type gplr ( cgains[iplr] );
		const complex_type gpll ( cgains[ipll] );

		const complex_type gqrr ( cgains[iqrr] );
		const complex_type gqrl ( cgains[iqrl] );
		const complex_type gqlr ( cgains[iqlr] );
		const complex_type gqll ( cgains[iqll] );

		// the following long expressions come from sympy
		// see :math_gradient.py:
		// see :math_gradient.code:

		// model forward depends on pb2corr
		complex_type model;

		/*
		We do this over polarbaseline loop, so that we keep track of all the baselines
		*/
		if ( pb2corr == 0 ) {

const complex_type Dgprr_coeff_dpqrr = -gqrl*conj(mrl) - gqrr*conj(mrr) ; 
const complex_type Dgprl_coeff_dpqrr = -gqrl*conj(mll) - gqrr*conj(mlr) ; 
const complex_type Dgqrr_coeff_dqprr = -gprl*mlr - gprr*mrr ; 
const complex_type Dgqrl_coeff_dqprr = -gprl*mll - gprr*mrl ; 
const complex_type Dgprr_coeff_gprr = gqll*mrl*conj(gqll)*conj(mrl) + gqll*mrr*conj(gqlr)*conj(mrl) + gqlr*mrl*conj(gqll)*conj(mrr) + gqlr*mrr*conj(gqlr)*conj(mrr) + gqrl*mrl*conj(gqrl)*conj(mrl) + gqrl*mrr*conj(gqrr)*conj(mrl) + gqrr*mrl*conj(gqrl)*conj(mrr) + gqrr*mrr*conj(gqrr)*conj(mrr) ; 
const complex_type Dgprl_coeff_gprr = gqll*mrl*conj(gqll)*conj(mll) + gqll*mrr*conj(gqlr)*conj(mll) + gqlr*mrl*conj(gqll)*conj(mlr) + gqlr*mrr*conj(gqlr)*conj(mlr) + gqrl*mrl*conj(gqrl)*conj(mll) + gqrl*mrr*conj(gqrr)*conj(mll) + gqrr*mrl*conj(gqrl)*conj(mlr) + gqrr*mrr*conj(gqrr)*conj(mlr) ; 
const complex_type Dgqrr_coeff_gqrr = gpll*mlr*conj(gpll)*conj(mlr) + gpll*mlr*conj(gplr)*conj(mrr) + gplr*mrr*conj(gpll)*conj(mlr) + gplr*mrr*conj(gplr)*conj(mrr) + gprl*mlr*conj(gprl)*conj(mlr) + gprl*mlr*conj(gprr)*conj(mrr) + gprr*mrr*conj(gprl)*conj(mlr) + gprr*mrr*conj(gprr)*conj(mrr) ; 
const complex_type Dgqrl_coeff_gqrr = gpll*mll*conj(gpll)*conj(mlr) + gpll*mll*conj(gplr)*conj(mrr) + gplr*mrl*conj(gpll)*conj(mlr) + gplr*mrl*conj(gplr)*conj(mrr) + gprl*mll*conj(gprl)*conj(mlr) + gprl*mll*conj(gprr)*conj(mrr) + gprr*mrl*conj(gprl)*conj(mlr) + gprr*mrl*conj(gprr)*conj(mrr) ; 

			grad[iprr] += Dgprr_coeff_dpqrr*data + Dgprr_coeff_gprr*gprr;
			grad[iprl] += Dgprl_coeff_dpqrr*data + Dgprl_coeff_gprr*gprr;

			grad[iqrr] += Dgqrr_coeff_dqprr*cata + Dgqrr_coeff_gqrr*gqrr;
			grad[iqrl] += Dgqrl_coeff_dqprr*cata + Dgqrl_coeff_gqrr*gqrr;

model = gprl*mll*conj(gqrl) + gprl*mlr*conj(gqrr) + gprr*mrl*conj(gqrl) + gprr*mrr*conj(gqrr) ;

		} // rr
		else if ( pb2corr == 1 ) {

const complex_type Dgprr_coeff_dpqrl = -gqll*conj(mrl) - gqlr*conj(mrr) ; 
const complex_type Dgprl_coeff_dpqrl = -gqll*conj(mll) - gqlr*conj(mlr) ; 
const complex_type Dgqlr_coeff_dqprl = -gprl*mlr - gprr*mrr ; 
const complex_type Dgqll_coeff_dqprl = -gprl*mll - gprr*mrl ; 
const complex_type Dgprr_coeff_gprl = gqll*mll*conj(gqll)*conj(mrl) + gqll*mlr*conj(gqlr)*conj(mrl) + gqlr*mll*conj(gqll)*conj(mrr) + gqlr*mlr*conj(gqlr)*conj(mrr) + gqrl*mll*conj(gqrl)*conj(mrl) + gqrl*mlr*conj(gqrr)*conj(mrl) + gqrr*mll*conj(gqrl)*conj(mrr) + gqrr*mlr*conj(gqrr)*conj(mrr) ; 
const complex_type Dgprl_coeff_gprl = gqll*mll*conj(gqll)*conj(mll) + gqll*mlr*conj(gqlr)*conj(mll) + gqlr*mll*conj(gqll)*conj(mlr) + gqlr*mlr*conj(gqlr)*conj(mlr) + gqrl*mll*conj(gqrl)*conj(mll) + gqrl*mlr*conj(gqrr)*conj(mll) + gqrr*mll*conj(gqrl)*conj(mlr) + gqrr*mlr*conj(gqrr)*conj(mlr) ; 
const complex_type Dgqrr_coeff_gqrl = gpll*mlr*conj(gpll)*conj(mll) + gpll*mlr*conj(gplr)*conj(mrl) + gplr*mrr*conj(gpll)*conj(mll) + gplr*mrr*conj(gplr)*conj(mrl) + gprl*mlr*conj(gprl)*conj(mll) + gprl*mlr*conj(gprr)*conj(mrl) + gprr*mrr*conj(gprl)*conj(mll) + gprr*mrr*conj(gprr)*conj(mrl) ; 
const complex_type Dgqrl_coeff_gqrl = gpll*mll*conj(gpll)*conj(mll) + gpll*mll*conj(gplr)*conj(mrl) + gplr*mrl*conj(gpll)*conj(mll) + gplr*mrl*conj(gplr)*conj(mrl) + gprl*mll*conj(gprl)*conj(mll) + gprl*mll*conj(gprr)*conj(mrl) + gprr*mrl*conj(gprl)*conj(mll) + gprr*mrl*conj(gprr)*conj(mrl) ; 

			grad[iprr] += Dgprr_coeff_dpqrl*data + Dgprr_coeff_gprl*gprl;
			grad[iprl] += Dgprl_coeff_dpqrl*data + Dgprl_coeff_gprl*gprl;

			grad[iqlr] += Dgqlr_coeff_dqprl*cata;
			grad[iqll] += Dgqll_coeff_dqprl*cata;

			grad[iqrr] += Dgqrr_coeff_gqrl*gqrl;
			grad[iqrl] += Dgqrl_coeff_gqrl*gqrl;

model = gprl*mll*conj(gqll) + gprl*mlr*conj(gqlr) + gprr*mrl*conj(gqll) + gprr*mrr*conj(gqlr) ;

		} // rl
		else if ( pb2corr == 2 ) {

const complex_type Dgplr_coeff_dpqlr = -gqrl*conj(mrl) - gqrr*conj(mrr) ; 
const complex_type Dgpll_coeff_dpqlr = -gqrl*conj(mll) - gqrr*conj(mlr) ; 
const complex_type Dgqrr_coeff_dqplr = -gpll*mlr - gplr*mrr ; 
const complex_type Dgqrl_coeff_dqplr = -gpll*mll - gplr*mrl ; 
const complex_type Dgplr_coeff_gplr = gqll*mrl*conj(gqll)*conj(mrl) + gqll*mrr*conj(gqlr)*conj(mrl) + gqlr*mrl*conj(gqll)*conj(mrr) + gqlr*mrr*conj(gqlr)*conj(mrr) + gqrl*mrl*conj(gqrl)*conj(mrl) + gqrl*mrr*conj(gqrr)*conj(mrl) + gqrr*mrl*conj(gqrl)*conj(mrr) + gqrr*mrr*conj(gqrr)*conj(mrr) ; 
const complex_type Dgpll_coeff_gplr = gqll*mrl*conj(gqll)*conj(mll) + gqll*mrr*conj(gqlr)*conj(mll) + gqlr*mrl*conj(gqll)*conj(mlr) + gqlr*mrr*conj(gqlr)*conj(mlr) + gqrl*mrl*conj(gqrl)*conj(mll) + gqrl*mrr*conj(gqrr)*conj(mll) + gqrr*mrl*conj(gqrl)*conj(mlr) + gqrr*mrr*conj(gqrr)*conj(mlr) ; 
const complex_type Dgqlr_coeff_gqlr = gpll*mlr*conj(gpll)*conj(mlr) + gpll*mlr*conj(gplr)*conj(mrr) + gplr*mrr*conj(gpll)*conj(mlr) + gplr*mrr*conj(gplr)*conj(mrr) + gprl*mlr*conj(gprl)*conj(mlr) + gprl*mlr*conj(gprr)*conj(mrr) + gprr*mrr*conj(gprl)*conj(mlr) + gprr*mrr*conj(gprr)*conj(mrr) ; 
const complex_type Dgqll_coeff_gqlr = gpll*mll*conj(gpll)*conj(mlr) + gpll*mll*conj(gplr)*conj(mrr) + gplr*mrl*conj(gpll)*conj(mlr) + gplr*mrl*conj(gplr)*conj(mrr) + gprl*mll*conj(gprl)*conj(mlr) + gprl*mll*conj(gprr)*conj(mrr) + gprr*mrl*conj(gprl)*conj(mlr) + gprr*mrl*conj(gprr)*conj(mrr) ; 

			grad[iplr] += Dgplr_coeff_dpqlr*data + Dgplr_coeff_gplr*gplr;
			grad[ipll] += Dgpll_coeff_dpqlr*data + Dgpll_coeff_gplr*gplr;

			grad[iqrr] += Dgqrr_coeff_dqplr*cata;
			grad[iqrl] += Dgqrl_coeff_dqplr*cata;

			grad[iqlr] += Dgqlr_coeff_gqlr*gqlr;
			grad[iqll] += Dgqll_coeff_gqlr*gqlr;

model = gpll*mll*conj(gqrl) + gpll*mlr*conj(gqrr) + gplr*mrl*conj(gqrl) + gplr*mrr*conj(gqrr) ;

		} // lr
		else if ( pb2corr == 3 ) {

const complex_type Dgplr_coeff_dpqll = -gqll*conj(mrl) - gqlr*conj(mrr) ; 
const complex_type Dgpll_coeff_dpqll = -gqll*conj(mll) - gqlr*conj(mlr) ; 
const complex_type Dgqlr_coeff_dqpll = -gpll*mlr - gplr*mrr ; 
const complex_type Dgqll_coeff_dqpll = -gpll*mll - gplr*mrl ; 
const complex_type Dgplr_coeff_gpll = gqll*mll*conj(gqll)*conj(mrl) + gqll*mlr*conj(gqlr)*conj(mrl) + gqlr*mll*conj(gqll)*conj(mrr) + gqlr*mlr*conj(gqlr)*conj(mrr) + gqrl*mll*conj(gqrl)*conj(mrl) + gqrl*mlr*conj(gqrr)*conj(mrl) + gqrr*mll*conj(gqrl)*conj(mrr) + gqrr*mlr*conj(gqrr)*conj(mrr) ; 
const complex_type Dgpll_coeff_gpll = gqll*mll*conj(gqll)*conj(mll) + gqll*mlr*conj(gqlr)*conj(mll) + gqlr*mll*conj(gqll)*conj(mlr) + gqlr*mlr*conj(gqlr)*conj(mlr) + gqrl*mll*conj(gqrl)*conj(mll) + gqrl*mlr*conj(gqrr)*conj(mll) + gqrr*mll*conj(gqrl)*conj(mlr) + gqrr*mlr*conj(gqrr)*conj(mlr) ; 
const complex_type Dgqlr_coeff_gqll = gpll*mlr*conj(gpll)*conj(mll) + gpll*mlr*conj(gplr)*conj(mrl) + gplr*mrr*conj(gpll)*conj(mll) + gplr*mrr*conj(gplr)*conj(mrl) + gprl*mlr*conj(gprl)*conj(mll) + gprl*mlr*conj(gprr)*conj(mrl) + gprr*mrr*conj(gprl)*conj(mll) + gprr*mrr*conj(gprr)*conj(mrl) ; 
const complex_type Dgqll_coeff_gqll = gpll*mll*conj(gpll)*conj(mll) + gpll*mll*conj(gplr)*conj(mrl) + gplr*mrl*conj(gpll)*conj(mll) + gplr*mrl*conj(gplr)*conj(mrl) + gprl*mll*conj(gprl)*conj(mll) + gprl*mll*conj(gprr)*conj(mrl) + gprr*mrl*conj(gprl)*conj(mll) + gprr*mrl*conj(gprr)*conj(mrl) ; 

			grad[iplr] += Dgplr_coeff_dpqll*data + Dgplr_coeff_gpll*gpll;
			grad[ipll] += Dgpll_coeff_dpqll*data + Dgpll_coeff_gpll*gpll;

			grad[iqlr] += Dgqlr_coeff_dqpll*cata + Dgqlr_coeff_gqll*gqll;
			grad[iqll] += Dgqll_coeff_dqpll*cata + Dgqll_coeff_gqll*gqll;

model = gpll*mll*conj(gqll) + gpll*mlr*conj(gqlr) + gplr*mrl*conj(gqll) + gplr*mrr*conj(gqlr) ;

		} // ll
			
		// update cost
		cost += std::norm ( data - model );
		//std::cout << " iterationcost=" << cost << " ";

	} // iterate over polar baselines
	
	//std::cout << " rgrads=";

	/* load complex grad into real and imaginary parts */
	for ( int igain = 0; igain < ngains; igain++ ) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain + 0] = gg.real();
		rgrad[2*igain + 1] = gg.imag();
		/* we are actually computing derivative wrt conjugate, so lets unconjugate */
		//std::cout << rgrad[2*igain + 0] << " " << rgrad[2*igain + 1] << " ";
		//std::cout << gg << " ";
	}

	//std::cout << " full_jones_cost=" << cost << std::endl; 

	return cost;
}

LBFGS::real_type LBFGS::diag_jones (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * GREF means gains are reference to the first antenna
	 *
	 * n  = (nant-1)*2 + 1 = 2*nant - 1
	 *
	 * ng = 2*nant
	 *
	 * Usual layout is like this:
	 * layout = R1 I1 R2 I2 R3 I3 R4 I4 ....  Rng     Ing
	 * gainidx= 0     1     2     3     ....  ng-1
	 * index  = 0  1  2  3  4  5  6     ....  2*ng-2  2*ng-1
	 * antidx = 0           1           ....  nant-1
	 *
	 * We are setting I1 = 0 and removing it from the array
	 * layout = R1 R2 I2 R3 I3 R4 I4 ....  Rng     Ing
	 * gainidx= 0  1     2     3     ....  ng-1
	 * index  = 0  1  2  3  4  5  6  ....  2*ng-3  2*ng-2
	 * antidx = 0        1                 nant-1
	 * 
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const data_t *pkg = reinterpret_cast<const data_t*>(instance);

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    ngains ( 2 * pkg->nantennas );

	/* populate complex gains vector */
	vc_type  gains ( ngains, complex_type(0.0f, 0.0f) );
#ifdef GREF
	/* this exponentiation is parameter transformation to ensure realpart is positive */
	gains[0]         = complex_type ( std::exp(rgains[0]), 0.0f );
	for ( int igain = 1; igain < ngains; igain++ ) {
		gains[igain]   = complex_type ( rgains[2*igain-1], rgains[2*igain] );
	}
#else
	for ( int igain = 0; igain < ngains; igain++ ) {
		gains[igain]   = complex_type ( rgains[2*igain+0], rgains[2*igain+1] );
	}
#endif

	/* create complex grad vector */
	vc_type  grad ( ngains, complex_type(0.0f, 0.0f) );

	/* iterate over the polar baselines */
	for ( int ibl = 0; ibl < pkg->npolarbaselines; ibl++ ) {

		// fetch the antenna index
		const int iant1 ( pkg->iant1[ibl] );
		const int iant2 ( pkg->iant2[ibl] );

		// fetch the pb2corr
		const int pb2corr    = pkg->pb2corr [ ibl ];

		// fetch complex data
		const complex_type data ( pkg->data[ibl] );

		// fetch the par corrected model
		const complex_type mrr ( pkg->par_model_rr[ibl] );
		const complex_type mrl ( pkg->par_model_rl[ibl] );
		const complex_type mlr ( pkg->par_model_lr[ibl] );
		const complex_type mll ( pkg->par_model_ll[ibl] );

		// set the gain indices
		const int iprr ( 2*iant1 + 0 );
		const int ipll ( 2*iant1 + 1 );

		const int iqrr ( 2*iant2 + 0 );
		const int iqll ( 2*iant2 + 1 );

		// fetch full gains for both antennas
		const complex_type gprr ( gains[iprr] );
		const complex_type gpll ( gains[ipll] );

		const complex_type gqrr ( gains[iqrr] );
		const complex_type gqll ( gains[iqll] );

		// the following long expressions come from sympy
		// see :math_gradient_parallel.py:
		// see :math_gradient_parallel.code:
		// see :math_gradient_parallel.pdf:

		// model forward depends on pb2corr
		complex_type model;

		// do on every pb2corr
		if ( pb2corr == 0 ) {

const complex_type Dgprr_coeff_dpqrr = -gqrr*conj(mrr) ; 
const complex_type Dgqrr_coeff_dqprr = -gprr*mrr ; 
const complex_type Dgprr_coeff_gprr = gqll*mrl*conj(gqll)*conj(mrl) + gqrr*mrr*conj(gqrr)*conj(mrr) ; 
const complex_type Dgqrr_coeff_gqrr = gpll*mlr*conj(gpll)*conj(mlr) + gprr*mrr*conj(gprr)*conj(mrr) ; 

			grad[iprr] += Dgprr_coeff_gprr*gprr + Dgprr_coeff_dpqrr*data;
			grad[iqrr] += Dgqrr_coeff_gqrr*gqrr + Dgqrr_coeff_dqprr*conj(data);

			model = gprr * mrr * conj(gqrr);

		} // rr
		else if ( pb2corr == 1 ) {

const complex_type Dgprr_coeff_dpqrl = -gqll*conj(mrl) ; 
const complex_type Dgqll_coeff_dqprl = -gprr*mrl ; 

			grad[iprr] += Dgprr_coeff_dpqrl*data; 
			grad[iqll] += Dgqll_coeff_dqprl*conj(data);

			model = gprr * mrl * conj(gqll);

		} // rl
		else if ( pb2corr == 2 ) {

const complex_type Dgpll_coeff_dpqlr = -gqrr*conj(mlr) ; 
const complex_type Dgqrr_coeff_dqplr = -gpll*mlr ; 

			grad[ipll] +=  Dgpll_coeff_dpqlr*data;
			grad[iqrr] +=  Dgqrr_coeff_dqplr*conj(data);

			model = gpll * mlr * conj(gqrr);

		} // lr
		else if ( pb2corr == 3 ) {

const complex_type Dgpll_coeff_dpqll = -gqll*conj(mll) ; 
const complex_type Dgqll_coeff_dqpll = -gpll*mll ; 
const complex_type Dgpll_coeff_gpll = gqll*mll*conj(gqll)*conj(mll) + gqrr*mlr*conj(gqrr)*conj(mlr) ; 
const complex_type Dgqll_coeff_gqll = gpll*mll*conj(gpll)*conj(mll) + gprr*mrl*conj(gprr)*conj(mrl) ; 

			grad[ipll] += Dgpll_coeff_gpll*gpll + Dgpll_coeff_dpqll*data;
			grad[iqll] += Dgqll_coeff_gqll*gqll + Dgqll_coeff_dqpll*conj(data);

			model = gpll * mll * conj(gqll);

		} // ll
			
		// update cost
		cost += std::norm ( data - model );

	} // iterate over polar baselines

#ifdef GREF
	/* load complex grad into real and imaginary parts */
	/* chain rule */
	rgrad[0]        = grad[0].real() * std::exp (rgains[0]);
	for (int igain = 1; igain < ngains; igain++) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain-1] = gg.real();
		rgrad[2*igain] = gg.imag();
	}
#else
	/* load complex grad into real and imaginary parts */
	for ( int igain = 0; igain < ngains; igain++ ) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain + 0] = gg.real();
		rgrad[2*igain + 1] = gg.imag();
	}
#endif

	return cost;
}

int LBFGS::progress_reporter (void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls ) {

	/* get data_t* ptr out of instance */
	data_t *pkg = reinterpret_cast<data_t*>(instance);

	/* save cost and gnorm */
	pkg->cost   = fx;
	pkg->gnorm  = gnorm;
	pkg->niter++;

#ifdef DPRINT
	std::cout << " iteration=" << k << " cost=" << fx << " gnorm=" << gnorm << std::endl;
#endif

	/* return 0 always */
	return 0;
}

int LBFGS::c_progress_reporter (void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls ) {

	/* get data_t* ptr out of instance */
	c_data_t *pkg = static_cast<c_data_t*>(instance);

	/* save cost and gnorm */
	pkg->cost   = fx;
	pkg->gnorm  = gnorm;
	pkg->niter++;

	const real_type iup ( std::exp(x[n-1]) );

#ifdef DPRINT
	std::cout << " iteration=" << k << " cost=" << fx << " gnorm=" << gnorm << " Iup=" << iup << std::endl;
#endif

	/* return 0 always */
	return 0;
}

LBFGS::real_type LBFGS::Solver::solve_full_jones (data_t& pkg) {

	real_type final_cost (0.0f);

	void* vpkg  = static_cast<void*>(&pkg);

	initialize_full_jones ();

  rcode = lbfgs(npar, xpar, &final_cost, full_jones, progress_reporter, vpkg, &param);

  cost   = final_cost;
  gnorm  = pkg.gnorm;
  niter  = pkg.niter;

#ifdef DPRINT
	std::cout << " rcode=" << rcode << std::endl;
#endif

  return final_cost;
}

LBFGS::real_type LBFGS::Solver::solve_diag_jones (data_t& pkg) {

	real_type final_cost (0.0f);

	void* vpkg  = static_cast<void*>(&pkg);

	initialize_diag_jones ();

  //rcode = lbfgs(npar, xpar, &final_cost, diag_jones, progress_reporter, vpkg, &param);
  rcode = lbfgs(npar, xpar, &final_cost, diag_jones_normalized, progress_reporter, vpkg, &param);

  cost   = final_cost;
  gnorm  = pkg.gnorm;
  niter  = pkg.niter;

#ifdef DPRINT
	std::cout << " rcode=" << rcode << std::endl;
#endif

  return final_cost;
}

LBFGS::real_type LBFGS::Solver::solve_full_jones (c_data_t& pkg) {

	real_type final_cost (0.0f);

	void* vpkg  = static_cast<void*>(&pkg);

	initialize_full_jones ();

	/* initialize unpolarized intensity */
	xpar[npar-1] = 3.8;

  rcode = lbfgs(npar, xpar, &final_cost, combined_jones, c_progress_reporter, vpkg, &param);
  //rcode = lbfgs(npar, xpar, &final_cost, combined_jones, nullptr, vpkg, &param);

  cost   = final_cost;
  gnorm  = pkg.gnorm;
  niter  = pkg.niter;

#ifdef DPRINT
	std::cout << " rcode=" << rcode << std::endl;
#endif

  return final_cost;
}

#ifdef GREF
LBFGS::real_type LBFGS::diag_jones_normalized (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * GREF means gains are reference to the first antenna
	 *
	 * n  = (nant-1)*2 + 1 = 2*nant - 1
	 *
	 * ng = 2*nant
	 *
	 * Usual layout is like this:
	 * layout = R1 I1 R2 I2 R3 I3 R4 I4 ....  Rng     Ing
	 * gainidx= 0     1     2     3     ....  ng-1
	 * index  = 0  1  2  3  4  5  6     ....  2*ng-2  2*ng-1
	 * antidx = 0           1           ....  nant-1
	 *
	 * We are setting I1 = 0 and removing it from the array
	 * layout = R1 R2 I2 R3 I3 R4 I4 ....  Rng     Ing
	 * gainidx= 0  1     2     3     ....  ng-1
	 * index  = 0  1  2  3  4  5  6  ....  2*ng-3  2*ng-2
	 * antidx = 0        1                 nant-1
	 * 
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const data_t *pkg = reinterpret_cast<const data_t*>(instance);

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    ngains ( 2 * pkg->nantennas );

	/* populate complex gains vector */
	vc_type  gains ( ngains, complex_type(0.0f, 0.0f) );
	/* this exponentiation is parameter transformation to ensure realpart is positive */
	gains[0]         = complex_type ( std::exp(rgains[0]), 0.0f );
	for ( int igain = 1; igain < ngains; igain++ ) {
		gains[igain]   = complex_type ( rgains[2*igain-1], rgains[2*igain] );
	}

	/* create complex grad vector */
	vc_type  grad ( ngains, complex_type(0.0f, 0.0f) );

	/* iterate over the polar baselines */
	for ( int ibl = 0; ibl < pkg->npolarbaselines; ibl++ ) {

		// fetch the antenna index
		const int iant1 ( pkg->iant1[ibl] );
		const int iant2 ( pkg->iant2[ibl] );

		// fetch the pb2corr
		const int pb2corr    = pkg->pb2corr [ ibl ];

		// fetch complex data
		const complex_type data ( pkg->data[ibl] );

		// fetch the par corrected model
		const complex_type mrr ( pkg->par_model_rr[ibl] );
		const complex_type mrl ( pkg->par_model_rl[ibl] );
		const complex_type mlr ( pkg->par_model_lr[ibl] );
		const complex_type mll ( pkg->par_model_ll[ibl] );

		// set the gain indices
		const int iprr ( 2*iant1 + 0 );
		const int ipll ( 2*iant1 + 1 );

		const int iqrr ( 2*iant2 + 0 );
		const int iqll ( 2*iant2 + 1 );

		// fetch full gains for both antennas
		const complex_type gprr ( gains[iprr] );
		const complex_type gpll ( gains[ipll] );

		const complex_type gqrr ( gains[iqrr] );
		const complex_type gqll ( gains[iqll] );

		// the following long expressions come from sympy
		// see :math_gradient_parallel.py:
		// see :math_gradient_parallel.code:
		// see :math_gradient_parallel.pdf:

		// model forward depends on pb2corr
		complex_type model;

		// do on every pb2corr
		if ( pb2corr == 0 ) {

const complex_type Dgprr_coeff_dpqrr = -gqrr/mrr ; 
const complex_type Dgqrr_coeff_dqprr = -gprr/conj(mrr) ; 
const complex_type Dgprr_coeff_gprr = gqll*conj(gqll) + gqrr*conj(gqrr) ; 
const complex_type Dgqrr_coeff_gqrr = gpll*conj(gpll) + gprr*conj(gprr) ; 

			grad[iprr] += Dgprr_coeff_gprr*gprr + Dgprr_coeff_dpqrr*data;
			grad[iqrr] += Dgqrr_coeff_gqrr*gqrr + Dgqrr_coeff_dqprr*conj(data);

			model = gprr * mrr * conj(gqrr);

			cost += std::norm ( ( data - model ) / mrr );

		} // rr
		else if ( pb2corr == 1 ) {

const complex_type Dgprr_coeff_dpqrl = -gqll/mrl ; 
const complex_type Dgqll_coeff_dqprl = -gprr/conj(mrl) ; 

			grad[iprr] += Dgprr_coeff_dpqrl*data; 
			grad[iqll] += Dgqll_coeff_dqprl*conj(data);

			model = gprr * mrl * conj(gqll);

			cost += std::norm ( ( data - model ) / mrl );

		} // rl
		else if ( pb2corr == 2 ) {

const complex_type Dgpll_coeff_dpqlr = -gqrr/mlr ; 
const complex_type Dgqrr_coeff_dqplr = -gpll/conj(mlr) ; 

			grad[ipll] +=  Dgpll_coeff_dpqlr*data;
			grad[iqrr] +=  Dgqrr_coeff_dqplr*conj(data);

			model = gpll * mlr * conj(gqrr);

			cost += std::norm ( ( data - model ) / mlr );

		} // lr
		else if ( pb2corr == 3 ) {

const complex_type Dgpll_coeff_dpqll = -gqll/mll ; 
const complex_type Dgqll_coeff_dqpll = -gpll/conj(mll) ; 
const complex_type Dgpll_coeff_gpll = gqll*conj(gqll) + gqrr*conj(gqrr) ; 
const complex_type Dgqll_coeff_gqll = gpll*conj(gpll) + gprr*conj(gprr) ; 

			grad[ipll] += Dgpll_coeff_gpll*gpll + Dgpll_coeff_dpqll*data;
			grad[iqll] += Dgqll_coeff_gqll*gqll + Dgqll_coeff_dqpll*conj(data);

			model = gpll * mll * conj(gqll);

			cost += std::norm ( ( data - model ) / mll );

		} // ll
			
	} // iterate over polar baselines

	/* load complex grad into real and imaginary parts */
	/* chain rule */
	rgrad[0]        = grad[0].real() * std::exp (rgains[0]);
	for (int igain = 1; igain < ngains; igain++) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain-1] = gg.real();
		rgrad[2*igain] = gg.imag();
	}

	return cost;
}
#endif
