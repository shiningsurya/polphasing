#include "LBFGSsuperiterative.hpp"

//#define DPRINT

#ifdef DPRINT
#include <iostream>
#endif

LBFGS::real_type LBFGS::unpolarized_rr (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * Solving for rr gain, populates the derivative of rr that appears in {rr,rl,lr} baselines.
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = reinterpret_cast<const uata_t*>(instance);

	/* model I */
	const real_type MI ( pkg->modeli );
	const real_type MI2 ( MI * MI );

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    ngains ( pkg->nantennas );

	/* populate complex gains vector */
	vc_type  cgains ( ngains, complex_type(0.0f, 0.0f) );
	for ( int igain = 0; igain < ngains; igain++ ) {
		cgains[igain]   = complex_type ( rgains[2*igain], rgains[2*igain+1] );
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

		// early continue on ll baseline
		if ( pb2corr == 3 ) continue;

		// fetch complex data
		const complex_type data ( pkg->data[ibl] );
		const complex_type cata ( conj(data) );

		// fetch the parallactic angle skip
		const complex_type z1 ( pkg->par_z1[ibl] );
		const complex_type z2 ( pkg->par_z2[ibl] );

		// set the gain indices
		const int iprr ( iant1 );
		const int iprl ( 4*iant1 + 1 );
		const int iplr ( 4*iant1 + 2 );
		const int ipll ( 4*iant1 + 3 );

		const int iqrr ( iant2 );
		const int iqrl ( 4*iant2 + 1 );
		const int iqlr ( 4*iant2 + 2 );
		const int iqll ( 4*iant2 + 3 );

		// fetch full gains for both antennas
		const complex_type gprr ( cgains[iprr] );
		const complex_type gprl ( pkg->cgains[iprl] );
		const complex_type gplr ( pkg->cgains[iplr] );
		const complex_type gpll ( pkg->cgains[ipll] );

		const complex_type gqrr ( cgains[iqrr] );
		const complex_type gqrl ( pkg->cgains[iqrl] );
		const complex_type gqlr ( pkg->cgains[iqlr] );
		const complex_type gqll ( pkg->cgains[iqll] );

		// the following long expressions come from sympy
		// see :math_gradient_full_unity.py:
		// see :math_gradient_full_unity.code:

		// model forward depends on pb2corr
		complex_type model;

		/*
		We do this over polarbaseline loop, so that we keep track of all the baselines
		*/
		if ( pb2corr == 0 ) {

const complex_type Dgprr_coeff_dpqrr = -MI*gqrr*z1*conj(z2) ; 
const complex_type Dgprr_coeff_gprr = MI2*gqlr*conj(gqlr) + MI2*gqrr*conj(gqrr) ; 
const complex_type Dgqrr_coeff_dqprr = -MI*gprr*z2*conj(z1) ; 
const complex_type Dgqrr_coeff_gqrr = MI2*gplr*conj(gplr) + MI2*gprr*conj(gprr) ; 

			grad[iprr] += Dgprr_coeff_dpqrr*data + Dgprr_coeff_gprr*gprr;

			grad[iqrr] += Dgqrr_coeff_dqprr*cata + Dgqrr_coeff_gqrr*gqrr;

model = gprl*z1*conj(gqrl)*conj(z2) + gprr*z2*conj(gqrr)*conj(z1) ;

		} // rr
		else if ( pb2corr == 1 ) {

const complex_type Dgprr_coeff_dpqrl = -MI*gqlr*z1*conj(z2) ; 
const complex_type Dgprr_coeff_gprl = MI2*gqlr*(z1 * z1)*conj(gqll)*(conj(z2) * conj(z2)) + MI2*gqrr*(z1 * z1)*conj(gqrl)*(conj(z2) * conj(z2)) ; 
const complex_type Dgqrr_coeff_gqrl = MI2*gplr*(z2 * z2)*conj(gpll)*(conj(z1) * conj(z1)) + MI2*gprr*(z2 * z2)*conj(gprl)*(conj(z1) * conj(z1)) ; 

			grad[iprr] += Dgprr_coeff_dpqrl*data + Dgprr_coeff_gprl*gprl;

			grad[iqrr] += Dgqrr_coeff_gqrl*gqrl;

model = gprl*z1*conj(gqll)*conj(z2) + gprr*z2*conj(gqlr)*conj(z1) ;

		} // rl
		else if ( pb2corr == 2 ) {

const complex_type Dgqrr_coeff_dqplr = -MI*gplr*z2*conj(z1) ; 

			grad[iqrr] += Dgqrr_coeff_dqplr*cata;

model = gpll*z1*conj(gqrl)*conj(z2) + gplr*z2*conj(gqrr)*conj(z1) ;

		} // lr
			
		// modeli
		model *= MI;
			
		// update cost
		cost += std::norm ( data - model );

	} // iterate over polar baselines
	
	/* load complex grad into real and imaginary parts */
	for ( int igain = 0; igain < ngains; igain++ ) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain + 0] = gg.real();
		rgrad[2*igain + 1] = gg.imag();
	}

	return cost;
}

LBFGS::real_type LBFGS::unpolarized_rl (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * Solving for rl gain, populates the derivative of rl that appears in {rr,rl,lr} baselines.
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = reinterpret_cast<const uata_t*>(instance);

	/* model I */
	const real_type MI ( pkg->modeli );
	const real_type MI2 ( MI * MI );

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    ngains ( pkg->nantennas );

	/* populate complex gains vector */
	vc_type  cgains ( ngains, complex_type(0.0f, 0.0f) );
	for ( int igain = 0; igain < ngains; igain++ ) {
		cgains[igain]   = complex_type ( rgains[2*igain], rgains[2*igain+1] );
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

		// early continue on ll baseline
		if ( pb2corr == 3 ) continue;

		// fetch complex data
		const complex_type data ( pkg->data[ibl] );
		const complex_type cata ( conj(data) );

		// fetch the parallactic angle skip
		const complex_type z1 ( pkg->par_z1[ibl] );
		const complex_type z2 ( pkg->par_z2[ibl] );

		// set the gain indices
		const int iprr ( 4*iant1 + 0 );
		const int iprl ( iant1 );
		const int iplr ( 4*iant1 + 2 );
		const int ipll ( 4*iant1 + 3 );

		const int iqrr ( 4*iant2 + 0 );
		const int iqrl ( iant2 );
		const int iqlr ( 4*iant2 + 2 );
		const int iqll ( 4*iant2 + 3 );

		// fetch full gains for both antennas
		const complex_type gprr ( pkg->cgains[iprr] );
		const complex_type gprl ( cgains[iprl] );
		const complex_type gplr ( pkg->cgains[iplr] );
		const complex_type gpll ( pkg->cgains[ipll] );

		const complex_type gqrr ( pkg->cgains[iqrr] );
		const complex_type gqrl ( cgains[iqrl] );
		const complex_type gqlr ( pkg->cgains[iqlr] );
		const complex_type gqll ( pkg->cgains[iqll] );

		// the following long expressions come from sympy
		// see :math_gradient_full_unity.py:
		// see :math_gradient_full_unity.code:

		// model forward depends on pb2corr
		complex_type model;

		/*
		We do this over polarbaseline loop, so that we keep track of all the baselines
		*/
		if ( pb2corr == 0 ) {

const complex_type Dgprl_coeff_dpqrr = -MI*gqrl*z2*conj(z1) ; 
const complex_type Dgprl_coeff_gprr = MI2*gqll*(z2 * z2)*conj(gqlr)*(conj(z1) * conj(z1)) + MI2*gqrl*(z2 * z2)*conj(gqrr)*(conj(z1) * conj(z1)) ; 
const complex_type Dgqrl_coeff_dqprr = -MI*gprl*z1*conj(z2) ; 
const complex_type Dgqrl_coeff_gqrr = MI2*gpll*(z1 * z1)*conj(gplr)*(conj(z2) * conj(z2)) + MI2*gprl*(z1 * z1)*conj(gprr)*(conj(z2) * conj(z2)) ; 

			grad[iprl] += Dgprl_coeff_dpqrr*data + Dgprl_coeff_gprr*gprr;

			grad[iqrl] += Dgqrl_coeff_dqprr*cata + Dgqrl_coeff_gqrr*gqrr;

model = gprl*z1*conj(gqrl)*conj(z2) + gprr*z2*conj(gqrr)*conj(z1) ;

		} // rr
		else if ( pb2corr == 1 ) {

const complex_type Dgprl_coeff_dpqrl = -MI*gqll*z2*conj(z1) ; 
const complex_type Dgprl_coeff_gprl = MI2*gqll*conj(gqll) + MI2*gqrl*conj(gqrl) ; 
const complex_type Dgqrl_coeff_gqrl = MI2*gpll*conj(gpll) + MI2*gprl*conj(gprl) ; 

			grad[iprl] += Dgprl_coeff_dpqrl*data + Dgprl_coeff_gprl*gprl;

			grad[iqrl] += Dgqrl_coeff_gqrl*gqrl;

model = gprl*z1*conj(gqll)*conj(z2) + gprr*z2*conj(gqlr)*conj(z1) ;

		} // rl
		else if ( pb2corr == 2 ) {

const complex_type Dgqrl_coeff_dqplr = -MI*gpll*z1*conj(z2) ; 

			grad[iqrl] += Dgqrl_coeff_dqplr*cata;

model = gpll*z1*conj(gqrl)*conj(z2) + gplr*z2*conj(gqrr)*conj(z1) ;

		} // lr

		// scaling model
		model *= MI;
			
		// update cost
		cost += std::norm ( data - model );

	} // iterate over polar baselines
	
	/* load complex grad into real and imaginary parts */
	for ( int igain = 0; igain < ngains; igain++ ) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain + 0] = gg.real();
		rgrad[2*igain + 1] = gg.imag();
	}

	return cost;
}

LBFGS::real_type LBFGS::unpolarized_lr (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * Solving for lr gain, populates the derivative of lr that appears in {rr,rl,lr} baselines.
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = reinterpret_cast<const uata_t*>(instance);

	/* model I */
	const real_type MI ( pkg->modeli );
	const real_type MI2 ( MI * MI );

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    ngains ( pkg->nantennas );

	/* populate complex gains vector */
	vc_type  cgains ( ngains, complex_type(0.0f, 0.0f) );
	for ( int igain = 0; igain < ngains; igain++ ) {
		cgains[igain]   = complex_type ( rgains[2*igain], rgains[2*igain+1] );
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

		// early continue on rr baseline
		if ( pb2corr == 0 ) continue;

		// fetch complex data
		const complex_type data ( pkg->data[ibl] );
		const complex_type cata ( conj(data) );

		// fetch the parallactic angle skip
		const complex_type z1 ( pkg->par_z1[ibl] );
		const complex_type z2 ( pkg->par_z2[ibl] );

		// set the gain indices
		const int iprr ( 4*iant1 + 0 );
		const int iprl ( 4*iant1 + 1 );
		const int iplr ( iant1 );
		const int ipll ( 4*iant1 + 3 );

		const int iqrr ( 4*iant2 + 0 );
		const int iqrl ( 4*iant2 + 1 );
		const int iqlr ( iant2 );
		const int iqll ( 4*iant2 + 3 );

		// fetch full gains for both antennas
		const complex_type gprr ( pkg->cgains[iprr] );
		const complex_type gprl ( pkg->cgains[iprl] );
		const complex_type gplr ( cgains[iplr] );
		const complex_type gpll ( pkg->cgains[ipll] );

		const complex_type gqrr ( pkg->cgains[iqrr] );
		const complex_type gqrl ( pkg->cgains[iqrl] );
		const complex_type gqlr ( cgains[iqlr] );
		const complex_type gqll ( pkg->cgains[iqll] );

		// the following long expressions come from sympy
		// see :math_gradient_full_unity.py:
		// see :math_gradient_full_unity.code:

		// model forward depends on pb2corr
		complex_type model;

		/*
		We do this over polarbaseline loop, so that we keep track of all the baselines
		*/
		if ( pb2corr == 1 ) {

const complex_type Dgqlr_coeff_dqprl = -MI*gprr*z2*conj(z1) ; 

			grad[iqlr] += Dgqlr_coeff_dqprl*cata;

model = gprl*z1*conj(gqll)*conj(z2) + gprr*z2*conj(gqlr)*conj(z1) ;

		} // rl
		else if ( pb2corr == 2 ) {

const complex_type Dgplr_coeff_dpqlr = -MI*gqrr*z1*conj(z2) ; 
const complex_type Dgplr_coeff_gplr = MI2*gqlr*conj(gqlr) + MI2*gqrr*conj(gqrr) ; 
const complex_type Dgqlr_coeff_gqlr = MI2*gplr*conj(gplr) + MI2*gprr*conj(gprr) ; 

			grad[iplr] += Dgplr_coeff_dpqlr*data + Dgplr_coeff_gplr*gplr;

			grad[iqlr] += Dgqlr_coeff_gqlr*gqlr;

model = gpll*z1*conj(gqrl)*conj(z2) + gplr*z2*conj(gqrr)*conj(z1) ;

		} // lr
		else if ( pb2corr == 3 ) {

const complex_type Dgplr_coeff_dpqll = -MI*gqlr*z1*conj(z2) ; 
const complex_type Dgplr_coeff_gpll = MI2*gqlr*(z1 * z1)*conj(gqll)*(conj(z2) * conj(z2)) + MI2*gqrr*(z1 * z1)*conj(gqrl)*(conj(z2) * conj(z2)) ; 
const complex_type Dgqlr_coeff_dqpll = -MI*gplr*z2*conj(z1) ; 
const complex_type Dgqlr_coeff_gqll = MI2*gplr*(z2 * z2)*conj(gpll)*(conj(z1) * conj(z1)) + MI2*gprr*(z2 * z2)*conj(gprl)*(conj(z1) * conj(z1)) ; 

			grad[iplr] += Dgplr_coeff_dpqll*data + Dgplr_coeff_gpll*gpll;

			grad[iqlr] += Dgqlr_coeff_dqpll*cata + Dgqlr_coeff_gqll*gqll;

model = gpll*z1*conj(gqll)*conj(z2) + gplr*z2*conj(gqlr)*conj(z1) ;

		} // ll
			
		// modeli
		model *= MI;
			
		//std::cout << " iterationcost=" << cost << " ";
		// update cost
		cost += std::norm ( data - model );

	} // iterate over polar baselines
	
	/* load complex grad into real and imaginary parts */
	for ( int igain = 0; igain < ngains; igain++ ) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain + 0] = gg.real();
		rgrad[2*igain + 1] = gg.imag();
	}

	return cost;
}

LBFGS::real_type LBFGS::unpolarized_ll (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * Solving for ll gain, populates the derivative of ll that appears in {rl,lr,ll} baselines.
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = reinterpret_cast<const uata_t*>(instance);

	/* model I */
	const real_type MI ( pkg->modeli );
	const real_type MI2 ( MI * MI );

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    ngains ( pkg->nantennas );

	/* populate complex gains vector */
	vc_type  cgains ( ngains, complex_type(0.0f, 0.0f) );
	for ( int igain = 0; igain < ngains; igain++ ) {
		cgains[igain]   = complex_type ( rgains[2*igain], rgains[2*igain+1] );
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

		// early continue with rr
		if (pb2corr == 0) continue;

		// fetch complex data
		const complex_type data ( pkg->data[ibl] );
		const complex_type cata ( conj(data) );

		// fetch the parallactic angle skip
		const complex_type z1 ( pkg->par_z1[ibl] );
		const complex_type z2 ( pkg->par_z2[ibl] );

		// set the gain indices
		const int iprr ( 4*iant1 + 0 );
		const int iprl ( 4*iant1 + 1 );
		const int iplr ( 4*iant1 + 2 );
		const int ipll ( iant1 );

		const int iqrr ( 4*iant2 + 0 );
		const int iqrl ( 4*iant2 + 1 );
		const int iqlr ( 4*iant2 + 2 );
		const int iqll ( iant2 );

		// fetch full gains for both antennas
		const complex_type gprr ( pkg->cgains[iprr] );
		const complex_type gprl ( pkg->cgains[iprl] );
		const complex_type gplr ( pkg->cgains[iplr] );
		const complex_type gpll ( cgains[ipll] );

		const complex_type gqrr ( pkg->cgains[iqrr] );
		const complex_type gqrl ( pkg->cgains[iqrl] );
		const complex_type gqlr ( pkg->cgains[iqlr] );
		const complex_type gqll ( cgains[iqll] );

		// the following long expressions come from sympy
		// see :math_gradient_full_unity.py:
		// see :math_gradient_full_unity.code:

		// model forward depends on pb2corr
		complex_type model;

		/*
		We do this over polarbaseline loop, so that we keep track of all the baselines
		*/
		if ( pb2corr == 1 ) {

const complex_type Dgqll_coeff_dqprl = -MI*gprl*z1*conj(z2) ; 

			grad[iqll] += Dgqll_coeff_dqprl*cata;

model = gprl*z1*conj(gqll)*conj(z2) + gprr*z2*conj(gqlr)*conj(z1) ;

		} // rl
		else if ( pb2corr == 2 ) {

const complex_type Dgpll_coeff_dpqlr = -MI*gqrl*z2*conj(z1) ; 
const complex_type Dgpll_coeff_gplr = MI2*gqll*(z2 * z2)*conj(gqlr)*(conj(z1) * conj(z1)) + MI2*gqrl*(z2 * z2)*conj(gqrr)*(conj(z1) * conj(z1)) ; 
const complex_type Dgqll_coeff_gqlr = MI2*gpll*(z1 * z1)*conj(gplr)*(conj(z2) * conj(z2)) + MI2*gprl*(z1 * z1)*conj(gprr)*(conj(z2) * conj(z2)) ; 

			grad[ipll] += Dgpll_coeff_dpqlr*data + Dgpll_coeff_gplr*gplr;

			grad[iqll] += Dgqll_coeff_gqlr*gqlr;

model = gpll*z1*conj(gqrl)*conj(z2) + gplr*z2*conj(gqrr)*conj(z1) ;

		} // lr
		else if ( pb2corr == 3 ) {

const complex_type Dgpll_coeff_dpqll = -MI*gqll*z2*conj(z1) ; 
const complex_type Dgpll_coeff_gpll = MI2*gqll*conj(gqll) + MI2*gqrl*conj(gqrl) ; 
const complex_type Dgqll_coeff_dqpll = -MI*gpll*z1*conj(z2) ; 
const complex_type Dgqll_coeff_gqll = MI2*gpll*conj(gpll) + MI2*gprl*conj(gprl) ; 

			grad[ipll] += Dgpll_coeff_dpqll*data + Dgpll_coeff_gpll*gpll;

			grad[iqll] += Dgqll_coeff_dqpll*cata + Dgqll_coeff_gqll*gqll;

model = gpll*z1*conj(gqll)*conj(z2) + gplr*z2*conj(gqlr)*conj(z1) ;

		} // ll
			
		// modeli
		model *= MI;

		// update cost
		cost += std::norm ( data - model );

	} // iterate over polar baselines
	
	/* load complex grad into real and imaginary parts */
	for ( int igain = 0; igain < ngains; igain++ ) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain + 0] = gg.real();
		rgrad[2*igain + 1] = gg.imag();
	}

	return cost;
}

LBFGS::cost_t LBFGS::IterativeUnpolarizedSolver::solve (uata_t& pkg, int n ) {

	cost_t    final_cost;

	void* vpkg  = static_cast<void*>(&pkg);

	for (int __i = 0; __i < n; __i++) {
		/* load gains from pkg into xpars */
		loadgains_rr ( pkg );
		loadgains_ll ( pkg );

		/* solve for rr and ll */
		rcode_rr = lbfgs(n_rr, xpar_rr, &cost_rr, unpolarized_rr, NULL, vpkg, &param);
		rcode_ll = lbfgs(n_ll, xpar_ll, &cost_ll, unpolarized_ll, NULL, vpkg, &param);

		/* save solved gains into pkg.cgains */
		savegains_rr ( pkg );
		savegains_ll ( pkg );

		/* solve for rl and lr */
		rcode_rl = lbfgs(n_rl, xpar_rl, &cost_rl, unpolarized_rl, NULL, vpkg, &param);
		rcode_lr = lbfgs(n_lr, xpar_lr, &cost_lr, unpolarized_lr, NULL, vpkg, &param);

		/* save solved gains into pkg.cgains */
		savegains_rl ( pkg );
		savegains_lr ( pkg );

		final_cost[0]  = cost_rr;
		final_cost[1]  = cost_rl;
		final_cost[2]  = cost_lr;
		final_cost[3]  = cost_ll;

#ifdef DPRINT
	std::cout << "iter=" << __i << " costs rr=" << cost_rr << " rl=" << cost_rl << " lr=" << cost_lr << " ll=" << cost_ll << std::endl;
#endif
		
	} // solving iterations

  return final_cost;
}
