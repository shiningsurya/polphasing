#include "LBFGSclean.hpp"

#ifdef DPRINT
#include <iostream>
#endif

using std::conj;

/*
 * Iterative Objective functions
*/

/*
 * Solve for rr
 *
 * gains are forevery antenna {grr, grl}
*/
lbfgsfloatval_t LBFGS::rr_unpolarized (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = static_cast<const uata_t*>(instance);

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    nant ( pkg->nantennas );
	const int    ngains ( 2 * pkg->nantennas );

	/* populate complex gains vector */
	vc_type  cgains ( ngains, complex_type(0.0f, 0.0f) );
	for ( int igain = 0; igain < ngains; igain++ ) {
		cgains[igain]   = complex_type ( rgains[2*igain  ], rgains[2*igain+1] );
	}

	/* create complex grad vector */
	vc_type  grad ( ngains, complex_type(0.0f, 0.0f) );

	/* iterate over the polar baselines */
	for ( int ibl = 0; ibl < pkg->npolarbaselines; ibl++ ) {

		// fetch the pb2corr
		const int pb2corr    = pkg->pb2corr [ ibl ];

		// we are only interested in rr
		if (pb2corr != 0) continue;

		// fetch the antenna index
		const int iant1 ( pkg->iant1[ibl] );
		const int iant2 ( pkg->iant2[ibl] );

		// set the gain indices
		const int iprr ( 2*iant1 + 0 );
		const int iprl ( 2*iant1 + 1 );
		const int iqrr ( 2*iant2 + 0 );
		const int iqrl ( 2*iant2 + 1 );

		// fetch complex data
		const complex_type data ( pkg->data[ibl] );
		const complex_type cata ( conj(data) );

		// fetch the parallactic angle skip
		const complex_type z1 ( pkg->par_z1[ibl] );
		const complex_type z2 ( pkg->par_z2[ibl] );

		// fetch gains
		const complex_type gprr ( cgains[iprr] );
		const complex_type gprl ( cgains[iprl] );
		const complex_type gqrr ( cgains[iqrr] );
		const complex_type gqrl ( cgains[iqrl] );

		// the following long expressions come from sympy
		// see :math_gradient_full_unity_iterative.py:
		// see :math_gradient_full_unity_iterative.code:

		// model forward depends on pb2corr
		complex_type model ( 0.0f, 0.0f );

		/*
			We do this over polarbaseline loop, so that we keep track of all the baselines
			pb2corr == 0
		*/

const complex_type Dgprr_coeff_dpqrr = -gqrr*z1*conj(z2) ; 
const complex_type Dgprr_coeff_gprr = gqrr*z1*z2*conj(gqrr)*conj(z1)*conj(z2) ; 
const complex_type Dgprr_coeff_gprl = gqrr*(z1 * z1)*conj(gqrl)*(conj(z2) * conj(z2)) ; 

const complex_type Dgqrr_coeff_dqprr = -gprr*z2*conj(z1) ; 
const complex_type Dgqrr_coeff_gqrl = gprr*(z2 * z2)*conj(gprl)*(conj(z1) * conj(z1)) ; 
const complex_type Dgqrr_coeff_gqrr = gprr*z1*z2*conj(gprr)*conj(z1)*conj(z2) ; 

const complex_type Dgprl_coeff_dpqrr = -gqrl*z2*conj(z1) ; 
const complex_type Dgprl_coeff_gprr = gqrl*(z2 * z2)*conj(gqrr)*(conj(z1) * conj(z1)) ; 
const complex_type Dgprl_coeff_gprl = gqrl*z1*z2*conj(gqrl)*conj(z1)*conj(z2) ; 

const complex_type Dgqrl_coeff_dqprr = -gprl*z1*conj(z2) ; 
const complex_type Dgqrl_coeff_gqrr = gprl*(z1 * z1)*conj(gprr)*(conj(z2) * conj(z2)) ; 
const complex_type Dgqrl_coeff_gqrl = gprl*z1*z2*conj(gprl)*conj(z1)*conj(z2) ; 

			grad[iprr] += Dgprr_coeff_dpqrr*data + Dgprr_coeff_gprr*gprr + Dgprr_coeff_gprl*gprl;

			grad[iqrr] += Dgqrr_coeff_dqprr*cata + Dgqrr_coeff_gqrr*gqrr + Dgqrr_coeff_gqrl*gqrl;

			grad[iprl] += Dgprl_coeff_dpqrr*data + Dgprl_coeff_gprr*gprr + Dgprl_coeff_gprl*gprl;

			grad[iqrl] += Dgqrl_coeff_dqprr*cata + Dgqrl_coeff_gqrr*gqrr + Dgqrl_coeff_gqrl*gqrl;

model = gprl*z1*conj(gqrl)*conj(z2) + gprr*z2*conj(gqrr)*conj(z1) ;

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

lbfgsfloatval_t LBFGS::ll_unpolarized (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = static_cast<const uata_t*>(instance);

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    nant ( pkg->nantennas );
	const int    ngains ( 2 * pkg->nantennas );

	/* populate complex gains vector */
	vc_type  cgains ( ngains, complex_type(0.0f, 0.0f) );
	for ( int igain = 0; igain < ngains; igain++ ) {
		cgains[igain]   = complex_type ( rgains[2*igain  ], rgains[2*igain+1] );
	}

	/* create complex grad vector */
	vc_type  grad ( ngains, complex_type(0.0f, 0.0f) );

	/* iterate over the polar baselines */
	for ( int ibl = 0; ibl < pkg->npolarbaselines; ibl++ ) {

		// fetch the pb2corr
		const int pb2corr    = pkg->pb2corr [ ibl ];

		// we are only interested in ll
		if (pb2corr != 3) continue;

		// fetch the antenna index
		const int iant1 ( pkg->iant1[ibl] );
		const int iant2 ( pkg->iant2[ibl] );

		// set the gain indices
		const int iplr ( 2*iant1 + 0 );
		const int ipll ( 2*iant1 + 1 );
		const int iqlr ( 2*iant2 + 0 );
		const int iqll ( 2*iant2 + 1 );

		// fetch complex data
		const complex_type data ( pkg->data[ibl] );
		const complex_type cata ( conj(data) );

		// fetch the parallactic angle skip
		const complex_type z1 ( pkg->par_z1[ibl] );
		const complex_type z2 ( pkg->par_z2[ibl] );

		// fetch gains
		const complex_type gplr ( cgains[iplr] );
		const complex_type gpll ( cgains[ipll] );
		const complex_type gqlr ( cgains[iqlr] );
		const complex_type gqll ( cgains[iqll] );

		// the following long expressions come from sympy
		// see :math_gradient_full_unity_iterative.py:
		// see :math_gradient_full_unity_iterative.code:

		// model forward depends on pb2corr
		complex_type model ( 0.0f, 0.0f );

		/*
			We do this over polarbaseline loop, so that we keep track of all the baselines

			this is pb2corr==3
		*/
const complex_type Dgplr_coeff_dpqll = -gqlr*z1*conj(z2) ; 
const complex_type Dgplr_coeff_gpll = gqlr*(z1 * z1)*conj(gqll)*(conj(z2) * conj(z2)) ; 
const complex_type Dgplr_coeff_gplr = gqlr*z1*z2*conj(gqlr)*conj(z1)*conj(z2) ; 

const complex_type Dgqlr_coeff_dqpll = -gplr*z2*conj(z1) ; 
const complex_type Dgqlr_coeff_gqll = gplr*(z2 * z2)*conj(gpll)*(conj(z1) * conj(z1)) ; 
const complex_type Dgqlr_coeff_gqlr = gplr*z1*z2*conj(gplr)*conj(z1)*conj(z2) ; 

const complex_type Dgpll_coeff_dpqll = -gqll*z2*conj(z1) ; 
const complex_type Dgpll_coeff_gplr = gqll*(z2 * z2)*conj(gqlr)*(conj(z1) * conj(z1)) ; 
const complex_type Dgpll_coeff_gpll = gqll*z1*z2*conj(gqll)*conj(z1)*conj(z2) ; 

const complex_type Dgqll_coeff_dqpll = -gpll*z1*conj(z2) ; 
const complex_type Dgqll_coeff_gqlr = gpll*(z1 * z1)*conj(gplr)*(conj(z2) * conj(z2)) ; 
const complex_type Dgqll_coeff_gqll = gpll*z1*z2*conj(gpll)*conj(z1)*conj(z2) ; 

			grad[iplr] += Dgplr_coeff_dpqll*data + Dgplr_coeff_gpll*gpll + Dgplr_coeff_gplr*gplr;
			grad[ipll] += Dgpll_coeff_dpqll*data + Dgpll_coeff_gpll*gpll + Dgpll_coeff_gpll*gpll;

			grad[iqlr] += Dgqlr_coeff_dqpll*cata + Dgqlr_coeff_gqll*gqll + Dgqlr_coeff_gqlr*gqlr;
			grad[iqll] += Dgqll_coeff_dqpll*cata + Dgqll_coeff_gqll*gqll + Dgqll_coeff_gqlr*gqlr;

model = gpll*z1*conj(gqll)*conj(z2) + gplr*z2*conj(gqlr)*conj(z1) ;
			
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

lbfgsfloatval_t LBFGS::rllr_unpolarized (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = static_cast<const uata_t*>(instance);

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    nant ( pkg->nantennas );
	const int    ngains ( 4 * pkg->nantennas );

	/* populate complex gains vector */
	vc_type  cgains ( ngains, complex_type(0.0f, 0.0f) );
	for ( int igain = 0; igain < ngains; igain++ ) {
		cgains[igain]   = complex_type ( rgains[2*igain  ], rgains[2*igain+1] );
	}

	/* create complex grad vector */
	vc_type  grad ( ngains, complex_type(0.0f, 0.0f) );

	/* iterate over the polar baselines */
	for ( int ibl = 0; ibl < pkg->npolarbaselines; ibl++ ) {

		// fetch the pb2corr
		const int pb2corr    = pkg->pb2corr [ ibl ];

		// we are only interested in ll
		if (pb2corr != 1 || pb2corr != 2) continue;

		// fetch the antenna index
		const int iant1 ( pkg->iant1[ibl] );
		const int iant2 ( pkg->iant2[ibl] );

		// set the gain indices
		const int iprr ( 4*iant1 + 0 );
		const int iprl ( 4*iant1 + 1 );
		const int iplr ( 4*iant1 + 2 );
		const int ipll ( 4*iant1 + 3 );

		const int iqrr ( 4*iant2 + 0 );
		const int iqrl ( 4*iant2 + 1 );
		const int iqlr ( 4*iant2 + 2 );
		const int iqll ( 4*iant2 + 3 );

		// fetch complex data
		const complex_type data ( pkg->data[ibl] );
		const complex_type cata ( conj(data) );

		// fetch the parallactic angle skip
		const complex_type z1 ( pkg->par_z1[ibl] );
		const complex_type z2 ( pkg->par_z2[ibl] );

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
		// see :math_gradient_full_unity.py:
		// see :math_gradient_full_unity.code:

		// model forward depends on pb2corr
		complex_type model ( 0.0f, 0.0f );

		/*
			We do this over polarbaseline loop, so that we keep track of all the baselines
			which at this point in code is either 1 or 2.
		*/
		if ( pb2corr == 1 ) {

const complex_type Dgprr_coeff_dpqrl = -gqlr*z1*conj(z2) ; 
const complex_type Dgprr_coeff_gprl = gqlr*(z1 * z1)*conj(gqll)*(conj(z2) * conj(z2)) ; 
const complex_type Dgprr_coeff_gprr = gqlr*z1*z2*conj(gqlr)*conj(z1)*conj(z2) ; 

const complex_type Dgprl_coeff_dpqrl = -gqll*z2*conj(z1) ; 
const complex_type Dgprl_coeff_gprr = gqll*(z2 * z2)*conj(gqlr)*(conj(z1) * conj(z1)) ; 
const complex_type Dgprl_coeff_gprl = gqll*z1*z2*conj(gqll)*conj(z1)*conj(z2) ; 

const complex_type Dgqlr_coeff_dqprl = -gprr*z2*conj(z1) ; 
const complex_type Dgqlr_coeff_gqll = gprr*(z2 * z2)*conj(gprl)*(conj(z1) * conj(z1)) ; 
const complex_type Dgqlr_coeff_gqlr = gprr*z1*z2*conj(gprr)*conj(z1)*conj(z2) ; 

const complex_type Dgqll_coeff_dqprl = -gprl*z1*conj(z2) ; 
const complex_type Dgqll_coeff_gqlr = gprl*(z1 * z1)*conj(gprr)*(conj(z2) * conj(z2)) ; 
const complex_type Dgqll_coeff_gqll = gprl*z1*z2*conj(gprl)*conj(z1)*conj(z2) ; 

			grad[iprr] += Dgprr_coeff_dpqrl*data + Dgprr_coeff_gprl*gprl + Dgprr_coeff_gprr*gprr;
			grad[iprl] += Dgprl_coeff_dpqrl*data + Dgprl_coeff_gprl*gprl + Dgprl_coeff_gprl*gprl;

			grad[iqlr] += Dgqlr_coeff_dqprl*cata + Dgqlr_coeff_gqll*gqll + Dgqlr_coeff_gqlr*gqlr;
			grad[iqll] += Dgqll_coeff_dqprl*cata + Dgqll_coeff_gqll*gqll + Dgqll_coeff_gqlr*gqlr;

model = gprl*z1*conj(gqll)*conj(z2) + gprr*z2*conj(gqlr)*conj(z1) ;

		} // rl
		else if ( pb2corr == 2 ) {

const complex_type Dgqrr_coeff_dqplr = -gplr*z2*conj(z1) ; 
const complex_type Dgqrr_coeff_gqrl = gplr*(z2 * z2)*conj(gpll)*(conj(z1) * conj(z1)) ; 
const complex_type Dgqrr_coeff_gqrr = gplr*z1*z2*conj(gplr)*conj(z1)*conj(z2) ; 

const complex_type Dgqrl_coeff_dqplr = -gpll*z1*conj(z2) ; 
const complex_type Dgqrl_coeff_gqrr = gpll*(z1 * z1)*conj(gplr)*(conj(z2) * conj(z2)) ; 
const complex_type Dgqrl_coeff_gqrl = gpll*z1*z2*conj(gpll)*conj(z1)*conj(z2) ; 

const complex_type Dgplr_coeff_dpqlr = -gqrr*z1*conj(z2) ; 
const complex_type Dgplr_coeff_gpll = gqrr*(z1 * z1)*conj(gqrl)*(conj(z2) * conj(z2)) ; 
const complex_type Dgplr_coeff_gplr = gqrr*z1*z2*conj(gqrr)*conj(z1)*conj(z2) ; 

const complex_type Dgpll_coeff_dpqlr = -gqrl*z2*conj(z1) ; 
const complex_type Dgpll_coeff_gplr = gqrl*(z2 * z2)*conj(gqrr)*(conj(z1) * conj(z1)) ; 
const complex_type Dgpll_coeff_gpll = gqrl*z1*z2*conj(gqrl)*conj(z1)*conj(z2) ; 

			grad[iplr] += Dgplr_coeff_dpqlr*data + Dgplr_coeff_gplr*gplr + Dgplr_coeff_gpll*gpll;
			grad[ipll] += Dgpll_coeff_dpqlr*data + Dgpll_coeff_gplr*gplr + Dgpll_coeff_gpll*gpll;

			grad[iqrr] += Dgqrr_coeff_dqplr*cata + Dgqrr_coeff_gqrr*gqrr + Dgqrr_coeff_gqrl*gqrl;
			grad[iqrl] += Dgqrl_coeff_dqplr*cata + Dgqrl_coeff_gqrr*gqrr + Dgqrl_coeff_gqrl*gqrl;

model = gpll*z1*conj(gqrl)*conj(z2) + gplr*z2*conj(gqrr)*conj(z1) ;

		} // lr
			
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

LBFGS::real_type LBFGS::IterativeUnpolarizedSolver::solve (uata_t& pkg) {
	/*
	 * One iteration consists of 
	 * - `rr` and `ll` done independently
	 * - `rllr`
	*/

	real_type final_cost (0.0f);

	void* vpkg  = static_cast<void*>(&pkg);

	for (int __i = 0; __i < 3; __i++) {
		/* load gains from pkg into xpars */
		loadgains_rr ( pkg );
		loadgains_ll ( pkg );

		/* solve for rr and ll */
		rcode_rr = lbfgs(n_rr, xpar_rr, &cost_rr, rr_unpolarized, NULL, vpkg, &param);
		rcode_ll = lbfgs(n_ll, xpar_ll, &cost_ll, ll_unpolarized, NULL, vpkg, &param);

		/* save solved gains into pkg.cgains */
		savegains_rr ( pkg );
		savegains_ll ( pkg );

		/* load gains from pkg into xpar_full */
		loadgains_full ( pkg );

		/* solve for rl and lr */
		rcode_full = lbfgs(n_full, xpar_full, &cost_full, rllr_unpolarized, NULL, vpkg, &param);

		/* save solved gains into pkg.cgains */
		savegains_full ( pkg );

		final_cost   = cost_rr + cost_ll + cost_full;
#ifdef DPRINT
	std::cout << " costs rr=" << cost_rr << " ll=" << cost_ll << " full=" << cost_full << std::endl;
#endif
		
	} // solving iterations

  //gnorm  = pkg.gnorm;
  //niter  = pkg.niter;

  return final_cost;
}
