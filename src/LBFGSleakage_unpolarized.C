#include "LBFGSclean.hpp"

/*
 * Unpolarized leakage
*/
#ifdef DPRINT
#include <iostream>
#endif

LBFGS::real_type LBFGS::UnpolarizedLeakageSolver::solve_leakage_unpolarized (uata_t& pkg) {
	/*
	 * First solve parallel and then solve only for leakages
	 * keeping the parallel solved gains as fixed
	*/

	real_type final_cost (0.0f);

	void* vpkg  = static_cast<void*>(&pkg);

	{
		// single unpolarized solving 
		SingleUnpolarizedSolver single_r ( nantennas );
		SingleUnpolarizedSolver single_l ( nantennas );

		single_r.solve_single_unpolarized<0> ( pkg );
		single_l.solve_single_unpolarized<3> ( pkg );

		final_cost += single_r.cost;
		final_cost += single_l.cost;

		// load parallel solved gains into uata_t
		// xpar      = {RI RI RI RI}
		//              01 23 45 67
		for (int iant = 0; iant < nantennas; iant++) {
			// rr
			pkg.cgains[4*iant + 0] = complex_type ( single_r.xpar[2*iant + 0], single_r.xpar[2*iant + 1] );

			// ll
			pkg.cgains[4*iant + 3] = complex_type ( single_l.xpar[2*iant + 0], single_l.xpar[2*iant + 1] );
		}
	}

	// full solve 
	real_type cost_full (0.0f);
	rcode   = lbfgs(n, xpar, &cost_full, leakage_unpolarized, NULL, vpkg, &param);

  final_cost  += cost_full;
  //gnorm  = pkg.gnorm;
  //niter  = pkg.niter;

#ifdef DPRINT
	std::cout << " rcode=" << rcode << std::endl;
#endif

	// load cross gains into cgains
	for (int iant = 0; iant < nantennas; iant++) {
		// rl
		pkg.cgains[4*iant + 1] = complex_type ( xpar[4*iant + 0], xpar[4*iant + 1] );

		// lr
		pkg.cgains[4*iant + 2] = complex_type ( xpar[4*iant + 2], xpar[4*iant + 3] );
	}

  return final_cost;
}

LBFGS::real_type LBFGS::leakage_unpolarized (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * Each antenna has four gains. LBFGS requires real parameters. We express the four gains are
	 * __ rgains layout 
	 * R I    R I .....
	 * | one antenna |
	 * |rl, lr|
	 *
	 * We fix the parallel gains as is, and only solve for cross gains.
	 *
	 * The parallel gains must come from unity_data_t;
	 *
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = reinterpret_cast<const uata_t*>(instance);

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    nant   ( pkg->nantennas );
	const int    ngains ( 4 * pkg->nantennas );

	/* populate complex gains vector */
	/* this is full gains */
	vc_type  cgains ( pkg->cgains );

	for ( int iant = 0; iant < nant; iant++ ) {
		// rl
		cgains[4*iant + 1]   = complex_type ( rgains[4*iant+0], rgains[4*iant+1] );
		// lr
		cgains[4*iant + 2]   = complex_type ( rgains[4*iant+2], rgains[4*iant+3] );
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

		// fetch the parallactic angle skip
		const complex_type z1 ( pkg->par_z1[ibl] );
		const complex_type z2 ( pkg->par_z2[ibl] );

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
		// see :math_gradient_full_unity.py:
		// see :math_gradient_full_unity.code:

		// model forward depends on pb2corr
		complex_type model;

		/*
		We do this over polarbaseline loop, so that we keep track of all the baselines
		*/
		if ( pb2corr == 0 ) {

const complex_type Dgprr_coeff_dpqrr = -gqrr*z1*conj(z2) ; 
const complex_type Dgprr_coeff_gprr = gqlr*z1*z2*conj(gqlr)*conj(z1)*conj(z2) + gqrr*z1*z2*conj(gqrr)*conj(z1)*conj(z2) ; 
const complex_type Dgqrr_coeff_dqprr = -gprr*z2*conj(z1) ; 
const complex_type Dgqrr_coeff_gqrr = gplr*z1*z2*conj(gplr)*conj(z1)*conj(z2) + gprr*z1*z2*conj(gprr)*conj(z1)*conj(z2) ; 
const complex_type Dgprl_coeff_dpqrr = -gqrl*z2*conj(z1) ; 
const complex_type Dgprl_coeff_gprr = gqll*(z2 * z2)*conj(gqlr)*(conj(z1) * conj(z1)) + gqrl*(z2 * z2)*conj(gqrr)*(conj(z1) * conj(z1)) ; 
const complex_type Dgqrl_coeff_dqprr = -gprl*z1*conj(z2) ; 
const complex_type Dgqrl_coeff_gqrr = gpll*(z1 * z1)*conj(gplr)*(conj(z2) * conj(z2)) + gprl*(z1 * z1)*conj(gprr)*(conj(z2) * conj(z2)) ; 

			grad[iprr] += Dgprr_coeff_dpqrr*data + Dgprr_coeff_gprr*gprr;
			grad[iprl] += Dgprl_coeff_dpqrr*data + Dgprl_coeff_gprr*gprr;

			grad[iqrr] += Dgqrr_coeff_dqprr*cata + Dgqrr_coeff_gqrr*gqrr;
			grad[iqrl] += Dgqrl_coeff_dqprr*cata + Dgqrl_coeff_gqrr*gqrr;

model = gprl*z1*conj(gqrl)*conj(z2) + gprr*z2*conj(gqrr)*conj(z1) ;

		} // rr
		else if ( pb2corr == 1 ) {

const complex_type Dgprr_coeff_dpqrl = -gqlr*z1*conj(z2) ; 
const complex_type Dgprr_coeff_gprl = gqlr*(z1 * z1)*conj(gqll)*(conj(z2) * conj(z2)) + gqrr*(z1 * z1)*conj(gqrl)*(conj(z2) * conj(z2)) ; 
const complex_type Dgqrr_coeff_gqrl = gplr*(z2 * z2)*conj(gpll)*(conj(z1) * conj(z1)) + gprr*(z2 * z2)*conj(gprl)*(conj(z1) * conj(z1)) ; 
const complex_type Dgprl_coeff_dpqrl = -gqll*z2*conj(z1) ; 
const complex_type Dgprl_coeff_gprl = gqll*z1*z2*conj(gqll)*conj(z1)*conj(z2) + gqrl*z1*z2*conj(gqrl)*conj(z1)*conj(z2) ; 
const complex_type Dgqrl_coeff_gqrl = gpll*z1*z2*conj(gpll)*conj(z1)*conj(z2) + gprl*z1*z2*conj(gprl)*conj(z1)*conj(z2) ; 
const complex_type Dgqlr_coeff_dqprl = -gprr*z2*conj(z1) ; 
const complex_type Dgqll_coeff_dqprl = -gprl*z1*conj(z2) ; 

			grad[iprr] += Dgprr_coeff_dpqrl*data + Dgprr_coeff_gprl*gprl;
			grad[iprl] += Dgprl_coeff_dpqrl*data + Dgprl_coeff_gprl*gprl;

			grad[iqlr] += Dgqlr_coeff_dqprl*cata;
			grad[iqll] += Dgqll_coeff_dqprl*cata;

			grad[iqrr] += Dgqrr_coeff_gqrl*gqrl;
			grad[iqrl] += Dgqrl_coeff_gqrl*gqrl;

model = gprl*z1*conj(gqll)*conj(z2) + gprr*z2*conj(gqlr)*conj(z1) ;

		} // rl
		else if ( pb2corr == 2 ) {

const complex_type Dgqrr_coeff_dqplr = -gplr*z2*conj(z1) ; 
const complex_type Dgqrl_coeff_dqplr = -gpll*z1*conj(z2) ; 
const complex_type Dgplr_coeff_dpqlr = -gqrr*z1*conj(z2) ; 
const complex_type Dgplr_coeff_gplr = gqlr*z1*z2*conj(gqlr)*conj(z1)*conj(z2) + gqrr*z1*z2*conj(gqrr)*conj(z1)*conj(z2) ; 
const complex_type Dgqlr_coeff_gqlr = gplr*z1*z2*conj(gplr)*conj(z1)*conj(z2) + gprr*z1*z2*conj(gprr)*conj(z1)*conj(z2) ; 
const complex_type Dgpll_coeff_dpqlr = -gqrl*z2*conj(z1) ; 
const complex_type Dgpll_coeff_gplr = gqll*(z2 * z2)*conj(gqlr)*(conj(z1) * conj(z1)) + gqrl*(z2 * z2)*conj(gqrr)*(conj(z1) * conj(z1)) ; 
const complex_type Dgqll_coeff_gqlr = gpll*(z1 * z1)*conj(gplr)*(conj(z2) * conj(z2)) + gprl*(z1 * z1)*conj(gprr)*(conj(z2) * conj(z2)) ; 

			grad[iplr] += Dgplr_coeff_dpqlr*data + Dgplr_coeff_gplr*gplr;
			grad[ipll] += Dgpll_coeff_dpqlr*data + Dgpll_coeff_gplr*gplr;

			grad[iqrr] += Dgqrr_coeff_dqplr*cata;
			grad[iqrl] += Dgqrl_coeff_dqplr*cata;

			grad[iqlr] += Dgqlr_coeff_gqlr*gqlr;
			grad[iqll] += Dgqll_coeff_gqlr*gqlr;

model = gpll*z1*conj(gqrl)*conj(z2) + gplr*z2*conj(gqrr)*conj(z1) ;

		} // lr
		else if ( pb2corr == 3 ) {

const complex_type Dgplr_coeff_dpqll = -gqlr*z1*conj(z2) ; 
const complex_type Dgplr_coeff_gpll = gqlr*(z1 * z1)*conj(gqll)*(conj(z2) * conj(z2)) + gqrr*(z1 * z1)*conj(gqrl)*(conj(z2) * conj(z2)) ; 
const complex_type Dgqlr_coeff_dqpll = -gplr*z2*conj(z1) ; 
const complex_type Dgqlr_coeff_gqll = gplr*(z2 * z2)*conj(gpll)*(conj(z1) * conj(z1)) + gprr*(z2 * z2)*conj(gprl)*(conj(z1) * conj(z1)) ; 
const complex_type Dgpll_coeff_dpqll = -gqll*z2*conj(z1) ; 
const complex_type Dgpll_coeff_gpll = gqll*z1*z2*conj(gqll)*conj(z1)*conj(z2) + gqrl*z1*z2*conj(gqrl)*conj(z1)*conj(z2) ; 
const complex_type Dgqll_coeff_dqpll = -gpll*z1*conj(z2) ; 
const complex_type Dgqll_coeff_gqll = gpll*z1*z2*conj(gpll)*conj(z1)*conj(z2) + gprl*z1*z2*conj(gprl)*conj(z1)*conj(z2) ; 

			grad[iplr] += Dgplr_coeff_dpqll*data + Dgplr_coeff_gpll*gpll;
			grad[ipll] += Dgpll_coeff_dpqll*data + Dgpll_coeff_gpll*gpll;

			grad[iqlr] += Dgqlr_coeff_dqpll*cata + Dgqlr_coeff_gqll*gqll;
			grad[iqll] += Dgqll_coeff_dqpll*cata + Dgqll_coeff_gqll*gqll;

model = gpll*z1*conj(gqll)*conj(z2) + gplr*z2*conj(gqlr)*conj(z1) ;

		} // ll
			
		// update cost
		cost += std::norm ( data - model );

	} // iterate over polar baselines

	
	// we ignore grad[rr,ll] as we are only
	// fitting for rl and lr gains
	for ( int iant = 0; iant < nant; iant++ ) {
		// rl
		const complex_type grl ( grad[4*iant + 1] );
		rgrad[4*iant + 0] = grl.real();
		rgrad[4*iant + 1] = grl.imag();
		// lr
		const complex_type glr ( grad[4*iant + 2] );
		rgrad[4*iant + 2] = glr.real();
		rgrad[4*iant + 3] = glr.imag();
	}

	return cost;
}

