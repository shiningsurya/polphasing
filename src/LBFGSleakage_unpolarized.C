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
		// rr
		const complex_type rr ( pkg.cgains[4*iant + 0] );

		// ll
		const complex_type ll ( pkg.cgains[4*iant + 3] );

		// drl
		const complex_type drl ( xpar[4*iant + 0], xpar[4*iant + 1] );

		// dlr
		const complex_type dlr ( xpar[4*iant + 2], xpar[4*iant + 3] );

		// rl
		pkg.cgains[4*iant + 1] = rr * drl;
		// lr
		pkg.cgains[4*iant + 2] = ll * dlr;

		// save leakages
		pkg.lgains[2*iant + 0] = drl;
		pkg.lgains[2*iant + 1] = dlr;

	}

  return final_cost;
}

LBFGS::real_type LBFGS::leakage_unpolarized (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * We solve for leakage terms - drl and dlr.
	 *
	 * Our Jones matrix is
	 * | grr   0  | | 1   drl |
	 * | 0    gll | | dlr  1  |
	 * =
	 * | grr  grr*drl |
	 * | gll*dlr  grr |
	 *
	 * The parallel gains come from unity_data_t;
	 *
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = reinterpret_cast<const uata_t*>(instance);

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    nant   ( pkg->nantennas );
	const int    ngains ( 2 * pkg->nantennas );

	/* populate complex gains vector */
	/* this is full gains */
	vc_type  cgains ( pkg->cgains );

	for ( int iant = 0; iant < nant; iant++ ) {

		// xpar is leakages
		const complex_type drl ( rgains[4*iant + 0], rgains[4*iant + 1] );
		const complex_type dlr ( rgains[4*iant + 2], rgains[4*iant + 3] );
		// rl
		// rr * drl
		cgains[4*iant + 1]   = cgains[4*iant + 0] * drl;
		// lr
		// ll * dlr
		cgains[4*iant + 2]   = cgains[4*iant + 3] * dlr;
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
		const int ipll ( 4*iant1 + 3 );

		const int iqrr ( 4*iant2 + 0 );
		const int iqll ( 4*iant2 + 3 );
		
		// leakage indices
		const int iprl ( 2*iant1 + 0 );
		const int iplr ( 2*iant1 + 1 );

		const int iqrl ( 2*iant2 + 0 );
		const int iqlr ( 2*iant2 + 1 );

		// fetch gains for both antennas
		const complex_type gprr ( cgains[iprr] );
		const complex_type gpll ( cgains[ipll] );

		const complex_type gqrr ( cgains[iqrr] );
		const complex_type gqll ( cgains[iqll] );

		// fetch leakages
		const complex_type dprl ( rgains[4*iant1 + 0], rgains[4*iant1 + 1] );
		const complex_type dplr ( rgains[4*iant1 + 2], rgains[4*iant1 + 3] );

		const complex_type dqrl ( rgains[4*iant2 + 0], rgains[4*iant2 + 1] );
		const complex_type dqlr ( rgains[4*iant2 + 2], rgains[4*iant2 + 3] );

		// the following long expressions come from sympy
		// see :math_gradient_leakage_unpolarized.py:
		// see :math_gradient_leakage_unpolarized.code:

		// model forward depends on pb2corr
		complex_type model;

		/*
		We do this over polarbaseline loop, so that we keep track of all the baselines
		*/
		if ( pb2corr == 0 ) {

const complex_type Ddprl_coeff_dpqrr = -dqrl*gqrr*z2*conj(gprr)*conj(z1) ; 
const complex_type Ddqrl_coeff_dqprr = -dprl*gprr*z1*conj(gqrr)*conj(z2) ; 

			grad[iprl] += Ddprl_coeff_dpqrr*data;
			grad[iqrl] += Ddqrl_coeff_dqprr*cata;

model = dprl*gprr*z1*conj(dqrl)*conj(gqrr)*conj(z2) + gprr*z2*conj(gqrr)*conj(z1) ;

		} // rr
		else if ( pb2corr == 1 ) {

const complex_type Ddprl_coeff_dpqrl = -gqll*z2*conj(gprr)*conj(z1) ; 
const complex_type Ddprl_coeff_dqrl = gprr*gqrr*(z2 * z2)*conj(gprr)*conj(gqrr)*(conj(z1) * conj(z1)) ; 
const complex_type Ddprl_coeff_dprl = dqrl*gprr*gqrr*z1*z2*conj(dqrl)*conj(gprr)*conj(gqrr)*conj(z1)*conj(z2) + gprr*gqll*z1*z2*conj(gprr)*conj(gqll)*conj(z1)*conj(z2) ; 

const complex_type Ddqrl_coeff_dprl = gprr*gqrr*(z1 * z1)*conj(gprr)*conj(gqrr)*(conj(z2) * conj(z2)) ; 
const complex_type Ddqrl_coeff_dqrl = dprl*gprr*gqrr*z1*z2*conj(dprl)*conj(gprr)*conj(gqrr)*conj(z1)*conj(z2) + gpll*gqrr*z1*z2*conj(gpll)*conj(gqrr)*conj(z1)*conj(z2) ; 

const complex_type Ddqlr_coeff_dqprl = -gprr*z2*conj(gqll)*conj(z1) ; 

const complex_type Ddprl_coeff_1 = gprr*gqll*(z2 * z2)*conj(dqlr)*conj(gprr)*conj(gqll)*(conj(z1) * conj(z1)) ; 
const complex_type Ddqrl_coeff_1 = gpll*gqrr*(z1 * z1)*conj(dplr)*conj(gpll)*conj(gqrr)*(conj(z2) * conj(z2)) ; 
			
			grad[iprl] += Ddprl_coeff_dpqrl*data + Ddprl_coeff_dqrl*dqrl + Ddprl_coeff_dprl*dprl;
			grad[iqrl] += Ddqrl_coeff_dprl*dprl + Ddqrl_coeff_dqrl*dqrl;

			grad[iqlr] += Ddqlr_coeff_dqprl*cata;

			grad[iprl] += Ddprl_coeff_1;
			grad[iqrl] += Ddqrl_coeff_1;

model = dprl*gprr*z1*conj(gqll)*conj(z2) + gprr*z2*conj(dqlr)*conj(gqll)*conj(z1) ;
		} // rl
		else if ( pb2corr == 2 ) {

const complex_type Ddqrl_coeff_dqplr = -gpll*z1*conj(gqrr)*conj(z2) ; 

const complex_type Ddplr_coeff_dpqlr = -gqrr*z1*conj(gpll)*conj(z2) ; 
const complex_type Ddplr_coeff_dqlr = gpll*gqll*(z1 * z1)*conj(gpll)*conj(gqll)*(conj(z2) * conj(z2)) ; 
const complex_type Ddplr_coeff_dplr = dqlr*gpll*gqll*z1*z2*conj(dqlr)*conj(gpll)*conj(gqll)*conj(z1)*conj(z2) + gpll*gqrr*z1*z2*conj(gpll)*conj(gqrr)*conj(z1)*conj(z2) ; 

const complex_type Ddqlr_coeff_dplr = gpll*gqll*(z2 * z2)*conj(gpll)*conj(gqll)*(conj(z1) * conj(z1)) ; 
const complex_type Ddqlr_coeff_dqlr = dplr*gpll*gqll*z1*z2*conj(dplr)*conj(gpll)*conj(gqll)*conj(z1)*conj(z2) + gprr*gqll*z1*z2*conj(gprr)*conj(gqll)*conj(z1)*conj(z2) ; 

const complex_type Ddplr_coeff_1 = gpll*gqrr*(z1 * z1)*conj(dqrl)*conj(gpll)*conj(gqrr)*(conj(z2) * conj(z2)) ; 
const complex_type Ddqlr_coeff_1 = gprr*gqll*(z2 * z2)*conj(dprl)*conj(gprr)*conj(gqll)*(conj(z1) * conj(z1)) ; 

			grad[iplr] += Ddplr_coeff_dpqlr*data + Ddplr_coeff_dqlr*dqlr + Ddplr_coeff_dplr*dplr;
			grad[iqlr] += Ddqlr_coeff_dplr*dplr + Ddqlr_coeff_dqlr*dqlr;

			grad[iqrl] += Ddqrl_coeff_dqplr*cata;

			grad[iplr] += Ddplr_coeff_1;
			grad[iqlr] += Ddqlr_coeff_1;

model = dplr*gpll*z2*conj(gqrr)*conj(z1) + gpll*z1*conj(dqrl)*conj(gqrr)*conj(z2) ;

		} // lr
		else if ( pb2corr == 3 ) {

const complex_type Ddplr_coeff_dpqll = -dqlr*gqll*z1*conj(gpll)*conj(z2) ; 
const complex_type Ddqlr_coeff_dqpll = -dplr*gpll*z2*conj(gqll)*conj(z1) ; 

			grad[iplr] += Ddplr_coeff_dpqll*data;
			grad[iqlr] += Ddqlr_coeff_dqpll*cata;

model = dplr*gpll*z2*conj(dqlr)*conj(gqll)*conj(z1) + gpll*z1*conj(gqll)*conj(z2) ;

		} // ll
			
		// update cost
		cost += std::norm ( data - model );

	} // iterate over polar baselines

	
	// leakage gradient
	for ( int iant = 0; iant < nant; iant++ ) {
		// rl
		const complex_type grl ( grad[2*iant + 0] );
		rgrad[4*iant + 0] = grl.real();
		rgrad[4*iant + 1] = grl.imag();
		// lr
		const complex_type glr ( grad[2*iant + 1] );
		rgrad[4*iant + 2] = glr.real();
		rgrad[4*iant + 3] = glr.imag();
	}

	return cost;
}

