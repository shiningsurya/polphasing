#include "LBFGSclean.hpp"

#ifdef DPRINT
#include <iostream>
#endif

LBFGS::real_type LBFGS::PolarizedSolver::solve_diag_polarized (data_t& pkg) {
	/*
	 * First solve parallel and then solve full
	 * using the parallel solve as initial for parallel gains
	 *
	 * While ensuring the phases of the parallel gains of the reference antenna are zero
	*/

	real_type cost (0.0f);

	void* vpkg  = static_cast<void*>(&pkg);

	// parallel solve 
	initialize_parallel ();
  rcode_para = lbfgs(n_para, xpar_para, &cost, parallel_polarized, NULL, vpkg, &param);

	gnorm_para  = pkg.gnorm;
	niter_para  = pkg.niter;

#ifdef DPRINT
	std::cout << " rcode=" << rcode_para << std::endl;
#endif

  return cost;
}

LBFGS::real_type LBFGS::parallel_polarized (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
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
		const complex_type mrr ( pkg->parleak_model_rr[ibl] );
		const complex_type mrl ( pkg->parleak_model_rl[ibl] );
		const complex_type mlr ( pkg->parleak_model_lr[ibl] );
		const complex_type mll ( pkg->parleak_model_ll[ibl] );

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
