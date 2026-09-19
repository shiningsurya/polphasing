#include "LBFGSclean.hpp"
/*
 *
 * All phases in the interferometer are relative.
 * Everything is relative. 
 *
 * - Unpolarized full jones with unity amplitude
 *   - Set parallel phases of reference antenna to zero.
 *   - Need to do stepwise:
 *   		(1) Parallel rr and ll
 *   		(2) everything using above as initial solution
 *   		// see if doing everything together maintains 
 *   		- Leakage  rl and lr
*/

//#define DPRINT

#ifdef DPRINT
#include <iostream>
#endif

using std::conj;

int LBFGS::progress_reporter (void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls ) {

	/* get data_t* ptr out of instance */
	uata_t *pkg = reinterpret_cast<uata_t*>(instance);

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

LBFGS::real_type LBFGS::Solver::solve_full_jones_unity_amp (uata_t& pkg) {

	real_type final_cost (0.0f);

	void* vpkg  = static_cast<void*>(&pkg);

	initialize_full_jones ();

  rcode = lbfgs(npar, xpar, &final_cost, full_jones_unity_amp, progress_reporter, vpkg, &param);

  cost   = final_cost;
  gnorm  = pkg.gnorm;
  niter  = pkg.niter;

#ifdef DPRINT
	std::cout << " rcode=" << rcode << std::endl;
#endif

  return final_cost;
}

LBFGS::real_type LBFGS::full_jones_unity_amp (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * phases of the parallel gains of the reference antenna (first antenna) is zero.
	 * ==> real part is always positive
	 * This constraint built into the gradient computation (see GREF)
	 *
	 * Each antenna has four gains. LBFGS requires real parameters. We express the four gains are
	 * __ rgains layout 
	 * R I R I R I R I .....
	 * | one antenna |
	 * |rr,rl, lr, ll|
	 *
	 * GREF layout
	 * R R I R I R | R I R I R I R I | ....
	 * rr
	 *
	 *
	 *
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = reinterpret_cast<const uata_t*>(instance);

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    ngains ( 4 * pkg->nantennas );

	/* populate complex gains vector */
	/* despite GREF, the gains are kept as complex to be uniform */
	vc_type  cgains ( ngains, complex_type(0.0f, 0.0f) );
	// rr
	cgains[0]   = complex_type ( std::exp(rgains[0]), 0.0f );
	// rl
	cgains[1]   = complex_type ( rgains[1], rgains[2] );
	// lr
	cgains[2]   = complex_type ( rgains[3], rgains[4] );
	// ll
	cgains[3]   = complex_type ( std::exp(rgains[5]), 0.0f );
	for ( int igain = 4; igain < ngains; igain++ ) {
		cgains[igain]   = complex_type ( rgains[2*igain-2], rgains[2*igain-1] );
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
			
		//std::cout << " iterationcost=" << cost << " ";
		// update cost
		cost += std::norm ( data - model );

	} // iterate over polar baselines
	
	//std::cout << " rgrads=";

	/* load complex grad into real and imaginary parts */
	// chain rule because of GREF
	// rr
	rgrad[0] = grad[0].real() * std::exp(rgains[0]);
	// rl
	rgrad[1] = grad[1].real();
	rgrad[2] = grad[1].imag();
	// lr
	rgrad[3] = grad[2].real();
	rgrad[4] = grad[2].imag();
	// ll
	rgrad[5] = grad[3].real() * std::exp(rgains[5]);
	// for the rest of grad
	for ( int igain = 4; igain < ngains; igain++ ) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain - 2] = gg.real();
		rgrad[2*igain - 1] = gg.imag();
	}

	//std::cout << " full_jones_cost=" << cost << std::endl; 

	return cost;
}

#if 0
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

  rcode = lbfgs(npar, xpar, &final_cost, diag_jones, progress_reporter, vpkg, &param);

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
#endif
