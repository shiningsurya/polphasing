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
	initialize ();
  rcode  = lbfgs(n, xpar, &cost, parallel_polarized, NULL, vpkg, &param);

	gnorm  = pkg.gnorm;
	niter  = pkg.niter;

#ifdef DPRINT
	std::cout << " rcode=" << rcode_para << std::endl;
#endif

  return cost;
}

LBFGS::real_type LBFGS::parallel_polarized (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * xpar has 2*complex gains for antennas + crosshand phase
	 *
	 * RI RI RI RI ... RI | xphase
	 * 
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const data_t *pkg = reinterpret_cast<const data_t*>(instance);

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    ngains ( 2 * pkg->nantennas );
	const int    pindex ( 4 * pkg->nantennas );

	/* populate complex gains vector */
	vc_type  gains ( ngains, complex_type(0.0f, 0.0f) );
	for ( int igain = 0; igain < ngains; igain++ ) {
		gains[igain]   = complex_type ( rgains[2*igain], rgains[2*igain+1] );
	}

	const real_type xphase ( rgains[pindex] );

	const complex_type zp  ( std::cos(xphase), std::sin(xphase) );
	const complex_type czp ( std::cos(xphase), -std::sin(xphase) );
	const complex_type zp2 ( std::cos(2.0f*xphase), std::sin(2.0f*xphase) );
	const complex_type czp2( std::cos(2.0f*xphase), std::sin(-2.0f*xphase) );

	// purely imaginary I
	const complex_type I ( 0.0f, 1.0f );

	/* create grad vector */
	vc_type  grad ( ngains, complex_type(0.0f, 0.0f) );
	// keeping complex because math is complex
	// but the gradient will be purely real
	complex_type xpgrad ( 0.0f, 0.0f );

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
		const int iprr ( 2*iant1 + 0 );
		const int ipll ( 2*iant1 + 1 );

		const int iqrr ( 2*iant2 + 0 );
		const int iqll ( 2*iant2 + 1 );

		// fetch the leakage terms
		const complex_type pdrl  ( pkg->pleakrl[ibl] );
		const complex_type pdlr  ( pkg->pleaklr[ibl] );

		const complex_type qdrl  ( pkg->qleakrl[ibl] );
		const complex_type qdlr  ( pkg->qleaklr[ibl] );

		// fetch full gains for both antennas
		const complex_type gprr ( gains[iprr] );
		const complex_type gpll ( gains[ipll] );

		const complex_type gqrr ( gains[iqrr] );
		const complex_type gqll ( gains[iqll] );

		// the following long expressions come from sympy
		// see :math_gradient_parallel_polarized_leakage.py:
		// see :math_gradient_parallel_polarized_leakage.code:

		// model forward depends on pb2corr
		complex_type model;

		// do on every pb2corr
		if ( pb2corr == 0 ) {

const complex_type Dgprr_coeff_dpqrr = -gqrr*qdrl*conj(mll)*conj(pdrl) - gqrr*qdrl*czp*conj(mrl) - gqrr*zp*conj(mlr)*conj(pdrl) - gqrr*conj(mrr) ; 
const complex_type Dgqrr_coeff_dqprr = -gprr*mll*pdrl*conj(qdrl) - gprr*mlr*pdrl*czp - gprr*mrl*zp*conj(qdrl) - gprr*mrr ; 

const complex_type Dgprr_coeff_gprr = gqll*mll*pdrl*qdlr*zp*conj(gqll)*conj(mlr)*conj(pdrl) + gqll*mll*pdrl*qdlr*conj(gqll)*conj(mrr) + gqll*mll*pdrl*conj(gqll)*conj(mll)*conj(pdrl) + gqll*mll*pdrl*czp*conj(gqll)*conj(mrl) + gqll*mlr*pdrl*qdlr*conj(gqll)*conj(mlr)*conj(pdrl)*conj(qdlr) + gqll*mlr*pdrl*qdlr*czp*conj(gqll)*conj(mrr)*conj(qdlr) + gqll*mlr*pdrl*czp*conj(gqll)*conj(mll)*conj(pdrl)*conj(qdlr) + gqll*mlr*pdrl*czp2*conj(gqll)*conj(mrl)*conj(qdlr) + gqll*mrl*qdlr*zp2*conj(gqll)*conj(mlr)*conj(pdrl) + gqll*mrl*qdlr*zp*conj(gqll)*conj(mrr) + gqll*mrl*zp*conj(gqll)*conj(mll)*conj(pdrl) + gqll*mrl*conj(gqll)*conj(mrl) + gqll*mrr*qdlr*zp*conj(gqll)*conj(mlr)*conj(pdrl)*conj(qdlr) + gqll*mrr*qdlr*conj(gqll)*conj(mrr)*conj(qdlr) + gqll*mrr*conj(gqll)*conj(mll)*conj(pdrl)*conj(qdlr) + gqll*mrr*czp*conj(gqll)*conj(mrl)*conj(qdlr) + gqrr*mll*pdrl*qdrl*conj(gqrr)*conj(mll)*conj(pdrl)*conj(qdrl) + gqrr*mll*pdrl*qdrl*czp*conj(gqrr)*conj(mrl)*conj(qdrl) + gqrr*mll*pdrl*zp*conj(gqrr)*conj(mlr)*conj(pdrl)*conj(qdrl) + gqrr*mll*pdrl*conj(gqrr)*conj(mrr)*conj(qdrl) + gqrr*mlr*pdrl*qdrl*czp*conj(gqrr)*conj(mll)*conj(pdrl) + gqrr*mlr*pdrl*qdrl*czp2*conj(gqrr)*conj(mrl) + gqrr*mlr*pdrl*conj(gqrr)*conj(mlr)*conj(pdrl) + gqrr*mlr*pdrl*czp*conj(gqrr)*conj(mrr) + gqrr*mrl*qdrl*zp*conj(gqrr)*conj(mll)*conj(pdrl)*conj(qdrl) + gqrr*mrl*qdrl*conj(gqrr)*conj(mrl)*conj(qdrl) + gqrr*mrl*zp2*conj(gqrr)*conj(mlr)*conj(pdrl)*conj(qdrl) + gqrr*mrl*zp*conj(gqrr)*conj(mrr)*conj(qdrl) + gqrr*mrr*qdrl*conj(gqrr)*conj(mll)*conj(pdrl) + gqrr*mrr*qdrl*czp*conj(gqrr)*conj(mrl) + gqrr*mrr*zp*conj(gqrr)*conj(mlr)*conj(pdrl) + gqrr*mrr*conj(gqrr)*conj(mrr) ; 

const complex_type Dgqrr_coeff_gqrr = gpll*mll*qdrl*conj(gpll)*conj(mll)*conj(qdrl) + gpll*mll*qdrl*czp*conj(gpll)*conj(mrl)*conj(pdlr)*conj(qdrl) + gpll*mll*zp*conj(gpll)*conj(mlr)*conj(qdrl) + gpll*mll*conj(gpll)*conj(mrr)*conj(pdlr)*conj(qdrl) + gpll*mlr*qdrl*czp*conj(gpll)*conj(mll) + gpll*mlr*qdrl*czp2*conj(gpll)*conj(mrl)*conj(pdlr) + gpll*mlr*conj(gpll)*conj(mlr) + gpll*mlr*czp*conj(gpll)*conj(mrr)*conj(pdlr) + gpll*mrl*pdlr*qdrl*zp*conj(gpll)*conj(mll)*conj(qdrl) + gpll*mrl*pdlr*qdrl*conj(gpll)*conj(mrl)*conj(pdlr)*conj(qdrl) + gpll*mrl*pdlr*zp2*conj(gpll)*conj(mlr)*conj(qdrl) + gpll*mrl*pdlr*zp*conj(gpll)*conj(mrr)*conj(pdlr)*conj(qdrl) + gpll*mrr*pdlr*qdrl*conj(gpll)*conj(mll) + gpll*mrr*pdlr*qdrl*czp*conj(gpll)*conj(mrl)*conj(pdlr) + gpll*mrr*pdlr*zp*conj(gpll)*conj(mlr) + gpll*mrr*pdlr*conj(gpll)*conj(mrr)*conj(pdlr) + gprr*mll*pdrl*qdrl*conj(gprr)*conj(mll)*conj(pdrl)*conj(qdrl) + gprr*mll*pdrl*qdrl*czp*conj(gprr)*conj(mrl)*conj(qdrl) + gprr*mll*pdrl*zp*conj(gprr)*conj(mlr)*conj(pdrl)*conj(qdrl) + gprr*mll*pdrl*conj(gprr)*conj(mrr)*conj(qdrl) + gprr*mlr*pdrl*qdrl*czp*conj(gprr)*conj(mll)*conj(pdrl) + gprr*mlr*pdrl*qdrl*czp2*conj(gprr)*conj(mrl) + gprr*mlr*pdrl*conj(gprr)*conj(mlr)*conj(pdrl) + gprr*mlr*pdrl*czp*conj(gprr)*conj(mrr) + gprr*mrl*qdrl*zp*conj(gprr)*conj(mll)*conj(pdrl)*conj(qdrl) + gprr*mrl*qdrl*conj(gprr)*conj(mrl)*conj(qdrl) + gprr*mrl*zp2*conj(gprr)*conj(mlr)*conj(pdrl)*conj(qdrl) + gprr*mrl*zp*conj(gprr)*conj(mrr)*conj(qdrl) + gprr*mrr*qdrl*conj(gprr)*conj(mll)*conj(pdrl) + gprr*mrr*qdrl*czp*conj(gprr)*conj(mrl) + gprr*mrr*zp*conj(gprr)*conj(mlr)*conj(pdrl) + gprr*mrr*conj(gprr)*conj(mrr) ; 


const complex_type Dphi_coeff_dpqrr = I*gqrr*qdrl*czp*conj(gprr)*conj(mrl) - I*gqrr*zp*conj(gprr)*conj(mlr)*conj(pdrl) ; 
const complex_type Dphi_coeff_dqprr = I*gprr*mlr*pdrl*czp*conj(gqrr) - I*gprr*mrl*zp*conj(gqrr)*conj(qdrl) ; 
const complex_type Dphi_coeff_zp = I*gpll*gqll*mll*qdlr*conj(gpll)*conj(gqll)*conj(mlr) + I*gpll*gqll*mrl*pdlr*qdlr*conj(gpll)*conj(gqll)*conj(mrr)*conj(pdlr) + I*gpll*gqll*mrl*pdlr*conj(gpll)*conj(gqll)*conj(mll) + I*gpll*gqll*mrr*pdlr*qdlr*conj(gpll)*conj(gqll)*conj(mlr)*conj(qdlr) + I*gpll*gqrr*mll*conj(gpll)*conj(gqrr)*conj(mlr)*conj(qdrl) + I*gpll*gqrr*mrl*pdlr*qdrl*conj(gpll)*conj(gqrr)*conj(mll)*conj(qdrl) + I*gpll*gqrr*mrl*pdlr*conj(gpll)*conj(gqrr)*conj(mrr)*conj(pdlr)*conj(qdrl) + I*gpll*gqrr*mrr*pdlr*conj(gpll)*conj(gqrr)*conj(mlr) + I*gprr*gqll*mll*pdrl*qdlr*conj(gprr)*conj(gqll)*conj(mlr)*conj(pdrl) + I*gprr*gqll*mrl*qdlr*conj(gprr)*conj(gqll)*conj(mrr) + I*gprr*gqll*mrl*conj(gprr)*conj(gqll)*conj(mll)*conj(pdrl) + I*gprr*gqll*mrr*qdlr*conj(gprr)*conj(gqll)*conj(mlr)*conj(pdrl)*conj(qdlr) + I*gprr*gqrr*mll*pdrl*conj(gprr)*conj(gqrr)*conj(mlr)*conj(pdrl)*conj(qdrl) + I*gprr*gqrr*mrl*qdrl*conj(gprr)*conj(gqrr)*conj(mll)*conj(pdrl)*conj(qdrl) + I*gprr*gqrr*mrl*conj(gprr)*conj(gqrr)*conj(mrr)*conj(qdrl) + I*gprr*gqrr*mrr*conj(gprr)*conj(gqrr)*conj(mlr)*conj(pdrl) ; 

			xpgrad  += Dphi_coeff_dpqrr*data;
			xpgrad  += Dphi_coeff_dqprr*cata;
			xpgrad  += Dphi_coeff_zp*zp;

			grad[iprr] += Dgprr_coeff_gprr*gprr + Dgprr_coeff_dpqrr*data;
			grad[iqrr] += Dgqrr_coeff_gqrr*gqrr + Dgqrr_coeff_dqprr*cata;

			model = gprr*mll*pdrl*conj(gqrr)*conj(qdrl) + gprr*mlr*pdrl*czp*conj(gqrr) + gprr*mrl*zp*conj(gqrr)*conj(qdrl) + gprr*mrr*conj(gqrr);

		} // rr
		else if ( pb2corr == 1 ) {

const complex_type Dgprr_coeff_dpqrl = -gqll*qdlr*zp*conj(mlr)*conj(pdrl) - gqll*qdlr*conj(mrr) - gqll*conj(mll)*conj(pdrl) - gqll*czp*conj(mrl) ; 
const complex_type Dgqll_coeff_dqprl = -gprr*mll*pdrl - gprr*mlr*pdrl*czp*conj(qdlr) - gprr*mrl*zp - gprr*mrr*conj(qdlr) ; 

const complex_type Dphi_coeff_dpqrl = -I*gqll*qdlr*zp*conj(gprr)*conj(mlr)*conj(pdrl) + I*gqll*czp*conj(gprr)*conj(mrl) ; 
const complex_type Dphi_coeff_dqprl = I*gprr*mlr*pdrl*czp*conj(gqll)*conj(qdlr) - I*gprr*mrl*zp*conj(gqll) ; 
const complex_type Dphi_coeff_zp2 = 2.0f*I*gpll*gqll*mrl*pdlr*qdlr*conj(gpll)*conj(gqll)*conj(mlr) + 2.0f*I*gpll*gqrr*mrl*pdlr*conj(gpll)*conj(gqrr)*conj(mlr)*conj(qdrl) + 2.0f*I*gprr*gqll*mrl*qdlr*conj(gprr)*conj(gqll)*conj(mlr)*conj(pdrl) + 2.0f*I*gprr*gqrr*mrl*conj(gprr)*conj(gqrr)*conj(mlr)*conj(pdrl)*conj(qdrl) ; 

			xpgrad += Dphi_coeff_dpqrl*data;
			xpgrad += Dphi_coeff_dqprl*cata;
			xpgrad += Dphi_coeff_zp2*zp2;

			grad[iprr] += Dgprr_coeff_dpqrl*data; 
			grad[iqll] += Dgqll_coeff_dqprl*cata;

			model = gprr*mll*pdrl*conj(gqll) + gprr*mlr*pdrl*czp*conj(gqll)*conj(qdlr) + gprr*mrl*zp*conj(gqll) + gprr*mrr*conj(gqll)*conj(qdlr);


		} // rl
		else if ( pb2corr == 2 ) {

const complex_type Dgpll_coeff_dpqlr = -gqrr*qdrl*conj(mll) - gqrr*qdrl*czp*conj(mrl)*conj(pdlr) - gqrr*zp*conj(mlr) - gqrr*conj(mrr)*conj(pdlr) ; 
const complex_type Dgqrr_coeff_dqplr = -gpll*mll*conj(qdrl) - gpll*mlr*czp - gpll*mrl*pdlr*zp*conj(qdrl) - gpll*mrr*pdlr ; 

const complex_type Dphi_coeff_dpqlr = I*gqrr*qdrl*czp*conj(gpll)*conj(mrl)*conj(pdlr) - I*gqrr*zp*conj(gpll)*conj(mlr) ; 
const complex_type Dphi_coeff_dqplr = I*gpll*mlr*czp*conj(gqrr) - I*gpll*mrl*pdlr*zp*conj(gqrr)*conj(qdrl) ; 
const complex_type Dphi_coeff_czp2 = -2.0f*I*gpll*gqll*mlr*conj(gpll)*conj(gqll)*conj(mrl)*conj(pdlr)*conj(qdlr) - 2.0f*I*gpll*gqrr*mlr*qdrl*conj(gpll)*conj(gqrr)*conj(mrl)*conj(pdlr) - 2.0f*I*gprr*gqll*mlr*pdrl*conj(gprr)*conj(gqll)*conj(mrl)*conj(qdlr) - 2.0f*I*gprr*gqrr*mlr*pdrl*qdrl*conj(gprr)*conj(gqrr)*conj(mrl) ; 

			xpgrad += Dphi_coeff_dpqlr*data;
			xpgrad += Dphi_coeff_dqplr*cata;
			xpgrad += Dphi_coeff_czp2*czp2;

			grad[ipll] +=  Dgpll_coeff_dpqlr*data;
			grad[iqrr] +=  Dgqrr_coeff_dqplr*cata;

			model = gpll*mll*conj(gqrr)*conj(qdrl) + gpll*mlr*czp*conj(gqrr) + gpll*mrl*pdlr*zp*conj(gqrr)*conj(qdrl) + gpll*mrr*pdlr*conj(gqrr);


		} // lr
		else if ( pb2corr == 3 ) {

const complex_type Dgpll_coeff_dpqll = -gqll*qdlr*zp*conj(mlr) - gqll*qdlr*conj(mrr)*conj(pdlr) - gqll*conj(mll) - gqll*czp*conj(mrl)*conj(pdlr) ; 
const complex_type Dgqll_coeff_dqpll = -gpll*mll - gpll*mlr*czp*conj(qdlr) - gpll*mrl*pdlr*zp - gpll*mrr*pdlr*conj(qdlr) ; 

const complex_type Dgpll_coeff_gpll = gqll*mll*qdlr*zp*conj(gqll)*conj(mlr) + gqll*mll*qdlr*conj(gqll)*conj(mrr)*conj(pdlr) + gqll*mll*conj(gqll)*conj(mll) + gqll*mll*czp*conj(gqll)*conj(mrl)*conj(pdlr) + gqll*mlr*qdlr*conj(gqll)*conj(mlr)*conj(qdlr) + gqll*mlr*qdlr*czp*conj(gqll)*conj(mrr)*conj(pdlr)*conj(qdlr) + gqll*mlr*czp*conj(gqll)*conj(mll)*conj(qdlr) + gqll*mlr*czp2*conj(gqll)*conj(mrl)*conj(pdlr)*conj(qdlr) + gqll*mrl*pdlr*qdlr*zp2*conj(gqll)*conj(mlr) + gqll*mrl*pdlr*qdlr*zp*conj(gqll)*conj(mrr)*conj(pdlr) + gqll*mrl*pdlr*zp*conj(gqll)*conj(mll) + gqll*mrl*pdlr*conj(gqll)*conj(mrl)*conj(pdlr) + gqll*mrr*pdlr*qdlr*zp*conj(gqll)*conj(mlr)*conj(qdlr) + gqll*mrr*pdlr*qdlr*conj(gqll)*conj(mrr)*conj(pdlr)*conj(qdlr) + gqll*mrr*pdlr*conj(gqll)*conj(mll)*conj(qdlr) + gqll*mrr*pdlr*czp*conj(gqll)*conj(mrl)*conj(pdlr)*conj(qdlr) + gqrr*mll*qdrl*conj(gqrr)*conj(mll)*conj(qdrl) + gqrr*mll*qdrl*czp*conj(gqrr)*conj(mrl)*conj(pdlr)*conj(qdrl) + gqrr*mll*zp*conj(gqrr)*conj(mlr)*conj(qdrl) + gqrr*mll*conj(gqrr)*conj(mrr)*conj(pdlr)*conj(qdrl) + gqrr*mlr*qdrl*czp*conj(gqrr)*conj(mll) + gqrr*mlr*qdrl*czp2*conj(gqrr)*conj(mrl)*conj(pdlr) + gqrr*mlr*conj(gqrr)*conj(mlr) + gqrr*mlr*czp*conj(gqrr)*conj(mrr)*conj(pdlr) + gqrr*mrl*pdlr*qdrl*zp*conj(gqrr)*conj(mll)*conj(qdrl) + gqrr*mrl*pdlr*qdrl*conj(gqrr)*conj(mrl)*conj(pdlr)*conj(qdrl) + gqrr*mrl*pdlr*zp2*conj(gqrr)*conj(mlr)*conj(qdrl) + gqrr*mrl*pdlr*zp*conj(gqrr)*conj(mrr)*conj(pdlr)*conj(qdrl) + gqrr*mrr*pdlr*qdrl*conj(gqrr)*conj(mll) + gqrr*mrr*pdlr*qdrl*czp*conj(gqrr)*conj(mrl)*conj(pdlr) + gqrr*mrr*pdlr*zp*conj(gqrr)*conj(mlr) + gqrr*mrr*pdlr*conj(gqrr)*conj(mrr)*conj(pdlr) ; 
const complex_type Dgqll_coeff_gqll = gpll*mll*qdlr*zp*conj(gpll)*conj(mlr) + gpll*mll*qdlr*conj(gpll)*conj(mrr)*conj(pdlr) + gpll*mll*conj(gpll)*conj(mll) + gpll*mll*czp*conj(gpll)*conj(mrl)*conj(pdlr) + gpll*mlr*qdlr*conj(gpll)*conj(mlr)*conj(qdlr) + gpll*mlr*qdlr*czp*conj(gpll)*conj(mrr)*conj(pdlr)*conj(qdlr) + gpll*mlr*czp*conj(gpll)*conj(mll)*conj(qdlr) + gpll*mlr*czp2*conj(gpll)*conj(mrl)*conj(pdlr)*conj(qdlr) + gpll*mrl*pdlr*qdlr*zp2*conj(gpll)*conj(mlr) + gpll*mrl*pdlr*qdlr*zp*conj(gpll)*conj(mrr)*conj(pdlr) + gpll*mrl*pdlr*zp*conj(gpll)*conj(mll) + gpll*mrl*pdlr*conj(gpll)*conj(mrl)*conj(pdlr) + gpll*mrr*pdlr*qdlr*zp*conj(gpll)*conj(mlr)*conj(qdlr) + gpll*mrr*pdlr*qdlr*conj(gpll)*conj(mrr)*conj(pdlr)*conj(qdlr) + gpll*mrr*pdlr*conj(gpll)*conj(mll)*conj(qdlr) + gpll*mrr*pdlr*czp*conj(gpll)*conj(mrl)*conj(pdlr)*conj(qdlr) + gprr*mll*pdrl*qdlr*zp*conj(gprr)*conj(mlr)*conj(pdrl) + gprr*mll*pdrl*qdlr*conj(gprr)*conj(mrr) + gprr*mll*pdrl*conj(gprr)*conj(mll)*conj(pdrl) + gprr*mll*pdrl*czp*conj(gprr)*conj(mrl) + gprr*mlr*pdrl*qdlr*conj(gprr)*conj(mlr)*conj(pdrl)*conj(qdlr) + gprr*mlr*pdrl*qdlr*czp*conj(gprr)*conj(mrr)*conj(qdlr) + gprr*mlr*pdrl*czp*conj(gprr)*conj(mll)*conj(pdrl)*conj(qdlr) + gprr*mlr*pdrl*czp2*conj(gprr)*conj(mrl)*conj(qdlr) + gprr*mrl*qdlr*zp2*conj(gprr)*conj(mlr)*conj(pdrl) + gprr*mrl*qdlr*zp*conj(gprr)*conj(mrr) + gprr*mrl*zp*conj(gprr)*conj(mll)*conj(pdrl) + gprr*mrl*conj(gprr)*conj(mrl) + gprr*mrr*qdlr*zp*conj(gprr)*conj(mlr)*conj(pdrl)*conj(qdlr) + gprr*mrr*qdlr*conj(gprr)*conj(mrr)*conj(qdlr) + gprr*mrr*conj(gprr)*conj(mll)*conj(pdrl)*conj(qdlr) + gprr*mrr*czp*conj(gprr)*conj(mrl)*conj(qdlr) ; 

const complex_type Dphi_coeff_dpqll = -I*gqll*qdlr*zp*conj(gpll)*conj(mlr) + I*gqll*czp*conj(gpll)*conj(mrl)*conj(pdlr) ; 
const complex_type Dphi_coeff_dqpll = I*gpll*mlr*czp*conj(gqll)*conj(qdlr) - I*gpll*mrl*pdlr*zp*conj(gqll) ; 
const complex_type Dphi_coeff_czp = -I*gpll*gqll*mll*conj(gpll)*conj(gqll)*conj(mrl)*conj(pdlr) - I*gpll*gqll*mlr*qdlr*conj(gpll)*conj(gqll)*conj(mrr)*conj(pdlr)*conj(qdlr) - I*gpll*gqll*mlr*conj(gpll)*conj(gqll)*conj(mll)*conj(qdlr) - I*gpll*gqll*mrr*pdlr*conj(gpll)*conj(gqll)*conj(mrl)*conj(pdlr)*conj(qdlr) - I*gpll*gqrr*mll*qdrl*conj(gpll)*conj(gqrr)*conj(mrl)*conj(pdlr)*conj(qdrl) - I*gpll*gqrr*mlr*qdrl*conj(gpll)*conj(gqrr)*conj(mll) - I*gpll*gqrr*mlr*conj(gpll)*conj(gqrr)*conj(mrr)*conj(pdlr) - I*gpll*gqrr*mrr*pdlr*qdrl*conj(gpll)*conj(gqrr)*conj(mrl)*conj(pdlr) - I*gprr*gqll*mll*pdrl*conj(gprr)*conj(gqll)*conj(mrl) - I*gprr*gqll*mlr*pdrl*qdlr*conj(gprr)*conj(gqll)*conj(mrr)*conj(qdlr) - I*gprr*gqll*mlr*pdrl*conj(gprr)*conj(gqll)*conj(mll)*conj(pdrl)*conj(qdlr) - I*gprr*gqll*mrr*conj(gprr)*conj(gqll)*conj(mrl)*conj(qdlr) - I*gprr*gqrr*mll*pdrl*qdrl*conj(gprr)*conj(gqrr)*conj(mrl)*conj(qdrl) - I*gprr*gqrr*mlr*pdrl*qdrl*conj(gprr)*conj(gqrr)*conj(mll)*conj(pdrl) - I*gprr*gqrr*mlr*pdrl*conj(gprr)*conj(gqrr)*conj(mrr) - I*gprr*gqrr*mrr*qdrl*conj(gprr)*conj(gqrr)*conj(mrl) ; 

			xpgrad += Dphi_coeff_dpqll*data;
			xpgrad += Dphi_coeff_dqpll*cata;
			xpgrad += Dphi_coeff_czp*czp;

			grad[ipll] += Dgpll_coeff_gpll*gpll + Dgpll_coeff_dpqll*data;
			grad[iqll] += Dgqll_coeff_gqll*gqll + Dgqll_coeff_dqpll*conj(data);

			model = gpll*mll*conj(gqll) + gpll*mlr*czp*conj(gqll)*conj(qdlr) + gpll*mrl*pdlr*zp*conj(gqll) + gpll*mrr*pdlr*conj(gqll)*conj(qdlr);

		} // ll
			
		// update cost
		cost += std::norm ( data - model );

	} // iterate over polar baselines

	/* load complex grad into real and imaginary parts */
	/* chain rule */
	for (int igain = 0; igain < ngains; igain++) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain]   = gg.real();
		rgrad[2*igain+1] = gg.imag();
	}
	rgrad[pindex]      = xpgrad.real(); 
	// xpgrad will be purely real

	return cost;
}
