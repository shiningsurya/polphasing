#include "LBFGSpolarizediterative.hpp"

//#define DPRINT

#ifdef DPRINT
#include <iostream>
#endif

LBFGS::real_type LBFGS::diagpolarized_rr (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * Solving for rr gain, populates the derivative of rr that appears in {rr,rl,lr} baselines.
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const data_t *pkg = reinterpret_cast<const data_t*>(instance);

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

		// set the gain indices
		const int iprr ( iant1 );
		const int iqrr ( iant2 );

		// fetch full gains for both antennas
		const complex_type gprr ( cgains[iprr] );
		const complex_type gpll ( pkg->cgains[2*iant1 + 1] );
		const complex_type gprl ( pkg->leakrl[iant1] * gprr );
		const complex_type gplr ( pkg->leaklr[iant1] * gpll );

		const complex_type gqrr ( cgains[iqrr] );
		const complex_type gqll ( pkg->cgains[2*iant2 + 1] );
		const complex_type gqrl ( pkg->leakrl[iant2] * gqrr );
		const complex_type gqlr ( pkg->leaklr[iant2] * gqll );

		// the following long expressions come from sympy
		// see :math_gradient_full_polarized.py:
		// see :math_gradient_full_polarized.code:

		// model forward depends on pb2corr
		complex_type model;

		/* model */
		const complex_type mrr ( pkg->par_model_rr[ibl] );
		const complex_type mrl ( pkg->par_model_rl[ibl] );
		const complex_type mlr ( pkg->par_model_lr[ibl] );
		const complex_type mll ( pkg->par_model_ll[ibl] );

		/*
		We do this over polarbaseline loop, so that we keep track of all the baselines
		*/
		if ( pb2corr == 0 ) {

const complex_type Dgprr_coeff_dpqrr = -gqrl*conj(mrl) - gqrr*conj(mrr) ; 
const complex_type Dgprr_coeff_gprr = gqll*mrl*conj(gqll)*conj(mrl) + gqll*mrr*conj(gqlr)*conj(mrl) + gqlr*mrl*conj(gqll)*conj(mrr) + gqlr*mrr*conj(gqlr)*conj(mrr) + gqrl*mrl*conj(gqrl)*conj(mrl) + gqrl*mrr*conj(gqrr)*conj(mrl) + gqrr*mrl*conj(gqrl)*conj(mrr) + gqrr*mrr*conj(gqrr)*conj(mrr) ; 

const complex_type Dgqrr_coeff_dqprr = -gprl*mlr - gprr*mrr ; 
const complex_type Dgqrr_coeff_gqrr = gpll*mlr*conj(gpll)*conj(mlr) + gpll*mlr*conj(gplr)*conj(mrr) + gplr*mrr*conj(gpll)*conj(mlr) + gplr*mrr*conj(gplr)*conj(mrr) + gprl*mlr*conj(gprl)*conj(mlr) + gprl*mlr*conj(gprr)*conj(mrr) + gprr*mrr*conj(gprl)*conj(mlr) + gprr*mrr*conj(gprr)*conj(mrr) ; 

			grad[iprr] += Dgprr_coeff_dpqrr*data + Dgprr_coeff_gprr*gprr;

			grad[iqrr] += Dgqrr_coeff_dqprr*cata + Dgqrr_coeff_gqrr*gqrr;

model = gprl*mll*conj(gqrl) + gprl*mlr*conj(gqrr) + gprr*mrl*conj(gqrl) + gprr*mrr*conj(gqrr) ;

		} // rr
		else if ( pb2corr == 1 ) {

const complex_type Dgprr_coeff_dpqrl = -gqll*conj(mrl) - gqlr*conj(mrr) ; 
const complex_type Dgprr_coeff_gprl = gqll*mll*conj(gqll)*conj(mrl) + gqll*mlr*conj(gqlr)*conj(mrl) + gqlr*mll*conj(gqll)*conj(mrr) + gqlr*mlr*conj(gqlr)*conj(mrr) + gqrl*mll*conj(gqrl)*conj(mrl) + gqrl*mlr*conj(gqrr)*conj(mrl) + gqrr*mll*conj(gqrl)*conj(mrr) + gqrr*mlr*conj(gqrr)*conj(mrr) ; 

const complex_type Dgqrr_coeff_gqrl = gpll*mlr*conj(gpll)*conj(mll) + gpll*mlr*conj(gplr)*conj(mrl) + gplr*mrr*conj(gpll)*conj(mll) + gplr*mrr*conj(gplr)*conj(mrl) + gprl*mlr*conj(gprl)*conj(mll) + gprl*mlr*conj(gprr)*conj(mrl) + gprr*mrr*conj(gprl)*conj(mll) + gprr*mrr*conj(gprr)*conj(mrl) ; 

			grad[iprr] += Dgprr_coeff_dpqrl*data + Dgprr_coeff_gprl*gprl;

			grad[iqrr] += Dgqrr_coeff_gqrl*gqrl;

model = gprl*mll*conj(gqll) + gprl*mlr*conj(gqlr) + gprr*mrl*conj(gqll) + gprr*mrr*conj(gqlr) ;

		} // rl
		else if ( pb2corr == 2 ) {

const complex_type Dgqrr_coeff_dqplr = -gpll*mlr - gplr*mrr ; 

			grad[iqrr] += Dgqrr_coeff_dqplr*cata;

model = gpll*mll*conj(gqrl) + gpll*mlr*conj(gqrr) + gplr*mrl*conj(gqrl) + gplr*mrr*conj(gqrr) ;

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

LBFGS::real_type LBFGS::diagpolarized_ll (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * Solving for ll gain, populates the derivative of rr that appears in {rr,rl,lr} baselines.
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const data_t *pkg = reinterpret_cast<const data_t*>(instance);

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

		// set the gain indices
		const int ipll ( iant1 );
		const int iqll ( iant2 );

		// fetch full gains for both antennas
		const complex_type gprr ( pkg->cgains[2*iant1 + 0] );
		const complex_type gpll ( cgains[ipll] );
		const complex_type gprl ( pkg->leakrl[iant1] * gprr );
		const complex_type gplr ( pkg->leaklr[iant1] * gpll );

		const complex_type gqrr ( pkg->cgains[2*iant2 + 0] );
		const complex_type gqll ( cgains[iqll] );
		const complex_type gqrl ( pkg->leakrl[iant2] * gqrr );
		const complex_type gqlr ( pkg->leaklr[iant2] * gqll );

		// the following long expressions come from sympy
		// see :math_gradient_full_polarized.py:
		// see :math_gradient_full_polarized.code:

		// model forward depends on pb2corr
		complex_type model;

		/* model */
		const complex_type mrr ( pkg->par_model_rr[ibl] );
		const complex_type mrl ( pkg->par_model_rl[ibl] );
		const complex_type mlr ( pkg->par_model_lr[ibl] );
		const complex_type mll ( pkg->par_model_ll[ibl] );

		/*
		We do this over polarbaseline loop, so that we keep track of all the baselines
		*/
		if ( pb2corr == 1 ) {

const complex_type Dgqll_coeff_dqprl = -gprl*mll - gprr*mrl ; 

			grad[iqll] += Dgqll_coeff_dqprl*cata;

model = gprl*mll*conj(gqll) + gprl*mlr*conj(gqlr) + gprr*mrl*conj(gqll) + gprr*mrr*conj(gqlr) ;

		} // rl
		else if ( pb2corr == 2 ) {

const complex_type Dgpll_coeff_dpqlr = -gqrl*conj(mll) - gqrr*conj(mlr) ; 
const complex_type Dgpll_coeff_gplr = gqll*mrl*conj(gqll)*conj(mll) + gqll*mrr*conj(gqlr)*conj(mll) + gqlr*mrl*conj(gqll)*conj(mlr) + gqlr*mrr*conj(gqlr)*conj(mlr) + gqrl*mrl*conj(gqrl)*conj(mll) + gqrl*mrr*conj(gqrr)*conj(mll) + gqrr*mrl*conj(gqrl)*conj(mlr) + gqrr*mrr*conj(gqrr)*conj(mlr) ; 
const complex_type Dgqll_coeff_gqlr = gpll*mll*conj(gpll)*conj(mlr) + gpll*mll*conj(gplr)*conj(mrr) + gplr*mrl*conj(gpll)*conj(mlr) + gplr*mrl*conj(gplr)*conj(mrr) + gprl*mll*conj(gprl)*conj(mlr) + gprl*mll*conj(gprr)*conj(mrr) + gprr*mrl*conj(gprl)*conj(mlr) + gprr*mrl*conj(gprr)*conj(mrr) ; 
			
			grad[ipll] += Dgpll_coeff_dpqlr*data + Dgpll_coeff_gplr*gplr;

			grad[iqll] += Dgqll_coeff_gqlr*gqlr;

model = gpll*mll*conj(gqrl) + gpll*mlr*conj(gqrr) + gplr*mrl*conj(gqrl) + gplr*mrr*conj(gqrr) ;

		} // lr
		else if ( pb2corr == 3 ) {

const complex_type Dgpll_coeff_dpqll = -gqll*conj(mll) - gqlr*conj(mlr) ; 
const complex_type Dgpll_coeff_gpll = gqll*mll*conj(gqll)*conj(mll) + gqll*mlr*conj(gqlr)*conj(mll) + gqlr*mll*conj(gqll)*conj(mlr) + gqlr*mlr*conj(gqlr)*conj(mlr) + gqrl*mll*conj(gqrl)*conj(mll) + gqrl*mlr*conj(gqrr)*conj(mll) + gqrr*mll*conj(gqrl)*conj(mlr) + gqrr*mlr*conj(gqrr)*conj(mlr) ; 

const complex_type Dgqll_coeff_dqpll = -gpll*mll - gplr*mrl ; 
const complex_type Dgqll_coeff_gqll = gpll*mll*conj(gpll)*conj(mll) + gpll*mll*conj(gplr)*conj(mrl) + gplr*mrl*conj(gpll)*conj(mll) + gplr*mrl*conj(gplr)*conj(mrl) + gprl*mll*conj(gprl)*conj(mll) + gprl*mll*conj(gprr)*conj(mrl) + gprr*mrl*conj(gprl)*conj(mll) + gprr*mrl*conj(gprr)*conj(mrl) ; 

			grad[ipll] += Dgpll_coeff_dpqll*data + Dgpll_coeff_gpll*gpll;
			
			grad[iqll] += Dgqll_coeff_dqpll*cata + Dgqll_coeff_gqll*gqll;

model = gpll*mll*conj(gqll) + gpll*mlr*conj(gqlr) + gplr*mrl*conj(gqll) + gplr*mrr*conj(gqlr) ;

		} // ll
			
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

LBFGS::diagcost_t LBFGS::IterativePolarizedDiagonalSolver::solve (data_t& pkg, int n ) {

	diagcost_t    final_cost;

	void* vpkg  = static_cast<void*>(&pkg);

	for (int __i = 0; __i < n; __i++) {
		/* load gains from pkg into xpars */
		loadgains ( pkg );

		/* solve for rr and ll */
		rcode_rr = lbfgs(n_rr, xpar_rr, &cost_rr, diagpolarized_rr, NULL, vpkg, &param);
		rcode_ll = lbfgs(n_ll, xpar_ll, &cost_ll, diagpolarized_ll, NULL, vpkg, &param);

		/* save solved gains into pkg.cgains */
		savegains ( pkg );

		final_cost[0]  = cost_rr;
		final_cost[1]  = cost_ll;

#ifdef DPRINT
	std::cout << "iter=" << __i << " costs rr=" << cost_rr << " ll=" << cost_ll << std::endl;
#endif
		
	} // solving iterations

  return final_cost;
}
