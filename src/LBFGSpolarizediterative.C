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
	for ( int ibl = 0; ibl < pkg->nbaselines; ibl++ ) {

		// fetch the antenna index
		const int iant1 ( pkg->iant1[ibl] );
		const int iant2 ( pkg->iant2[ibl] );

		// polarbaseline index
		const int pbl ( 4*ibl );

		// fetch complex data
		const complex_type dpqrr ( pkg->data[ pbl    ] );
		const complex_type dpqrl ( pkg->data[ pbl + 1] );
		const complex_type dpqlr ( pkg->data[ pbl + 2] );
		const complex_type dpqll ( pkg->data[ pbl + 3] );

		const complex_type cdpqrr ( conj ( dpqrr ) );
		const complex_type cdpqrl ( conj ( dpqrl ) );
		const complex_type cdpqlr ( conj ( dpqlr ) );
		const complex_type cdpqll ( conj ( dpqll ) );

		// set the gain indices
		const int iprr ( iant1 );
		const int iqrr ( iant2 );

		// fetch full gains for both antennas
		const complex_type gprr ( cgains[iprr] );
		const complex_type gpll ( pkg->cgains[2*iant1 + 1] );
		const complex_type lprl ( pkg->leakrl[iant1] );
		const complex_type lplr ( pkg->leaklr[iant1] );

		const complex_type cgprr ( conj ( gprr ) );
		const complex_type cgpll ( conj ( gpll ) );
		const complex_type clprl ( conj ( lprl ) );
		const complex_type clplr ( conj ( lplr ) );

		const complex_type gqrr ( cgains[iqrr] );
		const complex_type gqll ( pkg->cgains[2*iant2 + 1] );
		const complex_type lqrl ( pkg->leakrl[iant2] );
		const complex_type lqlr ( pkg->leaklr[iant2] );

		const complex_type cgqrr ( conj ( gqrr ) );
		const complex_type cgqll ( conj ( gqll ) );
		const complex_type clqrl ( conj ( lqrl ) );
		const complex_type clqlr ( conj ( lqlr ) );

		/* model */
		const complex_type mrr ( pkg->par_model_rr[ibl] );
		const complex_type mrl ( pkg->par_model_rl[ibl] );
		const complex_type mlr ( pkg->par_model_lr[ibl] );
		const complex_type mll ( pkg->par_model_ll[ibl] );

		const complex_type cmrr ( conj ( mrr ) );
		const complex_type cmrl ( conj ( mrl ) );
		const complex_type cmlr ( conj ( mlr ) );
		const complex_type cmll ( conj ( mll ) );

		// the following long expressions come from sympy
		// see :math_gradient_parallel_polarized.py:
		// see :math_gradient_parallel_polarized.stdout:
		
const complex_type Dgprr_dpqrr_gqrr = -lqrl*clprl*cmll - lqrl*cmrl - clprl*cmlr - cmrr ; 
const complex_type Dgprr_dpqrl_gqll = -lqlr*clprl*cmlr - lqlr*cmrr - clprl*cmll - cmrl ; 
const complex_type Dgprr_gprr_ngqll = lprl*lqlr*mll*clprl*cmlr + lprl*lqlr*mll*cmrr + lprl*lqlr*mlr*clprl*clqlr*cmlr + lprl*lqlr*mlr*clqlr*cmrr + lprl*mll*clprl*cmll + lprl*mll*cmrl + lprl*mlr*clprl*clqlr*cmll + lprl*mlr*clqlr*cmrl + lqlr*mrl*clprl*cmlr + lqlr*mrl*cmrr + lqlr*mrr*clprl*clqlr*cmlr + lqlr*mrr*clqlr*cmrr + mrl*clprl*cmll + mrl*cmrl + mrr*clprl*clqlr*cmll + mrr*clqlr*cmrl ; 
const complex_type Dgprr_gprr_ngqrr = lprl*lqrl*mll*clprl*clqrl*cmll + lprl*lqrl*mll*clqrl*cmrl + lprl*lqrl*mlr*clprl*cmll + lprl*lqrl*mlr*cmrl + lprl*mll*clprl*clqrl*cmlr + lprl*mll*clqrl*cmrr + lprl*mlr*clprl*cmlr + lprl*mlr*cmrr + lqrl*mrl*clprl*clqrl*cmll + lqrl*mrl*clqrl*cmrl + lqrl*mrr*clprl*cmll + lqrl*mrr*cmrl + mrl*clprl*clqrl*cmlr + mrl*clqrl*cmrr + mrr*clprl*cmlr + mrr*cmrr ; 
		
		grad[iprr] += Dgprr_dpqrr_gqrr*gqrr*dpqrr + Dgprr_dpqrl_gqll*gqll*dpqrl;
		grad[iprr] += Dgprr_gprr_ngqll*gqll*cgqll*gprr + Dgprr_gprr_ngqrr*gqrr*cgqrr*gprr;

const complex_type Dgqrr_cdpqrr_gprr = -lprl*mll*clqrl - lprl*mlr - mrl*clqrl - mrr ; 
const complex_type Dgqrr_cdpqlr_gpll = -lplr*mrl*clqrl - lplr*mrr - mll*clqrl - mlr ; 
const complex_type Dgqrr_gqrr_ngpll = lplr*lqrl*mrl*clplr*clqrl*cmrl + lplr*lqrl*mrl*clqrl*cmll + lplr*lqrl*mrr*clplr*cmrl + lplr*lqrl*mrr*cmll + lplr*mrl*clplr*clqrl*cmrr + lplr*mrl*clqrl*cmlr + lplr*mrr*clplr*cmrr + lplr*mrr*cmlr + lqrl*mll*clplr*clqrl*cmrl + lqrl*mll*clqrl*cmll + lqrl*mlr*clplr*cmrl + lqrl*mlr*cmll + mll*clplr*clqrl*cmrr + mll*clqrl*cmlr + mlr*clplr*cmrr + mlr*cmlr ; 
const complex_type Dgqrr_gqrr_ngprr = lprl*lqrl*mll*clprl*clqrl*cmll + lprl*lqrl*mll*clqrl*cmrl + lprl*lqrl*mlr*clprl*cmll + lprl*lqrl*mlr*cmrl + lprl*mll*clprl*clqrl*cmlr + lprl*mll*clqrl*cmrr + lprl*mlr*clprl*cmlr + lprl*mlr*cmrr + lqrl*mrl*clprl*clqrl*cmll + lqrl*mrl*clqrl*cmrl + lqrl*mrr*clprl*cmll + lqrl*mrr*cmrl + mrl*clprl*clqrl*cmlr + mrl*clqrl*cmrr + mrr*clprl*cmlr + mrr*cmrr ; 
		
		grad[iqrr] += Dgqrr_cdpqrr_gprr*gprr*cdpqrr + Dgqrr_cdpqlr_gpll*gpll*cdpqlr;
		grad[iqrr] += Dgqrr_gqrr_ngpll*gpll*cgpll*gqrr + Dgqrr_gqrr_ngprr*gprr*cgprr*gqrr;
			
		// update cost
		// rr
		cost += std::norm ( dpqrr - gprr*lprl*mll*cgqrr*clqrl - gprr*lprl*mlr*cgqrr - gprr*mrl*cgqrr*clqrl - gprr*mrr*cgqrr );
		// rl
		cost += std::norm ( dpqrl - gprr*lprl*mll*cgqll - gprr*lprl*mlr*cgqll*clqlr - gprr*mrl*cgqll - gprr*mrr*cgqll*clqlr );
		// lr
		cost += std::norm ( dpqlr - gpll*lplr*mrl*cgqrr*clqrl - gpll*lplr*mrr*cgqrr - gpll*mll*cgqrr*clqrl - gpll*mlr*cgqrr );
		// ll
		//cost += std::norm ( dpqll - gpll*lplr*mrl*cgqll - gpll*lplr*mrr*cgqll*clqlr - gpll*mll*cgqll - gpll*mlr*cgqll*clqlr );

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
	 * Solving for ll gain, populates the derivative of ll that appears in {rl,lr,ll} baselines.
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
	for ( int ibl = 0; ibl < pkg->nbaselines; ibl++ ) {

		// fetch the antenna index
		const int iant1 ( pkg->iant1[ibl] );
		const int iant2 ( pkg->iant2[ibl] );

		// polarbaseline index
		const int pbl ( 4*ibl );

		// fetch complex data
		const complex_type dpqrr ( pkg->data[ pbl    ] );
		const complex_type dpqrl ( pkg->data[ pbl + 1] );
		const complex_type dpqlr ( pkg->data[ pbl + 2] );
		const complex_type dpqll ( pkg->data[ pbl + 3] );

		const complex_type cdpqrr ( conj ( dpqrr ) );
		const complex_type cdpqrl ( conj ( dpqrl ) );
		const complex_type cdpqlr ( conj ( dpqlr ) );
		const complex_type cdpqll ( conj ( dpqll ) );

		// set the gain indices
		const int ipll ( iant1 );
		const int iqll ( iant2 );

		// fetch full gains for both antennas
		const complex_type gprr ( pkg->cgains[2*iant1] );
		const complex_type gpll ( cgains[ipll] );
		const complex_type lprl ( pkg->leakrl[iant1] );
		const complex_type lplr ( pkg->leaklr[iant1] );

		const complex_type cgprr ( conj ( gprr ) );
		const complex_type cgpll ( conj ( gpll ) );
		const complex_type clprl ( conj ( lprl ) );
		const complex_type clplr ( conj ( lplr ) );

		const complex_type gqrr ( pkg->cgains[2*iant2] );
		const complex_type gqll ( cgains[iqll] );
		const complex_type lqrl ( pkg->leakrl[iant2] );
		const complex_type lqlr ( pkg->leaklr[iant2] );

		const complex_type cgqrr ( conj ( gqrr ) );
		const complex_type cgqll ( conj ( gqll ) );
		const complex_type clqrl ( conj ( lqrl ) );
		const complex_type clqlr ( conj ( lqlr ) );

		/* model */
		const complex_type mrr ( pkg->par_model_rr[ibl] );
		const complex_type mrl ( pkg->par_model_rl[ibl] );
		const complex_type mlr ( pkg->par_model_lr[ibl] );
		const complex_type mll ( pkg->par_model_ll[ibl] );

		const complex_type cmrr ( conj ( mrr ) );
		const complex_type cmrl ( conj ( mrl ) );
		const complex_type cmlr ( conj ( mlr ) );
		const complex_type cmll ( conj ( mll ) );

		// the following long expressions come from sympy
		// see :math_gradient_parallel_polarized.py:
		// see :math_gradient_parallel_polarized.stdout:
		
const complex_type Dgpll_dpqlr_gqrr = -lqrl*clplr*cmrl - lqrl*cmll - clplr*cmrr - cmlr ; 
const complex_type Dgpll_dpqll_gqll = -lqlr*clplr*cmrr - lqlr*cmlr - clplr*cmrl - cmll ; 
const complex_type Dgpll_gpll_ngqll = lplr*lqlr*mrl*clplr*cmrr + lplr*lqlr*mrl*cmlr + lplr*lqlr*mrr*clplr*clqlr*cmrr + lplr*lqlr*mrr*clqlr*cmlr + lplr*mrl*clplr*cmrl + lplr*mrl*cmll + lplr*mrr*clplr*clqlr*cmrl + lplr*mrr*clqlr*cmll + lqlr*mll*clplr*cmrr + lqlr*mll*cmlr + lqlr*mlr*clplr*clqlr*cmrr + lqlr*mlr*clqlr*cmlr + mll*clplr*cmrl + mll*cmll + mlr*clplr*clqlr*cmrl + mlr*clqlr*cmll ; 
const complex_type Dgpll_gpll_ngqrr = lplr*lqrl*mrl*clplr*clqrl*cmrl + lplr*lqrl*mrl*clqrl*cmll + lplr*lqrl*mrr*clplr*cmrl + lplr*lqrl*mrr*cmll + lplr*mrl*clplr*clqrl*cmrr + lplr*mrl*clqrl*cmlr + lplr*mrr*clplr*cmrr + lplr*mrr*cmlr + lqrl*mll*clplr*clqrl*cmrl + lqrl*mll*clqrl*cmll + lqrl*mlr*clplr*cmrl + lqrl*mlr*cmll + mll*clplr*clqrl*cmrr + mll*clqrl*cmlr + mlr*clplr*cmrr + mlr*cmlr ; 

		grad[ipll] += Dgpll_dpqlr_gqrr*gqrr*dpqlr + Dgpll_dpqll_gqll*gqll*dpqll;
		grad[ipll] += Dgpll_gpll_ngqll*gqll*cgqll*gpll + Dgpll_gpll_ngqrr*gqrr*cgqrr*gpll;

const complex_type Dgqll_cdpqrl_gprr = -lprl*mll - lprl*mlr*clqlr - mrl - mrr*clqlr ; 
const complex_type Dgqll_cdpqll_gpll = -lplr*mrl - lplr*mrr*clqlr - mll - mlr*clqlr ; 
const complex_type Dgqll_gqll_ngpll = lplr*lqlr*mrl*clplr*cmrr + lplr*lqlr*mrl*cmlr + lplr*lqlr*mrr*clplr*clqlr*cmrr + lplr*lqlr*mrr*clqlr*cmlr + lplr*mrl*clplr*cmrl + lplr*mrl*cmll + lplr*mrr*clplr*clqlr*cmrl + lplr*mrr*clqlr*cmll + lqlr*mll*clplr*cmrr + lqlr*mll*cmlr + lqlr*mlr*clplr*clqlr*cmrr + lqlr*mlr*clqlr*cmlr + mll*clplr*cmrl + mll*cmll + mlr*clplr*clqlr*cmrl + mlr*clqlr*cmll ; 
const complex_type Dgqll_gqll_ngprr = lprl*lqlr*mll*clprl*cmlr + lprl*lqlr*mll*cmrr + lprl*lqlr*mlr*clprl*clqlr*cmlr + lprl*lqlr*mlr*clqlr*cmrr + lprl*mll*clprl*cmll + lprl*mll*cmrl + lprl*mlr*clprl*clqlr*cmll + lprl*mlr*clqlr*cmrl + lqlr*mrl*clprl*cmlr + lqlr*mrl*cmrr + lqlr*mrr*clprl*clqlr*cmlr + lqlr*mrr*clqlr*cmrr + mrl*clprl*cmll + mrl*cmrl + mrr*clprl*clqlr*cmll + mrr*clqlr*cmrl ; 

		grad[iqll] += Dgqll_cdpqrl_gprr*gprr*cdpqrl + Dgqll_cdpqll_gpll*gpll*cdpqll;
		grad[iqll] += Dgqll_gqll_ngpll*gpll*cgpll*gqll + Dgqll_gqll_ngprr*gprr*cgprr*gqll;
			
		// update cost
		// rr
		//cost += std::norm ( dpqrr - gprr*lprl*mll*cgqrr*clqrl - gprr*lprl*mlr*cgqrr - gprr*mrl*cgqrr*clqrl - gprr*mrr*cgqrr );
		// rl
		cost += std::norm ( dpqrl - gprr*lprl*mll*cgqll - gprr*lprl*mlr*cgqll*clqlr - gprr*mrl*cgqll - gprr*mrr*cgqll*clqlr );
		// lr
		cost += std::norm ( dpqlr - gpll*lplr*mrl*cgqrr*clqrl - gpll*lplr*mrr*cgqrr - gpll*mll*cgqrr*clqrl - gpll*mlr*cgqrr );
		// ll
		cost += std::norm ( dpqll - gpll*lplr*mrl*cgqll - gpll*lplr*mrr*cgqll*clqlr - gpll*mll*cgqll - gpll*mlr*cgqll*clqlr );

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
