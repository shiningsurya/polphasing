#include "LBFGS.hpp"

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
		// see :math_gradient.stdout:

		// model forward depends on pb2corr
		complex_type model;

		/*
		We do this over polarbaseline loop, so that we keep track of all the baselines
		*/
		if ( pb2corr == 0 ) {

			const complex_type Dgprr_coeff_gprr = gqll*mrl*conj(gqll)*conj(mrl) + gqll*mrr*conj(gqrl)*conj(mrl) + gqlr*mrl*conj(gqlr)*conj(mrl) + gqlr*mrr*conj(gqrr)*conj(mrl) + gqrl*mrl*conj(gqll)*conj(mrr) + gqrl*mrr*conj(gqrl)*conj(mrr) + gqrr*mrl*conj(gqlr)*conj(mrr) + gqrr*mrr*conj(gqrr)*conj(mrr) ; 
			const complex_type Dgprl_coeff_gprr = gqll*mrl*conj(gqll)*conj(mll) + gqll*mrr*conj(gqrl)*conj(mll) + gqlr*mrl*conj(gqlr)*conj(mll) + gqlr*mrr*conj(gqrr)*conj(mll) + gqrl*mrl*conj(gqll)*conj(mlr) + gqrl*mrr*conj(gqrl)*conj(mlr) + gqrr*mrl*conj(gqlr)*conj(mlr) + gqrr*mrr*conj(gqrr)*conj(mlr) ; 

			const complex_type Dgqrr_coeff_gqrr = gpll*mlr*conj(gpll)*conj(mlr) + gpll*mlr*conj(gplr)*conj(mrr) + gplr*mrr*conj(gpll)*conj(mlr) + gplr*mrr*conj(gplr)*conj(mrr) + gprl*mlr*conj(gprl)*conj(mlr) + gprl*mlr*conj(gprr)*conj(mrr) + gprr*mrr*conj(gprl)*conj(mlr) + gprr*mrr*conj(gprr)*conj(mrr) ; 
			const complex_type Dgqlr_coeff_gqrr = gpll*mll*conj(gpll)*conj(mlr) + gpll*mll*conj(gplr)*conj(mrr) + gplr*mrl*conj(gpll)*conj(mlr) + gplr*mrl*conj(gplr)*conj(mrr) + gprl*mll*conj(gprl)*conj(mlr) + gprl*mll*conj(gprr)*conj(mrr) + gprr*mrl*conj(gprl)*conj(mlr) + gprr*mrl*conj(gprr)*conj(mrr) ; 

			const complex_type Dgprr_coeff_dpqrr = -gqlr*conj(mrl) - gqrr*conj(mrr) ; 
			const complex_type Dgprl_coeff_dpqrr = -gqlr*conj(mll) - gqrr*conj(mlr) ; 

			const complex_type Dgqrr_coeff_dqprr = -gprl*mlr - gprr*mrr ; 
			const complex_type Dgqlr_coeff_dqprr = -gprl*mll - gprr*mrl ; 

			grad[iprr] += Dgprr_coeff_dpqrr*data + Dgprr_coeff_gprr*gprr;
			grad[iprl] += Dgprl_coeff_dpqrr*data + Dgprl_coeff_gprr*gprr;

			grad[iqrr] += Dgqrr_coeff_dqprr*cata + Dgqrr_coeff_gqrr*gqrr;
			grad[iqlr] += Dgqlr_coeff_dqprr*cata + Dgqlr_coeff_gqrr*gqrr;

			// gpra   mab   conj(gqbr)
			model = 
				(gprr * mrr * conj(gqrr)) +
				(gprr * mrl * conj(gqlr)) +
				(gprl * mlr * conj(gqrr)) +
				(gprl * mll * conj(gqlr)) ;

		} // rr
		else if ( pb2corr == 1 ) {

			const complex_type Dgprr_coeff_gprl = gqll*mll*conj(gqll)*conj(mrl) + gqll*mlr*conj(gqrl)*conj(mrl) + gqlr*mll*conj(gqlr)*conj(mrl) + gqlr*mlr*conj(gqrr)*conj(mrl) + gqrl*mll*conj(gqll)*conj(mrr) + gqrl*mlr*conj(gqrl)*conj(mrr) + gqrr*mll*conj(gqlr)*conj(mrr) + gqrr*mlr*conj(gqrr)*conj(mrr) ; 
			const complex_type Dgprl_coeff_gprl = gqll*mll*conj(gqll)*conj(mll) + gqll*mlr*conj(gqrl)*conj(mll) + gqlr*mll*conj(gqlr)*conj(mll) + gqlr*mlr*conj(gqrr)*conj(mll) + gqrl*mll*conj(gqll)*conj(mlr) + gqrl*mlr*conj(gqrl)*conj(mlr) + gqrr*mll*conj(gqlr)*conj(mlr) + gqrr*mlr*conj(gqrr)*conj(mlr) ; 

			const complex_type Dgprr_coeff_dpqrl = -gqll*conj(mrl) - gqrl*conj(mrr) ; 
			const complex_type Dgprl_coeff_dpqrl = -gqll*conj(mll) - gqrl*conj(mlr) ; 

			const complex_type Dgqrl_coeff_gqrl = gpll*mlr*conj(gpll)*conj(mlr) + gpll*mlr*conj(gplr)*conj(mrr) + gplr*mrr*conj(gpll)*conj(mlr) + gplr*mrr*conj(gplr)*conj(mrr) + gprl*mlr*conj(gprl)*conj(mlr) + gprl*mlr*conj(gprr)*conj(mrr) + gprr*mrr*conj(gprl)*conj(mlr) + gprr*mrr*conj(gprr)*conj(mrr) ; 
			const complex_type Dgqll_coeff_gqrl = gpll*mll*conj(gpll)*conj(mlr) + gpll*mll*conj(gplr)*conj(mrr) + gplr*mrl*conj(gpll)*conj(mlr) + gplr*mrl*conj(gplr)*conj(mrr) + gprl*mll*conj(gprl)*conj(mlr) + gprl*mll*conj(gprr)*conj(mrr) + gprr*mrl*conj(gprl)*conj(mlr) + gprr*mrl*conj(gprr)*conj(mrr) ; 

			const complex_type Dgqrl_coeff_dqprl = -gprl*mlr - gprr*mrr ; 
			const complex_type Dgqll_coeff_dqprl = -gprl*mll - gprr*mrl ; 

			grad[iprr] += Dgprr_coeff_dpqrl*data + Dgprr_coeff_gprl*gprl;
			grad[iprl] += Dgprl_coeff_dpqrl*data + Dgprl_coeff_gprl*gprl;

			grad[iqrl] += Dgqrl_coeff_dqprl*cata + Dgqrl_coeff_gqrl*gqrl;
			grad[iqll] += Dgqll_coeff_dqprl*cata + Dgqll_coeff_gqrl*gqrl;

			// gpra   mab   conj(gqbl)
			model = 
				(gprr * mrr * conj(gqrl)) +
				(gprr * mrl * conj(gqll)) +
				(gprl * mlr * conj(gqrl)) +
				(gprl * mll * conj(gqll)) ;

		} // rl
		else if ( pb2corr == 2 ) {

			const complex_type Dgplr_coeff_gplr = gqll*mrl*conj(gqll)*conj(mrl) + gqll*mrr*conj(gqrl)*conj(mrl) + gqlr*mrl*conj(gqlr)*conj(mrl) + gqlr*mrr*conj(gqrr)*conj(mrl) + gqrl*mrl*conj(gqll)*conj(mrr) + gqrl*mrr*conj(gqrl)*conj(mrr) + gqrr*mrl*conj(gqlr)*conj(mrr) + gqrr*mrr*conj(gqrr)*conj(mrr) ; 
			const complex_type Dgpll_coeff_gplr = gqll*mrl*conj(gqll)*conj(mll) + gqll*mrr*conj(gqrl)*conj(mll) + gqlr*mrl*conj(gqlr)*conj(mll) + gqlr*mrr*conj(gqrr)*conj(mll) + gqrl*mrl*conj(gqll)*conj(mlr) + gqrl*mrr*conj(gqrl)*conj(mlr) + gqrr*mrl*conj(gqlr)*conj(mlr) + gqrr*mrr*conj(gqrr)*conj(mlr) ; 

			const complex_type Dgplr_coeff_dpqlr = -gqlr*conj(mrl) - gqrr*conj(mrr) ; 
			const complex_type Dgpll_coeff_dpqlr = -gqlr*conj(mll) - gqrr*conj(mlr) ; 

			const complex_type Dgqrr_coeff_gqlr = gpll*mlr*conj(gpll)*conj(mll) + gpll*mlr*conj(gplr)*conj(mrl) + gplr*mrr*conj(gpll)*conj(mll) + gplr*mrr*conj(gplr)*conj(mrl) + gprl*mlr*conj(gprl)*conj(mll) + gprl*mlr*conj(gprr)*conj(mrl) + gprr*mrr*conj(gprl)*conj(mll) + gprr*mrr*conj(gprr)*conj(mrl) ; 
			const complex_type Dgqlr_coeff_gqlr = gpll*mll*conj(gpll)*conj(mll) + gpll*mll*conj(gplr)*conj(mrl) + gplr*mrl*conj(gpll)*conj(mll) + gplr*mrl*conj(gplr)*conj(mrl) + gprl*mll*conj(gprl)*conj(mll) + gprl*mll*conj(gprr)*conj(mrl) + gprr*mrl*conj(gprl)*conj(mll) + gprr*mrl*conj(gprr)*conj(mrl) ; 

			const complex_type Dgqrr_coeff_dqplr = -gpll*mlr - gplr*mrr ; 
			const complex_type Dgqlr_coeff_dqplr = -gpll*mll - gplr*mrl ; 

			grad[iplr] += Dgplr_coeff_dpqlr*data + Dgplr_coeff_gplr*gplr;
			grad[ipll] += Dgpll_coeff_dpqlr*data + Dgpll_coeff_gplr*gplr;

			grad[iqlr] += Dgqlr_coeff_dqplr*cata + Dgqlr_coeff_gqlr*gqlr;
			grad[iqrr] += Dgqrr_coeff_dqplr*cata + Dgqrr_coeff_gqlr*gqlr;

			// gpla   mab   conj(gqbr)
			model = 
				(gplr * mrr * conj(gqrr)) +
				(gplr * mrl * conj(gqlr)) +
				(gpll * mlr * conj(gqrr)) +
				(gpll * mll * conj(gqlr)) ;

		} // lr
		else if ( pb2corr == 3 ) {

			const complex_type Dgplr_coeff_gpll = gqll*mll*conj(gqll)*conj(mrl) + gqll*mlr*conj(gqrl)*conj(mrl) + gqlr*mll*conj(gqlr)*conj(mrl) + gqlr*mlr*conj(gqrr)*conj(mrl) + gqrl*mll*conj(gqll)*conj(mrr) + gqrl*mlr*conj(gqrl)*conj(mrr) + gqrr*mll*conj(gqlr)*conj(mrr) + gqrr*mlr*conj(gqrr)*conj(mrr) ; 
			const complex_type Dgpll_coeff_gpll = gqll*mll*conj(gqll)*conj(mll) + gqll*mlr*conj(gqrl)*conj(mll) + gqlr*mll*conj(gqlr)*conj(mll) + gqlr*mlr*conj(gqrr)*conj(mll) + gqrl*mll*conj(gqll)*conj(mlr) + gqrl*mlr*conj(gqrl)*conj(mlr) + gqrr*mll*conj(gqlr)*conj(mlr) + gqrr*mlr*conj(gqrr)*conj(mlr) ; 

			const complex_type Dgplr_coeff_dpqll = -gqll*conj(mrl) - gqrl*conj(mrr) ; 
			const complex_type Dgpll_coeff_dpqll = -gqll*conj(mll) - gqrl*conj(mlr) ; 

			const complex_type Dgqrl_coeff_gqll = gpll*mlr*conj(gpll)*conj(mll) + gpll*mlr*conj(gplr)*conj(mrl) + gplr*mrr*conj(gpll)*conj(mll) + gplr*mrr*conj(gplr)*conj(mrl) + gprl*mlr*conj(gprl)*conj(mll) + gprl*mlr*conj(gprr)*conj(mrl) + gprr*mrr*conj(gprl)*conj(mll) + gprr*mrr*conj(gprr)*conj(mrl) ; 
			const complex_type Dgqll_coeff_gqll = gpll*mll*conj(gpll)*conj(mll) + gpll*mll*conj(gplr)*conj(mrl) + gplr*mrl*conj(gpll)*conj(mll) + gplr*mrl*conj(gplr)*conj(mrl) + gprl*mll*conj(gprl)*conj(mll) + gprl*mll*conj(gprr)*conj(mrl) + gprr*mrl*conj(gprl)*conj(mll) + gprr*mrl*conj(gprr)*conj(mrl) ; 

			const complex_type Dgqrl_coeff_dqpll = -gpll*mlr - gplr*mrr ; 
			const complex_type Dgqll_coeff_dqpll = -gpll*mll - gplr*mrl ; 

			grad[iplr] += Dgplr_coeff_dpqll*data + Dgplr_coeff_gpll*gpll;
			grad[ipll] += Dgpll_coeff_dpqll*data + Dgpll_coeff_gpll*gpll;

			grad[iqrl] += Dgqrl_coeff_dqpll*cata + Dgqrl_coeff_gqll*gqll;
			grad[iqll] += Dgqll_coeff_dqpll*cata + Dgqll_coeff_gqll*gqll;

			// gpla   mab   conj(gqbl)
			model = 
				(gplr * mrr * conj(gqrl)) +
				(gplr * mrl * conj(gqll)) +
				(gpll * mlr * conj(gqrl)) +
				(gpll * mll * conj(gqll)) ;

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

LBFGS::real_type LBFGS::diag_jones (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const data_t *pkg = reinterpret_cast<const data_t*>(instance);

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    ngains ( 2 * pkg->nantennas );

	/* populate complex gains vector */
	vc_type  gains ( ngains, complex_type(0.0f, 0.0f) );
	for ( int igain = 0; igain < ngains; igain++ ) {
		gains[igain]   = complex_type ( rgains[2*igain+0], rgains[2*igain+1] );
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
		// see :math_gradient_parallel.stdout:
		// see :math_gradient_parallel.pdf:

		/*
		 * We follow the strategy in :fullgdsolver::gradient: impl. 
		 * We iterate through polarbaselines and update correspondingly.
		 *
		 */

		// model forward depends on pb2corr
		complex_type model;

		// do on every pb2corr
		if ( pb2corr == 0 ) {

			const complex_type Dgprr_coeff_dpqrr = -gqrr*conj(mrr) ; 
			const complex_type Dgqrr_coeff_dqprr = -gprr*mrr ; 

			const complex_type Dgprr_coeff_gprr = gqrr*mrr*conj(gqrr)*conj(mrr) ; 
			const complex_type Dgqrr_coeff_gqrr = gprr*mrr*conj(gprr)*conj(mrr) ; 

			grad[iprr] += Dgprr_coeff_gprr*gprr + Dgprr_coeff_dpqrr*data;
			grad[iqrr] += Dgqrr_coeff_gqrr*gqrr + Dgqrr_coeff_dqprr*conj(data);

			model = gprr * mrr * conj(gqrr);

		} // rr
		else if ( pb2corr == 1 ) {

			const complex_type Dgprr_coeff_dpqrl = -gqll*conj(mrl) ; 
			const complex_type Dgqll_coeff_dqprl = -gprr*mrl ; 

			const complex_type Dgprr_coeff_gprr = gqll*mrl*conj(gqll)*conj(mrl) ; 
			const complex_type Dgqll_coeff_gqll = gprr*mrl*conj(gprr)*conj(mrl) ; 

			grad[iprr] += Dgprr_coeff_gprr*gprr + Dgprr_coeff_dpqrl*data; 
			grad[iqll] += Dgqll_coeff_gqll*gqll + Dgqll_coeff_dqprl*conj(data);

			model = gprr * mrl * conj(gqll);

		} // rl
		else if ( pb2corr == 2 ) {

			const complex_type Dgpll_coeff_dpqlr = -gqrr*conj(mlr) ; 
			const complex_type Dgqrr_coeff_dqplr = -gpll*mlr ; 

			const complex_type Dgpll_coeff_gpll = gqrr*mlr*conj(gqrr)*conj(mlr) ; 
			const complex_type Dgqrr_coeff_gqrr = gpll*mlr*conj(gpll)*conj(mlr) ; 

			grad[ipll] +=  Dgpll_coeff_gpll*gpll + Dgpll_coeff_dpqlr*data;
			grad[iqrr] +=  Dgqrr_coeff_gqrr*gqrr + Dgqrr_coeff_dqplr*conj(data);

			model = gpll * mlr * conj(gqrr);

		} // lr
		else if ( pb2corr == 3 ) {

			const complex_type Dgpll_coeff_dpqll = -gqll*conj(mll) ; 
			const complex_type Dgqll_coeff_dqpll = -gpll*mll ; 

			const complex_type Dgpll_coeff_gpll = gqll*mll*conj(gqll)*conj(mll) ; 
			const complex_type Dgqll_coeff_gqll = gpll*mll*conj(gpll)*conj(mll) ; 

			grad[ipll] += Dgpll_coeff_gpll*gpll + Dgpll_coeff_dpqll*data;
			grad[iqll] += Dgqll_coeff_gqll*gqll + Dgqll_coeff_dqpll*conj(data);

			model = gpll * mll * conj(gqll);

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

int LBFGS::progress_reporter (void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls ) {

	/* get data_t* ptr out of instance */
	data_t *pkg = reinterpret_cast<data_t*>(instance);

	/* save cost and gnorm */
	pkg->cost   = fx;
	pkg->gnorm  = gnorm;
	pkg->niter++;

	/* return 0 always */
	return 0;
}

// for now this solves full jones
LBFGS::real_type LBFGS::Solver::operator() (data_t& pkg) {

	real_type final_cost (0.0f);

	void* vpkg  = static_cast<void*>(&pkg);

  rcode = lbfgs(npar, xpar, &final_cost, full_jones, progress_reporter, vpkg, &param);

  cost   = final_cost;
  niter  = pkg.niter;
  gnorm  = pkg.gnorm;

  return final_cost;
}
