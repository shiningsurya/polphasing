#include "gdsolver.hpp"

GDSolver::real_type GDSolver::norm ( const vc_type& g ) {
	real_type rnorm ( 0.0f );

	for ( const auto& ig : g ) rnorm += std::norm(ig);

	return rnorm;
}

GDSolver::real_type GDSolver::cost ( const data_t& pkg, const vc_type& gains ) {
	real_type cost ( 0.0f );

	for (int ibl = 0; ibl < npolarbaselines; ibl++) {

		// fetch the antenna index
		const int iant1 ( pkg.iant1[ibl] );
		const int iant2 ( pkg.iant2[ibl] );

		// fetch the pb2corr
		const int pb2corr    = pkg.pb2corr [ ibl ];

		// fetch complex data
		const complex_type data ( pkg.data[ibl] );

		// fetch the par corrected model
		const complex_type mrr ( pkg.par_model_rr[ibl] );
		const complex_type mrl ( pkg.par_model_rl[ibl] );
		const complex_type mlr ( pkg.par_model_lr[ibl] );
		const complex_type mll ( pkg.par_model_ll[ibl] );

		// fetch full gains for both antennas
		const complex_type gprr ( gains[2*iant1 + 0] );
		const complex_type gpll ( gains[2*iant1 + 1] );

		const complex_type gqrr ( gains[2*iant2 + 0] );
		const complex_type gqll ( gains[2*iant2 + 1] );

		// model forward depends on pb2corr
		complex_type model;

		if      ( pb2corr == 0 ) model = gprr * mrr * conj(gqrr);
		else if ( pb2corr == 1 ) model = gprr * mrl * conj(gqll);
		else if ( pb2corr == 2 ) model = gpll * mlr * conj(gqrr);
		else if ( pb2corr == 3 ) model = gpll * mll * conj(gqll);

		// update cost
		cost += std::norm ( data - model );

	} // for every polarbaseline

	return cost;
}

GDSolver::real_type GDSolver::solve ( const data_t& pkg, vc_type& gains ) {
	rcode = 0;
	niter = 0;

	real_type rcost (0.0f);

	// EMA of square of norm of gradient
	real_type ema_gnorm ( 0.0f );

	// EMAs of cost 
	real_type ema_cost_fast ( 0.0f );
	real_type ema_cost_slow ( 0.0f );

	Adam     apple ( ngains, 0.05f, 0.90f, 0.99f, 1000 );
	vc_type  grad ( ngains, complex_type(0.0f, 0.0f) );

	for ( int iter = 0; iter < max_iterations; iter++ ) {

		// find cost before iteration
		const real_type old_cost = cost ( pkg, gains );

		// zero out before solving
		std::fill ( grad.begin(), grad.end(), complex_type(0.0f, 0.0f) );

		// find gradient
		gradient ( pkg, gains, grad );

		// gradient norm
		const real_type gnorm = norm ( grad );

		// EMA of gnorm
		ema_gnorm  = betag*ema_gnorm + (1.0f - betag)*gnorm;

		// use ADAM to update gains
		// in place updation
		apple ( grad, gains );

		// find cost after iteration
		const real_type new_cost = cost ( pkg, gains );

		// EMAs of new cost
		ema_cost_fast = beta_cost_fast*ema_cost_fast + (1.0f - beta_cost_fast)*new_cost;
		ema_cost_slow = beta_cost_slow*ema_cost_slow + (1.0f - beta_cost_slow)*new_cost;

#ifdef CHANDEBUG
		std::cout << iter << " " << new_cost << " " << gnorm << " " << ema_gnorm << " " << ema_cost_slow << " " << ema_cost_fast << std::endl;
#endif 

		/*
		 * We do not have any validation dataset to measure validating error.
		 * We cannot set a threshold on the error as a termination condition, 
		 * because we do not know how the error would be. 
		 *
		 * Instead, we put a termination condition on the norm of the gradient.
		 * (precisely, the square of the norm of the gradient).
		 * Because when the gradient vanishes, we know we are the minimum point.
		 *
		 * Instead of directly using the gnorm which is noisy and does not really show the trend,
		 * we use exponential moving average with a suitable beta (betag)
		 * and set the condition as ema(gnorm) < 0.1
		 *
		 * This is a very stringent condition. It would probably be better to relax it.
		 *
		*/

#ifndef CHANDEBUG
		// termination condition
		if ( ema_gnorm <= delta ) {
			rcode  = 1;
			rcost  = new_cost;
			break;
		}
		if ( std::abs ( ema_cost_fast - ema_cost_slow ) <= gamma ) {
			rcode  = 2;
			rcost  = new_cost;
			break;
		}
#endif

		niter++;
	}

	return rcost;
}

int GDSolver::gradient ( const data_t& pkg, const vc_type& gains, vc_type& grad ) {
	/*
	 * This function computes gradient and saves in grad
	 * both should be 2*nantennas
	 */

	/* iterate over the polar baselines */
	for ( int ibl = 0; ibl < npolarbaselines; ibl++ ) {

		// fetch the antenna index
		const int iant1 ( pkg.iant1[ibl] );
		const int iant2 ( pkg.iant2[ibl] );

		// fetch the pb2corr
		const int pb2corr    = pkg.pb2corr [ ibl ];

		// fetch complex data
		const complex_type data ( pkg.data[ibl] );

		// fetch the par corrected model
		const complex_type mrr ( pkg.par_model_rr[ibl] );
		const complex_type mrl ( pkg.par_model_rl[ibl] );
		const complex_type mlr ( pkg.par_model_lr[ibl] );
		const complex_type mll ( pkg.par_model_ll[ibl] );

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

		// do on every pb2corr
		if ( pb2corr == 0 ) {

			const complex_type Dgprr_coeff_dpqrr = -gqrr*conj(mrr) ; 
			const complex_type Dgqrr_coeff_dqprr = -gprr*mrr ; 

			const complex_type Dgprr_coeff_gprr = gqrr*mrr*conj(gqrr)*conj(mrr) ; 
			const complex_type Dgqrr_coeff_gqrr = gprr*mrr*conj(gprr)*conj(mrr) ; 

			grad[iprr] += Dgprr_coeff_gprr*gprr + Dgprr_coeff_dpqrr*data;
			grad[iqrr] += Dgqrr_coeff_gqrr*gqrr + Dgqrr_coeff_dqprr*conj(data);

		} // rr
		else if ( pb2corr == 1 ) {

			const complex_type Dgprr_coeff_dpqrl = -gqll*conj(mrl) ; 
			const complex_type Dgqll_coeff_dqprl = -gprr*mrl ; 

			const complex_type Dgprr_coeff_gprr = gqll*mrl*conj(gqll)*conj(mrl) ; 
			const complex_type Dgqll_coeff_gqll = gprr*mrl*conj(gprr)*conj(mrl) ; 

			grad[iprr] += Dgprr_coeff_gprr*gprr + Dgprr_coeff_dpqrl*data; 
			grad[iqll] += Dgqll_coeff_gqll*gqll + Dgqll_coeff_dqprl*conj(data);

		} // rl
		else if ( pb2corr == 2 ) {

			const complex_type Dgpll_coeff_dpqlr = -gqrr*conj(mlr) ; 
			const complex_type Dgqrr_coeff_dqplr = -gpll*mlr ; 

			const complex_type Dgpll_coeff_gpll = gqrr*mlr*conj(gqrr)*conj(mlr) ; 
			const complex_type Dgqrr_coeff_gqrr = gpll*mlr*conj(gpll)*conj(mlr) ; 

			grad[ipll] +=  Dgpll_coeff_gpll*gpll + Dgpll_coeff_dpqlr*data;
			grad[iqrr] +=  Dgqrr_coeff_gqrr*gqrr + Dgqrr_coeff_dqplr*conj(data);

		} // lr
		else if ( pb2corr == 3 ) {

			const complex_type Dgpll_coeff_dpqll = -gqll*conj(mll) ; 
			const complex_type Dgqll_coeff_dqpll = -gpll*mll ; 

			const complex_type Dgpll_coeff_gpll = gqll*mll*conj(gqll)*conj(mll) ; 
			const complex_type Dgqll_coeff_gqll = gpll*mll*conj(gpll)*conj(mll) ; 

			grad[ipll] += Dgpll_coeff_gpll*gpll + Dgpll_coeff_dpqll*data;
			grad[iqll] += Dgqll_coeff_gqll*gqll + Dgqll_coeff_dqpll*conj(data);

		} // ll

	} // iterate over polar baselines

	return 0;
}
