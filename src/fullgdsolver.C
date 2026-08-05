#include "fullgdsolver.hpp"

using std::conj;

/*
 * i need a way to selectively use gains
 * Like, 
 * only use parallel hand gains (rr and ll)
 *
 * i can enforce gains to be zero.
 * But that does not make the corresponding gradient to be zero?
 * It makes it, but our gradient equations are condensed and derived when they are non zero.
*/

int FullGDSolver::gradient ( const solve_data_t& pkg, const vc_type& gains, vc_type& grad ) {
	/*
	 * This function computes gradient and saves in grad
	 * both should be 4*nantennas
	 *
	 * grad should be 4*nantennas;
	 */

	/* zero out gradient */
	// 20260727: the caller should zero out the gradient
	// std::fill ( grad.begin(), grad.end(), complex_type(0.0f, 0.0f) );

	/* iterate over the polar baselines */
	for ( int ibl = 0; ibl < npolarbaselines; ibl++ ) {

		// fetch the antenna index
		const int iant1 ( pkg.iant1[ibl] );
		const int iant2 ( pkg.iant2[ibl] );

		// fetch the pb2corr
		const int pb2corr    = pkg.pb2corr [ ibl ];

		// fetch complex data
		const complex_type data ( pkg.data[ibl] );
		const complex_type cata ( conj(data) );

		// fetch the par corrected model
		const complex_type mrr ( pkg.par_model_rr[ibl] );
		const complex_type mrl ( pkg.par_model_rl[ibl] );
		const complex_type mlr ( pkg.par_model_lr[ibl] );
		const complex_type mll ( pkg.par_model_ll[ibl] );

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
		const complex_type gprr ( gains[iprr] );
		const complex_type gprl ( gains[iprl] );
		const complex_type gplr ( gains[iplr] );
		const complex_type gpll ( gains[ipll] );

		const complex_type gqrr ( gains[iqrr] );
		const complex_type gqrl ( gains[iqrl] );
		const complex_type gqlr ( gains[iqlr] );
		const complex_type gqll ( gains[iqll] );

		// the following long expressions come from sympy
		// see :math_gradient.py:
		// see :math_gradient.stdout:

		/*
		We do this over polarbaseline loop, so that we keep track of all the baselines
		*/

		// do on every pb2corr
		if ( pb2corr == 0 ) {

const complex_type Dgprr_coeff_dpqrr = -gqrl*conj(mrl) - gqrr*conj(mrr) ; 
const complex_type Dgprl_coeff_dpqrr = -gqrl*conj(mll) - gqrr*conj(mlr) ; 
const complex_type Dgqrr_coeff_dqprr = -gprl*mlr - gprr*mrr ; 
const complex_type Dgqrl_coeff_dqprr = -gprl*mll - gprr*mrl ; 
const complex_type Dgprr_coeff_gprr = gqll*mrl*conj(gqll)*conj(mrl) + gqll*mrr*conj(gqlr)*conj(mrl) + gqlr*mrl*conj(gqll)*conj(mrr) + gqlr*mrr*conj(gqlr)*conj(mrr) + gqrl*mrl*conj(gqrl)*conj(mrl) + gqrl*mrr*conj(gqrr)*conj(mrl) + gqrr*mrl*conj(gqrl)*conj(mrr) + gqrr*mrr*conj(gqrr)*conj(mrr) ; 
const complex_type Dgprl_coeff_gprr = gqll*mrl*conj(gqll)*conj(mll) + gqll*mrr*conj(gqlr)*conj(mll) + gqlr*mrl*conj(gqll)*conj(mlr) + gqlr*mrr*conj(gqlr)*conj(mlr) + gqrl*mrl*conj(gqrl)*conj(mll) + gqrl*mrr*conj(gqrr)*conj(mll) + gqrr*mrl*conj(gqrl)*conj(mlr) + gqrr*mrr*conj(gqrr)*conj(mlr) ; 
const complex_type Dgqrr_coeff_gqrr = gpll*mlr*conj(gpll)*conj(mlr) + gpll*mlr*conj(gplr)*conj(mrr) + gplr*mrr*conj(gpll)*conj(mlr) + gplr*mrr*conj(gplr)*conj(mrr) + gprl*mlr*conj(gprl)*conj(mlr) + gprl*mlr*conj(gprr)*conj(mrr) + gprr*mrr*conj(gprl)*conj(mlr) + gprr*mrr*conj(gprr)*conj(mrr) ; 
const complex_type Dgqrl_coeff_gqrr = gpll*mll*conj(gpll)*conj(mlr) + gpll*mll*conj(gplr)*conj(mrr) + gplr*mrl*conj(gpll)*conj(mlr) + gplr*mrl*conj(gplr)*conj(mrr) + gprl*mll*conj(gprl)*conj(mlr) + gprl*mll*conj(gprr)*conj(mrr) + gprr*mrl*conj(gprl)*conj(mlr) + gprr*mrl*conj(gprr)*conj(mrr) ; 

			grad[iprr] += Dgprr_coeff_dpqrr*data + Dgprr_coeff_gprr*gprr;
			grad[iprl] += Dgprl_coeff_dpqrr*data + Dgprl_coeff_gprr*gprr;

			grad[iqrr] += Dgqrr_coeff_dqprr*cata + Dgqrr_coeff_gqrr*gqrr;
			grad[iqrl] += Dgqrl_coeff_dqprr*cata + Dgqrl_coeff_gqrr*gqrr;

//model = gprl*mll*conj(gqrl) + gprl*mlr*conj(gqrr) + gprr*mrl*conj(gqrl) + gprr*mrr*conj(gqrr) ;

		} // rr
		else if ( pb2corr == 1 ) {

const complex_type Dgprr_coeff_dpqrl = -gqll*conj(mrl) - gqlr*conj(mrr) ; 
const complex_type Dgprl_coeff_dpqrl = -gqll*conj(mll) - gqlr*conj(mlr) ; 
const complex_type Dgqlr_coeff_dqprl = -gprl*mlr - gprr*mrr ; 
const complex_type Dgqll_coeff_dqprl = -gprl*mll - gprr*mrl ; 
const complex_type Dgprr_coeff_gprl = gqll*mll*conj(gqll)*conj(mrl) + gqll*mlr*conj(gqlr)*conj(mrl) + gqlr*mll*conj(gqll)*conj(mrr) + gqlr*mlr*conj(gqlr)*conj(mrr) + gqrl*mll*conj(gqrl)*conj(mrl) + gqrl*mlr*conj(gqrr)*conj(mrl) + gqrr*mll*conj(gqrl)*conj(mrr) + gqrr*mlr*conj(gqrr)*conj(mrr) ; 
const complex_type Dgprl_coeff_gprl = gqll*mll*conj(gqll)*conj(mll) + gqll*mlr*conj(gqlr)*conj(mll) + gqlr*mll*conj(gqll)*conj(mlr) + gqlr*mlr*conj(gqlr)*conj(mlr) + gqrl*mll*conj(gqrl)*conj(mll) + gqrl*mlr*conj(gqrr)*conj(mll) + gqrr*mll*conj(gqrl)*conj(mlr) + gqrr*mlr*conj(gqrr)*conj(mlr) ; 
const complex_type Dgqrr_coeff_gqrl = gpll*mlr*conj(gpll)*conj(mll) + gpll*mlr*conj(gplr)*conj(mrl) + gplr*mrr*conj(gpll)*conj(mll) + gplr*mrr*conj(gplr)*conj(mrl) + gprl*mlr*conj(gprl)*conj(mll) + gprl*mlr*conj(gprr)*conj(mrl) + gprr*mrr*conj(gprl)*conj(mll) + gprr*mrr*conj(gprr)*conj(mrl) ; 
const complex_type Dgqrl_coeff_gqrl = gpll*mll*conj(gpll)*conj(mll) + gpll*mll*conj(gplr)*conj(mrl) + gplr*mrl*conj(gpll)*conj(mll) + gplr*mrl*conj(gplr)*conj(mrl) + gprl*mll*conj(gprl)*conj(mll) + gprl*mll*conj(gprr)*conj(mrl) + gprr*mrl*conj(gprl)*conj(mll) + gprr*mrl*conj(gprr)*conj(mrl) ; 

			grad[iprr] += Dgprr_coeff_dpqrl*data + Dgprr_coeff_gprl*gprl;
			grad[iprl] += Dgprl_coeff_dpqrl*data + Dgprl_coeff_gprl*gprl;

			grad[iqlr] += Dgqlr_coeff_dqprl*cata;
			grad[iqll] += Dgqll_coeff_dqprl*cata;

			grad[iqrr] += Dgqrr_coeff_gqrl*gqrl;
			grad[iqrl] += Dgqrl_coeff_gqrl*gqrl;

//model = gprl*mll*conj(gqll) + gprl*mlr*conj(gqlr) + gprr*mrl*conj(gqll) + gprr*mrr*conj(gqlr) ;

		} // rl
		else if ( pb2corr == 2 ) {

const complex_type Dgplr_coeff_dpqlr = -gqrl*conj(mrl) - gqrr*conj(mrr) ; 
const complex_type Dgpll_coeff_dpqlr = -gqrl*conj(mll) - gqrr*conj(mlr) ; 
const complex_type Dgqrr_coeff_dqplr = -gpll*mlr - gplr*mrr ; 
const complex_type Dgqrl_coeff_dqplr = -gpll*mll - gplr*mrl ; 
const complex_type Dgplr_coeff_gplr = gqll*mrl*conj(gqll)*conj(mrl) + gqll*mrr*conj(gqlr)*conj(mrl) + gqlr*mrl*conj(gqll)*conj(mrr) + gqlr*mrr*conj(gqlr)*conj(mrr) + gqrl*mrl*conj(gqrl)*conj(mrl) + gqrl*mrr*conj(gqrr)*conj(mrl) + gqrr*mrl*conj(gqrl)*conj(mrr) + gqrr*mrr*conj(gqrr)*conj(mrr) ; 
const complex_type Dgpll_coeff_gplr = gqll*mrl*conj(gqll)*conj(mll) + gqll*mrr*conj(gqlr)*conj(mll) + gqlr*mrl*conj(gqll)*conj(mlr) + gqlr*mrr*conj(gqlr)*conj(mlr) + gqrl*mrl*conj(gqrl)*conj(mll) + gqrl*mrr*conj(gqrr)*conj(mll) + gqrr*mrl*conj(gqrl)*conj(mlr) + gqrr*mrr*conj(gqrr)*conj(mlr) ; 
const complex_type Dgqlr_coeff_gqlr = gpll*mlr*conj(gpll)*conj(mlr) + gpll*mlr*conj(gplr)*conj(mrr) + gplr*mrr*conj(gpll)*conj(mlr) + gplr*mrr*conj(gplr)*conj(mrr) + gprl*mlr*conj(gprl)*conj(mlr) + gprl*mlr*conj(gprr)*conj(mrr) + gprr*mrr*conj(gprl)*conj(mlr) + gprr*mrr*conj(gprr)*conj(mrr) ; 
const complex_type Dgqll_coeff_gqlr = gpll*mll*conj(gpll)*conj(mlr) + gpll*mll*conj(gplr)*conj(mrr) + gplr*mrl*conj(gpll)*conj(mlr) + gplr*mrl*conj(gplr)*conj(mrr) + gprl*mll*conj(gprl)*conj(mlr) + gprl*mll*conj(gprr)*conj(mrr) + gprr*mrl*conj(gprl)*conj(mlr) + gprr*mrl*conj(gprr)*conj(mrr) ; 

			grad[iplr] += Dgplr_coeff_dpqlr*data + Dgplr_coeff_gplr*gplr;
			grad[ipll] += Dgpll_coeff_dpqlr*data + Dgpll_coeff_gplr*gplr;

			grad[iqrr] += Dgqrr_coeff_dqplr*cata;
			grad[iqrl] += Dgqrl_coeff_dqplr*cata;

			grad[iqlr] += Dgqlr_coeff_gqlr*gqlr;
			grad[iqll] += Dgqll_coeff_gqlr*gqlr;

//model = gpll*mll*conj(gqrl) + gpll*mlr*conj(gqrr) + gplr*mrl*conj(gqrl) + gplr*mrr*conj(gqrr) ;

		} // lr
		else if ( pb2corr == 3 ) {

const complex_type Dgplr_coeff_dpqll = -gqll*conj(mrl) - gqlr*conj(mrr) ; 
const complex_type Dgpll_coeff_dpqll = -gqll*conj(mll) - gqlr*conj(mlr) ; 
const complex_type Dgqlr_coeff_dqpll = -gpll*mlr - gplr*mrr ; 
const complex_type Dgqll_coeff_dqpll = -gpll*mll - gplr*mrl ; 
const complex_type Dgplr_coeff_gpll = gqll*mll*conj(gqll)*conj(mrl) + gqll*mlr*conj(gqlr)*conj(mrl) + gqlr*mll*conj(gqll)*conj(mrr) + gqlr*mlr*conj(gqlr)*conj(mrr) + gqrl*mll*conj(gqrl)*conj(mrl) + gqrl*mlr*conj(gqrr)*conj(mrl) + gqrr*mll*conj(gqrl)*conj(mrr) + gqrr*mlr*conj(gqrr)*conj(mrr) ; 
const complex_type Dgpll_coeff_gpll = gqll*mll*conj(gqll)*conj(mll) + gqll*mlr*conj(gqlr)*conj(mll) + gqlr*mll*conj(gqll)*conj(mlr) + gqlr*mlr*conj(gqlr)*conj(mlr) + gqrl*mll*conj(gqrl)*conj(mll) + gqrl*mlr*conj(gqrr)*conj(mll) + gqrr*mll*conj(gqrl)*conj(mlr) + gqrr*mlr*conj(gqrr)*conj(mlr) ; 
const complex_type Dgqlr_coeff_gqll = gpll*mlr*conj(gpll)*conj(mll) + gpll*mlr*conj(gplr)*conj(mrl) + gplr*mrr*conj(gpll)*conj(mll) + gplr*mrr*conj(gplr)*conj(mrl) + gprl*mlr*conj(gprl)*conj(mll) + gprl*mlr*conj(gprr)*conj(mrl) + gprr*mrr*conj(gprl)*conj(mll) + gprr*mrr*conj(gprr)*conj(mrl) ; 
const complex_type Dgqll_coeff_gqll = gpll*mll*conj(gpll)*conj(mll) + gpll*mll*conj(gplr)*conj(mrl) + gplr*mrl*conj(gpll)*conj(mll) + gplr*mrl*conj(gplr)*conj(mrl) + gprl*mll*conj(gprl)*conj(mll) + gprl*mll*conj(gprr)*conj(mrl) + gprr*mrl*conj(gprl)*conj(mll) + gprr*mrl*conj(gprr)*conj(mrl) ; 

			grad[iplr] += Dgplr_coeff_dpqll*data + Dgplr_coeff_gpll*gpll;
			grad[ipll] += Dgpll_coeff_dpqll*data + Dgpll_coeff_gpll*gpll;

			grad[iqlr] += Dgqlr_coeff_dqpll*cata + Dgqlr_coeff_gqll*gqll;
			grad[iqll] += Dgqll_coeff_dqpll*cata + Dgqll_coeff_gqll*gqll;

//model = gpll*mll*conj(gqll) + gpll*mlr*conj(gqlr) + gplr*mrl*conj(gqll) + gplr*mrr*conj(gqlr) ;

		} // ll

	} // iterate over polar baselines

	return 0;
}

int FullGDSolver::gradient ( const solve_model_t& pkg, 
		const complex_type& mrr, 
		const complex_type& mrl, 
		const complex_type& mlr, 
		const complex_type& mll, 
		vc_type& grad ) {
	/*
	 * This function computes gradient and saves in grad
	 *
	 * grad should be 4;
	 */

	/* zero out gradient */
	// 20260727: the caller should zero out the gradient
	// std::fill ( grad.begin(), grad.end(), complex_type(0.0f, 0.0f) );

	/* iterate over the polar baselines */
	for ( int ibl = 0; ibl < npolarbaselines; ibl++ ) {

		// fetch the antenna index
		const int iant1 ( pkg.iant1[ibl] );
		const int iant2 ( pkg.iant2[ibl] );

		// fetch the pb2corr
		const int pb2corr    = pkg.pb2corr [ ibl ];

		// fetch complex data
		const complex_type data ( pkg.data[ibl] );

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
		const complex_type gprr ( pkg.gains[iprr] );
		const complex_type gprl ( pkg.gains[iprl] );
		const complex_type gplr ( pkg.gains[iplr] );
		const complex_type gpll ( pkg.gains[ipll] );

		const complex_type gqrr ( pkg.gains[iqrr] );
		const complex_type gqrl ( pkg.gains[iqrl] );
		const complex_type gqlr ( pkg.gains[iqlr] );
		const complex_type gqll ( pkg.gains[iqll] );

		// the following long expressions come from sympy
		// see :math_selfcal.py:
		// see :math_selfcal.stdout:

		const int irr ( 0 );
		const int irl ( 1 );
		const int ilr ( 2 );
		const int ill ( 3 );

		// do on every pb2corr
		if ( pb2corr == 0 ) {

			const complex_type Drr_coeff_mrr = gplr*gqrl*conj(gplr)*conj(gqrl) + gplr*gqrr*conj(gplr)*conj(gqrr) + gprr*gqrl*conj(gprr)*conj(gqrl) + gprr*gqrr*conj(gprr)*conj(gqrr) ; 
			const complex_type Drl_coeff_mrr = gplr*gqll*conj(gplr)*conj(gqrl) + gplr*gqlr*conj(gplr)*conj(gqrr) + gprr*gqll*conj(gprr)*conj(gqrl) + gprr*gqlr*conj(gprr)*conj(gqrr) ; 
			const complex_type Dlr_coeff_mrr = gplr*gqrl*conj(gpll)*conj(gqrl) + gplr*gqrr*conj(gpll)*conj(gqrr) + gprr*gqrl*conj(gprl)*conj(gqrl) + gprr*gqrr*conj(gprl)*conj(gqrr) ; 
			const complex_type Dll_coeff_mrr = gplr*gqll*conj(gpll)*conj(gqrl) + gplr*gqlr*conj(gpll)*conj(gqrr) + gprr*gqll*conj(gprl)*conj(gqrl) + gprr*gqlr*conj(gprl)*conj(gqrr) ; 

			const complex_type Drr_coeff_drr = -gqrr*conj(gprr) ; 
			const complex_type Drl_coeff_drr = -gqlr*conj(gprr) ; 
			const complex_type Dlr_coeff_drr = -gqrr*conj(gprl) ; 
			const complex_type Dll_coeff_drr = -gqlr*conj(gprl) ; 
			
			grad[irr] += Drr_coeff_mrr*mrr  + Drr_coeff_drr*data;
			grad[irl] += Drl_coeff_mrr*mrr  + Drl_coeff_drr*data;
			grad[ilr] += Dlr_coeff_mrr*mrr  + Dlr_coeff_drr*data;
			grad[ill] += Dll_coeff_mrr*mrr  + Dll_coeff_drr*data;

		} // rr
		else if ( pb2corr == 1 ) {

			// mrl
			const complex_type Drr_coeff_mrl = gplr*gqrl*conj(gplr)*conj(gqll) + gplr*gqrr*conj(gplr)*conj(gqlr) + gprr*gqrl*conj(gprr)*conj(gqll) + gprr*gqrr*conj(gprr)*conj(gqlr) ; 
			const complex_type Drl_coeff_mrl = gplr*gqll*conj(gplr)*conj(gqll) + gplr*gqlr*conj(gplr)*conj(gqlr) + gprr*gqll*conj(gprr)*conj(gqll) + gprr*gqlr*conj(gprr)*conj(gqlr) ; 
			const complex_type Dlr_coeff_mrl = gplr*gqrl*conj(gpll)*conj(gqll) + gplr*gqrr*conj(gpll)*conj(gqlr) + gprr*gqrl*conj(gprl)*conj(gqll) + gprr*gqrr*conj(gprl)*conj(gqlr) ; 
			const complex_type Dll_coeff_mrl = gplr*gqll*conj(gpll)*conj(gqll) + gplr*gqlr*conj(gpll)*conj(gqlr) + gprr*gqll*conj(gprl)*conj(gqll) + gprr*gqlr*conj(gprl)*conj(gqlr) ; 
			
			// drl
			const complex_type Drr_coeff_drl = -gqrl*conj(gprr) ; 
			const complex_type Drl_coeff_drl = -gqll*conj(gprr) ; 
			const complex_type Dlr_coeff_drl = -gqrl*conj(gprl) ; 
			const complex_type Dll_coeff_drl = -gqll*conj(gprl) ; 

			// accumulate
			grad[irr] += Drr_coeff_mrl*mrl + Drr_coeff_drl*data;
			grad[irl] += Drl_coeff_mrl*mrl + Drl_coeff_drl*data;
			grad[ilr] += Dlr_coeff_mrl*mrl + Dlr_coeff_drl*data;
			grad[ill] += Dll_coeff_mrl*mrl + Dll_coeff_drl*data;

		} // rl
		else if ( pb2corr == 2 ) {

			// mlr
			const complex_type Drr_coeff_mlr = gpll*gqrl*conj(gplr)*conj(gqrl) + gpll*gqrr*conj(gplr)*conj(gqrr) + gprl*gqrl*conj(gprr)*conj(gqrl) + gprl*gqrr*conj(gprr)*conj(gqrr) ; 
			const complex_type Drl_coeff_mlr = gpll*gqll*conj(gplr)*conj(gqrl) + gpll*gqlr*conj(gplr)*conj(gqrr) + gprl*gqll*conj(gprr)*conj(gqrl) + gprl*gqlr*conj(gprr)*conj(gqrr) ; 
			const complex_type Dlr_coeff_mlr = gpll*gqrl*conj(gpll)*conj(gqrl) + gpll*gqrr*conj(gpll)*conj(gqrr) + gprl*gqrl*conj(gprl)*conj(gqrl) + gprl*gqrr*conj(gprl)*conj(gqrr) ; 
			const complex_type Dll_coeff_mlr = gpll*gqll*conj(gpll)*conj(gqrl) + gpll*gqlr*conj(gpll)*conj(gqrr) + gprl*gqll*conj(gprl)*conj(gqrl) + gprl*gqlr*conj(gprl)*conj(gqrr) ; 
			// dlr
			const complex_type Drr_coeff_dlr = -gqrr*conj(gplr) ; 
			const complex_type Drl_coeff_dlr = -gqlr*conj(gplr) ; 
			const complex_type Dlr_coeff_dlr = -gqrr*conj(gpll) ; 
			const complex_type Dll_coeff_dlr = -gqlr*conj(gpll) ; 

			grad[irr] += Drr_coeff_mlr*mlr + Drr_coeff_dlr*data;
			grad[irl] += Drl_coeff_mlr*mlr + Drl_coeff_dlr*data;
			grad[ilr] += Dlr_coeff_mlr*mlr + Dlr_coeff_dlr*data;
			grad[ill] += Dll_coeff_mlr*mlr + Dll_coeff_dlr*data;

		} // lr
		else if ( pb2corr == 3 ) {

			// mll
			const complex_type Drr_coeff_mll = gpll*gqrl*conj(gplr)*conj(gqll) + gpll*gqrr*conj(gplr)*conj(gqlr) + gprl*gqrl*conj(gprr)*conj(gqll) + gprl*gqrr*conj(gprr)*conj(gqlr) ; 
			const complex_type Drl_coeff_mll = gpll*gqll*conj(gplr)*conj(gqll) + gpll*gqlr*conj(gplr)*conj(gqlr) + gprl*gqll*conj(gprr)*conj(gqll) + gprl*gqlr*conj(gprr)*conj(gqlr) ; 
			const complex_type Dlr_coeff_mll = gpll*gqrl*conj(gpll)*conj(gqll) + gpll*gqrr*conj(gpll)*conj(gqlr) + gprl*gqrl*conj(gprl)*conj(gqll) + gprl*gqrr*conj(gprl)*conj(gqlr) ; 
			const complex_type Dll_coeff_mll = gpll*gqll*conj(gpll)*conj(gqll) + gpll*gqlr*conj(gpll)*conj(gqlr) + gprl*gqll*conj(gprl)*conj(gqll) + gprl*gqlr*conj(gprl)*conj(gqlr) ; 
			// dll
			const complex_type Drr_coeff_dll = -gqrl*conj(gplr) ; 
			const complex_type Drl_coeff_dll = -gqll*conj(gplr) ; 
			const complex_type Dlr_coeff_dll = -gqrl*conj(gpll) ; 
			const complex_type Dll_coeff_dll = -gqll*conj(gpll) ; 
			// accumulate
			grad[irr] += Drr_coeff_mll*mll + Drr_coeff_dll*data;
			grad[irl] += Drl_coeff_mll*mll + Drl_coeff_dll*data;
			grad[ilr] += Dlr_coeff_mll*mll + Dlr_coeff_dll*data;
			grad[ill] += Dll_coeff_mll*mll + Dll_coeff_dll*data;

		} // ll

	} // iterate over polar baselines

	return 0;
}

FullGDSolver::real_type FullGDSolver::norm ( const vc_type& g ) {
	real_type rnorm ( 0.0f );

	for ( const auto& ig : g ) rnorm += std::norm(ig);

	return rnorm;
}

int FullGDSolver::iterate ( const solve_data_t& pkg, vc_type& gains ) {
	/*
	 * This function will compute new gains using the old gains.
	 * It will not update the gains. 
	 * That will be done by the caller as it needs rate logic
	 *
	 * gains are old gains,
	 * new_gains contain the new gains
	 * both should be 4*nantennas
	 */

  // XXX this is obselete

	const complex_type mrr ( pkg.mrr );
	const complex_type mrl ( pkg.mrl );
	const complex_type mlr ( pkg.mlr );
	const complex_type mll ( pkg.mll );

	vc_type    eqn_drr ( nantennas*3, complex_type(0.0f, 0.0f) );
	vc_type    eqn_drl ( nantennas*3, complex_type(0.0f, 0.0f) );
	vc_type    eqn_dlr ( nantennas*3, complex_type(0.0f, 0.0f) );
	vc_type    eqn_dll ( nantennas*3, complex_type(0.0f, 0.0f) );

	for ( int ibl = 0; ibl < npolarbaselines; ibl++ ) {

		// fetch the antenna index
		const int iant1 ( pkg.iant1[ibl] );
		const int iant2 ( pkg.iant2[ibl] );

		// fetch the pb2corr
		const int pb2corr    = pkg.pb2corr [ ibl ];

		// fetch complex data
		const complex_type data ( pkg.data[ibl] );

		// fetch full gains for both antennas
		const complex_type gprr ( gains[4*iant1 + 0] );
		const complex_type gprl ( gains[4*iant1 + 1] );
		const complex_type gplr ( gains[4*iant1 + 2] );
		const complex_type gpll ( gains[4*iant1 + 3] );

		const complex_type gqrr ( gains[4*iant2 + 0] );
		const complex_type gqrl ( gains[4*iant2 + 1] );
		const complex_type gqlr ( gains[4*iant2 + 2] );
		const complex_type gqll ( gains[4*iant2 + 3] );


		// the following long expressions come from sympy
		// see :math_leakages.py:
		// gain coefficients are as is
		// but the data coefficients have been sign flipped
		// helps with formulating linear system later

		/*
		 * It is much easier to do over baselines, but we are doing it over polarbaselines
		 * if we compute the coefficients over every polar baseline, we will 4x counting
		 * so we do the following:
		 *
		 * to evaluate the coefficients of gp,gq, we do not need data. 
		 * We still do this over polarbaseline loop, so that we keep track of all the baselines
		 *
		 * case 'rr'
		 * 	accumulate equ_drr(gprr, one constant term), equ_drl(gprr, one constant term)
		 * case 'rl'
		 * 	accumulate equ_drr(gprl, one constant term), equ_drl(gprl, one constant term)
		 * case 'lr'
		 * 	accumulate equ_dll(gplr, one constant term), equ_dlr(gplr, one constant term)
		 * case 'll'
		 * 	accumulate equ_dll(gpll, one constant term), equ_dlr(gpll, one constant term)
		 *
		 * 	every pb2corr, we compute two gain expressions and two data expressions
		 *
		 * 	variable name convention:
		 * 	{equation_tag}_coeff_{coefficient of what?}
		 * 	{equation_tag}_const_{while solving for what?}
		*/

		/*
		 * For every polar baselines, 
		 * we have to do for both the antennas
		 *
		 * visibility is <antp x conj(antq)>
		 *
		 * our update equation is sum over all baselines
		 * our equation is solving for antp
		 * 
		 * so to solve for antq, we swap p<-->q
		 * 
		 * only place where the swap affects is in visibility
		 * it is complex conjugate.
		 *
		*/

		// do on every pb2corr
		if ( pb2corr == 0 ) {
			// ant1 case
			const complex_type drr_coeff_gprr = 
				gqll*mrl*conj(gqll)*conj(mrl) + 
				gqll*mrr*conj(gqrl)*conj(mrl) + 
				gqlr*mrl*conj(gqlr)*conj(mrl) + 
				gqlr*mrr*conj(gqrr)*conj(mrl) + 
				gqrl*mrl*conj(gqll)*conj(mrr) + 
				gqrl*mrr*conj(gqrl)*conj(mrr) + 
				gqrr*mrl*conj(gqlr)*conj(mrr) + 
				gqrr*mrr*conj(gqrr)*conj(mrr) ;

			const complex_type drl_coeff_gprr =
				gqll*mrl*conj(gqll)*conj(mll) + 
				gqll*mrr*conj(gqrl)*conj(mll) + 
				gqlr*mrl*conj(gqlr)*conj(mll) + 
				gqlr*mrr*conj(gqrr)*conj(mll) + 
				gqrl*mrl*conj(gqll)*conj(mlr) + 
				gqrl*mrr*conj(gqrl)*conj(mlr) + 
				gqrr*mrl*conj(gqlr)*conj(mlr) + 
				gqrr*mrr*conj(gqrr)*conj(mlr) ;

			const complex_type drr_const_gp = data * (gqlr*conj(mrl) + gqrr*conj(mrr));
			const complex_type drl_const_gp = data * (gqlr*conj(mll) + gqrr*conj(mlr));

			// ant2 case
			const complex_type drr_coeff_gqrr = 
				gpll*mrl*conj(gpll)*conj(mrl) + 
				gpll*mrr*conj(gprl)*conj(mrl) + 
				gplr*mrl*conj(gplr)*conj(mrl) + 
				gplr*mrr*conj(gprr)*conj(mrl) + 
				gprl*mrl*conj(gpll)*conj(mrr) + 
				gprl*mrr*conj(gprl)*conj(mrr) + 
				gprr*mrl*conj(gplr)*conj(mrr) + 
				gprr*mrr*conj(gprr)*conj(mrr) ;

			const complex_type drl_coeff_gqrr =
				gpll*mrl*conj(gpll)*conj(mll) + 
				gpll*mrr*conj(gprl)*conj(mll) + 
				gplr*mrl*conj(gplr)*conj(mll) + 
				gplr*mrr*conj(gprr)*conj(mll) + 
				gprl*mrl*conj(gpll)*conj(mlr) + 
				gprl*mrr*conj(gprl)*conj(mlr) + 
				gprr*mrl*conj(gplr)*conj(mlr) + 
				gprr*mrr*conj(gprr)*conj(mlr) ;

			const complex_type drr_const_gq = conj(data) * (gplr*conj(mrl) + gprr*conj(mrr));
			const complex_type drl_const_gq = conj(data) * (gplr*conj(mll) + gprr*conj(mlr));

			// updation because we need to sum over
			// all other baselines
			eqn_drr [3*iant1 + 0]  += drr_coeff_gprr;
			eqn_drr [3*iant1 + 2]  += drr_const_gp;

			eqn_drr [3*iant2 + 0]  += drr_coeff_gqrr;
			eqn_drr [3*iant2 + 2]  += drr_const_gq;

			eqn_drl [3*iant1 + 0]  += drl_coeff_gprr;
			eqn_drl [3*iant1 + 2]  += drl_const_gp;

			eqn_drl [3*iant2 + 0]  += drl_coeff_gqrr;
			eqn_drl [3*iant2 + 2]  += drl_const_gq;

		} // rr
		else if ( pb2corr == 1 ) {

			const complex_type drr_coeff_gprl = 
				gqll*mll*conj(gqll)*conj(mrl) +
				gqll*mlr*conj(gqrl)*conj(mrl) +
				gqlr*mll*conj(gqlr)*conj(mrl) +
				gqlr*mlr*conj(gqrr)*conj(mrl) +
				gqrl*mll*conj(gqll)*conj(mrr) +
				gqrl*mlr*conj(gqrl)*conj(mrr) +
				gqrr*mll*conj(gqlr)*conj(mrr) +
				gqrr*mlr*conj(gqrr)*conj(mrr) ;

			const complex_type drl_coeff_gprl =
				gqll*mll*conj(gqll)*conj(mll) + 
				gqll*mlr*conj(gqrl)*conj(mll) + 
				gqlr*mll*conj(gqlr)*conj(mll) + 
				gqlr*mlr*conj(gqrr)*conj(mll) + 
				gqrl*mll*conj(gqll)*conj(mlr) + 
				gqrl*mlr*conj(gqrl)*conj(mlr) + 
				gqrr*mll*conj(gqlr)*conj(mlr) + 
				gqrr*mlr*conj(gqrr)*conj(mlr) ;

			const complex_type drr_const_gp = data * (gqll*conj(mrl) + gqrl*conj(mrr));
			const complex_type drl_const_gp = data * (gqll*conj(mll) + gqrl*conj(mlr));

			const complex_type drr_coeff_gqrl = 
				gpll*mll*conj(gpll)*conj(mrl) +
				gpll*mlr*conj(gprl)*conj(mrl) +
				gplr*mll*conj(gplr)*conj(mrl) +
				gplr*mlr*conj(gprr)*conj(mrl) +
				gprl*mll*conj(gpll)*conj(mrr) +
				gprl*mlr*conj(gprl)*conj(mrr) +
				gprr*mll*conj(gplr)*conj(mrr) +
				gprr*mlr*conj(gprr)*conj(mrr) ;

			const complex_type drl_coeff_gqrl =
				gpll*mll*conj(gpll)*conj(mll) + 
				gpll*mlr*conj(gprl)*conj(mll) + 
				gplr*mll*conj(gplr)*conj(mll) + 
				gplr*mlr*conj(gprr)*conj(mll) + 
				gprl*mll*conj(gpll)*conj(mlr) + 
				gprl*mlr*conj(gprl)*conj(mlr) + 
				gprr*mll*conj(gplr)*conj(mlr) + 
				gprr*mlr*conj(gprr)*conj(mlr) ;

			const complex_type drr_const_gq = conj(data) * (gpll*conj(mrl) + gprl*conj(mrr));
			const complex_type drl_const_gq = conj(data) * (gpll*conj(mll) + gprl*conj(mlr));

			eqn_drr [3*iant1 + 1]  += drr_coeff_gprl;
			eqn_drr [3*iant1 + 2]  += drr_const_gp;

			eqn_drr [3*iant2 + 1]  += drr_coeff_gqrl;
			eqn_drr [3*iant2 + 2]  += drr_const_gq;

			eqn_drl [3*iant1 + 1]  += drl_coeff_gprl;
			eqn_drl [3*iant1 + 2]  += drl_const_gp;

			eqn_drl [3*iant2 + 1]  += drl_coeff_gqrl;
			eqn_drl [3*iant2 + 2]  += drl_const_gq;

		} // rl
		else if ( pb2corr == 2 ) {

			const complex_type dlr_coeff_gplr = 
				gqll*mrl*conj(gqll)*conj(mrl) + 
				gqll*mrr*conj(gqrl)*conj(mrl) + 
				gqlr*mrl*conj(gqlr)*conj(mrl) + 
				gqlr*mrr*conj(gqrr)*conj(mrl) + 
				gqrl*mrl*conj(gqll)*conj(mrr) + 
				gqrl*mrr*conj(gqrl)*conj(mrr) + 
				gqrr*mrl*conj(gqlr)*conj(mrr) + 
				gqrr*mrr*conj(gqrr)*conj(mrr) ;

			const complex_type dll_coeff_gplr = 
				gqll*mrl*conj(gqll)*conj(mll) + 
				gqll*mrr*conj(gqrl)*conj(mll) + 
				gqlr*mrl*conj(gqlr)*conj(mll) + 
				gqlr*mrr*conj(gqrr)*conj(mll) + 
				gqrl*mrl*conj(gqll)*conj(mlr) + 
				gqrl*mrr*conj(gqrl)*conj(mlr) + 
				gqrr*mrl*conj(gqlr)*conj(mlr) + 
				gqrr*mrr*conj(gqrr)*conj(mlr) ;

			const complex_type dlr_const_gp  = data * (gqlr*conj(mrl) + gqrr*conj(mrr));
			const complex_type dll_const_gp  = data * (gqlr*conj(mll) + gqrr*conj(mlr));

			const complex_type dlr_coeff_gqlr = 
				gpll*mrl*conj(gpll)*conj(mrl) + 
				gpll*mrr*conj(gprl)*conj(mrl) + 
				gplr*mrl*conj(gplr)*conj(mrl) + 
				gplr*mrr*conj(gprr)*conj(mrl) + 
				gprl*mrl*conj(gpll)*conj(mrr) + 
				gprl*mrr*conj(gprl)*conj(mrr) + 
				gprr*mrl*conj(gplr)*conj(mrr) + 
				gprr*mrr*conj(gprr)*conj(mrr) ;

			const complex_type dll_coeff_gqlr = 
				gpll*mrl*conj(gpll)*conj(mll) + 
				gpll*mrr*conj(gprl)*conj(mll) + 
				gplr*mrl*conj(gplr)*conj(mll) + 
				gplr*mrr*conj(gprr)*conj(mll) + 
				gprl*mrl*conj(gpll)*conj(mlr) + 
				gprl*mrr*conj(gprl)*conj(mlr) + 
				gprr*mrl*conj(gplr)*conj(mlr) + 
				gprr*mrr*conj(gprr)*conj(mlr) ;

			const complex_type dlr_const_gq  = conj(data) * (gplr*conj(mrl) + gprr*conj(mrr));
			const complex_type dll_const_gq  = conj(data) * (gplr*conj(mll) + gprr*conj(mlr));

			eqn_dlr[3*iant1 + 0]  += dlr_coeff_gplr;
			eqn_dlr[3*iant1 + 2]  += dlr_const_gp;

			eqn_dll[3*iant1 + 0]  += dll_coeff_gplr;
			eqn_dll[3*iant1 + 2]  += dll_const_gp;

			eqn_dlr[3*iant2 + 0]  += dlr_coeff_gqlr;
			eqn_dlr[3*iant2 + 2]  += dlr_const_gq;

			eqn_dll[3*iant2 + 0]  += dll_coeff_gqlr;
			eqn_dll[3*iant2 + 2]  += dll_const_gq;

		} // lr
		else if ( pb2corr == 3 ) {
			
			const complex_type dlr_coeff_gpll =
				gqll*mll*conj(gqll)*conj(mrl) + 
				gqll*mlr*conj(gqrl)*conj(mrl) + 
				gqlr*mll*conj(gqlr)*conj(mrl) + 
				gqlr*mlr*conj(gqrr)*conj(mrl) + 
				gqrl*mll*conj(gqll)*conj(mrr) + 
				gqrl*mlr*conj(gqrl)*conj(mrr) + 
				gqrr*mll*conj(gqlr)*conj(mrr) + 
				gqrr*mlr*conj(gqrr)*conj(mrr) ;

			const complex_type dll_coeff_gpll = 
				gqll*mll*conj(gqll)*conj(mll) + 
				gqll*mlr*conj(gqrl)*conj(mll) + 
				gqlr*mll*conj(gqlr)*conj(mll) + 
				gqlr*mlr*conj(gqrr)*conj(mll) + 
				gqrl*mll*conj(gqll)*conj(mlr) + 
				gqrl*mlr*conj(gqrl)*conj(mlr) + 
				gqrr*mll*conj(gqlr)*conj(mlr) + 
				gqrr*mlr*conj(gqrr)*conj(mlr) ;

			const complex_type dlr_coeff_gqll =
				gpll*mll*conj(gpll)*conj(mrl) + 
				gpll*mlr*conj(gprl)*conj(mrl) + 
				gplr*mll*conj(gplr)*conj(mrl) + 
				gplr*mlr*conj(gprr)*conj(mrl) + 
				gprl*mll*conj(gpll)*conj(mrr) + 
				gprl*mlr*conj(gprl)*conj(mrr) + 
				gprr*mll*conj(gplr)*conj(mrr) + 
				gprr*mlr*conj(gprr)*conj(mrr) ;

			const complex_type dll_coeff_gqll = 
				gpll*mll*conj(gpll)*conj(mll) + 
				gpll*mlr*conj(gprl)*conj(mll) + 
				gplr*mll*conj(gplr)*conj(mll) + 
				gplr*mlr*conj(gprr)*conj(mll) + 
				gprl*mll*conj(gpll)*conj(mlr) + 
				gprl*mlr*conj(gprl)*conj(mlr) + 
				gprr*mll*conj(gplr)*conj(mlr) + 
				gprr*mlr*conj(gprr)*conj(mlr) ;

			const complex_type dlr_const_gp  = data * (gqll*conj(mrl) + gqrl*conj(mrr));
			const complex_type dll_const_gp  = data * (gqll*conj(mll) + gqrl*conj(mlr));

			const complex_type dlr_const_gq  = conj(data) * (gpll*conj(mrl) + gprl*conj(mrr));
			const complex_type dll_const_gq  = conj(data) * (gpll*conj(mll) + gprl*conj(mlr));

			eqn_dlr[3*iant1 + 1] += dlr_coeff_gpll;
			eqn_dll[3*iant1 + 1] += dll_coeff_gpll;

			eqn_dlr[3*iant1 + 2] += dlr_const_gp;
			eqn_dll[3*iant1 + 2] += dll_const_gp;

			eqn_dlr[3*iant2 + 1] += dlr_coeff_gqll;
			eqn_dll[3*iant2 + 1] += dll_coeff_gqll;

			eqn_dlr[3*iant2 + 2] += dlr_const_gq;
			eqn_dll[3*iant2 + 2] += dll_const_gq;

		} // ll

	} // iterate over polar baselines

	// solve for new gains
	// the variable naming is defined in fullgdsolver.hpp:64
	for (int iant = 0; iant < nantennas; iant++) {

		const complex_type a1 ( eqn_drr [ 3*iant + 0 ] ); 
		const complex_type b1 ( eqn_drr [ 3*iant + 1 ] ); 
		const complex_type c1 ( eqn_drr [ 3*iant + 2 ] ); 

		const complex_type a2 ( eqn_drl [ 3*iant + 0 ] ); 
		const complex_type b2 ( eqn_drl [ 3*iant + 1 ] ); 
		const complex_type c2 ( eqn_drl [ 3*iant + 2 ] ); 

		const complex_type a3 ( eqn_dlr [ 3*iant + 0 ] ); 
		const complex_type b3 ( eqn_dlr [ 3*iant + 1 ] ); 
		const complex_type c3 ( eqn_dlr [ 3*iant + 2 ] ); 

		const complex_type a4 ( eqn_dll [ 3*iant + 0 ] ); 
		const complex_type b4 ( eqn_dll [ 3*iant + 1 ] ); 
		const complex_type c4 ( eqn_dll [ 3*iant + 2 ] ); 

		const complex_type det12 ( a1*b2 - a2*b1 );
		const complex_type det34 ( a3*b4 - a4*b3 );

		const complex_type new_grr ( ( b2*c1 - b1*c2 ) / det12 );
		const complex_type new_grl ( ( c2*a1 - c1*a2 ) / det12 );

		const complex_type new_glr ( ( b4*c3 - b3*c4 ) / det34 );
		const complex_type new_gll ( ( c4*a3 - c3*a4 ) / det34 );

		gains [ 4*iant + 0 ] = new_grr;
		gains [ 4*iant + 1 ] = new_grl;
		gains [ 4*iant + 2 ] = new_glr;
		gains [ 4*iant + 3 ] = new_gll;
	} // for every ant
	
	return 0;
}

/*
 * XXX maybe in future iterations, the cost computation can be done with :iteration:
 *
 * It puts a lot of math
*/
FullGDSolver::real_type FullGDSolver::cost ( const solve_data_t& pkg, const vc_type& gains ) {

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
		const complex_type gprr ( gains[4*iant1 + 0] );
		const complex_type gprl ( gains[4*iant1 + 1] );
		const complex_type gplr ( gains[4*iant1 + 2] );
		const complex_type gpll ( gains[4*iant1 + 3] );

		const complex_type gqrr ( gains[4*iant2 + 0] );
		const complex_type gqrl ( gains[4*iant2 + 1] );
		const complex_type gqlr ( gains[4*iant2 + 2] );
		const complex_type gqll ( gains[4*iant2 + 3] );

		// model forward depends on pb2corr
		complex_type model;

		if ( pb2corr == 0 ) {

			// gpra   mab   conj(gqbr)
			model = 
				(gprr * mrr * conj(gqrr)) +
				(gprr * mrl * conj(gqlr)) +
				(gprl * mlr * conj(gqrr)) +
				(gprl * mll * conj(gqlr)) ;

		} // rr 
		else if ( pb2corr == 1 ) {

			// gpra   mab   conj(gqbl)
			model = 
				(gprr * mrr * conj(gqrl)) +
				(gprr * mrl * conj(gqll)) +
				(gprl * mlr * conj(gqrl)) +
				(gprl * mll * conj(gqll)) ;

		} // rl
		else if ( pb2corr == 2 ) {

			// gpla   mab   conj(gqbr)
			model = 
				(gplr * mrr * conj(gqrr)) +
				(gplr * mrl * conj(gqlr)) +
				(gpll * mlr * conj(gqrr)) +
				(gpll * mll * conj(gqlr)) ;

		} // lr
		else if ( pb2corr == 3 ) {

			// gpla   mab   conj(gqbl)
			model = 
				(gplr * mrr * conj(gqrl)) +
				(gplr * mrl * conj(gqll)) +
				(gpll * mlr * conj(gqrl)) +
				(gpll * mll * conj(gqll)) ;

		} // ll
		
		// update cost
		cost += std::norm ( data - model );

	} // for every polarbaseline

	return cost;
}

FullGDSolver::real_type FullGDSolver::cost ( const solve_model_t& pkg, 
		const complex_type& mrr, const complex_type& mrl, const complex_type& mlr, const complex_type& mll ) {

	real_type cost ( 0.0f );

	for (int ibl = 0; ibl < npolarbaselines; ibl++) {

		// fetch the antenna index
		const int iant1 ( pkg.iant1[ibl] );
		const int iant2 ( pkg.iant2[ibl] );

		// fetch the pb2corr
		const int pb2corr    = pkg.pb2corr [ ibl ];

		// fetch complex data
		const complex_type data ( pkg.data[ibl] );

		// fetch full gains for both antennas
		const complex_type gprr ( pkg.gains[4*iant1 + 0] );
		const complex_type gprl ( pkg.gains[4*iant1 + 1] );
		const complex_type gplr ( pkg.gains[4*iant1 + 2] );
		const complex_type gpll ( pkg.gains[4*iant1 + 3] );

		const complex_type gqrr ( pkg.gains[4*iant2 + 0] );
		const complex_type gqrl ( pkg.gains[4*iant2 + 1] );
		const complex_type gqlr ( pkg.gains[4*iant2 + 2] );
		const complex_type gqll ( pkg.gains[4*iant2 + 3] );

		// model forward depends on pb2corr
		complex_type model;

		if ( pb2corr == 0 ) {

			// gpra   mab   conj(gqbr)
			model = 
				(gprr * mrr * conj(gqrr)) +
				(gprr * mrl * conj(gqlr)) +
				(gprl * mlr * conj(gqrr)) +
				(gprl * mll * conj(gqlr)) ;

		} // rr 
		else if ( pb2corr == 1 ) {

			// gpra   mab   conj(gqbl)
			model = 
				(gprr * mrr * conj(gqrl)) +
				(gprr * mrl * conj(gqll)) +
				(gprl * mlr * conj(gqrl)) +
				(gprl * mll * conj(gqll)) ;

		} // rl
		else if ( pb2corr == 2 ) {

			// gpla   mab   conj(gqbr)
			model = 
				(gplr * mrr * conj(gqrr)) +
				(gplr * mrl * conj(gqlr)) +
				(gpll * mlr * conj(gqrr)) +
				(gpll * mll * conj(gqlr)) ;

		} // lr
		else if ( pb2corr == 3 ) {

			// gpla   mab   conj(gqbl)
			model = 
				(gplr * mrr * conj(gqrl)) +
				(gplr * mrl * conj(gqll)) +
				(gpll * mlr * conj(gqrl)) +
				(gpll * mll * conj(gqll)) ;

		} // ll
		
		// update cost
		cost += std::norm ( data - model );

	} // for every polarbaseline

	return cost;
}

FullGDSolver::real_type FullGDSolver::solve ( const solve_data_t& pkg, vc_type& gains ) {
	rcode = 0;
	niter = 0;
	gnorm = 0.0f;

	real_type rcost (0.0f);

	// EMA of square of norm of gradient
	real_type ema_gnorm ( 0.0f );
	real_type ema_gnorm_slow ( 0.0f );

	Adam     apple ( ngains, 0.01f, 0.90f, 0.95f, 1000 );
	vc_type  grad ( ngains, complex_type(0.0f, 0.0f) );

	for ( int iter = 0; iter < max_iterations; iter++ ) {

		// find cost before iteration
		const real_type old_cost = cost ( pkg, gains );

		// zero out before solving
		std::fill ( grad.begin(), grad.end(), complex_type(0.0f, 0.0f) );

		// find gradient
		gradient ( pkg, gains, grad );

		// gradient norm
		gnorm = norm ( grad );

		// EMA of gnorm
		ema_gnorm      = betag*gnorm + (1.0f - betag)*ema_gnorm;
		ema_gnorm_slow = beta_gnorm_slow*gnorm + (1.0f - beta_gnorm_slow)*ema_gnorm_slow;

		// use ADAM to update gains
		// in place updation
		/*
		 * Instead of using one fixed alpha throughout the iterations, 
		 * let us use Adam strategy to update the ``learning rate''. 
		 * We will also pick one for every `gain`. 
		 * So that we get maximum granularity.
		 *
		 * This and more is in Adam.
		*/
		apple ( grad, gains );

		// find cost after iteration
		const real_type new_cost = cost ( pkg, gains );

		//std::cout << iter << " " << new_cost << " " << gnorm << " " << ema_gnorm << " " << ema_gnorm_slow << std::endl;
#ifdef CHANDEBUG
		std::cout << iter << " " << new_cost << " " << gnorm << " " << ema_gnorm << " " << ema_gnorm_slow << std::endl;
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
		// we do not want cost plateau condition.
		if ( std::abs ( ema_gnorm - ema_gnorm_slow ) <= gamma ) {
			rcode  = 2;
			rcost  = new_cost;
		  break;
		}
#endif

		niter++;
	}

	return rcost;
}

FullGDSolver::real_type FullGDSolver::solve ( const solve_model_t& pkg, complex_type& mrr, complex_type& mrl, complex_type& mlr, complex_type& mll ) {
	rcode = 0;
	niter = 0;

	real_type rcost (0.0f);

	// EMA of square of norm of gradient
	real_type ema_gnorm ( 0.0f );

	// EMAs of cost 
	real_type ema_cost_fast ( 0.0f );
	real_type ema_cost_slow ( 0.0f );

	const int npar ( 4 );

	Adam     apple ( npar, 0.01f, 0.90f, 0.99f, 1000 );
	vc_type  grad ( npar, complex_type(0.0f, 0.0f) );
	vc_type  updated_m ( npar, complex_type(0.0f, 0.0f) );

	for ( int iter = 0; iter < max_iterations; iter++ ) {

		// find cost before iteration
		const real_type old_cost = cost ( pkg, mrr, mrl, mlr, mll );

		// zero out the gradient
		std::fill ( grad.begin(), grad.end(), complex_type(0.0f, 0.0f) );

		// find gradient
		gradient ( pkg, mrr, mrl, mlr, mll, grad );

		// gradient norm
		const real_type gnorm = norm ( grad );

		// EMA of gnorm
		ema_gnorm  = betag*ema_gnorm + (1.0f - betag)*gnorm;

		// use ADAM to update gains
		// in place updation
		/*
		 * Instead of using one fixed alpha throughout the iterations, 
		 * let us use Adam strategy to update the ``learning rate''. 
		 * We will also pick one for every `gain`. 
		 * So that we get maximum granularity.
		 *
		 * This and more is in Adam.
		*/
		updated_m[0] = mrr; updated_m[1] = mrl;
		updated_m[2] = mlr; updated_m[3] = mll;
		apple ( grad, updated_m );
		mrr = updated_m[0]; mrl = updated_m[1];
		mlr = updated_m[2]; mll = updated_m[3];

		// find cost after iteration
		const real_type new_cost = cost ( pkg, mrr, mrl, mlr, mll );

		// EMAs of new cost
		//ema_cost_fast = beta_cost_fast*ema_cost_fast + (1.0f - beta_cost_fast)*new_cost;
		//ema_cost_slow = beta_cost_slow*ema_cost_slow + (1.0f - beta_cost_slow)*new_cost;

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
		if ( std::abs ( ema_cost_fast - ema_cost_slow) <= gamma ) {
			rcode  = 2;
			rcost  = new_cost;
			break;
		}
#endif

		niter++;
	}

	return rcost;
}

FullGDSolver::real_type FullGDSolver::solve ( const solve_data_t& pkg1, const solve_data_t& pkg2, vc_type& gains ) {
	rcode = 0;
	niter = 0;

	real_type rcost (0.0f);

	// EMA of square of norm of gradient
	real_type ema_gnorm ( 0.0f );

	// EMAs of cost 
	real_type ema_cost_fast ( 0.0f );
	real_type ema_cost_slow ( 0.0f );

	Adam     apple ( ngains, 0.01f, 0.90f, 0.99f, 1000 );
	vc_type  grad ( ngains, complex_type(0.0f, 0.0f) );

	for ( int iter = 0; iter < max_iterations; iter++ ) {

		// find cost before iteration
		const real_type old_cost = cost ( pkg1, gains ) + cost ( pkg2, gains );

		// zero out before solving
		std::fill ( grad.begin(), grad.end(), complex_type(0.0f, 0.0f) );

		// find gradient
		gradient ( pkg1, gains, grad );
		gradient ( pkg2, gains, grad );

		// gradient norm
		const real_type gnorm = norm ( grad );

		// EMA of gnorm
		ema_gnorm  = betag*ema_gnorm + (1.0f - betag)*gnorm;

		// use ADAM to update gains
		// in place updation
		/*
		 * Instead of using one fixed alpha throughout the iterations, 
		 * let us use Adam strategy to update the ``learning rate''. 
		 * We will also pick one for every `gain`. 
		 * So that we get maximum granularity.
		 *
		 * This and more is in Adam.
		*/
		apple ( grad, gains );

		// find cost after iteration
		const real_type new_cost = cost ( pkg1, gains ) + cost ( pkg2, gains );

		// EMAs of new cost
		//ema_cost_fast = beta_cost_fast*ema_cost_fast + (1.0f - beta_cost_fast)*new_cost;
		//ema_cost_slow = beta_cost_slow*ema_cost_slow + (1.0f - beta_cost_slow)*new_cost;

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
		if ( std::abs ( ema_cost_fast - ema_cost_slow) <= gamma ) {
			rcode  = 2;
			rcost  = new_cost;
			break;
		}
#endif

		niter++;
	}

	return rcost;
}
