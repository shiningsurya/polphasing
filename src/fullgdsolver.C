#include "fullgdsolver.hpp"

using std::conj;

int FullGDSolver::gradient ( const solve_data_t& pkg, const vc_type& gains, vc_type& grad ) {
	/*
	 * This function computes gradient and saves in grad
	 * both should be 4*nantennas
	 *
	 * grad should be 4*nantennas;
	 */

	/* load model parameters */
	const complex_type mrr ( pkg.mrr );
	const complex_type mrl ( pkg.mrl );
	const complex_type mlr ( pkg.mlr );
	const complex_type mll ( pkg.mll );

	/* zero out gradient */
	std::fill ( grad.begin(), grad.end(), complex_type(0.0f, 0.0f) );

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
		 * It is much easier to do over baselines, but we are doing it over polarbaselines
		 * if we compute the coefficients over every polar baseline, we will overcount 4x
		 * so we do the following:
		 *
		 * We do this over polarbaseline loop, so that we keep track of all the baselines
		 *
		 * small d denotes data. Big D denotes derivative.
		 * There are four gradient terms per antenna.
		 * The ant1 is denoted by p. ant2 by q.
		 *
		 * The four gradient terms are denoted by Dg{pq}{rr,rl,lr,ll}.
		 *
		 * for every polarbaseline:
		 *
		 * case 'rr'
		 * 	accumulate coefficients of drrpq and gprr in Dgprr and Dgprl, 
		 * case 'rl'
		 * 	accumulate coefficients of drlpq and gprl in Dgprr and Dgprl, 
		 * case 'lr'
		 * 	accumulate coefficients of dlrpq and gplr in Dgplr and Dgpll, 
		 * case 'll'
		 * 	accumulate coefficients of dllpq and gpll in Dgplr and Dgpll, 
		 *
		 * 	Each of the gradient term has four terms in it.
		 * 	Two are data terms, and two are gain terms.
		 *
		 * 	Every pb2corr, we compute two gain expressions and two data expressions
		 * 	So every pb2corr we update one gain and one data expression in two gradient term.
		 *
		 * 	variable name convention:
		 * 	gradient_tags = {Dgprr, Dgprl, Dgplr, Dgpll}
		 * 	coeff_tag     = { drr, drl, dlr, dll } | {gprr, gprl, gplr, gpll}
		 * 	d has pq baseline which is assumed.
		 *
		 * 	{gradient_tag}_coeff_{coeff_tag}
		 *
		 * With this naming convention, 
		 * index should match gradient_tag
		 * coefficient tag should match coefficient
		*/

		/*
		 * For every polar baselines, 
		 * we have to do for both the antennas - ant1 and ant2
		 *
		 * Solving for one gain is summing over all the baselines. 
		 * The data layout is that ant1 < ant2. 
		 * So when iterating, we compute ant1 and ant2.
		 * For ant2, we will take conjugate.
		 *
		 * visibility is <antp x conj(antq)>
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
			const complex_type Dgprr_coeff_gprr = 
				gqll*mrl*conj(gqll)*conj(mrl) + 
				gqll*mrr*conj(gqrl)*conj(mrl) + 
				gqlr*mrl*conj(gqlr)*conj(mrl) + 
				gqlr*mrr*conj(gqrr)*conj(mrl) + 
				gqrl*mrl*conj(gqll)*conj(mrr) + 
				gqrl*mrr*conj(gqrl)*conj(mrr) + 
				gqrr*mrl*conj(gqlr)*conj(mrr) + 
				gqrr*mrr*conj(gqrr)*conj(mrr) ;

			const complex_type Dgprl_coeff_gprr =
				gqll*mrl*conj(gqll)*conj(mll) + 
				gqll*mrr*conj(gqrl)*conj(mll) +
				gqlr*mrl*conj(gqlr)*conj(mll) +
				gqlr*mrr*conj(gqrr)*conj(mll) +
				gqrl*mrl*conj(gqll)*conj(mlr) +
				gqrl*mrr*conj(gqrl)*conj(mlr) +
				gqrr*mrl*conj(gqlr)*conj(mlr) +
				gqrr*mrr*conj(gqrr)*conj(mlr) ;

			const complex_type Dgprr_coeff_drr = -gqlr*conj(mrl) - gqrr*conj(mrr);

			const complex_type Dgprl_coeff_drr = -gqlr*conj(mll) - gqrr*conj(mlr);

			const complex_type Dgqrr_coeff_gqrr = 
				gpll*mrl*conj(gpll)*conj(mrl) + 
				gpll*mrr*conj(gprl)*conj(mrl) + 
				gplr*mrl*conj(gplr)*conj(mrl) + 
				gplr*mrr*conj(gprr)*conj(mrl) + 
				gprl*mrl*conj(gpll)*conj(mrr) + 
				gprl*mrr*conj(gprl)*conj(mrr) + 
				gprr*mrl*conj(gplr)*conj(mrr) + 
				gprr*mrr*conj(gprr)*conj(mrr) ;

			const complex_type Dgqrl_coeff_gqrr =
				gpll*mrl*conj(gpll)*conj(mll) + 
				gpll*mrr*conj(gprl)*conj(mll) +
				gplr*mrl*conj(gplr)*conj(mll) +
				gplr*mrr*conj(gprr)*conj(mll) +
				gprl*mrl*conj(gpll)*conj(mlr) +
				gprl*mrr*conj(gprl)*conj(mlr) +
				gprr*mrl*conj(gplr)*conj(mlr) +
				gprr*mrr*conj(gprr)*conj(mlr) ;

			const complex_type Dgqrr_coeff_drr = -gplr*conj(mrl) - gprr*conj(mrr);

			const complex_type Dgqrl_coeff_drr = -gplr*conj(mll) - gprr*conj(mlr);

			/* update gradiant */
			// index should match gradient_tag
			// coefficient tag should match coefficient
			grad[iprr] += Dgprr_coeff_gprr*gprr + Dgprr_coeff_drr*data; 
			grad[iprl] += Dgprl_coeff_gprr*gprr + Dgprl_coeff_drr*data;

			/* update gradiant */
			// index should match gradient_tag
			// coefficient tag should match coefficient
			grad[iqrr] += Dgqrr_coeff_gqrr*gqrr + Dgqrr_coeff_drr*conj(data); 
			grad[iqrl] += Dgqrl_coeff_gqrr*gqrr + Dgqrl_coeff_drr*conj(data);

		} // rr
		else if ( pb2corr == 1 ) {
			// ant1

			const complex_type Dgprr_coeff_gprl = 
				gqll*mll*conj(gqll)*conj(mrl) +
				gqll*mlr*conj(gqrl)*conj(mrl) +
				gqlr*mll*conj(gqlr)*conj(mrl) +
				gqlr*mlr*conj(gqrr)*conj(mrl) +
				gqrl*mll*conj(gqll)*conj(mrr) +
				gqrl*mlr*conj(gqrl)*conj(mrr) +
				gqrr*mll*conj(gqlr)*conj(mrr) +
				gqrr*mlr*conj(gqrr)*conj(mrr) ;

			const complex_type Dgprl_coeff_gprl = 
				gqll*mll*conj(gqll)*conj(mll) +
				gqll*mlr*conj(gqrl)*conj(mll) +
				gqlr*mll*conj(gqlr)*conj(mll) +
				gqlr*mlr*conj(gqrr)*conj(mll) +
				gqrl*mll*conj(gqll)*conj(mlr) +
				gqrl*mlr*conj(gqrl)*conj(mlr) +
				gqrr*mll*conj(gqlr)*conj(mlr) +
				gqrr*mlr*conj(gqrr)*conj(mlr) ;

			const complex_type Dgprr_coeff_drl = -gqll*conj(mrl) - gqrl*conj(mrr);
			const complex_type Dgprl_coeff_drl = -gqll*conj(mll) - gqrl*conj(mlr);

			// just swap p<-->q
			const complex_type Dgqrr_coeff_gqrl = 
				gpll*mll*conj(gpll)*conj(mrl) +
				gpll*mlr*conj(gprl)*conj(mrl) +
				gplr*mll*conj(gplr)*conj(mrl) +
				gplr*mlr*conj(gprr)*conj(mrl) +
				gprl*mll*conj(gpll)*conj(mrr) +
				gprl*mlr*conj(gprl)*conj(mrr) +
				gprr*mll*conj(gplr)*conj(mrr) +
				gprr*mlr*conj(gprr)*conj(mrr) ;

			const complex_type Dgqrl_coeff_gqrl = 
				gpll*mll*conj(gpll)*conj(mll) +
				gpll*mlr*conj(gprl)*conj(mll) +
				gplr*mll*conj(gplr)*conj(mll) +
				gplr*mlr*conj(gprr)*conj(mll) +
				gprl*mll*conj(gpll)*conj(mlr) +
				gprl*mlr*conj(gprl)*conj(mlr) +
				gprr*mll*conj(gplr)*conj(mlr) +
				gprr*mlr*conj(gprr)*conj(mlr) ;

			const complex_type Dgqrr_coeff_drl = -gpll*conj(mrl) - gprl*conj(mrr);
			const complex_type Dgqrl_coeff_drl = -gpll*conj(mll) - gprl*conj(mlr);

			/* update gradiant */
			// index should match gradient_tag
			// coefficient tag should match coefficient
			grad[iprr] += Dgprr_coeff_gprl*gprl + Dgprr_coeff_drl*data; 
			grad[iprl] += Dgprl_coeff_gprl*gprl + Dgprl_coeff_drl*data;

			grad[iqrr] += Dgqrr_coeff_gqrl*gqrl + Dgqrr_coeff_drl*conj(data); 
			grad[iqrl] += Dgqrl_coeff_gqrl*gqrl + Dgqrl_coeff_drl*conj(data);

		} // rl
		else if ( pb2corr == 2 ) {

			const complex_type Dgplr_coeff_gplr = 
				gqll*mrl*conj(gqll)*conj(mrl) +
				gqll*mrr*conj(gqrl)*conj(mrl) +
				gqlr*mrl*conj(gqlr)*conj(mrl) +
				gqlr*mrr*conj(gqrr)*conj(mrl) +
				gqrl*mrl*conj(gqll)*conj(mrr) +
				gqrl*mrr*conj(gqrl)*conj(mrr) +
				gqrr*mrl*conj(gqlr)*conj(mrr) +
				gqrr*mrr*conj(gqrr)*conj(mrr) ;

			const complex_type Dgpll_coeff_gplr = 
				gqll*mrl*conj(gqll)*conj(mll) +
				gqll*mrr*conj(gqrl)*conj(mll) +
				gqlr*mrl*conj(gqlr)*conj(mll) +
				gqlr*mrr*conj(gqrr)*conj(mll) +
				gqrl*mrl*conj(gqll)*conj(mlr) +
				gqrl*mrr*conj(gqrl)*conj(mlr) +
				gqrr*mrl*conj(gqlr)*conj(mlr) +
				gqrr*mrr*conj(gqrr)*conj(mlr) ;

			const complex_type Dgplr_coeff_dlr = -gqlr*conj(mrl) - gqrr*conj(mrr);

			const complex_type Dgpll_coeff_dlr = -gqlr*conj(mll) - gqrr*conj(mlr);

			const complex_type Dgqlr_coeff_gqlr = 
				gpll*mrl*conj(gpll)*conj(mrl) +
				gpll*mrr*conj(gprl)*conj(mrl) +
				gplr*mrl*conj(gplr)*conj(mrl) +
				gplr*mrr*conj(gprr)*conj(mrl) +
				gprl*mrl*conj(gpll)*conj(mrr) +
				gprl*mrr*conj(gprl)*conj(mrr) +
				gprr*mrl*conj(gplr)*conj(mrr) +
				gprr*mrr*conj(gprr)*conj(mrr) ;

			const complex_type Dgqll_coeff_gqlr = 
				gpll*mrl*conj(gpll)*conj(mll) +
				gpll*mrr*conj(gprl)*conj(mll) +
				gplr*mrl*conj(gplr)*conj(mll) +
				gplr*mrr*conj(gprr)*conj(mll) +
				gprl*mrl*conj(gpll)*conj(mlr) +
				gprl*mrr*conj(gprl)*conj(mlr) +
				gprr*mrl*conj(gplr)*conj(mlr) +
				gprr*mrr*conj(gprr)*conj(mlr) ;

			const complex_type Dgqlr_coeff_dlr = -gplr*conj(mrl) - gprr*conj(mrr);

			const complex_type Dgqll_coeff_dlr = -gplr*conj(mll) - gprr*conj(mlr);

			/* update gradiant */
			// index should match gradient_tag
			// coefficient tag should match coefficient
			grad[iplr] += Dgplr_coeff_gplr*gplr + Dgplr_coeff_dlr*data; 
			grad[ipll] += Dgpll_coeff_gplr*gplr + Dgpll_coeff_dlr*data;

			grad[iqlr] += Dgqlr_coeff_gqlr*gqlr + Dgqlr_coeff_dlr*conj(data); 
			grad[iqll] += Dgqll_coeff_gqlr*gqlr + Dgqll_coeff_dlr*conj(data);

		} // lr
		else if ( pb2corr == 3 ) {

			const complex_type Dgplr_coeff_gpll = 
				gqll*mll*conj(gqll)*conj(mrl) +
				gqll*mlr*conj(gqrl)*conj(mrl) +
				gqlr*mll*conj(gqlr)*conj(mrl) +
				gqlr*mlr*conj(gqrr)*conj(mrl) +
				gqrl*mll*conj(gqll)*conj(mrr) +
				gqrl*mlr*conj(gqrl)*conj(mrr) +
				gqrr*mll*conj(gqlr)*conj(mrr) +
				gqrr*mlr*conj(gqrr)*conj(mrr) ;

			const complex_type Dgpll_coeff_gpll = 
				gqll*mll*conj(gqll)*conj(mll) +
				gqll*mlr*conj(gqrl)*conj(mll) +
				gqlr*mll*conj(gqlr)*conj(mll) +
				gqlr*mlr*conj(gqrr)*conj(mll) +
				gqrl*mll*conj(gqll)*conj(mlr) +
				gqrl*mlr*conj(gqrl)*conj(mlr) +
				gqrr*mll*conj(gqlr)*conj(mlr) +
				gqrr*mlr*conj(gqrr)*conj(mlr) ;

			const complex_type Dgplr_coeff_dll = -gqll*conj(mrl) - gqrl*conj(mrr);

			const complex_type Dgpll_coeff_dll = -gqll*conj(mll) - gqrl*conj(mlr);
			
			const complex_type Dgqlr_coeff_gqll = 
				gpll*mll*conj(gpll)*conj(mrl) +
				gpll*mlr*conj(gprl)*conj(mrl) +
				gplr*mll*conj(gplr)*conj(mrl) +
				gplr*mlr*conj(gprr)*conj(mrl) +
				gprl*mll*conj(gpll)*conj(mrr) +
				gprl*mlr*conj(gprl)*conj(mrr) +
				gprr*mll*conj(gplr)*conj(mrr) +
				gprr*mlr*conj(gprr)*conj(mrr) ;

			const complex_type Dgqll_coeff_gqll = 
				gpll*mll*conj(gpll)*conj(mll) +
				gpll*mlr*conj(gprl)*conj(mll) +
				gplr*mll*conj(gplr)*conj(mll) +
				gplr*mlr*conj(gprr)*conj(mll) +
				gprl*mll*conj(gpll)*conj(mlr) +
				gprl*mlr*conj(gprl)*conj(mlr) +
				gprr*mll*conj(gplr)*conj(mlr) +
				gprr*mlr*conj(gprr)*conj(mlr) ;

			const complex_type Dgqlr_coeff_dll = -gpll*conj(mrl) - gprl*conj(mrr);

			const complex_type Dgqll_coeff_dll = -gpll*conj(mll) - gprl*conj(mlr);
			
			/* update gradiant */
			// index should match gradient_tag
			// coefficient tag should match coefficient
			grad[iplr] += Dgplr_coeff_gpll*gpll + Dgplr_coeff_dll*data;
			grad[ipll] += Dgpll_coeff_gpll*gpll + Dgpll_coeff_dll*data;

			grad[iqlr] += Dgqlr_coeff_gqll*gqll + Dgqlr_coeff_dll*conj(data);
			grad[iqll] += Dgqll_coeff_gqll*gqll + Dgqll_coeff_dll*conj(data);

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

	const complex_type mrr ( pkg.mrr );
	const complex_type mrl ( pkg.mrl );
	const complex_type mlr ( pkg.mlr );
	const complex_type mll ( pkg.mll );

	for (int ibl = 0; ibl < npolarbaselines; ibl++) {

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
		const real_type old_cost = cost ( pkg, gains );

		// find gradient
		gradient ( pkg, gains, grad );

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
		const real_type new_cost = cost ( pkg, gains );

		// EMAs of new cost
		ema_cost_fast = beta_cost_fast*ema_cost_fast + (1.0f - beta_cost_fast)*new_cost;
		ema_cost_slow = beta_cost_slow*ema_cost_slow + (1.0f - beta_cost_slow)*new_cost;

		std::cout << iter << " " << new_cost << " " << gnorm << " " << ema_gnorm << " " << ema_cost_slow << " " << ema_cost_fast << std::endl;

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

#if 0
		// termination condition
		// If the change in the cost is not a lot!
		if ( std::abs(old_cost - new_cost) <= delta ) {
			rcode  = 1;
			rcost  = new_cost;
			break;
		}
		// If the norm of the gradient is vanishing
		if ( gnorm  <= gamma ) {
			rcode  = 2;
			rcost  = new_cost;
			break;
		}
#endif
		// do not terminate on gainconvergence
		// only terminate if cost converges

		niter++;
	}

	return rcost;
}
