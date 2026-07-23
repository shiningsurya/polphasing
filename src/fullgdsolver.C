#include "fullgdsolver.hpp"

using std::conj;

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

	for ( int iter = 0; iter < max_iterations; iter++ ) {

		// find cost before iteration
		real_type old_cost = cost ( pkg, gains );

		// make copy of gains
		vc_type solutions ( gains );

		// iterate once
		iterate ( pkg, solutions );

		// update gains
		// ngains is 4 x nantennas
		for ( int igain = 0; igain < ngains; igain++ ) {
			const complex_type og = gains [ igain ];
			const complex_type ng = solutions [ igain ];

			// lerp with alpha
			gains [ igain ] = (1.0f - alpha)*og + alpha*ng;
		}

		real_type new_cost = cost ( pkg, gains );

		// termination condition
		if ( std::abs(old_cost - new_cost) <= delta ) {
			rcode = 1;
			rcost  = new_cost;
			break;
		}
		// do not terminate on gainconvergence
		// only terminate if cost converges

		niter++;
	}

	return rcost;
}
