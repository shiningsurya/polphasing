#include "LBFGSclean.hpp"

/*
 * Single unpolarized
*/
#ifdef DPRINT
#include <iostream>
#endif

template<int forcorr>
LBFGS::real_type LBFGS::SingleUnpolarizedSolver::solve_single_unpolarized (uata_t& pkg) {
	/*
	 *
	*/

	real_type cost (0.0f);

	void* vpkg  = static_cast<void*>(&pkg);

	// parallel solve 
	initialize ();
  rcode       = lbfgs(n, xpar, &cost, single_unpolarized<forcorr>, NULL, vpkg, &param);

	gnorm  = pkg.gnorm;
	niter  = pkg.niter;

#ifdef DPRINT
	std::cout << " rcode=" << rcode_para << std::endl;
#endif

	/* normalize the gains */
	/* phase of the first antenna is forced to be zero */
	complex_type gref ( xpar[0], xpar[1] );
	gref /= std::abs(gref);
	for ( int iant = 0; iant < nantennas; iant++ ) {
		const int ir ( 2*iant );
		const int ii ( 2*iant+1 );
		// read
		const complex_type cg ( xpar[ir], xpar[ii] );
		const complex_type gg ( cg / gref );
		xpar[ir]  = gg.real();
		xpar[ii]  = gg.imag();
	}

  return cost;
}

template<int forcorr>
LBFGS::real_type LBFGS::single_unpolarized (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * We assume the unpolarized source has I of unity.
	 * We only solve for a single hand, either r and l.
	 *
	 * forcorr sets the pb2corr. 0 ==> rr, 3 ==> ll
	 *
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = static_cast<const uata_t*>(instance);

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    nant ( pkg->nantennas );
	const int    ngains ( pkg->nantennas );

	/* populate complex gains vector */
	/* despite GREF, the gains are kept as complex to be uniform */
	vc_type  cgains ( ngains, complex_type(0.0f, 0.0f) );
	for ( int igain = 0; igain < ngains; igain++ ) {
		cgains[igain]   = complex_type ( rgains[2*igain  ], rgains[2*igain+1] );
	}

	/* create complex grad vector */
	vc_type  grad ( ngains, complex_type(0.0f, 0.0f) );

	/* iterate over the polar baselines */
	for ( int ibl = 0; ibl < pkg->npolarbaselines; ibl++ ) {

		// fetch the pb2corr
		const int pb2corr    = pkg->pb2corr [ ibl ];
		// we are only interested in rr and ll
		if (pb2corr != forcorr) continue;

		// fetch the antenna index
		const int iant1 ( pkg->iant1[ibl] );
		const int iant2 ( pkg->iant2[ibl] );

		// fetch complex data
		const complex_type data ( pkg->data[ibl] );
		const complex_type cata ( conj(data) );

		// fetch the parallactic angle skip
		const complex_type z1 ( pkg->par_z1[ibl] );
		const complex_type z2 ( pkg->par_z2[ibl] );

		// set the gain indices
		const int ip ( iant1 );
		const int iq ( iant2 );

		// fetch gains
		const complex_type gp ( cgains[ip] );
		const complex_type gq ( cgains[iq] );

		// the following long expressions come from sympy
		// see :math_gradient_parallel_unity.py:
		// see :math_gradient_parallel_unity.code:

		// model forward depends on pb2corr
		complex_type model;

		/*
		We do this over polarbaseline loop, so that we keep track of all the baselines

		XXX It would be ideal to have if constexpr (forcorr == 0)
		but we are sticking with c++11 which does not have that.

		So we still have to keep the if logic based on pb2corr
		which at this point in code is either 0 or 3.
		*/
		if ( pb2corr == 0 ) {

const complex_type Dgprr_coeff_dpqrr = -gq*z1*conj(z2) ; 
const complex_type Dgprr_coeff_gp = gq*z1*z2*conj(gq)*conj(z1)*conj(z2) ; 
const complex_type Dgqrr_coeff_dqprr = -gp*z2*conj(z1) ; 
const complex_type Dgqrr_coeff_gq = gp*z1*z2*conj(gp)*conj(z1)*conj(z2) ; 

			grad[ip] += Dgprr_coeff_dpqrr*data + Dgprr_coeff_gp*gp;
			grad[iq] += Dgqrr_coeff_dqprr*cata + Dgqrr_coeff_gq*gq;

model = gp*z2*conj(gq)*conj(z1) ;

		} // rr
		else if ( pb2corr == 3 ) {

const complex_type Dgpll_coeff_dpqll = -gq*z2*conj(z1) ; 
const complex_type Dgpll_coeff_gp = gq*z1*z2*conj(gq)*conj(z1)*conj(z2) ; 
const complex_type Dgqll_coeff_dqpll = -gp*z1*conj(z2) ; 
const complex_type Dgqll_coeff_gq = gp*z1*z2*conj(gp)*conj(z1)*conj(z2) ; 
			
			grad[ip] += Dgpll_coeff_dpqll*data + Dgpll_coeff_gp*gp;
			grad[iq] += Dgqll_coeff_dqpll*cata + Dgqll_coeff_gq*gq;

model = gp*z1*conj(gq)*conj(z2) ;
		} // ll
			
		// update cost
		cost += std::norm ( data - model );

	} // iterate over polar baselines
	

	/* load complex grad into real and imaginary parts */
	// chain rule because of GREF
	// 20260917: wirtinger derivative to real/imag
	// need a factor of two
	// for the rest of grad
	for ( int igain = 0; igain < ngains; igain++ ) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain + 0] = 1.0f * gg.real();
		rgrad[2*igain + 1] = 1.0f * gg.imag();
	}

	return cost;
}

// explicit instantiation
template LBFGS::real_type LBFGS::SingleUnpolarizedSolver::solve_single_unpolarized<0>(LBFGS::uata_t&);
template LBFGS::real_type LBFGS::SingleUnpolarizedSolver::solve_single_unpolarized<3>(LBFGS::uata_t&);
