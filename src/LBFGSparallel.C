#include "LBFGSclean.hpp"

/*
 * parallel
*/
#ifdef DPRINT
#include <iostream>
#endif


LBFGS::real_type LBFGS::parallel_unpolarized (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * We assume the unpolarized source has I of unity.
	 *
	 * Will only solve rr and ll.
	 *
	 * phases of the parallel gains of the reference antenna (first antenna) is zero.
	 * ==> real part is always positive
	 * This constraint built into the gradient computation (see GREF)
	 *
	 * Each antenna has two gains, rr and ll. LBFGS requires real parameters. We express the two gains are
	 * __ rgains layout 
	 * R I         R I  
	 * | one antenna |
	 * |rr,rl, lr, ll|
	 *
	 * GREF layout
	 * R R | R I R I R I R I | ....
	 * rr
	 *
	 * n = nantennas*2*2 - 2
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = reinterpret_cast<const uata_t*>(instance);

	/* zero out gradient */
	std::fill ( rgrad, rgrad + n, 0.0f );

	const int    ngains ( 2 * pkg->nantennas );

	/* populate complex gains vector */
	/* despite GREF, the gains are kept as complex to be uniform */
	vc_type  cgains ( ngains, complex_type(0.0f, 0.0f) );
	// rr
	cgains[0]   = complex_type ( std::exp(rgains[0]), 0.0f );
	// ll
	cgains[1]   = complex_type ( std::exp(rgains[1]), 0.0f );
	for ( int igain = 2; igain < ngains; igain++ ) {
		cgains[igain]   = complex_type ( rgains[2*igain-2], rgains[2*igain-1] );
	}

	//std::cout << " cgains=";
	//for ( int ig = 0; ig < ngains; ig++ ) {
		//std::cout <<  cgains[ig] << " ";
	//}
	//std::cout << std::endl;

	/* create complex grad vector */
	vc_type  grad ( ngains, complex_type(0.0f, 0.0f) );

	/* iterate over the polar baselines */
	for ( int ibl = 0; ibl < pkg->npolarbaselines; ibl++ ) {

		// fetch the antenna index
		const int iant1 ( pkg->iant1[ibl] );
		const int iant2 ( pkg->iant2[ibl] );

		// fetch the pb2corr
		const int pb2corr    = pkg->pb2corr [ ibl ];
		// we are only interested in rr and ll
		if (pb2corr != 0 && pb2corr != 3) continue;

		// fetch complex data
		const complex_type data ( pkg->data[ibl] );
		const complex_type cata ( conj(data) );

		// fetch the parallactic angle skip
		const complex_type z1 ( pkg->par_z1[ibl] );
		const complex_type z2 ( pkg->par_z2[ibl] );

		// set the gain indices
		const int iprr ( 2*iant1 + 0 );
		const int ipll ( 2*iant1 + 1 );

		const int iqrr ( 2*iant2 + 0 );
		const int iqll ( 2*iant2 + 1 );

		// fetch full gains for both antennas
		const complex_type gprr ( cgains[iprr] );
		const complex_type gpll ( cgains[ipll] );

		const complex_type gqrr ( cgains[iqrr] );
		const complex_type gqll ( cgains[iqll] );

		// the following long expressions come from sympy
		// see :math_gradient_parallel_unity.py:
		// see :math_gradient_parallel_unity.code:
		//if ( iant1 == 0 ) {
			//std::cout << gprr << " " << gpll << " " << gqrr << " " << gqll << " " << data << std::endl;
		//}

		// model forward depends on pb2corr
		complex_type model;

		/*
		We do this over polarbaseline loop, so that we keep track of all the baselines
		*/
		if ( pb2corr == 0 ) {

const complex_type Dgprr_coeff_dpqrr = -gqrr*z1*conj(z2) ; 
const complex_type Dgprr_coeff_gprr = gqrr*z1*z2*conj(gqrr)*conj(z1)*conj(z2) ; 
const complex_type Dgqrr_coeff_dqprr = -gprr*z2*conj(z1) ; 
const complex_type Dgqrr_coeff_gqrr = gprr*z1*z2*conj(gprr)*conj(z1)*conj(z2) ; 

			grad[iprr] += Dgprr_coeff_dpqrr*data + Dgprr_coeff_gprr*gprr;
			grad[iqrr] += Dgqrr_coeff_dqprr*cata + Dgqrr_coeff_gqrr*gqrr;

			//if (iant1 == 0) {
				//std::cout << "[pb2corr=0] " << Dgprr_coeff_dpqrr << " " << Dgprr_coeff_gprr << " " << grad[iprr] << std::endl;;
			//}

model = gprr*z2*conj(gqrr)*conj(z1) ;

		} // rr
		else if ( pb2corr == 3 ) {

const complex_type Dgpll_coeff_dpqll = -gqll*z2*conj(z1) ; 
const complex_type Dgpll_coeff_gpll = gqll*z1*z2*conj(gqll)*conj(z1)*conj(z2) ; 
const complex_type Dgqll_coeff_dqpll = -gpll*z1*conj(z2) ; 
const complex_type Dgqll_coeff_gqll = gpll*z1*z2*conj(gpll)*conj(z1)*conj(z2) ; 
			
			grad[ipll] += Dgpll_coeff_dpqll*data + Dgpll_coeff_gpll*gpll;
			grad[iqll] += Dgqll_coeff_dqpll*cata + Dgqll_coeff_gqll*gqll;

			//if (iant1 == 0) {
				//std::cout << "[pb2corr=3] " << Dgpll_coeff_dpqll << " " << Dgpll_coeff_gpll << " " << grad[ipll] << std::endl;;
			//}

model = gpll*z1*conj(gqll)*conj(z2) ;
		} // ll
			
		//std::cout << " iterationcost=" << cost << " grad @ " << iprr << " " << ipll << " " << grad[iprr] << " " << grad[iqll] <<   std::endl;
		// update cost
		cost += std::norm ( data - model );

	} // iterate over polar baselines
	

	/* load complex grad into real and imaginary parts */
	// chain rule because of GREF
	// rr
	// 20260917: wirtinger derivative to real/imag
	// need a factor of two
	rgrad[0] = 2.0f * grad[0].real() * std::exp(rgains[0]);
	// ll
	rgrad[1] = 2.0f * grad[1].real() * std::exp(rgains[1]);
	// for the rest of grad
	for ( int igain = 2; igain < ngains; igain++ ) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain - 2] = 2.0f * gg.real();
		rgrad[2*igain - 1] = 2.0f * gg.imag();
	}

	//std::cout << " full_jones_cost=" << cost << std::endl; 
	
	//std::cout << " rgrads=";
	//for ( int ig = 0; ig < n; ig++ ) {
		//std::cout <<  rgrad[ig] << " ";
	//}
	//std::cout << std::endl;

	return cost;
}
