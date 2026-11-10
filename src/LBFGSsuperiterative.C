#include "LBFGSsuperiterative.hpp"

//#define DPRINT

#ifdef DPRINT
#include <iostream>
#endif

LBFGS::real_type LBFGS::unpolarized_rr (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * Solving for rr gain, populates the derivative of rr that appears in {rr,rl,lr} baselines.
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = reinterpret_cast<const uata_t*>(instance);

	/* model I */
	const real_type MI ( pkg->modeli );
	const real_type MI2 ( MI * MI );

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

		const complex_type dqprr ( conj ( dpqrr ) );
		const complex_type dqprl ( conj ( dpqrl ) );
		const complex_type dqplr ( conj ( dpqlr ) );
		const complex_type dqpll ( conj ( dpqll ) );


		// fetch the parallactic angle skip
		const complex_type za ( pkg->par_za[ibl] );
		const complex_type zb ( pkg->par_zb[ibl] );

		const complex_type za2 ( za * za );
		const complex_type zb2 ( zb * zb );

		// set the gain indices
		const int iprr ( iant1 );
		const int iprl ( 4*iant1 + 1 );
		const int iplr ( 4*iant1 + 2 );
		const int ipll ( 4*iant1 + 3 );

		const int iqrr ( iant2 );
		const int iqrl ( 4*iant2 + 1 );
		const int iqlr ( 4*iant2 + 2 );
		const int iqll ( 4*iant2 + 3 );

		// fetch full gains for both antennas
		const complex_type gprr ( cgains[iprr] );
		const complex_type gprl ( pkg->cgains[iprl] );
		const complex_type gplr ( pkg->cgains[iplr] );
		const complex_type gpll ( pkg->cgains[ipll] );

		const complex_type cgprr ( conj ( gprr ) );
		const complex_type cgprl ( conj ( gprl ) );
		const complex_type cgplr ( conj ( gplr ) );
		const complex_type cgpll ( conj ( gpll ) );

		const complex_type gqrr ( cgains[iqrr] );
		const complex_type gqrl ( pkg->cgains[iqrl] );
		const complex_type gqlr ( pkg->cgains[iqlr] );
		const complex_type gqll ( pkg->cgains[iqll] );

		const complex_type cgqrr ( conj ( gqrr ) );
		const complex_type cgqrl ( conj ( gqrl ) );
		const complex_type cgqlr ( conj ( gqlr ) );
		const complex_type cgqll ( conj ( gqll ) );

		// the following long expressions come from sympy
		// see :math_gradient_full_unity.py:
		// see :math_gradient_full_unity.code:
	
const complex_type Dgprr_coeff_gprr = MI2*(gqlr*cgqlr + gqrr*cgqrr) ; 
const complex_type Dgprr_coeff_dpqrr = -MI*gqrr*za ; 
const complex_type Dgprr_coeff_dpqrl = -MI*gqlr*za ; 
const complex_type Dgprr_coeff_gprl = MI2*za2*(gqlr*cgqll + gqrr*cgqrl) ; 

const complex_type Dgqrr_coeff_gqrr = MI2*(gplr*cgplr + gprr*cgprr) ; 
const complex_type Dgqrr_coeff_dqplr = -MI*gplr*zb ; 
const complex_type Dgqrr_coeff_dqprr = -MI*gprr*zb ; 
const complex_type Dgqrr_coeff_gqrl = MI2*zb2*(gplr*cgpll + gprr*cgprl) ; 

		grad[iprr] += Dgprr_coeff_gprr*gprr + Dgprr_coeff_dpqrr*dpqrr + Dgprr_coeff_dpqrl*dpqrl + Dgprr_coeff_gprl*gprl;

		grad[iqrr] += Dgqrr_coeff_gqrr*gqrr + Dgqrr_coeff_dqplr*dqplr + Dgqrr_coeff_dqprr*dqprr + Dgqrr_coeff_gqrl*gqrl;

		// update cost
		// only rr, rl, lr
		cost += std::norm ( dpqrr - MI*gprl*za*cgqrl - MI*gprr*zb*cgqrr );
		cost += std::norm ( dpqrl - MI*gprl*za*cgqll - MI*gprr*zb*cgqlr );
		cost += std::norm ( dpqlr - MI*gpll*za*cgqrl - MI*gplr*zb*cgqrr );
		//cost += std::norm ( dpqll - MI*gpll*za*cgqll - MI*gplr*zb*cgqlr );

	} // iterate over baselines
	
	/* load complex grad into real and imaginary parts */
	for ( int igain = 0; igain < ngains; igain++ ) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain + 0] = gg.real();
		rgrad[2*igain + 1] = gg.imag();
	}

	return cost;
}

LBFGS::real_type LBFGS::unpolarized_rl (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * Solving for rl gain, populates the derivative of rl that appears in {rr,rl,lr} baselines.
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = reinterpret_cast<const uata_t*>(instance);

	/* model I */
	const real_type MI ( pkg->modeli );
	const real_type MI2 ( MI * MI );

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

		const complex_type dqprr ( conj ( dpqrr ) );
		const complex_type dqprl ( conj ( dpqrl ) );
		const complex_type dqplr ( conj ( dpqlr ) );
		const complex_type dqpll ( conj ( dpqll ) );


		// fetch the parallactic angle skip
		const complex_type za ( pkg->par_za[ibl] );
		const complex_type zb ( pkg->par_zb[ibl] );

		const complex_type za2 ( za * za );
		const complex_type zb2 ( zb * zb );

		// set the gain indices
		const int iprr ( 4*iant1 );
		const int iprl ( iant1 );
		const int iplr ( 4*iant1 + 2 );
		const int ipll ( 4*iant1 + 3 );

		const int iqrr ( 4*iant2 );
		const int iqrl ( iant2 );
		const int iqlr ( 4*iant2 + 2 );
		const int iqll ( 4*iant2 + 3 );

		// fetch full gains for both antennas
		const complex_type gprr ( pkg->cgains[iprr] );
		const complex_type gprl ( cgains[iprl] );
		const complex_type gplr ( pkg->cgains[iplr] );
		const complex_type gpll ( pkg->cgains[ipll] );

		const complex_type cgprr ( conj ( gprr ) );
		const complex_type cgprl ( conj ( gprl ) );
		const complex_type cgplr ( conj ( gplr ) );
		const complex_type cgpll ( conj ( gpll ) );

		const complex_type gqrr ( pkg->cgains[iqrr] );
		const complex_type gqrl ( cgains[iqrl] );
		const complex_type gqlr ( pkg->cgains[iqlr] );
		const complex_type gqll ( pkg->cgains[iqll] );

		const complex_type cgqrr ( conj ( gqrr ) );
		const complex_type cgqrl ( conj ( gqrl ) );
		const complex_type cgqlr ( conj ( gqlr ) );
		const complex_type cgqll ( conj ( gqll ) );

		// the following long expressions come from sympy
		// see :math_gradient_full_unity.py:
		// see :math_gradient_full_unity.code:
	
const complex_type Dgprl_coeff_gprl = MI2*(gqll*cgqll + gqrl*cgqrl) ; 
const complex_type Dgprl_coeff_dpqrr = -MI*gqrl*zb ; 
const complex_type Dgprl_coeff_dpqrl = -MI*gqll*zb ; 
const complex_type Dgprl_coeff_gprr = MI2*zb2*(gqll*cgqlr + gqrl*cgqrr) ; 

const complex_type Dgqrl_coeff_gqrl = MI2*(gpll*cgpll + gprl*cgprl) ; 
const complex_type Dgqrl_coeff_dqplr = -MI*gpll*za ; 
const complex_type Dgqrl_coeff_dqprr = -MI*gprl*za ; 
const complex_type Dgqrl_coeff_gqrr = MI2*za2*(gpll*cgplr + gprl*cgprr) ; 

		grad[iprl] += Dgprl_coeff_gprl*gprl + Dgprl_coeff_dpqrr*dpqrr + Dgprl_coeff_dpqrl*dpqrl + Dgprl_coeff_gprr*gprr ;

		grad[iqrl] += Dgqrl_coeff_gqrl*gqrl + Dgqrl_coeff_dqplr*dqplr + Dgqrl_coeff_dqprr*dqprr + Dgqrl_coeff_gqrr*gqrr ;

		// update cost
		// only rr, rl, lr
		cost += std::norm ( dpqrr - MI*gprl*za*cgqrl - MI*gprr*zb*cgqrr );
		cost += std::norm ( dpqrl - MI*gprl*za*cgqll - MI*gprr*zb*cgqlr );
		cost += std::norm ( dpqlr - MI*gpll*za*cgqrl - MI*gplr*zb*cgqrr );
		//cost += std::norm ( dpqll - MI*gpll*za*cgqll - MI*gplr*zb*cgqlr );

	} // iterate over baselines
	
	/* load complex grad into real and imaginary parts */
	for ( int igain = 0; igain < ngains; igain++ ) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain + 0] = gg.real();
		rgrad[2*igain + 1] = gg.imag();
	}

	return cost;
}

LBFGS::real_type LBFGS::unpolarized_lr (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * Solving for rl gain, populates the derivative of rl that appears in {rr,rl,lr} baselines.
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = reinterpret_cast<const uata_t*>(instance);

	/* model I */
	const real_type MI ( pkg->modeli );
	const real_type MI2 ( MI * MI );

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

		const complex_type dqprr ( conj ( dpqrr ) );
		const complex_type dqprl ( conj ( dpqrl ) );
		const complex_type dqplr ( conj ( dpqlr ) );
		const complex_type dqpll ( conj ( dpqll ) );


		// fetch the parallactic angle skip
		const complex_type za ( pkg->par_za[ibl] );
		const complex_type zb ( pkg->par_zb[ibl] );

		const complex_type za2 ( za * za );
		const complex_type zb2 ( zb * zb );

		// set the gain indices
		const int iprr ( 4*iant1 );
		const int iprl ( 4*iant1 + 1 );
		const int iplr ( iant1 );
		const int ipll ( 4*iant1 + 3 );

		const int iqrr ( 4*iant2 );
		const int iqrl ( 4*iant2 + 1 );
		const int iqlr ( iant2 );
		const int iqll ( 4*iant2 + 3 );

		// fetch full gains for both antennas
		const complex_type gprr ( pkg->cgains[iprr] );
		const complex_type gprl ( pkg->cgains[iprl] );
		const complex_type gplr ( cgains[iplr] );
		const complex_type gpll ( pkg->cgains[ipll] );

		const complex_type cgprr ( conj ( gprr ) );
		const complex_type cgprl ( conj ( gprl ) );
		const complex_type cgplr ( conj ( gplr ) );
		const complex_type cgpll ( conj ( gpll ) );

		const complex_type gqrr ( pkg->cgains[iqrr] );
		const complex_type gqrl ( pkg->cgains[iqrl] );
		const complex_type gqlr ( cgains[iqlr] );
		const complex_type gqll ( pkg->cgains[iqll] );

		const complex_type cgqrr ( conj ( gqrr ) );
		const complex_type cgqrl ( conj ( gqrl ) );
		const complex_type cgqlr ( conj ( gqlr ) );
		const complex_type cgqll ( conj ( gqll ) );

		// the following long expressions come from sympy
		// see :math_gradient_full_unity.py:
		// see :math_gradient_full_unity.code:
		//
const complex_type Dgplr_coeff_gplr = MI2*(gqlr*cgqlr + gqrr*cgqrr) ; 
const complex_type Dgplr_coeff_dpqlr = -MI*gqrr*za ; 
const complex_type Dgplr_coeff_dpqll = -MI*gqlr*za ; 
const complex_type Dgplr_coeff_gpll = MI2*za2*(gqlr*cgqll + gqrr*cgqrl) ; 

const complex_type Dgqlr_coeff_gqlr = MI2*(gplr*cgplr + gprr*cgprr) ; 
const complex_type Dgqlr_coeff_dqpll = -MI*gplr*zb ; 
const complex_type Dgqlr_coeff_dqprl = -MI*gprr*zb ; 
const complex_type Dgqlr_coeff_gqll = MI2*zb2*(gplr*cgpll + gprr*cgprl) ; 

		grad[iplr] += Dgplr_coeff_gplr*gplr + Dgplr_coeff_dpqlr*dpqlr + Dgplr_coeff_dpqll*dpqll + Dgplr_coeff_gpll*gpll ;

		grad[iqlr] += Dgqlr_coeff_gqlr*gqlr + Dgqlr_coeff_dqpll*dqpll + Dgqlr_coeff_dqprl*dqprl + Dgqlr_coeff_gqll*gqll ;
	
		// update cost
		// only rl, lr, ll
		//cost += std::norm ( dpqrr - MI*gprl*za*cgqrl - MI*gprr*zb*cgqrr );
		cost += std::norm ( dpqrl - MI*gprl*za*cgqll - MI*gprr*zb*cgqlr );
		cost += std::norm ( dpqlr - MI*gpll*za*cgqrl - MI*gplr*zb*cgqrr );
		cost += std::norm ( dpqll - MI*gpll*za*cgqll - MI*gplr*zb*cgqlr );

	} // iterate over baselines
	
	/* load complex grad into real and imaginary parts */
	for ( int igain = 0; igain < ngains; igain++ ) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain + 0] = gg.real();
		rgrad[2*igain + 1] = gg.imag();
	}

	return cost;
}

LBFGS::real_type LBFGS::unpolarized_ll (void *instance, const lbfgsfloatval_t *rgains, lbfgsfloatval_t *rgrad, const int n, const lbfgsfloatval_t step) {
	/*
	 * Solving for rl gain, populates the derivative of rl that appears in {rr,rl,lr} baselines.
	*/

	/* return this */
	real_type cost ( 0.0f );

	/* get data_t* ptr out of instance */
	const uata_t *pkg = reinterpret_cast<const uata_t*>(instance);

	/* model I */
	const real_type MI ( pkg->modeli );
	const real_type MI2 ( MI * MI );

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

		const complex_type dqprr ( conj ( dpqrr ) );
		const complex_type dqprl ( conj ( dpqrl ) );
		const complex_type dqplr ( conj ( dpqlr ) );
		const complex_type dqpll ( conj ( dpqll ) );


		// fetch the parallactic angle skip
		const complex_type za ( pkg->par_za[ibl] );
		const complex_type zb ( pkg->par_zb[ibl] );

		const complex_type za2 ( za * za );
		const complex_type zb2 ( zb * zb );

		// set the gain indices
		const int iprr ( 4*iant1 );
		const int iprl ( 4*iant1 + 1 );
		const int iplr ( 4*iant1 + 2 );
		const int ipll ( iant1 );

		const int iqrr ( 4*iant2 );
		const int iqrl ( 4*iant2 + 1 );
		const int iqlr ( 4*iant2 + 2 );
		const int iqll ( iant2 );

		// fetch full gains for both antennas
		const complex_type gprr ( pkg->cgains[iprr] );
		const complex_type gprl ( pkg->cgains[iprl] );
		const complex_type gplr ( pkg->cgains[iplr] );
		const complex_type gpll ( cgains[ipll] );

		const complex_type cgprr ( conj ( gprr ) );
		const complex_type cgprl ( conj ( gprl ) );
		const complex_type cgplr ( conj ( gplr ) );
		const complex_type cgpll ( conj ( gpll ) );

		const complex_type gqrr ( pkg->cgains[iqrr] );
		const complex_type gqrl ( pkg->cgains[iqrl] );
		const complex_type gqlr ( pkg->cgains[iqlr] );
		const complex_type gqll ( cgains[iqll] );

		const complex_type cgqrr ( conj ( gqrr ) );
		const complex_type cgqrl ( conj ( gqrl ) );
		const complex_type cgqlr ( conj ( gqlr ) );
		const complex_type cgqll ( conj ( gqll ) );

		// the following long expressions come from sympy
		// see :math_gradient_full_unity.py:
		// see :math_gradient_full_unity.code:
	
const complex_type Dgpll_coeff_gpll = MI2*(gqll*cgqll + gqrl*cgqrl) ; 
const complex_type Dgpll_coeff_dpqlr = -MI*gqrl*zb ; 
const complex_type Dgpll_coeff_dpqll = -MI*gqll*zb ; 
const complex_type Dgpll_coeff_gplr = MI2*zb2*(gqll*cgqlr + gqrl*cgqrr) ; 

const complex_type Dgqll_coeff_gqll = MI2*(gpll*cgpll + gprl*cgprl) ; 
const complex_type Dgqll_coeff_dqpll = -MI*gpll*za ; 
const complex_type Dgqll_coeff_dqprl = -MI*gprl*za ; 
const complex_type Dgqll_coeff_gqlr = MI2*za2*(gpll*cgplr + gprl*cgprr) ; 

		grad[ipll] += Dgpll_coeff_gpll*gpll + Dgpll_coeff_dpqlr*dpqlr + Dgpll_coeff_dpqll*dpqll + Dgpll_coeff_gplr*gplr ;

		grad[iqll] += Dgqll_coeff_gqll*gqll + Dgqll_coeff_dqpll*dqpll + Dgqll_coeff_dqprl*dqprl + Dgqll_coeff_gqlr*gqlr ;

		// update cost
		// only rr, rl, lr
		//cost += std::norm ( dpqrr - MI*gprl*za*cgqrl - MI*gprr*zb*cgqrr );
		cost += std::norm ( dpqrl - MI*gprl*za*cgqll - MI*gprr*zb*cgqlr );
		cost += std::norm ( dpqlr - MI*gpll*za*cgqrl - MI*gplr*zb*cgqrr );
		cost += std::norm ( dpqll - MI*gpll*za*cgqll - MI*gplr*zb*cgqlr );

	} // iterate over baselines
	
	/* load complex grad into real and imaginary parts */
	for ( int igain = 0; igain < ngains; igain++ ) {
		const complex_type gg ( grad[igain] );
		rgrad[2*igain + 0] = gg.real();
		rgrad[2*igain + 1] = gg.imag();
	}

	return cost;
}

LBFGS::cost_t LBFGS::IterativeUnpolarizedSolver::solve (uata_t& pkg, int n ) {

	cost_t    final_cost;

	void* vpkg  = static_cast<void*>(&pkg);

	for (int __i = 0; __i < n; __i++) {
		/* load gains from pkg into xpars */
		loadgains ( pkg );

		/* solve for rr and ll */
		rcode_rr = lbfgs(n_rr, xpar_rr, &cost_rr, unpolarized_rr, NULL, vpkg, &param);
		rcode_ll = lbfgs(n_ll, xpar_ll, &cost_ll, unpolarized_ll, NULL, vpkg, &param);

		/* save solved gains into pkg.cgains */
		savegains_rrll ( pkg );

		/* solve for rl and lr */
		rcode_rl = lbfgs(n_rl, xpar_rl, &cost_rl, unpolarized_rl, NULL, vpkg, &param);
		rcode_lr = lbfgs(n_lr, xpar_lr, &cost_lr, unpolarized_lr, NULL, vpkg, &param);

		/* save solved gains into pkg.cgains */
		savegains_rllr ( pkg );

		final_cost[0]  = cost_rr;
		final_cost[1]  = cost_rl;
		final_cost[2]  = cost_lr;
		final_cost[3]  = cost_ll;

#ifdef DPRINT
	std::cout << "iter=" << __i << " costs rr=" << cost_rr << " rl=" << cost_rl << " lr=" << cost_lr << " ll=" << cost_ll << std::endl;
#endif
		
	} // solving iterations

  return final_cost;
}
