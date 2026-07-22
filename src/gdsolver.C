#include "gdsolver.hpp"


/**
 * variable names are same as in :gradient_descent.py:
 *
 * usol is vector<complex<float>> ( 2*nant )
 * g_nr and g_dr also same size
 **/
GDSolver::real_type GDSolver::iterate(const ptrdata_t& pkg, vc_type& usol) {

	std::fill ( g_nr.begin(), g_nr.end(), complex_type(0.0f, 0.0f) );
	std::fill ( g_dr.begin(), g_dr.end(), 0.0f );

	real_type cost (0.0f);

	for ( int im = 0; im < m; im++ ) {

		const complex_type& data  = pkg->data [ im ];
		const complex_type& model = pkg->model[ im ];

		const int& b1    = pkg->index_b1 [ im ];
		const int& b2    = pkg->index_b2 [ im ];

		const complex_type& g1 = usol [ b1 ];
		const complex_type& g2 = usol [ b2 ];

		/* compute cost */
		cost += std::norm ( data - ( g1 * model * std::conj(g2) ) );

		const real_type    d_dr ( std::norm(model) * std::norm ( g2 ) );
		const complex_type d_nr ( data * std::conj ( model ) * g2 );

		const real_type    i_dr ( std::norm(model) * std::norm ( g1 ) );
		const complex_type i_nr ( std::conj(data) * std::conj ( model ) * g1 );

		g_nr[b1] += d_nr;
		g_dr[b1] += d_dr;

		g_nr[b2] += i_nr;
		g_dr[b2] += i_dr;

	} // polar baselines
	
	// update usol in place
	// we do lerp with alpha
	for ( int ia = 0; ia < n; ia++ ) {

		const complex_type gia ( g_nr[ia] / g_dr[ia] );

		usol[ia] = (1.0f - alpha)*usol[ia] + alpha*gia;
	} // vsol
	
	return cost;

} // iteration


GDSolver::real_type GDSolver::solve ( const ptrdata_t& pkg, vc_type& solution ) {

	rcode = 0;

	real_type last_cost ( 0.0f );

	for (int iter = 0; iter < max_iterations; iter++ ){

		// iterate once
		real_type cost = iterate ( pkg, solution );

		// note that this cost is one behind the update
		if ( last_cost - cost <= delta ) {
			rcode = 1;
			break;
		}

		last_cost = cost;

	} // iterations


	return last_cost;
}
