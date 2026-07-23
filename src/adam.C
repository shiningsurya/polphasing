#include "adam.hpp"


/*
 * Update the oldgains in place.
*/
int Adam::operator() (const float oldcost, const float newcost, vc_t& oldgains, const vc_t& newgains ) {

	using std::conj;
	using std::norm;
	using std::sqrt;

	/* compute gradient */
	complex_t   gradient;

	const float diff_cost ( newcost - oldcost );

	for ( int ipar = 0; ipar < npar; ipar++ ) {

		/*
		 * grad = Delta Cost / conjugate(Delta gain)
		 * The conjugation here is important for correct updation later.
		*/
		gradient  = diff_cost / conj ( newgains[ipar] - oldgains[ipar] );

		/* momemtum update */
		const complex_t mt = (beta1 * last_mt[ipar]) + ((1.0f - beta1)*gradient);

		/* RMSprop update */
		const real_t    vt = (beta2 * last_vt[ipar]) + ((1.0f - beta2)*norm(gradient));

		/* correct the bias */
		const complex_t hmt = mt / ( 1.0f - rbeta1 ); 
		const real_t    hvt = vt / ( 1.0f - rbeta2 ); 

		/* update the gain */
		const real_t    _dr_term ( eps + sqrt(hvt) );
		const complex_t update_term = hmt * alpha / _dr_term;
		oldgains[ipar] += update_term;

		/* record mt, vt */
		last_mt[ipar]   = mt;
		last_vt[ipar]   = vt;

		/* save running product */
		rbeta1 *= beta1;
		rbeta2 *= beta2;
	} // gains
	
	return 0;
}
