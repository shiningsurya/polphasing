#include "adam.hpp"


/*
 * Update the gains in place.
*/
int Adam::operator() (const vc_t& gradient, vc_t& gains ) {

	using std::conj;
	using std::norm;
	using std::sqrt;

	for ( int ipar = 0; ipar < npar; ipar++ ) {

		/* The conjugation here is required for the math */
		const complex_t grad = conj(gradient[ipar]);

		/* momemtum update */
		const complex_t mt = (beta1 * last_mt[ipar]) + ((1.0f - beta1)*grad);

		/* RMSprop update */
		const real_t    vt = (beta2 * last_vt[ipar]) + ((1.0f - beta2)*norm(grad));

		/* correct the bias */
		const complex_t hmt = mt / ( 1.0f - rbeta1 ); 
		const real_t    hvt = vt / ( 1.0f - rbeta2 ); 

		/* update the gain */
		const real_t    _dr_term ( eps + sqrt(hvt) );
		const complex_t update_term = hmt * alpha / _dr_term;

		/* are we sure this has to be minus? */
		gains[ipar] -= update_term;

		/* record mt, vt */
		last_mt[ipar]   = mt;
		last_vt[ipar]   = vt;

	} // gains
		
	/* save running product */
	/* after iterating over gains */
	rbeta1 *= beta1;
	rbeta2 *= beta2;

	return 0;
}
