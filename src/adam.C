#include "adam.hpp"


/*
 * Update the gains in place.
*/
int Adam::operator() (const vc_t& gradient, vc_t& gains ) {

	using std::conj;
	using std::norm;
	using std::sqrt;

	for ( int ipar = 0; ipar < npar; ipar++ ) {

		/* The conjugation here is required for the math 
		 *
		 * We are already computing the conjugated gradient.
		 * 
		 */
		//const complex_t grad = conj(gradient[ipar]);
		const complex_t grad = gradient[ipar];

		/* momemtum update */
		const complex_t mt = (betam * last_mt[ipar]) + ((1.0f - betam)*grad);

		/* RMSprop update */
		const real_t    vt = (betav * last_vt[ipar]) + ((1.0f - betav)*norm(grad));

		/* correct the bias */
		const complex_t hmt = mt / ( 1.0f - rbetam ); 
		const real_t    hvt = vt / ( 1.0f - rbetav ); 

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
	rbetam *= betam;
	rbetav *= betav;

	return 0;
}
