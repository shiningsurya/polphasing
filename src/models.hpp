#pragma once

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <complex>

/* fmt library */
//#define FMT_HEADER_ONLY
//#include "fmt/base.h"

/*
 * All interfaces to the models
 *
 * should i put them in a namespace?
*/

namespace models {
	using real_type    = float;
	using complex_type = std::complex<real_type>;
	using vf_type      = std::vector<real_type>;
	using vc_type      = std::vector<complex_type>;
	using model_type   = std::array<complex_type,4>;

	struct model_data_t {
		std::string name;

		real_type  fedge;
		real_type  fbw;
		int        nchans;

		vf_type  freqs;
		vf_type  stokes_i;
		vf_type  stokes_q;
		vf_type  stokes_u;

		model_data_t ( int nchan, real_type _fbw, real_type _fedge ) : 
			nchans ( nchan ), fbw ( _fbw ), fedge (_fedge),
			freqs (nchan),
			stokes_i (nchan), stokes_q (nchan), stokes_u (nchan) {}
	};

	struct model_vis_t {
		int        nchans;

		vc_type  rr;
		vc_type  rl;
		vc_type  lr;
		vc_type  ll;

		model_vis_t ( int nchan) : 
			nchans ( nchan ),
			rr (nchan), rl (nchan), lr (nchan), ll(nchan) {}
	};

	model_data_t read_model_file ( const std::string& str );
	/*
	 * 
	 * Model file format
	 *
	 * freq         stokes_i     stokes_q    stokes_u
	 * 9.6          7.4          7.4         7.4
	 *
	*/

	int write_model_file ( const model_data_t& model, const std::string& );

	int write_model_file ( const model_vis_t& model, const std::string& );

	model_type parallactic_correction ( const real_type, const real_type, const complex_type, const complex_type, const complex_type, const complex_type );

}; /* models namespace */

std::ostream& operator<< (std::ostream& os, const models::complex_type& g );
