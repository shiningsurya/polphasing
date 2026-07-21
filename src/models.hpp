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


}; /* models namespace */
