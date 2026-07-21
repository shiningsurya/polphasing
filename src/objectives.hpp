#pragma once
/*
 * Polarimetric phasing objective functions
 */

#include <vector>
#include <array>
#include <complex>

namespace polphasing {
	using real_type    = float; 
	using complex_type = std::complex<real_type>; 
	using vi_type      = std::vector<int>;
	using vf_type      = std::vector<real_type>;
	using vc_type      = std::vector<complex_type>;

	/*
	 * Single channel data package
	 * We pass this as first argument to fcn
	 *
	 * gains are solve parameters so they come from fcn.
	 * We need model and indices mapping.
	 *
	 * Which residual (m) maps to which (n)
	 */
	struct polphasing_data_t {
		/* complex arrays */
		vc_type data;
		vc_type model;

		/* index mapping */
		vi_type index_b1;
		vi_type index_b2;

		/* polarbaseline to correlation */
		// FITIQU
		// vi_type pb2corr;

		polphasing_data_t ( int ndata ) : 
			data(ndata), model(ndata), 
			index_b1(ndata), index_b2(ndata)  
			// FITIQU
			// pb2corr(ndata) 
		{}

	};



	/* given RR RL LR LL */
	int full_polar_fcn (void*, int, int, const real_type*, real_type*, real_type*, int, int);

	/* when only RR and LL */
	// XXX need to implement
	//int parallel_fcn (void*, int, int, const real_type*, real_type*, real_type*, int, int);

}; /* namespace polphasing */

