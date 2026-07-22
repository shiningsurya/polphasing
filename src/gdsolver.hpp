#pragma once
/*
 * Gradient descent solving
 *
 * similar to rantsol
 *
*/

#include <iostream>
#include <numeric>
#include <algorithm>
#include <array>
#include <vector>
#include <memory>

#include "objectives.hpp"

class GDSolver {
	public:
		using real_type  = float;
		using complex_type = std::complex<real_type>;
		using vr_type    = std::vector<real_type>;
		using vc_type    = std::vector<complex_type>;
		using data_t     = polphasing::data_t;
		using ptrdata_t  = polphasing::ptrdata_t;

	private:
		/*
		 * Interpolation between old solution and new solution
		*/
		static constexpr real_type alpha = 0.40;
		/*
		 * Change in SSE observed
		*/
		static constexpr real_type delta = 0.1;

		//static constexpr complex_type zero_complex = complex_type( 0.0f, 0.0f );

		static constexpr int max_iterations = 1000;

		/* one iteration */
		real_type iterate(const ptrdata_t& pkg, vc_type& usol);


		int m; // number of polar baselines
		int n; // number of antbands
		
		vc_type g_nr;
		vr_type g_dr;
	
	public:
		int rcode; 
		int niter;

		GDSolver ( int _m, int _n ) : m (_m), n (_n), g_nr (n), g_dr(n), rcode(-1), niter(0) {}

		real_type solve ( const ptrdata_t& pkg, vc_type& initial_solution );


};
