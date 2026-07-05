#pragma once
/*
 * Everything comes here
 *
 * class interface?
 */

#include <vector>
#include <string>
#include <fstream>

#include "objectives.hpp"

/* not to be confused with namespace::polphasing */
class PolPhasing {
	using real_type    = float; 
	using complex_type = std::complex<real_type>; 
	using vi_type      = std::vector<int>;
	using vf_type      = std::vector<real_type>;
	using vc_type      = std::vector<complex_type>;

	using data_t = polphasing_data_t;
	

	private:
		/* have it here */
		int nbaselines, nchannels;
		int ndata;
		int npar;

		/* data baselines,channels  */
		vc_type    data;
		vc_type    model;

		/* Stokes model channels  */
		vf_type    Imodel, Qmodel, Umodel;


		int solve_chan ( int );

	public:

		int solve ();

		int write_solutions ( const std::string& );

};
