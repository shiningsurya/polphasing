#pragma once

//#define DEBUG

#include "lbfgs.h"

#include <complex>
#include <array>
#include <vector>

namespace LBFGS {

	/* data package */
	using real_type    = lbfgsfloatval_t; 
	using complex_type = std::complex<real_type>; 
	using vi_type      = std::vector<int>;
	using vr_type      = std::vector<real_type>;
	using vc_type      = std::vector<complex_type>;

	/* fitting package */
	struct data_t {
		/*
		* This fitting is done per channel. 
		*/

		const int                npolarbaselines;
		const int                nbaselines;
		const int                nantennas;

		/*
		* Indices of antenna1 and antenna2.
		* size of nbaselines
		*/
		vi_type                  iant1;
		vi_type                  iant2;
 
		/*
		 * correlation type.
		 * rr:0, rl:1, lr:2, ll:3
		 * size npolarbaselines
		 * XXX do we still need this?
		*/
		vi_type                  pb2corr;

		/*
		* Complex data
		* size of npolarbaselines
		*/
		vc_type                  data;

		/*
		* gain solutions
		*
		* parallel gains and leakages
		* cgains = 4*nantennas
		*/
		vc_type                  cgains;

		/*
		* Par applied model data
		* size nbaselines
		*/
		vc_type                  par_model_rr;
		vc_type                  par_model_rl;
		vc_type                  par_model_lr;
		vc_type                  par_model_ll;

		data_t ( 
				int npbl, int nant
				) : npolarbaselines(npbl), nbaselines(npbl/4), nantennas(nant),
			iant1(nbaselines), iant2(nbaselines), pb2corr(npbl), data(npbl),
			cgains(4*nantennas, complex_type(0.0, 0.0)),
			par_model_rr (nbaselines), par_model_rl (nbaselines),
			par_model_lr (nbaselines), par_model_ll (nbaselines) {}
	}; 


} // namespace LBFGS


