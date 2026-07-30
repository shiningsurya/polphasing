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

#include "adam.hpp"

class GDSolver {
	public:
		using real_type    = float;
		using complex_type = std::complex<real_type>;
		using vi_type      = std::vector<int>;
		using vr_type      = std::vector<real_type>;
		using vc_type      = std::vector<complex_type>;

		struct data_t {

			const complex_type       mrr, mrl, mlr, mll;

			/* complex arrays */
			vc_type    data;

			/* index mapping */
			vi_type    iant1;
			vi_type    iant2;

			/* polarbaseline to correlation */
			vi_type    pb2corr;

			// parallactic angle corrected model
			// over polarbaselines
			vc_type    par_model_rr;
			vc_type    par_model_rl;
			vc_type    par_model_lr;
			vc_type    par_model_ll;

			data_t ( int ndata, const complex_type _mrr, const complex_type _mrl, const complex_type _mlr, const complex_type _mll ) : 
				mrr (_mrr), mrl (_mrl), mlr (_mlr), mll (_mll),
				data(ndata), 
				par_model_rr(ndata), par_model_rl(ndata), par_model_lr(ndata), par_model_ll(ndata), 
				iant1 (ndata), iant2 (ndata), pb2corr (ndata)
			{}

		};

	private:
		/*
		 * Termination if ema(norm(gradient)) < delta
		*/
		static constexpr real_type delta = 0.01;
		/* EMA beta parameter of norm(gradient) */
		/* Default as Adam */
		static constexpr real_type betag = 0.95;
		/* Fast and slow EMA beta parameter for cost */
		// higher beta fast changing
		static constexpr real_type beta_cost_fast = 0.9;
		static constexpr real_type beta_cost_slow = 0.6;
		/*
		 * if the difference between the fast_ema and slow_ema is <= gamma,
		 * terminate
		*/
		static constexpr real_type gamma = 0.01;

		static constexpr int max_iterations = 100000;


		int gradient ( const data_t& pkg, const vc_type& gains, vc_type& grad );

	public:
		int npolarbaselines;
		int nantennas;
		int ngains;

		int rcode; 
		int niter;

		GDSolver( int _npolarbaselines, int _nantennas ) : 
			npolarbaselines(_npolarbaselines),
			nantennas (_nantennas), ngains (2*_nantennas),
			rcode (0), niter(0) {}

		real_type solve ( const data_t& pkg, vc_type& gains );

		real_type cost ( const data_t& pkg, const vc_type& gains );

		real_type norm ( const vc_type& gains );

};
