#pragma once
/*
 * Gradient descent solving
 * for full Jones matrix
 *
 * four Complex gains per antenna
 *
 * This is run per channel
 * each channel is independent
 *
*/

#include <iostream>
#include <numeric>
#include <memory>
#include <complex>
#include <algorithm>
#include <array>
#include <vector>
#include <map>

#include "adam.hpp"

class FullGDSolver {
	public:
		using real_type    = float; 
		using complex_type = std::complex<real_type>; 
		using vi_type      = std::vector<int>;
		using vr_type      = std::vector<real_type>;
		using vc_type      = std::vector<complex_type>;

		using antname_t    = std::array<char,4>;
		using names_t      = std::vector<antname_t>;

		struct data_t {
			/*
			 * This fitting is done per channel. 
			 *
			 * ant2idx = {ant:idx}
			 * We need to map each antenna to some index
			 * which will help us fetch gains and equation indices
			 *
			 * data  = (ant1, ant2, band1, band2, [pbcorr], complex)
			 * size  = number of polar baselines
			 *
			 * model = (rr, rl, lr, ll)
			 * size  = 4
			 *
			 * gains = for every ant ( grr, grl, glr, gll ) 
			 * 				 each of which is complex
			 * 				 {ant: (grr, grl, glr, gll)}
			 * size  = nantennas * 4
			 *
			 * We get updated gains after solving system of four equations 
			 * in four unknown. (they are 2 in 2, but same thing)
			 * Each coefficient is complex and stored in keyvalue pairs
			 * where keys are antennas and values are 3-array of complex
			 * eqn_drr = {ant: (a1, b1, c1)} 
			 * eqn_drl = {ant: (a2, b2, c2)} 
			 * eqn_dlr = {ant: (a3, b3, c3)} 
			 * eqn_dll = {ant: (a4, b4, c4)} 
			 * size (each) = nantennas*3
			 *
			 * see math_leakages.{pdf,tex,py}
			 * 
			 * a1*gp_rr + b1*gp_rl = c1
			 * a2*gp_rr + b2*gp_rl = c2
			 * a3*gp_lr + b3*gp_ll = c3
			 * a4*gp_lr + b4*gp_ll = c4
			 *
			 * this numbering scheme is comforting
			 *
			 * We solve first two and last two separately.
			 * It makes our lives easier.
			 *
			 * det12, and det34 are the determinants resp.
			 *
			 * | a1  b1 | | rr | = | c1 |
			 * | a2  b2 | | rl |   | c2 |
			 *
			 * | rr | = (det12)**-1 * | b2  -b1 | | c1 | 
			 * | rl |                 | -a2  a1 | | c2 |
			 *
			 * = (det-term-here) | b2*c1 - b1*c2 | 
			 *                   | c2*a1 - c1*a2 |
			 *
			 * | a3  b3 | | lr | = | c3 |
			 * | a4  b4 | | ll |   | c4 |
			 *
			 * | lr | = (det34)**-1 * |  b4 -b3 | | c3 |
			 * | ll |                 | -a4  a3 | | c4 |
			 * = (det-term-here) | b4c3 - b3c4 |
			 *                   | c4a3 - c3a4 |
			 *
			 * implementation details:
			 * instead of an associative container for gains and eqns.
			 * We will use vectors and use ant2idx to fetch.
			 * This would be faster than any tree lookup.
			 */

			const complex_type       mrr, mrl, mlr, mll;

			//names_t                  ant1;
			//names_t                  ant2;
			vi_type                  iant1;
			vi_type                  iant2;

			vi_type                  pb2corr;

			vc_type                  data;

			// parallactic angle corrected model
			// over polarbaselines
			vc_type                  par_model_rr;
			vc_type                  par_model_rl;
			vc_type                  par_model_lr;
			vc_type                  par_model_ll;

			// why keep ant2idx
			//const std::map<antname_t,int>& ant2idx;

			data_t ( 
					int npbl, 
					const complex_type _mrr,
					const complex_type _mrl,
					const complex_type _mlr,
					const complex_type _mll
			) : 
				iant1(npbl), iant2(npbl), pb2corr(npbl), data(npbl),
				mrr(_mrr), mrl(_mrl), mlr(_mlr), mll(_mll),
				par_model_rr (npbl), par_model_rl (npbl), par_model_lr (npbl), par_model_ll (npbl) {}
		}; 

		struct model_t {
			/*
			 * When we are solving for model given gains
			*/

			//names_t                  ant1;
			//names_t                  ant2;
			vi_type                  iant1;
			vi_type                  iant2;

			vi_type                  pb2corr;

			vc_type                  data;

			vc_type                  gains;

			// why keep ant2idx
			//const std::map<antname_t,int>& ant2idx;

			model_t ( 
					int npbl, 
					int ngains
			) : 
				iant1(npbl), iant2(npbl), pb2corr(npbl), data(npbl),
				gains ( ngains ) {}
		}; 

		using solve_data_t    = struct data_t;
		using solve_ptrdata_t = std::unique_ptr<solve_data_t>;
		using solve_model_t   = struct model_t;
	
	private:
		/*
		 * Termination if ema(norm(gradient)) < delta
		*/
		static constexpr real_type delta = 0.1;
		/* EMA beta parameter of norm(gradient) */
		/* Default as Adam */
		static constexpr real_type betag = 0.99;
		/* Fast and slow EMA beta parameter for cost */
		// higher beta fast changing
		static constexpr real_type beta_gnorm_slow = 0.75;
		/*
		 * if the difference between the fast_ema and slow_ema is <= gamma,
		 * terminate
		*/
		static constexpr real_type gamma = 0.0001;
		/*
		 * minimum norm of the gradient vector
		 * This is probably arbitrary
		 * obselete
		*/
		// static constexpr real_type gamma = 100;

		//static constexpr complex_type zero_complex = complex_type( 0.0f, 0.0f );

		static constexpr int max_iterations = 10000;

		// iteration method
		int iterate( const solve_data_t& pkg, vc_type& gains );

	public:

		// gradient method
		int gradient ( const solve_data_t& pkg, const vc_type& gains, vc_type& grad );
		int gradient ( const solve_model_t& pkg, 
				const complex_type& mrr, 
				const complex_type& mrl, 
				const complex_type& mlr, 
				const complex_type& mll, 
				vc_type& grad );
	

		int npolarbaselines;
		int nantennas;
		int ngains;

		int rcode;
		int niter;

		real_type gnorm;

		FullGDSolver( int _npolarbaselines, int _nantennas ) : 
			npolarbaselines(_npolarbaselines),
			nantennas (_nantennas), ngains (4*_nantennas),
			rcode (0), niter(0) {}

		real_type solve ( const solve_data_t& pkg, vc_type& gains );
		real_type solve ( const solve_data_t& pkg, const solve_data_t& qkg, vc_type& gains );

		real_type solve ( const solve_model_t& pkg, complex_type& mrr, complex_type& mrl, complex_type& mlr, complex_type& mll );

		real_type cost ( const solve_data_t& pkg, const vc_type& gains );
		real_type cost ( const solve_model_t& pkg, 
		const complex_type& mrr, const complex_type& mrl, const complex_type& mlr, const complex_type& mll );

		real_type norm ( const vc_type& gains );

}; // FullGDSolve

