#pragma once

#define DEBUG

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
		* cgains = 2*nantennas
		* leak{rl,lr} = nantennas
		*/
		vc_type                  cgains;
		vc_type                  leakrl;
		vc_type                  leaklr;

		/*
		* Par applied model data
		* size nbaselines
		*/
		vc_type                  par_model_rr;
		vc_type                  par_model_rl;
		vc_type                  par_model_lr;
		vc_type                  par_model_ll;

		/*
		 * Precomputed
		* size of npolarbaselines
		*/
		vc_type                  precompute_prr;
		vc_type                  precompute_pll;
		vc_type                  precompute_qrr;
		vc_type                  precompute_qll;

		/* state */
		real_type                cost;
		real_type                gnorm;
		int                      niter;

		data_t ( 
				int npbl, int nant
				) : npolarbaselines(npbl), nbaselines(npbl/4), nantennas(nant),
			iant1(nbaselines), iant2(nbaselines), pb2corr(npbl), data(npbl),
			cgains(2*nantennas, complex_type(1.0, 0.0)),
			leakrl (nantennas), leaklr(nantennas),
			par_model_rr (nbaselines), par_model_rl (nbaselines),
			par_model_lr (nbaselines), par_model_ll (nbaselines), 
			precompute_prr (npbl), precompute_pll(npbl),
			precompute_qrr (npbl), precompute_qll(npbl),
			cost (0.0f), gnorm(0.0f), niter(0) {}
	}; 

	/*
	 * Objective functions
	 * compute cost and gradient in one pass.
	 * instance is a pointer to data_t
	 *
	*/
	lbfgsfloatval_t diagpolarized_rr (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	lbfgsfloatval_t diagpolarized_ll (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	lbfgsfloatval_t diagpolarized_all(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

	/*
	 * Solver interfaces
	*/
	using diagcost_t  = std::array<real_type,2>; 
	using cost_t      = std::array<real_type,4>; 

	struct IterativePolarizedDiagonalSolver {
		/*
		 * Solves 
		 * 	(a) Srr and Sll independently
		 * Repeats until convergence or until two or three iterations
		 * We need to repeat it because mrl and mlr terms contain expressions of grr * conj(gll) and conjugate
		 *
		*/
		/* number of variables/parameters */
		const int             nantennas;
		const int             n_rr; 
		const int             n_ll; 
		const int            n_all; 

		/* all parameters of the solver */
		lbfgs_parameter_t    param;


		/* xpar, grad */
		lbfgsfloatval_t      *xpar_rr;
		lbfgsfloatval_t       cost_rr;
		lbfgsfloatval_t      gnorm_rr;
		lbfgsfloatval_t      rcode_rr;
		int                  niter_rr;

		lbfgsfloatval_t      *xpar_ll;
		lbfgsfloatval_t       cost_ll;
		lbfgsfloatval_t      gnorm_ll;
		lbfgsfloatval_t      rcode_ll;
		int                  niter_ll;

		lbfgsfloatval_t     *xpar_all;
		lbfgsfloatval_t      cost_all;
		lbfgsfloatval_t     gnorm_all;
		lbfgsfloatval_t     rcode_all;
		int                 niter_all;

		// ctor
		IterativePolarizedDiagonalSolver (int _nant, int n_hessian_corrections = 24, int max_iterations = 1000) : 
			nantennas(_nant), 
			n_rr ( nantennas*2 ),  n_ll ( nantennas*2 ), n_all ( nantennas*4 ),
			niter_rr(0), niter_ll(0), niter_all(0)
		{
			/* load default first*/
			lbfgs_parameter_init(&param);

			/* number of corrections to the hessian matrix */
			param.m                 = n_hessian_corrections;
			/* max iterations */
			param.max_iterations    = max_iterations;
			/* max linesearch */
			param.max_linesearch    = 24;

			xpar_rr   = lbfgs_malloc ( n_rr );
			xpar_ll   = lbfgs_malloc ( n_ll );
			xpar_all  = lbfgs_malloc ( n_all );

			/* initialize xpar */
			std::fill ( xpar_rr, xpar_rr + n_rr, 0.0f ); 
			std::fill ( xpar_ll, xpar_ll + n_ll, 0.0f ); 
			std::fill ( xpar_all, xpar_all + n_all, 0.0f ); 
		}

		void loadgains ( const data_t &pkg ) {
			/*
			 * Load unity_data_t::cgains into xpar_rr
			 *
			 * xpar_rr only has grr for every antenna
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _grr ( pkg.cgains[2*iant + 0] );
				const complex_type _gll ( pkg.cgains[2*iant + 1] );

				xpar_rr[2*iant + 0] = _grr.real();
				xpar_rr[2*iant + 1] = _grr.imag();

				xpar_ll[2*iant + 0] = _gll.real();
				xpar_ll[2*iant + 1] = _gll.imag();
			}
		}

		void savegains ( data_t &pkg ) {
			/*
			 * Save xpar_rr,ll into unity_data_t::cgains
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _grr ( xpar_rr[2*iant + 0], xpar_rr[2*iant + 1] );
				const complex_type _gll ( xpar_ll[2*iant + 0], xpar_ll[2*iant + 1] );

				pkg.cgains[2*iant + 0] = _grr;
				pkg.cgains[2*iant + 1] = _gll;
			}
		}

		// dtor
		~IterativePolarizedDiagonalSolver() {
			if (xpar_rr) lbfgs_free ( xpar_rr );
			if (xpar_ll) lbfgs_free ( xpar_ll );
			if (xpar_all) lbfgs_free ( xpar_all );
		}

		/* solving methods */
		diagcost_t solve (data_t& pkg, int n = 3);
		real_type allsolve (data_t& pkg);

	}; // solver
		 //
	void precompute (data_t* pkg);

} // namespace LBFGS

