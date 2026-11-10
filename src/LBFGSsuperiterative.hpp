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

	struct unity_data_t {
		/*
		* This fitting is done per channel. 
		*/

		/*
		* These are brightness matrix elements.
		* That come from the model.
		*/
		const int                npolarbaselines;
		const int                nantennas;

		/*
		* Indices of antenna1 and antenna2.
		*/
		vi_type                  iant1;
		vi_type                  iant2;
 
		/*
		 * correlation type.
		 * rr:0, rl:1, lr:2, ll:3
		*/
		vi_type                  pb2corr;

		/*
		* Complex data
		*/
		vc_type                  data;

		/*
		 * full gains
		 * size=4*nantennas
		 *
		 * leakage gains
		 * size=2*nantennas
		*/
		vc_type                  cgains;
		vc_type                  lgains;

		/*
		* These are parallactic jones elements
		*/
		vc_type                  par_z1;
		vc_type                  par_z2;

		/* model only total intensity */
		real_type                modeli;

		/* state */
		real_type                cost;
		real_type                gnorm;
		int                      niter;

		unity_data_t ( 
				int npbl, int nant
				) : npolarbaselines(npbl), nantennas(nant),
			iant1(npbl), iant2(npbl), pb2corr(npbl), data(npbl),
			cgains (4*nantennas, complex_type(0.0f, 0.0f)),
			lgains (2*nantennas, complex_type(0.0f, 0.0f)),
			par_z1(npbl), par_z2(npbl), 
		 cost (0.0f), gnorm(0.0f), niter(0) {
		 	 for (int ia = 0; ia < nantennas; ia++) {
		 	 	 cgains[4*ia + 0]  = complex_type (1.0, 0.0);
		 	 	 cgains[4*ia + 3]  = complex_type (1.0, 0.0);
			 } 
		 }
	}; 

	using uata_t     = struct unity_data_t;

	/*
	 * Objective functions
	 * compute cost and gradient in one pass.
	 * instance is a pointer to data_t
	 *
	*/
	lbfgsfloatval_t unpolarized_rr (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	lbfgsfloatval_t unpolarized_rl (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	lbfgsfloatval_t unpolarized_lr (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	lbfgsfloatval_t unpolarized_ll (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

	int progress_reporter (void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls );

	/*
	 * Solver interfaces
	*/
	using cost_t  = std::array<real_type,4>; 
	
	struct IterativeUnpolarizedSolver {
		/*
		 * Solves 
		 * 	(a) Srr and Sll independently
		 * 	(b) Srl and Slr jointly
		 * Repeats until convergence or until two or three iterations
		 *
		 * At each step, unity_data_t::cgains are updated
		 *
		 * GREF condition is manually ensured by SingleUnpolarizedSolver
		 *
		 * For each antenna, 
		 * xpar layout is 
		 * R I R I R I R I
		 * 0 1 2 3 4 5 6 7
		 * rr |rl |lr |ll
		 *
		*/
		/* number of variables/parameters */
		const int             nantennas;
		const int             n_rr; 
		const int             n_rl; 
		const int             n_lr; 
		const int             n_ll; 

		/* all parameters of the solver */
		lbfgs_parameter_t    param;


		/* xpar, grad */
		lbfgsfloatval_t      *xpar_rr;
		lbfgsfloatval_t       cost_rr;
		lbfgsfloatval_t      gnorm_rr;
		lbfgsfloatval_t      rcode_rr;
		int                  niter_rr;

		lbfgsfloatval_t      *xpar_rl;
		lbfgsfloatval_t       cost_rl;
		lbfgsfloatval_t      gnorm_rl;
		lbfgsfloatval_t      rcode_rl;
		int                  niter_rl;

		lbfgsfloatval_t      *xpar_lr;
		lbfgsfloatval_t       cost_lr;
		lbfgsfloatval_t      gnorm_lr;
		lbfgsfloatval_t      rcode_lr;
		int                  niter_lr;

		lbfgsfloatval_t      *xpar_ll;
		lbfgsfloatval_t       cost_ll;
		lbfgsfloatval_t      gnorm_ll;
		lbfgsfloatval_t      rcode_ll;
		int                  niter_ll;

		// ctor
		IterativeUnpolarizedSolver (int _nant, int n_hessian_corrections = 32, int max_iterations = 1000) : 
			nantennas(_nant), 
			n_rr ( nantennas*2 ), n_rl ( nantennas*2 ), n_lr ( nantennas*2 ), n_ll ( nantennas*2 ),
			niter_rr(0), niter_ll(0), niter_rl(0), niter_lr(0)
		{
			/* load default first*/
			lbfgs_parameter_init(&param);

			/* number of corrections to the hessian matrix */
			param.m                 = n_hessian_corrections;
			/* max iterations */
			param.max_iterations    = max_iterations;
			/* max linesearch */
			param.max_linesearch    = 64;

			xpar_rr   = lbfgs_malloc ( n_rr );
			xpar_rl   = lbfgs_malloc ( n_rl );
			xpar_lr   = lbfgs_malloc ( n_lr );
			xpar_ll   = lbfgs_malloc ( n_ll );

			/* initialize xpar */
			std::fill ( xpar_rr, xpar_rr + n_rr, 0.0f ); 
			std::fill ( xpar_rl, xpar_rl + n_rl, 0.0f ); 
			std::fill ( xpar_lr, xpar_lr + n_lr, 0.0f ); 
			std::fill ( xpar_ll, xpar_ll + n_ll, 0.0f ); 
		}

		void loadgains_rr ( const uata_t &pkg ) {
			/*
			 * Load unity_data_t::cgains into xpar_rr
			 *
			 * xpar_rr only has grr for every antenna
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _grr ( pkg.cgains[4*iant + 0] );

				xpar_rr[2*iant + 0] = _grr.real();
				xpar_rr[2*iant + 1] = _grr.imag();
			}
		}

		void loadgains_rl ( const uata_t &pkg ) {
			/*
			 * Load unity_data_t::cgains into xpar_rl
			 *
			 * xpar_rl only has grl for every antenna
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _grl ( pkg.cgains[4*iant + 1] );

				xpar_rl[2*iant + 0] = _grl.real();
				xpar_rl[2*iant + 1] = _grl.imag();
			}
		}

		void loadgains_lr ( const uata_t &pkg ) {
			/*
			 * Load unity_data_t::cgains into xpar_lr
			 *
			 * xpar_rl only has grl for every antenna
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _glr ( pkg.cgains[4*iant + 2] );

				xpar_lr[2*iant + 0] = _glr.real();
				xpar_lr[2*iant + 1] = _glr.imag();
			}
		}

		void loadgains_ll ( const uata_t &pkg ) {
			/*
			 * Load unity_data_t::cgains into xpar_ll
			 *
			 * xpar_ll only has gll for every antenna
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _gll ( pkg.cgains[4*iant + 3] );

				xpar_ll[2*iant + 0] = _gll.real();
				xpar_ll[2*iant + 1] = _gll.imag();
			}
		}

		void savegains_rr ( uata_t &pkg ) {
			/*
			 * Save xpar_rr into unity_data_t::cgains
			 *
			 * xpar_rr only has grr and drl terms for every antenna
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _grr ( xpar_rr[2*iant + 0], xpar_rr[2*iant + 1] );

				pkg.cgains[4*iant + 0] = _grr;
			}
		}

		void savegains_rl ( uata_t &pkg ) {
			/*
			 * Save xpar_rl into unity_data_t::cgains
			 *
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _grl ( xpar_rl[2*iant + 0], xpar_rl[2*iant + 1] );

				pkg.cgains[4*iant + 1] = _grl;
			}
		}

		void savegains_lr ( uata_t &pkg ) {
			/*
			 * Save xpar_lr into unity_data_t::cgains
			 *
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _glr ( xpar_lr[2*iant + 0], xpar_lr[2*iant + 1] );

				pkg.cgains[4*iant + 2] = _glr;
			}
		}

		void savegains_ll ( uata_t &pkg ) {
			/*
			 * Save xpar_ll into unity_data_t::cgains
			 *
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _gll ( xpar_ll[2*iant + 0], xpar_ll[2*iant + 1] );

				pkg.cgains[4*iant + 3] = _gll;
			}
		}

		// dtor
		~IterativeUnpolarizedSolver() {

			if (xpar_rr) lbfgs_free ( xpar_rr );
			if (xpar_rl) lbfgs_free ( xpar_rl );
			if (xpar_lr) lbfgs_free ( xpar_lr );
			if (xpar_ll) lbfgs_free ( xpar_ll );
		}

		/* solving methods */
		cost_t solve (uata_t& pkg, int n = 3);

	}; // solver
	 
}; // namespace
