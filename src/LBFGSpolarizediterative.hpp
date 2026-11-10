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
		*/
		vc_type                  par_model_rr;
		vc_type                  par_model_rl;
		vc_type                  par_model_lr;
		vc_type                  par_model_ll;

		/* state */
		real_type                cost;
		real_type                gnorm;
		int                      niter;

		data_t ( 
				int npbl, int nant
				) : npolarbaselines(npbl), nantennas(nant),
			iant1(npbl), iant2(npbl), pb2corr(npbl), data(npbl),
			cgains(4*nantennas, complex_type(0.0, 0.0)),
			leakrl (nantennas), leaklr(nantennas),
			par_model_rr (npbl), par_model_rl (npbl),
			par_model_lr (npbl), par_model_ll (npbl), 
		 cost (0.0f), gnorm(0.0f), niter(0) {

		 	 for (int ia = 0; ia < nantennas; ia++) {
		 	 	 cgains[4*ia  ] = complex_type(1.0, 0.0);
		 	 	 cgains[4*ia+3] = complex_type(1.0, 0.0);
			 }
		 }
	}; 

	/*
	 * Objective functions
	 * compute cost and gradient in one pass.
	 * instance is a pointer to data_t
	 *
	*/
	lbfgsfloatval_t diagpolarized_rr (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	lbfgsfloatval_t diagpolarized_ll (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

	lbfgsfloatval_t polarized_rr (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	lbfgsfloatval_t polarized_rl (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	lbfgsfloatval_t polarized_lr (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	lbfgsfloatval_t polarized_ll (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

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

		// ctor
		IterativePolarizedDiagonalSolver (int _nant, int n_hessian_corrections = 32, int max_iterations = 1000) : 
			nantennas(_nant), 
			n_rr ( nantennas*2 ),  n_ll ( nantennas*2 ),
			niter_rr(0), niter_ll(0)
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
			xpar_ll   = lbfgs_malloc ( n_ll );

			/* initialize xpar */
			std::fill ( xpar_rr, xpar_rr + n_rr, 0.0f ); 
			std::fill ( xpar_ll, xpar_ll + n_ll, 0.0f ); 
		}

		void loadgains ( const data_t &pkg ) {
			/*
			 * Load unity_data_t::cgains into xpar_rr
			 *
			 * xpar_rr only has grr for every antenna
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _grr ( pkg.cgains[4*iant + 0] );
				const complex_type _gll ( pkg.cgains[4*iant + 3] );

				xpar_rr[2*iant + 0] = _grr.real();
				xpar_rr[2*iant + 1] = _grr.imag();

				xpar_ll[2*iant + 0] = _gll.real();
				xpar_ll[2*iant + 1] = _gll.imag();
			}
		}

		void savegains ( data_t &pkg ) {
			/*
			 * Save xpar_rr into unity_data_t::cgains
			 *
			 * xpar_rr only has grr and drl terms for every antenna
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _grr ( xpar_rr[2*iant + 0], xpar_rr[2*iant + 1] );
				const complex_type _gll ( xpar_ll[2*iant + 0], xpar_ll[2*iant + 1] );

				pkg.cgains[4*iant + 0] = _grr;
				pkg.cgains[4*iant + 3] = _gll;
			}
		}

		// dtor
		~IterativePolarizedDiagonalSolver() {

			if (xpar_rr) lbfgs_free ( xpar_rr );
			if (xpar_ll) lbfgs_free ( xpar_ll );
		}

		/* solving methods */
		diagcost_t solve (data_t& pkg, int n = 3);

	}; // solver
		 //
	struct IterativePolarizedSolver {
		/*
		 * Solves 
		 * 	(a) Grr and Gll 
		 * 	(b) Grl and Glr
		 * Repeats until convergence or until two or three iterations
		 * We need to repeat it because mrl and mlr terms contain expressions of grr * conj(gll) and conjugate
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
		IterativePolarizedSolver (int _nant, int n_hessian_corrections = 32, int max_iterations = 1000) : 
			nantennas(_nant), 
			n_rr ( nantennas*2 ),  
			n_rl ( nantennas*2 ),  
			n_lr ( nantennas*2 ),
			n_ll ( nantennas*2 ),
			niter_rr(0), 
			niter_rl(0), 
			niter_lr(0), 
			niter_ll(0)
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

		void loadgains ( const data_t &pkg ) {
			/*
			 * Load unity_data_t::cgains into xpar_rr
			 *
			 * xpar_rr only has grr for every antenna
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _grr ( pkg.cgains[4*iant + 0] );
				const complex_type _grl ( pkg.cgains[4*iant + 1] );
				const complex_type _glr ( pkg.cgains[4*iant + 2] );
				const complex_type _gll ( pkg.cgains[4*iant + 3] );

				xpar_rr[2*iant + 0] = _grr.real();
				xpar_rr[2*iant + 1] = _grr.imag();

				xpar_rl[2*iant + 0] = _grl.real();
				xpar_rl[2*iant + 1] = _grl.imag();

				xpar_lr[2*iant + 0] = _glr.real();
				xpar_lr[2*iant + 1] = _glr.imag();

				xpar_ll[2*iant + 0] = _gll.real();
				xpar_ll[2*iant + 1] = _gll.imag();
			}
		}

		void savegains_odg ( data_t &pkg ) {
			/*
			 * Save xpar_rr into unity_data_t::cgains
			 *
			 * xpar_rr only has grr and drl terms for every antenna
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _grl ( xpar_rl[2*iant + 0], xpar_rl[2*iant + 1] );
				const complex_type _glr ( xpar_lr[2*iant + 0], xpar_lr[2*iant + 1] );

				pkg.cgains[4*iant + 1] = _grl;
				pkg.cgains[4*iant + 2] = _glr;
			}
		}

		void savegains_par ( data_t &pkg ) {
			/*
			 * Save xpar_rr into unity_data_t::cgains
			 *
			 * xpar_rr only has grr and drl terms for every antenna
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _grr ( xpar_rr[2*iant + 0], xpar_rr[2*iant + 1] );
				const complex_type _gll ( xpar_ll[2*iant + 0], xpar_ll[2*iant + 1] );

				pkg.cgains[4*iant + 0] = _grr;
				pkg.cgains[4*iant + 3] = _gll;
			}
		}

		// dtor
		~IterativePolarizedSolver() {

			if (xpar_rr) lbfgs_free ( xpar_rr );
			if (xpar_rl) lbfgs_free ( xpar_rl );
			if (xpar_lr) lbfgs_free ( xpar_lr );
			if (xpar_ll) lbfgs_free ( xpar_ll );
		}

		/* solving methods */
		cost_t solve (data_t& pkg, int n = 3);

	}; // solver


} // namespace LBFGS

