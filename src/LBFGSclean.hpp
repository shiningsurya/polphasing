#pragma once

#define DEBUG

#include "lbfgs.h"

#include <complex>
#include <vector>


namespace LBFGS {

	/* type of Jones being solved */
	enum class jones_t {
		DIAG_COMPLEX,
		COMPLEX
	};

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
		*/
		vc_type                  cgains;

		/*
		* These are parallactic jones elements
		*/
		vc_type                  par_z1;
		vc_type                  par_z2;

		/* state */
		real_type                cost;
		real_type                gnorm;
		int                      niter;

		unity_data_t ( 
				int npbl, int nant
				) : npolarbaselines(npbl), nantennas(nant),
			iant1(npbl), iant2(npbl), pb2corr(npbl), data(npbl),
			cgains (4*nantennas, complex_type(0.0f, 0.0f)),
			par_z1(npbl), par_z2(npbl), 
		 cost (0.0f), gnorm(0.0f), niter(0) {}
	}; 

	struct diagleak_data_t {
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

		vc_type                  parleak_model_rr;
		vc_type                  parleak_model_rl;
		vc_type                  parleak_model_lr;
		vc_type                  parleak_model_ll;

		/* state */
		real_type                cost;
		real_type                gnorm;
		int                      niter;

		diagleak_data_t ( 
				int npbl, int nant
				) : npolarbaselines(npbl), nantennas(nant),
			iant1(npbl), iant2(npbl), pb2corr(npbl), data(npbl),
			parleak_model_rr (npbl), parleak_model_rl (npbl),
			parleak_model_lr (npbl), parleak_model_ll (npbl), 
		 cost (0.0f), gnorm(0.0f), niter(0) {}
	}; 

	using uata_t     = struct unity_data_t;
	using data_t     = struct diagleak_data_t;

	/*
	 * Objective functions
	 * compute cost and gradient in one pass.
	 * instance is a pointer to data_t
	 *
	 * {parallel,full}_unpolarized
	 * assume source is unpolarized with unity fluxden
	*/
	template<int forcorr>
	lbfgsfloatval_t single_unpolarized (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	lbfgsfloatval_t parallel_unpolarized (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	lbfgsfloatval_t leakage_unpolarized (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	lbfgsfloatval_t full_unpolarized (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

	lbfgsfloatval_t parallel_polarized (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	//lbfgsfloatval_t combined_jones (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

	//lbfgsfloatval_t diag_jones (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

	int progress_reporter (void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls );

	/*
	 * Solver interfaces
	*/
	struct SingleUnpolarizedSolver {
		/*
		 * Unpolarized solving for any single hand
		 *
		 * There is no absolute phase in an interferometer.
		 * So we can arbitrarily set the phase of the gain of the first antenna to zero.
		 * However, instead of setting a constraint while solving,
		 * we enforce this condition after solving.
		 *
		 * Earlier i used to parameterize the gain as exp, so that it is always positive.
		 * But it fails.
		 *
		 *
		*/

		/* all parameters of the solver */
		lbfgs_parameter_t    param;

		/* number of variables/parameters */
		const int             nantennas;
		const int             n;

		/* xpar, grad */
		int                  rcode;
		lbfgsfloatval_t       cost;
		lbfgsfloatval_t      gnorm;
		lbfgsfloatval_t      *xpar;
		int                  niter;

		// ctor
		SingleUnpolarizedSolver (int _nant, int n_hessian_corrections = 32, int max_iterations = 1000) : 
			nantennas(_nant), n ( nantennas*2 ), 
			niter(0) 
		{
			/* load default first*/
			lbfgs_parameter_init(&param);

			/* number of corrections to the hessian matrix */
			param.m                 = n_hessian_corrections;
			/* max iterations */
			param.max_iterations    = max_iterations;
			/* max linesearch */
			param.max_linesearch    = 64;

			/* GREF constraint */
			xpar       = lbfgs_malloc ( n );

			/* initialize xpar */
			std::fill ( xpar, xpar + n, 0.0f ); 
		}

		void initialize () {
			/*
			 * All the real-parts are set to 1.0
			 *
			*/
			for ( int ipar = 0; ipar < n; ipar+=2 ) xpar[ipar] = 1.0f;
		}

		// dtor
		~SingleUnpolarizedSolver() {
			if (xpar) lbfgs_free ( xpar);
		}

		/* solving methods */
		template<int forcorr>
		real_type solve_single_unpolarized (uata_t& pkg);

	}; // solver
		 //
	struct UnpolarizedLeakageSolver {
		/*
		 * Uses SingleUnpolarizedSolver to solve for r and l.
		 * And then solves only for leakage terms.
		 *
		 * GREF condition is ensured by SingleUnpolarizedSolver
		 *
		 * For each antenna, 
		 * we are solving for gains rl and lr.
		 * xpar layout is 
		 * R I R I 
		 * 0 1 2 3 
		 * rl |lr 
		*/

		/* all parameters of the solver */
		lbfgs_parameter_t    param;

		/* number of variables/parameters */
		const int             nantennas;
		const int             n;

		/* xpar, grad */
		int                  rcode;
		lbfgsfloatval_t       cost;
		lbfgsfloatval_t      gnorm;
		lbfgsfloatval_t      *xpar;
		int                  niter;

		// ctor
		UnpolarizedLeakageSolver (int _nant, int n_hessian_corrections = 32, int max_iterations = 1000) : 
			nantennas(_nant), n(nantennas*2*2),
			niter(0) 
		{
			/* load default first*/
			lbfgs_parameter_init(&param);

			/* number of corrections to the hessian matrix */
			param.m                 = n_hessian_corrections;
			/* max iterations */
			param.max_iterations    = max_iterations;
			/* max linesearch */
			param.max_linesearch    = 64;

			/* initialize xpar */
			xpar        = lbfgs_malloc ( n );
			std::fill ( xpar, xpar + n, 0.0f ); 
		}

		// dtor
		~UnpolarizedLeakageSolver() {
			if (xpar) lbfgs_free ( xpar );
		}

		/* solving methods */
		real_type solve_leakage_unpolarized (uata_t& pkg);

	}; // solver
		
	struct UnpolarizedSolver {
		/*
		 * Uses SingleUnpolarizedSolver to solve for r and l.
		 * And then solves for full Jones terms.
		 *
		 * GREF condition is manually ensured by SingleUnpolarizedSolver
		 *
		 * For each antenna, 
		 * xpar_full layout is 
		 * R I R I R I R I
		 * 0 1 2 3 4 5 6 7
		 * rr |rl |lr |ll
		 *
		 * XXX after checking if parallel gains are valid, remove _para
		*/

		/* all parameters of the solver */
		lbfgs_parameter_t    param;

		/* number of variables/parameters */
		const int             nantennas;
		const int             n_para;
		const int             n_full;

		/* xpar, grad */
		int                  rcode_para;
		lbfgsfloatval_t       cost_para;
		lbfgsfloatval_t      gnorm_para;
		lbfgsfloatval_t      *xpar_para;
		int                  niter_para;

		int                  rcode_full;
		lbfgsfloatval_t       cost_full;
		lbfgsfloatval_t      gnorm_full;
		lbfgsfloatval_t      *xpar_full;
		int                  niter_full;

		// ctor
		UnpolarizedSolver (int _nant, int n_hessian_corrections = 32, int max_iterations = 1000) : 
			nantennas(_nant), n_para ( nantennas*2*2 ), n_full(nantennas*4*2),
			niter_para(0), niter_full(0) 
		{
			/* load default first*/
			lbfgs_parameter_init(&param);

			/* number of corrections to the hessian matrix */
			param.m                 = n_hessian_corrections;
			/* max iterations */
			param.max_iterations    = max_iterations;
			/* max linesearch */
			param.max_linesearch    = 64;

			/* GREF constraint */
			xpar_para      = lbfgs_malloc ( n_para );
			xpar_full      = lbfgs_malloc ( n_full );

			/* initialize xpar */
			std::fill ( xpar_para, xpar_para + n_para, 0.0f ); 
			std::fill ( xpar_full, xpar_full + n_full, 0.0f ); 
		}

		void initialize_parallel () {
			/*
			 * Only the parallel gains are set to unity with zero imaginary.
			 * Which in case of diag_jones, is every gain
			 *
			 * Because of GREF, the layout is
			 * R R | R I R I|
			 * ant | ant    |
			 *
			 * GREF in unpolarized case requires us to completely eliminate crosshand phase
			 *
			*/
			//xpar_para[0] = 5.5f;
			//xpar_para[1] = 5.5f;
			for ( int ipar = 0; ipar < n_para; ipar+=2 ) xpar_para[ipar] = 5.0f;
		}

		void initialize_full () {
			/*
			 * initialize using parallel solution
			 *
			 * it is upto the caller to ensure parallel solution is solved
			 *
			 * parallel solution is 
			 * R R | RI RI | 
			 * ant |  ant  |
			 * 0 1 | 2  4  | 6
			 *
			 * 4iant - 2, 4iant -1, 4iant, 4iant+1
			 *
			 * full solution is 
			 * R RI RI R | RI RI RI RI |
			 *  ant      |     ant     |
			 * 0 12 34 5 | 6  8  10 12 | 14
			 *
			 * 8iant-2, 8iant-1, 8iant, 8iant+1, 8iant+2, 8iant+3, 8iant+4, 8iant+5
			 *
			 *  full[0] = para[0]
			 *  full[5] = para[1]
			 *
			 *  full[6,7]   = para[2,3]
			 *  full[12,13] = para[4,5] 
			 *
			*/
			// GREF constraint
			//xpar_full[0] = xpar_para[0];
			//xpar_full[5] = xpar_para[1];

			// loop over antennas
			// starting from 1
			for (int iant = 0; iant < nantennas; iant++) {
				// rr
				xpar_full[8*iant + 0] = xpar_para[4*iant + 0];
				xpar_full[8*iant + 1] = xpar_para[4*iant + 1];

				// ll
				xpar_full[8*iant + 7] = xpar_para[4*iant + 2];
				xpar_full[8*iant + 8] = xpar_para[4*iant + 3];
			}
		}

		// dtor
		~UnpolarizedSolver() {

			if (xpar_para) lbfgs_free ( xpar_para );
			if (xpar_full) lbfgs_free ( xpar_full );

		}

		/* solving methods */
		real_type solve_full_unpolarized (uata_t& pkg);

	}; // solver
	 
	struct PolarizedSolver {
		/*
		 * Internally manages two solvers
		 * - one for parallel solving
		 * - one for full solving
		*/

		/* all parameters of the solver */
		lbfgs_parameter_t    param;

		/* number of variables/parameters */
		const int             nantennas;
		const int             n_para;

		/* xpar, grad */
		int                  rcode_para;
		lbfgsfloatval_t       cost_para;
		lbfgsfloatval_t      gnorm_para;
		lbfgsfloatval_t      *xpar_para;
		int                  niter_para;

		// ctor
		PolarizedSolver (int _nant, int n_hessian_corrections = 32, int max_iterations = 1000) : 
			nantennas(_nant), n_para ( nantennas*2*2 - 1 ),
			niter_para(0)
		{
			/* load default first*/
			lbfgs_parameter_init(&param);

			/* number of corrections to the hessian matrix */
			param.m                 = n_hessian_corrections;
			/* max iterations */
			param.max_iterations    = max_iterations;
			/* max linesearch */
			param.max_linesearch    = 64;

			/* GREF constraint */
			xpar_para      = lbfgs_malloc ( n_para );

			/* initialize xpar */
			std::fill ( xpar_para, xpar_para + n_para, 0.0f ); 
		}

		void initialize_parallel () {
			/*
			 * Only the parallel gains are set to unity with zero imaginary.
			 * Which in case of diag_jones, is every gain
			 *
			 * Because of GREF, the layout is
			 * R R | R I R I|
			 * ant | ant    |
			 *
			 * GREF in unpolarized case requires us to completely eliminate crosshand phase
			 *
			*/
			xpar_para[0] = 1.5f;
			//xpar_para[1] = 5.5f;
			for ( int ipar = 1; ipar < n_para; ipar+=2 ) xpar_para[ipar] = 5.0f;
		}

		// dtor
		~PolarizedSolver() {

			if (xpar_para) lbfgs_free ( xpar_para );

		}

		/* solving methods */
		real_type solve_diag_polarized (data_t& pkg);

	}; // solver

}; // namespace
