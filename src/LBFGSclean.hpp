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
		 	 	
		 	 for (int iant = 0; iant < nantennas; iant++) {
		 	 	 cgains[4*iant + 0] = complex_type(1.0f, 0.0f);
		 	 	 cgains[4*iant + 1] = complex_type(0.0f, 0.0f);
		 	 	 cgains[4*iant + 2] = complex_type(0.0f, 0.0f);
		 	 	 cgains[4*iant + 3] = complex_type(1.0f, 0.0f);
			 } // initialize cgains
		 }
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

		/*
		* Leakage solutions
		*/
		vc_type                  pleakrl;
		vc_type                  pleaklr;
		vc_type                  qleakrl;
		vc_type                  qleaklr;

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

		diagleak_data_t ( 
				int npbl, int nant
				) : npolarbaselines(npbl), nantennas(nant),
			iant1(npbl), iant2(npbl), pb2corr(npbl), data(npbl),
			pleakrl (npbl), pleaklr(npbl),
			qleakrl (npbl), qleaklr(npbl),
			par_model_rr (npbl), par_model_rl (npbl),
			par_model_lr (npbl), par_model_ll (npbl), 
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
	
	/*
	 * Objective functions of IterativeSolving
	*/
	lbfgsfloatval_t rr_unpolarized (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	lbfgsfloatval_t ll_unpolarized (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	lbfgsfloatval_t rllr_unpolarized (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

	//lbfgsfloatval_t diag_jones (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

	int progress_reporter (void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls );
	int progress_reporter_unity (void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls );

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
			param.max_linesearch    = 128;

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

		int                  rcode;
		lbfgsfloatval_t       cost;
		lbfgsfloatval_t      gnorm;
		lbfgsfloatval_t      *xpar_full;
		int                  niter;

		// ctor
		UnpolarizedSolver (int _nant, int n_hessian_corrections = 16, int max_iterations = 1000) : 
			nantennas(_nant), n_para ( nantennas*2*2 ), n_full(nantennas*4*2),
			niter_para(0), niter(0) 
		{
			/* load default first*/
			lbfgs_parameter_init(&param);

			/* gnorm */
			param.epsilon           = 1E-9;
			/* fnorm */
			param.delta             = 1E-9;
			/* xtol */
			param.xtol              = 1E-9;

			/* number of corrections to the hessian matrix */
			param.m                 = n_hessian_corrections;
			/* max iterations */
			param.max_iterations    = max_iterations;
			/* max linesearch */
			param.max_linesearch    = 128;
			/* linesearch algo */
			param.linesearch        = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
			/* C||x|| */
			//param.orthantwise_c     = 0.1;
			//param.orthantwise_start = 0;
			//param.orthantwise_end   = n_full;

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
		real_type solve_parallel_unpolarized_single (uata_t& pkg);
		real_type solve_parallel_unpolarized (uata_t& pkg);
		real_type solve_leakage_unpolarized (uata_t& pkg);
		real_type solve_full_unpolarized (uata_t& pkg);

	}; // solver
		 //
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

		/* all parameters of the solver */
		lbfgs_parameter_t    param;

		/* number of variables/parameters */
		const int             nantennas;
		const int             n_rr, n_ll;
		const int             n_rllr;

		/* xpar, grad */
		lbfgsfloatval_t      *xpar_rr;
		lbfgsfloatval_t      cost_rr;
		lbfgsfloatval_t      gnorm_rr;
		lbfgsfloatval_t      rcode_rr;
		int                  niter_rr;

		lbfgsfloatval_t      *xpar_ll;
		lbfgsfloatval_t      cost_ll;
		lbfgsfloatval_t      gnorm_ll;
		lbfgsfloatval_t      rcode_ll;
		int                  niter_ll;

		lbfgsfloatval_t      *xpar_rllr;
		lbfgsfloatval_t      cost_rllr;
		lbfgsfloatval_t      gnorm_rllr;
		lbfgsfloatval_t      rcode_rllr;
		int                  niter_rllr;

		/*
		 * Enforce the phases of the reference antenna to be zero for rr and ll gains
		 * imaginary part is zero and real part is purely positive.
		 * The real part positive enforcement is done with parameterizing exp(..)
		 * This condition is termed as GREF.
		 *
		 * For each {rr,ll}, there are completely 2*2*nantennas parameters. 
		 * For each {rr,ll}, there are {parallel gain, leakage gain}, both of which are complex.
		 * parallel gain of the reference antenna is zero.
		 * free parameters = 2*2*nantennas - 1
		 *
		*/

		// ctor
		IterativeUnpolarizedSolver (int _nant, int n_hessian_corrections = 200, int max_iterations = 1000) : 
			nantennas(_nant), n_rr ( nantennas*2*2 - 1 ), n_ll ( nantennas*2*2 - 1 ), n_rllr(nantennas*2*2),
			niter_rr(0), niter_ll(0), niter_rllr(0) 
		{
			/* load default first*/
			lbfgs_parameter_init(&param);

			/* number of corrections to the hessian matrix */
			param.m                 = n_hessian_corrections;
			/* max iterations */
			param.max_iterations    = max_iterations;
			/* gnorm */
			param.epsilon           = 1E-9;
			/* fnorm */
			param.delta             = 1E-9;
			/* xtol */
			param.xtol              = 1E-9;
			/* max linesearch */
			param.max_linesearch    = 128;
			/* linesearch algo */
			param.linesearch        = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;

			xpar_rr   = lbfgs_malloc ( n_rr );
			xpar_ll   = lbfgs_malloc ( n_ll );
			xpar_rllr = lbfgs_malloc ( n_rllr );

			/* initialize xpar */
			std::fill ( xpar_rr, xpar_rr + n_rr, 0.0f ); 
			std::fill ( xpar_ll, xpar_ll + n_ll, 0.0f ); 
			std::fill ( xpar_rllr, xpar_rllr + n_rllr, 0.0f ); 
		}

		void loadgains_rllr ( const uata_t &pkg ) {
			/*
			 * Load unity_data_t::cgains into xpar_rllr
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _grl ( pkg.cgains[4*iant + 1] );
				const complex_type _glr ( pkg.cgains[4*iant + 2] );

				xpar_rllr[4*iant + 0] = _grl.real();
				xpar_rllr[4*iant + 1] = _grl.imag();

				xpar_rllr[4*iant + 2] = _glr.real();
				xpar_rllr[4*iant + 3] = _glr.imag();

			}
		}

		void loadgains_rr ( const uata_t &pkg ) {
			/*
			 * Load unity_data_t::cgains into xpar_rr
			 *
			 * xpar_rr only has grr and drl terms for every antenna
			*/

			{
				// GREF
				const complex_type _grr ( pkg.cgains[0] );
				const complex_type _grl ( pkg.cgains[1] );

				xpar_rr[0]  = std::log (_grr.real());
				xpar_rr[1]  = _grl.real();
				xpar_rr[2]  = _grl.imag();
			}

			for (int iant = 1; iant < nantennas; iant++) {

				const complex_type _grr ( pkg.cgains[4*iant + 0] );
				const complex_type _grl ( pkg.cgains[4*iant + 1] );

				xpar_rr[4*iant - 1] = _grr.real();
				xpar_rr[4*iant + 0] = _grr.imag();
				xpar_rr[4*iant + 1] = _grl.real();
				xpar_rr[4*iant + 2] = _grl.imag();
			}
		}

		void loadgains_ll ( const uata_t &pkg ) {
			/*
			 * Load unity_data_t::cgains into xpar_ll
			 *
			 * xpar_ll only has gll and dlr terms for every antenna
			*/

			{
				// GREF
				const complex_type _glr ( pkg.cgains[2] );
				const complex_type _gll ( pkg.cgains[3] );

				xpar_ll[0]  = _glr.real();
				xpar_ll[1]  = _glr.imag();
				xpar_ll[2]  = std::log(_gll.real());
			}

			for (int iant = 1; iant < nantennas; iant++) {

				const complex_type _glr ( pkg.cgains[4*iant + 2] );
				const complex_type _gll ( pkg.cgains[4*iant + 3] );

				xpar_ll[4*iant - 1] = _glr.real();
				xpar_ll[4*iant + 0] = _glr.imag();
				xpar_ll[4*iant + 1] = _gll.real();
				xpar_ll[4*iant + 2] = _gll.imag();
			}
		}

		void savegains_rllr ( uata_t &pkg ) {
			/*
			 * Save xpar_full into unity_data_t::cgains
			*/

			for (int iant = 0; iant < nantennas; iant++) {

				const complex_type _grl ( xpar_rllr[4*iant + 0], xpar_rllr[4*iant + 1] );
				const complex_type _glr ( xpar_rllr[4*iant + 2], xpar_rllr[4*iant + 3] );

				pkg.cgains[4*iant + 1] = _grl;
				pkg.cgains[4*iant + 2] = _glr;
			}
		}

		void savegains_rr ( uata_t &pkg ) {
			/*
			 * Save xpar_rr into unity_data_t::cgains
			 *
			 * xpar_rr only has grr and drl terms for every antenna
			*/

			{
				// GREF
				const complex_type _grr ( std::exp(xpar_rr[0]), 0.0f );
				const complex_type _grl ( xpar_rr[1], xpar_rr[2] );

				pkg.cgains[0] = _grr;
				pkg.cgains[1] = _grl;
			}

			for (int iant = 1; iant < nantennas; iant++) {

				const complex_type _grr ( xpar_rr[4*iant - 1], xpar_rr[4*iant + 0] );
				const complex_type _grl ( xpar_rr[4*iant + 1], xpar_rr[4*iant + 2] );

				pkg.cgains[4*iant + 0] = _grr;
				pkg.cgains[4*iant + 1] = _grl;
			}
		}

		void savegains_ll ( uata_t &pkg ) {
			/*
			 * Save xpar_ll into unity_data_t::cgains
			 *
			 * xpar_ll only has gll and dlr terms for every antenna
			*/

			{
				// GREF
				const complex_type _glr ( xpar_ll[0], xpar_ll[1] );
				const complex_type _gll ( std::exp(xpar_ll[2]), 0.0f );

				pkg.cgains[2] = _glr;
				pkg.cgains[3] = _gll;
			}

			for (int iant = 1; iant < nantennas; iant++) {

				const complex_type _glr ( xpar_ll[4*iant - 1], xpar_ll[4*iant + 0] );
				const complex_type _gll ( xpar_ll[4*iant + 1], xpar_ll[4*iant + 2] );

				pkg.cgains[4*iant + 2] = _glr;
				pkg.cgains[4*iant + 3] = _gll;
			}
		}

		// dtor
		~IterativeUnpolarizedSolver() {

			if (xpar_rr) lbfgs_free ( xpar_rr );
			if (xpar_ll) lbfgs_free ( xpar_ll );
			if (xpar_rllr) lbfgs_free ( xpar_rllr );
		}

		/* solving methods */
		real_type solve (uata_t& pkg);

	}; // solver
	 
	struct PolarizedSolver {
		/*
		 * - Solves for diag_complex gains and crosshand phase
		 *   We need to solve for both, because we consider leakages as well,
		 *   which makes things non-communitative
		 *
		 *   n = (2*2*nantennas) + 1
		 *   each antenna has two gains and each gain is complex
		 *   add one because of crosshand phase
		 *   initialize crosshand phase from 0
		*/

		/* all parameters of the solver */
		lbfgs_parameter_t    param;

		/* number of variables/parameters */
		const int             nantennas;
		const int             n;
		const int             pindex;

		/* xpar, grad */
		int                  rcode;
		lbfgsfloatval_t       cost;
		lbfgsfloatval_t      gnorm;
		lbfgsfloatval_t      *xpar;
		int                  niter;

		// ctor
		PolarizedSolver (int _nant, int n_hessian_corrections = 32, int max_iterations = 1000) : 
			nantennas(_nant), n( nantennas*2*2 - 1 ), pindex(4*nantennas - 2),
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
			xpar              = lbfgs_malloc ( n );
			std::fill ( xpar, xpar + n, 0.0f ); 
		}

		void initialize () {
			/*
			 * xpar has gains and crosshand phase
			 *
			 * All complex gains are initialized to (1,0)
			 * crosshand phase is initialized to 0.0
			 *
			 * GGREF makes the layout
			 * R R RI RI RI .... RI | xphase
			 *
			*/
			xpar[0] = 1.0f;
			xpar[1] = 1.0f;
			for ( int iant = 1; iant < nantennas; iant++) xpar [4*iant - 2] = 1.0f;
			xpar[pindex] = 0.0f;
		}

		// dtor
		~PolarizedSolver() {

			if (xpar) lbfgs_free ( xpar );

		}

		/* solving methods */
		real_type solve_diag_polarized (data_t& pkg);

	}; // solver

}; // namespace
