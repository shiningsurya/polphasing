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
		const complex_type       mrr, mrl, mlr, mll;

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
		* These are brightness matrix elements that have
		* been parallactic angle rotation applied.
		* The P-jones have been applied.
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
				int npbl, int nant,
				const complex_type _mrr,
				const complex_type _mrl,
				const complex_type _mlr,
				const complex_type _mll
				) : npolarbaselines(npbl), nantennas(nant),
			iant1(npbl), iant2(npbl), pb2corr(npbl), data(npbl),
			mrr(_mrr), mrl(_mrl), mlr(_mlr), mll(_mll),
			par_model_rr (npbl), par_model_rl (npbl), par_model_lr (npbl), par_model_ll (npbl),
		 cost (0.0f), gnorm(0.0f), niter(0) {}
	}; 

	struct combined_data_t {
		/*
		* This fitting is done per channel. 
		*
		* This struct contains both unpol (uol) and pol (pol)
		* datasets
		*
		* For unpol, we also fit for I. 
		* Because, we can. 
		* Because we are not always gauranteed fluxmodels of unpolarized calibrators
		* Our flux reference is from polarized calibrator, because we absolutely need IQU model
		*
		* Unless explicitly mentioned, assume the variables are for pol
		*
		* For unpolarized source, since we are also fitting for I,
		* we are conveniently ignoring the parang correction that is needed.
		*/

		/*
		* These are brightness matrix elements.
		* That come from the model.
		*/
		const int                npolarbaselines;
		const int                nantennas;
		const complex_type       polmrr, polmrl, polmlr, polmll;
		/* these must correspond to pol */

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
		vc_type                  pol_data;
		vc_type                  uol_data;

		/*
		* These are brightness matrix elements that have
		* been parallactic angle rotation applied.
		* The P-jones have been applied.
		*/
		vc_type                  polpar_model_rr;
		vc_type                  polpar_model_rl;
		vc_type                  polpar_model_lr;
		vc_type                  polpar_model_ll;

		/* state */
		real_type                cost;
		real_type                gnorm;
		int                      niter;

		combined_data_t ( 
				int npbl, int nant,
				const complex_type _mrr,
				const complex_type _mrl,
				const complex_type _mlr,
				const complex_type _mll
				) : npolarbaselines(npbl), nantennas(nant),
			iant1(npbl), iant2(npbl), pb2corr(npbl), pol_data(npbl), uol_data(npbl),
			polmrr(_mrr), polmrl(_mrl), polmlr(_mlr), polmll(_mll),
			polpar_model_rr (npbl), polpar_model_rl (npbl), polpar_model_lr (npbl), polpar_model_ll (npbl),
		 cost (0.0f), gnorm(0.0f), niter(0) {}
	}; 

	using data_t       = struct data_t;
	using c_data_t     = struct combined_data_t;

	/*
	 * Objective functions
	 * compute cost and gradient in one pass.
	 * instance is a pointer to data_t
	*/
	lbfgsfloatval_t full_jones (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

	lbfgsfloatval_t combined_jones (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

	lbfgsfloatval_t diag_jones (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

	int progress_reporter (void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls );
	int c_progress_reporter (void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls );

	/*
	 * Solver interface
	*/
	struct Solver {

		/* all parameters of the solver */
		lbfgs_parameter_t    param;

		/* number of variables/parameters */
		const int             npar;

		/* xpar, grad */
		int                  rcode;
		lbfgsfloatval_t       cost;
		lbfgsfloatval_t      gnorm;
		lbfgsfloatval_t      *xpar;
		lbfgsfloatval_t      *grad;
		int                  niter;

		// ctor
		Solver (int _npar, int n_hessian_corrections = 16, int max_iterations = 1000) : 
			npar(_npar),
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

			xpar      = lbfgs_malloc ( npar );
			grad      = lbfgs_malloc ( npar );

			/* initialize xpar */
			std::fill ( xpar, xpar + npar, 0.0f ); 
		}

		/* initializes xpar */
		void initialize_full_jones () {
			/*
			 * Only the parallel gains are set to unity with zero imaginary.
			*/
			for ( int ipar = 0; ipar < npar; ipar+=8 ) xpar[ipar] = 1.0f;
			for ( int ipar = 6; ipar < npar; ipar+=8 ) xpar[ipar] = 1.0f;
		}
		void initialize_diag_jones () {
			/*
			 * Only the parallel gains are set to unity with zero imaginary.
			 * Which in case of diag_jones, is every gain
			 *
			 * R I R I R I R I
			 *
			 * In case of GREF, the layout is
			 * R R I R I R I
			*/
#ifdef GREF
			xpar[0] = 1.0f;
			for ( int ipar = 1; ipar < npar; ipar+=2 ) xpar[ipar] = 1.0f;
#else
			for ( int ipar = 0; ipar < npar; ipar+=2 ) xpar[ipar] = 1.0f;
#endif
		}

		// dtor
		~Solver() {

			if (xpar) lbfgs_free ( xpar );
			if (grad) lbfgs_free ( grad );
		}

		/* solving methods */
		real_type solve_full_jones (data_t& pkg);
		real_type solve_full_jones (c_data_t& pkg );
		real_type solve_diag_jones (data_t& pkg);

	}; // solver

}; // namespace
