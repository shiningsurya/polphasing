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
			par_model_rr (npbl), par_model_rl (npbl), par_model_lr (npbl), par_model_ll (npbl) {}
	}; 

	using data_t       = struct data_t;

	/*
	 * Objective functions
	 * compute cost and gradient in one pass.
	 * instance is a pointer to data_t
	*/
	lbfgsfloatval_t full_jones (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

	lbfgsfloatval_t diag_jones (void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

	int progress_reporter (void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls );

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
			/* linesearch strong wolfe */
			// param.linesearch        = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;

			xpar      = lbfgs_malloc ( npar );
			grad      = lbfgs_malloc ( npar );

			/* initialize xpar */
			std::fill ( xpar, xpar + npar, 0.0f ); 
			for ( int ipar = 0; ipar < npar; ipar+=8 ) xpar[ipar] = 1.0f;
			for ( int ipar = 6; ipar < npar; ipar+=8 ) xpar[ipar] = 1.0f;
		}

		// dtor
		~Solver() {

			if (xpar) lbfgs_free ( xpar );
			if (grad) lbfgs_free ( grad );
		}

		//template<typename jones_t>
		// for now only implement full jones
		real_type operator()(data_t& pkg);

	}; // solver

}; // namespace
