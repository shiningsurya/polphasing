#pragma once
/*
 * My interface to lmder
 *
 * see if i can port example in C
 *
 * Compiling only what i need:
 *
 *  g++ -g -D__cminpack_float__ my_lmder.C lmder.c dpmpar.c enorm.c qrfac.c lmpar.c qrsolv.c  -o mylmder
 *
 *  Idea being i take only these c files to give a solver
 *  use lute to provide LTA file interface
 *  Write the fcn function
 *
 *  Need a csv file interface?
 *  Need a binary serialization?
 *  Maybe just write the structs
 *
 */

#include <iostream>
#include <numeric>
#include <algorithm>
#include <array>
#include <vector>


extern "C" {
#	include "cminpack.h"
}

#include "objectives.hpp"

constexpr int MAX_FUNCTION_EVALUATIONS = 4000;


enum class LMSolverType {
	PARALLEL_HANDS,
	FULL_POLAR
};

template<LMSolverType stype>
class LMSolver {
	public:
		using real_type  = float;
		using vr_type    = std::vector<real_type>;
		using data_t     = polphasing::data_t;
		using ptrdata_t  = polphasing::ptrdata_t;

	private:
		/* tolerances */
		/* rel. error in sum of squares of errors desired */
		static constexpr real_type ftol = 1.49012e-8;
		/* rel. error in approx. sum of squares of errors desired */
		static constexpr real_type xtol = 1.49012e-8;
		/* orthogonality between residuals and columns of jacobian matrix */
		static constexpr real_type gtol = 0;

		/* integer flags/counters */
		int m; /* number of residuals/data points */
		int n; /* number of variables/parameters */
		int iflag;
		/* assert n < m*/
		int ldfjac; /* dimension of jacobian (why isn't it always m?)*/
		/* assert ldfjac >= m */
		/* maximum function evaluations */
		static constexpr int maxfev = MAX_FUNCTION_EVALUATIONS; 
		/* we want lmder to scale parameters internally */
		int mode; 
		/* printing flag */
		int nprint;

		/* tuning parameters */
		real_type factor;

		/* diag (size=(n)) */
		vr_type diag;

		/* work arrays */
		/* following must be size=n*/
		std::vector<int> ipvt;
		vr_type qtf;
		vr_type wa1, wa2, wa3;
		/* following must be size=m*/
		vr_type wa4;

	public:
		/* evaluation counters */
		int nfev, njev;
		/* output flag */ 
		int info;

		/* input vector (size=n) */
		/* initial solution --> final estimate of solution */
		vr_type isolution;
		/* residuals vector (size=m) */
		vr_type residuals;
		/* jacobian matrix as vector (shape=(m,n)) */
		vr_type jacobian;

		/* ctor */
		LMSolver (int ndata, int npar, int _nprint = 0, real_type _factor = 100.) : n(npar), m(ndata), ldfjac(m),
			iflag (0),
			/* work arrays */
			ipvt(n), qtf(n), wa1(n), wa2(n), wa3(n), wa4(m),
			/* main vectors */
			/* isolution one to begin with */
			isolution(n, 1.), residuals(m, 0.), jacobian (m*n, 0.),
			/* rest of the vectors */
			diag(n),
			/* initialize counters */
			nfev(0), njev(0), 
			/* tuning */
			factor ( _factor ), mode (1),
			/* printing */
			nprint ( _nprint ),
			/* initialize output flag */
			info ( -1 ) 
		{}

		/* default dtor */
		~LMSolver () = default;

		/* i do not want to copy */
		LMSolver (const LMSolver &other) = delete;

		/* compute sum of squared errors */
		real_type get_sse () const {
			real_type ret = 0.0;
			for ( const real_type &res : residuals ) ret += ( res * res );
			return ret;
		}
		
		/* reset */
		int reset () {
			/* fill zero instead of creating a new object */

			/* few variables*/
			iflag   = 0;
			//nfev    = 0;
			//njev    = 0;
			info    = -1;

			/* vectors */
			//std::fill ( isolution.begin(), isolution.end(), 1.0f );
			std::fill ( residuals.begin(), residuals.end(), 0.0f );
			std::fill ( jacobian.begin(), jacobian.end(), 0.0f );

			std::fill ( diag.begin(), diag.end(), 0.0f );
			std::fill ( qtf.begin(), qtf.end(), 0.0f );

			std::fill ( wa1.begin(), wa1.end(), 0.0f );
			std::fill ( wa2.begin(), wa2.end(), 0.0f );
			std::fill ( wa3.begin(), wa3.end(), 0.0f );
			std::fill ( wa4.begin(), wa4.end(), 0.0f );

			std::fill ( ipvt.begin(), ipvt.end(), 0 );

			return 0;
		}

		/* main method - solve */
		int solve ( const ptrdata_t& pkg );
		int solve ( const ptrdata_t& pkg, vr_type& initial_solution ) {

			/* copy given initial solution to isolution*/
			std::copy ( initial_solution.cbegin(), initial_solution.cend(), isolution.begin() );

			/* run solve */
			solve ( pkg );

			/* copy isolution back into initial_solution */
			std::copy ( isolution.cbegin(), isolution.cend(), initial_solution.begin() );

			return info;
		}

		/* model pass */
		int forward ( const ptrdata_t& pkg, int iflag );
	
		/* need to define in class inline because of templates */
		friend std::ostream& operator<< (std::ostream& os, const LMSolver<stype>& p) {
			switch ( p.info ) {
				case 0: 
					os << "Improper input parameters" << std::endl;
					break;
				case 1: 
					os << "ftol condition reached" << std::endl;
					break;
				case 2: 
					os << "xtol condition reached" << std::endl;
					break;
				case 3: 
					os << "ftol and xtol conditions reached" << std::endl;
					break;
				case 4: 
					os << "gtol condition reached" << std::endl;
					break;
				case 5: 
					os << "maxfev reached" << std::endl;
					break;
				case 6: 
				case 7: 
				case 8: 
					os << "ftol, xtol or gtol is too small" << std::endl;
					break;
				case -1:
					os << "Something is off" << std::endl;
					break;
			}

			os << " isolution = ";
			for (int ipar = 0; ipar < p.n; ipar++) {
				os << p.isolution[ipar] << " ";
			}
			os <<  std::endl;

			os << " nfev=" << p.nfev << " njev=" << p.njev << std::endl;

			return os;
		}
};
