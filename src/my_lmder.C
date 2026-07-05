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

//#define __cminpack_float__

#include <iostream>
#include <numeric>
#include <algorithm>
#include <array>
#include <vector>

using std::vector;
using std::array;

using real_type = float;

/* simple test */
constexpr int npar  = 2; /* npar */
constexpr int ndata = 5;  /* ndata */
constexpr array<real_type,ndata> xx {1., 2., 3., 4., 5.};
constexpr array<real_type,ndata> yy {120., 220., 320., 420., 520.};

extern "C" {
#	include "cminpack.h"
}

constexpr int MAX_FUNCTION_EVALUATIONS = 1000;

int fcn (void*, int, int, const real_type*, real_type*, real_type*, int, int);

class Problem {
	public:
		using vr_type = std::vector<real_type>;
	private:
		/* tolerances */
		/* rel. error in sum of squares of errors desired */
		static constexpr real_type ftol = 1E-9;
		/* rel. error in approx. sum of squares of errors desired */
		static constexpr real_type xtol = 1E-9;
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
		/* output flag */ 
		int info;
		/* evaluation counters */
		int nfev, njev;

		/* tuning parameters */
		real_type factor;

		/* input vector (size=n) */
		/* initial solution --> final estimate of solution */
		vr_type isolution;
		/* residuals vector (size=m) */
		vr_type residuals;
		/* jacobian matrix as vector (shape=(m,n)) */
		vr_type jacobian;
		/* diag (size=(n)) */
		vr_type diag;

		/* work arrays */
		/* following must be size=n*/
		std::vector<int> ipvt;
		vr_type qtf;
		vr_type wa1, wa2, wa3;
		/* following must be size=m*/
		vr_type wa4;

		/*XXX why have fcn within the case?*/
		//int fcn (void*, int, int, const real_type*, real_type*, real_type*, int, int);

	public:
		/* ctor */
		Problem (int ndata, int npar, int _nprint = 0, real_type _factor = 100.) : n(npar), m(ndata), ldfjac(m),
			/* work arrays */
			ipvt(n), qtf(n), wa1(n), wa2(n), wa3(n), wa4(m),
			/* main vectors */
			isolution(n), residuals(m), jacobian (m*n),
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

		/* default dtor since nothing is dynamically allocated */
		~Problem() = default;

		/* i do not want to copy */
		Problem (const Problem &other) = delete;

		/* compute sum of squared errors */
		real_type get_sse () const {
			return std::transform_reduce(
					residuals.cbegin(), residuals.cend(),
					0.0,
					std::plus<real_type>(),
					[] (real_type x) { return x*x; }
			);
		}

		/* main method - solve */
		int solve ( const vr_type& data );
		int solve ( const vr_type& data, vr_type& initial_solution );
	
		/* write */
		friend std::ostream& operator<< (std::ostream&, const Problem&);
};

int Problem::solve ( const vr_type& data, vr_type& initial_solution ) {
	/* input size check */
	if ( data.size() != m || initial_solution.size() != n ) {
		info = 0;
		return -1;
	}

	/* copy given initial solution to isolution*/
	std::copy ( initial_solution.cbegin(), initial_solution.cend(), isolution.begin() );

	/* run solve */
	solve ( data );

	/* copy isolution back into initial_solution */
	std::copy ( isolution.cbegin(), isolution.cend(), initial_solution.begin() );

	return info;
}

int Problem::solve ( const vr_type& data ) {
	/* input size check */
	if ( data.size() != m ) {
		std::cout << "bad size" << std::endl;
		info = 0;
		return -1;
	}

	/* call lmder */
	info = __cminpack_func__ (lmder) (
		fcn,
		(void*) data.data(),
		m, n, 
		isolution.data(),
		residuals.data(),
		jacobian.data(),
		ldfjac,
		ftol, xtol, gtol,
		maxfev,
		diag.data(), mode, 
		factor, nprint, 
		&nfev, &njev,
		ipvt.data(), qtf.data(),
		wa1.data(), wa2.data(),
		wa3.data(), wa4.data()
	);

	return info;
}

std::ostream& operator<< (std::ostream& os, const Problem& p) {

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

/* 
 * for the linear model 
 *
 * model ( x | bias, slope ) = bias + slope*x
 * res ( x | bias, slope )   = data - bias - slope*x
 *
 * d res / dbias  = -1.0
 * d res / dslope = -x
 *
 */
int fcn( void* data, int m, int n, const real_type *x, real_type *fvec, real_type *fjac, int ldfjac, int iflag ) {
	/* outside of optimizer */
	/* for now use the constexpr data and model */
	/* m is ndata, n is npar */

	if ( iflag == 0 ) {
		//std::cout << " iflag=0 ";
		//for (int in = 0; in < n; in++) {
			//std::cout << x[in] << " ";
		//}
		//std::cout << std::endl;

	} // printing
	else if ( iflag == 1 ) {

		//std::cout << " iflag=1 fvec=";
		for (int im = 0; im < m; im++) {
			/* model prediction */
			real_type model = x[0] + ( x[1] * xx[im] );
			/* error = DATA - MODEL */
			fvec [ im ] = yy[im] - model;
			//std::cout << "(" << model << "," << yy[im] << ")" << " ";
		} /* every data point */
		//std::cout << std::endl;
		//std::cout << " iflag=1 x=";
		//for (int in = 0; in < n; in++) {
			//std::cout << x[in] << " ";
		//}
		//std::cout << std::endl;

	} /* function evaluation */
	else if ( iflag == 2 ) {
		for (int im = 0; im < m; im++) {
			/* npar, ndata */
			fjac [ ldfjac*0 + im ] = - 1.0;
			fjac [ ldfjac*1 + im ] = - xx[im];
		} /* every data point */
		//for (int in = 0; in < n; in++) {
			//std::cout << x[in] << " ";
		//}
		//std::cout << std::endl;
		
	} /* jacobian evaluation */

	return iflag;
}



int main() {
	Problem test(ndata, npar, 1);
	Problem::vr_type isol ( npar, 2.0 );
	//Problem::vr_type isol {20, 120};
	Problem::vr_type data ( ndata, 1.0 );
	std::cout << test << std::endl;

	//std::cout << "before solving SSE=" << test.get_sse() << std::endl;
	//std::cout << test << std::endl;
	test.solve ( data, isol );
	std::cout << test << std::endl;
	std::cout << "after solving SSE=" << test.get_sse() << std::endl;

	std::cout << isol[0] << " " << isol[1] << std::endl;

	return 0;
}

