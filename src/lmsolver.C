#include "lmsolver.hpp"

template<>
int LMSolver<LMSolverType::FULL_POLAR>::solve ( const ptrdata_t& pkg ) {

	/* call lmder */
	info = __cminpack_func__ (lmder) (
		polphasing::full_polar_fcn,
		static_cast<void*>(pkg.get()),
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

#if 0
	/* repeat the same after manual reset */
	reset ();
	/* using the found solution as initial solution */

	info = __cminpack_func__ (lmder) (
		polphasing::full_polar_fcn,
		static_cast<void*>(pkg.get()),
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
#endif

	return info;
}

template<>
int LMSolver<LMSolverType::FULL_POLAR>::forward ( const ptrdata_t& pkg, int i ) {

	/* call objective function */
	return polphasing::full_polar_fcn (
		static_cast<void*>(pkg.get()),
		m, n, 
		isolution.data(),
		residuals.data(),
		jacobian.data(),
		ldfjac,
		i
	);
}
