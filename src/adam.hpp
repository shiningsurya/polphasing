#pragma once

#include <complex>
#include <vector>
#include <cmath>


class Adam {
	using real_t     = float;
	using complex_t  = std::complex<real_t>;
	using vr_t = std::vector<real_t>;
	using vc_t = std::vector<complex_t>;

	// small number
	static constexpr float eps = 1E-6;

	// number of parameters
	const int  npar;

	// hyperparameters
	// alpha <- learning rate
	const float alpha; 
	// beta1 is weight for momentum
	const float beta1;
	// beta2 is weight for variance
	const float beta2;

	// to correct the initialization bias
	// we need to power beta?**iter
	// to optimize the operation
	// we save the running computation
	// saves us from using std::pow
	float rbeta1;
	float rbeta2;
	// these are initialized to one

	// time counter
	int iter;
	const int max_iterations;

	// record mt and vt
	vc_t   last_mt;
	vr_t   last_vt;
	
	public:
		Adam(
			int _npar, 
			float _alpha = 0.01,
			float _beta1 = 0.99,
			float _beta2 = 0.99,
			int _max_iterations = 1000
		) : 
			npar(_npar), 
			alpha(_alpha), beta1(_beta1), beta2(_beta2), 
			iter(0), rbeta1(_beta1), rbeta2(_beta2),
			max_iterations ( _max_iterations ),
			last_mt (_npar, complex_t(0.0f,0.0f)), last_vt (_npar, 0.0f) {}

		int operator() (const float oldcost, const float newcost, vc_t& oldgains, const vc_t& newgains );
		/*
		 * Update the oldgains in place.
		*/
};
