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
	static constexpr float eps = 1E-7;

	// number of parameters
	const int  npar;

	// hyperparameters
	// alpha <- learning rate
	const float alpha; 
	// beta1 is weight for momentum
	const float betam;
	// beta2 is weight for variance
	const float betav;

	// to correct the initialization bias
	// we need to power beta?**iter
	// to optimize the operation
	// we save the running computation
	// saves us from using std::pow
	float rbetam;
	float rbetav;
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
			float _alpha = 0.50,
			float _betam = 0.99,
			float _betav = 0.99,
			int _max_iterations = 1000
		) : 
			npar(_npar), 
			alpha(_alpha), betam(_betam), betav(_betav), 
			iter(0), rbetam(_betam), rbetav(_betav),
			max_iterations ( _max_iterations ),
			last_mt (_npar, complex_t(0.0f,0.0f)), last_vt (_npar, 0.0f) {}

		int operator() (const vc_t& gradient, vc_t& gains );
		/*
		 * Update the oldgains in place.
		*/
};
