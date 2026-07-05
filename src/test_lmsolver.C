#include "lmsolver.hpp"

int main() {

	/* full polar */
	using FullPolarLMSolver = LMSolver<LMSolverType::FULL_POLAR>;
	/* parallel polar */
	//using FullPolarLMSolver = LMSolver<LMSolver_type::PARALLEL_HANDS>;
	
	
	int nbaselines   = 1012;
	int nchannels    = 2048;
	int nantennas    = 22;
	int ndata        = 2 * nbaselines; /* complex -> real,imag */
	int npar         = 4 * nantennas; /* complex(R) and complex(L) gains */

	int ichan        = 500;
	std::array<std::complex<float>,4> model_chan;
	/* RR RL LR LL */
	model_chan[0]    = std::complex<float>(1.0, 0.0);
	model_chan[1]    = std::complex<float>(0.5, 0.5);
	model_chan[2]    = std::complex<float>(0.5,-0.5);
	model_chan[3]    = std::complex<float>(1.0, 0.0);

	int nbldata      = nbaselines * nchannels;
	std::vector<std::complex<float>>  fulldata ( nbldata );
	
	FullPolarLMSolver::ptrdata_t  pkg ( new FullPolarLMSolver::data_t ( nbaselines ) );

	/* read data */
	{
		std::ifstream ifs( "bl.data", std::ios::binary );
		ifs.read ( reinterpret_cast<char*>(fulldata.data()), nbldata * sizeof(std::complex<float>) );
	}
	/* populate pkg */
	for (int ibl = 0; ibl < nbaselines; ibl++ ) {
		/* baseline,channel layout */
		pkg->data  [ ibl ] = fulldata [ ichan + nchannels*ibl ];

		/* need to check further */
		pkg->model [ ibl ] = model_chan [ ibl % 4 ]; 

		/* indices */
	}

	FullPolarLMSolver          test(ndata, npar, 1);
	FullPolarLMSolver::vr_type isol ( npar, 1.0 );

	test.solve ( pkg, isol );

	std::cout << test << std::endl;
	std::cout << "after solving SSE=" << test.get_sse() << std::endl;

	return 0;
}
