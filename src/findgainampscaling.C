#include <iostream>
#include <cmath>

#include <unistd.h>

#include <string>

#include "lta_file.hpp"
#include "gaintable.hpp"
#include "models.hpp"
#include "ants.hpp"

#ifdef CUSTOM
#include "custom_file.hpp"
#endif

#include "measuring_xphase.hpp"

using real_type = LBFGS::real_type;
using complex_type = LBFGS::complex_type;

#ifdef TIMING
#include <chrono>
#endif

//#define CHANDEBUG

double get_median ( const std::vector<double>& vrr ) {
	// copy
	std::vector<double> rr ( vrr );

	// sort
	std::sort ( rr.begin(), rr.end() );

	// take mid
	const int nn ( rr.size () / 2 );

	return rr [ nn ];
}

std::ostream& operator<< (std::ostream& os, const LBFGS::complex_type& g ) {
	/* these flags are set before priting */
	//os << std::fixed << std::setprecision(3) << std::showpos;
	os << g.real() << g.imag() << 'j';
	return os;
}


#ifdef TIMING
auto start  = std::chrono::high_resolution_clock::now();
auto end    = std::chrono::high_resolution_clock::now();
#endif 

void print_help () {
	std::cout << "findgainampscaling" << std::endl;
	std::cout << "  Finds the median scaling that brings unapplied vis(pa) and applied vis(pa) to be in the same range. This scaling is then used to write the gwb solutions" << std::endl;
	std::cout << std::endl;
	std::cout << " findgainampscaling [ARGUMENTS] LTA_FILE" << std::endl;
	std::cout << "    -h Print help" << std::endl;
	std::cout << "    -s <scan> scan number of the LTA file used to generate complex gains" << std::endl;
	std::cout << "    -t <tag> tag/stem of the full (rr,rl,lr,ll) complex gains, " << std::endl;
	std::cout << std::endl;
}

int main(int argc, char *argv[]) {

	/* hello getopt, my old friend */
	int opt;
	int cal_scan_number;
	std::string ofile ("gainampscaling.out");
	std::string tag;
	std::string lta_path;

	if ( argc < 2 ) {
		print_help ();
		exit (EXIT_SUCCESS);
	}

	while ( (opt = getopt ( argc, argv, "hs:t:o:" )) != -1 ) {
		switch (opt) {
			case 'h':
				print_help ();
				exit (EXIT_SUCCESS);
				break;
			case 's':
				cal_scan_number = atoi ( optarg );
				break;
			case 't':
				tag  = optarg;
				break;
			case 'o':
				ofile  = optarg;
		} // switch
	} // getopt
	if ( optind >= argc ) {
		print_help ();
		exit (EXIT_SUCCESS);
	}
	/* lta file*/
	lta_path  = argv[optind];
	optind++;

	/* other files */
	const std::string path_gain_rr    = tag + std::string("_rr.gains");
	const std::string path_gain_rl    = tag + std::string("_rl.gains");
	const std::string path_gain_lr    = tag + std::string("_lr.gains");
	const std::string path_gain_ll    = tag + std::string("_ll.gains");

	std::cout << "[inputs] lta="   << lta_path << std::endl;
	std::cout << "[inputs] tag="   << tag  << std::endl;
	std::cout << "[inputs] outfile="   << ofile << std::endl;
	std::cout << "[inputs] scan="  << cal_scan_number << std::endl;

	/***************************************/
	/*      READ LTA FILE                  */
	/***************************************/
#ifdef CUSTOM
	std::cout << "[customfile] This code is modified to run with customfile input." << std::endl;
	std::cout << "[customfile] Given LTA file is not read." << std::endl;
	custom_file lta_file;
#else
	LTA lta_file ( lta_path );
#endif

	int nchannels    = lta_file.nchannels;
	float fbw        = lta_file.fbw;
	float fedge      = lta_file.fedge;

	const int npolarbaselines (lta_file.nbaselines);
	const int nbaselines ( npolarbaselines / 4 );

	/***************************************/
	/*      READ GAIN FILE                 */
	/***************************************/
	gaintable::gaintable_t g_rr  = gaintable::read_complex_solutions ( path_gain_rr );
	gaintable::gaintable_t g_rl  = gaintable::read_complex_solutions ( path_gain_rl );
	gaintable::gaintable_t g_lr  = gaintable::read_complex_solutions ( path_gain_lr );
	gaintable::gaintable_t g_ll  = gaintable::read_complex_solutions ( path_gain_ll );

	/***************************************/
	/*      AVERAGE LTA SCAN               */
	/***************************************/
	/* maybe directly read complex<float>  */
  /* NOTICE:avgbldata contains self-terms*/
	LTA::vf_type avgbldata ( 2 * npolarbaselines * nchannels, 0.0f );

	const LTA::scan_t cal_scan = lta_file.get_scan ( cal_scan_number );

	lta_file.time_average  ( cal_scan_number, avgbldata );

	/***************************************/
	/*        SOLVER RUN                   */
	/***************************************/

	/* main loop */
	std::cout << " Starting main solving loop" << std::endl;

	auto total_start = std::chrono::high_resolution_clock::now();

	LBFGS::vr_type      frr_un ( nchannels );
	LBFGS::vr_type      frr_ap ( nchannels );

	LBFGS::vr_type      fll_un ( nchannels );
	LBFGS::vr_type      fll_ap ( nchannels );

	LBFGS::vr_type      ratio_rr ( nchannels );
	LBFGS::vr_type      ratio_ll ( nchannels );

	constexpr std::array<char,2> bandmap{'r','l'};
	std::ofstream of ( ofile );

	of << "chan ant1 ant2 correlation unapplied applied" << std::endl;
	of << std::fixed << std::setprecision(3) << std::showpos;
	
#ifdef CHANDEBUG
	/* channel 458 has high loss.*/
	for (int ichan = 501; ichan < 502; ichan++) {
#else
//	#pragma omp parallel for num_threads(4) 
	for (int ichan = 0; ichan < nchannels; ichan++) {
#endif

		LBFGS::real_type   rr_unapplied (0.0);
		LBFGS::real_type   ll_unapplied (0.0);
		LBFGS::real_type   rr_applied (0.0);
		LBFGS::real_type   ll_applied (0.0);

		/* initialize data */
		for ( int ii = 0; ii < npolarbaselines; ii++ ) {

			const int ib ( ii / 4 );

			/* index in (baseline, channel) complex<float> */
			const int _i     = ichan + nchannels * ii;

			/* get baseline object */
			const baseline_t&  _bl = lta_file.baselines [ ii ];

			/* get antennas */
			const auto& ant1  = _bl.ant1;
			const auto& ant2  = _bl.ant2;

			/* get bands */
			const auto& band1 = _bl.band1;
			const auto& band2 = _bl.band2;

			/* ID correlation */
			const int pb2corr = band1*2 + band2;

			if ( pb2corr != 0 && pb2corr != 3 ) continue;

			/* copy data */
			const float _real ( avgbldata[2*_i] );
			const float _imag ( avgbldata[2*_i + 1] );
			const LBFGS::complex_type _data ( _real, _imag );

			/* get gains */
			const LBFGS::complex_type g1rr ( g_rr[ant1][ichan] );
			const LBFGS::complex_type g1rl ( g_rl[ant1][ichan] );
			const LBFGS::complex_type g1lr ( g_lr[ant1][ichan] );
			const LBFGS::complex_type g1ll ( g_ll[ant1][ichan] );

			const LBFGS::complex_type g2rr ( g_rr[ant2][ichan] );
			const LBFGS::complex_type g2rl ( g_rl[ant2][ichan] );
			const LBFGS::complex_type g2lr ( g_lr[ant2][ichan] );
			const LBFGS::complex_type g2ll ( g_ll[ant2][ichan] );

			const LBFGS::complex_type det1 ( 1.0 - (g1rl/g1rr)*(g1lr/g1ll) );
			const LBFGS::complex_type det2 ( 1.0 - (g2rl/g2rr)*(g2lr/g2ll) );

			LBFGS::complex_type _iata;
			if ( pb2corr == 0 ) {
					_iata  =  _data / g1rr / std::conj(g2rr) / det1 / std::conj(det2) ; 
			} else if ( pb2corr == 3 ) {
					_iata  =  _data / g1ll / std::conj(g2ll) / det1 / std::conj(det2) ; 
			}

			// write to file
			of << ichan << " " << ant1[0] << ant1[1] << ant1[2] << " " << ant2[0] << ant2[1] << ant2[2] << " " << bandmap[band1] <<  bandmap[band2] << " " << _data << " " << _iata << std::endl;

			// if self term
			if ( ant1 == ant2 ) {
				if ( pb2corr == 0 ) {
					rr_unapplied += _data.real();

					rr_applied   += _iata.real();
				} // rr
				else if ( pb2corr == 3 ) {
					ll_unapplied += _data.real();

					ll_applied   += _iata.real();
				} // ll
			} // self 
			else {
				if ( pb2corr == 0 ) {
					rr_unapplied += 2.0 * _data.real();

					rr_applied   += 2.0 * _iata.real();
				} // rr
				else if ( pb2corr == 3 ) {
					ll_unapplied += 2.0 * _data.real();

					ll_applied   += 2.0 * _iata.real();
				} // ll
			} // nonself baseline

		} /* polarbaselines */

		frr_un[ichan] = rr_unapplied;
		fll_un[ichan] = ll_unapplied;

		frr_ap[ichan] = rr_applied;
		fll_ap[ichan] = ll_applied;

		if ( std::isnan ( rr_applied ) ) rr_applied = rr_unapplied;
		if ( std::isnan ( ll_applied ) ) ll_applied = ll_unapplied;

		ratio_rr[ichan] = rr_unapplied/rr_applied;
		ratio_ll[ichan] = ll_unapplied/ll_applied;

	} // channel
	
	// find median of ratio
	const LBFGS::real_type scale_rr ( get_median ( ratio_rr ) );
	const LBFGS::real_type scale_ll ( get_median ( ratio_ll ) );

	std::cout << " scaling rr=" << scale_rr << " ll=" << scale_ll << std::endl;

	return 0;
}


