/*
 * Average scan and write selfcorr 
 *
 */

#include<iostream>

#include <unistd.h>

#include <fstream>
#include <complex>

#include "lta_file.hpp"
#include "gaintable.hpp"

#ifdef CUSTOM
#include "custom_file.hpp"
#endif

/* overload operator<< for my formats */
std::ostream& operator<< (std::ostream& os, const gaintable::antname_t& a ) {
	//os << a[0] << a[1] << a[2] << a[3];
	// ant names are only three characters anyway
	os << a[0] << a[1] << a[2];
	return os;
}
std::ostream& operator<< (std::ostream& os, const gaintable::gain_type& g ) {
	/* these flags are set before priting */
	//os << std::fixed << std::setprecision(3) << std::showpos;
	os << g.real() << g.imag() << 'j';
	return os;
}

void print_help () {
	std::cout << "write_self" << std::endl;
	std::cout << "  Write baselines in ant, chan, rr, ll, rl " << std::endl;
	std::cout << std::endl;
	std::cout << " write_baselines [ARGUMENTS] LTA_FILE" << std::endl;
	std::cout << "    -h Print help" << std::endl;
	std::cout << "    -s <scan> scan number of the LTA file" << std::endl;
	std::cout << "    -o <ofile> output file to save baseline data" << std::endl;
	std::cout << std::endl;
}


int main(int argc, char *argv[]) {

	/* hello getopt, my old friend */
	int opt;
	int cal_scan_number;
	std::string ofile ("baselines.out");
	std::string lta_path;

	if ( argc < 2 ) {
		print_help ();
		exit (EXIT_SUCCESS);
	}

	while ( (opt = getopt ( argc, argv, "hs:o:" )) != -1 ) {
		switch (opt) {
			case 'h':
				print_help ();
				exit (EXIT_SUCCESS);
				break;
			case 's':
				cal_scan_number = atoi ( optarg );
				break;
			case 'o':
				ofile  = optarg;
				break;
		} // switch
	} // getopt
	if ( optind >= argc ) {
		print_help ();
		exit (EXIT_SUCCESS);
	}
	/* lta file*/
	lta_path  = argv[optind];
	optind++;

	std::cout << "[inputs] lta=" << lta_path << std::endl;
	std::cout << "[inputs] ofile=" << ofile << std::endl;
	std::cout << "[inputs] scan=" << cal_scan_number << std::endl;

	/***************************************/
	/*      READ LTA FILE                  */
	/***************************************/
#ifdef CUSTOM
	custom_file lta_file;
#else
	LTA lta_file ( lta_path );
#endif

	int nbaselines   = lta_file.nbaselines;
	int nchannels    = lta_file.nchannels;
	float fbw        = lta_file.fbw;
	float fedge      = lta_file.fedge;

	/***************************************/
	/*      AVERAGE LTA SCAN               */
	/***************************************/
	/* maybe directly read complex<float>  */
  /* NOTICE:avgbldata contains self-terms*/
	LTA::vf_type avgbldata ( 2 * nbaselines * nchannels, 0.0f );
	lta_file.time_average  ( cal_scan_number, avgbldata );

	/***************************************/
	/*      WRITE DIRECTLY TO FILE         */
	/***************************************/
	// to map band{0,1} -> band{r,l}
	constexpr std::array<char,2> bandmap{'r','l'};
	std::ofstream of ( ofile );

	of << "ant chan corr real" << std::endl;
	of << std::fixed << std::setprecision(3) << std::showpos;

	for (int ichan = 0; ichan < nchannels; ichan++) {
		for (int ib = 0; ib < nbaselines; ib++) {

			/* get baseline object */
			const baseline_t&  _bl = lta_file.baselines [ ib ];

			const auto& ant1  = _bl.ant1;
			const auto& ant2  = _bl.ant2;

			if ( ant1 != ant2 ) continue;

			const auto& band1 = _bl.band1;
			const auto& band2 = _bl.band2;

			/* get data */
			/* index in (baseline, channel) complex<float> */
			const int _i      = ichan + nchannels * ib;
			const float real ( avgbldata[2*_i] );
			const float imag ( avgbldata[2*_i+1] );

			/*
			 * selfcorr are special, rr,ll are purely real
			 * rl,lr are complex conjugate pair,
			 * so we apply a trick
			 *
			 * rl is real(rl) and lr is imag(rl)
			*/

			// write to file
			if ( band1 == 0 && band2 == 0 ) {
				of << ant1 << " " << ichan << " " << bandmap[band1] <<  bandmap[band2] << " " << real << std::endl;
			} 
			else if ( band1 == 0 && band2 == 1 ) {
				of << ant1 << " " << ichan << " " << bandmap[band1] <<  bandmap[band2] << " " << real << std::endl;
			}
			else if ( band1 == 1 && band2 == 0 ) {
				of << ant1 << " " << ichan << " " << bandmap[band1] <<  bandmap[band2] << " " << imag << std::endl;
			}
			else if ( band1 == 1 && band2 == 1 ) {
				of << ant1 << " " << ichan << " " << bandmap[band1] <<  bandmap[band2] << " " << real << std::endl;
			}

		} // baseline
	} // channel


	return 0;
}

