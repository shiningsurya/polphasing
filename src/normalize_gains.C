/*
 * read gaintables
 * normalize
 * and write gaintables
 */

#include<iostream>

#include <unistd.h>

#include <fstream>
#include <complex>

#include "gaintable.hpp"


void print_help () {
	std::cout << "normalize_gains" << std::endl;
	std::cout << "  Normalizes gains " << std::endl;
	std::cout << std::endl;
	std::cout << " normalize_gains [ARGUMENTS] " << std::endl;
	std::cout << "    -h Print help" << std::endl;
	std::cout << "    -t <tag> tag of the gain solutions, will read {rr,rl,lr,ll}" << std::endl;
	std::cout << " For now, only uses C00 as reference antenna " << std::endl;
	std::cout << " multiplies from right " << std::endl;
	std::cout << std::endl;
}

int main(int argc, char *argv[]) {

	/* hello getopt, my old friend */
	int opt;

	std::string tag;
	std::string oag;

	if ( argc < 2 ) {
		print_help ();
		exit (EXIT_SUCCESS);
	}

	while ( (opt = getopt ( argc, argv, "ht:" )) != -1 ) {
		switch (opt) {
			case 'h':
				print_help ();
				exit (EXIT_SUCCESS);
				break;
			case 't':
				tag  = optarg;
				break;
		} // switch
	} // getopt
	
	oag  = tag + std::string("_normalized");

	const auto& antlist ( gaintable::sol_ant_order );
	const auto refant ( antlist[0] );
		
	std::cout << "[inputs] tag=" << tag << std::endl;

	const std::string read_file_rr  = tag + std::string("_rr.gains");
	const std::string read_file_rl  = tag + std::string("_rl.gains");
	const std::string read_file_lr  = tag + std::string("_lr.gains");
	const std::string read_file_ll  = tag + std::string("_ll.gains");

	const std::string write_file_rr  = oag + std::string("_rr.gains");
	const std::string write_file_rl  = oag + std::string("_rl.gains");
	const std::string write_file_lr  = oag + std::string("_lr.gains");
	const std::string write_file_ll  = oag + std::string("_ll.gains");

	/***************************************/
	/*      READ GAIN FILES                */
	/***************************************/

	gaintable::gaintable_t grr  = gaintable::read_complex_solutions ( read_file_rr );
	gaintable::gaintable_t grl  = gaintable::read_complex_solutions ( read_file_rl );
	gaintable::gaintable_t glr  = gaintable::read_complex_solutions ( read_file_lr );
	gaintable::gaintable_t gll  = gaintable::read_complex_solutions ( read_file_ll );

	const int nchannels ( grr[refant].size() );

	gaintable::gaintable_t    hrr = gaintable::prepare_gaintables ( nchannels );
	gaintable::gaintable_t    hrl = gaintable::prepare_gaintables ( nchannels );
	gaintable::gaintable_t    hlr = gaintable::prepare_gaintables ( nchannels );
	gaintable::gaintable_t    hll = gaintable::prepare_gaintables ( nchannels );

	/***************************************/
	/*      NORMALIZE GAINS                */
	/***************************************/

	for (int ichan = 0; ichan < nchannels; ichan++) {
	//for (int ichan = 500; ichan < 501; ichan++) {

		/* make normalizer object */
		const gaintable::normalizer_t  normer ( 
				grr[refant][ichan],
				grl[refant][ichan],
				glr[refant][ichan],
				gll[refant][ichan]
		);

		/* apply on every antenna including reference antenna */
		for ( const auto& iant : antlist ) {

			const auto& ret = normer (  
				grr[iant][ichan],
				grl[iant][ichan],
				glr[iant][ichan],
				gll[iant][ichan]
			);

			hrr[iant][ichan]  = ret[0];
			hrl[iant][ichan]  = ret[1];
			hlr[iant][ichan]  = ret[2];
			hll[iant][ichan]  = ret[3];

		} // antenna

	} // nchan

	/***************************************/
	/*        WRITE GAINTABLES             */
	/***************************************/
	gaintable::write_complex_solutions ( hrr, write_file_rr );
	gaintable::write_complex_solutions ( hrl, write_file_rl );
	gaintable::write_complex_solutions ( hlr, write_file_lr );
	gaintable::write_complex_solutions ( hll, write_file_ll );

	return 0;
}
