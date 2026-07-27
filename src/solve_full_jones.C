#include <iostream>
#include <cmath>

#include <unistd.h>

#include <string>

#include "lta_file.hpp"
#include "gaintable.hpp"
#include "models.hpp"
#include "logsolve.hpp"
#include "fullgdsolver.hpp"

#ifdef TIMING
#include <chrono>
#endif

//#define CHANDEBUG


#ifdef TIMING
auto start  = std::chrono::high_resolution_clock::now();
auto end    = std::chrono::high_resolution_clock::now();
#endif 

void print_help () {
	std::cout << "solve_full_jones" << std::endl;
	std::cout << "  Solving full Jones matrix for each antenna using IQU models of two sources using multithreading " << std::endl;
	std::cout << std::endl;
	std::cout << " fullpolphase [ARGUMENTS] LTA_FILE" << std::endl;
	std::cout << "    -h Print help" << std::endl;
	std::cout << "    -t <tag> tag/stem with which to save log and complex gains" << std::endl;
	std::cout << "    -i <scan> scan number of the first source in the LTA_FILE" << std::endl;
	std::cout << "    -j <scan> scan number of the second source in the LTA_FILE" << std::endl;
	std::cout << "    -m <model> path to model file of the first source" << std::endl;
	std::cout << "    -n <model> path to model file of the second source" << std::endl;
	std::cout << "Notes:" << std::endl;
	std::cout << "   1. One of the source is expected to be unpolarized. " << std::endl;
	std::cout << "      Not strictly required but is a real life constraint. " << std::endl;
	std::cout << std::endl;
}

int main(int argc, char *argv[]) {

	/* hello getopt, my old friend */
	int         opt;
	int         first_scan_number, second_scan_number;
	std::string first_model_path, second_model_path;
	std::string tag;
	std::string lta_path;

	if ( argc < 2 ) {
		print_help ();
		exit (EXIT_SUCCESS);
	}

	while ( (opt = getopt ( argc, argv, "ht:i:j:m:n:" )) != -1 ) {
		switch (opt) {
			case 'h':
				print_help ();
				exit (EXIT_SUCCESS);
				break;
			case 'i':
				first_scan_number  = atoi ( optarg );
				break;
			case 'j':
				second_scan_number = atoi ( optarg );
				break;
			case 't':
				tag  = optarg;
				break;
			case 'm':
				first_model_path  = optarg;
				break;
			case 'n':
				second_model_path = optarg;
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

	/* other files */
	const std::string save_file_rr  = tag + std::string("_rr.gains");
	const std::string save_file_rl  = tag + std::string("_rl.gains");
	const std::string save_file_lr  = tag + std::string("_lr.gains");
	const std::string save_file_ll  = tag + std::string("_ll.gains");
	const std::string log_file      = tag + std::string(".log");

	std::cout << "[inputs] lta=" << lta_path << " tag=" << tag << std::endl;
	std::cout << "[inputs] first model scan=" << first_scan_number << " model=" << first_model_path << std::endl;
	std::cout << "[inputs] second model scan=" << second_scan_number << " model=" << second_model_path << std::endl;

	/***************************************/
	/*      READ LTA FILE                  */
	/***************************************/
	LTA lta_file ( lta_path );

	int nbaselines   = lta_file.nbaselines;
	int nchannels    = lta_file.nchannels;
	float fbw        = lta_file.fbw;
	float fedge      = lta_file.fedge;

	/***************************************/
	/*      READ MODEL FILE                */
	/***************************************/
	models::model_data_t first_model  = models::read_model_file ( first_model_path );
	models::model_data_t second_model = models::read_model_file ( second_model_path );

	/***************************************/
	/*      AVERAGE LTA SCAN               */
	/***************************************/
	/* maybe directly read complex<float>  */
  /* NOTICE:avgbldata contains self-terms*/
	LTA::vf_type first_avgbldata  ( 2 * nbaselines * nchannels, 0. );
	LTA::vf_type second_avgbldata ( 2 * nbaselines * nchannels, 0. );

	lta_file.time_average  ( first_scan_number, first_avgbldata );
	lta_file.time_average  ( second_scan_number, second_avgbldata );

	/***************************************/
	/*        SOLVER PREPARE               */
	/***************************************/

	/* Following is for mapping */
	/* (1) Mapping antenna to index  -- ant2idx */

	/**/
	int rant    = 0;
	std::map<LTA::antname_t,int>  ant2idx;
	/* contains noself_baseline index */
	std::vector<int>     noself_baselines;

	/*
	 * This loop only computes ant2idx
	 * and collects noself baselines 
	 */
	for (int ib = 0; ib < nbaselines; ib++) {

		/* get baseline object */
		const baseline_t&  _bl = lta_file.baselines [ ib ];

		/* building ant2idx */
		const auto& ant1  = _bl.ant1;
		const auto& band1 = _bl.band1;
		if ( ant2idx.find ( ant1 ) == ant2idx.end() ) {
			ant2idx.emplace ( ant1, rant );
			rant++;
		}
		const auto& ant2  = _bl.ant2;
		const auto& band2 = _bl.band2;
		if ( ant2idx.find ( ant2 ) == ant2idx.end() ) {
			ant2idx.emplace ( ant2, rant );
			rant++;
		}

		/* we want noself baselines */
		if ( ant1 != ant2 ) noself_baselines.push_back ( ib );

	} /* polar baselines */

	/***************************************/
	/*        SOLVER RUN                   */
	/***************************************/
	/*
	 * Need to ignore self-terms
	 *
	 * given nantennas, 
	 * nbaselines will be 0.5*nantennas*(nantennas+1)*4
	 *
	 * noself_nbaselines will be 0.5*nantennas*(nantennas-1)*4
	 */
	int nantennas    = ant2idx.size(); 
	int n_noself_baselines  = noself_baselines.size();
	if ( n_noself_baselines != (0.5 * nantennas * ( nantennas - 1 ) * 4) ) {
		throw std::runtime_error ("baseline count mismatch");
	}

	/* logging */
	logging::log_t            logger ( nchannels );

	/* gain tables */
	gaintable::gaintable_t    solved_gains_rr = gaintable::prepare_gaintables ( nchannels );
	gaintable::gaintable_t    solved_gains_rl = gaintable::prepare_gaintables ( nchannels );
	gaintable::gaintable_t    solved_gains_lr = gaintable::prepare_gaintables ( nchannels );
	gaintable::gaintable_t    solved_gains_ll = gaintable::prepare_gaintables ( nchannels );

	/* main loop */
	std::cout << " Starting main solving loop" << std::endl;

	auto total_start = std::chrono::high_resolution_clock::now();
	
#ifdef CHANDEBUG
	/* channel 458 has high loss.*/
	for (int ichan = 458; ichan < 459; ichan++) {
#else
	#pragma omp parallel for num_threads(4) 
	for (int ichan = 0; ichan < nchannels; ichan++) {
#endif
		/* when parallelizing inside loop */
		/* this will be doing a lot of mallocs */

		/* read stokes IQU for ichan */
		const models::real_type    first_stokes_i ( first_model.stokes_i[ichan] );
		const models::real_type    first_stokes_q ( first_model.stokes_q[ichan] );
		const models::real_type    first_stokes_u ( first_model.stokes_u[ichan] );

		const models::real_type    second_stokes_i ( second_model.stokes_i[ichan] );
		const models::real_type    second_stokes_q ( second_model.stokes_q[ichan] );
		const models::real_type    second_stokes_u ( second_model.stokes_u[ichan] );

		/* populate rr, rl, lr, ll */
		const models::complex_type first_model_rr ( first_stokes_i, 0.0f );
		const models::complex_type first_model_rl ( first_stokes_q, first_stokes_u );
		const models::complex_type first_model_lr ( first_stokes_q,-first_stokes_u );
		const models::complex_type first_model_ll ( first_stokes_i, 0.0f );

		const models::complex_type second_model_rr ( second_stokes_i, 0.0f );
		const models::complex_type second_model_rl ( second_stokes_q, second_stokes_u );
		const models::complex_type second_model_lr ( second_stokes_q,-second_stokes_u );
		const models::complex_type second_model_ll ( second_stokes_i, 0.0f );

		/* data package */
		/*
		 * There is some amount of data duplication between first_pkg and second_pkg. 
		 * The index_b1, index_b2 are the same. 
		 * But it is so coupled with everything that i do not want to decouple now.
		*/
		FullGDSolver::solve_data_t  first_pkg ( 
				n_noself_baselines, 
				first_model_rr, first_model_rl, first_model_lr, first_model_ll
		);
		FullGDSolver::solve_data_t  second_pkg ( 
				n_noself_baselines, 
				second_model_rr, second_model_rl, second_model_lr, second_model_ll
		);

		/* initialize data */
		for ( int ii = 0; ii < n_noself_baselines; ii++ ) {

			/* ensures this is index of noself baseline */
			/* which is consistent with lta_file        */
			/* ib only for reading */
			const int ib     = noself_baselines[ii];

			/* index in (baseline, channel) complex<float> */
			const int _i     = ichan + nchannels * ib;

			/* get baseline object */
			const baseline_t&  _bl = lta_file.baselines [ ib ];

			/* every baseline is polar baseline */
			/* although redundant, we write a model column */
			/* because it makes future computations straightforward */
			/*
			 * (0,0) = rr -> 0 
			 * (0,1) = rl -> 1 
			 * (1,0) = lr -> 2 
			 * (1,1) = ll -> 3
			 *
			 * 20260721: 
			 * Actually, see comment before.
			 * This order is really important. 
			 */

			/* get antennas */
			const auto& ant1  = _bl.ant1;
			const auto& ant2  = _bl.ant2;

			/* get bands */
			const auto& band1 = _bl.band1;
			const auto& band2 = _bl.band2;

			/* ID correlation */
			const int pb2corr = _bl.band1*2 + _bl.band2;

			/* index_b1 b2 */
			const auto& iant1 = ant2idx.at(ant1);
			const auto& iant2 = ant2idx.at(ant2);

			/* copy data */
			const float _first_real ( first_avgbldata[2*_i] );
			const float _first_imag ( first_avgbldata[2*_i + 1] );

			const float _second_real ( second_avgbldata[2*_i] );
			const float _second_imag ( second_avgbldata[2*_i + 1] );

			// data
			first_pkg.data  [ ii ]    = FullGDSolver::complex_type (  _first_real,  _first_imag );
			second_pkg.data [ ii ]    = FullGDSolver::complex_type (  _second_real,  _second_imag );

			// antenna indices
			first_pkg.iant1 [ ii ]    = second_pkg.iant1[ ii ] = iant1;
			first_pkg.iant2 [ ii ]    = second_pkg.iant2[ ii ] = iant2;

			// correlation index
			first_pkg.pb2corr  [ ii ] = second_pkg.pb2corr [ ii ] = pb2corr;
		} /* baselines */

		/* perform solving */
#ifdef TIMING
		start  = std::chrono::high_resolution_clock::now();
#endif
		FullGDSolver               solver (n_noself_baselines, nantennas);
		FullGDSolver::vc_type      isol ( solver.ngains, FullGDSolver::complex_type (1.0f, 1.0f) );
		auto cost = solver.solve ( first_pkg, second_pkg, isol );

#ifdef TIMING
		end   = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double,std::milli> duration = end - start;
		logger.time_chan [ ichan ] = duration.count();
#endif

		logger.sse_chan [ ichan ]  = cost;
		logger.nfev [ ichan ]      = solver.niter;
		logger.info [ ichan ]      = solver.rcode;

		/* save into gain table */
		for (auto _i = ant2idx.begin(); _i != ant2idx.end(); ++_i) {

			const auto& iant = _i->first;
			const auto& idx  = _i->second;

			const FullGDSolver::complex_type rr ( isol[4*idx + 0] );
			const FullGDSolver::complex_type rl ( isol[4*idx + 1] );
			const FullGDSolver::complex_type lr ( isol[4*idx + 2] );
			const FullGDSolver::complex_type ll ( isol[4*idx + 3] );

			solved_gains_rr[iant][ichan]  = rr;
			solved_gains_rl[iant][ichan]  = rl;
			solved_gains_lr[iant][ichan]  = lr;
			solved_gains_ll[iant][ichan]  = ll;

		} /* ant2idx */

	} // channel

	/***************************************/
	/*        WRITE GAINTABLES             */
	/***************************************/
	gaintable::write_complex_solutions ( solved_gains_rr, save_file_rr );
	gaintable::write_complex_solutions ( solved_gains_rl, save_file_rl );
	gaintable::write_complex_solutions ( solved_gains_lr, save_file_lr );
	gaintable::write_complex_solutions ( solved_gains_ll, save_file_ll );

	/***************************************/
	/*        WRITE LOG                    */
	/***************************************/
	logging::write_log ( logger, log_file );

	auto total_end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> total_duration = total_end - total_start;

	std::cout << std::setprecision(3) << std::endl << " Total solving took " << total_duration.count() << " seconds ..." << std::endl;

	return 0;
}

