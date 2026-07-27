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
	std::cout << "solve_model" << std::endl;
	std::cout << "  Solving IQU model using solved Jones matrices of each antenna and the data used to generate those solutions." << std::endl;
	std::cout << std::endl;
	std::cout << " solve_model [ARGUMENTS] LTA_FILE" << std::endl;
	std::cout << "    -h Print help" << std::endl;
	std::cout << "    -s <scan> scan number of the LTA file" << std::endl;
	std::cout << "    -t <tag> tag/stem with which solutions are saved. " << std::endl;
	std::cout << std::endl;
}

int main(int argc, char *argv[]) {

	/* hello getopt, my old friend */
	int opt;
	int cal_scan_number;
	std::string tag;
	std::string lta_path;

	if ( argc < 2 ) {
		print_help ();
		exit (EXIT_SUCCESS);
	}

	while ( (opt = getopt ( argc, argv, "hs:t:" )) != -1 ) {
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
	const std::string save_model  = tag + std::string("_solved_model.txt");
	const std::string log_file    = tag + std::string("_solved_model.log");

	std::cout << "[inputs] lta=" << lta_path << std::endl;
	std::cout << "[inputs] tag=" << tag << " scan=" << cal_scan_number << std::endl;

	/***************************************/
	/*      READ SOLUTIONS                 */
	/***************************************/
	const std::string file_solved_rr  = tag + std::string("_rr.gains");
	const std::string file_solved_rl  = tag + std::string("_rl.gains");
	const std::string file_solved_lr  = tag + std::string("_lr.gains");
	const std::string file_solved_ll  = tag + std::string("_ll.gains");

	const gaintable::gaintable_t rrgains = gaintable::read_complex_solutions ( file_solved_rr );
	const gaintable::gaintable_t rlgains = gaintable::read_complex_solutions ( file_solved_rl );
	const gaintable::gaintable_t lrgains = gaintable::read_complex_solutions ( file_solved_lr );
	const gaintable::gaintable_t llgains = gaintable::read_complex_solutions ( file_solved_ll );
	/***************************************/
	/*      READ LTA FILE                  */
	/***************************************/
	LTA lta_file ( lta_path );

	int nbaselines   = lta_file.nbaselines;
	int nchannels    = lta_file.nchannels;
	float fbw        = lta_file.fbw;
	float fedge      = lta_file.fedge;

	/***************************************/
	/*      AVERAGE LTA SCAN               */
	/***************************************/
	/* maybe directly read complex<float>  */
  /* NOTICE:avgbldata contains self-terms*/
	LTA::vf_type avgbldata ( 2 * nbaselines * nchannels, 0. );
	lta_file.time_average  ( cal_scan_number, avgbldata );

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
	 *
	 * ndata = 2 * noself_nbaselines
	 * npar  = 4 * nantennas
	 *
	 */
	const int nantennas     = ant2idx.size(); 
	const int                 ngains ( nantennas * 4 );
	int n_noself_baselines  = noself_baselines.size();
	if ( n_noself_baselines != (0.5 * nantennas * ( nantennas - 1 ) * 4) ) {
		throw std::runtime_error ("baseline count mismatch");
	}

	/* logging */
	logging::log_t   logger ( nchannels );

	/* solved model */
	models::model_vis_t solved_model ( nchannels );

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

		/* data package */
		FullGDSolver::solve_model_t  pkg ( 
				n_noself_baselines, 
				ngains
		);

		/* read gains into gains vector */
		for ( auto _i = ant2idx.begin(); _i != ant2idx.end(); ++_i ) {

			const auto& ant  = _i->first;
			const auto& idx  = _i->second;

			pkg.gains[4*idx + 0] = rrgains.at(ant)[ichan];
			pkg.gains[4*idx + 1] = rlgains.at(ant)[ichan];
			pkg.gains[4*idx + 2] = lrgains.at(ant)[ichan];
			pkg.gains[4*idx + 3] = llgains.at(ant)[ichan];
		}

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
			const float _real ( avgbldata[2*_i] );
			const float _imag ( avgbldata[2*_i + 1] );

			// data
			pkg.data [ ii ]     = FullGDSolver::complex_type (  _real,  _imag );

			// antenna indices
			pkg.iant1 [ ii ]    = iant1;
			pkg.iant2 [ ii ]    = iant2;

			// correlation index
			pkg.pb2corr  [ ii ] = pb2corr;
		} /* baselines */

		/* perform solving */
#ifdef TIMING
		start  = std::chrono::high_resolution_clock::now();
#endif
		FullGDSolver               solver (n_noself_baselines, nantennas);
		FullGDSolver::complex_type mrr, mrl, mlr, mll;
		auto cost = solver.solve ( pkg, mrr, mrl, mlr, mll );

#ifdef TIMING
		end   = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double,std::milli> duration = end - start;
		logger.time_chan [ ichan ] = duration.count();
#endif

		logger.sse_chan [ ichan ]  = cost;
		logger.nfev [ ichan ]      = solver.niter;
		logger.info [ ichan ]      = solver.rcode;

		/* save into model table */
		solved_model.rr[ichan]     = mrr;
		solved_model.rl[ichan]     = mrl;
		solved_model.lr[ichan]     = mlr;
		solved_model.ll[ichan]     = mll;

	} // channel

	/***************************************/
	/*        WRITE SOLVED MODEL           */
	/***************************************/
	models::write_model_file ( solved_model, save_model );

	/***************************************/
	/*        WRITE LOG                    */
	/***************************************/
	logging::write_log ( logger, log_file );

	auto total_end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> total_duration = total_end - total_start;

	std::cout << std::setprecision(3) << std::endl << " Total solving took " << total_duration.count() << " seconds ..." << std::endl;

	return 0;
}
