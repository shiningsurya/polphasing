#include <iostream>
#include <cmath>

#include <unistd.h>

#include <string>

#include "lta_file.hpp"
#include "gaintable.hpp"
#include "models.hpp"
#include "logsolve.hpp"
#include "LBFGS.hpp"
#include "ants.hpp"

#ifdef CUSTOM
#include "custom_file.hpp"
#endif

#ifdef TIMING
#include <chrono>
#endif

//#define CHANDEBUG


#ifdef TIMING
auto start  = std::chrono::high_resolution_clock::now();
auto end    = std::chrono::high_resolution_clock::now();
#endif 

void print_help () {
	std::cout << "fullpolphase" << std::endl;
	std::cout << "  Solving full Jones matrix for each antenna using IQU model using multithreading " << std::endl;
	std::cout << std::endl;
	std::cout << " fullpolphase [ARGUMENTS] LTA_FILE" << std::endl;
	std::cout << "    -h Print help" << std::endl;
	std::cout << "    -s <scan> scan number of the LTA file" << std::endl;
	std::cout << "    -t <tag> tag/stem with which to save log and complex gains" << std::endl;
	std::cout << "    -m <model> path to model file" << std::endl;
	std::cout << std::endl;
}

int main(int argc, char *argv[]) {

	/* hello getopt, my old friend */
	int opt;
	int cal_scan_number;
	std::string tag;
	std::string model_path;
	std::string lta_path;

	if ( argc < 2 ) {
		print_help ();
		exit (EXIT_SUCCESS);
	}

	while ( (opt = getopt ( argc, argv, "hs:t:m:" )) != -1 ) {
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
			case 'm':
				model_path = optarg;
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

	std::cout << "[inputs] lta="   << lta_path << std::endl;
	std::cout << "[inputs] model=" << model_path << std::endl;
	std::cout << "[inputs] tag="   << tag  << std::endl;
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

	int nbaselines   = lta_file.nbaselines;
	int nchannels    = lta_file.nchannels;
	float fbw        = lta_file.fbw;
	float fedge      = lta_file.fedge;

	/***************************************/
	/*      READ MODEL FILE                */
	/***************************************/
	models::model_data_t calmodel = models::read_model_file ( model_path );

	/***************************************/
	/*      AVERAGE LTA SCAN               */
	/***************************************/
	/* maybe directly read complex<float>  */
  /* NOTICE:avgbldata contains self-terms*/
	LTA::vf_type avgbldata ( 2 * nbaselines * nchannels, 0. );

	const LTA::scan_t cal_scan = lta_file.get_scan ( cal_scan_number );

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
	/*        MEASURE PAR ANGLES           */ 
	/***************************************/

	/*
	 * We need coordinates of each antenna. 
	 * > we can get them from casa and keep them as static.
	 *
	 * We need source coordinates, LST of observation
	 * > we have to get them from lta_file
	 *
	 * --------
	 *  We put the par angles in map ant2par. 
	 *  <antname_t,real_type>
	*/

	const ants::ant2par_t antpar = ants::calculate_parallactic_angle ( cal_scan.mjd, cal_scan.ra, cal_scan.dec );

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
	int nantennas    = ant2idx.size(); 
	int n_noself_baselines  = noself_baselines.size();
	if ( n_noself_baselines != (0.5 * nantennas * ( nantennas - 1 ) * 4) ) {
		throw std::runtime_error ("baseline count mismatch");
	}

	/* logging */
	logging::log_t   logger ( nchannels );

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
	for (int ichan = 501; ichan < 502; ichan++) {
#else
	#pragma omp parallel for num_threads(4) 
	for (int ichan = 0; ichan < nchannels; ichan++) {
#endif
		/* when parallelizing inside loop */
		/* this will be doing a lot of mallocs */

		/* read stokes IQU for ichan */
		const models::real_type    stokes_i ( calmodel.stokes_i[ichan] );
		const models::real_type    stokes_q ( calmodel.stokes_q[ichan] );
		const models::real_type    stokes_u ( calmodel.stokes_u[ichan] );
		const models::real_type    stokes_l ( std::sqrt ( stokes_q*stokes_q + stokes_u*stokes_u ) );

		/* populate rr, rl, lr, ll */
		const models::complex_type model_rr ( stokes_i, 0.0f );
		const models::complex_type model_rl ( stokes_q, stokes_u );
		const models::complex_type model_lr ( stokes_q,-stokes_u );
		const models::complex_type model_ll ( stokes_i, 0.0f );

		/* data package */
		LBFGS::data_t  pkg ( 
				n_noself_baselines, 
				nantennas,
				model_rr, model_rl, model_lr, model_ll
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

			/*
			 * (0,0) = rr -> 0 
			 * (0,1) = rl -> 1 
			 * (1,0) = lr -> 2 
			 * (1,1) = ll -> 3
			 */

			/* get antennas */
			const auto& ant1  = _bl.ant1;
			const auto& ant2  = _bl.ant2;

			/* get bands */
			const auto& band1 = _bl.band1;
			const auto& band2 = _bl.band2;


			/* ID correlation */
			const int pb2corr = band1*2 + band2;

			/* index_b1 b2 */
			const auto& iant1 = ant2idx.at(ant1);
			const auto& iant2 = ant2idx.at(ant2);

			/* copy data */
			const float _real ( avgbldata[2*_i] );
			const float _imag ( avgbldata[2*_i + 1] );

			// data
			pkg.data [ ii ]     = LBFGS::complex_type (  _real,  _imag );

			// parallactic angle correct model
			// find parallactic angle
			const float _par1 ( antpar.at(ant1) );
			const float _par2 ( antpar.at(ant2) );

			// perform correction
			const auto& _par_model = models::parallactic_correction ( _par1, _par2, model_rr, model_rl, model_lr, model_ll );
			//const auto& _par_model = models::parallactic_correction ( -_par1, -_par2, model_rr, model_rl, model_lr, model_ll );

			pkg.par_model_rr [ ii ] = _par_model[0];
			pkg.par_model_rl [ ii ] = _par_model[1];
			pkg.par_model_lr [ ii ] = _par_model[2];
			pkg.par_model_ll [ ii ] = _par_model[3];

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

		/* full jones has four complex gains per antenna */
		/* LBFGS separates real and imaginary so double it*/
		const int npar        ( nantennas * 4 * 2 );
		LBFGS::Solver                 solver (npar);
		auto cost = solver.solve_full_jones ( pkg );

#ifdef CHANDEBUG 
		{
			std::ofstream of("pkgchandebug.bin", std::ios::binary);
of.write (reinterpret_cast<const char*>(&pkg.npolarbaselines), sizeof(int));
of.write (reinterpret_cast<const char*>(&pkg.nantennas), sizeof(int));
of.write (reinterpret_cast<const char*>(&pkg.mrr), sizeof(LBFGS::complex_type));
of.write (reinterpret_cast<const char*>(&pkg.mrl), sizeof(LBFGS::complex_type));
of.write (reinterpret_cast<const char*>(&pkg.mlr), sizeof(LBFGS::complex_type));
of.write (reinterpret_cast<const char*>(&pkg.mll), sizeof(LBFGS::complex_type));
			/* data */
of.write (reinterpret_cast<const char*>(pkg.data.data()), pkg.data.size()*sizeof(LBFGS::complex_type));
			/* models */
of.write (reinterpret_cast<const char*>(pkg.par_model_rr.data()), pkg.par_model_rr.size()*sizeof(LBFGS::complex_type));
of.write (reinterpret_cast<const char*>(pkg.par_model_rl.data()), pkg.par_model_rl.size()*sizeof(LBFGS::complex_type));
of.write (reinterpret_cast<const char*>(pkg.par_model_lr.data()), pkg.par_model_lr.size()*sizeof(LBFGS::complex_type));
of.write (reinterpret_cast<const char*>(pkg.par_model_ll.data()), pkg.par_model_ll.size()*sizeof(LBFGS::complex_type));
of.write (reinterpret_cast<const char*>(pkg.iant1.data()), pkg.iant1.size()*sizeof(int));
of.write (reinterpret_cast<const char*>(pkg.iant2.data()), pkg.iant2.size()*sizeof(int));
of.write (reinterpret_cast<const char*>(pkg.pb2corr.data()), pkg.pb2corr.size()*sizeof(int));
of.write (reinterpret_cast<const char*>(solver.xpar), npar*sizeof(float));
		}
#endif

#ifdef TIMING
		end   = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double,std::milli> duration = end - start;
		logger.time   [ ichan ]  = duration.count();
#endif

		logger.sse    [ ichan ]  = cost;
		logger.nfev   [ ichan ]  = solver.niter;
		logger.info   [ ichan ]  = solver.rcode;
		logger.gnorm  [ ichan ]  = solver.gnorm;

		/* save into gain table */
		for (auto _i = ant2idx.begin(); _i != ant2idx.end(); ++_i) {

			const auto& iant = _i->first;
			const auto& idx  = _i->second;

			const LBFGS::complex_type  rr ( solver.xpar[8*idx + 0], solver.xpar[8*idx + 1] );
			const LBFGS::complex_type  rl ( solver.xpar[8*idx + 2], solver.xpar[8*idx + 3] );
			const LBFGS::complex_type  lr ( solver.xpar[8*idx + 4], solver.xpar[8*idx + 5] );
			const LBFGS::complex_type  ll ( solver.xpar[8*idx + 6], solver.xpar[8*idx + 7] );

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

