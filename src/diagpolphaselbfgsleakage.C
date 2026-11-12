#include <iostream>
#include <cmath>

#include <unistd.h>

#include <string>

#include "lta_file.hpp"
#include "gaintable.hpp"
#include "selftable.hpp"
#include "models.hpp"
#include "logsolve.hpp"
#include "ants.hpp"

#ifdef CUSTOM
#include "custom_file.hpp"
#endif

#ifdef TIMING
#include <chrono>
#endif

#include "LBFGSpolarizediterative.hpp"

//#define CHANDEBUG


#ifdef TIMING
auto start  = std::chrono::high_resolution_clock::now();
auto end    = std::chrono::high_resolution_clock::now();
#endif 

void print_help () {
	std::cout << "diagpolphaseleakage" << std::endl;
	std::cout << "  Solving diagonal Jones matrix for each antenna using IQU model using multithreading " << std::endl;
	std::cout << std::endl;
	std::cout << " diagpolphaseleakage [ARGUMENTS] LTA_FILE" << std::endl;
	std::cout << "    -h Print help" << std::endl;
	std::cout << "    -s <scan> scan number of the LTA file" << std::endl;
	std::cout << "    -t <tag> tag/stem with which to save log and complex gains" << std::endl;
	std::cout << "    -m <model> path to model file" << std::endl;
	std::cout << "    -l <tag> tag/stem of the full/leakage complex gains, solved using unpolarized scan " << std::endl;
	std::cout << std::endl;
}

int main(int argc, char *argv[]) {

	/* hello getopt, my old friend */
	int opt;
	int cal_scan_number;
	std::string tag;
	std::string leakage_tag;
	std::string model_path;
	std::string lta_path;

	if ( argc < 2 ) {
		print_help ();
		exit (EXIT_SUCCESS);
	}

	while ( (opt = getopt ( argc, argv, "hs:t:m:l:" )) != -1 ) {
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
			case 'l':
				leakage_tag  = optarg;
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
	const std::string path_leakage_r  = leakage_tag + std::string("_drl.gains");
	const std::string path_leakage_l  = leakage_tag + std::string("_dlr.gains");

	const std::string save_file_rr  = tag + std::string("_rr.gains");
	const std::string save_file_rl  = tag + std::string("_rl.gains");
	const std::string save_file_lr  = tag + std::string("_lr.gains");
	const std::string save_file_ll  = tag + std::string("_ll.gains");
	const std::string save_file_xp  = tag + std::string("_xp.gains");
	const std::string save_file_gr  = tag + std::string("_gr.gains");
	const std::string save_file_gl  = tag + std::string("_gl.gains");
	const std::string log_file      = tag + std::string(".log");

	const std::string save_self_rr  = tag + std::string("_rr.self");
	const std::string save_self_ll  = tag + std::string("_ll.self");

	const std::string save_file_par_r  = tag + std::string("_pr.gains");
	const std::string save_file_par_l  = tag + std::string("_pl.gains");

	std::cout << "[inputs] lta="   << lta_path << std::endl;
	std::cout << "[inputs] model=" << model_path << std::endl;
	std::cout << "[inputs] tag="   << tag  << std::endl;
	std::cout << "[inputs] leakage_tag="   << leakage_tag  << std::endl;
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
	/*      READ LEAKAGE FILE              */
	/***************************************/
	gaintable::gaintable_t leakage_r  = gaintable::read_complex_solutions ( path_leakage_r );
	gaintable::gaintable_t leakage_l  = gaintable::read_complex_solutions ( path_leakage_l );

	/***************************************/
	/*      AVERAGE LTA SCAN               */
	/***************************************/
	/* maybe directly read complex<float>  */
  /* NOTICE:avgbldata contains self-terms*/
	LTA::vf_type avgbldata ( 2 * nbaselines * nchannels, 0.0f );

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
	/* only self_baseline index */
	std::vector<int>     self_baselines;
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

		/* save self baselines */
		if ( ant1 == ant2 ) self_baselines.push_back ( ib );

	} /* polar baselines */

	const int nantennas   (ant2idx.size()); 

	/***************************************/
	/*        WRITE SELF CORRELATIONS      */ 
	/***************************************/
	selftable::selftable_t    self_rr = selftable::prepare_selftables ( nchannels );
	selftable::selftable_t    self_ll = selftable::prepare_selftables ( nchannels );

	const int n_self ( self_baselines.size() );
	if ( n_self != nantennas * 4 ) {
		std::cout << self_baselines.size() << " " << nantennas << std::endl;
		throw std::runtime_error ("self baseline count mismatch");
	}
	for (int ii = 0; ii < n_self; ii++) {

			const int ib     = self_baselines[ii];

			/*
			 * (0,0) = rr -> 0 
			 * (0,1) = rl -> 1 
			 * (1,0) = lr -> 2 
			 * (1,1) = ll -> 3
			*/
			/* get baseline object */
			const baseline_t&  _bl = lta_file.baselines [ ib ];

			/* get antennas */
			const auto& ant1  = _bl.ant1;
			const auto& ant2  = _bl.ant2;

			if ( ant1 != ant2 ) throw std::runtime_error("expected selfcorr");

			/* get bands */
			const auto& band1 = _bl.band1;
			const auto& band2 = _bl.band2;

			/* ID correlation */
			const int pb2corr = band1*2 + band2;

			for (int ichan = 0; ichan < nchannels; ichan++) {
					/* index in (baseline, channel) complex<float> */
					const int _i     = ichan + nchannels * ib;

					/* copy data */
					// only real because rr and ll of selfcorr are purely real
					const float _real ( avgbldata[2*_i] );

					// data
					if ( pb2corr == 0 ) self_rr[ant1][ichan] = _real;
					else if ( pb2corr == 3 ) self_ll[ant1][ichan] = _real;

			} // channel

	} // selfcorr
	
	// write selfcorr
	selftable::write_selftables ( self_rr, save_self_rr );
	selftable::write_selftables ( self_ll, save_self_ll );

	/***************************************/
	/*        MEASURE PAR ANGLES           */ 
	/***************************************/

	/*
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

	gaintable::gaintable_t    solved_gains_xp = gaintable::prepare_gaintables ( nchannels );

	gaintable::gaintable_t    solved_gains_gr = gaintable::prepare_gaintables ( nchannels );
	gaintable::gaintable_t    solved_gains_gl = gaintable::prepare_gaintables ( nchannels );

	gaintable::gaintable_t    used_par_r = gaintable::prepare_gaintables ( nchannels );
	gaintable::gaintable_t    used_par_l = gaintable::prepare_gaintables ( nchannels );

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
				nantennas
		);

		/* initialize gains */
		for (auto _i = ant2idx.begin(); _i != ant2idx.end(); ++_i) {

			const auto& iant = _i->first;
			const auto& idx  = _i->second;

			// load leakages
			pkg.leakrl [ idx ]  = leakage_r[iant][ichan];
			pkg.leaklr [ idx ]  = leakage_l[iant][ichan];

		}

		/* initialize data */
		for ( int ii = 0; ii < n_noself_baselines; ii++ ) {

			/* ensures this is index of noself baseline */
			/* which is consistent with lta_file        */
			/* ib only for reading */
			const int ipb     = noself_baselines[ii];

			const int ib ( ii / 4 );

			/* index in (baseline, channel) complex<float> */
			const int _i     = ichan + nchannels * ipb;

			/* get baseline object */
			const baseline_t&  _bl = lta_file.baselines [ ipb ];

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
			// correlation index
			pkg.pb2corr  [ ii ] = pb2corr;

			// parallactic angle correct model
			// find parallactic angle
			const float _par1 ( antpar.at(ant1) );
			const float _par2 ( antpar.at(ant2) );

			// perform parang correction
			const auto& _par_model = models::parallactic_correction ( 
					// angles
					_par1, _par2, 
					// input model
					model_rr, model_rl, model_lr, model_ll 
			);

			pkg.par_model_rr [ ib ] = _par_model[0];
			pkg.par_model_rl [ ib ] = _par_model[1];
			pkg.par_model_lr [ ib ] = _par_model[2];
			pkg.par_model_ll [ ib ] = _par_model[3];

			// antenna indices
			pkg.iant1 [ ib ]    = iant1;
			pkg.iant2 [ ib ]    = iant2;

		} /* baselines */

		/* perform solving */
#ifdef TIMING
		start  = std::chrono::high_resolution_clock::now();
#endif

		LBFGS::IterativePolarizedDiagonalSolver     solver (nantennas);
		LBFGS::precompute ( &pkg );
		solver.allsolve ( pkg );
		auto cost = solver.solve ( pkg, 8 );

#ifdef TIMING
		end   = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double,std::milli> duration = end - start;
		logger.time   [ ichan ]  = duration.count();
#endif

		logger.sse    [ ichan ]  = cost[0] + cost[1];
		logger.sse_rr [ ichan ]  = cost[0];
		logger.sse_ll [ ichan ]  = cost[1];
		//logger.nfev   [ ichan ]  = solver.niter;
		//logger.info   [ ichan ]  = solver.rcode;
		//logger.gnorm  [ ichan ]  = solver.gnorm;

		/* save into gain table */
		//const LBFGS::real_type par_xphase ( solver.xpar[4*nantennas - 2] );
		//const LBFGS::real_type xphase ( 2.0f * std::atan ( par_xphase ) );
		/*
		 * Should i zerophase by dividing only the grr and gll by reference grr?
		*/
		for (auto _i = ant2idx.begin(); _i != ant2idx.end(); ++_i) {

			const auto& iant = _i->first;
			const auto& idx  = _i->second;

			// fetch gains
			const LBFGS::complex_type grr ( pkg.cgains[2*idx] );
			const LBFGS::complex_type gll ( pkg.cgains[2*idx + 1] );

			/*
			 * You can phase zero the gll, to get nice phases refant phase zero solutions
			*/

			// fetch leakage terms
			const LBFGS::complex_type drl ( leakage_r[iant][ichan] );
			const LBFGS::complex_type dlr ( leakage_l[iant][ichan] );

			const LBFGS::complex_type xp  ( grr / gll );

			// rr
			const LBFGS::complex_type rr ( grr );
			// rl
			const LBFGS::complex_type rl ( grr * drl );
			// lr
			const LBFGS::complex_type lr ( gll * dlr );
			// ll
			const LBFGS::complex_type ll ( gll );

			solved_gains_rr[iant][ichan]  = rr;
			solved_gains_rl[iant][ichan]  = rl;
			solved_gains_lr[iant][ichan]  = lr;
			solved_gains_ll[iant][ichan]  = ll;

			solved_gains_xp[iant][ichan]  = xp;

			solved_gains_gr[iant][ichan]  = grr;
			solved_gains_gl[iant][ichan]  = gll;

			const float _par ( antpar.at(iant) );
			used_par_r[iant][ichan]  = LBFGS::complex_type ( std::cos ( _par ), -1.0f*std::sin( _par ) );
			used_par_l[iant][ichan]  = LBFGS::complex_type ( std::cos ( _par ),  1.0f*std::sin( _par ) );

		} /* ant2idx */

	} // channel

	/***************************************/
	/*        WRITE GAINTABLES             */
	/***************************************/
	gaintable::write_complex_solutions ( solved_gains_rr, save_file_rr );
	gaintable::write_complex_solutions ( solved_gains_rl, save_file_rl );
	gaintable::write_complex_solutions ( solved_gains_lr, save_file_lr );
	gaintable::write_complex_solutions ( solved_gains_ll, save_file_ll );

	gaintable::write_complex_solutions ( solved_gains_xp, save_file_xp );

	gaintable::write_complex_solutions ( solved_gains_gr, save_file_gr );
	gaintable::write_complex_solutions ( solved_gains_gl, save_file_gl );

	gaintable::write_complex_solutions ( used_par_r, save_file_par_r );
	gaintable::write_complex_solutions ( used_par_l, save_file_par_l );

	/***************************************/
	/*        WRITE LOG                    */
	/***************************************/
	logging::write_log ( logger, log_file );

	auto total_end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> total_duration = total_end - total_start;

	std::cout << std::setprecision(3) << std::endl << " Total solving took " << total_duration.count() << " seconds ..." << std::endl;

	return 0;
}



