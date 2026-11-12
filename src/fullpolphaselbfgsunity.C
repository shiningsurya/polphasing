#include <iostream>
#include <cmath>

#include <unistd.h>

#include <string>

#include "lta_file.hpp"
#include "gaintable.hpp"
#include "models.hpp"
#include "logsolve.hpp"
#include "ants.hpp"

//#include "LBFGSclean.hpp"
#include "LBFGSsuperiterative.hpp"

#ifdef CUSTOM
#include "custom_file.hpp"
#endif

#ifdef TIMING
#include <chrono>
#endif

//#define CHANDEBUG

//#define IPRINT

#ifdef TIMING
auto start  = std::chrono::high_resolution_clock::now();
auto end    = std::chrono::high_resolution_clock::now();
#endif 

void print_help () {
	std::cout << "fullpolphaseunity" << std::endl;
	std::cout << "  Solving full Jones matrix for each antenna using IQU model using multithreading " << std::endl;
	std::cout << std::endl;
	std::cout << " fullpolphase [ARGUMENTS] LTA_FILE" << std::endl;
	std::cout << "    -h Print help" << std::endl;
	std::cout << "    -s <scan> scan number of the LTA file" << std::endl;
	std::cout << "    -m <model> path to model file" << std::endl;
	std::cout << "    -t <tag> tag/stem with which to save log and complex gains" << std::endl;
	std::cout << std::endl;
}

int main(int argc, char *argv[]) {

	/* hello getopt, my old friend */
	int opt;
	int cal_scan_number;
	std::string tag;
	std::string lta_path;
	std::string model_path;

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

	const std::string save_file_leakage_r  = tag + std::string("_drl.gains");
	const std::string save_file_leakage_l  = tag + std::string("_dlr.gains");

	const std::string save_file_par_r  = tag + std::string("_pr.gains");
	const std::string save_file_par_l  = tag + std::string("_pl.gains");

	std::cout << "[inputs] lta="   << lta_path << std::endl;
	std::cout << "[inputs] tag="   << tag  << std::endl;
	std::cout << "[inputs] scan="  << cal_scan_number << std::endl;
	std::cout << "[inputs] model=" << model_path << std::endl;

	/***************************************/
	/*      READ MODEL FILE                */
	/***************************************/
	models::model_data_t calmodel = models::read_model_file ( model_path );

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
	/*      AVERAGE LTA SCAN               */
	/***************************************/
	/* maybe directly read complex<float>  */
  /* NOTICE:avgbldata contains self-terms*/
	LTA::vf_type avgbldata ( 2 * nbaselines * nchannels, 0. );

	const LTA::scan_t cal_scan = lta_file.get_scan ( cal_scan_number );

	lta_file.time_average  ( cal_scan_number, avgbldata );

	/***************************************/
	/*        UVCUT  PREPARE               */
	/***************************************/

	const ants::baseline_uv buv = ants::calculate_uvwave ( cal_scan.mjd, cal_scan.ra, cal_scan.dec, calmodel.fedge );

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

	gaintable::gaintable_t    solved_leakage_drl = gaintable::prepare_gaintables ( nchannels );
	gaintable::gaintable_t    solved_leakage_dlr = gaintable::prepare_gaintables ( nchannels );

	gaintable::gaintable_t    used_par_r = gaintable::prepare_gaintables ( nchannels );
	gaintable::gaintable_t    used_par_l = gaintable::prepare_gaintables ( nchannels );

	/* parangle is same for all channels */

	/* main loop */
	std::cout << " Starting main solving loop" << std::endl;

	auto total_start = std::chrono::high_resolution_clock::now();
	
#ifdef CHANDEBUG
	/* channel 458 has high loss.*/
	for (int ichan = 500; ichan < 501; ichan++) {
#else
	//#pragma omp parallel for num_threads(4) 
	for (int ichan = 0; ichan < nchannels; ichan++) {
#endif
		/* when parallelizing inside loop */
		/* this will be doing a lot of mallocs */

		/* data package */
		LBFGS::uata_t  pkg ( 
				n_noself_baselines, 
				nantennas,
				calmodel.stokes_i[ichan]
		);

		/* initialize data */
		for ( int ii = 0; ii < n_noself_baselines; ii++ ) {

			/* ensures this is index of noself baseline */
			/* which is consistent with lta_file        */
			/* ipb only for reading */
			const int ipb    = noself_baselines[ii];

			/* actual baseline */
			const int ib  ( ii / 4 );

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

			// correlation index
			pkg.pb2corr  [ ii ] = pb2corr;

			// data
			pkg.data [ ii ]   = LBFGS::complex_type (  _real,  _imag );

			// parallactic angle correct model
			// find parallactic angle
			
			const float _par1 ( antpar.at(ant1) );
			const float _par2 ( antpar.at(ant2) );

			// z1, z2
			pkg.par_za [ ib ]   = LBFGS::complex_type ( std::cos ( _par1-_par2 ), std::sin ( _par1-_par2 ) );
			pkg.par_zb [ ib ]   = LBFGS::complex_type ( std::cos ( _par2-_par1 ), std::sin ( _par2-_par1 ) );

			// antenna indices
			pkg.iant1 [ ib ]    = iant1;
			pkg.iant2 [ ib ]    = iant2;

		} /* baselines */

		/* perform solving */
#ifdef TIMING
		start  = std::chrono::high_resolution_clock::now();
#endif

		LBFGS::IterativeUnpolarizedSolver    solver (nantennas);
		solver.allsolve ( pkg );
		auto cost = solver.solve ( pkg, 8 );

#ifdef TIMING
		end   = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double,std::milli> duration = end - start;
		logger.time   [ ichan ]  = duration.count();
#endif

		logger.sse    [ ichan ]  = cost[0] + cost[1] + cost[2] + cost[3];
		logger.sse_rr [ ichan ]  = cost[0];
		logger.sse_rl [ ichan ]  = cost[1];
		logger.sse_lr [ ichan ]  = cost[2];
		logger.sse_ll [ ichan ]  = cost[3];
		//logger.nfev   [ ichan ]  = solver.niter;
		//logger.info   [ ichan ]  = solver.rcode;
		//logger.gnorm  [ ichan ]  = solver.gnorm;

		/* save into gain table */
		/* gainnormalizer 
		 * with first antenna gains
		*/
		const gaintable::normalizer_t    gnormer ( pkg.cgains[0], pkg.cgains[1], pkg.cgains[2], pkg.cgains[3] );
		for (auto _i = ant2idx.begin(); _i != ant2idx.end(); ++_i) {

			const auto& iant = _i->first;
			const auto& idx  = _i->second;

			const LBFGS::complex_type rr ( pkg.cgains[4*idx + 0] );
			const LBFGS::complex_type rl ( pkg.cgains[4*idx + 1] );
			const LBFGS::complex_type lr ( pkg.cgains[4*idx + 2] );
			const LBFGS::complex_type ll ( pkg.cgains[4*idx + 3] );

			/* normalizing */
			const auto ng  ( gnormer ( rr, rl, lr, ll ) );

			//const LBFGS::complex_type drl ( pkg.lgains[2*idx + 0] );
			//const LBFGS::complex_type dlr ( pkg.lgains[2*idx + 1] );

			const LBFGS::complex_type drl ( ng[1] / ng[0] );
			const LBFGS::complex_type dlr ( ng[2] / ng[3] );

			solved_gains_rr[iant][ichan]  = ng[0];
			solved_gains_rl[iant][ichan]  = ng[1];
			solved_gains_lr[iant][ichan]  = ng[2];
			solved_gains_ll[iant][ichan]  = ng[3];

			solved_leakage_drl[iant][ichan] = drl;
			solved_leakage_dlr[iant][ichan] = dlr;

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

	gaintable::write_complex_solutions ( solved_leakage_drl, save_file_leakage_r );
	gaintable::write_complex_solutions ( solved_leakage_dlr, save_file_leakage_l );

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


