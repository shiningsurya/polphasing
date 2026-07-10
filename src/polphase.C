#include <iostream>

#include <unistd.h>

#include <string>

#include "lmsolver.hpp"
#include "lta_file.hpp"
#include "gaintable.hpp"
#include "models.hpp"
#include "logsolve.hpp"

#ifdef TIMING
#include <chrono>
#endif

//#include "cxxopts.hpp"

/*
Arguments:
	LTA_file
	model_file
	scan_number

	future work:
		scan_number_def=last scan with object from model_file
		model_file_def=a directory contains all the model files
		code decides on its own
*/

// inputs
//const int cal_scan_number ( 5 );
//const std::string model_file ("/home/shining/shit/gmrt_phase_test/polphasing/code/polphasing/models/3C138_2048_97656.25_550000000_model.txt");
//const std::string lta_file   ("/tmp/lta/relta_out.lta");
//const std::string save_file_r("solvedO2_pfreset_r.gains");
//const std::string save_file_l("solvedO2_pfreset_l.gains");
//const std::string log_file   ("solvedO2_pfreset.log");

using FullPolarLMSolver = LMSolver<LMSolverType::FULL_POLAR>;
/***************************************/
/* NEED TO REMOVE SELF-TERMS           */
/***************************************/

#ifdef TIMING
auto start  = std::chrono::high_resolution_clock::now();
auto end    = std::chrono::high_resolution_clock::now();
#endif 

void print_help () {
	std::cout << "polphase" << std::endl;
	std::cout << "  Solving complex gains for each antenna using IQU model" << std::endl;
	std::cout << std::endl;
	std::cout << " polphase [ARGUMENTS] LTA_FILE" << std::endl;
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
				model_path = strdup ( optarg );
				break;
		} // switch
	} // getopt
	/* lta file*/
	lta_path  = argv[optind];
	optind++;

	/* other files */
	const std::string save_file_r  = tag + std::string("_r.gains");
	const std::string save_file_l  = tag + std::string("_l.gains");
	const std::string log_file     = tag + std::string(".log");

	std::cout << "[inputs] lta=" << lta_path << " model=" << model_path << std::endl;
	std::cout << "[inputs] tag=" << tag << " scan=" << cal_scan_number << std::endl;

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
	models::model_data_t calmodel = models::read_model_file ( model_path );

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
	/* (2) polarbaseline to correlation product -- pb2cor  */
	/* (3) (antenna, band) to gain in solution vector -- index_b{1,2} */
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

	/*
	 * ant2idx : _ _ _ _ _ _ .... _ (nant)
	 * 2*ant2idx + band
	 *  :   iant0 iant1 jant0 jant1 ..... zant0 zant1 ( 2*nant )
	 * 
	 * 2*(2*ant2idx + band) = 4*ant2idx + 2*band
	 *  :   iant0_r iant0_i iant1_r iant1_i
	 *
	 */

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

	int ndata        = 2 * n_noself_baselines; /* complex -> real,imag */
	int npar         = 4 * nantennas; /* complex(R) and complex(L) gains */

	/* logging */
	logging::log_t   logger ( nchannels );

	/* gain tables */
	gaintable::gaintable_t    solved_gains_r = gaintable::prepare_gaintables ( nchannels );
	gaintable::gaintable_t    solved_gains_l = gaintable::prepare_gaintables ( nchannels );

	/* data package */
	FullPolarLMSolver::ptrdata_t  pkg ( new FullPolarLMSolver::data_t ( n_noself_baselines ) );
	std::vector<int>     pb2corr (n_noself_baselines);

	/*
	 * This loop prepares pb2corr
	 */

	for ( int ii = 0; ii < n_noself_baselines; ii++ ) {

		int ib    = noself_baselines[ii];

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
		 */
		pb2corr[ii]    = _bl.band1*2 + _bl.band2;

		/* get antennas */
		const auto& ant1  = _bl.ant1;
		const auto& ant2  = _bl.ant2;

		/* get bands */
		const auto& band1 = _bl.band1;
		const auto& band2 = _bl.band2;

		/* index_b1 b2 */
		const auto& iant1 = ant2idx.at(ant1);
		const auto& iant2 = ant2idx.at(ant2);

		int ix1   = iant1*2 + band1;
		int ix2   = iant2*2 + band2;

		pkg->index_b1 [ ii ] = ix1;
		pkg->index_b2 [ ii ] = ix2;

	} // no self baselines 


	/* main loop */
	std::cout << " Starting main solving loop" << std::endl;

	auto total_start = std::chrono::high_resolution_clock::now();
	
  // openmp parallelizing the whole thing 
  // takes 440s or 7 minutes ish
  #pragma omp parallel for num_threads(4)
	for (int ichan = 0; ichan < nchannels; ichan++) {
	//for (int ichan = 0; ichan < 32; ichan++) {
		/* initialize solver */
		FullPolarLMSolver          solver (ndata, npar);
		const FullPolarLMSolver::vr_type& isol = solver.isolution;
	/* testing */
	//for (int ichan = 500; ichan < 501; ichan++) {
		//if (ichan % 128 == 0) std::cout << ichan << " ";

		/* read stokes IQU for ichan */
		models::real_type    stokes_i ( calmodel.stokes_i[ichan] );
		models::real_type    stokes_q ( calmodel.stokes_q[ichan] );
		models::real_type    stokes_u ( calmodel.stokes_u[ichan] );

		/* populate rr, rl, lr, ll */
		models::complex_type model_rr ( stokes_i, 0.0f );
		models::complex_type model_rl ( stokes_q, stokes_u );
		models::complex_type model_lr ( stokes_q,-stokes_u );
		models::complex_type model_ll ( stokes_i, 0.0f );

		/* initialize data */
		for ( int ii = 0; ii < n_noself_baselines; ii++ ) {

			/* ensures this is index of noself baseline */
			/* which is consistent with lta_file        */
			/* ib only for reading */
			int ib     = noself_baselines[ii];

			/* index in (baseline, channel) complex<float> */
			int _i     = ichan + nchannels * ib;

			/*
			 * NEED TO IGNORE SELF-TERMS
			 * nbaselines <- noself_nbaselines
			 *
			 * noself_baselines = vector<int>
			 * then iterate over this vector<int>
			 */

			/* copy data */
			pkg->data [ ii ] = std::complex<float>( avgbldata[2*_i], avgbldata[2*_i + 1] );

			/* copy model using polarbaseline to correlation product mapping */
			if      ( pb2corr[ii] == 0 ) pkg->model[ii] = model_rr;
			else if ( pb2corr[ii] == 1 ) pkg->model[ii] = model_rl;
			else if ( pb2corr[ii] == 2 ) pkg->model[ii] = model_lr;
			else if ( pb2corr[ii] == 3 ) pkg->model[ii] = model_ll;

		} /* baselines */

#if 0	
		/* save pkg and exit */
		{
			std::ofstream bof("avgbldata.bin", std::ios::binary);
			bof.write (reinterpret_cast<const char*>(avgbldata.data()), avgbldata.size()*sizeof(float));

			std::ofstream nof("noself_baselines.bin", std::ios::binary);
			nof.write (reinterpret_cast<const char*>(noself_baselines.data()), noself_baselines.size()*sizeof(int));

			std::ofstream of("pkgchan500.bin", std::ios::binary);

			/* data */
			of.write (reinterpret_cast<const char*>(pkg->data.data()), pkg->data.size()*sizeof(polphasing::complex_type));
			// data is read correctly

			/* model */
			of.write (reinterpret_cast<const char*>(pkg->model.data()), pkg->model.size()*sizeof(polphasing::complex_type));
			// model is read correctly

			/* idxb1  */
			of.write (reinterpret_cast<const char*>(pkg->index_b1.data()), pkg->index_b1.size()*sizeof(int));

			/* idxb2  */
			of.write (reinterpret_cast<const char*>(pkg->index_b2.data()), pkg->index_b2.size()*sizeof(int));
			// they seem okay too!
		}
#endif

		 solver.reset ();

		/* perform solving */
#ifdef TIMING
		start  = std::chrono::high_resolution_clock::now();
#endif
		solver.solve ( pkg );

#if 0	
		{
			/* before writing, update residual and jacobian by call fcn */
			/* with 1 to compute residuals */
			test.forward ( pkg, 1 );
			/* with 2 to compute jacobian */
			test.forward ( pkg, 2 );

			std::ofstream of("solverstate500.bin", std::ios::binary);
			of.write (reinterpret_cast<const char*>(test.isolution.data()), test.isolution.size()*sizeof(float));

			of.write (reinterpret_cast<const char*>(test.residuals.data()), test.residuals.size()*sizeof(float));

			of.write (reinterpret_cast<const char*>(test.jacobian.data()), test.jacobian.size()*sizeof(float));
		}
#endif

#ifdef TIMING
		end   = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double,std::milli> duration = end - start;
		//std::cout << " duration = " << duration.count() << " ms" << std::endl;
		logger.time_chan [ ichan ] = duration.count();
#endif
		//std::cout << "after solving SSE=" << test.get_sse() << std::endl;
		logger.sse_chan [ ichan ]  = solver.get_sse();

		logger.nfev [ ichan ]      = solver.nfev;
		logger.njev [ ichan ]      = solver.njev;

		logger.info [ ichan ]      = solver.info;

		/* save into gain table */
		for (auto _i = ant2idx.begin(); _i != ant2idx.end(); ++_i) {

			const auto& iant = _i->first;
			const auto& idx  = _i->second;

			polphasing::complex_type rg ( isol[4*idx + 0], isol[4*idx + 1] );
			polphasing::complex_type lg ( isol[4*idx + 2], isol[4*idx + 3] );

			solved_gains_r[iant][ichan]  = rg;
			solved_gains_l[iant][ichan]  = lg;

		} /* ant2idx */

	} // channel

	/***************************************/
	/*        WRITE GAINTABLES             */
	/***************************************/
	gaintable::write_complex_solutions ( solved_gains_r, save_file_r );
	gaintable::write_complex_solutions ( solved_gains_l, save_file_l );

	/***************************************/
	/*        WRITE LOG                    */
	/***************************************/
	logging::write_log ( logger, log_file );

	auto total_end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> total_duration = total_end - total_start;

	std::cout << std::setprecision(3) << std::endl << " Total solving took " << total_duration.count() << " seconds ..." << std::endl;

	return 0;
}
