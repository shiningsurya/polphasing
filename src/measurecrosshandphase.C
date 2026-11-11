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


#ifdef TIMING
auto start  = std::chrono::high_resolution_clock::now();
auto end    = std::chrono::high_resolution_clock::now();
#endif 

void print_help () {
	std::cout << "measurecrosshandphase" << std::endl;
	std::cout << "  Measures crosshand phase using IQU model and solved complex gains" << std::endl;
	std::cout << std::endl;
	std::cout << " measurecrosshandphase [ARGUMENTS] LTA_FILE" << std::endl;
	std::cout << "    -h Print help" << std::endl;
	std::cout << "    -s <scan> scan number of the LTA file" << std::endl;
	std::cout << "    -m <model> path to model file" << std::endl;
	std::cout << "    -t <tag> tag/stem with which to save crosshandphase gains" << std::endl;
	std::cout << "    -l <tag> tag/stem of the full (rr,rl,lr,ll) complex gains, solved using unpolarized scan " << std::endl;
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
	const std::string path_gain_rr    = leakage_tag + std::string("_rr.gains");
	const std::string path_gain_rl    = leakage_tag + std::string("_rl.gains");
	const std::string path_gain_lr    = leakage_tag + std::string("_lr.gains");
	const std::string path_gain_ll    = leakage_tag + std::string("_ll.gains");

	const std::string save_file_xp  = tag + std::string("_xp.gains");

	std::cout << "[inputs] lta="   << lta_path << std::endl;
	std::cout << "[inputs] model=" << model_path << std::endl;
	std::cout << "[inputs] tag="   << tag  << std::endl;
	std::cout << "[inputs] solved gains="   << leakage_tag  << std::endl;
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

	/* gain tables */
	gaintable::gaintable_t    solved_gains_xp = gaintable::prepare_gaintables ( nchannels );

	/* main loop */
	std::cout << " Starting main solving loop" << std::endl;

	auto total_start = std::chrono::high_resolution_clock::now();
	
#ifdef CHANDEBUG
	/* channel 458 has high loss.*/
	for (int ichan = 501; ichan < 502; ichan++) {
#else
//	#pragma omp parallel for num_threads(4) 
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

			// load gains
			pkg.cgains[4*idx]     = g_rr[iant][ichan];
			pkg.cgains[4*idx+1]   = g_rl[iant][ichan];
			pkg.cgains[4*idx+2]   = g_lr[iant][ichan];
			pkg.cgains[4*idx+3]   = g_ll[iant][ichan];

		}

		/* initialize data */
		for ( int ii = 0; ii < n_noself_baselines; ii++ ) {

			/* ensures this is index of noself baseline */
			/* which is consistent with lta_file        */
			/* ib only for reading */
			const int ipb     = noself_baselines[ii];

			const int ib ( ipb / 4 );

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

			// antenna indices
			pkg.iant1 [ ib ]    = iant1;
			pkg.iant2 [ ib ]    = iant2;

			/* copy data */
			const float _real ( avgbldata[2*_i] );
			const float _imag ( avgbldata[2*_i + 1] );

			const int ji ( 4*ib + pb2corr );

			// data
			pkg.data [ ji ]     = LBFGS::complex_type (  _real,  _imag );
			// correlation index
			pkg.pb2corr  [ ji ] = pb2corr;

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

		} /* polarbaselines */

		/* summations */
		complex_type arl (0.0, 0.0), brl (0.0, 0.0);
		complex_type alr (0.0, 0.0), blr (0.0, 0.0);

		/* invert data */
		for (int ibl = 0; ibl < pkg.nbaselines; ibl++) {

			// fetch the antenna index
			const int iant1 ( pkg.iant1[ibl] );
			const int iant2 ( pkg.iant2[ibl] );

			// polarbaseline index
			const int pbl ( 4*ibl );

			// fetch complex data
			const complex_type dpqrr ( pkg.data[ pbl    ] );
			const complex_type dpqrl ( pkg.data[ pbl + 1] );
			const complex_type dpqlr ( pkg.data[ pbl + 2] );
			const complex_type dpqll ( pkg.data[ pbl + 3] );

			const complex_type gprr ( pkg.cgains[4*iant1    ] );
			const complex_type gprl ( pkg.cgains[4*iant1 + 1] );
			const complex_type gplr ( pkg.cgains[4*iant1 + 2] );
			const complex_type gpll ( pkg.cgains[4*iant1 + 3] );

			const complex_type cgqrr ( conj ( pkg.cgains[4*iant2    ] ) );
			const complex_type cgqrl ( conj ( pkg.cgains[4*iant2 + 1] ) );
			const complex_type cgqlr ( conj ( pkg.cgains[4*iant2 + 2] ) );
			const complex_type cgqll ( conj ( pkg.cgains[4*iant2 + 3] ) );

			/* model */
			const complex_type mrr ( pkg.par_model_rr[ibl] );
			const complex_type mrl ( pkg.par_model_rl[ibl] );
			const complex_type mlr ( pkg.par_model_lr[ibl] );
			const complex_type mll ( pkg.par_model_ll[ibl] );

			const complex_type cmrr ( conj ( mrr ) );
			const complex_type cmrl ( conj ( mrl ) );
			const complex_type cmlr ( conj ( mlr ) );
			const complex_type cmll ( conj ( mll ) );

 const complex_type x0 ( cgqrl ) ;
 const complex_type x1 ( gpll*gprr ) ;
 const complex_type x2 ( cgqll ) ;
 const complex_type x3 ( cgqrr ) ;
 const complex_type x4 ( x2*x3 ) ;
 const complex_type x5 ( cgqlr ) ;
 const complex_type x6 ( x0*x5 ) ;
 const complex_type x7 ( gplr*gprl ) ;
 const complex_type x8 ( 1.0/(x1*x4 - x1*x6 - x4*x7 + x6*x7) ) ;
 const complex_type x9 ( gprl*x8 ) ;
 const complex_type x10 ( gpll*x8 ) ;
 const complex_type x11 ( gprr*x8 ) ;
 const complex_type x12 ( gplr*x8 ) ;

 			// rr
 			//pkg.cata[pbl] = dpqll*x0*x9 - dpqlr*x2*x9 - dpqrl*x0*x10 + dpqrr*x10*x2 ;
 			//const complex_type irr ( dpqll*x0*x9 - dpqlr*x2*x9 - dpqrl*x0*x10 + dpqrr*x10*x2 ) ;

			// rl
			const complex_type irl ( -dpqll*x3*x9 + dpqlr*gprl*x5*x8 + dpqrl*gpll*x3*x8 - dpqrr*x10*x5 ) ;

			// lr
			const complex_type ilr ( -dpqll*x0*x11 + dpqlr*gprr*x2*x8 + dpqrl*gplr*x0*x8 - dpqrr*x12*x2 ) ;

			// ll
			//const complex_type ill ( dpqll*x11*x3 - dpqlr*x11*x5 - dpqrl*x12*x3 + dpqrr*x12*x5 ) ;

			// accumulation
			arl += irl * cmrl;
			brl += conj ( irl ) * mrl;

			alr += ilr * cmlr;
			blr += conj ( ilr ) * mlr;

		} /* baselines */

		const complex_type xp2 ( ( arl + blr ) / ( alr + brl ) );

		const real_type phi ( 0.5 * std::arg ( xp2 ) );

		const complex_type xp ( std::cos ( phi ), std::sin ( phi ) );

		for (auto _i = ant2idx.begin(); _i != ant2idx.end(); ++_i) {

			const auto& iant = _i->first;
			const auto& idx  = _i->second;

			solved_gains_xp[iant][ichan]  = xp;

		} /* ant2idx */

	} // channel

	/***************************************/
	/*        WRITE GAINTABLES             */
	/***************************************/
	gaintable::write_complex_solutions ( solved_gains_xp, save_file_xp );

	return 0;
}
