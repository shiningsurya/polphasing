#include <iostream>
#include <fstream>

#include "objectives.hpp"
#include "lta_file.hpp"

using data_t = polphasing::polphasing_data_t;

int main() {
	LTA test ( "/tmp/lta/relta_out.lta" );
	std::cout << test;

	LTA::vf_type bl ( 2 * test.nbaselines * test.nchannels, 0. );

	test.time_average ( 5, bl );

	data_t  pkg ( test.nbaselines );

	int ichan     = 500;
	float model_i = 1.0;
	float model_q = 0.5;
	float model_u = 0.5;

	std::complex<float> model_rr ( model_i, 0.0f );
	std::complex<float> model_rl ( model_q, model_u );
	std::complex<float> model_lr ( model_q,-model_u );
	std::complex<float> model_ll ( model_i, 0.0f );

	/*
	 * Running list
	 */
	std::map<LTA::antname_t,int> ant2idx;
	int ir = 0;

	for ( int ib = 0; ib < test.nbaselines; ib++ ) {

		int _i    = ichan + test.nchannels*ib;

		auto  _bl = test.baselines [ ib ];

		pkg.data [ ib ] = std::complex<float>( bl[2*_i], bl[2*_i + 1] );

		/* every baseline is polar baseline */
		/* although redundant, we write a model column */
		/* because it makes future computations straightforward */
		if ( _bl.band1 == 0 ) {
			if ( _bl.band2 == 0 ) {
				/* RR */
				pkg.model [ib] = model_rr;
			} else {
				/* RL */
				pkg.model [ib] = model_rl;
			}
		} else {
			if ( _bl.band2 == 0 ) {
				/* LR */
				pkg.model [ib] = model_lr;
			} else {
				/* LL */
				pkg.model [ib] = model_ll;
			}
		}

		/* populate index_b{1,2} */
		/* this is done on the fly */
		auto ant1  = _bl.ant1;
		auto it1   = ant2idx.find(ant1);
		if ( it1 == ant2idx.end() ) {
			ant2idx.emplace ( ant1, ir );
			pkg.index_b1 [ib] = ir;
			ir++;
		} 
		else {
			pkg.index_b1 [ib] = it1->second;
		}

		auto ant2  = _bl.ant2;
		auto it2   = ant2idx.find(ant2);
		if ( it2 == ant2idx.end() ) {
			ant2idx.emplace ( ant2, ir );
			pkg.index_b2 [ib] = ir;
			ir++;
		} 
		else {
			pkg.index_b2 [ib] = it2->second;
		}

	}

	/* strong sanity check */
	if ( ir != test.nantennas ) {
		throw std::runtime_error ("Antenna mapping mismatch.");
	}

	std::cout << " break here " << std::endl;

	{
		std::ofstream of("pkg.bin", std::ios::binary);

		/* data */
		of.write (reinterpret_cast<const char*>(pkg.data.data()), pkg.data.size()*sizeof(polphasing::complex_type));

		/* model */
		of.write (reinterpret_cast<const char*>(pkg.model.data()), pkg.model.size()*sizeof(polphasing::complex_type));
		
		/* idxb1  */
		of.write (reinterpret_cast<const char*>(pkg.index_b1.data()), pkg.index_b1.size()*sizeof(int));

		/* idxb2  */
		of.write (reinterpret_cast<const char*>(pkg.index_b2.data()), pkg.index_b1.size()*sizeof(int));
	}

	//{
		//std::ofstream of("bl.data", std::ios::binary);
		//of.write ( reinterpret_cast<const char*>(bl.data()), sizeof(float)*2*test.nbaselines*test.nchannels );
	//}
	//{
		//std::ofstream of("ant2idx.data")
	//}


	return 0;
}
