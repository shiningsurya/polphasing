#pragma once
/*
 * C++ interface to lta_file
 *
 * reuse code from lute as much as possible
 *
 * Expects consistency between header and corresponding data
 * What if some antenna drops in between scans?
 *
 */

#include <iostream>
#include <array>
#include <vector>
#include <map>
#include <string>
#include <numeric>
#include <algorithm>

#include <cstring>


//extern "C" {
#	include "mylta.h"
//}


struct baseline_t {
	using antname_t = std::array<char,4>;
	int  band1,    band2;
	antname_t ant1,  ant2;
	char samp1[8], samp2[8];

	baseline_t ( const std::string& sant1, const std::string& sant2, const int _band1, const int _band2 ) {
		band1 = _band1;
		band2 = _band2;

		for (int i = 0; i < 4; i++) {
			ant1[i] = sant1[i];
			ant2[i] = sant2[i];
		}
	}

	/* ctor using BaselineType from lta.h */
	baseline_t ( const BaselineType &btype ) {
		band1 = btype.band[0];
		band2 = btype.band[1];

		for (int i = 0; i < 4; i++) {
			ant1[i] = btype.antname[0][i];
			ant2[i] = btype.antname[1][i];
		}

		/* wasted some time in trying to manage 
		 * char[4] to array<char,4>
		 */

		//std::copy ( static_cast<const char*>(&btype.antname[0]), static_cast<const char*>(&btype.antname[0])+4, ant1.data() );
		//std::copy ( static_cast<const char*>(&btype.antname[1]), static_cast<const char*>(&btype.antname[1])+4, ant1.data() );
		
		//std::copy ( &btype.antname[0], &btype.antname[0]+4, ant1.data() );
		//std::copy ( &btype.antname[1], &btype.antname[1]+4, ant2.data() );
		//std::copy ( std::begin(btype.antname[0]), std::end(btype.antname[0]), ant1.begin() );
		//std::copy ( std::begin(btype.antname[1]), std::end(btype.antname[1]), ant2.begin() );
		//strncpy ( ant1, btype.antname[0], 4 );
		//strncpy ( ant2, btype.antname[1], 4 );

		strncpy ( samp1, btype.bandname[0], 8 );
		strncpy ( samp2, btype.bandname[1], 8 );
	}
};



using base_t    = struct baseline_t;

class LTA {
	public:
		using recl_type = unsigned long;

		using antname_t = std::array<char,4>;
		//using antname_t = char[4];

		using vb_type   = std::vector<char>;
		using vi_type   = std::vector<int>;
		using vf_type   = std::vector<float>;

		struct scan_t {
			const std::string source;
			// coordinates in radians
			const double      ra;
			const double      dec;
			const double      mjd;

			scan_t ( const std::string& _source, double _ra, double _dec, double _mjd ) :
				source(_source), ra(_ra), dec(_dec), mjd(_mjd) {}
		};

	private:
		std::string     filepath;

		/* c-like interface */
		FILE       *fp;

		/* LTA related things */
		LtaInfo    linfo;
		recl_type  recl;

		/* one row in LTA file */
		vb_type    drow;

	
	public:
		/* antenna, pol, baseline ordering */
		int        nantennas;
		int        nsamplers;
		int        nbaselines;
		int        nchannels;

		std::vector<antname_t>  all_antnames;
		std::vector<base_t>     baselines;

		/* frequency axis */
		float      fedge;
		float      fbw;
		vf_type    freqs_MHz;
		/* scan table */
		int nscans() const { return linfo.scans; };

		/* time average one scan */
		void time_average ( int iscan, vf_type& inout );

		/* get properties of :iscan: scan */
		scan_t get_scan ( int iscan );

	public:
		/* ctor dtor */
		LTA (const std::string &) ;
		~LTA();

		friend std::ostream& operator<< (std::ostream&, const LTA &);

};
