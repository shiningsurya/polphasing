#pragma once

#include <string>
#include <vector>
#include <fstream>

#include "lta_file.hpp"

struct custom_file {
	using vf_type = std::vector<float>;
	/**
	 * defining everything already 
	 * because this is custom file
	 *
	 * methods and members are defined
	 * to make this drop-in replacement for lta_file
	**/
	static constexpr int     nbaselines {1984};
	static constexpr int     nchannels {2048};
	static constexpr float   fbw { 97656.25 };
	static constexpr float   fedge { 550000000 };
	static constexpr double  mjd {61122.454626195824};

	using scan_t   = LTA::scan_t;
	const scan_t  scan { "3C138", 1.3887374,  0.28956375, mjd };

	//const std::string  custom_tag {"/tmp/hc/3C138_unphased.custom"};
	const std::string  custom_tag {"/wd/POLTEST/TST3243/HarshaCorrelator/hc/3C138_unphased.custom"};

	/* this is just dummy */
	scan_t get_scan( int iscan = 0 ) { return scan; }

	std::vector<base_t>     baselines;

	/* time average one scan */
	/* iscan here is dummy */
	void time_average ( int iscan, vf_type& inout );

	/* populate everything in ctor */
	custom_file();

};
