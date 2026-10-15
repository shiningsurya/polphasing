#pragma once

#include <iostream>
#include <array>
#include <vector>
#include <string>
#include <complex>
#include <map>
#include <fstream>
#include <iomanip>

#include "gaintable.hpp"

namespace selftable {
	using real_type   = float;
	using antname_t   = std::array<char,4>;
	using self_t      = std::vector<real_type>;
	using selftable_t = std::map< antname_t, self_t>;

	static const std::vector<antname_t> sol_ant_order {
			antname_t {"C00"},
			antname_t {"C01"},
			antname_t {"C02"},
			antname_t {"C03"},
			antname_t {"C04"},
			antname_t {"C05"},
			antname_t {"C06"},
			antname_t {"C08"},
			antname_t {"C09"},
			antname_t {"C10"},
			antname_t {"C11"},
			antname_t {"C12"},
			antname_t {"C13"},
			antname_t {"C14"},
			antname_t {"E02"},
			antname_t {"E03"},
			antname_t {"E04"},
			antname_t {"E05"},
			antname_t {"E06"},
			antname_t {"S01"},
			antname_t {"S02"},
			antname_t {"S03"},
			antname_t {"S04"},
			antname_t {"S06"},
			antname_t {"W01"},
			antname_t {"W02"},
			antname_t {"W03"},
			antname_t {"W04"},
			antname_t {"W05"},
			antname_t {"W06"},
			antname_t {"C07"},
			antname_t {"S05"} };

	/* initialize gaintables */
	selftable_t prepare_selftables ( int nchan );

	selftable_t read_selftables ( const std::string& infile );

	int write_selftables ( const selftable_t& st, const std::string& outfile );

}; /* gaintable */


