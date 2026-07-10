#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iomanip>

/* fmt library */
//#define FMT_HEADER_ONLY
//#include "fmt/base.h"
//#include "fmt/core.h"
//#include "fmt/format.h"
//#include "fmt/ostream.h"

namespace logging {
	using real_type = float;
	using vr_type   = std::vector<real_type>;
	using vi_type   = std::vector<int>;

	struct log_t {
	
		vr_type  sse_chan;
		vr_type  time_chan;
		vi_type  nfev;
		vi_type  njev;
		vi_type  info;

		log_t ( int nchan ) : 
			sse_chan ( nchan ), time_chan (nchan), nfev (nchan), njev(nchan), info(nchan) {}

	};

	int write_log ( const log_t& log, const std::string& outfile );

}; // logging
