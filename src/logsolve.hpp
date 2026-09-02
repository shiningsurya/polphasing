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
	
		vr_type  sse;
		vr_type  time;
		vr_type  gnorm;
		vr_type  fitted_i;
		vi_type  nfev;
		vi_type  info;

		log_t ( int nchan ) : 
			sse ( nchan ), time (nchan), gnorm (nchan), fitted_i(nchan), nfev (nchan), info(nchan) {}

	};

	int write_log ( const log_t& log, const std::string& outfile );

}; // logging
