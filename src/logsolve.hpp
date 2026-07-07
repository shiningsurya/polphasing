#pragma once

#include <vector>
#include <string>
#include <fstream>

/* fmt library */
#define FMT_HEADER_ONLY
#include "fmt/base.h"
#include "fmt/core.h"
#include "fmt/format.h"
#include "fmt/ostream.h"

namespace logging {
	using real_type = float;
	using vr_type   = std::vector<real_type>;

	struct log_t {
	
		vr_type  sse_chan;
		vr_type  time_chan;

		log_t ( int nchan ) : 
			sse_chan ( nchan ), time_chan (nchan) {}

	};

	int write_log ( const log_t& log, const std::string& outfile );

}; // logging
