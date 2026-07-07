#include "logsolve.hpp"

int logging::write_log (const log_t& gt, const std::string& outfile) {
	std::ofstream of ( outfile );

	/* assume nchan is the same */
	int nchan = gt.at(sol_ant_order[0]).size();

	fmt::print (of, "{: ^14} {: ^14}\n", "sse", "time_ms");

	for (int ichan = 0; ichan < nchan; ichan++) {
		fmt::print (of, "{:6.3f} {:6.3f}\n", log.sse_chan[ichan], log.time_chan[ichan]);
	}

	return 0;
}
