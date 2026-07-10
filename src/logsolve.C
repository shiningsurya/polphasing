#include "logsolve.hpp"

int logging::write_log (const log_t& log, const std::string& outfile) {
	std::ofstream of ( outfile );

	/* assume nchan is the same */
	int nchan = log.sse_chan.size();

	//fmt::print (of, "{: ^14} {: ^14} {: ^14} {: ^14} {: ^14}\n", "sse", "time_ms", "nfev", "njev", "info");
	of << "sse" << " " << "time_ms" << " " << "nfev" << " " << "njev" << " " << "info" << std::endl;

	of << std::fixed << std::setprecision(3) << std::setfill('0') << std::setw(3);
	

	for (int ichan = 0; ichan < nchan; ichan++) {
		//fmt::print (of, "{:6.3f} {:6.3f} {:03} {:03} {:03}\n", log.sse_chan[ichan], log.time_chan[ichan], log.nfev[ichan], log.njev[ichan], log.info[ichan]);
		of << log.sse_chan[ichan] << " " << log.time_chan[ichan] << " " << log.nfev[ichan] << " " << log.njev[ichan] << " " << log.info[ichan] << std::endl;
	}

	return 0;
}
