#include "logsolve.hpp"

int logging::write_log (const log_t& log, const std::string& outfile) {
	std::ofstream of ( outfile );

	/* assume nchan is the same */
	int nchan = log.sse.size();

	//fmt::print (of, "{: ^14} {: ^14} {: ^14} {: ^14} {: ^14}\n", "sse", "time_ms", "nfev", "njev", "info");
	of << "time_ms sse gnorm fitted_I nfev info" << std::endl;

	of << std::fixed << std::setprecision(3) << std::setfill('0') << std::setw(3);
	

	for (int ichan = 0; ichan < nchan; ichan++) {
		of << log.time[ichan] << " " << log.sse[ichan]  << " " << log.gnorm[ichan] << " " << log.fitted_i[ichan] << " " << log.nfev[ichan] << " " << log.info[ichan] << std::endl;
	}

	return 0;
}
