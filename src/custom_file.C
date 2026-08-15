#include "custom_file.hpp"

#include <iostream>

#include <sstream>

custom_file::custom_file() {

	const std::string chdr = custom_tag + std::string(".htxt");

	std::cout << chdr << std::endl;

	static const std::string _rr ("rr");
	static const std::string _rl ("rl");
	static const std::string _lr ("lr");
	static const std::string _ll ("ll");

/**
 * nbaselines rows
 * <ant1> <ant2> <corr>
**/
	std::ifstream ifs (chdr);
	std::vector<std::string>  ltoks;
	std::string line;

	while(std::getline(ifs,line)) {
		/* read line by line */
		ltoks.clear();
		//std::cout << line << std::endl;
		{
			std::string tok;
			std::stringstream ss ( line );
			while ( ss >> tok ) ltoks.push_back ( tok );
		}
		if ( ltoks.size() != 3 ) throw std::runtime_error("Found != 3 entries in the custom_file header.");

		/* emplace_back ant1, ant2, band1, band2 */
		if ( ltoks[2] == _rr ) {
			baselines.emplace_back ( ltoks[0], ltoks[1], 0, 0 );
		}
		else if ( ltoks[2] == _rl ) {
			baselines.emplace_back ( ltoks[0], ltoks[1], 0, 1 );
		}
		else if ( ltoks[2] == _lr ) {
			baselines.emplace_back ( ltoks[0], ltoks[1], 1, 0 );
		}
		else if ( ltoks[2] == _ll ) {
			baselines.emplace_back ( ltoks[0], ltoks[1], 1, 1 );
		}

	} // polar baselines

} // ctor
	
void custom_file::time_average ( int iscan, vf_type& inout ) {
	
	const std::string cdat = custom_tag + std::string(".raw");

	std::ifstream ifs (cdat, std::ios::binary);

	size_t  s = nbaselines*nchannels*2*sizeof(float);

	ifs.read ( reinterpret_cast<char*>(inout.data()), s );

	if ( ifs.gcount() != s ) {
		throw std::runtime_error("incomplete read here.");
	}

}
