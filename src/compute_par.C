#include <iostream>

#include "ants.hpp"


int main() {

	// manual testing now
	// will write an CLI later
	//const double mjd ( 61122.78 );
	const double mjd ( 61250.29252716616 );

	// 3C138
	const double source_ra  ( ants::D2R*79.5687917 );
	const double source_dec ( ants::D2R*16.5907806 );

	std::cout << "At MJD " << mjd << std::endl;

	// calculate
	
	auto ret  = ants::calculate_parallactic_angle ( mjd, source_ra, source_dec );

	for (auto it = ret.begin(); it != ret.end(); ++it) {
		std::cout << it->first << " " << ants::R2D*it->second << std::endl;
	}
	
	return 0;
}
