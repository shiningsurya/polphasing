#pragma once

#include <cmath>
#include <array>
#include <map>

#include "gaintable.hpp"

/*
 * Defines antenna positions
 *
*/

namespace ants {

	using antname_t   = std::array<char,4>;

	static constexpr double  R2D ( 180 / M_PI );
	static constexpr double  D2R ( M_PI / 180 );
	static constexpr double  TAU ( 2.0 * M_PI );
	static constexpr double  S2R ( M_PI / 12.0 / 3600.0 );

	struct antenna_position_t {
		/*
		 * In radians
		*/
		double   longitude;
		double   latitude;

		antenna_position_t ( double _long_deg, double _lat_deg ) : 
			longitude (D2R*_long_deg), 
			latitude  (D2R*_lat_deg)
		{}
	};

	using antpos_t    = struct antenna_position_t;

	// LONGITUDE_DEG	74.05656115769753
	// LATITUDE_DEG	19.09300278307055
	// taken from :gmrt_coordinates.py:
	static const antname_t GMRT_ARR {"GMR"};
	static const antpos_t  GMRT_POS (74.05656115769753, 19.09300278307055);

	static const std::map<antname_t,antpos_t> ant2pos {
		{ GMRT_ARR, GMRT_POS },
		{ antname_t{"C00"}, antpos_t(74.052103, 19.090998) },
		{ antname_t{"C01"}, antpos_t(74.051102, 19.091846) },
		{ antname_t{"C02"}, antpos_t(74.050371, 19.093129) },
		{ antname_t{"C03"}, antpos_t(74.049865, 19.095370) },
		{ antname_t{"C04"}, antpos_t(74.049359, 19.095832) },
		{ antname_t{"C05"}, antpos_t(74.049825, 19.090767) },
		{ antname_t{"C06"}, antpos_t(74.049638, 19.091278) },
		{ antname_t{"C08"}, antpos_t(74.049909, 19.088808) },
		{ antname_t{"C09"}, antpos_t(74.050036, 19.091676) },
		{ antname_t{"C10"}, antpos_t(74.048192, 19.088427) },
		{ antname_t{"C11"}, antpos_t(74.047859, 19.092015) },
		{ antname_t{"C12"}, antpos_t(74.048915, 19.087037) },
		{ antname_t{"C13"}, antpos_t(74.045332, 19.085103) },
		{ antname_t{"C14"}, antpos_t(74.047239, 19.088940) },
		{ antname_t{"E02"}, antpos_t(74.060896, 19.093556) },
		{ antname_t{"E03"}, antpos_t(74.068770, 19.097195) },
		{ antname_t{"E04"}, antpos_t(74.080130, 19.096643) },
		{ antname_t{"E05"}, antpos_t(74.087706, 19.093641) },
		{ antname_t{"E06"}, antpos_t(74.096171, 19.098869) },
		{ antname_t{"S01"}, antpos_t(74.043412, 19.066610) },
		{ antname_t{"S02"}, antpos_t(74.036138, 19.056457) },
		{ antname_t{"S03"}, antpos_t(74.031285, 19.035723) },
		{ antname_t{"S04"}, antpos_t(74.024778, 19.011272) },
		{ antname_t{"S06"}, antpos_t(74.007453, 18.976143) },
		{ antname_t{"W01"}, antpos_t(74.048053, 19.102865) },
		{ antname_t{"W02"}, antpos_t(74.046683, 19.114451) },
		{ antname_t{"W03"}, antpos_t(74.045847, 19.133480) },
		{ antname_t{"W04"}, antpos_t(74.047829, 19.157845) },
		{ antname_t{"W05"}, antpos_t(74.053937, 19.185155) },
		{ antname_t{"W06"}, antpos_t(74.049341, 19.203848) },
		{ antname_t{"C07"}, antpos_t(74.049341, 19.203848) },
		{ antname_t{"S05"}, antpos_t(74.049341, 19.203848) },
	//};
	};
/*
* C07 and S05 are dummy antennas. They are forever flagged. The antsys.hdr shows dummy coordiantes. so i am also giving dummy coordinates here.
*/

	using ant2par_t   = std::map<antname_t,float>;

	ant2par_t calculate_parallactic_angle ( double tmjd, double source_ra_rad, double source_dec_rad );

	double calculate_lst ( double tmjd, antpos_t ant );

}; // namespace 
	
std::ostream& operator<< (std::ostream& os, const ants::antname_t& a );
