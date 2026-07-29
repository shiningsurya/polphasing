#pragma once

#include <string>
#include <utility>
#include <cmath>
#include <map>

namespace sources {
	/*
	 * Because GMRT source coordinates are in a precessing system, 
	 * those given in the headers of data products cannot be directly used.
	 *
	 * We provide the coordinates here.
	*/

	// always radians
	using coords_t  = std::pair<double,double>;
	using source_t  = std::string;
	using sources_t = std::map<source_t,coords_t>;

	static constexpr double  D2R ( M_PI / 180 );

	static const sources_t table {
		{"3C138",{ 79.5687917*D2R, 16.5907806*D2R }},
		{"3C147",{ 84.6812917*D2R, 49.8285556*D2R }},
		{"3C48",{ 24.4220417*D2R, 33.1597417*D2R }},
	};

}; // namespace
