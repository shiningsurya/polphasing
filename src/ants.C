#include "ants.hpp"

double ants::calculate_lst ( double mjd, antpos_t ant ) {
	/*
	 * This code is shamelessly taken from slalib 
	*/

	const double tu ( ( mjd - 51544.5 ) / 36525.0 );

	const double gmst ( 
			std::fmod ( mjd, 1.0 ) * TAU + 
			( 
			 24110.54841 + ( 8640184.812866 +  ( 0.093104 - 6.2e-6 * tu ) * tu ) * tu 
			) * S2R 
	);

	const double lst (std::fmod(gmst + ant.longitude, TAU));

	return lst;
}

ants::ant2par_t ants::calculate_parallactic_angle ( double tmjd, double source_ra_rad, double source_dec_rad ) {

	// return this
	ant2par_t  ret;

	/* get lst */
	/* get lst using phase center coordinates */
	const double lst = calculate_lst ( tmjd, GMRT_POS );

	// hourangle
	const double  ha ( lst - source_ra_rad );

	//std::cout <<  tmjd << " " << lst << " " << ha << std::endl;

	/* iterate over ant2pos and compute angle for every antenna */
	for (auto it = ant2pos.begin(); it != ant2pos.end(); ++it) {

		// fetch antenna name and its position
		const antname_t _ant = it->first;
		const antpos_t  _pos = it->second;

		const double ant_lat_rad  = _pos.latitude;

		const double q       = std::atan2 ( 
							std::sin(ha), 
							std::tan(ant_lat_rad)*std::cos(source_dec_rad) - 
							std::sin(source_dec_rad)*std::cos(ha)
		);

		ret[_ant]  = q;

	}

	return ret;
}

ants::baseline_uv ants::calculate_uvwave ( double tmjd, double source_ra_rad, double source_dec_rad ){

	/*return this*/
	baseline_uv ret;

	/* get hourangle */
	const double lst  = calculate_lst ( tmjd, GMRT_POS );
	const double ha   = lst - source_ra_rad;

	/* precompute cos and sin */
	const double sinha ( std::sin ( ha ) );
	const double cosha ( std::cos ( ha ) );
	const double sindec( std::sin ( source_dec_rad ) );
	const double cosdec( std::cos ( source_dec_rad ) );

	/* iterate over and save */
	for ( auto _it = bl2xyz.begin(); _it != bl2xyz.end(); _it++ ) {

		const auto& _bl   = _it->first;
		const auto& _xyz  = _it->second;

		const double lu ( sinha*_xyz[0] + cosha*_xyz[1] );
		const double lv ( -sindec*cosha*_xyz[0] + sindec*sinha*_xyz[1] + cosdec*_xyz[2] );

		ret[_bl]  = std::sqrt ( lu*lu + lv*lv );

	} // iteration

	return ret;
}
