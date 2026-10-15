#include "selftable.hpp"

/* overload operator<< for my formats */
inline std::ostream& operator<< (std::ostream& os, const selftable::antname_t& a ) {
	os << a[0] << a[1] << a[2] << a[3];
	return os;
}

int selftable::write_selftables (const selftable_t& gt, const std::string& outfile) {

	std::ofstream of ( outfile );

	/* assume nchan is the same */
	int nchan = gt.at(sol_ant_order[0]).size();

	of << 
		sol_ant_order[0] << " " << 
		sol_ant_order[1] << " " << 
		sol_ant_order[2] << " " << 
		sol_ant_order[3] << " " << 
		sol_ant_order[4] << " " << 
		sol_ant_order[5] << " " << 
		sol_ant_order[6] << " " << 
		sol_ant_order[7] << " " << 
		sol_ant_order[8] << " " << 
		sol_ant_order[9] << " " << 
		sol_ant_order[10] << " " << 
		sol_ant_order[11] << " " << 
		sol_ant_order[12] << " " << 
		sol_ant_order[13] << " " << 
		sol_ant_order[14] << " " << 
		sol_ant_order[15] << " " << 
		sol_ant_order[16] << " " << 
		sol_ant_order[17] << " " << 
		sol_ant_order[18] << " " << 
		sol_ant_order[19] << " " << 
		sol_ant_order[20] << " " << 
		sol_ant_order[21] << " " << 
		sol_ant_order[22] << " " << 
		sol_ant_order[23] << " " << 
		sol_ant_order[24] << " " << 
		sol_ant_order[25] << " " << 
		sol_ant_order[26] << " " << 
		sol_ant_order[27] << " " << 
		sol_ant_order[28] << " " << 
		sol_ant_order[29] << " " << 
		sol_ant_order[30] << " " << 
		sol_ant_order[31] << " " << 
		std::endl;

	//of << std::setprecision(3);
	of << std::fixed << std::setprecision(3) << std::showpos;

	/* XXX always sign here */
	// [+-]6.3f
	for (int ichan = 0; ichan < nchan; ichan++) {
		of << 
			gt.at(sol_ant_order[0])[ichan] << " " << 
			gt.at(sol_ant_order[1])[ichan] << " " <<
			gt.at(sol_ant_order[2])[ichan] <<" " <<
			gt.at(sol_ant_order[3])[ichan] <<" " <<
			gt.at(sol_ant_order[4])[ichan] <<" " <<
			gt.at(sol_ant_order[5])[ichan] <<" " <<
			gt.at(sol_ant_order[6])[ichan] <<" " <<
			gt.at(sol_ant_order[7])[ichan] <<" " <<
			gt.at(sol_ant_order[8])[ichan] <<" " <<
			gt.at(sol_ant_order[9])[ichan] <<" " <<
			gt.at(sol_ant_order[10])[ichan] <<" " <<
			gt.at(sol_ant_order[11])[ichan] <<" " <<
			gt.at(sol_ant_order[12])[ichan] <<" " <<
			gt.at(sol_ant_order[13])[ichan] <<" " <<
			gt.at(sol_ant_order[14])[ichan] <<" " <<
			gt.at(sol_ant_order[15])[ichan] <<" " <<
			gt.at(sol_ant_order[16])[ichan] <<" " <<
			gt.at(sol_ant_order[17])[ichan] <<" " <<
			gt.at(sol_ant_order[18])[ichan] <<" " <<
			gt.at(sol_ant_order[19])[ichan] <<" " <<
			gt.at(sol_ant_order[20])[ichan] <<" " <<
			gt.at(sol_ant_order[21])[ichan] <<" " <<
			gt.at(sol_ant_order[22])[ichan] <<" " <<
			gt.at(sol_ant_order[23])[ichan] <<" " <<
			gt.at(sol_ant_order[24])[ichan] <<" " <<
			gt.at(sol_ant_order[25])[ichan] <<" " <<
			gt.at(sol_ant_order[26])[ichan] <<" " <<
			gt.at(sol_ant_order[27])[ichan] <<" " <<
			gt.at(sol_ant_order[28])[ichan] <<" " <<
			gt.at(sol_ant_order[29])[ichan] <<" " <<
			gt.at(sol_ant_order[30])[ichan] <<" " <<
			gt.at(sol_ant_order[31])[ichan] <<" " <<
			std::endl;
	}

	return 0;
}

selftable::selftable_t selftable::prepare_selftables ( int nchan ) {
	selftable_t ret;

	for (const antname_t& iant : sol_ant_order) {
		ret[iant] = self_t ( nchan );
	}

	return ret;
}

selftable::selftable_t selftable::read_selftables ( const std::string& infile ) {
	/* antname : vector<complex>  */
	selftable_t               ret;
	/* to store the order of the gains */
	std::vector<antname_t>    antnames;

	std::string               line;
	std::ifstream             ifs ( infile );
	/* tokens from every line */
	self_t                    gtoks;

	/* read header, populate antenna names  */
	std::getline ( ifs, line );
	{
		std::string tok;
		std::stringstream ss ( line );
		while ( ss >> tok ) {
			//const antname_t _ant { tok[0], tok[1], tok[2], ' ' };
			antname_t _ant;
			tok.copy ( _ant.data(), 4 );
			antnames.push_back ( _ant );
		}
	}
	const int nant   = antnames.size();

	/* read line by line, populate channel gains */
	while ( std::getline (ifs, line) ) {
		/* clear previous reads */
		gtoks.clear ();

		{
			std::string tok;
			std::stringstream ss ( line );
			while ( ss >> tok ) {
				/* parse the real number */
					const real_type _real ( std::stof ( tok ) );

					gtoks.push_back ( _real );

			} 
			
		}
		// sanity check
		if ( gtoks.size() != nant ) {
			throw std::runtime_error ("self file not consistent.");
		}

		// load into table
		for ( int iant = 0; iant < nant; iant++ ) {
			ret[ antnames[ iant ] ].push_back ( gtoks[iant] );
		}
	}

	/* return */
	return ret;
}
