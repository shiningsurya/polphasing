#include "gaintable.hpp"

/* overload operator<< for my formats */
std::ostream& operator<< (std::ostream& os, const gaintable::antname_t& a ) {
	os << a[0] << a[1] << a[2] << a[3];
	return os;
}
std::ostream& operator<< (std::ostream& os, const gaintable::gain_type& g ) {
	/* these flags are set before priting */
	//os << std::fixed << std::setprecision(3) << std::showpos;
	os << g.real() << g.imag() << 'j';
	return os;
}


int gaintable::write_complex_solutions (const gaintable_t& gt, const std::string& outfile) {
	std::ofstream of ( outfile );

	/* assume nchan is the same */
	int nchan = gt.at(sol_ant_order[0]).size();

/*
	fmt::print ( of, "{: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14} {: ^14}\n",
			sol_ant_order[0],
			sol_ant_order[1],
			sol_ant_order[2],
			sol_ant_order[3],
			sol_ant_order[4],
			sol_ant_order[5],
			sol_ant_order[6],
			sol_ant_order[7],
			sol_ant_order[8],
			sol_ant_order[9],
			sol_ant_order[10],
			sol_ant_order[11],
			sol_ant_order[12],
			sol_ant_order[13],
			sol_ant_order[14],
			sol_ant_order[15],
			sol_ant_order[16],
			sol_ant_order[17],
			sol_ant_order[18],
			sol_ant_order[19],
			sol_ant_order[20],
			sol_ant_order[21],
			sol_ant_order[22],
			sol_ant_order[23],
			sol_ant_order[24],
			sol_ant_order[25],
			sol_ant_order[26],
			sol_ant_order[27],
			sol_ant_order[28],
			sol_ant_order[29],
			sol_ant_order[30],
			sol_ant_order[31] 
	);
*/

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
	// [+-]6.3f[+-j]6.3f 
	for (int ichan = 0; ichan < nchan; ichan++) {
		/*
		fmt::print (of,
			"{:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f} {:+6.3f}\n",
			gt.at(sol_ant_order[0])[ichan],
			gt.at(sol_ant_order[1])[ichan],
			gt.at(sol_ant_order[2])[ichan],
			gt.at(sol_ant_order[3])[ichan],
			gt.at(sol_ant_order[4])[ichan],
			gt.at(sol_ant_order[5])[ichan],
			gt.at(sol_ant_order[6])[ichan],
			gt.at(sol_ant_order[7])[ichan],
			gt.at(sol_ant_order[8])[ichan],
			gt.at(sol_ant_order[9])[ichan],
			gt.at(sol_ant_order[10])[ichan],
			gt.at(sol_ant_order[11])[ichan],
			gt.at(sol_ant_order[12])[ichan],
			gt.at(sol_ant_order[13])[ichan],
			gt.at(sol_ant_order[14])[ichan],
			gt.at(sol_ant_order[15])[ichan],
			gt.at(sol_ant_order[16])[ichan],
			gt.at(sol_ant_order[17])[ichan],
			gt.at(sol_ant_order[18])[ichan],
			gt.at(sol_ant_order[19])[ichan],
			gt.at(sol_ant_order[20])[ichan],
			gt.at(sol_ant_order[21])[ichan],
			gt.at(sol_ant_order[22])[ichan],
			gt.at(sol_ant_order[23])[ichan],
			gt.at(sol_ant_order[24])[ichan],
			gt.at(sol_ant_order[25])[ichan],
			gt.at(sol_ant_order[26])[ichan],
			gt.at(sol_ant_order[27])[ichan],
			gt.at(sol_ant_order[28])[ichan],
			gt.at(sol_ant_order[29])[ichan],
			gt.at(sol_ant_order[30])[ichan],
			gt.at(sol_ant_order[31])[ichan]
		);
		*/

		// this will unfortunately print as (real,imag)
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

gaintable::gaintable_t gaintable::prepare_gaintables ( int nchan ) {
	gaintable_t ret;

	for (const antname_t& iant : sol_ant_order) {
		ret[iant] = gains_t ( nchan );
	}

	return ret;
}
