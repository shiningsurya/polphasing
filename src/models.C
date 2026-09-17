#include "models.hpp"

#include <iostream>
#include <sstream>

/*
 *
 * Because i am defining the same complex_type in multiple namespaces
 * and since the operator<< for the complex_type in global namespace
 * this models.C should appear after gaintable.C while compilation
*/

models::model_data_t models::read_model_file ( const std::string& filepath ) {

	/* parse the filename  */
	std::string filename;
	auto start_filename = filepath.find_last_of ( "/" );
	if ( start_filename == std::string::npos ) {
		/* filepath is just filename */
		filename = filepath;
	}
	else {
		filename = filepath.substr ( start_filename + 1 );
	}

	/* tokenize the filename */
	std::vector<std::string> toks;
	{
		std::string tok;
		std::stringstream ss (filename);
		while ( std::getline ( ss, tok, '_' ) ) toks.push_back ( tok );
	}

	/* there should be 5 tokens */
	// {CAL}_{NCHAN}_{FBW}_{FEDGE}_model.txt
	if ( toks.size() != 5 ) throw std::runtime_error ("Model filename is not according to convention.");
	/* parse the tokens */
	std::string name ( toks[0] );
	int         nchan ( std::stoi ( toks[1] ) );
	real_type   fbw   ( std::stof ( toks[2] ) );
	real_type   fedge ( std::stof ( toks[3] ) );

	//fmt::print ("{}\n{} nchan={:d} fbw={:.3f} fedge={:.3f}\n", filename, name, nchan, fbw*1E-3, fedge*1E-6);

	/* prepare the model */
	model_data_t model ( nchan, fbw, fedge );

	/* read the file */
	std::string   line;
	std::vector<std::string>  ltoks;
	std::ifstream ifs (filepath);

	static const std::string s_freqs ("freqs");
	static const std::string s_i ("stokes_i");
	static const std::string s_q ("stokes_q");
	static const std::string s_u ("stokes_u");

	/* read header */
	ltoks.clear();
	std::getline ( ifs, line );
	{
		std::string tok;
		std::stringstream ss ( line );
		while ( ss >> tok ) ltoks.push_back ( tok );
	}
	
	if ( ltoks.size() != 4 ) throw std::runtime_error("Found != 4 entries in header of model file.");
	if ( ltoks[0] != s_freqs || ltoks[1] != s_i || ltoks[2] != s_q || ltoks[3] != s_u ) {
		throw std::runtime_error ("Header not matching.");
	}

	for (int ichan = 0; ichan < nchan; ichan++) {
		/* read line by line */
		ltoks.clear();
		std::getline ( ifs, line );
		{
			std::string tok;
			std::stringstream ss ( line );
			while ( ss >> tok ) ltoks.push_back ( tok );
		}
		if ( ltoks.size() != 4 ) throw std::runtime_error("Found != 4 entries in the model file.");

		/* save */
		model.freqs[ichan]    = std::stof ( ltoks[0] );
		model.stokes_i[ichan] = std::stof ( ltoks[1] );
		model.stokes_q[ichan] = std::stof ( ltoks[2] );
		model.stokes_u[ichan] = std::stof ( ltoks[3] );
	}

	return model;
}


int models::write_model_file ( const model_data_t& model, const std::string& outfile ) {
	std::ofstream of ( outfile );
	/*freqs I Q U*/

	int nchan = model.freqs.size();

	//fmt::print ( "{: <9} {: <9} {: <9} {: <9}\n", "freqs", "stokes_i", "stokes_q", "stokes_u" );
	of << "freqs" << " " << "stokes_i" << " " << "stokes_q" << " " << "stokes_u" << std::endl;

	std::cout << std::setprecision(3);

	for ( int ichan = 0; ichan < nchan; ichan++ ) {
		/*
		fmt::print ( "{: <9.6f} {: <+9.4f} {: <+9.4f} {: <+9.4f}\n", 
				model.freqs[ichan], model.stokes_i[ichan], 
				model.stokes_q[ichan], model.stokes_u[ichan]
		);
		*/
		of << model.freqs[ichan] << " " << model.stokes_i[ichan] << " " << model.stokes_q[ichan] << " " << model.stokes_u[ichan] << std::endl;
	}

	return 0;
}

int models::write_model_file ( const model_vis_t& model, const std::string& outfile ) {
	std::ofstream of ( outfile );
	/*freqs I Q U*/

	int nchan = model.nchans;

	//fmt::print ( "{: <9} {: <9} {: <9} {: <9}\n", "freqs", "stokes_i", "stokes_q", "stokes_u" );
	of << "rr" << " " << "rl" << " " << "lr" << " " << "ll" << std::endl;

	of << std::fixed << std::setprecision(3) << std::showpos;

	for ( int ichan = 0; ichan < nchan; ichan++ ) {
		of << model.rr[ichan] << " " << model.rl[ichan] << " " << model.lr[ichan] << " " << model.ll[ichan] << std::endl;
	}

	return 0;
}


models::model_type models::parallactic_correction ( 
		const real_type pa1, const real_type pa2, 
		const complex_type mrr, const complex_type mrl, const complex_type mlr, const complex_type mll
	) {

	/*
	 * The following is parallactic angle correction in linear basis.
	 * We need to do in circular basis.
	// rr
	ret[0] = mrr*c1*c2 - mll*s1*s2 + mrl*s2*c1 - mlr*s1*c2;

	// rl
	ret[1] = -mll*s1*c2 + mlr*s1*s2 + mrl*c1*c2 - mrr*s2*c1;

	// lr
	ret[2] = mll*s2*c1 + mlr*c1*c2 + mrl*s1*s2 + mrr*s1*c2;

	// ll
	ret[3] = mll*c1*c2 - mlr*s2*c1 + mrl*s1*c2 - mrr*s1*s2;
	*/

	// return this
	model_type  ret;

	// precompute the angles
	// The order matters.
	const real_type cd ( std::cos(pa1 - pa2) );
	const real_type cs ( std::cos(pa1 + pa2) );
	const real_type sd ( std::sin(pa1 - pa2) );
	const real_type ss ( std::sin(pa1 + pa2) );

	// rr
	ret[0]  = mrr * complex_type(cd,-sd);
	// rl
	ret[1]  = mrl * complex_type(cs,-ss);
	// lr
	ret[2]  = mlr * complex_type(cs, ss);
	// ll
	ret[3]  = mll * complex_type(cd, sd);

	return ret;
}

models::model_type models::parallactic_leakage_correction ( 
		// par angles
		const real_type parp, const real_type parq, 
		// leakage terms
		const complex_type lpr, const complex_type lpl, 
		const complex_type lqr, const complex_type lql, 
		// input model
		const model_type input) {

	using std::conj;

	// par complex numbers
	const complex_type zp ( std::cos(parp), std::sin(parp) );
	const complex_type zq ( std::cos(parq), std::sin(parq) );

	// unpack input
	const complex_type mrr ( input[0] );
	const complex_type mrl ( input[1] );
	const complex_type mlr ( input[2] );
	const complex_type mll ( input[3] );

	const complex_type orr = lpr*mll*zp*conj(lqr)*conj(zq) + lpr*mlr*zp*zq + mrl*conj(lqr)*conj(zp)*conj(zq) + mrr*zq*conj(zp);

	const complex_type orl = lpr*mll*zp*conj(zq) + lpr*mlr*zp*zq*conj(lql) + mrl*conj(zp)*conj(zq) + mrr*zq*conj(lql)*conj(zp);

	const complex_type olr = lpl*mrl*conj(lqr)*conj(zp)*conj(zq) + lpl*mrr*zq*conj(zp) + mll*zp*conj(lqr)*conj(zq) + mlr*zp*zq;

	const complex_type oll = lpl*mrl*conj(zp)*conj(zq) + lpl*mrr*zq*conj(lql)*conj(zp) + mll*zp*conj(zq) + mlr*zp*zq*conj(lql);

	return model_type { orr, orl, olr, oll };

}
