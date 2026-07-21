#include "models.hpp"

#include <iostream>
#include <sstream>

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
